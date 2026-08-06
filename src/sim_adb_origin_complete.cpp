#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <vector>
using namespace Rcpp;

// partial convolution using fft
arma::vec convolve_fft(
    const arma::cx_vec& fx, 
    const arma::vec& y, 
    int n, 
    double eps) {
  int nz = 2 * y.n_elem;  
  arma::vec y_ext = y;
  y_ext.resize(nz);  
  arma::cx_vec fy = fft(y_ext);   
  arma::cx_vec fz = fx % fy;      
  arma::cx_vec z = ifft(fz); 
  return real(z.head(n)) * eps;
}

// Extinction prob - corrected version
// [[Rcpp::export]]
arma::mat get_X(
    double rho, 
    NumericVector a, 
    NumericVector b, 
    NumericVector d,
    NumericMatrix Xi_a,
    NumericMatrix Xi_s,
    arma::vec t, 
    double dx, 
    int maxit,
    double tol) {
  int ntype = a.size();
  int ntime = t.size();
  NumericVector t_nv = wrap(t);
  
  arma::mat X00(ntime, ntype); 
  arma::cx_mat F_t(2*ntime, ntype);
  
  for (int i = 0; i < ntype; i++){
    NumericVector pdf_nv = Rcpp::dgamma(t_nv, b(i), a(i));
    NumericVector cdf_nv = Rcpp::pgamma(t_nv, b(i), a(i));
    arma::vec pdf = as<arma::vec>(pdf_nv);
    arma::vec cdf = as<arma::vec>(cdf_nv);
    
    // P0 initial: probability of no sampled descendants
    // This includes both: 
    // 1. Lineage survives to present but is unsampled: (1-rho) * (1-cdf)
    // 2. Lineage dies before present: d * cdf
    X00.col(i) = (1.0 - rho) * (1.0 - cdf) + d(i) * cdf;
    
    arma::vec f_t = arma::join_cols(pdf, arma::zeros(pdf.n_elem));
    F_t.col(i) = fft(f_t);
  }
  
  arma::mat X0 = X00;
  arma::mat Xi = X00;
  double err = 1;
  int it = 0;
  
  while (err > tol && it < maxit) {
    for (int j=0; j < ntype; j++) {
      arma::vec sumX(ntime);
      for (int k=0; k < ntype; k++){
        sumX += Xi_a(j,k) * (X0.col(j) % X0.col(k)) +
          Xi_s(j,k) * (X0.col(k) % X0.col(k));
      }
      arma::vec I = convolve_fft(F_t.col(j), sumX, ntime, dx);
      Xi.col(j) = X00.col(j) + (1.0 - d(j)) * I; 
    }
    err = arma::norm(X0-Xi, 2);
    X0 = Xi;
    it++;
  }
  if (it == maxit) {
    Rcpp::warning("max iterations reached with error: %f", err);
  }
  
  return X0;
}

// helper: sample child types
std::pair<int,int> sample_child_types(int parent_type,
                                      const NumericMatrix& Xi_as,
                                      const NumericMatrix& Xi_s,
                                      int ntype) {
  std::vector<double> probs;
  std::vector<std::pair<int,int>> combos;
  
  for (int j = 0; j < ntype; j++) {
    for (int k = 0; k < ntype; k++) {
      if (Xi_as(parent_type, j) > 0 && j != k) {
        probs.push_back(Xi_as(parent_type, j));
        combos.push_back({j, k});
      }
      if (Xi_s(parent_type, k) > 0 && j == k) {
        probs.push_back(Xi_s(parent_type, k));
        combos.push_back({k, k});
      }
    }
  }
  
  double total = 0; 
  for (double p : probs) total += p;
  double u = R::runif(0, 1) * total;
  double cum = 0;
  for (int i = 0; i < (int)probs.size(); i++) {
    cum += probs[i];
    if (u <= cum) return combos[i];
  }
  return combos.back();
}

// [[Rcpp::export]]
List sim_adb_origin_loop_cpp(double origin_time,
                             NumericVector a,
                             NumericVector b,
                             NumericVector d,
                             double rho,
                             NumericMatrix Xi_a,
                             NumericMatrix Xi_s,
                             int origin_type = 0,
                             int m = 1024,
                             int maxit = 100,
                             double tol = 1e-6) {
  
  int ntype = a.size();
  
  int n_divisions = 0;
  int n_deaths    = 0;
  int n_censored  = 0;   // children born after the present (height < 0)
  
  // Pre-compute P0 over [0, origin_time]
  double dx = origin_time / (m - 1.0);
  arma::vec t_seq = arma::linspace(0.0, origin_time, m);
  arma::mat P0 = get_X(rho, a, b, d, Xi_a, Xi_s, t_seq, dx, maxit, tol);
  
  // Helper: look up P0 at a given time since origin with linear interpolation
  auto lookup_p0 = [&](double height, int type) -> double {
    if (rho >= 1.0) return 0.0;
    double pos = height / dx;
    int lo = (int)std::floor(pos);
    int hi = lo + 1;
    lo = std::max(0, std::min(lo, m - 1));
    hi = std::max(0, std::min(hi, m - 1));
    double frac = pos - std::floor(pos);
    return (1.0 - frac) * P0(lo, type) + frac * P0(hi, type);
  };
  
  // Simulation loop 
  int max_nodes = 1 << 20;
  std::vector<int>    v_id(max_nodes), v_type(max_nodes),
  v_parent(max_nodes, NA_INTEGER),
  v_left(max_nodes, NA_INTEGER),
  v_right(max_nodes, NA_INTEGER),
  v_status(max_nodes);
  std::vector<double> v_height(max_nodes);
  std::vector<double> v_edge_lengths;
  std::vector<int>    v_edges_from, v_edges_to;
  
  double root_lifetime = R::rgamma(b[origin_type], a[origin_type]);
  double root_height   = origin_time - root_lifetime;
  
  v_id[0]=1; v_type[0]=origin_type; v_height[0]=root_height;
  v_parent[0]=NA_INTEGER; v_left[0]=NA_INTEGER; v_right[0]=NA_INTEGER;
  v_status[0]=1;
  
  int event_counter = 1;
  int n_nodes = 1;
  std::vector<int> event_stack;
  event_stack.push_back(0);
  
  while (!event_stack.empty()) {
    int idx = event_stack.back();
    event_stack.pop_back();
    
    // For rho = 1, we need explicit death events
    if (rho >= 1.0) {
      // Explicit death check (original behavior)
      if (R::runif(0,1) < d[v_type[idx]]) {
        v_status[idx] = 0;  // Dies
        n_deaths++;
        continue;
      }
      // Otherwise, it divides (handled below)
    } else {
      // For rho < 1, use P0 for pruning
      double p0_current = lookup_p0(v_height[idx], v_type[idx]);
      
      if (R::runif(0,1) < p0_current) {
        v_status[idx] = 0;  // No sampled descendants
        n_deaths++;
        continue;
      }
      // Otherwise, it divides
    }
    
    // If not extinct, it must divide (since death is already accounted for in P0)
    v_status[idx] = 2;
    n_divisions++;
    
    if (n_nodes + 2 > max_nodes) {
      max_nodes *= 2;
      v_id.resize(max_nodes);    v_type.resize(max_nodes);
      v_parent.resize(max_nodes, NA_INTEGER);
      v_left.resize(max_nodes, NA_INTEGER);
      v_right.resize(max_nodes, NA_INTEGER);
      v_status.resize(max_nodes); v_height.resize(max_nodes);
    }
    
    int left_id  = event_counter + 1;
    int right_id = event_counter + 2;
    event_counter += 2;
    
    std::pair<int,int> child_types = sample_child_types(v_type[idx], Xi_a, Xi_s, ntype);
    int lt = child_types.first;
    int rt = child_types.second;
    
    // --- left child ---
    double left_lifetime = R::rgamma(b[lt], a[lt]);
    double left_height   = v_height[idx] - left_lifetime;
    bool   left_censored = left_height < 0;
    if (left_censored) { 
      left_lifetime = v_height[idx]; 
      left_height = 0.0; 
      n_censored++;
    }
    
    int li = n_nodes++;
    v_id[li]=left_id; v_type[li]=lt; v_height[li]=left_height;
    v_parent[li]=v_id[idx]; v_left[li]=NA_INTEGER; v_right[li]=NA_INTEGER; 
    v_status[li]=1;
    v_edges_from.push_back(v_id[idx]); v_edges_to.push_back(left_id);
    v_edge_lengths.push_back(left_lifetime);
    
    // Only process if not censored (i.e., born before present)
    if (!left_censored) {
      // We'll process this child when it's popped from the stack
      // Its extinction probability will be evaluated then
      event_stack.push_back(li);
    }
    // If censored, it's a tip at present time - we'll keep it as alive (status=1)
    // It will be marked as a tip in the final tree
    
    // --- right child ---
    double right_lifetime = R::rgamma(b[rt], a[rt]);
    double right_height   = v_height[idx] - right_lifetime;
    bool   right_censored = right_height < 0;
    if (right_censored) { 
      right_lifetime = v_height[idx]; 
      right_height = 0.0; 
      n_censored++;
    }
    
    int ri = n_nodes++;
    v_id[ri]=right_id; v_type[ri]=rt; v_height[ri]=right_height;
    v_parent[ri]=v_id[idx]; v_left[ri]=NA_INTEGER; v_right[ri]=NA_INTEGER; 
    v_status[ri]=1;
    v_edges_from.push_back(v_id[idx]); v_edges_to.push_back(right_id);
    v_edge_lengths.push_back(right_lifetime);
    
    if (!right_censored) {
      event_stack.push_back(ri);
    }
    
    v_left[idx]  = left_id;
    v_right[idx] = right_id;
  }
  
  int N = n_nodes;
  return List::create(
    _["id"]          = IntegerVector(v_id.begin(),          v_id.begin()+N),
    _["height"]      = NumericVector(v_height.begin(),      v_height.begin()+N),
    _["type"]        = IntegerVector(v_type.begin(),        v_type.begin()+N),
    _["parent"]      = IntegerVector(v_parent.begin(),      v_parent.begin()+N),
    _["leftchild"]   = IntegerVector(v_left.begin(),        v_left.begin()+N),
    _["rightchild"]  = IntegerVector(v_right.begin(),       v_right.begin()+N),
    _["status"]      = IntegerVector(v_status.begin(),      v_status.begin()+N),
    _["edges_from"]  = IntegerVector(v_edges_from.begin(),  v_edges_from.end()),
    _["edges_to"]    = IntegerVector(v_edges_to.begin(),    v_edges_to.end()),
    _["edge_lengths"]= NumericVector(v_edge_lengths.begin(),v_edge_lengths.end()),
    _["root_edge"]   = root_lifetime,
    _["n_divisions"] = n_divisions,
    _["n_deaths"]    = n_deaths,
    _["n_censored"]  = n_censored
  );
}