#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <vector>
using namespace Rcpp;

// partial convolution using fft
// used in iterative alg
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
  return real(z.head(n)) * eps;   // take real part & scale 
}

// Extinction prob
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
  // Initialize
  arma::mat X00(ntime, ntype); 
  arma::cx_mat F_t(2*ntime, ntype);
  for (int i = 0; i < ntype; i++){
    // get pdf/cdf for each type
    NumericVector pdf_nv = Rcpp::dgamma(t_nv, b(i), a(i));
    NumericVector cdf_nv = Rcpp::pgamma(t_nv, b(i), a(i));
    arma::vec pdf = as<arma::vec>(pdf_nv);
    arma::vec cdf = as<arma::vec>(cdf_nv);
    // Initialise LHS for each type
    X00.col(i) = (1.0 - rho) * (1.0 - cdf) + //lives to present, unsampled
      d(i) * cdf; //death event
    arma::vec f_t = arma::join_cols(pdf, arma::zeros(pdf.n_elem));
    F_t.col(i) = fft(f_t);    // FFT of f_t for each type
  }
  
  arma::mat X0 = X00; //previous iterate
  arma::mat Xi = X00; //current iterate
  double err = 1;
  int it = 0;
  
  // solve integral equation:
  // prob of branching but no descendants in sample
  while (err > tol && it < maxit) {
    // iterate for each type
    for (int j=0; j < ntype; j++) {
      arma::vec sumX(ntime);
      // sumX: 
      for (int k=0; k < ntype; k++){
        sumX += Xi_a(j,k) * (X0.col(j) % X0.col(k))+
          Xi_s(j,k) * (X0.col(k) % X0.col(k));
      }
      arma::vec I = convolve_fft(F_t.col(j), sumX, ntime, dx);
      Xi.col(j) = X00.col(j) + (1.-d(j))*I; 
    }
    err = arma::norm(X0-Xi, 2); // euclidean distance
    X0 = Xi;
    it++;
  }
  if (it == maxit) {
    Rcpp::warning("max iterations reached with error: %f", err);
  }
  
  return X0;
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
                             int m = 500,
                             int maxit = 100,
                             double tol = 1e-6) {
  
  // Pre-compute P0 over [0, origin_time]
  double dx = origin_time / (m - 1.0);
  arma::vec t_seq = arma::linspace(0.0, origin_time, m);
  arma::mat P0 = get_X(rho, a, b, d, Xi_a, Xi_s, t_seq, dx, maxit, tol);
  // P0 is (m x ntype); P0(i, type) = extinction prob at time t_seq(i)
  
  // Helper: look up P0 at a given absolute time for a given type
  // t here is absolute time from origin (i.e. height in the tree)
  auto lookup_p0 = [&](double t, int type) -> double {
    // t_seq runs 0..origin_time; find nearest index
    int idx = (int)std::round(t / dx);
    idx = std::max(0, std::min(idx, m - 1));
    return P0(idx, type);
  };
  
  // Simulation loop 
  int max_nodes = 1 << 20; // adjust as needed
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
    
    if (R::runif(0, 1) < d[v_type[idx]]) {
      v_status[idx] = 0;
    } else {
      v_status[idx] = 2;
      
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
      int lt = origin_type, rt = origin_type; // single-type
      
      // --- left child ---
      double left_lifetime = R::rgamma(b[lt], a[lt]);
      double left_height   = v_height[idx] - left_lifetime;
      bool   left_censored = left_height < 0;
      if (left_censored) { left_lifetime = v_height[idx]; left_height = 0.0; }
      
      int li = n_nodes++;
      v_id[li]=left_id; v_type[li]=lt; v_height[li]=left_height;
      v_parent[li]=v_id[idx]; v_left[li]=NA_INTEGER; v_right[li]=NA_INTEGER; v_status[li]=1;
      v_edges_from.push_back(v_id[idx]); v_edges_to.push_back(left_id);
      v_edge_lengths.push_back(left_lifetime);
      
      // if not censored, flip coin weighted by P0 to decide if lineage goes extinct
      if (!left_censored) {
        double p0_left = lookup_p0(left_height, lt);
        if (R::runif(0, 1) < p0_left) {
          v_status[li] = 0; // mark as extinct, don't enqueue
        } else {
          event_stack.push_back(li);
        }
      }
      
      // --- right child ---
      double right_lifetime = R::rgamma(b[rt], a[rt]);
      double right_height   = v_height[idx] - right_lifetime;
      bool   right_censored = right_height < 0;
      if (right_censored) { right_lifetime = v_height[idx]; right_height = 0.0; }
      
      int ri = n_nodes++;
      v_id[ri]=right_id; v_type[ri]=rt; v_height[ri]=right_height;
      v_parent[ri]=v_id[idx]; v_left[ri]=NA_INTEGER; v_right[ri]=NA_INTEGER; v_status[ri]=1;
      v_edges_from.push_back(v_id[idx]); v_edges_to.push_back(right_id);
      v_edge_lengths.push_back(right_lifetime);
      
      if (!right_censored) {
        double p0_right = lookup_p0(right_height, rt);
        if (R::runif(0, 1) < p0_right) {
          v_status[ri] = 0;
        } else {
          event_stack.push_back(ri);
        }
      }
      
      v_left[idx]  = left_id;
      v_right[idx] = right_id;
    }
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
    _["root_edge"]   = root_lifetime
  );
}