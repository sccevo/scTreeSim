#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <vector>
using namespace Rcpp;


// helper: sample child types
std::pair<int,int> sample_child_types(int parent_type,
                                      const NumericMatrix& Xi_as,
                                      const NumericMatrix& Xi_s,
                                      int ntype) {
  double r = R::runif(0, 1);
  double cum_prob = 0.0;
  
  for (int i = 0; i < ntype; i++) {
    // symmetric: both children type i
    cum_prob += Xi_s(parent_type, i);
    if (r < cum_prob) {
      return {i, i};
    }
    // asymmetric: one child stays parent_type, other becomes type i
    cum_prob += Xi_as(parent_type, i);
    if (r < cum_prob) {
      if (R::runif(0, 1) < 0.5) {
        return {parent_type, i};
      } else {
        return {i, parent_type};
      }
    }
  }
  // fallback (should not reach here if matrices sum to 1)
  return {parent_type, parent_type};
}

// [[Rcpp::export]]
List sim_adb_origin_loop_cpp(double origin_time,
                             NumericVector a,
                             NumericVector b,
                             NumericVector d,
                             NumericMatrix Xi_as,
                             NumericMatrix Xi_s,
                             int origin_type = 0) {
  
  int ntype = a.size();
  
  // Simulation loop 
  int max_nodes = 1024;
  std::vector<int>    v_id(max_nodes), v_type(max_nodes),
  v_parent(max_nodes, NA_INTEGER),
  v_left(max_nodes, NA_INTEGER),
  v_right(max_nodes, NA_INTEGER),
  v_status(max_nodes);
  std::vector<double> v_height(max_nodes);
  std::vector<double> v_edge_lengths;
  std::vector<int>    v_edges_from, v_edges_to;
  
  // drawing root note parameters
  double root_lifetime = R::rgamma(b[origin_type], a[origin_type]);
  double root_height   = origin_time - root_lifetime;
  
  v_id[0]=1; v_type[0]=origin_type; v_height[0]=root_height;
  v_parent[0]=NA_INTEGER; v_left[0]=NA_INTEGER; v_right[0]=NA_INTEGER;
  v_status[0]=1;
  
  int event_counter = 1;
  int n_nodes = 1;
  
  // event stack: store indices into v_ arrays (0-based)
  std::vector<int> event_stack;
  event_stack.push_back(0);
  
  while (!event_stack.empty()) {
    int idx = event_stack.back();
    event_stack.pop_back();
    

    // Explicit death check (original behavior)
    if (R::runif(0,1) < d[v_type[idx]]) {
      v_status[idx] = 0;  // Dies
      continue;
    }
    
    // If not extinct, it must divide (since death is already accounted for in P0)
    //creates two new children
    v_status[idx] = 2;
    
    //resize storage array if you run out of space
    if (n_nodes + 2 > max_nodes) {
      max_nodes *= 2;
      v_id.resize(max_nodes);    v_type.resize(max_nodes);
      v_parent.resize(max_nodes, NA_INTEGER);
      v_left.resize(max_nodes, NA_INTEGER);
      v_right.resize(max_nodes, NA_INTEGER);
      v_status.resize(max_nodes); v_height.resize(max_nodes);
    }
    
    // assigns ids to the new children
    int left_id  = event_counter + 1;
    int right_id = event_counter + 2;
    event_counter += 2;
    
    int lt, rt;
    // single type case
    if (ntype == 1) {
      lt = origin_type;
      rt = origin_type;
    } else {
    // multi type case
      auto child_types = sample_child_types(v_type[idx], Xi_as, Xi_s, ntype);
      lt = child_types.first;
      rt = child_types.second;
    }
    
    // sample lifetimes and properties for each new child
    // --- left child ---
    double left_lifetime = R::rgamma(b[lt], a[lt]);
    double left_height   = v_height[idx] - left_lifetime;
    bool   left_censored = left_height < 0;
    if (left_censored) { 
      left_lifetime = v_height[idx]; 
      left_height = 0.0;
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
    _["root_edge"]   = root_lifetime
  );
}