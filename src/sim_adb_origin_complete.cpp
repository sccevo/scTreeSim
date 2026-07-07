// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List sim_adb_origin_loop_cpp(double origin_time,
                             NumericVector a,
                             NumericVector b,
                             NumericVector d,
                             int origin_type = 0) {
  
  int max_nodes = 1024;
  std::vector<int>    v_id(max_nodes), v_type(max_nodes),
  v_parent(max_nodes, NA_INTEGER),
  v_left(max_nodes, NA_INTEGER),
  v_right(max_nodes, NA_INTEGER),
  v_status(max_nodes);
  std::vector<double> v_height(max_nodes);
  std::vector<double> v_edge_lengths; // collected during simulation
  std::vector<int>    v_edges_from, v_edges_to;
  
  // root node
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
    
    if (R::runif(0, 1) < d[v_type[idx]]) {
      v_status[idx] = 0;  // dies
    } else {
      v_status[idx] = 2;  // divides
      
      // grow storage if needed
      if (n_nodes + 2 > max_nodes) {
        max_nodes *= 2;
        v_id.resize(max_nodes);      v_type.resize(max_nodes);
        v_parent.resize(max_nodes, NA_INTEGER);
        v_left.resize(max_nodes, NA_INTEGER);
        v_right.resize(max_nodes, NA_INTEGER);
        v_status.resize(max_nodes);  v_height.resize(max_nodes);
      }
      
      int left_id  = event_counter + 1;
      int right_id = event_counter + 2;
      event_counter += 2;
      
      // single-type: children inherit origin_type
      int lt = origin_type, rt = origin_type;
      
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
      if (!left_censored) event_stack.push_back(li);
      
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
      if (!right_censored) event_stack.push_back(ri);
      
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