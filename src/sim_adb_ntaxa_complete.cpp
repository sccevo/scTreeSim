// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <queue>
#include <vector>
using namespace Rcpp;


// helper: sample child types
std::pair<int,int> sample_child_types_ntaxa(int parent_type,
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


struct Node {
  int    id, type, parent, leftchild, rightchild, status;
  double height;
};


// min-heap comparator by height
struct CmpHeight {
  bool operator()(const Node& a, const Node& b) { return a.height > b.height; }
};


// [[Rcpp::export]]
List sim_adb_loop_cpp(int ntaxa,
                      NumericVector a,   // scale params (length = n_types)
                      NumericVector b,   // shape params
                      NumericVector d,   // death probs
                      NumericMatrix Xi_as,
                      NumericMatrix Xi_s,
                      int origin_type = 0) {

  int ntype = a.size();

  //Simulation loop
  int max_nodes = ntaxa * 4 + 10;   // upper bound; resize if needed
  std::vector<int>    v_id(max_nodes), v_type(max_nodes),
  v_parent(max_nodes, NA_INTEGER),
  v_left(max_nodes, NA_INTEGER),
  v_right(max_nodes, NA_INTEGER),
  v_status(max_nodes);
  std::vector<double> v_height(max_nodes);

  // root
  double root_edge = R::rgamma(b[origin_type], a[origin_type]);
  v_id[0]=1; v_type[0]=origin_type; v_height[0]=root_edge;
  v_parent[0]=NA_INTEGER; v_left[0]=NA_INTEGER; v_right[0]=NA_INTEGER;
  v_status[0]=1;

  int event_counter = 1;
  int living = 1;
  int n_nodes = 1;

  std::priority_queue<Node, std::vector<Node>, CmpHeight> pq;
  Node root_node{1, origin_type, NA_INTEGER, NA_INTEGER, NA_INTEGER, 1, root_edge};
  pq.push(root_node);

  while (living < ntaxa && !pq.empty()) {
    Node ev = pq.top(); pq.pop();
    int idx = ev.id - 1;   // 0-based index

    if (R::runif(0,1) < d[ev.type]) {
      v_status[idx] = 0;
      living--;
    } else {
      v_status[idx] = 2;
      living++;

      //resize storage array if you run out of space
      if (n_nodes + 2 > max_nodes) {
        max_nodes *= 2;
        v_id.resize(max_nodes); v_type.resize(max_nodes);
        v_parent.resize(max_nodes, NA_INTEGER); v_left.resize(max_nodes, NA_INTEGER);
        v_right.resize(max_nodes, NA_INTEGER); v_status.resize(max_nodes);
        v_height.resize(max_nodes);
      }

      // assigns ids to the new children
      int left_id  = event_counter + 1;
      int right_id = event_counter + 2;
      event_counter += 2;

      int lt, rt;
      // single-type: both children same type as origin
      if (ntype == 1) {
        lt = origin_type;
        rt = origin_type;
      } else {
        // multi-type: sample child types based on transition matrices
        auto child_types = sample_child_types_ntaxa(ev.type, Xi_as, Xi_s, ntype);
        lt = child_types.first;
        rt = child_types.second;
      }

      double lh = ev.height + R::rgamma(b[lt], a[lt]);
      double rh = ev.height + R::rgamma(b[rt], a[rt]);

      // sample lifetimes and properties for each new child
      // --- left child ---
      int li = n_nodes++;
      v_id[li]=left_id; v_type[li]=lt; v_height[li]=lh;
      v_parent[li]=ev.id; v_left[li]=NA_INTEGER; v_right[li]=NA_INTEGER; v_status[li]=1;

      // right child
      int ri = n_nodes++;
      v_id[ri]=right_id; v_type[ri]=rt; v_height[ri]=rh;
      v_parent[ri]=ev.id; v_left[ri]=NA_INTEGER; v_right[ri]=NA_INTEGER; v_status[ri]=1;

      v_left[idx]  = left_id;
      v_right[idx] = right_id;

      pq.push({left_id,  lt, ev.id, NA_INTEGER, NA_INTEGER, 1, lh});
      pq.push({right_id, rt, ev.id, NA_INTEGER, NA_INTEGER, 1, rh});
    }
  }

  // trim to actual used nodes
  int N = n_nodes;
  return List::create(
    _["id"]         = IntegerVector(v_id.begin(),     v_id.begin()+N),
    _["height"]     = NumericVector(v_height.begin(), v_height.begin()+N),
    _["type"]       = IntegerVector(v_type.begin(),   v_type.begin()+N),
    _["parent"]     = IntegerVector(v_parent.begin(), v_parent.begin()+N),
    _["leftchild"]  = IntegerVector(v_left.begin(),   v_left.begin()+N),
    _["rightchild"] = IntegerVector(v_right.begin(),  v_right.begin()+N),
    _["status"]     = IntegerVector(v_status.begin(), v_status.begin()+N),
    _["root_edge"]  = root_edge
  );
}