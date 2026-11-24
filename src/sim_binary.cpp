#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstring>
using namespace Rcpp;

// Engine of sim_binary.
//
// States: 0/1, one byte each in a sites x nodes raw matrix.
// On a branch every site flips with the closed-form probability p[state] of exp(Qd), as does the
// error at the tips. Sparse flips skip ahead by geometric gaps of p_max and accept with
// p[state] / p_max, so the cost scales with the number of flips.


// flip each of the n sites of st with probability p[st[j]]
static void flip_sites(unsigned char* st, int n, const double p[2]) {
  const double p_max = std::max(p[0], p[1]);
  if (!(p_max > 0)) return;
  if (p_max > 0.25) {
    for (int j = 0; j < n; j++) {
      if (unif_rand() < p[st[j]]) st[j] = 1 - st[j];
    }
    return;
  }
  // gaps between candidate sites are Geometric(p_max): floor(Exp(1) / -log(1 - p_max))
  const double rate = -std::log1p(-p_max);
  for (double j = std::floor(exp_rand() / rate); j < n; j += 1 + std::floor(exp_rand() / rate)) {
    unsigned char& s = st[(int) j];
    if (p[s] == p_max || unif_rand() * p_max < p[s]) s = 1 - s;
  }
}


// [[Rcpp::export]]
RawVector sim_binary_cpp(IntegerVector parent, IntegerVector child, NumericVector distance,
                         IntegerVector origin_state, double lambda, int n_tips, double alpha, double beta) {
  // parent, child: branches in preorder (1-based, parent 0 = the origin)
  // distance: distance of the branch ending in each node
  // alpha, beta: false negative and false positive rates at the tips (nodes 1..n_tips)
  const int n_branches = child.size();
  const int n_sites = origin_state.size();
  const int n_nodes = distance.size();

  // P(flip | state) on a branch of distance d: the other state's equilibrium frequency times
  // 1 - exp(-r d), with r = 1 / (2 pi0 pi1) for one expected substitution per unit d
  const double pi0 = lambda / (lambda + 1), pi1 = 1 / (lambda + 1);
  const double r = 1 / (2 * pi0 * pi1);

  RawVector out((size_t) n_sites * n_nodes);
  out.attr("dim") = IntegerVector::create(n_sites, n_nodes);
  unsigned char* states = RAW(out);

  std::vector<unsigned char> origin(n_sites);
  for (int j = 0; j < n_sites; j++) origin[j] = (unsigned char) origin_state[j];

  for (int b = 0; b < n_branches; b++) {
    const int node = child[b] - 1;
    unsigned char* st = states + (size_t) node * n_sites;
    const unsigned char* pst = parent[b] == 0 ? origin.data() : states + (size_t) (parent[b] - 1) * n_sites;
    std::memcpy(st, pst, n_sites);

    const double d = distance[node];
    if (!(d > 0)) continue;
    const double change = -std::expm1(-r * d);
    const double p[2] = {pi1 * change, pi0 * change};
    flip_sites(st, n_sites, p);
  }

  // error is a measurement effect at the tips, after evolution
  if (alpha > 0 || beta > 0) {
    const double p[2] = {beta, alpha};
    for (int i = 0; i < n_tips; i++) flip_sites(states + (size_t) i * n_sites, n_sites, p);
  }
  return out;
}


// [[Rcpp::export]]
SEXP binary_format_cpp(RawVector states, IntegerVector columns, bool as_string) {
  // the selected columns of the sites x nodes `states`, as one "0"/"1" string each or as an
  // integer columns x sites matrix
  const int n_sites = Rf_nrows(states);
  const int n = columns.size();
  auto node_states = [&](int i) { return RAW(states) + (size_t) (columns[i] - 1) * n_sites; };

  if (as_string) {
    CharacterVector out(n);
    std::string buf(n_sites, '0');
    for (int i = 0; i < n; i++) {
      const unsigned char* st = node_states(i);
      for (int j = 0; j < n_sites; j++) buf[j] = '0' + st[j];
      SET_STRING_ELT(out, i, Rf_mkCharLen(buf.data(), buf.size()));
    }
    return out;
  }

  IntegerMatrix out(n, n_sites);
  for (int i = 0; i < n; i++) {
    const unsigned char* st = node_states(i);
    for (int j = 0; j < n_sites; j++) out(i, j) = st[j];
  }
  return out;
}
