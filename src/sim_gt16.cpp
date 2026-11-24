#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]
#include <R_ext/Random.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
using namespace Rcpp;

// Engine of sim_gt16_nt_diploid.
//
// States: GT16 genotypes 1..16, one byte each in a sites x nodes raw matrix.
// Branches are simulated by uniformization: Poisson(n_sites * lambda * d) candidate events, with
// lambda the maximum total rate, each at a random site in state s, firing with probability
// total[s] / lambda.


// [[Rcpp::export]]
RawVector sim_gt16_cpp(IntegerVector parent, IntegerVector child, NumericVector distance,
                       IntegerVector origin_state, IntegerMatrix dest, NumericMatrix cumulative) {
  // parent, child: branches in preorder (1-based, parent 0 = the origin)
  // distance: distance of the branch ending in each node
  // dest, cumulative: event table of .gt16_events()
  const int n_branches = child.size();
  const int n_sites = origin_state.size();
  const int n_nodes = distance.size();
  const int n_events = dest.ncol();

  // row-major copies of the event table
  std::vector<double> cum(16 * n_events);
  std::vector<unsigned char> next_state(16 * n_events);
  double lambda = 0.0;
  for (int s = 0; s < 16; s++) {
    for (int k = 0; k < n_events; k++) {
      cum[s * n_events + k] = cumulative(s, k);
      next_state[s * n_events + k] = (unsigned char) dest(s, k);
    }
    lambda = std::max(lambda, cumulative(s, n_events - 1));
  }

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
    if (!(d > 0) || !(lambda > 0)) continue;
    const size_t n_candidates = (size_t) R::rpois(n_sites * lambda * d);
    for (size_t c = 0; c < n_candidates; c++) {
      const int j = (int) R_unif_index(n_sites);
      const int s = st[j] - 1;
      // one draw decides whether the candidate fires and which event
      const double u = unif_rand() * lambda;
      const double* row = &cum[s * n_events];
      if (!(u < row[n_events - 1])) continue;
      const int k = std::upper_bound(row, row + n_events, u) - row;
      st[j] = next_state[s * n_events + std::min(k, n_events - 1)];
    }
  }
  return out;
}


// [[Rcpp::export]]
RawVector add_error_gt16_cpp(RawVector states, IntegerVector nodes, NumericMatrix error_matrix) {
  // observed genotypes of `nodes`, each site drawn from its true genotype's row of the error matrix
  const int n_sites = Rf_nrows(states);
  std::vector<double> cum(16 * 16);
  for (int s = 0; s < 16; s++) {
    double total = 0.0;
    for (int x = 0; x < 16; x++) cum[s * 16 + x] = (total += error_matrix(s, x));
  }

  RawVector out((size_t) n_sites * nodes.size());
  out.attr("dim") = IntegerVector::create(n_sites, nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    const unsigned char* st = RAW(states) + (size_t) (nodes[i] - 1) * n_sites;
    unsigned char* o = RAW(out) + (size_t) i * n_sites;
    for (int j = 0; j < n_sites; j++) {
      const double* row = &cum[(st[j] - 1) * 16];
      const double u = unif_rand() * row[15];
      const int x = std::upper_bound(row, row + 16, u) - row;
      o[j] = (unsigned char) (std::min(x, 15) + 1);
    }
  }
  return out;
}


// [[Rcpp::export]]
List gt16_format_cpp(RawVector states, IntegerVector columns, bool genotype) {
  // one string per column: genotypes joined by "_" ("AG_AC"), or the two allele sequences ("AA", "GC")
  static const char nucleotides[] = "ACGT";
  const int n_sites = Rf_nrows(states);
  const int n = columns.size();

  auto node_states = [&](int i) { return RAW(states) + (size_t) (columns[i] - 1) * n_sites; };

  if (genotype) {
    CharacterVector out(n);
    std::string buf(n_sites > 0 ? 3 * n_sites - 1 : 0, '_');
    for (int i = 0; i < n; i++) {
      const unsigned char* st = node_states(i);
      for (int j = 0; j < n_sites; j++) {
        buf[3 * j] = nucleotides[(st[j] - 1) / 4];
        buf[3 * j + 1] = nucleotides[(st[j] - 1) % 4];
      }
      SET_STRING_ELT(out, i, Rf_mkCharLen(buf.data(), buf.size()));
    }
    return List::create(out);
  }

  CharacterVector allele_1(n), allele_2(n);
  std::string a1(n_sites, 'A'), a2(n_sites, 'A');
  for (int i = 0; i < n; i++) {
    const unsigned char* st = node_states(i);
    for (int j = 0; j < n_sites; j++) {
      a1[j] = nucleotides[(st[j] - 1) / 4];
      a2[j] = nucleotides[(st[j] - 1) % 4];
    }
    SET_STRING_ELT(allele_1, i, Rf_mkCharLen(a1.data(), a1.size()));
    SET_STRING_ELT(allele_2, i, Rf_mkCharLen(a2.data(), a2.size()));
  }
  return List::create(allele_1, allele_2);
}
