#ifndef SCTREESIM_SAMPLE_CHILD_TYPES_H
#define SCTREESIM_SAMPLE_CHILD_TYPES_H

#include <Rcpp.h>
#include <utility>

// Sample the types of the two children of a dividing particle.
// Mirrored in R by sample_types() (R/tree_utils.R).
inline std::pair<int,int> sample_child_types(int parent_type,
                                             const Rcpp::NumericMatrix& asym_trans_prob,
                                             const Rcpp::NumericMatrix& sym_trans_prob,
                                             int ntype) {
  // total probability of the row, accumulated in the same order as the loop
  // below, so that the scaled draw r always falls below the final cum_prob:
  // outcomes are only ever chosen within the support (no fallback needed)
  double total = 0.0;
  for (int i = 0; i < ntype; i++) {
    total += sym_trans_prob(parent_type, i);
    total += asym_trans_prob(parent_type, i);
  }
  double r = R::runif(0, 1) * total;
  double cum_prob = 0.0;

  for (int i = 0; i < ntype; i++) {
    // symmetric: both children type i
    cum_prob += sym_trans_prob(parent_type, i);
    if (r < cum_prob) {
      return {i, i};
    }
    // asymmetric: one child stays parent_type, other becomes type i
    cum_prob += asym_trans_prob(parent_type, i);
    if (r < cum_prob) {
      if (R::runif(0, 1) < 0.5) {
        return {parent_type, i};
      } else {
        return {i, parent_type};
      }
    }
  }
  // unreachable: r < total and the last positive-probability outcome brings
  // cum_prob to exactly total
  Rcpp::stop("Failed to sample child types; check the transition matrices.");
}

#endif
