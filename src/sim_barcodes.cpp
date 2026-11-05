#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Engine of all barcode simulators (sim_barcode_generic, sim_barcode_seq, sim_barcode_nonseq).
//
// On a branch, editing and silencing are independent Poisson processes, and silencing overwrites
// every edit, so the state at the end of the branch only depends on the expected numbers of events
// (the integrated rates), not on when they happen:
// * sequential: the barcode is silenced with probability 1 - exp(-L_s); otherwise
//   min(Poisson(L_e), number of unedited sites) sites are edited, left to right;
// * independent sites: each site is silenced with probability 1 - exp(-L_s); otherwise an unedited
//   site is edited with probability 1 - exp(-L_e).
// Profiles without a closed-form integral (rate_fn) are simulated by thinning instead.
//
// States: 0 = unedited, 1..E = edit outcome, E + 1 = silenced or dropped out.


// One rate (editing or silencing) on every branch; the branch ending in node i (0-based) is branch i.
// Closed form: scale[i] = type multiplier x integral of the time profile over the branch,
// so a base rate b gives b * scale[i] expected events.
// Thinning (g given): scale[i] = type multiplier; the branch is cut into pieces
// [start[k], end[k]) for k in ptr[i]..ptr[i + 1] - 1, on which the profile is at most max[k].
struct RateProcess {
  NumericVector scale;
  bool thinning;
  IntegerVector ptr;
  NumericVector start, end, max;
  Function g;

  explicit RateProcess(const List& p)
    : scale(as<NumericVector>(p["scale"])), thinning(!Rf_isNull(p["g"])),
      g(thinning ? as<Function>(p["g"]) : Function("identity")) {
    if (thinning) {
      ptr = as<IntegerVector>(p["ptr"]);
      start = as<NumericVector>(p["start"]);
      end = as<NumericVector>(p["end"]);
      max = as<NumericVector>(p["max"]);
    }
  }

  // number of events with base rate `base` on branch i, at most `cap`
  int count(double base, int i, int cap) {
    if (cap <= 0) return 0;
    double mult = base * scale[i];
    if (!(mult > 0)) return 0;
    if (!thinning) {
      if (cap == 1) return unif_rand() < -std::expm1(-mult);
      return (int) std::min(R::rpois(mult), (double) cap);
    }
    int n = 0;
    for (int k = ptr[i]; k < ptr[i + 1]; k++) {
      double bound = mult * max[k];
      if (!(bound > 0)) continue;
      double t = start[k];
      while (true) {
        t += exp_rand() / bound;
        if (t >= end[k]) break;
        double accept = as<double>(g(t)) / max[k];
        if (accept > 1 + 1e-8) stop("A rate profile exceeds its stated maximum at time %g.", t);
        if (accept >= 1 || unif_rand() < accept) {
          if (++n >= cap) return n;
        }
      }
    }
    return n;
  }
};


// [[Rcpp::export]]
IntegerVector sim_barcodes_cpp(IntegerVector parent, IntegerVector child, IntegerVector tips,
                               int n_barcodes, int n_sites, bool sequential,
                               NumericVector edit_base, List edit_process,
                               NumericVector silencing_base, List silencing_process,
                               NumericVector dropout_base, NumericVector dropout_mult,
                               NumericVector edit_probs) {
  // parent, child: branches in preorder (1-based node ids, parent 0 = the origin, for the root edge)
  // tips: tip nodes (1-based), dropout_mult: type multiplier of the dropout probability of each tip
  const int n_nodes = child.size();
  const int n_outcomes = edit_probs.size();
  const int missing = n_outcomes + 1;

  RateProcess edit(edit_process), silencing(silencing_process);

  std::vector<double> cum(n_outcomes);
  double total = 0.0;
  for (int k = 0; k < n_outcomes; k++) cum[k] = (total += edit_probs[k]);
  // first outcome whose cumulative probability exceeds the draw (zero-probability outcomes never are)
  auto draw_outcome = [&]() {
    int k = std::upper_bound(cum.begin(), cum.end(), unif_rand() * total) - cum.begin();
    return std::min(k, n_outcomes - 1) + 1;
  };

  // states of one barcode: node-major, sites contiguous
  std::vector<int> state((size_t) n_nodes * n_sites);
  IntegerVector out((size_t) n_nodes * n_sites * n_barcodes);
  out.attr("dim") = IntegerVector::create(n_nodes, n_sites, n_barcodes);

  for (int bc = 0; bc < n_barcodes; bc++) {
    const double e_base = edit_base[bc], s_base = silencing_base[bc];

    for (int b = 0; b < n_nodes; b++) {
      const int node = child[b] - 1;
      int* st = &state[(size_t) node * n_sites];
      if (parent[b] == 0) {
        std::fill(st, st + n_sites, 0);
      } else {
        const int* pst = &state[(size_t) (parent[b] - 1) * n_sites];
        std::copy(pst, pst + n_sites, st);
      }

      if (sequential) {
        if (st[0] == missing) continue;  // silenced upstream
        if (silencing.count(s_base, node, 1)) {
          std::fill(st, st + n_sites, missing);
          continue;
        }
        int first_free = 0;
        while (first_free < n_sites && st[first_free] != 0) first_free++;
        int n_edits = edit.count(e_base, node, n_sites - first_free);
        for (int j = 0; j < n_edits; j++) st[first_free + j] = draw_outcome();
      } else {
        for (int s = 0; s < n_sites; s++) {
          if (st[s] == missing) continue;
          if (silencing.count(s_base, node, 1)) {
            st[s] = missing;
          } else if (st[s] == 0 && edit.count(e_base, node, 1)) {
            st[s] = draw_outcome();
          }
        }
      }
    }

    // dropout at the tips: the whole barcode if sequential, otherwise each site
    for (int i = 0; i < tips.size(); i++) {
      double p_drop = dropout_base[bc] * dropout_mult[i];
      if (!(p_drop > 0)) continue;
      int* st = &state[(size_t) (tips[i] - 1) * n_sites];
      if (sequential) {
        if (st[0] != missing && unif_rand() < p_drop) std::fill(st, st + n_sites, missing);
      } else {
        for (int s = 0; s < n_sites; s++) {
          if (st[s] != missing && unif_rand() < p_drop) st[s] = missing;
        }
      }
    }

    // nodes x sites x barcodes array
    int* o = &out[(size_t) bc * n_nodes * n_sites];
    for (int node = 0; node < n_nodes; node++) {
      for (int s = 0; s < n_sites; s++) o[(size_t) s * n_nodes + node] = state[(size_t) node * n_sites + s];
    }
  }
  return out;
}
