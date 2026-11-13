#' Binary sequence evolution along a tree
#'
#' Simulates one or more binary sequences along a fixed tree under a 2-state substitution model,
#' with optional false positive/negative error at the tips. Each site evolves independently with
#' transition matrix \eqn{P(d) = \exp(Qd)}, where \eqn{d} is the evolutionary distance of a
#' branch under the clock model and \eqn{Q} is scaled to one expected substitution per site per
#' unit of distance.
#'
#' The binary model is adapted from the BEAST2 package phylonco:
#' Chen, K., Moravec, J. C., Gavryushkin, A., Welch, D., & Drummond, A. J. (2022).
#' Accounting for errors in data improves divergence time estimates in single-cell cancer evolution.
#' Molecular Biology and Evolution, 39(8), msac143.
#'
#' @param tree a treedata object
#' @param n_sites number of binary sites per sequence
#' @param n_sequences number of independent sequences per node
#' @param lambda rate of back-mutation 1 -> 0 relative to mutation 0 -> 1, i.e.
#'   \eqn{Q = [-1, 1; \lambda, -\lambda]} before scaling; length 1 or `n_sequences`
#' @param origin_seq optional vector of 0s and 1s of length `n_sites`, the sequence at the origin;
#'   `NULL` (default) = drawn from the equilibrium of each sequence
#' @param clock a clock model, e.g. [clock_strict()], drawn once and shared by all sequences;
#'   `NULL` = use the `distance` column of `tree` (see [clock_tree()])
#' @param alpha false negative rate (P(observed 0 | true 1))
#' @param beta false positive rate (P(observed 1 | true 0))
#' @param write_as_string if `TRUE` (default), each sequence is returned as one string per node;
#'   if `FALSE`, as a node x site data frame
#'
#' @return if `write_as_string = TRUE`, a data frame with one row per node and one column per
#'   sequence: `node`, `binary_1`, ..., `binary_<n_sequences>`. Each entry is a string of `n_sites`
#'   characters "0"/"1".
#'   If `write_as_string = FALSE`, a named list (`binary_1`, ...) with one data frame per sequence,
#'   with one row per node and columns `node`, `site_1`, ..., `site_<n_sites>` (integer 0/1).
#' @seealso [clock_strict()], [clock_relaxed()], [spike()] for clock models, [clock_tree()] to store
#'   and reuse a clock draw, and `vignette("simulate_binary")`
#' @export
sim_binary <- function(tree, n_sites, n_sequences = 1, lambda = 1, origin_seq = NULL,
                       clock = clock_strict(), alpha = 0, beta = 0, write_as_string = TRUE) {
  stopifnot(methods::is(tree, "treedata"))
  stopifnot(n_sites >= 1L, n_sequences >= 1L)
  stopifnot(is.numeric(lambda), length(lambda) %in% c(1, n_sequences), all(lambda > 0))
  stopifnot(is.null(clock) || inherits(clock, "clock_model"))
  stopifnot(alpha >= 0, alpha <= 1, beta >= 0, beta <= 1)
  stopifnot(is.logical(write_as_string), length(write_as_string) == 1, !is.na(write_as_string))
  lambda <- rep_len(lambda, n_sequences)

  if (!is.null(origin_seq)) {
    stopifnot(length(origin_seq) == n_sites, all(origin_seq %in% c(0L, 1L)))
  }

  # one clock draw for all sequences
  if (!is.null(clock)) tree <- clock_tree(tree, clock)
  d <- .node_distances(tree)

  # edges in cladewise order, so every parent is simulated before its children
  edges <- ape::reorder.phylo(tree@phylo, order = "cladewise")$edge
  n_tips <- ape::Ntip(tree@phylo)
  root <- n_tips + 1
  nodes <- tibble::as_tibble(tree)$node

  seqs <- vector("list", n_sequences)
  for (k in seq_len(n_sequences)) {
    lam <- lambda[k]

    # node x site matrix of states 0/1; the origin sequence evolves along the root edge
    states <- matrix(NA_integer_, length(d), n_sites)
    start <- if (is.null(origin_seq)) stats::rbinom(n_sites, 1, .binary_pi(lam)[2]) else origin_seq
    states[root, ] <- .mutate_binary_seq(as.integer(start), lam, d[root])
    for (i in seq_len(nrow(edges))) {
      states[edges[i, 2], ] <- .mutate_binary_seq(states[edges[i, 1], ], lam, d[edges[i, 2]])
    }

    # error is a measurement effect at the tips, after evolution
    # (the true tip states could also be returned here, as in sim_gt16_nt_diploid)
    if (alpha > 0 || beta > 0) {
      tips <- seq_len(n_tips)
      states[tips, ] <- .add_error(states[tips, , drop = FALSE], alpha, beta)
    }

    seqs[[k]] <- states[nodes, , drop = FALSE]
  }
  names(seqs) <- paste0("binary_", seq_len(n_sequences))

  if (write_as_string) {
    return(data.frame(node = nodes, lapply(seqs, apply, 1, paste, collapse = "")))
  }
  lapply(seqs, function(m) {
    colnames(m) <- paste0("site_", seq_len(n_sites))
    data.frame(node = nodes, m)
  })
}

# equilibrium frequencies of states 0 and 1
.binary_pi <- function(lambda) {
  c(lambda / (lambda + 1), 1 / (lambda + 1))
}

# Probability that a site in state s (0 or 1) is in the other state after distance d (vectorised
# over s and d). Closed form of exp(Qd): the other state's equilibrium frequency times
# 1 - exp(-r d), with r = q01 + q10 = 1 / (2 pi0 pi1) for one expected substitution per unit d.
.binary_flip_prob <- function(lambda, s, d) {
  pi <- .binary_pi(lambda)
  r <- 1 / (2 * pi[1] * pi[2])
  other <- 1L - s
  pi[other + 1L] * (1 - exp(-r * d))
}

# flip each site independently with its flip probability
.mutate_binary_seq <- function(seq, lambda, d) {
  flip <- stats::runif(length(seq)) < .binary_flip_prob(lambda, seq, d)
  seq[flip] <- 1L - seq[flip]
  seq
}

# Observation error on true states (vector or matrix of 0/1, each site independent):
# a true 1 reads 0 with probability alpha (false negative), a true 0 reads 1 with probability
# beta (false positive).
.add_error <- function(states, alpha, beta) {
  flip <- stats::runif(length(states)) < c(beta, alpha)[states + 1L]
  states[flip] <- 1L - states[flip]
  states
}
