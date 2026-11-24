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

  # edges in cladewise order, so every parent is simulated before its children; the root edge
  # comes first, from the origin (parent 0)
  edges <- ape::reorder.phylo(tree@phylo, order = "cladewise")$edge
  n_tips <- ape::Ntip(tree@phylo)
  root <- n_tips + 1L
  nodes <- tibble::as_tibble(tree)$node

  # one sites x nodes raw matrix of states 0/1 per sequence (see src/sim_binary.cpp), with the
  # error at the tips (the true tip states could also be returned, as in sim_gt16_nt_diploid)
  seqs <- lapply(lambda, function(lam) {
    start <- if (is.null(origin_seq)) stats::rbinom(n_sites, 1, .binary_pi(lam)[2]) else origin_seq
    sim_binary_cpp(c(0L, edges[, 1]), c(root, edges[, 2]), d, as.integer(start), lam, n_tips, alpha, beta)
  })
  names(seqs) <- paste0("binary_", seq_len(n_sequences))

  if (write_as_string) {
    return(data.frame(node = nodes, lapply(seqs, binary_format_cpp, nodes, TRUE)))
  }
  lapply(seqs, function(states) {
    m <- binary_format_cpp(states, nodes, FALSE)
    colnames(m) <- paste0("site_", seq_len(n_sites))
    data.frame(node = nodes, m)
  })
}

# equilibrium frequencies of states 0 and 1
.binary_pi <- function(lambda) {
  c(lambda / (lambda + 1), 1 / (lambda + 1))
}
