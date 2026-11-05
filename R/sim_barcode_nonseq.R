# This code is adapted from Sophie Seidel's TideTree implementation and Antoine Zwaans's implementation with dropout.
# The original authors and their respective licenses are retained where applicable.
# source: https://github.com/azwaans/tidetree-dropout/blob/main/src/tidetree/simulation/SimulatedAlignment.java

#' Simulator of lineage barcodes with independent targets
#'
#' Evolves sequences according to the edit-and-silencing model from TiDeTree:
#' Seidel, S., Stadler, T. (2022) TiDeTree: a Bayesian phylogenetic framework to estimate single-cell trees and population dynamic parameters from genetic lineage tracing data. Proc. R. Soc. B, 289 (1986): 20221844. https://doi.org/10.1098/rspb.2022.1844
#'
#' Each target site evolves as a continuous-time Markov chain with three fates: 
#' staying unedited, transitioning to one of E edit outcomes, or being silenced.
#' Silencing acts at a constant rate at all times; only editing is restricted to the editing window.
#' At the end of the process, targets may drop out.
#'
#' @param tree a treedata object of a sampled tree: ultrametric, i.e. all tips
#'   are sampled cells at present (e.g. from [sim_adb_origin_samp()])
#' @param n_targets number of independent target sites
#' @param edit_rate scalar, or length-n_targets vector giving a per-target rate of any edit per time unit
#' @param edit_probs numeric vector of length E, summing to 1: the relative frequencies of the edit outcomes 1,...,E
#' @param edit_height start of the editing window, as a time before the present;
#'   defaults to the origin of the tree (the root plus the root edge, if any)
#' @param edit_duration length of the editing window; targets can be edited
#'   between \code{edit_height} and \code{edit_height - edit_duration} before the present;
#'   defaults to \code{edit_height}, i.e. editing from the start of the window until the present
#' @param silencing_rate scalar, or length-n_targets vector giving a per-target silencing rate;
#'   0 (default) = no silencing
#' @param dropout_prob scalar, or length-n_targets vector giving a per-target dropout probability; 0 (default) = no dropout
#' @param missing_state value marking silenced/dropout targets,
#'   defaults to theninteger `E + 1` (`E = length(edit_probs)`). Either a whole number or a string
#'   (e.g. "-"), different from the unedited state 0 and the edit outcomes 1,...,E.
#'   The target columns are integer if `missing_state` is numeric, and character otherwise.
#'
#' @return a data frame with one row per node in the tree and one column
#'   per target: \code{node}, \code{site_1}, ..., \code{site_<n_targets>}, values
#'   0 = unedited, 1..E = edit outcome, `missing_state` = silenced/dropout
#' @family barcode simulators
#' @export
sim_barcode_nonseq <- function(tree, n_targets, edit_rate, edit_probs,
                               edit_height = NULL, edit_duration = NULL,
                               silencing_rate = 0, dropout_prob = 0,
                               missing_state = length(edit_probs) + 1L) {

  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }
  times <- .barcode_tree_times(tree)

  # default editing window: from the origin (or the root) until the present
  if (is.null(edit_height)) edit_height <- times$origin
  if (is.null(edit_duration)) edit_duration <- edit_height
  stopifnot(
    n_targets >= 1,
    length(edit_rate) %in% c(1, n_targets), all(edit_rate >= 0),
    length(silencing_rate) %in% c(1, n_targets), all(silencing_rate >= 0),
    length(dropout_prob) %in% c(1, n_targets), all(dropout_prob >= 0 & dropout_prob <= 1),
    all(edit_probs >= 0),
    length(edit_height) == 1, length(edit_duration) == 1, edit_duration >= 0
  )
  if (!isTRUE(all.equal(sum(edit_probs), 1))) {
    stop("`edit_probs` must sum to 1.", call. = FALSE)
  }
  if (length(missing_state) != 1 || is.na(missing_state) ||
      as.character(missing_state) %in% as.character(0:length(edit_probs))) {
    stop("`missing_state` must be a single value different from 0 and the edit outcomes 1,...,E.", call. = FALSE)
  }
  n_targets <- as.integer(n_targets)
  E <- length(edit_probs)

  # the generic simulator with one barcode of a single site per target, edited inside the window
  # (as time since the origin) and silenced at a constant rate
  window_start <- times$origin - edit_height
  edit <- barcode_rate(edit_rate, time = rate_window(window_start, window_start + edit_duration))
  codes <- .sim_barcode_codes(times, n_targets, 1L, .as_barcode_rate(edit, n_targets, "edit_rate"),
                              .as_barcode_rate(silencing_rate, n_targets, "silencing_rate"),
                              .as_barcode_rate(dropout_prob, n_targets, "dropout_prob"),
                              edit_probs, sequential = FALSE)

  # nodes x targets; the silenced code E + 1 is relabelled to missing_state
  state_mat <- codes[, 1, , drop = TRUE]
  dim(state_mat) <- dim(codes)[c(1, 3)]
  if (is.character(missing_state)) {
    state_mat[] <- as.character(state_mat)
    state_mat[state_mat == as.character(E + 1L)] <- missing_state
  } else {
    state_mat[state_mat == E + 1L] <- as.integer(missing_state)
  }

  out <- data.frame(node = seq_len(nrow(state_mat)))
  for (site in seq_len(n_targets)) out[[paste0("site_", site)]] <- state_mat[, site]
  out
}
