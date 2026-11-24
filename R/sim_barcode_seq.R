# This code is adapted from Sophie Seidel's and Antoine Zwaans's SciPhy implementation with missing data.
# The original authors and their respective licenses are retained where applicable.
# source: https://github.com/azwaans/SciPhy/blob/missing_data/src/sciphy/evolution/simulation/SimulatedSciPhyAlignmentHeritableMissingBcodes.java

#' Simulator of sequentially-edited lineage barcodes
#'
#' Based on the SciPhy model:
#' Seidel, S., Zwaans, A., Regalado, S. et al. (2026) SciPhy: A Bayesian phylogenetic framework using sequential genetic lineage tracing data. Nat Commun, 17, 7398. https://doi.org/10.1038/s41467-026-73377-6
#'
#' Simulates barcodes under a sequential (left-to-right) editing process
#' (unedited sites are filled in position order at rate `edit_rate`,
#' as a Poisson process, absorbing once all `n_sites` sites are edited), 
#' with an additional whole-barcode silencing state competing against editing.
#'
#' @param tree a treedata object of a sampled tree: ultrametric, i.e. all tips
#'   are sampled cells at present (e.g. from [sim_adb_origin_samp()])
#' @param n_barcodes number of barcodes per cell
#' @param n_sites number of target sites per barcode
#' @param edit_rate scalar or length-`n_barcodes` editing rate(s) per barcode
#' @param edit_probs numeric vector of length E, summing to 1: the relative
#'   frequencies of the edit outcomes 1,...,E
#' @param silencing_rate scalar or length-`n_barcodes` rate(s) at which an entire
#'   barcode is silenced, racing against the editing process
#' @param dropout_prob scalar or length-`n_barcodes` probability that a
#'   barcode drops out at a tip
#' @param missing_state value marking silenced/dropout barcodes, 
#'   defaults to the integer `E + 1` (`E = length(edit_probs)`).
#'   Either a whole number or a string (e.g. "-"), different from the unedited
#'   state 0 and the edit outcomes 1, ..., E. 
#' @param write_as_string if `TRUE` (default), each barcode is returned as one
#'   "_"-separated string per node; if `FALSE`, as a node x site data frame
#'
#' @return if `write_as_string = TRUE`, a data frame with one row per node in the tree and one column
#'   per barcode: \code{node}, \code{barcode_1}, ..., \code{barcode_<n_barcodes>},
#'   each entry a string of `n_sites` positions separated by "_" (e.g. "0_2_1_0_0"), 
#'   where 0 = unedited, 1...E = edit outcome and `missing_state` = silenced/dropout.
#'   If `write_as_string = FALSE`, a named list (\code{barcode_1}, ...) with one data frame per barcode, 
#'   with one row per node and columns \code{node}, \code{site_1}, ..., \code{site_<n_sites>}.
#' @family barcode simulators
#' @export
sim_barcode_seq <- function(tree, n_barcodes, n_sites, edit_rate, edit_probs,
                            silencing_rate = 0, dropout_prob = 0, missing_state = length(edit_probs) + 1L, write_as_string = TRUE) {

  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }

  stopifnot(
    n_barcodes >= 1, n_sites >= 1,
    length(edit_rate) %in% c(1, n_barcodes), all(edit_rate >= 0),
    length(silencing_rate) %in% c(1, n_barcodes), all(silencing_rate >= 0),
    length(dropout_prob) %in% c(1, n_barcodes), all(dropout_prob >= 0 & dropout_prob <= 1),
    all(edit_probs >= 0)
  )

  # allow a single shared rate to stand in for all barcodes
  if (length(edit_rate) == 1) {
    edit_rate <- rep(edit_rate, n_barcodes)
  }
  if (length(silencing_rate) == 1) {
    silencing_rate <- rep(silencing_rate, n_barcodes)
  }
  if (length(dropout_prob) == 1) {
    dropout_prob <- rep(dropout_prob, n_barcodes)
  }
  
  if (!isTRUE(all.equal(sum(edit_probs), 1))) {
    stop("`edit_probs` must sum to 1.", call. = FALSE)
  }
  
  if (length(missing_state) != 1 || is.na(missing_state) ||
      as.character(missing_state) %in% as.character(0:length(edit_probs))) {
    stop("`missing_state` must be a single value different from 0 and the edit outcomes 1,...,E.", call. = FALSE)
  }
  numeric_states <- is.numeric(missing_state)
  E <- length(edit_probs)

  # the generic simulator with sequential editing and constant rates
  codes <- .sim_barcode_codes(.barcode_tree_times(tree), n_barcodes, n_sites, barcode_rate(edit_rate),
                              barcode_rate(silencing_rate), barcode_rate(dropout_prob), edit_probs,
                              sequential = TRUE)
  .format_barcodes(codes, c("0", seq_len(E), as.character(missing_state)), write_as_string,
                   integer_missing = if (numeric_states) missing_state)
}
