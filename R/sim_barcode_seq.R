#' Simulator of sequentially-edited barcodes
#'
#' Based on the SciPhy model:
#' Seidel, S., Zwaans, A., Regalado, S. et al. SciPhy: A Bayesian phylogenetic framework using sequential genetic lineage tracing data. Nat Commun 17, 7398 (2026). https://doi.org/10.1038/s41467-026-73377-6
#'
#' R version of the code has been adapted from Nicola Mulberry's implementation
#' source: https://github.com/sccevo/scTreeSim/blob/barcode_simulator/R/typewriter_barcodes.R
#' The original authors and their respective licenses are retained where applicable.
#'
#' Simulates barcodes under a sequential/left-to-right editing process
#' (unedited sites are filled in position order at rate `edit_rate`, as a
#' Poisson process, absorbing once all `n_sites` sites are edited), with an
#' additional whole-barcode silencing state competing against editing.
#'
#' @param tree a treedata object of a sampled tree: ultrametric, i.e. all tips
#'   are sampled cells at present (e.g. from [sim_adb_origin_samp()])
#' @param n_sites number of target sites per barcode
#' @param edit_rate scalar or length-`n_barcodes` editing rate(s) per barcode
#' @param silencing_rate scalar or length-`n_barcodes` rate(s) at which an entire
#'   barcode is silenced, racing against the editing process
#' @param dropout_prob scalar or length-`n_barcodes` probability that a
#'   barcode drops out at a tip
#' @param n_barcodes number of barcodes (tapes) per cell
#' @param chars vector of unique characters/state labels to be inserted
#' @param char_probs probabilities of inserting each of `chars` (same length as
#'   `chars`); `NULL` (default) for a uniform draw
#'
#' @return a data frame with one row per node in the tree and one column
#'   per barcode: \code{node}, \code{barcode_1}, ..., \code{barcode_<n_barcodes>},
#'   each entry a string of `n_sites` positions separated by "_"
#'   (e.g. "0_A_B_0_0"), where "0" = unedited and "-" = silenced/dropout
#'   (the silenced state; every site reads "-" once a barcode is silenced
#'   or dropped out)
#' @export
sim_barcode_seq <- function(tree, n_sites, edit_rate, silencing_rate = 0,
                            dropout_prob = 0, n_barcodes, chars, char_probs = NULL) {

  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }

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
  stopifnot(length(edit_rate) == n_barcodes)
  stopifnot(length(silencing_rate) == n_barcodes)
  stopifnot(length(dropout_prob) == n_barcodes)
  if (anyNA(dropout_prob) || any(dropout_prob < 0 | dropout_prob > 1)) {
    stop("`dropout_prob` must contain probabilities in [0, 1].", call. = FALSE)
  }

  # NULL: uniform draw over insertion characters
  if (!is.null(char_probs) &&
      (length(char_probs) != length(chars) || anyNA(char_probs) || any(char_probs < 0) ||
       !isTRUE(all.equal(sum(char_probs), 1)))) {
    stop("`char_probs` must be NULL or non-negative probabilities summing to 1, one per element of `chars`.",
         call. = FALSE)
  }

  missing_state <- "-"

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  # heights: time-since-present (0 = present, root = age of tree)
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  heights <- max(depth_from_root) - depth_from_root

  # visit parent before child so state can propagate down the tree
  order_df <- tree_df[order(-heights[tree_df$node]), ]

  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]

  # tips are nodes 1..Ntip (all sampled cells, since the tree is ultrametric)
  tip_nodes <- seq_along(tree@phylo$tip.label)

  bc_cols <- vector("list", n_barcodes)
  for (bc in seq_len(n_barcodes)) {
    # state_at: barcode string per node; silenced_at: whether that barcode
    # has already been silenced (needed since silencing must persist down
    # every descendant branch once it happens)
    state_at <- list()
    silenced_at <- list()
    state_at[[as.character(root)]] <- rep("0", n_sites)
    silenced_at[[as.character(root)]] <- FALSE

    # evolve along the root/origin edge first, if the tree has one
    if (!is.null(root_edge) && root_edge > 0) {
      res <- .evolve_branch_seq(
        state_at[[as.character(root)]], silenced_at[[as.character(root)]],
        origin_height - heights[root], edit_rate[bc], silencing_rate[bc], chars, char_probs, missing_state
      )
      state_at[[as.character(root)]] <- res$state
      silenced_at[[as.character(root)]] <- res$silenced
    }

    # then walk every remaining branch, parent state -> child state
    for (i in seq_len(nrow(order_df))) {
      node <- order_df$node[i]
      if (node == root) next
      parent <- order_df$parent[i]
      branch_length <- heights[parent] - heights[node]

      res <- .evolve_branch_seq(
        state_at[[as.character(parent)]], silenced_at[[as.character(parent)]],
        branch_length, edit_rate[bc], silencing_rate[bc], chars, char_probs, missing_state
      )
      state_at[[as.character(node)]] <- res$state
      silenced_at[[as.character(node)]] <- res$silenced
    }

    # apply dropout after simulation,
    # skipped for barcodes already fully silenced (nothing left to mask)
    if (dropout_prob[bc] > 0) {
      for (nd in tip_nodes) {
        key <- as.character(nd)
        if (!silenced_at[[key]] && stats::runif(1) < dropout_prob[bc]) {
          state_at[[key]][] <- missing_state
        }
      }
    }

    bc_cols[[bc]] <- vapply(
      tree_df$node,
      function(nd) paste(state_at[[as.character(nd)]], collapse = "_"),
      character(1)
    )
  }

  out <- data.frame(node = tree_df$node)
  for (bc in seq_len(n_barcodes)) out[[paste0("barcode_", bc)]] <- bc_cols[[bc]]
  out
}


# Evolve one barcode copy along one branch. At each event, draw whether it
# is an edit or silencing
.evolve_branch_seq <- function(parent_state, parent_silenced, branch_length,
                               edit_rate, silencing_rate, chars, char_probs, missing_state) {
  state <- parent_state

  # already silenced upstream, or no time on this branch: nothing to do
  if (parent_silenced || branch_length <= 0) {
    return(list(state = state, silenced = parent_silenced))
  }

  t <- 0

  repeat {
    # once every site is edited, only silencing can still happen
    saturated <- all(state != "0")
    current_edit_rate <- if (saturated) 0 else edit_rate
    event_rate <- current_edit_rate + silencing_rate

    if (event_rate <= 0) break  # no editing left to do and no silencing possible

    # draw time to the next event (edit or silencing) on the combined clock
    t <- t + stats::rexp(1, event_rate)
    if (t > branch_length) break  # no more events fit in this branch

    # which of the two competing events fired, weighted by relative rate
    is_silencing <- stats::runif(1) < (silencing_rate / event_rate)
    if (is_silencing) {
      state[] <- missing_state
      return(list(state = state, silenced = TRUE))
    }

    # edit event: fill the leftmost still-unedited site
    # (index chars via sample.int: sample(chars) would draw from 1:chars for a single number)
    pos <- min(which(state == "0"))
    state[pos] <- chars[sample.int(length(chars), size = 1, prob = char_probs)]
  }

  list(state = state, silenced = FALSE)
}
