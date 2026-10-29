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
  # states are handled as strings internally
  missing_state <- as.character(missing_state)

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  # heights: time-since-present (0 = present, root = age of tree)
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  heights <- max(depth_from_root) - depth_from_root

  # visit parent before child so state can propagate down the tree; use the
  # topology rather than heights, which tie on zero-length branches
  child_order <- ape::reorder.phylo(tree@phylo, order = "cladewise")$edge[, 2]
  order_df <- tree_df[match(child_order, tree_df$node), ]

  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]
  if (!is.null(tree@phylo$origin) && !isTRUE(all.equal(tree@phylo$origin, origin_height))) {
    warning("`tree@phylo$origin` differs from the tree height plus root edge; the latter is used.", call. = FALSE)
  }

  # tips are nodes 1..Ntip (all sampled cells, since the tree is ultrametric)
  tip_nodes <- seq_along(tree@phylo$tip.label)

  bc_states <- vector("list", n_barcodes)
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
        origin_height - heights[root], edit_rate[bc], silencing_rate[bc], edit_probs, missing_state
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
        branch_length, edit_rate[bc], silencing_rate[bc], edit_probs, missing_state
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

    # nodes x sites character matrix
    # (explicit matrix(): vapply drops to a vector when n_sites = 1)
    bc_states[[bc]] <- matrix(
      unlist(lapply(tree_df$node, function(nd) state_at[[as.character(nd)]])),
      ncol = n_sites, byrow = TRUE
    )
  }

  if (write_as_string) {
    out <- data.frame(node = tree_df$node)
    for (bc in seq_len(n_barcodes)) {
      out[[paste0("barcode_", bc)]] <- apply(bc_states[[bc]], 1, paste, collapse = "_")
    }
    return(out)
  }

  out <- lapply(bc_states, function(m) {
    if (numeric_states) storage.mode(m) <- "integer"
    colnames(m) <- paste0("site_", seq_len(n_sites))
    data.frame(node = tree_df$node, m, stringsAsFactors = FALSE)
  })
  names(out) <- paste0("barcode_", seq_len(n_barcodes))
  out
}


# Evolve one barcode copy along one branch: 
# silencing happens at most once (and overwrites everything); 
# otherwise a Poisson number of edits fills the free sites left to right
.evolve_branch_seq <- function(parent_state, parent_silenced, branch_length,
                               edit_rate, silencing_rate, edit_probs, missing_state) {
  state <- parent_state

  # already silenced upstream, or no time on this branch: nothing to do
  if (parent_silenced || branch_length <= 0) {
    return(list(state = state, silenced = parent_silenced))
  }

  # silencing: independent of editing and overwrites any edits on this branch,
  # so only whether it occurs matters, not when
  if (silencing_rate > 0 &&
      stats::runif(1) < 1 - exp(-silencing_rate * branch_length)) {
    state[] <- missing_state
    return(list(state = state, silenced = TRUE))
  }

  # editing: number of events is Poisson; sites fill left to right,
  # so events beyond the number of free sites are discarded
  free <- which(state == "0")
  n_edits <- min(stats::rpois(1, edit_rate * branch_length), length(free))
  if (n_edits > 0) {
    state[free[seq_len(n_edits)]] <-
      as.character(sample.int(length(edit_probs), size = n_edits, replace = TRUE, prob = edit_probs))
  }

  list(state = state, silenced = FALSE)
}
