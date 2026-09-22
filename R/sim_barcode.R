#' Simulator of Sequentially-Edited (TypeWriter-style) Barcodes
#' 
#' Parts of this code contains adaptations of code from Antoine Zwaans's implementation 
#' of sequential barcode simulator for SciPhy software with dropout features
#' source: https://github.com/azwaans/tidetree-dropout/blob/main/src/tidetree/simulation/SimulatedAlignment.java
#' 
#' R version of the code has been adapted from Nicola Mulberry's implementation
#' source: https://github.com/sccevo/scTreeSim/blob/barcode_simulator/R/typewriter_barcodes.R
#' The original authors and their respective licenses are retained where applicable.
#'
#' Simulates barcodes under a sequential/left-to-right editing process
#' (unedited sites are filled in position order at rate lambda, as a
#' Poisson process, absorbing once all k sites are edited), with an
#' additional whole-barcode silencing state competing against editing.
#'
#' @param tree a treedata object
#' @param k number of target sites per barcode
#' @param lambda_vec array of barcode-specific editing rates (length m)
#' @param silencing_rate_vec scalar or length-m rate(s) at which an entire
#'   barcode is silenced, racing against the editing process
#' @param dropout_p_vec scalar or length-m probability that a
#'   barcode copy drops out at a tip
#' @param m number of tapes per cell
#' @param chars array of unique characters/state labels to be inserted
#'
#' @return a data frame with one row per node in the tree and one column
#'   per barcode copy: \code{node}, \code{barcode_1}, ..., \code{barcode_m},
#'   each entry a length-k string with positions separated by "_"
#'   (e.g. "0_A_B_0_0"), where "0" = unedited and "-" = silenced/dropout
#'   (the silenced state; every site reads "-" once a barcode is silenced
#'   or dropped out)
#' @export
typewriter_barcodes <- function(tree, k, lambda_vec, silencing_rate_vec = (0),
                                dropout_p_vec = (0), m = length(lambda_vec), chars) {

  stopifnot(methods::is(tree, "treedata"))

  # allow a single shared rate to stand in for all m barcode copies
  if (length(silencing_rate_vec) == 1) {
    silencing_rate_vec <- rep(silencing_rate_vec[1], m)
  }
  if (length(dropout_p_vec) == 1) {
    dropout_p_vec <- rep(dropout_p_vec[1], m)
  }
  stopifnot(length(lambda_vec) == m)
  stopifnot(length(silencing_rate_vec) == m)
  stopifnot(length(dropout_p_vec) == m)

  silenced_state <- "-"
  # uniform draw over insertion characters
  sample_p <- rep(1 / length(chars), length(chars))

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

  tip_nodes <- tree_df$node[tree_df$status == 1]

  bc_cols <- vector("list", m)
  for (bc in seq_len(m)) {
    lambda <- lambda_vec[bc]
    silencing_rate <- silencing_rate_vec[bc]
    dropout_p <- dropout_p_vec[bc]

    # state_at: barcode string per node; silenced_at: whether that barcode
    # has already been silenced (needed since silencing must persist down
    # every descendant branch once it happens)
    state_at <- list()
    silenced_at <- list()
    state_at[[as.character(root)]] <- rep("0", k)
    silenced_at[[as.character(root)]] <- FALSE

    # evolve along the root/origin edge first, if the tree has one
    if (!is.null(root_edge) && root_edge > 0) {
      res <- .evolve_typewriter_branch(
        state_at[[as.character(root)]], silenced_at[[as.character(root)]],
        origin_height - heights[root], lambda, silencing_rate, chars, sample_p, silenced_state
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

      res <- .evolve_typewriter_branch(
        state_at[[as.character(parent)]], silenced_at[[as.character(parent)]],
        branch_length, lambda, silencing_rate, chars, sample_p, silenced_state
      )
      state_at[[as.character(node)]] <- res$state
      silenced_at[[as.character(node)]] <- res$silenced
    }

    # apply dropout after simulation,
    # skipped for barcodes already fully silenced (nothing left to mask)
    if (dropout_p > 0) {
      for (nd in tip_nodes) {
        key <- as.character(nd)
        if (!silenced_at[[key]] && stats::runif(1) < dropout_p) {
          state_at[[key]][] <- silenced_state
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
  for (bc in seq_len(m)) out[[paste0("barcode_", bc)]] <- bc_cols[[bc]]
  out
}

# Evolve one barcode copy along one branch. At each event, draw whether it
# is an edit or silencing

.evolve_typewriter_branch <- function(parent_state, parent_silenced, branch_length,
                                      lambda, silencing_rate, chars, sample_p, silenced_state) {
  state <- parent_state
  k <- length(state)

  # already silenced upstream, or no time on this branch: nothing to do
  if (parent_silenced || branch_length <= 0) {
    return(list(state = state, silenced = parent_silenced))
  }

  t <- 0

  repeat {
    # once every site is edited, only silencing can still happen
    saturated <- all(state != "0")
    edit_rate <- if (saturated) 0 else lambda
    event_rate <- edit_rate + silencing_rate

    if (event_rate <= 0) break  # no editing left to do and no silencing possible

    # draw time to the next event (edit or silencing) on the combined clock
    t <- t + stats::rexp(1, event_rate)
    if (t > branch_length) break  # no more events fit in this branch

    # which of the two competing events fired, weighted by relative rate
    is_silencing <- stats::runif(1) < (silencing_rate / event_rate)
    if (is_silencing) {
      state[] <- silenced_state
      return(list(state = state, silenced = TRUE))
    }

    # edit event: fill the leftmost still-unedited site
    pos <- min(which(state == "0"))
    state[pos] <- sample(chars, size = 1, prob = sample_p)
  }

  list(state = state, silenced = FALSE)
}