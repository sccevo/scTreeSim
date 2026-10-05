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
#' @param n_sites number of sites in the barcode
#' @param edit_rate scalar, or length-n_sites vector giving a per-site rate of any edit per time unit
#' @param edit_probs numeric vector of length E, summing to 1: the relative frequencies of the edit outcomes 1,...,E
#' @param silencing_rate scalar, or length-n_sites vector giving a per-site silencing rate
#' @param edit_height start of the editing window, as a time before the present;
#'   defaults to the origin time of the tree (\code{tree@phylo$origin}) if available, otherwise to the tree height
#' @param edit_duration length of the editing window; sites can be edited
#'   between \code{edit_height} and \code{edit_height - edit_duration} before the present;
#'   defaults to \code{edit_height}, i.e. editing from the start of the window until the present
#' @param dropout_prob scalar, or length-n_sites vector giving a per-site dropout probability; 0 = no dropout
#'
#' @return a data frame with one row per node in the tree and one column
#'   per site: \code{node}, \code{site_1}, ..., \code{site_k}, values
#'   0 = unedited, 1..E = edit outcome, E+1 = silenced/dropout
#' @export
sim_barcode_nonseq <- function(tree, n_sites, edit_rate, edit_probs, silencing_rate,
                               edit_height = NULL, edit_duration = NULL, dropout_prob = 0) {

  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }
  # default editing window: from the origin (or the root) until the present
  if (is.null(edit_height)) {
    edit_height <- if (!is.null(tree@phylo$origin)) tree@phylo$origin else max(ape::node.depth.edgelength(tree@phylo))
  }
  if (is.null(edit_duration)) edit_duration <- edit_height
  stopifnot(
    is.numeric(n_sites), length(n_sites) == 1, n_sites >= 1, n_sites == round(n_sites),
    is.numeric(edit_rate), !anyNA(edit_rate), all(edit_rate >= 0),
    is.numeric(edit_probs), length(edit_probs) >= 1, !anyNA(edit_probs), all(edit_probs >= 0),
    isTRUE(all.equal(sum(edit_probs), 1)),
    is.numeric(silencing_rate), !anyNA(silencing_rate), all(silencing_rate >= 0),
    is.numeric(edit_height), length(edit_height) == 1, !is.na(edit_height),
    is.numeric(edit_duration), length(edit_duration) == 1, !is.na(edit_duration), edit_duration >= 0,
    is.numeric(dropout_prob), !anyNA(dropout_prob), all(dropout_prob >= 0 & dropout_prob <= 1)
  )
  n_sites <- as.integer(n_sites)
  E <- length(edit_probs)
  nstates <- E + 2L
  silenced_state <- nstates

  # allow a single shared rate to stand in for all n_sites sites
  if (length(edit_rate) == 1) edit_rate <- rep(edit_rate, n_sites)
  if (length(silencing_rate) == 1) silencing_rate <- rep(silencing_rate, n_sites)
  if (length(dropout_prob) == 1) dropout_prob <- rep(dropout_prob, n_sites)
  stopifnot(length(edit_rate) == n_sites, length(silencing_rate) == n_sites, length(dropout_prob) == n_sites)

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  # heights: time-since-present (0 = present, root = age of tree)
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  tree_age <- max(depth_from_root)
  heights <- tree_age - depth_from_root

  # visit parent before child so state can propagate down the tree; use the
  # topology rather than heights, which tie on zero-length branches
  child_order <- ape::reorder.phylo(tree@phylo, order = "cladewise")$edge[, 2]
  order_df <- tree_df[match(child_order, tree_df$node), ]

  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]

  # tips are nodes 1..Ntip (all sampled cells, since the tree is ultrametric)
  tip_nodes <- seq_along(tree@phylo$tip.label)

  # state_at: n_sites-length integer state vector per node
  # 1-indexed internally
  # converted back to 0-indexed only in the final output
  state_at <- list()
  state_at[[as.character(root)]] <- rep(1L, n_sites)

  # evolve along the root/origin edge first, if the tree has one
  if (!is.null(root_edge) && root_edge > 0) {
    state_at[[as.character(root)]] <- .evolve_branch(
      state_at[[as.character(root)]], origin_height, heights[root],
      edit_rate, edit_probs, silencing_rate, edit_height, edit_duration
    )
  }

  # then walk every remaining branch, parent state -> child state
  for (i in seq_len(nrow(order_df))) {
    node <- order_df$node[i]
    if (node == root) next
    parent <- order_df$parent[i]
    state_at[[as.character(node)]] <- .evolve_branch(
      state_at[[as.character(parent)]], heights[parent], heights[node],
      edit_rate, edit_probs, silencing_rate, edit_height, edit_duration
    )
  }

  # dropout applied after simulation
  for (nd in tip_nodes) {
    key <- as.character(nd)
    mask <- stats::runif(n_sites) < dropout_prob
    state_at[[key]][mask] <- silenced_state
  }

  out <- data.frame(node = tree_df$node)
  state_mat <- t(vapply(tree_df$node, function(nd) state_at[[as.character(nd)]] - 1L, integer(n_sites)))
  for (site in seq_len(n_sites)) out[[paste0("site_", site)]] <- state_mat[, site]
  out
}


# Closed-form transition probability matrix for one time segment of length delta,
# for a single site's own edit and silencing rates.
# row = from-state, col = to-state
# 1-indexed: 1 = unedited, 2:(E+1) = edited outcomes, E+2 = silenced.
.edit_silencing_transition_probs <- function(edit_rate, edit_probs, silencing_rate, delta, in_edit_window) {
  E <- length(edit_probs)
  nstates <- E + 2L
  silenced <- nstates

  P <- matrix(0, nstates, nstates)
  exp_loss <- exp(-delta * silencing_rate)

  # baseline outcomes can be stay or silenced, applies to every state
  # (silencing_rate is a scalar here: one site's rate)
  diag(P) <- exp_loss
  P[, silenced] <- 1 - exp_loss
  P[silenced, ] <- 0
  P[silenced, silenced] <- 1  # silenced is absorbing

  if (in_edit_window && edit_rate > 0) {
    # stay unedited, move to one of the E edit outcomes, or silence
    p_stay <- exp(-delta * (silencing_rate + edit_rate))
    P[1, 1] <- p_stay
    P[1, 1 + seq_len(E)] <- edit_probs * (exp_loss - p_stay)
    P[1, silenced] <- 1 - exp_loss
  }
  P
}


# Draw the next state for each site independently, using that site's own
# edit and silencing rates (edit_probs are shared across sites).
# Sites sharing both rates share a transition matrix, which is built once.
.evolve_segment <- function(state, edit_rate, edit_probs, silencing_rate, delta, in_window) {
  if (delta <= 0) return(state)
  new_state <- state
  rates <- unique(data.frame(edit = edit_rate, silencing = silencing_rate))
  for (r in seq_len(nrow(rates))) {
    sites <- which(edit_rate == rates$edit[r] & silencing_rate == rates$silencing[r])
    P <- .edit_silencing_transition_probs(rates$edit[r], edit_probs, rates$silencing[r], delta, in_window)
    for (site in sites) {
      new_state[site] <- sample.int(nrow(P), 1, prob = P[state[site], ])
    }
  }
  new_state
}


# Split one branch at the editing-window boundaries so every segment
# passed to .evolve_segment is either fully inside or fully outside the window,
# then evolve through each segment in turn.
.evolve_branch <- function(parent_state, parent_height, child_height,
                           edit_rate, edit_probs, silencing_rate, edit_height, edit_duration) {
  window_top <- edit_height
  window_bot <- edit_height - edit_duration

  bounds <- sort(unique(c(parent_height, child_height, window_top, window_bot)), decreasing = TRUE)
  bounds <- bounds[bounds <= parent_height & bounds >= child_height]

  state <- parent_state
  for (i in seq_len(length(bounds) - 1)) {
    seg_top <- bounds[i]; seg_bot <- bounds[i + 1]
    delta <- seg_top - seg_bot
    in_window <- (seg_bot >= window_bot) && (seg_bot < window_top)
    state <- .evolve_segment(state, edit_rate, edit_probs, silencing_rate, delta, in_window)
  }
  state
}
