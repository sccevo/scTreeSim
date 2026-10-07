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
#' @export
sim_barcode_nonseq <- function(tree, n_targets, edit_rate, edit_probs,
                               edit_height = NULL, edit_duration = NULL,
                               silencing_rate = 0, dropout_prob = 0,
                               missing_state = length(edit_probs) + 1L) {

  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }
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
  if (!is.null(tree@phylo$origin) && !isTRUE(all.equal(tree@phylo$origin, origin_height))) {
    warning("`tree@phylo$origin` differs from the tree height plus root edge; the latter is used.", call. = FALSE)
  }

  # default editing window: from the origin (or the root) until the present
  if (is.null(edit_height)) edit_height <- origin_height
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
  nstates <- E + 2L
  silenced_state <- nstates

  # allow a single shared rate to stand in for all n_targets sites
  if (length(edit_rate) == 1) edit_rate <- rep(edit_rate, n_targets)
  if (length(silencing_rate) == 1) silencing_rate <- rep(silencing_rate, n_targets)
  if (length(dropout_prob) == 1) dropout_prob <- rep(dropout_prob, n_targets)

  # tips are nodes 1..Ntip (all sampled cells, since the tree is ultrametric)
  tip_nodes <- seq_along(tree@phylo$tip.label)

  # state_at: n_targets-length integer state vector per node
  # 1-indexed internally
  # converted back to 0-indexed only in the final output
  state_at <- list()
  state_at[[as.character(root)]] <- rep(1L, n_targets)

  # evolve along the root/origin edge first, if the tree has one
  if (!is.null(root_edge) && root_edge > 0) {
    state_at[[as.character(root)]] <- .evolve_branch_nonseq(
      state_at[[as.character(root)]], origin_height, heights[root],
      edit_rate, edit_probs, silencing_rate, edit_height, edit_duration
    )
  }

  # then walk every remaining branch, parent state -> child state
  for (i in seq_len(nrow(order_df))) {
    node <- order_df$node[i]
    if (node == root) next
    parent <- order_df$parent[i]
    state_at[[as.character(node)]] <- .evolve_branch_nonseq(
      state_at[[as.character(parent)]], heights[parent], heights[node],
      edit_rate, edit_probs, silencing_rate, edit_height, edit_duration
    )
  }

  # dropout applied after simulation
  for (nd in tip_nodes) {
    key <- as.character(nd)
    mask <- stats::runif(n_targets) < dropout_prob
    state_at[[key]][mask] <- silenced_state
  }

  # nodes x targets matrix, 0-indexed (explicit matrix(): vapply drops to a vector when n_targets = 1)
  state_mat <- matrix(
    unlist(lapply(tree_df$node, function(nd) state_at[[as.character(nd)]] - 1L)),
    ncol = n_targets, byrow = TRUE
  )
  # the internal silenced code E + 1 is relabelled to missing_state
  if (is.character(missing_state)) {
    state_mat[] <- as.character(state_mat)
    state_mat[state_mat == as.character(E + 1L)] <- missing_state
  } else {
    state_mat[state_mat == E + 1L] <- as.integer(missing_state)
  }

  out <- data.frame(node = tree_df$node)
  for (site in seq_len(n_targets)) out[[paste0("site_", site)]] <- state_mat[, site]
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
.evolve_segment_nonseq <- function(state, edit_rate, edit_probs, silencing_rate, delta, in_window) {
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
# passed to .evolve_segment_nonseq is either fully inside or fully outside the window,
# then evolve through each segment in turn.
.evolve_branch_nonseq <- function(parent_state, parent_height, child_height,
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
    state <- .evolve_segment_nonseq(state, edit_rate, edit_probs, silencing_rate, delta, in_window)
  }
  state
}
