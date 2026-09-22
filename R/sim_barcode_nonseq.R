#' Non-sequential (TiDeTree-style) Barcode Simulator
#'
#' This code contains adaptations of code from Sophie Seidel's TideTree implementation
#' and Antoine Zwaans's implementation with dropout
#' The original authors and their respective licenses are retained where applicable.
#' source: https://github.com/azwaans/tidetree-dropout/blob/main/src/tidetree/simulation/SimulatedAlignment.java
#'  
#' Simulates a single barcode with k independent sites under the
#' Edit-and-Silencing model: each site evolves as a continuous-time Markov
#' chain with three fates: staying unedited, transitioning to one of E
#' edit outcomes, or being silenced.
#'
#' @param tree a treedata object; node identity is taken from tree@data$node
#' @param k number of sites in the barcode
#' @param edit_rates numeric vector of length E, one rate per distinct edit
#'   outcome
#' @param silencing_rate scalar, or length-k vector giving a per-site
#'   silencing rate
#' @param edit_height time (from origin) at which the editing window ends
#' @param edit_duration length of the editing window
#' @param dropout_p scalar, or length-k vector giving a per-site dropout
#'   probability; 0 = no dropout
#'
#' @return a data frame with one row per node in the tree and one column
#'   per site: \code{node}, \code{site_1}, ..., \code{site_k}, values
#'   0 = unedited, 1..E = edit outcome, E+1 = silenced/dropout
#' @export
sim_barcode_nonseq <- function(tree, k, edit_rates, silencing_rate,
                               edit_height, edit_duration, dropout_p = 0) {

  stopifnot(methods::is(tree, "treedata"))
  E <- length(edit_rates)
  # 1 = unedited, 2:(E+1) = edited outcomes, E+2 = silenced/dropout
  nstates <- E + 2L
  silenced_state <- nstates

  # allow a single shared rate to stand in for all k sites
  if (length(silencing_rate) == 1) silencing_rate <- rep(silencing_rate, k)
  if (length(dropout_p) == 1) dropout_p <- rep(dropout_p, k)
  stopifnot(length(silencing_rate) == k, length(dropout_p) == k)

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  # heights: time-since-present (0 = present, root = age of tree)
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  tree_age <- max(depth_from_root)
  heights <- tree_age - depth_from_root

  # visit parent before child so state can propagate down the tree
  order_df <- tree_df[order(-heights[tree_df$node]), ]

  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]

  tip_nodes <- tree_df$node[tree_df$status == 1]

  # state_at: k-length integer state vector per node
  # 1-indexed internally
  # converted back to 0-indexed only in the final output
  state_at <- list()
  state_at[[as.character(root)]] <- rep(1L, k)

  # evolve along the root/origin edge first, if the tree has one
  if (!is.null(root_edge) && root_edge > 0) {
    state_at[[as.character(root)]] <- .evolve_branch(
      state_at[[as.character(root)]], origin_height, heights[root],
      edit_rates, silencing_rate, edit_height, edit_duration
    )
  }

  # then walk every remaining branch, parent state -> child state
  for (i in seq_len(nrow(order_df))) {
    node <- order_df$node[i]
    if (node == root) next
    parent <- order_df$parent[i]
    state_at[[as.character(node)]] <- .evolve_branch(
      state_at[[as.character(parent)]], heights[parent], heights[node],
      edit_rates, silencing_rate, edit_height, edit_duration
    )
  }

  # dropout applied after simulation
  for (nd in tip_nodes) {
    key <- as.character(nd)
    mask <- stats::runif(k) < dropout_p
    state_at[[key]][mask] <- silenced_state
  }

  out <- data.frame(node = tree_df$node)
  state_mat <- t(vapply(tree_df$node, function(nd) state_at[[as.character(nd)]] - 1L, integer(k)))
  for (site in seq_len(k)) out[[paste0("site_", site)]] <- state_mat[, site]
  out
}

# Closed-form transition probability matrix for one time segment of length delta,
# for a single site's own silencing rate.
# row = from-state, col = to-state
# 1-indexed: 1 = unedited, 2:(E+1) = edited outcomes, E+2 = silenced.
.edit_silencing_transition_probs <- function(edit_rates, silencing_rate, delta, in_edit_window) {
  E <- length(edit_rates)
  nstates <- E + 2L
  silenced <- nstates

  P <- matrix(0, nstates, nstates)
  exp_loss <- exp(-delta * silencing_rate)

  # baseline outcomees can be stay or silenced, applies to every state
  diag(P) <- exp_loss
  P[, silenced] <- 1 - exp_loss
  P[silenced, ] <- 0
  P[silenced, silenced] <- 1  # silenced is absorbing

  edit_rate_sum <- sum(edit_rates)
  if (in_edit_window && edit_rate_sum > 0) {
    # stay unedited, move to one of the E edit outcomes, or silence
    total_rate <- silencing_rate + edit_rate_sum
    p_stay <- exp(-delta * total_rate)
    P[1, 1] <- p_stay
    for (i in seq_len(E)) {
      P[1, i + 1] <- (edit_rates[i] * exp_loss - edit_rates[i] * p_stay) / edit_rate_sum
    }
    P[1, silenced] <- 1 - exp_loss
  }
  P
}

# Draw the next state for each site independently, using that site's own
# silencing rate (edit_rates are shared across sites).
.evolve_segment <- function(state, edit_rates, silencing_rate, delta, in_window) {
  if (delta <= 0) return(state)
  k <- length(state)
  new_state <- integer(k)
  for (site in seq_len(k)) {
    P <- .edit_silencing_transition_probs(edit_rates, silencing_rate[site], delta, in_window)
    new_state[site] <- sample.int(nrow(P), 1, prob = P[state[site], ])
  }
  new_state
}

# Split one branch at the editing-window boundaries so every segment
# passed to .evolve_segment is either fully inside or fully outside the
# window, then evolve through each segment in turn.
.evolve_branch <- function(parent_state, parent_height, child_height,
                           edit_rates, silencing_rate, edit_height, edit_duration) {
  window_top <- edit_height
  window_bot <- edit_height - edit_duration

  bounds <- sort(unique(c(parent_height, child_height, window_top, window_bot)), decreasing = TRUE)
  bounds <- bounds[bounds <= parent_height & bounds >= child_height]

  state <- parent_state
  for (i in seq_len(length(bounds) - 1)) {
    seg_top <- bounds[i]; seg_bot <- bounds[i + 1]
    delta <- seg_top - seg_bot
    in_window <- (seg_bot >= window_bot) && (seg_bot < window_top)
    state <- .evolve_segment(state, edit_rates, silencing_rate, delta, in_window)
  }
  state
}