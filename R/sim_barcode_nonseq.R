#' Non-sequential (TiDeTree-style) Barcode Simulator
#'
#' Simulates barcodes under the Edit-and-Silencing model: each of k target
#' sites independently evolves under a continuous-time Markov chain with
#' three fates : staying unedited, transitioning to one of E
#' edit outcomes, or being silenced
#'
#' @param tree a treedata object (as returned by sim_adb_origin_samp);
#'   node identity is taken from tree@data$node
#' @param k number of target sites per barcode
#' @param edit_rates numeric vector of length E, one rate per distinct edit
#'   outcome
#' @param silencing_rate scalar silencing rate, active across the whole tree
#' @param edit_height time (from origin) at which the editing window ends
#'   (i.e. the more recent boundary, "duration between onset of edit and
#'   sampling of the cells")
#' @param edit_duration length of the editing window; the window spans
#'   (edit_height - edit_duration, edit_height]
#' @param dropout_p probability an individual tip site is masked to the
#'   sentinel (missing/silenced) state after simulation; 0 = no dropout
#' @param m number of independent barcode copies per cell
#'
#' @return a data frame with one row per node in the tree and one column
#'   per barcode copy: \code{node}, \code{barcode_1}, ..., \code{barcode_m},
#'   each entry a length-k integer vector-as-string (sites separated by
#'   "_"), values 0 = unedited, 1..E = edit outcome, E+1 = silenced/dropout
#' @export
sim_barcode_nonseq <- function(tree, k, edit_rates, silencing_rate,
                                edit_height, edit_duration,
                                dropout_p = 0, m = 1) {

  stopifnot(methods::is(tree, "treedata"))
  E <- length(edit_rates)
  nstates <- E + 2       # 1 = unedited, 2:(E+1) = edited, E+2 = silenced/dropout
  silenced_state <- nstates

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  # Heights are time-since-present (0 = present, root = age of tree, maximum
  # at origin). ape::node.depth.edgelength returns time-since-root (0 at
  # root, increasing toward tips), indexed by ape's standard node numbering
  # (1:Ntip = tips, (Ntip+1):(Ntip+Nnode) = internal) -- the same convention
  # already used for tree_df$node elsewhere in this package, so position i
  # in the returned vector corresponds to node id i.
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  tree_age <- max(depth_from_root)
  heights <- tree_age - depth_from_root   # heights[node_id] = time-since-present

  # process nodes in decreasing height (root/oldest first, tips/present last)
  order_df <- tree_df[order(-heights[tree_df$node]), ]

  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]

  tip_nodes <- tree_df$node[tree_df$status == 1]

  bc_cols <- vector("list", m)
  for (bc in seq_len(m)) {
    state_at <- list()
    state_at[[as.character(root)]] <- rep(1L, k)  # 1 = unedited (1-indexed)

    # root/origin edge, if present: evolve from origin_height down to the root's height
    if (!is.null(root_edge) && root_edge > 0) {
      state_at[[as.character(root)]] <- .evolve_branch(
        state_at[[as.character(root)]], origin_height, heights[root],
        edit_rates, silencing_rate, edit_height, edit_duration
      )
    }

    for (i in seq_len(nrow(order_df))) {
      node <- order_df$node[i]
      if (node == root) next
      parent <- order_df$parent[i]
      state_at[[as.character(node)]] <- .evolve_branch(
        state_at[[as.character(parent)]], heights[parent], heights[node],
        edit_rates, silencing_rate, edit_height, edit_duration
      )
    }

    if (dropout_p > 0) {
      for (nd in tip_nodes) {
        key <- as.character(nd)
        mask <- stats::runif(k) < dropout_p
        state_at[[key]][mask] <- silenced_state
      }
    }

    bc_cols[[bc]] <- vapply(
      tree_df$node,
      function(nd) paste(state_at[[as.character(nd)]] - 1L, collapse = "_"), # back to 0-indexed for output
      character(1)
    )
  }

  out <- data.frame(node = tree_df$node)
  for (bc in seq_len(m)) out[[paste0("barcode_", bc)]] <- bc_cols[[bc]]
  out
}

# Closed-form transition probability matrix for the Edit-and-Silencing model
# over one (already window-clipped) time segment of length `delta`.
# row = from-state, col = to-state, 1-indexed:
#   1 = unedited, 2:(E+1) = edited outcomes, E+2 = silenced
.edit_silencing_transition_probs <- function(edit_rates, silencing_rate, delta, in_edit_window) {
  E <- length(edit_rates)
  nstates <- E + 2
  silenced <- nstates

  P <- matrix(0, nstates, nstates)
  exp_loss <- exp(-delta * silencing_rate)

  diag(P) <- exp_loss
  P[, silenced] <- 1 - exp_loss
  P[silenced, ] <- 0
  P[silenced, silenced] <- 1

  edit_rate_sum <- sum(edit_rates)
  if (in_edit_window && edit_rate_sum > 0) {
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

# Draw next states for all k sites at once, given a shared transition matrix.
.evolve_segment <- function(state, edit_rates, silencing_rate, delta, in_window) {
  if (delta <= 0) return(state)
  P <- .edit_silencing_transition_probs(edit_rates, silencing_rate, delta, in_window)
  vapply(state, function(s) sample.int(nrow(P), 1, prob = P[s, ]), integer(1))
}

# Walk one branch, splitting at the editing-window boundaries so each
# segment passed to .evolve_segment is either fully inside or fully
# outside the window (mirrors SimulatedAlignment's helper-node splitting).
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
