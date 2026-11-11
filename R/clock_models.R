# Clock models convert branch lengths (time) into evolutionary distances (expected substitutions
# per site). A clock holds
# - rate: substitutions per site per unit of time,
# - multipliers: a function of the node table giving each branch a mean-1 factor on the rate,
# - spike: an optional extra distance per branch, not scaled by time.
# Distance of the branch ending in a node = multiplier * rate * branch length (+ spike).

.clock_model <- function(rate, multipliers, spike = NULL) {
  structure(list(rate = rate, multipliers = multipliers, spike = spike), class = "clock_model")
}

#' Strict clock
#'
#' All branches evolve at the same `rate`, so the evolutionary distance of a branch is
#' `rate * branch length`, with branch lengths in time.
#'
#' @param rate substitutions per site per unit of time, non-negative
#'
#' @return an object of class `clock_model`, for the `clock` argument of [sim_binary()]
#' @family clock models
#' @export
clock_strict <- function(rate = 1) {
  stopifnot(is.numeric(rate), length(rate) == 1, is.finite(rate), rate >= 0)
  .clock_model(rate, function(tree_df) rep(1, nrow(tree_df)))
}

#' Evolutionary distances on a tree
#'
#' `clock_tree()` draws the evolutionary distance (expected substitutions per site) of every branch
#' under a clock model and stores it in a `distance` column of the tree data, one value per node for
#' the branch ending in it (the root's value belongs to the root edge). Branch lengths stay in time
#' units. Random clocks give a new draw on every call; pass the result to [sim_binary()] with
#' `clock = NULL` to simulate on a given draw.
#'
#' `distance_tree()` returns the tree with the distances as branch lengths
#'
#' @param tree a treedata object
#' @param clock a clock model, e.g. [clock_strict()]
#'
#' @return a treedata object
#' @family clock models
#' @export
clock_tree <- function(tree, clock = clock_strict()) {
  stopifnot(methods::is(tree, "treedata"), inherits(clock, "clock_model"))
  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()

  # time length of the branch ending in each node; the root takes the root edge (0 if none)
  lengths <- tree_df$branch.length
  is_root <- tree_df$parent == tree_df$node
  lengths[is_root] <- if (is.null(tree@phylo$root.edge)) 0 else tree@phylo$root.edge

  d <- numeric(nrow(tree_df))
  d[tree_df$node] <- clock$multipliers(tree_df) * clock$rate * lengths

  # the column is per node, so it goes stale if branches are merged later (collapsing, pruning)
  if (nrow(tree@data) == 0) tree@data <- tibble::tibble(node = seq_along(d))
  tree@data$distance <- d[tree@data$node]
  tree
}

#' @rdname clock_tree
#' @export
distance_tree <- function(tree) {
  d <- .node_distances(tree)
  # each edge takes the distance of its child node
  phy <- tree@phylo
  phy$edge.length <- d[phy$edge[, 2]]
  if (!is.null(phy$root.edge)) phy$root.edge <- d[ape::Ntip(phy) + 1]
  tree@phylo <- phy
  tree
}

# the `distance` column as a vector indexed by node id
.node_distances <- function(tree) {
  stopifnot(methods::is(tree, "treedata"))
  d <- rep(NA_real_, ape::Ntip(tree@phylo) + tree@phylo$Nnode)
  if ("distance" %in% names(tree@data)) d[tree@data$node] <- tree@data$distance
  if (anyNA(d)) stop("`tree` needs a `distance` for every node; see clock_tree().", call. = FALSE)
  d
}
