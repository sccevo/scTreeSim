# Clock models convert branch lengths (time) into evolutionary distances (expected substitutions
# per site). A clock holds
# - rate: substitutions per site per unit of time,
# - multipliers: function(branches) returning a mean-1 factor on the rate for each branch,
# - spike: an optional extra distance per branch, not scaled by time.
# `branches` is the node table of the tree, one row per node for the branch ending in it (the root's
# row is the root edge), with columns `node`, `parent`, `length` (time) and the node data (e.g. `type`).
# Distance of a branch = multiplier * rate * length (+ spike).

.clock_model <- function(rate, multipliers, spike = NULL, ...) {
  structure(list(rate = rate, multipliers = multipliers, spike = spike, ...), class = "clock_model")
}

.check_scalar <- function(x, lower = 0, strict = FALSE) {
  name <- deparse(substitute(x))
  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) && (if (strict) x > lower else x >= lower)
  if (!ok) {
    stop("`", name, "` must be a finite number ", if (strict) ">" else ">=", " ", lower, ".", call. = FALSE)
  }
}

#' Strict clock
#'
#' All branches evolve at the same `rate`, so the evolutionary distance of a branch is
#' `rate * branch length`, with branch lengths in time.
#'
#' @param rate substitutions per site per unit of time, non-negative
#'
#' @return an object of class `clock_model`, see [clock_tree()]
#' @family clock models
#' @export
clock_strict <- function(rate = 1) {
  .check_scalar(rate)
  .clock_model(rate, function(branches) rep(1, nrow(branches)))
}

#' Uncorrelated relaxed clock
#'
#' Every branch, including the root edge, evolves at `multiplier * rate`, with multipliers drawn
#' independently per branch from a distribution with mean 1. Each call of [clock_tree()] is a new draw.
#'
#' * `"lognormal"`: `exp(sd * z - sd^2 / 2)` with `z ~ N(0, 1)`, the BEAST2 uncorrelated lognormal clock
#'   (`LogNormalDistributionModel` with `M = 1`, `S = sd`, `meanInRealSpace = true`), drawn
#'   continuously and not normalised. `sd = 0` is the strict clock.
#' * `"gamma"`: `Gamma(shape, rate = shape)`, coefficient of variation `1 / sqrt(shape)`.
#' * `"exponential"`: `Exp(1)`, the gamma with `shape = 1`.
#'
#' @param rate mean substitutions per site per unit of time, non-negative
#' @param dist distribution of the multipliers
#' @param sd standard deviation of the log multiplier, non-negative (`"lognormal"` only)
#' @param shape shape of the multiplier, positive (`"gamma"` only)
#'
#' @return an object of class `clock_model`, see [clock_tree()]
#' @family clock models
#' @export
clock_relaxed <- function(rate = 1, dist = c("lognormal", "gamma", "exponential"), sd = NULL, shape = NULL) {
  .check_scalar(rate)
  dist <- match.arg(dist)

  # each distribution takes its own parameter and no other
  needed <- switch(dist, lognormal = "sd", gamma = "shape", exponential = character(0))
  given <- c("sd", "shape")[c(!is.null(sd), !is.null(shape))]
  if (!setequal(given, needed)) {
    stop("The ", dist, " clock takes ",
         switch(dist, lognormal = "`sd` only", gamma = "`shape` only", exponential = "neither `sd` nor `shape`"),
         ".", call. = FALSE)
  }

  multipliers <- switch(dist,
    lognormal = {
      .check_scalar(sd)
      function(branches) exp(sd * stats::rnorm(nrow(branches)) - sd^2 / 2)
    },
    gamma = {
      .check_scalar(shape, strict = TRUE)
      function(branches) stats::rgamma(nrow(branches), shape = shape, rate = shape)
    },
    exponential = function(branches) stats::rexp(nrow(branches))
  )
  .clock_model(rate, multipliers, dist = dist, sd = sd, shape = shape)
}

#' Evolutionary distances on a tree
#'
#' `clock_tree()` draws the evolutionary distance (expected substitutions per site) of every branch
#' under a clock model and stores it in a `distance` column of the tree data, one value per node for
#' the branch ending in it (the root's value belongs to the root edge). Branch lengths stay in time
#' units. Random clocks give a new draw on every call, so storing the distances lets several
#' simulations share one draw.
#'
#' `distance_tree()` returns the tree with the distances as branch lengths.
#'
#' @param tree a treedata object
#' @param clock a clock model, e.g. [clock_strict()]
#'
#' @return a treedata object
#' @family clock models
#' @export
clock_tree <- function(tree, clock = clock_strict()) {
  stopifnot(methods::is(tree, "treedata"), inherits(clock, "clock_model"))
  branches <- as.data.frame(tibble::as_tibble(tree))
  root_edge <- if (is.null(tree@phylo$root.edge)) 0 else tree@phylo$root.edge
  branches$length <- ifelse(branches$parent == branches$node, root_edge, branches$branch.length)

  d <- numeric(nrow(branches))
  d[branches$node] <- clock$multipliers(branches) * clock$rate * branches$length

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
