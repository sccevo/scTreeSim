# Clock models convert branch lengths (time) into evolutionary distances (expected substitutions
# per site). A clock holds
# - rate: substitutions per site per unit of time,
# - multipliers: function(branches) returning a mean-1 factor on the rate for each branch,
# - spike: see spike(), an optional extra distance per branch
# `branches` is the node table of the tree, one row per node for the branch ending in it (the root's
# row is the root edge), with columns `node`, `parent`, `length` (time) and the node data (e.g. `type`).
# Distance of a branch = multiplier * rate * length (+ spike).

.clock_model <- function(rate, multipliers, spike = NULL, ...) {
  if (!is.null(spike) && !inherits(spike, "clock_spike")) {
    stop("`spike` must be NULL or created by spike().", call. = FALSE)
  }
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
#' @param spike optional division-triggered spikes, see [spike()]
#'
#' @return an object of class `clock_model`, see [clock_tree()]
#' @family clock models
#' @export
clock_strict <- function(rate = 1, spike = NULL) {
  .check_scalar(rate)
  .clock_model(rate, function(branches) rep(1, nrow(branches)), spike = spike)
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
#' @inheritParams clock_strict
#'
#' @return an object of class `clock_model`, see [clock_tree()]
#' @family clock models
#' @export
clock_relaxed <- function(rate = 1, dist = c("lognormal", "gamma", "exponential"), sd = NULL, shape = NULL,
                          spike = NULL) {
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
  .clock_model(rate, multipliers, spike = spike, dist = dist, sd = sd, shape = shape)
}


#' Division-triggered spikes
#'
#' A spike is a burst of evolution at a cell division. Each branch whose parent node divided
#' (`status == 2`, or any internal node without a `status` column) gets an independent spike,
#' added to its distance and not scaled by the clock. The root edge gets none.
#'
#' * `shape = Inf` (default): fixed spike `mean` per division. Across `l` sites this gives about
#'   `Poisson(l * mean)` substitutions per daughter branch, as for a genome-wide count of
#'   `mu = l * mean` new mutations per offspring (Dieselhorst & Berg).
#' * finite `shape`: `mean * Gamma(shape, rate = shape)` per division, coefficient of variation
#'   `1 / sqrt(shape)` (MultiTypeSpikeModel, Ciuffi et al.), so substitution counts are
#'   over-dispersed (variance `mu + mu^2 / shape`).
#'
#' Collapsed trees hide divisions on their branches. An optional `hidden_divisions` column of the
#' tree data (counts per node, for the branch ending in it) adds that many spikes to a branch,
#' including the root edge if the tree has one. Leave it out on trees with single-child nodes (`collapse = FALSE`),
#' which show these divisions as nodes.
#'
#' @param mean spike per division and daughter branch in expected substitutions per site,
#'   non-negative
#' @param shape shape of the gamma-distributed spike, positive; `Inf` for a fixed spike
#' @param by_type optional non-negative multipliers of `mean` per type of the dividing (parent)
#'   node, indexed by type + 1 (requires a `type` column); `0` = no spikes for that type, `NA`
#'   only for types that never divide
#'
#' @return an object of class `clock_spike`, for the `spike` argument of the clock models
#' @family clock models
#' @export
spike <- function(mean, shape = Inf, by_type = NULL) {
  .check_scalar(mean)
  if (!(is.numeric(shape) && length(shape) == 1 && !is.na(shape) && shape > 0)) {
    stop("`shape` must be a positive number or Inf.", call. = FALSE)
  }
  if (!is.null(by_type) && !(is.numeric(by_type) && length(by_type) >= 1 && all(by_type >= 0, na.rm = TRUE) &&
                             !all(is.na(by_type)))) {
    stop("`by_type` must be NULL or a vector of non-negative numbers (NA for unused types).", call. = FALSE)
  }

  distances <- function(branches) {
    n_events <- .division_counts(branches)
    parent_row <- match(branches$parent, branches$node)
    size <- mean * .spike_type_multipliers(by_type, branches, parent_row, n_events > 0)
    if (is.infinite(shape)) size * n_events else size * stats::rgamma(nrow(branches), n_events * shape, shape)
  }
  structure(list(mean = mean, shape = shape, by_type = by_type, distances = distances), class = "clock_spike")
}

# Number of divisions on each branch: 1 if its parent node divided (0 for the root edge), plus the
# divisions hidden on the branch (optional `hidden_divisions` column). Sums of N gamma spikes are
# Gamma(N * shape, shape), so the draw in spike() needs no change.
.division_counts <- function(branches) {
  is_root <- branches$parent == branches$node
  divided <- if ("status" %in% names(branches)) {
    branches$node[branches$status %in% 2]
  } else {
    unique(branches$parent[!is_root])  # every internal node
  }
  explicit <- as.numeric(!is_root & branches$parent %in% divided)
  if (is.null(branches$hidden_divisions)) explicit else explicit + branches$hidden_divisions
}

# by_type multiplier of each branch from the type of its parent (1 without by_type, 0 if not spiked);
# NA for missing or unknown types, caught by the check in clock_tree()
.spike_type_multipliers <- function(by_type, branches, parent_row, spiked) {
  if (is.null(by_type)) return(rep(1, nrow(branches)))
  if (is.null(branches$type)) return(rep(NA_real_, nrow(branches)))
  mult <- by_type[match(branches$type[parent_row], seq_along(by_type) - 1)]
  mult[!spiked] <- 0
  mult
}


#' Evolutionary distances on a tree
#'
#' `clock_tree()` draws the evolutionary distance (expected substitutions per site) of every branch
#' under a clock model and stores it in a `distance` column of the tree data, one value per node for
#' the branch ending in it (the root's value belongs to the root edge). Branch lengths stay in time
#' units. Spikes of the clock (see [spike()]) are included. Random clocks and spikes give a new draw on every call, so storing the distances lets several
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
  if (!is.null(clock$spike)) d[branches$node] <- d[branches$node] + clock$spike$distances(branches)
  # no root edge, no evolution on it (hidden divisions on the root row are dropped)
  if (is.null(tree@phylo$root.edge)) d[ape::Ntip(tree@phylo) + 1L] <- 0
  if (anyNA(d) || any(d < 0)) {
    stop("The clock gave missing or negative distances; check the tree's `type` and ",
         "`hidden_divisions` against the spike's `by_type`.", call. = FALSE)
  }

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
