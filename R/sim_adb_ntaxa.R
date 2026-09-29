# NOTE (for later consideration): only a single sampling_prob is supported here.
# Options for type-dependent sampling (sampling_prob per type) with exactly ntaxa sampled tips:
# (a) mark each particle as sampled at birth with probability sampling_prob[type] (types do
#     not change during a lifetime) and stop the simulation when ntaxa sampled
#     particles are alive; exact analogue of the ntaxa stopping rule, but requires
#     changes to the C++ loop (and changes the behaviour for a single sampling_prob)
# (b) simulate until ceiling(ntaxa / min(sampling_prob)) particles are alive, then keep exactly
#     ntaxa tips drawn with weights sampling_prob[type]; simple, but the type composition of
#     the sample only approximately follows sampling_prob
# (c) sample tips independently with probability sampling_prob[type] and re-simulate until
#     exactly ntaxa tips are sampled; exact, but may need many attempts
#' Simulator of a phylogeny from an Age-Dependent Branching process for a fixed number of sampled particles
#' @param ntaxa number of sampled particles - tips in the phylogeny (at least 2)
#' @param scale vector of scale parameters per type (alternatively, give `mean_lifetime`)
#' @param shape vector of shape parameters per type
#' @param mean_lifetime vector of mean lifetimes per type (= scale * shape),
#'   alternative to `scale`
#' @param death_prob vector of death probabilities per type
#' @param sampling_prob sampling probability
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param asym_trans_prob matrix of asymmetric type transition probabilities
#' @param sym_trans_prob matrix of symmetric type transition probabilities
#' @param collapse if TRUE, nodes with a single descendant after pruning are
#'   removed (node-typed tree), otherwise they are kept (branch-typed tree)
#' @return a treedata object with `ntaxa` tips, or `NULL` (with a message)
#'   if too many particles died
#' @export
sim_adb_ntaxa_samp <- function(ntaxa, scale = NULL, shape, death_prob = 0, sampling_prob = 1, origin_type = 0,
                               asym_trans_prob = matrix(0), sym_trans_prob = matrix(1), collapse = TRUE,
                               mean_lifetime = NULL) {
  scale <- .resolve_scale(scale, shape, mean_lifetime)
  # the tree parameters are validated by sim_adb_ntaxa_complete_fast
  if (length(sampling_prob) != 1 || !is.numeric(sampling_prob) || is.na(sampling_prob) ||
      sampling_prob <= 0 || sampling_prob > 1) {
    stop("`sampling_prob` must be a single sampling probability in (0, 1].", call. = FALSE)
  }

  # estimate the number of taxa in the full tree (the C++ loop takes an integer)
  nfull <- ceiling(ntaxa / sampling_prob)

  # simulate full tree
  tree <- sim_adb_ntaxa_complete_fast(ntaxa = nfull, scale = scale, shape = shape, death_prob = death_prob,
                                      origin_type = origin_type,
                                      asym_trans_prob = asym_trans_prob, sym_trans_prob = sym_trans_prob)
  if (is.null(tree)) {
    return(NULL)
  }

  prune_tree(obj = tree, ntips = ntaxa, collapse = collapse)
}


# NOTE (for later consideration): the tree is truncated at stopping_time =
# min(height of alive particles), i.e. at the NEXT event after ntaxa particles
# are alive. This time depends on the residual lifetimes of the alive particles,
# so their pendant branches are systematically too long (by about one waiting
# time between events; strongest for small ntaxa). Truncating at the LAST
# division (when the ntaxa-th particle is born) would avoid this: a one-sample
# log-rank test of the lifetimes rejects the current truncation (p < 0.01 at
# ntaxa = 20 and 100) but not the last-division truncation.
# See also: Hartmann, Wong & Stadler (2010) on sampling trees with n taxa (SSA vs GSA).
#' Simulator of the complete Age-Dependent Branching Process (up to a fixed number of living particles)
#' @param ntaxa number of living particles at which the process is stopped (at least 2)
#' @param scale vector of scale parameters per type (alternatively, give `mean_lifetime`)
#' @param shape vector of shape parameters per type
#' @param mean_lifetime vector of mean lifetimes per type (= scale * shape),
#'   alternative to `scale`
#' @param death_prob vector of death probabilities per type
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param asym_trans_prob matrix of asymmetric type transition probabilities
#' @param sym_trans_prob matrix of symmetric type transition probabilities
#' @return a treedata object including dead particles, with `@data` columns
#'   `node`, `status` (0 = dead, 1 = alive, 2 = divided) and `type`, or
#'   `NULL` (with a message) if too many particles died
#' @export
sim_adb_ntaxa_complete_fast <- function(ntaxa, scale = NULL, shape, death_prob, origin_type = 0,
                                        asym_trans_prob = matrix(0), sym_trans_prob = matrix(1),
                                        mean_lifetime = NULL) {
  scale <- .resolve_scale(scale, shape, mean_lifetime)
  .check_adb_params(scale, shape, death_prob, origin_type, asym_trans_prob, sym_trans_prob)
  raw <- sim_adb_loop_cpp(ntaxa, scale, shape, death_prob, asym_trans_prob, sym_trans_prob, origin_type)
  nodes <- as.data.frame(raw[c("id", "height", "type", "parent", "status")])

  if (sum(nodes$status == 1) < ntaxa) {
    message("Too many particles died! Try another seed or decrease death probability.")
    return(NULL)
  }

  # truncate to the minimal height of living particles
  stopping_time <- min(nodes$height[nodes$status == 1])
  if (any(nodes$status == 2)) {
    # ensure that truncation occurs after any cell division
    stopifnot(stopping_time > max(nodes$height[nodes$status == 2]))
  }
  nodes$height[nodes$status == 1] <- stopping_time

  # branch lengths post-hoc (heights: time since origin)
  nodes$edge_length <- nodes$height - nodes$height[match(nodes$parent, nodes$id)]
  .nodes_to_treedata(nodes, root_edge = raw$root_edge, origin = stopping_time)
}


# R loop (slower than Rcpp), kept internally as a reference implementation of
# sim_adb_ntaxa_complete_fast
sim_adb_ntaxa_complete <- function(ntaxa, scale = NULL, shape, death_prob, origin_type = 0,
                                   asym_trans_prob = matrix(0), sym_trans_prob = matrix(1),
                                   mean_lifetime = NULL) {
  scale <- .resolve_scale(scale, shape, mean_lifetime)
  .check_adb_params(scale, shape, death_prob, origin_type, asym_trans_prob, sym_trans_prob)

  # initialize
  nodes <- data.frame(
    id = integer(0),
    height = numeric(0),
    type = integer(0),
    parent = integer(0),
    leftchild = integer(0),
    rightchild = integer(0),
    status = integer(0) # 0 = dead, 1 = alive, 2 = divided
  )

  # sample the lifetime of the first particle
  root_edge <- rgamma(1, shape = shape[origin_type + 1], scale = scale[origin_type + 1])
  nodes <- dplyr::bind_rows(nodes, c(id = 1, height = root_edge, type = origin_type,
                                     parent = NA, leftchild = NA, rightchild = NA, status = 1))
  events <- nodes
  event_counter <- 1
  living <- 1

  while (living < ntaxa && nrow(events) > 0) {
    if (living %% 10000 == 0) {
      message(living, " cells alive")
    }
    # look at one event
    event <- as.list(events[1, ])
    events <- events[-1, ]

    if (runif(1) < death_prob[event$type + 1]) {
      # particle dies
      nodes[event$id, "status"] <- 0
      living <- living - 1
    } else {
      # particle divides, create two new children
      left_id <- event_counter + 1
      right_id <- event_counter + 2
      event_counter <- event_counter + 2
      nodes[event$id, "status"] <- 2
      living <- living + 1

      if (ncol(sym_trans_prob) == 1) {
        # single-type case
        children_types <- rep(origin_type, 2)
      } else {
        # multi-type case: sample types
        children_types <- sample_types(parent_type = event$type, asym_trans_prob = asym_trans_prob,
                                       sym_trans_prob = sym_trans_prob)
      }

      # sample lifetimes and add new nodes
      left_type <- children_types[1]
      left_lifetime <- rgamma(1, shape = shape[left_type + 1], scale = scale[left_type + 1])
      left_node <- c(id = left_id, height = event$height + left_lifetime, type = left_type,
                     parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      nodes <- dplyr::bind_rows(nodes, left_node)

      right_type <- children_types[2]
      right_lifetime <- rgamma(1, shape = shape[right_type + 1], scale = scale[right_type + 1])
      right_node <- c(id = right_id, height = event$height + right_lifetime, type = right_type,
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      nodes <- dplyr::bind_rows(nodes, right_node)

      # to be processed, sort by height
      events <- dplyr::bind_rows(events, left_node, right_node)
      events <- events[order(events$height), ]

      # add child relationships
      nodes[event$id, c("leftchild", "rightchild")] <- c(left_id, right_id)
    }
  }

  if (sum(nodes$status == 1) < ntaxa) {
    message("Too many particles died! Try another seed or decrease death probability.")
    return(NULL)
  }

  # truncate to the minimal height of living particles
  stopping_time <- min(nodes$height[nodes$status == 1])
  if (any(nodes$status == 2)) {
    # ensure that truncation occurs after any cell division
    stopifnot(stopping_time > max(nodes$height[nodes$status == 2]))
  }
  nodes$height[nodes$status == 1] <- stopping_time

  # branch lengths post-hoc (heights: time since origin)
  nodes$edge_length <- nodes$height - nodes$height[match(nodes$parent, nodes$id)]
  .nodes_to_treedata(nodes, root_edge = root_edge, origin = stopping_time)
}
