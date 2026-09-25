#' Simulator of a phylogeny from an Age-Dependent Branching process for a fixed number of sampled particles
#' @param ntaxa number of sampled particles - tips in the phylogeny (at least 2)
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param rho sampling probability
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param Xi_as matrix of asymmetric type transition probabilities
#' @param Xi_s matrix of symmetric type transition probabilities
#' @param collapse if TRUE, nodes with a single descendant after pruning are
#'   removed (node-typed tree), otherwise they are kept (branch-typed tree)
#' @return a treedata object with `ntaxa` tips, or `NULL` (with a message)
#'   if too many particles died
#' @export
sim_adb_ntaxa_samp <- function(ntaxa, a, b, d = 0, rho = 1, origin_type = 0,
                               Xi_as = matrix(0), Xi_s = matrix(1), collapse = TRUE) {
  # the tree parameters are validated by sim_adb_ntaxa_complete_fast
  if (length(rho) != 1 || !is.numeric(rho) || rho <= 0 || rho > 1) {
    stop("`rho` must be a single sampling probability in (0, 1].", call. = FALSE)
  }

  # estimate the number of taxa in the full tree (the C++ loop takes an integer)
  nfull <- ceiling(ntaxa / rho)

  # simulate full tree
  tree <- sim_adb_ntaxa_complete_fast(ntaxa = nfull, a = a, b = b, d = d, origin_type = origin_type,
                                      Xi_as = Xi_as, Xi_s = Xi_s)
  if (is.null(tree)) {
    return(NULL)
  }

  prune_tree(obj = tree, ntips = ntaxa, collapse = collapse)
}


#' Simulator of the complete Age-Dependent Branching Process (up to a fixed number of living particles)
#' @param ntaxa number of living particles at which the process is stopped (at least 2)
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param Xi_as matrix of asymmetric type transition probabilities
#' @param Xi_s matrix of symmetric type transition probabilities
#' @return a treedata object including dead particles, with `@data` columns
#'   `node`, `status` (0 = dead, 1 = alive, 2 = divided) and `type`, or
#'   `NULL` (with a message) if too many particles died
#' @export
sim_adb_ntaxa_complete_fast <- function(ntaxa, a, b, d, origin_type = 0,
                                        Xi_as = matrix(0), Xi_s = matrix(1)) {
  .check_adb_params(a, b, d, origin_type, Xi_as, Xi_s)
  raw <- sim_adb_loop_cpp(ntaxa, a, b, d, Xi_as, Xi_s, origin_type)
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
sim_adb_ntaxa_complete <- function(ntaxa, a, b, d, origin_type = 0,
                                   Xi_as = matrix(0), Xi_s = matrix(1)) {
  .check_adb_params(a, b, d, origin_type, Xi_as, Xi_s)

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
  root_edge <- rgamma(1, shape = b[origin_type + 1], scale = a[origin_type + 1])
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

    if (runif(1) < d[event$type + 1]) {
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

      if (ncol(Xi_s) == 1) {
        # single-type case
        children_types <- rep(origin_type, 2)
      } else {
        # multi-type case: sample types
        children_types <- sample_types(parent_type = event$type, Xi_as = Xi_as, Xi_s = Xi_s)
      }

      # sample lifetimes and add new nodes
      left_type <- children_types[1]
      left_lifetime <- rgamma(1, shape = b[left_type + 1], scale = a[left_type + 1])
      left_node <- c(id = left_id, height = event$height + left_lifetime, type = left_type,
                     parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      nodes <- dplyr::bind_rows(nodes, left_node)

      right_type <- children_types[2]
      right_lifetime <- rgamma(1, shape = b[right_type + 1], scale = a[right_type + 1])
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
