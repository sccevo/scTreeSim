#' Simulator of a phylogeny from an Age-Dependent Branching Process for a fixed time interval (since origin)
#' @param origin_time time of birth of the initial particle
#' @param scale vector of scale parameters per type (alternatively, give `mean_lifetime`)
#' @param shape vector of shape parameters per type
#' @param mean_lifetime vector of mean lifetimes per type (= scale * shape),
#'   alternative to `scale`
#' @param death_prob vector of death probabilities per type
#' @param sampling_prob sampling probability at present: a single value for all types,
#'   or a vector with one value per type (values of 0 allowed, but not all)
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param asym_trans_prob matrix of asymmetric type transition probabilities
#' @param sym_trans_prob matrix of symmetric type transition probabilities
#' @param min_taxa minimum number of tips in the phylogeny
#' @param collapse if TRUE, nodes with a single descendant after pruning are
#'   removed (node-typed tree), otherwise they are kept (branch-typed tree)
#' @return a treedata object, or `NULL` (with a message) if fewer than
#'   `min_taxa` tips are alive or sampled
#' @export
sim_adb_origin_samp <- function(origin_time, scale = NULL, shape, death_prob = 0, sampling_prob = 1, origin_type = 0,
                                asym_trans_prob = matrix(0), sym_trans_prob = matrix(1), min_taxa = 2, collapse = TRUE,
                                mean_lifetime = NULL) {
  scale <- .resolve_scale(scale, shape, mean_lifetime)
  # the tree parameters are validated by sim_adb_origin_complete_fast
  if (!is.numeric(sampling_prob) || !(length(sampling_prob) %in% c(1, length(shape))) || anyNA(sampling_prob) ||
      any(sampling_prob < 0 | sampling_prob > 1) || all(sampling_prob == 0)) {
    stop("`sampling_prob` must be a single sampling probability in (0, 1], or one sampling probability ",
         "in [0, 1] per type (", length(shape), ", not all 0).", call. = FALSE)
  }

  # simulate full tree
  tree <- sim_adb_origin_complete_fast(origin_time = origin_time, scale = scale, shape = shape, death_prob = death_prob,
                                       origin_type = origin_type, asym_trans_prob = asym_trans_prob,
                                       sym_trans_prob = sym_trans_prob, min_taxa = min_taxa)
  if (is.null(tree)) {
    return(NULL)
  }

  prune_tree(obj = tree, sampling_prob = sampling_prob, min_taxa = min_taxa, collapse = collapse)
}


#' R wrapper for a faster C++ implementation of complete Age-Dependent Branching Process
#' for a fixed time interval (since origin)
#' @param origin_time time of birth of the initial particle
#' @param scale vector of scale parameters per type (alternatively, give `mean_lifetime`)
#' @param shape vector of shape parameters per type
#' @param mean_lifetime vector of mean lifetimes per type (= scale * shape),
#'   alternative to `scale`
#' @param death_prob vector of death probabilities per type
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param asym_trans_prob matrix of asymmetric type transition probabilities
#' @param sym_trans_prob matrix of symmetric type transition probabilities
#' @param min_taxa minimum number of tips in the tree
#' @return a treedata object including dead particles, with `@data` columns
#'   `node`, `status` (0 = dead, 1 = alive, 2 = divided) and `type`, or
#'   `NULL` (with a message) if fewer than `min_taxa` particles are alive
#' @export
sim_adb_origin_complete_fast <- function(origin_time, scale = NULL, shape, death_prob, origin_type = 0,
                                         asym_trans_prob = matrix(0), sym_trans_prob = matrix(1), min_taxa = 2,
                                         mean_lifetime = NULL) {
  scale <- .resolve_scale(scale, shape, mean_lifetime)
  .check_adb_params(scale, shape, death_prob, origin_type, asym_trans_prob, sym_trans_prob)
  raw <- sim_adb_origin_loop_cpp(origin_time, scale, shape, death_prob, asym_trans_prob, sym_trans_prob, origin_type)
  nodes <- as.data.frame(raw[c("id", "height", "type", "parent", "status")])

  # check number of tips
  if (sum(nodes$status == 1) < min_taxa) {
    message("The simulated tree has too few tips. Try another seed.")
    return(NULL)
  }

  # branch lengths post-hoc (heights: time before present)
  nodes$edge_length <- nodes$height[match(nodes$parent, nodes$id)] - nodes$height
  .nodes_to_treedata(nodes, root_edge = raw$root_edge, origin = origin_time)
}


# R loop (slower than Rcpp), kept internally as a reference implementation of
# sim_adb_origin_complete_fast
sim_adb_origin_complete <- function(origin_time, scale = NULL, shape, death_prob, origin_type = 0,
                                    asym_trans_prob = matrix(0), sym_trans_prob = matrix(1), min_taxa = 2,
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
  # censor the root lifetime: if it outlives the origin interval, the root is
  # a single tip alive at present and is not processed as an event
  root_censored <- root_edge > origin_time
  if (root_censored) root_edge <- origin_time
  nodes <- dplyr::bind_rows(nodes, c(id = 1, height = origin_time - root_edge, type = origin_type,
                                     parent = NA, leftchild = NA, rightchild = NA, status = 1))
  events <- if (root_censored) nodes[0, ] else nodes
  event_counter <- 1

  while (nrow(events) > 0) {
    # look at one event
    event <- as.list(events[1, ])
    events <- events[-1, ]

    if (runif(1) < death_prob[event$type + 1]) {
      # particle dies
      nodes[event$id, "status"] <- 0
    } else {
      # particle divides, create two new children
      left_id <- event_counter + 1
      right_id <- event_counter + 2
      event_counter <- event_counter + 2
      nodes[event$id, "status"] <- 2

      if (ncol(sym_trans_prob) == 1) {
        # single-type case
        children_types <- rep(origin_type, 2)
      } else {
        # multi-type case: sample types
        children_types <- sample_types(parent_type = event$type, asym_trans_prob = asym_trans_prob,
                                       sym_trans_prob = sym_trans_prob)
      }

      # sample lifetimes and add new nodes; children that outlive the present
      # are censored at height 0 (tips), all others are processed later
      left_type <- children_types[1]
      left_lifetime <- rgamma(1, shape = shape[left_type + 1], scale = scale[left_type + 1])
      left_censored <- event$height - left_lifetime < 0
      if (left_censored) left_lifetime <- event$height
      left_node <- c(id = left_id, height = event$height - left_lifetime, type = left_type,
                     parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      if (!left_censored) events <- dplyr::bind_rows(events, left_node)

      right_type <- children_types[2]
      right_lifetime <- rgamma(1, shape = shape[right_type + 1], scale = scale[right_type + 1])
      right_censored <- event$height - right_lifetime < 0
      if (right_censored) right_lifetime <- event$height
      right_node <- c(id = right_id, height = event$height - right_lifetime, type = right_type,
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      if (!right_censored) events <- dplyr::bind_rows(events, right_node)

      # add child relationships
      nodes <- dplyr::bind_rows(nodes, left_node, right_node)
      nodes[event$id, c("leftchild", "rightchild")] <- c(left_id, right_id)
    }
  }

  # check number of tips
  if (sum(nodes$status == 1) < min_taxa) {
    message("The simulated tree has too few tips. Try another seed.")
    return(NULL)
  }

  # branch lengths post-hoc (heights: time before present)
  nodes$edge_length <- nodes$height[match(nodes$parent, nodes$id)] - nodes$height
  .nodes_to_treedata(nodes, root_edge = root_edge, origin = origin_time)
}
