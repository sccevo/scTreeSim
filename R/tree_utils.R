#' Prune tree
#'
#' Removes dead particles and then samples tips at present, either with
#' probability `rho` each or as a fixed number `ntips`.
#' @param obj treedata object
#' @param rho sampling probability (at present only)
#' @param ntips number of sampled tips
#' @param min_tips minimum number of tips in the pruned tree
#' @param collapse if TRUE, stores the standard structure (node-typed tree), otherwise keeps nodes with one descendant (branch-typed tree)
#' @return the pruned treedata object, or `NULL` (with a message) if too few tips remain
#' @keywords internal
prune_tree <- function(obj, rho = NA, ntips = NA, min_tips = 2, collapse = TRUE) {
  if (is.na(rho) == is.na(ntips)) {
    stop("Please provide either the sampling probability or the desired number of tips in pruned tree.", call. = FALSE)
  }

  # store origin
  origin <- obj@phylo$origin

  # prune dead particles (dead particles are always tips)
  dead_nodes <- obj@data$node[obj@data$status == 0]
  # treeio::drop.tip emits "Invalid edge matrix for <phylo>. A <tbl_df> is returned." --> suppressed
  obj <- suppressMessages(treeio::drop.tip(obj, obj@phylo$tip.label[dead_nodes], collapse.singles = collapse))

  if (is.null(obj)) {
    message("All particles died! Try another seed.")
    return(NULL)
  }

  # prune unsampled particles
  tips <- obj@phylo$tip.label
  if (!is.na(rho)) {
    sampled_tips <- tips[sample(c(TRUE, FALSE), length(tips), prob = c(rho, 1 - rho), replace = TRUE)]
  } else {
    sampled_tips <- sample(tips, ntips)
  }

  if (length(sampled_tips) < min_tips) {
    message("Not enough tips sampled! Try another seed.")
    return(NULL)
  }

  obj <- suppressMessages(treeio::drop.tip(obj, setdiff(tips, sampled_tips), collapse.singles = collapse))

  # add origin and re-calculate root.edge
  obj@phylo$root.edge <- origin - max(ape::node.depth.edgelength(obj@phylo))
  obj@phylo$origin <- origin

  obj
}


#' Assemble a treedata object from the particles of a simulated ADB process
#'
#' Shared by all tree simulators (R and C++ loops). Nodes are numbered as
#' in ape: alive tips 1..n_alive, then dead tips, then internal nodes
#' (divided particles). Tip labels equal the tip node numbers.
#' @param nodes data frame with one row per particle and columns `id`,
#'   `parent` (NA for the root), `type`, `status` (0 = dead, 1 = alive,
#'   2 = divided) and `edge_length` (length of the branch leading to the
#'   particle; ignored for the root)
#' @param root_edge length of the root edge (lifetime of the first particle)
#' @param origin time from the origin to the present (stored in `phylo$origin`)
#' @return a treedata object with `@data` columns `node`, `status`, `type`
#' @noRd
.nodes_to_treedata <- function(nodes, root_edge, origin) {
  is_alive <- nodes$status == 1
  is_dead <- nodes$status == 0
  is_divided <- nodes$status == 2
  n_alive <- sum(is_alive)
  n_dead <- sum(is_dead)
  n_tip <- n_alive + n_dead
  n_node <- sum(is_divided)

  # node numbers in ape convention
  label <- integer(nrow(nodes))
  label[is_alive] <- seq_len(n_alive)
  label[is_dead] <- n_alive + seq_len(n_dead)
  label[is_divided] <- n_tip + seq_len(n_node)

  # one edge per non-root particle: parent -> particle
  child_rows <- which(!is.na(nodes$parent))
  edge <- cbind(label[match(nodes$parent[child_rows], nodes$id)], label[child_rows])

  phylo_tree <- list(edge = edge, edge.length = nodes$edge_length[child_rows],
                     Nnode = n_node, tip.label = as.character(seq_len(n_tip)))
  class(phylo_tree) <- "phylo"

  tree <- treeio::as.treedata(phylo_tree)
  tree@phylo$root.edge <- root_edge
  tree@phylo$origin <- origin
  node_order <- order(label)
  tree@data <- tibble::tibble(
    node = label[node_order],
    status = as.integer(nodes$status[node_order]),
    type = as.integer(nodes$type[node_order])
  )
  tree
}


#' Validate the parameters of the Age-Dependent Branching process
#'
#' Called by every tree simulator before any sampling, so that invalid
#' indices never reach the C++ code.
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param Xi_as matrix of asymmetric type transition probabilities
#' @param Xi_s matrix of symmetric type transition probabilities
#' @return `NULL`, invisibly; throws an error if any parameter is invalid
#' @noRd
.check_adb_params <- function(a, b, d, origin_type, Xi_as, Xi_s) {
  ntypes <- length(a)
  if (ntypes < 1 || !is.numeric(a) || any(a <= 0)) {
    stop("`a` must be a non-empty numeric vector of positive scale parameters.", call. = FALSE)
  }
  if (length(b) != ntypes || !is.numeric(b) || any(b <= 0)) {
    stop("`b` must be a numeric vector of positive shape parameters with one entry per type (", ntypes, ").", call. = FALSE)
  }
  if (length(d) != ntypes || !is.numeric(d) || any(d < 0 | d >= 1)) {
    stop("`d` must be a numeric vector of death probabilities in [0, 1) with one entry per type (", ntypes, ").", call. = FALSE)
  }
  if (length(origin_type) != 1 || !(origin_type %in% (seq_len(ntypes) - 1L))) {
    stop("`origin_type` must be a single value in 0, ..., ", ntypes - 1, ".", call. = FALSE)
  }
  if (!is.matrix(Xi_as) || !is.matrix(Xi_s) ||
      !identical(dim(Xi_as), c(ntypes, ntypes)) || !identical(dim(Xi_s), c(ntypes, ntypes))) {
    stop("`Xi_as` and `Xi_s` must both be ", ntypes, " x ", ntypes, " matrices.", call. = FALSE)
  }
  if (any(Xi_as < 0) || any(Xi_s < 0)) {
    stop("`Xi_as` and `Xi_s` must not contain negative probabilities.", call. = FALSE)
  }
  row_sums <- rowSums(Xi_as) + rowSums(Xi_s)
  if (!isTRUE(all.equal(row_sums, rep(1, ntypes), check.attributes = FALSE))) {
    stop("The transition probabilities in each row of `Xi_as` + `Xi_s` must sum to 1.", call. = FALSE)
  }
  invisible(NULL)
}


#' Helper function for sampling types of offspring upon division
#' @param parent_type ID of parent type
#' @param Xi_as matrix of asymmetric type transition probabilities
#' @param Xi_s matrix of symmetric type transition probabilities
#' @return integer vector of length 2 with the types of the two children
#' @keywords internal
#' @noRd
sample_types <- function(parent_type, Xi_as, Xi_s) {
  # mirrors sample_child_types in the C++ code
  ntypes <- ncol(Xi_s)
  row <- parent_type + 1

  # total probability of the row, accumulated in the same order as the loop
  # below, so that the scaled draw r always falls below the final cum_prob:
  # outcomes are only ever chosen within the support (no fallback needed,
  # even in case of rounding errors when summing probabilities)
  total <- 0
  for (i in seq_len(ntypes)) {
    total <- total + Xi_s[row, i]
    total <- total + Xi_as[row, i]
  }
  r <- runif(1) * total
  
  cum_prob <- 0
  for (i in seq(0, ntypes - 1)) {
    # symmetric: both children type i
    cum_prob <- cum_prob + Xi_s[row, i + 1]
    if (r < cum_prob) {
      return(c(i, i))
    }
    # asymmetric: one child stays parent_type, other becomes type i
    cum_prob <- cum_prob + Xi_as[row, i + 1]
    if (r < cum_prob) {
      if (runif(1) < 0.5) {
        return(c(parent_type, i))
      } else {
        return(c(i, parent_type))
      }
    }
  }
  # unreachable: r < total and the last positive-probability outcome brings
  # cum_prob to exactly total
  stop("Failed to sample child types; check the transition matrices.", call. = FALSE)
}
