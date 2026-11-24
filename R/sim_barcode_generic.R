# This simulator of lineage barcodes extends the TiDeTree and SciPhy models, i.e.,
# follows a generic editing process with time- and type-dependent editing and silencing rates and dropout probabilties. 
# sim_barcode_nonseq and sim_barcode_seq are special cases.

#' Simulator of generic lineage barcodes
#'
#' @param tree a treedata object of a sampled tree: ultrametric, i.e. all tips
#'   are sampled cells at present (e.g. from [sim_adb_origin_samp()])
#' @param n_barcodes number of barcodes per cell
#' @param n_sites number of target sites per barcode
#' @param edit_rate scalar or length-`n_barcodes` editing rate(s) per barcode
#'   (a number is constant in time). Time- and type-dependent rates are described by
#'   [barcode_rate()], with time profiles from [rate_profiles].
#' @param edit_probs numeric vector of length E, summing to 1: the relative
#'   frequencies of the edit outcomes. If named, the names are the labels of the
#'   outcomes written to the barcodes (e.g. `c(A = 0.7, T = 0.3)`); otherwise
#'   the outcomes are labelled 1,...,E.
#' @param silencing_rate scalar or length-`n_barcodes` silencing rate(s) (or [barcode_rate()]), racing
#'   against the editing process: the rate at which the entire barcode is silenced
#'   if `sequential = TRUE`, and the rate of each site if `sequential = FALSE`
#' @param dropout_prob scalar or length-`n_barcodes` probability (or [barcode_rate()] with `by_type` only,
#'   which may be `NA` for types not occurring at tips) that a
#'   barcode (if `sequential = TRUE`) or each site (if `sequential = FALSE`) drops out at a tip
#' @param missing_state string marking silenced/dropout barcodes (default "-"),
#'   different from the unedited state "0" and the outcome labels
#' @param sequential if `TRUE` (default), sites are edited left to right (only the leftmost
#'   unedited site is edited, at rate `edit_rate`), and silencing and dropout affect the
#'   whole barcode (SciPhy-like). If `FALSE`, every unedited site is edited independently at rate
#'   `edit_rate`, and silencing and dropout act on each site independently (TiDeTree-like).
#' @param write_as_string if `TRUE` (default), each barcode is returned as one
#'   "_"-separated string per node; if `FALSE`, as a node x site data frame
#'
#' @return if `write_as_string = TRUE`, a data frame with one row per node in the tree and one column
#'   per barcode: \code{node}, \code{barcode_1}, ..., \code{barcode_<n_barcodes>},
#'   each entry a string of `n_sites` positions separated by "_"
#'   (e.g. "0_A_B_0_0"), where "0" = unedited and `missing_state` = silenced/dropout
#'   (every site reads `missing_state` once a barcode is silenced or dropped out).
#'   If `write_as_string = FALSE`, a named list (\code{barcode_1}, ...) with one data frame per barcode,
#'   with one row per node and (character) columns \code{node}, \code{site_1}, ..., \code{site_<n_sites>}.
#' @family barcode simulators
#' @seealso [barcode_rate()] and [rate_profiles] for time- and type-dependent rates
#' @export
sim_barcode_generic <- function(tree, n_barcodes, n_sites, edit_rate, edit_probs,
                                silencing_rate = 0, dropout_prob = 0, missing_state = "-", write_as_string = TRUE,
                                sequential = TRUE) {

  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }

  stopifnot(
    is.logical(write_as_string), length(write_as_string) == 1, !is.na(write_as_string),
    is.logical(sequential), length(sequential) == 1, !is.na(sequential),
    length(n_barcodes) == 1, n_barcodes >= 1, length(n_sites) == 1, n_sites >= 1
  )

  # rates: a number, a per-barcode vector, or a barcode_rate() with a time profile
  edit_rate <- .as_barcode_rate(edit_rate, n_barcodes, "edit_rate")
  silencing_rate <- .as_barcode_rate(silencing_rate, n_barcodes, "silencing_rate")
  dropout_prob <- .as_barcode_rate(dropout_prob, n_barcodes, "dropout_prob")
  if (dropout_prob$time_given) {
    stop("`dropout_prob` is a probability at the tips and has no time profile.", call. = FALSE)
  }
  if (any(dropout_prob$base > 1)) stop("`dropout_prob` must not exceed 1.", call. = FALSE)

  if (!is.numeric(edit_probs) || length(edit_probs) < 1 || anyNA(edit_probs) || any(edit_probs < 0) ||
      !isTRUE(all.equal(sum(edit_probs), 1))) {
    stop("`edit_probs` must be non-negative numbers summing to 1.", call. = FALSE)
  }
  # outcome labels: names of edit_probs, otherwise 1..E
  labels <- if (is.null(names(edit_probs))) as.character(seq_along(edit_probs)) else names(edit_probs)
  if (anyNA(labels) || any(labels == "") || anyDuplicated(labels) || any(grepl("_", labels, fixed = TRUE))) {
    stop("Names of `edit_probs` must be unique, non-empty and must not contain \"_\".", call. = FALSE)
  }
  if (length(missing_state) != 1 || is.na(missing_state) ||
      as.character(missing_state) %in% c("0", labels) || grepl("_", missing_state, fixed = TRUE)) {
    stop("`missing_state` must be a single string different from \"0\" and the edit outcomes, without \"_\".",
         call. = FALSE)
  }
  missing_state <- as.character(missing_state)

  times <- .barcode_tree_times(tree)

  # type of each node (by node number), if any rate depends on it; the branch ending in a node
  # has that node's type (the root edge: the root's type)
  # tips are nodes 1..Ntip (all sampled cells, since the tree is ultrametric)
  node_type <- .check_node_types(.node_data(tree), edit_rate, silencing_rate, dropout_prob, times$tips)
  if (!is.null(node_type) && !is.null(dropout_prob$by_type) &&
      any(outer(dropout_prob$base, dropout_prob$by_type[unique(node_type[times$tips]) + 1]) > 1)) {
    stop("`dropout_prob` times `by_type` must not exceed 1.", call. = FALSE)
  }

  codes <- .sim_barcode_codes(times, n_barcodes, n_sites, edit_rate, silencing_rate, dropout_prob,
                              edit_probs, sequential, unname(node_type))
  .format_barcodes(codes, c("0", labels, missing_state), write_as_string)
}


# Nodes 1..N of a treedata object with their `type` (NA if missing), if the tree has types
.node_data <- function(tree) {
  node <- seq_len(length(tree@phylo$tip.label) + tree@phylo$Nnode)
  if (!"type" %in% names(tree@data)) return(data.frame(node = node))
  data.frame(node = node, type = tree@data$type[match(node, tree@data$node)])
}


# Branches of a tree for the barcode simulators: the time since the origin at the start (t0) and end (t1)
# of the branch ending in each node, the branches in preorder (parent before child, by the topology
# rather than by heights, which tie on zero-length branches) and the tips.
# The root edge is the branch ending in the root; the origin is the tree height plus the root edge.
.barcode_tree_times <- function(tree) {
  phy <- tree@phylo
  root <- length(phy$tip.label) + 1L
  root_edge <- if (is.null(phy$root.edge)) 0 else phy$root.edge
  t1 <- ape::node.depth.edgelength(phy) + root_edge
  origin <- max(t1)
  if (!is.null(phy$origin) && !isTRUE(all.equal(phy$origin, origin))) {
    warning("`tree@phylo$origin` differs from the tree height plus root edge; the latter is used.", call. = FALSE)
  }
  edge <- ape::reorder.phylo(phy, order = "cladewise")$edge
  t0 <- numeric(length(t1))
  t0[edge[, 2]] <- t1[edge[, 1]]
  t0[root] <- min(0, root_edge)
  list(t0 = t0, t1 = t1, origin = origin, parent = c(0L, edge[, 1]), child = c(root, edge[, 2]),
       tips = seq_along(phy$tip.label))
}


# Simulates the barcodes of every node with the C++ engine (src/sim_barcodes.cpp), shared by all
# barcode simulators. Rates are barcode_rate() objects with one base rate per barcode; node_type is
# the type of each node (NULL: no type-dependent rates).
# Returns an integer array nodes x sites x barcodes of states:
# 0 = unedited, 1..E = edit outcome, E + 1 = silenced or dropped out.
.sim_barcode_codes <- function(times, n_barcodes, n_sites, edit_rate, silencing_rate, dropout_prob,
                               edit_probs, sequential, node_type = NULL) {
  multiplier <- function(rate, nodes) {
    if (is.null(rate$by_type)) rep(1, length(nodes)) else rate$by_type[node_type[nodes] + 1]
  }
  nodes <- seq_along(times$t1)
  sim_barcodes_cpp(
    times$parent, times$child, times$tips, as.integer(n_barcodes), as.integer(n_sites), sequential,
    edit_rate$base, .rate_process(edit_rate, times$t0, times$t1, multiplier(edit_rate, nodes)),
    silencing_rate$base, .rate_process(silencing_rate, times$t0, times$t1, multiplier(silencing_rate, nodes)),
    dropout_prob$base, multiplier(dropout_prob, times$tips), as.numeric(edit_probs)
  )
}


# Output of sim_barcode_generic() and sim_barcode_seq() from the state array of .sim_barcode_codes():
# `states` are the labels of the states 0, 1..E, E + 1. If `integer_missing` is a number, the
# data frames are integer, with the missing state relabelled to it.
.format_barcodes <- function(codes, states, write_as_string, integer_missing = NULL) {
  n_nodes <- dim(codes)[1]
  n_sites <- dim(codes)[2]
  n_barcodes <- dim(codes)[3]
  barcode_names <- paste0("barcode_", seq_len(n_barcodes))
  if (write_as_string) {
    strings <- lapply(seq_len(n_barcodes), function(bc) {
      do.call(paste, c(lapply(seq_len(n_sites), function(s) states[codes[, s, bc] + 1L]), sep = "_"))
    })
    names(strings) <- barcode_names
    return(data.frame(node = seq_len(n_nodes), strings))
  }
  out <- lapply(seq_len(n_barcodes), function(bc) {
    if (is.null(integer_missing)) {
      m <- matrix(states[codes[, , bc] + 1L], n_nodes, n_sites)
    } else {
      m <- codes[, , bc, drop = FALSE]
      dim(m) <- c(n_nodes, n_sites)
      m[m == length(states) - 1L] <- as.integer(integer_missing)
    }
    colnames(m) <- paste0("site_", seq_len(n_sites))
    data.frame(node = seq_len(n_nodes), m, stringsAsFactors = FALSE)
  })
  names(out) <- barcode_names
  out
}
