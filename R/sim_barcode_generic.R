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
#' @param edit_probs numeric vector of length E, summing to 1: the relative
#'   frequencies of the edit outcomes. If named, the names are the labels of the
#'   outcomes written to the barcodes (e.g. `c(A = 0.7, T = 0.3)`); otherwise
#'   the outcomes are labelled 1,...,E.
#' @param silencing_rate scalar or length-`n_barcodes` rate(s) at which an entire
#'   barcode is silenced, racing against the editing process
#' @param dropout_prob scalar or length-`n_barcodes` probability that a
#'   barcode drops out at a tip
#' @param missing_state string marking silenced/dropout barcodes (default "-"),
#'   different from the unedited state "0" and the outcome labels
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
#' @export
sim_barcode_generic <- function(tree, n_barcodes, n_sites, edit_rate, edit_probs,
                                silencing_rate = 0, dropout_prob = 0, missing_state = "-", write_as_string = TRUE) {

  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }

  stopifnot(
    is.logical(write_as_string), length(write_as_string) == 1, !is.na(write_as_string),
    length(n_barcodes) == 1, n_barcodes >= 1, length(n_sites) == 1, n_sites >= 1,
    is.numeric(edit_rate), !anyNA(edit_rate), length(edit_rate) %in% c(1, n_barcodes), all(edit_rate >= 0),
    is.numeric(silencing_rate), !anyNA(silencing_rate), length(silencing_rate) %in% c(1, n_barcodes),
    all(silencing_rate >= 0),
    is.numeric(dropout_prob), !anyNA(dropout_prob), length(dropout_prob) %in% c(1, n_barcodes),
    all(dropout_prob >= 0 & dropout_prob <= 1)
  )

  # allow a single shared rate to stand in for all barcodes
  if (length(edit_rate) == 1) edit_rate <- rep(edit_rate, n_barcodes)
  if (length(silencing_rate) == 1) silencing_rate <- rep(silencing_rate, n_barcodes)
  if (length(dropout_prob) == 1) dropout_prob <- rep(dropout_prob, n_barcodes)

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

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  # visit parent before child so state can propagate down the tree; use the
  # topology rather than heights, which tie on zero-length branches
  child_order <- ape::reorder.phylo(tree@phylo, order = "cladewise")$edge[, 2]
  order_df <- tree_df[match(child_order, tree_df$node), ]

  # branch lengths come from the treedata; the root edge from the phylo object
  root_edge <- if (is.null(tree@phylo$root.edge)) 0 else tree@phylo$root.edge
  origin <- max(ape::node.depth.edgelength(tree@phylo)) + root_edge
  if (!is.null(tree@phylo$origin) && !isTRUE(all.equal(tree@phylo$origin, origin))) {
    warning("`tree@phylo$origin` differs from the tree height plus root edge; the latter is used.", call. = FALSE)
  }

  # tips are nodes 1..Ntip (all sampled cells, since the tree is ultrametric)
  tip_nodes <- seq_along(tree@phylo$tip.label)

  bc_states <- vector("list", n_barcodes)
  for (bc in seq_len(n_barcodes)) {
    # state_at: barcode state per node; silenced_at: whether that barcode
    # has already been silenced (needed since silencing must persist down
    # every descendant branch once it happens)
    state_at <- list()
    silenced_at <- list()
    state_at[[as.character(root)]] <- rep("0", n_sites)
    silenced_at[[as.character(root)]] <- FALSE

    # evolve along the root/origin edge first, if the tree has one
    if (root_edge > 0) {
      res <- .evolve_branch_generic(
        state_at[[as.character(root)]], silenced_at[[as.character(root)]],
        root_edge, edit_rate[bc], silencing_rate[bc], edit_probs, labels, missing_state
      )
      state_at[[as.character(root)]] <- res$state
      silenced_at[[as.character(root)]] <- res$silenced
    }

    # then walk every remaining branch, parent state -> child state
    for (i in seq_len(nrow(order_df))) {
      node <- order_df$node[i]
      if (node == root) next
      parent <- order_df$parent[i]

      res <- .evolve_branch_generic(
        state_at[[as.character(parent)]], silenced_at[[as.character(parent)]],
        order_df$branch.length[i], edit_rate[bc], silencing_rate[bc], edit_probs, labels, missing_state
      )
      state_at[[as.character(node)]] <- res$state
      silenced_at[[as.character(node)]] <- res$silenced
    }

    # apply dropout after simulation,
    # skipped for barcodes already fully silenced (nothing left to mask)
    if (dropout_prob[bc] > 0) {
      for (nd in tip_nodes) {
        key <- as.character(nd)
        if (!silenced_at[[key]] && stats::runif(1) < dropout_prob[bc]) {
          state_at[[key]][] <- missing_state
        }
      }
    }

    # nodes x sites character matrix
    # (explicit matrix(): vapply drops to a vector when n_sites = 1)
    bc_states[[bc]] <- matrix(
      unlist(lapply(tree_df$node, function(nd) state_at[[as.character(nd)]])),
      ncol = n_sites, byrow = TRUE
    )
  }

  if (write_as_string) {
    out <- data.frame(node = tree_df$node)
    for (bc in seq_len(n_barcodes)) {
      out[[paste0("barcode_", bc)]] <- apply(bc_states[[bc]], 1, paste, collapse = "_")
    }
    return(out)
  }

  out <- lapply(bc_states, function(m) {
    colnames(m) <- paste0("site_", seq_len(n_sites))
    data.frame(node = tree_df$node, m, stringsAsFactors = FALSE)
  })
  names(out) <- paste0("barcode_", seq_len(n_barcodes))
  out
}


# Evolve one barcode copy along one branch. At each event, draw whether it
# is an edit or silencing
.evolve_branch_generic <- function(parent_state, parent_silenced, branch_length,
                                   edit_rate, silencing_rate, edit_probs, labels, missing_state) {
  state <- parent_state

  # already silenced upstream, or no time on this branch: nothing to do
  if (parent_silenced || branch_length <= 0) {
    return(list(state = state, silenced = parent_silenced))
  }

  t <- 0

  repeat {
    # once every site is edited, only silencing can still happen
    saturated <- all(state != "0")
    current_edit_rate <- if (saturated) 0 else edit_rate
    event_rate <- current_edit_rate + silencing_rate

    if (event_rate <= 0) break  # no editing left to do and no silencing possible

    # draw time to the next event (edit or silencing) on the combined clock
    t <- t + stats::rexp(1, event_rate)
    if (t > branch_length) break  # no more events fit in this branch

    # which of the two competing events fired, weighted by relative rate
    is_silencing <- stats::runif(1) < (silencing_rate / event_rate)
    if (is_silencing) {
      state[] <- missing_state
      return(list(state = state, silenced = TRUE))
    }

    # edit event: fill the leftmost still-unedited site
    # (index labels via sample.int: sample() on a single number would draw from 1:n)
    pos <- min(which(state == "0"))
    state[pos] <- labels[sample.int(length(labels), size = 1, prob = edit_probs)]
  }

  list(state = state, silenced = FALSE)
}
