#' Simulator of Sequentially-Edited (TypeWriter-style) Barcodes
#'
#' Simulates barcodes under a sequential/left-to-right editing process
#' (unedited sites are filled in position order at rate lambda, as a
#' Poisson process, absorbing once all k sites are edited), with an
#' additional whole-barcode silencing state competing against editing.

#'
#' @param tree a treedata object (as returned by sim_adb_origin_samp);
#'   node identity is taken from tree@data$node
#' @param k number of target sites per barcode
#' @param lambda_vec array of barcode-specific editing rates (length m)
#' @param silencing_rate scalar rate at which an entire barcode is
#'   silenced, racing against the editing process as an alternative
#'   "character" draw; 0 = no silencing
#' @param dropout_p probability an individual tip site is masked to the
#'   silenced state after simulation; 0 = no dropout
#' @param m number of tapes per cell
#' @param chars array of unique characters/state labels to be inserted
#'
#' @return a data frame with one row per node in the tree and one column
#'   per barcode copy: \code{node}, \code{barcode_1}, ..., \code{barcode_m},
#'   each entry a length-k string with positions separated by "_"
#'   (e.g. "0_A_B_0_0"), where "0" = unedited and "-" = silenced/dropout
#'   (the silenced state; every site reads "-" once a barcode is silenced)
#' @export
typewriter_barcodes <- function(tree, k, lambda_vec, silencing_rate = 0,
                                dropout_p = 0, m = length(lambda_vec), chars) {
  
  stopifnot(methods::is(tree, "treedata"))
  stopifnot(length(lambda_vec) == m)
  
  silenced_state <- "-"
  # silencing is drawn from the same event race as ordinary edits: each
  # "character" competes at its own rate. chars/sample_p describe the
  # edit-outcome states only; silenced_state is a further, separate
  # competing outcome with its own rate (silencing_rate), not folded into
  # sample_p, since it applies to the WHOLE barcode rather than one site.
  sample_p <- rep(1 / length(chars), length(chars))
  
  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)
  
  # heights: time-since-present (0 = present, root = age of tree)
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  heights <- max(depth_from_root) - depth_from_root
  
  order_df <- tree_df[order(-heights[tree_df$node]), ]
  
  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]
  
  tip_nodes <- tree_df$node[tree_df$status == 1]
  
  bc_cols <- vector("list", m)
  for (bc in seq_len(m)) {
    lambda <- lambda_vec[bc]
    
    state_at <- list()
    silenced_at <- list()
    state_at[[as.character(root)]] <- rep("0", k)
    silenced_at[[as.character(root)]] <- FALSE
    
    if (!is.null(root_edge) && root_edge > 0) {
      res <- .evolve_typewriter_branch(
        state_at[[as.character(root)]], silenced_at[[as.character(root)]],
        origin_height - heights[root], lambda, silencing_rate, chars, sample_p, silenced_state
      )
      state_at[[as.character(root)]] <- res$state
      silenced_at[[as.character(root)]] <- res$silenced
    }
    
    for (i in seq_len(nrow(order_df))) {
      node <- order_df$node[i]
      if (node == root) next
      parent <- order_df$parent[i]
      branch_length <- heights[parent] - heights[node]
      
      res <- .evolve_typewriter_branch(
        state_at[[as.character(parent)]], silenced_at[[as.character(parent)]],
        branch_length, lambda, silencing_rate, chars, sample_p, silenced_state
      )
      state_at[[as.character(node)]] <- res$state
      silenced_at[[as.character(node)]] <- res$silenced
    }
    
    if (dropout_p > 0) {
      for (nd in tip_nodes) {
        key <- as.character(nd)
        if (!silenced_at[[key]]) {
          mask <- stats::runif(k) < dropout_p
          state_at[[key]][mask] <- silenced_state
        }
      }
    }
    
    bc_cols[[bc]] <- vapply(
      tree_df$node,
      function(nd) paste(state_at[[as.character(nd)]], collapse = "_"),
      character(1)
    )
  }
  
  out <- data.frame(node = tree_df$node)
  for (bc in seq_len(m)) out[[paste0("barcode_", bc)]] <- bc_cols[[bc]]
  out
}

# Evolve one barcode copy along one branch. Editing and silencing are two
# outcomes of the same competing-event race: at each event, draw whether it
# is an edit (rate lambda, applies to the leftmost unedited site) or
# silencing (rate silencing_rate, applies to the WHOLE barcode at once,
# overwriting every site -- edited and unedited alike -- to silenced_state
# and ending evolution for this branch and all descendants).
.evolve_typewriter_branch <- function(parent_state, parent_silenced, branch_length,
                                      lambda, silencing_rate, chars, sample_p, silenced_state) {
  state <- parent_state
  k <- length(state)
  
  if (parent_silenced || branch_length <= 0) {
    return(list(state = state, silenced = parent_silenced))
  }
  
  t <- 0
  
  repeat {
    saturated <- all(state != "0")
    edit_rate <- if (saturated) 0 else lambda
    event_rate <- edit_rate + silencing_rate
    
    if (event_rate <= 0) break  # no editing left to do and no silencing possible
    
    t <- t + stats::rexp(1, event_rate)
    if (t > branch_length) break
    
    is_silencing <- stats::runif(1) < (silencing_rate / event_rate)
    if (is_silencing) {
      state[] <- silenced_state
      return(list(state = state, silenced = TRUE))
    }
    
    pos <- min(which(state == "0"))
    state[pos] <- sample(chars, size = 1, prob = sample_p)
  }
  
  list(state = state, silenced = FALSE)
}