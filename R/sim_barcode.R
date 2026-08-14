#' Simulator of Sequentially-Edited (TypeWriter-style) Barcodes
#'
#' @param tree a phylo object (e.g. tree@phylo from a treedata object)
#' @param k number of target sites per barcode
#' @param lambda_vec array of barcode-specific editing rates (length m)
#' @param m number of barcodes per cell
#' @param chars array of unique characters to be inserted
#'
#' @return a data frame with columns \code{node}, \code{barcode_id}, \code{barcode}
#'   -- one row per node (tips and internal) per barcode copy. \code{node}
#'   matches the same node numbering used by tree@data in your ADB output,
#'   so this joins directly onto the gene-expression output and onto
#'   tree@data by \code{node}.
#' @export
typewriter_barcodes <- function(tree, k, lambda_vec, m, chars) {
  root <- paste0(rep("0", k), collapse = "")
  sample_p <- rep(1 / length(chars), length(chars))
  res <- lapply(seq_along(lambda_vec), function(i) {
    sim_typewriter_barcode(lambda_vec[i], chars, sample_p, tree, root, barcode_id = i)
  })
  dplyr::bind_rows(res)
}

sim_typewriter_barcode <- function(rate, chars, sampling, tree, root, barcode_id) {
  new_root <- colouring(root, tree$root.edge, rate, sampling, chars)
  mut_tree <- ape::rTraitMult(
    tree,
    colouring,
    root.value = new_root,
    ancestor = TRUE,
    lambda = rate,
    sample_p = sampling,
    chars = chars
  )

  node <- as.integer(rownames(mut_tree))
  if (any(is.na(node))) {
    # ape sometimes returns node ids as plain sequential row order instead
    # of row names -- fall back to 1:(Ntip+Nnode) in that case
    ntip <- length(tree$tip.label)
    node <- seq_len(ntip + tree$Nnode)
  }
  data.frame(
    node = node,
    barcode_id = barcode_id,
    barcode = as.character(mut_tree$x1)
  )
}

colouring <- function(x, l, lambda, sample_p, chars) {
  p1 <- stringr::str_split(x, "")[[1]]
  t <- 0
  M <- length(p1)
  while (t <= l) {
    time_to_mut <- stats::rexp(1, lambda)
    t <- t + time_to_mut
    if (p1[M] == "0" && t <= l) {
      pos <- min(which(p1 == "0"))
      new <- sample(chars, size = 1, prob = sample_p)
      p1[pos] <- new
    }
  }
  return(paste(p1, collapse = ""))
}