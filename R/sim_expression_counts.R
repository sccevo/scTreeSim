#' Negative binomial counts from underlying (log) expression values
#'
#' Simulates single-cell counts from the underlying expression values, for
#' example the values at the tips of a tree simulated with
#' \code{\link{sim_expression_ou}}. The count of each cell is drawn from a
#' negative binomial distribution with mean \eqn{s_i \exp(x_i)} and
#' dispersion \code{phi}, i.e. variance \eqn{\mu + \phi \mu^2}, where
#' \eqn{s_i} is the size factor of the cell. With \code{phi = 0} the model
#' reduces to a Poisson distribution.
#'
#' @param x named numeric vector of underlying expression values on the log
#'   scale, one per cell (e.g. the tip nodes of the tree). The names identify
#'   the cells and must be unique.
#' @param phi negative binomial dispersion (non-negative), length 1 or the
#'   same length as \code{x}
#' @param size_factor optional named numeric vector of positive size factors
#'   (library size scaling), one per cell, with the same names as \code{x}
#'   (in any order). If \code{NULL} (default), all cells have size factor 1.
#'
#' @return A named numeric vector of counts (whole numbers), one per cell, in
#'   the order and with the names of \code{x}.
#'
#' @export
sim_expression_counts <- function(x, phi, size_factor = NULL) {

  if (!is.numeric(x) || length(x) == 0 || anyNA(x)) {
    stop("`x` must be a non-empty numeric vector without missing values.", call. = FALSE)
  }
  if (is.null(names(x)) || anyNA(names(x)) || anyDuplicated(names(x)) > 0) {
    stop("`x` must be named with unique cell names.", call. = FALSE)
  }
  if (!is.numeric(phi) || anyNA(phi) || any(phi < 0) ||
      !length(phi) %in% c(1, length(x))) {
    stop("`phi` must be a non-negative numeric vector of length 1 or length(x).",
         call. = FALSE)
  }

  # without size factors, every cell is scaled equally
  if (is.null(size_factor)) {
    size_factor <- stats::setNames(rep(1, length(x)), names(x))
  }
  if (is.null(names(size_factor))) {
    stop("`size_factor` must be named.", call. = FALSE)
  }
  if (anyDuplicated(names(size_factor)) > 0 ||
      !setequal(names(x), names(size_factor))) {
    stop("The names of `size_factor` must match the names of `x`.", call. = FALSE)
  }
  if (!is.numeric(size_factor) || anyNA(size_factor) || any(size_factor <= 0)) {
    stop("`size_factor` must be a numeric vector of positive values.", call. = FALSE)
  }

  # size_factor is matched to x by name, so the order of its entries is free
  counts <- stats::rnbinom(length(x),
                           size = 1 / phi,
                           mu = size_factor[names(x)] * exp(x))

  return(stats::setNames(counts, names(x)))
}


#' Negative binomial counts for multiple genes at the tips of a tree
#'
#' Wrapper around \code{\link{sim_expression_counts}} that simulates the counts
#' of every gene in the output of \code{\link{sim_expression_matrix}}. Only the
#' sampled tips (status 1 in \code{tree@data}) are cells, so the internal and
#' dead nodes are dropped. Each gene is simulated independently, with its own
#' dispersion, and all genes share the size factor of each cell.
#'
#' @param expr_matrix data frame with columns node, type and one numeric
#'   column per gene with the underlying expression values on the log scale,
#'   as returned by \code{\link{sim_expression_matrix}}
#' @param tree the treedata object that \code{expr_matrix} was simulated on,
#'   used to identify the sampled tips
#' @param phi negative binomial dispersion (non-negative) of each gene: a
#'   numeric vector of length 1 (the same value for every gene) or with one
#'   value per gene
#' @param size_factor optional named numeric vector of positive size factors,
#'   one per cell and shared by all genes. The names are the node ids of the
#'   tips (in any order). If \code{NULL} (default), all cells have size factor 1.
#'
#' @return A data frame with one row per sampled tip and columns node, type,
#'   followed by the gene columns of \code{expr_matrix} (same names), with the
#'   simulated counts (whole numbers) of that gene in that cell.
#'
#' @export
sim_expression_counts_matrix <- function(expr_matrix, tree, phi, size_factor = NULL) {

  stopifnot(methods::is(tree, "treedata"))
  if (!is.data.frame(expr_matrix) || !all(c("node", "type") %in% colnames(expr_matrix))) {
    stop("`expr_matrix` must be a data frame with columns 'node' and 'type' ",
         "and one column per gene, as returned by sim_expression_matrix().",
         call. = FALSE)
  }
  gene_cols <- setdiff(colnames(expr_matrix), c("node", "type"))
  if (length(gene_cols) == 0) {
    stop("`expr_matrix` has no gene columns.", call. = FALSE)
  }
  if (!all(vapply(expr_matrix[gene_cols], is.numeric, logical(1)))) {
    stop("All gene columns of `expr_matrix` must be numeric.", call. = FALSE)
  }

  # the cells are the sampled tips; they need an expression value
  tip_nodes <- tree@data$node[tree@data$status == 1]
  if (!all(tip_nodes %in% expr_matrix$node)) {
    stop("`expr_matrix` has no row for some tips of `tree`; was it simulated ",
         "on this tree?", call. = FALSE)
  }
  tips <- expr_matrix[expr_matrix$node %in% tip_nodes, ]

  if (!is.numeric(phi) || anyNA(phi) || any(phi < 0) ||
      !length(phi) %in% c(1, length(gene_cols))) {
    stop("`phi` must be a non-negative numeric vector of length 1 or one value ",
         "per gene (", length(gene_cols), ").", call. = FALSE)
  }
  phi <- rep_len(phi, length(gene_cols))

  # size factors are cell-level, so a single vector serves all genes
  if (!is.null(size_factor) &&
      (is.null(names(size_factor)) ||
       !setequal(names(size_factor), as.character(tips$node)))) {
    stop("`size_factor` must be named by the node ids of the tips of `tree`.",
         call. = FALSE)
  }

  result <- tips[, c("node", "type")]
  for (g in seq_along(gene_cols)) {
    x <- stats::setNames(tips[[gene_cols[g]]], as.character(tips$node))
    result[[gene_cols[g]]] <- unname(sim_expression_counts(x,
                                                           phi = phi[g],
                                                           size_factor = size_factor))
  }
  rownames(result) <- NULL

  return(result)
}
