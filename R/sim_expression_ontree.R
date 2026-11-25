#' Ornstein-Uhlenbeck Model for Gene Expression along a tree
#'
#' Code originally written by Dr. Laura Tomas Lopez, Computational Evolution Group.
#' 
#' Simulates a single gene expression trait along a fixed tree under an
#' Ornstein-Uhlenbeck (OU) process, with a type-dependent optimum (theta).
#' Setting \code{alpha = 0} reduces the model to Brownian motion.
#'
#' @param tree a treedata object
#' @param alpha strength of selection (mean-reversion rate)
#' @param sigmasq variance of the stochastic driving noise
#' @param root_value expression value at the root
#' @param theta data frame with columns type and theta, giving
#'   the OU optimum for each type in the tree.
#'
#' @return A data frame with columns node, type, expr, with
#'   one row per node in the tree (tips and internal nodes), giving the
#'   simulated expression value at that node.
#'
#' @export
sim_expression_ou <- function(tree,
                                       alpha = 0.5,
                                       sigmasq = 0.2,
                                       root_value = 0,
                                       theta = NULL # data frame with columns 'type' and 'theta'
) {

  if (is.null(theta) || !all(c("type", "theta") %in% colnames(theta))) {
    message("Error: 'theta' must be a data frame with columns 'type' and 'theta'.")
    return(NULL)
  }

  # every type present in the tree needs a matching optimum, or downstream
  # rows for that type would silently evaluate to NA
  tree_types <- unique(tree@data$type)
  if (!all(tree_types %in% theta$type)) {
    missing_types <- setdiff(tree_types, theta$type)
    message(
      "Error: theta is missing an entry for type(s): ",
      paste(missing_types, collapse = ", "),
      ". Provide theta as a data frame with columns 'type' and 'theta' ",
      "covering every type present in tree@data$type."
    )
    return(NULL)
  }

  sigma <- sqrt(sigmasq)

  evolved_expression <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- evolved_expression %>% dplyr::filter(parent == node) %>% .$node

  # expr is filled in progressively as the tree is walked; only the root
  # is known at the start, everything else begins as NA
  evolved_expression$expr <- NA
  evolved_expression$expr[which(evolved_expression$node == root)] <- root_value

  # draw each node's own stochastic noise term up front, once, since it
  # only depends on that node's branch length and not on the parent's
  # value (alpha = 0 reduces the OU noise term to ordinary Brownian motion)
  evolved_expression_noise <- evolved_expression %>%
    dplyr::filter(node != root) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(noise = dplyr::case_when(
      alpha == 0 ~ stats::rnorm(
        1,
        mean = 0,
        sd = sigma * sqrt(branch.length)
      ),
      TRUE ~ stats::rnorm(
        1,
        mean = 0,
        sd = sigma * sqrt((1 - exp(-2 * alpha * branch.length)) / (2 * alpha))
      )
    )) %>%
    dplyr::ungroup() %>%
    dplyr::select(node, noise)

  evolved_expression <- evolved_expression %>%
    dplyr::left_join(evolved_expression_noise, by = "node") %>%
    dplyr::left_join(theta, by = "type")

  # the root has no incoming branch, so it carries no noise of its own
  evolved_expression$noise[which(evolved_expression$node == root)] <- 0

  # propagate expr from parent to child one tree level at a time: each
  # pass joins in every node's parent value and fills in expr wherever the
  # parent is already known, repeating until every node has been reached
  while (anyNA(evolved_expression$expr)) {

    evolved_expression <- evolved_expression %>%
      dplyr::select(-dplyr::any_of("parent_expr")) %>% # to allow overwriting
      dplyr::left_join(
        evolved_expression %>% dplyr::select(parent = node, parent_expr = expr),
        by = "parent"
      ) %>%
      dplyr::mutate(expr = dplyr::case_when(
        !is.na(expr) ~ expr,
        is.na(parent_expr) ~ expr,
        TRUE ~ theta + (parent_expr - theta) * exp(-alpha * branch.length) + noise
      ))
  }

  evolved_expression <- evolved_expression %>% dplyr::select(node, type, expr)

  return(evolved_expression)
}

#' Ornstein-Uhlenbeck Model for the expression of multiple genes along a tree
#'
#' Simulates several genes independently on the same tree, each with its own
#' OU parameters. It is the model of \code{\link{sim_expression_ou}}, simulated
#' for all the genes at once: the tree is traversed a single time and the genes
#' are handled as columns of a matrix, which is much faster than calling
#' \code{\link{sim_expression_ou}} once per gene. The random numbers are drawn
#' in a different order, so for a given seed the values are not the same as
#' those of \code{\link{sim_expression_ou}}, but they follow the same
#' distribution.
#'
#' The OU parameters \code{alpha}, \code{sigmasq}, \code{root_value} and
#' \code{theta}, and the optional \code{gene_names}, give one value per gene.
#' Each must have length 1 (the same value is used for every gene) or length
#' \code{N}, where \code{N} is the number of genes. For \code{theta}, the
#' length is the number of gene columns (all columns except \code{type}).
#' \code{N} is the length of the parameters that have length greater than 1.
#' If all of them have length 1, \code{N} is 1, unless \code{n_genes} is
#' given, in which case it sets \code{N} and the length-1 parameters are
#' shared by all the genes. If some parameter has length greater than 1,
#' \code{n_genes} must equal that length.
#'
#' @inheritParams sim_expression_ou
#' @param alpha numeric vector (length 1 or \code{N}): strength of selection
#'   for each gene
#' @param sigmasq numeric vector (length 1 or \code{N}): variance of the
#'   stochastic driving noise for each gene
#' @param root_value numeric vector (length 1 or \code{N}): expression value
#'   at the root for each gene
#' @param theta data frame with a column \code{type} and one numeric column
#'   per gene (\code{N} columns) giving the OU optimum of that gene for each
#'   type in the tree. With a single gene column (any name), the same optima are 
#'   used for every gene.
#' @param n_genes number of genes \code{N}. If \code{NULL} (default), it is
#'   inferred from the length of the parameters. It can only differ from the
#'   length of the parameters when all of them have length 1.
#' @param gene_names optional character vector of unique gene names used as
#'   column names (length 1 or \code{N}; a single name is only valid for
#'   \code{N = 1}). Like the other parameters, if it has length \code{N} it
#'   sets the number of genes. If \code{NULL} (default), genes are named by
#'   their id (\code{1}, ..., \code{N}).
#'
#' @return A data frame with one row per node in the tree (tips and internal
#'   nodes) and columns node, type, followed by one column per gene containing
#'   the simulated expression of that gene at that node. The gene columns are
#'   named by \code{gene_names}, or by gene id (\code{1}, ..., \code{N}) by
#'   default.
#'
#' @export
sim_expression_matrix <- function(tree,
                                  alpha = 0.5,
                                  sigmasq = 0.2,
                                  root_value = 0,
                                  theta = NULL,
                                  n_genes = NULL,
                                  gene_names = NULL) {

  if (!is.data.frame(theta) || !"type" %in% colnames(theta) || ncol(theta) < 2) {
    stop("`theta` must be a data frame with a column 'type' and one column ",
         "per gene (or a single gene column shared by all genes).", call. = FALSE)
  }
  gene_cols <- setdiff(colnames(theta), "type")
  if (!all(vapply(theta[gene_cols], is.numeric, logical(1)))) {
    stop("All gene columns of `theta` must be numeric.", call. = FALSE)
  }

  for (nm in c("alpha", "sigmasq", "root_value")) {
    x <- get(nm)
    if (!is.numeric(x) || length(x) == 0 || anyNA(x)) {
      stop("`", nm, "` must be a numeric vector without missing values.", call. = FALSE)
    }
  }

  if (!is.null(n_genes) &&
      (length(n_genes) != 1 || !is.numeric(n_genes) || is.na(n_genes) ||
       n_genes < 1 || n_genes != round(n_genes))) {
    stop("`n_genes` must be NULL or a single positive integer.", call. = FALSE)
  }

  # every parameter must have length 1 (shared by all genes) or length N
  par_lengths <- c(alpha = length(alpha), sigmasq = length(sigmasq),
                   root_value = length(root_value), theta = length(gene_cols))
  if (!is.null(gene_names)) {
    par_lengths <- c(par_lengths, gene_names = length(gene_names))
  }

  # N is given by the parameters of length > 1; n_genes can only set it when
  # all of them have length 1
  n_params <- max(par_lengths)
  bad <- par_lengths != 1 & par_lengths != n_params
  if (any(bad)) {
    stop("All parameters must have length 1 or the number of genes (",
         n_params, "), but ",
         paste0("`", names(par_lengths)[bad], "` has length ", par_lengths[bad],
                collapse = ", "),
         ".", call. = FALSE)
  }
  if (is.null(n_genes)) {
    n_genes <- n_params
  } else {
    n_genes <- as.integer(n_genes)
    if (n_params > 1 && n_genes != n_params) {
      stop("`n_genes` (", n_genes, ") must equal the length of the parameters (",
           n_params, "); it can only differ when all parameters have length 1.",
           call. = FALSE)
    }
  }

  if (!is.null(gene_names)) {
    gene_names <- as.character(gene_names)
    if (anyNA(gene_names) || any(gene_names == "")) {
      stop("`gene_names` must be non-empty and without missing values.",
           call. = FALSE)
    }
    # a length-1 gene_names is not recycled: names must be unique
    if (length(gene_names) != n_genes) {
      stop("`gene_names` must be unique; a single name cannot be shared by ",
           n_genes, " genes.", call. = FALSE)
    }
    if (anyDuplicated(gene_names) > 0) {
      stop("`gene_names` must be unique.", call. = FALSE)
    }
    if (any(gene_names %in% c("node", "type"))) {
      stop("`gene_names` cannot contain 'node' or 'type', which are reserved ",
           "column names.", call. = FALSE)
    }
  } else {
    gene_names <- as.character(seq_len(n_genes))
  }

  # recycle length-1 parameters to one value per gene
  alpha <- rep_len(alpha, n_genes)
  sigmasq <- rep_len(sigmasq, n_genes)
  root_value <- rep_len(root_value, n_genes)
  gene_cols <- rep_len(gene_cols, n_genes)

  # tree table (one row per node), in the same order as sim_expression_ou
  nodes <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- nodes$node[nodes$parent == nodes$node]

  # every type present in the tree needs a matching optimum (as in sim_expression_ou)
  tree_types <- unique(tree@data$type)
  if (!all(tree_types %in% theta$type)) {
    missing_types <- setdiff(tree_types, theta$type)
    message(
      "Error: theta is missing an entry for type(s): ",
      paste(missing_types, collapse = ", "),
      ". Provide theta as a data frame with a column 'type' covering every ",
      "type present in tree@data$type."
    )
    return(NULL)
  }

  n_nodes <- nrow(nodes)
  parent_idx <- match(nodes$parent, nodes$node)
  is_root <- nodes$node == root
  theta_mat <- as.matrix(theta[gene_cols])[match(nodes$type, theta$type), , drop = FALSE]

  # level of each node (number of branches from the root), so that all nodes of
  # a level can be updated at once, for every gene
  level <- rep(NA_integer_, n_nodes)
  level[is_root] <- 0L
  while (anyNA(level)) {
    level <- ifelse(is.na(level) & !is.na(level[parent_idx]), level[parent_idx] + 1L, level)
  }

  # noise of the incoming branch of every node and gene: N(0, sd^2) with
  # sd^2 = sigmasq * (1 - exp(-2 alpha t)) / (2 alpha), or sigmasq * t for
  # alpha = 0 (Brownian motion). alpha = 0 is replaced by 1 in the OU formula
  # to avoid 0/0; those entries are overwritten by the Brownian one
  bl <- nodes$branch.length
  bl[is_root] <- 0
  alpha_mat <- matrix(alpha, n_nodes, n_genes, byrow = TRUE)
  sigmasq_mat <- matrix(sigmasq, n_nodes, n_genes, byrow = TRUE)
  alpha_safe <- ifelse(alpha_mat == 0, 1, alpha_mat)
  var_ou <- sigmasq_mat * (1 - exp(-2 * alpha_safe * bl)) / (2 * alpha_safe)
  noise_var <- ifelse(alpha_mat == 0, sigmasq_mat * bl, var_ou)
  noise <- matrix(stats::rnorm(n_nodes * n_genes), n_nodes, n_genes) * sqrt(noise_var)

  # propagate expression from the root, level by level: the value of a node is
  # the OU transition from its parent's value over the branch, plus its noise
  expr <- matrix(NA_real_, n_nodes, n_genes)
  expr[is_root, ] <- matrix(root_value, 1, n_genes)
  for (lv in seq_len(max(level))) {
    idx <- which(level == lv)
    decay <- exp(-outer(bl[idx], alpha))
    expr[idx, ] <- theta_mat[idx, , drop = FALSE] +
      (expr[parent_idx[idx], , drop = FALSE] - theta_mat[idx, , drop = FALSE]) * decay +
      noise[idx, , drop = FALSE]
  }
  colnames(expr) <- gene_names

  result <- data.frame(node = nodes$node, type = nodes$type, expr, check.names = FALSE)

  return(result)
}
