#' Ornstein-Uhlenbeck Model for Gene Expression along a tree
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