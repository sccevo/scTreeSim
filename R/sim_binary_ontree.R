#' Binary substitution model along a tree
#'
#' Simulates one or more binary character sequences along a fixed tree under
#' a general binary substitution model. Each site evolves independently with
#' transition matrix \eqn{P(t) = \exp(Qt)}, where \eqn{Q} is parameterized by
#' \code{lambda} and normalized to its equilibrium frequencies.
#'
#' @param tree a treedata object (as returned by e.g. \code{sim_adb_origin_samp});
#'   node identity is taken from \code{tree@data$node}
#' @param l sequence length (number of binary sites)
#' @param lambda rate parameter for the Q matrix; equilibrium frequency of state
#'   \code{"0"} is \code{lambda / (lambda + 1)}
#' @param m number of independent binary sequences per node
#' @param lambda_vec optional vector of length \code{m} giving a separate
#'   \code{lambda} for each sequence; if \code{NULL}, all copies use \code{lambda}
#' @param root_seq optional integer vector of length \code{l} with values
#'   \code{1 = "0"}, \code{2 = "1"} at the root; if \code{NULL}, sites are
#'   drawn from the equilibrium distribution
#'
#' @return A data frame with one row per node and one column per binary copy:
#'   \code{node}, \code{binary_1}, ..., \code{binary_m}. Each entry is a
#'   length-\code{l} string with sites separated by \code{"_"} (e.g.
#'   \code{"0_1_0_1"}).
#' @export
sim_binary_ontree <- function(tree, l, lambda = 1, m = 1,
                              lambda_vec = NULL, root_seq = NULL) {

  stopifnot(methods::is(tree, "treedata"))
  stopifnot(l >= 1L, lambda > 0)

  if (is.null(lambda_vec)) {
    lambda_vec <- rep(lambda, m)
  }
  stopifnot(length(lambda_vec) == m, all(lambda_vec > 0))

  if (!is.null(root_seq)) {
    stopifnot(length(root_seq) == l, all(root_seq %in% c(1L, 2L)))
  }

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  heights <- max(depth_from_root) - depth_from_root

  order_df <- tree_df[order(-heights[tree_df$node]), ]

  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]

  binary_cols <- vector("list", m)
  for (bc in seq_len(m)) {
    lam <- lambda_vec[bc]
    pi <- .binary_pi(lam)

    state_at <- list()
    state_at[[as.character(root)]] <- if (is.null(root_seq)) {
      .generate_binary_seq(l, pi)
    } else {
      as.integer(root_seq)
    }

    if (!is.null(root_edge) && root_edge > 0) {
      P <- .binary_transition_probs(lam, origin_height - heights[root])
      state_at[[as.character(root)]] <- .mutate_binary_seq(
        state_at[[as.character(root)]], P
      )
    }

    for (i in seq_len(nrow(order_df))) {
      node <- order_df$node[i]
      if (node == root) next
      parent <- order_df$parent[i]
      branch_length <- heights[parent] - heights[node]

      P <- .binary_transition_probs(lam, branch_length)
      state_at[[as.character(node)]] <- .mutate_binary_seq(
        state_at[[as.character(parent)]], P
      )
    }

    binary_cols[[bc]] <- vapply(
      tree_df$node,
      function(nd) paste(.translate_binary_seq(state_at[[as.character(nd)]]), collapse = "_"),
      character(1)
    )
  }

  out <- data.frame(node = tree_df$node)
  for (bc in seq_len(m)) out[[paste0("binary_", bc)]] <- binary_cols[[bc]]
  out
}

.binary_states <- c("0", "1")

.binary_pi <- function(lambda) {
  c(lambda / (lambda + 1), 1 / (lambda + 1))
}

# Closed-form P(t) = exp(Qt) for the normalized 2-state CTMC; equivalent to
# expm(get_q(lambda) * t) in the reference implementation.
.binary_transition_probs <- function(lambda, t) {
  if (t <= 0) {
    return(diag(2))
  }

  pi <- .binary_pi(lambda)
  beta <- 1 / (pi[1] + pi[2] * lambda)
  rate <- beta * (1 + lambda)
  decay <- exp(-rate * t)

  P <- matrix(0, 2, 2)
  for (i in seq_len(2)) {
    for (j in seq_len(2)) {
      P[i, j] <- pi[j] + decay * ((i == j) - pi[j])
    }
  }
  P
}

.generate_binary_seq <- function(l, pi) {
  sample.int(2, l, replace = TRUE, prob = pi)
}

.mutate_binary_seq <- function(seq, P) {
  vapply(seq, function(s) sample.int(2, 1, prob = P[s, ]), integer(1))
}

.translate_binary_seq <- function(seq) {
  .binary_states[seq]
}
