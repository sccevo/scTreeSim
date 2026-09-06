#' Binary substitution model along a tree, with error and gamma spikes
#'
#' Simulates one or more binary character sequences along a fixed tree under
#' a general binary substitution model. Each site evolves independently with
#' transition matrix \eqn{P(t) = \exp(Qt)}. Optionally adds a gamma-distributed
#' burst of extra mutation at branches where a division occurred (status == 2
#' in tree@data, including internal nodes later collapsed into a single
#' surviving branch), and applies false positive/negative error at the tips.
#'
#' @param tree a treedata object
#' @param l sequence length (number of binary sites)
#' @param lambda rate parameter for the Q matrix
#' @param m number of independent binary sequences per node
#' @param lambda_vec optional vector of length m, one lambda per copy
#' @param root_seq optional integer vector of length l, 1 = "0", 2 = "1"
#' @param clock_rate substitution rate applied per unit branch length
#' @param gamma_spike if TRUE, add extra mutation at branches where a
#'   division occurred
#' @param spike_mean mean of the gamma-distributed spike size
#' @param spike_shape shape of the gamma-distributed spike size
#' @param alpha false negative rate (P(observed 0 | true 1))
#' @param beta false positive rate (P(observed 1 | true 0))
#'
#' @return A data frame with one row per node and one column per binary copy:
#'   \code{node}, \code{binary_1}, ..., \code{binary_m}. Each entry is a
#'   length-l string with sites separated by "_".
#' @export
sim_binary_ontree <- function(tree, l, lambda = 1, m = 1,
                              lambda_vec = NULL, root_seq = NULL,
                              clock_rate = 1,
                              gamma_spike = FALSE, spike_mean = 0.01, spike_shape = 2,
                              alpha = 0, beta = 0) {
  
  stopifnot(methods::is(tree, "treedata"))
  stopifnot(l >= 1L, lambda > 0, clock_rate > 0)
  stopifnot(spike_mean >= 0, spike_shape > 0)
  stopifnot(alpha >= 0, alpha <= 1, beta >= 0, beta <= 1)
  
  if (is.null(lambda_vec)) lambda_vec <- rep(lambda, m)
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
  
  status_at <- stats::setNames(tree_df$status, tree_df$node)
  tip_nodes <- tree_df$node[tree_df$status == 1]
  
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
      dist <- .evolutionary_distance(
        root_edge, clock_rate, division_occurred = status_at[[as.character(root)]] == 2,
        gamma_spike, spike_mean, spike_shape
      )
      P <- .binary_transition_probs(lam, dist)
      state_at[[as.character(root)]] <- .mutate_binary_seq(state_at[[as.character(root)]], P)
    }
    
    for (i in seq_len(nrow(order_df))) {
      node <- order_df$node[i]
      if (node == root) next
      parent <- order_df$parent[i]
      branch_length <- heights[parent] - heights[node]
      
      dist <- .evolutionary_distance(
        branch_length, clock_rate, division_occurred = status_at[[as.character(node)]] == 2,
        gamma_spike, spike_mean, spike_shape
      )
      P <- .binary_transition_probs(lam, dist)
      state_at[[as.character(node)]] <- .mutate_binary_seq(state_at[[as.character(parent)]], P)
    }
    
    if (alpha > 0 || beta > 0) {
      for (nd in tip_nodes) {
        key <- as.character(nd)
        state_at[[key]] <- .add_error(state_at[[key]], alpha, beta)
      }
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

# P(t) = exp(Qt) for the normalized 2-state CTMC
.binary_transition_probs <- function(lambda, t) {
  if (t <= 0) return(diag(2))
  
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

# clock-driven branch length, plus one gamma-distributed spike if a
# division occurred on this branch
.evolutionary_distance <- function(branch_length, clock_rate, division_occurred,
                                   gamma_spike, spike_mean, spike_shape) {
  gradual <- clock_rate * branch_length
  spike <- if (gamma_spike && division_occurred) {
    spike_mean * stats::rgamma(1, shape = spike_shape, rate = spike_shape)
  } else {
    0
  }
  gradual + spike
}

# false positive/negative error matrix: column = true state, row = observed state
.error_matrix <- function(alpha, beta) {
  matrix(c(1 - beta, alpha, beta, 1 - alpha), nrow = 2, byrow = TRUE)
}

.add_error <- function(seq, alpha, beta) {
  M <- .error_matrix(alpha, beta)
  vapply(seq, function(s) sample.int(2, 1, prob = M[, s]), integer(1))
}