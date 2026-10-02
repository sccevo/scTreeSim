#' GT16 diploid nucleotide SNV simulator
#'
#' The GT16 substitution model and its error model, and the respective code
#' logic, have been adapted from the BEAST phylonco paper and software as
#' cited below.
#' Paper: Chen, K., Moravec, J. C., Gavryushkin, A., Welch, D., & Drummond, A. J. (2022).
#' Accounting for errors in data improves divergence time estimates in single-cell cancer evolution.
#' Molecular biology and evolution, 39(8), msac143.
#' Beast phylonco software: A BEAST2 package for single-cell phylogenetic analysis of cancer evolution.
#' source: https://github.com/bioDS/beast-phylonco-paper/tree/main
#' The GT16 model itself is due to Kozlov et al. (2022), CellPhy, Genome Biology 23, 37.
#'

#'
#' @param tree a treedata object
#' @param sequences the root genotype, given as the two aligned allele
#'   sequences of length `l`, e.g. `c("CGTAAC", "CGTATC")`; the two strings
#'   are paired site by site into the GT16 states (here CC, GG, TT, AA, AT,
#'   CC). Supply this or `l`, not both
#' @param rates the six nucleotide exchangeabilities, in the order
#'   \eqn{r_{AC}, r_{AG}, r_{AT}, r_{CG}, r_{CT}, r_{GT}}; a common scaling
#'   of all six is removed by the rate normalization
#' @param pi equilibrium frequencies: either 16 genotype frequencies in the
#'   state order above, or 4 nucleotide frequencies (in the order A, C, G, T)
#'   from which the genotype frequencies are formed as
#'   \eqn{\pi_{ab} = \pi_a \pi_b}; `NULL` (default) for a uniform 1/16
#' @param l number of SNV sites, used in place of `sequences` to draw the
#'   root genotypes from `pi` instead of fixing them
#' @param clock_rate substitution rate applied per unit branch length
#' @param epsilon combined amplification and sequencing error probability
#' @param delta allelic dropout probability
#'
#' @return a data frame with one row per node in the tree and columns
#'   \code{node}, \code{genotype} (the true simulated sequence) and
#'   \code{observed} (the sequence after the error model), each entry a
#'   string of `l` genotypes separated by "_" (e.g. "AA_CG_TT"). The error
#'   model applies to the sampled tips (\code{status == 1}), so
#'   \code{observed} equals \code{genotype} at every other node.
#' @export
sim_gt16_nt_diploid <- function(tree, sequences = NULL, rates = rep(1, 6), pi = NULL,
                                l = NULL, clock_rate = 1,
                                epsilon = 0, delta = 0) {

  stopifnot(methods::is(tree, "treedata"))
  stopifnot(clock_rate > 0)

  # the root genotypes are either read off the two aligned sequences or,
  # failing that, drawn from the equilibrium distribution over `l` sites
  if (is.null(sequences) == is.null(l)) {
    stop("Supply either `sequences` (the two aligned root sequences) or `l` (the number of sites), not both.",
         call. = FALSE)
  }
  root_state <- if (is.null(sequences)) NULL else .gt16_pair_sequences(sequences)
  l <- if (is.null(root_state)) as.integer(l) else length(root_state)
  stopifnot(l >= 1L)

  if (length(rates) != 6 || anyNA(rates) || any(rates < 0) || sum(rates) <= 0) {
    stop("`rates` must be six non-negative exchangeabilities (AC, AG, AT, CG, CT, GT), not all zero.",
         call. = FALSE)
  }
  pi <- .gt16_pi(pi)
  if (epsilon < 0 || epsilon > 1 || delta < 0 || delta > 1) {
    stop("`epsilon` and `delta` must be probabilities in [0, 1].", call. = FALSE)
  }

  events <- .gt16_events(rates, pi)

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)

  # heights are measured as time since the present
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  heights <- max(depth_from_root) - depth_from_root
  order_df <- tree_df[order(-heights[tree_df$node]), ]

  root_edge <- tree@phylo$root.edge

  tip_nodes <- tree_df$node[tree_df$status == 1]

  state_at <- list()
  state_at[[as.character(root)]] <- if (is.null(root_state)) {
    sample.int(length(.gt16_states), l, replace = TRUE, prob = pi)
  } else {
    root_state
  }

  # evolve along the origin-to-root edge first, if the tree has one
  if (!is.null(root_edge) && root_edge > 0) {
    state_at[[as.character(root)]] <- .evolve_branch_gt16(
      state_at[[as.character(root)]], clock_rate * root_edge, events
    )
  }

  # then walk every remaining branch, parent genotype -> child genotype, so
  # that a node's genotypes are inherited before any further substitution
  for (i in seq_len(nrow(order_df))) {
    node <- order_df$node[i]
    if (node == root) next
    parent <- order_df$parent[i]
    branch_length <- heights[parent] - heights[node]

    state_at[[as.character(node)]] <- .evolve_branch_gt16(
      state_at[[as.character(parent)]], clock_rate * branch_length, events
    )
  }

  # sequencing is a measurement effect, so the observed genotypes are a
  # separate copy, derived from the finished true genotypes at the tips
  observed_at <- state_at
  if (epsilon > 0 || delta > 0) {
    error_matrix <- .gt16_error_matrix(epsilon, delta)
    for (nd in tip_nodes) {
      key <- as.character(nd)
      observed_at[[key]] <- .add_gt16_error(state_at[[key]], error_matrix)
    }
  }

  collapse_seq <- function(states) paste(.gt16_states[states], collapse = "_")
  data.frame(
    node = tree_df$node,
    genotype = vapply(tree_df$node, function(nd) collapse_seq(state_at[[as.character(nd)]]), character(1)),
    observed = vapply(tree_df$node, function(nd) collapse_seq(observed_at[[as.character(nd)]]), character(1)),
    stringsAsFactors = FALSE
  )
}

.gt16_nucleotides <- c("A", "C", "G", "T")

# the 16 ordered diploid states: AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
.gt16_states <- paste0(rep(.gt16_nucleotides, each = 4), rep(.gt16_nucleotides, times = 4))

# allele indices of each of the 16 states, used throughout to decide which
# copy a candidate change would affect
.gt16_allele1 <- rep(seq_along(.gt16_nucleotides), each = 4)
.gt16_allele2 <- rep(seq_along(.gt16_nucleotides), times = 4)

# accept either 16 genotype frequencies or 4 nucleotide frequencies, in
# which case the genotype frequencies are the products pi_a * pi_b
.gt16_pi <- function(pi) {
  if (is.null(pi)) return(rep(1 / 16, 16))
  if (anyNA(pi) || any(pi < 0) || !isTRUE(all.equal(sum(pi), 1))) {
    stop("`pi` must be non-negative frequencies summing to 1.", call. = FALSE)
  }
  # nucleotide frequencies: pi_ab = pi_a * pi_b, in the GT16 state order
  if (length(pi) == 4) return(rep(pi, each = 4) * rep(pi, times = 4))
  if (length(pi) != 16) {
    stop("`pi` must have length 16 (genotype frequencies) or 4 (nucleotide frequencies).", call. = FALSE)
  }
  as.vector(pi)
}

# Pair the two aligned allele sequences site by site into GT16 states, so
# that c("CGTAAC", "CGTATC") becomes CC, GG, TT, AA, AT, CC. The first
# sequence supplies the first allele of every state and the second the
# second, which is what makes the states ordered (phased).
.gt16_pair_sequences <- function(sequences) {
  if (!is.character(sequences) || length(sequences) != 2 || anyNA(sequences)) {
    stop("`sequences` must be a character vector of the two aligned allele sequences, e.g. c(\"CGTAAC\", \"CGTATC\").",
         call. = FALSE)
  }

  alleles <- lapply(strsplit(toupper(sequences), "", fixed = TRUE), match, table = .gt16_nucleotides)
  if (length(alleles[[1]]) != length(alleles[[2]])) {
    stop("The two sequences in `sequences` must be aligned, i.e. of equal length.", call. = FALSE)
  }
  if (length(alleles[[1]]) == 0 || anyNA(unlist(alleles))) {
    stop("`sequences` must be non-empty and contain only the nucleotides A, C, G and T.", call. = FALSE)
  }

  (alleles[[1]] - 1L) * 4L + alleles[[2]]
}

# the 4x4 symmetric nucleotide exchangeabilities, from the six GTR-like
# parameters in the order AC, AG, AT, CG, CT, GT
.gt16_nuc_exchangeability <- function(rates) {
  R <- matrix(0, 4, 4, dimnames = list(.gt16_nucleotides, .gt16_nucleotides))
  R[1, 2] <- R[2, 1] <- rates[1]  # A <-> C
  R[1, 3] <- R[3, 1] <- rates[2]  # A <-> G
  R[1, 4] <- R[4, 1] <- rates[3]  # A <-> T
  R[2, 3] <- R[3, 2] <- rates[4]  # C <-> G
  R[2, 4] <- R[4, 2] <- rates[5]  # C <-> T
  R[3, 4] <- R[4, 3] <- rates[6]  # G <-> T
  R
}

# The full normalized 16x16 GT16 rate matrix. Only used to build the event
# rates once per call (and as the reference the tests check against); the
# simulation itself never exponentiates it.
.gt16_rate_matrix <- function(rates, pi) {
  R <- .gt16_nuc_exchangeability(rates)
  n <- length(.gt16_states)

  Q <- matrix(0, n, n, dimnames = list(.gt16_states, .gt16_states))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      # exactly one allele may differ; a change at both copies at once
      # (e.g. AG -> CT) is not an instantaneous event
      exch <- if (.gt16_allele1[i] == .gt16_allele1[j]) {
        R[.gt16_allele2[i], .gt16_allele2[j]]
      } else if (.gt16_allele2[i] == .gt16_allele2[j]) {
        R[.gt16_allele1[i], .gt16_allele1[j]]
      } else {
        0
      }
      Q[i, j] <- exch * pi[j]
    }
  }

  # normalize so that one unit of branch length is one expected substitution
  d <- -1 * rowSums(Q)
  diag(Q) <- d
  total_rate <- as.vector(pi %*% d)
  if (total_rate >= 0) {
    # every state carrying equilibrium mass is a dead end, so there is no
    # timescale to normalize against
    stop("`rates` and `pi` give a degenerate GT16 rate matrix with no substitutions.", call. = FALSE)
  }
  (-1 / total_rate) * Q
}

# The six states reachable from each state by a single-allele substitution:
# three for the first allele copy, three for the second.
.gt16_neighbours <- function() {
  n <- length(.gt16_states)
  dest <- matrix(0L, n, 6L)
  for (i in seq_len(n)) {
    a1 <- .gt16_allele1[i]
    a2 <- .gt16_allele2[i]
    dest[i, ] <- c(
      (setdiff(1:4, a1) - 1L) * 4L + a2,  # first allele substituted
      (a1 - 1L) * 4L + setdiff(1:4, a2)   # second allele substituted
    )
  }
  dest
}

# Event rates taken directly from the off-diagonal entries of the normalized
# GT16 Q, so the mechanistic simulation and the matrix model cannot drift
# apart. `cumulative` is carried so that drawing an event is one comparison
# per site, and its last column doubles as the total outgoing rate (sharing
# the same floating-point sum keeps the draw inside the six events).
.gt16_events <- function(rates, pi) {
  Q <- .gt16_rate_matrix(rates, pi)
  dest <- .gt16_neighbours()
  n <- nrow(dest)
  rate <- matrix(Q[cbind(rep(seq_len(n), ncol(dest)), as.vector(dest))], n, ncol(dest))
  cumulative <- t(apply(rate, 1, cumsum))
  list(Q = Q, dest = dest, rate = rate, cumulative = cumulative,
       total = cumulative[, ncol(cumulative)])
}

# Gillespie simulation of the GT16 chain along one branch, run for all sites
# at once: sites still "in flight" draw their next waiting time together,
# those whose next event falls past the end of the branch drop out, and the
# rest take one single-allele substitution and go round again. Sites that
# never have an event simply inherit the parent genotype.
.evolve_branch_gt16 <- function(state, branch_time, events) {
  if (branch_time <= 0) return(state)

  elapsed <- numeric(length(state))
  active <- seq_along(state)

  while (length(active) > 0) {
    s <- state[active]
    # a state with no outgoing rate gives an infinite waiting time and so
    # simply leaves the loop below
    elapsed[active] <- elapsed[active] + stats::rexp(length(active), events$total[s])

    fits <- elapsed[active] <= branch_time
    active <- active[fits]
    if (length(active) == 0) break

    # which of the six events fired, with probability proportional to its rate
    s <- s[fits]
    u <- stats::runif(length(s)) * events$total[s]
    k <- rowSums(events$cumulative[s, , drop = FALSE] < u) + 1L
    state[active] <- events$dest[cbind(s, k)]
  }

  state
}

# P(observed | true) for the GT16 error model of Chen et al. (2022), eq. (2),
# with combined amplification/sequencing error epsilon and allelic dropout
# delta. Rows are the true genotype, columns the observed one.
.gt16_error_matrix <- function(epsilon, delta) {
  n <- length(.gt16_states)
  M <- matrix(0, n, n, dimnames = list(.gt16_states, .gt16_states))

  for (y in seq_len(n)) {
    a <- .gt16_allele1[y]
    b <- .gt16_allele2[y]
    for (x in seq_len(n)) {
      o1 <- .gt16_allele1[x]
      o2 <- .gt16_allele2[x]

      M[y, x] <- if (a == b) {
        # true genotype is homozygous aa
        if (o1 == a && o2 == a) {
          1 - epsilon + 0.5 * delta * epsilon                     # P(aa|aa)
        } else if (o1 == a || o2 == a) {
          (1 - delta) * epsilon / 6                               # P(ab|aa), P(ba|aa)
        } else if (o1 == o2) {
          delta * epsilon / 6                                     # P(bb|aa)
        } else {
          0
        }
      } else {
        # true genotype is heterozygous ab
        if (o1 == a && o2 == b) {
          (1 - delta) * (1 - epsilon)                             # P(ab|ab)
        } else if (o1 == o2 && (o1 == a || o1 == b)) {
          0.5 * delta + epsilon / 6 - delta * epsilon / 3         # P(aa|ab), P(bb|ab)
        } else if (o1 == o2) {
          delta * epsilon / 6                                     # P(cc|ab)
        } else if (o1 == a || o2 == b) {
          (1 - delta) * epsilon / 6                               # P(ac|ab), P(cb|ab)
        } else {
          0
        }
      }
    }
  }
  M
}

# applies the error model independently to each site of an already-simulated
# sequence, leaving the true sequence it was given untouched
.add_gt16_error <- function(state, error_matrix) {
  vapply(state, function(s) sample.int(nrow(error_matrix), 1, prob = error_matrix[s, ]), integer(1))
}
