.gt16_nucleotides <- c("A", "C", "G", "T")

# the 16 ordered diploid states: AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
.gt16_states <- paste0(rep(.gt16_nucleotides, each = 4), rep(.gt16_nucleotides, times = 4))

# the nucleotide on each allele copy of each state
.gt16_allele1 <- rep(seq_along(.gt16_nucleotides), each = 4)
.gt16_allele2 <- rep(seq_along(.gt16_nucleotides), times = 4)

#' GT16 diploid nucleotide SNV simulator
#' @details
#' Simulates diploid single-nucleotide variants along a fixed tree under the
#' GT16 model. The 16 states are the ordered pairs of nucleotides
#' \code{AA, AC, AG, AT, CA, CC, ..., TT}
#' Rather than forming \eqn{P(t) = \exp(Qt)} for every
#' branch, the chain is simulated event by event.
#'
#' The GT16 model is due to Kozlov et al. (2022), CellPhy, Genome Biology 23,
#' 37. The error model and the code logic here are adapted from Chen, K.,
#' Moravec, J. C., Gavryushkin, A., Welch, D., & Drummond, A. J. (2022),
#' Accounting for errors in data improves divergence time estimates in
#' single-cell cancer evolution, Molecular Biology and Evolution 39(8),
#' msac143, and the BEAST2 phylonco package
#' (https://github.com/bioDS/beast-phylonco-paper).
#'
#' @param tree a treedata object
#' @param sequences the root genotype, given as the two aligned allele
#'   sequences of length `l`; they are paired site by site, so
#'   `c("CGTAAC", "CGTATC")` is the six-site genotype CC, GG, TT, AA, AT, CC.
#'   Supply this or `l`
#' @param rates the six nucleotide exchangeabilities, in the order
#'   \eqn{r_{AC}, r_{AG}, r_{AT}, r_{CG}, r_{CT}, r_{GT}}; a common scaling of
#'   all six is removed by the rate normalization
#' @param pi equilibrium frequencies: either 16 genotype frequencies in the
#'   state order above, or 4 nucleotide frequencies (A, C, G, T) from which the
#'   genotype frequencies are formed as \eqn{\pi_{ab} = \pi_a \pi_b}; `NULL`
#'   (default) for a uniform 1/16
#' @param l number of SNV sites. Used in place of `sequences` to draw the root
#'   genotypes from `pi` instead of fixing them; alongside `sequences` it is
#'   optional, and must then give their length
#' @param clock_rate substitution rate applied per unit branch length
#' @param epsilon combined amplification and sequencing error probability
#' @param delta allelic dropout probability
#'
#' @return a data frame with one row per node in the tree and columns
#'   \code{node}, \code{genotype} (the true simulated sequence) and
#'   \code{observed} (the sequence after the error model), each entry a string
#'   of `l` genotypes separated by "_" (e.g. "AA_CG_TT"). The error model
#'   applies to the sampled tips (\code{status == 1}), so \code{observed}
#'   equals \code{genotype} at every other node.
#' @export
sim_gt16_nt_diploid <- function(tree, sequences = NULL, rates = rep(1, 6), pi = NULL,
                                l = NULL, clock_rate = 1,
                                epsilon = 0, delta = 0) {

  stopifnot(methods::is(tree, "treedata"))
  stopifnot(clock_rate > 0)

  # the root genotypes are either read off the two aligned sequences or,
  # failing that, drawn from the equilibrium distribution over `l` sites
  if (is.null(sequences) && is.null(l)) {
    stop("Supply `sequences` (the two aligned root sequences) or `l` (the number of sites).", call. = FALSE)
  }
  root_state <- if (is.null(sequences)) NULL else .gt16_pair_sequences(sequences)
  if (!is.null(root_state) && !is.null(l) && as.integer(l) != length(root_state)) {
    # `l` is redundant alongside `sequences`, but if given it has to agree
    stop("`l` (", as.integer(l), ") must match the length of `sequences` (", length(root_state), ").",
         call. = FALSE)
  }
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

  # edges in cladewise order, so a parent is always evolved before its
  # children; sorting nodes on height does not guarantee that, because a
  # zero-length edge leaves parent and child tied
  phylo <- ape::reorder.phylo(tree@phylo, "cladewise")

  state_at <- list()
  state_at[[as.character(root)]] <- if (is.null(root_state)) {
    sample.int(length(.gt16_states), l, replace = TRUE, prob = pi)
  } else {
    root_state
  }

  # evolve along the origin-to-root edge first, if the tree has one
  root_edge <- tree@phylo$root.edge
  if (!is.null(root_edge) && root_edge > 0) {
    state_at[[as.character(root)]] <- .evolve_branch_gt16(
      state_at[[as.character(root)]], clock_rate * root_edge, events
    )
  }

  # then walk every remaining branch, parent genotype -> child genotype
  # (the root never appears as a child, so it is covered by the block above)
  for (i in seq_len(nrow(phylo$edge))) {
    parent <- phylo$edge[i, 1]
    node <- phylo$edge[i, 2]

    state_at[[as.character(node)]] <- .evolve_branch_gt16(
      state_at[[as.character(parent)]], clock_rate * phylo$edge.length[i], events
    )
  }

  # sequencing is a measurement effect, so the observed genotypes are a
  # separate copy, derived site by site from the finished true genotypes
  observed_at <- state_at
  if (epsilon > 0 || delta > 0) {
    error_matrix <- .gt16_error_matrix(epsilon, delta)
    for (node in tree_df$node[tree_df$status == 1]) {
      key <- as.character(node)
      observed_at[[key]] <- vapply(
        state_at[[key]],
        function(s) sample.int(nrow(error_matrix), 1, prob = error_matrix[s, ]),
        integer(1)
      )
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

# accept either 16 genotype frequencies or 4 nucleotide frequencies, in which
# case the genotype frequencies are the products pi_a * pi_b
.gt16_pi <- function(pi) {
  if (is.null(pi)) return(rep(1 / 16, 16))
  if (anyNA(pi) || any(pi < 0) || !isTRUE(all.equal(sum(pi), 1))) {
    stop("`pi` must be non-negative frequencies summing to 1.", call. = FALSE)
  }
  if (length(pi) == 4) return(rep(pi, each = 4) * rep(pi, times = 4))
  if (length(pi) != 16) {
    stop("`pi` must have length 16 (genotype frequencies) or 4 (nucleotide frequencies).", call. = FALSE)
  }
  as.vector(pi)
}

# Pair the two aligned allele sequences site by site into GT16 states, so that
# c("CGTAAC", "CGTATC") becomes CC, GG, TT, AA, AT, CC. The first sequence
# supplies the first allele of every state and the second the second, which is
# what makes the states ordered (phased).
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

# The normalized 16x16 GT16 rate matrix. Built once per call to derive the
# event rates (and used as the reference the tests check against); the
# simulation itself never exponentiates it.
.gt16_rate_matrix <- function(rates, pi) {
  # symmetric nucleotide exchangeabilities; the lower triangle taken in column
  # order is exactly the input order AC, AG, AT, CG, CT, GT
  R <- matrix(0, 4, 4, dimnames = list(.gt16_nucleotides, .gt16_nucleotides))
  R[lower.tri(R)] <- rates
  R <- R + t(R)

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
  diag(Q) <- -rowSums(Q)
  total_rate <- as.vector(pi %*% diag(Q))
  if (total_rate >= 0) {
    # every state carrying equilibrium mass is a dead end, so there is no
    # timescale to normalize against
    stop("`rates` and `pi` give a degenerate GT16 rate matrix with no substitutions.", call. = FALSE)
  }
  (-1 / total_rate) * Q
}

# The event table driving the simulation: for each state, the six states
# reachable by a single-allele substitution (`dest`) and their rates, read
# straight off the off-diagonals of the normalized Q so that the mechanistic
# simulation and the matrix model cannot drift apart. `cumulative` makes
# drawing an event one comparison per site, and its last column doubles as the
# total outgoing rate (sharing the same floating-point sum keeps the draw
# inside the six events).
.gt16_events <- function(rates, pi) {
  Q <- .gt16_rate_matrix(rates, pi)
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

  rate <- matrix(Q[cbind(rep(seq_len(n), 6L), as.vector(dest))], n, 6L)
  cumulative <- t(apply(rate, 1, cumsum))
  list(Q = Q, dest = dest, rate = rate, cumulative = cumulative,
       total = cumulative[, 6L])
}

# Gillespie simulation of the GT16 chain along one branch, run for all sites at
# once: sites still "in flight" draw their next waiting time together, those
# whose next event falls past the end of the branch drop out, and the rest take
# one single-allele substitution and go round again. Sites that never have an
# event simply inherit the parent genotype.
.evolve_branch_gt16 <- function(state, branch_time, events) {
  if (branch_time <= 0) return(state)

  elapsed <- numeric(length(state))
  active <- seq_along(state)

  while (length(active) > 0) {
    s <- state[active]
    # a state with no outgoing rate gets an infinite waiting time and so simply
    # drops out below
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
