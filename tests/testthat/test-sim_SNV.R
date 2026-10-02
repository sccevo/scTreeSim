# GT16 diploid nucleotide SNVs

withr::with_seed(1, {
  tree <- sim_adb_origin_samp(5, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
})
tree_df <- tibble::as_tibble(tree)
split_sites <- function(sequences) strsplit(sequences, "_", fixed = TRUE)

# a deliberately uneven parameterisation, so that symmetry cannot hide a
# mis-indexed rate or frequency
rates <- c(0.5, 1.5, 0.7, 1.2, 2.3, 0.9)   # AC AG AT CG CT GT
pi <- (1:16) / sum(1:16)

# exp(Qt) by symmetrising the reversible Q, used only as the reference the
# mechanistic simulator is checked against
transition_probs <- function(Q, t, pi) {
  d <- sqrt(pi)
  S <- diag(d) %*% Q %*% diag(1 / d)
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  diag(1 / d) %*% e$vectors %*% diag(exp(e$values * t)) %*% t(e$vectors) %*% diag(d)
}

test_that("the 16 GT16 states are the ordered nucleotide pairs", {
  expect_equal(scTreeSim:::.gt16_states,
               c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                 "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"))
  # AG and GA are distinct states, i.e. the model is phased
  expect_equal(length(unique(scTreeSim:::.gt16_states)), 16)
})

test_that("only one allele changes in an individual event", {
  dest <- scTreeSim:::.gt16_neighbours()
  states <- scTreeSim:::.gt16_states

  expect_equal(dim(dest), c(16L, 6L))
  for (i in seq_len(16)) {
    from <- strsplit(states[i], "")[[1]]
    to <- strsplit(states[dest[i, ]], "")
    expect_equal(length(unique(dest[i, ])), 6L)
    expect_true(all(vapply(to, function(g) sum(g != from), integer(1)) == 1L))
  }

  # the worked example from the documentation
  expect_setequal(states[dest[match("AG", states), ]],
                  c("CG", "GG", "TG", "AA", "AC", "AT"))

  # a two-allele change is not an instantaneous event anywhere in Q
  Q <- scTreeSim:::.gt16_rate_matrix(rates, pi)
  expect_equal(Q["AG", "CT"], 0)
  for (i in seq_len(16)) {
    off <- setdiff(seq_len(16), c(i, dest[i, ]))
    expect_true(all(Q[i, off] == 0))
  }
})

test_that("event rates equal the off-diagonal entries of the normalized GT16 Q", {
  Q <- scTreeSim:::.gt16_rate_matrix(rates, pi)
  events <- scTreeSim:::.gt16_events(rates, pi)
  states <- scTreeSim:::.gt16_states

  for (i in seq_len(16)) {
    expect_equal(events$rate[i, ], Q[i, events$dest[i, ]], ignore_attr = TRUE)
    expect_equal(events$total[i], -Q[i, i])
  }

  # the normalization gives exactly one expected substitution per unit time
  expect_equal(sum(pi * events$total), 1)

  # and the off-diagonals follow Q_ij = R_ij * pi_j, up to that common factor:
  # AG -> AC is a G->C change of the second allele, AG -> AA a G->A change
  expect_equal(Q["AG", "AC"] / Q["AG", "AA"],
               (rates[4] * pi[match("AC", states)]) / (rates[2] * pi[match("AA", states)]))
  # pi is the stationary distribution of Q
  expect_equal(as.vector(pi %*% Q), rep(0, 16), tolerance = 1e-12)
  # and Q is reversible with respect to pi
  expect_equal(pi * Q, t(pi * Q), ignore_attr = TRUE)
})

test_that("a zero-length branch leaves every genotype unchanged", {
  events <- scTreeSim:::.gt16_events(rates, pi)
  state <- withr::with_seed(1, sample.int(16, 50, replace = TRUE))

  expect_equal(scTreeSim:::.evolve_branch_gt16(state, 0, events), state)
  expect_equal(scTreeSim:::.evolve_branch_gt16(state, -1, events), state)

  # over a very short branch almost every site is inherited untouched
  short <- withr::with_seed(2, scTreeSim:::.evolve_branch_gt16(state, 1e-4, events))
  expect_true(mean(short == state) > 0.95)
})

test_that("multiple substitutions can occur on a single branch", {
  events <- scTreeSim:::.gt16_events(rates, pi)
  # AG -> CT needs at least two events, since both alleles differ
  ends <- withr::with_seed(3, scTreeSim:::.evolve_branch_gt16(rep(match("AG", scTreeSim:::.gt16_states), 2000), 2, events))
  both_changed <- scTreeSim:::.gt16_allele1[ends] != 1L & scTreeSim:::.gt16_allele2[ends] != 3L
  expect_true(any(both_changed))
  expect_true(any(scTreeSim:::.gt16_states[ends] == "CT"))
})

test_that("mechanistic simulation reproduces the GT16 transition probabilities", {
  events <- scTreeSim:::.gt16_events(rates, pi)
  n <- 20000

  for (from in c("AA", "AG")) {
    for (t in c(0.2, 1.5)) {
      start <- match(from, scTreeSim:::.gt16_states)
      ends <- withr::with_seed(4, scTreeSim:::.evolve_branch_gt16(rep(start, n), t, events))
      observed <- table(factor(ends, levels = seq_len(16)))
      expected <- transition_probs(events$Q, t, pi)[start, ]

      # drop states too rare for the chi-squared approximation, renormalising
      keep <- expected * n > 10
      expect_gt(stats::chisq.test(observed[keep], p = expected[keep] / sum(expected[keep]))$p.value, 0.001)
    }
  }
})

test_that("the two aligned sequences are paired site by site", {
  pair <- scTreeSim:::.gt16_pair_sequences(c("CGTAAC", "CGTATC"))
  expect_equal(scTreeSim:::.gt16_states[pair], c("CC", "GG", "TT", "AA", "AT", "CC"))

  # the first sequence gives the first allele, so the states stay ordered
  expect_equal(scTreeSim:::.gt16_states[scTreeSim:::.gt16_pair_sequences(c("AG", "GA"))], c("AG", "GA"))
  expect_equal(scTreeSim:::.gt16_pair_sequences(c("acgt", "ACGT")),
               scTreeSim:::.gt16_pair_sequences(c("ACGT", "ACGT")))

  expect_error(scTreeSim:::.gt16_pair_sequences("ACGT"), "two aligned allele sequences")
  expect_error(scTreeSim:::.gt16_pair_sequences(c("ACGT", "ACG")), "must be aligned")
  expect_error(scTreeSim:::.gt16_pair_sequences(c("ACGN", "ACGT")), "only the nucleotides")
  expect_error(scTreeSim:::.gt16_pair_sequences(c("", "")), "non-empty")
})

test_that("a root given as two sequences starts the simulation from that genotype", {
  out <- sim_gt16_nt_diploid(tree, sequences = c("CGTAAC", "CGTATC"), rates = rates, pi = pi,
                             clock_rate = 1e-12)
  expect_true(all(out$genotype == "CC_GG_TT_AA_AT_CC"))

  # stating the length alongside the sequences is allowed, and changes nothing
  expect_equal(sim_gt16_nt_diploid(tree, sequences = c("CGTAAC", "CGTATC"), rates = rates, pi = pi,
                                   l = 6, clock_rate = 1e-12), out)
})

test_that("sim_gt16_nt_diploid returns one genotype sequence of l sites per node", {
  withr::local_seed(5)
  out <- sim_gt16_nt_diploid(tree, l = 8, rates = rates, pi = pi)

  expect_named(out, c("node", "genotype", "observed"))
  expect_equal(out$node, tree_df$node)
  sites <- split_sites(out$genotype)
  expect_true(all(lengths(sites) == 8))
  expect_true(all(unlist(sites) %in% scTreeSim:::.gt16_states))
  # without an error model the observations are the true genotypes
  expect_equal(out$observed, out$genotype)
})

test_that("genotypes are inherited from the parent and a slow clock keeps them", {
  withr::local_seed(6)
  out <- sim_gt16_nt_diploid(tree, sequences = c(strrep("A", 200), strrep("G", 200)),
                             rates = rates, pi = pi, clock_rate = 1e-5)
  sites <- split_sites(out$genotype)

  expect_true(all(vapply(sites, function(s) mean(s == "AG"), numeric(1)) > 0.95))

  # with no root edge and a zero clock every node keeps the root sequence
  frozen <- sim_gt16_nt_diploid(tree, sequences = c("ACTGA", "AGTAC"), rates = rates, pi = pi,
                                clock_rate = 1e-12)
  expect_true(all(frozen$genotype == "AA_CG_TT_GA_AC"))
})

test_that("input is validated", {
  expect_error(sim_gt16_nt_diploid(tree, l = 3, rates = rep(1, 5)), "`rates` must be six")
  expect_error(sim_gt16_nt_diploid(tree, l = 3, rates = rep(0, 6)), "`rates` must be six")
  expect_error(sim_gt16_nt_diploid(tree, l = 3, pi = rep(0.5, 16)), "summing to 1")
  expect_error(sim_gt16_nt_diploid(tree, l = 3, pi = rep(0.2, 5)), "must have length 16")
  expect_error(sim_gt16_nt_diploid(tree, l = 3, epsilon = 1.2), "must be probabilities")
  expect_error(sim_gt16_nt_diploid(tree), "Supply `sequences`")
  expect_error(sim_gt16_nt_diploid(tree, sequences = c("ACG", "AC")), "must be aligned")
  expect_error(sim_gt16_nt_diploid(tree, sequences = c("AC", "GT"), l = 3), "must match the length")
  # all equilibrium mass on one state leaves no substitutions to normalize by
  expect_error(sim_gt16_nt_diploid(tree, l = 3, pi = c(1, rep(0, 15))), "degenerate GT16 rate matrix")
})

test_that("nucleotide frequencies expand to the product genotype frequencies", {
  nuc <- c(0.1, 0.2, 0.3, 0.4)
  expanded <- scTreeSim:::.gt16_pi(nuc)
  expect_equal(sum(expanded), 1)
  expect_equal(expanded[match("AG", scTreeSim:::.gt16_states)], nuc[1] * nuc[3])
  expect_equal(expanded[match("TC", scTreeSim:::.gt16_states)], nuc[4] * nuc[2])
})

# --- error model (Chen et al. 2022, eq. 2) ------------------------------------

test_that("the GT16 error matrix follows equation (2) and is a proper distribution", {
  epsilon <- 0.3
  delta <- 0.2
  M <- scTreeSim:::.gt16_error_matrix(epsilon, delta)

  expect_equal(rowSums(M), rep(1, 16), ignore_attr = TRUE)
  expect_true(all(M >= 0))

  # true homozygous AA
  expect_equal(M["AA", "AA"], 1 - epsilon + 0.5 * delta * epsilon)
  expect_equal(M["AA", "AC"], (1 - delta) * epsilon / 6)
  expect_equal(M["AA", "CA"], M["AA", "AC"])
  expect_equal(M["AA", "CC"], delta * epsilon / 6)
  expect_equal(M["AA", "CG"], 0)

  # true heterozygous AC
  expect_equal(M["AC", "AC"], (1 - delta) * (1 - epsilon))
  expect_equal(M["AC", "AA"], 0.5 * delta + epsilon / 6 - delta * epsilon / 3)
  expect_equal(M["AC", "CC"], M["AC", "AA"])
  expect_equal(M["AC", "GG"], delta * epsilon / 6)
  expect_equal(M["AC", "AG"], (1 - delta) * epsilon / 6)   # P(ac|ab)
  expect_equal(M["AC", "GC"], M["AC", "AG"])               # P(cb|ab)
  expect_equal(M["AC", "CA"], 0)                           # phased: not an error outcome

  # no error collapses to the identity
  expect_equal(scTreeSim:::.gt16_error_matrix(0, 0), diag(16), ignore_attr = TRUE)
})

test_that("the error model changes observations but not the true genotypes", {
  args <- list(tree, sequences = c(strrep("A", 100), strrep("A", 100)), rates = rates, pi = pi)
  with_error <- withr::with_seed(7, do.call(sim_gt16_nt_diploid, c(args, epsilon = 0.3, delta = 0.2)))
  no_error <- withr::with_seed(7, do.call(sim_gt16_nt_diploid, c(args, epsilon = 0, delta = 0)))

  # the evolutionary simulation is untouched by the error parameters
  expect_equal(with_error$genotype, no_error$genotype)

  is_tip <- with_error$node %in% tree_df$node[tree_df$status == 1]
  expect_true(any(with_error$observed[is_tip] != with_error$genotype[is_tip]))
  # errors are a measurement effect at the sequenced cells only
  expect_equal(with_error$observed[!is_tip], with_error$genotype[!is_tip])
})

test_that("observed genotypes at the tips follow the error distribution", {
  epsilon <- 0.3
  delta <- 0.2
  out <- withr::with_seed(8, sim_gt16_nt_diploid(tree, sequences = c(strrep("A", 3000), strrep("C", 3000)),
                                                 rates = rates, pi = pi, clock_rate = 1e-12,
                                                 epsilon = epsilon, delta = delta))
  tip <- out[out$node %in% tree_df$node[tree_df$status == 1], ][1, ]
  expect_equal(unique(split_sites(tip$genotype)[[1]]), "AC")

  observed <- table(factor(split_sites(tip$observed)[[1]], levels = scTreeSim:::.gt16_states))
  expected <- scTreeSim:::.gt16_error_matrix(epsilon, delta)["AC", ]
  keep <- expected * 3000 > 10
  expect_gt(stats::chisq.test(observed[keep], p = expected[keep] / sum(expected[keep]))$p.value, 0.001)
})
