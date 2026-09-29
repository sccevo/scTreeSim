# Lifetimes of simulated particles (branch lengths of divided/dead particles
# and the root edge) follow the per-type gamma lifetime distribution.
# See helper-lifetimes.R for why a censoring-aware (log-rank) test is used.
#
# Origin simulator: death_prob = 0 and min_taxa = 1, so no tree is discarded. Discarding
# trees without enough alive tips would condition on survival, which depends
# on the lifetimes and would bias the test.
# Ntaxa simulator: d > 0 is fine; whether ntaxa particles are reached depends
# only on the sequence of death/division decisions, not on the lifetimes.

a <- c(0.5, 1)  # scale per type: mean lifetime 1 for both types
b <- c(2, 1)    # shape per type
r <- 0.2
Xi_as <- matrix(c(0, 0, r, 0), 2)     # type 0 -> asymmetric division into (0, 1)
Xi_s <- matrix(c(1 - r, 0, 0, 1), 2)

test_that("origin simulator: lifetimes follow the lifetime distribution (single type)", {
  withr::local_seed(1)
  trees <- simulate_trees(150, function() sim_adb_origin_complete_fast(5, scale = a[1], shape = b[1], death_prob = 0, min_taxa = 1))
  lifetimes <- do.call(rbind, lapply(trees, lifetime_data, censor_alive_at = "tip"))

  expect_gt(sum(lifetimes$event), 5000)
  expect_gt(logrank_p(lifetimes, a[1], b[1]), 0.001)
  # the test has power against a wrong lifetime distribution
  expect_lt(logrank_p(lifetimes, 1.25 * a[1], b[1]), 0.001)
})

test_that("origin simulator: lifetimes follow the per-type lifetime distributions (multi-type)", {
  withr::local_seed(2)
  trees <- simulate_trees(150, function() {
    sim_adb_origin_complete_fast(5, scale = a, shape = b, death_prob = c(0, 0), asym_trans_prob = Xi_as, sym_trans_prob = Xi_s, min_taxa = 1)
  })
  lifetimes <- do.call(rbind, lapply(trees, lifetime_data, censor_alive_at = "tip"))

  expect_gt(sum(lifetimes$event & lifetimes$type == 1), 1000)
  expect_gt(logrank_p(lifetimes, a, b), 0.001)
  # type 1 has exponential lifetimes (shape 1); testing it against type 0's distribution must fail
  expect_lt(logrank_p(lifetimes, a[c(1, 1)], b[c(1, 1)]), 0.001)
})

test_that("ntaxa simulator: lifetimes follow the lifetime distribution (single type)", {
  withr::local_seed(4)
  trees <- simulate_trees(150, function() sim_adb_ntaxa_complete_fast(50, scale = a[1], shape = b[1], death_prob = 0.1))
  lifetimes <- do.call(rbind, lapply(trees, lifetime_data, censor_alive_at = "last_division"))

  expect_gt(sum(lifetimes$event), 5000)
  expect_gt(logrank_p(lifetimes, a[1], b[1]), 0.001)
  expect_lt(logrank_p(lifetimes, 1.25 * a[1], b[1]), 0.001)
})

test_that("ntaxa simulator: lifetimes follow the per-type lifetime distributions (multi-type)", {
  withr::local_seed(5)
  trees <- simulate_trees(150, function() {
    sim_adb_ntaxa_complete_fast(50, scale = a, shape = b, death_prob = c(0.1, 0.1), asym_trans_prob = Xi_as, sym_trans_prob = Xi_s)
  })
  lifetimes <- do.call(rbind, lapply(trees, lifetime_data, censor_alive_at = "last_division"))

  expect_gt(sum(lifetimes$event & lifetimes$type == 1), 1000)
  expect_gt(logrank_p(lifetimes, a, b), 0.001)
  expect_lt(logrank_p(lifetimes, a[c(1, 1)], b[c(1, 1)]), 0.001)
})

test_that("ntaxa simulator: the root edge is a complete lifetime", {
  # the root always divides when ntaxa >= 2, and reaching ntaxa does not depend
  # on the root's lifetime, so root edges are a plain sample of lifetimes
  withr::local_seed(6)
  trees <- simulate_trees(300, function() sim_adb_ntaxa_complete_fast(10, scale = a[1], shape = b[1], death_prob = 0.1))
  root_edges <- vapply(trees, function(tree) tree@phylo$root.edge, numeric(1))

  expect_gt(stats::ks.test(root_edges, "pgamma", shape = b[1], scale = a[1])$p.value, 0.001)
  expect_lt(stats::ks.test(root_edges, "pgamma", shape = b[1], scale = 1.5 * a[1])$p.value, 0.001)
})
