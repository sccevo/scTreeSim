# Binary substitution model

# Small tree from helper-barcodes.R: origin -0.5- root(4) -0.5- node 5 -1- a(1), b(2); root -1.5- c(3).
# Distances from the origin: root 0.5, node 5 1.0, tips 2.0; between tips: a-b 2, a-c 3.
tree <- small_barcode_tree()

# node x site matrix of states ("0"/"1") for one copy
site_matrix <- function(col) do.call(rbind, strsplit(col, "", fixed = TRUE))

# closed form for the normalised 2-state model: pi = (lambda, 1) / (lambda + 1),
# relaxation rate r = (1 + lambda)^2 / (2 * lambda) (1 expected substitution per unit distance)
pi_one <- function(lambda) 1 / (lambda + 1)
relax_rate <- function(lambda) (1 + lambda)^2 / (2 * lambda)
p_01 <- function(lambda, d) pi_one(lambda) * (1 - exp(-relax_rate(lambda) * d))
p_mismatch <- function(lambda, d) {
  2 * pi_one(lambda) * (1 - pi_one(lambda)) * (1 - exp(-relax_rate(lambda) * d))
}

test_that("sim_binary returns one string of l sites per node and copy", {
  withr::local_seed(1)
  out <- sim_binary(tree, n_sites = 7, n_sequences = 2)

  expect_named(out, c("node", "binary_1", "binary_2"))
  expect_equal(out$node, tibble::as_tibble(tree)$node)
  strings <- c(out$binary_1, out$binary_2)
  expect_true(all(nchar(strings) == 7))
  expect_true(all(grepl("^[01]+$", strings)))
})

test_that("the root sequence is drawn from the equilibrium of each copy's lambda", {
  withr::local_seed(2)
  l <- 4000
  out <- sim_binary(tree, n_sites = l, n_sequences = 2, lambda = c(0.5, 3))

  for (k in 1:2) {
    p <- pi_one(c(0.5, 3)[k])
    root_freq <- mean(site_matrix(out[[paste0("binary_", k)]])[4, ] == "1")
    expect_within(root_freq, p, tolerance_4se(p, l))
  }
})

test_that("a fixed root sequence diverges as the closed-form P(d), with d = rate * distance from origin", {
  withr::local_seed(3)
  l <- 4000
  lambda <- 2
  rate <- 0.4
  out <- sim_binary(tree, n_sites = l, lambda = lambda, origin_seq = rep(0L, l), clock = clock_strict(rate))
  states <- site_matrix(out$binary_1)

  # the origin-to-root edge is evolved too
  dist_from_origin <- c(2, 2, 2, 0.5, 1)
  for (nd in 1:5) {
    p <- p_01(lambda, rate * dist_from_origin[nd])
    expect_within(mean(states[nd, ] == "1"), p, tolerance_4se(p, l))
  }
})

test_that("tip divergence at equilibrium follows the closed form of the path distance", {
  withr::local_seed(4)
  l <- 4000
  lambda <- 0.5
  out <- sim_binary(tree, n_sites = l, lambda = lambda, clock = clock_strict(0.3))
  states <- site_matrix(out$binary_1)

  p_ab <- p_mismatch(lambda, 0.3 * 2)
  p_ac <- p_mismatch(lambda, 0.3 * 3)
  expect_within(mean(states[1, ] != states[2, ]), p_ab, tolerance_4se(p_ab, l))
  expect_within(mean(states[1, ] != states[3, ]), p_ac, tolerance_4se(p_ac, l))
})

test_that("clock = NULL simulates on the distance column of the tree", {
  withr::local_seed(5)
  l <- 4000
  lambda <- 2
  ct <- clock_tree(tree, clock_strict(0.4))
  out <- sim_binary(ct, n_sites = l, lambda = lambda, origin_seq = rep(0L, l), clock = NULL)
  p <- p_01(lambda, 0.4 * 2)
  expect_within(mean(site_matrix(out$binary_1)[1, ] == "1"), p, tolerance_4se(p, l))

  expect_error(sim_binary(tree, n_sites = 5, clock = NULL), "distance")

  # a random clock is drawn once, so passing it equals drawing it with clock_tree first
  clk <- clock_relaxed(sd = 0.8)
  direct <- withr::with_seed(1, sim_binary(tree, n_sites = 20, n_sequences = 2, clock = clk))
  stored <- withr::with_seed(1, sim_binary(clock_tree(tree, clk), n_sites = 20, n_sequences = 2, clock = NULL))
  expect_equal(direct, stored)
})

test_that("every node of a sampled tree inherits its parent's sequence when nothing evolves", {
  withr::local_seed(6)
  adb <- sim_adb_origin_samp(10, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
  origin_seq <- sample(0:1, 20, replace = TRUE)
  out <- sim_binary(adb, n_sites = 20, origin_seq = origin_seq, clock = clock_strict(0))

  expect_equal(nrow(out), nrow(tibble::as_tibble(adb)))
  expect_true(all(out$binary_1 == paste(origin_seq, collapse = "")))
})

test_that("alpha and beta flip tip states only", {
  withr::local_seed(7)
  l <- 4000
  alpha <- 0.2
  beta <- 0.1
  # no evolution, so the true state everywhere is the root state
  sim <- function(origin_state) {
    out <- sim_binary(tree, n_sites = l, origin_seq = rep(origin_state, l), clock = clock_strict(0), alpha = alpha, beta = beta)
    site_matrix(out$binary_1)
  }

  zeros <- sim(0L)
  expect_true(all(zeros[4:5, ] == "0"))
  for (tip in 1:3) expect_within(mean(zeros[tip, ] == "1"), beta, tolerance_4se(beta, l))

  ones <- sim(1L)
  expect_true(all(ones[4:5, ] == "1"))
  for (tip in 1:3) expect_within(mean(ones[tip, ] == "0"), alpha, tolerance_4se(alpha, l))
})

test_that("invalid inputs are rejected", {
  expect_error(sim_binary(tree@phylo, n_sites = 5))
  expect_error(sim_binary(tree, n_sites = 0))
  expect_error(sim_binary(tree, n_sites = 5, lambda = 0))
  expect_error(sim_binary(tree, n_sites = 5, n_sequences = 2, lambda = c(1, 2, 3)))
  expect_error(sim_binary(tree, n_sites = 5, origin_seq = c(0L, 1L)))
  expect_error(sim_binary(tree, n_sites = 2, origin_seq = c(1L, 2L)))
  expect_error(sim_binary(tree, n_sites = 5, clock = 1))
  expect_error(sim_binary(tree, n_sites = 5, alpha = 1.5))
  expect_error(sim_binary(tree, n_sites = 5, beta = -0.1))
})
