# Binary substitution model

# Small tree from helper-barcodes.R: origin -0.5- root(4) -0.5- node 5 -1- a(1), b(2); root -1.5- c(3).
# Distances from the origin: root 0.5, node 5 1.0, tips 2.0; between tips: a-b 2, a-c 3.
tree <- small_barcode_tree()

# node x site matrix of states ("0"/"1") for one copy
site_matrix <- function(col) do.call(rbind, strsplit(col, "_", fixed = TRUE))

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
  out <- sim_binary(tree, l = 7, m = 2)

  expect_named(out, c("node", "binary_1", "binary_2"))
  expect_equal(out$node, tibble::as_tibble(tree)$node)
  sites <- strsplit(c(out$binary_1, out$binary_2), "_", fixed = TRUE)
  expect_true(all(lengths(sites) == 7))
  expect_true(all(unlist(sites) %in% c("0", "1")))
})

test_that("the root sequence is drawn from the equilibrium of each copy's lambda", {
  withr::local_seed(2)
  l <- 4000
  out <- sim_binary(tree, l = l, m = 2, lambda_vec = c(0.5, 3))

  for (k in 1:2) {
    p <- pi_one(c(0.5, 3)[k])
    root_freq <- mean(site_matrix(out[[paste0("binary_", k)]])[4, ] == "1")
    expect_within(root_freq, p, tolerance_4se(p, l))
  }
})

test_that("a fixed root sequence diverges as the closed-form P(d), with d = clock_rate * distance from origin", {
  withr::local_seed(3)
  l <- 4000
  lambda <- 2
  rate <- 0.4
  out <- sim_binary(tree, l = l, lambda = lambda, root_seq = rep(1L, l), clock_rate = rate)
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
  out <- sim_binary(tree, l = l, lambda = lambda, clock_rate = 0.3)
  states <- site_matrix(out$binary_1)

  p_ab <- p_mismatch(lambda, 0.3 * 2)
  p_ac <- p_mismatch(lambda, 0.3 * 3)
  expect_within(mean(states[1, ] != states[2, ]), p_ab, tolerance_4se(p_ab, l))
  expect_within(mean(states[1, ] != states[3, ]), p_ac, tolerance_4se(p_ac, l))
})

test_that("alpha and beta flip tip states only", {
  withr::local_seed(5)
  l <- 4000
  alpha <- 0.2
  beta <- 0.1
  # effectively no evolution, so the true state everywhere is the root state
  sim <- function(root_state) {
    out <- sim_binary(tree, l = l, root_seq = rep(root_state, l), clock_rate = 1e-12, alpha = alpha, beta = beta)
    site_matrix(out$binary_1)
  }

  zeros <- sim(1L)
  expect_true(all(zeros[4:5, ] == "0"))
  for (tip in 1:3) expect_within(mean(zeros[tip, ] == "1"), beta, tolerance_4se(beta, l))

  ones <- sim(2L)
  expect_true(all(ones[4:5, ] == "1"))
  for (tip in 1:3) expect_within(mean(ones[tip, ] == "0"), alpha, tolerance_4se(alpha, l))
})

test_that("invalid inputs are rejected", {
  expect_error(sim_binary(tree@phylo, l = 5))
  expect_error(sim_binary(tree, l = 0))
  expect_error(sim_binary(tree, l = 5, lambda = 0))
  expect_error(sim_binary(tree, l = 5, m = 2, lambda_vec = 1))
  expect_error(sim_binary(tree, l = 5, root_seq = c(1L, 2L)))
  expect_error(sim_binary(tree, l = 2, root_seq = c(0L, 1L)))
  expect_error(sim_binary(tree, l = 5, alpha = 1.5))
  expect_error(sim_binary(tree, l = 5, beta = -0.1))
})
