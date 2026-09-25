# Sampling of the children's types upon division

# parent type 0 never self-renews: symmetric to 0 has probability 0
Xi_s <- matrix(c(0, 0, 0, 0.3, 1, 0, 0.1, 0, 1), 3)
Xi_as <- matrix(c(0, 0, 0, 0.4, 0, 0, 0.2, 0, 0), 3)

test_that("sample_types draws outcomes with the probabilities given by Xi_s and Xi_as", {
  withr::local_seed(1)
  draws <- replicate(10000, paste(sample_types(0, Xi_as, Xi_s), collapse = "_"))

  # symmetric outcomes (i, i); asymmetric outcomes (0, i) or (i, 0) with equal probability
  expected <- c("1_1" = 0.3, "2_2" = 0.1, "0_1" = 0.2, "1_0" = 0.2, "0_2" = 0.1, "2_0" = 0.1)
  expect_setequal(unique(draws), names(expected)) # the impossible outcome "0_0" never occurs
  observed <- table(factor(draws, levels = names(expected)))
  expect_gt(stats::chisq.test(observed, p = expected)$p.value, 0.001)
})

test_that("sample_types always returns valid types when rows sum to 1 up to rounding error", {
  withr::local_seed(2)
  Xi_s_rounded <- Xi_s
  Xi_s_rounded[1, 3] <- Xi_s_rounded[1, 3] - 1e-12
  draws <- replicate(1000, sample_types(0, Xi_as, Xi_s_rounded))

  expect_true(all(draws %in% 0:2))
  expect_false(any(draws[1, ] == 0 & draws[2, ] == 0))
})

test_that("sample_types is deterministic for a certain outcome", {
  expect_equal(sample_types(1, Xi_as, Xi_s), c(1, 1))
  expect_equal(sample_types(2, Xi_as, Xi_s), c(2, 2))
})
