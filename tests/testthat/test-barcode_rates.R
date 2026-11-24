# Time profiles and barcode_rate() of the generic barcode simulator

test_that("rate profiles give the multiplier, breakpoints and a valid bound", {
  w <- rate_window(1, 2)
  expect_equal(w$g(c(0.5, 1, 1.5, 2)), c(0, 1, 1, 0))
  expect_equal(w$breaks, c(1, 2))
  expect_equal(w$max(0, 1), 0)
  expect_equal(w$max(1, 2), 1)
  d <- rate_exp_decay(half_life = 2)
  expect_equal(d$g(c(0, 2, 4)), c(1, 0.5, 0.25))
  expect_equal(d$max(1, 3), d$g(1))
  pw <- rate_piecewise(c(1, 3), c(0, 2, 1))
  expect_equal(pw$g(c(0.5, 1, 2, 3, 5)), c(0, 2, 2, 1, 1))
  expect_error(rate_piecewise(c(1, 3), c(1, 2)))
  expect_error(rate_window(2, 1))
  expect_error(rate_exp_decay(0))
  expect_error(barcode_rate(-1))
  expect_error(barcode_rate(1, time = "a"))
})

test_that("rate_product multiplies profiles, pools breakpoints and bounds each piece", {
  prod_rate <- rate_product(rate_window(1, 2), rate_exp_decay(half_life = 0.5))
  expect_equal(prod_rate$breaks, c(1, 2))
  expect_equal(prod_rate$g(c(0.5, 1, 1.5, 2)), c(0, 0.25, 0.125, 0))
  expect_equal(prod_rate$max(1, 2), 0.25)  # the decay at the start of the piece
  expect_equal(prod_rate$max(0, 1), 0)
  expect_equal(rate_product(rate_window(1, 3), rate_window(2, 4))$breaks, 1:4)
  expect_error(rate_product())
  expect_error(rate_product(1))
})

test_that("closed-form integrals of the profiles match numerical integration", {
  t0 <- c(-1, 0, 0.3, 1.2, 2, 5)
  t1 <- c(0.5, 2.5, 1.7, 1.2, 4, 5)
  numerical <- function(g, a, b) {
    if (b <= a) return(0)
    stats::integrate(Vectorize(g), a, b, subdivisions = 1000L, rel.tol = 1e-10)$value
  }
  profiles <- list(
    rate_constant(), rate_window(1, 2), rate_window(1, 1), rate_exp_decay(0.7),
    rate_piecewise(c(1, 3), c(0, 2, 1)),
    rate_product(rate_window(0.5, 3), rate_exp_decay(1.5)),
    rate_product(rate_piecewise(1, c(2, 0.5)), rate_exp_decay(1), rate_exp_decay(2), rate_window(0, 2.2))
  )
  for (p in profiles) {
    expected <- mapply(function(a, b) numerical(p$g, a, b), t0, t1)
    expect_equal(scTreeSim:::.profile_integral(p$pexp, t0, t1), expected, tolerance = 1e-6)
  }
  # rate_fn has no closed form, also not inside a product
  expect_null(rate_fn(function(t) 1, max = 1)$pexp)
  expect_null(rate_product(rate_window(0, 1), rate_fn(function(t) 1, max = 1))$pexp)
})

test_that("rate_fn (thinning) and the equivalent closed-form profile give the same distribution", {
  skip_on_cran()
  tr <- small_barcode_tree()
  n_rep <- 3000
  closed <- barcode_rate(1.2, rate_product(rate_window(0.4, 1.7), rate_exp_decay(1)))
  thinned <- barcode_rate(1.2, rate_product(rate_window(0.4, 1.7), rate_fn(function(t) 2^-t, max = 1)))
  run <- function(rate, sequential) {
    replicate_states(sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 2, edit_probs = c(A = 0.5, B = 0.5),
                                         edit_rate = rate, silencing_rate = 0.2, sequential = sequential))
  }
  withr::local_seed(41)
  for (sequential in c(TRUE, FALSE)) {
    expect_gt(homogeneity_p(joint_counts(run(closed, sequential), 1:2),
                            joint_counts(run(thinned, sequential), 1:2)), 0.001)
  }
})

test_that("a rate_fn exceeding its stated maximum is an error", {
  tr <- small_barcode_tree()
  expect_error(
    sim_barcode_generic(tr, n_barcodes = 20, n_sites = 1, edit_probs = 1,
                        edit_rate = barcode_rate(5, rate_fn(function(t) 3, max = 1))),
    "exceeds"
  )
})
