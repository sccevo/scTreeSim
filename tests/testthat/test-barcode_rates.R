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

test_that("a rate_fn exceeding its stated maximum is an error", {
  tr <- small_barcode_tree()
  expect_error(
    sim_barcode_generic(tr, n_barcodes = 20, n_sites = 1, edit_probs = 1,
                        edit_rate = barcode_rate(5, rate_fn(function(t) 3, max = 1))),
    "exceeds"
  )
})
