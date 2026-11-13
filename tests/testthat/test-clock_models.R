# Clock models and per-branch evolutionary distances

# Small tree from helper-barcodes.R: nodes a(1), b(2), c(3), root(4), internal node 5.
# Time lengths of the branches ending in nodes 1..5 (root: root edge).
tree <- small_barcode_tree()
lengths <- c(1, 1, 1.5, 0.5, 0.5)

# distances drawn by a clock, indexed by node
distances <- function(clock, tr = tree) .node_distances(clock_tree(tr, clock))

# `n_rep` draws of the per-branch multipliers of a clock, pooled
draw_multipliers <- function(clock, n_rep = 4000) {
  as.vector(replicate(n_rep, distances(clock) / (clock$rate * lengths)))
}

# ---- clock_strict ----

test_that("clock_strict validates its rate", {
  expect_s3_class(clock_strict(0.3), "clock_model")
  expect_equal(clock_strict()$rate, 1)
  for (bad in list(-1, NA_real_, Inf, c(1, 2), "1")) expect_error(clock_strict(bad), "rate")
})

test_that("strict distances are rate * branch length", {
  expect_equal(distances(clock_strict()), lengths)
  expect_equal(distances(clock_strict(0.4)), 0.4 * lengths)
  expect_equal(distances(clock_strict(0)), rep(0, 5))
})

# ---- clock_relaxed ----

test_that("clock_relaxed takes exactly the parameter of its distribution", {
  expect_equal(clock_relaxed(0.3, sd = 0.5)[c("rate", "dist", "sd")], list(rate = 0.3, dist = "lognormal", sd = 0.5))
  expect_equal(clock_relaxed(dist = "gamma", shape = 2)[c("dist", "shape")], list(dist = "gamma", shape = 2))
  expect_equal(clock_relaxed(dist = "exponential")$dist, "exponential")

  expect_error(clock_relaxed(), "`sd` only")
  expect_error(clock_relaxed(sd = 0.5, shape = 2), "`sd` only")
  expect_error(clock_relaxed(dist = "gamma"), "`shape` only")
  expect_error(clock_relaxed(dist = "gamma", shape = 2, sd = 0.5), "`shape` only")
  expect_error(clock_relaxed(dist = "exponential", sd = 1), "neither")
  expect_error(clock_relaxed(dist = "cauchy", sd = 0.5))

  expect_error(clock_relaxed(-1, sd = 0.5), "rate")
  for (bad in list(-0.1, NA_real_, Inf, c(0.1, 0.2), "0.5")) expect_error(clock_relaxed(sd = bad), "sd")
  for (bad in list(0, -1, Inf)) expect_error(clock_relaxed(dist = "gamma", shape = bad), "shape")
})

test_that("relaxed distances are rate * branch length * one iid multiplier per branch", {
  sd <- 0.5
  expected <- list(
    lognormal   = function() exp(sd * stats::rnorm(5) - sd^2 / 2),
    gamma       = function() stats::rgamma(5, shape = 3, rate = 3),
    exponential = function() stats::rexp(5)
  )
  clocks <- list(
    lognormal   = clock_relaxed(0.4, sd = sd),
    gamma       = clock_relaxed(0.4, dist = "gamma", shape = 3),
    exponential = clock_relaxed(0.4, dist = "exponential")
  )
  for (dist in names(clocks)) {
    expect_equal(withr::with_seed(1, distances(clocks[[dist]])),
                 withr::with_seed(1, 0.4 * lengths * expected[[dist]]()), label = dist)
  }

  expect_equal(distances(clock_relaxed(0.4, sd = 0)), distances(clock_strict(0.4)))
  withr::local_seed(2)
  expect_false(isTRUE(all.equal(distances(clocks$gamma), distances(clocks$gamma))))  # new draw per call
})

test_that("lognormal multipliers have mean 1 and log-scale sd `sd`", {
  sd <- 0.5
  m <- withr::with_seed(3, draw_multipliers(clock_relaxed(sd = sd)))
  n <- length(m)
  expect_within(mean(m), 1, 4 * sqrt((exp(sd^2) - 1) / n))
  expect_within(mean(log(m)), -sd^2 / 2, 4 * sd / sqrt(n))
  expect_within(stats::sd(log(m)), sd, 4 * sd / sqrt(2 * n))
})

test_that("gamma multipliers have mean 1 and variance 1 / shape", {
  shape <- 4
  m <- withr::with_seed(4, draw_multipliers(clock_relaxed(dist = "gamma", shape = shape)))
  n <- length(m)
  m4 <- 3 * (shape + 2) / shape^3  # 4th central moment of Gamma(shape, rate = shape)
  expect_within(mean(m), 1, 4 * sqrt(1 / shape / n))
  expect_within(stats::var(m), 1 / shape, 4 * sqrt((m4 - 1 / shape^2) / n))
})

test_that("exponential multipliers have mean 1 and P(m < 1) = 1 - exp(-1)", {
  m <- withr::with_seed(5, draw_multipliers(clock_relaxed(dist = "exponential")))
  n <- length(m)
  p <- 1 - exp(-1)
  expect_within(mean(m), 1, 4 / sqrt(n))
  expect_within(mean(m < 1), p, tolerance_4se(p, n))
})

# ---- clock_tree, distance_tree ----

test_that("the root distance is the root edge, or 0 without one", {
  no_edge <- tree
  no_edge@phylo$root.edge <- NULL
  expect_equal(distances(clock_strict(2), no_edge), 2 * c(1, 1, 1.5, 0, 0.5))
})

test_that("clock_tree keeps branch lengths in time and the other node data", {
  ct <- clock_tree(tree, clock_strict(0.4))
  expect_equal(ct@phylo, tree@phylo)
  expect_equal(ct@data$status, tree@data$status)
  expect_false("distance" %in% names(tree@data))
  expect_equal(distances(clock_strict(), ct), lengths)  # an existing column is replaced
})

test_that("clock_tree works on a tree without node data", {
  expect_equal(distances(clock_strict(2), treeio::as.treedata(tree@phylo)), 2 * lengths)
})

test_that("clock_tree checks its inputs", {
  expect_error(clock_tree(tree@phylo, clock_strict()))
  expect_error(clock_tree(tree, 1))
})

test_that("distance_tree uses the distances as branch lengths", {
  phy <- distance_tree(clock_tree(tree, clock_strict(0.4)))@phylo
  expect_equal(phy$edge.length, 0.4 * lengths[phy$edge[, 2]])
  expect_equal(phy$root.edge, 0.4 * 0.5)

  no_edge <- tree
  no_edge@phylo$root.edge <- NULL
  expect_null(distance_tree(clock_tree(no_edge))@phylo$root.edge)
})

test_that("a complete distance column is required", {
  expect_error(distance_tree(tree), "distance")
  ct <- clock_tree(tree)
  ct@data$distance[2] <- NA
  expect_error(distance_tree(ct), "distance")
})

# ---- spike ----

# number of spikes per branch of the small tree: both internal nodes (root 4, node 5) divide, so
# every branch except the root edge gets one
daughter <- c(1, 1, 1, 0, 1)

with_data <- function(...) {
  tr <- tree
  cols <- list(...)
  for (col in names(cols)) tr@data[[col]] <- cols[[col]]
  tr
}

test_that("spike validates its inputs", {
  s <- spike(0.1)
  expect_s3_class(s, "clock_spike")
  expect_equal(s[c("mean", "shape", "by_type")], list(mean = 0.1, shape = Inf, by_type = NULL))

  for (bad in list(-1, NA_real_, Inf, c(1, 2), "1")) expect_error(spike(bad), "mean")
  for (bad in list(0, -1, NA_real_, c(1, 2), "1")) expect_error(spike(0.1, shape = bad), "shape")
  for (bad in list(-1, numeric(0), "1", NA)) expect_error(spike(0.1, by_type = bad), "by_type")

  expect_error(clock_strict(spike = 0.1), "spike")
  expect_identical(clock_relaxed(dist = "exponential", spike = s)$spike, s)
})

test_that("a fixed spike adds `mean` per division on the branches below it, plus hidden divisions", {
  spk <- clock_strict(0, spike = spike(0.2))
  expect_equal(distances(clock_strict(0.4, spike = spike(0.2))), 0.4 * lengths + 0.2 * daughter)

  # without `status` every internal node divides; with it only status 2
  expect_equal(distances(spk, treeio::as.treedata(tree@phylo)), 0.2 * daughter)
  expect_equal(distances(spk, with_data(status = c(1, 1, 1, 2, 1))), 0.2 * c(0, 0, 1, 0, 1))

  # hidden divisions add to the count, also on the root edge
  hidden <- c(0, 2, 1, 3, 0)
  expect_equal(distances(spk, with_data(hidden_divisions = hidden)), 0.2 * (daughter + hidden))
})

test_that("gamma spikes are one Gamma(N * shape, shape) per branch, drawn after and not scaled by the multipliers", {
  hidden <- c(0, 2, 1, 3, 0)
  clk <- clock_relaxed(0.4, dist = "gamma", shape = 3, spike = spike(0.1, shape = 2))
  expected <- function() {
    m <- stats::rgamma(5, shape = 3, rate = 3)
    0.4 * lengths * m + 0.1 * stats::rgamma(5, shape = 2 * (daughter + hidden), rate = 2)
  }
  expect_equal(withr::with_seed(1, distances(clk, with_data(hidden_divisions = hidden))),
               withr::with_seed(1, expected()))
})

test_that("gamma spikes have mean `mean` and variance mean^2 / shape", {
  mu <- 0.1
  shape <- 2
  clk <- clock_strict(0, spike = spike(mu, shape = shape))
  s <- withr::with_seed(6, as.vector(replicate(4000, distances(clk)[daughter == 1])))
  n <- length(s)
  v <- mu^2 / shape
  m4 <- 3 * (shape + 2) * mu^4 / shape^2  # 4th central moment of mu * Gamma(shape, rate = shape)
  expect_within(mean(s), mu, 4 * sqrt(v / n))
  expect_within(stats::var(s), v, 4 * sqrt((m4 - v^2) / n))
})

test_that("by_type multiplies the spike by the type of the dividing (parent) node", {
  # root 4 type 0 is the parent of 3, 5 and the root edge; node 5 type 1 the parent of 1, 2
  types <- c(0L, 1L, 1L, 0L, 1L)
  clk <- clock_strict(0, spike = spike(0.1, by_type = c(1, 3)))
  expect_equal(distances(clk, with_data(type = types)), 0.1 * c(3, 3, 1, 0, 1))
  expect_equal(distances(clk, with_data(type = types, hidden_divisions = c(0, 2, 1, 3, 0))),
               0.1 * c(3, 3 * 3, 1 * 2, 1 * 3, 1))

  # types that never divide need no entry (NA or beyond by_type)
  expect_equal(distances(clock_strict(0, spike = spike(0.1, by_type = c(2, NA))),
                         with_data(type = c(1L, 1L, 5L, 0L, 0L))), 0.2 * daughter)
})

test_that("clock_tree refuses missing or negative distances", {
  clk <- clock_strict(0, spike = spike(0.1, by_type = c(1, 2)))
  msg <- "missing or negative distances"
  expect_error(clock_tree(treeio::as.treedata(tree@phylo), clk), msg)            # no type column
  expect_error(clock_tree(with_data(type = c(0L, 0L, 0L, 0L, NA)), clk), msg)    # dividing node without type
  expect_error(clock_tree(with_data(type = c(0L, 0L, 0L, 0L, 2L)), clk), msg)    # dividing type beyond by_type
  expect_error(clock_tree(with_data(type = c(0L, 0L, 0L, 0L, 1L)),
                          clock_strict(0, spike = spike(0.1, by_type = c(1, NA)))), msg)
  expect_error(clock_tree(with_data(hidden_divisions = c(0, -2, 0, 0, 0)), clock_strict(0, spike = spike(0.1))), msg)
})
