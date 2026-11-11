# Clock models and per-branch evolutionary distances

# Small tree from helper-barcodes.R: nodes a(1), b(2), c(3), root(4), internal node 5.
# Branch lengths ending in nodes 1..5: 1, 1, 1.5, root edge 0.5, 0.5.
tree <- small_barcode_tree()
branch_lengths <- c(1, 1, 1.5, 0.5, 0.5)

test_that("clock_strict stores its rate and validates it", {
  clk <- clock_strict(0.3)
  expect_s3_class(clk, "clock_model")
  expect_equal(clk$rate, 0.3)
  expect_null(clk$spike)
  expect_equal(clock_strict()$rate, 1)
  expect_equal(clock_strict(0)$rate, 0)  # spikes only, no gradual evolution

  expect_error(clock_strict(-1))
  expect_error(clock_strict(NA_real_))
  expect_error(clock_strict(Inf))
  expect_error(clock_strict(c(1, 2)))
  expect_error(clock_strict("1"))
})

# distances written by clock_tree, as a vector indexed by node id
clocked <- function(tr, clock) .node_distances(clock_tree(tr, clock))

test_that("strict clock distances are rate times branch length, indexed by node", {
  expect_equal(clocked(tree, clock_strict()), branch_lengths)
  expect_equal(clocked(tree, clock_strict(0.4)), 0.4 * branch_lengths)
  expect_equal(clocked(tree, clock_strict(0)), rep(0, 5))
})

test_that("the root distance is the root edge, or 0 without one", {
  no_edge <- tree
  no_edge@phylo$root.edge <- NULL
  expect_equal(clocked(no_edge, clock_strict(2)), 2 * c(1, 1, 1.5, 0, 0.5))
})

test_that("clock_tree keeps branch lengths in time and the other node data", {
  ct <- clock_tree(tree, clock_strict(0.4))
  expect_equal(ct@phylo, tree@phylo)
  expect_equal(ct@data$status, tree@data$status)
  expect_false("distance" %in% names(tree@data))  # the input is not modified

  # default strict clock at rate 1, and an existing column is replaced
  expect_equal(clocked(ct, clock_strict()), branch_lengths)
})

test_that("clock_tree works on a tree without node data", {
  expect_equal(clocked(treeio::as.treedata(tree@phylo), clock_strict(2)), 2 * branch_lengths)
})

test_that("clock_tree checks its inputs", {
  expect_error(clock_tree(tree@phylo, clock_strict()))
  expect_error(clock_tree(tree, 1))
})

test_that("distance_tree uses the distances as branch lengths", {
  dt <- distance_tree(clock_tree(tree, clock_strict(0.4)))
  phy <- dt@phylo
  expect_equal(phy$edge.length, 0.4 * branch_lengths[phy$edge[, 2]])
  expect_equal(phy$root.edge, 0.4 * 0.5)
  # tips at distance 0.4 * 2 from the origin
  expect_equal(ape::node.depth.edgelength(phy)[1:3] + phy$root.edge, rep(0.8, 3))

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
