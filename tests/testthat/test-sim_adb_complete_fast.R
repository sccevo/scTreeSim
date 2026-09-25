# Structure of complete trees from the C++ simulators

expect_complete_tree <- function(tree, origin) {
  ph <- tree@phylo
  n_tip <- length(ph$tip.label)
  alive <- tree@data$node[tree@data$status == 1]

  expect_s4_class(tree, "treedata")
  expect_equal(ph$Nnode, n_tip - 1)                               # binary tree
  expect_true(all(ph$edge.length >= 0))                           # positive branch lengths
  expect_equal(ph$tip.label, as.character(seq_len(n_tip)))        # tip labels = node numbers
  expect_equal(alive, seq_along(alive))                           # alive tips come first
  expect_true(all(tree@data$status[tree@data$node > n_tip] == 2)) # internal nodes divided
  # alive tips end at the present
  expect_equal(ape::node.depth.edgelength(ph)[alive] + ph$root.edge, rep(origin, length(alive)))
}

test_that("sim_adb_ntaxa_complete_fast returns a complete tree with ntaxa alive tips", {
  withr::local_seed(1)
  tree <- sim_adb_ntaxa_complete_fast(20, a = 1, b = 1, d = 0.2)

  expect_complete_tree(tree, origin = tree@phylo$origin)
  expect_equal(sum(tree@data$status == 1), 20)
  expect_true(any(tree@data$status == 0))
})

test_that("sim_adb_origin_complete_fast returns a complete tree ending at origin_time", {
  withr::local_seed(1)
  tree <- sim_adb_origin_complete_fast(5, a = 1, b = 1, d = 0.2)

  expect_complete_tree(tree, origin = 5)
  expect_equal(tree@phylo$origin, 5)
  expect_true(any(tree@data$status == 0))
})

test_that("large trees are complete trees (storage is resized in the C++ loops)", {
  withr::local_seed(1)
  # ntaxa: initial capacity 4 * ntaxa + 10 nodes, exceeded when many particles die
  tree <- simulate_trees(10, function() sim_adb_ntaxa_complete_fast(50, a = 1, b = 1, d = 0.35))[[1]]
  expect_gt(nrow(tree@data), 4 * 50 + 10)
  expect_complete_tree(tree, origin = tree@phylo$origin)

  # origin: initial capacity 1024 nodes
  tree <- sim_adb_origin_complete_fast(10, a = 1, b = 1, d = 0)
  expect_gt(nrow(tree@data), 1024)
  expect_complete_tree(tree, origin = 10)
})

test_that("failed simulations return NULL with a message", {
  withr::local_seed(1)
  expect_message(expect_null(sim_adb_ntaxa_complete_fast(20, a = 1, b = 1, d = 0.9)), "Too many particles died")
  expect_message(expect_null(sim_adb_origin_complete_fast(5, a = 1, b = 1, d = 0, min_tips = 1e6)), "too few tips")
  # the root outlives the origin interval: a single alive tip, never divides
  expect_message(expect_null(sim_adb_origin_complete_fast(1, a = 1000, b = 1, d = 0)), "too few tips")
})
