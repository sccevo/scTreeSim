# Pruning of dead and unsampled particles, and the sampled-tree simulators

# expected structure of a pruned tree: only alive tips, ending at the present
expect_pruned_tree <- function(tree, origin) {
  ph <- tree@phylo
  tips <- seq_along(ph$tip.label)
  expect_true(all(tree@data$status[match(tips, tree@data$node)] == 1))
  expect_equal(ph$origin, origin)
  expect_equal(ape::node.depth.edgelength(ph)[tips] + ph$root.edge, rep(origin, length(tips)))
}

withr::with_seed(1, {
  complete_tree <- sim_adb_ntaxa_complete_fast(20, scale = 1, shape = 1, death_prob = 0.2) # contains dead tips
})
origin <- complete_tree@phylo$origin

test_that("prune_tree removes dead particles and samples a fixed number of tips", {
  withr::local_seed(1)
  tree <- prune_tree(complete_tree, ntips = 8)

  expect_length(tree@phylo$tip.label, 8)
  expect_equal(tree@phylo$Nnode, 7)
  expect_pruned_tree(tree, origin)
})

test_that("prune_tree samples tips with probability sampling_prob", {
  withr::local_seed(1)
  expect_length(prune_tree(complete_tree, sampling_prob = 1)@phylo$tip.label, 20)
  expect_pruned_tree(prune_tree(complete_tree, sampling_prob = 0.5), origin)
})

test_that("prune_tree samples tips with a type-dependent sampling_prob", {
  withr::local_seed(3)
  Xi_as <- matrix(c(0, 0, 0.3, 0), 2)
  Xi_s <- matrix(c(0.7, 0, 0, 1), 2)
  multitype_tree <- sim_adb_ntaxa_complete_fast(30, scale = c(1, 1), shape = c(1, 1), death_prob = c(0.2, 0.2), asym_trans_prob = Xi_as, sym_trans_prob = Xi_s)
  alive_type0 <- sum(multitype_tree@data$status == 1 & multitype_tree@data$type == 0)
  tree <- prune_tree(multitype_tree, sampling_prob = c(1, 0))
  tip_types <- tree@data$type[match(seq_along(tree@phylo$tip.label), tree@data$node)]

  expect_true(all(tip_types == 0))
  expect_length(tip_types, alive_type0)
  expect_pruned_tree(tree, multitype_tree@phylo$origin)
})

test_that("prune_tree keeps nodes with a single descendant if collapse = FALSE", {
  withr::local_seed(1)
  tree <- prune_tree(complete_tree, sampling_prob = 1, collapse = FALSE)

  expect_gt(tree@phylo$Nnode, length(tree@phylo$tip.label) - 1)
  expect_pruned_tree(tree, origin)
})

test_that("prune_tree handles too few tips and invalid arguments", {
  withr::local_seed(1)
  expect_message(expect_null(prune_tree(complete_tree, ntips = 1, min_taxa = 2)), "Not enough tips")
  expect_error(prune_tree(complete_tree), "either the sampling probability or the desired number")
  expect_error(prune_tree(complete_tree, sampling_prob = 0.5, ntips = 5), "either the sampling probability or the desired number")
})

test_that("sim_adb_ntaxa_samp and sim_adb_origin_samp return pruned trees", {
  withr::local_seed(1)
  tree <- sim_adb_ntaxa_samp(10, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
  expect_length(tree@phylo$tip.label, 10)
  expect_pruned_tree(tree, tree@phylo$origin)

  tree <- sim_adb_origin_samp(5, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
  expect_pruned_tree(tree, 5)
})
