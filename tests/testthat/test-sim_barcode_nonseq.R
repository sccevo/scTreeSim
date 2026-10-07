# Non-sequentially edited barcodes (TiDeTree's edit-and-silencing model)

withr::with_seed(1, {
  tree <- sim_adb_origin_samp(5, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
})
tree_df <- tibble::as_tibble(tree)
tip_nodes <- tree_df$node[tree_df$status == 1]
tree_age <- max(ape::node.depth.edgelength(tree@phylo))

out_tips <- function(tr) which(tibble::as_tibble(tr)$status == 1)

test_that("transition probability rows sum to 1 and silenced is absorbing", {
  for (in_window in c(TRUE, FALSE)) {
    P <- .edit_silencing_transition_probs(1.3, c(0.25, 0.75), 0.1, 0.7, in_window)
    expect_equal(rowSums(P), rep(1, 4))
    expect_equal(P[4, ], c(0, 0, 0, 1))
  }
})

test_that("sim_barcode_nonseq returns one row per node and n_targets site columns", {
  withr::local_seed(1)
  out <- sim_barcode_nonseq(tree, n_targets = 4, edit_rate = 1.3, edit_probs = c(0.25, 0.75), silencing_rate = 0.05,
                            edit_height = tree_age, edit_duration = tree_age)
  expect_named(out, c("node", paste0("site_", 1:4)))
  expect_equal(out$node, tree_df$node)
  expect_true(all(unlist(out[-1]) %in% 0:3))
})

test_that("no edits happen without edit rates or outside the editing window", {
  withr::local_seed(2)
  out <- sim_barcode_nonseq(tree, n_targets = 5, edit_rate = 0, edit_probs = c(0.5, 0.5), silencing_rate = 0,
                            edit_height = tree_age, edit_duration = tree_age)
  expect_true(all(unlist(out[-1]) == 0))

  out <- sim_barcode_nonseq(tree, n_targets = 5, edit_rate = 5, edit_probs = 1, silencing_rate = 0,
                            edit_height = tree_age + 10, edit_duration = 1)
  expect_true(all(unlist(out[-1]) == 0))
})

test_that("dropout_prob = 1 silences all tips only", {
  withr::local_seed(3)
  out <- sim_barcode_nonseq(tree, n_targets = 3, edit_rate = 1, edit_probs = 1, silencing_rate = 0,
                            edit_height = tree_age, edit_duration = tree_age, dropout_prob = 1)
  tips <- out$node %in% tip_nodes
  expect_true(all(unlist(out[tips, -1]) == 2))
  expect_false(any(unlist(out[!tips, -1]) == 2))
})

test_that("states are inherited: edits and silencing are never undone", {
  withr::local_seed(4)
  out <- sim_barcode_nonseq(tree, n_targets = 5, edit_rate = 3, edit_probs = c(1/3, 2/3), silencing_rate = 0.2,
                            edit_height = tree_age, edit_duration = tree_age)
  m <- as.matrix(out[-1])
  pm <- m[match(tree_df$parent, out$node), , drop = FALSE]
  expect_true(all(m[pm != 0] == pm[pm != 0] | m[pm != 0] == 3))
})

test_that("edit frequency on a single branch matches the closed form", {
  withr::local_seed(5)
  phy <- ape::read.tree(text = "(a:1,b:1);")
  tr <- treeio::as.treedata(phy)
  tr@data <- tibble::tibble(node = 1:3, status = c(1, 1, 2))
  R <- 2; p <- c(0.25, 0.75); s <- 0.3
  out <- sim_barcode_nonseq(tr, n_targets = 4000, edit_rate = R, edit_probs = p, silencing_rate = s,
                            edit_height = 1, edit_duration = 1)
  a <- unlist(out[out$node == 1, -1])
  expected <- c(exp(-(R + s)), p * (exp(-s) - exp(-(R + s))), 1 - exp(-s))
  expect_equal(as.numeric(table(factor(a, 0:3))) / 4000, expected, tolerance = 0.03)
})

test_that("zero-length branches do not break the parent-before-child order", {
  withr::local_seed(6)
  phy <- ape::read.tree(text = "((a:1,b:1):0,c:1);")
  tr <- treeio::as.treedata(phy)
  tr@data <- tibble::tibble(node = 1:5, status = c(1, 1, 1, 2, 2))
  out <- sim_barcode_nonseq(tr, n_targets = 3, edit_rate = 1, edit_probs = 1, silencing_rate = 0.1,
                            edit_height = 1, edit_duration = 1)
  expect_equal(nrow(out), 5)
})

test_that("edit_rate can be site-specific", {
  withr::local_seed(7)
  out <- sim_barcode_nonseq(tree, n_targets = 2, edit_rate = c(0, 50), edit_probs = 1, silencing_rate = 0)
  expect_true(all(out$site_1 == 0))
  expect_true(all(out$site_2[out$node %in% tip_nodes] == 1))
})

test_that("editing window defaults to the origin, or to the tree height", {
  # large rate: everything is edited iff the window covers the tree
  run <- function(tr) sim_barcode_nonseq(tr, n_targets = 3, edit_rate = 100, edit_probs = 1, silencing_rate = 0)
  expect_true(all(unlist(run(tree)[out_tips(tree), -1]) == 1))
  tree_no_origin <- tree
  tree_no_origin@phylo$origin <- NULL
  expect_true(all(unlist(run(tree_no_origin)[out_tips(tree_no_origin), -1]) == 1))
})

test_that("the default editing window starts at the origin, including the root edge", {
  withr::local_seed(8)
  tree_root_edge <- tree
  tree_root_edge@phylo$origin <- NULL
  expect_false(is.null(tree_root_edge@phylo$root.edge))
  # a rate this large edits every target on the root edge alone
  out <- sim_barcode_nonseq(tree_root_edge, n_targets = 3, edit_rate = 1e3, edit_probs = 1)
  expect_true(all(unlist(out[out$node == tree_df$node[tree_df$parent == tree_df$node], -1]) == 1))
})

test_that("missing_state defaults to E + 1 and its type sets the column type", {
  args <- list(tree, n_targets = 3, edit_rate = 1, edit_probs = c(0.5, 0.5), dropout_prob = 1)
  tips <- tree_df$node %in% tip_nodes
  out_default <- withr::with_seed(10, do.call(sim_barcode_nonseq, args))
  out_numeric <- withr::with_seed(10, do.call(sim_barcode_nonseq, c(args, missing_state = -1)))
  out_string <- withr::with_seed(10, do.call(sim_barcode_nonseq, c(args, missing_state = "-")))

  expect_true(all(vapply(out_default[-1], is.integer, logical(1))))
  expect_true(all(unlist(out_default[tips, -1]) == 3L))
  expect_true(all(unlist(out_numeric[tips, -1]) == -1))
  expect_true(all(vapply(out_string[-1], is.character, logical(1))))
  expect_true(all(unlist(out_string[tips, -1]) == "-"))
  # relabelling only: the same draws otherwise
  expect_equal(as.matrix(out_string[-1]) == "-", as.matrix(out_default[-1]) == 3L, ignore_attr = TRUE)
})

test_that("a mismatching tree origin raises a warning", {
  tree_bad_origin <- tree
  tree_bad_origin@phylo$origin <- tree@phylo$origin + 1
  expect_warning(sim_barcode_nonseq(tree_bad_origin, n_targets = 3, edit_rate = 1, edit_probs = 1), "origin")
  expect_no_warning(sim_barcode_nonseq(tree, n_targets = 3, edit_rate = 1, edit_probs = 1))
})
