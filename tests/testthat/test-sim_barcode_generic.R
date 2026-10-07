# Generic barcodes (exponential waiting times)

withr::with_seed(1, {
  tree <- sim_adb_origin_samp(5, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
})
tree_df <- tibble::as_tibble(tree)
split_sites <- function(barcodes) strsplit(barcodes, "_", fixed = TRUE)

test_that("sim_barcode_generic returns one barcode string of n_sites per node and barcode", {
  withr::local_seed(1)
  out <- sim_barcode_generic(tree, n_barcodes = 2, n_sites = 5, edit_rate = c(0.3, 1),
                             edit_probs = c(0.5, 0.5), silencing_rate = 0.05)

  expect_named(out, c("node", "barcode_1", "barcode_2"))
  expect_equal(out$node, tree_df$node)
  sites <- split_sites(c(out$barcode_1, out$barcode_2))
  expect_true(all(lengths(sites) == 5))
  expect_true(all(unlist(sites) %in% c("0", "1", "2", "-")))
})

test_that("write_as_string = FALSE returns one node x site data frame per barcode", {
  set <- function(...) {
    withr::local_seed(8)
    sim_barcode_generic(tree, n_barcodes = 2, n_sites = 4, edit_rate = c(0.3, 1), edit_probs = c(A = 0.5, B = 0.5),
                        silencing_rate = 0.05, dropout_prob = 0.2, ...)
  }
  str_out <- set()
  df_out <- set(write_as_string = FALSE)

  expect_named(df_out, c("barcode_1", "barcode_2"))
  for (bc in names(df_out)) {
    expect_named(df_out[[bc]], c("node", paste0("site_", 1:4)))
    expect_equal(df_out[[bc]]$node, tree_df$node)
    # same draws, so the string form is the pasted data frame
    expect_equal(apply(df_out[[bc]][-1], 1, paste, collapse = "_"), str_out[[bc]], ignore_attr = TRUE)
  }

  # single site: still a data frame, not a vector
  one <- sim_barcode_generic(tree, 1, 1, edit_rate = 1, edit_probs = 1, write_as_string = FALSE)
  expect_named(one$barcode_1, c("node", "site_1"))
  expect_equal(nrow(one$barcode_1), nrow(tree_df))
})

test_that("names of edit_probs are the outcome labels", {
  withr::local_seed(2)
  out <- sim_barcode_generic(tree, 1, 5, edit_rate = 1, edit_probs = c(A = 0.7, B = 0.3), missing_state = "x",
                             silencing_rate = 0.1, dropout_prob = 0.3)
  states <- unique(unlist(split_sites(out$barcode_1)))
  expect_true(all(states %in% c("0", "A", "B", "x")))
  expect_true(all(c("A", "B") %in% states))

  # unnamed: 1..E
  out <- sim_barcode_generic(tree, 1, 5, edit_rate = 1, edit_probs = c(0.2, 0.3, 0.5))
  expect_true(all(unlist(split_sites(out$barcode_1)) %in% as.character(0:3)))
})

test_that("edits fill sites left to right and are inherited, and silencing is inherited", {
  withr::local_seed(3)
  out <- sim_barcode_generic(tree, 1, 5, edit_rate = 1, edit_probs = c(A = 0.5, B = 0.5), silencing_rate = 0.1)
  sites <- split_sites(out$barcode_1)
  parent_sites <- sites[match(tree_df$parent, out$node)]

  # no unedited site before an edited one
  expect_true(all(vapply(sites, function(s) !is.unsorted(s == "0"), logical(1))))
  for (i in seq_along(sites)) {
    if (all(parent_sites[[i]] == "-")) {
      expect_true(all(sites[[i]] == "-"))
    } else if (!all(sites[[i]] == "-")) {
      edited <- parent_sites[[i]] != "0"
      expect_equal(sites[[i]][edited], parent_sites[[i]][edited])
    }
  }
})

test_that("dropout masks whole barcodes at the tips only", {
  withr::local_seed(4)
  out <- sim_barcode_generic(tree, 1, 3, edit_rate = 1, edit_probs = 1, dropout_prob = 1)
  tips <- out$node %in% seq_along(tree@phylo$tip.label)
  expect_true(all(out$barcode_1[tips] == "-_-_-"))
  expect_false(any(grepl("-", out$barcode_1[!tips])))
})

test_that("zero-length branches do not break the parent-before-child order", {
  withr::local_seed(5)
  phy <- ape::read.tree(text = "((a:1,b:1):0,c:1);")
  tr <- treeio::as.treedata(phy)
  tr@data <- tibble::tibble(node = 1:5, status = c(1, 1, 1, 2, 2))
  out <- sim_barcode_generic(tr, 1, 3, edit_rate = 1, edit_probs = 1, silencing_rate = 0.1)
  expect_equal(nrow(out), 5)
  expect_equal(out$barcode_1[4], out$barcode_1[5])  # zero-length branch: no change
})

test_that("the root edge is evolved too", {
  withr::local_seed(6)
  phy <- ape::read.tree(text = "(a:1,b:1);")
  phy$root.edge <- 100
  tr <- treeio::as.treedata(phy)
  tr@data <- tibble::tibble(node = 1:3, status = c(1, 1, 2))
  out <- sim_barcode_generic(tr, 1, 2, edit_rate = 1, edit_probs = 1)
  expect_true(all(out$barcode_1 == "1_1"))
})

test_that("invalid arguments are rejected", {
  expect_error(sim_barcode_generic(tree, 1, 3, edit_rate = -1, edit_probs = 1))
  expect_error(sim_barcode_generic(tree, 2, 3, edit_rate = c(1, 1, 1), edit_probs = 1))
  expect_error(sim_barcode_generic(tree, 1, 3, edit_rate = 1, edit_probs = c(0.5, 0.6)), "summing to 1")
  expect_error(sim_barcode_generic(tree, 1, 3, edit_rate = 1, edit_probs = c(A = 0.5, A = 0.5)), "unique")
  expect_error(sim_barcode_generic(tree, 1, 3, edit_rate = 1, edit_probs = c(`a_b` = 1)), "_")
  expect_error(sim_barcode_generic(tree, 1, 3, edit_rate = 1, edit_probs = c(A = 1), missing_state = "A"), "missing_state")
  expect_error(sim_barcode_generic(tree, 1, 3, edit_rate = 1, edit_probs = 1, dropout_prob = 2))
  expect_error(sim_barcode_generic(tree, 1, 3, edit_rate = 1, edit_probs = 1, write_as_string = NA))
})
