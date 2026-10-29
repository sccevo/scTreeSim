# Sequentially-edited barcodes (SciPhy model)

withr::with_seed(1, {
  tree <- sim_adb_origin_samp(5, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
})
parent <- tibble::as_tibble(tree)$parent
n_tips <- length(tree@phylo$tip.label)
split_sites <- function(barcodes) strsplit(barcodes, "_", fixed = TRUE)

test_that("sim_barcode_seq returns one barcode string of n_sites per node and barcode", {
  withr::local_seed(1)
  out <- sim_barcode_seq(tree, n_barcodes = 2, n_sites = 5, edit_rate = c(0.3, 1),
                         edit_probs = c(0.5, 0.5), silencing_rate = 0.05)

  expect_named(out, c("node", "barcode_1", "barcode_2"))
  expect_equal(out$node, tibble::as_tibble(tree)$node)
  sites <- split_sites(c(out$barcode_1, out$barcode_2))
  expect_true(all(lengths(sites) == 5))
  expect_true(all(unlist(sites) %in% c("0", "1", "2", "3")))
})

test_that("edits fill sites left to right and are inherited, and silencing is inherited", {
  withr::local_seed(2)
  out <- sim_barcode_seq(tree, n_barcodes = 1, n_sites = 5, edit_rate = 1, edit_probs = c(0.5, 0.5), silencing_rate = 0.1)
  sites <- split_sites(out$barcode_1)
  parent_sites <- sites[match(parent, out$node)]

  # no unedited site before an edited one
  expect_true(all(vapply(sites, function(s) !is.unsorted(s == "0"), logical(1))))
  for (i in seq_along(sites)) {
    if (all(parent_sites[[i]] == "3")) {
      expect_true(all(sites[[i]] == "3"))
    } else if (!all(sites[[i]] == "3")) {
      edited <- parent_sites[[i]] != "0"
      expect_equal(sites[[i]][edited], parent_sites[[i]][edited])
    }
  }
})

test_that("write_as_string = FALSE returns one node x site data frame per barcode", {
  args <- list(tree, n_barcodes = 2, n_sites = 4, edit_rate = 1, edit_probs = c(0.5, 0.5),
               silencing_rate = 0.1, dropout_prob = 0.3)
  as_string <- withr::with_seed(7, do.call(sim_barcode_seq, args))
  as_list <- withr::with_seed(7, do.call(sim_barcode_seq, c(args, write_as_string = FALSE)))

  expect_named(as_list, c("barcode_1", "barcode_2"))
  for (bc in names(as_list)) {
    df <- as_list[[bc]]
    expect_s3_class(df, "data.frame")
    expect_named(df, c("node", paste0("site_", 1:4)))
    expect_equal(df$node, as_string$node)
    expect_equal(apply(df[-1], 1, paste, collapse = "_"), as_string[[bc]])
  }
  # one site and one node still give a data frame
  one <- sim_barcode_seq(tree, n_barcodes = 1, n_sites = 1, edit_rate = 1, edit_probs = 1, write_as_string = FALSE)
  expect_equal(dim(one$barcode_1), c(nrow(as_string), 2L))
})

test_that("missing_state defaults to E + 1 and its type sets the column type", {
  args <- list(tree, n_barcodes = 1, n_sites = 3, edit_rate = 1, edit_probs = c(0.5, 0.5),
               dropout_prob = 1, write_as_string = FALSE)
  out_default <- do.call(sim_barcode_seq, args)$barcode_1
  out_numeric <- do.call(sim_barcode_seq, c(args, missing_state = -1))$barcode_1
  out_string <- do.call(sim_barcode_seq, c(args, missing_state = "X"))$barcode_1
  is_tip <- out_default$node <= n_tips

  expect_true(all(vapply(out_default[-1], is.integer, logical(1))))
  expect_true(all(unlist(out_default[is_tip, -1]) == 3L))
  expect_true(all(vapply(out_numeric[-1], is.integer, logical(1))))
  expect_true(all(unlist(out_numeric[is_tip, -1]) == -1L))
  expect_true(all(vapply(out_string[-1], is.character, logical(1))))
  expect_true(all(unlist(out_string[is_tip, -1]) == "X"))

  as_string <- do.call(sim_barcode_seq, c(args[names(args) != "write_as_string"], missing_state = "X"))
  expect_true(all(as_string$barcode_1[is_tip] == "X_X_X"))
})

test_that("parents are visited before children on zero-length branches with unordered node numbers", {
  # internal node 6 is a child of node 7, but has a lower number; both branches have length 0
  phy <- structure(list(
    edge = matrix(c(5L, 7L, 7L, 6L, 6L, 1L, 6L, 2L, 7L, 3L, 5L, 4L), ncol = 2, byrow = TRUE),
    edge.length = c(0, 0, 1, 1, 1, 1), tip.label = c("A", "B", "C", "D"), Nnode = 3L
  ), class = "phylo")
  out <- withr::with_seed(8, sim_barcode_seq(treeio::as.treedata(phy), n_barcodes = 1, n_sites = 3,
                                             edit_rate = 2, edit_probs = c(0.5, 0.5)))
  expect_equal(out$node, 1:7)
  # nodes joined by zero-length branches share the same barcode
  expect_equal(out$barcode_1[c(5, 6, 7)], rep(out$barcode_1[5], 3))
})

test_that("dropout masks whole barcodes at the tips only", {
  args <- list(tree, n_barcodes = 1, n_sites = 5, edit_rate = 1, edit_probs = c(0.5, 0.5))
  with_dropout <- withr::with_seed(3, do.call(sim_barcode_seq, c(args, dropout_prob = 1)))
  without_dropout <- withr::with_seed(3, do.call(sim_barcode_seq, c(args, dropout_prob = 0)))
  is_tip <- with_dropout$node <= n_tips

  expect_true(all(with_dropout$barcode_1[is_tip] == "3_3_3_3_3"))
  expect_equal(with_dropout$barcode_1[!is_tip], without_dropout$barcode_1[!is_tip])
})

test_that("edit outcomes are 1..E and follow edit_probs", {
  withr::local_seed(4)
  edited <- function(out) sort(unique(setdiff(unlist(split_sites(out$barcode_1)), c("0", "-"))))

  expect_equal(edited(sim_barcode_seq(tree, n_barcodes = 1, n_sites = 3, edit_rate = 5, edit_probs = 1,
                                  missing_state = "-")), "1")
  expect_equal(edited(sim_barcode_seq(tree, n_barcodes = 1, n_sites = 3, edit_rate = 5, edit_probs = c(0, 1, 0),
                                  missing_state = "-")), "2")
})


# Distributions: barcodes are independent given the tree (tips within a barcode are
# not, due to shared ancestry), so statistics are computed per barcode and tested
# across barcodes.
n_sites <- 5
edit_rate <- 0.4
silencing_rate <- 0.1
dropout_prob <- 0.2
edit_probs <- c(0.6, 0.3, 0.1)
many_barcodes <- withr::with_seed(5, {
  sim_barcode_seq(tree, n_barcodes = 500, n_sites, edit_rate, edit_probs, silencing_rate, dropout_prob)
})
barcode_sites <- lapply(many_barcodes[-1], function(b) do.call(rbind, split_sites(b))) # nodes x sites per barcode

test_that("numbers of missing and unedited sites at the tips match their expectations", {
  # every tip is observed after time obs_time since the origin (ultrametric tree); a barcode
  # is missing if silenced (rate silencing_rate) or dropped out; otherwise silencing
  # and editing are independent, and min(Poisson(edit_rate * obs_time), n_sites) sites are edited
  obs_time <- tree@phylo$origin
  p_missing <- 1 - exp(-silencing_rate * obs_time) * (1 - dropout_prob)
  expected_missing <- n_sites * p_missing
  expected_unedited <- (1 - p_missing) * sum((n_sites - 0:(n_sites - 1)) * stats::dpois(0:(n_sites - 1), edit_rate * obs_time))

  tip_sites <- lapply(barcode_sites, function(s) s[seq_len(n_tips), , drop = FALSE])
  missing <- vapply(tip_sites, function(s) mean(rowSums(s == "4")), numeric(1))
  unedited <- vapply(tip_sites, function(s) mean(rowSums(s == "0")), numeric(1))

  expect_gt(stats::t.test(missing, mu = expected_missing)$p.value, 0.001)
  expect_gt(stats::t.test(unedited, mu = expected_unedited)$p.value, 0.001)
})

test_that("edit outcomes follow edit_probs", {
  # an edit is inherited by all descendants, so count each edit once: at the node
  # where it first appears (site edited, but unedited in the parent; root: unedited
  # at the origin). These are independent draws from edit_probs.
  new_edit_chars <- unlist(lapply(barcode_sites, function(s) {
    parent_sites <- s[match(parent, many_barcodes$node), , drop = FALSE]
    parent_sites[parent == many_barcodes$node, ] <- "0"
    s[!(s %in% c("0", "4")) & parent_sites == "0"]
  }))
  observed <- table(factor(new_edit_chars, levels = seq_along(edit_probs)))

  expect_gt(sum(observed), 1000)
  expect_gt(stats::chisq.test(observed, p = edit_probs)$p.value, 0.001)
})

test_that("sim_barcode_seq requires an ultrametric (sampled) tree", {
  complete_tree <- withr::with_seed(1, sim_adb_ntaxa_complete_fast(20, scale = 1, shape = 1, death_prob = 0.2))
  expect_error(sim_barcode_seq(complete_tree, n_barcodes = 1, n_sites = 3, edit_rate = 1, edit_probs = 1),
               "must be ultrametric")
})

test_that("a mismatching tree origin raises a warning", {
  tree_bad_origin <- tree
  tree_bad_origin@phylo$origin <- tree@phylo$origin + 1
  expect_warning(sim_barcode_seq(tree_bad_origin, n_barcodes = 1, n_sites = 3, edit_rate = 1, edit_probs = 1),
                 "origin")
  expect_no_warning(sim_barcode_seq(tree, n_barcodes = 1, n_sites = 3, edit_rate = 1, edit_probs = 1))
})
