# Sequentially-edited barcodes

withr::with_seed(1, {
  tree <- sim_adb_origin_samp(5, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
})
parent <- tibble::as_tibble(tree)$parent
n_tips <- length(tree@phylo$tip.label)
split_sites <- function(barcodes) strsplit(barcodes, "_", fixed = TRUE)

test_that("sim_barcode_seq returns one barcode string of n_sites per node and barcode", {
  withr::local_seed(1)
  out <- sim_barcode_seq(tree, n_sites = 5, edit_rate = c(0.3, 1), silencing_rate = 0.05,
                         n_barcodes = 2, chars = c("A", "B"))

  expect_named(out, c("node", "barcode_1", "barcode_2"))
  expect_equal(out$node, tibble::as_tibble(tree)$node)
  sites <- split_sites(c(out$barcode_1, out$barcode_2))
  expect_true(all(lengths(sites) == 5))
  expect_true(all(unlist(sites) %in% c("0", "A", "B", "-")))
})

test_that("edits fill sites left to right and are inherited, and silencing is inherited", {
  withr::local_seed(2)
  out <- sim_barcode_seq(tree, n_sites = 5, edit_rate = 1, silencing_rate = 0.1, n_barcodes = 1, chars = c("A", "B"))
  sites <- split_sites(out$barcode_1)
  parent_sites <- sites[match(parent, out$node)]

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
  args <- list(tree, n_sites = 5, edit_rate = 1, n_barcodes = 1, chars = c("A", "B"))
  with_dropout <- withr::with_seed(3, do.call(sim_barcode_seq, c(args, dropout_prob = 1)))
  without_dropout <- withr::with_seed(3, do.call(sim_barcode_seq, c(args, dropout_prob = 0)))
  is_tip <- with_dropout$node <= n_tips

  expect_true(all(with_dropout$barcode_1[is_tip] == "-_-_-_-_-"))
  expect_equal(with_dropout$barcode_1[!is_tip], without_dropout$barcode_1[!is_tip])
})

test_that("inserted characters follow chars and char_probs", {
  withr::local_seed(4)
  edited_chars <- function(out) setdiff(unlist(split_sites(out$barcode_1)), c("0", "-"))

  # a single number is used as a character, not as 1:n
  expect_equal(edited_chars(sim_barcode_seq(tree, n_sites = 3, edit_rate = 5, n_barcodes = 1, chars = 7)), "7")
  expect_equal(edited_chars(sim_barcode_seq(tree, n_sites = 3, edit_rate = 5, n_barcodes = 1, chars = c("A", "B"),
                                            char_probs = c(1, 0))), "A")
  expect_error(sim_barcode_seq(tree, n_sites = 3, edit_rate = 5, n_barcodes = 1, chars = c("A", "B"),
                               char_probs = 1), "`char_probs` must be")
  expect_error(sim_barcode_seq(tree, n_sites = 3, edit_rate = 5, n_barcodes = 1, chars = c("A", "B"),
                               char_probs = c(0.5, 0.6)), "`char_probs` must be")
  expect_error(sim_barcode_seq(tree, n_sites = 3, edit_rate = 5, dropout_prob = 1.2, n_barcodes = 1, chars = "A"),
               "`dropout_prob` must contain probabilities")
})

# Distributions: barcodes are independent given the tree (tips within a barcode are
# not, due to shared ancestry), so statistics are computed per barcode and tested
# across barcodes.
n_sites <- 5
edit_rate <- 0.4
silencing_rate <- 0.1
dropout_prob <- 0.2
chars <- c("A", "B", "C")
char_probs <- c(0.6, 0.3, 0.1)
many_barcodes <- withr::with_seed(5, {
  sim_barcode_seq(tree, n_sites, edit_rate, silencing_rate, dropout_prob, n_barcodes = 500, chars, char_probs)
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
  missing <- vapply(tip_sites, function(s) mean(rowSums(s == "-")), numeric(1))
  unedited <- vapply(tip_sites, function(s) mean(rowSums(s == "0")), numeric(1))

  expect_gt(stats::t.test(missing, mu = expected_missing)$p.value, 0.001)
  expect_gt(stats::t.test(unedited, mu = expected_unedited)$p.value, 0.001)
})

test_that("inserted characters follow char_probs", {
  # an edit is inherited by all descendants, so count each edit once: at the node
  # where it first appears (site edited, but unedited in the parent; root: unedited
  # at the origin). These are independent draws from char_probs.
  new_edit_chars <- unlist(lapply(barcode_sites, function(s) {
    parent_sites <- s[match(parent, many_barcodes$node), , drop = FALSE]
    parent_sites[parent == many_barcodes$node, ] <- "0"
    s[!(s %in% c("0", "-")) & parent_sites == "0"]
  }))
  observed <- table(factor(new_edit_chars, levels = chars))

  expect_gt(sum(observed), 1000)
  expect_gt(stats::chisq.test(observed, p = char_probs)$p.value, 0.001)
})

test_that("sim_barcode_seq requires an ultrametric (sampled) tree", {
  complete_tree <- withr::with_seed(1, sim_adb_ntaxa_complete_fast(20, scale = 1, shape = 1, death_prob = 0.2))
  expect_error(sim_barcode_seq(complete_tree, n_sites = 3, edit_rate = 1, n_barcodes = 1, chars = "A"),
               "must be ultrametric")
})
