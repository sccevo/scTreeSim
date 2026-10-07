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


# sequential switch and equivalence with the other simulators ------------------------------

n_rep <- 4000
r <- 1.2; s <- 0.3; p <- c(0.25, 0.75)

test_that("single-site tip states match the closed form, for both settings of sequential", {
  skip_on_cran()
  tr <- small_barcode_tree()
  expected <- single_site_probs(r, p, s)
  for (sequential in c(TRUE, FALSE)) {
    withr::local_seed(10)
    out <- sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_rate = r, edit_probs = p,
                               silencing_rate = s, sequential = sequential)
    states <- replicate_states(out)
    for (tip in 1:3) {
      observed <- as.numeric(table(factor(states[tip, ], names(expected)))) / n_rep
      expect_within(observed, expected, tolerance_4se(expected, n_rep))
    }
  }
})

test_that("single-site generic matches sim_barcode_seq and sim_barcode_nonseq jointly at sibling tips", {
  skip_on_cran()
  tr <- small_barcode_tree()
  withr::local_seed(11)
  generic <- replicate_states(sim_barcode_generic(tr, n_rep, 1, r, p, silencing_rate = s))
  generic_indep <- replicate_states(sim_barcode_generic(tr, n_rep, 1, r, p, silencing_rate = s, sequential = FALSE))
  seq <- replicate_states(sim_barcode_seq(tr, n_rep, 1, r, p, silencing_rate = s, missing_state = "-"))
  nonseq <- replicate_states(sim_barcode_nonseq(tr, n_rep, r, p, silencing_rate = s, missing_state = "-"))

  ref <- joint_counts(nonseq, 1:2)
  expect_gt(homogeneity_p(joint_counts(generic, 1:2), ref), 0.001)
  expect_gt(homogeneity_p(joint_counts(generic_indep, 1:2), ref), 0.001)
  expect_gt(homogeneity_p(joint_counts(seq, 1:2), ref), 0.001)
})

test_that("independent sites in one barcode match independent barcodes of one site (and nonseq)", {
  skip_on_cran()
  tr <- small_barcode_tree()
  withr::local_seed(12)
  # one barcode with n_rep independent sites: states over sites follow the single-site distribution
  sites <- sim_barcode_generic(tr, 1, n_rep, r, p, silencing_rate = s, sequential = FALSE, write_as_string = FALSE)
  sites <- as.matrix(sites$barcode_1[-1])  # nodes x sites
  expected <- single_site_probs(r, p, s)
  observed <- as.numeric(table(factor(sites[1, ], names(expected)))) / n_rep
  expect_within(observed, expected, tolerance_4se(expected, n_rep))

  nonseq <- replicate_states(sim_barcode_nonseq(tr, n_rep, r, p, silencing_rate = s, missing_state = "-"))
  expect_gt(homogeneity_p(joint_counts(sites, 1:2), joint_counts(nonseq, 1:2)), 0.001)
})

test_that("sequential = FALSE edits sites in any order, sequential = TRUE strictly left to right", {
  tr <- small_barcode_tree()
  withr::local_seed(13)
  tips <- 1:3
  split_tips <- function(sequential) {
    out <- sim_barcode_generic(tr, 200, 6, edit_rate = 0.5, edit_probs = c(A = 0.5, B = 0.5),
                               silencing_rate = 0.1, sequential = sequential)
    unlist(lapply(out[tips, -1], function(x) strsplit(x, "_", fixed = TRUE)), recursive = FALSE)
  }
  unedited_before_edited <- function(s) {
    edited <- which(!s %in% c("0", "-"))
    length(edited) > 0 && any(s[seq_len(max(edited))] == "0")
  }
  expect_false(any(vapply(split_tips(TRUE), unedited_before_edited, logical(1))))
  expect_true(any(vapply(split_tips(FALSE), unedited_before_edited, logical(1))))
})

test_that("silencing acts on the whole barcode if sequential, on single sites otherwise", {
  tr <- small_barcode_tree()
  withr::local_seed(14)
  mixed <- function(sequential) {
    out <- sim_barcode_generic(tr, 200, 6, edit_rate = 1, edit_probs = 1, silencing_rate = 0.4,
                               sequential = sequential)
    sites <- unlist(lapply(out[1:3, -1], function(x) strsplit(x, "_", fixed = TRUE)), recursive = FALSE)
    vapply(sites, function(s) any(s == "-") && !all(s == "-"), logical(1))
  }
  expect_false(any(mixed(TRUE)))
  expect_true(any(mixed(FALSE)))
})

test_that("dropout acts on the whole barcode if sequential, on single sites otherwise", {
  tr <- small_barcode_tree()
  withr::local_seed(15)
  mixed <- function(sequential) {
    out <- sim_barcode_generic(tr, 200, 6, edit_rate = 0, edit_probs = 1, dropout_prob = 0.5,
                               sequential = sequential)
    tips <- unlist(lapply(out[1:3, -1], function(x) strsplit(x, "_", fixed = TRUE)), recursive = FALSE)
    inner <- unlist(lapply(out[4:5, -1], function(x) strsplit(x, "_", fixed = TRUE)), recursive = FALSE)
    list(mixed = vapply(tips, function(s) any(s == "-") && !all(s == "-"), logical(1)),
         inner_missing = any(unlist(inner) == "-"))
  }
  seq_res <- mixed(TRUE)
  ind_res <- mixed(FALSE)
  expect_false(any(seq_res$mixed))
  expect_true(any(ind_res$mixed))
  # dropout is applied at the tips only
  expect_false(seq_res$inner_missing)
  expect_false(ind_res$inner_missing)

  # rate of per-site dropout
  out <- sim_barcode_generic(tr, 1, 4000, edit_rate = 0, edit_probs = 1, dropout_prob = 0.3,
                             sequential = FALSE, write_as_string = FALSE)
  frac <- mean(out$barcode_1[1, -1] == "-")
  expect_within(frac, 0.3, tolerance_4se(0.3, 4000))
})


# time-varying rates ---------------------------------------------

# fraction of replicate single-site barcodes that are edited (not "0" and not "-") at the tips
edited_fraction <- function(out) mean(!replicate_states(out)[1:3, ] %in% c("0", "-"))

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

test_that("a constant profile via barcode_rate() or rate_fn() leaves the simulator unchanged", {
  tr <- small_barcode_tree()
  args <- list(tr, n_barcodes = 5, n_sites = 3, edit_probs = c(A = 0.5, B = 0.5), silencing_rate = 0.1)
  withr::local_seed(3)
  plain <- do.call(sim_barcode_generic, c(args, edit_rate = 0.8))
  withr::local_seed(3)
  const <- do.call(sim_barcode_generic, c(args, list(edit_rate = barcode_rate(0.8, rate_constant()))))
  expect_identical(const, plain)
  # a vector base rate with a profile is accepted per barcode
  out <- sim_barcode_generic(tr, n_barcodes = 2, n_sites = 2, edit_rate = barcode_rate(c(0, 5), rate_constant()),
                             edit_probs = 1)
  expect_true(all(out$barcode_1 == "0_0"))
  expect_error(sim_barcode_generic(tr, n_barcodes = 2, n_sites = 2, edit_probs = 1,
                                   edit_rate = barcode_rate(c(1, 2, 3))), "length")
})

test_that("a window outside the tree gives no edits, and a covering window matches a constant rate", {
  skip_on_cran()
  tr <- small_barcode_tree()
  n_rep <- 3000
  run <- function(rate) {
    sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_probs = 1, edit_rate = rate)
  }
  withr::local_seed(21)
  expect_equal(edited_fraction(run(barcode_rate(5, rate_window(10, 20)))), 0)
  expect_equal(edited_fraction(run(barcode_rate(5, rate_window(-5, 0)))), 0)
  covering <- edited_fraction(run(barcode_rate(0.7, rate_window(-1, 10))))
  expected <- 1 - exp(-0.7 * small_barcode_time)
  expect_within(covering, expected, tolerance_4se(expected, 3 * n_rep))
})

test_that("window gives the closed form of the time spent inside it, matching nonseq's edit window", {
  skip_on_cran()
  tr <- small_barcode_tree()
  n_rep <- 4000
  r <- 1; start <- 0.8; end <- 1.5   # time since origin; origin - height = start
  # nonseq: times before the present, tips at height 0, origin at height 2
  edit_height <- small_barcode_time - start
  edit_duration <- end - start
  expected <- 1 - exp(-r * (end - start))
  withr::local_seed(22)
  gen <- sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_probs = 1,
                             edit_rate = barcode_rate(r, rate_window(start, end)))
  non <- sim_barcode_nonseq(tr, n_targets = n_rep, edit_rate = r, edit_probs = 1,
                            edit_height = edit_height, edit_duration = edit_duration)
  tol <- tolerance_4se(expected, 3 * n_rep)
  expect_within(edited_fraction(gen), expected, tol)
  expect_within(mean(replicate_states(non)[1:3, ] != "0"), expected, tol)
  expect_gt(homogeneity_p(joint_counts(replicate_states(gen), 1:2),
                          joint_counts(replicate_states(non), 1:2)), 0.001)
})

test_that("exponentially decaying editing follows the integrated rate and edits less late than a constant rate", {
  skip_on_cran()
  tr <- small_barcode_tree()
  n_rep <- 4000
  r <- 1.5; half_life <- 0.5
  # P(edited by T) = 1 - exp(-integral of r 2^(-t/h)); tips at time 2
  cum <- r * half_life / log(2) * (1 - 2^(-small_barcode_time / half_life))
  expected <- 1 - exp(-cum)
  withr::local_seed(23)
  decaying <- sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_probs = 1,
                                  edit_rate = barcode_rate(r, rate_exp_decay(half_life)))
  constant <- sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_probs = 1, edit_rate = r)
  expect_within(edited_fraction(decaying), expected, tolerance_4se(expected, 3 * n_rep))
  expect_lt(edited_fraction(decaying), edited_fraction(constant))
})

test_that("a piecewise profile and rate_fn agree with the closed form; silencing is time-varying too", {
  skip_on_cran()
  tr <- small_barcode_tree()
  n_rep <- 4000
  # rate 0 until 1, then 2, so the integrated rate to the tips (t = 2) is 2
  expected <- 1 - exp(-2)
  pw <- barcode_rate(1, rate_piecewise(1, c(0, 2)))
  fn <- barcode_rate(1, rate_fn(function(t) if (t < 1) 0 else 2, max = 2))
  for (rate in list(pw, fn)) {
    withr::local_seed(24)
    out <- sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_probs = 1, edit_rate = rate)
    expect_within(edited_fraction(out), expected, tolerance_4se(expected, 3 * n_rep))
  }
  # silencing only in the window [0, 1): P(silenced) = 1 - exp(-s)
  s <- 1.2
  withr::local_seed(25)
  out <- sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_probs = 1, edit_rate = 0,
                             silencing_rate = barcode_rate(s, rate_window(0, 1)))
  silenced <- mean(replicate_states(out)[1:3, ] == "-")
  expect_within(silenced, 1 - exp(-s), tolerance_4se(1 - exp(-s), 3 * n_rep))
})

test_that("a rate_fn exceeding its stated maximum is an error", {
  tr <- small_barcode_tree()
  expect_error(
    sim_barcode_generic(tr, n_barcodes = 20, n_sites = 1, edit_probs = 1,
                        edit_rate = barcode_rate(5, rate_fn(function(t) 3, max = 1))),
    "exceeds"
  )
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

test_that("decay inside a window follows the integrated rate", {
  skip_on_cran()
  tr <- small_barcode_tree()
  n_rep <- 4000
  r <- 2; h <- 0.5; start <- 0.8; end <- 1.5
  # window [0.8, 1.5) lies inside the tree and straddles the internal node at t = 1;
  # decay is measured from the origin, as for any profile
  cum <- r * h / log(2) * (2^(-start / h) - 2^(-end / h))
  expected <- 1 - exp(-cum)
  withr::local_seed(26)
  out <- sim_barcode_generic(tr, n_barcodes = n_rep, n_sites = 1, edit_probs = 1,
                             edit_rate = barcode_rate(r, rate_product(rate_window(start, end), rate_exp_decay(h))))
  expect_within(edited_fraction(out), expected, tolerance_4se(expected, 3 * n_rep))
})
