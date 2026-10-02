# GESTALT-like barcodes

# example barcode from GABI's XML
xml_barcode <- c(
  "CG",
  "GATACGATACGCGCACGCTATGG",
  "AGTC",
  "GACACGACTCGCGCATACGATGG",
  "AGTC",
  "GATAGTATGCGTATACGCTATGG",
  "AGTC",
  "GATATGCATAGCGCATGCTATGG",
  "GAAAAAAAAAAAAAAA"
)

test_that(".gestalt_layout reproduces manually computed positions", {
  layout <- .gestalt_layout(xml_barcode, cut_site = 6, crucial_pos = c(6, 6))

  expect_equal(layout$n, 4)
  expect_equal(layout$length, 122)
  expect_equal(layout$targets$target, 0:3)
  expect_equal(layout$targets$cut_pos, c(19, 46, 73, 100))
  expect_equal(layout$targets$prot_left, c(13, 40, 67, 94))
  expect_equal(layout$targets$prot_right, c(24, 51, 78, 105))
  expect_equal(layout$targets$Lmax, c(19, 26, 26, 26))
  expect_equal(layout$targets$Llong, c(20, 22, 22, 22))
  expect_equal(layout$targets$Rmax, c(26, 26, 26, 21))
  expect_equal(layout$targets$Rlong, c(22, 22, 22, 22))
})

test_that(".gestalt_layout handles a single target: boundary +1 rule and clipping", {
  # prefix "AAAA" (4), target "TTTTTTTT" (8), suffix "GGGG" (4); length 16, cut_site 3
  # cut_pos = 4 + 8 - 3 = 9
  layout <- .gestalt_layout(c("AAAA", "TTTTTTTT", "GGGG"), cut_site = 3, crucial_pos = c(2, 2))

  expect_equal(layout$n, 1)
  expect_equal(layout$length, 16)
  expect_equal(layout$targets$cut_pos, 9)
  expect_equal(layout$targets$prot_left, 7)
  expect_equal(layout$targets$prot_right, 10)
  # no neighbour on either side: long-trim min = max trim + 1 (boundary rule)
  expect_equal(layout$targets$Lmax, 9)
  expect_equal(layout$targets$Llong, 10)
  expect_equal(layout$targets$Rmax, 6)
  expect_equal(layout$targets$Rlong, 7)

  # a wide crucial_pos clips the protected region to the barcode ends
  layout_clipped <- .gestalt_layout(c("AAAA", "TTTTTTTT", "GGGG"), cut_site = 3, crucial_pos = c(20, 20))
  expect_equal(layout_clipped$targets$prot_left, 0)
  expect_equal(layout_clipped$targets$prot_right, 15)
})


tract_cols <- c("min_deac", "min_cut", "max_cut", "max_deac", "left_long", "right_long", "focal")

test_that(".gestalt_tracts finds 25 tracts when all 4 targets are active", {
  tracts <- .gestalt_tracts(active = 0:3, n_targets = 4)

  expect_named(tracts, tract_cols)
  expect_equal(nrow(tracts), 25)
  expect_equal(nrow(unique(tracts[, c("min_deac", "min_cut", "max_cut", "max_deac")])), 25)
  # focal <=> min_cut == max_cut; left/right long <=> deac differs from cut
  expect_equal(tracts$focal, tracts$min_cut == tracts$max_cut)
  expect_equal(tracts$left_long, tracts$min_deac != tracts$min_cut)
  expect_equal(tracts$right_long, tracts$max_deac != tracts$max_cut)
  # no left-long trim starting at target 0, no right-long trim ending at target 3 (barcode ends)
  expect_false(any(tracts$min_cut == 0 & tracts$left_long))
  expect_false(any(tracts$max_cut == 3 & tracts$right_long))
})

test_that(".gestalt_tracts: only target 1 (of 4) active gives its 4 focal tracts", {
  tracts <- .gestalt_tracts(active = 1, n_targets = 4)

  expect_equal(nrow(tracts), 4)
  expect_true(all(tracts$min_cut == 1 & tracts$max_cut == 1 & tracts$focal))
  expect_equal(sort(tracts$min_deac), c(0, 0, 1, 1))
  expect_equal(sort(tracts$max_deac), c(1, 1, 2, 2))
})

test_that(".gestalt_tracts: a single target has no long trims at all", {
  tracts <- .gestalt_tracts(active = 0, n_targets = 1)

  expect_equal(nrow(tracts), 1)
  expect_equal(unlist(tracts[tract_cols], use.names = FALSE),
               c(0, 0, 0, 0, FALSE, FALSE, TRUE))
})


# helper function for extracting a target tract from list
tract_row <- function(tracts, min_deac, min_cut, max_cut, max_deac) {
  which(tracts$min_deac == min_deac & tracts$min_cut == min_cut &
          tracts$max_cut == max_cut & tracts$max_deac == max_deac)
}

test_that(".gestalt_hazards reproduces hand-computed values for individual tracts", {
  tracts <- .gestalt_tracts(active = 0:3, n_targets = 4)
  hazards <- .gestalt_hazards(tracts, cut_rates = c(1, 2, 3, 4),
                               long_trim_factors = c(0.1, 0.2), double_cut_weight = 0.05)

  # focal at target 1, left long trim only: lambda_1 * f_L^1 * f_R^0
  expect_equal(hazards[tract_row(tracts, 0, 1, 1, 1)], 2 * 0.1)
  # focal at target 2, both long trims: lambda_2 * f_L * f_R
  expect_equal(hazards[tract_row(tracts, 1, 2, 2, 3)], 3 * 0.1 * 0.2)
  # double cut targets 0-1, no long trims: (lambda_0 + lambda_1) * w
  expect_equal(hazards[tract_row(tracts, 0, 0, 1, 1)], (1 + 2) * 0.05)
  # double cut targets 1-2, both long trims: (lambda_1 + lambda_2) * w * f_L * f_R
  expect_equal(hazards[tract_row(tracts, 0, 1, 2, 3)], (2 + 3) * 0.05 * 0.1 * 0.2)
})


# use example layout and tracts for the following tests
xml_layout <- .gestalt_layout(xml_barcode, cut_site = 6, crucial_pos = c(6, 6))
xml_tracts <- .gestalt_tracts(active = 0:3, n_targets = 4)

# example parameters set from XML (means converted from log to natural scale)
repair_args <- list(
  insert_zero_prob = 0.6, insert_mean = 0.69,
  trim_zero_probs = matrix(c(0.081, 0.52, 0.095, 0.31), 2),
  trim_short_mean = c(1.49, 2.72), trim_long_mean = c(8.04, 3.58)
)

# helper function to generate multiple repair outcomes per tract
repair_many <- function(tract, args, n = 200) {
  do.call(rbind, replicate(n, do.call(.gestalt_repair, c(list(tract = tract, layout = xml_layout), args)),
                           simplify = FALSE))
}

test_that(".gestalt_repair: short deletions stay below the long-trim minimum, long deletions stay within bounds", {
  withr::local_seed(1)
  # focal at target 1: left long (forced), right short (optional)
  tract <- xml_tracts[tract_row(xml_tracts, 0, 1, 1, 1), ]
  Llong1 <- xml_layout$targets$Llong[2]
  Lmax1 <- xml_layout$targets$Lmax[2]
  Rlong1 <- xml_layout$targets$Rlong[2]

  draws <- repair_many(tract, repair_args)

  expect_true(all(draws$left_del >= Llong1 & draws$left_del <= Lmax1))
  expect_true(all(draws$right_del == 0 | draws$right_del <= Rlong1 - 1))
  expect_true(any(draws$right_del > 0))  # the optional short trim does occur sometimes
})

test_that(".gestalt_repair: a long side always gets a deletion", {
  withr::local_seed(2)
  # focal at target 3: left long (forced); right has no long option (barcode
  # end), only the optional short trim, using the full 1..Rmax range (boundary rule)
  tract <- xml_tracts[tract_row(xml_tracts, 2, 3, 3, 3), ]
  Rmax3 <- xml_layout$targets$Rmax[4]
  draws <- repair_many(tract, repair_args, n = 50)
  expect_true(all(draws$left_del > 0))
  expect_true(all(draws$right_del >= 0 & draws$right_del <= Rmax3))
})

test_that(".gestalt_repair: a focal cut with no long trim always leaves a trace", {
  withr::local_seed(3)
  # focal at target 1, no long trims; near-zero occurrence probabilities exercise the redraw loop
  tract <- xml_tracts[tract_row(xml_tracts, 1, 1, 1, 1), ]
  args <- list(insert_zero_prob = 0.99, insert_mean = 0.5,
               trim_zero_probs = matrix(0.99, 2, 2),
               trim_short_mean = c(0.5, 0.5), trim_long_mean = c(0.5, 0.5))
  draws <- repair_many(tract, args, n = 100)
  expect_true(all(draws$left_del > 0 | draws$right_del > 0 | nchar(draws$insert) > 0))
})

test_that(".gestalt_repair: a double cut can leave no side deletion and no insertion", {
  withr::local_seed(4)
  # double cut targets 0-1, no long trims
  tract <- xml_tracts[tract_row(xml_tracts, 0, 0, 1, 1), ]
  args <- list(insert_zero_prob = 1, insert_mean = 0.5,
               trim_zero_probs = matrix(1, 2, 2),
               trim_short_mean = c(0.5, 0.5), trim_long_mean = c(0.5, 0.5))
  draws <- repair_many(tract, args, n = 20)
  expect_true(all(draws$left_del == 0 & draws$right_del == 0 & draws$insert == ""))
})

test_that(".gestalt_repair: insertions are lowercase acgt of length 1 + Poisson(mean)", {
  withr::local_seed(5)
  tract <- xml_tracts[tract_row(xml_tracts, 0, 0, 1, 1), ]
  args <- list(insert_zero_prob = 0, insert_mean = 3,
               trim_zero_probs = matrix(1, 2, 2),
               trim_short_mean = c(0.5, 0.5), trim_long_mean = c(0.5, 0.5))
  draws <- repair_many(tract, args)
  expect_true(all(nchar(draws$insert) >= 1))
  expect_true(all(grepl("^[acgt]+$", draws$insert)))
  expect_gt(mean(nchar(draws$insert)), 1)  # mean length = 1 + insert_mean = 4
})


# resets the allele state to the original barcode
empty_allele <- function() list(deleted = rep(FALSE, xml_layout$length), insertions = list())

test_that(".gestalt_apply: focal cut deletes a contiguous span around its own cut site", {
  # focal at target 1 (cut_pos 46): left_del 5, right_del 3, no insertion
  tract <- xml_tracts[tract_row(xml_tracts, 0, 1, 1, 1), ]
  indel <- data.frame(left_del = 5, right_del = 3, insert = "", stringsAsFactors = FALSE)

  allele <- .gestalt_apply(empty_allele(), tract, indel, xml_layout)

  expect_equal(which(allele$deleted), 42:49)  # 0-based 41:48 -> R index 42:49
  expect_equal(length(allele$insertions), 0)
})

test_that(".gestalt_apply: double cut deletes everything between the two cuts, plus outer sides", {
  # double cut targets 1-2 (cut_pos 46, 73): left_del 2, right_del 0, insert "cg"
  tract <- xml_tracts[tract_row(xml_tracts, 1, 1, 2, 2), ]
  indel <- data.frame(left_del = 2, right_del = 0, insert = "cg", stringsAsFactors = FALSE)

  allele <- .gestalt_apply(empty_allele(), tract, indel, xml_layout)

  expect_equal(which(allele$deleted), 45:73)  # 0-based 44:72 -> R index 45:73
  expect_equal(allele$insertions[["46"]], "cg")
})

test_that(".gestalt_apply: a double cut with no outer deletion still deletes the middle", {
  # double cut targets 0-1 (cut_pos 19, 46), no outer deletion, no insertion
  tract <- xml_tracts[tract_row(xml_tracts, 0, 0, 1, 1), ]
  indel <- data.frame(left_del = 0, right_del = 0, insert = "", stringsAsFactors = FALSE)

  allele <- .gestalt_apply(empty_allele(), tract, indel, xml_layout)

  expect_equal(which(allele$deleted), 20:46)  # 0-based 19:45 (between the two cuts) -> R index 20:46
})

test_that(".gestalt_apply: a later double cut removes an insertion it now spans", {
  # first, an earlier focal cut at target 1 (cut_pos 46) leaves an insertion there
  tract1 <- xml_tracts[tract_row(xml_tracts, 1, 1, 1, 1), ]
  indel1 <- data.frame(left_del = 0, right_del = 0, insert = "aaa", stringsAsFactors = FALSE)
  allele <- .gestalt_apply(empty_allele(), tract1, indel1, xml_layout)
  expect_equal(allele$insertions[["46"]], "aaa")

  # later, a double cut spanning targets 0 and 2 engulfs target 1's insertion
  tract2 <- xml_tracts[tract_row(xml_tracts, 0, 0, 2, 2), ]
  indel2 <- data.frame(left_del = 0, right_del = 0, insert = "", stringsAsFactors = FALSE)
  allele2 <- .gestalt_apply(allele, tract2, indel2, xml_layout)

  expect_equal(length(allele2$insertions), 0)
  expect_true(all(allele2$deleted[20:73]))  # 0-based 19:72 -> R index 20:73
})

test_that(".gestalt_apply preserves deletions and insertions already present in the allele", {
  tract <- xml_tracts[tract_row(xml_tracts, 0, 0, 0, 0), ]
  indel <- data.frame(left_del = 3, right_del = 0, insert = "", stringsAsFactors = FALSE)

  prior <- empty_allele()
  prior$deleted[100] <- TRUE
  prior$insertions[["73"]] <- "tt"

  allele <- .gestalt_apply(prior, tract, indel, xml_layout)

  expect_true(allele$deleted[100])
  expect_equal(allele$insertions[["73"]], "tt")
  expect_true(all(allele$deleted[17:19]))  # new deletion also present (0-based 16:18)
})


# remaining parameters from example
xml_params <- c(
  list(cut_rates = rep(1, 4), long_trim_factors = c(0.1, 0.1), double_cut_weight = 0.05, clock_rate = 1),
  repair_args
)

# function with aggregated parameters
evolve_branch <- function(state, branch_length, params = xml_params) {
  do.call(.gestalt_evolve_branch, c(list(state = state, branch_length = branch_length, layout = xml_layout), params))
}

# state at origin
full_active_state <- function() list(active = 0:3, allele = empty_allele())

# Target status check: 
# a target's protected region is untouched while active, and shows an edit once deactivated; 
# an edit is either a deletion overlapping the region or an insertion anchored inside it
check_target_invariant <- function(active, allele, layout) {
  targets <- layout$targets
  ins_pos <- as.integer(names(allele$insertions))
  for (t in targets$target) {
    region <- targets$prot_left[t + 1]:targets$prot_right[t + 1]  # 0-based
    touched <- any(allele$deleted[region + 1]) || any(ins_pos %in% region)
    if (touched == (t %in% active)) return(FALSE)
  }
  TRUE
}

test_that(".gestalt_evolve_branch: a zero-length branch leaves the state unchanged", {
  state <- full_active_state()
  out <- evolve_branch(state, branch_length = 0)
  expect_equal(out, state)
})

test_that(".gestalt_evolve_branch: no active targets leaves the state unchanged, however long the branch", {
  state <- list(active = integer(0), allele = empty_allele())
  out <- evolve_branch(state, branch_length = 100)
  expect_equal(out, state)
})

test_that(".gestalt_evolve_branch: active targets keep an untouched protected region, deactivated ones show an edit", {
  withr::local_seed(1)
  for (i in 1:30) {
    out <- evolve_branch(full_active_state(), branch_length = 3)
    expect_true(check_target_invariant(out$active, out$allele, xml_layout))
  }
})

test_that(".gestalt_evolve_branch: time to the first event on a branch ~ Exp(clock_rate * sum(hazards))", {
  withr::local_seed(2)
  tracts <- .gestalt_tracts(0:3, 4)
  H <- sum(.gestalt_hazards(tracts, xml_params$cut_rates, xml_params$long_trim_factors, xml_params$double_cut_weight))
  clock_rate <- 2
  branch_length <- 0.01
  n_rep <- 4000

  params <- xml_params
  params$clock_rate <- clock_rate
  events <- vapply(seq_len(n_rep), function(i) {
    out <- evolve_branch(full_active_state(), branch_length, params)
    length(out$active) < 4
  }, logical(1))

  expected_p <- 1 - exp(-clock_rate * H * branch_length)
  expect_gt(stats::prop.test(sum(events), n_rep, p = expected_p)$p.value, 0.001)
})


# rendering
test_that(".gestalt_render places insertions at the left cut and keeps segments space-separated", {
  barcode <- c("AAAA", "TTTTTTTT", "GGGG")
  allele <- list(deleted = rep(FALSE, 16), insertions = list())
  allele$deleted[1:2] <- TRUE  # 0-based positions 0-1 -> R indices 1-2
  allele$insertions[["4"]] <- "cg"  # inserted right before position 4 (first base of the target)

  expect_equal(.gestalt_render(allele, barcode),
               paste("--AA", paste0("cg", strrep("T", 8)), "GGGG"))
})

test_that(".gestalt_render: an unedited allele reproduces the original barcode", {
  barcode <- c("AAAA", "TTTTTTTT", "GGGG")
  allele <- list(deleted = rep(FALSE, 16), insertions = list())
  expect_equal(.gestalt_render(allele, barcode), paste(barcode, collapse = " "))
})


# simulate tree for barcode evolution along lineages
withr::with_seed(10, {
  adb_tree <- sim_adb_origin_samp(5, scale = 1, shape = 1, death_prob = 0.1, sampling_prob = 0.5)
})

# join segments into a barcode for the full simulation
barcode <- paste(xml_barcode, collapse = " ")

test_that("sim_barcode_gestalt returns one edited sequence per node in the expected format", {
  withr::local_seed(11)
  out <- sim_barcode_gestalt(adb_tree, barcode, cut_rates = rep(0.3, 4))

  expect_named(out, c("node", "sequence"))
  expect_equal(out$node, tibble::as_tibble(adb_tree)$node)
  tokens <- strsplit(out$sequence, " ", fixed = TRUE)
  expect_true(all(lengths(tokens) == length(xml_barcode)))
  expect_true(all(grepl("^[ACGTacgt -]*$", out$sequence)))
})

test_that("sim_barcode_gestalt is reproducible with set.seed", {
  out1 <- withr::with_seed(42, sim_barcode_gestalt(adb_tree, barcode, cut_rates = rep(0.3, 4)))
  out2 <- withr::with_seed(42, sim_barcode_gestalt(adb_tree, barcode, cut_rates = rep(0.3, 4)))
  expect_identical(out1, out2)
})

test_that("sim_barcode_gestalt requires an ultrametric (sampled) tree", {
  complete_tree <- withr::with_seed(1, sim_adb_ntaxa_complete_fast(20, scale = 1, shape = 1, death_prob = 0.2))
  expect_error(sim_barcode_gestalt(complete_tree, barcode, cut_rates = rep(0.3, 4)), "must be ultrametric")
})

test_that("sim_barcode_gestalt validates cut_rates length against the number of targets", {
  expect_error(sim_barcode_gestalt(adb_tree, barcode, cut_rates = rep(0.3, 3)), "cut_rates")
})

test_that("sim_barcode_gestalt: a single cut_rate is recycled to all targets", {
  out_single <- withr::with_seed(7, sim_barcode_gestalt(adb_tree, barcode, cut_rates = 0.3))
  out_full <- withr::with_seed(7, sim_barcode_gestalt(adb_tree, barcode, cut_rates = rep(0.3, 4)))
  expect_identical(out_single, out_full)
})

test_that("sim_barcode_gestalt simulates a tree with zero-length edges (a parent is simulated before its children)", {
  # tips t1-t3 hang from their parent by zero-length edges, so parent and child have the same height
  tree <- treeio::as.treedata(ape::read.tree(text = "(((t1:0,t2:0,t3:0):1,t4:1):1,t5:2);"))
  tree_df <- as.data.frame(tibble::as_tibble(tree))
  for (seed in 1:3) {
    out <- withr::with_seed(seed, sim_barcode_gestalt(tree, barcode, cut_rates = rep(2, 4)))
    expect_false(any(grepl("NA", out$sequence)))
    zero <- tree_df$branch.length == 0 & !is.na(tree_df$branch.length) & tree_df$parent != tree_df$node
    expect_identical(out$sequence[match(tree_df$node[zero], out$node)], out$sequence[match(tree_df$parent[zero], out$node)])
  }
})
test_that("sim_barcode_gestalt evolves along the root edge like an ordinary branch", {
  withr::local_seed(20)
  tree_big_root <- adb_tree
  tree_big_root@phylo$root.edge <- 5  # long root edge; with the high cut rate below, an edit during it is near-certain

  out <- sim_barcode_gestalt(tree_big_root, barcode, cut_rates = rep(5, 4))
  tree_df <- tibble::as_tibble(tree_big_root)
  root <- tree_df$node[tree_df$parent == tree_df$node]

  expect_false(out$sequence[out$node == root] == barcode)
})


# alignment and event format
# toy barcode: prefix "CG" (2bp), target 0 "GATACG" (6bp, cut_pos 6), spacer "AG" (2bp),
# target 1 "GACACG" (6bp, cut_pos 14), suffix "GA" (2bp); cut_site = 2
toy_barcode_segs <- c("CG", "GATACG", "AG", "GACACG", "GA")
toy_barcode <- paste(toy_barcode_segs, collapse = " ")

toy_sequences <- c(
  "CG GATA-- -- GACACG GA",     # right-deletion at target 0 (pos 6-9), target 1 untouched
  "CG GATA-- -- GACAttCG GA",   # same target-0 deletion; pure insertion "tt" at target 1's cut
  "CG GATACG AG GACAg-- GA"     # target 0 untouched; insertion "g" + right-deletion at target 1
)

toy_aligned <- c(
  "CG GATA-- -- GACA..CG GA",
  "CG GATA-- -- GACAttCG GA",
  "CG GATACG AG GACAg.-- GA"
)

test_that("align_sequences left-justifies inserts and pads with '.' so all sequences have equal length", {
  aligned <- align_sequences(toy_sequences, toy_barcode, cut_site = 2)

  expect_equal(aligned, toy_aligned)
  expect_length(unique(nchar(aligned)), 1)
})

test_that("align_sequences infers the cut from the insertions when cut_site is not given", {
  expect_equal(align_sequences(toy_sequences, toy_barcode), toy_aligned)
})

test_that("align_sequences warns when cut_site is inconsistent with the insertions and aligns at the insertions", {
  expect_warning(aligned <- align_sequences(toy_sequences, toy_barcode, cut_site = 3), "cut_site")
  expect_equal(aligned, toy_aligned)
})

test_that("align_sequences stops when insertions at one target are at different positions", {
  sequences <- c("CG GATACG AG GACAttCG GA", "CG GATACG AG GACACtG GA")
  expect_error(align_sequences(sequences, toy_barcode), "different positions")
})

test_that("align_sequences adds no padding to a target with no insertion in any sequence", {
  aligned <- align_sequences(toy_sequences, toy_barcode, cut_site = 2)
  target0 <- vapply(strsplit(aligned, " ", fixed = TRUE), `[`, character(1), 2)
  expect_equal(nchar(target0), rep(6, 3))  # original target-0 length, unchanged
})

test_that("convert_event_format reports events with positions in the unedited barcode coordinates", {
  events <- convert_event_format(toy_sequences, toy_barcode, cut_site = 2)

  expect_equal(events, c(
    "6_4_0_0_",
    "6_4_0_0_,14_0_1_1_tt",
    "14_2_1_1_g"
  ))
})

test_that("convert_event_format: inserted bases do not shift the positions of later events", {
  sequences <- c(
    "CG GATAgg-- AG GAC--- GA",
    "CG GATA-- AG GAC--- GA"
  )
  events <- convert_event_format(sequences, toy_barcode, cut_site = 2)
  expect_equal(events[1], "6_2_0_0_gg,13_3_1_1_")
  expect_equal(events[2], "6_2_0_0_,13_3_1_1_")
})

test_that("convert_event_format: an insertion between a left and a right deletion is one event", {
  events <- convert_event_format("CG GA--tt-- AG GACACG GA", toy_barcode, cut_site = 2)
  expect_equal(events, "4_4_0_0_tt")
})

test_that("convert_event_format: a deletion ending exactly at a target's cut is credited to it", {
  events <- convert_event_format("CG GATACG AG ----CG GA", toy_barcode, cut_site = 2)
  expect_equal(events, "10_4_1_1_")
})

test_that("convert_event_format reports no events for an unedited sequence", {
  events <- convert_event_format(c(toy_barcode, toy_sequences[2]), toy_barcode, cut_site = 2)
  expect_equal(events[1], "")
})

test_that("convert_event_format does not depend on the other sequences", {
  expect_equal(
    convert_event_format(toy_sequences, toy_barcode, cut_site = 2)[1],
    convert_event_format(toy_sequences[1], toy_barcode, cut_site = 2)
  )
})


test_that("convert_event_to_indel_format copies an event to every affected target and pads with NONE", {
  events <- "138_2_0_0_,217_90_3_6_,323_50_7_8_"
  expect_equal(
    convert_event_to_indel_format(events, 10),
    "2D+138,NONE,NONE,90D+217,90D+217,90D+217,90D+217,50D+323,50D+323,NONE"
  )
})

test_that("convert_event_to_indel_format concatenates several events at one target with '&'", {
  expect_equal(
    convert_event_to_indel_format("243_1_4_4_,248_13_4_4_", 6),
    "NONE,NONE,NONE,NONE,1D+243&13D+248,NONE"
  )
})

test_that("convert_event_to_indel_format writes insertions with their length and bases", {
  expect_equal(
    convert_event_to_indel_format(c("140_0_0_0_TACTG", "10_3_0_0_ac"), 2),
    c("5I+140+TACTG,NONE", "3D+10&2I+10+ac,NONE")
  )
})

test_that("convert_event_to_indel_format gives NONE for every target for no events", {
  expect_equal(convert_event_to_indel_format(c("", "6_4_0_0_"), 2), c("NONE,NONE", "4D+6,NONE"))
})

test_that("convert_indel_format converts sequences to the per-target indel format", {
  expect_equal(
    convert_indel_format(toy_sequences, toy_barcode, cut_site = 2),
    c("4D+6,NONE", "4D+6,2I+14+tt", "NONE,2D+14&1I+14+g")
  )
})


# another example, not simulated sequences: 
# the insertions do not sit at the cut implied by cut_site = 6, so align_sequences() warns and aligns at the position the insertions occur at
toy_barcode <- "CG GATACGATACGCGCACGCTATGG AGTC GATACGATACGCGCACGCTATGG AGTC GATACGATACGCGCACGCTATGG AGTC GATACGATACGCGCACGCTATGG GAAAAAAAAAAAAAAA"
toy_sequences <- c(
  "CG GATACGATACGCGCACGCTATGG AGTC GATACGATACGCG---------- AGTC GATACGATACGCGCACGCTATGG AGTC GATACGATACGCGC--tatGCTATGG GAAAAAAAAAAAAAAA",
  "CG GATACGAa-------------TGG AGTC GATACGATACGCG---------- AGTC GATACGATACGCGCAtaaggtc------GG AGTC GATACGATACGCGC--tatGCTATGG GAAAAAAAAAAAAAAA"
)

test_that("align_sequences reproduces example", {
  expect_warning(aligned <- align_sequences(toy_sequences, toy_barcode, cut_site = 6), "cut_site")
  expect_equal(aligned, align_sequences(toy_sequences, toy_barcode))

  expect_equal(aligned, c(
    "CG GATACGA.TACGCGCACGCTATGG AGTC GATACGATACGCG---------- AGTC GATACGATACGCGCA.......CGCTATGG AGTC GATACGATACGCGC--tatGCTATGG GAAAAAAAAAAAAAAA",
    "CG GATACGAa-------------TGG AGTC GATACGATACGCG---------- AGTC GATACGATACGCGCAtaaggtc------GG AGTC GATACGATACGCGC--tatGCTATGG GAAAAAAAAAAAAAAA"
  ))
})
