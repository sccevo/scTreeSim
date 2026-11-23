# Interface to CellCoal. The tests that run CellCoal are skipped when it is not installed (cellcoal_install()).

quickstart <- "-n1 -s6 -l1000 -u1e-6 -e10000 -1 -6 -#1"

skip_if_no_cellcoal <- function() {
  skip_on_cran()
  skip_if_not(cellcoal_available(), "CellCoal is not installed")
}

# ---- splitting a command line ----

test_that("a command line is split like a shell, without expanding anything", {
  expect_equal(.cellcoal_split(quickstart), c("-n1", "-s6", "-l1000", "-u1e-6", "-e10000", "-1", "-6", "-#1"))
  expect_equal(.cellcoal_split("  -n1\t -P0.2 0.02 "), c("-n1", "-P0.2", "0.02"))
  expect_equal(.cellcoal_split("-T'my tree.nwk' --ado \"mean=0.1,block=3000\""),
               c("-Tmy tree.nwk", "--ado", "mean=0.1,block=3000"))
  expect_equal(.cellcoal_split("-Tmy\\ tree.nwk; $HOME"), c("-Tmy tree.nwk;", "$HOME"))
  expect_equal(.cellcoal_split(""), character())
  expect_error(.cellcoal_split("-T'tree.nwk"), "Unclosed quote")
})

# ---- finding CellCoal ----

test_that("a path set by the user must exist", {
  withr::local_options(scTreeSim.cellcoal = file.path(tempdir(), "no-cellcoal"))
  expect_error(cellcoal_path(), "CellCoal is not at")
  expect_false(cellcoal_available())
})

test_that("the docs list the manual", {
  expect_output(links <- cellcoal_docs(), "manual")
  expect_true(all(c("repository", "manual", "paper") %in% names(links)))
})

# ---- running CellCoal ----

test_that("the quickstart runs from a command line and from a vector of arguments", {
  skip_if_no_cellcoal()
  expect_match(cellcoal_version(), "^[0-9]+\\.[0-9]+")

  run <- cellcoal_run(quickstart)
  expect_s3_class(run, "cellcoal_run")
  expect_true(run$ok)
  expect_equal(run$status, 0L)
  expect_true(all(c("log", "run.json", "snv_gen", "trees") %in% run$files))
  expect_equal(run$results, file.path(run$dir, "results"))
  expect_equal(tail(run$args, 1), "-oresults")

  # 6 cells and the healthy outgroup
  tree <- ape::read.tree(file.path(run$results, "trees"))
  expect_equal(ape::Ntip(tree), 7)
  expect_true("outgcell" %in% tree$tip.label)

  # the same seed gives the same tree; the program name in front is dropped
  again <- cellcoal_run(c("cellcoal", .cellcoal_split(quickstart)))
  expect_equal(readLines(file.path(again$results, "trees")), readLines(file.path(run$results, "trees")))
  expect_output(print(run), "Results:")
})

test_that("the run record is read", {
  skip_if_no_cellcoal()
  skip_if_not_installed("jsonlite")
  run <- cellcoal_run(quickstart)
  expect_equal(run$run$seed, 1)
  expect_equal(run$run$cells, 6)
})

test_that("scTreeSim sets the results folder", {
  skip_if_no_cellcoal()
  expect_warning(run <- cellcoal_run(paste(quickstart, "-oMine")), "results folder")
  expect_true("trees" %in% run$files)
  expect_false(dir.exists(file.path(run$dir, "Mine")))

  expect_warning(run <- cellcoal_run(c("-n1", "-s4", "-#1", "--output", "dir=Mine,trees")), "results folder")
  expect_true("trees" %in% run$files)
})

test_that("a failed run stops with the error of CellCoal, or returns it", {
  skip_if_no_cellcoal()
  expect_error(cellcoal_run("-n1 --bogus"), "Unknown option: --bogus")
  run <- cellcoal_run("-n1 --bogus", error = FALSE)
  expect_false(run$ok)
  expect_match(run$errors, "Unknown option")
  expect_true(is.na(run$results))
})

test_that("the run folder must be new or empty", {
  skip_if_no_cellcoal()
  dir <- withr::local_tempdir()
  writeLines("x", file.path(dir, "keep.txt"))
  expect_error(cellcoal_run(quickstart, dir = dir), "not empty")
  expect_true(file.exists(file.path(dir, "keep.txt")))
})

test_that("check validates the options without simulating", {
  skip_if_no_cellcoal()
  check <- cellcoal_check(quickstart)
  expect_true(check$ok)
  expect_length(check$files, 0)
  expect_output(print(check), "the options are valid")

  check <- cellcoal_check("-n1 --bogus")
  expect_false(check$ok)
  expect_match(check$errors, "Unknown option")
})

test_that("a user tree in another folder is copied into the run folder", {
  skip_if_no_cellcoal()
  tree_dir <- withr::local_tempdir()
  tree_file <- file.path(tree_dir, "user.nwk")
  writeLines("((D:0.03,B:0.04):0.02,(C:0.06,A:0.07):0.05);", tree_file)

  for (args in list(c("-n1", "-l100", "-6", "-#1", paste0("-T", tree_file)),
                    c("-n1", "-l100", "-#1", "--tree", paste0("file=", tree_file), "--output", "trees"))) {
    run <- cellcoal_run(args)
    expect_true(file.exists(file.path(run$dir, "user.nwk")))
    expect_true(any(args != run$args[seq_along(args)]))
    tree <- ape::read.tree(file.path(run$results, "trees"))
    expect_setequal(tree$tip.label, c("A", "B", "C", "D", "outgcell"))
  }
})

test_that("a config is read from a file or written from a list, and the arguments override it", {
  skip_if_no_cellcoal()
  skip_if_not_installed("yaml")
  skip_if_not_installed("jsonlite")

  run <- cellcoal_run(config = list(replicates = 1, cells = 5, sites = 200, print_trees = TRUE, seed = 2))
  expect_true(file.exists(file.path(run$dir, "cellcoal_config.yaml")))
  expect_equal(ape::Ntip(ape::read.tree(file.path(run$results, "trees"))), 6)

  config <- withr::local_tempfile(fileext = ".yaml")
  writeLines(c("replicates: 1", "cells: 5", "sites: 200", "print_trees: true", "seed: 2"), config)
  run <- cellcoal_run("-s8", config = config)
  expect_equal(run$run$cells, 8)

  # files named in a config are passed in `files`
  tree_file <- file.path(withr::local_tempdir(), "t.nwk")
  writeLines("((A:0.001,B:0.001):0.001,C:0.002);", tree_file)
  run <- cellcoal_run(config = list(replicates = 1, sites = 1000, tree_file = "t.nwk", print_trees = TRUE, seed = 1),
                      files = tree_file)
  expect_setequal(ape::read.tree(file.path(run$results, "trees"))$tip.label, c("A", "B", "C", "outgcell"))
  expect_false(any(grepl("\\.0$", readLines(file.path(run$dir, "cellcoal_config.yaml")))))
})

test_that("an error of the simulation is reported", {
  skip_if_no_cellcoal()
  tree_file <- file.path(withr::local_tempdir(), "long.nwk")
  writeLines("((A:1,B:1):1,C:2);", tree_file)
  run <- cellcoal_run(c("-n1", "-l100", "-#1", paste0("-T", tree_file)), error = FALSE)
  expect_false(run$ok)
  expect_match(run$errors, "infinite sites model \\(ISM\\) has been violated")
})
