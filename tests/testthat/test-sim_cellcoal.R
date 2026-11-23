# CellCoal on a tree of scTreeSim. The runs are skipped when CellCoal is not installed (helper-cellcoal.R).

# Small tree from helper-barcodes.R: ((a:1, b:1):0.5, c:1.5) with a root edge of 0.5, time lengths
tree <- small_barcode_tree()

# the tree CellCoal should read: distances per haplotype, no root edge
haplotype_tree <- function(rate) {
  phy <- tree@phylo
  phy$edge.length <- rate * phy$edge.length / 2
  phy$root.edge <- NULL
  phy
}

# ---- the tree ----

test_that("CellCoal reads the tree with half the distances of the clock", {
  skip_if_no_cellcoal()
  dir <- withr::local_tempdir()
  run <- sim_cellcoal(tree, "-n1 -l1000 -6 --no-outgroup -#1", clock = clock_strict(0.01), dir = file.path(dir, "run"))
  expect_true(ape::all.equal.phylo(cellcoal_tree(run), haplotype_tree(0.01)))
  expect_equal(run$cells, data.frame(cell = c("a", "b", "c"), node = 1:3))
  expect_equal(run$tree@data$distance, 0.01 * c(1, 1, 1.5, 0.5, 0.5))

  # with the outgroup, the tips are the cells and outgcell
  run <- sim_cellcoal(tree, "-n1 -l1000 -6 -#1", clock = clock_strict(0.01), dir = file.path(dir, "outgroup"))
  expect_setequal(cellcoal_tree(run)$tip.label, c("a", "b", "c", "outgcell"))
})

test_that("a stored clock draw is used with clock = NULL", {
  skip_if_no_cellcoal()
  dir <- withr::local_tempdir()
  clocked <- clock_tree(tree, clock_relaxed(0.01, sd = 0.5))
  run <- sim_cellcoal(clocked, "-n1 -l1000 -6 --no-outgroup -#1", clock = NULL, dir = file.path(dir, "run"))
  expected <- tree@phylo
  expected$edge.length <- .node_distances(clocked)[expected$edge[, 2]] / 2
  expected$root.edge <- NULL
  expect_true(ape::all.equal.phylo(cellcoal_tree(run), expected))
})

# ---- the trunk ----

test_that("the root edge is the trunk of CellCoal", {
  skip_if_no_cellcoal()
  dir <- withr::local_tempdir()
  # the trunk is -k times the mean depth of the tips: 0.5 / mean(1.5, 1.5, 1.5)
  run <- sim_cellcoal(tree, "-n1 -l100 -1 -#1", clock = clock_strict(0.001), dir = file.path(dir, "root"))
  expect_equal(as.numeric(sub("^-k", "", grep("^-k", run$args, value = TRUE))), 0.5 / 1.5)

  no_root <- tree
  no_root@phylo$root.edge <- NULL
  run <- sim_cellcoal(no_root, "-n1 -l100 -1 -#1", clock = clock_strict(0.001), dir = file.path(dir, "no_root"))
  expect_true("-k0" %in% run$args)

  run <- sim_cellcoal(tree, "-n1 -l100 -1 -k0.2 -#1", clock = clock_strict(0.001), root_edge = FALSE,
                      dir = file.path(dir, "cellcoal"))
  expect_equal(grep("^-k", run$args, value = TRUE), "-k0.2")
})

test_that("a branch of distance d carries d * sites mutations on average, the root edge included", {
  skip_if_no_cellcoal()
  skip_on_ci()
  dir <- withr::local_tempdir()
  rate <- 0.002
  n_sites <- 10000
  n_rep <- 200
  # total time 4.5 (with the root edge 0.5): 90 mutations expected, standard error about 0.7
  run <- sim_cellcoal(tree, sprintf("-n%d -l%d -1 --no-outgroup -#7 -y0", n_rep, n_sites),
                      clock = clock_strict(rate), dir = file.path(dir, "run"))
  snvs <- vapply(seq_len(n_rep), function(i) length(attr(cellcoal_genotypes(run, i), "sites")), numeric(1))
  expect_equal(mean(snvs), rate * 4.5 * n_sites, tolerance = 0.05)
})

# ---- checks ----

test_that("the tree must hold sampled cells with labels that CellCoal reads", {
  dead <- tree
  dead@data$status[1] <- 0
  expect_error(sim_cellcoal(dead), "dead cells")

  bad <- tree
  bad@phylo$tip.label[1] <- "cell a"
  expect_error(sim_cellcoal(bad), "tip labels")
  bad@phylo$tip.label[1] <- "outgcell"
  expect_error(sim_cellcoal(bad), "outgroup")
  bad@phylo$tip.label[1] <- "b"
  expect_error(sim_cellcoal(bad), "unique")
  expect_error(sim_cellcoal(tree@phylo), "treedata")
})

test_that("the options must not give another tree or trunk", {
  expect_error(sim_cellcoal(tree, "-n1 -Tother.nwk"), "remove `-T`")
  expect_error(sim_cellcoal(tree, c("--tree", "file=other.nwk")), "remove `-T`")
  expect_error(sim_cellcoal(tree, config = list(tree_file = "other.nwk")), "remove `-T`")
  expect_error(sim_cellcoal(tree, "-n1 -k0.2"), "remove `-k`")
  expect_error(sim_cellcoal(tree, c("--tree", "transform=0.2")), "remove `-k`")

  skip_if_no_cellcoal()
  dir <- withr::local_tempdir()
  expect_warning(sim_cellcoal(tree, "-n1 -l100 -u1e-6 -1 -#1", clock = clock_strict(0.001),
                              dir = file.path(dir, "run")), "no effect")
})
