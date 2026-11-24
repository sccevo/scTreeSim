# a small simulation list with the structure documented in plot_sim_results()
make_sim <- function(n_gene = 12, n_cell = 20) {
  set.seed(1)
  gene_ids <- paste0("gene_", seq_len(n_gene))
  cell_ids <- paste0("cell_", seq_len(n_cell))
  counts <- matrix(stats::rnbinom(n_gene * n_cell, size = 2, mu = 3), n_gene, n_cell,
                   dimnames = list(gene_ids, cell_ids))
  # genes with very low counts, so that every plot has data
  counts[1, ] <- 0
  counts[2, ] <- c(rep(1, 5), rep(0, n_cell - 5))
  is_bm <- rep(c(TRUE, FALSE), length.out = n_gene)
  truth <- data.frame(
    gene = gene_ids,
    class = ifelse(is_bm, "BM", "OU"),
    class_eff = ifelse(is_bm, "BM", rep(c("OU1", "OUM"), length.out = n_gene)),
    log_g = stats::rnorm(n_gene),
    alpha = ifelse(is_bm, 0, stats::runif(n_gene, 0.1, 2)),
    sigma2 = stats::runif(n_gene, 0.1, 1),
    eta = ifelse(is_bm, 0, stats::runif(n_gene)),
    phi = stats::runif(n_gene, 0.1, 1),
    psi1 = ifelse(is_bm, 0.5, stats::runif(n_gene))
  )
  cell_total <- colSums(counts)
  list(counts = counts, truth = truth,
       cell_size = stats::setNames(stats::runif(n_cell, 0.5, 1.5), cell_ids),
       cell_size_obs = cell_total / mean(cell_total),
       cell_total = cell_total, gene_total = rowSums(counts))
}

test_that("plot_sim_tree returns a ggplot for every layout", {
  skip_if_not_installed("ggtree")
  tree <- make_tree()
  for (layout in c("rectangular", "slanted", "circular", "equal_angle")) {
    p <- plot_sim_tree(tree, layout = layout, show_tips = TRUE)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_s3_class(plot_sim_tree(tree, show_axis = FALSE), "ggplot")
})

test_that("show_tips adds one label layer, with the requested size and no legend entry", {
  skip_if_not_installed("ggtree")
  tree <- make_tree()
  without <- plot_sim_tree(tree)
  with <- plot_sim_tree(tree, show_tips = TRUE, tip_size = 4)
  expect_equal(length(with$layers), length(without$layers) + 1)
  tip_layer <- with$layers[[length(with$layers)]]
  expect_equal(tip_layer$aes_params$size, 4)
  expect_false(isTRUE(tip_layer$show.legend))
  expect_no_error(ggplot2::ggplot_build(with))
})

test_that("plot_sim_tree validates its input", {
  skip_if_not_installed("ggtree")
  expect_error(plot_sim_tree("tree"))
})

test_that("plot_prior_draws returns a faceted histogram and validates its input", {
  skip_if_not_installed("ggplot2")
  p <- plot_prior_draws(list(alpha = stats::rexp(100), sigma2 = stats::rexp(100, 2)))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(levels(p$data$parameter), c("alpha", "sigma2"))

  expect_error(plot_prior_draws(1:10), "named list")
  expect_error(plot_prior_draws(list(1:10)), "named list")
  expect_error(plot_prior_draws(list(a = "x")), "named list")
})

test_that("plot_prior_draws marks the mean of the draws and the bounds of the priors", {
  skip_if_not_installed("ggplot2")
  priors <- list(alpha = sim_prior("exponential", rate = 1, lower = 0.1),
                 sigma2 = sim_prior("uniform", min = 0, max = 2, lower = 0.5, upper = 1.5))
  set.seed(1)
  draws <- lapply(priors, function(p) rsim_prior(500, p))
  lines_of <- function(p) Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)[[1]]$data
  with_priors <- plot_prior_draws(draws, priors)
  expect_no_error(ggplot2::ggplot_build(with_priors))
  lines <- lines_of(with_priors)
  expect_equal(lines$x[lines$what == "mean"], unname(vapply(draws, mean, numeric(1))))
  expect_equal(lines$x[lines$what == "lower bound"], c(0.1, 0.5))   # one finite lower bound per parameter
  expect_equal(lines$x[lines$what == "upper bound"], 1.5)           # only sigma2 has an upper bound
  expect_equal(levels(with_priors$data$parameter), c("alpha", "sigma2"))
  # without the priors, only the mean
  no_priors <- lines_of(plot_prior_draws(draws))
  expect_equal(as.character(unique(no_priors$what)), "mean")
  expect_error(plot_prior_draws(draws, list(alpha = priors$alpha)), "sim_prior")
  expect_error(plot_prior_draws(draws, list(alpha = 1, sigma2 = 2)), "sim_prior")
})

test_that("the colours of the model classes are fixed and distinct", {
  cols <- .model_colours()
  expect_named(cols, c("BM", "OU1", "OUM"))
  expect_length(unique(cols), 3)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
})

test_that("plot_sim_results returns all the plots, and they can be built", {
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("viridis")
  plots <- plot_sim_results(make_sim())

  expect_named(plots, c("oumparams", "modelclass", "genesize", "cellsize", "countsmatrix",
                        "countspergene", "countspercell", "lowcounts", "zeros_gene", "zeros_cell"))
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("plot_sim_results shows an evenly spaced subset of genes without using the RNG", {
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("viridis")
  sim <- make_sim(n_gene = 30)
  set.seed(9); a <- stats::runif(1)
  set.seed(9)
  plots <- plot_sim_results(sim, max_genes_heatmap = 10)
  expect_equal(stats::runif(1), a)
  expect_equal(length(unique(plots$countsmatrix$data$gene)), 10)
})

test_that("the combined panels are arranged figures", {
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("viridis")
  plots <- plot_sim_results(make_sim())
  expect_s3_class(plot_sim_counts_panel(plots), "ggplot")
  expect_s3_class(plot_sim_sizes_panel(plots), "ggplot")
  # the size panel has the plot of the genes with few counts too, and draws
  panel <- plot_sim_sizes_panel(plots)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(print(panel))
  expect_match(plots$lowcounts$labels$x, "< 20 counts")
})

test_that(".require_pkgs names the missing packages", {
  expect_error(.require_pkgs("notARealPackage123"), "notARealPackage123")
  expect_silent(.require_pkgs("stats"))
})

test_that("launch_gui finds the app in the package", {
  expect_true(file.exists(system.file("app", "app.R", package = "scTreeSim")))
  expect_true(is.function(launch_gui))
})

# ---- plots of the layers ----

test_that("plot_sim_snv draws the genotypes in both formats, with errors and as a distance tree", {
  skip_if_not_installed("ggtree")
  tree <- clock_tree(make_tree(), clock_strict(0.2))
  set.seed(1)
  for (format in c("genotype", "alleles")) {
    snv <- sim_gt16_nt_diploid(tree, n_sites = 8, clock = NULL, epsilon = 0.2, delta = 0.3, format = format)
    for (args in list(list(), list(errors = TRUE), list(distance = TRUE), list(errors = TRUE, distance = TRUE))) {
      p <- do.call(plot_sim_snv, c(list(tree, snv), args))
      expect_s3_class(p, "ggplot")
      expect_no_error(ggplot2::ggplot_build(p))
    }
  }
  # one tile of 16 colours per cell and site
  for (args in list(list(genotypes = TRUE), list(genotypes = TRUE, errors = TRUE), list(genotypes = TRUE, distance = TRUE))) {
    p <- do.call(plot_sim_snv, c(list(tree, snv), args))
    expect_no_error(ggplot2::ggplot_build(p))
  }
  geno_plot <- plot_sim_snv(tree, snv, genotypes = TRUE)
  is_tile <- vapply(geno_plot$layers, function(l) inherits(l$geom, "GeomTile"), logical(1))
  expect_equal(nrow(geno_plot$layers[is_tile][[1]]$data), length(tree@phylo$tip.label) * 8)
  expect_length(unique(.genotype_cols), 16)
  # the error marks add a point layer
  with_errors <- plot_sim_snv(tree, snv, errors = TRUE)
  expect_gt(length(with_errors$layers), length(plot_sim_snv(tree, snv)$layers))
})

test_that("plot_sim_snv validates its input and needs distances for the distance tree", {
  skip_if_not_installed("ggtree")
  tree <- make_tree()
  expect_error(plot_sim_snv(tree, data.frame(node = 1)))
  snv <- sim_gt16_nt_diploid(clock_tree(tree), n_sites = 4, clock = NULL)
  expect_error(plot_sim_snv(tree, snv, distance = TRUE), "distance")
})

test_that(".error_kinds classifies the changes of the genotypes", {
  true <- list(allele1 = matrix(c("A", "A", "A", "A"), 1), allele2 = matrix(c("A", "C", "A", "C"), 1))
  obs  <- list(allele1 = matrix(c("A", "A", "G", "C"), 1), allele2 = matrix(c("A", "A", "C", "C"), 1))
  kinds <- .error_kinds(true, obs)
  # site 1 unchanged; site 2 AC read as AA (het -> hom); site 3 AA read as GC (both alleles);
  # site 4 AC read as CC (het -> hom)
  expect_equal(kinds$allele1[1, ], c(NA, NA, "both alleles", "het -> hom"))
  expect_equal(kinds$allele2[1, ], c(NA, "het -> hom", "both alleles", NA))
})

test_that("plot_sim_binary draws one or several sequences, with errors and as a distance tree", {
  skip_if_not_installed("ggtree")
  tree <- clock_tree(make_tree(), clock_strict(0.3))
  set.seed(2)
  one <- sim_binary(tree, n_sites = 20, clock = NULL, alpha = 0.1, beta = 0.1)
  set.seed(2)
  true <- sim_binary(tree, n_sites = 20, clock = NULL)
  several <- sim_binary(tree, n_sites = 20, n_sequences = 3, clock = NULL)

  for (p in list(plot_sim_binary(tree, one), plot_sim_binary(tree, one, errors = true),
                 plot_sim_binary(tree, one, distance = TRUE), plot_sim_binary(tree, several))) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  # one frame per sequence
  frames <- plot_sim_binary(tree, several)$layers
  expect_equal(nrow(frames[[length(frames)]]$data), 3)

  expect_error(plot_sim_binary(tree, data.frame(node = 1)), "sim_binary")
})

test_that("plot_sim_barcodes draws the output of the three barcode simulators", {
  skip_if_not_installed("ggtree")
  tree <- make_tree()
  probs <- c(0.5, 0.3, 0.2)
  set.seed(3)
  seq_bc <- sim_barcode_seq(tree, 2, 6, 1, probs, silencing_rate = 0.1, dropout_prob = 0.1)
  nonseq <- sim_barcode_nonseq(tree, 8, 1, probs, silencing_rate = 0.1, dropout_prob = 0.1)
  generic <- sim_barcode_generic(tree, 2, 6, barcode_rate(1, rate_window(0, 2)), probs,
                                 silencing_rate = 0.1, dropout_prob = 0.1)

  for (p in list(plot_sim_barcodes(tree, seq_bc, n_edits = 3, missing_state = 4L),
                 plot_sim_barcodes(tree, nonseq, n_edits = 3, missing_state = 4L),
                 plot_sim_barcodes(tree, generic, n_edits = 3, missing_state = "-"))) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that(".barcode_matrices gives one tips x sites matrix per barcode and marks missing sites", {
  tree <- make_tree()
  n_tips <- length(tree@phylo$tip.label)
  bc <- data.frame(node = seq_len(n_tips + 3), barcode_1 = "0_1_-", barcode_2 = "2_-_-")
  mats <- .barcode_matrices(tree, bc, "-")
  expect_length(mats, 2)
  expect_equal(dim(mats[[1]]), c(n_tips, 3))
  expect_equal(rownames(mats[[1]]), tree@phylo$tip.label)
  expect_equal(unname(mats[[1]][1, ]), c("0", "1", "missing"))

  # one column per site (sim_barcode_nonseq)
  nonseq <- data.frame(node = seq_len(n_tips + 3), site_1 = 1L, site_2 = 4L)
  expect_equal(unname(.barcode_matrices(tree, nonseq, 4L)[[1]][1, ]), c("1", "missing"))
})

# ---- plots of the true parameters of the genes ----

test_that("plot_sim_params works without counts, and with or without phi", {
  skip_if_not_installed("ggplot2")
  truth <- make_sim()$truth

  with_phi <- plot_sim_params(truth)
  expect_named(with_phi, c("oumparams", "modelclass"))
  expect_setequal(levels(with_phi$oumparams$data$param), c("alpha", "sigma2", "eta", "phi", "psi1"))

  # right after the expression is simulated there is no phi yet
  without_phi <- plot_sim_params(truth[, setdiff(names(truth), "phi")])
  expect_setequal(levels(without_phi$oumparams$data$param), c("alpha", "sigma2", "eta", "psi1"))
  for (p in c(with_phi, without_phi)) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("plot_sim_results includes the plots of plot_sim_params", {
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("viridis")
  sim <- make_sim()
  plots <- plot_sim_results(sim)
  params <- plot_sim_params(sim$truth)
  expect_equal(plots$oumparams$data, params$oumparams$data)
  expect_equal(plots$modelclass$data, params$modelclass$data)
})

# ---- tree with the expression of the tips ----

# the tile layer of a plot, and the y positions of the tips of its tree
tile_layer <- function(p) {
  is_tile <- vapply(p$layers, function(l) inherits(l$geom, "GeomTile"), logical(1))
  p$layers[[which(is_tile)]]$data
}

test_that("plot_sim_expression draws one tile row per tip, aligned with the tree, and one column per gene", {
  skip_if_not_installed("ggtree")
  tree <- make_tree()
  n_tips <- length(tree@phylo$tip.label)
  expr <- sim_expression_matrix(tree, alpha = c(0.5, 0, 1), theta = theta_shared(tree),
                                gene_names = c("a", "b", "c"))

  p <- plot_sim_expression(expr, tree)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  tiles <- tile_layer(p)
  expect_equal(nrow(tiles), 3 * n_tips)
  # the rows are the tips of the tree: same y positions, and only the tips (not the internal nodes)
  tip_y <- p$data$y[p$data$isTip]
  expect_setequal(unique(tiles$y), tip_y)
  expect_equal(length(unique(tiles$y)), n_tips)
  # three columns, side by side, to the right of the tree
  expect_equal(length(unique(tiles$x)), 3)
  expect_gt(min(tiles$x), max(p$data$x))

  # the values of a column are the expression of that gene at the tips, row by row
  tip_label <- p$data$label[p$data$isTip]
  first_col <- tiles[tiles$x == min(tiles$x), ]
  cell_nodes <- match(tip_label[match(first_col$y, tip_y)], tree@phylo$tip.label)
  expect_equal(first_col$expression, expr$a[match(cell_nodes, expr$node)])
})

test_that("plot_sim_expression groups the genes by class, with a gap, and shows a subset of many genes", {
  skip_if_not_installed("ggtree")
  tree <- make_tree()
  expr <- sim_expression_matrix(tree, theta = theta_shared(tree), n_genes = 30)
  genes <- setdiff(names(expr), c("node", "type"))
  truth <- data.frame(gene = genes, class_eff = rep(c("BM", "OU1", "OUM"), length.out = 30))

  grouped <- plot_sim_expression(expr, tree, truth, gap = 2)
  expect_no_error(ggplot2::ggplot_build(grouped))
  cols <- sort(unique(tile_layer(grouped)$x))
  expect_equal(length(cols), 30)
  # two gaps of two empty columns: the steps between the columns show them
  steps <- round(diff(cols) / min(diff(cols)))
  expect_equal(sum(steps == 3), 2)
  # a frame and a label per class
  labels <- unlist(lapply(grouped$layers, function(l) if (inherits(l$geom, "GeomText")) l$aes_params$label))
  expect_true(all(c("BM", "OU1", "OUM") %in% labels))

  subset <- plot_sim_expression(expr, tree, max_genes = 10)
  expect_equal(length(unique(tile_layer(subset)$x)), 10)
  expect_match(subset$labels$caption, "subset of 10")
})

test_that("plot_sim_expression validates its input", {
  skip_if_not_installed("ggtree")
  tree <- make_tree()
  expr <- sim_expression_matrix(tree, theta = theta_shared(tree))
  expect_error(plot_sim_expression(expr[, c("node", "type")], tree), "no gene columns")
  expect_error(plot_sim_expression(expr[expr$node != 1, ], tree), "no row for some tips")
  expect_error(plot_sim_expression(as.matrix(expr), tree))
  expect_error(plot_sim_expression(expr, "tree"))
})

test_that("trees with more than two types get one colour per type", {
  skip_if_not_installed("ggtree")
  expect_length(unique(.type_colours(as.character(0:9))), 10)
  expect_equal(unname(.type_colours(c("0", "1"))), c("#2C7BB6", "#D7191C"))
  tree <- make_tree()
  tree@data$type <- (seq_len(nrow(tree@data)) - 1L) %% 4L
  p <- plot_sim_tree(tree)
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplot_build(plot_sim_snv(clock_tree(tree, clock_strict(0.2)),
                                                    sim_gt16_nt_diploid(clock_tree(tree, clock_strict(0.2)), n_sites = 4, clock = NULL))))
})

test_that("plot_sim_params draws the proportions of the optima of any number of types", {
  skip_if_not_installed("ggplot2")
  truth <- make_sim(n_gene = 30)$truth
  truth$psi_0 <- 0.3
  truth$psi_1 <- 0.3
  truth$psi1 <- NULL
  p <- plot_sim_params(truth)
  expect_true(all(c("psi_0", "psi_1") %in% levels(p$oumparams$data$param)))
})

# ---- many cells: subsample, automatic layout, UMAP ----

test_that("the plots with a row per cell draw a random subsample, with the tree pruned to it, and say so", {
  skip_if_not_installed("ggtree")
  tree <- make_tree3(40)
  ct <- clock_tree(tree, clock_strict(0.1))
  n_total <- ape::Ntip(ct@phylo)
  snv <- sim_gt16_nt_diploid(ct, n_sites = 5, clock = NULL, epsilon = 0.1, delta = 0.1)
  tile_rows <- function(p) nrow(Filter(function(l) inherits(l$geom, "GeomTile"), p$layers)[[1]]$data)
  for (args in list(list(), list(distance = TRUE), list(errors = TRUE), list(genotypes = TRUE))) {
    p <- do.call(plot_sim_snv, c(list(ct, snv, max_cells = 12), args))
    expect_no_error(ggplot2::ggplot_build(p))
    expect_match(ggplot2::get_labs(p)$caption, sprintf("Random subsample of 12 of %d cells", n_total))
    expect_equal(tile_rows(p), 12 * 5 * (if (isTRUE(args$genotypes)) 1 else 2))
  }
  # all the cells: no subsample, no note
  full <- plot_sim_snv(ct, snv, max_cells = n_total)
  expect_null(ggplot2::get_labs(full)$caption)
  expect_equal(tile_rows(full), n_total * 5 * 2)
  # the same seed gives the same cells, another seed other cells
  cells <- function(seed) plot_sim_snv(ct, snv, max_cells = 12, seed = seed)$data$label[plot_sim_snv(ct, snv, max_cells = 12, seed = seed)$data$isTip]
  expect_identical(cells(1), cells(1))
  expect_false(identical(cells(1), cells(2)))
})

test_that("the subsample keeps the data of the cells that are shown", {
  tree <- make_tree3(40)
  ct <- clock_tree(tree, clock_strict(0.1))
  sub <- .subsample_cells(ct, 10, seed = 3)
  expect_equal(ape::Ntip(sub$tree@phylo), 10)
  expect_equal(sub$tree@phylo$tip.label, as.character(sub$keep))
  expect_equal(sub$tree@data$type[sub$tree@data$node <= 10], ct@data$type[match(sub$keep, ct@data$node)])
  expect_true("distance" %in% names(sub$tree@data))
  # every kept node has its parent: the pruned tree is connected up to a single root
  edge <- sub$tree@phylo$edge
  expect_equal(length(setdiff(edge[, 1], edge[, 2])), 1)
  expect_equal(nrow(sub$tree@data), 10 + sub$tree@phylo$Nnode)
  for (seed in 1:30) {   # many subsamples: none loses a branch (a regression of the climb to the root)
    pruned <- .subsample_cells(ct, 10, seed)$tree@phylo
    expect_equal(length(setdiff(pruned$edge[, 1], pruned$edge[, 2])), 1)
    expect_equal(nrow(pruned$edge), ape::Ntip(pruned) + pruned$Nnode - 1)
  }
  expect_equal(sub$tree@phylo$root.edge, ct@phylo$root.edge)
  # the table of a layer is renumbered like the tips of the pruned tree
  snv <- sim_gt16_nt_diploid(ct, n_sites = 4, clock = NULL)
  full <- .allele_matrices(ct, snv, "true")
  pruned <- .allele_matrices(sub$tree, lapply(snv, .rekey_cells, keep = sub$keep), "true")
  expect_equal(unname(pruned$allele1), unname(full$allele1[sub$keep, ]))
  expect_equal(rownames(pruned$allele1), rownames(full$allele1)[sub$keep])
  # the random numbers of the session are not touched
  set.seed(10); a <- stats::runif(1)
  set.seed(10); .subsample_cells(ct, 10, seed = 99); expect_equal(stats::runif(1), a)
  expect_identical(.subsample_cells(ct, 1000, 1)$tree, ct)
})

test_that("the binary, barcode and expression plots subsample the cells too", {
  skip_if_not_installed("ggtree")
  tree <- make_tree3(40)
  ct <- clock_tree(tree, clock_strict(0.2))
  set.seed(1)
  bin <- sim_binary(ct, n_sites = 8, clock = NULL, alpha = 0.1, beta = 0.1)
  true <- sim_binary(ct, n_sites = 8, clock = NULL)
  bc <- sim_barcode_seq(tree, n_barcodes = 1, n_sites = 4, edit_rate = 1, edit_probs = c(0.5, 0.5))
  genes <- sim_gene_params(tree, n_genes = 6)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  plots <- list(plot_sim_binary(ct, bin, max_cells = 9), plot_sim_binary(ct, bin, errors = true, max_cells = 9),
                plot_sim_barcodes(tree, bc, n_edits = 2, max_cells = 9), plot_sim_expression(expr, tree, genes, max_cells = 9))
  for (p in plots) {
    expect_no_error(ggplot2::ggplot_build(p))
    expect_match(ggplot2::get_labs(p)$caption, "Random subsample of 9 of")
  }
})

test_that("plot_sim_tree chooses a circular layout for big trees and hides the labels of many tips", {
  skip_if_not_installed("ggtree")
  tree <- make_tree3(40)
  circular <- plot_sim_tree(tree, max_tips_rectangular = 10)
  expect_s3_class(circular$coordinates, "CoordPolar")
  expect_false(inherits(plot_sim_tree(tree)$coordinates, "CoordPolar"))
  expect_s3_class(plot_sim_tree(tree, layout = "rectangular", max_tips_rectangular = 10)$coordinates, "CoordCartesian")
  hidden <- plot_sim_tree(tree, show_tips = TRUE, max_tip_labels = 10)
  expect_match(ggplot2::get_labs(hidden)$caption, "tip labels are not shown")
  expect_null(ggplot2::get_labs(plot_sim_tree(tree, show_tips = TRUE))$caption)
})

test_that("plot_sim_params also draws the optima of the types, and the counts matrix is subsampled in cells", {
  skip_if_not_installed("ggplot2")
  tree <- make_tree3(30)
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 40, pi_ou = 0.8)
  p <- plot_sim_params(genes)
  expect_named(p, c("oumparams", "modelclass", "thetas"))
  expect_equal(levels(p$thetas$data$type), paste("Type", 0:2))
  expect_true(all(p$thetas$data$class_eff != "BM"))
  expect_no_error(ggplot2::ggplot_build(p$thetas))
  expect_named(plot_sim_params(make_sim()$truth), c("oumparams", "modelclass"))   # no optima in the table
  skip_if_not_installed("ggpubr"); skip_if_not_installed("viridis")
  sim <- make_sim(n_cell = 60)
  big <- plot_sim_results(sim, max_cells_heatmap = 20)
  expect_equal(length(unique(big$countsmatrix$data$cell)), 20)
  expect_match(big$countsmatrix$labels$x, "evenly spaced subset of 20")
})

test_that("plot_sim_umap draws one point per cell, coloured by cell type, and is reproducible", {
  skip_if_not_installed("uwot")
  tree <- make_tree3(40)
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 30, pi_ou = 1)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  p <- plot_sim_umap(expr, tree)
  expect_equal(nrow(p$data), sum(tree@data$status == 1))
  expect_equal(levels(p$data$type), as.character(sort(unique(tree@data$type))))
  expect_equal(as.character(p$data$type), as.character(tree@data$type[match(as.integer(p$data$cell), tree@data$node)]))
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(plot_sim_umap(expr, tree)$data, p$data)               # the default seed
  expect_false(isTRUE(all.equal(plot_sim_umap(expr, tree, seed = 2)$data$UMAP1, p$data$UMAP1)))
  set.seed(5); a <- stats::runif(1)
  set.seed(5); plot_sim_umap(expr, tree); expect_equal(stats::runif(1), a)   # the session's random numbers are not touched
  expect_equal(levels(p$data$kind), c("sampled tip", "extinct tip", "internal node"))
  expect_true(all(p$data$kind == "sampled tip"))               # by default only the tip cells
  expect_error(plot_sim_umap(expr[1:3], tree), "at least 2 genes")
  expect_no_error(plot_sim_umap(expr, tree, max_genes = 5))
})

test_that("plot_sim_umap can include the internal nodes and the extinct tips, marked by kind", {
  skip_if_not_installed("uwot")
  tree <- make_tree3(40)
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 30, pi_ou = 1)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  all_nodes <- plot_sim_umap(expr, tree, tips_only = FALSE)
  expect_equal(nrow(all_nodes$data), nrow(tree@data))
  expect_equal(sum(all_nodes$data$kind == "internal node"), sum(tree@data$status == 2))
  expect_equal(sum(all_nodes$data$kind == "sampled tip"), sum(tree@data$status == 1))
  expect_equal(sum(all_nodes$data$kind == "extinct tip"), sum(tree@data$status == 0))
  expect_true(all(grepl("^node_", all_nodes$data$cell[all_nodes$data$kind == "internal node"])))
  expect_no_error(ggplot2::ggplot_build(all_nodes))
  expect_match(ggplot2::get_labs(all_nodes)$caption, "internal nodes")
  # only the tips: the sampled tips, as before
  tips <- plot_sim_umap(expr, tree)
  expect_equal(nrow(tips$data), sum(tree@data$status == 1))
  expect_match(ggplot2::get_labs(tips)$caption, "the sampled tips")
  # the types of the internal nodes are those of the tree
  internal <- all_nodes$data[all_nodes$data$kind == "internal node", ]
  ids <- as.integer(sub("node_", "", internal$cell))
  expect_equal(as.character(internal$type), as.character(tree@data$type[match(ids, tree@data$node)]))
})

test_that("the heatmap of the layers keeps a share of the plot however many sites there are", {
  skip_if_not_installed("ggtree")
  tree <- make_tree3(30)
  ct <- clock_tree(tree, clock_strict(0.1))
  heat_width <- function(p) {
    tiles <- Filter(function(l) inherits(l$geom, "GeomTile"), p$layers)[[1]]$data
    diff(range(tiles$x)) / diff(range(p$data$x))
  }
  many <- sim_gt16_nt_diploid(ct, n_sites = 300, clock = NULL)
  few <- sim_gt16_nt_diploid(ct, n_sites = 10, clock = NULL)
  expect_lte(heat_width(plot_sim_snv(ct, many)), 1.02)
  expect_lte(heat_width(plot_sim_snv(ct, many, max_heat_width = 0.5)), 0.51)
  # few sites keep the width per site: the centres of 10 sites of 0.015 of the tree (two tiles each) span 9.5 sites
  expect_equal(heat_width(plot_sim_snv(ct, few)), 0.015 * 9.5, tolerance = 0.02)
  # the binary trait and the barcodes cap their width too
  bin <- sim_binary(ct, n_sites = 600, clock = NULL)
  pb <- ggplot2::ggplot_build(plot_sim_binary(ct, bin, max_heat_width = 0.6))
  expect_no_error(pb)
  bc <- sim_barcode_seq(tree, n_barcodes = 2, n_sites = 30, edit_rate = 1, edit_probs = c(0.5, 0.5))
  expect_no_error(ggplot2::ggplot_build(plot_sim_barcodes(tree, bc, n_edits = 2, max_heat_width = 0.4)))
})

# ---- consistent titles: capitals, and Greek letters for the parameters ----

test_that("the names of the parameters are written as plotmath, and the other names are quoted", {
  expect_equal(.param_math(c("alpha", "sigma2", "eta", "phi", "psi_0", "psi1", "theta_12", "other")),
               c("alpha", "sigma^2", "eta", "phi", "psi[0]", "psi[1]", "theta[12]", "\"other\""))
  expect_equal(.cap(c("sampled tip", "edit 1", "A")), c("Sampled tip", "Edit 1", "A"))
  # every label can be parsed
  expect_no_error(lapply(.param_math(c("alpha", "sigma2", "psi_3", "theta_0", "x y")), function(x) parse(text = x)))
})

test_that("the strips and the legends of the plots use Greek letters and start with a capital", {
  skip_if_not_installed("ggplot2")
  skip_if_not(exists("get_strip_labels", asNamespace("ggplot2")), "needs a newer ggplot2")
  tree <- make_tree3(30)
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 40, pi_ou = 0.8)
  p <- plot_sim_params(genes)
  expect_equal(levels(p$oumparams$data$label),
               c("alpha", "sigma^2", "eta", "psi[0]", "psi[1]"))   # plotmath: two shares for three types
  expect_equal(as.character(ggplot2::get_strip_labels(p$thetas)$facets$type), c("Type 0", "Type 1", "Type 2"))
  expect_match(paste(deparse(p$thetas$labels$x), collapse = ""), "theta")   # an expression with the letter
  expect_equal(p$modelclass$labels$x, "Model class")

  priors <- list(alpha = sim_prior("exponential", rate = 1, lower = 0.1),
                 sigma2 = sim_prior("exponential", rate = 1, upper = 2))
  set.seed(2)
  pp <- plot_prior_draws(lapply(priors, function(x) rsim_prior(200, x)), priors)
  expect_match(levels(pp$data$title), "^(alpha|sigma\\^2) ~ \"\\(mean ")
  expect_equal(as.character(ggplot2::get_guide_data(pp, "colour")$.label), c("Mean", "Lower bound", "Upper bound"))

  skip_if_not_installed("ggtree")
  ct <- clock_tree(tree, clock_strict(0.1))
  snv <- sim_gt16_nt_diploid(ct, n_sites = 5, clock = NULL, epsilon = 0.2, delta = 0.2)
  ps <- plot_sim_snv(ct, snv, errors = TRUE)
  expect_equal(as.character(ggplot2::get_guide_data(ps, "shape")$.label),
               c("One allele", "Heterozygous to homozygous", "Both alleles"))
  expect_equal(ps$labels$x, "Time before present")
  expect_match(ggplot2::get_labs(plot_sim_snv(ct, snv, genotypes = TRUE, errors = TRUE))$caption, "^Points: ")
  bin <- sim_binary(ct, n_sites = 8, clock = NULL, alpha = 0.2, beta = 0.2)
  true <- sim_binary(ct, n_sites = 8, clock = NULL)
  pb <- plot_sim_binary(ct, bin, errors = true)
  labels <- as.character(ggplot2::get_guide_data(pb, "fill")$.label)
  expect_true(all(labels %in% c("0", "1", "False positive", "False negative")))
  bc <- sim_barcode_seq(tree, n_barcodes = 1, n_sites = 4, edit_rate = 1, edit_probs = c(0.5, 0.5))
  pc <- plot_sim_barcodes(tree, bc, n_edits = 2)
  expect_true(all(grepl("^[A-Z]", as.character(ggplot2::get_guide_data(pc, "fill")$.label))))
  expect_equal(ggplot2::get_labs(plot_sim_tree(tree))$colour, "Type")
})
