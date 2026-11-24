# the summaries are what the app reports as text; they return the numbers too

make_layers <- function() {
  tree <- make_tree3(30)
  ct <- clock_tree(tree, clock_strict(0.2))
  set.seed(1)
  list(tree = tree, ct = ct,
       snv = sim_gt16_nt_diploid(ct, n_sites = 12, clock = NULL, epsilon = 0.1, delta = 0.1),
       bin = sim_binary(ct, n_sites = 20, clock = NULL, alpha = 0.1, beta = 0.1),
       bin_true = sim_binary(ct, n_sites = 20, clock = NULL))
}

test_that("summarize_tree reports the size, height and types of the tree", {
  tree <- make_tree3(30)
  s <- summarize_tree(tree)
  expect_s3_class(s, "sim_summary")
  expect_equal(s$stats$n_tips, ape::Ntip(tree@phylo))
  expect_equal(s$stats$height, max(ape::node.depth.edgelength(tree@phylo)))
  expect_equal(sum(s$stats$tips_by_type), sum(tree@data$status == 1))
  expect_match(format(s), "^Tips: [0-9]+   Internal nodes: [0-9]+   Height: [0-9.]+\nTips by type: 0 = ")
  expect_output(print(s), "Tips by type")
  expect_error(summarize_tree(1))
})

test_that("summarize_snv reports heterozygosity and the effect of the error, in both formats", {
  l <- make_layers()
  s <- summarize_snv(l$ct, l$snv)
  expect_equal(s$stats$n_sites, 12)
  expect_equal(s$stats$n_cells, ape::Ntip(l$tree@phylo))
  expect_true(all(c("heterozygous", "variable_sites", "changed_by_error") %in% names(s$stats)))
  expect_match(format(s), "Heterozygous sites in the cells \\(true\\): [0-9.]+%")
  expect_match(format(s), "Genotypes changed by the error: [0-9.]+%")
  # no error: nothing changes
  none <- sim_gt16_nt_diploid(l$ct, n_sites = 12, clock = NULL)
  expect_equal(summarize_snv(l$ct, none)$stats$changed_by_error, 0)
  # the other format gives the same numbers
  set.seed(1)
  alleles <- sim_gt16_nt_diploid(l$ct, n_sites = 12, clock = NULL, epsilon = 0.1, delta = 0.1, format = "alleles")
  set.seed(1)
  genotypes <- sim_gt16_nt_diploid(l$ct, n_sites = 12, clock = NULL, epsilon = 0.1, delta = 0.1)
  expect_equal(summarize_snv(l$ct, alleles)$stats, summarize_snv(l$ct, genotypes)$stats)
})

test_that("summarize_binary reports the sites, and the flips of the error when the truth is given", {
  l <- make_layers()
  s <- summarize_binary(l$ct, l$bin)
  expect_equal(s$stats$n_sites, 20)
  expect_equal(s$stats$n_sequences, 1)
  expect_null(s$stats$flipped_by_error)
  expect_false(grepl("flipped", format(s)))
  with_truth <- summarize_binary(l$ct, l$bin, true = l$bin_true)
  expect_gte(with_truth$stats$flipped_by_error, 0)
  expect_match(format(with_truth), "Sites flipped by the error")
  expect_error(summarize_binary(l$ct, data.frame(node = 1)), "binary_1")
})

test_that("summarize_barcodes reports the fractions of unedited, edited and missing sites", {
  tree <- make_tree3(30)
  set.seed(2)
  bc <- sim_barcode_seq(tree, n_barcodes = 2, n_sites = 4, edit_rate = 1, edit_probs = c(0.5, 0.3, 0.2),
                        silencing_rate = 0.05)
  s <- summarize_barcodes(tree, bc, missing_state = 4)
  expect_named(s$stats$barcodes, c("barcode_1", "barcode_2"))
  for (b in s$stats$barcodes) expect_equal(b$unedited + b$edited + b$missing, 1)
  expect_match(format(s), "barcode_1: [0-9.]+% unedited, [0-9.]+% edited, [0-9.]+% missing")
  set.seed(2)
  nonseq <- sim_barcode_nonseq(tree, n_targets = 5, edit_rate = 1, edit_probs = c(0.5, 0.5))
  expect_match(format(summarize_barcodes(tree, nonseq, missing_state = 3)), "targets: ")
})

test_that("summarize_expression and summarize_counts report the genes and the cells", {
  tree <- make_tree3(30)
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 40)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  s <- summarize_expression(expr, tree, genes)
  expect_equal(s$stats$n_genes, 40)
  expect_equal(s$stats$n_cells, sum(tree@data$status == 1))
  expect_equal(sum(s$stats$class), 40)
  expect_match(format(s), "Effective class \\(OU1 / OUM by gap\\)")

  cp <- sim_count_params(expr, tree, genes, target_mean_count = 2)
  cdf <- sim_expression_counts_matrix(cp$expr, tree, phi = cp$phi, size_factor = cp$size_factor)
  counts <- t(as.matrix(cdf[genes$gene]))
  dimnames(counts) <- list(genes$gene, tree@phylo$tip.label[cdf$node])
  cs <- summarize_counts(counts, cp$gene_params, size_factor = cp$size_factor, shift = cp$shift)
  expect_equal(cs$stats$n_cells, ncol(counts))
  expect_equal(cs$stats$mean_count, mean(counts))
  expect_equal(cs$stats$fraction_zeros, mean(counts == 0))
  expect_equal(cs$stats$shift, cp$shift)
  expect_gt(cs$stats$cor_cell_size, 0.5)
  expect_match(format(cs), "Level shift applied")
  # without the true size factors and the shift, those lines are left out
  bare <- summarize_counts(counts, cp$gene_params)
  expect_false(grepl("Level shift|true cell size", format(bare)))
  expect_error(summarize_counts(counts[-1, ], cp$gene_params))
})
