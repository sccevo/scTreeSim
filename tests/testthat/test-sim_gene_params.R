# ---- gene_prior() and rgene_prior() ----

test_that("gene_prior describes the four families, by name or by position", {
  expect_s3_class(gene_prior("exponential", rate = 2), "gene_prior")
  expect_equal(gene_prior("lognormal", -2, 1)$params, list(meanlog = -2, sdlog = 1))
  expect_equal(gene_prior("normal", mean = 1, sd = 0.5)$params, list(mean = 1, sd = 0.5))
  expect_equal(gene_prior("uniform", 0, 3, lower = 1)$lower, 1)
  expect_output(print(gene_prior("exponential", rate = 2, lower = 0.1)), "exponential\\(rate = 2\\).*\\[0.1, Inf\\]")
})

test_that("gene_prior validates the parameters and the bounds", {
  expect_error(gene_prior("exponential"), "rate")
  expect_error(gene_prior("exponential", rate = -1), "positive")
  expect_error(gene_prior("exponential", sd = 1), "rate")
  expect_error(gene_prior("normal", mean = 0, sd = 0), "standard deviation")
  expect_error(gene_prior("uniform", min = 2, max = 1), "minimum")
  expect_error(gene_prior("exponential", rate = c(1, 2)), "rate")
  expect_error(gene_prior("exponential", rate = 1, lower = 2, upper = 1), "lower")
  expect_error(gene_prior("gamma"), "should be one of")
})

test_that("rgene_prior draws inside the bounds and is reproducible", {
  prior <- gene_prior("normal", mean = 0, sd = 1, lower = 0, upper = 1)
  set.seed(1)
  x <- rgene_prior(500, prior)
  expect_length(x, 500)
  expect_true(all(x >= 0 & x <= 1))
  set.seed(1)
  expect_identical(rgene_prior(500, prior), x)
  # untruncated draws are those of the base R generator
  set.seed(2); a <- rgene_prior(5, gene_prior("lognormal", meanlog = -1, sdlog = 2))
  set.seed(2); expect_identical(a, stats::rlnorm(5, -1, 2))
  expect_error(rgene_prior(10, gene_prior("uniform", min = 0, max = 1, lower = 5)), "incompatible")
  expect_error(rgene_prior(3, list(rate = 1)))
})

# ---- sim_gene_params() ----

test_that("sim_gene_params gives the same genes for the same seed", {
  tree <- make_tree3()
  set.seed(5); a <- sim_gene_params(tree, n_genes = 40)
  set.seed(5); b <- sim_gene_params(tree, n_genes = 40)
  expect_identical(a, b)
  set.seed(6); expect_false(identical(a, sim_gene_params(tree, n_genes = 40)))
})

test_that("sim_gene_params returns the table of true parameters for any number of types", {
  for (n_types in 1:3) {
    tree <- make_tree3()
    tree@data$type <- tree@data$type %% n_types
    set.seed(1)
    genes <- sim_gene_params(tree, n_genes = 60, n_types = n_types)
    expect_equal(nrow(genes), 60)
    expect_equal(genes$gene, paste0("gene_", 1:60))
    expect_true(all(paste0("theta_", seq_len(n_types) - 1) %in% names(genes)))
    expect_equal(sum(grepl("^psi_", names(genes))), n_types - 1)
    expect_setequal(unique(genes$class), c("BM", "OU"))
    expect_true(all(genes$class_eff %in% c("BM", "OU1", "OUM")))
    expect_true(all(genes$alpha[genes$class == "BM"] == 0))
    # the optima of the types of an OU gene deviate from the baseline by amounts that sum to zero
    theta <- as.matrix(genes[paste0("theta_", seq_len(n_types) - 1)])
    expect_equal(unname(rowMeans(theta)), genes$theta_base)
    if (n_types == 1) expect_false("OUM" %in% genes$class_eff)
  }
})

test_that("sim_gene_params respects the bounds of the priors and the probability of OU", {
  tree <- make_tree3()
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 300, pi_ou = 1,
                           alpha = gene_prior("exponential", rate = 1, lower = 0.5),
                           sigma2 = gene_prior("uniform", min = 0.1, max = 0.4))
  expect_true(all(genes$class == "OU"))
  expect_true(all(genes$alpha >= 0.5))
  expect_true(all(genes$sigma2 >= 0.1 & genes$sigma2 <= 0.4))
  expect_true(all(genes$sd_stat == sqrt(genes$sigma2 / (2 * genes$alpha))))
  expect_true(all(sim_gene_params(tree, n_genes = 20, pi_ou = 0)$class == "BM"))
  # eta = 0: every type has the same optimum
  flat <- sim_gene_params(tree, n_genes = 20, pi_ou = 1, eta = gene_prior("uniform", min = 0, max = 1e-6))
  expect_true(all(flat$theta_gap < 1e-5))
})

test_that("the mean of the baseline level is mean_count_gene", {
  tree <- make_tree3()
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 20000, mean_count_gene = 3, log_g_sd = 0.5)
  expect_equal(mean(exp(genes$log_g)), 3, tolerance = 0.03)
  expect_equal(stats::sd(genes$log_g), 0.5, tolerance = 0.02)
})

test_that("sim_gene_params validates its input", {
  tree <- make_tree3()
  expect_error(sim_gene_params(1), "treedata")
  expect_error(sim_gene_params(tree, n_genes = 0), "n_genes")
  expect_error(sim_gene_params(tree, pi_ou = 2), "pi_ou")
  expect_error(sim_gene_params(tree, alpha = 1), "gene_prior")
  expect_error(sim_gene_params(tree, psi_concentration = 0), "psi_concentration")
  expect_error(sim_gene_params(tree, n_types = 2), "types outside")
  expect_error(sim_gene_params(tree, mean_count_gene = -1), "mean_count_gene")
})

test_that("the parameters of sim_gene_params drive sim_expression_matrix", {
  tree <- make_tree3()
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 8)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  expect_equal(names(expr), c("node", "type", genes$gene))
  expect_equal(nrow(expr), nrow(tree@data))
})

# ---- gene_params_from_table() and gene_optima() ----

test_that("gene_params_from_table builds the table of true parameters of the given genes", {
  tree <- make_tree3()
  tab <- data.frame(gene = c("a", "b", "c"), alpha = c(0, 1.5, 2), sigma2 = c(0.5, 0.4, 0.3),
                    theta_0 = c(NA, 0, -1), theta_1 = c(NA, 0, 1), theta_2 = c(NA, 0, 2))
  genes <- gene_params_from_table(tab, tree)
  expect_equal(genes$gene, c("a", "b", "c"))
  expect_equal(genes$class, c("BM", "OU", "OU"))
  expect_equal(genes$class_eff[2], "OU1")
  expect_equal(genes$class_eff[3], "OUM")
  expect_equal(genes$theta_gap, c(0, 0, 3))
  expect_equal(genes$sd_stat[2], sqrt(0.4 / 3))
  expect_equal(genes$theta_base, c(0, 0, 2 / 3))
  # the default root: the optimum of the type of the root for OU genes, 0 for BM genes
  root_type <- tree@data$type[tree@data$node == ape::Ntip(tree@phylo) + 1]
  expect_equal(genes$x_root, c(0, 0, c(-1, 1, 2)[root_type + 1]))
  expect_equal(gene_params_from_table(transform(tab, root_value = c(1, 2, 3)), tree)$x_root, c(1, 2, 3))
  # no gene names: gene_1, ...
  expect_equal(gene_params_from_table(tab[-1], tree)$gene, paste0("gene_", 1:3))
})

test_that("gene_params_from_table validates the table", {
  tree <- make_tree3()
  ok <- data.frame(alpha = c(0, 1), sigma2 = c(0.5, 0.5), theta_0 = c(0, 1), theta_1 = c(0, -1), theta_2 = c(0, 0))
  expect_error(gene_params_from_table(ok["alpha"], tree), "sigma2")
  expect_error(gene_params_from_table(transform(ok, sigma2 = c(0.5, 0)), tree), "sigma2 must")
  expect_error(gene_params_from_table(transform(ok, alpha = c(-1, 1)), tree), "alpha must")
  expect_error(gene_params_from_table(ok[c("alpha", "sigma2", "theta_0")], tree), "theta_1")
  expect_error(gene_params_from_table(transform(ok, theta_1 = c(0, NA)), tree), "missing for some OU")
  expect_error(gene_params_from_table(transform(ok, gene = c("a", "a")), tree), "unique")
  expect_error(gene_params_from_table(transform(ok, root_value = c(1, NA)), tree), "root_value")
  expect_error(gene_params_from_table(ok[0, ], tree), "no rows")
  # only Brownian motions: the optima are not needed
  bm <- gene_params_from_table(data.frame(alpha = c(0, 0), sigma2 = c(1, 2)), tree)
  expect_equal(bm$theta_2, c(0, 0))
  expect_equal(gene_params_from_table(transform(ok, sigmasq = sigma2, sigma2 = NULL), tree)$sigma2, c(0.5, 0.5))
})

test_that("gene_optima arranges the optima as the theta of sim_expression_matrix", {
  tree <- make_tree3()
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 4)
  theta <- gene_optima(genes)
  expect_equal(names(theta), c("type", genes$gene))
  expect_equal(theta$type, 0:2)
  expect_equal(unname(unlist(theta[2, -1])), genes$theta_1)
  expect_error(gene_optima(data.frame(gene = "a")), "theta_0")
})

# ---- sim_count_params() ----

test_that("sim_count_params draws the dispersion and the size factors", {
  tree <- make_tree3()
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 15)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  set.seed(2); cp <- sim_count_params(expr, tree, genes)
  set.seed(2); expect_identical(sim_count_params(expr, tree, genes), cp)
  tips <- tree@data$node[tree@data$status == 1]
  expect_equal(names(cp$size_factor), as.character(tips))
  expect_equal(mean(cp$size_factor), 1)
  expect_length(cp$phi, 15)
  expect_true(all(cp$phi >= 1e-8))
  expect_identical(cp$expr, expr)            # no calibration: nothing is shifted
  expect_equal(cp$shift, 0)
  expect_equal(cp$gene_params$phi, cp$phi)
  # the draws feed the counts simulation
  counts <- sim_expression_counts_matrix(cp$expr, tree, phi = cp$phi, size_factor = cp$size_factor)
  expect_equal(nrow(counts), length(tips))
})

test_that("sim_count_params calibrates the mean count exactly", {
  tree <- make_tree3()
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 15)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  set.seed(2)
  cp <- sim_count_params(expr, tree, genes, target_mean_count = 2.5)
  tips <- tree@data$node[tree@data$status == 1]
  x_tip <- as.matrix(cp$expr[match(tips, cp$expr$node), genes$gene])
  expect_equal(mean(exp(t(x_tip)) * matrix(cp$size_factor, 15, length(tips), byrow = TRUE)), 2.5)
  # the constant is added to the log-scale parameters too
  expect_equal(cp$gene_params$log_g, genes$log_g + cp$shift)
  expect_equal(cp$gene_params$theta_1, genes$theta_1 + cp$shift)
  expect_equal(cp$gene_params$alpha, genes$alpha)
})

test_that("sim_count_params can make the dispersion decrease with the level of the gene", {
  tree <- make_tree3()
  set.seed(1)
  genes <- sim_gene_params(tree, n_genes = 200, log_g_sd = 2)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  cp <- sim_count_params(expr, tree, genes, phi_mode = "trend", phi_sd = 0.1)
  expect_lt(stats::cor(genes$log_g, log(cp$phi), method = "spearman"), -0.8)
})

test_that("sim_count_params validates its input", {
  tree <- make_tree3()
  genes <- sim_gene_params(tree, n_genes = 3)
  expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2, root_value = genes$x_root,
                                theta = gene_optima(genes), gene_names = genes$gene)
  expect_error(sim_count_params(expr, tree, genes, phi = 1), "gene_prior")
  expect_error(sim_count_params(expr, tree, genes, cell_size_sd = -1), "cell_size_sd")
  expect_error(sim_count_params(expr, tree, genes, target_mean_count = 0), "target_mean_count")
  expect_error(sim_count_params(expr[-3], tree, genes), "column for every gene")
  expect_error(sim_count_params(expr, tree, genes[names(genes) != "log_g"], phi_mode = "trend"), "log_g")
})
