test_that("output has node, type and one column per gene, one row per node", {
  tree <- make_tree()
  th <- data.frame(type = sort(unique(tree@data$type)),
                   g1 = 0, g2 = 1, g3 = 2)
  res <- sim_expression_matrix(tree, alpha = c(0.1, 0.5, 1), theta = th)

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("node", "type", "1", "2", "3"))
  expect_equal(nrow(res), nrow(tree@data))
  expect_setequal(res$node, tree@data$node)
  expect_false(anyNA(res[, c("1", "2", "3")]))
})

test_that("length-1 parameters are recycled across genes", {
  tree <- make_tree()
  res <- sim_expression_matrix(tree, alpha = c(0.1, 0.5), theta = theta_shared(tree))
  expect_equal(ncol(res), 2 + 2)
  expect_equal(colnames(res)[-(1:2)], c("1", "2"))
})

test_that("n_genes is inferred from the parameters, and is 1 if all have length 1", {
  tree <- make_tree()
  th <- theta_shared(tree)
  expect_equal(ncol(sim_expression_matrix(tree, theta = th)), 3)
  expect_equal(ncol(sim_expression_matrix(tree, sigmasq = c(0.1, 0.2, 0.3, 0.4),
                                          theta = th)), 2 + 4)
  expect_equal(ncol(sim_expression_matrix(tree, root_value = 1:5, theta = th)), 2 + 5)
})

test_that("explicit n_genes recycles length-1 parameters", {
  tree <- make_tree()
  res <- sim_expression_matrix(tree, theta = theta_shared(tree), n_genes = 4)
  expect_equal(colnames(res), c("node", "type", "1", "2", "3", "4"))
  # genes are simulated independently, so the columns differ
  expect_false(isTRUE(all.equal(res[["1"]], res[["2"]])))
})

test_that("n_genes must equal the parameter length unless all parameters have length 1", {
  tree <- make_tree()
  th <- theta_shared(tree)
  # a parameter of length 2 fixes N = 2: n_genes can only repeat it
  expect_error(sim_expression_matrix(tree, alpha = c(1, 2), theta = th, n_genes = 3),
               "`n_genes` \\(3\\) must equal the length of the parameters \\(2\\)")
  expect_error(sim_expression_matrix(tree, alpha = c(1, 2), theta = th, n_genes = 1),
               "`n_genes` \\(1\\) must equal")
  expect_equal(ncol(sim_expression_matrix(tree, alpha = c(1, 2), theta = th, n_genes = 2)), 4)
  # all length 1: n_genes sets N
  expect_equal(ncol(sim_expression_matrix(tree, theta = th, n_genes = 5)), 2 + 5)
  expect_equal(ncol(sim_expression_matrix(tree, theta = th, n_genes = 1)), 2 + 1)
})

test_that("incompatible parameter lengths are rejected, naming the parameter", {
  tree <- make_tree()
  th <- theta_shared(tree)
  expect_error(sim_expression_matrix(tree, alpha = c(1, 2), sigmasq = c(1, 2, 3),
                                     theta = th),
               "`alpha` has length 2")
  th3 <- data.frame(type = th$type, a = 0, b = 1, c = 2)
  # N is inferred as 3, so it is `alpha` that does not fit
  expect_error(sim_expression_matrix(tree, alpha = c(1, 2), theta = th3),
               "`alpha` has length 2")
  # with an explicit n_genes, the 3 theta columns fix N = 3
  expect_error(sim_expression_matrix(tree, theta = th3, n_genes = 2),
               "`n_genes` \\(2\\) must equal the length of the parameters \\(3\\)")
})

test_that("invalid n_genes is rejected", {
  tree <- make_tree()
  th <- theta_shared(tree)
  for (bad in list(0, -1, 1.5, NA, c(2, 3), "2")) {
    expect_error(sim_expression_matrix(tree, theta = th, n_genes = bad), "n_genes")
  }
})

test_that("invalid numeric parameters are rejected", {
  tree <- make_tree()
  th <- theta_shared(tree)
  expect_error(sim_expression_matrix(tree, alpha = "a", theta = th), "`alpha`")
  expect_error(sim_expression_matrix(tree, sigmasq = NA_real_, theta = th), "`sigmasq`")
  expect_error(sim_expression_matrix(tree, root_value = numeric(0), theta = th),
               "`root_value`")
})

test_that("theta must be a data frame with type and numeric gene columns", {
  tree <- make_tree()
  expect_error(sim_expression_matrix(tree), "`theta`")
  expect_error(sim_expression_matrix(tree, theta = c(0, 1)), "`theta`")
  expect_error(sim_expression_matrix(tree, theta = data.frame(type = 0L)), "`theta`")
  expect_error(sim_expression_matrix(tree, theta = data.frame(a = 0, b = 1)), "`theta`")
  expect_error(sim_expression_matrix(tree, theta = data.frame(type = 0:1, g = c("a", "b"))),
               "numeric")
})

test_that("a theta missing a type present in the tree returns NULL with a message", {
  tree <- make_tree()
  th <- data.frame(type = 0L, theta = 0)
  expect_message(res <- sim_expression_matrix(tree, theta = th), "missing an entry")
  expect_null(res)
})

test_that("each gene uses its own parameters and matches sim_expression_ou", {
  tree <- make_tree()
  types <- sort(unique(tree@data$type))
  th <- data.frame(type = types, g1 = 0, g2 = 5)

  # gene 2 has a very different optimum and no noise: it stays at its optimum
  # for every type
  res <- sim_expression_matrix(tree, alpha = 1, sigmasq = c(0.2, 1e-12),
                               root_value = c(0, 5), theta = th)
  expect_equal(res[["2"]], rep(5, nrow(res)), tolerance = 1e-4)

  # a single gene reproduces sim_expression_ou for the same seed
  set.seed(7)
  single <- sim_expression_ou(tree, alpha = 0.5, sigmasq = 0.2, root_value = 1,
                              theta = data.frame(type = types, theta = 2))
  set.seed(7)
  multi <- sim_expression_matrix(tree, alpha = 0.5, sigmasq = 0.2, root_value = 1,
                                 theta = data.frame(type = types, theta = 2))
  expect_equal(multi$node, single$node)
  expect_equal(multi$type, single$type)
  expect_equal(multi[["1"]], single$expr)
})

test_that("type-specific optima are applied per gene", {
  tree <- make_tree()
  types <- sort(unique(tree@data$type))
  skip_if(length(types) < 2, "tree has a single type")
  th <- data.frame(type = types, g1 = 0, g2 = seq_along(types) * 10)
  res <- sim_expression_matrix(tree, alpha = 5, sigmasq = 1e-6, root_value = 0, theta = th)
  # with strong mean reversion and tiny noise, long-lived nodes sit near the
  # optimum of their type for gene 2 (the order of magnitude is enough)
  expect_true(all(abs(res[["1"]]) < 1))
  expect_gt(mean(res[["2"]][res$type == max(types)]),
            mean(res[["2"]][res$type == min(types)]))
})

test_that("results are reproducible with a seed", {
  tree <- make_tree()
  th <- theta_shared(tree)
  set.seed(3); a <- sim_expression_matrix(tree, alpha = c(0.1, 0.5), theta = th)
  set.seed(3); b <- sim_expression_matrix(tree, alpha = c(0.1, 0.5), theta = th)
  expect_equal(a, b)
})

test_that("gene_names replaces the default gene ids", {
  tree <- make_tree()
  th <- theta_shared(tree)
  res <- sim_expression_matrix(tree, alpha = c(0.1, 0.5, 1), theta = th,
                               gene_names = c("Sox2", "Pax6", "Gata1"))
  expect_equal(colnames(res), c("node", "type", "Sox2", "Pax6", "Gata1"))

  # also works with n_genes and recycled parameters
  res <- sim_expression_matrix(tree, theta = th, n_genes = 2, gene_names = c("a", "b"))
  expect_equal(colnames(res), c("node", "type", "a", "b"))

  # names only relabel the columns: values are the same as with default names
  set.seed(5); named <- sim_expression_matrix(tree, theta = th, n_genes = 2,
                                              gene_names = c("a", "b"))
  set.seed(5); default <- sim_expression_matrix(tree, theta = th, n_genes = 2)
  expect_equal(unname(as.list(named[, c("a", "b")])), unname(as.list(default[, c("1", "2")])))
})

test_that("numeric gene_names are converted to column names", {
  tree <- make_tree()
  res <- sim_expression_matrix(tree, theta = theta_shared(tree), n_genes = 2,
                               gene_names = c(10, 20))
  expect_equal(colnames(res), c("node", "type", "10", "20"))
})

test_that("invalid gene_names are rejected", {
  tree <- make_tree()
  th <- theta_shared(tree)
  # a single name cannot be shared by several genes
  expect_error(sim_expression_matrix(tree, theta = th, n_genes = 2, gene_names = "a"),
               "single name cannot be shared by 2 genes")
  expect_error(sim_expression_matrix(tree, alpha = c(1, 2, 3), theta = th,
                                     gene_names = c("a", "b")),
               "number of genes \\(3\\).*`gene_names` has length 2")
  expect_error(sim_expression_matrix(tree, theta = th, n_genes = 2,
                                     gene_names = c("a", "a")), "unique")
  expect_error(sim_expression_matrix(tree, theta = th, n_genes = 2,
                                     gene_names = c("a", NA)), "missing")
  expect_error(sim_expression_matrix(tree, theta = th, n_genes = 2,
                                     gene_names = c("a", "")), "non-empty")
  expect_error(sim_expression_matrix(tree, theta = th, n_genes = 2,
                                     gene_names = c("a", "type")), "reserved")
})

test_that("gene_names of length N sets the number of genes like any other parameter", {
  tree <- make_tree()
  th <- theta_shared(tree)
  # all other parameters have length 1
  res <- sim_expression_matrix(tree, theta = th, gene_names = c("a", "b", "c"))
  expect_equal(colnames(res), c("node", "type", "a", "b", "c"))

  # agrees with the other parameters, and with n_genes
  res <- sim_expression_matrix(tree, alpha = c(0.1, 0.2), theta = th,
                               gene_names = c("a", "b"))
  expect_equal(colnames(res), c("node", "type", "a", "b"))
  res <- sim_expression_matrix(tree, theta = th, gene_names = c("a", "b"), n_genes = 2)
  expect_equal(colnames(res), c("node", "type", "a", "b"))

  # conflicts with n_genes when it sets N
  expect_error(sim_expression_matrix(tree, theta = th, gene_names = c("a", "b"),
                                     n_genes = 3),
               "`n_genes` \\(3\\) must equal the length of the parameters \\(2\\)")

  # a single name is valid for a single gene
  res <- sim_expression_matrix(tree, theta = th, gene_names = "Sox2")
  expect_equal(colnames(res), c("node", "type", "Sox2"))
})
