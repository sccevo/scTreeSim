test_that("sim_expression_counts returns named non-negative whole-number counts in the order of x", {
  x <- c(a = 1, b = 2, c = 0)
  set.seed(1)
  res <- sim_expression_counts(x, phi = 0.5)
  expect_named(res, c("a", "b", "c"))
  expect_true(all(res >= 0))
  expect_equal(unname(res), round(unname(res)))
})

test_that("sim_expression_counts has the expected mean and dispersion", {
  n <- 20000
  x <- stats::setNames(rep(log(10), n), paste0("c", seq_len(n)))
  set.seed(1)
  nb <- sim_expression_counts(x, phi = 0.5)
  expect_equal(mean(nb), 10, tolerance = 0.05)
  # variance = mu + phi * mu^2 = 60
  expect_equal(stats::var(nb), 60, tolerance = 0.1)

  # phi = 0 is a Poisson: variance = mean
  pois <- sim_expression_counts(x, phi = 0)
  expect_equal(mean(pois), 10, tolerance = 0.05)
  expect_equal(stats::var(pois), 10, tolerance = 0.1)
})

test_that("size factors scale the counts and are matched by name", {
  n <- 5000
  nm <- paste0("c", seq_len(n))
  x <- stats::setNames(rep(log(5), n), nm)
  sf <- stats::setNames(rep(c(1, 4), length.out = n), nm)
  set.seed(1)
  res <- sim_expression_counts(x, phi = 0.1, size_factor = sf)
  expect_equal(mean(res[sf == 1]), 5, tolerance = 0.1)
  expect_equal(mean(res[sf == 4]), 20, tolerance = 0.1)

  # the order of size_factor does not matter
  set.seed(2); a <- sim_expression_counts(x[1:10], phi = 0.1, size_factor = sf[1:10])
  set.seed(2); b <- sim_expression_counts(x[1:10], phi = 0.1, size_factor = rev(sf[1:10]))
  expect_equal(a, b)
})

test_that("sim_expression_counts does not touch the random seed and is reproducible", {
  x <- c(a = 1, b = 2, c = 0)
  set.seed(3); a <- sim_expression_counts(x, phi = 0.5)
  set.seed(3); b <- sim_expression_counts(x, phi = 0.5)
  expect_equal(a, b)
  expect_false("seed" %in% names(formals(sim_expression_counts)))
})

test_that("sim_expression_counts validates its input", {
  x <- c(a = 1, b = 2)
  expect_error(sim_expression_counts(c(1, 2), phi = 0.1), "named")
  expect_error(sim_expression_counts(c(a = 1, a = 2), phi = 0.1), "unique")
  expect_error(sim_expression_counts(c(a = 1, b = NA), phi = 0.1), "`x`")
  expect_error(sim_expression_counts("a", phi = 0.1), "`x`")
  expect_error(sim_expression_counts(x, phi = -1), "`phi`")
  expect_error(sim_expression_counts(x, phi = c(1, 2, 3)), "`phi`")
  expect_error(sim_expression_counts(x, phi = NA_real_), "`phi`")
  expect_error(sim_expression_counts(x, phi = 0.1, size_factor = c(1, 1)), "named")
  expect_error(sim_expression_counts(x, phi = 0.1, size_factor = c(a = 1, z = 1)), "match")
  expect_error(sim_expression_counts(x, phi = 0.1, size_factor = c(a = 1, a = 1)), "match")
  expect_error(sim_expression_counts(x, phi = 0.1, size_factor = c(a = 1, b = 0)),
               "positive")
  # phi may be given per cell
  expect_length(sim_expression_counts(x, phi = c(0.1, 0.2)), 2)
})

# ---- sim_expression_counts_matrix ----

make_expr <- function(tree, ...) {
  sim_expression_matrix(tree, theta = theta_shared(tree), ...)
}

test_that("counts matrix has node, type and the gene columns, one row per tip", {
  tree <- make_tree()
  expr <- make_expr(tree, n_genes = 3, gene_names = c("a", "b", "c"))
  set.seed(1)
  res <- sim_expression_counts_matrix(expr, tree, phi = 0.3)

  tips <- tree@data$node[tree@data$status == 1]
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("node", "type", "a", "b", "c"))
  expect_setequal(res$node, tips)
  expect_equal(nrow(res), length(tips))
  # type is carried over from the tips
  expect_equal(res$type, expr$type[match(res$node, expr$node)])
  # counts are non-negative integers
  counts <- as.matrix(res[, c("a", "b", "c")])
  expect_true(all(counts >= 0))
  expect_true(all(counts == round(counts)))
})

test_that("default gene ids 1..N are kept as column names", {
  tree <- make_tree()
  set.seed(1)
  res <- sim_expression_counts_matrix(make_expr(tree, n_genes = 2), tree, phi = 0.3)
  expect_equal(colnames(res), c("node", "type", "1", "2"))
})

test_that("counts matrix matches sim_expression_counts called once per gene", {
  tree <- make_tree()
  expr <- make_expr(tree, n_genes = 2)
  tips <- expr[expr$node %in% tree@data$node[tree@data$status == 1], ]
  cells <- as.character(tips$node)

  set.seed(10)
  res <- sim_expression_counts_matrix(expr, tree, phi = c(0.1, 0.5))
  set.seed(10)
  g1 <- sim_expression_counts(stats::setNames(tips[["1"]], cells), phi = 0.1)
  g2 <- sim_expression_counts(stats::setNames(tips[["2"]], cells), phi = 0.5)
  expect_equal(res[["1"]], unname(g1))
  expect_equal(res[["2"]], unname(g2))
})

test_that("phi has length 1 or one value per gene", {
  tree <- make_tree()
  expr <- make_expr(tree, n_genes = 3)
  expect_silent(sim_expression_counts_matrix(expr, tree, phi = 0.2))
  expect_silent(sim_expression_counts_matrix(expr, tree, phi = c(0.1, 0.2, 0.3)))
  expect_error(sim_expression_counts_matrix(expr, tree, phi = c(0.1, 0.2)), "`phi`")
  expect_error(sim_expression_counts_matrix(expr, tree, phi = -1), "`phi`")
  expect_error(sim_expression_counts_matrix(expr, tree, phi = NA_real_), "`phi`")
})

test_that("size factors are shared by all genes and matched by node id", {
  tree <- make_tree()
  expr <- make_expr(tree, n_genes = 2, root_value = 3)
  tips <- tree@data$node[tree@data$status == 1]
  # one cell gets a tiny size factor: its counts must be ~0 for every gene
  sf <- stats::setNames(rep(1, length(tips)), as.character(tips))
  sf[1] <- 1e-9
  set.seed(1)
  res <- sim_expression_counts_matrix(expr, tree, phi = 0.1, size_factor = sf)
  low <- res$node == tips[1]
  expect_equal(unlist(res[low, c("1", "2")], use.names = FALSE), c(0, 0))

  # order of size_factor is free
  set.seed(2); a <- sim_expression_counts_matrix(expr, tree, phi = 0.1, size_factor = sf)
  set.seed(2); b <- sim_expression_counts_matrix(expr, tree, phi = 0.1, size_factor = rev(sf))
  expect_equal(a, b)

  expect_error(sim_expression_counts_matrix(expr, tree, phi = 0.1,
                                            size_factor = unname(sf)), "size_factor")
  expect_error(sim_expression_counts_matrix(expr, tree, phi = 0.1,
                                            size_factor = sf[-1]), "size_factor")
})

test_that("counts matrix is reproducible with a seed", {
  tree <- make_tree()
  expr <- make_expr(tree, n_genes = 2)
  set.seed(4); a <- sim_expression_counts_matrix(expr, tree, phi = 0.3)
  set.seed(4); b <- sim_expression_counts_matrix(expr, tree, phi = 0.3)
  expect_equal(a, b)
})

test_that("counts matrix validates expr_matrix and tree", {
  tree <- make_tree()
  expr <- make_expr(tree, n_genes = 2)
  expect_error(sim_expression_counts_matrix(as.matrix(expr), tree, phi = 0.1), "expr_matrix")
  expect_error(sim_expression_counts_matrix(expr[, c("node", "1")], tree, phi = 0.1),
               "expr_matrix")
  expect_error(sim_expression_counts_matrix(expr[, c("node", "type")], tree, phi = 0.1),
               "no gene columns")
  bad <- expr
  bad[["1"]] <- as.character(bad[["1"]])
  expect_error(sim_expression_counts_matrix(bad, tree, phi = 0.1), "numeric")
  expect_error(sim_expression_counts_matrix(expr, "tree", phi = 0.1))
  # an expression matrix that lacks some tips of the tree
  tips <- tree@data$node[tree@data$status == 1]
  expect_error(sim_expression_counts_matrix(expr[expr$node != tips[1], ], tree, phi = 0.1),
               "no row for some tips")
})
