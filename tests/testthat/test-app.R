# Tests of the Shiny app (inst/app/app.R) with shiny::testServer(): the settings are set as inputs, the buttons are
# "pressed", and the results are compared with those of the script that the app writes.

app_file <- function() system.file("app", "app.R", package = "scTreeSim")

# the functions defined at the top of app.R (the code of the layers, the draws of the parameters, the table reader)
app_functions <- function() {
  env <- new.env(parent = globalenv())
  sys.source(app_file(), env)
  env
}

skip_if_no_app <- function() {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("ggtree")
  skip_if_not_installed("ggplot2")
  skip_if(app_file() == "", "the app is not installed")
}

# the settings of the six tabs, small enough to be fast
set_tree <- function(session, stop = "ntaxa", n_types = 3, normalize = TRUE) {
  session$setInputs(tree_stop = stop, n_tips = 25, origin_time = 3, tree_seed = 1, n_types = n_types, rho = 1,
                    keep_single = TRUE, normalize_height = normalize, tree_height = 1, sim_tree = 1)
}
set_layers <- function(session) {
  session$setInputs("snv_clock-type" = "strict", "snv_clock-rate" = 0.15, "snv_clock-spike" = FALSE,
                    snv_seed = 2, snv_n_sites = 10, snv_origin = "homozygous",
                    snv_r_AC = 1, snv_r_AG = 1, snv_r_AT = 1, snv_r_CG = 1, snv_r_CT = 1, snv_r_GT = 1,
                    snv_pi_A = 0.25, snv_pi_C = 0.25, snv_pi_G = 0.25, snv_pi_T = 0.25,
                    snv_epsilon = 0.1, snv_delta = 0.1)
  session$setInputs("bin_clock-type" = "lognormal", "bin_clock-rate" = 0.15, "bin_clock-sd" = 1, "bin_clock-spike" = TRUE,
                    "bin_clock-spike_mean" = 0.05, "bin_clock-spike_fixed" = FALSE, "bin_clock-spike_shape" = 2,
                    bin_seed = 3, bin_n_sites = 15, bin_n_sequences = 1, bin_lambda = 1, bin_origin = "random",
                    bin_alpha = 0.1, bin_beta = 0.1)
  session$setInputs(bc_seed = 4, bc_model = "generic", bc_probs_source = "fixed", bc_edit_probs = "0.5, 0.3, 0.2", bc_edit_rate = 1,
                    bc_silencing_rate = 0.01, bc_dropout_prob = 0.05, bc_n_barcodes = 2, bc_n_sites = 4,
                    bc_sequential = TRUE, bc_profile = "decay", bc_half_life = 1, bc_by_type = TRUE)
  session$setInputs(expr_seed = 1, expr_source = "priors", gap_cutoff = 1, n_gene = 12, pi_ou = 0.6,
                    psi_concentration = 0.5, mean_count_gene = 1, log_g_sd = 1.5,
                    "alpha_prior-type" = "Exponential", "alpha_prior-p1" = 1, "alpha_prior-lower" = 0.1, "alpha_prior-upper" = NA,
                    "sigma2_prior-type" = "Exponential", "sigma2_prior-p1" = 1, "sigma2_prior-lower" = 1e-8, "sigma2_prior-upper" = 2,
                    "eta_prior-type" = "Exponential", "eta_prior-p1" = 1, "eta_prior-lower" = 0, "eta_prior-upper" = NA,
                    "phi_prior-type" = "Lognormal", "phi_prior-p1" = -2, "phi_prior-p2" = 2, "phi_prior-lower" = 1e-8, "phi_prior-upper" = NA)
  session$setInputs(counts_seed = 5, cell_size_sd = 0.9, use_target = TRUE, target_mean_count = 1, phi_mode = "prior")
}

test_that("the script of the app reproduces what it showed, with any number of types", {
  skip_if_no_app()
  skip_on_cran()
  shiny::testServer(shiny::shinyAppFile(app_file()), {
    set_tree(session)
    tr <- tree_sim()
    expect_equal(sort(unique(tr@data$type)), 0:2)
    expect_equal(max(ape::node.depth.edgelength(tr@phylo)), 1)   # normalized
    set_layers(session)
    session$setInputs(sim_snv = 1, sim_bin = 1, sim_bc = 1, sim_expr = 1, sim_counts = 1)
    es <- expr_sim()
    expect_true(all(c("theta_0", "theta_1", "theta_2", "psi_0", "psi_1") %in% names(es$truth)))
    expect_equal(nrow(es$truth), 12)

    script <- script_text()
    env <- new.env()
    eval(parse(text = script), env)
    expect_equal(env$tree@phylo, tr@phylo)
    expect_identical(env$snv, snv_sim()$snv)
    expect_identical(env$bin, bin_sim()$bin)
    expect_identical(env$bin_true, bin_sim()$true)
    expect_identical(env$bc, bc_sim()$bc)
    expect_identical(env$gene_params, es$truth)
    expect_identical(env$expr, es$expr)
    expect_identical(env$count_params$expr, counts_sim()$expr)   # shifted by the level calibration
    expect_identical(env$count_params$gene_params, counts_sim()$truth)
    expect_identical(env$counts, counts_sim()$counts)
    # the seeds are in the script
    expect_match(script, "set.seed(1)", fixed = TRUE)
    expect_match(script, "set.seed(5)", fixed = TRUE)
    expect_match(script, "sim_adb_ntaxa_samp", fixed = TRUE)
    # the parameters of the genes are drawn in the script, not saved in it
    expect_match(script, "sim_gene_params(tree, n_genes = 12", fixed = TRUE)
    expect_match(script, 'alpha = sim_prior("exponential", rate = 1, lower = 0.1)', fixed = TRUE)
    expect_match(script, 'sigma2 = sim_prior("exponential", rate = 1, lower = 1e-08, upper = 2)', fixed = TRUE)
    expect_match(script, "sim_count_params(expr, tree, gene_params", fixed = TRUE)
    expect_lt(nchar(script), 6000)

    # the true parameters are tables with a row per gene / cell
    expect_equal(nrow(counts_sim()$truth), 12)
    expect_true("phi" %in% names(counts_sim()$truth))
    expect_equal(nrow(cell_sizes_table(counts_sim())), 25)
  })
})

test_that("the app simulates a tree at a fixed time and reads the parameters of the genes from a table", {
  skip_if_no_app()
  skip_on_cran()
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(gene = c("a", "b", "c"), alpha = c(0, 1.5, 2), sigma2 = c(0.5, 0.4, 0.3),
                              theta_0 = c(NA, 0, -1), theta_1 = c(NA, 0, 1)), path, row.names = FALSE)
  shiny::testServer(shiny::shinyAppFile(app_file()), {
    set_tree(session, stop = "time", n_types = 2, normalize = FALSE)
    session$setInputs(tree_seed = 2, origin_time = 4, sim_tree = 2)
    tr <- tree_sim()
    expect_match(tree_run()$code, "sim_adb_origin_samp(\n  origin_time = 4", fixed = TRUE)
    set_layers(session)
    session$setInputs(bc_probs_source = "prior", bc_n_outcomes = 4, bc_probs_concentration = 0.5, sim_bc = 1)
    bc <- bc_sim()
    expect_length(bc$edit_probs, 4)
    expect_equal(sum(bc$edit_probs), 1)
    expect_match(bc$code, "sim_dirichlet(4, concentration = 0.5)", fixed = TRUE)
    session$setInputs(expr_source = "table", expr_table = list(datapath = path, name = "genes.csv"), sim_expr = 1)
    es <- expr_sim()
    expect_equal(es$truth$gene, c("a", "b", "c"))
    expect_equal(es$truth$class_eff[1], "BM")
    expect_equal(es$truth$theta_1, c(0, 0, 1))
    expect_match(es$code, "genes.csv", fixed = TRUE)
    expect_match(es$code, "gene_params_from_table", fixed = TRUE)
    env <- new.env()
    eval(parse(text = script_text()), env)
    expect_identical(env$expr, es$expr)
    expect_identical(env$bc, bc_sim()$bc)
    expect_identical(env$edit_probs, bc$edit_probs)
  })
})

test_that("nothing is simulated until the user asks for a tree, and cells x genes is limited", {
  skip_if_no_app()
  shiny::testServer(shiny::shinyAppFile(app_file()), {
    expect_error(tree_run(), class = "shiny.silent.error")   # no tree at the start
    set_tree(session, n_types = 2)
    expect_s4_class(tree_sim(), "treedata")
    set_layers(session)
    session$setInputs(n_gene = 1e6, sim_expr = 1)
    expect_error(expr_sim(), class = "shiny.silent.error")   # 25 cells x 1,000,000 genes is too much
    # more tips than the limit
    session$setInputs(tree_stop = "ntaxa", n_tips = 20000, sim_tree = 2)
    expect_error(tree_run(), class = "shiny.silent.error")
  })
})

test_that("a new tree sends the layers back to the state in which nothing has been simulated", {
  skip_if_no_app()
  skip_on_cran()
  shiny::testServer(shiny::shinyAppFile(app_file()), {
    set_tree(session, n_types = 2)
    set_layers(session)
    session$setInputs(sim_snv = 1, sim_expr = 1)
    expect_true(layer_ready("sim_snv"))
    expect_true(layer_ready("sim_expr"))
    session$setInputs(tree_seed = 7, sim_tree = 2)
    expect_false(layer_ready("sim_snv"))
    expect_false(layer_ready("sim_expr"))
    session$setInputs(sim_snv = 2)   # pressing the button again simulates on the new tree
    expect_true(layer_ready("sim_snv"))
  })
})

test_that("the template of the uploads is a valid table, and files are read as data frames", {
  skip_if_no_app()
  app <- app_functions()
  tree <- make_tree3()
  path <- tempfile(fileext = ".csv")
  utils::write.csv(app$gene_table_template(3), path, row.names = FALSE)
  genes <- gene_params_from_table(app$read_gene_csv(path), tree, n_types = 3)
  expect_equal(genes$gene, c("gene_BM", "gene_OU1", "gene_OUM"))
  expect_equal(genes$class_eff, c("BM", "OU1", "OUM"))
  writeLines("not,a\n'table", path)
  expect_error(gene_params_from_table(app$read_gene_csv(path), tree, n_types = 3), "alpha and sigma2")
  # the code of a table reads back to the same table
  tab <- data.frame(gene = c("a", "b"), alpha = c(0, 1.5), sigma2 = c(0.5, 0.25), theta_0 = c(NA, 1 / 3))
  expect_identical(eval(parse(text = app$fmt_df(tab))), tab)
})

test_that("numbers are written as R code that reads back to the same value", {
  skip_if_no_app()
  app <- app_functions()
  x <- c(0.1, 1 / 3, 2 / 3 * pi, 1e-12, 123456789.123, -0.5, 5, Inf)
  expect_identical(eval(parse(text = app$fmt_vec(x))), x)
  m <- matrix(c(0.8, 0, 0.2, 1), 2, byrow = TRUE)
  expect_identical(eval(parse(text = app$fmt_matrix(m))), m)
  expect_identical(eval(parse(text = app$fmt_chr(c("a", "b\"c")))), c("a", "b\"c"))
})
