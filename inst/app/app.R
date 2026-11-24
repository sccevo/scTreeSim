# Shiny app of the scTreeSim package, started with scTreeSim::launch_gui():
#   1. Simulate tree         - age-dependent branching process with any number of types, stopped at a number of tips
#                              (sim_adb_ntaxa_samp()) or at a time (sim_adb_origin_samp()).
#   2. Simulate mutations    - diploid SNVs under the GT16 model (sim_gt16_nt_diploid()).
#   3. Simulate methylation  - binary sites (sim_binary()).
#   4. Simulate barcodes     - lineage barcodes (sim_barcode_seq(), sim_barcode_nonseq(), sim_barcode_generic()).
#   5. Simulate expression   - BM / OU1 / OUM process per gene (sim_expression_matrix()), with the parameters of the
#                              genes drawn from priors or read from a table.
#   6. Simulate counts       - negative-binomial counts per cell and gene (sim_expression_counts_matrix()).
#
# The simulations are done by the functions of the package and the figures by its plot_*() functions
# (R/plot_utils.R). This file contains the interface, the draws of the per-gene parameters from their priors, and the
# glue between the package functions.
#
# Every simulation runs from R code that is built from the settings at the moment the button is pressed, and the code is
# kept with the result. "Reproduce in R" writes the code of what is shown into a script, with the seeds.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(scTreeSim)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Limits for many cells. The plots with a row per cell draw a random subsample above MAX_CELLS_PLOT cells (the default of
# the plot functions of the package), the tree is drawn circular above MAX_TIPS_RECTANGULAR tips and without labels above
# MAX_TIP_LABELS tips, and the number of cells times the number of genes is limited by the memory.
MAX_TIPS <- 10000
MAX_CELLS_PLOT <- 500
MAX_TIPS_RECTANGULAR <- 1000
MAX_TIP_LABELS <- 300
MAX_CELLS_X_GENES <- 2e7

# ======================================================================================================
# R code of the simulations
# ======================================================================================================

# numbers as R code that reads back to exactly the same value
fmt_num <- function(x) {
  s <- sprintf("%.15g", x)
  bad <- is.finite(x) & suppressWarnings(as.numeric(s)) != x
  s[bad] <- sprintf("%.17g", x[bad])
  s[is.na(x)] <- "NA"
  s[!is.na(x) & x == Inf] <- "Inf"
  s[!is.na(x) & x == -Inf] <- "-Inf"
  s
}

# a vector as R code, in several lines when it is long
fmt_lines <- function(s, per_line) {
  if (length(s) == 1) return(s)
  rows <- split(s, ceiling(seq_along(s) / per_line))
  paste0("c(", paste(vapply(rows, paste, "", collapse = ", "), collapse = ",\n    "), ")")
}
fmt_vec <- function(x, per_line = 8) fmt_lines(fmt_num(x), per_line)
fmt_chr <- function(x, per_line = 5) fmt_lines(vapply(x, function(e) paste(deparse(e), collapse = ""), ""), per_line)
fmt_bool <- function(x) if (isTRUE(x)) "TRUE" else "FALSE"
fmt_matrix <- function(m) {
  rows <- vapply(seq_len(nrow(m)), function(i) paste(fmt_num(m[i, ]), collapse = ", "), "")
  sprintf("matrix(c(%s),\n    nrow = %d, byrow = TRUE)", paste(rows, collapse = ",\n             "), nrow(m))
}

# run the code of a layer, with the objects of the previous layers, and return its environment
eval_code <- function(code, vars = list()) {
  env <- list2env(vars, parent = globalenv())
  eval(parse(text = code), envir = env)
  env
}

# p: tree settings (see tree_run() in the server)
tree_code <- function(p) {
  fixed_time <- p$stop == "time"
  args <- c(if (fixed_time) sprintf("origin_time = %s", fmt_num(p$origin_time)) else sprintf("ntaxa = %s", fmt_num(p$ntaxa)),
            sprintf("scale = %s", fmt_vec(p$scale)), sprintf("shape = %s", fmt_vec(p$shape)),
            sprintf("death_prob = %s", fmt_vec(p$death)), sprintf("sampling_prob = %s", fmt_num(p$rho)),
            sprintf("origin_type = %s", fmt_num(p$origin_type)), sprintf("collapse = %s", fmt_bool(p$collapse)),
            sprintf("asym_trans_prob = %s", fmt_matrix(p$asym)), sprintf("sym_trans_prob = %s", fmt_matrix(p$sym)))
  code <- c(sprintf("set.seed(%s)", fmt_num(p$seed)),
            sprintf("tree <- %s(\n  %s)", if (fixed_time) "sim_adb_origin_samp" else "sim_adb_ntaxa_samp",
                    paste(args, collapse = ",\n  ")),
            "stopifnot(!is.null(tree))  # a simulation with too few tips returns NULL: try another seed")
  if (isTRUE(p$normalize)) {
    code <- c(code, sprintf("# normalize the height of the tree (time from the root to the tips) to %s", fmt_num(p$height)),
              sprintf("scale_factor <- %s / max(ape::node.depth.edgelength(tree@phylo))", fmt_num(p$height)),
              "tree@phylo$edge.length <- tree@phylo$edge.length * scale_factor",
              "tree@phylo$root.edge <- tree@phylo$root.edge * scale_factor",
              "tree@phylo$origin <- tree@phylo$origin * scale_factor")
  }
  paste(c(code, "# plot_sim_tree(tree)"), collapse = "\n")
}

snv_code <- function(p) {
  origin <- if (p$origin == "homozygous") {
    c(sprintf('x <- paste(sample(c("A", "C", "G", "T"), %s, replace = TRUE, prob = nuc_freq), collapse = "")', fmt_num(p$n_sites)),
      "origin_seq <- c(x, x)  # homozygous origin: the same random sequence on both alleles")
  } else "origin_seq <- NULL  # the origin is drawn from the equilibrium frequencies"
  paste(c(sprintf("set.seed(%s)", fmt_num(p$seed)),
          sprintf("clock <- %s", p$clock),
          "ct <- clock_tree(tree, clock)  # one draw of the clock, shared by all sites",
          sprintf("nuc_freq <- %s", fmt_vec(p$pi_raw)), "nuc_freq <- nuc_freq / sum(nuc_freq)",
          origin,
          sprintf("snv <- sim_gt16_nt_diploid(ct, n_sites = %s, rates = %s, pi = nuc_freq,\n                            origin_seq = origin_seq, clock = NULL, epsilon = %s, delta = %s)",
                  fmt_num(p$n_sites), fmt_vec(p$rates), fmt_num(p$epsilon), fmt_num(p$delta)),
          "# plot_sim_snv(ct, snv)"), collapse = "\n")
}

bin_code <- function(p) {
  call <- function(var, alpha, beta) {
    c(sprintf("set.seed(%s)  # the error is drawn after the sequence: the same seed gives the same true sequence", fmt_num(p$seed + 1)),
      sprintf("%s <- sim_binary(ct, n_sites = %s, n_sequences = %s, lambda = %s,\n                  origin_seq = origin_seq, clock = NULL, alpha = %s, beta = %s)",
              var, fmt_num(p$n_sites), fmt_num(p$n_sequences), fmt_num(p$lambda), fmt_num(alpha), fmt_num(beta)))
  }
  origin <- switch(p$origin, random = "origin_seq <- NULL  # drawn from the equilibrium", zeros = sprintf("origin_seq <- rep(0L, %s)", fmt_num(p$n_sites)),
                   ones = sprintf("origin_seq <- rep(1L, %s)", fmt_num(p$n_sites)))
  true <- if ((p$alpha > 0 || p$beta > 0) && p$n_sequences == 1) c("# the same sequence without the error", call("bin_true", 0, 0))
  paste(c(sprintf("set.seed(%s)", fmt_num(p$seed)), sprintf("clock <- %s", p$clock), "ct <- clock_tree(tree, clock)", origin,
          call("bin", p$alpha, p$beta), true, "# plot_sim_binary(ct, bin)"), collapse = "\n")
}

bc_code <- function(p) {
  common <- sprintf("edit_probs = edit_probs, silencing_rate = %s, dropout_prob = %s", fmt_num(p$silencing_rate), fmt_num(p$dropout_prob))
  call <- switch(p$model,
    seq = sprintf("sim_barcode_seq(tree, n_barcodes = %s, n_sites = %s,\n                       edit_rate = %s, %s)",
                  fmt_num(p$n_barcodes), fmt_num(p$n_sites), fmt_num(p$edit_rate), common),
    nonseq = sprintf("sim_barcode_nonseq(tree, n_targets = %s,\n                          edit_rate = %s, %s%s)",
                     fmt_num(p$n_targets), fmt_num(p$edit_rate), common,
                     if (isTRUE(p$window)) sprintf(",\n                          edit_height = %s, edit_duration = %s", fmt_num(p$edit_height), fmt_num(p$edit_duration)) else ""),
    generic = sprintf("sim_barcode_generic(tree, n_barcodes = %s, n_sites = %s,\n                          edit_rate = %s, %s,\n                          sequential = %s)",
                      fmt_num(p$n_barcodes), fmt_num(p$n_sites),
                      sprintf("barcode_rate(%s%s%s)", fmt_num(p$edit_rate),
                              if (!is.null(p$profile)) paste0(", time = ", p$profile) else "",
                              if (!is.null(p$by_type)) paste0(", by_type = ", fmt_vec(p$by_type)) else ""),
                      common, fmt_bool(p$sequential)))
  # the frequencies of the edit outcomes: fixed, or drawn from a Dirichlet distribution (with seed + 1)
  probs <- if (identical(p$probs_source, "prior")) {
    c("# The frequencies of the edit outcomes are drawn from a Dirichlet distribution: the same seed gives the same ones.",
      sprintf("set.seed(%s)", fmt_num(p$seed + 1)),
      sprintf("edit_probs <- sim_dirichlet(%s, concentration = %s)", fmt_num(p$n_outcomes), fmt_num(p$concentration)),
      sprintf("set.seed(%s)", fmt_num(p$seed)))
  } else {
    c(sprintf("set.seed(%s)", fmt_num(p$seed)), sprintf("edit_probs <- %s", fmt_vec(p$edit_probs)))
  }
  paste(c(probs, sprintf("bc <- %s", call), "# plot_sim_barcodes(tree, bc, n_edits = length(edit_probs))"), collapse = "\n")
}

# a prior (sim_prior()) as R code
prior_code <- function(p) {
  par <- paste(names(p$params), "=", vapply(p$params, fmt_num, ""), collapse = ", ")
  bounds <- paste0(if (is.finite(p$lower)) paste0(", lower = ", fmt_num(p$lower)),
                   if (is.finite(p$upper)) paste0(", upper = ", fmt_num(p$upper)))
  sprintf('sim_prior("%s", %s%s)', p$family, par, bounds)
}

# a data frame of numbers and strings as R code
fmt_df <- function(df) {
  cols <- lapply(df, function(x) if (is.character(x)) fmt_chr(x) else fmt_vec(as.numeric(x)))
  paste0("data.frame(\n", paste0("  ", names(cols), " = ", unlist(cols), collapse = ",\n"), ",\n  stringsAsFactors = FALSE)")
}

# p: source ("priors" or "table"), seed, n_types, gap_cutoff and, for the priors, n_genes, pi_ou, priors (alpha, sigma2, eta),
# psi_concentration, mean_count_gene, log_g_sd; for a table, table (as read) and table_name
expr_code <- function(p) {
  genes <- if (p$source == "priors") {
    c("# The parameters of the genes are drawn from their priors: the same seed gives the same genes.",
      "# To change some of them, save the table (write.csv(gene_params, ...)), edit it and read it back with",
      "# gene_params_from_table().",
      sprintf("set.seed(%s)", fmt_num(p$seed + 1)),
      sprintf("gene_params <- sim_gene_params(tree, n_genes = %s, n_types = %s, pi_ou = %s,\n                               alpha = %s,\n                               sigma2 = %s,\n                               eta = %s,\n                               psi_concentration = %s, mean_count_gene = %s, log_g_sd = %s,\n                               gap_cutoff = %s)",
              fmt_num(p$n_genes), fmt_num(p$n_types), fmt_num(p$pi_ou), prior_code(p$priors$alpha), prior_code(p$priors$sigma2),
              prior_code(p$priors$eta), fmt_num(p$psi_concentration), fmt_num(p$mean_count_gene), fmt_num(p$log_g_sd),
              fmt_num(p$gap_cutoff)))
  } else {
    c(sprintf("# The parameters of the genes come from the table that was uploaded in the app%s;", if (is.null(p$table_name)) "" else paste0(" (", p$table_name, ")")),
      '# to use the file itself, replace it by gene_table <- read.csv("<the file>").',
      paste0("gene_table <- ", fmt_df(p$table)),
      sprintf("gene_params <- gene_params_from_table(gene_table, tree, n_types = %s, gap_cutoff = %s)",
              fmt_num(p$n_types), fmt_num(p$gap_cutoff)))
  }
  paste(c(genes,
          sprintf("set.seed(%s)", fmt_num(p$seed)),
          "expr <- sim_expression_matrix(tree, alpha = gene_params$alpha, sigmasq = gene_params$sigma2,\n                               root_value = gene_params$x_root, theta = gene_optima(gene_params),\n                               gene_names = gene_params$gene)",
          "# plot_sim_expression(expr, tree, gene_params)",
          "# plot_sim_umap(expr, tree)"), collapse = "\n")
}

# p: seed, cell_size_sd, target_mean_count (or NULL), phi_mode, phi (a sim_prior, in "prior" mode), phi_a0, phi_a1, phi_sd
counts_code <- function(p) {
  args <- c(sprintf("cell_size_sd = %s", fmt_num(p$cell_size_sd)),
            if (!is.null(p$target_mean_count)) sprintf("target_mean_count = %s", fmt_num(p$target_mean_count)),
            if (p$phi_mode == "prior") sprintf("phi = %s", prior_code(p$phi))
            else sprintf('phi_mode = "trend", phi_a0 = %s, phi_a1 = %s, phi_sd = %s', fmt_num(p$phi_a0), fmt_num(p$phi_a1), fmt_num(p$phi_sd)))
  paste(c("# The dispersion of the genes and the size factors of the cells are drawn: the same seed gives the same values.",
          sprintf("set.seed(%s)", fmt_num(p$seed + 1)),
          sprintf("count_params <- sim_count_params(expr, tree, gene_params,\n                                 %s)", paste(args, collapse = ",\n                                 ")),
          "# count_params$gene_params are the true parameters of the genes, with phi (and shifted by the level calibration)",
          sprintf("set.seed(%s)", fmt_num(p$seed)),
          "counts_df <- sim_expression_counts_matrix(count_params$expr, tree, phi = count_params$phi,\n                                           size_factor = count_params$size_factor)",
          "counts <- t(as.matrix(counts_df[gene_params$gene]))  # genes x cells",
          "dimnames(counts) <- list(gene_params$gene, tree@phylo$tip.label[counts_df$node])"), collapse = "\n")
}

# ======================================================================================================
# Simulation of many genes
# ======================================================================================================

# The parameters of the genes are drawn by sim_gene_params() (from priors) or read by gene_params_from_table(), and the
# parameters of the counts by sim_count_params(): they are functions of the package, and the script of "Reproduce in R"
# calls them. The app builds the code of every simulation (see the top of the file) and runs it.
#
# 1. simulate_expression(): an independent BM / OU1 / OUM process per gene on the tree (sim_expression_matrix()).
#    The seed of the draws of the parameters is seed + 1; the evolution on the tree uses `seed`.
# 2. simulate_counts(): negative-binomial counts at the tips (the cells), from the result of step 1
#    (sim_expression_counts_matrix()). The seed of the draws of phi and of the size factors is seed + 1; the counts use `seed`.

# A table with the columns of the uploads, as an example (genes: a BM, an OU1 and an OUM process)
gene_table_template <- function(n_regime) {
  theta <- matrix(0, 3, n_regime)
  if (n_regime > 1) theta[3, ] <- seq(-1, 1, length.out = n_regime)
  tab <- data.frame(gene = c("gene_BM", "gene_OU1", "gene_OUM"), alpha = c(0, 1, 1), sigma2 = c(0.5, 0.5, 0.5))
  for (k in seq_len(n_regime)) tab[[paste0("theta_", k - 1L)]] <- theta[, k]
  tab
}

# the uploaded file, as a data frame
read_gene_csv <- function(path) {
  tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE),
           error = function(e) stop("The file could not be read as a CSV table: ", conditionMessage(e), call. = FALSE))
}

# tr: treedata object; n_regime: number of types of the tree.
# source: "priors" (the other arguments are those of sim_gene_params()) or "table" (`table`: the uploaded data frame).
# Returns list(tree, expr, truth, code), with the true parameters of the genes.
simulate_expression <- function(tr, n_regime, source = c("priors", "table"), table = NULL, table_name = NULL, seed = 1,
                                n_genes = 100, pi_ou = 0.5, gap_cutoff = 1,
                                mean_count_gene = 1, log_g_sd = 1.5, priors = NULL, psi_concentration = 0.5) {
  source <- match.arg(source)
  if (!is.null(table)) {
    keep_cols <- intersect(c("gene", "alpha", "sigma2", "sigmasq", "root_value", paste0("theta_", seq_len(n_regime) - 1L)), names(table))
    table <- table[keep_cols]
  }
  code <- expr_code(list(source = source, seed = seed, n_types = n_regime, gap_cutoff = gap_cutoff, n_genes = n_genes,
                         pi_ou = pi_ou, priors = priors, psi_concentration = psi_concentration,
                         mean_count_gene = mean_count_gene, log_g_sd = log_g_sd, table = table, table_name = table_name))
  env <- eval_code(code, list(tree = tr))
  if (is.null(env$expr)) stop("The expression simulation failed.")
  list(tree = tr, expr = env$expr, truth = env$gene_params, code = code, n_regime = n_regime)
}

# es                  result of simulate_expression()
# cell_size_sd, target_mean_count, phi_mode, phi, phi_a0, phi_a1, phi_sd   see sim_count_params()
# seed                the counts use `seed`; phi and the size factors are drawn with seed + 1
simulate_counts <- function(es, cell_size_sd = 0.9, target_mean_count = NULL,
                            phi = sim_prior("lognormal", meanlog = -2, sdlog = 2, lower = 1e-8),
                            phi_mode = c("prior", "trend"), phi_a0 = 0.1, phi_a1 = 0.5, phi_sd = 0.3, seed = 1) {
  phi_mode <- match.arg(phi_mode)
  tr <- es$tree
  code <- counts_code(list(seed = seed, cell_size_sd = cell_size_sd, target_mean_count = target_mean_count,
                           phi_mode = phi_mode, phi = phi, phi_a0 = phi_a0, phi_a1 = phi_a1, phi_sd = phi_sd))
  env <- eval_code(code, list(tree = tr, expr = es$expr, gene_params = es$truth))
  cp <- env$count_params                                  # the expression is shifted when the level is calibrated
  counts <- env$counts
  cell_names <- colnames(counts)

  # What would be observed with real data: totals per cell and gene, and the size factor derived from them
  cell_total    <- colSums(counts)
  gene_total    <- rowSums(counts)
  cell_size_obs <- cell_total / mean(cell_total)

  list(counts = counts, tree = tr, expr = cp$expr, truth = cp$gene_params,
       cell_size = setNames(unname(cp$size_factor[as.character(env$counts_df$node)]), cell_names),
       cell_size_obs = cell_size_obs,
       cell_total = cell_total, gene_total = gene_total, shift = cp$shift, code = code)
}

# ======================================================================================================
# Prior input module
# ======================================================================================================

# the families of sim_prior(), with the names shown in the app
PRIOR_FAMILIES <- c("Exponential" = "exponential", "Lognormal" = "lognormal", "Normal" = "normal", "Uniform" = "uniform")

# One prior: a family with its parameters and, optionally, a lower and an upper bound (empty: no bound). `lower` and
# `upper` are the defaults of the bounds (NA: no bound); `min_lower` is the smallest lower bound that makes sense
# (0 for parameters that cannot be negative).
priorUI <- function(id, label, selected = "Exponential", lower = NA, upper = NA, min_lower = NA) {
  ns <- NS(id)
  tagList(
    selectInput(ns("type"), label, choices = names(PRIOR_FAMILIES), selected = selected),
    uiOutput(ns("pars")),
    fluidRow(column(6, numericInput(ns("lower"), "Lower bound", value = lower, min = min_lower, step = 0.1)),
             column(6, numericInput(ns("upper"), "Upper bound", value = upper, step = 0.1)))
  )
}

# defaults: list with one numeric vector c(p1, p2) per prior family; min_lower as in priorUI().
# Returns a reactive with a sim_prior(), shared by all the parameters that are drawn from a prior
priorServer <- function(id, defaults, min_lower = -Inf) {
  moduleServer(id, function(input, output, session) {
    output$pars <- renderUI({
      ns <- session$ns
      d <- defaults[[input$type]]
      switch(input$type,
        Exponential = numericInput(ns("p1"), "Rate", value = d[1], min = 1e-6, step = 0.1),
        Lognormal   = fluidRow(column(6, numericInput(ns("p1"), "Mean (log scale)", value = d[1], step = 0.5)),
                               column(6, numericInput(ns("p2"), "SD (log scale)", value = d[2], min = 1e-6, step = 0.5))),
        Normal      = fluidRow(column(6, numericInput(ns("p1"), "Mean", value = d[1], step = 0.5)),
                               column(6, numericInput(ns("p2"), "SD", value = d[2], min = 1e-6, step = 0.5))),
        Uniform     = fluidRow(column(6, numericInput(ns("p1"), "Minimum", value = d[1], step = 0.5)),
                               column(6, numericInput(ns("p2"), "Maximum", value = d[2], step = 0.5)))
      )
    })
    # the parameter inputs live in a tab that may not have been opened yet: render them anyway so that the
    # simulation can read their defaults
    outputOptions(output, "pars", suspendWhenHidden = FALSE)
    reactive({
      req(input$type, input$p1)
      family <- PRIOR_FAMILIES[[input$type]]
      pars <- if (family == "exponential") list(input$p1) else { req(input$p2); list(input$p1, input$p2) }
      lower <- if (is.null(input$lower) || is.na(input$lower)) -Inf else input$lower
      upper <- if (is.null(input$upper) || is.na(input$upper)) Inf else input$upper
      tryCatch(do.call(sim_prior, c(list(family), pars, list(lower = max(lower, min_lower), upper = upper))),
               error = function(e) validate(need(FALSE, conditionMessage(e))))
    })
  })
}

# ======================================================================================================
# Clock module (mutations and methylation)
# ======================================================================================================

clockUI <- function(id, default_rate = 0.15) {
  ns <- NS(id)
  tagList(
    helpText("This tab draws its own clock, independent of the other tabs. The draw depends on the random seed of the tab: ",
             "two tabs with the same seed and the same clock settings get the same draw."),
    selectInput(ns("type"), "Clock model",
                choices = c("Strict" = "strict", "Relaxed, lognormal" = "lognormal",
                            "Relaxed, gamma" = "gamma", "Relaxed, exponential" = "exponential")),
    numericInput(ns("rate"), "Rate (expected substitutions per site per unit of time)", value = default_rate,
                 min = 0, step = 0.05),
    conditionalPanel("input.type == 'lognormal'", ns = ns,
                     numericInput(ns("sd"), "SD of the log multipliers of the rate", value = 1, min = 0, step = 0.1)),
    conditionalPanel("input.type == 'gamma'", ns = ns,
                     numericInput(ns("shape"), "Shape of the multipliers of the rate", value = 2, min = 0.01, step = 0.5)),
    checkboxInput(ns("spike"), "Division-triggered spikes", value = FALSE),
    conditionalPanel("input.spike", ns = ns,
      helpText("A spike is a burst of evolution at a division. It is only placed at divisions that are nodes of the tree: ",
               "keep 'Keep single-child nodes' ticked in the tree tab."),
      numericInput(ns("spike_mean"), "Mean spike (expected substitutions per site per division)", value = 0.05,
                   min = 0, step = 0.01),
      checkboxInput(ns("spike_fixed"), "Same size at every division (shape = Inf)", value = TRUE),
      conditionalPanel("!input.spike_fixed", ns = ns,
                       numericInput(ns("spike_shape"), "Shape of the gamma spikes (CV = 1 / sqrt(shape))", value = 2,
                                    min = 0.01, step = 0.5)))
  )
}

# returns a reactive with list(clock = the clock model (clock_strict() / clock_relaxed(), with or without spike()),
# code = the same call as R code)
clockServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactive({
      req(input$type, input$rate)
      spk <- NULL
      spk_code <- ""
      if (isTRUE(input$spike)) {
        req(input$spike_mean)
        shape <- if (isTRUE(input$spike_fixed)) Inf else { req(input$spike_shape); input$spike_shape }
        spk <- spike(mean = input$spike_mean, shape = shape)
        spk_code <- sprintf(", spike = spike(mean = %s, shape = %s)", fmt_num(input$spike_mean), fmt_num(shape))
      }
      rate <- fmt_num(input$rate)
      tryCatch(
        switch(input$type,
          strict      = list(clock = clock_strict(input$rate, spike = spk),
                             code = sprintf("clock_strict(%s%s)", rate, spk_code)),
          lognormal   = { req(input$sd)
                          list(clock = clock_relaxed(input$rate, "lognormal", sd = input$sd, spike = spk),
                               code = sprintf('clock_relaxed(%s, "lognormal", sd = %s%s)', rate, fmt_num(input$sd), spk_code)) },
          gamma       = { req(input$shape)
                          list(clock = clock_relaxed(input$rate, "gamma", shape = input$shape, spike = spk),
                               code = sprintf('clock_relaxed(%s, "gamma", shape = %s%s)', rate, fmt_num(input$shape), spk_code)) },
          exponential = list(clock = clock_relaxed(input$rate, "exponential", spike = spk),
                             code = sprintf('clock_relaxed(%s, "exponential"%s)', rate, spk_code))),
        error = function(e) validate(need(FALSE, conditionMessage(e))))
    })
  })
}

# probabilities typed as numbers separated by commas; NULL if they are not valid or do not sum to 1
parse_probs <- function(txt) {
  p <- suppressWarnings(as.numeric(strsplit(gsub("[[:space:]]", "", txt), ",", fixed = TRUE)[[1]]))
  if (length(p) == 0 || anyNA(p) || any(p < 0) || abs(sum(p) - 1) > 1e-6) NULL else p
}

# plot height (pixels) that grows with the number of cells that are drawn (at most MAX_CELLS_PLOT: above that the plots show
# a random subsample), for the heatmaps of the layers
layer_plot_height <- function(tree_fun, per_tip, min_height = 550) {
  function() {
    n_tips <- tryCatch(length(tree_fun()@phylo$tip.label), error = function(e) 0)
    max(min_height, per_tip * min(n_tips, MAX_CELLS_PLOT))
  }
}

# the notice above the plots that show a subsample of the cells
subsample_notice <- function(tree_fun) {
  renderUI({
    n_tips <- tryCatch(ape::Ntip(tree_fun()@phylo), error = function(e) 0)
    if (n_tips > MAX_CELLS_PLOT) {
      div(class = "alert alert-warning py-1 small",
          sprintf("The tree has %s cells: the plot shows a random subsample of %d of them, with the tree pruned to them. The simulation and the downloads use all the cells.",
                  format(n_tips, big.mark = ","), MAX_CELLS_PLOT))
    }
  })
}

# ======================================================================================================
# User interface
# ======================================================================================================

# default transition probabilities for n types: a type divides symmetrically (both daughters keep the type, 0.8) or
# asymmetrically (one daughter becomes the next type, 0.2); the last type only divides symmetrically
default_transitions <- function(n) {
  asym <- matrix(0, n, n)
  sym <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (i < n) { asym[i, i + 1] <- 0.2; sym[i, i] <- 0.8 } else sym[i, i] <- 1
  }
  list(asym = asym, sym = sym)
}

# a grid of inputs for a matrix of probabilities (rows: from, columns: to; types are numbered 0, 1, ...).
# `keep(id, default)` gives the value an input had before the grid was redrawn.
trans_matrix_ui <- function(prefix, title, default, keep) {
  n <- nrow(default)
  cell <- function(i, j) {
    id <- paste0(prefix, "_", i, j)
    numericInput(id, label = NULL, value = keep(id, default[i, j]), min = 0, max = 1, step = 0.05, width = "100%")
  }
  head_cells <- lapply(c("", paste("to", seq_len(n) - 1L)), function(x) div(class = "small text-muted text-center", x))
  rows <- unlist(lapply(seq_len(n), function(i) {
    c(list(div(class = "small text-muted", paste("from", i - 1L))), lapply(seq_len(n), function(j) cell(i, j)))
  }), recursive = FALSE)
  tagList(
    strong(title),
    div(style = sprintf("display: grid; grid-template-columns: 4em repeat(%d, 1fr); gap: 0 6px; align-items: center;", n),
        head_cells, rows)
  )
}

ui <- page_navbar(
  title = "scTreeSim",
  theme = bs_theme(version = 5),
  id = "main_tabs",
  header = tagList(
    tags$style(HTML(".shiny-output-error-validation { color: #b30000; font-weight: 500; }")),
    # The tabs and the buttons that need a tree are disabled until there is one (and the counts until there is
    # expression): the server sends the list of the ones to enable or disable. Everything but the tree starts disabled.
    tags$script(HTML("
      (function() {
        var treeTabs = ['mutations', 'methylation', 'barcodes', 'expression', 'counts'];
        var treeButtons = ['download_script', 'download_all'];
        var buttonsDisabled = true;
        function setDisabled(el, disabled) {
          if (!el) return;
          el.classList.toggle('disabled', disabled);
          if (disabled) { el.setAttribute('aria-disabled', 'true'); el.setAttribute('tabindex', '-1'); }
          else { el.removeAttribute('aria-disabled'); el.removeAttribute('tabindex'); }
        }
        function applyButtons() {
          treeButtons.forEach(function(id) { setDisabled(document.getElementById(id), buttonsDisabled); });
        }
        function apply(m) {
          [].concat(m.tabs || []).forEach(function(v) { setDisabled(document.querySelector('a.nav-link[data-value=\"' + v + '\"]'), m.disabled); });
          if (m.buttons) { buttonsDisabled = m.disabled; applyButtons(); }
        }
        document.addEventListener('DOMContentLoaded', function() { apply({tabs: treeTabs, buttons: treeButtons, disabled: true}); });
        // a disabled tab or button does nothing, whoever clicks it
        document.addEventListener('click', function(e) {
          var el = e.target.closest ? e.target.closest('a.nav-link.disabled, a.btn.disabled') : null;
          if (el) { e.preventDefault(); e.stopImmediatePropagation(); }
        }, true);
        // Shiny enables a download link again when it gives it its address: apply the state again then
        $(document).on('shiny:value', function(e) {
          if (treeButtons.indexOf(e.name) >= 0) setTimeout(applyButtons, 0);
        });
        if (window.Shiny && Shiny.addCustomMessageHandler) Shiny.addCustomMessageHandler('set_disabled', apply);
        else $(document).on('shiny:connected', function() { Shiny.addCustomMessageHandler('set_disabled', apply); });
      })();
    "))),

  # ---------------------------------------------------------------- Simulate tree
  nav_panel(
    "Simulate tree", value = "tree",
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        h5("scTreeSim parameters"),
        radioButtons("tree_stop", "The simulation stops at",
                     choices = c("A number of sampled tips: sim_adb_ntaxa_samp()" = "ntaxa",
                                 "A fixed time since the origin: sim_adb_origin_samp()" = "time")),
        conditionalPanel("input.tree_stop == 'ntaxa'",
          numericInput("n_tips", "Number of tips (cells), up to 10,000", value = 100, min = 10, max = 10000, step = 10),
          helpText("With many tips the plots of the tree and of the layers are drawn for a random subsample of the cells: the simulations themselves use all the cells.")),
        conditionalPanel("input.tree_stop == 'time'",
          numericInput("origin_time", "Time of the origin (from the first cell to the present)", value = 5, min = 0.1, step = 0.5),
          helpText("The number of tips is random: it grows quickly with the time and with the short lifetimes, so long times can give very large trees. A tree needs at least 2 sampled tips.")),
        numericInput("tree_seed", "Random seed (tree)", value = 1, step = 1),
        numericInput("n_types", "Number of types", value = 2, min = 1, max = 10, step = 1),
        helpText("The mean lifetime of a type is scale x shape, so the height of the tree follows from the lifetimes unless it is normalized below."),
        uiOutput("type_params"),
        sliderInput("rho", "Sampling probability (rho)", min = 0.05, max = 1, value = 1, step = 0.05),
        checkboxInput("keep_single", "Keep single-child nodes (collapse = FALSE)", value = TRUE),
        helpText("Keeps the divisions where a lineage died and the type changes that happen there. Without them a type change is hidden on a longer branch, which then takes the type of the node it ends in. This is needed for the type-dependent expression optima, the division-triggered spikes of the mutations and methylation, and the type-dependent barcode rates."),
        checkboxInput("normalize_height", "Normalize the tree height", value = FALSE),
        conditionalPanel("input.normalize_height",
          numericInput("tree_height", "Height of the tree (time from the root to the tips)", value = 1, min = 0.001, step = 0.5),
          helpText("All branch lengths, the root edge and the origin time are multiplied by the same factor. The clock rates of the other tabs are then per unit of this height, so their defaults give about the same amount of change in every tree.")),
        hr(),
        h6("Type transitions at division"),
        helpText("For each type i, the row i of the asymmetric matrix plus the row i of the symmetric matrix must sum to 1. Asymmetric: one daughter changes type; symmetric: both daughters get the new type."),
        uiOutput("trans_ui"),
        uiOutput("trans_check"),
        hr(),
        selectInput("tree_layout", "Layout", choices = c("auto", "rectangular", "slanted", "circular", "equal_angle")),
        helpText("auto: rectangular, and circular for trees with more than 1,000 tips."),
        checkboxInput("tree_tips", "Show tip labels (up to 300 tips)", value = FALSE),
        checkboxInput("tree_axis", "Show time axis (rectangular / slanted layouts)", value = TRUE),
        actionButton("sim_tree", "Simulate tree", class = "btn-primary")
      ),
      card(card_header("Simulated tree, coloured by type"),
           # the plot can be taller than the card (many tip labels): the card body scrolls
           card_body(fillable = FALSE, style = "overflow-y: auto;",
                     plotOutput("tree_plot", height = "auto", fill = FALSE),
                     verbatimTextOutput("tree_info")))
    )
  ),

  # ---------------------------------------------------------------- Simulate mutations
  nav_panel(
    "Simulate mutations", value = "mutations",
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        h5("Diploid SNVs (GT16 model)"),
        numericInput("snv_seed", "Random seed (mutations)", value = 2, step = 1),
        sliderInput("snv_n_sites", "Number of sites", min = 5, max = 200, value = 40, step = 5),
        radioButtons("snv_origin", "Genotype at the origin",
                     choices = c("Drawn from the equilibrium frequencies" = "random",
                                 "Homozygous: a random sequence on both alleles" = "homozygous")),
        h6("Substitution model"),
        helpText("Exchangeabilities of the six nucleotide pairs (a common scaling has no effect) and equilibrium nucleotide frequencies (rescaled to sum to 1)."),
        fluidRow(column(4, numericInput("snv_r_AC", "AC", value = 1, min = 0, step = 0.5)),
                 column(4, numericInput("snv_r_AG", "AG", value = 1, min = 0, step = 0.5)),
                 column(4, numericInput("snv_r_AT", "AT", value = 1, min = 0, step = 0.5))),
        fluidRow(column(4, numericInput("snv_r_CG", "CG", value = 1, min = 0, step = 0.5)),
                 column(4, numericInput("snv_r_CT", "CT", value = 1, min = 0, step = 0.5)),
                 column(4, numericInput("snv_r_GT", "GT", value = 1, min = 0, step = 0.5))),
        fluidRow(column(3, numericInput("snv_pi_A", "pi A", value = 0.25, min = 0, step = 0.05)),
                 column(3, numericInput("snv_pi_C", "pi C", value = 0.25, min = 0, step = 0.05)),
                 column(3, numericInput("snv_pi_G", "pi G", value = 0.25, min = 0, step = 0.05)),
                 column(3, numericInput("snv_pi_T", "pi T", value = 0.25, min = 0, step = 0.05))),
        hr(),
        h6("Clock"),
        clockUI("snv_clock"),
        hr(),
        h6("Sequencing error at the tips"),
        numericInput("snv_epsilon", "epsilon: amplification and sequencing error", value = 0, min = 0, max = 1, step = 0.05),
        numericInput("snv_delta", "delta: allelic dropout", value = 0, min = 0, max = 1, step = 0.05),
        hr(),
        h6("Display"),
        checkboxInput("snv_errors", "Show the observed genotypes and mark the errors", value = FALSE),
        helpText("The marks show where the error changed the genotype of a cell: they only appear when epsilon or delta is above 0."),
        checkboxInput("snv_genotypes", "Show the 16 genotypes in 16 colours (instead of two tiles of four nucleotides per site)", value = FALSE),
        checkboxInput("snv_distance", "Draw the distance tree (branch lengths in expected substitutions)", value = FALSE),
        actionButton("sim_snv", "Simulate mutations", class = "btn-primary"),
        downloadButton("download_snv", "Download (.rds)")
      ),
      card(card_header("Genotypes of the cells (one row per cell; two tiles per site: allele 1 left, allele 2 right)"),
           card_body(fillable = FALSE, style = "overflow-y: auto;",
                     uiOutput("snv_note"),
                     plotOutput("snv_plot", height = "auto", fill = FALSE),
                     verbatimTextOutput("snv_summary")))
    )
  ),

  # ---------------------------------------------------------------- Simulate methylation
  nav_panel(
    "Simulate methylation", value = "methylation",
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        h5("Binary sites"),
        numericInput("bin_seed", "Random seed (methylation)", value = 3, step = 1),
        sliderInput("bin_n_sites", "Number of sites per sequence", min = 10, max = 500, value = 100, step = 10),
        sliderInput("bin_n_sequences", "Number of independent sequences", min = 1, max = 5, value = 1, step = 1),
        numericInput("bin_lambda", "lambda: rate of 1 -> 0 relative to 0 -> 1", value = 1, min = 0.01, step = 0.1),
        radioButtons("bin_origin", "Sequence at the origin",
                     choices = c("Drawn from the equilibrium" = "random", "All 0" = "zeros", "All 1" = "ones")),
        hr(),
        h6("Clock"),
        clockUI("bin_clock"),
        hr(),
        h6("Error at the tips"),
        numericInput("bin_alpha", "alpha: false negative rate (1 read as 0)", value = 0, min = 0, max = 1, step = 0.05),
        numericInput("bin_beta", "beta: false positive rate (0 read as 1)", value = 0, min = 0, max = 1, step = 0.05),
        hr(),
        h6("Display"),
        checkboxInput("bin_errors", "Mark the errors (false positives and negatives; one sequence only)", value = FALSE),
        checkboxInput("bin_distance", "Draw the distance tree (branch lengths in expected substitutions)", value = FALSE),
        actionButton("sim_bin", "Simulate methylation", class = "btn-primary"),
        downloadButton("download_bin", "Download (.rds)")
      ),
      card(card_header("Binary sequences of the cells (black = 1, white = 0)"),
           card_body(fillable = FALSE, style = "overflow-y: auto;",
                     uiOutput("bin_note"),
                     plotOutput("bin_plot", height = "auto", fill = FALSE),
                     verbatimTextOutput("bin_summary")))
    )
  ),

  # ---------------------------------------------------------------- Simulate barcodes
  nav_panel(
    "Simulate barcodes", value = "barcodes",
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        h5("Lineage barcodes"),
        helpText("The barcodes evolve in the time of the tree: the rates below are per unit of time, so they are the clock of this tab. The relaxed clocks and spikes of the mutations and methylation tabs do not apply to barcodes."),
        numericInput("bc_seed", "Random seed (barcodes)", value = 4, step = 1),
        radioButtons("bc_model", "Simulator",
                     choices = c("Sequential editing, SciPhy-like: sim_barcode_seq()" = "seq",
                                 "Independent targets, TiDeTree-like: sim_barcode_nonseq()" = "nonseq",
                                 "Generic, with time and type dependent editing: sim_barcode_generic()" = "generic")),
        conditionalPanel("input.bc_model != 'nonseq'",
          sliderInput("bc_n_barcodes", "Number of barcodes per cell", min = 1, max = 5, value = 2, step = 1),
          sliderInput("bc_n_sites", "Number of target sites per barcode", min = 1, max = 30, value = 8, step = 1)),
        conditionalPanel("input.bc_model == 'nonseq'",
          sliderInput("bc_n_targets", "Number of independent target sites", min = 1, max = 60, value = 20, step = 1)),
        conditionalPanel("input.bc_model == 'generic'",
          checkboxInput("bc_sequential", "Sequential: sites are edited left to right; silencing and dropout affect the whole barcode", value = TRUE)),
        radioButtons("bc_probs_source", "Frequencies of the edit outcomes",
                     choices = c("Fixed values" = "fixed", "Drawn from a Dirichlet prior" = "prior")),
        conditionalPanel("input.bc_probs_source == 'fixed'",
          textInput("bc_edit_probs", "Frequencies (comma-separated, sum to 1)", value = "0.5, 0.3, 0.2")),
        conditionalPanel("input.bc_probs_source == 'prior'",
          numericInput("bc_n_outcomes", "Number of edit outcomes", value = 3, min = 1, max = 20, step = 1),
          numericInput("bc_probs_concentration", "Concentration of the Dirichlet", value = 1, min = 0.01, step = 0.5),
          helpText("The frequencies are drawn once, with the seed + 1. Concentration 1: every set of frequencies is equally likely; below 1: a few outcomes dominate; above 1: the outcomes are similar.")),
        numericInput("bc_edit_rate", "Editing rate (per unit of time)", value = 1, min = 0, step = 0.1),
        numericInput("bc_silencing_rate", "Silencing rate (per unit of time)", value = 0, min = 0, step = 0.01),
        numericInput("bc_dropout_prob", "Dropout probability at the tips", value = 0, min = 0, max = 1, step = 0.05),
        conditionalPanel("input.bc_model == 'nonseq'",
          checkboxInput("bc_window", "Restrict editing to a time window", value = FALSE),
          conditionalPanel("input.bc_window",
            helpText("Times are measured before the present. By default editing starts at the origin of the tree and lasts until the present."),
            numericInput("bc_edit_height", "Start of the editing window (time before present)", value = 3, min = 0, step = 0.5),
            numericInput("bc_edit_duration", "Length of the editing window", value = 1, min = 0, step = 0.5))),
        conditionalPanel("input.bc_model == 'generic'",
          hr(),
          h6("Editing rate: time and type dependence"),
          helpText("The editing rate is multiplied by a time profile (time measured from the origin) and by a factor per type."),
          selectInput("bc_profile", "Time profile",
                      choices = c("Constant" = "constant", "Editing window" = "window", "Exponential decay" = "decay")),
          conditionalPanel("input.bc_profile == 'window'",
            fluidRow(column(6, numericInput("bc_win_start", "Window start", value = 0, min = 0, step = 0.5)),
                     column(6, numericInput("bc_win_end", "Window end", value = 2, min = 0, step = 0.5)))),
          conditionalPanel("input.bc_profile == 'decay'",
            numericInput("bc_half_life", "Half-life of the rate", value = 1, min = 0.01, step = 0.5)),
          checkboxInput("bc_by_type", "Different editing rate per type", value = FALSE),
          conditionalPanel("input.bc_by_type", uiOutput("bc_mult_ui"))),
        hr(),
        actionButton("sim_bc", "Simulate barcodes", class = "btn-primary"),
        downloadButton("download_bc", "Download (.rds)")
      ),
      card(card_header("Barcodes of the cells (grey = unedited, colours = edit outcomes, pink = missing)"),
           card_body(fillable = FALSE, style = "overflow-y: auto;",
                     uiOutput("bc_note"),
                     plotOutput("bc_plot", height = "auto", fill = FALSE),
                     verbatimTextOutput("bc_summary")))
    )
  ),

  # ---------------------------------------------------------------- Gene expression: expression and counts
  nav_panel(
    "Gene expression", value = "expression",
    navset_tab(
      id = "expr_subtab",
      nav_panel(
        "Simulate expression", value = "expr",
        layout_sidebar(
          sidebar = sidebar(
            width = 420,
            h5("BM / OU1 / OUM process per gene"),
            numericInput("expr_seed", "Random seed (expression)", value = 1, step = 1),
            helpText("The expression evolves on the tree with this seed; the parameters of the genes are drawn from the priors with the seed + 1."),
            radioButtons("expr_source", "Parameters of the genes",
                         choices = c("Drawn from priors" = "priors", "Uploaded table: one row per gene" = "table")),
            conditionalPanel("input.expr_source == 'table'",
              fileInput("expr_table", "Table of gene parameters (.csv)", accept = c(".csv", "text/csv", "text/plain")),
              uiOutput("expr_table_help"),
              downloadButton("expr_template", "Download a template (.csv)"),
              hr()),
            conditionalPanel("input.expr_source == 'priors'",
              sliderInput("n_gene", "Number of genes (one independent process per gene)", min = 100, max = 20000, value = 1000, step = 100),
              helpText("Large numbers of genes take longer to simulate and plot; with many cells the number of genes is limited (cells x genes at most 20 million)."),
              sliderInput("pi_ou", "Probability that a gene is OU (otherwise BM)", min = 0, max = 1, value = 0.5, step = 0.05),
              helpText("Within OU genes, OU1 vs OUM is controlled by eta: eta = 0 gives OU1, larger eta separates the optima of the types."),
              hr(),
              h6("Gene baseline level"),
              helpText("The baseline level g_g of a gene is the level around which its expression evolves (theta_base = log g_g), so it belongs to the expression. On the count scale, the mean of g_g is:"),
              numericInput("mean_count_gene", "Mean of the gene baseline level g_g", value = 1, min = 0.001, step = 0.5),
              sliderInput("log_g_sd", "SD of log g_g across genes", min = 0, max = 4, value = 1.5, step = 0.1),
              hr(),
              h6("Priors of the parameters"),
              helpText("Every parameter is drawn from the chosen distribution and parameters, truncated to the bounds (leave a bound empty for none). The optima of the types, theta, follow from the baseline level, eta and psi: their distribution is shown in the plots."),
              priorUI("alpha_prior", "alpha (strength of the reversion, OU genes)", selected = "Exponential", lower = 0.1, min_lower = 0),
              helpText("The lower bound of alpha is the threshold from which a gene is OU rather than a Brownian motion."),
              priorUI("sigma2_prior", "sigma2 (diffusion)", selected = "Exponential", lower = 1e-8, upper = 2, min_lower = 0),
              priorUI("eta_prior", "eta (separation between the optima of the types, at least 0)", selected = "Exponential", lower = 0, min_lower = 0),
              sliderInput("psi_concentration", "psi: concentration of its Dirichlet prior (the share of each type in the deviations of the optima; < 1: optima far apart, > 1: close)",
                          min = 0.1, max = 5, value = 0.5, step = 0.1)),
            sliderInput("gap_cutoff", "Effective OUM class: largest difference between the optima > cutoff x stationary sd", min = 0.1, max = 5, value = 1, step = 0.1),
            uiOutput("collapsed_tree_note"),
            hr(),
            actionButton("sim_expr", "Simulate expression", class = "btn-primary"),
            downloadButton("download_expr", "Download (.rds)"),
            downloadButton("download_expr_truth", "True parameters (.csv)")
          ),
          navset_card_tab(
            nav_panel("Priors", card_body(fillable = FALSE, plotOutput("prior_preview", height = "320px", fill = FALSE))),
            nav_panel("Simulated expression",
                      # the heatmap grows with the number of cells: the sub-tab scrolls instead of squeezing the plots
                      card_body(fillable = FALSE, style = "overflow-y: auto;",
                                uiOutput("expr_note"),
                                plotOutput("expr_heatmap", height = "auto", fill = FALSE),
                                plotOutput("oum_params_plot", height = "380px", fill = FALSE),
                                plotOutput("thetas_plot", height = "320px", fill = FALSE),
                                plotOutput("class_plot", height = "320px", fill = FALSE))),
            nav_panel("UMAP", card_body(fillable = FALSE, plotOutput("expr_umap", height = "650px", fill = FALSE))),
            nav_panel("Summary", verbatimTextOutput("expr_summary"))
          )
        )
      ),
      nav_panel(
        "Simulate counts", value = "counts",
        layout_sidebar(
          sidebar = sidebar(
            width = 420,
            h5("Negative-binomial counts"),
            helpText("Simulates the counts at the tips from the expression of the previous sub-tab."),
            numericInput("counts_seed", "Random seed (counts)", value = 1, step = 1),
            h6("Gene level (counts per gene)"),
            checkboxInput("use_target", "Calibrate the global level to a target mean count", value = FALSE),
            conditionalPanel("input.use_target",
                             numericInput("target_mean_count", "Target mean count per gene and cell", value = 1, min = 0.001, step = 0.5)),
            h6("Sequencing depth (counts per cell)"),
            sliderInput("cell_size_sd", "SD (log scale) of the cell size factor", min = 0, max = 2, value = 0.9, step = 0.05),
            hr(),
            h6("Overdispersion (NB: var = mu + phi mu^2)"),
            radioButtons("phi_mode", "Overdispersion model",
                         choices = c("Prior, independent of the mean" = "prior", "Trend: decreases with the gene level" = "trend")),
            conditionalPanel("input.phi_mode == 'prior'",
                             priorUI("phi_prior", "Prior of phi", selected = "Lognormal", lower = 1e-8, min_lower = 0)),
            conditionalPanel("input.phi_mode == 'trend'",
                             helpText("phi ~ lognormal(log(a0 + a1 / g_g), sd), with g_g the gene baseline level."),
                             numericInput("phi_a0", "a0 (asymptotic overdispersion)", value = 0.1, min = 0, step = 0.05),
                             numericInput("phi_a1", "a1 (extra overdispersion at low levels)", value = 0.5, min = 0, step = 0.1),
                             numericInput("phi_sd", "sd (log scale) around the trend", value = 0.3, min = 0, step = 0.1)),
            hr(),
            actionButton("sim_counts", "Simulate counts", class = "btn-primary"),
            downloadButton("download_counts", "Download (.rds)"),
            downloadButton("download_counts_truth", "True gene parameters (.csv)"),
            downloadButton("download_cell_sizes", "True cell sizes (.csv)")
          ),
          navset_card_tab(
            nav_panel("Counts", plotOutput("counts_plots", height = "1000px")),
            nav_panel("Gene and cell sizes", plotOutput("sizes_plots", height = "420px"),
                      plotOutput("lowcounts_plot", height = "320px")),
            nav_panel("Summary", verbatimTextOutput("sim_summary"))
          )
        )
      )
    )
  ),

  nav_spacer(),
  nav_item(div(class = "d-flex gap-2",
               downloadButton("download_script", "Reproduce in R", class = "btn-outline-success btn-sm"),
               downloadButton("download_all", "Download all", class = "btn-success btn-sm")))
)

# ======================================================================================================
# Server
# ======================================================================================================

server <- function(input, output, session) {

  # ---- Tree ------------------------------------------------------------------------------------
  # The settings of every type (and the matrices of transition probabilities) are inputs that are drawn for the current
  # number of types. When the number changes the grid is drawn again and keeps what was typed for the types that remain.
  # The tree is also simulated when the page opens, before these inputs exist: then the defaults are used.
  n_types_ui <- reactive({
    k <- input$n_types
    if (is.null(k) || is.na(k)) 2L else max(1L, min(10L, as.integer(round(k))))
  })
  keep <- function(id, default) {
    v <- isolate(input[[id]])
    if (is.null(v) || length(v) != 1 || is.na(v)) default else v
  }

  output$type_params <- renderUI({
    k <- n_types_ui()
    tagList(
      lapply(seq_len(k), function(i) {
        tagList(
          strong(paste("Type", i - 1L)),
          fluidRow(
            column(4, numericInput(paste0("scale", i), "Scale", value = keep(paste0("scale", i), 0.2), min = 0.01, step = 0.1)),
            column(4, numericInput(paste0("shape", i), "Shape", value = keep(paste0("shape", i), 5), min = 0.01, step = 0.5)),
            column(4, numericInput(paste0("death", i), "Death prob.", value = keep(paste0("death", i), 0.05), min = 0, max = 0.99, step = 0.01))))
      }),
      selectInput("origin_type", "Type of the first cell (the origin)",
                  choices = stats::setNames(seq_len(k) - 1L, paste("type", seq_len(k) - 1L)),
                  selected = keep("origin_type", 0)))
  })

  output$trans_ui <- renderUI({
    d <- default_transitions(n_types_ui())
    # the defaults of the matrices depend on the number of types (the last type is the one that only divides
    # symmetrically), so a change of the number of types starts again from the defaults
    fresh <- function(id, default) default
    tagList(trans_matrix_ui("asym", "Asymmetric transition probabilities", d$asym, fresh),
            trans_matrix_ui("sym", "Symmetric transition probabilities", d$sym, fresh))
  })

  # the matrix of an input family (rows: from, columns: to), the defaults where an input does not exist (yet)
  read_matrix <- function(prefix, default) {
    m <- default
    for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) {
      v <- input[[paste0(prefix, "_", i, j)]]
      if (!is.null(v)) m[i, j] <- v
    }
    m
  }
  read_types <- function(prefix, n, default) {
    vapply(seq_len(n), function(i) { v <- input[[paste0(prefix, i)]]; if (is.null(v)) default else v }, numeric(1))
  }

  output$trans_check <- renderUI({
    d <- default_transitions(n_types_ui())
    asym <- read_matrix("asym", d$asym); sym <- read_matrix("sym", d$sym)
    if (anyNA(asym) || anyNA(sym)) return(NULL)
    sums <- rowSums(asym) + rowSums(sym)
    ok <- all(abs(sums - 1) < 1e-8)
    div(class = if (ok) "text-success" else "text-danger",
        sprintf("Row sums (asymmetric + symmetric), by type: %s %s", paste(sprintf("%.3f", sums), collapse = ", "),
                if (ok) "(OK)" else "(must all be 1)"))
  })

  # the tree, the R code that simulated it and its number of types
  tree_run <- eventReactive(input$sim_tree, {
    k <- n_types_ui()
    d <- default_transitions(k)
    asym <- read_matrix("asym", d$asym); sym <- read_matrix("sym", d$sym)
    scale <- read_types("scale", k, 0.2); shape <- read_types("shape", k, 5); death <- read_types("death", k, 0.05)
    origin_type <- as.integer(input$origin_type %||% 0L)
    by_time <- identical(input$tree_stop, "time")
    validate(
      need(!anyNA(c(asym, sym)), "All transition probabilities must be filled in."),
      need(all(abs(rowSums(asym) + rowSums(sym) - 1) < 1e-8), "For each type, the asymmetric and symmetric transition probabilities must sum to 1."),
      need(!anyNA(c(scale, shape, death, input$tree_seed, input$rho)), "All tree parameters must be filled in."),
      need(!by_time || (!is.na(input$origin_time) && input$origin_time > 0), "The time of the origin must be positive."),
      need(by_time || (!is.na(input$n_tips %||% 100) && (input$n_tips %||% 100) >= 2 && (input$n_tips %||% 100) <= MAX_TIPS),
           sprintf("The number of tips must be between 2 and %s.", format(MAX_TIPS, big.mark = ","))),
      need(!isTRUE(input$normalize_height) || (!is.na(input$tree_height) && input$tree_height > 0),
           "The tree height must be positive.")
    )
    p <- list(stop = if (by_time) "time" else "ntaxa", ntaxa = input$n_tips %||% 100, origin_time = input$origin_time,
              seed = input$tree_seed, scale = scale, shape = shape, death = death, rho = input$rho,
              origin_type = origin_type, collapse = !isTRUE(input$keep_single), asym = asym, sym = sym,
              normalize = isTRUE(input$normalize_height), height = input$tree_height)
    code <- tree_code(p)
    tr <- withProgress(message = "Simulating tree...", value = 0.5, {
      tryCatch(eval_code(code)$tree,
               error = function(e) { message("Tree simulation failed: ", conditionMessage(e)); NULL })
    })
    validate(need(!is.null(tr), if (by_time)
      "Tree simulation failed (for example, fewer than 2 tips were sampled). Try another seed, a longer time or a lower death probability."
      else "Tree simulation failed (for example, too many lineages died). Try another seed or a lower death probability."))
    list(tree = tr, code = code, n_types = k)
  })
  tree_sim <- reactive(tree_run()$tree)
  # number of types of the tree that is shown (or of the settings, while there is no tree)
  n_types_tree <- reactive(tryCatch(tree_run()$n_types, error = function(e) n_types_ui()))

  # tip labels need room: in the rectangular and slanted layouts the plot grows with the number of tips
  # (15 px per tip) so that the labels do not overlap; trees with many tips are drawn circular, without labels
  n_tips_tree <- reactive(tryCatch(ape::Ntip(tree_sim()@phylo), error = function(e) 0))
  tree_layout <- reactive({
    layout <- input$tree_layout %||% "auto"
    if (layout == "auto") layout <- if (n_tips_tree() > MAX_TIPS_RECTANGULAR) "circular" else "rectangular"
    layout
  })
  tips_in_rows <- reactive(isTRUE(input$tree_tips) && n_tips_tree() <= MAX_TIP_LABELS && tree_layout() %in% c("rectangular", "slanted"))

  output$tree_plot <- renderPlot({
    validate(need(isTRUE(input$sim_tree > 0), "Press 'Simulate tree' to simulate a tree."))
    plot_sim_tree(tree_sim(), layout = tree_layout(), show_tips = input$tree_tips, show_axis = input$tree_axis,
                  tip_size = if (tips_in_rows()) 3 else 2, max_tip_labels = MAX_TIP_LABELS)
  }, height = function() {
    if (tips_in_rows()) max(650, 15 * n_tips_tree()) else if (tree_layout() %in% c("circular", "equal_angle")) 800 else 650
  })

  output$tree_info <- renderText({
    validate(need(isTRUE(input$sim_tree > 0), ""))
    format(summarize_tree(tree_sim()))
  })

  # the optima of the expression depend on the type of each branch: warn if type changes can be hidden
  output$collapsed_tree_note <- renderUI({
    if (!isTRUE(input$keep_single)) {
      div(class = "text-warning",
          "The tree is collapsed: type changes at divisions with a dead lineage are hidden on its branches, ",
          "which take the type of the node they end in. Tick 'Keep single-child nodes' in the tree tab to keep them.")
    }
  })

  # ---- Mutations (SNVs) ------------------------------------------------------------------------
  # A layer is simulated on the tree. Its result only counts if its button was pressed after the last time the
  # tree was simulated: `base` holds the value of every button at that moment, so a new tree sends all layers back to
  # the state in which nothing has been simulated (the old results stay cached but are never shown). The counts are
  # simulated from the expression, so they need the expression of the current tree too.
  layer_buttons <- c("sim_snv", "sim_bin", "sim_bc", "sim_expr", "sim_counts")
  base <- do.call(reactiveValues, setNames(as.list(rep(0, length(layer_buttons))), layer_buttons))
  observeEvent(tree_sim(), {
    for (id in layer_buttons) base[[id]] <- if (is.null(input[[id]])) 0 else input[[id]]
  })
  layer_ready <- function(id) {
    isTRUE(input[[id]] > base[[id]]) && (id != "sim_counts" || layer_ready("sim_expr"))
  }
  need_press <- function(id, message) validate(need(layer_ready(id), message))
  fail <- function(what, e) validate(need(FALSE, paste0(what, " failed: ", conditionMessage(e))))

  # Download of the result of one layer, as an .rds file. The download links are in tabs that are not visible when the
  # page loads, and Shiny gives no URL to the output of a hidden element (the button stays disabled), so they must not
  # be suspended when hidden. Before the layer has been simulated, the download explains why there is nothing to get.
  layer_download <- function(id, prefix, pressed, layer, what) {
    output[[id]] <- downloadHandler(
      filename = function() paste0(prefix, format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
      content = function(file) {
        if (!isTRUE(pressed())) {
          showNotification(paste0("Nothing to download yet: ", what), type = "warning")
          stop("Nothing to download yet.")
        }
        saveRDS(layer(), file)
      })
    outputOptions(output, id, suspendWhenHidden = FALSE)
  }

  # Download of a table of true values, as a .csv file
  csv_download <- function(id, prefix, pressed, data, what) {
    output[[id]] <- downloadHandler(
      filename = function() paste0(prefix, format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        if (!isTRUE(pressed())) {
          showNotification(paste0("Nothing to download yet: ", what), type = "warning")
          stop("Nothing to download yet.")
        }
        utils::write.csv(data(), file, row.names = FALSE)
      })
    outputOptions(output, id, suspendWhenHidden = FALSE)
  }

  snv_clock <- clockServer("snv_clock")

  snv_sim <- eventReactive(input$sim_snv, {
    tr <- tree_sim()
    rates <- c(input$snv_r_AC, input$snv_r_AG, input$snv_r_AT, input$snv_r_CG, input$snv_r_CT, input$snv_r_GT)
    pi_raw <- c(input$snv_pi_A, input$snv_pi_C, input$snv_pi_G, input$snv_pi_T)
    validate(
      need(!anyNA(c(rates, pi_raw, input$snv_epsilon, input$snv_delta, input$snv_n_sites, input$snv_seed)),
           "All mutation parameters must be filled in."),
      need(all(rates >= 0) && sum(rates) > 0, "The exchangeabilities must be non-negative and not all zero."),
      need(all(pi_raw > 0), "The nucleotide frequencies must be positive."),
      need(input$snv_epsilon >= 0 && input$snv_epsilon <= 1 && input$snv_delta >= 0 && input$snv_delta <= 1,
           "epsilon and delta must be probabilities.")
    )
    clock <- snv_clock()
    code <- snv_code(list(seed = input$snv_seed, clock = clock$code, n_sites = input$snv_n_sites, rates = rates, pi_raw = pi_raw,
                          origin = input$snv_origin, epsilon = input$snv_epsilon, delta = input$snv_delta))
    withProgress(message = "Simulating mutations...", value = 0.5, {
      tryCatch({
        env <- eval_code(code, list(tree = tr))
        list(tree = env$ct, snv = env$snv, code = code)   # the tree with the distances of the clock draw
      }, error = function(e) fail("The mutation simulation", e))
    })
  })

  output$snv_note <- subsample_notice(tree_sim)
  output$snv_plot <- renderPlot({
    need_press("sim_snv", "Press 'Simulate mutations' to simulate SNVs on the tree.")
    s <- snv_sim()
    plot_sim_snv(s$tree, s$snv, errors = isTRUE(input$snv_errors), distance = isTRUE(input$snv_distance),
                 genotypes = isTRUE(input$snv_genotypes))
  }, height = layer_plot_height(tree_sim, per_tip = 12))

  output$snv_summary <- renderText({ need_press("sim_snv", ""); format(summarize_snv(snv_sim()$tree, snv_sim()$snv)) })
  layer_download("download_snv", "sim_snv_", function() layer_ready("sim_snv"), snv_sim, "press 'Simulate mutations' first.")

  # ---- Methylation (binary sites) --------------------------------------------------------------
  bin_clock <- clockServer("bin_clock")

  bin_sim <- eventReactive(input$sim_bin, {
    tr <- tree_sim()
    validate(
      need(!anyNA(c(input$bin_n_sites, input$bin_n_sequences, input$bin_lambda, input$bin_alpha, input$bin_beta, input$bin_seed)),
           "All methylation parameters must be filled in."),
      need(input$bin_lambda > 0, "lambda must be positive."),
      need(input$bin_alpha >= 0 && input$bin_alpha <= 1 && input$bin_beta >= 0 && input$bin_beta <= 1,
           "alpha and beta must be probabilities.")
    )
    clock <- bin_clock()
    # the true sequence can only be recovered for one sequence: with several, the error draws of one sequence shift the
    # random numbers of the next
    code <- bin_code(list(seed = input$bin_seed, clock = clock$code, n_sites = input$bin_n_sites,
                          n_sequences = input$bin_n_sequences, lambda = input$bin_lambda, origin = input$bin_origin,
                          alpha = input$bin_alpha, beta = input$bin_beta))
    withProgress(message = "Simulating methylation...", value = 0.5, {
      tryCatch({
        env <- eval_code(code, list(tree = tr))
        list(tree = env$ct, bin = env$bin, true = env$bin_true, code = code)
      }, error = function(e) fail("The methylation simulation", e))
    })
  })

  output$bin_note <- subsample_notice(tree_sim)
  output$bin_plot <- renderPlot({
    need_press("sim_bin", "Press 'Simulate methylation' to simulate binary sites on the tree.")
    s <- bin_sim()
    plot_sim_binary(s$tree, s$bin, errors = if (isTRUE(input$bin_errors)) s$true, distance = isTRUE(input$bin_distance))
  }, height = layer_plot_height(tree_sim, per_tip = 8))

  output$bin_summary <- renderText({ need_press("sim_bin", ""); { b <- bin_sim(); format(summarize_binary(b$tree, b$bin, b$true)) } })
  layer_download("download_bin", "sim_methylation_", function() layer_ready("sim_bin"), bin_sim, "press 'Simulate methylation' first.")

  # ---- Barcodes -------------------------------------------------------------------------------
  # the multipliers of the editing rate of the types of the tree (default: 1 for type 0, 2 for the others)
  output$bc_mult_ui <- renderUI({
    k <- n_types_tree()
    tagList(lapply(seq_len(k), function(i) {
      numericInput(paste0("bc_mult", i), paste("Factor, type", i - 1L), value = keep(paste0("bc_mult", i), if (i == 1) 1 else 2),
                   min = 0, step = 0.5)
    }))
  })

  bc_sim <- eventReactive(input$sim_bc, {
    tr <- tree_sim()
    by_prior <- identical(input$bc_probs_source, "prior")
    probs_fixed <- if (!by_prior) parse_probs(input$bc_edit_probs)
    n_outcomes <- input$bc_n_outcomes %||% 3
    concentration <- input$bc_probs_concentration %||% 1
    k <- n_types_tree()
    mult <- vapply(seq_len(k), function(i) input[[paste0("bc_mult", i)]] %||% (if (i == 1) 1 else 2), numeric(1))
    validate(
      need(by_prior || !is.null(probs_fixed), "The edit outcome frequencies must be numbers separated by commas, and sum to 1."),
      need(!by_prior || (!is.na(n_outcomes) && n_outcomes >= 1 && n_outcomes == round(n_outcomes) && !is.na(concentration) && concentration > 0),
           "The number of edit outcomes must be a positive whole number, and the concentration positive."),
      need(!anyNA(c(input$bc_edit_rate, input$bc_silencing_rate, input$bc_dropout_prob, input$bc_seed)),
           "All barcode parameters must be filled in."),
      need(input$bc_edit_rate >= 0 && input$bc_silencing_rate >= 0, "The rates must be non-negative."),
      need(input$bc_dropout_prob >= 0 && input$bc_dropout_prob <= 1, "The dropout probability must be in [0, 1]."),
      need(!isTRUE(input$bc_by_type) || (!anyNA(mult) && all(mult >= 0)), "The editing factors of the types must be non-negative numbers.")
    )
    mode <- input$bc_model
    profile <- if (mode == "generic") switch(input$bc_profile,
      constant = NULL,
      window = sprintf("rate_window(%s, %s)", fmt_num(input$bc_win_start), fmt_num(input$bc_win_end)),
      decay = sprintf("rate_exp_decay(%s)", fmt_num(input$bc_half_life)))
    code <- bc_code(list(seed = input$bc_seed, model = mode, probs_source = if (by_prior) "prior" else "fixed",
                         edit_probs = probs_fixed, n_outcomes = n_outcomes, concentration = concentration,
                         edit_rate = input$bc_edit_rate,
                         silencing_rate = input$bc_silencing_rate, dropout_prob = input$bc_dropout_prob,
                         n_barcodes = input$bc_n_barcodes, n_sites = input$bc_n_sites, n_targets = input$bc_n_targets,
                         window = isTRUE(input$bc_window), edit_height = input$bc_edit_height, edit_duration = input$bc_edit_duration,
                         profile = profile, by_type = if (isTRUE(input$bc_by_type)) mult, sequential = isTRUE(input$bc_sequential)))
    withProgress(message = "Simulating barcodes...", value = 0.5, {
      tryCatch({
        env <- eval_code(code, list(tree = tr))
        n_edits <- length(env$edit_probs)
        list(tree = tr, bc = env$bc, n_edits = n_edits, edit_probs = env$edit_probs,
             missing_state = if (mode == "generic") "-" else n_edits + 1L, code = code)
      }, error = function(e) fail("The barcode simulation", e))
    })
  })

  output$bc_note <- subsample_notice(tree_sim)
  output$bc_plot <- renderPlot({
    need_press("sim_bc", "Press 'Simulate barcodes' to simulate barcodes on the tree.")
    s <- bc_sim()
    plot_sim_barcodes(s$tree, s$bc, n_edits = s$n_edits, missing_state = s$missing_state)
  }, height = layer_plot_height(tree_sim, per_tip = 8))

  output$bc_summary <- renderText({ need_press("sim_bc", ""); { b <- bc_sim(); format(summarize_barcodes(b$tree, b$bc, b$missing_state)) } })
  layer_download("download_bc", "sim_barcodes_", function() layer_ready("sim_bc"), bc_sim, "press 'Simulate barcodes' first.")

  # ---- Priors ---------------------------------------------------------------------------------
  alpha_prior  <- priorServer("alpha_prior",  list(Exponential = c(1, NA), Lognormal = c(0, 1), Normal = c(1, 1), Uniform = c(0, 3)), min_lower = 0)
  sigma2_prior <- priorServer("sigma2_prior", list(Exponential = c(1, NA), Lognormal = c(0, 1), Normal = c(1, 0.5), Uniform = c(0, 2)), min_lower = 0)
  eta_prior    <- priorServer("eta_prior",    list(Exponential = c(1, NA), Lognormal = c(0, 1), Normal = c(1, 0.5), Uniform = c(0, 2)), min_lower = 0)
  phi_prior    <- priorServer("phi_prior",    list(Exponential = c(5, NA), Lognormal = c(-2, 2), Normal = c(0.5, 0.3), Uniform = c(0, 2)), min_lower = 0)

  # the priors of the parameters of the genes, with their bounds
  expr_priors <- reactive(list(alpha = alpha_prior(), sigma2 = sigma2_prior(), eta = eta_prior()))

  output$prior_preview <- renderPlot({
    set.seed(1)
    pr <- expr_priors()
    p <- tryCatch(plot_prior_draws(lapply(pr, function(p) rsim_prior(5000, p))),
                  error = function(e) conditionMessage(e))
    validate(need(!is.character(p), if (is.character(p)) p else ""))
    p
  })

  # ---- Expression --------------------------------------------------------------------------------
  # help on the format of the uploaded table, for the number of types of the tree
  output$expr_table_help <- renderUI({
    theta <- paste(paste0("theta_", seq_len(n_types_tree()) - 1L), collapse = ", ")
    helpText("One row per gene; the number of rows is the number of genes. Columns: ", code("alpha"),
             " (reversion strength; 0 is a Brownian motion), ", code("sigma2"), " (diffusion, > 0), ", code(theta),
             " (optimum of each type, log scale; can be empty for Brownian motions), and optionally ", code("gene"),
             " (names) and ", code("root_value"), " (default: the optimum of the type of the root for OU genes, 0 for Brownian motions).")
  })
  output$expr_template <- downloadHandler(
    filename = function() "gene_parameters_template.csv",
    content = function(file) utils::write.csv(gene_table_template(n_types_tree()), file, row.names = FALSE))
  outputOptions(output, "expr_template", suspendWhenHidden = FALSE)

  # the memory limits the number of cells times the number of genes
  check_size <- function(n_cells, n_genes) {
    validate(need(n_cells * n_genes <= MAX_CELLS_X_GENES,
                  sprintf("%s cells x %s genes is too much: at most %s (cells x genes). Use fewer genes or cells.",
                          format(n_cells, big.mark = ","), format(n_genes, big.mark = ","), format(MAX_CELLS_X_GENES, big.mark = ","))))
  }

  expr_sim <- eventReactive(input$sim_expr, {
    tr <- tree_sim()
    k <- n_types_tree()
    validate(need(!anyNA(c(input$expr_seed, input$gap_cutoff)), "All expression parameters must be filled in."))
    if (identical(input$expr_source, "table")) {
      f <- input$expr_table
      validate(need(!is.null(f), "Choose a CSV file with the parameters of the genes, or draw them from priors."))
      tab <- tryCatch({
        tab <- read_gene_csv(f$datapath)
        gene_params_from_table(tab, tr, k, input$gap_cutoff)   # validates the table, with a clear message
        tab
      }, error = function(e) validate(need(FALSE, conditionMessage(e))))
      check_size(ape::Ntip(tr@phylo), nrow(tab))
      withProgress(message = "Simulating expression...", value = 0.5, {
        tryCatch(
          simulate_expression(tr, k, "table", table = tab, table_name = f$name, seed = input$expr_seed, gap_cutoff = input$gap_cutoff),
          error = function(e) fail("The expression simulation", e))
      })
    } else {
      validate(need(!anyNA(c(input$n_gene, input$mean_count_gene, input$log_g_sd)), "All expression parameters must be filled in."))
      check_size(ape::Ntip(tr@phylo), input$n_gene)
      pr <- expr_priors()
      withProgress(message = "Simulating expression...", value = 0.5, {
        tryCatch(
          simulate_expression(
            tr, k, "priors", seed = input$expr_seed, n_genes = input$n_gene,
            pi_ou = input$pi_ou, gap_cutoff = input$gap_cutoff,
            mean_count_gene = input$mean_count_gene, log_g_sd = input$log_g_sd,
            priors = pr, psi_concentration = input$psi_concentration),
          error = function(e) fail("The expression simulation", e))
      })
    }
  })

  expr_params_plots <- reactive({
    need_press("sim_expr", "Press 'Simulate expression' to simulate the expression of the genes on the tree.")
    plot_sim_params(expr_sim()$truth)
  })
  # underlying expression of every gene (columns, grouped by class) in every cell (rows, aligned with the tips of the tree)
  output$expr_heatmap <- renderPlot({
    need_press("sim_expr", "Press 'Simulate expression' to simulate the expression of the genes on the tree.")
    es <- expr_sim()
    plot_sim_expression(es$expr, es$tree, es$truth)
  }, height = layer_plot_height(function() expr_sim()$tree, per_tip = 8))
  output$expr_note <- subsample_notice(function() expr_sim()$tree)
  output$oum_params_plot <- renderPlot({ expr_params_plots()$oumparams })
  output$thetas_plot     <- renderPlot({
    p <- expr_params_plots()$thetas
    validate(need(!is.null(p), "No OU genes: there are no optima to show."))
    p
  })
  output$class_plot      <- renderPlot({ expr_params_plots()$modelclass })

  # UMAP of the cells from their latent expression, coloured by cell type (the seed of the plot function is its default)
  output$expr_umap <- renderPlot({
    need_press("sim_expr", "Press 'Simulate expression' to simulate the expression of the genes on the tree.")
    es <- expr_sim()
    withProgress(message = "Computing the UMAP...", value = 0.5, plot_sim_umap(es$expr, es$tree))
  })

  output$expr_summary <- renderText({
    need_press("sim_expr", "Press 'Simulate expression' to simulate the expression of the genes on the tree.")
    es <- expr_sim()
    format(summarize_expression(es$expr, es$tree, es$truth))
  })
  layer_download("download_expr", "sim_expression_", function() layer_ready("sim_expr"), expr_sim, "press 'Simulate expression' first.")
  # the true parameters of every gene: drawn from the priors or read from the table (plus the derived classes)
  csv_download("download_expr_truth", "sim_expression_true_parameters_", function() layer_ready("sim_expr"),
               function() expr_sim()$truth, "press 'Simulate expression' first.")

  # ---- Counts ------------------------------------------------------------------------------------
  counts_sim <- eventReactive(input$sim_counts, {
    need_press("sim_expr", "Simulate the expression first (sub-tab 'Simulate expression'): the counts are simulated from it.")
    es <- expr_sim()
    validate(
      need(!anyNA(c(input$counts_seed, input$cell_size_sd)), "All count parameters must be filled in.")
    )
    phi <- if (input$phi_mode == "prior") phi_prior()
    withProgress(message = "Simulating counts...", value = 0.5, {
      tryCatch(
        simulate_counts(
          es, cell_size_sd = input$cell_size_sd,
          target_mean_count = if (isTRUE(input$use_target)) input$target_mean_count else NULL,
          phi = phi, phi_mode = input$phi_mode, phi_a0 = input$phi_a0, phi_a1 = input$phi_a1, phi_sd = input$phi_sd,
          seed = input$counts_seed),
        error = function(e) fail("The counts simulation", e))
    })
  })

  plots <- reactive({
    s <- counts_sim()
    withProgress(message = "Building plots...", value = 0.5, plot_sim_results(s))
  })

  need_counts <- function() {
    need_press("sim_expr", "Simulate the expression first (sub-tab 'Simulate expression'): the counts are simulated from it.")
    need_press("sim_counts", "Press 'Simulate counts' to simulate the counts at the tips from the expression.")
  }

  output$counts_plots <- renderPlot({
    need_counts()
    plot_sim_counts_panel(plots())
  })

  output$sizes_plots <- renderPlot({
    need_counts()
    plot_sim_sizes_panel(plots())
  })

  output$lowcounts_plot <- renderPlot({ need_counts(); plots()$lowcounts })

  output$sim_summary <- renderText({ need_counts(); { cs <- counts_sim(); format(summarize_counts(cs$counts, cs$truth, cs$cell_size, cs$shift)) } })

  layer_download("download_counts", "sim_counts_", function() layer_ready("sim_counts"), counts_sim, "press 'Simulate counts' first.")
  # true values of the counts: the overdispersion and the (log-scale) parameters of every gene, and the size of every cell
  cell_sizes_table <- function(s) {
    data.frame(cell = names(s$cell_size), true_size_factor = unname(s$cell_size),
               observed_size_factor = unname(s$cell_size_obs[names(s$cell_size)]),
               total_counts = unname(s$cell_total[names(s$cell_size)]))
  }
  csv_download("download_counts_truth", "sim_counts_true_gene_parameters_", function() layer_ready("sim_counts"),
               function() counts_sim()$truth, "press 'Simulate counts' first.")
  csv_download("download_cell_sizes", "sim_counts_true_cell_sizes_", function() layer_ready("sim_counts"),
               function() cell_sizes_table(counts_sim()), "press 'Simulate counts' first.")

  # ---- Tabs ---------------------------------------------------------------------------------------
  # Nothing but the tree can be used until there is a tree, and the counts need the expression. The page starts with
  # everything disabled (see the header of the UI); here the server enables and disables the tabs and the buttons.
  observe({
    has_tree <- !is.null(tryCatch(tree_run()$tree, error = function(e) NULL))
    session$sendCustomMessage("set_disabled", list(tabs = c("mutations", "methylation", "barcodes", "expression"),
                                                   buttons = c("download_script", "download_all"), disabled = !has_tree))
  })
  observe({
    has_expression <- isTRUE(tryCatch(layer_ready("sim_expr"), error = function(e) FALSE))
    session$sendCustomMessage("set_disabled", list(tabs = "counts", disabled = !has_expression))
    # the counts sub-tab cannot stay selected when it is disabled (for example after a new tree)
    if (!has_expression) bslib::nav_select("expr_subtab", "expr", session = session)
  })

  # ---- Script ------------------------------------------------------------------------------------
  # The R code of what is shown: the tree and every layer that has been simulated on it. Every layer ran from this code
  # (see the top of the file), so the script gives the same results with the same version of the package.
  script_text <- function() {
    tree <- tree_run()
    layer <- function(id, title, fun) {
      if (!layer_ready(id)) return(NULL)
      code <- tryCatch(fun()$code, error = function(e) NULL)
      if (!is.null(code)) c(paste0("# ---- ", title, " ----"), code, "")
    }
    header <- c(
      "# R script written by the scTreeSim app. It reproduces the simulation that was shown when it was downloaded:",
      "# the same functions of the package, with the same settings and seeds. Run it from top to bottom, as every section",
      "# uses the objects of the sections before it.",
      sprintf("# scTreeSim %s, %s, written %s", as.character(utils::packageVersion("scTreeSim")), R.version.string,
              format(Sys.time(), "%Y-%m-%d %H:%M")),
      "",
      "library(scTreeSim)",
      'RNGkind("Mersenne-Twister", "Inversion", "Rejection")  # the random number generator of the app',
      "")
    paste(c(header, "# ---- 1. Tree ----", tree$code, "",
            layer("sim_snv", "2. Mutations", snv_sim), layer("sim_bin", "3. Methylation", bin_sim),
            layer("sim_bc", "4. Barcodes", bc_sim), layer("sim_expr", "5. Expression", expr_sim),
            layer("sim_counts", "6. Counts", counts_sim)), collapse = "\n")
  }
  output$download_script <- downloadHandler(
    filename = function() paste0("scTreeSim_reproduce_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
    content = function(file) {
      text <- tryCatch(script_text(), error = function(e) {
        showNotification("There is no tree to reproduce: simulate the tree first.", type = "warning")
        stop("No tree.")
      })
      writeLines(text, file)
    })
  outputOptions(output, "download_script", suspendWhenHidden = FALSE)

  # ---- All results ------------------------------------------------------------------------------
  # one .rds with a list: the tree and the result of every layer that has been simulated (NULL for the others), the
  # true parameters of the expression and the counts as tables (also inside the results), and the script.
  # Every layer keeps the tree it was simulated on (with the clock distances for mutations and methylation).
  output$download_all <- downloadHandler(
    filename = function() paste0("sim_all_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
    content = function(file) {
      result <- function(id, layer) if (layer_ready(id)) tryCatch(layer(), error = function(e) NULL)
      expression <- result("sim_expr", expr_sim)
      counts <- result("sim_counts", counts_sim)
      saveRDS(list(tree = tryCatch(tree_sim(), error = function(e) NULL),
                   mutations = result("sim_snv", snv_sim),
                   methylation = result("sim_bin", bin_sim),
                   barcodes = result("sim_bc", bc_sim),
                   expression = expression,
                   counts = counts,
                   true = list(expression_genes = expression$truth,
                               counts_genes = counts$truth,
                               counts_cells = if (!is.null(counts)) cell_sizes_table(counts)),
                   script = tryCatch(script_text(), error = function(e) NULL)),
              file)
    }
  )
  outputOptions(output, "download_all", suspendWhenHidden = FALSE)
}

shinyApp(ui, server)
