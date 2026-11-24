# Summary statistics of the simulated layers. Every function returns the numbers as a list (`stats`) and the text
# that the Shiny app shows (`text`), as an object of class "sim_summary": print() and format() give the text.

.sim_summary <- function(stats, text) {
  structure(list(stats = stats, text = text), class = "sim_summary")
}

#' @export
print.sim_summary <- function(x, ...) {
  cat(x$text, sep = "\n")
  invisible(x)
}

#' @export
format.sim_summary <- function(x, ...) {
  paste(x$text, collapse = "\n")
}

#' Summary statistics of the simulated layers
#'
#' Quick sanity checks of a simulated tree and of the layers on top of it, which
#' are also what the Shiny app (\code{\link{launch_gui}}) reports as text. Every
#' function returns an object of class \code{sim_summary}: a list with the numbers
#' (\code{stats}) and the text (\code{text}). \code{print()} shows the text and
#' \code{format()} gives it as one string.
#'
#' @param tree a treedata object, from \code{\link{sim_adb_ntaxa_samp}} or
#'   \code{\link{sim_adb_origin_samp}} (with the distances of a clock for the
#'   mutations and the methylation: see \code{\link{clock_tree}})
#' @param snv the output of \code{\link{sim_gt16_nt_diploid}}: a list with the data
#'   frames \code{true} and \code{observed}, in either format
#' @param bin the output of \code{\link{sim_binary}}
#' @param true optionally, the output of \code{\link{sim_binary}} without error
#'   (same seed, \code{alpha = 0} and \code{beta = 0}), to report the sites flipped
#'   by the error
#' @param bc the output of \code{\link{sim_barcode_seq}},
#'   \code{\link{sim_barcode_nonseq}} or \code{\link{sim_barcode_generic}}
#' @param missing_state the state of a missing barcode site (\code{"-"} for
#'   \code{sim_barcode_generic}, the number of edit outcomes plus one for the others)
#' @param expr the expression of the genes on the tree, from
#'   \code{\link{sim_expression_matrix}}
#' @param gene_params the table of true parameters of the genes, from
#'   \code{\link{sim_gene_params}} or \code{\link{gene_params_from_table}} (for the
#'   counts, the table that \code{\link{sim_count_params}} returns, which has
#'   \code{phi})
#' @param counts the counts, a matrix with one row per gene and one column per cell
#'   (as in the app: \code{t()} of the output of
#'   \code{\link{sim_expression_counts_matrix}} without its \code{node} and
#'   \code{type} columns)
#' @param size_factor optionally, the true size factors of the cells, to compare
#'   with those that the totals of the counts give
#' @param shift optionally, the level shift of the calibration of
#'   \code{\link{sim_count_params}}
#'
#' @return an object of class \code{sim_summary}
#' @name sim_summary
#' @examples
#' set.seed(1)
#' tree <- sim_adb_ntaxa_samp(20, scale = c(0.2, 0.2), shape = c(5, 5), death_prob = c(0.05, 0.05),
#'                            asym_trans_prob = matrix(c(0, 0, 0.2, 0), 2),
#'                            sym_trans_prob = matrix(c(0.8, 0, 0, 1), 2), collapse = FALSE)
#' summarize_tree(tree)
#' summarize_tree(tree)$stats$height
NULL

#' @rdname sim_summary
#' @export
summarize_tree <- function(tree) {
  stopifnot(methods::is(tree, "treedata"))
  tips <- tree@data$status == 1
  by_type <- table(tree@data$type[tips])
  stats <- list(n_tips = length(tree@phylo$tip.label), n_internal = tree@phylo$Nnode,
                height = max(ape::node.depth.edgelength(tree@phylo)), tips_by_type = by_type)
  text <- c(sprintf("Tips: %d   Internal nodes: %d   Height: %.3f", stats$n_tips, stats$n_internal, stats$height),
            sprintf("Tips by type: %s", paste(names(by_type), by_type, sep = " = ", collapse = ", ")))
  .sim_summary(stats, text)
}

#' @rdname sim_summary
#' @export
summarize_snv <- function(tree, snv) {
  true <- .allele_matrices(tree, snv, "true")
  obs <- .allele_matrices(tree, snv, "observed")
  het <- true$allele1 != true$allele2
  genotype <- matrix(paste0(true$allele1, true$allele2), nrow(true$allele1))
  variable <- sum(apply(genotype, 2, function(x) length(unique(x)) > 1))
  changed <- true$allele1 != obs$allele1 | true$allele2 != obs$allele2
  stats <- list(n_cells = nrow(genotype), n_sites = ncol(genotype), heterozygous = mean(het),
                variable_sites = variable, changed_by_error = mean(changed))
  text <- c(sprintf("Cells: %d   Sites: %d", stats$n_cells, stats$n_sites),
            sprintf("Heterozygous sites in the cells (true): %.1f%%", 100 * stats$heterozygous),
            sprintf("Sites that differ between cells (true): %d of %d", variable, stats$n_sites),
            sprintf("Genotypes changed by the error: %.1f%%", 100 * stats$changed_by_error))
  .sim_summary(stats, text)
}

#' @rdname sim_summary
#' @export
summarize_binary <- function(tree, bin, true = NULL) {
  tips <- seq_len(ape::Ntip(tree@phylo))
  cols <- grep("^binary_", names(bin), value = TRUE)
  if (length(cols) == 0) stop("`bin` has no columns binary_1, binary_2, ...", call. = FALSE)
  text <- sprintf("Cells: %d   Sequences: %d   Sites per sequence: %d", length(tips), length(cols), nchar(bin[[cols[1]]][1]))
  sequences <- list()
  for (col in cols) {
    m <- do.call(rbind, strsplit(bin[[col]][match(tips, bin$node)], "", fixed = TRUE))
    sequences[[col]] <- list(proportion_1 = mean(m == "1"), variable_sites = sum(apply(m, 2, function(x) length(unique(x)) > 1)))
    text <- c(text, sprintf("%s: %.1f%% of the sites are 1; %d sites differ between cells", col,
                            100 * sequences[[col]]$proportion_1, sequences[[col]]$variable_sites))
  }
  flipped <- NULL
  if (!is.null(true)) {
    m <- do.call(rbind, strsplit(bin[[cols[1]]][match(tips, bin$node)], "", fixed = TRUE))
    t <- do.call(rbind, strsplit(true[[cols[1]]][match(tips, true$node)], "", fixed = TRUE))
    flipped <- mean(m != t)
    text <- c(text, sprintf("Sites flipped by the error: %.1f%%", 100 * flipped))
  }
  .sim_summary(list(n_cells = length(tips), n_sequences = length(cols), n_sites = nchar(bin[[cols[1]]][1]),
                    sequences = sequences, flipped_by_error = flipped), text)
}

#' @rdname sim_summary
#' @export
summarize_barcodes <- function(tree, bc, missing_state = "-") {
  tips <- seq_len(ape::Ntip(tree@phylo))
  bc <- bc[match(tips, bc$node), , drop = FALSE]
  cols <- grep("^barcode_", names(bc), value = TRUE)
  mats <- if (length(cols) > 0) {
    stats::setNames(lapply(cols, function(col) do.call(rbind, strsplit(bc[[col]], "_", fixed = TRUE))), cols)
  } else {
    list(targets = as.matrix(bc[grep("^site_", names(bc))]))
  }
  text <- sprintf("Cells: %d", length(tips))
  barcodes <- list()
  for (nm in names(mats)) {
    m <- matrix(as.character(mats[[nm]]), nrow(mats[[nm]]))
    barcodes[[nm]] <- list(unedited = mean(m == "0"), edited = mean(m != "0" & m != as.character(missing_state)),
                           missing = mean(m == as.character(missing_state)))
    text <- c(text, sprintf("%s: %.1f%% unedited, %.1f%% edited, %.1f%% missing (silenced or dropout)", nm,
                            100 * barcodes[[nm]]$unedited, 100 * barcodes[[nm]]$edited, 100 * barcodes[[nm]]$missing))
  }
  .sim_summary(list(n_cells = length(tips), barcodes = barcodes), text)
}

#' @rdname sim_summary
#' @export
summarize_expression <- function(expr, tree, gene_params) {
  stopifnot(methods::is(tree, "treedata"), is.data.frame(expr), is.data.frame(gene_params))
  tips <- tree@data$node[tree@data$status == 1]
  x <- as.matrix(expr[match(tips, expr$node), gene_params$gene])   # cells x genes, log scale
  sd_genes <- apply(x, 2, stats::sd)
  by_class <- tapply(sd_genes, gene_params$class_eff, mean)
  stats <- list(n_genes = nrow(gene_params), n_cells = nrow(x), class = table(gene_params$class),
                class_eff = table(gene_params$class_eff), mean = mean(x), sd = stats::sd(x),
                median_sd_gene = stats::median(sd_genes), mean_sd_by_class_eff = by_class)
  tab_text <- function(tab) paste(names(tab), tab, sep = " = ", collapse = ", ")
  text <- c(sprintf("Genes: %d   Cells: %d", stats$n_genes, stats$n_cells),
            sprintf("Simulated class (BM vs OU): %s", tab_text(stats$class)),
            sprintf("Effective class (OU1 / OUM by gap): %s", tab_text(stats$class_eff)),
            sprintf("Expression at the tips (log scale): mean %.2f, sd %.2f", stats$mean, stats$sd),
            sprintf("Median sd of a gene across cells: %.2f", stats$median_sd_gene),
            "Standard deviation across cells of the genes, by effective class:",
            paste0("  ", utils::capture.output(print(by_class))))
  .sim_summary(stats, text)
}

#' @rdname sim_summary
#' @export
summarize_counts <- function(counts, gene_params, size_factor = NULL, shift = NULL) {
  stopifnot(is.matrix(counts), is.data.frame(gene_params), nrow(counts) == nrow(gene_params))
  cell_total <- colSums(counts)
  gene_total <- rowSums(counts)
  size_obs <- cell_total / mean(cell_total)
  cor_size <- NULL
  if (!is.null(size_factor)) {
    if (!is.null(names(size_factor)) && all(colnames(counts) %in% names(size_factor))) size_factor <- size_factor[colnames(counts)]
    cor_size <- stats::cor(log(size_factor), log(size_obs))
  }
  cor_gene <- stats::cor(gene_params$log_g, log(gene_total + 1))
  by_class <- tapply(rowMeans(counts), gene_params$class_eff, mean)
  stats <- list(n_genes = nrow(counts), n_cells = ncol(counts), class = table(gene_params$class),
                class_eff = table(gene_params$class_eff), cor_cell_size = cor_size, cor_gene_level = cor_gene,
                fraction_zeros = mean(counts == 0), mean_count = mean(counts), max_count = max(counts), shift = shift,
                cell_total = summary(cell_total), mean_count_by_class_eff = by_class)
  tab_text <- function(tab) paste(names(tab), tab, sep = " = ", collapse = ", ")
  text <- c(sprintf("Genes: %d   Cells: %d", stats$n_genes, stats$n_cells),
            sprintf("Simulated class (BM vs OU): %s", tab_text(stats$class)),
            sprintf("Effective class (OU1 / OUM by gap): %s", tab_text(stats$class_eff)),
            if (!is.null(cor_size)) sprintf("cor(log true cell size, log realised cell size): %.3f", cor_size),
            sprintf("cor(log g_g, log total counts of the gene): %.3f", cor_gene),
            sprintf("Fraction of zeros: %.3f", stats$fraction_zeros),
            sprintf("Mean count per gene and cell: %.3f", stats$mean_count),
            sprintf("Maximum count: %s", format(stats$max_count, big.mark = ",")),
            if (!is.null(shift)) sprintf("Level shift applied (target mean count): %.3f", shift),
            "Total counts per cell:",
            paste0("  ", utils::capture.output(print(stats$cell_total))),
            "Mean count per gene, by effective class:",
            paste0("  ", utils::capture.output(print(by_class))))
  .sim_summary(stats, text)
}
