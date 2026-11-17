# Plotting functions for simulated trees, expression and counts. They are used by
# the Shiny app (see launch_gui()) and can also be called directly. The plotting
# packages are optional (Suggests): each function checks that it can load them.

# stop with an informative message if some suggested packages are not installed
.require_pkgs <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("This function needs the package(s) ",
         paste0("'", missing_pkgs, "'", collapse = ", "),
         ", which are not installed.", call. = FALSE)
  }
}

# colours of the gene models: Brownian motion, single-optimum OU, multi-optima OU
.model_colours <- function() {
  cols <- if (requireNamespace("scico", quietly = TRUE)) {
    scico::scico(3, palette = "batlow")
  } else {
    c("#001959", "#7D7D3A", "#FACCFA")
  }
  stats::setNames(cols, c("BM", "OU1", "OUM"))
}

# One colour per type: the first two are the blue and red of the two-type trees, then a colour-blind friendly
# palette, and an evenly spaced palette for more than eight types
.type_colours <- function(types) {
  base <- c("#2C7BB6", "#D7191C", "#FDAE61", "#7B3294", "#1A9641", "#8C510A", "#4DAC26", "#666666")
  cols <- if (length(types) <= length(base)) base[seq_along(types)] else grDevices::hcl.colors(length(types), "Dark 3")
  stats::setNames(cols, types)
}

#' Plot a simulated tree coloured by type
#'
#' The type of a branch is the type of its child node.
#'
#' @param tree a treedata object with a \code{type} column in \code{tree@data},
#'   e.g. from \code{\link{sim_adb_ntaxa_samp}}
#' @param layout tree layout, passed to \code{ggtree::ggtree()}
#'   (e.g. "rectangular", "slanted", "circular", "equal_angle")
#' @param show_tips if TRUE, label the tips
#' @param tip_size size of the tip labels (in mm); with many tips, use a small
#'   size or a tall figure so that the labels do not overlap
#' @param show_axis if TRUE, add a time axis (rectangular and slanted layouts only)
#'
#' @return a ggplot object
#' @export
plot_sim_tree <- function(tree, layout = "rectangular", show_tips = FALSE, show_axis = TRUE,
                          tip_size = 2.5) {
  .require_pkgs(c("ggplot2", "ggtree"))
  stopifnot(methods::is(tree, "treedata"))

  types <- levels(factor(tree@data$type))
  cols <- .type_colours(types)

  p <- ggtree::ggtree(tree, ggplot2::aes(color = factor(.data$type)), layout = layout) +
    ggplot2::scale_color_manual(values = cols, name = "Type") +
    ggplot2::theme(legend.position = "top")
  if (show_tips) {
    p <- p + ggtree::geom_tiplab(size = tip_size, show.legend = FALSE)
  }
  if (show_axis && layout %in% c("rectangular", "slanted")) {
    p <- p + ggtree::theme_tree2() + ggplot2::theme(legend.position = "top")
  }
  return(p)
}

#' Histograms of draws from parameter priors
#'
#' @param draws named list of numeric vectors, one per parameter (the names are
#'   used as facet labels)
#'
#' @return a ggplot object
#' @export
plot_prior_draws <- function(draws) {
  .require_pkgs("ggplot2")
  if (!is.list(draws) || is.null(names(draws)) || !all(vapply(draws, is.numeric, logical(1)))) {
    stop("`draws` must be a named list of numeric vectors.", call. = FALSE)
  }

  df <- do.call(rbind, lapply(names(draws), function(nm) {
    data.frame(parameter = nm, value = draws[[nm]])
  }))
  df$parameter <- factor(df$parameter, levels = names(draws))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 50, fill = "lightgray", colour = "white") +
    ggplot2::facet_wrap(~parameter, scales = "free", nrow = 1) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank()) +
    ggplot2::xlab("Value") + ggplot2::ylab("Count")
}

#' Plots of the true parameters of the genes
#'
#' Histograms of the process parameters of the simulated genes (\code{alpha},
#' \code{sigma2}, \code{eta}, the proportions \code{psi1}, \code{psi_0}, \code{psi_1}, ... and \code{phi}, whichever are present) and a bar plot of the
#' model classes. They only need the table of true parameters, so they can be drawn after the
#' expression is simulated, before the counts.
#'
#' @param truth data frame with one row per gene and columns \code{class} (BM or OU),
#'   \code{class_eff} (BM, OU1 or OUM), \code{alpha}, \code{sigma2}, \code{eta}, and optionally
#'   \code{phi} and the proportions of the optima of the types (columns named \code{psi1}, or
#'   \code{psi_0}, \code{psi_1}, ..., one per type but the last)
#'
#' @return a named list of two ggplot objects: \code{oumparams} and \code{modelclass}
#' @export
plot_sim_params <- function(truth) {
  .require_pkgs("ggplot2")
  truth <- as.data.frame(truth)
  model_names <- .model_colours()
  top_legend <- ggplot2::theme(panel.background = ggplot2::element_blank(),
                               legend.position = "top")

  # true classes
  p_modelclasseffcount <- ggplot2::ggplot(truth,
                                          ggplot2::aes(x = .data$class, fill = .data$class_eff)) +
    ggplot2::geom_bar() +
    ggplot2::scale_fill_manual(values = model_names, name = "Effective model") +
    top_legend +
    ggplot2::xlab("Model class") + ggplot2::ylab("Number of genes")

  # true process parameters (alpha, eta and psi are not defined for BM genes,
  # psi only matters for OUM genes)
  psi <- grep("^psi_?[0-9]+$", names(truth), value = TRUE)
  params <- c(intersect(c("alpha", "sigma2", "eta", "phi"), names(truth)), psi)
  long <- do.call(rbind, lapply(params, function(p) {
    data.frame(class = truth$class, class_eff = truth$class_eff, param = p, value = truth[[p]])
  }))
  long$param <- factor(long$param, levels = params)
  long <- long[!(long$class == "BM" & long$param %in% c("alpha", "eta", psi)), ]
  long <- long[long$class_eff == "OUM" | !(long$param %in% psi), ]
  below_q <- stats::ave(long$value, long$param,
                        FUN = function(v) v <= stats::quantile(v, 0.995, na.rm = TRUE)) == 1
  long <- long[below_q, ]
  p_oumparams <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, fill = .data$class_eff)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::facet_wrap(~param, scales = "free", nrow = 1) +
    ggplot2::scale_fill_manual(values = model_names, name = "Model") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   legend.position = "top") +
    ggplot2::xlab("Value") + ggplot2::ylab("Number of genes")

  return(list(oumparams = p_oumparams, modelclass = p_modelclasseffcount))
}

#' Tree with the simulated expression of the genes in the cells
#'
#' Plots a tree and, next to it, the underlying (latent) expression at the tips, simulated by
#' \code{\link{sim_expression_matrix}}, as a heatmap of tiles with one row per cell, aligned with the tip of
#' the tree, and one column per gene. If the table of true parameters is given, the genes are grouped by their
#' effective model class (BM, OU1, OUM), with a gap between the groups.
#'
#' @inheritParams plot_sim_snv
#' @param expr the output of \code{\link{sim_expression_matrix}}: a data frame with the columns
#'   \code{node}, \code{type} and one column per gene
#' @param tree the treedata object the expression was simulated on; only its tips (the cells) are drawn
#' @param truth optional data frame with one row per gene and the columns \code{gene} (the gene
#'   names, as in \code{expr}) and \code{class_eff} (BM, OU1 or OUM). If given, the genes are grouped by
#'   class.
#' @param max_genes maximum number of genes shown; if there are more, an evenly spaced subset is drawn
#' @param gap number of empty columns between the classes of genes
#' @param heat_width width of the heatmap as a fraction of the width of the tree
#'
#' @return a ggplot object
#' @export
plot_sim_expression <- function(expr, tree, truth = NULL, max_genes = 200, gap = 2, heat_width = 0.8) {
  .require_pkgs(c("ggplot2", "ggtree"))
  stopifnot(methods::is(tree, "treedata"), is.data.frame(expr), all(c("node", "type") %in% names(expr)))
  gene_cols <- setdiff(names(expr), c("node", "type"))
  if (length(gene_cols) == 0) {
    stop("`expr` has no gene columns.", call. = FALSE)
  }

  # the tips (cells), as in the tree plot
  phy <- tree@phylo
  tips <- seq_len(ape::Ntip(phy))
  if (!all(tips %in% expr$node)) {
    stop("`expr` has no row for some tips of `tree`; was it simulated on this tree?", call. = FALSE)
  }

  genes <- if (length(gene_cols) > max_genes) {
    gene_cols[unique(round(seq(1, length(gene_cols), length.out = max_genes)))]
  } else {
    gene_cols
  }
  m <- as.matrix(expr[match(tips, expr$node), genes, drop = FALSE])  # cells x genes
  rownames(m) <- phy$tip.label[tips]

  # columns: genes, grouped by class with a gap between the groups
  if (is.null(truth)) {
    groups <- list(all = seq_along(genes))
  } else {
    class_eff <- factor(truth$class_eff[match(genes, truth$gene)], levels = c("BM", "OU1", "OUM"))
    ord <- order(class_eff)
    genes <- genes[ord]; m <- m[, ord, drop = FALSE]; class_eff <- class_eff[ord]
    groups <- split(seq_along(genes), droplevels(class_eff))
  }
  group_size <- vapply(groups, length, 1L)
  first <- c(0, cumsum(group_size + gap))[seq_along(groups)] + 1
  pos <- numeric(length(genes))
  for (k in seq_along(groups)) {
    pos[groups[[k]]] <- first[k] + seq_along(groups[[k]]) - 1
  }
  last <- first + group_size - 1

  p <- .plot_layer_tree(tree, pad = 0.02)
  tip_data <- p$data[p$data$isTip, ]
  y <- tip_data$y[match(rownames(m), tip_data$label)]
  tree_width <- diff(range(p$data$x))
  w <- heat_width * tree_width / max(pos)  # width of a column
  x0 <- max(p$data$x) + 0.04 * tree_width

  tiles <- data.frame(x = x0 + rep(pos, each = nrow(m)) * w, y = rep(y, ncol(m)),
                      expression = as.vector(m))
  frames <- data.frame(xmin = x0 + (first - 0.5) * w, xmax = x0 + (last + 0.5) * w,
                       ymin = min(y) - 0.5, ymax = max(y) + 0.5)

  p <- p +
    ggplot2::geom_tile(data = tiles, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$expression),
                       width = w, height = 1, inherit.aes = FALSE) +
    ggplot2::geom_rect(data = frames,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
                       fill = NA, color = "grey40", linewidth = 0.3, inherit.aes = FALSE) +
    ggplot2::scale_fill_viridis_c(name = "expression (log scale)") +
    ggplot2::labs(caption = sprintf("Genes: %s%d (columns), cells: %d (rows, aligned with the tips)",
                                    if (length(genes) < length(gene_cols)) "subset of " else "",
                                    length(genes), nrow(m)))
  if (!is.null(truth)) {
    p <- p + ggplot2::annotate("text", x = (frames$xmin + frames$xmax) / 2, y = max(y) + 1,
                               label = names(groups), vjust = 0, size = 3, color = "grey30")
  }
  return(p)
}

#' Result plots of a simulation of expression and counts
#'
#' Builds all the diagnostic plots of a simulation with many genes: counts
#' matrix, counts per gene and per cell, zeros, gene and cell sizes, and the
#' true process parameters.
#'
#' @param sim a list describing the simulation, with elements
#'   \itemize{
#'     \item \code{counts}: counts matrix (genes x cells) with gene ids as row
#'       names and cell names as column names
#'     \item \code{truth}: data frame with one row per gene (in the order of
#'       \code{counts}) and columns \code{gene}, \code{class} (BM or OU),
#'       \code{class_eff} (BM, OU1 or OUM), \code{log_g}, \code{alpha},
#'       \code{sigma2}, \code{eta}, \code{phi} and \code{psi1}
#'     \item \code{cell_size}: named numeric vector with the simulated size
#'       factor of each cell
#'     \item \code{cell_size_obs}: named numeric vector with the size factor
#'       realised in the counts (total counts of the cell, divided by their mean)
#'     \item \code{cell_total}, \code{gene_total}: named numeric vectors with
#'       the total counts per cell and per gene
#'   }
#' @param max_genes_heatmap maximum number of genes shown in the counts matrix;
#'   if there are more, an evenly spaced subset is shown
#'
#' @return a named list of ggplot objects: \code{oumparams} and \code{modelclass}
#'   (see \code{\link{plot_sim_params}}), \code{genesize}, \code{cellsize},
#'   \code{countsmatrix}, \code{countspergene}, \code{countspercell},
#'   \code{lowcounts}, \code{zeros_gene} and \code{zeros_cell}
#' @export
plot_sim_results <- function(sim, max_genes_heatmap = 300) {
  .require_pkgs(c("ggplot2", "ggpubr", "viridis"))

  counts <- sim$counts
  truth <- as.data.frame(sim$truth)
  model_names <- .model_colours()
  blank <- ggplot2::theme(panel.background = ggplot2::element_blank())
  no_legend <- ggplot2::theme(panel.background = ggplot2::element_blank(),
                              legend.position = "none")
  top_legend <- ggplot2::theme(panel.background = ggplot2::element_blank(),
                               legend.position = "top")

  # proportion of zeros per gene / per cell
  zeros_gene <- data.frame(gene = rownames(counts), propzeros = rowMeans(counts == 0))
  zeros_gene$class_eff <- truth$class_eff[match(zeros_gene$gene, truth$gene)]
  p_propzeros_pergene <- ggplot2::ggplot(zeros_gene,
                                         ggplot2::aes(x = .data$propzeros, fill = .data$class_eff)) +
    ggplot2::geom_histogram(bins = 30) +
    no_legend +
    ggplot2::xlab("Proportion of cells with zero counts") + ggplot2::ylab("Number of genes") +
    ggplot2::scale_fill_manual(values = model_names, name = "Model") +
    ggplot2::coord_cartesian(xlim = c(0, 1))

  zeros_cell <- data.frame(cell = colnames(counts), propzeros = colMeans(counts == 0))
  p_propzeros_percell <- ggplot2::ggplot(zeros_cell, ggplot2::aes(x = .data$propzeros)) +
    ggplot2::geom_histogram(bins = 30) +
    blank +
    ggplot2::xlab("Proportion of genes with zero counts") + ggplot2::ylab("Number of cells") +
    ggplot2::coord_cartesian(xlim = c(0, 1))

  # realised counts per cell and per gene
  cell_df <- data.frame(cell = names(sim$cell_total), totalcounts = sim$cell_total)
  p_countspercell <- ggplot2::ggplot(cell_df, ggplot2::aes(x = .data$totalcounts)) +
    ggplot2::geom_histogram(bins = 100, fill = "lightgray") +
    ggplot2::xlab("Counts per cell") + ggplot2::ylab("Number of cells") +
    blank

  gene_total_df <- data.frame(gene = names(sim$gene_total), totalcounts = sim$gene_total)
  gene_total_df$class_eff <- truth$class_eff[match(gene_total_df$gene, truth$gene)]

  below_q <- gene_total_df$totalcounts <= stats::quantile(gene_total_df$totalcounts, 0.995, na.rm = TRUE)
  p_countspergene <- ggplot2::ggplot(gene_total_df[below_q, ],
                                     ggplot2::aes(x = .data$totalcounts, fill = .data$class_eff)) +
    ggplot2::geom_histogram(bins = 60) +
    ggplot2::xlab("Counts per gene (up to the 99.5th percentile)") + ggplot2::ylab("Number of genes") +
    no_legend +
    ggplot2::scale_fill_manual(values = model_names, name = "Model")

  p_lowcounts <- ggplot2::ggplot(gene_total_df[gene_total_df$totalcounts < 20, ],
                                 ggplot2::aes(x = .data$totalcounts, fill = .data$class_eff)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::xlab("Realised gene size (genes with fewer than 20 counts)") +
    ggplot2::ylab("Number of genes") +
    top_legend +
    ggplot2::scale_fill_manual(values = model_names, name = "Model")

  # counts matrix (evenly spaced subset of genes for speed)
  genes_show <- if (nrow(counts) > max_genes_heatmap) {
    unique(round(seq(1, nrow(counts), length.out = max_genes_heatmap)))
  } else {
    seq_len(nrow(counts))
  }
  hm <- as.data.frame(as.table(pmin(counts[genes_show, , drop = FALSE], 100)))
  names(hm) <- c("gene", "cell", "counts")
  p_countsmatrix <- ggplot2::ggplot(hm, ggplot2::aes(x = .data$cell, y = .data$gene, fill = .data$counts)) +
    ggplot2::geom_raster() +
    viridis::scale_fill_viridis(option = "A", trans = "sqrt", name = "Counts (capped at 100)") +
    ggplot2::xlab("Cells") +
    ggplot2::ylab(sprintf("Genes (subset of %d)", length(genes_show))) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "top")

  # simulated vs realised size factors
  cell_size_df <- data.frame(cell_size_sim = sim$cell_size,
                             cell_size_obs = sim$cell_size_obs[names(sim$cell_size)])
  p_cellsize <- ggplot2::ggplot(cell_size_df,
                                ggplot2::aes(x = .data$cell_size_sim, y = .data$cell_size_obs)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggpubr::stat_cor() +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::xlab("Simulated cell size") + ggplot2::ylab("Realised cell size") +
    blank

  gene_size_sim <- exp(truth$log_g)
  gene_size_df <- data.frame(
    gene_size_sim = gene_size_sim / mean(gene_size_sim),
    gene_size_obs = (sim$gene_total / mean(sim$gene_total))[truth$gene]
  )
  p_genesize <- ggplot2::ggplot(gene_size_df,
                                ggplot2::aes(x = .data$gene_size_sim, y = .data$gene_size_obs)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggpubr::stat_cor() +
    ggplot2::xlab("Simulated gene size") + ggplot2::ylab("Realised gene size") +
    blank

  return(c(plot_sim_params(truth),
           list(genesize = p_genesize, cellsize = p_cellsize,
              countsmatrix = p_countsmatrix, countspergene = p_countspergene,
              countspercell = p_countspercell, lowcounts = p_lowcounts,
              zeros_gene = p_propzeros_pergene, zeros_cell = p_propzeros_percell)))
}

#' Combined figure of the counts plots
#'
#' Arranges the counts matrix, counts per gene and per cell, and the
#' proportions of zeros from \code{\link{plot_sim_results}} in one figure.
#'
#' @param plots the list returned by \code{\link{plot_sim_results}}
#'
#' @return a ggplot object (arranged figure)
#' @export
plot_sim_counts_panel <- function(plots) {
  .require_pkgs("ggpubr")
  ggpubr::ggarrange(
    ggpubr::ggarrange(plots$countsmatrix, plots$countspergene, plots$countspercell, nrow = 1),
    ggpubr::ggarrange(plots$zeros_gene, plots$zeros_cell, nrow = 1),
    ncol = 1
  )
}

#' Combined figure of the gene and cell size plots
#'
#' Arranges the simulated versus realised gene and cell sizes from
#' \code{\link{plot_sim_results}} side by side.
#'
#' @inheritParams plot_sim_counts_panel
#'
#' @return a ggplot object (arranged figure)
#' @export
plot_sim_sizes_panel <- function(plots) {
  .require_pkgs("ggpubr")
  ggpubr::ggarrange(plots$genesize, plots$cellsize, nrow = 1)
}


# ---- trees with the data of the layers (mutations, methylation, barcodes) ----

# colours of the heatmaps
.nucleotide_cols <- c(A = "#009E73", C = "#0072B2", G = "#F0E442", T = "#D55E00")
# the 16 genotypes: the hue is the nucleotide of allele 1, the shade (dark to light) the nucleotide of allele 2
.genotype_cols <- local({
  shade <- function(col, amount) {  # amount < 0 darkens, > 0 lightens
    rgb <- grDevices::col2rgb(col)[, 1] / 255
    target <- if (amount < 0) 0 else 1
    out <- rgb + (target - rgb) * abs(amount)
    grDevices::rgb(out[1], out[2], out[3])
  }
  amounts <- c(A = -0.45, C = -0.12, G = 0.25, T = 0.55)
  cols <- unlist(lapply(names(.nucleotide_cols), function(a1) {
    vapply(amounts, function(a) shade(.nucleotide_cols[[a1]], a), "")
  }))
  names(cols) <- paste0(rep(names(.nucleotide_cols), each = 4), rep(names(.nucleotide_cols), 4))
  cols
})
.error_shapes <- c("one allele" = 16, "het -> hom" = 1, "both alleles" = 4)
.binary_cols <- c("0" = "white", "1" = "black", "false positive" = "#D55E00", "false negative" = "#009E73")
.unedited_col <- "#F2F2F2"
.missing_col <- "#CC79A7"

# Tree for the plots of the layers, coloured by type, with the root edge. The time tree has a time axis
# (time before present, 0 = present) and marks for the origin and the present. The distance tree
# (`distance = TRUE`, needs the `distance` column of clock_tree()) has branch lengths in expected
# substitutions per site, from the origin.
.plot_layer_tree <- function(tree, distance = FALSE, pad = 0.05) {
  stopifnot(methods::is(tree, "treedata"))
  if (distance) {
    tree <- distance_tree(tree)
  }
  phy <- tree@phylo
  has_type <- "type" %in% names(tree@data)
  types <- levels(factor(tree@data$type))
  cols <- .type_colours(types)

  p <- if (has_type) {
    ggtree::ggtree(tree, ggplot2::aes(color = factor(.data$type)))
  } else {
    ggtree::ggtree(tree, color = "black")
  }
  root_type <- if (has_type) as.character(tree@data$type[tree@data$node == ape::Ntip(phy) + 1]) else NULL
  root_col <- if (has_type) cols[[root_type]] else "black"

  if (distance) {
    # ggtree puts the root at 0: shift so that distances are measured from the origin
    if (!is.null(phy$root.edge)) {
      p$data$x <- p$data$x + phy$root.edge
    }
    if (!is.null(phy$root.edge)) p <- p + ggtree::geom_rootedge(color = root_col)
    if (has_type) p <- p + ggplot2::scale_color_manual(values = cols, name = "type")
    return(p +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, pad))) +
      ggtree::theme_tree2() +
      ggplot2::labs(x = "expected substitutions per site from the origin") +
      ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(16, 8, 4, 4)))
  }

  p <- ggtree::revts(p)  # present = 0
  if (!is.null(phy$root.edge)) p <- p + ggtree::geom_rootedge(color = root_col)
  if (has_type) p <- p + ggplot2::scale_color_manual(values = cols, name = "type")
  origin <- if (!is.null(phy$origin)) phy$origin else max(ape::node.depth.edgelength(phy)) + phy$root.edge
  breaks <- pretty(c(-origin, 0))
  breaks <- breaks[breaks >= -origin & breaks <= 0]
  p +
    ggplot2::geom_vline(xintercept = c(-origin, 0), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::annotate("text", x = c(-origin, 0), y = Inf, label = c("origin", "present"),
                      vjust = -0.4, hjust = c(0, 1), size = 3, color = "grey50") +
    ggplot2::scale_x_continuous(breaks = breaks, labels = abs, expand = ggplot2::expansion(mult = c(0.02, pad))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggtree::theme_tree2() +
    ggplot2::labs(x = "time before present") +
    ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(16, 8, 4, 4))
}

# one character matrix (tips x sites, rownames = tip labels) per entry of a list of tip x site matrices
.name_tip_matrix <- function(tree, m) {
  m[] <- as.character(m)
  dimnames(m) <- list(tree@phylo$tip.label, paste0("site_", seq_len(ncol(m))))
  m
}

# several matrices side by side, with `gap` empty columns between them, as a data frame of factors
.combine_blocks <- function(mats, prefix, gap, levels) {
  blocks <- Map(function(m, k) {
    colnames(m) <- paste0(prefix, k, "_", colnames(m))  # unique column names
    if (k < length(mats)) {
      m <- cbind(m, matrix(NA_character_, nrow(m), gap,
                           dimnames = list(rownames(m), paste0("gap_", k, "_", seq_len(gap)))))
    }
    m
  }, mats, seq_along(mats))
  df <- as.data.frame(do.call(cbind, blocks))
  df[] <- lapply(df, factor, levels = levels)
  df
}

#' Tree with the SNV genotypes of the tips
#'
#' Plots a tree and, next to it, the diploid genotypes of the tips simulated by
#' \code{\link{sim_gt16_nt_diploid}} as a heatmap with one row per cell. Every site has two
#' tiles side by side (allele 1 left, allele 2 right) with one colour per nucleotide, and a white
#' line separates the sites. A site is heterozygous where the two tiles of a cell differ.
#'
#' @param tree the treedata object the SNVs were simulated on. If it has a \code{distance} column
#'   (see \code{\link{clock_tree}}), \code{distance = TRUE} draws the tree with the distances as
#'   branch lengths.
#' @param snv the output of \code{\link{sim_gt16_nt_diploid}}: a list with the data frames
#'   \code{true} and \code{observed}, in either format
#' @param errors if TRUE, the observed genotypes are drawn and the alleles that the sequencing
#'   error changed are marked, with one point shape per kind of error. If FALSE (default), the
#'   true genotypes are drawn.
#' @param distance if TRUE, draw the distance tree instead of the time tree
#' @param site_width width of a site as a fraction of the width of the tree
#' @param genotypes if TRUE, draw one tile per cell and site, coloured by the genotype (16 colours:
#'   the hue is the nucleotide of allele 1, the shade the nucleotide of allele 2), instead of two
#'   tiles of four colours. With \code{errors = TRUE}, the genotypes changed by the error are marked
#'   with a point.
#'
#' @return a ggplot object
#' @export
plot_sim_snv <- function(tree, snv, errors = FALSE, distance = FALSE, site_width = 0.015, genotypes = FALSE) {
  .require_pkgs(c("ggplot2", "ggtree"))
  stopifnot(is.list(snv), all(c("true", "observed") %in% names(snv)))

  shown <- .allele_matrices(tree, snv, if (errors) "observed" else "true")
  p <- .plot_layer_tree(tree, distance = distance, pad = 0.02)

  tips <- p$data[p$data$isTip, ]
  y <- tips$y[match(rownames(shown$allele1), tips$label)]
  tree_width <- diff(range(p$data$x))
  w <- site_width * tree_width
  x0 <- max(p$data$x) + 0.04 * tree_width
  n_sites <- ncol(shown$allele1)

  if (genotypes) return(.plot_snv_genotypes(p, shown, tree, snv, errors, y, x0, w))

  # long format: the two alleles of a site side by side (allele 1 left, allele 2 right), each half a site wide
  left <- x0 + w / 2
  cells <- function(m, dx) {
    data.frame(x = left + (rep(seq_len(n_sites), each = nrow(m)) - 1) * w + dx * w, y = rep(y, n_sites),
               value = as.vector(m))
  }
  dx <- c(allele1 = 0.25, allele2 = 0.75)
  tiles <- rbind(cells(shown$allele1, dx[["allele1"]]), cells(shown$allele2, dx[["allele2"]]))
  tiles$value <- factor(tiles$value, levels = names(.nucleotide_cols))
  top <- max(y) + 0.5
  bottom <- min(y) - 0.5

  p <- p +
    ggplot2::geom_tile(data = tiles, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       width = w / 2, height = 0.9, inherit.aes = FALSE) +
    ggplot2::annotate("segment", x = left + seq_len(n_sites - 1) * w, xend = left + seq_len(n_sites - 1) * w,
                      y = bottom, yend = top, color = "black", linewidth = 0.6) +  # between sites
    ggplot2::annotate("rect", xmin = left, xmax = left + n_sites * w, ymin = bottom, ymax = top,
                      fill = NA, color = "grey40", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = .nucleotide_cols, name = "nucleotide", drop = FALSE,
                               guide = ggplot2::guide_legend(override.aes = list(shape = NA)))

  if (errors) {
    kinds <- .error_kinds(.allele_matrices(tree, snv, "true"), shown)
    marks <- rbind(cells(kinds$allele1, dx[["allele1"]]), cells(kinds$allele2, dx[["allele2"]]))
    marks <- marks[!is.na(marks$value), ]
    marks$value <- factor(marks$value, levels = names(.error_shapes))
    p <- p +
      ggplot2::geom_point(data = marks, ggplot2::aes(x = .data$x, y = .data$y, shape = .data$value),
                          size = 1.1, color = "black", inherit.aes = FALSE, show.legend = TRUE) +
      ggplot2::scale_shape_manual(values = .error_shapes, name = "error", limits = names(.error_shapes),
                                  drop = FALSE)
  }
  return(p)
}

# plot_sim_snv() with one tile of 16 colours per cell and site
.plot_snv_genotypes <- function(p, shown, tree, snv, errors, y, x0, w) {
  n_sites <- ncol(shown$allele1)
  long <- function(m) data.frame(x = x0 + rep(seq_len(n_sites), each = nrow(m)) * w, y = rep(y, n_sites),
                                 value = as.vector(m))
  tiles <- long(matrix(paste0(shown$allele1, shown$allele2), nrow(shown$allele1)))
  tiles$value <- factor(tiles$value, levels = names(.genotype_cols))
  p <- p +
    ggplot2::geom_tile(data = tiles, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       width = w, height = 0.9, inherit.aes = FALSE) +
    ggplot2::annotate("rect", xmin = x0 + w / 2, xmax = x0 + (n_sites + 0.5) * w,
                      ymin = min(y) - 0.5, ymax = max(y) + 0.5, fill = NA, color = "grey40", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = .genotype_cols, name = "genotype", drop = FALSE,
                               guide = ggplot2::guide_legend(nrow = 2, override.aes = list(shape = NA)))
  if (errors) {
    true <- .allele_matrices(tree, snv, "true")
    changed <- true$allele1 != shown$allele1 | true$allele2 != shown$allele2
    marks <- long(changed)
    marks <- marks[marks$value == "TRUE", ]
    p <- p +
      ggplot2::geom_point(data = marks, ggplot2::aes(x = .data$x, y = .data$y), shape = 16, size = 1.1,
                          color = "black", inherit.aes = FALSE) +
      ggplot2::labs(caption = "points: genotypes changed by the error")
  }
  p
}

# Genotypes of the tips as two tips x sites character matrices (allele 1, allele 2), from the `true` or
# `observed` part of the output of sim_gt16_nt_diploid(), in either format
.allele_matrices <- function(tree, snv, which = c("true", "observed")) {
  df <- snv[[match.arg(which)]]
  df <- df[match(seq_along(tree@phylo$tip.label), df$node), , drop = FALSE]
  chars <- function(x) do.call(rbind, strsplit(x, "", fixed = TRUE))
  m <- if ("genotype" %in% names(df)) {
    g <- do.call(rbind, strsplit(df$genotype, "_", fixed = TRUE))
    list(allele1 = substr(g, 1, 1), allele2 = substr(g, 2, 2))
  } else {
    list(allele1 = chars(df$allele_1), allele2 = chars(df$allele_2))
  }
  lapply(m, function(a) .name_tip_matrix(tree, matrix(a, nrow = nrow(df))))
}

# Kind of error of every allele: NA (unchanged), "one allele" (one allele read as another nucleotide),
# "het -> hom" (a heterozygous site read as homozygous, as after allelic dropout) or "both alleles"
.error_kinds <- function(true, observed) {
  changed1 <- true$allele1 != observed$allele1
  changed2 <- true$allele2 != observed$allele2
  het_to_hom <- true$allele1 != true$allele2 & observed$allele1 == observed$allele2
  site_kind <- ifelse(changed1 & changed2, "both alleles",
                      ifelse(het_to_hom, "het -> hom", "one allele"))
  list(allele1 = ifelse(changed1, site_kind, NA_character_),
       allele2 = ifelse(changed2, site_kind, NA_character_))
}

#' Tree with the binary sequences of the tips
#'
#' Plots a tree and, next to it, the binary sequences (e.g. methylation) of the tips simulated by
#' \code{\link{sim_binary}} as a black (1) and white (0) heatmap with one row per tip. Several
#' sequences are drawn as blocks side by side.
#'
#' @inheritParams plot_sim_snv
#' @param bin the output of \code{\link{sim_binary}} with \code{write_as_string = TRUE}: a data
#'   frame with the columns \code{node}, \code{binary_1}, ...
#' @param errors optional data frame with the true states, in the format of \code{bin}, e.g. the
#'   output of \code{\link{sim_binary}} with the same seed and without error. The sites that the
#'   error flipped are marked as false positives and false negatives.
#' @param gap number of empty columns between the sequences
#' @param site_width width of a site as a fraction of the width of the tree
#'
#' @return a ggplot object
#' @export
plot_sim_binary <- function(tree, bin, errors = NULL, distance = FALSE, gap = 2, site_width = 0.006) {
  .require_pkgs(c("ggplot2", "ggtree"))
  seq_cols <- grep("^binary_", names(bin), value = TRUE)
  if (length(seq_cols) == 0) {
    stop("`bin` must be the data frame returned by sim_binary() with `write_as_string = TRUE`.",
         call. = FALSE)
  }
  tip_matrices <- function(df) {
    df <- df[match(seq_along(tree@phylo$tip.label), df$node), , drop = FALSE]
    lapply(seq_cols, function(col) .name_tip_matrix(tree, do.call(rbind, strsplit(df[[col]], "", fixed = TRUE))))
  }
  mats <- tip_matrices(bin)
  if (!is.null(errors)) {
    truth <- tip_matrices(errors)
    mats <- Map(function(m, t) {
      m[t == "0" & m == "1"] <- "false positive"
      m[t == "1" & m == "0"] <- "false negative"
      m
    }, mats, truth)
  }
  n_cols <- sum(vapply(mats, ncol, 1L)) + gap * (length(mats) - 1)
  df <- .combine_blocks(mats, "binary_", gap, names(.binary_cols))

  p <- .plot_layer_tree(tree, distance = distance, pad = 0.02)
  tree_width <- diff(range(p$data$x))
  offset <- 0.04 * tree_width
  width <- n_cols * site_width
  w <- width * tree_width / n_cols  # gheatmap places column j at max(x) + offset + j * w
  # frame around every block
  first <- c(1, 1 + cumsum(vapply(mats, ncol, 1L) + gap))[seq_along(mats)]
  last <- first + vapply(mats, ncol, 1L) - 1
  frames <- data.frame(xmin = max(p$data$x) + offset + (first - 0.5) * w,
                       xmax = max(p$data$x) + offset + (last + 0.5) * w,
                       ymin = 0.5, ymax = nrow(mats[[1]]) + 0.5)
  shown <- names(.binary_cols)[names(.binary_cols) %in% unlist(mats)]
  suppressMessages(
    ggtree::gheatmap(p, df, offset = offset, width = width, colnames = FALSE, color = NA) +
      ggplot2::geom_rect(data = frames,
                         ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
                         fill = NA, color = "grey40", linewidth = 0.3, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = .binary_cols, name = "state", na.translate = FALSE, breaks = shown,
                                 guide = ggplot2::guide_legend(override.aes = list(color = "grey40")))
  )
}

#' Tree with the lineage barcodes of the tips
#'
#' Plots a tree and, next to it, the barcodes of the tips as colour-coded heatmaps, one block per
#' barcode: unedited sites are light grey, every edit outcome has its own colour, and missing
#' sites (silenced or dropout) are pink.
#'
#' @inheritParams plot_sim_snv
#' @param bc the output of \code{\link{sim_barcode_seq}} or \code{\link{sim_barcode_generic}} with
#'   \code{write_as_string = TRUE} (columns \code{barcode_1}, ...), or of
#'   \code{\link{sim_barcode_nonseq}} (columns \code{site_1}, ...)
#' @param n_edits number of edit outcomes E (\code{length(edit_probs)}), so that the colours do
#'   not depend on which outcomes happen to occur
#' @param missing_state the value that marks silenced and dropout barcodes (sites)
#' @param gap number of empty columns between the barcodes
#' @param site_width width of a site as a fraction of the width of the tree
#'
#' @return a ggplot object
#' @export
plot_sim_barcodes <- function(tree, bc, n_edits, missing_state = "-", gap = 2, site_width = 0.025) {
  .require_pkgs(c("ggplot2", "ggtree"))
  mats <- .barcode_matrices(tree, bc, missing_state)
  cols <- c(stats::setNames(c(.unedited_col, grDevices::hcl.colors(n_edits, "Viridis", rev = TRUE)),
                            c("0", seq_len(n_edits))),
            missing = .missing_col)
  n_cols <- sum(vapply(mats, ncol, 1L)) + gap * (length(mats) - 1)
  df <- .combine_blocks(mats, "barcode_", gap, names(cols))

  p <- .plot_layer_tree(tree, pad = n_cols * site_width * 0.5 + 0.1)
  tree_width <- diff(range(p$data$x))
  # legend: the states that occur, with fixed colours and names
  labels <- stats::setNames(c("unedited", paste("edit", seq_len(n_edits)), "missing"), names(cols))
  shown <- names(cols)[names(cols) %in% unlist(mats)]
  suppressMessages(
    ggtree::gheatmap(p, df, offset = 0.1 * tree_width, width = n_cols * site_width,
                     colnames = FALSE, color = "white") +
      ggplot2::scale_fill_manual(values = cols, name = "state", na.translate = FALSE,
                                 breaks = shown, labels = labels[shown])
  )
}

# Barcode states of the tips as one matrix per barcode (tips x sites, rownames = tip labels); `bc` has
# one column per barcode (`barcode_1`, ...: strings with sites separated by "_") or one column per site
# (`site_1`, ...)
.barcode_matrices <- function(tree, bc, missing_state = "-") {
  bc <- bc[match(seq_along(tree@phylo$tip.label), bc$node), , drop = FALSE]
  cols <- grep("^barcode_", names(bc), value = TRUE)
  mats <- if (length(cols) > 0) {
    lapply(cols, function(col) do.call(rbind, strsplit(bc[[col]], "_", fixed = TRUE)))
  } else {
    list(as.matrix(bc[grep("^site_", names(bc))]))
  }
  lapply(mats, function(m) {
    m <- .name_tip_matrix(tree, m)
    m[m == as.character(missing_state)] <- "missing"
    m
  })
}
