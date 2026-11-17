# Plotting helpers shared by the vignettes (not part of the package).

ggplot2::set_theme(ggplot2::theme_classic())

# color palette
palette <- c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
# tree colors
type_cols <- setNames(palette, as.character(c(0:4)))
status_cols <- c("alive" = "#009E73", "dead" = "#999999")
# barcode colors
unedited_col <- "#F2F2F2"
missing_col <- "#CC79A7"

# type of the root node (= the origin type if the tree keeps its type changes)
root_type <- function(tree) {
  root <- length(tree@phylo$tip.label) + 1
  tree@data$type[tree@data$node == root]
}

# Tree with a time axis (time before present, 0 = present) and marks for the origin and the present.
#   color: "none", "type" (branches by type; root edge by the root type), "status"
#          (tips by alive / dead) or "distance" (branches by the `distance` column of a
#          clock_tree() result; `limits` fixes the colour scale, e.g. to share it between plots)
#   pad: extra space on the right (fraction of the plot), e.g. for the barcodes of `plot_barcodes()`
plot_tree <- function(tree, color = c("none", "type", "status", "distance"), pad = 0.05,
                      tip_size = 1, show_nodes = FALSE, limits = NULL) {
  color <- match.arg(color)
  origin <- tree@phylo$origin
  has_type <- "type" %in% names(tree@data)

  p <- if (color == "type" && has_type) {
    ggtree::ggtree(tree, ggplot2::aes(color = factor(type)))
  } else if (color == "distance") {
    ggtree::ggtree(tree, ggplot2::aes(color = distance), linewidth = 0.9)
  } else {
    ggtree::ggtree(tree, color = "black")
  }
  p <- ggtree::revts(p)  # present = 0

  if (color == "distance") {
    # the root edge as a segment, so that it takes its colour from the same scale
    if (!is.null(tree@phylo$root.edge)) {
      root <- p$data[p$data$node == p$data$parent, ]
      p <- p + ggplot2::geom_segment(
        data = data.frame(x = root$x - tree@phylo$root.edge, xend = root$x, y = root$y, distance = root$distance),
        ggplot2::aes(x = x, xend = xend, y = y, yend = y, color = distance), linewidth = 0.9, inherit.aes = FALSE
      )
    }
    p <- p + ggplot2::scale_color_viridis_c(name = "distance", limits = limits, end = 0.9)
  } else {
    root_col <- if (color == "type" && has_type) type_cols[[as.character(root_type(tree))]] else "black"
    if (!is.null(tree@phylo$root.edge)) p <- p + ggtree::geom_rootedge(color = root_col)
  }

  if (color == "type" && has_type) {
    p <- p + ggplot2::scale_color_manual(values = type_cols, name = "type")
  } else if (color == "status") {
    p <- p +
      ggtree::geom_tippoint(ggplot2::aes(subset = isTip & status == 1), color = status_cols[["alive"]], size = tip_size) +
      ggtree::geom_tippoint(ggplot2::aes(subset = isTip & status == 0), color = status_cols[["dead"]], size = tip_size)
  }
  
  if (show_nodes) {
    p <- p + ggtree::geom_nodepoint(size = tip_size / 2)
  } 

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
    ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(16, 8, 4, 4),
                   plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 14)))
}

# Barcode states of the tips as one matrix per barcode (tips x sites, rownames = tip labels).
# `bc` has one column per barcode (`barcode_1`, ..., strings with sites separated by "_";
# sim_barcode_seq / _generic) or one column per site (`site_1`, ...; sim_barcode_nonseq).
barcode_matrices <- function(tree, bc, missing_state = "-") {
  tips <- seq_along(tree@phylo$tip.label)
  bc <- bc[match(tips, bc$node), , drop = FALSE]
  cols <- grep("^barcode_", names(bc), value = TRUE)
  mats <- if (length(cols) > 0) {
    lapply(cols, function(col) do.call(rbind, strsplit(bc[[col]], "_", fixed = TRUE)))
  } else {
    list(as.matrix(bc[grep("^site_", names(bc))]))
  }
  lapply(mats, function(m) {
    m[] <- as.character(m)
    m[m == as.character(missing_state)] <- "missing"
    dimnames(m) <- list(tree@phylo$tip.label, paste0("site_", seq_len(ncol(m))))
    m
  })
}

# colours of the barcode states: unedited, edit outcomes 1, ..., E, missing
state_cols <- function(n_edits) {
  c(stats::setNames(c(unedited_col, grDevices::hcl.colors(n_edits, "Viridis", rev = TRUE)), c("0", seq_len(n_edits))),
    missing = missing_col)
}

# Tree with the barcodes of the tips as colour-coded heatmaps, one block per barcode:
# unedited sites are light grey, edit outcomes have their own colour, missing sites (silenced
# or dropout) are pink. `n_edits` = number of edit outcomes E (`length(edit_probs)`), so that
# the colours do not depend on which outcomes happen to occur. `gap` = number of empty columns
# between the barcodes, `site_width` = width of a column.
plot_barcodes <- function(tree, bc, n_edits, color = "type", missing_state = "-", site_width = 0.025, gap = 2, ...) {
  mats <- barcode_matrices(tree, bc, missing_state)
  cols <- state_cols(n_edits)
  n_cols <- sum(vapply(mats, ncol, 1L)) + gap * (length(mats) - 1)
  # one heatmap with an empty column between the barcodes (a single call: one fill scale)
  blocks <- Map(function(m, k) {
    colnames(m) <- paste0("barcode_", k, "_", colnames(m))  # unique column names
    if (k < length(mats)) m <- cbind(m, matrix(NA_character_, nrow(m), gap, dimnames = list(rownames(m), paste0("gap_", k, "_", seq_len(gap)))))
    m
  }, mats, seq_along(mats))
  df <- as.data.frame(do.call(cbind, blocks))
  df[] <- lapply(df, factor, levels = names(cols))
  p <- plot_tree(tree, color = color, pad = n_cols * site_width * 0.5 + 0.1, ...)
  tree_width <- diff(range(p$data$x))
  # legend: the states that occur, with fixed colours and names
  labels <- stats::setNames(c("unedited", paste("edit", seq_len(n_edits)), "missing"), names(cols))
  shown <- names(cols)[names(cols) %in% unlist(mats)]
  # gheatmap adds its own y and fill scales, which the fill scale below replaces (messages)
  suppressMessages(
    ggtree::gheatmap(p, df, offset = 0.1 * tree_width, width = n_cols * site_width,
                     colnames = FALSE, color = "white") +
      ggplot2::scale_fill_manual(values = cols, name = "state", na.translate = FALSE,
                                 breaks = shown, labels = labels[shown])
  )
}


# ---- binary sequences ----

# error colours differ from the type colours (palette[1:2])
binary_cols <- c("0" = "white", "1" = "black", "false positive" = palette[[4]], "false negative" = palette[[3]])

# Tree with branch lengths in expected substitutions per site (distance_tree() of a clock_tree() result),
# coloured by type as in plot_tree(), with the origin at 0. No present mark: the tips are not aligned
# in distance.
plot_distance_tree <- function(tree, color = c("type", "none"), pad = 0.05) {
  color <- match.arg(color)
  dt <- distance_tree(tree)
  by_type <- color == "type" && "type" %in% names(tree@data)
  p <- if (by_type) ggtree::ggtree(dt, ggplot2::aes(color = factor(type))) else ggtree::ggtree(dt, color = "black")
  # ggtree puts the root at 0; shift so that distances are measured from the origin
  if (!is.null(dt@phylo$root.edge)) p$data$x <- p$data$x + dt@phylo$root.edge
  root_col <- if (by_type) type_cols[[as.character(root_type(tree))]] else "black"
  if (!is.null(dt@phylo$root.edge)) p <- p + ggtree::geom_rootedge(color = root_col)
  if (by_type) p <- p + ggplot2::scale_color_manual(values = type_cols, name = "type")
  p +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, pad))) +
    ggtree::theme_tree2() +
    ggplot2::labs(x = "expected substitutions per site from the origin") +
    ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(16, 8, 4, 4))
}

# Binary states of the tips as a tips x sites character matrix (rownames = tip labels), from the
# output of sim_binary(): the data frame of strings (column `binary_1`) or one data frame of the
# `write_as_string = FALSE` list (columns `site_1`, ...).
binary_matrix <- function(tree, bin) {
  tips <- seq_along(tree@phylo$tip.label)
  bin <- bin[match(tips, bin$node), , drop = FALSE]
  m <- if ("binary_1" %in% names(bin)) {
    do.call(rbind, strsplit(bin$binary_1, "", fixed = TRUE))
  } else {
    as.matrix(bin[grep("^site_", names(bin))])
  }
  m[] <- as.character(m)
  dimnames(m) <- list(tree@phylo$tip.label, paste0("site_", seq_len(ncol(m))))
  m
}

# Tree with the binary sequences of the tips as a black (1) / white (0) heatmap, one row per tip.
# `errors`: the true states (sim_binary() output without error), to mark the sites that the error
# flipped. `distance = TRUE` draws the distance tree instead of the time tree.
plot_binary <- function(tree, bin, errors = NULL, distance = FALSE, color = "type", site_width = 0.006) {
  m <- binary_matrix(tree, bin)
  if (!is.null(errors)) {
    truth <- binary_matrix(tree, errors)
    m[truth == "0" & m == "1"] <- "false positive"
    m[truth == "1" & m == "0"] <- "false negative"
  }
  df <- as.data.frame(m)
  df[] <- lapply(df, factor, levels = names(binary_cols))
  p <- if (distance) plot_distance_tree(tree, color = color, pad = 0.02) else plot_tree(tree, color = color, pad = 0.02)
  tree_width <- diff(range(p$data$x))
  offset <- 0.04 * tree_width
  width <- ncol(m) * site_width
  # frame around the heatmap (gheatmap places column j at max(x) + offset + j * w, w = width * tree_width / ncol)
  w <- width * tree_width / ncol(m)
  frame <- c(xmin = max(p$data$x) + offset + w / 2, xmax = max(p$data$x) + offset + (ncol(m) + 0.5) * w,
             ymin = 0.5, ymax = nrow(m) + 0.5)
  shown <- names(binary_cols)[names(binary_cols) %in% m]
  suppressMessages(
    ggtree::gheatmap(p, df, offset = offset, width = width, colnames = FALSE, color = NA) +
      ggplot2::annotate("rect", xmin = frame[["xmin"]], xmax = frame[["xmax"]], ymin = frame[["ymin"]],
                        ymax = frame[["ymax"]], fill = NA, color = "grey40", linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = binary_cols, name = "state", na.translate = FALSE, breaks = shown,
                                 guide = ggplot2::guide_legend(override.aes = list(color = "grey40")))
  )
}


# ---- diploid SNVs (GT16) ----

# one colour per nucleotide (Okabe-Ito, distinct from the type colours palette[1:2])
nucleotide_cols <- c(A = "#009E73", C = "#0072B2", G = "#F0E442", T = "#D55E00")
# point shapes of the error kinds (all black)
error_shapes <- c("one allele" = 16, "het -> hom" = 1, "both alleles" = 4)

# Genotypes of the tips as two tips x sites character matrices (allele 1, allele 2; rownames = tip
# labels), from the `true` or `observed` part of the output of sim_gt16_nt_diploid(), in either
# format ("genotype": "AA_CG_...", or "alleles": `allele_1`, `allele_2`).
allele_matrices <- function(tree, snv, which = c("true", "observed")) {
  df <- snv[[match.arg(which)]]
  df <- df[match(seq_along(tree@phylo$tip.label), df$node), , drop = FALSE]
  chars <- function(x) do.call(rbind, strsplit(x, "", fixed = TRUE))
  m <- if ("genotype" %in% names(df)) {
    g <- do.call(rbind, strsplit(df$genotype, "_", fixed = TRUE))
    list(allele1 = substr(g, 1, 1), allele2 = substr(g, 2, 2))
  } else {
    list(allele1 = chars(df$allele_1), allele2 = chars(df$allele_2))
  }
  lapply(m, function(a) {
    a <- matrix(a, nrow = nrow(df))
    dimnames(a) <- list(tree@phylo$tip.label, paste0("site_", seq_len(ncol(a))))
    a
  })
}

# Kind of error of every allele cell, from the true and observed allele_matrices():
# NA (unchanged), "one allele" (one allele read as another nucleotide), "het -> hom" (a heterozygous
# site read as homozygous for one of its alleles, as after allelic dropout) or "both alleles".
error_kinds <- function(true, observed) {
  changed1 <- true$allele1 != observed$allele1
  changed2 <- true$allele2 != observed$allele2
  het_to_hom <- true$allele1 != true$allele2 & observed$allele1 == observed$allele2
  site_kind <- ifelse(changed1 & changed2, "both alleles",
                      ifelse(het_to_hom, "het -> hom", "one allele"))
  list(allele1 = ifelse(changed1, site_kind, NA_character_),
       allele2 = ifelse(changed2, site_kind, NA_character_))
}

# Tree with the diploid genotypes of the tips as a heatmap with two rows per tip (allele 1 above,
# allele 2 below), one colour per nucleotide, from the output of sim_gt16_nt_diploid(). By default
# the true genotypes are drawn; `errors = TRUE` draws the observed genotypes and marks the alleles
# that the error changed, with one point shape per kind of error. `distance = TRUE` draws the
# distance tree.
plot_alleles <- function(tree, snv, errors = FALSE, distance = FALSE, color = "type", site_width = 0.015) {
  shown <- allele_matrices(tree, snv, if (errors) "observed" else "true")
  p <- if (distance) plot_distance_tree(tree, color = color, pad = 0.02) else plot_tree(tree, color = color, pad = 0.02)

  tips <- p$data[p$data$isTip, ]
  y <- tips$y[match(rownames(shown$allele1), tips$label)]
  tree_width <- diff(range(p$data$x))
  w <- site_width * tree_width
  x0 <- max(p$data$x) + 0.04 * tree_width
  n_sites <- ncol(shown$allele1)

  # long format: allele 1 slightly above the tip, allele 2 slightly below
  cells <- function(m, dy) {
    data.frame(x = x0 + rep(seq_len(n_sites), each = nrow(m)) * w, y = rep(y, n_sites) + dy,
               value = as.vector(m))
  }
  dy <- c(allele1 = 0.22, allele2 = -0.22)
  tiles <- rbind(cells(shown$allele1, dy[["allele1"]]), cells(shown$allele2, dy[["allele2"]]))
  tiles$value <- factor(tiles$value, levels = names(nucleotide_cols))

  p <- p +
    ggplot2::geom_tile(data = tiles, ggplot2::aes(x = x, y = y, fill = value), width = w, height = 0.42,
                       inherit.aes = FALSE) +
    ggplot2::annotate("rect", xmin = x0 + w / 2, xmax = x0 + (n_sites + 0.5) * w,
                      ymin = min(y) - 0.5, ymax = max(y) + 0.5, fill = NA, color = "grey40", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = nucleotide_cols, name = "nucleotide", drop = FALSE,
                               guide = ggplot2::guide_legend(override.aes = list(shape = NA)))

  if (errors) {
    kinds <- error_kinds(allele_matrices(tree, snv, "true"), shown)
    marks <- rbind(cells(kinds$allele1, dy[["allele1"]]), cells(kinds$allele2, dy[["allele2"]]))
    marks <- marks[!is.na(marks$value), ]
    marks$value <- factor(marks$value, levels = names(error_shapes))
    p <- p +
      ggplot2::geom_point(data = marks, ggplot2::aes(x = x, y = y, shape = value), size = 1.1,
                          color = "black", inherit.aes = FALSE, show.legend = TRUE) +
      # all kinds in the legend, so that the legends of several plots can be collected
      ggplot2::scale_shape_manual(values = error_shapes, name = "error", limits = names(error_shapes),
                                  drop = FALSE)
  }
  p
}
