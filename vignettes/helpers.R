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
#   color: "none", "type" (branches by type; root edge by the root type) or "status"
#          (tips by alive / dead)
#   pad: extra space on the right (fraction of the plot), e.g. for the barcodes of `plot_barcodes()`
plot_tree <- function(tree, color = c("none", "type", "status"), pad = 0.05,
                      tip_size = 1, show_nodes = FALSE) {
  color <- match.arg(color)
  origin <- tree@phylo$origin
  has_type <- "type" %in% names(tree@data)

  p <- if (color == "type" && has_type) {
    ggtree::ggtree(tree, ggplot2::aes(color = factor(type)))
  } else {
    ggtree::ggtree(tree, color = "black")
  }
  p <- ggtree::revts(p)  # present = 0

  root_col <- if (color == "type" && has_type) type_cols[[as.character(root_type(tree))]] else "black"
  if (!is.null(tree@phylo$root.edge)) p <- p + ggtree::geom_rootedge(color = root_col)

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
    ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(16, 8, 4, 4))
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
