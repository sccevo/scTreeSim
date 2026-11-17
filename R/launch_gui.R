#' Launch the scTreeSim Shiny app
#'
#' Opens a graphical interface to simulate a tree and, on top of it, the layers
#' of data of this package, using its simulation functions: SNV mutations
#' (\code{\link{sim_gt16_nt_diploid}}), binary methylation
#' (\code{\link{sim_binary}}), lineage barcodes
#' (\code{\link{sim_barcode_seq}}, \code{\link{sim_barcode_nonseq}} and
#' \code{\link{sim_barcode_generic}}), the expression of many genes
#' (\code{\link{sim_expression_matrix}}) and the single-cell counts at the tips
#' (\code{\link{sim_expression_counts_matrix}}). Every layer is a tab, in
#' this order, and uses the tree of the first tab. The tree has any number of types
#' and is simulated up to a number of tips (\code{\link{sim_adb_ntaxa_samp}}) or up
#' to a time (\code{\link{sim_adb_origin_samp}}). The parameters of the genes are
#' drawn from priors, or read from a table with one row per gene.
#'
#' The button "Reproduce in R" downloads an R script with the code and the seeds
#' of the simulation that is shown, which uses the functions of this package. The
#' button "Download all" saves every result, the true parameters of the genes and
#' the script in one \code{.rds} file; the tabs of the expression and the counts
#' also download their true parameters as tables.
#'
#' The app needs the suggested packages \pkg{shiny}, \pkg{bslib},
#' \pkg{ggplot2}, \pkg{ggtree}, \pkg{ggpubr} and \pkg{viridis}.
#'
#' @param ... further arguments passed to \code{shiny::runApp()}, e.g.
#'   \code{port} or \code{launch.browser}
#'
#' @return Called for its side effect of running the app; returns the value of
#'   \code{shiny::runApp()} invisibly.
#' @export
launch_gui <- function(...) {
  .require_pkgs(c("shiny", "bslib", "ggplot2", "ggtree", "ggpubr", "viridis"))

  app_dir <- system.file("app", package = "scTreeSim")
  if (app_dir == "") {
    stop("The Shiny app was not found in the installed package; try reinstalling scTreeSim.",
         call. = FALSE)
  }

  invisible(shiny::runApp(app_dir, ...))
}
