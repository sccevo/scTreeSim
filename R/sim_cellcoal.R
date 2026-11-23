# CellCoal on a tree of scTreeSim: the tree is the genealogy of the cells (CellCoal -T) and CellCoal simulates
# the mutations and the single-cell experiment on it. CellCoal reads a branch length as the expected number
# of mutations per site of each of the two haplotypes of a cell, and ignores the root edge of the tree (it
# adds its own trunk, of -k times the mean depth of the tips, and the healthy outgroup).

# name of the tree file in the run folder
.cellcoal_tree_file <- "sctreesim_tree.nwk"

#' CellCoal on a tree
#'
#' Runs CellCoal (see [cellcoal_run()]) with `tree` as the genealogy of the cells: CellCoal simulates the
#' somatic mutations, the single-cell amplification and the sequencing on it, with the options in `args`
#' and `config`. CellCoal adds a healthy outgroup cell (`outgcell`), unless `--no-outgroup`.
#'
#' The branch lengths are the evolutionary distances of the clock model, in expected substitutions per
#' (diploid) site as in [sim_gt16_nt_diploid()]: a branch of distance \eqn{d} carries on average \eqn{d L}
#' mutations over \eqn{L} sites. CellCoal counts mutations per site of each haplotype, so the tree it reads
#' has branch lengths \eqn{d / 2}. The mutation rate of CellCoal (`-u`) has no effect on a given tree.
#'
#' CellCoal ignores the root edge of a tree and puts a trunk of its own between the healthy root and the
#' MRCA of the cells, of `-k` times the mean distance from the MRCA to the tips. With `root_edge = TRUE`,
#' `-k` is set so that the trunk is the root edge of `tree` (no trunk if the tree has no root edge); with
#' `root_edge = FALSE`, the trunk is that of CellCoal (`-k`, default 0.5). Every replicate (`-n`) uses the
#' same tree. Options of CellCoal for its own genealogies (population size, growth, demographic periods,
#' demes) do not apply.
#'
#' @param tree a treedata object, with the sampled cells only (no dead cells, `status` 0): e.g. from
#'   [sim_adb_ntaxa_samp()]. Its tip labels name the cells in CellCoal.
#' @param args,config,files,dir,timeout,echo,error options of CellCoal and of the run, see [cellcoal_run()].
#'   They must not give a tree (`-T`, `--tree file=`, `tree_file`), nor `-k` if `root_edge = TRUE`.
#' @param clock a clock model, e.g. [clock_strict()]; `NULL` = use the `distance` column of `tree` (see
#'   [clock_tree()])
#' @param root_edge if `TRUE` (default), the root edge of `tree` is the trunk of CellCoal
#'
#' @return a run of CellCoal (see [cellcoal_run()]; read its results with [cellcoal_read()]), with also
#'   * `tree`: `tree` with the `distance` of every branch (see [clock_tree()]),
#'   * `cells`: a data frame with columns `cell` (the name in CellCoal) and `node` (the tip in `tree`).
#' @seealso [clock_strict()], [clock_relaxed()], [spike()] for clock models, [clock_tree()] to store
#'   and reuse a clock draw
#' @family CellCoal interface
#' @export
sim_cellcoal <- function(tree, args = character(), config = NULL, clock = clock_strict(), root_edge = TRUE,
                         files = NULL, dir = NULL, timeout = Inf, echo = FALSE, error = TRUE) {
  stopifnot(methods::is(tree, "treedata"))
  stopifnot(is.null(clock) || inherits(clock, "clock_model"))
  stopifnot(is.logical(root_edge), length(root_edge) == 1, !is.na(root_edge))
  args <- .cellcoal_args(args)
  .cellcoal_check_tree_options(args, config, root_edge)

  phy <- tree@phylo
  n_tips <- ape::Ntip(phy)
  if (n_tips < 2) stop("`tree` must have at least 2 cells.", call. = FALSE)
  if ("status" %in% names(tree@data) && any(tree@data$status[tree@data$node <= n_tips] == 0)) {
    stop("`tree` has dead cells (`status` 0). Simulate it with sampling, e.g. sim_adb_ntaxa_samp().", call. = FALSE)
  }
  .cellcoal_check_labels(phy$tip.label)

  # one clock draw for all replicates
  if (!is.null(clock)) tree <- clock_tree(tree, clock)
  d <- .node_distances(tree)

  # branch lengths per haplotype, without the root edge and the node labels
  phy$edge.length <- d[phy$edge[, 2]] / 2
  phy$root.edge <- NULL
  phy$node.label <- NULL
  tree_dir <- tempfile("cellcoal-tree-")
  dir.create(tree_dir)
  on.exit(unlink(tree_dir, recursive = TRUE), add = TRUE)
  tree_file <- file.path(tree_dir, .cellcoal_tree_file)
  ape::write.tree(phy, tree_file, digits = 15)

  tree_args <- paste0("-T", .cellcoal_tree_file)
  if (root_edge) {
    # the trunk of CellCoal is -k times the mean depth of the tips
    depth <- mean(ape::node.depth.edgelength(phy)[seq_len(n_tips)])
    trunk <- d[n_tips + 1L] / 2
    if (trunk > 0 && depth == 0) {
      stop("The tips of `tree` are at distance 0 from its root, so CellCoal cannot make a trunk of the root edge; ",
           "use `root_edge = FALSE`.", call. = FALSE)
    }
    tree_args <- c(tree_args, paste0("-k", if (trunk == 0) 0 else format(trunk / depth, digits = 15)))
  }

  run <- cellcoal_run(c(tree_args, args), config = config, files = c(files, tree_file), dir = dir,
                      timeout = timeout, echo = echo, error = error)
  run$tree <- tree
  run$cells <- data.frame(cell = phy$tip.label, node = seq_len(n_tips))
  run
}

# the tree and its trunk are set by sim_cellcoal(); a mutation rate has no effect on a given tree
.cellcoal_check_tree_options <- function(args, config, root_edge) {
  # the value of a group of options (--tree file=...,transform=...) is the argument after it
  group <- c("", utils::head(args, -1))
  config_keys <- if (is.list(config)) names(config) else character()
  if (any(grepl("^-T", args) | (group == "--tree" & grepl("(^|,)file=", args))) ||
      any(config_keys %in% c("tree_file", "T"))) {
    stop("The tree of CellCoal is `tree`: remove `-T` (or `--tree file=`, `tree_file`).", call. = FALSE)
  }
  if (root_edge && (any(grepl("^-k", args) | (group == "--tree" & grepl("(^|,)transform=", args))) ||
                    any(config_keys %in% c("transforming_branch_length", "k")))) {
    stop("With `root_edge = TRUE`, the trunk of CellCoal is the root edge of `tree`: remove `-k` ",
         "(or `--tree transform=`, `transforming_branch_length`), or use `root_edge = FALSE`.", call. = FALSE)
  }
  if (any(grepl("^-u", args) | (group == "--mutation" & grepl("(^|,)rate=", args))) ||
      any(config_keys %in% c("mutation_rate", "u"))) {
    warning("The mutation rate of CellCoal has no effect on a given tree: the clock model sets the distances.",
            call. = FALSE)
  }
}

# tip labels that CellCoal reads as cell names
.cellcoal_check_labels <- function(labels) {
  bad <- labels[grepl("[\\s(),:;'\\[\\]]", labels, perl = TRUE) | nchar(labels) > 119 | !nzchar(labels)]
  if (length(bad) > 0) {
    stop("CellCoal cannot read the tip labels ", paste0("\"", utils::head(bad, 3), "\"", collapse = ", "),
         ": use names of at most 119 characters without spaces or ( ) , : ; ' [ ].", call. = FALSE)
  }
  if (anyDuplicated(labels) > 0) stop("The tip labels of `tree` must be unique.", call. = FALSE)
  if ("outgcell" %in% labels) stop("`outgcell` is the outgroup of CellCoal; rename that tip.", call. = FALSE)
}
