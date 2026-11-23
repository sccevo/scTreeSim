# Readers of the results of CellCoal (see cellcoal.R). A results folder holds either all the replicates in
# one file, each in a block that starts with a line "[#k]" (snv_gen, snv_hap, ado_cells, vcf; the trees
# file has one line per replicate), or, with -v, one file per replicate in a folder of its own
# (snv_genotypes_dir/snv_gen.0001, ...). A replicate without variable sites has no block.

# folder of the files of each replicate (-v)
.cellcoal_replicate_dirs <- c(trees = "trees_dir", snv_gen = "snv_genotypes_dir", snv_hap = "snv_haplotypes_dir",
                              ado_cells = "ado_cells_dir", vcf = "vcf_dir")

.cellcoal_results_dir <- function(run) {
  results <- if (inherits(run, "cellcoal_run")) run$results else run
  if (!is.character(results) || length(results) != 1 || is.na(results) || !dir.exists(results)) {
    stop("`run` must be a run of cellcoal_run() that wrote results, or a results folder of CellCoal.", call. = FALSE)
  }
  results
}

# the file of the results, either of one replicate (-v) or of all, possibly compressed (-z); NA if none
.cellcoal_file <- function(results, file, replicate = 1) {
  own <- file.path(results, .cellcoal_replicate_dirs[[file]], sprintf("%s.%04d", file, as.integer(replicate)))
  candidates <- c(own, paste0(own, ".gz"), file.path(results, c(file, paste0(file, ".gz"))))
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) NA_character_ else found[1]
}

.cellcoal_has <- function(results, file) !is.na(.cellcoal_file(results, file))

# the lines of one replicate of a file of the results
.cellcoal_lines <- function(results, file, replicate) {
  stopifnot(is.numeric(replicate), length(replicate) == 1, replicate >= 1, replicate == round(replicate))
  path <- .cellcoal_file(results, file, replicate)
  if (is.na(path)) {
    stop("The results have no file `", file, "` (is it among the outputs of the run? see cellcoal_docs()).",
         call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  # a file of one replicate is all of it
  if (dirname(path) != results) return(lines)

  if (file == "trees") {
    trees <- lines[nzchar(trimws(lines))]
    if (replicate > length(trees)) stop("The results have ", length(trees), " replicates.", call. = FALSE)
    return(trees[replicate])
  }
  starts <- grep("^\\[#[0-9]+\\]$", trimws(lines))
  numbers <- as.integer(gsub("\\D", "", lines[starts]))
  # the number of replicates heads the file, except for ado_cells, which has every replicate
  n_replicates <- if (grepl("^[0-9]+$", trimws(lines[1]))) as.integer(lines[1]) else max(numbers, 0L)
  if (replicate > n_replicates) stop("The results have ", n_replicates, " replicates.", call. = FALSE)
  k <- match(replicate, numbers)
  if (is.na(k)) return(character())
  end <- if (k < length(starts)) starts[k + 1] - 1 else length(lines)
  if (end <= starts[k]) character() else lines[(starts[k] + 1):end]
}

# rows of a Phylip-like file (snv_gen, snv_hap): "rows sites", the positions of the sites, then one row per
# name with its values
.cellcoal_phylip <- function(lines) {
  lines <- trimws(lines[nzchar(trimws(lines))])
  if (length(lines) == 0) return(list(sites = integer(), names = character(), values = list()))
  dims <- as.integer(strsplit(lines[1], "\\s+")[[1]])
  sites <- as.integer(strsplit(lines[2], "\\s+")[[1]])
  rows <- strsplit(lines[2 + seq_len(dims[1])], "\\s+")
  if (length(sites) != dims[2] || anyNA(sites)) {
    stop("Could not read the positions of the ", dims[2], " sites of the results.", call. = FALSE)
  }
  list(sites = sites, names = vapply(rows, `[`, "", 1), values = lapply(rows, `[`, -1))
}

# one column site_<position> per site
.cellcoal_site_columns <- function(values, sites) {
  values <- matrix(values, ncol = length(sites))
  colnames(values) <- paste0("site_", sites)
  as.data.frame(values, stringsAsFactors = FALSE)
}


#' Read the results of CellCoal
#'
#' Read the files written by a run of [cellcoal_run()] (or any results folder of CellCoal), one replicate at
#' a time. Each reader needs its output among the options of the run: `-6` for the tree, `-1` for the
#' genotypes, `-2` for the haplotypes, read counts (`-C`) for the VCF and the allelic dropout of the cells.
#' `cellcoal_read()` reads all of them that the run wrote. The cells are named as in CellCoal (`cell0001`,
#' ..., and the healthy outgroup `outgcell`) or as the tips of a user tree.
#'
#' * `cellcoal_tree()`: the genealogy, with branch lengths in expected mutations per site.
#' * `cellcoal_genotypes()`: the genotypes of the variable sites (file `snv_gen`): with allelic dropout and
#'   genotyping error, or, with read counts, the maximum likelihood genotypes. `?` is an allele lost by
#'   dropout and `N` a missing genotype; with the binary alphabet, the alleles are `0` and `1`.
#' * `cellcoal_haplotypes()`: the maternal and paternal haplotypes of the variable sites (file `snv_hap`).
#' * `cellcoal_ado()`: the allelic dropout rate of every cell (file `ado_cells`).
#' * `cellcoal_vcf()`: fields of the VCF, e.g. `GT` (called genotype), `DP` (depth), `RC` (read counts of
#'   A, C, G, T), `NG` (genotype in the library, after dropout) and `TG` (true genotype); the header of the
#'   VCF and the manual describe all of them. Compressed VCFs (`-z`) are read; BCFs (`--bcf`) are not.
#'
#' A replicate without variable sites gives no rows (genotypes, haplotypes) or no sites (VCF).
#'
#' @param run a run of [cellcoal_run()], or the path to a results folder of CellCoal
#' @param replicate number of the replicate
#' @param write_as_string if `TRUE` (default), the genotypes (haplotypes) of a cell are one string; if
#'   `FALSE`, one column per site
#' @param fields FORMAT fields of the VCF to read
#'
#' @return
#' * `cellcoal_tree()`: a `phylo` object.
#' * `cellcoal_genotypes()`: a data frame with one row per cell and columns `cell` and `genotype`, the
#'   genotypes separated by "_" (e.g. "AC_TT_CG", as `sim_gt16_nt_diploid(format = "genotype")`), or, if
#'   `write_as_string = FALSE`, `cell` and `site_<position>` with one genotype per site. The positions of the
#'   sites are in `attr(x, "sites")`.
#' * `cellcoal_haplotypes()`: a data frame with one row per cell and columns `cell`, `allele_1` (maternal)
#'   and `allele_2` (paternal), one string each (as `sim_gt16_nt_diploid(format = "alleles")`), or, if
#'   `write_as_string = FALSE`, two rows per cell with columns `cell`, `allele` (1 maternal, 2 paternal) and
#'   `site_<position>`. The positions of the sites are in `attr(x, "sites")`.
#' * `cellcoal_ado()`: a data frame with columns `cell` and `ado_rate`.
#' * `cellcoal_vcf()`: a list with `sites`, a data frame with one row per site (`chrom`, `position`, `id`,
#'   `ref`, `alt`, `info`), and one data frame per field, with one row per cell and columns `cell` and
#'   `site_<position>`; numbers are numeric and missing values `NA`.
#' * `cellcoal_read()`: a list with the results the run wrote, among `tree`, `genotypes`, `haplotypes`,
#'   `ado` and `vcf` (with the default fields).
#' @family CellCoal interface
#' @export
cellcoal_read <- function(run, replicate = 1) {
  results <- .cellcoal_results_dir(run)
  readers <- list(tree = list("trees", cellcoal_tree), genotypes = list("snv_gen", cellcoal_genotypes),
                  haplotypes = list("snv_hap", cellcoal_haplotypes), ado = list("ado_cells", cellcoal_ado),
                  vcf = list("vcf", cellcoal_vcf))
  present <- vapply(readers, function(reader) .cellcoal_has(results, reader[[1]]), logical(1))
  lapply(readers[present], function(reader) reader[[2]](results, replicate))
}

#' @rdname cellcoal_read
#' @export
cellcoal_tree <- function(run, replicate = 1) {
  lines <- .cellcoal_lines(.cellcoal_results_dir(run), "trees", replicate)
  ape::read.tree(text = paste(lines[nzchar(trimws(lines))], collapse = ""))
}

#' @rdname cellcoal_read
#' @export
cellcoal_genotypes <- function(run, replicate = 1, write_as_string = TRUE) {
  stopifnot(is.logical(write_as_string), length(write_as_string) == 1, !is.na(write_as_string))
  snv <- .cellcoal_phylip(.cellcoal_lines(.cellcoal_results_dir(run), "snv_gen", replicate))

  genotypes <- if (write_as_string) {
    data.frame(cell = snv$names, genotype = vapply(snv$values, paste, "", collapse = "_"))
  } else {
    # one row of genotypes per cell
    values <- if (length(snv$values) > 0) do.call(rbind, snv$values) else character()
    data.frame(cell = snv$names, .cellcoal_site_columns(values, snv$sites))
  }
  attr(genotypes, "sites") <- snv$sites
  genotypes
}

#' @rdname cellcoal_read
#' @export
cellcoal_haplotypes <- function(run, replicate = 1, write_as_string = TRUE) {
  stopifnot(is.logical(write_as_string), length(write_as_string) == 1, !is.na(write_as_string))
  snv <- .cellcoal_phylip(.cellcoal_lines(.cellcoal_results_dir(run), "snv_hap", replicate))
  # rows <cell>m (maternal) and <cell>p (paternal), one string each
  cells <- sub("[mp]$", "", snv$names)
  haplotypes <- vapply(snv$values, paste, "", collapse = "")
  maternal <- grepl("m$", snv$names)

  out <- if (write_as_string) {
    cell <- unique(cells)
    data.frame(cell = cell, allele_1 = haplotypes[maternal][match(cell, cells[maternal])],
               allele_2 = haplotypes[!maternal][match(cell, cells[!maternal])])
  } else {
    sites <- do.call(rbind, strsplit(haplotypes, ""))
    data.frame(cell = cells, allele = ifelse(maternal, 1L, 2L),
               .cellcoal_site_columns(if (is.null(sites)) character() else sites, snv$sites))
  }
  attr(out, "sites") <- snv$sites
  out
}

#' @rdname cellcoal_read
#' @export
cellcoal_ado <- function(run, replicate = 1) {
  lines <- .cellcoal_lines(.cellcoal_results_dir(run), "ado_cells", replicate)
  lines <- lines[nzchar(trimws(lines))]
  utils::read.table(text = lines, header = TRUE, sep = "\t", col.names = c("cell", "ado_rate"),
                    colClasses = c("character", "numeric"))
}

#' @rdname cellcoal_read
#' @export
cellcoal_vcf <- function(run, replicate = 1, fields = c("GT", "DP", "TG")) {
  stopifnot(is.character(fields), length(fields) >= 1)
  lines <- .cellcoal_lines(.cellcoal_results_dir(run), "vcf", replicate)
  meta <- grep("^##", lines, value = TRUE)
  header <- grep("^#CHROM", lines, value = TRUE)
  if (length(header) == 0) stop("The VCF of replicate ", replicate, " has no header line.", call. = FALSE)
  cells <- strsplit(header[1], "\t", fixed = TRUE)[[1]][-(1:9)]
  rows <- strsplit(grep("^#", lines, value = TRUE, invert = TRUE), "\t", fixed = TRUE)
  rows <- rows[lengths(rows) > 0]

  # Number and Type of every FORMAT field, from the header
  format_lines <- grep("^##FORMAT=<ID=", meta, value = TRUE)
  ids <- sub("^##FORMAT=<ID=([^,]+),.*", "\\1", format_lines)
  types <- stats::setNames(sub(".*Type=([^,>]+).*", "\\1", format_lines), ids)
  numbers <- stats::setNames(sub(".*Number=([^,>]+).*", "\\1", format_lines), ids)
  unknown <- setdiff(fields, ids)
  if (length(unknown) > 0) {
    stop("The VCF has no field ", paste(unknown, collapse = ", "), "; it has ", paste(ids, collapse = ", "), ".",
         call. = FALSE)
  }

  column <- function(i) vapply(rows, `[`, "", i)
  positions <- as.integer(column(2))
  sites <- data.frame(chrom = column(1), position = positions, id = column(3), ref = column(4), alt = column(5),
                      info = column(8))
  if (length(rows) == 0) sites <- sites[0, ]
  site_names <- make.unique(paste0("site_", positions))

  # sites x cells entries "GT:DP:...", split into their fields by the FORMAT of their site
  keys <- strsplit(column(9), ":", fixed = TRUE)
  entries <- strsplit(unlist(lapply(rows, `[`, -(1:9))), ":", fixed = TRUE)
  site_of_entry <- rep(seq_along(rows), each = length(cells))

  values <- lapply(fields, function(field) {
    index <- vapply(keys, match, 0L, x = field)
    value <- vapply(seq_along(entries), function(k) entries[[k]][index[site_of_entry[k]]], "")
    value[value %in% c(".", "") | is.na(value)] <- NA
    if (numbers[[field]] == "1" && types[[field]] == "Integer") value <- as.integer(value)
    if (numbers[[field]] == "1" && types[[field]] == "Float") value <- as.numeric(value)
    # cells x sites
    value <- matrix(value, nrow = length(cells), ncol = length(rows), dimnames = list(NULL, site_names))
    data.frame(cell = cells, value, check.names = FALSE)
  })
  c(list(sites = sites), stats::setNames(values, fields))
}
