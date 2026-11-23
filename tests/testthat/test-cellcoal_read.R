# Readers of the results of CellCoal. The runs are skipped when CellCoal is not installed (helper-cellcoal.R).

genotype_args <- "-n3 -s5 -l500 -u1e-6 -1 -2 -6 -#2"
readcount_args <- "-n1 -s6 -l500 -u1e-6 -C10 -E0.001 -P0.1 0.01 -1 -2 -6 -#3"

# ---- tree, genotypes and haplotypes ----

test_that("the tree has the cells and the outgroup", {
  skip_if_no_cellcoal()
  run <- local_cellcoal_run(genotype_args)
  tree <- cellcoal_tree(run)
  expect_s3_class(tree, "phylo")
  expect_setequal(tree$tip.label, c(sprintf("cell%04d", 1:5), "outgcell"))
  expect_equal(ape::write.tree(tree), readLines(file.path(run$results, "trees"))[1])
})

test_that("genotypes as strings and as columns agree", {
  skip_if_no_cellcoal()
  run <- local_cellcoal_run(genotype_args)
  strings <- cellcoal_genotypes(run)
  columns <- cellcoal_genotypes(run, write_as_string = FALSE)
  sites <- attr(strings, "sites")

  expect_equal(names(strings), c("cell", "genotype"))
  expect_equal(strings$cell, c(sprintf("cell%04d", 1:5), "outgcell"))
  expect_equal(names(columns), c("cell", paste0("site_", sites)))
  expect_equal(attr(columns, "sites"), sites)
  expect_equal(strings$genotype, apply(as.matrix(columns[, -1]), 1, paste, collapse = "_"))
  expect_true(all(grepl("^[ACGT]{2}$", as.matrix(columns[, -1]))))
})

test_that("without errors, a genotype is the maternal and the paternal allele", {
  skip_if_no_cellcoal()
  run <- local_cellcoal_run(genotype_args)
  genotypes <- cellcoal_genotypes(run, write_as_string = FALSE)
  haplotypes <- cellcoal_haplotypes(run)
  expect_equal(names(haplotypes), c("cell", "allele_1", "allele_2"))
  expect_equal(haplotypes$cell, genotypes$cell)
  expect_equal(nchar(haplotypes$allele_1), rep(length(attr(genotypes, "sites")), nrow(genotypes)))

  paired <- t(mapply(function(m, p) paste0(strsplit(m, "")[[1]], strsplit(p, "")[[1]]),
                     haplotypes$allele_1, haplotypes$allele_2, USE.NAMES = FALSE))
  expect_equal(unname(as.matrix(genotypes[, -1])), paired)

  by_site <- cellcoal_haplotypes(run, write_as_string = FALSE)
  expect_equal(nrow(by_site), 2 * nrow(haplotypes))
  expect_equal(by_site$allele, rep(1:2, nrow(haplotypes)))
  expect_equal(paste(by_site[1, -(1:2)], collapse = ""), haplotypes$allele_1[1])
})

# ---- replicates ----

test_that("every replicate is read, from one file or from a file per replicate", {
  skip_if_no_cellcoal()
  run <- local_cellcoal_run(genotype_args)
  per_replicate <- local_cellcoal_run(paste(genotype_args, "-v"))
  expect_true(any(startsWith(per_replicate$files, "snv_genotypes_dir/")))

  for (replicate in 1:3) {
    expect_equal(cellcoal_genotypes(per_replicate, replicate), cellcoal_genotypes(run, replicate))
    expect_equal(cellcoal_tree(per_replicate, replicate), cellcoal_tree(run, replicate))
  }
  expect_false(identical(cellcoal_tree(run, 1), cellcoal_tree(run, 2)))
  expect_error(cellcoal_genotypes(run, 4), "3 replicates")
  expect_error(cellcoal_tree(run, 4), "3 replicates")
})

test_that("a replicate without variable sites has no rows", {
  skip_if_no_cellcoal()
  run <- local_cellcoal_run("-n1 -s3 -l10 -u1e-9 -1 -2 -6 -#2")
  expect_equal(nrow(cellcoal_genotypes(run)), 0)
  expect_equal(nrow(cellcoal_haplotypes(run)), 0)
  expect_equal(ape::Ntip(cellcoal_tree(run)), 4)
})

# ---- read counts ----

test_that("the VCF has the sites of the genotypes and typed fields", {
  skip_if_no_cellcoal()
  run <- local_cellcoal_run(readcount_args)
  vcf <- cellcoal_vcf(run, fields = c("GT", "DP", "RC", "TG"))
  genotypes <- cellcoal_genotypes(run)

  expect_equal(names(vcf), c("sites", "GT", "DP", "RC", "TG"))
  expect_equal(vcf$sites$position, attr(genotypes, "sites"))
  expect_equal(vcf$DP$cell, genotypes$cell)
  expect_type(vcf$DP[[2]], "integer")
  expect_type(vcf$RC[[2]], "character")
  # the depth is the sum of the read counts of A, C, G and T
  counts <- vapply(strsplit(as.matrix(vcf$RC[, -1]), ","), function(x) sum(as.integer(x)), 0)
  expect_equal(counts, as.vector(as.matrix(vcf$DP[, -1])))
  expect_error(cellcoal_vcf(run, fields = "XX"), "no field XX")
})

test_that("the dropout rate of every cell is read", {
  skip_if_no_cellcoal()
  run <- local_cellcoal_run(readcount_args)
  ado <- cellcoal_ado(run)
  expect_equal(names(ado), c("cell", "ado_rate"))
  expect_equal(ado$cell, c(sprintf("cell%04d", 1:6), "outgcell"))
  expect_true(all(ado$ado_rate >= 0 & ado$ado_rate <= 1))
})

# ---- all the results ----

test_that("cellcoal_read reads what the run wrote", {
  skip_if_no_cellcoal()
  expect_named(cellcoal_read(local_cellcoal_run(genotype_args)), c("tree", "genotypes", "haplotypes"))
  run <- local_cellcoal_run(readcount_args)
  results <- cellcoal_read(run)
  expect_named(results, c("tree", "genotypes", "haplotypes", "ado", "vcf"))
  expect_equal(results$vcf, cellcoal_vcf(run))
  expect_equal(results$ado, cellcoal_ado(run))
})

test_that("the readers need a run with results", {
  expect_error(cellcoal_tree(file.path(tempdir(), "no-results")), "results folder")
  expect_error(cellcoal_tree(NA_character_), "results folder")
  skip_if_no_cellcoal()
  run <- local_cellcoal_run("-n1 -s4 -l100 -#1 -1")
  expect_error(cellcoal_tree(run), "no file `trees`")
})
