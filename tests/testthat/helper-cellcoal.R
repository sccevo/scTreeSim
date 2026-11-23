# CellCoal is an external program: the tests that run it are skipped when it is not installed
# (cellcoal_install())
skip_if_no_cellcoal <- function() {
  skip_on_cran()
  skip_if_not(cellcoal_available(), "CellCoal is not installed")
}

# a finished run of CellCoal in a new temporary folder, removed at the end of the calling test
local_cellcoal_run <- function(args, ..., env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  cellcoal_run(args, dir = file.path(dir, "run"), ...)
}
