# Interface to CellCoal (Posada 2020, https://github.com/dapogon/cellcoal), an external program (C, GPL-2)
# that simulates the genealogy, somatic mutations and single-cell sequencing of a sample of cells.
# scTreeSim does not bundle it: cellcoal_install() downloads a release binary, cellcoal_path() finds it and
# cellcoal_run() runs it in a folder of its own and reports what it wrote.
#
# What the code relies on (checked with CellCoal 2.7.0):
# - single-letter options take their value attached (-n10, -Ttree.nwk), and file and folder names can be
#   at most 119 characters, so CellCoal runs inside the run folder with short relative names;
# - the last -o wins over an earlier -o, `--output dir=` and `results_folder` in a config file, and a
#   results folder that exists is deleted;
# - options given after --config override the file;
# - errors are lines starting with "!!!" and give a non-zero exit status.

.cellcoal_version <- "2.7.0"
.cellcoal_repo <- "https://github.com/dapogon/cellcoal"
.cellcoal_results <- "results"

# folder of the binaries installed by cellcoal_install()
.cellcoal_home <- function() file.path(tools::R_user_dir("scTreeSim", "data"), "cellcoal")

# name of the release asset of this computer, e.g. "macos-arm64"
.cellcoal_platform <- function() {
  info <- Sys.info()
  arch <- switch(tolower(info[["machine"]]), arm64 = , aarch64 = "arm64",
                 x86_64 = , `x86-64` = , amd64 = "x86_64", tolower(info[["machine"]]))
  os <- switch(info[["sysname"]], Darwin = "macos", Linux = "linux", Windows = "windows", tolower(info[["sysname"]]))
  platform <- paste0(os, "-", arch)
  if (!platform %in% c("macos-arm64", "macos-x86_64", "linux-x86_64", "windows-x86_64")) {
    stop("There is no CellCoal release for ", platform, ". Build it from source (", .cellcoal_repo,
         ") and set options(scTreeSim.cellcoal = \"/path/to/cellcoal\").", call. = FALSE)
  }
  platform
}

.cellcoal_exe <- function(folder) {
  file.path(folder, if (.Platform$OS.type == "windows") "cellcoal.exe" else "cellcoal")
}

# binaries installed by cellcoal_install(), newest version first
.cellcoal_installed <- function() {
  versions <- list.files(.cellcoal_home())
  versions <- versions[!is.na(numeric_version(versions, strict = FALSE))]
  versions <- versions[order(numeric_version(versions), decreasing = TRUE)]
  exes <- vapply(versions, function(v) {
    .cellcoal_exe(file.path(.cellcoal_home(), v, paste0("cellcoal-", v, "-", .cellcoal_platform())))
  }, character(1))
  unname(exes[file.exists(exes)])
}


#' CellCoal documentation
#'
#' Links to the documentation of CellCoal: the repository, the manual (all the options, models and output
#' files), the list of changes, the releases and the paper. If CellCoal was installed with
#' [cellcoal_install()], the manual of the installed version is also listed (`local_manual`).
#'
#' @param open open the manual in the browser
#'
#' @return the links, a named character vector (invisibly)
#' @family CellCoal
#' @export
cellcoal_docs <- function(open = FALSE) {
  links <- c(
    repository = .cellcoal_repo,
    manual = "https://dapogon.github.io/cellcoal/cellcoal.manual.v2.0.html",
    manual_pdf = paste0(.cellcoal_repo, "/blob/master/docs/cellcoal.manual.v2.0.pdf"),
    changes = paste0(.cellcoal_repo, "/blob/master/docs/CHANGELOG.md"),
    releases = paste0(.cellcoal_repo, "/releases"),
    paper = "https://academic.oup.com/mbe/article/37/5/1535/5728646"
  )
  installed <- .cellcoal_installed()
  local_manual <- file.path(dirname(installed[1]), "cellcoal.manual.v2.0.html")
  if (length(installed) > 0 && file.exists(local_manual)) links[["local_manual"]] <- local_manual

  cat(sprintf("%-13s %s", paste0(names(links), ":"), links), sep = "\n")
  cat("Citation: Posada D. (2020) CellCoal: coalescent simulation of single-cell sequencing samples.",
      "Molecular Biology and Evolution 37: 1535-1542.\n")
  if (open) utils::browseURL(links[["manual"]])
  invisible(links)
}


#' Install CellCoal
#'
#' Downloads a release of CellCoal for this computer (Linux x86_64, macOS arm64 or x86_64, Windows x86_64)
#' from GitHub, checks it against its SHA-256 checksum and unpacks it in the user data folder of scTreeSim
#' (`tools::R_user_dir("scTreeSim", "data")`), where [cellcoal_path()] finds it. The release also contains
#' the manual and examples of parameter and config files.
#'
#' CellCoal can also be installed with Homebrew (`brew tap dapogon/cellcoal https://github.com/dapogon/cellcoal`
#' then `brew install cellcoal`) or built from source; see [cellcoal_docs()].
#'
#' @param version CellCoal version
#' @param force install again if this version is already installed
#'
#' @return the path to the binary (invisibly)
#' @family CellCoal
#' @export
cellcoal_install <- function(version = .cellcoal_version, force = FALSE) {
  name <- paste0("cellcoal-", version, "-", .cellcoal_platform())
  windows <- grepl("windows", name)
  archive <- paste0(name, if (windows) ".zip" else ".tar.gz")
  dest <- file.path(.cellcoal_home(), version)
  exe <- .cellcoal_exe(file.path(dest, name))
  if (file.exists(exe) && !force) {
    message("CellCoal ", version, " is already installed: ", exe)
    return(invisible(exe))
  }

  url <- paste0(.cellcoal_repo, "/releases/download/v", version, "/", archive)
  tmp <- tempfile("cellcoal-download-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  file <- file.path(tmp, archive)
  for (suffix in c("", ".sha256")) {
    ok <- tryCatch(utils::download.file(paste0(url, suffix), paste0(file, suffix), mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) stop("Could not download ", url, suffix, ". Is the version right? See ", .cellcoal_repo, "/releases", call. = FALSE)
  }
  .check_sha256(file, paste0(file, ".sha256"))

  unlink(file.path(dest, name), recursive = TRUE)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  if (windows) utils::unzip(file, exdir = dest) else utils::untar(file, exdir = dest)
  if (!file.exists(exe)) stop("The release ", archive, " does not contain ", basename(exe), ".", call. = FALSE)
  Sys.chmod(exe, "755")

  message("CellCoal ", cellcoal_version(exe), " installed: ", exe)
  invisible(exe)
}

.check_sha256 <- function(file, sums) {
  expected <- strsplit(trimws(readLines(sums, warn = FALSE)[1]), "\\s+")[[1]][1]
  # tools::sha256sum() exists from R 4.5.0
  if (!exists("sha256sum", envir = asNamespace("tools"))) {
    warning("The checksum of the download was not checked (it needs R >= 4.5.0).", call. = FALSE)
    return(invisible())
  }
  actual <- unname(get("sha256sum", envir = asNamespace("tools"))(file))
  if (!identical(tolower(expected), tolower(actual))) {
    stop("The checksum of ", basename(file), " is wrong; the download is damaged. Try again.", call. = FALSE)
  }
}


#' Find CellCoal
#'
#' Looks for the CellCoal binary, in this order:
#' 1. the option `scTreeSim.cellcoal`, e.g. `options(scTreeSim.cellcoal = "/path/to/cellcoal")`,
#' 2. the environment variable `CELLCOAL`,
#' 3. the binaries installed by [cellcoal_install()], newest version first,
#' 4. `cellcoal` on the `PATH` (e.g. installed with Homebrew).
#'
#' `cellcoal_available()` says whether it is found, and `cellcoal_version()` gives the version of a binary.
#'
#' @param path path to a CellCoal binary
#'
#' @return `cellcoal_path()`: the path to the binary; `cellcoal_available()`: `TRUE` or `FALSE`;
#'   `cellcoal_version()`: the version, a string such as `"2.7.0"`
#' @family CellCoal
#' @export
cellcoal_path <- function() {
  # a path the user set must exist: do not fall back to another binary silently
  for (set in c(option = getOption("scTreeSim.cellcoal", ""), environment = Sys.getenv("CELLCOAL"))) {
    if (nzchar(set)) {
      if (!file.exists(set)) stop("CellCoal is not at ", set, " (set by the option scTreeSim.cellcoal or CELLCOAL).", call. = FALSE)
      return(normalizePath(set))
    }
  }
  found <- c(.cellcoal_installed(), unname(Sys.which("cellcoal")))
  found <- found[nzchar(found)]
  if (length(found) == 0) {
    stop("CellCoal was not found. Install it with cellcoal_install(), or set ",
         "options(scTreeSim.cellcoal = \"/path/to/cellcoal\").", call. = FALSE)
  }
  normalizePath(found[1])
}

#' @rdname cellcoal_path
#' @export
cellcoal_available <- function() {
  tryCatch(nzchar(cellcoal_path()), error = function(e) FALSE)
}

#' @rdname cellcoal_path
#' @export
cellcoal_version <- function(path = cellcoal_path()) {
  out <- processx::run(path, "--version", error_on_status = FALSE)
  version <- regmatches(out$stdout, regexpr("[0-9]+\\.[0-9]+(\\.[0-9]+)?", out$stdout))
  if (length(version) == 0) stop(path, " does not look like CellCoal: `--version` printed no version.", call. = FALSE)
  version
}


#' Run CellCoal
#'
#' Runs CellCoal with the options of the command line (`args`), of a config file (`config`) or both, in a
#' folder of its own, and returns what it did and the files it wrote. The options are those of CellCoal
#' itself, see [cellcoal_docs()]; `cellcoal_check()` validates them (CellCoal `--check`) and shows what the
#' run would do, without simulating.
#'
#' The run folder (`dir`) holds the config file and the input files, and CellCoal writes its results to its
#' subfolder `results`: scTreeSim sets the results folder, so `-o`, `--output dir=` and `results_folder`
#' are ignored, with a warning for `-o` and `--output dir=`.
#'
#' Files named in `args` (the tree of `-T`, the genome of `-U`, a parameter file of `-F`, a file after an
#' option, or a `key=file` of a group of options) are copied into the run folder, and `args` refers to their
#' copies. Files named inside a config file are not: pass them in `files` and name them in the config by
#' their file name only (e.g. `tree_file: tree.nwk`).
#'
#' @param args options of the command line: a character vector (`c("-n1", "-s6")`) or one string as typed
#'   in a terminal (`"-n1 -s6 -l1000 -1 -6 -#1"`); a leading program name (`cellcoal`, `./bin/cellcoal-2.7.0`)
#'   is dropped. Single-letter options take their value attached, as in CellCoal (`-n10`, not `-n 10`).
#' @param config a config file (YAML or JSON), or a named list of options written to one (keys as in
#'   `cellcoal --help`, e.g. `list(cells = 20, sites = 1000, print_trees = TRUE)`; needs the yaml package).
#'   Options in `args` override the config.
#' @param files more input files to copy into the run folder
#' @param dir run folder: it must not exist or must be empty. By default, a new folder in the temporary
#'   directory of the session.
#' @param timeout maximum time of the run in seconds
#' @param echo print the output of CellCoal while it runs
#' @param error stop with the error of CellCoal if the run fails; with `FALSE`, return the result with
#'   `ok = FALSE`
#'
#' @return an object of class `cellcoal_run`, a list with
#'   * `ok`: whether CellCoal succeeded,
#'   * `status`: the exit status,
#'   * `args`: the arguments CellCoal was run with,
#'   * `dir`, `results`: the run folder and the results folder (`NA` if it was not written),
#'   * `files`: the files in the results folder,
#'   * `errors`, `notes`: the errors (lines `!!!`) and notes (lines `NOTE:`) of CellCoal,
#'   * `run`: the record of the run (`results/run.json`, needs the jsonlite package), with the options,
#'     the seed and the version,
#'   * `stdout`, `stderr`: the output of CellCoal,
#'   * `version`: the version of CellCoal.
#' @family CellCoal
#' @export
#' @examples
#' \dontrun{
#' run <- cellcoal_run("-n1 -s6 -l1000 -u1e-6 -e10000 -1 -6 -#1")
#' run$files
#' readLines(file.path(run$results, "trees"))
#'
#' cellcoal_check(config = list(cells = 20, sites = 1000, coverage = 10))
#' }
cellcoal_run <- function(args = character(), config = NULL, files = NULL, dir = NULL, timeout = Inf,
                         echo = FALSE, error = TRUE) {
  .cellcoal_exec(args, config, files, dir, check = FALSE, timeout = timeout, echo = echo, error = error)
}

#' @rdname cellcoal_run
#' @export
cellcoal_check <- function(args = character(), config = NULL, files = NULL, error = FALSE) {
  .cellcoal_exec(args, config, files, dir = NULL, check = TRUE, timeout = 60, echo = FALSE, error = error)
}

.cellcoal_exec <- function(args, config, files, dir, check, timeout, echo, error) {
  path <- cellcoal_path()
  version <- .cellcoal_check_version(path)
  if (!is.character(args)) stop("`args` must be a character vector.", call. = FALSE)
  if (length(args) == 1) args <- .cellcoal_split(args)
  if (length(args) > 0 && grepl("^cellcoal", basename(args[1]), ignore.case = TRUE)) args <- args[-1]
  .cellcoal_warn_output(args, config)

  if (is.null(dir)) dir <- tempfile("cellcoal-")
  if (dir.exists(dir) && length(list.files(dir, all.files = TRUE, no.. = TRUE)) > 0) {
    stop("The run folder ", dir, " is not empty.", call. = FALSE)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir)

  staged <- .cellcoal_stage(files, dir)
  args <- .cellcoal_stage_args(args, dir, staged)
  if (!is.null(config)) config <- .cellcoal_config(config, dir, staged)

  # the config comes first so that `args` override it, and the results folder last so that it wins
  full_args <- c(if (!is.null(config)) c("--config", config), args, if (check) "--check",
                 paste0("-o", .cellcoal_results))
  out <- processx::run(path, full_args, wd = dir, error_on_status = FALSE, timeout = timeout, echo = echo)

  lines <- trimws(unlist(strsplit(c(out$stdout, out$stderr), "\r?\n")))
  results <- file.path(dir, .cellcoal_results)
  if (!dir.exists(results)) results <- NA_character_
  run_json <- file.path(results, "run.json")
  run <- if (!is.na(results) && file.exists(run_json) && requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::read_json(run_json, simplifyVector = TRUE)
  }
  timed_out <- isTRUE(out$timeout)
  result <- structure(list(
    ok = identical(out$status, 0L) && !timed_out,
    status = out$status,
    args = full_args,
    dir = dir,
    results = results,
    files = if (is.na(results)) character() else list.files(results, recursive = TRUE),
    errors = c(.cellcoal_errors(lines),
               if (timed_out) paste0("The run took longer than ", timeout, " seconds and was stopped.")),
    notes = grep("^NOTE:", lines, value = TRUE),
    run = run,
    check = check,
    stdout = out$stdout,
    stderr = out$stderr,
    version = version
  ), class = "cellcoal_run")

  if (error && !result$ok) {
    stop("CellCoal failed (exit status ", out$status, "):\n",
         paste(if (length(result$errors) > 0) result$errors else "no error message, see the output", collapse = "\n"),
         "\nRun folder: ", dir, call. = FALSE)
  }
  result
}

# errors of CellCoal: a line "!!! ... ERROR: ..." (options) or "ERROR: ..." (simulation), with the lines
# that follow it up to a blank line
.cellcoal_errors <- function(lines) {
  starts <- grep("^(!!!|ERROR)", lines)
  vapply(starts, function(start) {
    end <- start
    while (end < length(lines) && nzchar(lines[end + 1])) end <- end + 1
    sub("^!!!\\s*", "", paste(lines[start:end], collapse = " "))
  }, character(1))
}

# version of the binary, warning once per binary if it is older than the supported one
.cellcoal_versions <- new.env(parent = emptyenv())
.cellcoal_check_version <- function(path) {
  if (is.null(.cellcoal_versions[[path]])) {
    version <- cellcoal_version(path)
    if (numeric_version(version) < numeric_version(.cellcoal_version)) {
      warning("CellCoal ", version, " is older than ", .cellcoal_version, ", the version scTreeSim supports. ",
              "Update it with cellcoal_install().", call. = FALSE)
    }
    .cellcoal_versions[[path]] <- version
  }
  .cellcoal_versions[[path]]
}

# split a command line into arguments like a shell: whitespace separates, quotes group, backslash escapes
# (but nothing is expanded and no shell runs it)
.cellcoal_split <- function(command) {
  chars <- strsplit(command, "")[[1]]
  args <- character()
  current <- ""
  in_arg <- FALSE
  quote <- ""
  i <- 1
  while (i <= length(chars)) {
    ch <- chars[i]
    if (nzchar(quote)) {
      if (ch == quote) {
        quote <- ""
      } else if (ch == "\\" && quote == "\"" && i < length(chars) && chars[i + 1] %in% c("\"", "\\")) {
        i <- i + 1
        current <- paste0(current, chars[i])
      } else {
        current <- paste0(current, ch)
      }
    } else if (grepl("\\s", ch)) {
      if (in_arg) args <- c(args, current)
      current <- ""
      in_arg <- FALSE
    } else if (ch %in% c("\"", "'")) {
      quote <- ch
      in_arg <- TRUE
    } else if (ch == "\\" && i < length(chars)) {
      i <- i + 1
      current <- paste0(current, chars[i])
      in_arg <- TRUE
    } else {
      current <- paste0(current, ch)
      in_arg <- TRUE
    }
    i <- i + 1
  }
  if (nzchar(quote)) stop("Unclosed quote (", quote, ") in the CellCoal command.", call. = FALSE)
  if (in_arg) args <- c(args, current)
  args
}

.cellcoal_warn_output <- function(args, config) {
  # -oNAME, or dir=NAME in the group that follows --output
  after_output <- c(FALSE, utils::head(args, -1) == "--output")
  set <- grepl("^-o", args) | (after_output & grepl("(^|,)dir=", args))
  if (is.list(config) && any(names(config) %in% c("results_folder", "o"))) set <- c(set, TRUE)
  if (any(set)) {
    warning("scTreeSim sets the results folder of CellCoal (the folder `results` of the run folder); ",
            "the results folder given is ignored.", call. = FALSE)
  }
}

# copy input files into the run folder; returns their original paths, named by file name
.cellcoal_stage <- function(files, dir, staged = character()) {
  for (file in files) {
    if (!file.exists(file) || dir.exists(file)) stop("The input file ", file, " does not exist.", call. = FALSE)
    file <- normalizePath(file)
    name <- basename(file)
    if (name %in% names(staged)) {
      if (staged[[name]] != file) stop("Two input files are called ", name, "; rename one.", call. = FALSE)
      next
    }
    if (nchar(name) > 119) stop("The name of the input file ", name, " is longer than 119 characters.", call. = FALSE)
    file.copy(file, file.path(dir, name))
    staged[[name]] <- file
  }
  staged
}

# copy the files named in the arguments into the run folder and refer to them by file name
.cellcoal_stage_args <- function(args, dir, staged) {
  is_file <- function(x) nzchar(x) && file.exists(x) && !dir.exists(x)
  stage <- function(x) {
    staged <<- .cellcoal_stage(x, dir, staged)
    basename(x)
  }
  for (i in seq_along(args)) {
    arg <- args[i]
    if (grepl("^-[TUF].", arg) && is_file(substring(arg, 3))) {
      # single-letter options with a file: -T tree, -U genome, -F parameter file
      args[i] <- paste0(substring(arg, 1, 2), stage(substring(arg, 3)))
    } else if (!startsWith(arg, "-") && !grepl("=", arg) && is_file(arg)) {
      # the value of an option such as --config FILE or --fasta FILE
      args[i] <- stage(arg)
    } else if (grepl("=", arg)) {
      # key=file in a group of options (--tree file=tree.nwk,no-outgroup) or --option=file
      parts <- strsplit(arg, ",", fixed = TRUE)[[1]]
      for (j in seq_along(parts)) {
        key <- sub("=.*", "", parts[j])
        value <- sub("^[^=]*=", "", parts[j])
        if (grepl("=", parts[j]) && is_file(value)) parts[j] <- paste0(key, "=", stage(value))
      }
      args[i] <- paste(parts, collapse = ",")
    }
  }
  args
}

# a config file in the run folder, from a file or a list of options; returns its name
.cellcoal_config <- function(config, dir, staged) {
  if (is.character(config) && length(config) == 1) {
    .cellcoal_stage(config, dir, staged)
    return(basename(config))
  }
  if (!is.list(config) || is.null(names(config)) || any(!nzchar(names(config)))) {
    stop("`config` must be a config file or a named list of options.", call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Writing a config needs the yaml package.", call. = FALSE)
  config <- config[!names(config) %in% c("results_folder", "o")]
  # whole numbers as integers (yaml writes 1 as 1.0, which CellCoal passes on as -n1.0)
  config <- lapply(config, function(x) {
    if (is.double(x) && all(is.finite(x)) && all(x == round(x)) && all(abs(x) < .Machine$integer.max)) as.integer(x) else x
  })
  file <- "cellcoal_config.yaml"
  # logicals as true/false (the default of yaml is yes/no)
  yaml::write_yaml(config, file.path(dir, file), handlers = list(logical = yaml::verbatim_logical))
  file
}


#' @export
print.cellcoal_run <- function(x, ...) {
  cat("CellCoal", x$version, if (x$check) "check" else "run", if (x$ok) "(ok)" else paste0("(failed, exit status ", x$status, ")"), "\n")
  cat("Arguments:", paste(x$args, collapse = " "), "\n")
  cat("Run folder:", x$dir, "\n")
  if (length(x$errors) > 0) cat(paste("Error:", x$errors), sep = "\n")
  if (length(x$notes) > 0) cat(x$notes, sep = "\n")
  if (x$check && x$ok) {
    # the summary of --check: from "the options are valid" to the end
    lines <- strsplit(x$stdout, "\r?\n")[[1]]
    start <- grep("the options are valid", lines)
    if (length(start) > 0) cat(lines[start[1]:length(lines)], sep = "\n")
  } else if (!x$check) {
    cat("Results:", if (is.na(x$results)) "none" else x$results, "\n")
    if (length(x$files) > 0) cat("Files:", paste(x$files, collapse = ", "), "\n")
  }
  invisible(x)
}
