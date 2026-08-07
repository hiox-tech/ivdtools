#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(flag) flag %in% args
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[[i + 1L]]
}

if (!has_flag("--yes")) {
  stop("Installation changes an R library. Re-run with --yes only after explicit user authorization.", call. = FALSE)
}

repos <- value_after("--repos", "https://cloud.r-project.org")
library_path <- value_after("--library", NULL)
if (is.null(library_path)) {
  library_path <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(library_path)) library_path <- .libPaths()[1L]
}
library_path <- path.expand(library_path)
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(library_path) || file.access(library_path, 2L) != 0L) {
  stop("R library is not writable: ", library_path, call. = FALSE)
}
.libPaths(c(library_path, .libPaths()))

# ivdtools is published on CRAN; require at least the version that includes
# lob_lod_loq() (analytical sensitivity / EP17-A2 LoB-LoD-LoQ).
min_version <- "0.1.1"
dependencies <- c("ggplot2", "ggrepel", "minpack.lm", "nloptr", "nls2",
                  "nortest", "rlang", "VCA", "VFP", "jsonlite", "readxl",
                  "rmarkdown", "knitr")

missing <- dependencies[!vapply(dependencies, requireNamespace, quietly = TRUE,
                                FUN.VALUE = logical(1))]
if (length(missing)) {
  message("Installing dependencies from CRAN: ", paste(missing, collapse = ", "))
  utils::install.packages(missing, lib = library_path, repos = repos)
}

if (!requireNamespace("ivdtools", quietly = TRUE)) {
  message("Installing ivdtools from CRAN (>= ", min_version, ")")
  utils::install.packages("ivdtools", lib = library_path, repos = repos,
                          dependencies = FALSE)
}

if (!requireNamespace("ivdtools", quietly = TRUE))
  stop("ivdtools installation failed", call. = FALSE)
installed <- as.character(utils::packageVersion("ivdtools"))
if (package_version(installed) < min_version) {
  stop("Installed ivdtools version is ", installed, "; need >= ", min_version,
       ". Update from CRAN (or install into a library you control) before using this skill.",
       call. = FALSE)
}
cat("Installed ivdtools", installed, "in", library_path, "\n")
