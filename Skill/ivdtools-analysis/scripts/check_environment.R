#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  i <- match(flag, args)
  if (is.na(i)) return(NULL)
  if (i == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[[i + 1L]]
}

extra_lib <- value_after("--library")
if (!is.null(extra_lib)) .libPaths(c(path.expand(extra_lib), .libPaths()))

# ivdtools is published on CRAN; require at least the version that includes
# lob_lod_loq() (analytical sensitivity / EP17-A2 LoB-LoD-LoQ).
min_version <- "0.1.1"
required_imports <- c("ggplot2", "ggrepel", "minpack.lm", "nloptr", "nls2",
                      "nortest", "rlang", "VCA", "VFP")
helper_packages <- c("jsonlite", "readxl", "rmarkdown", "knitr")
problems <- character()

if (getRversion() < "4.1.0") {
  problems <- c(problems, paste0("R >= 4.1.0 is required; found ", getRversion()))
}

if (!requireNamespace("ivdtools", quietly = TRUE)) {
  problems <- c(problems, "ivdtools is not installed")
  ivd_version <- NA_character_
} else {
  ivd_version <- as.character(utils::packageVersion("ivdtools"))
  if (package_version(ivd_version) < min_version) {
    problems <- c(problems, paste0("ivdtools >= ", min_version,
                                   " is required; found ", ivd_version))
  }
}

missing_imports <- required_imports[!vapply(required_imports, requireNamespace,
                                            quietly = TRUE, FUN.VALUE = logical(1))]
missing_helpers <- helper_packages[!vapply(helper_packages, requireNamespace,
                                           quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_imports)) {
  problems <- c(problems, paste("Missing ivdtools dependencies:", paste(missing_imports, collapse = ", ")))
}
if (length(missing_helpers)) {
  problems <- c(problems, paste("Missing report/input helpers:", paste(missing_helpers, collapse = ", ")))
}

pandoc <- Sys.which("pandoc")
if (!nzchar(pandoc) && requireNamespace("rmarkdown", quietly = TRUE) &&
    !rmarkdown::pandoc_available()) {
  problems <- c(problems, "Pandoc is unavailable; HTML rendering will fail")
}

cat("R:", R.version.string, "\n")
cat("Library paths:\n", paste0("- ", .libPaths(), collapse = "\n"), "\n", sep = "")
cat("ivdtools:", ifelse(is.na(ivd_version), "not installed", ivd_version), "\n")
cat("Pandoc:", if (nzchar(pandoc)) pandoc else if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) rmarkdown::pandoc_version() else "not found", "\n")

if (length(problems)) {
  cat("STATUS: NOT READY\n", paste0("- ", problems, collapse = "\n"), "\n", sep = "")
  quit(status = 1L)
}
cat("STATUS: READY\n")
