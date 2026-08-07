#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[[i + 1L]]
}

input <- value_after("--input")
if (is.null(input)) stop("Usage: render_report.R --input report.Rmd [--output-dir DIR] [--output-file report.html]", call. = FALSE)
input <- normalizePath(input, mustWork = TRUE)
output_dir <- value_after("--output-dir", dirname(input))
output_file <- value_after("--output-file", "report.html")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Package rmarkdown is required", call. = FALSE)
if (!rmarkdown::pandoc_available()) stop("Pandoc is required", call. = FALSE)

rendered <- rmarkdown::render(
  input = input,
  output_file = output_file,
  output_dir = output_dir,
  envir = new.env(parent = globalenv()),
  clean = TRUE,
  quiet = FALSE,
  encoding = "UTF-8"
)
cat(normalizePath(rendered, mustWork = TRUE), "\n")
