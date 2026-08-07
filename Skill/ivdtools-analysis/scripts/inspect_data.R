#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[[i + 1L]]
}

input <- value_after("--input")
if (is.null(input)) stop("Usage: inspect_data.R --input FILE [--sheet SHEET] [--encoding UTF-8] [--output FILE|-]", call. = FALSE)
input <- normalizePath(input, mustWork = TRUE)
sheet <- value_after("--sheet")
encoding <- value_after("--encoding", "UTF-8")
output <- value_after("--output", "-")
max_levels <- as.integer(value_after("--max-levels", "20"))
if (is.na(max_levels) || max_levels < 1L) stop("--max-levels must be a positive integer", call. = FALSE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package jsonlite is required", call. = FALSE)

ext <- tolower(tools::file_ext(input))
if (ext %in% c("xls", "xlsx")) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Package readxl is required for Excel files", call. = FALSE)
  sheets <- readxl::excel_sheets(input)
  if (is.null(sheet)) {
    if (length(sheets) > 1L) stop("Workbook has multiple sheets; pass --sheet. Available: ", paste(sheets, collapse = ", "), call. = FALSE)
    sheet <- sheets[[1L]]
  }
  data <- as.data.frame(readxl::read_excel(input, sheet = sheet), check.names = FALSE)
} else if (ext == "rds") {
  data <- readRDS(input)
  if (!is.data.frame(data)) stop("RDS object must be a data frame", call. = FALSE)
} else if (ext == "csv") {
  data <- utils::read.csv(input, check.names = FALSE, fileEncoding = encoding)
} else if (ext %in% c("tsv", "txt")) {
  data <- utils::read.delim(input, check.names = FALSE, fileEncoding = encoding)
} else {
  stop("Unsupported extension: ", ext, call. = FALSE)
}

summarize_column <- function(x, name) {
  nonmissing <- x[!is.na(x)]
  unique_values <- unique(nonmissing)
  samples <- head(unique_values, 5L)
  id_like_name <- grepl(
    "(^id$|(^|[._ -])id($|[._ -])|sample.*id|样本.*(id|编号|号)|^编号$)",
    tolower(name), perl = TRUE
  )
  list(
    name = name,
    class = paste(class(x), collapse = "/"),
    missing = sum(is.na(x)),
    missing_percent = if (length(x)) round(100 * mean(is.na(x)), 2) else 0,
    unique = length(unique_values),
    candidate_binary = length(unique_values) == 2L,
    candidate_id = id_like_name && length(x) > 0L && !anyNA(x) &&
      length(unique_values) == length(x),
    levels = if (length(unique_values) <= max_levels) as.character(unique_values) else NULL,
    sample = as.character(samples),
    numeric_min = if (is.numeric(x) && length(nonmissing)) min(nonmissing) else NULL,
    numeric_max = if (is.numeric(x) && length(nonmissing)) max(nonmissing) else NULL
  )
}

result <- list(
  source = input,
  format = ext,
  sheet = sheet,
  rows = nrow(data),
  columns_count = ncol(data),
  duplicate_column_names = unique(names(data)[duplicated(names(data))]),
  duplicate_rows = sum(duplicated(data)),
  columns = unname(Map(summarize_column, data, names(data)))
)
json <- jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
if (identical(output, "-")) {
  cat(json, "\n")
} else {
  output <- path.expand(output)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  writeLines(json, output, useBytes = TRUE)
  cat("Wrote", output, "\n")
}
