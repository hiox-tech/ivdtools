# Internal helpers shared by analysis modules.

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but is not installed.", call. = FALSE)
  }
  invisible(TRUE)
}

.format_num <- function(x, digits = 4L) {
  formatC(x, format = "f", digits = digits)
}

.print_df <- function(df) {
  lines <- utils::capture.output(print.data.frame(df, row.names = FALSE))
  if (!length(lines)) {
    return(invisible(df))
  }
  cat(lines[1L], "\n", sep = "")
  cat(strrep("-", nchar(lines[1L])), "\n", sep = "")
  if (length(lines) > 1L) {
    cat(paste(lines[-1L], collapse = "\n"), "\n")
  }
  invisible(df)
}

.print_df_sections <- function(df, sections, indent = "  ") {
  if (!is.data.frame(df)) stop("`df` must be a data frame.", call. = FALSE)
  if (!is.list(sections) || !length(sections) || is.null(names(sections)) ||
      any(!nzchar(names(sections)))) {
    stop("`sections` must be a named list of column groups.", call. = FALSE)
  }
  used <- character()
  for (title in names(sections)) {
    columns <- unique(as.character(sections[[title]]))
    if (!length(columns)) next
    missing <- setdiff(columns, names(df))
    if (length(missing)) {
      stop("Unknown column(s) in section '", title, "': ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    cat(indent, title, "\n", sep = "")
    .print_df(df[, columns, drop = FALSE])
    used <- union(used, columns)
  }
  remaining <- setdiff(names(df), used)
  if (length(remaining)) {
    cat(indent, "Additional fields\n", sep = "")
    .print_df(df[, remaining, drop = FALSE])
  }
  invisible(df)
}

.print_value_text <- function(value) {
  if (is.factor(value)) value <- as.character(value)
  if (is.list(value)) value <- unlist(value, use.names = FALSE)
  if (!length(value)) return("")
  if (is.numeric(value)) {
    return(paste(format(value, digits = 8L, trim = TRUE), collapse = ", "))
  }
  paste(as.character(value), collapse = ", ")
}

.print_kv_sections <- function(df, sections, indent = "  ",
                               show_sections = FALSE) {
  if (!is.data.frame(df) || nrow(df) != 1L) {
    stop("`df` must be a one-row data frame.", call. = FALSE)
  }
  if (!is.list(sections) || !length(sections) || is.null(names(sections)) ||
      any(!nzchar(names(sections)))) {
    stop("`sections` must be a named list of fields.", call. = FALSE)
  }
  used <- character()
  section_fields <- list()
  for (title in names(sections)) {
    fields <- sections[[title]]
    if (is.null(fields) || !length(fields)) next
    columns <- as.character(unname(fields))
    labels <- names(fields)
    if (is.null(labels)) labels <- columns
    labels[is.na(labels) | !nzchar(labels)] <- columns[is.na(labels) | !nzchar(labels)]
    missing <- setdiff(columns, names(df))
    if (length(missing)) {
      stop("Unknown field(s) in section '", title, "': ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    section_fields[[title]] <- list(columns = columns, labels = labels)
    used <- union(used, columns)
  }
  remaining <- setdiff(names(df), used)
  if (length(remaining)) section_fields[["Additional fields"]] <- list(
    columns = remaining, labels = remaining)
  all_labels <- unlist(lapply(section_fields, `[[`, "labels"),
                       use.names = FALSE)
  common_width <- if (length(all_labels)) max(nchar(all_labels), 1L) else 1L
  for (title in names(section_fields)) {
    fields <- section_fields[[title]]
    if (isTRUE(show_sections)) cat(indent, title, "\n", sep = "")
    width <- if (isTRUE(show_sections)) max(nchar(fields$labels), 1L) else common_width
    field_indent <- if (isTRUE(show_sections)) "  " else ""
    for (i in seq_along(fields$columns)) {
      cat(indent, field_indent, sprintf("%-*s", width, fields$labels[i]), " : ",
          .print_value_text(df[[fields$columns[i]]][[1L]]), "\n", sep = "")
    }
  }
  invisible(df)
}

.print_kv_rows <- function(df, sections, row_title = NULL, indent = "  ",
                           show_sections = FALSE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  if (!is.list(sections) || !length(sections) || is.null(names(sections)) ||
      any(!nzchar(names(sections)))) {
    stop("`sections` must be a named list of fields.", call. = FALSE)
  }
  if (!is.null(row_title)) {
    if (!is.character(row_title) || !length(row_title) ||
        any(is.na(row_title)) || any(!nzchar(unname(row_title)))) {
      stop("`row_title` must contain existing field names.", call. = FALSE)
    }
    title_columns <- unname(row_title)
    title_labels <- names(row_title)
    if (is.null(title_labels)) title_labels <- title_columns
    title_labels[is.na(title_labels) | !nzchar(title_labels)] <-
      title_columns[is.na(title_labels) | !nzchar(title_labels)]
    missing <- setdiff(title_columns, names(df))
    if (length(missing)) {
      stop("Unknown row-title field(s): ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
  } else {
    title_columns <- character()
    title_labels <- character()
  }

  row_heading <- function(i) {
    if (!length(title_columns)) return(paste0("Row ", i))
    values <- vapply(title_columns, function(column) {
      .print_value_text(df[[column]][[i]])
    }, character(1))
    if (is.null(names(row_title))) return(paste(values, collapse = " | "))
    paste(paste0(title_labels, ": ", values), collapse = " | ")
  }

  used <- title_columns
  section_fields <- list()
  for (title in names(sections)) {
    fields <- sections[[title]]
    if (is.null(fields) || !length(fields)) next
    columns <- as.character(unname(fields))
    labels <- names(fields)
    if (is.null(labels)) labels <- columns
    labels[is.na(labels) | !nzchar(labels)] <-
      columns[is.na(labels) | !nzchar(labels)]
    missing <- setdiff(columns, names(df))
    if (length(missing)) {
      stop("Unknown field(s) in section '", title, "': ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    section_fields[[title]] <- list(columns = columns, labels = labels)
    used <- union(used, columns)
  }
  remaining <- setdiff(names(df), used)
  if (length(remaining)) section_fields[["Additional fields"]] <- list(
    columns = remaining, labels = remaining)
  all_labels <- unlist(lapply(section_fields, `[[`, "labels"),
                       use.names = FALSE)
  common_width <- if (length(all_labels)) max(nchar(all_labels), 1L) else 1L

  for (i in seq_len(nrow(df))) {
    cat(indent, row_heading(i), "\n", sep = "")
    for (title in names(section_fields)) {
      fields <- section_fields[[title]]
      if (isTRUE(show_sections)) cat(indent, "  ", title, "\n", sep = "")
      width <- if (isTRUE(show_sections)) max(nchar(fields$labels), 1L) else common_width
      for (j in seq_along(fields$columns)) {
        field_indent <- if (isTRUE(show_sections)) "    " else "  "
        cat(indent, field_indent, sprintf("%-*s", width, fields$labels[j]), " : ",
            .print_value_text(df[[fields$columns[j]]][[i]]), "\n", sep = "")
      }
    }
    if (i < nrow(df)) cat("\n")
  }
  invisible(df)
}
