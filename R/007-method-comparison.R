#' Method Comparison Regression (MCR) Analysis
#'
#' Description: Build mcr class objects, descriptive statistics, correlation/regression/Bland-Altman analysis,
#' outlier detection, bias analysis, and visualization.
#'
#' @importFrom stats cor.test lm coef vcov residuals qt qnorm sd
#'   median quantile
#' @importFrom utils capture.output
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_line
#'   geom_hline geom_ribbon labs theme_bw coord_fixed
#'   coord_equal scale_linetype_manual element_text margin
#'   annotate
#' @noRd
NULL

# Zero-safe divisor (preserves sign, avoids division by zero)
.safe_div <- function(x) {
  ifelse(abs(x) < .Machine$double.eps, .Machine$double.eps, x)
}

#' Split a data frame at fixed numeric cutpoints
#'
#' Split an input data frame into a named list of data frames before analysis.
#' Cutpoints are applied to one numeric variable using left-closed,
#' right-open intervals. For example, cutpoints `c(100, 500)` create
#' `[-Inf, 100)`, `[100, 500)`, and `[500, Inf)`.
#'
#' @param data A data frame.
#' @param variable A single character string naming the numeric variable used
#'   to split the rows.
#' @param cutpoints A non-empty, finite, strictly increasing numeric vector.
#' @param labels Optional unique names for the returned data frames. Its length
#'   must be `length(cutpoints) + 1`. By default the names are `segment_1`,
#'   `segment_2`, and so on.
#' @param na_action How non-finite splitting values are handled: `"error"`
#'   (default) stops; `"drop"` excludes them and records the excluded count.
#' @param drop_empty Whether empty intervals are removed. The default is
#'   `FALSE`, so the returned list preserves all prespecified intervals.
#'
#' @return A named list of data frames with class `cutpoint_split`. The
#'   `segment_info` attribute records interval bounds, labels, and row counts.
#'
#' @examples
#' parts <- split_by_cutpoints(
#'   ivd_mcr_example, variable = "ref", cutpoints = c(18, 26)
#' )
#' parts
#' parts$segment_1
#' fits <- lapply(parts, function(d) mcr(d, "sid", "test", "ref"))
#' @export
split_by_cutpoints <- function(data, variable, cutpoints, labels = NULL,
                               na_action = c("error", "drop"),
                               drop_empty = FALSE) {
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(variable) || length(variable) != 1L ||
      is.na(variable) || !nzchar(variable) || !variable %in% names(data))
    stop("'variable' must name one column in 'data'.", call. = FALSE)
  values <- data[[variable]]
  if (!is.numeric(values))
    stop("'variable' must identify a numeric column.", call. = FALSE)
  if (!is.numeric(cutpoints) || !length(cutpoints) ||
      any(!is.finite(cutpoints)))
    stop("'cutpoints' must be a non-empty finite numeric vector.",
         call. = FALSE)
  if (is.unsorted(cutpoints, strictly = TRUE))
    stop("'cutpoints' must be strictly increasing and unique.",
         call. = FALSE)
  if (!is.logical(drop_empty) || length(drop_empty) != 1L ||
      is.na(drop_empty))
    stop("'drop_empty' must be TRUE or FALSE.", call. = FALSE)

  na_action <- match.arg(na_action)
  n_segments <- length(cutpoints) + 1L
  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != n_segments ||
        anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
      stop("'labels' must contain one unique, non-empty label per interval.",
           call. = FALSE)
    list_names <- labels
  } else {
    list_names <- paste0("segment_", seq_len(n_segments))
  }

  finite <- is.finite(values)
  n_excluded <- sum(!finite)
  if (n_excluded && na_action == "error")
    stop("The splitting variable contains ", n_excluded,
         " non-finite value(s); use na_action = \"drop\" to exclude them.",
         call. = FALSE)
  if (n_excluded) {
    data <- data[finite, , drop = FALSE]
    values <- values[finite]
  }

  membership <- findInterval(values, cutpoints) + 1L
  lower <- c(-Inf, cutpoints)
  upper <- c(cutpoints, Inf)
  format_bound <- function(value) {
    if (is.infinite(value)) return(if (value < 0) "-Inf" else "Inf")
    format(value, trim = TRUE, scientific = FALSE)
  }
  interval <- vapply(seq_len(n_segments), function(i) {
    paste0("[", format_bound(lower[i]), ", ",
           format_bound(upper[i]), ")")
  }, character(1))
  counts <- tabulate(membership, nbins = n_segments)
  info <- data.frame(
    segment = seq_len(n_segments),
    name = list_names,
    interval = interval,
    lower = lower,
    upper = upper,
    n = counts,
    stringsAsFactors = FALSE
  )
  result <- lapply(seq_len(n_segments), function(i) {
    data[membership == i, , drop = FALSE]
  })
  names(result) <- list_names

  if (drop_empty) {
    keep <- counts > 0L
    result <- result[keep]
    info <- info[keep, , drop = FALSE]
    rownames(info) <- NULL
  }
  attr(result, "variable") <- variable
  attr(result, "cutpoints") <- cutpoints
  attr(result, "segment_info") <- info
  attr(result, "n_excluded") <- n_excluded
  attr(result, "call") <- match.call()
  class(result) <- c("cutpoint_split", "list")
  result
}

#' Print a fixed-cutpoint data split
#'
#' @param x A `cutpoint_split` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object, invisibly.
#' @export
print.cutpoint_split <- function(x, ...) {
  info <- attr(x, "segment_info")
  cat("\nData split by cutpoints\n")
  cat(sprintf("  Variable: %s\n", attr(x, "variable")))
  cat(sprintf("  Cutpoints: %s\n\n",
              paste(format(attr(x, "cutpoints"), trim = TRUE),
                    collapse = ", ")))
  if (nrow(info)) {
    width <- max(nchar(info$name), 9L)
    for (i in seq_len(nrow(info))) {
      cat(sprintf("  %-*s  %-20s  n = %d\n", width, info$name[i],
                  info$interval[i], info$n[i]))
    }
  } else {
    cat("  No non-empty intervals.\n")
  }
  n_excluded <- attr(x, "n_excluded")
  if (!is.null(n_excluded) && n_excluded > 0L)
    cat(sprintf("\n  Excluded non-finite values: %d\n", n_excluded))
  invisible(x)
}

#' Split a data frame by the levels of one categorical variable
#'
#' Split an input data frame into a named list of data frames before analysis.
#' Factor levels retain their declared order; character and logical values use
#' first-appearance order. The splitting column is retained in every returned
#' data frame.
#'
#' @param data A data frame.
#' @param variable A single character string naming the factor, character, or
#'   logical variable used to split the rows.
#' @param na_action How missing splitting values are handled: `"drop"`
#'   (default) excludes them and records the excluded count; `"error"` stops.
#' @param drop_empty Whether unused factor levels are removed. The default is
#'   `TRUE`. This argument has no effect for character or logical variables,
#'   whose groups are determined from observed non-missing values.
#'
#' @return A named list of data frames with class `factor_split`. List names
#'   are the original factor levels or categorical values. The `group_info`
#'   attribute records their order and row counts.
#'
#' @examples
#' d <- data.frame(
#'   batch = factor(c("Batch-A", "Batch-A", "Batch-B")),
#'   result = c(10.1, 9.9, 10.4)
#' )
#' parts <- split_by_factors(d, variable = "batch")
#' parts
#' parts[["Batch-A"]]
#' means <- lapply(parts, function(z) mean(z$result))
#' @export
split_by_factors <- function(data, variable,
                             na_action = c("drop", "error"),
                             drop_empty = TRUE) {
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(variable) || length(variable) != 1L ||
      is.na(variable) || !nzchar(variable) || !variable %in% names(data))
    stop("'variable' must name one column in 'data'.", call. = FALSE)
  values <- data[[variable]]
  if (!(is.factor(values) || is.character(values) || is.logical(values)))
    stop("'variable' must identify a factor, character, or logical column.",
         call. = FALSE)
  if (!is.logical(drop_empty) || length(drop_empty) != 1L ||
      is.na(drop_empty))
    stop("'drop_empty' must be TRUE or FALSE.", call. = FALSE)

  na_action <- match.arg(na_action)
  missing <- is.na(values)
  n_excluded <- sum(missing)
  if (n_excluded && na_action == "error")
    stop("The splitting variable contains ", n_excluded,
         " missing value(s); use na_action = \"drop\" to exclude them.",
         call. = FALSE)
  if (n_excluded) {
    data <- data[!missing, , drop = FALSE]
    values <- values[!missing]
  }

  observed <- as.character(values)
  if (is.factor(values)) {
    group_levels <- levels(values)
    if (drop_empty)
      group_levels <- group_levels[group_levels %in% observed]
  } else {
    group_levels <- unique(observed)
  }
  if (!length(group_levels))
    stop("No non-missing groups remain after applying the split rules.",
         call. = FALSE)

  counts <- vapply(group_levels, function(level) {
    sum(observed == level)
  }, integer(1))
  result <- lapply(group_levels, function(level) {
    data[observed == level, , drop = FALSE]
  })
  names(result) <- group_levels

  info <- data.frame(
    group = seq_along(group_levels),
    name = group_levels,
    level = group_levels,
    n = unname(counts),
    stringsAsFactors = FALSE
  )
  attr(result, "variable") <- variable
  attr(result, "group_info") <- info
  attr(result, "n_excluded") <- n_excluded
  attr(result, "call") <- match.call()
  class(result) <- c("factor_split", "list")
  result
}

#' Print a categorical data split
#'
#' @param x A `factor_split` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object, invisibly.
#' @export
print.factor_split <- function(x, ...) {
  info <- attr(x, "group_info")
  n_excluded <- attr(x, "n_excluded") %||% 0L
  cat("\nData split by factors\n")
  cat(sprintf("  Variable: %s\n", attr(x, "variable")))
  cat(sprintf("  Groups: %d\n", nrow(info)))
  cat(sprintf("  Missing/excluded: %d\n\n", n_excluded))
  display_name <- ifelse(nzchar(info$name), info$name, "<empty>")
  width <- max(nchar(display_name), 5L)
  for (i in seq_len(nrow(info))) {
    cat(sprintf("  %-*s  n = %d\n", width, display_name[i], info$n[i]))
  }
  invisible(x)
}

#' mcr constructor
#'
#' Create a method comparison regression (mcr) object containing one summarized
#' row per sample ID. Raw replicate rows must first be summarized with
#' [replicate_to_mean()] or by the user.
#'
#' @param data  Data frame with id, candidate, reference columns
#' @param id    ID variable name (character), used to identify samples
#' @param candidate Candidate/test method variable name (character), used as
#'   the response (Y) throughout regression and plotting.
#' @param reference Reference/comparative method variable name (character),
#'   used as the predictor (X) throughout regression and plotting.
#' @param weights Optional column name (character) in data containing non-negative finite
#'   numeric weights for weighted regression. Default NULL (no weighting).
#' @param candidate_sd,reference_sd Optional column names containing within-sample
#'   SDs for summarized replicate data. Supply both when `lambda = "replicates"`
#'   will be used with one mean pair per sample.
#' @param candidate_n,reference_n Optional column names containing replicate
#'   counts for summarized data. Supply both or neither. Without them, equal
#'   replicate counts are assumed and cancel from the variance ratio.
#'
#' @return Returns an S3 "mcr" object containing:
#'   \item{call}{Original function call}
#'   \item{data}{Complete data frame}
#'   \item{id / id_name}{ID vector and column name}
#'   \item{candidate / reference}{Test/reference method numeric vectors (complete pairs only)}
#'   \item{candidate_name / reference_name}{Method column names}
#'   \item{n}{Number of complete pairs}
#'   \item{complete_idx}{Row indices of complete pairs in the original data}
#'   and analysis result storage slots: $correlation, $regression,
#'   $outlier, $bland_altman
#'
#' @examples
#'   df <- data.frame(sid = 1:30,
#'                    test = rnorm(30, 50, 10),
#'                    ref  = rnorm(30, 50, 10))
#'   obj <- mcr(df, "sid", "test", "ref")
#' @export
mcr <- function(data, id, candidate, reference, weights = NULL,
                candidate_sd = NULL, reference_sd = NULL, candidate_n = NULL,
                reference_n = NULL) {

  # Input validation
  if (!is.data.frame(data))
    stop("'data' must be a data frame.")
  if (!is.character(id) || length(id) != 1L || !id %in% names(data))
    stop("'id' must be a single column name in 'data'.")
  if (!is.character(candidate) || length(candidate) != 1L || !candidate %in% names(data))
    stop("'candidate' must be a single column name in 'data'.")
  if (!is.character(reference) || length(reference) != 1L || !reference %in% names(data))
    stop("'reference' must be a single column name in 'data'.")

  # Extract data
  id_vec  <- data[[id]]
  cand    <- data[[candidate]]
  ref     <- data[[reference]]

  if (!is.numeric(cand))
    stop("'candidate' column must be numeric.")
  if (!is.numeric(ref))
    stop("'reference' column must be numeric.")

  summary_names <- list(candidate_sd = candidate_sd, reference_sd = reference_sd,
                        candidate_n = candidate_n, reference_n = reference_n)
  supplied <- !vapply(summary_names, is.null, logical(1))
  if (xor(supplied["candidate_sd"], supplied["reference_sd"]))
    stop("'candidate_sd' and 'reference_sd' must be supplied together.",
         call. = FALSE)
  if (xor(supplied["candidate_n"], supplied["reference_n"]))
    stop("'candidate_n' and 'reference_n' must be supplied together.",
         call. = FALSE)
  for (nm in names(summary_names)[supplied]) {
    column <- summary_names[[nm]]
    if (!is.character(column) || length(column) != 1L ||
        !column %in% names(data))
      stop("'", nm, "' must name a column in 'data'.", call. = FALSE)
    values <- data[[column]]
    if (!is.numeric(values) || any(!is.finite(values)))
      stop("'", nm, "' column must contain finite numeric values.",
           call. = FALSE)
    if (grepl("_sd$", nm) && any(values < 0))
      stop("'", nm, "' column must be non-negative.", call. = FALSE)
    if (grepl("_n$", nm) &&
        any(values < 2 | values != as.integer(values)))
      stop("'", nm, "' column must contain integer replicate counts >= 2.",
           call. = FALSE)
  }

  # Complete-pair analysis
  ok <- is.finite(cand) & is.finite(ref)
  n  <- sum(ok)
  n_total <- length(cand)
  n_miss  <- n_total - n

  if (n < 3L)
    stop("Need at least 3 complete pairs for method comparison, got ", n, ".")
  if (anyDuplicated(id_vec[!is.na(id_vec)]))
    stop("'mcr()' requires one summarized row per sample ID. Summarize raw ",
         "replicates first with replicate_to_mean() or provide precomputed means.",
         call. = FALSE)

  # Weights column validation
  if (!is.null(weights)) {
    if (!is.character(weights) || length(weights) != 1L || !weights %in% names(data))
      stop("'weights' must be a single column name in 'data' or NULL.", call. = FALSE)
    w_col <- data[[weights]]
    if (!is.numeric(w_col))
      stop("'weights' column must be numeric.", call. = FALSE)
    weights_vec <- w_col[ok]
    if (any(!is.finite(weights_vec) | weights_vec < 0))
      stop("'weights' column must contain non-negative finite values.", call. = FALSE)
  } else {
    weights_vec <- NULL
  }

  # Construct object
  # Print slot: data overview
  id_vec_all <- data[[id]]
  dups <- duplicated(id_vec_all) | duplicated(id_vec_all, fromLast = TRUE)
  n_dup <- sum(dups)
  dup_ids <- if (n_dup > 0) which(dups) else integer(0)

  miss_cand <- !is.finite(cand)
  miss_ref  <- !is.finite(ref)

  print_slot <- list(
    data_class   = class(data)[1L],
    n            = n,
    n_total      = n_total,
    n_miss       = n_miss,
    n_dup        = n_dup,
    var_names    = list(id = id, candidate = candidate, reference = reference),
    dup_ids      = dup_ids,
    miss_ids     = list(
      candidate  = which(miss_cand & !miss_ref),
      reference  = which(miss_ref & !miss_cand),
      both       = which(miss_cand & miss_ref)
    ),
    weights_name = weights
  )

  out <- list(
    call            = match.call(),
    data            = data,
    id              = id_vec[ok],
    id_name         = id,
    candidate       = cand[ok],
    reference       = ref[ok],
    candidate_name  = candidate,
    reference_name  = reference,
    n               = n,
    complete_idx    = which(ok),
    print           = print_slot,
    weights_name    = weights,
    weights         = weights_vec,
    replicate_columns = summary_names,
    candidate_sd    = if (supplied["candidate_sd"])
      data[[candidate_sd]][ok] else NULL,
    reference_sd    = if (supplied["reference_sd"])
      data[[reference_sd]][ok] else NULL,
    candidate_n     = if (supplied["candidate_n"])
      as.integer(data[[candidate_n]][ok]) else NULL,
    reference_n     = if (supplied["reference_n"])
      as.integer(data[[reference_n]][ok]) else NULL,
    describe        = NULL,
    correlation     = NULL,
    regression      = NULL,
    outlier         = NULL,
    bland_altman    = NULL,
    bias            = NULL
  )
  class(out) <- "mcr"

  invisible(out)
}


#' S3 print method
#'
#' Print data overview and analysis slot status based on \code{$print} slot.
#'
#' @param x \code{mcr} object
#' @param ... Additional arguments
#'
#' @return Invisible \code{mcr} object
#'
#' @export
#'
#' @examples
#'   df <- data.frame(sid = 1:30,
#'                    test = rnorm(30, 50, 10),
#'                    ref  = rnorm(30, 50, 10))
#'   obj <- mcr(df, "sid", "test", "ref")
#'   print(obj)
print.mcr <- function(x, ...) {

  if (!inherits(x, "mcr"))
    stop("Input must be an 'mcr' object.")

  p <- x$print %||% list(
    data_class = class(x$data)[1L],
    n          = x$n,
    n_total    = nrow(x$data),
    n_miss     = nrow(x$data) - x$n,
    n_dup      = 0L,
    var_names  = list(id = x$id_name, candidate = x$candidate_name,
                      reference = x$reference_name),
    dup_ids    = integer(0),
    miss_ids   = list(candidate = integer(0), reference = integer(0), both = integer(0))
  )

  cat("\nMethod Comparison\n")
  cat(sprintf("  Data class:        %s\n", p$data_class))
  cat(sprintf("  Total rows:        %d\n", p$n_total))
  cat(sprintf("  Complete pairs:    %d\n", p$n))
  if (p$n_miss > 0)
    cat(sprintf("  Missing/excluded:  %d\n", p$n_miss))
  cat(sprintf("  ID variable:       %s\n", p$var_names$id))
  cat(sprintf("  Reference (X):     %s\n", p$var_names$reference))
  cat(sprintf("  Candidate (Y):     %s\n", p$var_names$candidate))
  if (p$n_dup > 0)
    cat(sprintf("  Duplicate IDs:     %d (rows: %s)\n",
                p$n_dup, paste(p$dup_ids, collapse = ", ")))
  else
    cat("  Duplicate IDs:     none\n")

  if (!is.null(p$weights_name))
    cat(sprintf("  Weights column:   %s\n", p$weights_name))

  # Missing detail
  if (p$n_miss > 0) {
    cat("  Missing value detail:\n")
    if (length(p$miss_ids$candidate))
      cat(sprintf("    %s: %d row(s): %s\n",
                  p$var_names$candidate, length(p$miss_ids$candidate),
                  paste(p$miss_ids$candidate, collapse = ", ")))
    if (length(p$miss_ids$reference))
      cat(sprintf("    %s: %d row(s): %s\n",
                  p$var_names$reference, length(p$miss_ids$reference),
                  paste(p$miss_ids$reference, collapse = ", ")))
    if (length(p$miss_ids$both))
      cat(sprintf("    both: %d row(s): %s\n",
                  length(p$miss_ids$both),
                  paste(p$miss_ids$both, collapse = ", ")))
  }

  # Analysis slot status
  cat("\nAnalysis slots:\n")
  slots <- list(
    describe    = "Describe",
    correlation = "Correlation",
    regression  = "Regression",
    outlier     = "Outlier",
    bland_altman = "Bland-Altman",
    bias        = "Bias"
  )
  for (nm in names(slots)) {
    val <- x[[nm]]
    if (!is.null(val)) {
      extra <- ""
      if (nm == "correlation") extra <- sprintf(" (%s)", val$method)
      else if (nm == "regression") extra <- sprintf(" (%s)", val$method)
      else if (nm == "outlier") extra <- sprintf(" (%s)", val$method)
      else if (nm == "bland_altman") extra <- sprintf(" (%s)", val$type)
      cat(sprintf("    [x] %s%s\n", slots[[nm]], extra))
    } else {
      cat(sprintf("    [ ] %s\n", slots[[nm]]))
    }
  }

  invisible(x)
}


#' S3 summary method
#'
#' Summarize key information from all executed analyses in the \code{mcr} object
#' (descriptive statistics, correlation analysis, regression analysis,
#' outlier detection, Bland-Altman analysis, and bias analysis).
#'
#' @param object \code{mcr} object
#' @param ... Additional arguments
#'
#' @return Invisible \code{mcr} object
#'
#' @examples
#'   df <- data.frame(sid = 1:30,
#'                    test = rnorm(30, 50, 10),
#'                    ref  = rnorm(30, 50, 10))
#'   obj <- mcr(df, "sid", "test", "ref")
#'   summary(obj)
#'   obj <- correlation(obj)
#'   obj <- regression(obj, method = "ols")
#'   obj <- bland_altman(obj)
#'   summary(obj)
#' @export
summary.mcr <- function(object, ...) {

  if (!inherits(object, "mcr"))
    stop("Input must be an 'mcr' object.")

  cat("\nMethod Comparison Regression (mcr)\n")
  cat("-----------------------------------------------------\n")
  print.mcr(object)


  # Detail: Describe
  if (!is.null(object$describe)) print(object$describe)

  # Detail: Correlation
  if (!is.null(object$correlation)) print(object$correlation)

  # Detail: Regression
  if (!is.null(object$regression)) print(object$regression)

  # Detail: Outlier
  if (!is.null(object$outlier)) print(object$outlier)

  # Detail: Bland-Altman
  if (!is.null(object$bland_altman)) print(object$bland_altman)

  # Detail: Bias
  if (!is.null(object$bias)) print(object$bias)

  invisible(object)
}


#' Print method for describe results
#' @param x mcr_describe object
#' @param digits Number of decimal places, default 4
#' @param ... additional arguments
#' @return The unchanged \code{mcr_describe} object, invisibly.
#' @export
print.mcr_describe <- function(x, digits = NULL, ...) {
  d <- x
  if (is.null(digits)) digits <- 4L

  # Method comparison table + other numeric columns
  cat("\nDescriptive Statistics\n")

  # Three rows for method comparison
  tbl_rows <- list(
    list(name = sprintf("  %s (X)", d$reference_name), vals = d$reference_values),
    list(name = sprintf("  %s (Y)", d$candidate_name), vals = d$candidate_values),
    list(name = "  Diff (Y-X)",              vals = d$diff_values)
  )

  # Add other numeric columns to the main table
  cat_cols <- list()
  for (col in names(d$columns)) {
    ci <- d$columns[[col]]
    if (ci$type == "numeric") {
      tbl_rows[[length(tbl_rows) + 1]] <- list(name = sprintf("  %s", col), vals = d$numeric_extra[[col]])
    } else if (ci$type == "categorical") {
      cat_cols[[col]] <- ci
    }
  }

  # Table header
  sep75 <- paste(rep("-", 75), collapse = "")
  cat(sprintf("  %-18s %5s %8s %8s %8s %8s %8s %8s %8s  %s\n",
              "Variable", "n", "Mean", "SD", "Min", "Q1", "Median", "Q3", "Max", "NA"))
  cat(" ", sep75, "\n")

  for (m in tbl_rows) {
    v <- m$vals
    n <- length(v)
    n_na <- sum(is.na(v))
    if (n == 0) {
      cat(sprintf("  %-18s %5d %8s %8s %8s %8s %8s %8s %8s  %s\n",
                  m$name, 0L, "-", "-", "-", "-", "-", "-", "-", if (n_na > 0) n_na else ""))
    } else {
      ok <- v[!is.na(v)]
      if (length(ok) == 0) {
        cat(sprintf("  %-18s %5d %8s %8s %8s %8s %8s %8s %8s  %d\n",
                    m$name, n, "-", "-", "-", "-", "-", "-", "-", n_na))
      } else {
        na_str <- if (n_na > 0) sprintf("%d", n_na) else ""
        cat(sprintf("  %-18s %5d %8s %8s %8s %8s %8s %8s %8s  %s\n",
                    m$name, length(ok),
                    .format_num(mean(ok), digits),
                    .format_num(sd(ok),   digits),
                    .format_num(min(ok),  digits),
                    .format_num(quantile(ok, 0.25, names = FALSE), digits),
                    .format_num(median(ok), digits),
                    .format_num(quantile(ok, 0.75, names = FALSE), digits),
                    .format_num(max(ok),  digits),
                    na_str))
      }
    }
  }

  # Categorical table
  for (col in names(cat_cols)) {
    ci <- cat_cols[[col]]
    na_str <- if (ci$n_na > 0)
      sprintf(" (%d NA, %.1f%%)", ci$n_na, ci$pct_na) else ""

    cat(sprintf("\n  %s%s\n", col, na_str))
    cat(sprintf("  %-20s %5s  %7s\n", "Level", "Count", "Percent"))
    cat("  ", paste(rep("-", 38), collapse = ""), "\n")
    for (r in seq_len(nrow(ci$table))) {
      lv <- as.character(ci$table$Level[r])
      if (is.na(lv)) lv <- "<NA>"
      cat(sprintf("  %-20s %5d  %6.1f%%\n",
                  lv, ci$table$Freq[r], ci$table$Pct[r]))
    }
  }

  invisible(x)
}

#' Print method for correlation results
#' @param x mcr_correlation object
#' @param ... additional arguments
#' @return The unchanged \code{mcr_correlation} object, invisibly.
#' @export
print.mcr_correlation <- function(x, ...) {
  cl_pct <- if (!is.null(x$conf.int))
    sprintf(", %.0f%% CI [%.4f, %.4f]",
            attr(x$conf.int, "conf.level") * 100,
            x$conf.int[1L], x$conf.int[2L]) else ""
  cat("\nCorrelation\n")
  cat(sprintf("  Method: %s\n", x$method))
  cat(sprintf("  r = %.4f%s\n", x$estimate, cl_pct))
  cat(sprintf("  p-value: %s\n", format.pval(x$p.value, digits = 4)))
  cat(sprintf("  n = %d\n", x$n))
  invisible(x)
}

#' Print method for regression results
#' @param x mcr_regression object
#' @param ... additional arguments
#' @return The unchanged \code{mcr_regression} object, invisibly.
#' @export
print.mcr_regression <- function(x, ...) {
  cat("\nRegression\n")
  cat(sprintf("  Method: %s\n", x$method))
  if (!is.null(x$weights_spec))
    cat(sprintf("  Weights: %s\n", x$weights_spec))
  if (!is.null(x$weight_iterations))
    cat(sprintf("  Weight updates: %d\n", x$weight_iterations))
  if (!is.null(x$ci_method))
    cat(sprintf("  CI method: %s\n", x$ci_method))
  if (!is.null(x$lambda_metadata))
    cat(sprintf("  Lambda: %s (from %s)\n", .format_num(x$lambda),
                x$lambda_metadata$source))
  if (!is.null(x$bootstrap))
    cat(sprintf("  Bootstrap: %d effective / %d requested (%d failed)\n",
                x$bootstrap$effective, x$bootstrap$requested,
                x$bootstrap$failed))
  cat(sprintf("  Intercept = %s  (SE = %s)\n",
              .format_num(x$intercept), .format_num(x$intercept_se)))
  cat(sprintf("  Slope     = %s  (SE = %s)\n",
              .format_num(x$slope), .format_num(x$slope_se)))
  cl_pct <- x$conf.level * 100
  cat(sprintf("  Intercept %.0f%% CI: [%s, %s]\n", cl_pct,
              .format_num(x$intercept_ci[1L]), .format_num(x$intercept_ci[2L])))
  cat(sprintf("  Slope     %.0f%% CI: [%s, %s]\n", cl_pct,
              .format_num(x$slope_ci[1L]), .format_num(x$slope_ci[2L])))
  if (!is.na(x$slope_se) && is.finite(x$slope_se)) {
    t_slope1 <- (x$slope - 1) / x$slope_se
    p_slope1 <- 2 * pt(-abs(t_slope1), df = x$n - 2)
    ci_contains_1 <- x$slope_ci[1L] <= 1 && x$slope_ci[2L] >= 1
    cat(sprintf("    H0: slope = 1 -> t = %.4f, p = %s\n",
                t_slope1, format.pval(p_slope1, digits = 4)))
    if (p_slope1 >= 1 - x$conf.level) {
      cat("    slope not significantly different from 1",
          if (ci_contains_1) " (CI contains 1)" else "", "\n", sep = "")
    } else {
      cat("    slope significantly different from 1",
          if (!ci_contains_1) " (CI excludes 1)" else "", "\n", sep = "")
    }
  }
  if (!is.null(x$r_squared))
    cat(sprintf("  R-squared = %.4f\n", x$r_squared))
  cat(sprintf("  Sigma     = %s\n", .format_num(x$sigma)))
  cat(sprintf("  n = %d\n", x$n))
  invisible(x)
}

#' Print method for outlier results
#' @param x mcr_outlier object
#' @param ... additional arguments
#' @return The unchanged \code{mcr_outlier} object, invisibly.
#' @export
print.mcr_outlier <- function(x, ...) {
  cat("\nOutlier\n")
  cat(sprintf("  Method: %s\n", x$method))
  cat(sprintf("  Data:   %s\n", x$data_label))
  cat(sprintf("  Alpha:  %.2f", x$alpha))
  if (x$method == "iqr")
    cat(sprintf("  (coef = %.1f)", x$coef))
  cat("\n\n")
  if (x$n_out == 0L) {
    cat("  No outliers detected.\n")
  } else {
    cat(sprintf("  %d outlier(s) found:\n", x$n_out))
    cat(sprintf("  %-12s  %10s  %10s  %s\n", "Row", "Value", "Statistic", "Critical"))
    cat("  ", paste(rep("-", 50), collapse = ""), "\n")
    for (i in seq_len(x$n_out)) {
      s <- if (length(x$statistic) >= i) x$statistic[i] else x$statistic
      c <- if (length(x$critical) >= i) x$critical[i] else x$critical
      cat(sprintf("  %-12d  %10.4f  %10.4f  %.4f\n",
                  x$indices[i], x$values[i], s, c))
    }
  }

  invisible(x)
}

#' Print method for Bland-Altman results
#' @param x mcr_bland_altman object
#' @param ... additional arguments
#' @return The unchanged \code{mcr_bland_altman} object, invisibly.
#' @export
print.mcr_bland_altman <- function(x, ...) {
  agree_pct <- x$agree.level * 100
  conf_pct  <- x$conf.level * 100
  cat("\nBland-Altman\n")
  cat(sprintf("  Type:        %s\n", x$type))
  cat(sprintf("  Y-axis:      %s\n", x$y_label))
  cat(sprintf("  X-axis:      %s\n", x$x_label))
  if (!is.null(x$percent_denominator))
    cat(sprintf("  Percent base: %s\n", x$percent_denominator))
  if (!is.null(x$method)) cat(sprintf("  Method:       %s\n", x$method))
  cat(sprintf("  n = %d\n\n", x$n))
  center_label <- switch(x$method,
    nonparametric = "Median diff",
    hodges_lehmann = "HL pseudomedian",
    "Mean diff")
  cat(sprintf("  %-16s %s\n", paste0(center_label, ":"), .format_num(x$center)))
  if (!is.null(x$center_ci)) {
    if (identical(x$method, "nonparametric") &&
        !is.null(x$center_ci_coverage)) {
      cat(sprintf("  %.0f%% requested CI for %s (actual %.1f%%): [%s, %s]\n",
                  conf_pct, tolower(center_label), x$center_ci_coverage * 100,
                  .format_num(x$center_ci[1L]), .format_num(x$center_ci[2L])))
    } else {
      cat(sprintf("  %.0f%% CI for %s: [%s, %s]\n", conf_pct,
                  tolower(center_label), .format_num(x$center_ci[1L]),
                  .format_num(x$center_ci[2L])))
    }
  }
  cat(sprintf("  SD diff:         %s\n", .format_num(x$sd_diff)))
  cat(sprintf("  %.0f%% LoA: [%s, %s]\n", agree_pct,
              .format_num(x$loa[1L]), .format_num(x$loa[2L])))
  cat(sprintf("\n  %.0f%% CI for LoA:\n", conf_pct))
  cat(sprintf("    Lower LoA: [%s, %s]\n",
              .format_num(x$loa_ci$lower[1L]), .format_num(x$loa_ci$lower[2L])))
  cat(sprintf("    Upper LoA: [%s, %s]\n",
              .format_num(x$loa_ci$upper[1L]), .format_num(x$loa_ci$upper[2L])))
  invisible(x)
}

#' Print method for bias results
#' @param x mcr_bias object
#' @param ... additional arguments
#' @return The unchanged \code{mcr_bias} data frame, invisibly.
#' @export
print.mcr_bias <- function(x, ...) {
  cat("\nBias\n")
  cat(sprintf("  Method: %s\n", attr(x, "method") %||% "N/A"))
  interval <- attr(x, "interval") %||% ""
  level    <- attr(x, "level") %||% 0.95
  if (nchar(interval) > 0)
    cat(sprintf("  Interval: %s (%.0f%%)\n", interval, level * 100))
  else
    cat("  Interval: point estimate only\n")
  cat(sprintf("  MDL points: %d\n", nrow(x)))
  cat(sprintf("  Bias range: [%s, %s]\n",
              .format_num(min(x$bias)), .format_num(max(x$bias))))
  cat("\n")
  txt <- capture.output(print.data.frame(x, row.names = FALSE))
  cat(txt[1L], "\n")
  cat(paste(rep("-", nchar(txt[1L])), collapse = ""), "\n")
  cat(paste(txt[-1L], collapse = "\n"), "\n")
  invisible(x)
}
#'
#' Compute descriptive statistics and store results
#'
#' Computes data overview (rows/columns/ID duplicates), per-column summaries,
#' and a side-by-side comparison table of candidate, reference, and Diff.
#' Results are stored in \code{x$describe} (class \code{mcr_describe}) and
#' displayed via \code{print.mcr_describe()} (e.g., when calling \code{summary()}).
#'
#' @param x       mcr object
#' @param cols    Column name vector to analyze; by default, analyzes all columns except
#'                id, candidate, and reference. Supports exclusion with \code{-} prefix,
#'                e.g., cols = -c("age", "sex").
#' @param digits  Number of decimal places, default 4
#' @param ...     Additional arguments
#'
#' @return Updated \code{mcr} object (invisible), results stored in \code{x$describe}
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   describe(obj)
#' @export
describe.mcr <- function(x, cols = NULL, digits = 4L, ...) {

  if (!inherits(x, "mcr"))
    stop("Input must be an 'mcr' object.")

  data <- x$data
  n_col <- ncol(data)
  all_names <- names(data)

  # Determine which columns to analyze
  # By default, exclude id, candidate, reference columns
  exclude_default <- c(x$id_name, x$candidate_name, x$reference_name)

  cols_expr <- substitute(cols)
  negated <- FALSE
  if (is.call(cols_expr) &&
      length(cols_expr) == 2 && identical(cols_expr[[1]], as.name("-"))) {
    negated <- TRUE
    cols <- eval(cols_expr[[2]], parent.frame())
  }

  if (is.null(cols)) {
    cols <- setdiff(all_names, exclude_default)
  } else if (negated) {
    cols <- setdiff(all_names, cols)
  } else {
    include <- grep("^-", cols, value = TRUE, invert = TRUE)
    exclude <- grep("^-", cols, value = TRUE)
    exclude <- sub("^-", "", exclude)
    if (length(include) > 0) {
      bad <- setdiff(include, all_names)
      if (length(bad))
        stop("Column(s) not found: ", paste(bad, collapse = ", "))
      cols <- include
    } else {
      cols <- setdiff(all_names, exclude_default)
    }
    if (length(exclude) > 0) {
      bad <- setdiff(exclude, all_names)
      if (length(bad))
        stop("Column(s) to exclude not found: ", paste(bad, collapse = ", "))
      cols <- setdiff(cols, exclude)
    }
  }
  # Always exclude id
  cols <- setdiff(cols, x$id_name)

  # ID duplicate check (based on complete pairs)
  id_dup <- NULL
  id_vec <- x$id
  dups <- duplicated(id_vec) | duplicated(id_vec, fromLast = TRUE)
  n_dup <- sum(dups)
  n_unique <- length(unique(id_vec))
  id_dup <- list(n_dup = n_dup, n_unique = n_unique)

  # Per-column summary (based on original data frame)
  columns <- list()

  for (col in cols) {
    vec     <- data[[col]]
    n_total <- length(vec)
    n_na    <- sum(is.na(vec))
    pct_na  <- round(n_na / n_total * 100, 1)

    if (is.factor(vec) || is.character(vec)) {
      tbl   <- table(vec, useNA = "ifany")
      n_ok  <- sum(!is.na(vec))
      freqs <- as.data.frame(tbl, responseName = "Freq")
      names(freqs)[1] <- "Level"
      freqs$Pct <- round(freqs$Freq / n_ok * 100, 1)
      na_row <- which(is.na(freqs$Level))
      if (length(na_row))
        freqs$Pct[na_row] <- round(freqs$Freq[na_row] / n_total * 100, 1)

      columns[[col]] <- list(
        type     = "categorical",
        n_na     = n_na,
        pct_na   = pct_na,
        n_levels = nrow(freqs),
        table    = freqs
      )

    } else if (is.numeric(vec)) {
      ok  <- vec[!is.na(vec)]
      if (length(ok) == 0) {
        columns[[col]] <- list(
          type = "numeric", n_na = n_na, pct_na = pct_na,
          n = 0L, mean = NA_real_, sd = NA_real_,
          min = NA_real_, q25 = NA_real_, median = NA_real_,
          q75 = NA_real_, max = NA_real_
        )
      } else {
        columns[[col]] <- list(
          type   = "numeric",
          n_na   = n_na,
          pct_na = pct_na,
          n      = length(ok),
          mean   = mean(ok),
          sd     = sd(ok),
          min    = min(ok),
          q25    = quantile(ok, 0.25, names = FALSE),
          median = median(ok),
          q75    = quantile(ok, 0.75, names = FALSE),
          max    = max(ok)
        )
      }
    } else {
      columns[[col]] <- list(type = "other", n_na = n_na, pct_na = pct_na)
    }
  }

  # Gather complete-pair values for numeric extra columns (for print method)
  numeric_table_cols <- list()
  for (col in names(columns)) {
    ci <- columns[[col]]
    if (ci$type == "numeric") {
      numeric_table_cols[[col]] <- data[[col]][x$complete_idx]
    }
  }

  result <- structure(list(
    n_complete       = x$n,
    n_total          = nrow(data),
    n_col            = n_col,
    id_dup           = id_dup,
    columns          = columns,
    id_name          = x$id_name,
    candidate_name   = x$candidate_name,
    reference_name   = x$reference_name,
    candidate_values = x$candidate,
    reference_values = x$reference,
    diff_values      = x$candidate - x$reference,
    numeric_extra    = numeric_table_cols
  ), class = "mcr_describe")
  x$describe <- result
  print(result)
  invisible(x)
}


#' Correlation analysis (S3 method)
#'
#' Perform correlation analysis between the reference method (X) and candidate method (Y),
#' supporting Pearson, Spearman, and Kendall methods.
#'
#' @param x       mcr object
#' @param method  Correlation method: "pearson" (default), "spearman", or "kendall"
#' @param ...     Additional arguments passed to cor.test
#'
#' @return Updated mcr object (invisible), results stored in $correlation
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   correlation(obj)
#'   correlation(obj, method = "spearman")
#' @export
correlation.mcr <- function(x, method = c("pearson", "spearman", "kendall"), ...) {

  if (!inherits(x, "mcr"))
    stop("Input must be an 'mcr' object.")

  method <- match.arg(method)

  # Computation
  ct <- cor.test(x$candidate, x$reference, method = method, ...)

  # Storage
  result <- structure(list(
    method    = method,
    estimate  = as.numeric(ct$estimate),
    statistic = as.numeric(ct$statistic),
    parameter = as.numeric(ct$parameter),
    p.value   = ct$p.value,
    conf.int  = ct$conf.int,
    n         = x$n
  ), class = "mcr_correlation")
  x$correlation <- result

  print(result)
  invisible(x)
}


# Internal function: Deming regression (closed-form, Linnet 1993)
.deming_fit <- function(x, y, lambda = 1, conf.level = 0.95) {

  n <- length(x)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      is.na(lambda) || lambda <= 0)
    stop("'lambda' must be a single positive number.", call. = FALSE)
  mx <- mean(x); my <- mean(y)

  s_xx <- sum((x - mx)^2)
  s_yy <- sum((y - my)^2)
  s_xy <- sum((x - mx) * (y - my))

  if (s_xy == 0)
    stop("Zero covariance between x and y; cannot fit Deming regression.")

  # Deming slope & intercept (Linnet 1993, eq. 5)
  U <- s_yy - lambda * s_xx
  slope <- (U + sqrt(U^2 + 4 * lambda * s_xy^2)) / (2 * s_xy)
  intercept <- my - slope * mx

  # Jackknife SE
  jack_slope <- numeric(n)
  jack_int   <- numeric(n)
  for (i in seq_len(n)) {
    xi <- x[-i]; yi <- y[-i]
    mi_x <- mean(xi); mi_y <- mean(yi)
    sxx_i <- sum((xi - mi_x)^2)
    syy_i <- sum((yi - mi_y)^2)
    sxy_i <- sum((xi - mi_x) * (yi - mi_y))
    if (sxy_i == 0) next
    U_i <- syy_i - lambda * sxx_i
    jack_slope[i] <- (U_i + sqrt(U_i^2 + 4 * lambda * sxy_i^2)) / (2 * sxy_i)
    jack_int[i]   <- mi_y - jack_slope[i] * mi_x
  }
  slope_se <- sqrt((n - 1) / n * sum((jack_slope - mean(jack_slope, na.rm = TRUE))^2, na.rm = TRUE))
  int_se   <- sqrt((n - 1) / n * sum((jack_int   - mean(jack_int,   na.rm = TRUE))^2, na.rm = TRUE))
  cov_xy   <- (n - 1) / n * sum((jack_slope - mean(jack_slope, na.rm = TRUE)) *
                                 (jack_int   - mean(jack_int,   na.rm = TRUE)),
                                 na.rm = TRUE)

  # Residuals and sigma
  predicted <- intercept + slope * x
  residuals <- y - predicted
  sigma <- sqrt(sum(residuals^2) / (n - 2))

  t_val <- qt(1 - (1 - conf.level) / 2, n - 2)

  list(
    method        = "Deming",
    intercept     = intercept,
    slope         = slope,
    intercept_se  = int_se,
    slope_se      = slope_se,
    cov_xy        = cov_xy,
    intercept_ci  = c(intercept - t_val * int_se, intercept + t_val * int_se),
    slope_ci      = c(slope - t_val * slope_se, slope + t_val * slope_se),
    sigma         = sigma,
    residuals     = residuals,
    lambda        = lambda,
    n             = n,
    conf.level    = conf.level
  )
}


# Internal function: Weighted Deming regression
.wdeming_fit <- function(x, y, lambda = 1, weights = NULL, conf.level = 0.95) {

  n <- length(x)
  if (is.null(weights))
    weights <- rep(1, n)

  if (!is.numeric(weights) || length(weights) != n)
    stop("'weights' must be a numeric vector of length ", n,
         " (number of complete pairs), got ", length(weights), ".", call. = FALSE)
  if (any(!is.finite(weights) | weights < 0))
    stop("'weights' must be non-negative and finite.", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      is.na(lambda) || lambda <= 0)
    stop("'lambda' must be a single positive number.", call. = FALSE)

  w <- weights / sum(weights) * n  # Normalize
  mx <- sum(w * x) / n
  my <- sum(w * y) / n

  s_xx <- sum(w * (x - mx)^2)
  s_yy <- sum(w * (y - my)^2)
  s_xy <- sum(w * (x - mx) * (y - my))

  if (s_xy == 0)
    stop("Zero weighted covariance between x and y.")

  U <- s_yy - lambda * s_xx
  slope <- (U + sqrt(U^2 + 4 * lambda * s_xy^2)) / (2 * s_xy)
  intercept <- my - slope * mx

  # Jackknife SE
  jack_slope <- numeric(n)
  jack_int   <- numeric(n)
  for (i in seq_len(n)) {
    xi <- x[-i]; yi <- y[-i]; wi <- weights[-i]
    wi <- wi / sum(wi) * (n - 1)
    mi_x <- sum(wi * xi) / (n - 1)
    mi_y <- sum(wi * yi) / (n - 1)
    sxx_i <- sum(wi * (xi - mi_x)^2)
    syy_i <- sum(wi * (yi - mi_y)^2)
    sxy_i <- sum(wi * (xi - mi_x) * (yi - mi_y))
    if (sxy_i == 0) next
    U_i <- syy_i - lambda * sxx_i
    jack_slope[i] <- (U_i + sqrt(U_i^2 + 4 * lambda * sxy_i^2)) / (2 * sxy_i)
    jack_int[i]   <- mi_y - jack_slope[i] * mi_x
  }
  slope_se <- sqrt((n - 1) / n * sum((jack_slope - mean(jack_slope, na.rm = TRUE))^2, na.rm = TRUE))
  int_se   <- sqrt((n - 1) / n * sum((jack_int   - mean(jack_int,   na.rm = TRUE))^2, na.rm = TRUE))
  cov_xy   <- (n - 1) / n * sum((jack_slope - mean(jack_slope, na.rm = TRUE)) *
                                 (jack_int   - mean(jack_int,   na.rm = TRUE)),
                                 na.rm = TRUE)

  predicted <- intercept + slope * x
  residuals <- y - predicted
  sigma <- sqrt(sum(w * residuals^2) / (n - 2))

  t_val <- qt(1 - (1 - conf.level) / 2, n - 2)

  list(
    method        = "Weighted Deming",
    intercept     = intercept,
    slope         = slope,
    intercept_se  = int_se,
    slope_se      = slope_se,
    cov_xy        = cov_xy,
    intercept_ci  = c(intercept - t_val * int_se, intercept + t_val * int_se),
    slope_ci      = c(slope - t_val * slope_se, slope + t_val * slope_se),
    sigma         = sigma,
    residuals     = residuals,
    lambda        = lambda,
    n             = n,
    conf.level    = conf.level
  )
}


# Internal point-estimate helpers used by jackknife and bootstrap refits.
.deming_point <- function(x, y, lambda = 1, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  w <- weights / sum(weights)
  mx <- sum(w * x); my <- sum(w * y)
  sxx <- sum(w * (x - mx)^2)
  syy <- sum(w * (y - my)^2)
  sxy <- sum(w * (x - mx) * (y - my))
  if (!is.finite(sxy) || abs(sxy) <= .Machine$double.eps)
    stop("Zero covariance between x and y; cannot fit Deming regression.",
         call. = FALSE)
  u <- syy - lambda * sxx
  slope <- (u + sqrt(u^2 + 4 * lambda * sxy^2)) / (2 * sxy)
  c(intercept = my - slope * mx, slope = slope)
}

.constant_cv_point <- function(x, y, lambda, epsilon = 1e-10) {
  if (any(!is.finite(x) | !is.finite(y) | x <= 0 | y <= 0))
    stop("Constant-CV Deming requires strictly positive finite concentrations ",
         "for both methods.", call. = FALSE)
  .require_pkg("ppwdeming")
  # ppwdeming defines lambda as Var(X)/Var(Y); ivdtools consistently exposes
  # Var(Y)/Var(X), hence the reciprocal at this single package boundary.
  raw <- ppwdeming::WD_Linnet(
    X = x, Y = y, lambda = 1 / lambda, getCI = FALSE, epsilon = epsilon
  )
  if (any(!is.finite(c(raw$alpha, raw$beta))))
    stop("Constant-CV Deming did not return finite coefficients.", call. = FALSE)

  # Reconstruct Linnet's final latent readings and weights for diagnostics.
  latent_x <- x; latent_y <- y; iter <- 0L; rel_change <- Inf
  pkg_lambda <- 1 / lambda
  repeat {
    iter <- iter + 1L
    old_x <- latent_x; old_y <- latent_y
    w <- ((latent_x + latent_y) / 2)^(-2)
    sw <- sum(w); xb <- sum(w * x) / sw; yb <- sum(w * y) / sw
    u <- sum(w * (x - xb)^2); q <- sum(w * (y - yb)^2)
    p <- sum(w * (x - xb) * (y - yb))
    disc <- (u - pkg_lambda * q)^2 + 4 * pkg_lambda * p^2
    b <- (pkg_lambda * q - u + sqrt(disc)) / (2 * pkg_lambda * p)
    a <- yb - b * xb
    d <- y - a - b * x; den <- 1 + pkg_lambda * b^2
    latent_x <- x + pkg_lambda * b * d / den
    latent_y <- y - d / den
    rel_change <- sum((old_x - latent_x)^2 + (old_y - latent_y)^2) /
      sum(old_x^2 + old_y^2)
    if (!is.finite(rel_change) || rel_change <= epsilon || iter >= 10000L) break
  }
  list(
    coefficients = c(intercept = raw$alpha, slope = raw$beta),
    latent_reference = latent_x,
    latent_candidate = latent_y,
    latent_concentration = (latent_x + latent_y) / 2,
    weights = ((latent_x + latent_y) / 2)^(-2),
    convergence = list(converged = is.finite(rel_change) && rel_change <= epsilon,
                       iterations = iter, relative_change = rel_change,
                       tolerance = epsilon)
  )
}

.profile_variance_function <- function(model, component, label) {
  if (is.function(model)) {
    fun <- function(z) {
      value <- model(z)
      if (!is.numeric(value) || length(value) != length(z) ||
          any(!is.finite(value) | value <= 0))
        stop(label, " variance function must return one positive finite ",
             "variance per concentration.", call. = FALSE)
      as.numeric(value)
    }
    return(list(fun = fun, range = NULL, source = "function"))
  }
  if (!inherits(model, "precision") || is.null(model$profile$fits))
    stop(label, " variance model must be a function or a precision object ",
         "after profile().", call. = FALSE)
  fit <- model$profile$fits[[component]]
  if (is.null(fit) || is.null(fit$vfp) || is.null(fit$best_no))
    stop("Precision profile component '", component, "' is unavailable for ",
         label, ".", call. = FALSE)
  observed_range <- range(fit$data$mean, finite = TRUE)
  fun <- function(z) {
    prediction <- suppressWarnings(VFP::predict.VFP(
      fit$vfp, model.no = fit$best_no, newdata = z, type = "sd"
    ))
    value <- as.numeric(prediction$Fitted)^2
    if (length(value) != length(z) || any(!is.finite(value) | value <= 0))
      stop(label, " precision profile produced non-positive or non-finite ",
           "variance predictions.", call. = FALSE)
    value
  }
  list(fun = fun, range = observed_range,
       source = paste0("precision:", component))
}

.precision_profile_point <- function(x, y, reference_fun, candidate_fun,
                                     epsilon = 1e-8) {
  .require_pkg("ppwdeming")
  raw <- ppwdeming::PWD_known(
    X = x, Y = y,
    gfun = function(z, parms) reference_fun(z),
    hfun = function(z, parms) candidate_fun(z),
    gparms = numeric(0), hparms = numeric(0), epsilon = epsilon,
    getCI = FALSE
  )
  if (any(!is.finite(c(raw$alpha, raw$beta))))
    stop("Precision-profile Deming did not return finite coefficients.",
         call. = FALSE)
  list(
    coefficients = c(intercept = raw$alpha, slope = raw$beta),
    latent_concentration = raw$mu,
    fitted_candidate = raw$fity,
    reference_variance = reference_fun(raw$mu),
    candidate_variance = candidate_fun(raw$fity),
    scaled_residuals = raw$scalr,
    convergence = list(converged = all(is.finite(c(raw$mu, raw$fity, raw$scalr))),
                       criterion = raw$L, tolerance = epsilon)
  )
}

.replicate_lambda_data <- function(x) {
  if (!is.null(x$candidate_sd)) {
    cn <- x$candidate_n; rn <- x$reference_n
    if (xor(is.null(cn), is.null(rn)))
      stop("Both candidate and reference replicate counts are required.",
           call. = FALSE)
    if (!is.null(cn) && (length(unique(cn)) != 1L || length(unique(rn)) != 1L))
      stop("Replicate counts must be constant across samples for scalar ",
           "lambda; use weights = \"precision_profile\" for varying counts.",
           call. = FALSE)
    if (is.null(cn)) {
      candidate_var <- mean(x$candidate_sd^2)
      reference_var <- mean(x$reference_sd^2)
      mean_ratio <- 1
      counts <- list(candidate = NA_integer_, reference = NA_integer_)
      pooling <- "equal-replicate-count assumption"
    } else {
      candidate_var <- sum((cn - 1) * x$candidate_sd^2) / sum(cn - 1)
      reference_var <- sum((rn - 1) * x$reference_sd^2) / sum(rn - 1)
      mean_ratio <- unique(rn) / unique(cn)
      counts <- list(candidate = unique(cn), reference = unique(rn))
      pooling <- "degrees-of-freedom pooled"
    }
    lambda <- candidate_var / reference_var * mean_ratio
    if (!is.finite(lambda) || lambda <= 0)
      stop("Pooled candidate and reference replicate variances must both be ",
           "positive and finite to estimate lambda.", call. = FALSE)
    return(list(x = x$reference, y = x$candidate, id = x$id,
                complete_idx = x$complete_idx, lambda = lambda,
                metadata = list(source = "summarized SD data", pooling = pooling,
                  replicate_counts = counts,
                  pooled_variance = c(reference = reference_var,
                                      candidate = candidate_var))))
  }
  stop("lambda = \"replicates\" requires summarized mean data plus ",
       "candidate_sd and reference_sd columns. Summarize raw measurements ",
       "with replicate_to_mean() first.", call. = FALSE)
}

.with_seed_preserved <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed))
    stop("'seed' must be NULL or a single finite integer.", call. = FALSE)
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
          else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.fit_point_spec <- function(spec, index = seq_along(spec$x)) {
  xx <- spec$x[index]; yy <- spec$y[index]
  ww <- if (is.null(spec$weights)) NULL else spec$weights[index]
  if (spec$method == "ols") {
    fit <- if (is.null(ww)) stats::lm(yy ~ xx) else stats::lm(yy ~ xx, weights = ww)
    return(setNames(unname(stats::coef(fit)), c("intercept", "slope")))
  }
  if (spec$method == "wls") {
    if (identical(spec$weights_mode, "residual")) {
      fit <- stats::lm(yy ~ xx)
      for (iteration in seq_len(spec$weight_iterations)) {
        sd_fit <- stats::lm(abs(stats::residuals(fit)) ~ xx)
        fitted_sd <- as.numeric(stats::fitted(sd_fit))
        if (any(!is.finite(fitted_sd) | fitted_sd <= 0))
          stop("Residual SD refit failed.", call. = FALSE)
        fit <- stats::lm(yy ~ xx, weights = 1 / fitted_sd^2)
      }
    } else {
      if (identical(spec$weights_mode, "constant_cv")) ww <- 1 / xx^2
      fit <- stats::lm(yy ~ xx, weights = ww)
    }
    return(setNames(unname(stats::coef(fit)), c("intercept", "slope")))
  }
  if (spec$method == "deming") return(.deming_point(xx, yy, spec$lambda))
  if (spec$method == "wdeming") {
    if (identical(spec$weights_mode, "constant_cv"))
      return(.constant_cv_point(xx, yy, spec$lambda,
                                spec$epsilon)$coefficients)
    if (identical(spec$weights_mode, "precision_profile"))
      return(.precision_profile_point(xx, yy, spec$reference_fun,
        spec$candidate_fun, spec$epsilon)$coefficients)
    return(.deming_point(xx, yy, spec$lambda, ww))
  }
  if (spec$method == "pb") return(.pb_point(xx, yy)[c("intercept", "slope")])
  stop("Unknown regression fit specification.", call. = FALSE)
}

.resample_indices <- function(ids) {
  sample.int(length(ids), replace = TRUE)
}

.bootstrap_regression <- function(spec, replicates, seed = NULL) {
  draws <- .with_seed_preserved(seed, replicate(as.integer(replicates), {
    index <- .resample_indices(spec$id)
    tryCatch(.fit_point_spec(spec, index), error = function(e)
      c(intercept = NA_real_, slope = NA_real_))
  }))
  draws <- t(draws)
  colnames(draws) <- c("intercept", "slope")
  ok <- stats::complete.cases(draws) & apply(draws, 1L, function(z) all(is.finite(z)))
  valid <- draws[ok, , drop = FALSE]
  failed <- sum(!ok)
  if (nrow(valid) < max(100L, ceiling(replicates / 2)))
    stop("Too few successful bootstrap fits (", nrow(valid), " of ",
         replicates, ").", call. = FALSE)
  if (failed)
    warning(failed, " of ", replicates, " bootstrap fits failed and were omitted.",
            call. = FALSE)
  list(samples = valid, requested = as.integer(replicates),
       effective = nrow(valid), failed = failed, seed = seed)
}


# Internal function: Passing-Bablok Algorithm I point estimate
.pb_point <- function(x, y) {
  n <- length(x)
  pairs <- utils::combn(seq_len(n), 2L)
  dx <- x[pairs[2L, ]] - x[pairs[1L, ]]
  dy <- y[pairs[2L, ]] - y[pairs[1L, ]]
  duplicate_pair <- abs(dx) <= .Machine$double.eps &
    abs(dy) <= .Machine$double.eps
  slopes <- dy[!duplicate_pair] / dx[!duplicate_pair]

  if (!length(slopes))
    stop("No valid pairwise slopes for Passing-Bablok regression.")

  k <- sum(slopes < -1)
  slopes <- slopes[abs(slopes + 1) > .Machine$double.eps]
  slopes <- sort(slopes)
  n_slopes <- length(slopes)
  if (!n_slopes)
    stop("No valid pairwise slopes remain after excluding slopes equal to -1.")

  shifted_rank <- function(rank) ((rank + k - 1L) %% n_slopes) + 1L
  middle <- if (n_slopes %% 2L) {
    (n_slopes + 1L) / 2L
  } else {
    c(n_slopes / 2L, n_slopes / 2L + 1L)
  }
  slope <- mean(slopes[shifted_rank(middle)])
  c(intercept = stats::median(y - slope * x), slope = slope,
    k = k, n_slopes = n_slopes)
}


# Internal function: Passing-Bablok regression
.pb_fit <- function(x, y, conf.level = 0.95) {

  n <- length(x)
  if (n < 5L)
    stop("Passing-Bablok regression requires at least 5 complete pairs, got ", n, ".")

  # All pairwise slopes. Vertical pairs are represented by signed infinities;
  # only exact duplicate pairs are undefined and omitted (EP09c Appendix I).
  pairs <- utils::combn(seq_len(n), 2L)
  denom <- x[pairs[2L, ]] - x[pairs[1L, ]]
  delta_y <- y[pairs[2L, ]] - y[pairs[1L, ]]
  duplicate_pair <- abs(denom) <= .Machine$double.eps &
    abs(delta_y) <= .Machine$double.eps
  slopes <- delta_y[!duplicate_pair] / denom[!duplicate_pair]

  if (!length(slopes))
    stop("No valid pairwise slopes for Passing-Bablok regression.")

  # Algorithm I bias correction: omit slopes equal to -1, count slopes below
  # -1 as K, then shift all point and interval ranks by K.
  K <- sum(slopes < -1)
  slopes <- slopes[abs(slopes + 1) > .Machine$double.eps]
  slopes <- sort(slopes)
  m <- length(slopes)
  if (m == 0L)
    stop("No valid pairwise slopes remain after excluding slopes equal to -1.")
  shifted_rank <- function(rank) ((rank + K - 1L) %% m) + 1L

  # Regression estimates
  middle <- if (m %% 2L) (m + 1L) / 2L else c(m / 2L, m / 2L + 1L)
  slope_est <- mean(slopes[shifted_rank(middle)])
  intercept_est <- median(y - slope_est * x)

  # Jackknife parameter covariance is retained for prediction intervals.
  # Passing-Bablok confidence limits below continue to use the method's
  # rank-based construction; the jackknife quantities are used only when a
  # new-point confidence or prediction interval is requested.
  jack <- t(vapply(seq_len(n), function(i) {
    z <- tryCatch(.pb_point(x[-i], y[-i]), error = function(e) NULL)
    if (is.null(z)) c(intercept = NA_real_, slope = NA_real_)
    else z[c("intercept", "slope")]
  }, numeric(2)))
  jack_ok <- stats::complete.cases(jack)
  if (sum(jack_ok) >= 3L) {
    jack <- jack[jack_ok, , drop = FALSE]
    nj <- nrow(jack)
    jack_cov <- (nj - 1) / nj * stats::cov(jack) * (nj - 1)
    int_se <- sqrt(jack_cov["intercept", "intercept"])
    slope_se <- sqrt(jack_cov["slope", "slope"])
    cov_xy <- jack_cov["intercept", "slope"]
  } else {
    int_se <- slope_se <- cov_xy <- NA_real_
  }

  # Rank confidence interval from Algorithm I (EP09c I6-I11).
  z <- qnorm(1 - (1 - conf.level) / 2)
  c_gamma <- z * sqrt(n * (n - 1) * (2 * n + 5) / 18)
  rank_lower <- max(1L, as.integer(round((m - c_gamma) / 2)))
  rank_upper <- min(m, m - rank_lower + 1L)
  ci_lower <- slopes[shifted_rank(rank_lower)]
  ci_upper <- slopes[shifted_rank(rank_upper)]
  if (ci_lower > ci_upper) {
    tmp <- ci_lower; ci_lower <- ci_upper; ci_upper <- tmp
  }

  # Intercept confidence interval (based on slope CI bounds)
  int_lower <- median(y - ci_upper * x)
  int_upper <- median(y - ci_lower * x)

  # Residuals (for sigma)
  residuals <- y - (intercept_est + slope_est * x)
  sigma <- sqrt(sum(residuals^2) / (n - 2))

  list(
    method        = "Passing-Bablok",
    intercept     = intercept_est,
    slope         = slope_est,
    intercept_se  = int_se,
    slope_se      = slope_se,
    cov_xy        = cov_xy,
    intercept_ci  = c(int_lower, int_upper),
    slope_ci      = c(ci_lower, ci_upper),
    sigma         = sigma,
    residuals     = residuals,
    lambda        = NA_real_,
    k_adjustment  = K,
    n_slopes      = m,
    n             = n,
    conf.level    = conf.level
  )
}


# Internal function: resolve weights specification to numeric vector
.resolve_weights <- function(x, weights) {
  if (is.null(weights)) return(NULL)

  if (!is.character(weights) || length(weights) != 1L)
    stop("'weights' must be a character string or NULL.", call. = FALSE)

  # Built-in modes
  if (weights %in% c("equal", "1/y", "1/y^2", "1/x", "1/x^2",
                     "constant_cv")) {
    w <- switch(weights,
      "equal" = rep(1, x$n),
      "1/y"   = 1 / pmax(abs(x$candidate), .Machine$double.eps),
      "1/y^2" = 1 / pmax(x$candidate^2, .Machine$double.eps),
      "1/x"   = 1 / pmax(abs(x$reference), .Machine$double.eps),
      "1/x^2" = 1 / pmax(x$reference^2, .Machine$double.eps),
      "constant_cv" = {
        if (any(x$reference <= 0))
          stop("weights = \"constant_cv\" requires positive reference ",
               "concentrations.", call. = FALSE)
        1 / x$reference^2
      }
    )
    return(w)
  }

  # Column name in data
  if (weights %in% names(x$data)) {
    w <- x$data[[weights]][x$complete_idx]
    if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0))
      stop("'weights' column must contain non-negative finite numeric values.",
           call. = FALSE)
    return(w)
  }

  stop("'weights' must be one of 'equal', '1/y', '1/y^2', '1/x', '1/x^2', ",
       "'constant_cv', 'precision_profile', 'residual',\n",
       "  a column name in data, or NULL.", call. = FALSE)
}


#' Regression analysis (S3 method)
#'
#' Perform regression analysis with the reference/comparative method on the
#' X-axis and the candidate/test method on the Y-axis.
#' Supports five methods: Ordinary Least Squares (OLS), Weighted Least Squares (WLS),
#' Deming regression, Weighted Deming regression, and Passing-Bablok regression.
#'
#' @param x          mcr object
#' @param method     Regression method: "ols" (default), "wls", "deming", "wdeming", "pb"
#' @param conf.level Confidence level, default 0.95
#' @param lambda     Variance ratio for Deming regression (candidate/Y variance
#'   divided by reference/X variance), or `"replicates"` to estimate it from
#'   the summarized SD/n columns supplied to [mcr()]. Default 1.
#' @param weights    Weights specification, default NULL (unweighted).
#'   Character string for built-in modes: "equal", "1/y", "1/y^2", "1/x", "1/x^2",
#'   "constant_cv", "precision_profile", "residual", or a column name in data
#'   containing non-negative weights. Constant-CV means exactly `1/x^2` for
#'   WLS and Linnet's iterative fit for weighted Deming.
#'   The residual model follows EP09c Appendix H: fit absolute residuals on
#'   reference concentration and use the inverse squared fitted SD.
#'   For OLS, non-residual weights are passed to lm(). WLS requires a weight
#'   specification and is the only method supporting "residual". Weighted
#'   Deming requires non-residual weights. Deming does not support weights;
#'   Passing-Bablok ignores them with a warning.
#' @param variance_models For precision-profile Deming, a two-element list
#'   identifying `reference` and `candidate` variance models. Each is either a
#'   `precision` object after [profile()] or a function of concentration that
#'   returns positive finite variances.
#' @param profile_components Precision-profile components for reference and
#'   candidate methods; both default to `"total"`.
#' @param ci_method Regression-parameter interval method. `"auto"` uses
#'   analytic intervals for OLS/WLS, full-fit jackknife intervals for Deming
#'   methods, and the Appendix I rank interval for Passing-Bablok. Explicit
#'   Appendix K2 choices are `"bootstrap_percentile"` and
#'   `"bootstrap_standard"`.
#' @param bootstrap_replicates Number of paired Bootstrap resamples; at least
#'   5000 when a Bootstrap interval is requested.
#' @param seed Optional integer seed. The caller's random-number state is
#'   restored after a seeded analysis.
#' @param mdl Optional medical decision levels for which bias Bootstrap draws
#'   are retained with the regression result.
#' @param weight_iterations Number of residual-weight updates when
#'   \code{weights = "residual"}; default 4, as recommended by EP09c.
#' @param epsilon Convergence tolerance for constant-CV or precision-profile
#'   Deming. NULL uses the underlying algorithm's default.
#' @param ...        Additional arguments passed to lm (OLS/WLS only)
#'
#' @return Updated mcr object (invisible), results stored in $regression
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   regression(obj, method = "ols")
#'   regression(obj, method = "deming")
#' @export
regression.mcr <- function(x,
                            method = c("ols", "wls", "deming", "wdeming", "pb"),
                            conf.level = 0.95, lambda = 1,
                            weights = NULL, variance_models = NULL,
                            profile_components = c(reference = "total",
                                                   candidate = "total"),
                            ci_method = c("auto", "analytic", "jackknife", "rank",
                              "bootstrap_percentile", "bootstrap_standard"),
                            bootstrap_replicates = 5000L, seed = NULL, mdl = NULL,
                            weight_iterations = 4L,
                            epsilon = NULL, ...) {
  if (!inherits(x, "mcr")) stop("Input must be an 'mcr' object.")
  method <- match.arg(method); ci_method <- match.arg(ci_method)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(weight_iterations) || length(weight_iterations) != 1L ||
      !is.finite(weight_iterations) || weight_iterations < 1L ||
      weight_iterations != as.integer(weight_iterations))
    stop("'weight_iterations' must be a positive integer.", call. = FALSE)
  if (!is.null(mdl) && (!is.numeric(mdl) || any(!is.finite(mdl))))
    stop("'mdl' must be NULL or a finite numeric vector.", call. = FALSE)

  lambda_info <- NULL
  replicate_mode <- identical(lambda, "replicates")
  if (is.character(lambda) && !replicate_mode)
    stop("'lambda' must be a positive number or \"replicates\".", call. = FALSE)
  if (replicate_mode) {
    if (!method %in% c("deming", "wdeming"))
      stop("lambda = \"replicates\" is available only for Deming methods.",
           call. = FALSE)
    rep_data <- .replicate_lambda_data(x)
    all_ref <- rep_data$x; all_cand <- rep_data$y; all_ids <- rep_data$id
    all_complete_idx <- rep_data$complete_idx
    lambda <- rep_data$lambda; lambda_info <- rep_data$metadata
    lambda_info$estimate <- lambda
  } else {
    all_ref <- x$reference; all_cand <- x$candidate; all_ids <- x$id
    all_complete_idx <- x$complete_idx
  }
  if (method %in% c("deming", "wdeming") &&
      (!is.numeric(lambda) || length(lambda) != 1L ||
       !is.finite(lambda) || lambda <= 0))
    stop("'lambda' must be a single positive number.", call. = FALSE)

  minimum <- if (method == "pb") 5L else 3L
  if (length(all_ref) < minimum)
    stop("Need at least ", minimum, " complete pairs for regression, got ",
         length(all_ref), ".", call. = FALSE)
  ref_vals <- all_ref; cand_vals <- all_cand
  ids <- all_ids; complete_idx <- all_complete_idx
  n_selected <- length(ref_vals)

  residual_weighting <- identical(weights, "residual")
  constant_cv <- identical(weights, "constant_cv")
  precision_profile <- identical(weights, "precision_profile")
  if (residual_weighting && method != "wls")
    stop("weights = \"residual\" is supported only for method = \"wls\".",
         call. = FALSE)
  if (precision_profile && method != "wdeming")
    stop("weights = \"precision_profile\" requires method = \"wdeming\".",
         call. = FALSE)
  if (constant_cv && !method %in% c("wls", "wdeming"))
    stop("weights = \"constant_cv\" requires WLS or weighted Deming.",
         call. = FALSE)
  if (replicate_mode && !is.null(weights) &&
      !weights %in% c("constant_cv", "precision_profile"))
    stop("Column and fixed observation weights cannot be combined with ",
         "lambda = \"replicates\" after replicate aggregation.", call. = FALSE)

  w <- NULL
  if (!residual_weighting && !precision_profile && !constant_cv &&
      !is.null(weights)) {
    resolved <- .resolve_weights(x, weights)
    if (replicate_mode)
      stop("Observation weights cannot be mapped after replicate aggregation.",
           call. = FALSE)
    w <- resolved
  }
  if (constant_cv && method == "wls") {
    if (any(ref_vals <= 0))
      stop("Constant-CV WLS requires strictly positive reference concentrations.",
           call. = FALSE)
    w <- 1 / ref_vals^2
  }

  reference_profile <- candidate_profile <- NULL
  if (precision_profile) {
    if (!is.list(variance_models) || length(variance_models) != 2L)
      stop("'variance_models' must be a list with reference and candidate models.",
           call. = FALSE)
    if (!is.null(names(variance_models)) &&
        all(c("reference", "candidate") %in% names(variance_models))) {
      reference_model <- variance_models$reference
      candidate_model <- variance_models$candidate
    } else {
      reference_model <- variance_models[[1L]]
      candidate_model <- variance_models[[2L]]
    }
    if (length(profile_components) == 1L)
      profile_components <- rep(profile_components, 2L)
    if (is.null(names(profile_components)))
      names(profile_components) <- c("reference", "candidate")
    if (!all(c("reference", "candidate") %in% names(profile_components)))
      stop("'profile_components' must identify reference and candidate components.",
           call. = FALSE)
    reference_profile <- .profile_variance_function(
      reference_model, profile_components[["reference"]], "Reference")
    candidate_profile <- .profile_variance_function(
      candidate_model, profile_components[["candidate"]], "Candidate")
    if (!is.null(reference_profile$range) &&
        any(ref_vals < reference_profile$range[1L] |
            ref_vals > reference_profile$range[2L]))
      warning("Reference concentrations extend beyond the fitted precision-profile range.",
              call. = FALSE)
    if (!is.null(candidate_profile$range) &&
        any(cand_vals < candidate_profile$range[1L] |
            cand_vals > candidate_profile$range[2L]))
      warning("Candidate concentrations extend beyond the fitted precision-profile range.",
              call. = FALSE)
  }

  fit_spec <- list(
    method = method, x = ref_vals, y = cand_vals, id = ids, weights = w,
    weights_mode = weights, lambda = if (method %in% c("deming", "wdeming"))
      lambda else NA_real_, weight_iterations = as.integer(weight_iterations),
    reference_fun = reference_profile$fun %||% NULL,
    candidate_fun = candidate_profile$fun %||% NULL,
    epsilon = epsilon %||% if (constant_cv) 1e-10 else 1e-8
  )

  residual_history <- residual_sd_model <- diagnostic <- NULL
  if (method %in% c("ols", "wls")) {
    if (method == "wls" && is.null(w) && !residual_weighting)
      stop("'weights' must be provided for WLS.", call. = FALSE)
    if (residual_weighting) {
      fit <- stats::lm(cand_vals ~ ref_vals, ...)
      residual_history <- data.frame(iteration = 0L,
        intercept = unname(coef(fit)[1L]), slope = unname(coef(fit)[2L]))
      for (iteration in seq_len(as.integer(weight_iterations))) {
        residual_sd_model <- stats::lm(abs(stats::residuals(fit)) ~ ref_vals)
        fitted_sd <- as.numeric(stats::fitted(residual_sd_model))
        if (any(!is.finite(fitted_sd) | fitted_sd <= 0))
          stop("Residual SD model produced non-positive or non-finite values.",
               call. = FALSE)
        w <- 1 / fitted_sd^2; fit_spec$weights <- w
        fit <- stats::lm(cand_vals ~ ref_vals, weights = w, ...)
        residual_history <- rbind(residual_history, data.frame(
          iteration = iteration, intercept = unname(coef(fit)[1L]),
          slope = unname(coef(fit)[2L])))
      }
    } else {
      fit <- if (is.null(w)) stats::lm(cand_vals ~ ref_vals, ...)
        else stats::lm(cand_vals ~ ref_vals, weights = w, ...)
    }
    s <- summary(fit); coefs <- coef(s); covariance <- stats::vcov(fit)
    t_val <- stats::qt(1 - (1 - conf.level) / 2, stats::df.residual(fit))
    result <- list(method = toupper(method), intercept = coefs[1L, 1L],
      slope = coefs[2L, 1L], intercept_se = coefs[1L, 2L],
      slope_se = coefs[2L, 2L], cov_xy = covariance[1L, 2L],
      intercept_ci = coefs[1L, 1L] + c(-1, 1) * t_val * coefs[1L, 2L],
      slope_ci = coefs[2L, 1L] + c(-1, 1) * t_val * coefs[2L, 2L],
      sigma = s$sigma, r_squared = s$r.squared,
      residuals = stats::residuals(fit), lambda = NA_real_)
  } else if (method == "deming") {
    if (!is.null(weights))
      stop("Weights are not supported for unweighted Deming regression.",
           call. = FALSE)
    result <- .deming_fit(ref_vals, cand_vals, lambda, conf.level)
  } else if (method == "wdeming") {
    if (constant_cv) {
      diagnostic <- .constant_cv_point(ref_vals, cand_vals, lambda,
                                        fit_spec$epsilon)
      coefficients <- diagnostic$coefficients
      residuals <- cand_vals - coefficients[1L] - coefficients[2L] * ref_vals
      result <- list(method = "Constant-CV Deming", intercept = coefficients[1L],
        slope = coefficients[2L], sigma = sqrt(sum(diagnostic$weights * residuals^2) /
          sum(diagnostic$weights) * n_selected / (n_selected - 2L)),
        residuals = residuals, lambda = lambda,
        latent_concentration = diagnostic$latent_concentration,
        latent_reference = diagnostic$latent_reference,
        latent_candidate = diagnostic$latent_candidate,
        weights = diagnostic$weights, convergence = diagnostic$convergence)
    } else if (precision_profile) {
      diagnostic <- .precision_profile_point(ref_vals, cand_vals,
        reference_profile$fun, candidate_profile$fun, fit_spec$epsilon)
      coefficients <- diagnostic$coefficients
      residuals <- cand_vals - coefficients[1L] - coefficients[2L] * ref_vals
      result <- list(method = "Precision-profile Deming",
        intercept = coefficients[1L], slope = coefficients[2L],
        sigma = sqrt(sum(residuals^2) / (n_selected - 2L)),
        residuals = residuals, lambda = lambda,
        latent_concentration = diagnostic$latent_concentration,
        fitted_candidate = diagnostic$fitted_candidate,
        reference_variance = diagnostic$reference_variance,
        candidate_variance = diagnostic$candidate_variance,
        final_weights = list(reference = 1 / diagnostic$reference_variance,
                             candidate = 1 / diagnostic$candidate_variance),
        scaled_residuals = diagnostic$scaled_residuals,
        convergence = diagnostic$convergence,
        profile_sources = c(reference = reference_profile$source,
                            candidate = candidate_profile$source),
        profile_ranges = list(reference = reference_profile$range,
                              candidate = candidate_profile$range))
    } else {
      if (is.null(w)) stop("'weights' must be provided for Weighted Deming.",
                           call. = FALSE)
      result <- .wdeming_fit(ref_vals, cand_vals, lambda, w, conf.level)
    }
  } else {
    if (!is.null(weights))
      warning("Weights are ignored for Passing-Bablok regression.", call. = FALSE)
    result <- .pb_fit(ref_vals, cand_vals, conf.level)
  }

  resolved_ci <- if (ci_method == "auto") {
    if (method %in% c("ols", "wls")) "analytic"
    else if (method == "pb") "rank" else "jackknife"
  } else ci_method
  if (resolved_ci == "analytic" && !method %in% c("ols", "wls"))
    stop("ci_method = \"analytic\" is available only for OLS/WLS.",
         call. = FALSE)
  if (resolved_ci == "rank" && method != "pb")
    stop("ci_method = \"rank\" is available only for Passing-Bablok.",
         call. = FALSE)

  parameter_samples <- NULL; bootstrap <- NULL
  if (resolved_ci == "jackknife") {
    parameter_samples <- t(vapply(seq_len(n_selected), function(i)
      tryCatch(.fit_point_spec(fit_spec, setdiff(seq_len(n_selected), i)),
               error = function(e) c(intercept = NA_real_, slope = NA_real_)),
      numeric(2)))
    ok <- stats::complete.cases(parameter_samples)
    if (sum(ok) < 3L) stop("Too few successful jackknife fits.", call. = FALSE)
    parameter_samples <- parameter_samples[ok, , drop = FALSE]
    nj <- nrow(parameter_samples)
    center <- colMeans(parameter_samples)
    covariance <- (nj - 1) / nj * crossprod(sweep(parameter_samples, 2L, center))
    result$intercept_se <- sqrt(covariance[1L, 1L])
    result$slope_se <- sqrt(covariance[2L, 2L])
    result$cov_xy <- covariance[1L, 2L]
    t_val <- stats::qt(1 - (1 - conf.level) / 2, n_selected - 2L)
    result$intercept_ci <- result$intercept + c(-1, 1) * t_val * result$intercept_se
    result$slope_ci <- result$slope + c(-1, 1) * t_val * result$slope_se
  } else if (grepl("^bootstrap_", resolved_ci)) {
    if (!is.numeric(bootstrap_replicates) || length(bootstrap_replicates) != 1L ||
        !is.finite(bootstrap_replicates) || bootstrap_replicates < 5000L ||
        bootstrap_replicates != as.integer(bootstrap_replicates))
      stop("Bootstrap regression intervals require at least 5000 replicates.",
           call. = FALSE)
    bootstrap <- .bootstrap_regression(fit_spec, bootstrap_replicates, seed)
    bootstrap$method <- resolved_ci
    parameter_samples <- bootstrap$samples
    covariance <- stats::cov(parameter_samples)
    result$intercept_se <- sqrt(covariance[1L, 1L])
    result$slope_se <- sqrt(covariance[2L, 2L])
    result$cov_xy <- covariance[1L, 2L]
    alpha <- (1 - conf.level) / 2
    if (resolved_ci == "bootstrap_percentile") {
      result$intercept_ci <- as.numeric(stats::quantile(parameter_samples[, 1L],
        c(alpha, 1 - alpha), names = FALSE))
      result$slope_ci <- as.numeric(stats::quantile(parameter_samples[, 2L],
        c(alpha, 1 - alpha), names = FALSE))
    } else {
      t_val <- stats::qt(1 - alpha, n_selected - 2L)
      result$intercept_ci <- result$intercept + c(-1, 1) * t_val * result$intercept_se
      result$slope_ci <- result$slope + c(-1, 1) * t_val * result$slope_se
    }
    if (length(mdl))
      bootstrap$bias_samples <- outer(parameter_samples[, 2L], mdl, `*`) +
        parameter_samples[, 1L] -
        matrix(mdl, nrow(parameter_samples), length(mdl), byrow = TRUE)
  }

  result$weights_spec <- weights
  result$weights <- result[["weights", exact = TRUE]] %||% w
  result$residual_sd_coefficients <- if (residual_weighting)
    unname(coef(residual_sd_model)) else NULL
  result$residual_weight_history <- residual_history
  result$weight_iterations <- if (residual_weighting)
    as.integer(weight_iterations) else NULL
  result$n <- n_selected; result$conf.level <- conf.level
  result$ci_method <- resolved_ci; result$ci_method_requested <- ci_method
  result$parameter_samples <- parameter_samples
  result$bootstrap <- bootstrap
  result$bootstrap_replicates <- if (is.null(bootstrap)) NULL else bootstrap$requested
  result$seed <- if (is.null(bootstrap)) NULL else seed
  result$lambda_metadata <- lambda_info
  result$replicate_aggregation <- replicate_mode
  result$complete_idx <- complete_idx
  result$x_vals <- ref_vals; result$y_vals <- cand_vals
  result$fit_spec <- fit_spec

  x$regression <- structure(result, class = "mcr_regression")
  print(x$regression)
  invisible(x)
}


#' Outlier detection (S3 method)
#'
#' Perform outlier detection based on the difference or percent difference between two measurements.
#' Reuses four methods from ../precision/outliers.R.
#'
#' @param x      mcr object
#' @param method Outlier detection method: "grubbs" (default), "esd", "dixon", "iqr"
#' @param type   Data type to compute: "difference" (default) or "percent" (percent difference)
#' @param percent_denominator Denominator for percent differences: reference
#'   (default), pair mean, or candidate.
#' @param alpha  Significance level, default 0.05 (Grubbs/ESD/Dixon only)
#' @param coef   IQR multiplier, default 1.5 (IQR only)
#' @param ...    Additional arguments passed to esd_test / dixon_test
#'
#' @return Updated mcr object (invisible), results stored in $outlier
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   outlier(obj, method = "grubbs")
#'   outlier(obj, method = "iqr", type = "percent")
#' @export
outlier.mcr <- function(x, method = c("grubbs", "esd", "dixon", "iqr"),
                         type = c("difference", "percent"),
                         percent_denominator = c("reference", "mean", "candidate"),
                         alpha = 0.05, coef = 1.5, ...) {

  if (!inherits(x, "mcr"))
    stop("Input must be an 'mcr' object.")

  method <- match.arg(method)
  type   <- match.arg(type)
  percent_denominator <- match.arg(percent_denominator)

  if (x$n < 3L)
    stop("Need at least 3 complete pairs for outlier detection, got ", x$n, ".")

  # Compute difference/percent difference
  if (type == "difference") {
    d     <- x$candidate - x$reference
    label <- sprintf("Difference (%s - %s)", x$candidate_name, x$reference_name)
  } else {
    denominator <- switch(percent_denominator,
      reference = x$reference,
      mean = (x$candidate + x$reference) / 2,
      candidate = x$candidate)
    d <- (x$candidate - x$reference) / .safe_div(denominator) * 100
    denominator_label <- switch(percent_denominator,
      reference = x$reference_name, mean = "pair mean",
      candidate = x$candidate_name)
    label <- sprintf("Percent diff ((%s - %s) / %s * 100)",
                     x$candidate_name, x$reference_name, denominator_label)
  }

  # Call outliers_test (requires data.frame + column name)
  tmp_data <- data.frame(.mcr_diff = d, check.names = FALSE)
  raw <- switch(method,
    grubbs = outliers_test(tmp_data, ".mcr_diff", method = "grubbs", alpha = alpha),
    esd    = outliers_test(tmp_data, ".mcr_diff", method = "esd", alpha = alpha, ...),
    dixon  = outliers_test(tmp_data, ".mcr_diff", method = "dixon", alpha = alpha, ...),
    iqr    = outliers_test(tmp_data, ".mcr_diff", method = "iqr", coef = coef)
  )

  # Map back to original data row numbers
  n_out <- raw$n_outliers
  idx <- raw$indices
  if (n_out > 0 && length(idx) > 0) {
    orig_rows <- x$complete_idx[idx]
  } else {
    orig_rows <- integer(0)
    n_out <- 0L
  }

  # Storage: adapt to outliers_test return structure
  statistic <- if (method == "esd") raw$details$R %||% NA_real_
               else raw$details$statistic %||% NA_real_
  critical  <- raw$details$critical %||% NA_real_

  result <- structure(list(
    method     = method,
    data_type  = type,
    percent_denominator = if (type == "percent") percent_denominator else NULL,
    data_label = label,
    n_out      = n_out,
    indices    = orig_rows,
    values     = raw$values %||% numeric(0),
    statistic  = statistic,
    critical   = critical,
    alpha      = alpha,
    coef       = coef
  ), class = "mcr_outlier")
  x$outlier <- result

  print(result)
  invisible(x)
}


#' Bland-Altman analysis (S3 method)
#'
#' Calculate Limits of Agreement between two measurement methods,
#' supporting three y-axis metrics: difference, ratio, and percent difference.
#'
#' @param x           mcr object
#' @param type        Y-axis metric: "difference" (default),
#'                    "ratio" or "percent"
#' @param x_axis      X-axis: "mean" (average of the two methods, default), "candidate" (test method values),
#'                    "reference" (reference method values)
#' @param percent_denominator Denominator used when \code{type = "percent"}:
#'                    "reference" (default), "mean" (average of candidate and
#'                    reference), or "candidate". This is independent of
#'                    \code{x_axis}.
#' @param agree.level Agreement level for limits of agreement, default 0.95 (i.e., 95% LoA)
#' @param conf.level  Confidence level for LoA confidence intervals, default 0.95
#' @param method      Estimation method: "parametric" uses the mean, its t
#'                    confidence interval, the standard deviation, and parametric
#'                    limits of agreement; "nonparametric" uses the median, an
#'                    exact sign-test order-statistic interval, empirical limits,
#'                    and percentile Bootstrap confidence intervals for the limits.
#'                    "hodges_lehmann" uses the Walsh-average pseudomedian and
#'                    an uncorrected Wilcoxon/Tukey interval; its LoA remain
#'                    empirical quantiles with percentile Bootstrap intervals.
#' @param median_ci_method For the ordinary nonparametric median, either the
#'                    conservative order-statistic interval (default) or the
#'                    EP09c continuous-rank normal approximation with linear
#'                    interpolation between adjacent order statistics.
#' @param bootstrap_replicates Number of Bootstrap resamples for the nonparametric method,
#'                    default 5000. Must be at least 100.
#' @param seed         Optional integer random seed for reproducible Bootstrap intervals.
#' @param ...         Additional arguments
#'
#' @return Updated mcr object (invisible), with results stored in $bland_altman.
#' The result records the selected method, its centre estimate and
#' \code{center_ci}. Method-specific aliases are \code{mean_diff_ci} for the
#' parametric mean and \code{median_diff_ci} for the nonparametric median.
#' The original \code{mean_diff} and \code{sd_diff} fields are retained for
#' compatibility.
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   bland_altman(obj)
#'   bland_altman(obj, type = "ratio")
#'   bland_altman(obj, type = "percent", x_axis = "reference")
#'   bland_altman(obj, type = "percent", x_axis = "mean",
#'                 percent_denominator = "mean")
#' @export
bland_altman.mcr <- function(x, type = c("difference", "ratio", "percent"),
                              x_axis = c("mean", "candidate", "reference"),
                              percent_denominator = c("reference", "mean", "candidate"),
                              agree.level = 0.95, conf.level = 0.95,
                              method = c("parametric", "nonparametric",
                                         "hodges_lehmann"),
                              median_ci_method = c("order_statistic", "interpolated"),
                              bootstrap_replicates = 5000L, seed = NULL, ...) {

  if (!inherits(x, "mcr"))
    stop("Input must be an 'mcr' object.")

  type   <- match.arg(type)
  x_axis <- match.arg(x_axis)
  percent_denominator <- match.arg(percent_denominator)
  method <- match.arg(method)
  median_ci_method <- match.arg(median_ci_method)

  candidate <- x$candidate
  reference <- x$reference

  # Parameter range validation
  if (!is.numeric(agree.level) || length(agree.level) != 1L ||
      is.na(agree.level) || agree.level <= 0 || agree.level >= 1)
    stop("'agree.level' must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(bootstrap_replicates) || length(bootstrap_replicates) != 1L ||
      is.na(bootstrap_replicates) || bootstrap_replicates < 100 || bootstrap_replicates != as.integer(bootstrap_replicates))
    stop("'bootstrap_replicates' must be an integer of at least 100.", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed != as.integer(seed)))
    stop("'seed' must be NULL or a single finite integer.", call. = FALSE)

  # Compute x-axis
  means <- (candidate + reference) / 2
  x_vals <- switch(x_axis,
    mean      = means,
    candidate = candidate,
    reference = reference
  )

  # Compute y-axis
  percent_base <- switch(percent_denominator,
    reference = reference,
    mean      = means,
    candidate = candidate
  )
  y_vals <- switch(type,
    difference = candidate - reference,
    ratio      = candidate / .safe_div(reference),
    percent    = (candidate - reference) / .safe_div(percent_base) * 100
  )

  # Limits of agreement
  n <- length(y_vals)
  mean_diff <- mean(y_vals)
  sd_diff   <- sd(y_vals)
  center_ranks <- NULL
  if (method == "parametric") {
    center <- mean_diff; z_agree <- qnorm(1 - (1 - agree.level) / 2)
    t_val <- qt(1 - (1 - conf.level) / 2, n - 1)
    center_ci <- center + c(-1, 1) * t_val * sd_diff / sqrt(n)
    center_ci_coverage <- conf.level
    loa_lower <- center - z_agree * sd_diff; loa_upper <- center + z_agree * sd_diff
    se_loa <- sqrt(sd_diff^2 * (1 / n + z_agree^2 / (2 * (n - 1))))
    loa_lower_ci <- c(loa_lower - t_val * se_loa, loa_lower + t_val * se_loa)
    loa_upper_ci <- c(loa_upper - t_val * se_loa, loa_upper + t_val * se_loa)
  } else if (method == "nonparametric") {
    center <- median(y_vals)
    sorted_y <- sort(y_vals)
    if (median_ci_method == "order_statistic") {
      # Narrowest symmetric sign-test interval meeting the requested coverage.
      candidate_ranks <- seq_len(floor((n + 1L) / 2L))
      coverages <- 1 - 2 * stats::pbinom(candidate_ranks - 1L, n, 0.5)
      eligible <- candidate_ranks[coverages >= conf.level]
      if (length(eligible)) {
        lower_rank <- max(eligible); upper_rank <- n - lower_rank + 1L
        center_ci <- c(sorted_y[lower_rank], sorted_y[upper_rank])
        center_ci_coverage <- coverages[lower_rank]
        center_ranks <- c(lower = lower_rank, upper = upper_rank)
      } else {
        center_ci <- range(y_vals); center_ci_coverage <- 1 - 2^(1 - n)
        center_ranks <- c(lower = 1, upper = n)
        warning("The requested nonparametric median confidence level cannot be ",
                "attained with n = ", n, "; using the full observed range ",
                "(actual coverage ", signif(center_ci_coverage, 4), ").",
                call. = FALSE)
      }
    } else {
      z_ci <- stats::qnorm(1 - (1 - conf.level) / 2)
      continuous <- c(lower = (n + 1 - z_ci * sqrt(n)) / 2,
                      upper = (n + 1 + z_ci * sqrt(n)) / 2)
      interpolate_order <- function(rank) {
        rank <- min(max(rank, 1), n)
        lo <- floor(rank); hi <- ceiling(rank)
        if (lo == hi) return(sorted_y[lo])
        sorted_y[lo] + (rank - lo) * (sorted_y[hi] - sorted_y[lo])
      }
      center_ci <- vapply(continuous, interpolate_order, numeric(1))
      center_ci_coverage <- conf.level
      center_ranks <- continuous
    }
    probs <- c((1 - agree.level) / 2, 1 - (1 - agree.level) / 2)
    loa_lower <- as.numeric(quantile(y_vals, probs[1L], names = FALSE, type = 7, na.rm = TRUE))
    loa_upper <- as.numeric(quantile(y_vals, probs[2L], names = FALSE, type = 7, na.rm = TRUE))
  } else {
    hl_test <- suppressWarnings(stats::wilcox.test(
      y_vals, conf.int = TRUE, conf.level = conf.level,
      exact = FALSE, correct = FALSE
    ))
    # wilcox.test computes the Walsh-average pseudomedian and the Tukey
    # interval under the same no-continuity-correction normal approximation.
    center <- unname(hl_test$estimate)
    center_ci <- as.numeric(hl_test$conf.int)
    center_ci_coverage <- conf.level
    probs <- c((1 - agree.level) / 2, 1 - (1 - agree.level) / 2)
    loa_lower <- as.numeric(quantile(y_vals, probs[1L], names = FALSE, type = 7, na.rm = TRUE))
    loa_upper <- as.numeric(quantile(y_vals, probs[2L], names = FALSE, type = 7, na.rm = TRUE))
  }

  if (method != "parametric") {
    oldseed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
    if (!is.null(seed)) set.seed(as.integer(seed))
    boots <- replicate(as.integer(bootstrap_replicates), {
      yy <- sample(y_vals, n, replace = TRUE)
      quantile(yy, probs, names = FALSE, type = 7, na.rm = TRUE)
    })
    if (!is.null(seed)) { if (is.null(oldseed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", oldseed, envir = .GlobalEnv) }
    a <- (1 - conf.level) / 2
    loa_lower_ci <- as.numeric(quantile(boots[1L, ], c(a, 1 - a), names = FALSE, na.rm = TRUE))
    loa_upper_ci <- as.numeric(quantile(boots[2L, ], c(a, 1 - a), names = FALSE, na.rm = TRUE))
  }

  # Labels
  y_label <- switch(type,
    difference = sprintf("%s - %s", x$candidate_name, x$reference_name),
    ratio      = sprintf("%s / %s", x$candidate_name, x$reference_name),
    percent    = switch(percent_denominator,
      reference = sprintf("(%s - %s) / %s * 100",
                          x$candidate_name, x$reference_name,
                          x$reference_name),
      mean = sprintf("(%s - %s) / mean(%s, %s) * 100",
                     x$candidate_name, x$reference_name,
                     x$candidate_name, x$reference_name),
      candidate = sprintf("(%s - %s) / %s * 100",
                          x$candidate_name, x$reference_name,
                          x$candidate_name)
    )
  )
  x_label <- switch(x_axis,
    mean      = sprintf("Mean of %s and %s", x$candidate_name, x$reference_name),
    candidate = x$candidate_name,
    reference = x$reference_name
  )

  # Storage
  result <- list(
    type        = type,
    x_axis      = x_axis,
    percent_denominator = if (type == "percent")
      percent_denominator else NULL,
    n           = n,
    mean_diff   = mean_diff,
    sd_diff     = sd_diff,
    median_diff = median(y_vals),
    center      = center,
    center_estimand = switch(method, parametric = "mean",
      nonparametric = "median", hodges_lehmann = "Hodges-Lehmann pseudomedian"),
    center_ci   = setNames(center_ci, c("lower", "upper")),
    center_ci_coverage = center_ci_coverage,
    center_ci_ranks = center_ranks,
    mean_diff_ci = if (method == "parametric")
      setNames(center_ci, c("lower", "upper")) else NULL,
    median_diff_ci = if (method == "nonparametric")
      setNames(center_ci, c("lower", "upper")) else NULL,
    median_ci_method = if (method == "nonparametric") median_ci_method else NULL,
    hl_estimate = if (method == "hodges_lehmann") center else NULL,
    hl_ci = if (method == "hodges_lehmann")
      setNames(center_ci, c("lower", "upper")) else NULL,
    method      = method,
    bootstrap_replicates = if (method != "parametric") as.integer(bootstrap_replicates) else NULL,
    loa         = c(lower = loa_lower, upper = loa_upper),
    loa_ci      = list(
      lower = loa_lower_ci,
      upper = loa_upper_ci
    ),
    agree.level = agree.level,
    conf.level  = conf.level,
    complete_idx = x$complete_idx,
    y_label     = y_label,
    x_label     = x_label,
    x_vals      = x_vals,   # for plot
    y_vals      = y_vals    # for plot
  )
  x$bland_altman <- structure(result, class = "mcr_bland_altman")
  print(x$bland_altman)
  invisible(x)
}


# Resolve joint regression-parameter uncertainty for predict() and bias().
.regression_interval_samples <- function(rg, ci_method, bootstrap_replicates,
                                         seed) {
  resolved <- if (ci_method == "auto") {
    if (identical(rg$method, "Passing-Bablok")) "bootstrap_percentile"
    else rg$ci_method %||% "analytic"
  } else ci_method
  if (!grepl("^bootstrap_", resolved))
    return(list(method = resolved, samples = NULL, bootstrap = NULL))
  if (!is.numeric(bootstrap_replicates) || length(bootstrap_replicates) != 1L ||
      !is.finite(bootstrap_replicates) || bootstrap_replicates < 5000L ||
      bootstrap_replicates != as.integer(bootstrap_replicates))
    stop("Bootstrap intervals require at least 5000 replicates.", call. = FALSE)
  cached <- rg$bootstrap
  reusable <- !is.null(cached) && identical(cached$method, resolved) &&
    cached$requested >= bootstrap_replicates &&
    (is.null(seed) || identical(cached$seed, seed))
  boot <- if (reusable) cached else
    .bootstrap_regression(rg$fit_spec, bootstrap_replicates, seed)
  boot$method <- resolved
  list(method = resolved, samples = boot$samples, bootstrap = boot)
}


#' Predict (S3 method) -- predict based on regression results
#'
#' Predict values based on stored regression results.
#'
#' @param object     mcr object (must call regression() first)
#' @param newdata    Data frame containing input values. If provided,
#'                   reference or candidate names a column in it.
#' @param reference  Reference-method input for forward prediction. With
#'                   newdata, a column name; otherwise a numeric vector.
#'                   If NULL, uses the original reference values.
#' @param candidate  Candidate-method input for inverse prediction. With
#'                   newdata, a column name; otherwise a numeric vector.
#'                   If NULL, uses the original candidate values.
#' @param inverse    If FALSE, predict candidate from reference. If TRUE,
#'                   algebraically predict reference from candidate.
#' @param interval   Interval type: "none" (default), "confidence", "prediction", "both".
#'                   Only forward prediction (inverse=FALSE) supports intervals.
#' @param conf.level Confidence level, default 0.95
#' @param ci_method Parameter-uncertainty method. `"auto"` reuses the stored
#'   regression method except that Passing-Bablok uses percentile Bootstrap.
#' @param bootstrap_replicates Number of paired Bootstrap resamples; at least
#'   5000 when a Bootstrap method is used.
#' @param seed Optional integer seed.
#' @param ...        Additional arguments
#'
#' @return data.frame containing predicted values and (if requested) interval bounds
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   obj <- regression(obj, method = "ols")
#'   predict(obj, reference = c(40, 50, 60))
#'   predict(obj, reference = c(40, 50, 60), interval = "both")
#'   predict(obj, newdata = data.frame(xx = c(40, 50)), reference = "xx")
#'   predict(obj, inverse = TRUE)  # Inverse prediction on original data
#' @export
predict.mcr <- function(object, newdata = NULL,
                         reference = NULL, candidate = NULL,
                         inverse = FALSE,
                         interval = c("none", "confidence", "prediction", "both"),
                         conf.level = 0.95,
                         ci_method = c("auto", "analytic", "jackknife", "rank",
                           "bootstrap_percentile", "bootstrap_standard"),
                         bootstrap_replicates = 5000L, seed = NULL, ...) {

  if (!inherits(object, "mcr"))
    stop("Input must be an 'mcr' object.")
  if (is.null(object$regression))
    stop("Run regression() first before predict().")

  rg <- object$regression
  interval <- match.arg(interval)
  ci_method <- match.arg(ci_method)

  # Determine input values
  newdata_provided <- !is.null(newdata)

  if (newdata_provided) {
    # newdata mode: candidate/reference are column names
    if (!is.data.frame(newdata))
      stop("'newdata' must be a data.frame.")
  }

  # Inverse prediction (candidate to reference)
  if (inverse) {
    if (abs(rg$slope) < .Machine$double.eps)
      stop("Slope is near zero; cannot perform inverse prediction.")

    if (newdata_provided) {
      cand_col <- if (!is.null(candidate)) candidate else object$candidate_name
      if (!cand_col %in% names(newdata))
        stop("Column '", cand_col, "' not found in newdata.")
      cand_vals <- newdata[[cand_col]]
    } else {
      cand_vals <- if (!is.null(candidate)) candidate else object$candidate
    }

    ref_pred <- (cand_vals - rg$intercept) / rg$slope
    result <- data.frame(candidate = cand_vals, reference_pred = ref_pred)

    return(result)
  }

  # Forward prediction (reference to candidate)
  if (newdata_provided) {
    ref_col <- if (!is.null(reference)) reference else object$reference_name
    if (!ref_col %in% names(newdata))
      stop("Column '", ref_col, "' not found in newdata.")
    ref_vals <- newdata[[ref_col]]
  } else {
    ref_vals <- if (!is.null(reference)) reference else object$reference
  }

  if (!is.numeric(ref_vals))
    stop("Input values must be numeric.")

  pred <- rg$intercept + rg$slope * ref_vals
  result <- data.frame(reference = ref_vals, candidate_pred = pred)

  if (interval != "none") {
    n <- rg$n
    uncertainty <- .regression_interval_samples(
      rg, ci_method, bootstrap_replicates, seed)
    t_val <- qt(1 - (1 - conf.level) / 2, n - 2)

    if (!is.null(uncertainty$samples)) {
      draws <- outer(uncertainty$samples[, 2L], ref_vals, `*`) +
        uncertainty$samples[, 1L]
      se_fit <- apply(draws, 2L, stats::sd)
    } else {
      if (is.null(rg$slope_se) || !is.finite(rg$slope_se))
        stop("The regression object does not contain usable parameter uncertainty.",
             call. = FALSE)
      se_fit <- sqrt(pmax(0, rg$intercept_se^2 +
                     ref_vals^2 * rg$slope_se^2 +
                     2 * ref_vals * rg$cov_xy))
    }

    if (interval %in% c("confidence", "both")) {
      if (identical(uncertainty$method, "bootstrap_percentile")) {
        alpha <- (1 - conf.level) / 2
        limits <- vapply(seq_len(ncol(draws)), function(j)
          stats::quantile(draws[, j], c(alpha, 1 - alpha), names = FALSE),
          numeric(2))
        result$ci_lower <- limits[1L, ]; result$ci_upper <- limits[2L, ]
      } else {
        result$ci_lower <- pred - t_val * se_fit
        result$ci_upper <- pred + t_val * se_fit
      }
    }
    if (interval %in% c("prediction", "both")) {
      se_pred <- sqrt(se_fit^2 + rg$sigma^2)
      result$pi_lower <- pred - t_val * se_pred
      result$pi_upper <- pred + t_val * se_pred
    }
    attr(result, "ci_method") <- uncertainty$method
    attr(result, "bootstrap") <- uncertainty$bootstrap
  }

  result
}


#' Bias analysis (S3 method) -- estimate bias at given medical decision levels
#'
#' Compute bias at specified medical decision levels (MDL) based on stored regression results:
#'   bias = predicted_candidate - reference
#'        = (intercept + slope * x0) - x0
#'        = intercept + (slope - 1) * x0
#'
#' Bias follows the conventional candidate-minus-reference direction; a
#' positive value indicates that the candidate method is higher.
#'
#' @param x        mcr object (must call regression() first)
#' @param mdl      Medical decision levels (numeric vector) on the reference/X scale
#' @param level    Confidence level, default 0.95
#' @param interval Interval type: "" (point estimate, default), "confidence" (confidence interval),
#'                 "prediction" (prediction interval), "both" (confidence + prediction)
#' @param ci_method Parameter-uncertainty method. `"auto"` uses percentile
#'   Bootstrap for Passing-Bablok and otherwise reuses the regression method.
#' @param bootstrap_replicates Number of paired Bootstrap resamples; at least
#'   5000 when a Bootstrap method is used.
#' @param seed Optional integer seed.
#' @param ...      Additional arguments (currently unused).
#'
#' @return Updated mcr object (invisible), results stored in $bias
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   obj <- regression(obj, method = "ols")
#'   bias(obj, mdl = c(200, 400, 600))
#'   bias(obj, mdl = c(200, 400, 600), interval = "both")
#' @export
bias.mcr <- function(x, mdl, level = 0.95,
                      interval = c("", "confidence", "prediction", "both"),
                      ci_method = c("auto", "analytic", "jackknife", "rank",
                        "bootstrap_percentile", "bootstrap_standard"),
                      bootstrap_replicates = 5000L, seed = NULL,
                      ...) {

  if (!inherits(x, "mcr"))
    stop("Input must be an 'mcr' object.")
  if (is.null(x$regression))
    stop("Run regression() first before bias().")

  rg <- x$regression
  interval <- match.arg(interval)
  ci_method <- match.arg(ci_method)
  use_interval <- interval != ""

  if (missing(mdl) || is.null(mdl))
    stop("'mdl' must be provided (medical decision level(s)).")
  if (!is.numeric(mdl))
    stop("'mdl' must be numeric.")

  # Compute bias: predicted candidate - reference
  pred_candidate <- rg$intercept + rg$slope * mdl
  bias_val <- pred_candidate - mdl

  result <- data.frame(mdl = mdl, bias = bias_val)

  # Interval estimation
  uncertainty <- NULL
  if (use_interval) {
    n <- rg$n
    uncertainty <- .regression_interval_samples(
      rg, ci_method, bootstrap_replicates, seed)
    if (!is.null(uncertainty$samples)) {
      draws <- outer(uncertainty$samples[, 2L], mdl, `*`) +
        uncertainty$samples[, 1L] -
        matrix(mdl, nrow(uncertainty$samples), length(mdl), byrow = TRUE)
      se_bias <- apply(draws, 2L, stats::sd)
    } else {
      if (is.null(rg$slope_se) || !is.finite(rg$slope_se))
        stop("The regression object does not contain usable parameter uncertainty.",
             call. = FALSE)
      se_bias <- sqrt(pmax(0, rg$intercept_se^2 +
                      mdl^2 * rg$slope_se^2 + 2 * mdl * rg$cov_xy))
    }
    t_val <- qt(1 - (1 - level) / 2, n - 2)

    if (interval %in% c("confidence", "both")) {
      if (identical(uncertainty$method, "bootstrap_percentile")) {
        alpha <- (1 - level) / 2
        limits <- vapply(seq_len(ncol(draws)), function(j)
          stats::quantile(draws[, j], c(alpha, 1 - alpha), names = FALSE),
          numeric(2))
        result$ci_lower <- limits[1L, ]; result$ci_upper <- limits[2L, ]
      } else {
        result$ci_lower <- bias_val - t_val * se_bias
        result$ci_upper <- bias_val + t_val * se_bias
      }
    }
    if (interval %in% c("prediction", "both")) {
      se_pred <- sqrt(se_bias^2 + rg$sigma^2)
      result$pi_lower <- bias_val - t_val * se_pred
      result$pi_upper <- bias_val + t_val * se_pred
    }
  }

  # Store with metadata attributes
  result <- structure(result,
    class = c("mcr_bias", "data.frame"),
    method   = rg$method,
    interval = interval,
    level    = level,
    ci_method = uncertainty$method %||% NULL,
    bootstrap = uncertainty$bootstrap %||% NULL
  )
  x$bias <- result
  print(result)
  invisible(x)
}


#' Plot (S3 method) -- unified plotting interface
#'
#' Draw one of four plot types for an mcr object.
#'
#' @param x          mcr object
#' @param type       Plot type: "scatter" (scatter + identity line, default),
#'                   "regression" (regression line + confidence/prediction bands),
#'                   "bland_altman" (Bland-Altman plot),
#'                   "bias" (bias trend plot)
#' @param interval   Interval type (only for type="regression"): "none" (default),
#'                   "confidence", "prediction", "both"
#' @param conf.level Confidence level; if NULL, reads from analysis slot, otherwise uses this value
#' @param mdl        Medical decision levels (only for type="bias"); if NULL, tries to read
#'                   from x$bias$mdl
#' @param ...        Additional arguments
#'
#' @return Invisible ggplot object
#'
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   plot(obj)
#'   plot(obj, type = "scatter")
#'   obj <- regression(obj)
#'   plot(obj, type = "regression", interval = "confidence")
#'   obj <- bland_altman(obj)
#'   plot(obj, type = "bland_altman")
#'   obj <- bias(obj, mdl = c(200, 400))
#'   plot(obj, type = "bias")
#' @export
plot.mcr <- function(x,
                     type = c("scatter", "regression", "bland_altman", "bias"),
                     interval = c("none", "confidence", "prediction", "both"),
                     conf.level = NULL, mdl = NULL, ...) {

  if (!inherits(x, "mcr"))
    stop("Input must be an 'mcr' object.")

  type <- match.arg(type)
  interval <- match.arg(interval)
  .require_pkg("ggplot2")

  #
  # Scatter plot + identity line
  #
  if (type == "scatter") {
    df <- data.frame(xv = x$reference, yv = x$candidate)
    x_range <- range(x$candidate, x$reference, na.rm = TRUE)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$xv, y = .data$yv)) +
      ggplot2::geom_point(size = 2, na.rm = TRUE) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "gray50", linewidth = 0.5) +
      ggplot2::coord_fixed(xlim = x_range, ylim = x_range) +
      ggplot2::labs(x = x$reference_name, y = x$candidate_name,
                    title = "Scatter Plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    print(p)
    return(invisible(p))
  }

  #
  # Regression plot
  #
  if (type == "regression") {
    rg <- x$regression
    if (is.null(rg))
      stop("Run regression() first before plot(type = 'regression').", call. = FALSE)

    cl <- conf.level %||% rg$conf.level %||% 0.95

    # Regression line: candidate = intercept + slope * reference
    ref_data <- rg$x_vals %||% x$reference
    cand_data <- rg$y_vals %||% x$candidate
    ref_seq <- seq(min(ref_data), max(ref_data), length.out = 200)
    cand_pred <- rg$intercept + rg$slope * ref_seq

    df <- data.frame(xv = ref_data, yv = cand_data)
    plot_df <- data.frame(x = ref_seq, y = cand_pred)

    show_ci <- interval %in% c("confidence", "both")
    show_pi <- interval %in% c("prediction", "both")
    has_se  <- !is.na(rg$slope_se) && is.finite(rg$slope_se)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$xv, y = .data$yv)) +
      ggplot2::geom_point(size = 2, na.rm = TRUE) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "gray50", linewidth = 0.5)

    if ((show_ci || show_pi) && has_se) {
      n <- rg$n
      t_val <- qt(1 - (1 - cl) / 2, n - 2)
      se_ref <- sqrt(rg$intercept_se^2 +
                     ref_seq^2 * rg$slope_se^2 +
                     2 * ref_seq * rg$cov_xy)

      if (show_ci) {
        ci_lower <- cand_pred - t_val * se_ref
        ci_upper <- cand_pred + t_val * se_ref
        p <- p + ggplot2::geom_ribbon(
          data = data.frame(x = ref_seq, ymin = ci_lower, ymax = ci_upper),
          inherit.aes = FALSE,
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          fill = "steelblue", alpha = 0.2)
      }
      if (show_pi) {
        se_pred <- sqrt(se_ref^2 + rg$sigma^2)
        pi_lower <- cand_pred - t_val * se_pred
        pi_upper <- cand_pred + t_val * se_pred
        p <- p + ggplot2::geom_ribbon(
          data = data.frame(x = ref_seq, ymin = pi_lower, ymax = pi_upper),
          inherit.aes = FALSE,
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          fill = "gray40", alpha = 0.1)
      }
    }

    p <- p +
      ggplot2::geom_line(data = plot_df, inherit.aes = FALSE,
                         ggplot2::aes(x = x, y = y),
                         color = "steelblue", linewidth = 1) +
      ggplot2::labs(x = x$reference_name, y = x$candidate_name,
                    title = paste0(rg$method, " Regression")) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    print(p)
    return(invisible(p))
  }

  #
  # Bland-Altman plot
  #
  if (type == "bland_altman") {
    ba <- x$bland_altman
    if (is.null(ba))
      stop("Run bland_altman() first before plot(type = 'bland_altman').", call. = FALSE)

    x_vals <- ba$x_vals
    y_vals <- ba$y_vals
    mean_diff <- if (!is.null(ba$center)) ba$center else ba$mean_diff
    loa_lower <- ba$loa[1L]
    loa_upper <- ba$loa[2L]

    loa_lower_ci <- if (!is.null(ba$loa_ci$lower)) ba$loa_ci$lower else c(NA, NA)
    loa_upper_ci <- if (!is.null(ba$loa_ci$upper)) ba$loa_ci$upper else c(NA, NA)

    df <- data.frame(xv = x_vals, yv = y_vals)
    x_rng <- range(x_vals, na.rm = TRUE)

    ribbon_df <- data.frame(
      x = c(x_rng[1L], x_rng[2L]),
      loa_lower_lo = loa_lower_ci[1L],
      loa_lower_hi = loa_lower_ci[2L],
      loa_upper_lo = loa_upper_ci[1L],
      loa_upper_hi = loa_upper_ci[2L]
    )
    center_ci <- ba$center_ci %||% c(mean_diff, mean_diff)
    center_ribbon_df <- data.frame(
      x = c(x_rng[1L], x_rng[2L]),
      center_lo = unname(center_ci[1L]),
      center_hi = unname(center_ci[2L])
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$xv, y = .data$yv)) +
      ggplot2::geom_point(size = 2, na.rm = TRUE) +
      ggplot2::geom_ribbon(data = center_ribbon_df, inherit.aes = FALSE,
                           ggplot2::aes(x = .data$x, ymin = .data$center_lo,
                                        ymax = .data$center_hi),
                           fill = "darkorange", alpha = 0.16) +
      ggplot2::geom_ribbon(data = ribbon_df, inherit.aes = FALSE,
                           ggplot2::aes(x = x, ymin = loa_lower_lo, ymax = loa_lower_hi),
                           fill = "steelblue", alpha = 0.15) +
      ggplot2::geom_ribbon(data = ribbon_df, inherit.aes = FALSE,
                           ggplot2::aes(x = x, ymin = loa_upper_lo, ymax = loa_upper_hi),
                           fill = "steelblue", alpha = 0.15) +
      ggplot2::geom_hline(yintercept = mean_diff,
                           color = "steelblue", linewidth = 1) +
      ggplot2::geom_hline(yintercept = loa_lower,
                           linetype = "dashed", color = "steelblue", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = loa_upper,
                           linetype = "dashed", color = "steelblue", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = if (ba$type == "ratio") 1 else 0,
                           linetype = "dotted", color = "gray50", linewidth = 0.5) +
      ggplot2::labs(x = ba$x_label, y = ba$y_label,
                    title = "Bland-Altman Plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    label_y <- max(y_vals, na.rm = TRUE)
    center_name <- switch(ba$method, nonparametric = "Median",
                          hodges_lehmann = "HL pseudomedian", "Mean")
    center_conf_pct <- if (identical(ba$method, "nonparametric"))
      ba$center_ci_coverage * 100 else ba$conf.level * 100
    p <- p + ggplot2::annotate("text",
                                x = min(x_vals, na.rm = TRUE),
                                y = label_y,
                                label = sprintf(
                                  "%s = %s (%.0f%% CI %s to %s)",
                                  center_name, .format_num(mean_diff),
                                  center_conf_pct,
                                  .format_num(center_ci[1L]),
                                  .format_num(center_ci[2L])
                                ),
                                hjust = 0, vjust = 1.5, size = 3.5)
    print(p)
    return(invisible(p))
  }

  #
  # Bias trend plot
  #
  if (type == "bias") {
    rg <- x$regression
    if (is.null(rg))
      stop("Run regression() first before plot(type = 'bias').", call. = FALSE)

    cl <- conf.level %||% rg$conf.level %||% 0.95

    # Cover the entire data range
    x_range <- range(x$reference, na.rm = TRUE)
    x_seq <- seq(x_range[1L], x_range[2L], length.out = 200)
    pred_candidate_seq <- rg$intercept + rg$slope * x_seq
    bias_seq <- pred_candidate_seq - x_seq

    plot_df <- data.frame(x = x_seq, bias = bias_seq)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$bias)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                           color = "gray50", linewidth = 0.5) +
      ggplot2::geom_line(color = "steelblue", linewidth = 1)


    # Confidence / prediction bands
    show_ci <- interval %in% c("confidence", "both")
    show_pi <- interval %in% c("prediction", "both")
    has_se  <- !is.na(rg$slope_se) && is.finite(rg$slope_se)

    if ((show_ci || show_pi) && has_se) {
      n <- rg$n
      se_bias <- sqrt(rg$intercept_se^2 +
                      x_seq^2 * rg$slope_se^2 +
                      2 * x_seq * rg$cov_xy)
      t_val <- qt(1 - (1 - cl) / 2, n - 2)

      if (show_ci) {
        ci_lower <- bias_seq - t_val * se_bias
        ci_upper <- bias_seq + t_val * se_bias
        ribbon_ci <- data.frame(x = x_seq, ymin = ci_lower, ymax = ci_upper)
        p <- p + ggplot2::geom_ribbon(
          data = ribbon_ci,
          inherit.aes = FALSE,
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          fill = "steelblue", alpha = 0.2)
      }
      if (show_pi) {
        se_pred <- sqrt(se_bias^2 + rg$sigma^2)
        pi_lower <- bias_seq - t_val * se_pred
        pi_upper <- bias_seq + t_val * se_pred
        ribbon_pi <- data.frame(x = x_seq, ymin = pi_lower, ymax = pi_upper)
        p <- p + ggplot2::geom_ribbon(
          data = ribbon_pi,
          inherit.aes = FALSE,
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          fill = "gray40", alpha = 0.1)
      }
    }

    # Annotate MDL points: prefer x$bias$mdl, fall back to mdl argument
    mdl_vals <- x$bias$mdl %||% mdl
    if (!is.null(mdl_vals)) {
      bias_at_mdl <- (rg$intercept + rg$slope * mdl_vals) - mdl_vals
      point_df <- data.frame(x = mdl_vals, bias = bias_at_mdl)
      p <- p + ggplot2::geom_point(data = point_df, inherit.aes = FALSE,
                                    ggplot2::aes(x = x, y = bias),
                                    color = "darkorange", size = 3)
    }

    p <- p +
      ggplot2::labs(x = x$reference_name,
                    y = sprintf("Bias (%s - %s)", x$candidate_name, x$reference_name),
                    title = paste0("Bias Plot (", rg$method, ")")) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    print(p)
    return(invisible(p))
  }
}
