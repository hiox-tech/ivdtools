#' Method Qualitative Analysis -- 2x2 Fourfold Table (Candidate Method x Reference Method)
#'
#' Functionality: Raw data descriptive statistics, contingency frequency table generation, fourfold table formatted printing.
#'
#' @importFrom stats qnorm binom.test pchisq sd median quantile
#' @noRd
NULL

#' Descriptive Statistics: Raw Data Overview (S3 Method)
#'
#' Print variable description information for a fourfold_table object (from raw_to_table).
#' If the object comes from counts_to_table (no raw data), indicates unavailability.
#'
#' @param x   a fourfold_table object
#' @param ... reserved arguments
#'
#' @return Updated fourfold_table object, with $describe slot filled with description results
#'
#' @export
#' @examples
#'   df <- data.frame(
#'     new    = c("positive","positive","negative","negative","positive","negative"),
#'     gold   = c("positive","negative","positive","negative","positive","positive"),
#'     age    = c(45, 52, 38, 61, 47, 55),
#'     gender = c("M","F","F","M","M","F"),
#'     id     = c("S001","S002","S003","S004","S005","S006")
#'   )
#'   tab <- raw_to_table(df, "new", "gold", id = "id")
#'   tab <- describe(tab)
describe.fourfold_table <- function(x, ...) {
  if (!inherits(x, "fourfold_table"))
    stop("Input must be a fourfold_table object (from raw_to_table / counts_to_table)")

  # If constructed from counts, no raw data to describe
  if (!is.null(x$print$source) && x$print$source == "counts") {
    cat("\n  No description available for count-derived table.\n")
    cat("  Use raw_to_table() with raw data to enable data description.\n\n")
    x$describe <- list(source = "counts",
      note = "No description available for count-derived table.")
    return(invisible(x))
  }

  data <- x$print$data_raw
  if (is.null(data))
    stop("Raw data not found in object (required for describe).")

  n_row <- nrow(data)
  n_col <- ncol(data)

  # Identify candidate and reference method names (excluding response column)
  id_var    <- x$print$var_names$id

  # Determine which columns to analyze: only exclude id column
  exclude_cols <- if (!is.null(id_var)) id_var else character(0)
  cols <- setdiff(names(data), exclude_cols)

  # ID duplication information
  id_dup <- NULL
  if (!is.null(id_var)) {
    if (id_var %in% names(data)) {
      id_vec <- data[[id_var]]
      dups   <- duplicated(id_vec) | duplicated(id_vec, fromLast = TRUE)
      n_dup  <- sum(dups)
      n_unique <- length(unique(id_vec))
      id_dup <- list(n_dup = n_dup, n_unique = n_unique)
    }
  }

  # Per-column summary
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
      if (length(na_row)) freqs$Pct[na_row] <- round(freqs$Freq[na_row] / n_total * 100, 1)

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
          type   = "numeric",
          n_na   = n_na,
          pct_na = pct_na,
          mean   = NA_real_, sd = NA_real_, median = NA_real_,
          q25 = NA_real_, q75 = NA_real_, min = NA_real_, max = NA_real_
        )
      } else {
        columns[[col]] <- list(
          type   = "numeric",
          n_na   = n_na,
          pct_na = pct_na,
          mean   = mean(ok), sd = sd(ok), median = median(ok),
          q25    = quantile(ok, 0.25, names = FALSE),
          q75    = quantile(ok, 0.75, names = FALSE),
          min    = min(ok), max = max(ok)
        )
      }
    } else {
      columns[[col]] <- list(type = "other", n_na = n_na, pct_na = pct_na)
    }
  }

  # Missing value summary
  missing_cols <- list()
  for (col in names(columns)) {
    ci <- columns[[col]]
    if (ci$n_na > 0) {
      na_idx <- which(is.na(data[[col]]))
      missing_cols[[col]] <- list(n_na = ci$n_na, pct_na = ci$pct_na, indices = na_idx)
    }
  }

  describe_result <- list(
    n_row        = n_row,
    n_col        = n_col,
    id_dup       = id_dup,
    columns      = columns,
    missing_cols = missing_cols
  )
  class(describe_result) <- "describe_table"
  x$describe <- describe_result

  print(x$describe)

  invisible(x)
}

#' S3 print method for describe_table
#' @param x describe_table object
#' @param ... reserved arguments
#' @return Invisible describe_table object
#' @export
#' @examples
#'   df <- data.frame(
#'     new    = c("positive","positive","negative","negative","positive","negative"),
#'     gold   = c("positive","negative","positive","negative","positive","positive"),
#'     age    = c(45, 52, 38, 61, 47, 55),
#'     gender = c("M","F","F","M","M","F"),
#'     id     = c("S001","S002","S003","S004","S005","S006")
#'   )
#'   tab <- raw_to_table(df, "new", "gold", id = "id")
#'   tab <- describe(tab)
print.describe_table <- function(x, ...) {
  n_row <- x$n_row
  columns <- x$columns
  id_dup <- x$id_dup
  missing_cols <- x$missing_cols

  cat("\n")
  cat("Data Summary\n")
  cat(sprintf("  Rows:    %d\n", n_row))
  cat(sprintf("  Columns: %d\n", x$n_col))

  if (!is.null(id_dup)) {
    cat(sprintf("  ID:      %d unique", id_dup$n_unique))
    if (id_dup$n_dup > 0) {
      cat(sprintf(", %d duplicated (%.1f%%)\n", id_dup$n_dup,
                  id_dup$n_dup / n_row * 100))
    } else {
      cat(", no duplicates\n")
    }
  }

  if (length(missing_cols) > 0) {
    cat("\n  Missing values:\n")
    for (col in names(missing_cols)) {
      m <- missing_cols[[col]]
      cat(sprintf("    %s: %d / %d (%.1f%%)", col, m$n_na, n_row, m$pct_na))
      if (length(m$indices) <= 10)
        cat(sprintf("  row(s): %s", paste(m$indices, collapse = ", ")))
      else
        cat(sprintf("  first 10 rows: %s ...", paste(m$indices[1:10], collapse = ", ")))
      cat("\n")
    }
  }

  cat("\n")

  for (col in names(columns)) {
    ci <- columns[[col]]
    na_str <- if (ci$n_na > 0) sprintf(" (%d NA, %.1f%%)", ci$n_na, ci$pct_na) else ""

    if (ci$type == "categorical") {
      cat(sprintf("  %s%s\n", col, na_str))
      for (r in seq_len(nrow(ci$table))) {
        lv <- as.character(ci$table$Level[r])
        if (is.na(lv)) lv <- "<NA>"
        cat(sprintf("    %-12s %5d  (%5.1f%%)\n",
                    lv, ci$table$Freq[r], ci$table$Pct[r]))
      }
    } else if (ci$type == "numeric") {
      cat(sprintf("  %s%s\n", col, na_str))
      if (is.na(ci$mean)) {
        cat("    (all missing)\n")
      } else {
        cat(sprintf("    Mean (SD):       %.2f (%.2f)\n", ci$mean, ci$sd))
        cat(sprintf("    Median (Q1-Q3):  %.2f (%.2f - %.2f)\n", ci$median, ci$q25, ci$q75))
        cat(sprintf("    Range:           %.2f - %.2f\n", ci$min, ci$max))
      }
    } else {
      cat(sprintf("  %s (other)%s\n", col, na_str))
    }
  }
  invisible(x)
}


#' Core function: raw data -> frequency fourfold table (silent construction)
#'
#' @param data        Data frame, one observation per row
#' @param candidate   Candidate method variable name (character), e.g. "method_new"
#' @param reference   Reference method variable name (character), e.g. "method_standard"
#' @param id          ID variable name (character), optional, used for duplicate detection in describe
#' @param na.rm       Whether to remove NA, default TRUE
#' @param positive    Specify the positive label (character); if non-NULL, levels will be uniformly mapped to positive/negative.
#'                   If NULL, the last sorted level is treated as positive.
#'
#' @return A fourfold_table object (list with S3 class "fourfold_table"), containing:
#'         - $freq_raw        Raw frequencies (long format)
#'         - $table           Fourfold matrix (with margins)
#'         - $n               Total sample size
#'         - $print           Print slot (descriptive metadata)
#'         - $candidate_levels  Levels of the candidate method
#'         - $reference_levels  Levels of the reference method
#'         Use print() directly to display the fourfold table.
#'
#' @export
#' @examples
#'   df <- data.frame(
#'     new    = c("positive","positive","negative","negative","positive","negative"),
#'     gold   = c("positive","negative","positive","negative","positive","positive"),
#'     age    = c(45, 52, 38, 61, 47, 55),
#'     gender = c("M","F","F","M","M","F"),
#'     id     = c("S001","S002","S003","S004","S005","S006")
#'   )
#'   result <- raw_to_table(df, "new", "gold", id = "id")
raw_to_table <- function(data, candidate, reference, id = NULL, na.rm = TRUE,
                         positive = NULL) {

  # Input validation
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!candidate %in% names(data)) stop("Variable '", candidate, "' not found in the data frame (candidate)")
  if (!reference %in% names(data)) stop("Variable '", reference, "' not found in the data frame (reference)")
  if (!is.null(positive) && (!is.character(positive) || length(positive) != 1 || is.na(positive)))
    stop("'positive' must be a single string or NULL")

  # Extraction and factor conversion
  cand_vec <- data[[candidate]]
  ref_vec  <- data[[reference]]

  miss_cand <- is.na(cand_vec)
  miss_ref  <- is.na(ref_vec)
  has_na    <- miss_cand | miss_ref
  n_na      <- sum(has_na)

  # Missing details (for print slot)
  miss_ids <- list()
  if (n_na > 0) {
    if (na.rm) {
      cand_vec <- cand_vec[!has_na]
      ref_vec  <- ref_vec[!has_na]
      miss_ids$candidate <- which(miss_cand)
      miss_ids$reference <- which(miss_ref)
      miss_ids$both      <- which(miss_cand & miss_ref)
    } else {
      warning(
        n_na, " row(s) contain NA values and will be dropped by table(). ",
        "Set na.rm = TRUE to see details, or pre-process your data to remove NAs.",
        call. = FALSE
      )
    }
  }

  cand_vec <- factor(cand_vec)
  ref_vec  <- factor(ref_vec)

  if (nlevels(cand_vec) != 2 || nlevels(ref_vec) != 2)
    stop("Both candidate and reference must have exactly 2 distinct levels to form a 2x2 table",
         "\n  candidate '", candidate, "' has ", nlevels(cand_vec), " level(s)",
         "\n  reference '", reference, "' has ", nlevels(ref_vec), " level(s)")

  # Put the positive level first so downstream metrics consistently read
  # [1, 1] as TP. With no explicit label, match the package default used by
  # roc(): the last sorted level is treated as positive.
  if (!is.null(positive)) {
    for (nm in c("cand_vec", "ref_vec")) {
      vec <- get(nm)
      lv  <- levels(vec)
      if (!positive %in% lv)
        stop(sprintf("Positive label '%s' not found in data", positive))

      neg <- setdiff(lv, positive)
      new_lv <- c("positive", "negative")
      names(new_lv) <- c(positive, neg)
      vec <- factor(new_lv[as.character(vec)], levels = new_lv)
      assign(nm, vec)
    }
  } else {
    cand_vec <- factor(cand_vec, levels = rev(sort(levels(cand_vec))))
    ref_vec  <- factor(ref_vec, levels = rev(sort(levels(ref_vec))))
  }

  # Generate frequency table
  freq_raw <- as.data.frame(
    table(cand_vec, ref_vec, dnn = c(candidate, reference)),
    responseName = "Freq"
  )

  tbl        <- table(cand_vec, ref_vec, dnn = c(candidate, reference))
  tbl_margin <- addmargins(tbl)

  if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0))
    stop("Degenerate table: one or more margins are zero (insufficient data in some categories)")

  # print slot (descriptive metadata)
  n_total <- nrow(data)
  n <- sum(tbl)

  # Duplicate ID check
  n_dup <- 0L
  dup_ids <- integer(0L)
  if (!is.null(id)) {
    if (id %in% names(data)) {
      id_vec <- data[[id]]
      dups <- duplicated(id_vec) | duplicated(id_vec, fromLast = TRUE)
      n_dup <- sum(dups)
      if (n_dup > 0) dup_ids <- which(dups)
    }
  }

  print_info <- list(
    data_class  = class(data)[1L],
    n           = n,
    n_total     = n_total,
    n_miss      = n_total - n,
    n_dup       = n_dup,
    source      = "raw",
    var_names   = list(candidate = candidate, reference = reference, id = id),
    dup_ids     = dup_ids,
    miss_ids    = miss_ids,
    data_raw    = data
  )

  # Return (silent)
  out <- structure(list(
    freq_raw         = freq_raw,
    table            = tbl_margin,
    n                = n,
    print            = print_info,
    describe         = NULL,
    candidate_levels = levels(cand_vec),
    reference_levels = levels(ref_vec),
    diagnostics      = NULL,
    kappa            = NULL,
    mcnemar          = NULL
  ), class = "fourfold_table")

  invisible(out)
}


#' Construct a 2x2 contingency table from TP / FP / TN / FN frequencies (silent construction)
#'
#' When raw data are not available and only the four diagnostic test frequencies are known,
#' use this function to directly construct a result object with the same structure as
#' raw_to_table(), compatible with print() for displaying the fourfold table.
#'
#' @param tp               True positives (candidate=level1, reference=level1)
#' @param fp               False positives (candidate=level1, reference=level2)
#' @param tn               True negatives (candidate=level2, reference=level2)
#' @param fn               False negatives (candidate=level2, reference=level1)
#' @param candidate        Candidate method name (character)
#' @param reference        Reference method name (character)
#' @param candidate_levels Two levels of the candidate method, default c("positive", "negative")
#' @param reference_levels Two levels of the reference method, default c("positive", "negative")
#'
#' @return A fourfold_table object (S3 class "fourfold_table") with the same structure as raw_to_table():
#'         - $freq_raw         Long-format frequency data.frame
#'         - $table            3x3 matrix with margins
#'         - $n                Total sample size
#'         - $print            Print slot
#'         - $candidate_levels Levels of the candidate method
#'         - $reference_levels Levels of the reference method
#'         Use print() directly to display the fourfold table.
#'
#' @export
#' @examples
#'   result <- counts_to_table(tp = 45, fp = 12, tn = 80, fn = 5,
#'                           candidate = "new method", reference = "gold standard")
counts_to_table <- function(tp, fp, tn, fn,
                          candidate = "candidate",
                          reference = "reference",
                          candidate_levels = c("positive", "negative"),
                          reference_levels = c("positive", "negative")) {

  # Input validation
  counts <- list(tp = tp, fp = fp, tn = tn, fn = fn)
  for (nm in names(counts)) {
    val <- counts[[nm]]
    if (!is.numeric(val) || length(val) != 1 || is.na(val) || is.infinite(val))
      stop(nm, " must be a single finite numeric value")
    if (val < 0 || val != round(val))
      stop(nm, " must be a non-negative integer")
  }

  if (!is.character(candidate) || length(candidate) != 1)
    stop("'candidate' must be a single string")
  if (!is.character(reference) || length(reference) != 1)
    stop("'reference' must be a single string")
  if (length(candidate_levels) != 2)
    stop("'candidate_levels' must have exactly 2 levels")
  if (length(reference_levels) != 2)
    stop("'reference_levels' must have exactly 2 levels")

  # Construct frequency table
  cand_vec <- factor(
    rep(candidate_levels, each = 2),
    levels = candidate_levels
  )
  ref_vec <- factor(
    rep(reference_levels, times = 2),
    levels = reference_levels
  )

  tbl <- table(cand_vec, ref_vec, dnn = c(candidate, reference))
  tbl[] <- c(tp, fn, fp, tn)

  if (sum(tbl) == 0)
    stop("Table has zero total count (all four frequencies are 0)")
  if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0))
    stop("Degenerate table: one or more margins are zero (insufficient data in some categories)")

  freq_raw <- as.data.frame(tbl, responseName = "Freq")
  tbl_margin <- addmargins(tbl)
  n <- sum(c(tp, fp, tn, fn))

  # print slot (from counts, no raw data information)
  print_info <- list(
    data_class = "counts",
    n          = n,
    source     = "counts",
    var_names  = list(candidate = candidate, reference = reference),
    n_dup      = NA_integer_,
    dup_ids    = NULL,
    miss_ids   = NULL
  )

  out <- structure(list(
    freq_raw         = freq_raw,
    table            = tbl_margin,
    n                = n,
    print            = print_info,
    describe         = NULL,
    candidate_levels = candidate_levels,
    reference_levels = reference_levels,
    diagnostics      = NULL,
    kappa            = NULL,
    mcnemar          = NULL
  ), class = "fourfold_table")

  invisible(out)
}


#' Print helper: pad with spaces to specified width
#' @param s character string
#' @param width target width
#' @keywords internal
#' @noRd
pad_w <- function(s, width) {
  s <- as.character(s)
  cur <- nchar(s, type = "width")
  if (cur < width) paste0(s, strrep(" ", width - cur)) else s
}


#' S3 print method: format a \code{fourfold_table} object as a fourfold table
#'
#' Prints data source, total sample size, variable names, fourfold table (with margins), and analysis status.
#'
#' @param x \code{fourfold_table} object (result from \code{raw_to_table} / \code{counts_to_table})
#' @param margin Whether to display margins, default \code{TRUE}
#' @param ... reserved arguments
#'
#' @return Invisible \code{fourfold_table} object
#'
#' @export
#' @examples
#'   df <- data.frame(
#'     new    = c("positive","positive","negative","negative","positive","negative"),
#'     gold   = c("positive","negative","positive","negative","positive","positive"),
#'     age    = c(45, 52, 38, 61, 47, 55),
#'     gender = c("M","F","F","M","M","F"),
#'     id     = c("S001","S002","S003","S004","S005","S006")
#'   )
#'   tab <- raw_to_table(df, "new", "gold", id = "id")
#'   print(tab)
print.fourfold_table <- function(x, margin = TRUE, ...) {

  if (!inherits(x, "fourfold_table"))
    stop("Input must be a fourfold_table object (from raw_to_table / counts_to_table)")

  p <- x$print

  # $print slot overview
  cat("\nFourfold Table\n")

  if (p$source == "raw") {
    cat(sprintf("  %-18s: %s\n", "Source", "raw data"))
    cat(sprintf("  %-18s: %d\n", "Total n", p$n))
    if (p$n_miss > 0) cat(sprintf("  %-18s: %d\n", "Excluded (NA)", p$n_miss))

    # Duplicate IDs
    if (!is.null(p$var_names$id)) {
      if (p$n_dup > 0) {
        dup_msg <- sprintf("%d", p$n_dup)
        if (length(p$dup_ids) > 0)
          dup_msg <- sprintf("%d  row(s): %s", p$n_dup, paste(p$dup_ids, collapse = ", "))
        cat(sprintf("  %-18s: %s\n", "Duplicate IDs", dup_msg))
      } else {
        cat(sprintf("  %-18s: %s\n", "Duplicate IDs", "none"))
      }
    }
  } else {
    cat(sprintf("  %-18s: %s\n", "Source", "counts"))
    cat(sprintf("  %-18s: %d\n", "Total n", p$n))
  }

  cat(sprintf("  %-18s: %s\n", "Candidate", p$var_names$candidate))
  cat(sprintf("  %-18s: %s\n", "Reference", p$var_names$reference))

  # Fourfold table
  cand_name <- names(dimnames(x$table))[1]
  ref_name  <- names(dimnames(x$table))[2]

  mat     <- if (margin) x$table else x$table[1:2, 1:2]
  ref_lv  <- colnames(mat)
  cand_lv <- rownames(mat)

  col_w <- 8
  cand_label <- cand_name
  cand_w <- nchar(cand_label, type = "width")
  lv_w <- max(nchar(c(cand_lv, "sum"), type = "width"))

  cat("\n  2x2 Contingency Table\n")

  sep <- paste(rep("-", 2 + cand_w + 1 + lv_w + 2 + col_w * length(ref_lv)), collapse = "")
  cat(" ", sep, "\n")

  cat(" ", pad_w("", cand_w), " ", pad_w("", lv_w), " ", pad_w(ref_name, (col_w + 1) * length(ref_lv)), "\n", sep = "")

  cat(" ", pad_w("", cand_w), " ", pad_w("", lv_w), " ")
  for (lv in ref_lv) cat(pad_w(lv, col_w), " ")
  cat("\n")

  cat(" ", sep, "\n")

  for (i in seq_len(nrow(mat))) {
    left <- if (i == 1) cand_label else ""
    right <- if (i == nrow(mat)) "sum" else cand_lv[i]
    cat(" ", pad_w(left, cand_w), " ", pad_w(right, lv_w), " ")
    for (j in seq_len(ncol(mat))) {
      cat(pad_w(mat[i, j], col_w), " ")
    }
    cat("\n")
  }

  cat(" ", sep, "\n")

  # Analysis status
  has_desc <- !is.null(x$describe)
  has_diag <- !is.null(x$diagnostics)
  has_kap  <- !is.null(x$kappa)
  has_mcn  <- !is.null(x$mcnemar)

  if (has_desc || has_diag || has_kap || has_mcn) {
    cat("\nAnalyses performed:\n\n")
    if (has_desc) cat("  [x] Describe\n")  else cat("  [ ] Describe\n")
    if (has_diag) cat("  [x] Diagnostics\n") else cat("  [ ] Diagnostics\n")
    if (has_kap)  cat("  [x] Kappa\n")      else cat("  [ ] Kappa\n")
    if (has_mcn)  cat("  [x] McNemar\n")    else cat("  [ ] McNemar\n")
  }
  invisible(x)
}


#' Single proportion confidence interval (internal helper)
#'
#' Supports multiple methods: wald, wald-cc, wilson, wilson-cc, agresti-coull, jeffreys, clopper-pearson.
#' @param x number of successes
#' @param n number of trials
#' @param conf.level confidence level, default 0.95
#' @param method CI method, default "wilson"
#' @return named vector with lower and upper bounds
#' @keywords internal
#' @noRd
.conf_prop_ci <- function(x, n, conf.level = 0.95,
                          method = c("wilson", "wald", "wald-cc", "agresti-coull",
                                     "jeffreys", "wilson-cc", "clopper-pearson")) {

  method <- match.arg(method)
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  if (n < 0 || x < 0 || x > n)
    return(c(lower = NA_real_, upper = NA_real_))

  p <- x / n
  a <- 1 - conf.level
  z <- qnorm(1 - a / 2)

  lower <- upper <- NA_real_

  switch(method,
    wald = {
      se <- sqrt(p * (1 - p) / n)
      lower <- p - z * se
      upper <- p + z * se
    },
    `wald-cc` = {
      se <- sqrt(p * (1 - p) / n)
      lower <- max(0, p - z * se - 1 / (2 * n))
      upper <- min(1, p + z * se + 1 / (2 * n))
    },
    wilson = {
      denom <- 1 + z^2 / n
      center <- (p + z^2 / (2 * n)) / denom
      se <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
      lower <- center - se
      upper <- center + se
    },
    `wilson-cc` = {
      denom <- 1 + z^2 / n
      center <- (p + z^2 / (2 * n)) / denom
      se <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
      cc <- 1 / (2 * n)
      lower <- max(0, center - se - cc)
      upper <- min(1, center + se + cc)
    },
    `agresti-coull` = {
      n_tilde <- n + z^2
      p_tilde <- (x + z^2 / 2) / n_tilde
      se <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)
      lower <- p_tilde - z * se
      upper <- p_tilde + z * se
    },
    jeffreys = {
      if (x == 0) {
        lower <- 0
        upper <- qbeta(1 - a / 2, 0.5, n + 0.5)
      } else if (x == n) {
        lower <- qbeta(a / 2, n + 0.5, 0.5)
        upper <- 1
      } else {
        lower <- qbeta(a / 2,     x + 0.5, n - x + 0.5)
        upper <- qbeta(1 - a / 2, x + 0.5, n - x + 0.5)
      }
    },
    `clopper-pearson` = {
      if (x == 0) {
        lower <- 0
        upper <- qbeta(1 - a / 2, 1, n)
      } else if (x == n) {
        lower <- qbeta(a / 2, n, 1)
        upper <- 1
      } else {
        lower <- qbeta(a / 2,     x,     n - x + 1)
        upper <- qbeta(1 - a / 2, x + 1, n - x)
      }
    }
  )

  c(lower = max(0, lower), upper = min(1, upper))
}


#' Diagnostic performance metrics: calculate sensitivity, specificity, likelihood ratios, etc. with confidence intervals for a fourfold_table
#'
#' @param object      a fourfold_table object
#' @param conf.level  Confidence level, default 0.95
#' @param ci.method   Proportion confidence interval method, default "wilson".
#'                   Supports "wald", "wald-cc", "wilson", "wilson-cc",
#'                   "agresti-coull", "jeffreys", "clopper-pearson"
#' @param prevalence  User-specified prevalence. If provided, PPV/NPV will be calculated
#'                    based on this prevalence via Bayes' theorem (rather than the observed
#'                    prevalence in the data).
#' @param ...         Additional arguments passed to print
#'
#' @return Updated fourfold_table object, with diagnostic performance metrics stored in $diagnostics:
#'         - $tp / $fp / $fn / $tn / $n   Four-fold frequencies and total sample size
#'         - $sensitivity    Sensitivity with CI
#'         - $specificity    Specificity with CI
#'         - $ppv            Positive predictive value with CI
#'         - $npv            Negative predictive value with CI
#'         - $accuracy       Accuracy with CI
#'         - $prevalence     Prevalence with CI
#'         - $lr_positive    Positive likelihood ratio with CI
#'         - $lr_negative    Negative likelihood ratio with CI
#'         - $odds_ratio     Odds ratio with CI
#'         - $conf_level     Confidence level used
#'         - $ci_method      Method name used
#'         - $prev_specified Whether user specified prevalence
#'
#' @export
#' @examples
#'   tab <- raw_to_table(ivd_qualitative_example, "new", "gold",
#'                       positive = "positive", id = "id")
#'   tab <- diagnostics(tab)
#'   tab <- diagnostics(tab, prevalence = 0.3)
#'   tab <- diagnostics(tab, conf.level = 0.90)
#'   tab <- diagnostics(tab, ci.method = "clopper-pearson")
diagnostics.fourfold_table <- function(object, conf.level = 0.95,
                                       ci.method = "wilson",
                                       prevalence = NULL, ...) {

  if (!inherits(object, "fourfold_table"))
    stop("Input must be a fourfold_table object (from raw_to_table / counts_to_table)")
  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single value in (0, 1)")

  mat <- object$table
  if (nrow(mat) != 3 || ncol(mat) != 3)
    stop("'table' must be a 3x3 matrix with margins")

  ci.method <- match.arg(ci.method,
    c("wilson", "wald", "wald-cc", "agresti-coull",
      "jeffreys", "wilson-cc", "clopper-pearson"))

  a <- mat[1, 1]  # TP
  b <- mat[1, 2]  # FP
  c <- mat[2, 1]  # FN
  d <- mat[2, 2]  # TN
  n <- a + b + c + d

  # Proportion confidence intervals (generic method)
  prop_ci <- function(x, m) {
    .conf_prop_ci(x, m, conf.level = conf.level, method = ci.method)
  }

  # Sensitivity: TP / (TP + FN)
  sens <- a / (a + c)
  sens_ci <- prop_ci(a, a + c)

  # Specificity: TN / (FP + TN)
  spec <- d / (b + d)
  spec_ci <- prop_ci(d, b + d)

  # Positive predictive value / Negative predictive value
  if (!is.null(prevalence)) {
    stopifnot(is.numeric(prevalence), length(prevalence) == 1,
              prevalence > 0, prevalence < 1)
    p_prev <- prevalence

    ppv <- sens * p_prev / (sens * p_prev + (1 - spec) * (1 - p_prev))
    npv <- spec * (1 - p_prev) / ((1 - sens) * p_prev + spec * (1 - p_prev))

    sens_low <- unname(sens_ci["lower"])
    sens_up  <- unname(sens_ci["upper"])
    spec_low <- unname(spec_ci["lower"])
    spec_up  <- unname(spec_ci["upper"])

    ppv_lower <- if (is.na(sens_low) || is.na(spec_low)) NA_real_
      else sens_low * p_prev / (sens_low * p_prev + (1 - spec_low) * (1 - p_prev))
    ppv_upper <- if (is.na(sens_up) || is.na(spec_up)) NA_real_
      else sens_up * p_prev / (sens_up * p_prev + (1 - spec_up) * (1 - p_prev))

    npv_lower <- if (is.na(sens_low) || is.na(spec_low)) NA_real_
      else spec_low * (1 - p_prev) / ((1 - sens_low) * p_prev + spec_low * (1 - p_prev))
    npv_upper <- if (is.na(sens_up) || is.na(spec_up)) NA_real_
      else spec_up * (1 - p_prev) / ((1 - sens_up) * p_prev + spec_up * (1 - p_prev))

    ppv_ci <- c(lower = ppv_lower, upper = ppv_upper)
    npv_ci <- c(lower = npv_lower, upper = npv_upper)

    prev <- p_prev
    prev_ci <- c(lower = NA_real_, upper = NA_real_)
    prev_specified <- TRUE

  } else {
    ppv <- a / (a + b)
    ppv_ci <- prop_ci(a, a + b)
    npv <- d / (c + d)
    npv_ci <- prop_ci(d, c + d)
    prev <- (a + c) / n
    prev_ci <- prop_ci(a + c, n)
    prev_specified <- FALSE
  }

  # Accuracy: (TP + TN) / n
  acc <- (a + d) / n
  acc_ci <- prop_ci(a + d, n)

  # Likelihood ratio confidence intervals (log scale)
  q <- qnorm(1 - (1 - conf.level) / 2)

  # Positive likelihood ratio: sens / (1 - spec)
  lr_pos <- sens / (1 - spec)
  if (is.finite(lr_pos) && lr_pos > 0 && a > 0 && b > 0) {
    se_ln_lr_pos <- sqrt(1 / a - 1 / (a + c) + 1 / b - 1 / (b + d))
    lr_pos_ci <- exp(log(lr_pos) + c(-1, 1) * q * se_ln_lr_pos)
  } else {
    lr_pos_ci <- c(lower = NA_real_, upper = NA_real_)
  }
  names(lr_pos_ci) <- c("lower", "upper")

  # Negative likelihood ratio: (1 - sens) / spec
  lr_neg <- (1 - sens) / spec
  if (is.finite(lr_neg) && lr_neg > 0 && c > 0 && d > 0) {
    se_ln_lr_neg <- sqrt(1 / c - 1 / (a + c) + 1 / d - 1 / (b + d))
    lr_neg_ci <- exp(log(lr_neg) + c(-1, 1) * q * se_ln_lr_neg)
  } else {
    lr_neg_ci <- c(lower = NA_real_, upper = NA_real_)
  }
  names(lr_neg_ci) <- c("lower", "upper")

  # Odds ratio: (TP * TN) / (FP * FN)
  or <- (a * d) / (b * c)
  if (is.finite(or) && or > 0 && a > 0 && b > 0 && c > 0 && d > 0) {
    se_ln_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    or_ci <- exp(log(or) + c(-1, 1) * q * se_ln_or)
  } else {
    or_ci <- c(lower = NA_real_, upper = NA_real_)
  }
  names(or_ci) <- c("lower", "upper")

  result <- list(
    tp = a, fp = b, fn = c, tn = d, n = n,
    sensitivity  = c(est = sens, sens_ci),
    specificity  = c(est = spec, spec_ci),
    ppv          = c(est = ppv, ppv_ci),
    npv          = c(est = npv, npv_ci),
    accuracy     = c(est = acc, acc_ci),
    prevalence   = c(est = prev, prev_ci),
    lr_positive  = c(est = lr_pos, lr_pos_ci),
    lr_negative  = c(est = lr_neg, lr_neg_ci),
    odds_ratio   = c(est = or, or_ci),
    conf_level   = conf.level,
    ci_method    = ci.method,
    prev_specified = prev_specified
  )
  class(result) <- "diagnostics"
  object$diagnostics <- result

  # Auto print
  print(object$diagnostics)
  invisible(object)
}


#' S3 print method: format \code{diagnostics} output (diagnostic performance metrics)
#'
#' Prints four-fold table frequencies, sensitivity, specificity, PPV, NPV, accuracy, likelihood ratios, etc.
#'
#' @param x \code{diagnostics} object
#' @param digits Number of decimal places, default \code{3}
#' @param ... reserved arguments
#'
#' @return Invisible \code{diagnostics} object
#'
#' @export
#' @examples
#'   df <- data.frame(
#'     new    = c("positive","positive","negative","negative","positive","negative"),
#'     gold   = c("positive","negative","positive","negative","positive","positive"),
#'     age    = c(45, 52, 38, 61, 47, 55),
#'     gender = c("M","F","F","M","M","F"),
#'     id     = c("S001","S002","S003","S004","S005","S006")
#'   )
#'   tab <- raw_to_table(df, "new", "gold", positive = "positive", id = "id")
#'   tab <- diagnostics(tab)
print.diagnostics <- function(x, digits = 3, ...) {

  if (!inherits(x, "diagnostics"))
    stop("Input must be a diagnostics object (from diagnostics.fourfold_table)")

  # Title
  cat("\nDiagnostic Performance Summary\n")

  # Metrics list
  fmt_est_ci <- function(label, vec) {
    est <- vec["est"]
    lo  <- vec["lower"]
    up  <- vec["upper"]
    if (is.na(lo) || is.na(up)) {
      sprintf("  %-18s %s", label, format(round(est, digits), nsmall = digits))
    } else {
      sprintf("  %-18s %s  (%s, %s)",
              label,
              format(round(est, digits), nsmall = digits),
              format(round(lo, digits),  nsmall = digits),
              format(round(up, digits),  nsmall = digits))
    }
  }

  cat(fmt_est_ci("Sensitivity",     x$sensitivity),  "\n")
  cat(fmt_est_ci("Specificity",     x$specificity),  "\n")
  cat(fmt_est_ci("PPV",             x$ppv),          "\n")
  cat(fmt_est_ci("NPV",             x$npv),          "\n")
  cat(fmt_est_ci("Accuracy",        x$accuracy),     "\n")
  cat(fmt_est_ci("Prevalence",      x$prevalence),   "\n")
  cat(fmt_est_ci("LR+",             x$lr_positive),  "\n")
  cat(fmt_est_ci("LR-",             x$lr_negative),  "\n")
  cat(fmt_est_ci("Odds Ratio",      x$odds_ratio),   "\n")

  # Footnote comments (separated by blank line above)
  notes <- c()
  if (x$prev_specified) {
    notes <- c(notes, sprintf(
      "PPV & NPV are based on user-specified prevalence of %.3f.", x$prevalence["est"]))
  } else {
    notes <- c(notes, sprintf(
      "PPV & NPV are based on data prevalence of %.3f.", x$prevalence["est"]))
  }
  notes <- c(notes, sprintf(
    "Confidence intervals for proportions use '%s' method.", x$ci_method))
  cat("\n")
  for (n in notes) cat(sprintf("  # %s\n", n))
  invisible(x)
}


#' S3 summary method
#'
#' Summarizes all key analysis information already executed in a fourfold_table object,
#' formatted following summary.mcr() pattern (mcr.R).
#'
#' @param object a fourfold_table object
#' @param digits Number of decimal places, default 3
#' @param ...    Reserved arguments
#'
#' @return Invisible fourfold_table object
#'
#' @export
#' @examples
#'   tab <- counts_to_table(tp = 85, fp = 3, tn = 90, fn = 2,
#'                          candidate = "new method", reference = "gold standard")
#'   summary(tab)
#'   tab <- diagnostics(tab); tab <- kappa(tab); tab <- mcnemar(tab)
#'   summary(tab)
summary.fourfold_table <- function(object, digits = 3, ...) {

  if (!inherits(object, "fourfold_table"))
    stop("Input must be a fourfold_table object (from raw_to_table / counts_to_table)")

  cat("\nQualitative Analysis - summary\n")
  cat("-------------------------------------------------\n")

    # Full print output (fourfold table)
  print(object)

  # Analysis status
  has_desc <- !is.null(object$describe) && !is.null(object$describe$columns)
  has_diag <- !is.null(object$diagnostics)
  has_kap  <- !is.null(object$kappa)
  has_mcn  <- !is.null(object$mcnemar)

  # Each analysis output via its own print method
  if (has_desc) print(object$describe)
  if (has_diag) print(object$diagnostics)
  if (has_kap)  print(object$kappa)
  if (has_mcn)  print(object$mcnemar)

  invisible(object)
}


#' Cohen's Kappa / PABAK Agreement Analysis (2x2 table only)
#'
#' Calculate Cohen's Kappa coefficient or PABAK (Prevalence-Adjusted Bias-Adjusted Kappa)
#' and its confidence interval, used to evaluate agreement between two binary classification methods.
#'
#' When data prevalence is extreme, Cohen's Kappa may be low even when observed agreement is high.
#' In such cases, specify the \code{prevalence} parameter to compute PABAK as a correction.
#' PABAK = 2 x p_obs - 1, which assumes expected agreement = 0.5.
#'
#' @param x           a fourfold_table object
#' @param prevalence  Patient prevalence. If provided, calculates PABAK instead of Cohen's Kappa.
#' @param conf.level  Confidence level, default 0.95
#' @param ...         Additional arguments (currently unused).
#'
#' @return Updated fourfold_table object, with Kappa results stored in $kappa:
#'         - $kappa         Kappa / PABAK coefficient
#'         - $se            Standard error
#'         - $lower / $upper  Lower/upper confidence interval bounds
#'         - $p_obs         Observed agreement
#'         - $p_exp         Expected agreement
#'         - $n             Total sample size
#'         - $conf_level    Confidence level used
#'         - $method        Method name used
#'
#' @export
#' @examples
#'   tab <- counts_to_table(tp = 85, fp = 3, tn = 90, fn = 2,
#'                          candidate = "new method", reference = "gold standard")
#'   tab <- kappa(tab)                        # Cohen's Kappa
#'   tab <- kappa(tab, prevalence = 0.5)      # PABAK
kappa.fourfold_table <- function(x, prevalence = NULL, conf.level = 0.95, ...) {

  if (!inherits(x, "fourfold_table"))
    stop("x must be a fourfold_table object")
  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single value in (0, 1)")

  mat <- x$table[1:2, 1:2, drop = FALSE]

  a <- mat[1, 1]; b <- mat[1, 2]
  c <- mat[2, 1]; d <- mat[2, 2]
  n <- a + b + c + d
  stopifnot(n > 0)

  p_obs <- (a + d) / n   # Observed agreement

  if (!is.null(prevalence)) {
    # PABAK: correction for Kappa underestimation due to extreme prevalence
    stopifnot(is.numeric(prevalence), length(prevalence) == 1,
              prevalence > 0, prevalence < 1)
    kappa_val <- 2 * p_obs - 1
    var_kappa <- 4 * p_obs * (1 - p_obs) / n
    p_exp <- 0.5
    method_label <- "PABAK"
  } else {
    # Cohen's Kappa
    p_exp <- ((a + b) * (a + c) + (c + d) * (b + d)) / n^2
    kappa_val <- (p_obs - p_exp) / (1 - p_exp)

    p_mat <- mat / n
    p_row <- rowSums(p_mat)
    p_col <- colSums(p_mat)

    A <- sum(p_mat[1, 1] * (1 - p_exp - (p_row[1] + p_col[1]) * (1 - p_obs))^2,
             p_mat[2, 2] * (1 - p_exp - (p_row[2] + p_col[2]) * (1 - p_obs))^2)
    B <- (1 - p_obs)^2 * (p_mat[1, 2] * (p_row[1] + p_col[2])^2 +
                          p_mat[2, 1] * (p_row[2] + p_col[1])^2)
    C <- (p_obs * p_exp - 2 * p_exp + p_obs)^2
    var_kappa <- (A + B - C) / (n * (1 - p_exp)^4)
    method_label <- "Cohen's Kappa"
  }

  se    <- sqrt(max(0, var_kappa))
  q     <- qnorm(1 - (1 - conf.level) / 2)
  lower <- max(-1, kappa_val - q * se)
  upper <- min(1,  kappa_val + q * se)

  # Construct results
  cand_name <- x$print$var_names$candidate
  ref_name  <- x$print$var_names$reference

  result <- structure(list(
    kappa      = kappa_val,
    se         = se,
    lower      = lower,
    upper      = upper,
    p_obs      = p_obs,
    p_exp      = p_exp,
    n          = n,
    conf.level = conf.level,
    method     = method_label,
    prevalence = prevalence,
    cand_name  = cand_name,
    ref_name   = ref_name
  ), class = "kappa_table")

  x$kappa <- result

  # Auto print
  print(x$kappa)
  invisible(x)
}


#' S3 print method: format \code{kappa_table} output
#'
#' Prints Cohen's Kappa / PABAK coefficient, standard error, confidence interval, and agreement rates.
#'
#' @param x \code{kappa_table} object
#' @param digits Number of decimal places, default \code{4}
#' @param ... reserved arguments
#'
#' @return Invisible \code{kappa_table} object
#'
#' @export
#' @examples
#'   tab <- counts_to_table(tp = 85, fp = 3, tn = 90, fn = 2,
#'                          candidate = "new method", reference = "gold standard")
#'   tab <- kappa(tab)
print.kappa_table <- function(x, digits = 4, ...) {

  if (!inherits(x, "kappa_table"))
    stop("Input must be a kappa_table object")

  cl <- x$conf.level * 100

  cat("\n", x$method, "\n", sep = "")
  cat(sprintf("  %s  vs  %s\n", x$cand_name, x$ref_name))
  cat(sprintf("  n = %d\n\n", x$n))

  cat(sprintf("  %-20s  %.4f\n", "Kappa", x$kappa))
  cat(sprintf("  %-20s  %.4f\n", "SE", x$se))
  cat(sprintf("  %-20s  (%.4f, %.4f)\n", paste0(cl, "% CI"), x$lower, x$upper))
  cat(sprintf("  %-20s  %.4f\n", "Observed agreement", x$p_obs))
  cat(sprintf("  %-20s  %.4f\n", "Expected agreement", x$p_exp))

  if (x$method == "PABAK") {
    cat(sprintf("  # PABAK (prevalence = %.3f): 2 * p_obs - 1, assumes expected agreement = 0.5\n", x$prevalence))
  }

  invisible(x)
}


#' McNemar Test (for fourfold_table only)
#'
#' Perform McNemar's test on a 2x2 paired contingency table to evaluate whether
#' there is a systematic difference between two binary classification methods
#' (i.e., whether marginal probabilities are equal).
#'
#' When discordant pairs (b + c) are few (< 25), automatically uses the exact binomial test
#' (binom.test); otherwise uses McNemar's chi-squared test with continuity correction.
#'
#' @param x           a fourfold_table object
#' @param conf.level  Confidence level, default 0.95
#' @param ...         Additional arguments (currently unused).
#'
#' @return Updated fourfold_table object, with McNemar results stored in $mcnemar:
#'         - $chi_sq      McNemar chi-squared statistic (NA for exact test)
#'         - $df          Degrees of freedom (NA for exact test)
#'         - $p_value     Test p-value
#'         - $method      Test method used
#'         - $b           Discordant pairs (positive, negative) count
#'         - $c           Discordant pairs (negative, positive) count
#'         - $n_bc        Total discordant pairs
#'         - $conf.level  Confidence level used
#'
#' @export
#' @examples
#'   tab <- counts_to_table(tp = 45, fp = 12, tn = 80, fn = 5,
#'                          candidate = "new method", reference = "gold standard")
#'   tab <- mcnemar(tab)
mcnemar.fourfold_table <- function(x, conf.level = 0.95, ...) {

  if (!inherits(x, "fourfold_table"))
    stop("x must be a fourfold_table object")
  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single value in (0, 1)")

  mat <- x$table[1:2, 1:2, drop = FALSE]

  b <- mat[1, 2]  # FP
  c <- mat[2, 1]  # FN
  n_bc <- b + c

  cand_name <- x$print$var_names$candidate
  ref_name  <- x$print$var_names$reference

  # No discordant pairs -> perfect agreement, p-value = 1
  if (n_bc == 0) {
    p_value <- 1
    method  <- "No discordant pairs (perfect agreement)"
    chi_sq  <- NA_real_
    df      <- NA_integer_
  } else if (n_bc < 25) {
    bt <- binom.test(b, n_bc, p = 0.5, conf.level = conf.level)
    p_value <- bt$p.value
    method  <- "Exact binomial test"
    chi_sq  <- NA_real_
    df      <- NA_integer_
  } else {
    # McNemar chi-squared test with continuity correction
    chi_sq <- (abs(b - c) - 1)^2 / n_bc
    p_value <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
    method  <- "McNemar's chi-squared test with continuity correction"
    df      <- 1L
  }

  result <- structure(list(
    chi_sq     = chi_sq,
    df         = df,
    p_value    = p_value,
    method     = method,
    b          = b,
    c          = c,
    n_bc       = n_bc,
    conf.level = conf.level,
    cand_name  = cand_name,
    ref_name   = ref_name
  ), class = "mcnemar_table")

  x$mcnemar <- result

  # Auto print
  print(x$mcnemar)
  invisible(x)
}


#' S3 print method: format \code{mcnemar_table} output
#'
#' Prints McNemar test results, including discordant pair counts, test statistic, and p-value.
#'
#' @param x \code{mcnemar_table} object
#' @param ... reserved arguments
#'
#' @return Invisible \code{mcnemar_table} object
#'
#' @export
#' @examples
#'   tab <- counts_to_table(tp = 45, fp = 12, tn = 80, fn = 5,
#'                          candidate = "new method", reference = "gold standard")
#'   tab <- mcnemar(tab)
print.mcnemar_table <- function(x, ...) {

  if (!inherits(x, "mcnemar_table"))
    stop("Input must be a mcnemar_table object")

  cl <- x$conf.level * 100

  cat("\nMcNemar Test\n")
  cat(sprintf("  %s  vs  %s\n", x$cand_name, x$ref_name))
  cat(sprintf("  Discordant pairs:  b (pos, neg) = %d,  c (neg, pos) = %d\n", x$b, x$c))
  cat(sprintf("  n (b + c)        = %d\n\n", x$n_bc))

  if (grepl("perfect agreement", x$method)) {
    cat(sprintf("  Method: %s\n", x$method))
    cat(sprintf("  P-value:          %.4f\n", x$p_value))
  } else if (x$method == "Exact binomial test") {
    cat(sprintf("  Method: %s\n", x$method))
    cat(sprintf("  P-value:          %.4f\n", x$p_value))
  } else {
    cat(sprintf("  Method: %s\n", x$method))
    cat(sprintf("  McNemar chi-squared = %.3f, df = %d\n", x$chi_sq, x$df))
    cat(sprintf("  P-value:             %.4f\n", x$p_value))
  }

  sig <- if (x$p_value < 1 - x$conf.level) "Significant" else "Not significant"
  cat(sprintf("  Conclusion: %s at %.0f%% confidence level\n", sig, cl))
  invisible(x)
}
