#' Outlier detection
#'
#' Integrates Grubbs test, generalized ESD test, Dixon Q test, and IQR method
#' under a consistent API with an S3 print method.
#'
#' @param data  A data frame.
#' @param col  Column name to test (single string).
#' @param method  Detection method: "grubbs" (default), "esd", "dixon", "iqr".
#' @param ...  Additional arguments passed to internal methods:
#'   \describe{
#'     \item{\code{alpha}}{Significance level for grubbs / esd / dixon (default 0.05). dixon only accepts 0.01 and 0.05.}
#'     \item{\code{r}}{Maximum number of outliers for ESD (default \code{min(5, n-2)}).}
#'     \item{\code{type}}{Tail to test for Dixon: \code{"both"} (default), \code{"min"}, \code{"max"}.}
#'     \item{\code{coef}}{IQR multiplier (default 1.5).}
#'   }
#'
#' @return An S3 object of class \code{"outliers_test"} containing:
#'   \describe{
#'     \item{\code{method}}{Method name used.}
#'     \item{\code{data_name}}{Name of the input data.}
#'     \item{\code{col}}{Column name tested.}
#'     \item{\code{n_total}}{Total number of rows.}
#'     \item{\code{n}}{Number of finite observations used.}
#'     \item{\code{missing}}{Number of missing (non-finite) values.}
#'     \item{\code{parameters}}{List of method parameters.}
#'     \item{\code{n_outliers}}{Number of outliers detected.}
#'     \item{\code{indices}}{Indices of outliers in the original vector (1-based).}
#'     \item{\code{values}}{Outlier values.}
#'     \item{\code{details}}{Method-specific details (statistics, critical values, bounds, etc.).}
#'   }
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(x = c(rnorm(20), 10, -8))
#' outliers_test(df, "x", "grubbs")
#' outliers_test(df, "x", "esd", r = 3)
#' outliers_test(df, "x", "dixon")
#' outliers_test(df, "x", "iqr", coef = 2)
#'
#' @export
outliers_test <- function(data, col, method = c("grubbs", "esd", "dixon", "iqr"), ...) {

  # internal: Grubbs test
  .grubbs_test <- function(y, alpha) {
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
      stop("alpha must be a single finite value in (0, 1).")
    valid_idx <- which(is.finite(y))
    y <- y[valid_idx]
    n <- length(y)
    if (n < 3L) {
      return(list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
                  details = list(statistic = NA_real_, critical = NA_real_)))
    }
    mu <- mean(y)
    s  <- stats::sd(y)
    if (s == 0) {
      return(list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
                  details = list(statistic = 0, critical = NA_real_)))
    }
    g  <- max(abs(y - mu)) / s
    t  <- stats::qt(1 - alpha / (2 * n), n - 2)
    gcrit <- ((n - 1L) / sqrt(n)) * sqrt(t^2 / (n - 2 + t^2))

    if (g > gcrit) {
      idx <- which.max(abs(y - mu))
      orig_idx <- valid_idx[idx]
      list(n_outliers = 1L, indices = orig_idx, values = y[idx],
           details = list(statistic = g, critical = gcrit))
    } else {
      list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
           details = list(statistic = g, critical = gcrit))
    }
  }

  # internal: generalized ESD test
  .esd_test <- function(y, r, alpha) {
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
      stop("alpha must be a single finite value in (0, 1).")
    valid_idx <- which(is.finite(y))
    y_orig <- y[valid_idx]
    y <- y_orig
    idx_map <- valid_idx
    n <- length(y)
    if (n < 5L) {
      return(list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
                  details = list(R = numeric(0L), critical = numeric(0L))))
    }
    if (is.null(r)) r <- min(5L, n - 2L)
    if (length(r) != 1L || !is.numeric(r) || r < 1L || r != as.integer(r))
      stop("r must be a positive integer.")
    if (r > n - 2L) r <- n - 2L

    R_vals <- numeric(r)
    crit_vals <- numeric(r)
    idx_removed <- integer(r)

    for (i in seq_len(r)) {
      mu <- mean(y)
      s <- stats::sd(y)
      if (s == 0) break
      dev <- abs(y - mu)
      R_vals[i] <- max(dev) / s
      idx_local <- which.max(dev)
      p <- 1 - alpha / (2 * (n - i + 1))
      t <- stats::qt(p, n - i - 1L)
      crit_vals[i] <- ((n - i) * t) / sqrt((n - i - 1 + t^2) * (n - i + 1))
      idx_removed[i] <- idx_map[idx_local]
      idx_map <- idx_map[-idx_local]
      y <- y[-idx_local]
    }

    n_out <- if (any(R_vals > crit_vals, na.rm = TRUE)) {
      max(which(R_vals > crit_vals))
    } else {
      0L
    }
    list(n_outliers = n_out,
         indices = idx_removed[seq_len(n_out)],
         values  = y_orig[match(idx_removed[seq_len(n_out)], valid_idx)],
         details = list(R = R_vals, critical = crit_vals))
  }

  # internal: Dixon Q test
  .dixon_test <- function(y, alpha, type) {
    type <- match.arg(type, c("both", "min", "max"))
    alpha <- match.arg(as.character(alpha), c("0.01", "0.05"))
    alpha <- as.numeric(alpha)
    valid_idx <- which(is.finite(y))
    y <- y[valid_idx]
    n <- length(y)
    if (n < 3L || n > 30L) {
      return(list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
                  details = list(statistic = NA_real_, critical = NA_real_,
                                 type_detected = NA_character_)))
    }
    q_table <- matrix(c(
      3,  0.941, 0.988,  4,  0.765, 0.889,
      5,  0.642, 0.780,  6,  0.560, 0.698,
      7,  0.507, 0.637,  8,  0.468, 0.590,
      9,  0.437, 0.555,  10, 0.412, 0.527,
      11, 0.392, 0.503,  12, 0.376, 0.482,
      13, 0.361, 0.464,  14, 0.349, 0.449,
      15, 0.338, 0.436,  16, 0.329, 0.425,
      17, 0.320, 0.414,  18, 0.313, 0.405,
      19, 0.306, 0.397,  20, 0.300, 0.389,
      21, 0.295, 0.383,  22, 0.290, 0.377,
      23, 0.286, 0.372,  24, 0.282, 0.367,
      25, 0.278, 0.363,  26, 0.275, 0.359,
      27, 0.272, 0.356,  28, 0.269, 0.353,
      29, 0.266, 0.350,  30, 0.264, 0.347
    ), ncol = 3, byrow = TRUE)
    colnames(q_table) <- c("n", "Q95", "Q99")
    qcrit <- q_table[q_table[, "n"] == n, if (alpha <= 0.01) "Q99" else "Q95"]
    if (length(qcrit) == 0) qcrit <- NA_real_

    ys <- sort(y)
    range_y <- ys[n] - ys[1L]
    if (range_y == 0) {
      return(list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
                  details = list(statistic = 0, critical = qcrit, type_detected = NA_character_)))
    }
    Q_min <- (ys[2L] - ys[1L]) / range_y
    Q_max <- (ys[n] - ys[n - 1L]) / range_y

    if (type %in% c("min", "both") && Q_min > qcrit) {
      orig_idx <- valid_idx[which(y == ys[1L])[1L]]
      return(list(n_outliers = 1L, indices = orig_idx, values = y[which(y == ys[1L])[1L]],
                  details = list(statistic = Q_min, critical = qcrit, type_detected = "min")))
    }
    if (type %in% c("max", "both") && Q_max > qcrit) {
      orig_idx <- valid_idx[which(y == ys[n])[1L]]
      return(list(n_outliers = 1L, indices = orig_idx, values = y[which(y == ys[n])[1L]],
                  details = list(statistic = Q_max, critical = qcrit, type_detected = "max")))
    }
    list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
         details = list(statistic = max(Q_min, Q_max), critical = qcrit, type_detected = NA_character_))
  }

  # internal: IQR method
  .iqr_test <- function(y, coef) {
    if (length(coef) != 1L || !is.finite(coef) || coef < 0)
      stop("coef must be a single finite non-negative number.")
    valid_idx <- which(is.finite(y))
    y <- y[valid_idx]
    n <- length(y)
    if (n < 3L) {
      return(list(n_outliers = 0L, indices = integer(0L), values = numeric(0L),
                  details = list(lower_bound = NA_real_, upper_bound = NA_real_)))
    }
    qs    <- stats::quantile(y, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr   <- qs[2L] - qs[1L]
    lower <- qs[1L] - coef * iqr
    upper <- qs[2L] + coef * iqr
    local_idx <- which(y < lower | y > upper)
    orig_idx <- valid_idx[local_idx]
    list(n_outliers = length(orig_idx), indices = orig_idx, values = y[local_idx],
         details = list(lower_bound = lower, upper_bound = upper))
  }

  # dispatch
  method <- match.arg(method)
  col <- as.character(col)

  # input validation
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)

  if (!is.character(col) || length(col) != 1L)
    stop("'col' must be a single column name.", call. = FALSE)

  if (!col %in% names(data))
    stop("Column '", col, "' not found in 'data'.", call. = FALSE)

  if (!is.numeric(data[[col]]))
    stop("Column '", col, "' must be numeric.", call. = FALSE)

  data_name <- deparse(substitute(data))
  y <- data[[col]]

  dots <- list(...)
  alpha <- if (!is.null(dots$alpha)) dots$alpha else 0.05

  result <- switch(method,
    grubbs = .grubbs_test(y, alpha = alpha),
    esd    = .esd_test(y, alpha = alpha, r = dots$r),
    dixon  = .dixon_test(y, alpha = alpha,
                         type = if (!is.null(dots$type)) dots$type else "both"),
    iqr    = .iqr_test(y, coef = if (!is.null(dots$coef)) dots$coef else 1.5)
  )

  # build S3 object (internal fns no longer return parameters)
  n_total <- length(y)
  esd_r <- if (!is.null(dots$r)) dots$r else min(5L, max(1L, length(which(is.finite(y))) - 2L))
  params <- switch(method,
    grubbs = list(alpha = alpha),
    esd    = list(alpha = alpha, r = esd_r),
    dixon  = list(alpha = alpha, type = if (!is.null(dots$type)) dots$type else "both"),
    iqr    = list(coef = if (!is.null(dots$coef)) dots$coef else 1.5)
  )

  n <- length(which(is.finite(y)))
  n_missing <- sum(!is.finite(y))

  out <- list(
    method     = method,
    data_name  = data_name,
    col        = col,
    n_total    = n_total,
    n          = n,
    missing    = n_missing,
    parameters = params,
    n_outliers = result$n_outliers,
    indices    = result$indices,
    values     = result$values,
    details    = result$details
  )
  class(out) <- "outliers_test"
  out
}

# S3 print method

#' Print an outlier detection result
#'
#' @param x An \code{outliers_test} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.outliers_test <- function(x, ...) {
  method_labels <- c(
    grubbs = "Grubbs Outliers Test",
    esd    = "Generalized ESD Outliers Test",
    dixon  = "Dixon Q Outliers Test",
    iqr    = "IQR Outliers Test"
  )

  cat(sprintf("%s\n", method_labels[x$method]))
  cat(sprintf("  Data:   %s\n", x$data_name))
  cat(sprintf("  Column: %s (n = %d, missing = %d)\n", x$col, x$n,
              if (is.null(x$missing)) 0L else x$missing))

  # parameters
  param_str <- switch(x$method,
    grubbs = sprintf("alpha = %.2f", x$parameters$alpha),
    esd    = sprintf("alpha = %.2f, r = %d", x$parameters$alpha, x$parameters$r),
    dixon  = sprintf("alpha = %.2f, type = %s", x$parameters$alpha, x$parameters$type),
    iqr    = sprintf("coef = %.1f", x$parameters$coef)
  )
  cat(sprintf("  Parameters: %s\n", param_str))

  cat(sprintf("  Outliers: %d\n", x$n_outliers))

  if (x$n_outliers > 0L) {
    cat("\n")
    cat(sprintf("  %-10s  %10s", "Index", "Value"))

    # method-specific extra columns
    extra_cols <- character(0L)
    switch(x$method,
      grubbs = { extra_cols <- c("G stat", "Critical") },
      esd    = { extra_cols <- c("R stat", "Critical") },
      dixon  = { extra_cols <- c("Q stat", "Critical", "Tail") },
      iqr    = { extra_cols <- c() }
    )
    if (length(extra_cols) > 0L) {
      for (ec in extra_cols) {
        cat(sprintf("  %10s", ec))
      }
    }
    cat("\n")
    cat(sprintf("  %s\n", paste(rep("-", 10 + 11 * (1 + length(extra_cols))), collapse = "")))

    for (i in seq_len(x$n_outliers)) {
      idx <- x$indices[i]
      val <- x$values[i]
      cat(sprintf("  %-10d  %10.4f", idx, val))

      switch(x$method,
        grubbs = {
          cat(sprintf("  %10.4f  %10.4f", x$details$statistic, x$details$critical))
        },
        esd = {
          cat(sprintf("  %10.4f  %10.4f", x$details$R[i], x$details$critical[i]))
        },
        dixon = {
          cat(sprintf("  %10.4f  %10.4f  %10s",
            x$details$statistic, x$details$critical, x$details$type_detected))
        },
        iqr = {
          # IQR has no per-row extra columns
        }
      )
      cat("\n")
    }

    # IQR method shows fences
    if (x$method == "iqr") {
      cat(sprintf("\n  fences: [%.4f, %.4f]\n",
        x$details$lower_bound, x$details$upper_bound))
    }
  } else {
    cat("\n  No outliers detected.\n")
  }

  invisible(x)
}


#' Normality test
#'
#' Performs normality tests on a single numeric column of a data frame.
#' Plots are generated separately via \code{plot()}.
#'
#' @param data  A data frame.
#' @param col  Column name to test (single string).
#' @param method  Test method: "auto" (default, chooses by sample size),
#'   "shapiro"/"sw", "ad", "lillie", "cvm".
#' @param level  Confidence level (default 0.95, used for QQ confidence bands).
#'
#' @return An S3 object of class \code{"normal_test"} containing:
#'   \describe{
#'     \item{\code{call}}{Function call.}
#'     \item{\code{data_name}}{Data object name.}
#'     \item{\code{col}}{Column name.}
#'     \item{\code{method}}{User-specified method.}
#'     \item{\code{level}}{Confidence level.}
#'     \item{\code{n_total}}{Total number of rows.}
#'     \item{\code{n}}{Number of finite values (non-missing).}
#'     \item{\code{missing}}{Number of missing values.}
#'     \item{\code{test_method}}{Actual test method used.}
#'     \item{\code{statistic}}{Test statistic.}
#'     \item{\code{p_value}}{p-value.}
#'     \item{\code{y}}{Clean numeric vector, used by \code{plot()}.}
#'   }
#'
#' @examples
#'   df <- data.frame(x = rnorm(50))
#'   normal_test(df, "x")
#'   normal_test(df, "x", method = "ad")
#'   plot(normal_test(df, "x"))
#'   plot(normal_test(df, "x"), type = "his")
#'   plot(normal_test(df, "x"), type = c("qq", "his"))
#' @importFrom nortest ad.test lillie.test cvm.test
#' @export
normal_test <- function(data, col, method = c("auto", "shapiro", "sw", "ad", "lillie", "cvm"),
                        level = 0.95) {

  # internal helpers

  .select_method <- function(n) {
    if (n < 3L)       return(NA_character_)
    if (n <= 50L)     return("shapiro")
    if (n <= 5000L)   return("ad")
    return("lillie")
  }

  .shapiro_test <- function(y) {
    st <- tryCatch(stats::shapiro.test(y), error = function(e) NULL)
    if (is.null(st)) return(list(statistic = NA_real_, p_value = NA_real_))
    list(statistic = unname(st$statistic), p_value = st$p.value)
  }

  .ad_test <- function(y) {
    st <- tryCatch(nortest::ad.test(y), error = function(e) NULL)
    if (is.null(st)) return(list(statistic = NA_real_, p_value = NA_real_))
    list(statistic = unname(st$statistic), p_value = st$p.value)
  }

  .lillie_test <- function(y) {
    st <- tryCatch(nortest::lillie.test(y), error = function(e) NULL)
    if (is.null(st)) return(list(statistic = NA_real_, p_value = NA_real_))
    list(statistic = unname(st$statistic), p_value = st$p.value)
  }

  .cvm_test <- function(y) {
    st <- tryCatch(nortest::cvm.test(y), error = function(e) NULL)
    if (is.null(st)) return(list(statistic = NA_real_, p_value = NA_real_))
    list(statistic = unname(st$statistic), p_value = st$p.value)
  }

  # input validation
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)

  if (!is.character(col) || length(col) != 1L)
    stop("'col' must be a single column name.", call. = FALSE)

  if (!col %in% names(data))
    stop("Column '", col, "' not found in 'data'.", call. = FALSE)

  if (!is.numeric(data[[col]]))
    stop("Column '", col, "' must be numeric.", call. = FALSE)

  method <- match.arg(method)

  # level validation
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single numeric value between 0 and 1.", call. = FALSE)

  # data cleaning
  data_name <- deparse(substitute(data))
  y <- data[[col]]
  n_total <- length(y)
  miss <- is.na(y) | !is.finite(y)
  n_miss <- sum(miss)
  y_clean <- y[!miss]
  n <- length(y_clean)

  # run test
  statistic <- NA_real_
  p_value <- NA_real_
  test_method <- NA_character_

  if (n >= 3L && stats::sd(y_clean) > 0) {
    test_method <- if (method == "auto") .select_method(n) else method

    if (!is.na(test_method)) {
      result <- switch(test_method,
        shapiro = .shapiro_test(y_clean),
        sw      = .shapiro_test(y_clean),
        ad      = .ad_test(y_clean),
        lillie  = .lillie_test(y_clean),
        cvm     = .cvm_test(y_clean)
      )
      statistic <- result$statistic
      p_value <- result$p_value
    }
  }

  # build S3 object
  out <- list(
    call        = match.call(),
    data_name   = data_name,
    col         = col,
    method      = method,
    level       = level,
    n_total     = n_total,
    n           = n,
    missing     = n_miss,
    test_method = test_method,
    statistic   = statistic,
    p_value     = p_value,
    y           = y_clean
  )
  class(out) <- "normal_test"
  out
}


# S3 print method

#' Print a normality test result
#'
#' @param x A \code{normal_test} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.normal_test <- function(x, ...) {
  .format_num <- function(xx, digits = 4L) {
    formatC(xx, format = "f", digits = digits)
  }
  .method_label <- function(method, n) {
    if (is.na(method)) {
      if (n < 3L) return("Insufficient data")
      return("Constant data")
    }
    switch(method,
      shapiro = "Shapiro-Wilk", sw = "Shapiro-Wilk",
      ad = "Anderson-Darling", lillie = "Lilliefors", cvm = "Cramer-von Mises",
      method
    )
  }

  method_label <- .method_label(x$test_method, x$n)

  stat_name <- switch(x$test_method,
    shapiro = "W", sw = "W", ad = "A", lillie = "D", cvm = "W"
  )

  cat("\nNormality Test\n")
  cat(sprintf("  Data:   %s\n", x$data_name))
  cat(sprintf("  Column: %s (n = %d, missing = %d)\n", x$col, x$n, x$missing))
  cat(sprintf("  Method:     %s\n\n", method_label))

  if (!is.na(x$statistic) && !is.na(x$p_value)) {
    stat_str <- .format_num(x$statistic, 4L)
    pv_str <- .format_num(x$p_value, 4L)
    cat(sprintf("  %s = %s, p = %s\n", stat_name, stat_str, pv_str))
  } else {
    cat("  Test could not be performed.\n")
  }

  invisible(x)
}

# S3 plot method

#' Plot normality diagnostics
#'
#' @param x A \code{normal_test} object.
#' @param type One or both of \code{"qq"} and \code{"his"}.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly; requested plots are drawn as a
#'   side effect. For fewer than three observations no plot is drawn.
#' @examples
#' d <- data.frame(value = rnorm(30))
#' nt <- normal_test(d, "value")
#' plot(nt, type = c("qq", "his"))
#' @import ggplot2
#' @export
plot.normal_test <- function(x, type = c("qq", "his"), ...) {
  type <- match.arg(type, several.ok = TRUE)

  .method_label <- function(method, n) {
    if (is.na(method)) {
      if (n < 3L) return("Insufficient data")
      return("Constant data")
    }
    switch(method,
      shapiro = "Shapiro-Wilk", sw = "Shapiro-Wilk",
      ad = "Anderson-Darling", lillie = "Lilliefors", cvm = "Cramer-von Mises",
      method
    )
  }

  y <- x$y
  n <- length(y)
  if (n < 3L) {
    cat("Insufficient data for plotting.\n")
    return(invisible(x))
  }

  method_label <- .method_label(x$test_method, n)

  if ("qq" %in% type) {
    y_sorted <- sort(y)
    p <- (seq_len(n) - 0.5) / n
    theoretical <- stats::qnorm(p)

    q_y <- stats::quantile(y, c(0.25, 0.75), names = FALSE)
    q_t <- stats::qnorm(c(0.25, 0.75))
    slope <- diff(q_y) / diff(q_t)
    intercept <- q_y[1L] - slope * q_t[1L]
    predicted <- intercept + slope * theoretical

    z <- stats::qnorm(1 - (1 - x$level) / 2)
    se <- slope * sqrt(p * (1 - p) / n) / stats::dnorm(theoretical)
    upper <- predicted + z * se
    lower <- predicted - z * se

    qq_data <- data.frame(theoretical = theoretical, observed = y_sorted,
                          predicted = predicted, lower = lower, upper = upper,
                          stringsAsFactors = FALSE)

    ann_label <- if (!is.na(x$p_value))
      sprintf("%s p = %.4f", method_label, x$p_value)
    else
      "N too small"
    ann_df <- data.frame(label = ann_label, stringsAsFactors = FALSE)

    p_qq <- ggplot2::ggplot(qq_data, ggplot2::aes(x = .data$theoretical)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = 0.2, fill = "steelblue", na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = .data$lower),
        color = "steelblue", linewidth = 0.4, linetype = "dotted", na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = .data$upper),
        color = "steelblue", linewidth = 0.4, linetype = "dotted", na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = .data$predicted),
        color = "steelblue", linewidth = 0.6, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = .data$observed),
        size = 2, na.rm = TRUE) +
      ggplot2::geom_text(data = ann_df,
        ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
        hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE) +
      ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
           title = paste0("Normal Q-Q Plot -- ", x$col)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(p_qq)
  }

  if ("his" %in% type) {
    z <- (y - mean(y)) / stats::sd(y)
    plot_data <- data.frame(z = z, stringsAsFactors = FALSE)

    seq_z <- seq(-4, 4, length.out = 200)
    norm_curve <- data.frame(z = seq_z, density = stats::dnorm(seq_z))

    p_his <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$z)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = 30, fill = "steelblue", alpha = 0.4, color = "white", na.rm = TRUE) +
      ggplot2::geom_density(color = "steelblue", linewidth = 0.5,
        na.rm = TRUE) +
      ggplot2::geom_line(data = norm_curve,
        ggplot2::aes(x = .data$z, y = .data$density),
        color = "darkorange", linewidth = 0.5, na.rm = TRUE) +
      ggplot2::labs(x = "Standardized value (z)", y = "Density",
           title = paste0("Distribution Plot -- ", x$col)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(p_his)
  }

  invisible(x)
}


#' Estimate bias using reference materials
#'
#' Calculates measurement bias and its interval according to Section 5.3 of
#' YY/T 1789.2-2021. For each reference-material level, the uncertainty of the
#' measured mean is calculated as \code{s / sqrt(n)} and expanded by \code{k}.
#' This expanded uncertainty is then combined with the reference-material
#' uncertainty. The function reports results only and applies no acceptance
#' criterion.
#'
#' @param data A data frame containing replicate measurement results.
#' @param result Name of the numeric measurement-result column.
#' @param reference Reference-material assigned value. Supply a finite scalar,
#'   a numeric vector with one value per row or level, or the name of a numeric
#'   column in \code{data}.
#' @param reference_uncertainty Expanded uncertainty of the reference material,
#'   at the same coverage level as the uncertainty calculated with \code{k}.
#'   Supply a non-negative scalar, a numeric vector with one value per row or
#'   level, or the name of a numeric column in \code{data}.
#' @param level Optional name of the reference-material level column. When
#'   omitted, all rows are treated as one level.
#' @param k Positive coverage factor used to expand the standard uncertainty of
#'   the measured mean. The standard's example uses \code{k = 2}.
#'
#' @return An object of class \code{"reference_bias"}. Its \code{result}
#'   component contains the sample size, mean, standard deviation, uncertainty,
#'   bias, relative bias, and interval for each reference-material level.
#'
#' @references YY/T 1789.2-2021, Section 5.3 and Appendix A.
#' @examples
#' reference_data <- data.frame(
#'   level = rep(c("low", "high"), each = 6),
#'   result = c(9.8, 10.1, 10.0, 9.9, 10.2, 10.0,
#'              19.7, 20.1, 20.0, 19.9, 20.2, 20.1)
#' )
#' bias_result <- reference_bias(
#'   reference_data, result = "result",
#'   reference = c(10, 20),
#'   reference_uncertainty = c(0.2, 0.3),
#'   level = "level"
#' )
#' bias_result$result
#' @export
reference_bias <- function(data, result, reference, reference_uncertainty,
                           level = NULL, k = 2) {
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(result) || length(result) != 1L || !nzchar(result))
    stop("'result' must be a single column name.", call. = FALSE)
  if (!result %in% names(data))
    stop("Column '", result, "' not found in 'data'.", call. = FALSE)
  if (!is.numeric(data[[result]]))
    stop("Column '", result, "' must be numeric.", call. = FALSE)
  if (nrow(data) == 0L)
    stop("'data' must contain at least one row.", call. = FALSE)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0)
    stop("'k' must be a single positive finite number.", call. = FALSE)

  if (is.null(level)) {
    level_values <- rep("Overall", nrow(data))
    level_name <- NULL
  } else {
    if (!is.character(level) || length(level) != 1L || !nzchar(level))
      stop("'level' must be NULL or a single column name.", call. = FALSE)
    if (!level %in% names(data))
      stop("Column '", level, "' not found in 'data'.", call. = FALSE)
    if (anyNA(data[[level]]))
      stop("Column '", level, "' must not contain missing values.", call. = FALSE)
    level_values <- as.character(data[[level]])
    level_name <- level
  }
  level_order <- unique(level_values)

  .resolve_input <- function(expr, value, argument, nonnegative = FALSE) {
    expr_name <- if (is.symbol(expr)) as.character(expr) else NULL
    if (!is.null(expr_name) && expr_name %in% names(data)) {
      value <- data[[expr_name]]
    } else {
      value <- force(value)
      if (is.character(value) && length(value) == 1L && value %in% names(data))
        value <- data[[value]]
    }
    if (!is.numeric(value))
      stop("'", argument, "' must be numeric or name a numeric column.",
           call. = FALSE)
    if (length(value) == 1L) {
      value <- rep(value, nrow(data))
    } else if (length(value) == length(level_order)) {
      value <- value[match(level_values, level_order)]
    } else if (length(value) != nrow(data)) {
      stop("'", argument, "' must have length 1, the number of levels, or the number of rows.",
           call. = FALSE)
    }
    if (any(!is.finite(value)))
      stop("'", argument, "' must contain only finite values.", call. = FALSE)
    if (nonnegative && any(value < 0))
      stop("'", argument, "' must not contain negative values.", call. = FALSE)
    value
  }

  reference_values <- .resolve_input(
    substitute(reference), reference, "reference")
  uncertainty_values <- .resolve_input(
    substitute(reference_uncertainty), reference_uncertainty,
    "reference_uncertainty", nonnegative = TRUE)

  measured <- data[[result]]
  valid <- is.finite(measured)
  rows <- lapply(level_order, function(current_level) {
    in_level <- level_values == current_level
    ref <- unique(reference_values[in_level])
    ref_u <- unique(uncertainty_values[in_level])
    if (length(ref) != 1L)
      stop("'reference' must be constant within level '", current_level, "'.",
           call. = FALSE)
    if (length(ref_u) != 1L)
      stop("'reference_uncertainty' must be constant within level '",
           current_level, "'.", call. = FALSE)

    y <- measured[in_level & valid]
    n <- length(y)
    if (n < 2L)
      stop("Level '", current_level,
           "' must contain at least two finite measurement results.",
           call. = FALSE)
    measured_mean <- mean(y)
    measured_sd <- stats::sd(y)
    standard_u <- measured_sd / sqrt(n)
    expanded_u <- k * standard_u
    bias <- measured_mean - ref
    half_width <- sqrt(expanded_u^2 + ref_u^2)

    data.frame(
      level = current_level,
      n_total = sum(in_level),
      n = n,
      excluded = sum(in_level & !valid),
      mean = measured_mean,
      sd = measured_sd,
      standard_uncertainty = standard_u,
      expanded_uncertainty = expanded_u,
      reference = ref,
      reference_uncertainty = ref_u,
      bias = bias,
      relative_bias = if (ref == 0) NA_real_ else 100 * bias / ref,
      interval_half_width = half_width,
      lower = bias - half_width,
      upper = bias + half_width,
      stringsAsFactors = FALSE
    )
  })

  out <- list(
    call = match.call(),
    data_name = deparse(substitute(data)),
    result_name = result,
    level_name = level_name,
    k = k,
    n_total = nrow(data),
    n = sum(valid),
    excluded = sum(!valid),
    result = do.call(rbind, rows)
  )
  rownames(out$result) <- NULL
  class(out) <- "reference_bias"
  out
}


#' Print reference-material bias results
#'
#' @param x A \code{reference_bias} object.
#' @param ... Reserved arguments.
#'
#' @return \code{x}, invisibly.
#' @examples
#' d <- data.frame(result = c(9.8, 10.0, 10.1, 9.9, 10.2, 10.0))
#' print(reference_bias(d, "result", reference = 10,
#'                      reference_uncertainty = 0.2))
#' @export
print.reference_bias <- function(x, ...) {
  if (!inherits(x, "reference_bias"))
    stop("'x' must be a reference_bias object.", call. = FALSE)

  measurement_table <- x$result[c(
    "level", "n", "mean", "sd", "reference", "reference_uncertainty",
    "standard_uncertainty", "expanded_uncertainty"
  )]
  names(measurement_table) <- c(
    "Level", "n", "Mean", "SD", "Reference", "U_ref", "u_X", "U_X"
  )
  bias_table <- x$result[c(
    "level", "bias", "relative_bias", "interval_half_width", "lower", "upper"
  )]
  names(bias_table) <- c(
    "Level", "Bias", "Bias (%)", "Interval half-width", "Lower", "Upper"
  )

  cat("\nReference Material Bias\n")
  cat("  Standard: YY/T 1789.2-2021, Section 5.3\n")
  cat(sprintf("  Data:     %s\n", x$data_name))
  cat(sprintf("  Result:   %s\n", x$result_name))
  cat(sprintf("  Levels:   %d\n", nrow(x$result)))
  cat(sprintf("  Records:  %d valid, %d excluded\n", x$n, x$excluded))
  cat(sprintf("  Coverage factor: k = %s\n\n",
              format(x$k, digits = 6L, trim = TRUE)))
  cat("  Measurement and uncertainty\n")
  .print_kv_rows(measurement_table, list(
    "Measurement summary" = c(
      "n" = "n", "Mean" = "Mean", "SD" = "SD", "Reference" = "Reference"),
    "Uncertainty" = c("U_ref" = "U_ref", "u_X" = "u_X", "U_X" = "U_X")
  ), row_title = c("Level" = "Level"))
  cat("\n  Bias interval\n")
  .print_kv_rows(bias_table, list(
    "Bias estimate" = c("Bias" = "Bias", "Bias (%)" = "Bias (%)"),
    "Bias interval" = c(
      "Interval half-width" = "Interval half-width", "Lower" = "Lower",
      "Upper" = "Upper")
  ), row_title = c("Level" = "Level"))
  cat("\n  Interval: Bias +/- sqrt(U_X^2 + U_ref^2)\n")

  invisible(x)
}
