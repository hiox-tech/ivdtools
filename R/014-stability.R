# Shared helper functions

#' @noRd
check_column <- function(data, col, arg, numeric = FALSE) {
  if (!is.character(col) || length(col) != 1L || !col %in% names(data))
    stop("'", arg, "' must be a single column name in 'data'.", call. = FALSE)
  if (numeric && !is.numeric(data[[col]]))
    stop("Column '", col, "' supplied to '", arg, "' must be numeric.", call. = FALSE)
  invisible(TRUE)
}

#' @noRd
parse_units <- function(temp_unit = "C", time_unit = "d") {
  temp_unit <- match.arg(temp_unit, c("C", "K"))
  aliases <- c(m = "m", min = "m", h = "h", d = "d", month = "month", y = "y")
  if (length(time_unit) != 1L || !time_unit %in% names(aliases))
    stop("'time_unit' must be one of: m, min, h, d, month, y.", call. = FALSE)
  list(temp = temp_unit, time = unname(aliases[time_unit]))
}

#' @noRd
resolve_direction <- function(direction, slope) {
  direction <- match.arg(direction, c("auto", "decrease", "increase"))
  if (direction == "auto") {
    if (!is.finite(slope) || slope == 0) "conservative" else if (slope < 0) "decrease" else "increase"
  } else direction
}

#' @noRd
validate_limit <- function(limit) {
  if (is.null(limit)) return(NULL)
  if (!is.numeric(limit) || length(limit) != 1L || !is.finite(limit) || limit <= 0)
    stop("'limit' must be NULL or one positive finite number.", call. = FALSE)
  as.numeric(limit)
}

#' @noRd
clean_data <- function(data, cols) {
  ok <- rep(TRUE, nrow(data))
  for (nm in cols) {
    ok <- ok & !is.na(data[[nm]])
    if (is.numeric(data[[nm]])) ok <- ok & is.finite(data[[nm]])
  }
  list(data = data[ok, , drop = FALSE], removed = which(!ok))
}

#' @noRd
summary_stats <- function(data, group, value) {
  key <- interaction(data[group], drop = TRUE, lex.order = TRUE)
  spl <- split(seq_len(nrow(data)), key)
  rows <- lapply(spl, function(i) {
    vals <- data[[value]][i]
    z <- data[i[1L], group, drop = FALSE]
    z$n <- length(vals)
    z$mean <- mean(vals)
    z$sd <- if (length(vals) > 1L) stats::sd(vals) else NA_real_
    z$se <- if (length(vals) > 1L) z$sd / sqrt(length(vals)) else NA_real_
    z
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @noRd
self_predictions <- function(fit, times, conf.level, metric, direction) {
  cf <- stats::coef(fit); vc <- suppressWarnings(stats::vcov(fit)); df <- stats::df.residual(fit)
  b0 <- unname(cf[1L]); b1 <- unname(cf[2L])
  if (metric == "relative" && abs(b0) < .Machine$double.eps^0.5)
    stop("Regression intercept is zero or too close to zero for relative change.", call. = FALSE)
  alpha <- 1 - conf.level; tc <- stats::qt(1 - alpha, df)
  if (metric == "absolute") {
    est_rate <- b1; se_rate <- sqrt(vc[2L, 2L]); scale <- 1
  } else {
    est_rate <- b1 / b0
    var_rate <- vc[2L, 2L] / b0^2 + b1^2 * vc[1L, 1L] / b0^4 -
      2 * b1 * vc[1L, 2L] / b0^3
    se_rate <- sqrt(max(var_rate, 0)); scale <- 100
  }
  lower_rate <- est_rate - tc * se_rate; upper_rate <- est_rate + tc * se_rate
  est <- scale * est_rate * times
  lower <- scale * lower_rate * times; upper <- scale * upper_rate * times
  ci <- if (direction == "decrease") lower else if (direction == "increase") upper else
    ifelse(abs(lower) >= abs(upper), lower, upper)
  data.frame(time = times, estimate = est, ci_limit = ci,
             ci_lower = lower, ci_upper = upper, stringsAsFactors = FALSE)
}

#' @title stability_bias --- Node-wise bias analysis across storage conditions
#' @description Each storage condition is baseline-corrected to its own first
#'   time point. Absolute and relative biases are computed and compared
#'   across conditions at matched time points.
#' @param data A data frame.
#' @param condition Column name for storage condition (e.g. temperature).
#' @param time Column name for time point.
#' @param value Column name for measurement value.
#' @param reference_condition Reference condition for between-condition comparison;
#'   NULL uses the first occurring condition.
#' @param limit A positive number for symmetric allowable bias; NULL to skip.
#' @param bias_type \code{"relative"} (percent) or \code{"absolute"}.
#' @param time_unit Time unit: m, min, h, d, month, or y.
#' @return An S3 object of class \code{"stability_bias"}.
#' @export
#' @examples
#' set.seed(42)
#' df <- data.frame(cond = rep(c("2-8C","25C"), each = 8),
#'                  t = rep(c(0,1,3,7), 4),
#'                  val = c(rnorm(4,100,2), rnorm(4,98,2),
#'                          rnorm(4,100,2), rnorm(4,92,3)))
#' stability_bias(df, "cond", "t", "val", limit = 5)
stability_bias <- function(data, condition, time, value, reference_condition = NULL,
                           limit = NULL, bias_type = c("relative", "absolute"),
                           time_unit = "d") {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  check_column(data, condition, "condition")
  check_column(data, time, "time", numeric = TRUE)
  check_column(data, value, "value", numeric = TRUE)
  time_aliases <- c(m = "m", min = "m", h = "h", d = "d", month = "month", y = "y")
  if (length(time_unit) != 1L || !time_unit %in% names(time_aliases))
    stop("'time_unit' must be one of: m, min, h, d, month, y.", call. = FALSE)
  units <- list(time = unname(time_aliases[time_unit]))
  bias_type <- match.arg(bias_type)
  limit <- validate_limit(limit)

  # Remove rows with missing or infinite values in key columns
  cl <- clean_data(data, c(condition, time, value))
  d <- cl$data
  if (!nrow(d)) stop("No complete observations available.", call. = FALSE)
  d$.condition <- as.character(d[[condition]])
  d$.time <- d[[time]]
  d$.value <- d[[value]]
  cond_order <- unique(d$.condition)
  if (length(cond_order) < 1L) stop("No storage condition found.", call. = FALSE)
  reference_condition <- as.character(if (is.null(reference_condition)) cond_order[1L] else reference_condition)
  if (!reference_condition %in% cond_order)
    stop("'reference_condition' not found in condition column.", call. = FALSE)

  # Compute per-condition, per-time-point summary statistics
  sm <- summary_stats(d, c(".condition", ".time"), ".value")
  names(sm)[names(sm) == "mean"] <- "value"

  # Extract baseline (first time point) for each condition
  bases <- do.call(rbind, lapply(split(sm, sm$.condition), function(z) {
    z <- z[order(z$.time), , drop = FALSE]
    z[1L, c(".condition", ".time", "value", "n", "sd", "se"), drop = FALSE]
  }))
  rownames(bases) <- NULL
  names(bases)[names(bases) == ".time"] <- "baseline_time"
  names(bases)[names(bases) == "value"] <- "baseline"
  if (any(abs(bases$baseline) < .Machine$double.eps^0.5))
    stop("A baseline is zero or too close to zero; relative bias is undefined.", call. = FALSE)

  # Merge baselines back into summaries and compute biases vs baseline
  sm <- merge(sm, bases[c(".condition", "baseline_time", "baseline")], by = ".condition", sort = FALSE)
  sm$absolute_bias <- sm$value - sm$baseline
  sm$relative_bias <- 100 * sm$absolute_bias / sm$baseline
  sm <- sm[order(match(sm$.condition, cond_order), sm$.time), ]
  rownames(sm) <- NULL

  within_bias <- do.call(rbind, lapply(split(sm, sm$.condition), function(z) {
    z <- z[order(z$.time), ]
    if (nrow(z) < 2L) return(NULL)
    cmb <- utils::combn(seq_len(nrow(z)), 2L)
    do.call(rbind, lapply(seq_len(ncol(cmb)), function(j) {
      a <- z[cmb[1L, j], ]; b <- z[cmb[2L, j], ]
      data.frame(condition = a$.condition, time_from = a$.time, time_to = b$.time,
                 value_from = a$value, value_to = b$value,
                 absolute_bias = b$value - a$value,
                 relative_bias = if (abs(a$value) > .Machine$double.eps^0.5)
                   100 * (b$value - a$value) / a$value else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))
  if (is.null(within_bias)) within_bias <- data.frame()

  ref <- sm[sm$.condition == reference_condition, c(".time", "value")]
  names(ref)[2L] <- "reference_value"
  between_bias <- do.call(rbind, lapply(setdiff(cond_order, reference_condition), function(cd) {
    z <- sm[sm$.condition == cd, c(".time", "value")]
    m <- merge(ref, z, by = ".time", all = FALSE)
    if (!nrow(m)) return(NULL)
    data.frame(time = m$.time, reference_condition = reference_condition,
               test_condition = cd, reference_value = m$reference_value,
               test_value = m$value,
               absolute_bias = m$value - m$reference_value,
               relative_bias = ifelse(abs(m$reference_value) > .Machine$double.eps^0.5,
                 100 * (m$value - m$reference_value) / m$reference_value, NA_real_),
               stringsAsFactors = FALSE)
  }))
  if (is.null(between_bias)) between_bias <- data.frame()

  sm$limit_metric <- if (bias_type == "relative") sm$relative_bias else sm$absolute_bias
  sm$pass <- if (is.null(limit)) NA else abs(sm$limit_metric) <= limit
  out <- list(call = match.call(), data = data, summary = sm, baselines = bases,
              within_bias = within_bias, between_bias = between_bias,
              condition_name = condition, time_name = time, value_name = value,
              condition_order = cond_order, reference_condition = reference_condition,
              limit = limit, bias_type = bias_type, units = units,
              removed_rows = cl$removed,
              method = "Condition-specific first-node normalization")
  class(out) <- "stability_bias"
  out
}

#' Print stability-bias results
#'
#' @param x A \code{stability_bias} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.stability_bias <- function(x, ...) {
  num <- function(x, digits = 4L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), as.character(x))
  }
  print_removed <- function(rows) {
    if (length(rows))
      cat(sprintf("  Removed rows:    %d (%s)\n", length(rows), paste(rows, collapse = ", ")))
  }
  cat("\nStability Bias Analysis\n")
  cat(sprintf("  Conditions:      %d (%s)\n", length(x$condition_order), paste(x$condition_order, collapse = ", ")))
  cat(sprintf("  Reference:       %s\n", x$reference_condition))
  cat(sprintf("  Time unit:       %s\n", x$units$time))
  cat(sprintf("  Node comparisons:%d\n", nrow(x$within_bias)))
  if (is.null(x$limit)) cat("  Acceptance:      not evaluated\n") else
    cat(sprintf("  Acceptance:      +/- %s (%s); %d/%d nodes passed\n",
                num(x$limit), x$bias_type, sum(x$summary$pass), nrow(x$summary)))
  print_removed(x$removed_rows)
  cat("\n")

  # Within-condition bias table
  cat("  Within-condition bias (vs first node)\n")
  z <- x$summary
  bias_lab <- if (x$bias_type == "relative") "Bias(%)" else "Bias"
  status_lab <- if (is.null(x$limit)) " " else "Status"
  hdr <- sprintf("  %-8s %8s %10s %10s", "Cond", "Time", "Value", bias_lab)
  if (!is.null(x$limit)) hdr <- paste0(hdr, " ", status_lab)
  cat(hdr, "\n", sep = "")
  cat("  ", paste(rep("-", nchar(hdr) - 2), collapse = ""), "\n", sep = "")
  for (i in seq_len(nrow(z))) {
    line <- sprintf("  %-8s %8s %10s %10s",
                    z[i, ".condition"],
                    num(z[i, ".time"]),
                    num(z[i, "value"]),
                    num(z[i, "relative_bias"]))
    if (!is.null(x$limit))
      line <- paste0(line, sprintf(" %s", ifelse(isTRUE(z[i, "pass"]), "Pass", "Fail")))
    cat(line, "\n", sep = "")
  }
  cat("\n")

  # Between-condition comparison table (informational only)
  if (length(x$condition_order) > 1L && nrow(x$between_bias)) {
    bb <- x$between_bias
    cat(sprintf("  Between-condition bias (reference: %s)\n", x$reference_condition))
    hdr2 <- sprintf("  %8s %-8s %10s %8s", "Time", "Condition", "TestVal", "Bias(%)")
    cat(hdr2, "\n", sep = "")
    cat("  ", paste(rep("-", nchar(hdr2) - 2), collapse = ""), "\n", sep = "")
    for (i in seq_len(nrow(bb))) {
      cat(sprintf("  %8s %-8s %10s %10s\n",
                  num(bb[i, "time"]),
                  bb[i, "test_condition"],
                  num(bb[i, "test_value"]),
                  num(bb[i, "relative_bias"])))
    }
    cat("\n")
  }

  invisible(x)
}

#' Plot stability bias across storage conditions
#'
#' @param x A \code{stability_bias} object.
#' @param ... Reserved for the S3 generic.
#' @return A \code{ggplot} object, invisibly; the plot is also drawn.
#' @examples
#' d <- data.frame(
#'   condition = rep("2-8C", each = 8),
#'   time = rep(c(0, 1, 3, 7), 2),
#'   value = c(100, 99.8, 99.4, 98.9, 100.2, 99.9, 99.5, 98.8)
#' )
#' plot(stability_bias(d, "condition", "time", "value"))
#' @export
plot.stability_bias <- function(x, ...) {
  z <- x$summary
  z$.condition <- factor(z$.condition, levels = x$condition_order)
  p <- ggplot2::ggplot(z, ggplot2::aes(x = .time, y = relative_bias,
                                       color = .condition, group = .condition)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey35", linewidth = 0.45) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5)
  cross <- z[duplicated(z$.time) | duplicated(z$.time, fromLast = TRUE), ]
  if (nrow(cross)) {
    cross <- cross[order(cross$.time, suppressWarnings(as.numeric(as.character(cross$.condition))), cross$.condition), ]
    p <- p + ggplot2::geom_line(data = cross,
      ggplot2::aes(x = .time, y = relative_bias, group = .time),
      inherit.aes = FALSE, color = "grey55", linetype = "dashed", linewidth = 0.55)
  }
  if (!is.null(x$limit)) {
    if (x$bias_type == "relative") {
      p <- p + ggplot2::geom_hline(yintercept = c(-x$limit, x$limit),
                                   linetype = "dotdash", color = "#B2182B")
    } else {
      lim <- x$baselines
      lim$upper <- 100 * x$limit / abs(lim$baseline)
      lim$lower <- -lim$upper
      lim$.condition <- factor(lim$.condition, levels = x$condition_order)
      p <- p +
        ggplot2::geom_hline(data = lim, ggplot2::aes(yintercept = upper, color = .condition), linetype = "dotdash") +
        ggplot2::geom_hline(data = lim, ggplot2::aes(yintercept = lower, color = .condition), linetype = "dotdash")
    }
  }
  p <- p + ggplot2::labs(x = paste0(x$time_name, " (", x$units$time, ")"), y = "Change from first node (%)",
                  color = x$condition_name, title = "Stability bias by storage condition") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom",
                                          plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)
  invisible(p)
}

# stability_regression -- CLSI EP25 stability regression

#' @title stability_regression --- CLSI EP25 regression-based stability assessment
#' @description Performs simple linear regression on a single condition's
#'   measurements (self mode) or on paired biases between test and reference
#'   conditions (compare mode), with one-sided confidence limits to evaluate
#'   whether drift remains within allowable limits.
#' @param data A data frame.
#' @param time Column name for time point.
#' @param value Column name for measurement value.
#' @param condition Column name for storage condition; required for compare mode.
#' @param mode \code{"self"} (single-condition self-reference) or
#'   \code{"compare"} (test vs reference condition).
#' @param test_condition Test condition; auto-selected if NULL.
#' @param reference_condition Reference condition; NULL uses the first condition.
#' @param bias_type \code{"relative"} (percent) or \code{"absolute"}.
#' @param direction \code{"auto"}, \code{"decrease"}, or \code{"increase"}.
#' @param conf.level One-sided confidence level, default 0.95.
#' @param limit A positive number for symmetric allowable bias; NULL to skip.
#' @param time_unit Time unit: m, min, h, d, month, or y.
#' @return An S3 object of class \code{"stability_regression"}.
#' @export
#' @examples
#' set.seed(42)
#' df <- data.frame(cond = rep(c("2-8C","25C"), each = 8),
#'                  t = rep(c(0,1,3,7), 4),
#'                  val = c(rnorm(4,100,2), rnorm(4,98,2),
#'                          rnorm(4,100,2), rnorm(4,92,3)))
#' stability_regression(df, "t", "val", "cond", test_condition = "25C", limit = 5)
#'
#' set.seed(42)
#' df <- data.frame(cond = rep(c("Ref","Test"), each = 8),
#'                  t = rep(c(0,1,3,7), 4),
#'                  val = c(rnorm(4,100,2), rnorm(4,99,2),
#'                          rnorm(4,100,2), rnorm(4,94,3)))
#' stability_regression(df, "t", "val", "cond", mode = "compare",
#'                      test_condition = "Test", reference_condition = "Ref", limit = 5)
stability_regression <- function(data, time, value, condition = NULL,
                                 mode = c("self", "compare"), test_condition = NULL,
                                 reference_condition = NULL,
                                 bias_type = c("relative", "absolute"),
                                 direction = c("auto", "decrease", "increase"),
                                 conf.level = 0.95, limit = NULL, time_unit = "d") {
  # Internal: fit OLS or WLS (weighted by 1/SE^2) regression
  fit_model <- function(dat, x, y, se = NULL) {
    use_wls <- !is.null(se) && length(se) == nrow(dat) && all(is.finite(se) & se > 0)
    fit <- if (use_wls) stats::lm(stats::reformulate(x, y), data = dat, weights = 1 / se^2) else
      stats::lm(stats::reformulate(x, y), data = dat)
    list(fit = fit, method = if (use_wls) "WLS (1/SE^2)" else "OLS",
         fallback = !is.null(se) && !use_wls)
  }

  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  check_column(data, time, "time", numeric = TRUE)
  check_column(data, value, "value", numeric = TRUE)
  mode <- match.arg(mode); bias_type <- match.arg(bias_type); direction <- match.arg(direction)
  units <- parse_units("C", time_unit)
  limit <- validate_limit(limit)
  if (!is.numeric(conf.level) || length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be in (0, 1).", call. = FALSE)
  if (!is.null(condition)) check_column(data, condition, "condition")
  needed <- c(time, value, if (!is.null(condition)) condition)
  cl <- clean_data(data, needed); d <- cl$data
  d$.time <- d[[time]]; d$.value <- d[[value]]
  d$.condition <- if (is.null(condition)) "condition" else as.character(d[[condition]])

  # Self mode: regress measurement vs time on a single condition
  if (mode == "self") {
    choices <- unique(d$.condition)
    test_condition <- as.character(if (is.null(test_condition)) {
      if (length(choices) == 1L) choices else NA_character_
    } else test_condition)
    if (is.na(test_condition) || !test_condition %in% choices)
      stop("With multiple conditions, supply a valid 'test_condition'.", call. = FALSE)
    z <- d[d$.condition == test_condition, ]
    sm <- summary_stats(z, ".time", ".value")
    names(sm)[names(sm) == "mean"] <- "value"
    sm <- sm[order(sm$.time), ]; rownames(sm) <- NULL
    if (nrow(sm) < 4L) stop("At least 4 valid time points are required.", call. = FALSE)
    ff <- fit_model(sm, ".time", "value", sm$se)
    slope <- unname(stats::coef(ff$fit)[2L]); used_direction <- resolve_direction(direction, slope)
    pred <- self_predictions(ff$fit, sm$.time, conf.level, bias_type, used_direction)
    baseline_observed <- sm$value[1L]
    if (bias_type == "relative" && abs(baseline_observed) < .Machine$double.eps^0.5)
      stop("First-node mean is zero or too close to zero.", call. = FALSE)
    sm$observed_bias <- if (bias_type == "relative")
      100 * (sm$value - baseline_observed) / baseline_observed else sm$value - baseline_observed
    nodes <- merge(sm, pred, by.x = ".time", by.y = "time", sort = FALSE)
    nodes <- nodes[order(nodes$.time), ]
    ref_condition <- NULL
  } else {
    # Compare mode: regress paired bias vs time (test vs reference)
    if (is.null(condition)) stop("'condition' is required for compare mode.", call. = FALSE)
    choices <- unique(d$.condition)
    reference_condition <- as.character(if (is.null(reference_condition)) choices[1L] else reference_condition)
    test_condition <- as.character(if (is.null(test_condition)) setdiff(choices, reference_condition)[1L] else test_condition)
    if (!all(c(reference_condition, test_condition) %in% choices) || reference_condition == test_condition)
      stop("Supply distinct valid test and reference conditions.", call. = FALSE)
    sm <- summary_stats(d[d$.condition %in% c(reference_condition, test_condition), ],
                        c(".condition", ".time"), ".value")
    r <- sm[sm$.condition == reference_condition, c(".time", "mean", "se", "n")]
    t <- sm[sm$.condition == test_condition, c(".time", "mean", "se", "n")]
    names(r)[-1L] <- paste0(c("mean", "se", "n"), "_ref")
    names(t)[-1L] <- paste0(c("mean", "se", "n"), "_test")
    nodes <- merge(r, t, by = ".time")
    if (nrow(nodes) < 4L) stop("At least 4 matched time points are required.", call. = FALSE)
    if (bias_type == "relative" && any(abs(nodes$mean_ref) < .Machine$double.eps^0.5))
      stop("A reference mean is zero or too close to zero.", call. = FALSE)
    nodes$observed_bias <- if (bias_type == "relative")
      100 * (nodes$mean_test - nodes$mean_ref) / nodes$mean_ref else nodes$mean_test - nodes$mean_ref
    if (bias_type == "absolute") {
      nodes$bias_se <- sqrt(nodes$se_test^2 + nodes$se_ref^2)
    } else {
      nodes$bias_se <- 100 * sqrt(nodes$se_test^2 / nodes$mean_ref^2 +
        nodes$mean_test^2 * nodes$se_ref^2 / nodes$mean_ref^4)
    }
    ff <- fit_model(nodes, ".time", "observed_bias", nodes$bias_se)
    slope <- unname(stats::coef(ff$fit)[2L]); used_direction <- resolve_direction(direction, slope)
    pr <- stats::predict(ff$fit, newdata = nodes, se.fit = TRUE)
    tc <- stats::qt(conf.level, stats::df.residual(ff$fit))
    nodes$estimate <- as.numeric(pr$fit)
    nodes$ci_lower <- nodes$estimate - tc * as.numeric(pr$se.fit)
    nodes$ci_upper <- nodes$estimate + tc * as.numeric(pr$se.fit)
    nodes$ci_limit <- if (used_direction == "decrease") nodes$ci_lower else if (used_direction == "increase") nodes$ci_upper else
      ifelse(abs(nodes$ci_lower) >= abs(nodes$ci_upper), nodes$ci_lower, nodes$ci_upper)
    ref_condition <- reference_condition
  }

  # Evaluate point estimates and CIS against the allowable limit
  nodes$pass_estimate <- if (is.null(limit)) NA else abs(nodes$estimate) <= limit
  nodes$pass_ci <- if (is.null(limit)) NA else abs(nodes$ci_limit) <= limit

  # Fit diagnostics: R-squared, residual SD, Cook's distance
  fit_sum <- suppressWarnings(summary(ff$fit))
  cooks <- stats::cooks.distance(ff$fit)
  diagnostics <- list(r_squared = fit_sum$r.squared, adjusted_r_squared = fit_sum$adj.r.squared,
                      residual_sd = fit_sum$sigma, df = stats::df.residual(ff$fit),
                      cook_distance = cooks, influential = which(cooks > 1),
                      wls_fallback = ff$fallback)
  out <- list(call = match.call(), mode = mode, fit = ff$fit, fit_method = ff$method,
              nodes = nodes, bias_type = bias_type, direction = used_direction,
              conf.level = conf.level, limit = limit, units = units,
              time_name = time, value_name = value, condition_name = condition,
              test_condition = test_condition, reference_condition = ref_condition,
              diagnostics = diagnostics, removed_rows = cl$removed,
              observed_range = range(nodes$.time),
              method = if (mode == "self") "CLSI EP25 intercept-based drift (Appendix D)" else
                "Reference-anchored paired bias regression")
  class(out) <- "stability_regression"
  out
}

#' Print stability-regression results
#'
#' @param x A \code{stability_regression} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.stability_regression <- function(x, ...) {
  num <- function(x, digits = 4L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), as.character(x))
  }
  print_removed <- function(rows) {
    if (length(rows))
      cat(sprintf("  Removed rows:    %d (%s)\n", length(rows), paste(rows, collapse = ", ")))
  }
  cf <- stats::coef(x$fit)
  cat("\nStability Regression\n")
  cat(sprintf("  Mode:            %s\n", x$mode))
  cat(sprintf("  Test condition:  %s\n", x$test_condition))
  if (!is.null(x$reference_condition)) cat(sprintf("  Reference:       %s\n", x$reference_condition))
  cat(sprintf("  Fit:             %s; intercept = %s, slope = %s\n",
              x$fit_method, num(cf[1L]), num(cf[2L])))
  cat(sprintf("  Direction:       %s\n", x$direction))
  cat(sprintf("  One-sided CI:    %.1f%%\n", 100 * x$conf.level))
  cat(sprintf("  R-squared:       %s\n", num(x$diagnostics$r_squared)))
  if (length(x$diagnostics$influential))
    cat(sprintf("  Cook distance:   row(s) %s > 1; not removed\n", paste(x$diagnostics$influential, collapse = ", ")))
  if (x$diagnostics$wls_fallback) cat("  Note:            WLS unavailable; OLS used.\n")
  print_removed(x$removed_rows)
  cat("\n")

  z <- x$nodes
  hdr <- sprintf("  %8s %12s %10s %10s %6s", "Time", "ObsBias", "Estimate", "CI_limit", "Status")
  cat(hdr, "\n", sep = "")
  cat("  ", paste(rep("-", nchar(hdr) - 2), collapse = ""), "\n", sep = "")
  for (i in seq_len(nrow(z))) {
    ci <- z[i, "ci_limit"]
    star <- if (is.null(x$limit)) "" else {
      if (isTRUE(z[i, "pass_estimate"]) && isTRUE(z[i, "pass_ci"])) "" else "*"
    }
    cat(sprintf("  %8s %12s %10s %10s %6s\n",
                num(z[i, ".time"]),
                num(z[i, "observed_bias"]),
                num(z[i, "estimate"]),
                num(if (is.finite(ci)) ci else NA_real_),
                star))
  }
  cat("\n")
  invisible(x)
}

#' Plot a stability regression
#'
#' @param x A \code{stability_regression} object.
#' @param ... Reserved for the S3 generic.
#' @return A \code{ggplot} object, invisibly; the plot is also drawn.
#' @examples
#' d <- data.frame(
#'   condition = rep("2-8C", each = 8),
#'   time = rep(c(0, 1, 3, 7), 2),
#'   value = c(100, 99.8, 99.4, 98.9, 100.2, 99.9, 99.5, 98.8)
#' )
#' fit <- stability_regression(d, "time", "value", "condition")
#' plot(fit)
#' @export
plot.stability_regression <- function(x, ...) {
  z <- x$nodes
  p <- ggplot2::ggplot(z, ggplot2::aes(x = .time)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey40", linewidth = 0.45) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmin(estimate, ci_limit), ymax = pmax(estimate, ci_limit)),
                         fill = "#9ECAE1", alpha = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = estimate), color = "#2166AC", linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = observed_bias), size = 2.3)
  if (!is.null(x$limit))
    p <- p + ggplot2::geom_hline(yintercept = c(-x$limit, x$limit), color = "#B2182B", linetype = "dashed")
  p <- p + ggplot2::labs(x = paste0(x$time_name, " (", x$units$time, ")"),
                          y = if (x$bias_type == "relative") "Bias (%)" else "Absolute bias",
                          title = "Stability regression", subtitle = paste(x$fit_method, x$method, sep = " | ")) +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)
  invisible(p)
}

#' @title stability_time --- Predict time-to-limit from a stability regression
#' @description Given a fitted regression model, computes the time at which
#'   the point estimate or one-sided confidence limit first reaches the
#'   allowable bias. May extrapolate beyond the observed time range.
#' @param object A \code{stability_regression} object.
#' @param limit Symmetric allowable bias.
#' @param max_time Maximum search time; NULL defaults to 100x the maximum observed time.
#' @return An S3 object of class \code{"stability_time"}.
#' @export
#' @examples
#' set.seed(42)
#' df <- data.frame(cond = rep(c("2-8C","25C"), each = 8),
#'                  t = rep(c(0,1,3,7), 4),
#'                  val = c(rnorm(4,100,2), rnorm(4,98,2),
#'                          rnorm(4,100,2), rnorm(4,92,3)))
#' r <- stability_regression(df, "t", "val", "cond", test_condition = "25C", limit = 5)
#' stability_time(r)
#'
#' set.seed(42)
#' df <- data.frame(cond = rep(c("Ref","Test"), each = 8),
#'                  t = rep(c(0,1,3,7), 4),
#'                  val = c(rnorm(4,100,2), rnorm(4,99,2),
#'                          rnorm(4,100,2), rnorm(4,94,3)))
#' r2 <- stability_regression(df, "t", "val", "cond", mode = "compare",
#'                            test_condition = "Test", reference_condition = "Ref", limit = 5)
#' stability_time(r2, limit = 5)
stability_time <- function(object, limit = NULL, max_time = NULL) {
  # Internal: predict bias at a given time for compare mode (point estimate + one-sided CI)
  compare_at_time <- function(object, t) {
    nd <- data.frame(.time = t)
    pr <- stats::predict(object$fit, newdata = nd, se.fit = TRUE)
    tc <- stats::qt(object$conf.level, stats::df.residual(object$fit))
    est <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
    ci <- if (object$direction == "decrease") est - tc * se else if (object$direction == "increase") est + tc * se else
      if (abs(est - tc * se) >= abs(est + tc * se)) est - tc * se else est + tc * se
    c(estimate = est, ci = ci)
  }

  # Internal: find earliest time where |fun(t)| crosses target, using grid search + uniroot
  find_root <- function(fun, target, start, max_time) {
    grid <- seq(start, max_time, length.out = 501L)
    vals <- vapply(grid, function(t) abs(fun(t)) - target, numeric(1L))
    hit <- which(vals >= 0)[1L]
    if (is.na(hit)) return(Inf)
    if (hit == 1L) return(grid[1L])
    stats::uniroot(function(t) abs(fun(t)) - target, c(grid[hit - 1L], grid[hit]))$root
  }

  # -- Validate inputs --
  if (!inherits(object, "stability_regression")) stop("'object' must be a stability_regression object.", call. = FALSE)
  limit <- if (is.null(limit)) object$limit else validate_limit(limit)
  if (is.null(limit)) stop("Supply 'limit' or fit stability_regression with a limit.", call. = FALSE)
  obs_max <- max(object$observed_range)
  max_time <- if (is.null(max_time)) max(obs_max * 100, obs_max + 100) else max_time
  if (!is.numeric(max_time) || length(max_time) != 1L || max_time <= 0)
    stop("'max_time' must be one positive number.", call. = FALSE)

  # Choose the appropriate prediction function based on regression mode
  if (object$mode == "self") {
    f <- function(t) self_predictions(object$fit, t, object$conf.level,
                                        object$bias_type, object$direction)
    est_fun <- function(t) f(t)$estimate
    ci_fun <- function(t) f(t)$ci_limit
  } else {
    est_fun <- function(t) compare_at_time(object, t)["estimate"]
    ci_fun <- function(t) compare_at_time(object, t)["ci"]
  }

  # Solve for the limit-crossing time of the point estimate and the one-sided CI
  t_est <- find_root(est_fun, limit, 0, max_time)
  t_ci <- find_root(ci_fun, limit, 0, max_time)
  ans <- data.frame(target = c("estimate", "one_sided_ci"), time = c(t_est, t_ci),
                    extrapolated = c(t_est, t_ci) > obs_max,
                    within_search = is.finite(c(t_est, t_ci)), stringsAsFactors = FALSE)
  out <- list(call = match.call(), result = ans, limit = limit, bias_type = object$bias_type,
              direction = object$direction, time_unit = object$units$time,
              observed_max = obs_max, max_time = max_time, source = object)
  class(out) <- "stability_time"
  out
}

#' Print stability limit-time predictions
#'
#' @param x A \code{stability_time} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.stability_time <- function(x, ...) {
  num <- function(x, digits = 4L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), as.character(x))
  }
  cat("\nStability Limit Time\n")
  cat(sprintf("  Limit: %s (%s); direction: %s\n", num(x$limit), x$bias_type, x$direction))
  for (i in seq_len(nrow(x$result))) {
    r <- x$result[i, ]
    val <- if (is.finite(r$time)) paste0(num(r$time), " ", x$time_unit) else "not reached"
    note <- if (isTRUE(r$extrapolated)) " [EXTRAPOLATED]" else ""
    cat(sprintf("  %-14s %s%s\n", r$target, val, note))
  }
  cat("\n")
  invisible(x)
}

#' Plot stability limit-time predictions
#'
#' @param x A \code{stability_time} object.
#' @param ... Reserved for the S3 generic.
#' @return A \code{ggplot} object, invisibly; the plot is also drawn.
#' @examples
#' d <- data.frame(
#'   condition = rep("2-8C", each = 8),
#'   time = rep(c(0, 1, 3, 7), 2),
#'   value = c(100, 99.8, 99.4, 98.9, 100.2, 99.9, 99.5, 98.8)
#' )
#' fit <- stability_regression(d, "time", "value", "condition")
#' limit_time <- stability_time(fit, limit = 5)
#' plot(limit_time)
#' @export
plot.stability_time <- function(x, ...) {
  # Reuse the regression plot and overlay limit-crossing markers
  src <- x$source
  z <- src$nodes

  # Determine y-axis label based on bias type
  ylab <- if (src$bias_type == "relative") "Bias (%)" else "Absolute bias"

  # Build the regression plot (shared aesthetics with plot.stability_regression)
  p <- ggplot2::ggplot(z, ggplot2::aes(x = .time)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey40", linewidth = 0.45) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmin(estimate, ci_limit), ymax = pmax(estimate, ci_limit)),
                         fill = "#9ECAE1", alpha = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = estimate), color = "#2166AC", linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = observed_bias), size = 2.3) +
    ggplot2::geom_hline(yintercept = c(-x$limit, x$limit), color = "#B2182B", linetype = "dashed") +
    ggplot2::labs(x = paste0(src$time_name, " (", src$units$time, ")"),
                  y = ylab,
                  title = "Stability limit time predictions",
                  subtitle = paste0("Limit: ", x$limit, " (", x$bias_type, "), direction: ", x$direction)) +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # Add vertical markers at the crossing times where |bias| reaches the limit
  colors <- c("estimate" = "#2166AC", "one_sided_ci" = "#E08214")
  labels <- c("estimate" = "Est.", "one_sided_ci" = "CI")
  for (i in seq_len(nrow(x$result))) {
    r <- x$result[i, ]
    if (!is.finite(r$time)) next
    time_str <- formatC(r$time, format = "f", digits = 2)
    ext <- if (isTRUE(r$extrapolated)) " [EXT]" else ""
    label <- paste0(labels[r$target], " t=", time_str, ext)
    p <- p + ggplot2::geom_vline(xintercept = r$time, color = colors[r$target], linewidth = 0.7, linetype = "dashed") +
      ggplot2::annotate("text", x = r$time, y = max(abs(z$ci_limit), abs(z$estimate), abs(x$limit), na.rm = TRUE) * 0.95,
                         label = label, color = colors[r$target], hjust = 1.05, size = 3.2)
  }

  print(p)
  invisible(p)
}

#' @title mkt --- Mean Kinetic Temperature
#' @description Calculates the mean kinetic temperature from one or more
#'   temperature probes. Supports equal-interval and time-weighted
#'   trapezoidal methods.
#' @param data A data frame.
#' @param temp_cols One or more column names for temperature probes.
#' @param time Optional time column; accepts numeric, Date, or POSIXct.
#' @param temp_unit Temperature unit: \code{"C"} or \code{"K"}.
#' @param ea Activation energy in kJ/mol; default 83.144.
#' @return An S3 object of class \code{"mkt"}.
#' @export
#' @examples
#' temp_df <- data.frame(t1 = c(25,26,27), t2 = c(24,25,26), time = 1:3)
#' mkt(temp_df, c("t1","t2"), "time")
mkt <- function(data, temp_cols, time = NULL, temp_unit = "C", ea = 83.144) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(temp_cols) || !length(temp_cols) || any(!temp_cols %in% names(data)))
    stop("'temp_cols' must contain valid column names.", call. = FALSE)
  if (any(!vapply(data[temp_cols], is.numeric, logical(1L))))
    stop("All temperature columns must be numeric.", call. = FALSE)
  temp_unit <- match.arg(temp_unit, c("C", "K"))
  if (!is.numeric(ea) || length(ea) != 1L || !is.finite(ea) || ea <= 0)
    stop("'ea' must be one positive number in kJ/mol.", call. = FALSE)
  if (!is.null(time)) check_column(data, time, "time")
  R <- 8.3145e-3
  default_ea <- missing(ea)
  tt <- NULL

  # Parse time column if provided (supports numeric, Date, POSIXct)
  if (!is.null(time)) {
    raw <- data[[time]]
    tt <- if (inherits(raw, "POSIXt")) as.numeric(raw) else if (inherits(raw, "Date")) as.numeric(raw) else raw
    if (!is.numeric(tt) || any(!is.finite(tt))) stop("Time values must be finite numeric/date-time values.", call. = FALSE)
    if (is.unsorted(tt, strictly = TRUE)) stop("Time values must be strictly increasing and unique.", call. = FALSE)
  }
  # Internal: compute MKT for a single probe using equal-interval or trapezoidal method
  calc <- function(v) {
    kelvin <- if (temp_unit == "C") v + 273.15 else v
    ok <- is.finite(kelvin) & kelvin > 0
    if (sum(ok) < 2L) return(list(mkt = NA_real_, numerator = NA_real_, weight = 0, n = sum(ok)))
    q <- exp(-ea / (R * kelvin))
    if (is.null(tt)) {
      # Equal-interval MKT: simple mean of exp(-Ea/(RT))
      val <- mean(q[ok]); wt <- sum(ok)
    } else {
      # Time-weighted trapezoidal MKT: integrate over time intervals
      interval_ok <- ok[-length(ok)] & ok[-1L]
      dt <- diff(tt)
      val_num <- sum(dt[interval_ok] * (q[-length(q)][interval_ok] + q[-1L][interval_ok]) / 2)
      wt <- sum(dt[interval_ok])
      if (wt <= 0) return(list(mkt = NA_real_, numerator = NA_real_, weight = 0, n = sum(ok)))
      val <- val_num / wt
    }
    # Back-transform the averaged Boltzmann factor to temperature
    list(mkt = -ea / (R * log(val)), numerator = val * wt, weight = wt, n = sum(ok))
  }

  # Compute MKT for each probe separately
  res <- lapply(data[temp_cols], calc)
  tab <- data.frame(probe = temp_cols,
                    mkt_K = vapply(res, `[[`, numeric(1L), "mkt"),
                    n = vapply(res, `[[`, numeric(1L), "n"), stringsAsFactors = FALSE)
  tab$mkt_C <- tab$mkt_K - 273.15

  # Pool valid probes into a combined MKT (weighted by time or count)
  valid <- vapply(res, function(z) is.finite(z$numerator) && z$weight > 0, logical(1L))
  if (!any(valid)) stop("No probe has enough valid temperature exposure data.", call. = FALSE)
  pooled_q <- sum(vapply(res[valid], `[[`, numeric(1L), "numerator")) /
    sum(vapply(res[valid], `[[`, numeric(1L), "weight"))
  pooled_K <- -ea / (R * log(pooled_q))
  combined <- data.frame(probe = "Combined", mkt_K = pooled_K, n = sum(tab$n[valid]), mkt_C = pooled_K - 273.15)
  out <- list(call = match.call(), probes = tab, combined = combined, ea = ea,
              default_ea = default_ea, temp_unit = temp_unit, time_name = time,
              method = if (is.null(time)) "Equal-interval MKT" else "Time-weighted trapezoidal MKT",
              reference = "CLSI EP25 (2023), equations 1-2")
  class(out) <- "mkt"
  out
}

#' Print mean kinetic temperature results
#'
#' @param x A \code{mkt} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.mkt <- function(x, ...) {
  num <- function(x, digits = 4L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), as.character(x))
  }
  cat("\nMean Kinetic Temperature\n")
  cat(sprintf("  Ea: %s kJ/mol%s\n", num(x$ea, 3L), if (x$default_ea) " (default assumption)" else ""))
  cat(sprintf("  Method: %s\n\n", x$method))
  all <- rbind(x$probes, x$combined)
  for (i in seq_len(nrow(all)))
    cat(sprintf("  %-18s %9s K  %9s C\n", all$probe[i], num(all$mkt_K[i], 3L), num(all$mkt_C[i], 3L)))
  cat("\n")
  invisible(x)
}

#' @title arrhenius --- Arrhenius accelerated stability analysis
#' @description Estimates degradation rates at multiple elevated temperatures,
#'   fits the Arrhenius equation, and extrapolates to the target storage
#'   temperature to predict shelf life.
#' @param data A data frame.
#' @param temperature Column name for storage temperature (numeric).
#' @param time Column name for time point (numeric).
#' @param value Column name for measurement value (numeric).
#' @param target_temp Target storage temperature.
#' @param order \code{"auto"}, \code{"zero"}, or \code{"first"}.
#' @param direction \code{"auto"}, \code{"decrease"}, or \code{"increase"}.
#' @param limit A positive number for symmetric allowable bias.
#' @param bias_type \code{"relative"} (percent) or \code{"absolute"}.
#' @param temp_unit Temperature unit: \code{"C"} or \code{"K"}.
#' @param time_unit Time unit: m, min, h, d, month, or y.
#' @param conf.level Confidence level, default 0.95.
#' @return An S3 object of class \code{"arrhenius"}.
#' @export
#' @examples
#' temp_df <- data.frame(temp = rep(c(25,30,37), each = 3),
#'                       time = rep(c(0,1,3), 3),
#'                       val = c(100,98,95, 100,96,90, 100,94,85))
#' arrhenius(temp_df, "temp", "time", "val", target_temp = 5, limit = 10)
arrhenius <- function(data, temperature, time, value, target_temp,
                      order = c("auto", "zero", "first"),
                      direction = c("auto", "decrease", "increase"),
                      limit, bias_type = c("relative", "absolute"),
                      temp_unit = "C", time_unit = "d", conf.level = 0.95) {
  # Internal: fit zero-order or first-order kinetics at each temperature
  fit_kinetic_order <- function(sm, order, direction) {
    spl <- split(sm, sm$.temperature)
    fits <- lapply(spl, function(z) {
      z <- z[order(z$.time), ]
      if (nrow(z) < 3L) return(NULL)
      if (order == "first" && any(z$value <= 0)) return(NULL)
      yy <- if (order == "first") log(z$value) else z$value
      model_data <- data.frame(response = yy, time = z$.time)
      fit <- stats::lm(response ~ time, data = model_data)
      slope <- unname(stats::coef(fit)[2L])
      rate <- if (direction == "decrease") -slope else slope
      pred_raw <- if (order == "first") exp(stats::fitted(fit)) else stats::fitted(fit)
      base <- abs(z$value[1L]); nrm <- if (base > .Machine$double.eps^0.5) base else 1
      data.frame(temperature = z$.temperature[1L], rate = rate, slope = slope,
                 intercept = unname(stats::coef(fit)[1L]), p_value = stats::coef(summary(fit))[2L, 4L],
                 rmse_norm = sqrt(mean((z$value - pred_raw)^2)) / nrm,
                 n = nrow(z), stringsAsFactors = FALSE) -> row
      list(fit = fit, data = z, row = row, predicted = pred_raw)
    })
    fits <- fits[!vapply(fits, is.null, logical(1L))]
    rates <- if (length(fits)) do.call(rbind, lapply(fits, `[[`, "row")) else data.frame()
    list(order = order, fits = fits, rates = rates,
         score = if (nrow(rates)) sqrt(mean(rates$rmse_norm^2)) else Inf)
  }
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  check_column(data, temperature, "temperature", numeric = TRUE)
  check_column(data, time, "time", numeric = TRUE)
  check_column(data, value, "value", numeric = TRUE)
  order <- match.arg(order); direction <- match.arg(direction); bias_type <- match.arg(bias_type)
  units <- parse_units(temp_unit, time_unit); limit <- validate_limit(limit)
  if (!is.numeric(target_temp) || length(target_temp) != 1L || !is.finite(target_temp))
    stop("'target_temp' must be one finite number.", call. = FALSE)
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be in (0, 1).", call. = FALSE)

  # Clean and prepare data
  cl <- clean_data(data, c(temperature, time, value)); d <- cl$data
  d$.temperature <- d[[temperature]]; d$.time <- d[[time]]; d$.value <- d[[value]]
  temp_K <- if (temp_unit == "C") d$.temperature + 273.15 else d$.temperature
  target_K <- if (temp_unit == "C") target_temp + 273.15 else target_temp
  if (any(temp_K <= 0) || target_K <= 0) stop("All Kelvin temperatures must be positive.", call. = FALSE)
  sm <- summary_stats(d, c(".temperature", ".time"), ".value")
  names(sm)[names(sm) == "mean"] <- "value"
  if (length(unique(sm$.temperature)) < 3L) stop("At least 3 elevated temperatures are required.", call. = FALSE)

  # Infer degradation direction from the median slope across temperatures
  raw_slopes <- vapply(split(sm, sm$.temperature), function(z) stats::coef(stats::lm(value ~ .time, z))[2L], numeric(1L))
  used_direction <- resolve_direction(direction, stats::median(raw_slopes, na.rm = TRUE))
  if (used_direction == "conservative") stop("Cannot infer a degradation direction; supply 'direction'.", call. = FALSE)

  # Fit both kinetic orders and score each by normalized RMSE
  zero <- fit_kinetic_order(sm, "zero", used_direction)
  first <- fit_kinetic_order(sm, "first", used_direction)
  candidates <- list(zero = zero, first = first)
  if (order == "auto") {
    scores <- vapply(candidates, `[[`, numeric(1L), "score")
    selected_name <- if (isTRUE(all.equal(scores[1L], scores[2L]))) "first" else names(which.min(scores))
  } else selected_name <- order
  selected <- candidates[[selected_name]]
  rates <- selected$rates
  rates$valid_rate <- is.finite(rates$rate) & rates$rate > 0
  valid_rates <- rates[rates$valid_rate, ]
  if (nrow(valid_rates) < 3L)
    stop("Fewer than 3 temperatures have a positive degradation rate.", call. = FALSE)
  valid_rates$temp_K <- if (temp_unit == "C") valid_rates$temperature + 273.15 else valid_rates$temperature
  valid_rates$inv_temp <- 1 / valid_rates$temp_K

  # Arrhenius regression: ln(rate) ~ 1/T  (Eq. 7-13 in CLSI EP25)
  arr_fit <- stats::lm(log(rate) ~ inv_temp, data = valid_rates)
  nd <- data.frame(inv_temp = 1 / target_K)
  pr <- stats::predict(arr_fit, nd, interval = "confidence", level = conf.level)
  rate_target <- exp(pr[1L, "fit"]); rate_lo <- exp(pr[1L, "lwr"]); rate_hi <- exp(pr[1L, "upr"])
  base_rows <- do.call(rbind, lapply(split(sm, sm$.temperature), function(z) z[which.min(z$.time), ]))
  baseline <- mean(base_rows$value)

  # Internal: compute duration to reach the allowable limit at a given degradation rate k
  # Zero-order: t = change/k; First-order: t = -ln(1 - fraction)/k or ln(1 + fraction)/k
  duration_from_rate <- function(k) {
    if (selected_name == "zero") {
      change <- if (bias_type == "relative") abs(baseline) * limit / 100 else limit
      change / k
    } else if (bias_type == "relative") {
      frac <- limit / 100
      if (used_direction == "decrease") {
        if (frac >= 1) stop("Relative decrease limit must be < 100% for first-order kinetics.", call. = FALSE)
        -log(1 - frac) / k
      } else log(1 + frac) / k
    } else {
      if (used_direction == "decrease") {
        if (limit >= baseline) stop("Absolute decrease limit must be below the pooled baseline.", call. = FALSE)
        -log((baseline - limit) / baseline) / k
      } else log((baseline + limit) / baseline) / k
    }
  }
  duration <- duration_from_rate(rate_target)
  duration_ci <- sort(c(duration_from_rate(rate_hi), duration_from_rate(rate_lo)))
  R <- 8.3145e-3
  activation_energy <- -stats::coef(arr_fit)[2L] * R
  degradation <- vapply(split(sm, sm$.temperature), function(z) {
    z <- z[order(z$.time), ]; 100 * abs(tail(z$value, 1L) - z$value[1L]) / abs(z$value[1L])
  }, numeric(1L))

  # Collect warnings about data quality
  warnings <- character()
  if (any(valid_rates$p_value > 0.05)) warnings <- c(warnings, "One or more temperature-specific slopes are not significant.")
  if (any(degradation > 20)) warnings <- c(warnings, "Observed degradation exceeds 20% at one or more temperatures.")
  if (any(!rates$valid_rate)) warnings <- c(warnings, "One or more temperatures had a nonpositive degradation rate and were excluded.")
  out <- list(call = match.call(), selected_order = selected_name, requested_order = order,
              direction = used_direction, candidates = candidates, rates = rates,
              arrhenius_fit = arr_fit, target_temp = target_temp, target_K = target_K,
              target_rate = rate_target, target_rate_ci = c(rate_lo, rate_hi),
              activation_energy = activation_energy, baseline = baseline,
              limit = limit, bias_type = bias_type, duration = duration,
              duration_ci = duration_ci, conf.level = conf.level, units = units,
              warnings = warnings, removed_rows = cl$removed,
              observed = sm, temperature_name = temperature, time_name = time, value_name = value,
              method = "CLSI EP25 Arrhenius method (equations 7-13)")
  class(out) <- "arrhenius"
  out
}

#' Print Arrhenius stability-analysis results
#'
#' @param x An \code{arrhenius} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.arrhenius <- function(x, ...) {
  num <- function(x, digits = 4L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), as.character(x))
  }
  print_removed <- function(rows) {
    if (length(rows))
      cat(sprintf("  Removed rows:    %d (%s)\n", length(rows), paste(rows, collapse = ", ")))
  }
  cat("\nArrhenius Stability Analysis\n")
  cat(sprintf("  Kinetic order:   %s%s\n", x$selected_order, if (x$requested_order == "auto") " (auto-selected)" else ""))
  cat(sprintf("  Direction:       %s\n", x$direction))
  cat(sprintf("  Activation Ea:   %s kJ/mol\n", num(x$activation_energy, 3L)))
  cat(sprintf("  Target rate:     %s /%s\n", num(x$target_rate, 6L), x$units$time))
  cat(sprintf("  Limit time:      %s %s [%s, %s]\n", num(x$duration), x$units$time,
              num(x$duration_ci[1L]), num(x$duration_ci[2L])))
  if (length(x$warnings)) for (w in x$warnings) cat("  WARNING:         ", w, "\n", sep = "")
  print_removed(x$removed_rows)
  cat("\n")
  invisible(x)
}

#' Plot Arrhenius stability analysis
#'
#' @param x An \code{arrhenius} object.
#' @param ... Reserved for the S3 generic.
#' @return A named list of \code{ggplot} objects, invisibly. The kinetic-fit
#'   and Arrhenius-regression plots are also drawn.
#' @examples
#' d <- data.frame(
#'   temperature = rep(c(25, 30, 37), each = 3),
#'   time = rep(c(0, 1, 3), 3),
#'   value = c(100, 98, 95, 100, 96, 90, 100, 94, 85)
#' )
#' fit <- arrhenius(d, "temperature", "time", "value",
#'                  target_temp = 5, limit = 10)
#' plot(fit)
#' @export
plot.arrhenius <- function(x, ...) {
  sel <- x$candidates[[x$selected_order]]
  kin <- do.call(rbind, lapply(sel$fits, function(z) {
    data.frame(temperature = z$data$.temperature, time = z$data$.time,
               observed = z$data$value, fitted = z$predicted)
  }))
  p1 <- ggplot2::ggplot(kin, ggplot2::aes(x = time)) +
    ggplot2::geom_point(ggplot2::aes(y = observed)) +
    ggplot2::geom_line(ggplot2::aes(y = fitted), color = "#2166AC") +
    ggplot2::facet_wrap(~ temperature, scales = "free_y") +
    ggplot2::labs(x = paste0(x$time_name, " (", x$units$time, ")"), y = x$value_name,
                  title = paste("Kinetic fits:", x$selected_order)) + ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  rr <- x$rates[x$rates$valid_rate, ]
  rr$temp_K <- if (x$units$temp == "C") rr$temperature + 273.15 else rr$temperature
  rr$inv_temp1000 <- 1000 / rr$temp_K
  rr$log_rate <- log(rr$rate)
  grid <- data.frame(inv_temp = seq(min(1 / rr$temp_K), max(1 / rr$temp_K), length.out = 100L))
  grid$fit <- stats::predict(x$arrhenius_fit, grid)
  grid$inv_temp1000 <- 1000 * grid$inv_temp
  p2 <- ggplot2::ggplot(rr, ggplot2::aes(x = inv_temp1000, y = log_rate)) +
    ggplot2::geom_point(size = 2.5) + ggplot2::geom_line(data = grid, ggplot2::aes(y = fit), color = "#B2182B") +
    ggplot2::labs(x = "1000 / Temperature (K)", y = "ln(rate)", title = "Arrhenius regression") +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p1); print(p2)
  invisible(list(kinetics = p1, arrhenius = p2))
}

#' @title stability_plan --- CLSI EP25 Appendix A time-point planning
#' @description Based on the ratio of expected drift to allowable drift and
#'   the reproducibility variability, looks up the minimum number of time
#'   points required per regression series in the EP25 Appendix A tables.
#' @param allowable_drift Allowable drift (positive number).
#' @param variability Repeatability CV/SD (simple mode) or a data frame of
#'   variance components.
#' @param replicates Number of replicate measurements per time point.
#' @param expected_drift Expected drift; NULL defaults to half of allowable_drift.
#' @param power Statistical power: 0.8 or 0.9.
#' @param mode \code{"simple"} or \code{"components"}.
#' @param bias_type \code{"relative"} or \code{"absolute"}.
#' @return An S3 object of class \code{"stability_plan"}.
#' @export
#' @examples
#' stability_plan(4, 1.5, replicates = c(1,2,3))
#'
#' comp <- data.frame(source = c("repeatability","lot","operator"),
#'                    variability = c(1.5,1.0,0.8),
#'                    levels = c(1,3,2))
#' stability_plan(4, comp, replicates = c(1,2,3), mode = "components")
stability_plan <- function(allowable_drift, variability, replicates = 1L,
                           expected_drift = NULL, power = c(0.8, 0.9),
                           mode = c("simple", "components"),
                           bias_type = c("relative", "absolute")) {
  power <- match.arg(as.character(power), c("0.8", "0.9")); power <- as.numeric(power)
  mode <- match.arg(mode); bias_type <- match.arg(bias_type)

  # Internal: EP25 Appendix A lookup tables (Table A1 for 80% power, A2 for 90% power)
  # Rows: precision ratio (allowable / residual); Cols: drift ratio (expected / allowable)
  plan_tables <- function(power) {
    cols <- sprintf("%.1f", seq(0, 0.8, 0.1)); rows <- as.character(10:2)
    a80 <- rbind(
      c(4,4,4,4,4,5,6,9,20), c(4,4,4,4,5,5,7,11,25),
      c(4,4,4,4,5,6,8,14,32), c(4,4,4,5,6,7,11,18,NA),
      c(4,4,5,6,7,9,14,25,NA), c(5,5,6,7,9,13,20,35,NA),
      c(6,7,8,10,14,20,32,NA,NA), c(9,11,14,18,25,36,NA,NA,NA),
      c(20,24,30,NA,NA,NA,NA,NA,NA))
    a90 <- rbind(
      c(4,4,4,4,5,6,8,12,28), c(4,4,4,5,5,6,9,16,35),
      c(4,4,5,5,6,8,11,20,NA), c(4,5,5,6,7,10,15,25,NA),
      c(5,5,6,7,9,13,20,35,NA), c(6,6,8,9,13,18,28,NA,NA),
      c(7,9,11,14,19,28,NA,NA,NA), c(12,15,19,25,34,NA,NA,NA,NA),
      c(28,34,NA,NA,NA,NA,NA,NA,NA))
    tab <- if (power == 0.8) a80 else a90
    dimnames(tab) <- list(rows, cols)
    tab
  }
  if (!is.numeric(allowable_drift) || length(allowable_drift) != 1L || allowable_drift <= 0)
    stop("'allowable_drift' must be one positive number.", call. = FALSE)
  assumed_expected <- is.null(expected_drift)
  expected_drift <- if (is.null(expected_drift)) (0.5 * allowable_drift) else expected_drift
  if (!is.numeric(expected_drift) || length(expected_drift) != 1L || expected_drift < 0)
    stop("'expected_drift' must be one nonnegative number.", call. = FALSE)
  if (!is.numeric(replicates) || any(!is.finite(replicates) | replicates < 1 | replicates != floor(replicates)))
    stop("'replicates' must contain positive integers.", call. = FALSE)
  replicates <- unique(as.integer(replicates))
  if (expected_drift / allowable_drift > 0.8)
    warning("Expected/allowable drift ratio exceeds the EP25 table range (0.8).", call. = FALSE)

  # Internal: compute residual variability for a given number of replicates
  residual_for <- function(n) {
    if (mode == "simple") {
      if (!is.numeric(variability) || length(variability) != 1L || variability <= 0)
        stop("Simple mode requires one positive 'variability'.", call. = FALSE)
      variability / sqrt(n)
    } else {
      req <- c("source", "variability", "levels")
      if (!is.data.frame(variability) || any(!req %in% names(variability)))
        stop("Components mode requires columns: source, variability, levels.", call. = FALSE)
      v <- variability
      if (any(!is.finite(v$variability) | v$variability < 0) || any(!is.finite(v$levels) | v$levels <= 0))
        stop("Component variability and levels must be valid nonnegative/positive numbers.", call. = FALSE)
      idx <- tolower(as.character(v$source)) == "repeatability"
      if (sum(idx) != 1L) stop("Components must contain exactly one 'repeatability' row.", call. = FALSE)
      v$levels[idx] <- n
      sqrt(sum(v$variability^2 / v$levels))
    }
  }
  tab <- plan_tables(power)

  # Map drift ratio to table column, then look up minimum time points per replicate scenario
  drift_ratio <- expected_drift / allowable_drift
  drift_col <- ceiling((drift_ratio - 1e-12) * 10) / 10
  rows <- lapply(seq_along(replicates), function(i) {
    res <- residual_for(replicates[i])
    # EP25 Appendix A rounds the residual to 2 significant figures for table lookup.
    # The unrounded value is retained for audit trail.
    planning_res <- signif(res, 2L)
    precision_ratio <- allowable_drift / planning_res
    table_row <- if (precision_ratio >= 10) 10 else floor(precision_ratio + 1e-12)
    supported <- table_row >= 2 && drift_col >= 0 && drift_col <= 0.8
    nodes <- NA_integer_; reason <- ""
    if (supported) {
      nodes <- tab[as.character(table_row), sprintf("%.1f", drift_col)]
      if (is.na(nodes)) { supported <- FALSE; reason <- "No feasible entry in EP25 table." }
    } else reason <- "Ratios are outside the EP25 table range."
    data.frame(replicates = replicates[i],
               residual_variability = res, residual_for_lookup = planning_res,
               allowable_residual_ratio = precision_ratio,
               table_precision_row = if (table_row >= 2) table_row else NA_real_,
               expected_allowable_ratio = drift_ratio, table_drift_column = drift_col,
               min_time_points = nodes, supported = supported, reason = reason,
               stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  out <- list(call = match.call(), result = result, allowable_drift = allowable_drift,
              expected_drift = expected_drift, assumed_expected = assumed_expected,
              power = power, mode = mode, bias_type = bias_type,
              method = paste0("CLSI EP25 Appendix A Table A", if (power == 0.8) "1" else "2"))
  class(out) <- "stability_plan"
  out
}

#' Print a stability study time-point plan
#'
#' @param x A \code{stability_plan} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.stability_plan <- function(x, ...) {
  num <- function(x, digits = 4L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), as.character(x))
  }
  cat("\nStability Study Time-Point Plan\n")
  cat(sprintf("  Allowable drift: %s; expected drift: %s%s\n", num(x$allowable_drift),
              num(x$expected_drift), if (x$assumed_expected) " (default assumption)" else ""))
  cat(sprintf("  Success target:  %.0f%%; method: %s\n\n", 100 * x$power, x$method))
  cat(sprintf("  %10s %12s %12s %10s\n", "Replicates", "Residual", "Allow/Resid", "Min nodes"))
  cat("  ", paste(rep("-", 48), collapse = ""), "\n", sep = "")
  notes <- character()
  for (i in seq_len(nrow(x$result))) {
    z <- x$result[i, ]
    min_show <- if (z$supported) as.character(z$min_time_points) else "n/a"
    cat(sprintf("  %10d %12s %12s %10s\n", z$replicates, num(z$residual_variability),
                num(z$allowable_residual_ratio), min_show))
    if (!z$supported) notes <- c(notes, z$reason)
  }
  for (n in unique(notes)) cat("\n  n/a -- ", n, sep = "")
  cat("\n\n")
  invisible(x)
}
