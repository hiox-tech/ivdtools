#' Precision variance component analysis -- Sadler precision profile fitting based on VCA and VFP
#'
#' @importFrom stats aggregate sd qnorm dnorm quantile qchisq
#' @importFrom utils capture.output
#' @importFrom VCA anovaVCA vcovVC VCAinference
#' @importFrom VFP fit_vfp get_model predict.VFP
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline
#'   geom_histogram geom_density geom_ribbon geom_col geom_text
#'   facet_wrap facet_grid labs theme_bw scale_fill_brewer
#'   scale_x_continuous expansion after_stat theme element_text margin
#' @noRd
NULL

.safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

.safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) stats::sd(x) else NA_real_
}

# Check data balance
.check_balance <- function(data, group_vars) {
  if (!length(group_vars)) return(list(ratio = NA_real_, is_balanced = TRUE, warning = NULL))
  if (!nrow(data)) return(list(ratio = NA_real_, is_balanced = TRUE, warning = NULL))
  factor_vars <- group_vars[sapply(group_vars, function(v) {
    !is.numeric(data[[v]]) || is.factor(data[[v]])
  })]
  if (!length(factor_vars)) return(list(ratio = NA_real_, is_balanced = TRUE, warning = NULL))
  counts <- stats::aggregate(rep(1L, nrow(data)),
    by = lapply(factor_vars, function(v) as.factor(data[[v]])), FUN = sum)
  x <- counts$x
  ratio <- if (length(x) > 1L) max(x) / min(x) else 1L
  warn <- NULL
  if (ratio > 5) {
    warn <- paste0("Data is highly unbalanced (max/min cell ratio = ",
      .format_num(ratio, 2), "). Consider using REML estimation.")
  } else if (ratio > 2) {
    warn <- paste0("Data is moderately unbalanced (max/min cell ratio = ",
      .format_num(ratio, 2), ").")
  }
  list(ratio = ratio, is_balanced = ratio <= 2, warning = warn)
}

# Build data description table
.make_data_info <- function(data, response, by) {
  if (is.null(by)) {
    samples <- data.frame(
      sample = "all", n = sum(!is.na(data[[response]])),
      mean = .safe_mean(data[[response]]),
      sd = .safe_sd(data[[response]]),
      cv = .safe_sd(data[[response]]) / .safe_mean(data[[response]]) * 100,
      missing = sum(is.na(data[[response]])), stringsAsFactors = FALSE)
    return(list(samples = samples, levels = list(), missing = samples,
      balance = .check_balance(data, NULL)))
  }
  if (length(by) == 1L) {
    grp <- data[[by]]
  } else {
    grp <- interaction(data[by], drop = TRUE)
  }
  split_data <- split(data, grp, drop = TRUE)
  samples <- do.call(rbind, lapply(names(split_data), function(sname) {
    df <- split_data[[sname]]
    y <- df[[response]]
    data.frame(sample = sname,
      n = sum(!is.na(y)), mean = .safe_mean(y), sd = .safe_sd(y),
      cv = .safe_sd(y) / .safe_mean(y) * 100, missing = sum(is.na(y)),
      stringsAsFactors = FALSE)
  }))
  rownames(samples) <- NULL
  factor_cols <- setdiff(names(data), c(response, by))
  levels_info <- list()
  for (fc in factor_cols) {
    levels_info[[fc]] <- table(data[[fc]], useNA = "ifany")
  }
  list(samples = samples, levels = levels_info,
    missing = samples[, c("sample", "missing")],
    balance = .check_balance(data, factor_cols))
}


# Main analysis function

#' Precision variance component analysis
#'
#' For each sample (grouped by \code{by}), estimate variance components via
#' \code{VCA::anovaVCA}, compute SD / CV and their Satterthwaite confidence intervals.
#'
#' @param data  Data frame containing the response variable and design factor columns.
#' @param form  Formula parsed by VCA (e.g. \code{y ~ day/run/rep}).
#' @param by    Column name(s) for sample grouping, supports multiple columns.
#' @param NegVC Allow negative variance components? Default \code{FALSE}.
#' @param ...   Additional arguments passed to \code{VCA::anovaVCA}.
#' @return Returns an object of class \code{precision}.
#'   Use \code{\link{summary}}, \code{\link{profile}}, etc. to view results.
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   precision(VCAdata1, y ~ day/run, by = "sample")
#'
#'   # More complex: deeply nested factors with more sample groups
#'   data(realData, package = "VCA")
#'   precision(realData, y ~ lot/calibration/day/run, by = "PID")
precision <- function(data, form, by = NULL, NegVC = FALSE, ...) {
  stopifnot(is.data.frame(data), !is.null(form), inherits(form, "formula"))
  response <- as.character(form)[2L]
  stopifnot(response %in% names(data))
  formula_vars <- all.vars(form)
  missing_vars <- setdiff(formula_vars, names(data))
  if (length(missing_vars)) {
    stop(
      "Variable(s) in 'form' not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  # Split by sample
  if (is.null(by)) {
    split_list <- list(all = data)
  } else {
    if (is.character(by) && length(by) == 1L) {
      stopifnot(by %in% names(data))
      split_list <- split(data, data[[by]], drop = TRUE)
    } else if (is.character(by) && length(by) >= 2L) {
      for (nm in by) stopifnot(nm %in% names(data))
      grp <- interaction(data[by], drop = TRUE)
      split_list <- split(data, grp, drop = TRUE)
    } else {
      stop("'by' must be NULL or a character vector.")
    }
  }

  # Unified grouping key (single/multi-column by shared)
  by_key <- if (is.null(by)) {
    NULL
  } else if (length(by) == 1L) {
    as.character(data[[by]])
  } else {
    as.character(interaction(data[by], drop = TRUE))
  }

  # data_info (pure descriptive statistics, no VCA needed)
  data_info <- .make_data_info(data, response, by)

  # print slot (descriptive metadata)
  factor_cols <- setdiff(names(data), c(response, by))
  factor_level_info <- list()
  for (fc in factor_cols) {
    if (is.factor(data[[fc]]) || is.character(data[[fc]])) {
      factor_level_info[[fc]] <- table(data[[fc]], useNA = "ifany")
    }
  }
  balance_info <- .check_balance(data, factor_cols)

  print_info <- list(
    data_class    = class(data)[1L],
    n_total       = nrow(data),
    n_samples     = length(split_list),
    var_names     = list(response = response, by = by, form = deparse(form)),
    factor_levels = factor_level_info,
    balance       = balance_info
  )

  out <- list(call = match.call(), form = form, by = by, response = response,
    data = data, by_key = by_key, data_info = data_info, print = print_info,
    results = NULL, outlier = NULL, normal = NULL, vc = NULL, ci = NULL,
    profile = NULL)
  class(out) <- "precision"

  invisible(out)
}


#' @title Print a precision object overview
#' @description
#' Print the data overview of a \code{precision} object (data class, total rows, formula,
#' grouping variables, sample count, factor levels, balance warnings) and the analysis
#' slot status (\code{[x]} / \code{[ ]}).
#'
#' @param x \code{precision} object
#' @param ... Reserved arguments
#'
#' @return Invisible \code{precision} object
#'
#' @export
print.precision <- function(x, ...) {
  p <- x$print

  cat("Precision analysis object\n")
  cat(sprintf("  %-20s: %s\n", "Data class", p$data_class))
  cat(sprintf("  %-20s: %d\n", "Total rows", p$n_total))
  cat(sprintf("  %-20s: %s\n", "Formula", p$var_names$form))
  if (!is.null(p$var_names$by))
    cat(sprintf("  %-20s: %s\n", "Group by",
      paste(p$var_names$by, collapse = " + ")))
  cat(sprintf("  %-20s: %d\n", "Samples", p$n_samples))

  # Factor levels
  if (length(p$factor_levels)) {
    cat("\nFactor levels\n")
    for (nm in names(p$factor_levels))
      cat(sprintf("  %-12s: %s\n", nm,
        paste(names(p$factor_levels[[nm]]), " (", p$factor_levels[[nm]], ")", sep = "",
              collapse = ", ")))
  }

  # Balance
  if (!is.null(p$balance$warning))
    cat("\nNote:", p$balance$warning, "\n")

  # Analysis slot status
  cat("\nAnalysis status\n\n")
  ok <- function(val) if (is.null(val)) "[ ]" else "[x]"
  cat(sprintf("  %s outlier\n",      ok(x$outlier)))
  cat(sprintf("  %s normal\n",       ok(x$normal)))
  cat(sprintf("  %s variance\n",     ok(x$results)))
  cat(sprintf("  %s ci\n",           ok(x$ci)))
  cat(sprintf("  %s profile\n",      ok(x$profile)))

  invisible(x)
}


#' S3 plot method
#'
#' Draw one of five plot types for a \code{precision} object. \code{dot} (run-order scatter plot) and
#' \code{his} (histogram + N(0,1) density) require no prior analysis; \code{var} (variance component bar chart)
#' requires \code{variance()} to be run first; \code{qq} (normal Q-Q plot) requires \code{normal()} first;
#' \code{profile} (precision profile) requires \code{profile()} first.
#'
#' @param x \code{precision} object
#' @param type Plot type: \code{"dot"} (run-order scatter plot, default preferred),
#'   \code{"his"} (histogram), \code{"var"} (variance component plot),
#'   \code{"qq"} (normal Q-Q plot), \code{"profile"} (precision profile)
#' @param ... Reserved arguments
#'
#' @return Invisible \code{ggplot} object
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   plot(obj, type = "dot")
#'   plot(obj, type = "his")
plot.precision <- function(x, type = c("dot", "his", "var", "qq", "profile"), ...) {
  type <- match.arg(type)
  .require_pkg("ggplot2")

  if (type == "dot") {
    # Scatter plot: show y values in run order (no x$results needed)
    if (is.null(x$by)) {
      sample_names <- "all"
    } else {
      sample_names <- unique(x$by_key)
    }
    plot_data <- do.call(rbind, lapply(sample_names, function(sname) {
      rows <- if (is.null(x$by)) seq_len(nrow(x$data)) else which(x$by_key == sname)
      y <- x$data[[x$response]][rows]
      y <- y[is.finite(y)]
      data.frame(sample = sname, index = seq_along(y), y = y,
                 stringsAsFactors = FALSE)
    }))
    means <- do.call(rbind, lapply(sample_names, function(sname) {
      rows <- if (is.null(x$by)) seq_len(nrow(x$data)) else which(x$by_key == sname)
      y <- x$data[[x$response]][rows]
      data.frame(sample = sname, mean = mean(y, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$index, y = .data$y)) +
      ggplot2::geom_line(color = "steelblue", linewidth = 0.4, na.rm = TRUE) +
      ggplot2::geom_point(size = 1, na.rm = TRUE) +
      ggplot2::geom_hline(data = means, ggplot2::aes(yintercept = .data$mean),
        color = "darkorange", linewidth = 0.6, linetype = "dashed") +
      ggplot2::facet_wrap(~ sample, scales = "free_y") +
      ggplot2::labs(x = "Index", y = x$response,
           title = "Run-order scatter plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  } else if (type == "his") {
    # Histogram: standardized data + N(0,1) density curve (no x$results needed)
    if (is.null(x$by)) {
      sample_names <- "all"
    } else {
      sample_names <- unique(x$by_key)
    }
    plot_data <- do.call(rbind, lapply(sample_names, function(sname) {
      rows <- if (is.null(x$by)) seq_len(nrow(x$data)) else which(x$by_key == sname)
      y <- x$data[[x$response]][rows]
      y <- y[is.finite(y)]
      if (length(y) < 3L) return(NULL)
      z <- (y - mean(y)) / stats::sd(y)
      data.frame(sample = sname, z = z, stringsAsFactors = FALSE)
    }))

    seq_z <- seq(-4, 4, length.out = 200)
    norm_curve <- data.frame(z = seq_z, density = stats::dnorm(seq_z))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$z)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = 30, fill = "steelblue", alpha = 0.4, color = "white", na.rm = TRUE) +
      ggplot2::geom_density(color = "steelblue", linewidth = 0.5,
        na.rm = TRUE) +
      ggplot2::geom_line(data = norm_curve,
        ggplot2::aes(x = .data$z, y = .data$density),
        color = "darkorange", linewidth = 0.5, na.rm = TRUE) +
      ggplot2::facet_wrap(~ sample, scales = "free") +
      ggplot2::labs(x = "Standardized value (z)", y = "Density",
           title = "Histogram with N(0,1) density") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  } else if (type == "var") {
    # Variance component bar chart: %Total per sample
    if (is.null(x$results)) {
      cat("  [Info] Run variance() first to compute variance components.\n")
      return(invisible(x))
    }

    plot_data <- do.call(rbind, lapply(x$results, function(r) {
      if (is.null(r$anova)) return(NULL)
      an <- r$anova
      data.frame(sample = r$sample,
                 component = rownames(an),
                 pct_total = an[, "%Total"],
                 stringsAsFactors = FALSE)
    }))
    plot_data$component <- factor(plot_data$component,
      levels = rev(sort(unique(plot_data$component))))

    p <- ggplot2::ggplot(plot_data,
        ggplot2::aes(x = .data$sample, y = .data$pct_total,
                     fill = .data$component)) +
      ggplot2::geom_col(width = 0.6, color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::labs(x = "Sample", y = "%Total",
           title = "Variance components by sample",
           fill = "Component") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5))

  } else if (type == "qq") {
    # QQ plot (with confidence band) -- requires normal() first
    if (is.null(x$normal)) {
      cat("  [Info] Run normal() first to compute normality test results.\n")
      return(invisible(x))
    }

    qq_list <- lapply(names(x$normal), function(sname) {
      rows <- if (is.null(x$by)) seq_len(nrow(x$data)) else which(x$by_key == sname)
      y <- x$data[[x$response]][rows]
      y <- y[is.finite(y)]
      n <- length(y)
      if (n < 3L) return(NULL)

      y_sorted <- sort(y)
      p <- (seq_len(n) - 0.5) / n
      theoretical <- stats::qnorm(p)

      q_y <- stats::quantile(y, c(0.25, 0.75), names = FALSE)
      q_t <- stats::qnorm(c(0.25, 0.75))
      slope <- diff(q_y) / diff(q_t)
      intercept <- q_y[1L] - slope * q_t[1L]
      predicted <- intercept + slope * theoretical

      z <- stats::qnorm(0.975)
      se <- slope * sqrt(p * (1 - p) / n) / stats::dnorm(theoretical)
      upper <- predicted + z * se
      lower <- predicted - z * se

      data.frame(sample = sname, theoretical = theoretical,
                 observed = y_sorted, predicted = predicted,
                 lower = lower, upper = upper,
                 stringsAsFactors = FALSE)
    })
    qq_data <- do.call(rbind, qq_list)

    ann_data <- do.call(rbind, lapply(names(x$normal), function(sname) {
      pv <- x$normal[[sname]]$p_value
      lbl <- if (is.na(pv)) "N too small" else sprintf("Shapiro-Wilk p = %.4f", pv)
      data.frame(sample = sname, label = lbl, stringsAsFactors = FALSE)
    }))

    p <- ggplot2::ggplot(qq_data, ggplot2::aes(x = .data$theoretical)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = 0.2, fill = "steelblue", na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = .data$lower),
        color = "steelblue", linewidth = 0.4, linetype = "dotted", na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = .data$upper),
        color = "steelblue", linewidth = 0.4, linetype = "dotted", na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = .data$predicted),
        color = "steelblue", linewidth = 0.6, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = .data$observed), size = 1, na.rm = TRUE) +
      ggplot2::geom_text(data = ann_data,
        ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
        hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE) +
      ggplot2::facet_wrap(~ sample, scales = "free_y") +
      ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
           title = "Normal Q-Q plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  } else if (type == "profile") {
    # Precision profile -- requires profile() first
    if (is.null(x$profile) || !length(x$profile$fits)) {
      cat("  [Info] Run profile() first to fit precision profiles.\n")
      return(invisible(x))
    }

    fits <- x$profile$fits

    plot_list <- lapply(fits, function(f) {
      d <- f$data
      rbind(
        data.frame(mean = d$mean, sample = d$sample, component = d$component,
                   metric = "SD", value = d$SD, stringsAsFactors = FALSE),
        data.frame(mean = d$mean, sample = d$sample, component = d$component,
                   metric = "CV", value = d$CV, stringsAsFactors = FALSE)
      )
    })
    plot_data <- do.call(rbind, plot_list)
    plot_data$metric <- factor(plot_data$metric, levels = c("SD", "CV"))

    pred_list <- lapply(fits, function(f) {
      p <- f$pred
      if (is.null(p)) return(NULL)
      rbind(
        data.frame(mean = p$mean, component = p$component,
                   metric = "SD", fit = p$fit, lwr = p$lwr, upr = p$upr,
                   stringsAsFactors = FALSE),
        data.frame(mean = p$mean, component = p$component,
                   metric = "CV",
                   fit = p$fit / p$mean * 100,
                   lwr = p$lwr / p$mean * 100,
                   upr = p$upr / p$mean * 100,
                   stringsAsFactors = FALSE)
      )
    })
    pred_data <- do.call(rbind, pred_list)
    pred_data$metric <- factor(pred_data$metric, levels = c("SD", "CV"))

    eq_data <- do.call(rbind, lapply(names(fits), function(comp) {
      f <- fits[[comp]]
      template <- .MODEL_FORMULAS[as.character(f$best_no)]
      best_key <- paste0("model", sub("Model_", "", f$best))
      best_mod <- f$models[[best_key]]
      eq <- template
      eq <- gsub(" = ", " == ", eq, fixed = TRUE)
      if (!is.null(best_mod)) {
        s <- tryCatch(summary(best_mod), error = function(e) NULL)
        if (!is.null(s) && !is.null(s$coefficients)) {
          est <- s$coefficients[, "Estimate"]
          for (i in seq_along(est))
            eq <- gsub(paste0("b", i), formatC(est[i], digits = 4, format = "f"), eq, fixed = TRUE)
        }
      }
      data.frame(component = comp, metric = "SD", label = eq, stringsAsFactors = FALSE)
    }))
    eq_data$metric <- factor(eq_data$metric, levels = c("SD", "CV"))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$mean)) +
      ggplot2::geom_ribbon(data = pred_data,
        ggplot2::aes(y = .data$fit, ymin = .data$lwr, ymax = .data$upr),
        alpha = 0.2, fill = "steelblue", na.rm = TRUE) +
      ggplot2::geom_line(data = pred_data,
        ggplot2::aes(y = .data$fit),
        color = "steelblue", linewidth = 1, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = .data$value), size = 1.2, na.rm = TRUE) +
      ggplot2::geom_text(ggplot2::aes(y = .data$value, label = .data$sample),
        vjust = -0.8, size = 2.5, check_overlap = TRUE, na.rm = TRUE) +
      ggplot2::geom_text(data = eq_data,
        ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
        hjust = -0.08, vjust = 1.5, size = 2.5, parse = TRUE, inherit.aes = FALSE) +
      ggplot2::facet_grid(metric ~ component, scales = "free_y") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.25))) +
      ggplot2::labs(x = "Mean", y = NULL, title = "Precision profiles") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  print(p)
  invisible(p)
}


#' @title Print outlier detection results
#' @description
#' Prints the outlier detection method, count of outliers found, and the detail table
#' (row, sample, value, statistic, critical value).
#' @param x \code{precision_outlier} object
#' @param ... Reserved arguments
#' @return Invisibly returns \code{x}.
#' @export
print.precision_outlier <- function(x, ...) {
  cat(sprintf("Outlier detection -- %s\n", attr(x, "method")))
  if (x$n == 0L) {
    cat("  No outliers detected.\n")
    return(invisible(x))
  }
  cat(sprintf("  %d outlier(s) found:\n", x$n))
  .print_df(x$data[, c("row", "sample", "value", "stat", "critical")])
  invisible(x)
}

#' S3 outlier detection method
#'
#' Detect outliers in the response values for each sample, supporting Grubbs test and IQR method.
#' Underlying call to \code{outliers_test()} (from \code{outliers-and-normal.R}).
#'
#' @param x \code{precision} object
#' @param alpha Significance level, default \code{0.05} (Grubbs only)
#' @param method Detection method: \code{"grubbs"} (default) or \code{"iqr"}
#' @param ... Additional arguments passed to \code{outliers_test()}
#'
#' @return Updated \code{precision} object with results stored in \code{$outlier}
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   obj <- outlier(obj, method = "iqr")
outlier.precision <- function(x, alpha = 0.05, method = c("grubbs", "iqr"), ...) {
  method <- match.arg(method)
  data <- x$data; response <- x$response

  # Get sample names (independent of x$results)
  if (is.null(x$by)) {
    sample_names <- "all"
  } else {
    sample_names <- unique(x$by_key)
  }

  out_list <- list()
  for (sname in sample_names) {
    # Global row indices for the current sample
    global_rows <- if (is.null(x$by)) {
      seq_len(nrow(data))
    } else {
      which(x$by_key == sname)
    }
    df <- data[global_rows, , drop = FALSE]
    finite <- is.finite(df[[response]])
    group_rows <- global_rows[finite]
    y <- df[[response]][finite]
    if (length(y) < 3L) next

    # Use external outliers_test wrapper function
    res_test <- outliers_test(data.frame(y = y), col = "y",
      method = method, alpha = alpha, ...)

    if (method == "grubbs") {
      # Grubbs: detects at most 1 outlier
      if (!is.null(res_test$details$index) && is.finite(res_test$details$index)) {
        out_list[[length(out_list) + 1L]] <- data.frame(sample = sname,
          row = group_rows[res_test$details$index], value = res_test$values[1L],
          stat = res_test$details$statistic,
          critical = res_test$details$critical,
          method = "Grubbs", stringsAsFactors = FALSE)
      }
    } else {
      # IQR: may detect multiple outliers
      for (j in seq_along(res_test$indices))
        out_list[[length(out_list) + 1L]] <- data.frame(sample = sname,
          row = group_rows[res_test$indices[j]], value = res_test$values[j],
          stat = NA_real_, critical = NA_real_,
          method = "IQR", stringsAsFactors = FALSE)
    }
  }

  # Store as structured slot and print via S3 method
  method_label <- if (method == "grubbs") sprintf("Grubbs (alpha = %.2f)", alpha) else "IQR"
  if (!length(out_list)) {
    x$outlier <- structure(list(n = 0L, rows = integer(0L), data = data.frame()),
      method = method_label, class = "precision_outlier")
    print(x$outlier)
    return(invisible(x))
  }
  out_df <- do.call(rbind, out_list)
  x$outlier <- structure(list(n = nrow(out_df), rows = out_df$row, data = out_df),
    method = method_label, class = "precision_outlier")
  print(x$outlier)
  invisible(x)
}


#' @title Print normality test results
#' @description
#' Prints the selected normality test result per sample, including the actual
#' method, sample size, statistic name, statistic, and p-value.
#' @param x \code{precision_normal} object
#' @param ... Reserved arguments
#' @return Invisibly returns \code{x}.
#' @export
print.precision_normal <- function(x, ...) {
  method_label <- function(method) {
    switch(method,
      shapiro = "Shapiro-Wilk", sw = "Shapiro-Wilk",
      ad = "Anderson-Darling", lillie = "Lilliefors",
      cvm = "Cramer-von Mises", "Unknown")
  }
  stat_name <- function(method) {
    switch(method,
      shapiro = "W", sw = "W", ad = "A", lillie = "D", cvm = "W", "statistic")
  }

  methods <- vapply(x, function(item) {
    if (!is.list(item) || is.null(item$n)) return(NA_character_)
    item$test_method %||% NA_character_
  }, character(1))
  methods <- methods[!is.na(methods)]
  if (length(unique(methods)) == 1L) {
    cat(sprintf("Normality test -- %s\n", method_label(unique(methods))))
  } else {
    cat("Normality test\n")
  }
  cat(sprintf("  %-12s  %8s  %-18s  %10s  %s\n", "Sample", "n", "Method", "Statistic", "p-value"))
  cat(sprintf("  %s\n", paste(rep("-", 75), collapse = "")))
  for (nm in names(x)) {
    if (!is.list(x[[nm]]) || is.null(x[[nm]]$n)) next
    item <- x[[nm]]
    method <- item$test_method %||% if (!is.null(item$W)) "shapiro" else NA_character_
    stat <- item$statistic %||% item$W
    name <- item$statistic_name %||% stat_name(method)
    pv <- item$p_value
    nn <- item$n
    method_str <- if (is.na(method)) "--" else method_label(method)
    stat_str <- if (is.na(stat)) "--" else sprintf("%s=%.4f", name, stat)
    pv_str <- if (is.na(pv)) "  --" else sprintf("%.4f", pv)
    cat(sprintf("  %-12s  %8d  %-18s  %10s  %s\n", nm, nn, method_str, stat_str, pv_str))
  }
  invisible(x)
}

#' S3 normality test method
#'
#' Perform normality tests on the response values for each sample, underlying call to \code{normal_test()}.
#' Supports automatic test selection (Shapiro-Wilk, Anderson-Darling, Lilliefors, CVM).
#'
#' @param x \code{precision} object
#' @param method Test method: \code{"auto"} (default, automatic selection), \code{"shapiro"},
#'   \code{"sw"}, \code{"ad"}, \code{"lillie"}, \code{"cvm"}
#' @param level Confidence level, default \code{0.95}
#' @param ... Reserved arguments
#'
#' @return Updated \code{precision} object with results stored in \code{$normal}.
#'   Each sample entry contains \code{test_method}, \code{statistic_name},
#'   \code{statistic}, \code{p_value}, and \code{n}; \code{W} is retained as a
#'   legacy alias only for Shapiro-Wilk results.
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   obj <- normal(obj)
normal.precision <- function(x, method = c("auto", "shapiro", "sw", "ad", "lillie", "cvm"),
                             level = 0.95, ...) {
  method <- match.arg(method)
  data <- x$data; response <- x$response

  # Get sample names (independent of x$results)
  if (is.null(x$by)) {
    sample_names <- "all"
  } else {
    sample_names <- unique(x$by_key)
  }

  out <- list()
  for (sname in sample_names) {
    rows <- if (is.null(x$by)) seq_len(nrow(data)) else which(x$by_key == sname)
    y <- data[[response]][rows]
    y <- y[is.finite(y)]
    n <- length(y)

    if (n >= 3L && n <= 5000L && stats::sd(y) > 0) {
      nt <- normal_test(data.frame(y = y), col = "y", method = method, level = level)
      test_method <- nt$test_method
      statistic <- nt$statistic
      statistic_name <- switch(test_method,
        shapiro = "W", sw = "W", ad = "A", lillie = "D", cvm = "W",
        "statistic"
      )
      pv <- nt$p_value
    } else {
      test_method <- NA_character_
      statistic <- NA_real_
      statistic_name <- NA_character_
      pv <- NA_real_
    }
    out[[sname]] <- list(
      test_method = test_method,
      statistic_name = statistic_name,
      statistic = statistic,
      W = if (identical(test_method, "shapiro") || identical(test_method, "sw")) statistic else NA_real_,
      p_value = pv,
      n = n
    )
  }

  x$normal <- structure(out, class = "precision_normal")
  print(x$normal)
  invisible(x)
}


#' @title Print variance component table
#' @description
#' Prints the variance component summary table with columns: sample, component,
#' VC (variance component), \%Total, SD, and CV[\%].
#' @param x \code{precision_vc} object (a data frame)
#' @param ... Reserved arguments
#' @return Invisibly returns \code{x}.
#' @export
print.precision_vc <- function(x, ...) {
  cat("Variance components\n")
  .print_df(x)
  invisible(x)
}

#' S3 variance component estimation method
#'
#' For each sample, run VCA (\code{VCA::anovaVCA} for balanced designs / ANOVA,
#' or \code{VCA::remlVCA} for REML) to estimate variance components, compute SD,
#' CV%, and their percentage of total variance. Results are stored in the
#' \code{$results} and \code{$vc} slots.
#'
#' @param x \code{precision} object
#' @param digits Number of decimal places, default \code{4}
#' @param method Estimation method: \code{"auto"} (default) chooses per sample
#'   group \code{anovaVCA} when the design is balanced and \code{remlVCA}
#'   otherwise, using \code{VCA::isBalanced()} (VCA's own balancedness
#'   criterion); \code{"anova"} forces ANOVA, \code{"reml"} forces REML.
#' @param ... Additional arguments passed to \code{VCA::anovaVCA()} /
#'   \code{VCA::remlVCA()} (e.g. \code{NegVC}, \code{conf.level})
#'
#' @return Updated \code{precision} object, \code{$results} contains per-sample VCA results,
#'   \code{$vc} is the summary variance component table
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   obj <- variance(obj)
variance.precision <- function(x, digits = 4, method = c("auto", "anova", "reml"), ...) {
  .require_pkg("VCA")
  method <- match.arg(method)

  form <- x$form
  data <- x$data
  response <- x$response
  by <- x$by

  # Split by sample
  if (is.null(by)) {
    split_list <- list(all = data)
  } else {
    if (length(by) == 1L) {
      split_list <- split(data, data[[by]], drop = TRUE)
    } else {
      grp <- interaction(data[by], drop = TRUE)
      split_list <- split(data, grp, drop = TRUE)
    }
  }

  # Run VCA per group
  results <- lapply(names(split_list), function(sname) {
    df <- split_list[[sname]]
    y <- df[[response]]
    y_ok <- !is.na(y)
    df <- df[y_ok, , drop = FALSE]
    y <- y[y_ok]

    # resolve estimation method for this group (balancedness judged by VCA::isBalanced)
    resolved <- if (method == "auto") {
      if (tryCatch(VCA::isBalanced(form, df), error = function(e) TRUE)) "anova" else "reml"
    } else {
      method
    }
    fit_args <- list(...)
    if (resolved == "reml") {
      # remlVCA has no ANOVA-only arguments (NegVC, VarVC.method, MME)
      fit_args[c("NegVC", "VarVC.method", "MME")] <- NULL
    }
    fit_fun <- getExportedValue("VCA", if (resolved == "reml") "remlVCA" else "anovaVCA")
    vca_obj <- tryCatch(
      suppressMessages(do.call(fit_fun, c(list(form = form, Data = df), fit_args))),
      error = function(e) {
        warning("VCA::", resolved, " failed for sample '", sname, "': ", e$message)
        NULL
      })

    vc <- NULL; vcov <- NULL; an <- NULL; sd_comp <- NULL; cv_comp <- NULL

    if (!is.null(vca_obj)) {
      tab <- vca_obj$aov.tab
      if ("total" %in% rownames(tab))
        tab <- tab[setdiff(rownames(tab), "total"), , drop = FALSE]
      an <- as.data.frame(tab)
      vc <- tab[, "VC"]
      sd_comp <- tab[, "SD"]
      cv_comp <- tab[, "CV[%]"]
      names(vc) <- names(sd_comp) <- names(cv_comp) <- rownames(tab)

      vcov <- tryCatch(VCA::vcovVC(vca_obj), error = function(e) NULL)
      if (!is.null(vcov)) {
        idx <- intersect(rownames(vcov), names(vc))
        vcov <- vcov[idx, idx, drop = FALSE]
      }
    }

    list(sample = sname, n = sum(y_ok), mean = .safe_mean(y),
      sd = .safe_sd(y), cv = .safe_sd(y) / .safe_mean(y) * 100,
      anova = an, vc = vc, sd_comp = sd_comp, cv_comp = cv_comp,
      vcov = vcov, vca = vca_obj, method = resolved)
  })

  names(results) <- names(split_list)
  x$results <- results

  # Variance component table
  vc_rows <- list()
  for (r in results) {
    if (is.null(r$anova)) next
    an <- r$anova
    for (i in seq_len(nrow(an))) {
      vc_rows[[length(vc_rows) + 1L]] <- data.frame(sample = r$sample,
        component = rownames(an)[i], VC = an[i, "VC"],
        `%Total` = an[i, "%Total"], SD = an[i, "SD"], `CV[%]` = an[i, "CV[%]"],
        Method = r$method,
        check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  if (length(vc_rows)) {
    vc_df <- do.call(rbind, vc_rows)
    vc_df$VC <- .format_num(vc_df$VC, digits)
    vc_df$`%Total` <- .format_num(vc_df$`%Total`, 1)
    vc_df$SD <- .format_num(vc_df$SD, digits)
    vc_df$`CV[%]` <- .format_num(vc_df$`CV[%]`, digits)
    x$vc <- structure(vc_df, class = c("precision_vc", "data.frame"))
    print(x$vc)
  }

  invisible(x)
}

#' Backward-compatible alias for variance.precision
#'
#' @param x A `precision` object.
#' @param digits Number of decimal places, default `4`.
#' @param method Estimation method: `"auto"`, `"anova"`, or
#'   `"reml"`.
#' @param ... Additional arguments passed to the underlying VCA method.
#' @return The same result as \code{\link{variance.precision}}: an updated
#'   \code{precision} object, invisibly.
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   obj <- vc(obj)
#' @export
vc.precision <- variance.precision

.precision_term_key <- function(x) {
  x <- as.character(x)
  vapply(strsplit(x, ":", fixed = TRUE), function(parts) {
    paste(sort(parts), collapse = ":")
  }, character(1))
}

.precision_component_index <- function(available, requested) {
  match(.precision_term_key(requested), .precision_term_key(available))
}

.precision_combined_interval <- function(result, components, alpha = 0.05) {
  if (is.null(result$vc) || is.null(names(result$vc)))
    stop("Variance components are unavailable for sample '", result$sample,
         "'.", call. = FALSE)

  vc_index <- .precision_component_index(names(result$vc), components)
  if (anyNA(vc_index)) {
    stop("Cannot form the requested combined variance for sample '",
         result$sample, "'. Missing component(s): ",
         paste(components[is.na(vc_index)], collapse = ", "), ".",
         call. = FALSE)
  }
  estimate <- sum(as.numeric(result$vc[vc_index]))

  covariance <- result$vcov
  if (is.null(covariance) || is.null(rownames(covariance)) ||
      is.null(colnames(covariance))) {
    stop("The variance-component covariance matrix is unavailable for sample '",
         result$sample, "'.", call. = FALSE)
  }
  row_index <- .precision_component_index(rownames(covariance), components)
  col_index <- .precision_component_index(colnames(covariance), components)
  if (anyNA(row_index) || anyNA(col_index)) {
    stop("The variance-component covariance matrix does not contain all ",
         "components required for sample '", result$sample, "'.",
         call. = FALSE)
  }
  variance <- sum(covariance[row_index, col_index, drop = FALSE])

  if (!is.finite(estimate) || estimate < 0 ||
      !is.finite(variance) || variance <= 0) {
    return(c(estimate = sqrt(max(estimate, 0)), lower = NA_real_,
             upper = NA_real_, df = NA_real_))
  }

  df <- 2 * estimate^2 / variance
  lower_variance <- df * estimate / stats::qchisq(1 - alpha / 2, df)
  upper_variance <- df * estimate / stats::qchisq(alpha / 2, df)
  c(
    estimate = sqrt(estimate),
    lower = sqrt(max(lower_variance, 0)),
    upper = sqrt(max(upper_variance, 0)),
    df = df
  )
}

#' Generate an EP05 precision summary table
#'
#' @param x A precision object after \code{variance()} has been run.
#' @param day Name of the day (or equivalent) design variable.
#' @param run Name of the run variable, when present.
#' @param site Name of the site variable for a multisite design.
#' @param table One of \code{"simple"}, \code{"detailed"}, or \code{"CI"}.
#' @param digits Number of decimal places.
#' @param ... Reserved for future extensions.
#' @return A data.frame containing the requested report table.
#' @details Interaction terms are matched independently of the textual order
#'   of their factors. For a multisite design, within-laboratory precision is
#'   formed from the within-site components (for example, `site:day + error`),
#'   whereas reproducibility includes the between-site component. With
#'   `table = "CI"`, confidence intervals for combined components use the
#'   variance-component covariance matrix and a Satterthwaite approximation.
#' @examples
#' if (requireNamespace("VCA", quietly = TRUE)) {
#'   data(VCAdata1, package = "VCA")
#'   p <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   p <- variance(p)
#'   report(p, day = "day", run = "run", table = "simple")
#' }
#' @export
report.precision <- function(x, day, run = NULL, site = NULL,
                             table = c("simple", "detailed", "CI"), digits = 4L, ...) {
  if (!inherits(x, "precision")) stop("x must be a 'precision' object.", call. = FALSE)
  if (is.null(x$results)) {
    stop("Precision variance components have not been computed. Run variance(x) before report(x).", call. = FALSE)
  }
  table <- match.arg(table, c("simple", "detailed", "CI"))
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0)
    stop("'digits' must be a non-negative number.", call. = FALSE)
  rhs <- all.vars(x$form[[3L]])
  mapped <- unique(c(day, run, site))
  if (any(!is.character(mapped)) || any(!nzchar(mapped)))
    stop("day/run/site must be variable names.", call. = FALSE)
  extra <- setdiff(rhs, mapped)
  if (length(extra)) stop("Unmapped design variables on the right-hand side of the formula: ",
                          paste(extra, collapse = ", "),
                          ". Map each variable explicitly to day, run, or site.", call. = FALSE)
  multi <- !is.null(site)
  if (is.null(day) || (!multi && is.null(run)))
    stop("Specify 'day' explicitly. A single-site design also requires 'run'; ",
         "a multisite design requires 'site'.", call. = FALSE)
  term_labels <- attr(stats::terms(x$form), "term.labels")
  observed <- sort(.precision_term_key(term_labels))
  expected_raw <- if (!multi) c(day, paste(day, run, sep = ":")) else {
    z <- c(site, paste(site, day, sep = ":"))
    if (!is.null(run)) z <- c(z, paste(site, day, run, sep = ":"))
    z
  }
  expected <- sort(.precision_term_key(expected_raw))
  if (!setequal(observed, expected))
    stop("The formula does not match a supported EP05 nested design. Observed terms: ",
         paste(observed, collapse = ", "), "; expected terms: ",
         paste(expected, collapse = ", "), ".", call. = FALSE)

  getvc <- function(r, key) {
    if (is.null(r$vc) || is.null(names(r$vc))) return(NA_real_)
    index <- .precision_component_index(names(r$vc), key)
    if (is.na(index)) NA_real_ else as.numeric(r$vc[[index]])
  }
  cv <- function(sd, mu) ifelse(is.finite(mu) && mu > 0, 100 * sd / mu, NA_real_)
  run_key <- if (!is.null(run)) {
    if (multi) paste(site, day, run, sep = ":") else paste(day, run, sep = ":")
  } else NULL
  day_key <- if (multi) paste(site, day, sep = ":") else day
  site_key <- if (multi) site else NULL
  within_lab_keys <- c("error", run_key, day_key)
  reproducibility_keys <- c(within_lab_keys, site_key)

  rows <- lapply(x$results, function(r) {
    mu <- as.numeric(r$mean); n <- as.numeric(r$n)
    required_keys <- reproducibility_keys
    required_values <- vapply(required_keys, function(key) getvc(r, key), numeric(1))
    if (anyNA(required_values)) {
      stop("Variance components required by the EP05 mapping are missing for ",
           "sample '", r$sample, "': ",
           paste(required_keys[is.na(required_values)], collapse = ", "), ".",
           call. = FALSE)
    }
    repvc <- getvc(r, "error")
    brvc <- if (!is.null(run_key)) getvc(r, run_key) else NA_real_
    bdvc <- getvc(r, day_key)
    bsvc <- if (multi) getvc(r, site_key) else NA_real_
    wl <- sqrt(sum(vapply(within_lab_keys, function(key) getvc(r, key), numeric(1))))
    repro <- sqrt(sum(required_values))
    data.frame(`Sample Description` = r$sample, `Mean Value` = mu, N = n,
      `Repeatability SD` = sqrt(repvc), `%CV Repeatability` = cv(sqrt(repvc), mu),
      `Between-Run SD` = sqrt(brvc), `%CV Between-Run` = cv(sqrt(brvc), mu),
      `Between-Day SD` = sqrt(bdvc), `%CV Between-Day` = cv(sqrt(bdvc), mu),
      `Between-Site SD` = sqrt(bsvc), `%CV Between-Site` = cv(sqrt(bsvc), mu),
      `Within-Laboratory Precision SD` = wl, `%CV Within-Laboratory Precision` = cv(wl, mu),
      `Reproducibility SD` = repro, `%CV Reproducibility` = cv(repro, mu), check.names = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (table == "simple") {
    keep <- if (multi) c("Sample Description", "Mean Value", "Repeatability SD", "%CV Repeatability", "Reproducibility SD", "%CV Reproducibility") else c("Sample Description", "Mean Value", "Repeatability SD", "%CV Repeatability", "Within-Laboratory Precision SD", "%CV Within-Laboratory Precision")
  } else if (table == "detailed") {
    keep <- names(out)
  } else {
    if (is.null(x$ci) || !length(x$ci))
      stop("Confidence intervals have not been computed. Run ci(x) before requesting table = \"CI\".", call. = FALSE)
    component_map <- c(
      "Repeatability" = "error",
      "Between-Run" = if (!is.null(run_key)) run_key else NA_character_,
      "Between-Day" = day_key,
      "Between-Site" = if (multi) site_key else NA_character_,
      "Reproducibility" = "total"
    )
    for (ci_type in c("SD", "CV")) {
      ci_data <- x$ci[[ci_type]]
      if (is.null(ci_data) || !nrow(ci_data)) next
      ci_suffix <- if (ci_type == "SD") "SD" else "CV"
      ci_keys <- .precision_term_key(ci_data$component)
      for (nm in names(component_map)) {
        key <- component_map[[nm]]
        if (is.na(key)) next
        z <- ci_data[ci_keys == .precision_term_key(key), , drop = FALSE]
        if (nrow(z)) {
          sample_index <- match(out$`Sample Description`, z$sample)
          out[[paste(nm, ci_suffix, "Lower")]] <- z$lower[sample_index]
          out[[paste(nm, ci_suffix, "Upper")]] <- z$upper[sample_index]
        }
      }
    }

    if (multi) {
      within_lab_ci <- lapply(
        x$results, .precision_combined_interval,
        components = within_lab_keys, alpha = 0.05
      )
      within_lab_ci <- do.call(rbind, within_lab_ci)
      sample_order <- match(
        out$`Sample Description`,
        vapply(x$results, `[[`, character(1), "sample")
      )
      out[["Within-Laboratory Precision SD Lower"]] <-
        within_lab_ci[sample_order, "lower"]
      out[["Within-Laboratory Precision SD Upper"]] <-
        within_lab_ci[sample_order, "upper"]
      out[["Within-Laboratory Precision CV Lower"]] <-
        100 * within_lab_ci[sample_order, "lower"] / out$`Mean Value`
      out[["Within-Laboratory Precision CV Upper"]] <-
        100 * within_lab_ci[sample_order, "upper"] / out$`Mean Value`
    } else {
      for (ci_type in c("SD", "CV")) {
        ci_data <- x$ci[[ci_type]]
        if (is.null(ci_data) || !nrow(ci_data)) next
        z <- ci_data[ci_data$component == "total", , drop = FALSE]
        sample_index <- match(out$`Sample Description`, z$sample)
        ci_suffix <- if (ci_type == "SD") "SD" else "CV"
        out[[paste("Within-Laboratory Precision", ci_suffix, "Lower")]] <-
          z$lower[sample_index]
        out[[paste("Within-Laboratory Precision", ci_suffix, "Upper")]] <-
          z$upper[sample_index]
      }
    }
    keep <- names(out)
  }
  out <- out[, intersect(keep, names(out)), drop = FALSE]
  num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], function(z) round(z, digits))
  out
}


#' @title Print confidence intervals
#' @description
#' Prints confidence interval tables for SD and \%CV components, with columns:
#' sample, component, estimate, lower, upper.
#' @param x \code{precision_ci} object (a list of data frames, keyed by "SD" and "CV")
#' @param ... Reserved arguments
#' @return Invisibly returns \code{x}.
#' @export
print.precision_ci <- function(x, ...) {
  for (ci_type in names(x)) {
    ci_df <- x[[ci_type]]
    ci_label <- if (ci_type == "SD") "SD" else "%CV"
    cat(sprintf("Confidence intervals (%s)\n", ci_label))
    .print_df(ci_df)
    cat("\n")
  }
  invisible(x)
}

#' S3 confidence interval method
#'
#' Compute confidence intervals for each variance component (SD and %CV) based on
#' Satterthwaite approximation. Requires \code{variance()} to be run first to populate
#' the \code{$results} slot.
#'
#' @param x \code{precision} object
#' @param digits Number of decimal places, default \code{4}
#' @param ... Reserved arguments
#'
#' @return Updated \code{precision} object with results stored in \code{$ci}
#' @details This method reports intervals for the fitted model components and
#'   the total component. EP05 semantic combinations such as multisite
#'   within-laboratory precision require the `day`, `run`, and `site` mapping
#'   supplied to \code{\link{report.precision}} and are therefore constructed
#'   by `report(x, ..., table = "CI")`.
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   obj <- variance(obj)
#'   obj <- ci(obj)
ci.precision <- function(x, digits = 4, ...) {
  if (is.null(x$results)) {
    cat("  [Info] Run variance() first to compute variance components.\n")
    return(invisible(x))
  }

  ci_types <- c("SD", "CV")
  ci_tables <- list()
  for (ci_type in ci_types) {
    # Inline original .precision_ci: extract CI table per sample via VCAinference
    ci_df <- tryCatch({
      out_list <- list()
      ci_type_internal <- ci_type
      for (r in x$results) {
        if (is.null(r$vca)) next
        inf <- suppressMessages(VCA::VCAinference(r$vca, VarVC = TRUE, alpha = 0.05))
        ci_table <- if (ci_type_internal == "SD") inf$ConfInt$SD$TwoSided else inf$ConfInt$CV$TwoSided
        if (is.null(ci_table)) next
        ci_df_tmp <- as.data.frame(ci_table)
        ci_df_tmp$sample <- r$sample
        rn <- rownames(ci_table)
        ci_df_tmp$component <- rn
        est_vec <- if (ci_type_internal == "SD") r$sd_comp else r$cv_comp
        ci_df_tmp$est. <- est_vec[rn]
        if (any(rn == "total") && is.na(ci_df_tmp$est.[rn == "total"][1L])) {
          total_sd <- sqrt(sum(r$vc, na.rm = TRUE))
          ci_df_tmp$est.[rn == "total"] <- if (ci_type_internal == "SD") total_sd else total_sd / r$mean * 100
        }
        ci_df_tmp <- ci_df_tmp[, c("sample", "component", "est.", "LCL", "UCL")]
        names(ci_df_tmp) <- c("sample", "component", "estimate", "lower", "upper")
        out_list[[length(out_list) + 1L]] <- ci_df_tmp
      }
      df <- do.call(rbind, out_list)
      rownames(df) <- NULL
      df
    }, error = function(e) NULL)

    if (!is.null(ci_df) && nrow(ci_df)) {
      ci_df$estimate <- .format_num(ci_df$estimate, digits)
      ci_df$lower <- .format_num(ci_df$lower, digits)
      ci_df$upper <- .format_num(ci_df$upper, digits)
      ci_tables[[ci_type]] <- ci_df
    }
  }
  x$ci <- structure(ci_tables, class = "precision_ci")
  if (length(ci_tables)) print(x$ci)
  invisible(x)
}


#' S3 summary method
#'
#' Summarize all executed analyses in a \code{precision} object, including data description,
#' analysis status checklist, outlier detection, normality tests, variance components,
#' confidence intervals, and precision profile.
#'
#' @param object \code{precision} object
#' @param ... Reserved arguments
#'
#' @return Invisible \code{precision} object
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   summary(obj)
summary.precision <- function(object, ...) {
  if (!inherits(object, "precision"))
    stop("Input must be a 'precision' object.")

  cat("\nPrecision analysis -- summary\n")
  cat(paste(rep("-", 40), collapse = ""), "\n\n")

  # Data description
  print.precision(object)

  # Delegate to each slot's S3 print method
  if (!is.null(object$outlier)) {
    cat("\n")
    print(object$outlier)
  }
  if (!is.null(object$normal)) {
    cat("\n")
    print(object$normal)
  }
  if (!is.null(object$vc)) {
    cat("\n")
    print(object$vc)
  }
  if (!is.null(object$ci) && length(object$ci)) {
    cat("\n")
    print(object$ci)
  }
  if (!is.null(object$profile) && length(object$profile$fits)) {
    print(object$profile$fits)
  }

  cat("\n")
  invisible(object)
}


# Profile model definitions

.MODEL_FORMULAS <- c(
  "1"  = "sigma^2 = beta1",
  "2"  = "sigma^2 = beta1 * u^2",
  "3"  = "sigma^2 = beta1 + beta2 * u^2",
  "4"  = "sigma^2 = (beta1 + beta2 * u)^2",
  "5"  = "sigma^2 = beta1 + beta2 * u^K",
  "6"  = "sigma^2 = beta1 + beta2 * u + beta3 * u^J",
  "7"  = "sigma^2 = beta1 + beta2 * u^J",
  "8"  = "sigma^2 = (beta1 + beta2 * u)^J",
  "9"  = "sigma^2 = beta1 * u^J",
  "10" = "CV = beta1 * u^J")

.MODEL_TYPES <- c(
  "1"  = "Constant variance",         "2"  = "Constant CV",
  "3"  = "Mixed (const + prop var)",  "4"  = "Constrained power (const exponent)",
  "5"  = "Alternative constrained power", "6"  = "Unconstrained power (with min)",
  "7"  = "Alternative unconstrained power", "8"  = "Unconstrained power (Sadler default)",
  "9"  = "CLSI EP17-like",           "10" = "CLSI EP17 exact (log-log)")

# Internal Sadler precision profile fitting (10 candidate models)
.sadler <- function(df, type = c("sd", "cv"), model.no = 1:10, K = 2, quiet = TRUE) {
  stopifnot(is.data.frame(df), all(c("mean", "CV", "SD") %in% names(df)))
  type <- match.arg(type)
  df <- df[is.finite(df$mean) & (is.finite(df$SD) | is.finite(df$CV)), ]
  if (nrow(df) < 3L) { warning("Need at least 3 points for Sadler fitting."); return(list()) }
  .require_pkg("VFP")

  Mean <- df$mean; SD <- df$SD; CV <- df$CV
  VC <- if (type == "sd") SD^2 else (CV / 100 * Mean)^2
  DF <- if ("DF" %in% names(df)) df$DF else rep(max(1L, nrow(df) - 1L), nrow(df))
  mat <- data.frame(Mean = Mean, VC = VC, DF = DF, SD = SD, CV = CV)
  fit <- tryCatch(
    suppressMessages(VFP::fit_vfp(Data = mat, model.no = model.no, K = K, quiet = quiet,
      col.mean = "Mean", col.var = "VC", col.df = "DF", col.sd = "SD", col.cv = "CV")),
    error = function(e) return(list()))
  if (!length(fit)) return(list())
  num <- VFP::get_model(fit)
  best_name <- attr(num, "model")
  best_no <- as.integer(sub("Model_", "", best_name))

  xseq <- seq(min(Mean), max(Mean), length.out = 100L)
  pred <- VFP::predict.VFP(fit, model.no = best_no, newdata = xseq,
    type = if (type == "sd") "sd" else "cv")
  pred <- data.frame(mean = pred$Mean, fit = pred$Fitted, lwr = pred$LCL, upr = pred$UCL)

  list(best = best_name, best_no = best_no, aic = fit$AIC, models = fit$Models,
    data = df, pred = pred, type = type, vfp = fit)
}

# Helper: format model formula with coefficient estimates substituted
.format_model_str <- function(f) {
  template <- .MODEL_FORMULAS[as.character(f$best_no)]
  best_key <- paste0("model", sub("Model_", "", f$best))
  best_mod <- f$models[[best_key]]
  if (is.null(best_mod)) return(template)
  est <- tryCatch(stats::coef(best_mod), error = function(e) NULL)
  if (is.null(est) || !length(est)) return(template)
  eq <- template
  for (i in seq_along(est))
    eq <- gsub(paste0("beta", i), formatC(est[i], digits = 4, format = "f"), eq, fixed = TRUE)
  eq
}

# Helper: compute R^2 for the best model (on the variance scale)
.model_rsq <- function(f) {
  best_key <- paste0("model", sub("Model_", "", f$best))
  best_mod <- f$models[[best_key]]
  if (is.null(best_mod)) return(NA_real_)
  y_obs <- best_mod$y
  y_pred <- best_mod$fitted.values
  if (is.null(y_pred) || is.null(y_obs)) return(NA_real_)
  rss <- sum((y_obs - y_pred)^2, na.rm = TRUE)
  tss <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(tss) || tss == 0) return(NA_real_)
  1 - rss / tss
}




#' @title Print precision profile fits
#' @description
#' Prints the best Sadler model formula, AIC, and R-squared for each variance component.
#' @param x \code{precision_profile} object (a list of per-component fit results)
#' @param ... Reserved arguments
#' @return Invisibly returns \code{x}.
#' @export
print.precision_profile <- function(x, ...) {
  cat("Precision profiles -- Sadler\n")
  for (comp in names(x)) {
    f <- x[[comp]]
    eq <- .format_model_str(f)
    rsq <- .model_rsq(f)
    rsq_str <- if (is.finite(rsq)) sprintf("R^2 = %.4f", rsq) else "R^2 = --"
    cat(sprintf("  %s: %s  (AIC = %.2f, %s)\n",
      comp, eq, f$aic[f$best], rsq_str))
  }
  invisible(x)
}

#' S3 precision profile method
#'
#' Fit precision (SD) vs. mean across samples using the Sadler model selection algorithm
#' (10 candidate models) based on \code{VFP::fit_vfp()}.
#'
#' @param object \code{precision} object (requires \code{variance()} to be run first)
#' @param model.no Vector of candidate model numbers, default \code{1:10}
#' @param ... Additional arguments passed to \code{.sadler()} (e.g. \code{K})
#'
#' @return Updated \code{precision} object with results stored in \code{$profile}
#'
#' @export
#' @examples
#'   data(VCAdata1, package = "VCA")
#'   obj <- precision(VCAdata1, y ~ day/run, by = "sample")
#'   obj <- variance(obj)
#'   obj <- profile(obj, model.no = 1)
profile.precision <- function(object, model.no = 1:10, ...) {
  if (is.null(object$results)) {
    cat("  [Info] Run variance() first to compute variance components.\n")
    return(invisible(object))
  }

  res <- object$results

  # Parse all available variance components
  all_comps <- unique(unlist(lapply(res, function(r) names(r$vc))))
  all_comps <- all_comps[!is.na(all_comps) & nzchar(all_comps)]
  comps <- c("total", all_comps)

  # Fit each component (based on SD)
  fits <- list()
  for (comp in comps) {
    prof <- do.call(rbind, lapply(res, function(r) {
      if (comp == "total") {
        sd_val <- if (!is.null(r$vc)) sqrt(sum(r$vc, na.rm = TRUE)) else NA_real_
      } else {
        sd_val <- if (!is.null(r$sd_comp) && comp %in% names(r$sd_comp))
                    as.numeric(r$sd_comp[comp]) else NA_real_
      }
      if (!is.finite(sd_val)) return(NULL)
      data.frame(mean = r$mean, SD = sd_val, CV = sd_val / r$mean * 100,
                 sample = r$sample, component = comp, stringsAsFactors = FALSE)
    }))
    if (is.null(prof) || nrow(prof) < 3L) next

    fit <- NULL
    utils::capture.output(fit <- tryCatch(
      .sadler(prof, type = "sd", model.no = model.no, quiet = TRUE, ...),
      error = function(e) NULL))
    if (!length(fit) || is.null(fit$best)) next

    if (!is.null(fit$pred)) fit$pred$component <- comp
    fits[[comp]] <- fit
  }

  if (!length(fits)) {
    warning("Precision profile could not be fitted for any component.")
    object$profile <- list(fits = structure(fits, class = "precision_profile"))
    return(invisible(object))
  }

  object$profile <- list(fits = structure(fits, class = "precision_profile"))
  print(object$profile$fits)
  invisible(object)
}


# Standalone function: list Sadler model equations

#' List Sadler precision profile model equations
#'
#' @details List the 10 candidate Sadler precision profile model formulas and types
#'   used in the VFP package. Includes: Constant SD, Constant CV, Linear (variance),
#'   Power model, Exponential model, etc.
#'
#' @rawRd \arguments{This function has no arguments.}
#' @return Invisibly returns \code{NULL}. The model table is printed to the
#'   console as a side effect.
#' @export
#' @examples
#'   list_sadler()
list_sadler <- function() {
  cat("Sadler Precision Profile Models\n")
  cat("\n")
  for (i in seq_along(.MODEL_TYPES)) {
    nm <- names(.MODEL_TYPES)[i]
    cat(sprintf("  %2s  %-20s  %s\n", nm, .MODEL_TYPES[nm], .MODEL_FORMULAS[nm]))
  }
  invisible()
}
