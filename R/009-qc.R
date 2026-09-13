#' List Westgard QC rules
#'
#' Prints commonly used Westgard multi-rule QC rules and their meanings.
#' Rule names printed here can be passed to \code{qc_chart(rules = "...")}.
#'
#' @param brief Set to \code{TRUE} to suppress printing and return only
#'   the named character vector.
#'
#' @return A named character vector of rule descriptions (invisibly).
#'
#' @export
#'
#' @examples
#'   list_westgard()
list_westgard <- function(brief = FALSE) {
  rules <- c(
    "1-2s"  = "1 point exceeds +/-2s (warning)",
    "1-3s"  = "1 point exceeds +/-3s (out of control)",
    "2-2s"  = "2 consecutive points exceed +/-2s (same side)",
    "R-4s"  = "2 consecutive points differ by at least 4s",
    "4-1s"  = "4 consecutive points exceed +/-1s (same side)",
    "8x"    = "8 consecutive points on same side of mean",
    "10x"   = "10 consecutive points on same side of mean"
  )

  if (!brief) {
    cat("\nWestgard QC Rules\n")
    cat(" ", paste(rep("-", 50), collapse = ""), "\n")
    for (nm in names(rules)) {
      cat(sprintf("  %-8s%s\n", nm, rules[[nm]]))
    }
    cat(" ", paste(rep("-", 50), collapse = ""), "\n\n")
  }

  invisible(rules)
}

#' Levey-Jennings QC chart + Westgard rule evaluation
#'
#' Draws a Levey-Jennings QC chart and detects out-of-control points
#' based on selected Westgard rules. If target mean and SD are not
#' provided, they are estimated from the data.
#'
#' @param data  A data frame.
#' @param value Column name (numeric) containing measured values.
#' @param rules Comma-separated rule names, e.g. "1-2s,1-3s,10x".
#'              Rule names must match those printed by list_westgard().
#' @param mean  Target mean; NULL to estimate from the selected data column.
#' @param sd    Target SD; NULL to estimate from the selected data column.
#' @param group Optional grouping column name (e.g. QC lot), used for
#'              faceted plots.
#' @param run   Run-order column name (optional, default 1:n).
#'
#' @return An object of class \code{"qc_chart"} with violation results
#'   and plot data. Use \code{print()} to show summary and \code{plot()}
#'   to draw the Levey-Jennings chart.
#'
#' @examples
#'   set.seed(42)
#'   df <- data.frame(
#'     run   = 1:30,
#'     value = c(rnorm(27, 100, 5), 108, 115, 85)
#'   )
#'   res <- qc_chart(df, "value", rules = "1-2s,1-3s,10x")
#'   print(res)
#'   plot(res)
#' @export
qc_chart <- function(data, value, rules = "1-2s,1-3s,2-2s,R-4s,4-1s,10x",
                     mean = NULL, sd = NULL,
                     group = NULL, run = NULL) {

  # Helper functions

  `%||%` <- function(a, b) if (is.null(a)) b else a

  KNOWN_RULES <- c("1-2s", "1-3s", "2-2s", "R-4s", "4-1s", "8x", "10x")

  .apply_westgard_rules <- function(values, mean, sd, rule_names) {
    z <- (values - mean) / sd
    n <- length(values)
    result <- setNames(vector("list", length(rule_names)), rule_names)

    for (rule in rule_names) {
      idx <- integer(0L)

      if (rule == "1-2s") {
        idx <- which(abs(z) > 2)

      } else if (rule == "1-3s") {
        idx <- which(abs(z) > 3)

      } else if (rule == "2-2s") {
        hit <- rep(FALSE, n)
        if (n >= 2L) {
          for (i in seq_len(n - 1L)) {
            if ((z[i] > 2 && z[i + 1L] > 2) || (z[i] < -2 && z[i + 1L] < -2)) {
              hit[i] <- TRUE
              hit[i + 1L] <- TRUE
            }
          }
        }
        idx <- which(hit)

      } else if (rule == "R-4s") {
        hit <- rep(FALSE, n)
        if (n >= 2L) {
          for (i in seq_len(n - 1L)) {
            if (abs(z[i]) >= 2 && abs(z[i + 1L]) >= 2 &&
                abs(z[i] - z[i + 1L]) >= 4) {
              hit[i] <- TRUE
              hit[i + 1L] <- TRUE
            }
          }
        }
        idx <- which(hit)

      } else if (rule == "4-1s") {
        hit <- rep(FALSE, n)
        if (n >= 4L) {
          for (i in seq_len(n - 3L)) {
            win <- z[i:(i + 3L)]
            if (all(abs(win) > 1) && (all(win > 0) || all(win < 0))) {
              hit[i:(i + 3L)] <- TRUE
            }
          }
        }
        idx <- which(hit)

      } else if (rule == "8x") {
        hit <- rep(FALSE, n)
        if (n >= 8L) {
          for (i in seq_len(n - 7L)) {
            win <- z[i:(i + 7L)]
            if (all(win > 0) || all(win < 0)) {
              hit[i:(i + 7L)] <- TRUE
            }
          }
        }
        idx <- which(hit)

      } else if (rule == "10x") {
        hit <- rep(FALSE, n)
        if (n >= 10L) {
          for (i in seq_len(n - 9L)) {
            win <- z[i:(i + 9L)]
            if (all(win > 0) || all(win < 0)) {
              hit[i:(i + 9L)] <- TRUE
            }
          }
        }
        idx <- which(hit)
      }

      result[[rule]] <- idx
    }

    result
  }

  # Input validation
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(value) || length(value) != 1L || !value %in% names(data))
    stop("'value' must be a single column name in 'data'.", call. = FALSE)
  if (!is.numeric(data[[value]]))
    stop("Column '", value, "' must be numeric.", call. = FALSE)

  if (!is.null(group) && (!is.character(group) || length(group) != 1L || !group %in% names(data)))
    stop("'group' must be NULL or a single column name in 'data'.", call. = FALSE)
  if (!is.null(run) && (!is.character(run) || length(run) != 1L || !run %in% names(data)))
    stop("'run' must be NULL or a single column name in 'data'.", call. = FALSE)

  # Parse rule names
  if (!is.character(rules) || length(rules) != 1L)
    stop("'rules' must be a single comma-separated string.", call. = FALSE)

  rule_vec <- trimws(strsplit(rules, ",")[[1L]])
  rule_vec <- rule_vec[nchar(rule_vec) > 0]
  bad_rules <- setdiff(rule_vec, KNOWN_RULES)
  if (length(bad_rules) > 0)
    stop("Unknown rule(s): ", paste(bad_rules, collapse = ", "),
         ". See list_westgard() for valid rule names.", call. = FALSE)
  if (length(rule_vec) == 0L)
    stop("No valid rules specified.", call. = FALSE)

  # Extract data
  values_vec <- data[[value]]
  n_total    <- length(values_vec)
  run_vec    <- if (!is.null(run)) data[[run]] else seq_len(n_total)
  group_vec  <- if (!is.null(group)) data[[group]] else rep(1L, n_total)

  # Handle missing values
  ok <- !is.na(values_vec)
  n_complete <- sum(ok)
  n_miss     <- n_total - n_complete

  if (n_complete == 0L)
    stop("No complete cases in 'value' column.", call. = FALSE)

  values_clean <- values_vec[ok]
  run_clean    <- run_vec[ok]
  group_clean  <- group_vec[ok]

  # Contextual unique lines for remaining replacements

  # Mean / SD
  target_mean <- mean %||% mean(values_clean, na.rm = TRUE)
  target_sd   <- sd   %||% stats::sd(values_clean, na.rm = TRUE)

  if (target_sd <= 0)
    stop("SD must be positive.", call. = FALSE)

  # Rule detection
  violations <- .apply_westgard_rules(values_clean, target_mean, target_sd, rule_vec)

  # Global violation indices (deduplicated)
  all_violation_indices <- sort(unique(unlist(violations)))
  n_violations <- length(all_violation_indices)

  # Build per-violation rule labels (joined with "/" for display)
  pt_rules <- character(n_complete)
  for (rule_nm in rule_vec) {
    for (i in violations[[rule_nm]]) {
      pt_rules[i] <- if (nchar(pt_rules[i]) == 0) rule_nm else paste0(pt_rules[i], "/", rule_nm)
    }
  }

  # Severity classification
  if ("1-3s" %in% rule_vec) {
    is_oc <- seq_len(n_complete) %in% violations[["1-3s"]]
  } else {
    is_oc <- rep(FALSE, n_complete)
  }
  severity <- rep(NA_character_, n_complete)
  severity[all_violation_indices] <- "warning"
  severity[is_oc] <- "out_of_control"

  n_warn <- sum(severity == "warning", na.rm = TRUE)
  n_oc   <- sum(severity == "out_of_control", na.rm = TRUE)

  # Build plot_data
  plot_data <- data.frame(
    run    = run_clean,
    value  = values_clean,
    group  = factor(group_clean),
    severity = factor(severity, levels = c("warning", "out_of_control")),
    vio_label = ifelse(nchar(pt_rules) > 0, pt_rules, NA_character_),
    stringsAsFactors = FALSE
  )

  run_name   <- if (!is.null(run)) run else "Run"
  group_name <- group
  chart_title <- paste0("Levey-Jennings Chart (", value, ")")

  structure(list(
    value       = value,
    target_mean = target_mean,
    target_sd   = target_sd,
    n_total     = n_total,
    n_complete  = n_complete,
    n_miss      = n_miss,
    rule_vec    = rule_vec,
    violations  = violations,
    n_violations = n_violations,
    n_warn      = n_warn,
    n_oc        = n_oc,
    plot_data   = plot_data,
    run_name    = run_name,
    group_name  = group_name,
    title       = chart_title
  ), class = "qc_chart")
}

#' Print a Levey-Jennings QC chart analysis
#'
#' @param x A \code{qc_chart} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.qc_chart <- function(x, ...) {
  .format_num <- function(v) formatC(v, format = "f", digits = 4L)

  .format_index_list <- function(idx) {
    if (length(idx) == 0L) return("")
    idx <- sort(unique(as.integer(idx)))
    n <- length(idx)
    if (n == 1L) return(as.character(idx))
    parts <- character(0L)
    start <- idx[1L]
    end   <- idx[1L]
    for (i in 2L:n) {
      if (idx[i] == end + 1L) {
        end <- idx[i]
      } else {
        parts <- c(parts, if (start == end) as.character(start) else paste0(start, "-", end))
        start <- idx[i]
        end   <- idx[i]
      }
    }
    parts <- c(parts, if (start == end) as.character(start) else paste0(start, "-", end))
    paste(parts, collapse = ", ")
  }

  cat("\nQC Control Chart\n")
  cat(sprintf("  Value column:    %s\n", x$value))
  cat(sprintf("  Mean (target):   %s\n", .format_num(x$target_mean)))
  cat(sprintf("  SD (target):     %s\n", .format_num(x$target_sd)))
  cat(sprintf("  N (total):       %d  |  Complete: %d  |  Missing: %d\n",
              x$n_total, x$n_complete, x$n_miss))
  cat(sprintf("  Rules applied:   %s\n\n", paste(x$rule_vec, collapse = ", ")))

  if (x$n_violations == 0L) {
    cat("  Violations: none\n")
  } else {
    cat("  Violations:\n")
    cat(" ", paste(rep("-", 50), collapse = ""), "\n")
    for (rule_nm in x$rule_vec) {
      pts <- x$violations[[rule_nm]]
      if (length(pts) > 0) {
        pt_str <- .format_index_list(pts)
        cat(sprintf("    %-6s: points %s\n", rule_nm, pt_str))
      }
    }
    cat(sprintf("  -> %d violation(s): %d warning(s), %d out of control\n",
                x$n_violations, x$n_warn, x$n_oc))
  }
  cat("\n")
  invisible(x)
}


#' Plot a Levey-Jennings QC chart
#'
#' @param x A \code{qc_chart} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly; the chart is drawn as a side
#'   effect.
#' @importFrom ggplot2 ggplot aes geom_hline geom_line geom_point labs
#'   theme_bw theme element_text scale_color_manual facet_wrap element_rect
#' @importFrom ggrepel geom_text_repel
#' @examples
#' set.seed(42)
#' d <- data.frame(value = rnorm(30, 100, 5), run = 1:30)
#' plot(qc_chart(d, "value", run = "run"))
#' @export
plot.qc_chart <- function(x, ...) {
  sev_colors <- c("warning" = "#E67E22",
                  "out_of_control" = "#C0392B")

  pd <- x$plot_data
  tm <- x$target_mean
  ts <- x$target_sd

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$run, y = .data$value)) +
    ggplot2::geom_hline(yintercept = tm,            color = "gray30", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = tm + ts, color = "steelblue", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = tm - ts, color = "steelblue", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = tm + 2 * ts, color = "orange", linewidth = 0.5, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = tm - 2 * ts, color = "orange", linewidth = 0.5, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = tm + 3 * ts, color = "firebrick", linewidth = 0.6, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = tm - 3 * ts, color = "firebrick", linewidth = 0.6, linetype = "dashed") +
    ggplot2::geom_line(color = "gray50", linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_point(color = "gray30", size = 2, na.rm = TRUE) +
    ggplot2::labs(
      x     = x$run_name,
      y     = x$value,
      title = x$title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )

  # Highlight violation points
  if (x$n_violations > 0) {
    vio_data <- pd[!is.na(pd$severity), , drop = FALSE]
    if (nrow(vio_data) > 0) {
      p <- p +
        ggplot2::geom_point(data = vio_data,
                            ggplot2::aes(x = .data$run, y = .data$value,
                                         color = .data$severity),
                            size = 3, inherit.aes = FALSE) +
        ggplot2::scale_color_manual(
          values = sev_colors,
          labels = c("warning" = "Warning (1-2s)", "out_of_control" = "Out of Control (1-3s)"),
          name   = NULL,
          drop   = FALSE
        )

      has_label <- vio_data[!is.na(vio_data$vio_label), , drop = FALSE]
      if (nrow(has_label) > 0) {
        p <- p + ggrepel::geom_text_repel(
          data = has_label,
          ggplot2::aes(x = .data$run, y = .data$value,
                       label = .data$vio_label, color = .data$severity),
          size = 3, inherit.aes = FALSE,
          show.legend = FALSE,
          min.segment.length = 0.2, max.overlaps = 20
        )
      }
    }
  }

  # Facet wrap
  if (!is.null(x$group_name)) {
    p <- p + ggplot2::facet_wrap(~ group, scales = "free_x") +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey90"))
  }

  print(p)
  invisible(x)
}


#' Youden plot
#'
#' Draws a Youden plot comparing two measurement systems.
#' If target means and SDs are not provided, they are estimated
#' from the data.
#'
#' @param data    A data frame.
#' @param sample1 Column name for sample 1.
#' @param sample2 Column name for sample 2.
#' @param mean1   Target mean for sample 1; NULL to estimate.
#' @param mean2   Target mean for sample 2; NULL to estimate.
#' @param sd1     Target SD for sample 1; NULL to estimate.
#' @param sd2     Target SD for sample 2; NULL to estimate.
#'
#' @return An object of class \code{"youden_plot"} with summary statistics
#'   and plot data. Use \code{print()} to show summary and \code{plot()}
#'   to draw the Youden plot.
#'
#' @examples
#'   set.seed(42)
#'   df <- data.frame(
#'     m1 = rnorm(30, 100, 5),
#'     m2 = rnorm(30,  98, 5)
#'   )
#'   res <- youden_plot(df, "m1", "m2", mean1 = 100, mean2 = 100, sd1 = 5, sd2 = 5)
#'   print(res)
#'   plot(res)
#' @export
youden_plot <- function(data, sample1, sample2,
                        mean1 = NULL, mean2 = NULL,
                        sd1 = NULL, sd2 = NULL) {

  # Helper functions

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Input validation
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)

  for (nm in c(sample1, sample2)) {
    if (!is.character(nm) || length(nm) != 1L || !nm %in% names(data))
      stop("'", nm, "' must be a single column name in 'data'.", call. = FALSE)
    if (!is.numeric(data[[nm]]))
      stop("Column '", nm, "' must be numeric.", call. = FALSE)
  }

  # Extract data
  s1 <- data[[sample1]]
  s2 <- data[[sample2]]

  # Handle missing values
  ok <- !is.na(s1) & !is.na(s2)
  n_total <- nrow(data)
  n_complete <- sum(ok)
  n_miss <- n_total - n_complete

  if (n_complete == 0L)
    stop("No complete cases.", call. = FALSE)

  s1_clean <- s1[ok]
  s2_clean <- s2[ok]

  # Statistics
  target_mean1 <- mean1 %||% mean(s1_clean, na.rm = TRUE)
  target_mean2 <- mean2 %||% mean(s2_clean, na.rm = TRUE)
  target_sd1   <- sd1   %||% stats::sd(s1_clean, na.rm = TRUE)
  target_sd2   <- sd2   %||% stats::sd(s2_clean, na.rm = TRUE)
  cor_val      <- stats::cor(s1_clean, s2_clean, use = "complete.obs")

  # Check whether points are within the ±2s rectangle
  z1 <- (s1_clean - target_mean1) / target_sd1
  z2 <- (s2_clean - target_mean2) / target_sd2
  within_2s <- abs(z1) <= 2 & abs(z2) <= 2
  outside_idx <- which(!within_2s)
  n_outside   <- length(outside_idx)

  # Build plot_data
  plot_data <- data.frame(
    x         = s1_clean,
    y         = s2_clean,
    within_2s = within_2s,
    stringsAsFactors = FALSE
  )

  xlim <- range(s1_clean, target_mean1 + c(-3, 3) * target_sd1)
  ylim <- range(s2_clean, target_mean2 + c(-3, 3) * target_sd2)

  plot_title <- sprintf("Youden Plot (n = %d)", n_complete)

  outside_rows <- which(ok)[outside_idx]

  structure(list(
    sample1     = sample1,
    sample2     = sample2,
    mean1       = target_mean1,
    mean2       = target_mean2,
    sd1         = target_sd1,
    sd2         = target_sd2,
    correlation = cor_val,
    n_complete  = n_complete,
    n_miss      = n_miss,
    n_total     = n_total,
    n_outside   = n_outside,
    outside_rows = outside_rows,
    plot_data   = plot_data,
    xlim       = xlim,
    ylim       = ylim,
    title      = plot_title
  ), class = "youden_plot")
}

#'
#' Print a Youden plot analysis
#'
#' @param x A \code{youden_plot} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly.
#' @export
print.youden_plot <- function(x, ...) {
  .format_num <- function(v) formatC(v, format = "f", digits = 4L)

  cat("\nYouden Plot\n")
  cat(sprintf("  Sample 1 (%s):\n", x$sample1))
  cat(sprintf("    Mean:   %s\n", .format_num(x$mean1)))
  cat(sprintf("    SD:     %s\n", .format_num(x$sd1)))
  cat(sprintf("  Sample 2 (%s):\n", x$sample2))
  cat(sprintf("    Mean:   %s\n", .format_num(x$mean2)))
  cat(sprintf("    SD:     %s\n", .format_num(x$sd2)))
  cat(sprintf("  Correlation: %s\n", .format_num(x$correlation)))
  cat(sprintf("  N (total):   %d  |  Complete: %d  |  Missing: %d\n",
              x$n_total, x$n_complete, x$n_miss))

  if (x$n_outside > 0L) {
    rows <- paste(x$outside_rows, collapse = ", ")
    cat(sprintf("  Points outside +/-2SD: %d  (rows: %s)\n", x$n_outside, rows))
  } else {
    cat("  Points outside +/-2SD: 0\n")
  }

  cat("\n")
  invisible(x)
}

#' Plot a Youden analysis
#'
#' @param x A \code{youden_plot} object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged \code{x}, invisibly; the plot is drawn as a side
#'   effect.
#' @importFrom ggplot2 ggplot aes geom_abline annotate geom_vline
#'   geom_hline geom_point scale_color_manual labs coord_equal
#'   theme_bw theme element_text
#' @examples
#' set.seed(42)
#' d <- data.frame(m1 = rnorm(30, 100, 5), m2 = rnorm(30, 98, 5))
#' plot(youden_plot(d, "m1", "m2"))
#' @export
plot.youden_plot <- function(x, ...) {
  p <- ggplot2::ggplot(x$plot_data, ggplot2::aes(x = .data$x, y = .data$y)) +
    # Identity line y = x
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         color = "gray60", linewidth = 0.5, linetype = "dotted") +
    # ±2s rectangle
    ggplot2::annotate("rect",
             xmin = x$mean1 - 2 * x$sd1,
             xmax = x$mean1 + 2 * x$sd1,
             ymin = x$mean2 - 2 * x$sd2,
             ymax = x$mean2 + 2 * x$sd2,
             fill = "orange", alpha = 0.06, color = "orange", linewidth = 0.5, linetype = "dashed") +
    # ±1s rectangle
    ggplot2::annotate("rect",
             xmin = x$mean1 - x$sd1,
             xmax = x$mean1 + x$sd1,
             ymin = x$mean2 - x$sd2,
             ymax = x$mean2 + x$sd2,
             fill = "steelblue", alpha = 0.06, color = "steelblue", linewidth = 0.4, linetype = "dashed") +
    # Center crosshairs
    ggplot2::geom_vline(xintercept = x$mean1, color = "gray30", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = x$mean2, color = "gray30", linewidth = 0.6) +
    # Scatter points
    ggplot2::geom_point(
      ggplot2::aes(color = .data$within_2s),
      size = 2.5, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "gray30", "FALSE" = "firebrick"),
      labels = c("TRUE" = "Within +/-2s", "FALSE" = "Outside +/-2s"),
      name   = NULL
    ) +
    ggplot2::labs(
      x     = x$sample1,
      y     = x$sample2,
      title = x$title
    ) +

    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )

  print(p)
  invisible(x)
}
