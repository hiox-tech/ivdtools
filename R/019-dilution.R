# CLSI EP34 dilution, spiking, and high-dose hook-effect calculations.
#
# These functions calculate estimates and uncertainty. They intentionally do
# not apply acceptance limits or make establishment, validation, verification,
# or reportable-interval decisions.

utils::globalVariables(c(".sample", ".series", ".concentration", ".mean",
                         ".target", ".recovery", ".specimen_content"))

.ep34_col <- function(data, column, argument, numeric = FALSE,
                      allow_null = FALSE) {
  if (is.null(column) && allow_null) return(NULL)
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column) || !column %in% names(data)) {
    stop("`", argument, "` must be a single column name in `data`.",
         call. = FALSE)
  }
  if (numeric && !is.numeric(data[[column]])) {
    stop("Column '", column, "' supplied to `", argument,
         "` must be numeric.", call. = FALSE)
  }
  column
}

.ep34_conf_level <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop("`conf.level` must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(x)
}

.ep34_stats <- function(x) {
  n <- length(x)
  average <- mean(x)
  standard_deviation <- if (n > 1L) stats::sd(x) else NA_real_
  data.frame(
    n = n,
    mean = average,
    sd = standard_deviation,
    cv = if (n > 1L && abs(average) > .Machine$double.eps^0.5)
      100 * standard_deviation / abs(average) else NA_real_,
    row.names = NULL
  )
}

.ep34_mean_ci <- function(x, conf.level) {
  x <- x[is.finite(x)]
  n <- length(x)
  average <- if (n) mean(x) else NA_real_
  standard_deviation <- if (n > 1L) stats::sd(x) else NA_real_
  standard_error <- if (n > 1L) standard_deviation / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt((1 + conf.level) / 2, df = n - 1L) else
    NA_real_
  data.frame(
    n = n,
    mean = average,
    sd = standard_deviation,
    se = standard_error,
    df = if (n > 1L) n - 1L else NA_integer_,
    conf_low = average - critical * standard_error,
    conf_high = average + critical * standard_error,
    row.names = NULL
  )
}

.ep34_single_value <- function(x, label, sample) {
  values <- unique(x[is.finite(x)])
  if (length(values) != 1L) {
    stop("`", label, "` must have exactly one finite value for sample '",
         sample, "'.", call. = FALSE)
  }
  values
}

#' Calculate dilution recovery
#'
#' Summarizes replicate measurements and calculates target concentration,
#' dilution-corrected concentration, and recovery for each specimen and
#' dilution factor. Across-specimen mean recovery and its Student t confidence
#' interval are calculated for every dilution factor. No acceptance criterion
#' or suitability decision is applied.
#'
#' @param data Data frame in long format, with one row per measurement.
#' @param sample Sample-identification column.
#' @param result Numeric measurement-result column.
#' @param dilution_factor Numeric dilution-factor column. Undiluted material is
#'   represented by 1, a 1/2 dilution by 2, and a 1/20 dilution by 20.
#' @param assigned Optional column containing one assigned original-sample
#'   concentration per sample. When omitted, the mean undiluted result is used.
#' @param diluent_value Measurand concentration in the diluent. The default is 0.
#' @param conf.level Confidence level for the across-sample mean recovery.
#'
#' @return An object of class `dilution_recovery` containing `results`,
#'   `summary`, the retained data, and excluded row indices.
#' @references CLSI. Establishing and Verifying an Extended Measuring Interval
#'   Through Specimen Dilution and Spiking. EP34, 1st ed. 2018.
#' @examples
#' dilution_data <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 6),
#'   factor = rep(rep(c(1, 2, 4), each = 2), 2),
#'   result = c(99, 101, 49, 51, 24, 26,
#'              199, 201, 99, 101, 49, 51)
#' )
#' dilution <- dilution_recovery(
#'   dilution_data, "sample", "result", "factor"
#' )
#' dilution$summary
#' @export
dilution_recovery <- function(data, sample, result, dilution_factor,
                              assigned = NULL, diluent_value = 0,
                              conf.level = 0.95) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  sample <- .ep34_col(data, sample, "sample")
  result <- .ep34_col(data, result, "result", numeric = TRUE)
  dilution_factor <- .ep34_col(data, dilution_factor, "dilution_factor",
                               numeric = TRUE)
  assigned <- .ep34_col(data, assigned, "assigned", numeric = TRUE,
                        allow_null = TRUE)
  conf.level <- .ep34_conf_level(conf.level)
  if (!is.numeric(diluent_value) || length(diluent_value) != 1L ||
      is.na(diluent_value) || !is.finite(diluent_value)) {
    stop("`diluent_value` must be one finite numeric value.", call. = FALSE)
  }

  keep <- !is.na(data[[sample]]) & is.finite(data[[result]]) &
    is.finite(data[[dilution_factor]]) & data[[dilution_factor]] >= 1
  if (!is.null(assigned)) keep <- keep & is.finite(data[[assigned]])
  excluded <- which(!keep)
  d <- data[keep, , drop = FALSE]
  if (!nrow(d)) stop("No usable observations remain.", call. = FALSE)
  d$.sample <- as.character(d[[sample]])
  d$.result <- as.numeric(d[[result]])
  d$.factor <- as.numeric(d[[dilution_factor]])
  d$.assigned <- if (is.null(assigned)) NA_real_ else as.numeric(d[[assigned]])

  sample_levels <- unique(d$.sample)
  assigned_values <- setNames(numeric(length(sample_levels)), sample_levels)
  for (id in sample_levels) {
    z <- d[d$.sample == id, , drop = FALSE]
    if (is.null(assigned)) {
      neat <- z$.result[z$.factor == 1]
      if (!length(neat)) {
        stop("Sample '", id, "' has no undiluted result; supply `assigned`.",
             call. = FALSE)
      }
      assigned_values[id] <- mean(neat)
    } else {
      assigned_values[id] <- .ep34_single_value(z$.assigned, "assigned", id)
    }
  }

  groups <- split(seq_len(nrow(d)), interaction(d$.sample, d$.factor,
                                                 drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(index) {
    z <- d[index, , drop = FALSE]
    descriptive <- .ep34_stats(z$.result)
    id <- z$.sample[1L]
    factor <- z$.factor[1L]
    fraction <- 1 / factor
    original <- unname(assigned_values[id])
    target <- original * fraction + diluent_value * (1 - fraction)
    corrected <- (descriptive$mean - diluent_value * (1 - fraction)) /
      fraction
    data.frame(
      sample = id,
      dilution_factor = factor,
      specimen_content = 100 * fraction,
      n = descriptive$n,
      mean = descriptive$mean,
      sd = descriptive$sd,
      cv = descriptive$cv,
      assigned = original,
      target = target,
      corrected_concentration = corrected,
      recovery = if (abs(target) > .Machine$double.eps^0.5)
        100 * descriptive$mean / target else NA_real_,
      row.names = NULL,
      check.names = FALSE
    )
  })
  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  results <- results[order(results$dilution_factor, results$sample), , drop = FALSE]

  factors <- sort(unique(results$dilution_factor))
  summaries <- lapply(factors, function(factor) {
    z <- results[results$dilution_factor == factor, , drop = FALSE]
    ci <- .ep34_mean_ci(z$recovery, conf.level)
    data.frame(
      dilution_factor = factor,
      specimen_content = 100 / factor,
      n_samples = ci$n,
      mean_recovery = ci$mean,
      sd_recovery = ci$sd,
      se_recovery = ci$se,
      df = ci$df,
      conf_low = ci$conf_low,
      conf_high = ci$conf_high,
      row.names = NULL
    )
  })
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL

  structure(list(
    call = match.call(),
    columns = list(sample = sample, result = result,
                   dilution_factor = dilution_factor, assigned = assigned),
    diluent_value = as.numeric(diluent_value),
    conf.level = conf.level,
    results = results,
    summary = summary,
    data = d,
    excluded = excluded
  ), class = "dilution_recovery")
}

#' Print dilution-recovery results
#'
#' @param x A `dilution_recovery` object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.dilution_recovery <- function(x, ...) {
  cat("\nDilution recovery\n")
  cat(sprintf("  Samples: %d\n", length(unique(x$results$sample))))
  cat(sprintf("  Dilution factors: %s\n",
              paste(format(x$summary$dilution_factor), collapse = ", ")))
  cat(sprintf("  Confidence level: %.1f%%\n\n", 100 * x$conf.level))
  .print_kv_rows(x$summary, list(
    "Recovery summary" = c(
      "Specimen content" = "specimen_content", "Samples" = "n_samples",
      "Mean recovery" = "mean_recovery", "SD recovery" = "sd_recovery"),
    "Confidence interval" = c(
      "SE recovery" = "se_recovery", "Degrees of freedom" = "df",
      "CI lower" = "conf_low", "CI upper" = "conf_high")
  ), row_title = c("Dilution factor" = "dilution_factor"))
  cat("\n  Sample results\n\n")
  .print_kv_rows(x$results, list(
    "Measurement summary" = c(
      "Replicates" = "n", "Mean result" = "mean", "SD result" = "sd",
      "CV result" = "cv"),
    "Recovery result" = c(
      "Specimen content" = "specimen_content", "Assigned value" = "assigned",
      "Target concentration" = "target",
      "Corrected concentration" = "corrected_concentration",
      "Recovery" = "recovery")
  ), row_title = c("Sample" = "sample",
                   "Dilution factor" = "dilution_factor"))
  invisible(x)
}

#' Plot dilution-recovery results
#'
#' @param x A `dilution_recovery` object.
#' @param type `"recovery"`, `"summary"`, or `"observed"`.
#' @param ... Reserved arguments.
#' @return A ggplot object, invisibly.
#' @examples
#' dilution_data <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 6),
#'   factor = rep(rep(c(1, 2, 4), each = 2), 2),
#'   result = c(99, 101, 49, 51, 24, 26,
#'              199, 201, 99, 101, 49, 51)
#' )
#' dilution <- dilution_recovery(
#'   dilution_data, "sample", "result", "factor"
#' )
#' plot(dilution)
#' @export
plot.dilution_recovery <- function(x,
                                   type = c("recovery", "summary", "observed"),
                                   ...) {
  type <- match.arg(type)
  if (type == "recovery") {
    p <- ggplot2::ggplot(
      x$results,
      ggplot2::aes(x = .data$specimen_content, y = .data$recovery,
                   group = .data$sample, color = .data$sample)
    ) +
      ggplot2::geom_hline(yintercept = 100, linetype = 2, color = "grey40") +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::labs(x = "Specimen content (%)", y = "Recovery (%)",
                    color = "Sample", title = "Dilution recovery by sample")
  } else if (type == "summary") {
    p <- ggplot2::ggplot(
      x$summary,
      ggplot2::aes(x = .data$specimen_content, y = .data$mean_recovery)
    ) +
      ggplot2::geom_hline(yintercept = 100, linetype = 2, color = "grey40") +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
        width = 1.5, na.rm = TRUE
      ) +
      ggplot2::geom_line(color = "steelblue") +
      ggplot2::geom_point(color = "steelblue", size = 2.7) +
      ggplot2::labs(x = "Specimen content (%)", y = "Mean recovery (%)",
                    title = "Mean dilution recovery with confidence intervals")
  } else {
    p <- ggplot2::ggplot(
      x$results,
      ggplot2::aes(x = .data$target, y = .data$mean, color = .data$sample)
    ) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2,
                           color = "grey40") +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::labs(x = "Target concentration", y = "Observed mean",
                    color = "Sample", title = "Observed and target concentrations")
  }
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)
  invisible(p)
}

#' Calculate the concentration of a spiked sample
#'
#' Applies the mass-balance equation for mixing a sample and a spiking
#' material. Arguments are vectorized using ordinary R recycling.
#'
#' @param sample_concentration Concentration in the sample before spiking.
#' @param sample_volume Volume of sample used.
#' @param stock_concentration Concentration of the spiking material.
#' @param spike_volume Volume of spiking material added.
#'
#' @return A data frame containing total volume, spike fraction, concentration
#'   added by the spike, and final calculated concentration.
#' @references CLSI. Establishing and Verifying an Extended Measuring Interval
#'   Through Specimen Dilution and Spiking. EP34, 1st ed. 2018.
#' @examples
#' spike_concentration(
#'   sample_concentration = 10, sample_volume = 0.95,
#'   stock_concentration = 210, spike_volume = 0.05
#' )
#' @export
spike_concentration <- function(sample_concentration, sample_volume,
                                stock_concentration, spike_volume) {
  args <- list(sample_concentration = sample_concentration,
               sample_volume = sample_volume,
               stock_concentration = stock_concentration,
               spike_volume = spike_volume)
  if (any(!vapply(args, is.numeric, logical(1))) ||
      any(vapply(args, function(z) anyNA(z) || any(!is.finite(z)), logical(1)))) {
    stop("All concentration and volume inputs must be finite numeric values.",
         call. = FALSE)
  }
  n <- max(vapply(args, length, integer(1)))
  lengths <- vapply(args, length, integer(1))
  if (any(!lengths %in% c(1L, n))) {
    stop("Inputs must have length 1 or a common maximum length.", call. = FALSE)
  }
  args <- lapply(args, rep_len, length.out = n)
  if (any(args$sample_concentration < 0) ||
      any(args$stock_concentration < 0) ||
      any(args$sample_volume <= 0) || any(args$spike_volume <= 0)) {
    stop("Concentrations must be nonnegative and volumes must be positive.",
         call. = FALSE)
  }
  total <- args$sample_volume + args$spike_volume
  fraction <- args$spike_volume / total
  data.frame(
    sample_concentration = args$sample_concentration,
    sample_volume = args$sample_volume,
    stock_concentration = args$stock_concentration,
    spike_volume = args$spike_volume,
    total_volume = total,
    spike_fraction = fraction,
    added_concentration = args$stock_concentration * fraction,
    final_concentration =
      (args$sample_concentration * args$sample_volume +
       args$stock_concentration * args$spike_volume) / total,
    row.names = NULL
  ) |> structure(class = c("spike_concentration", "data.frame"))
}

#' Print spike-concentration calculations
#'
#' @param x A `spike_concentration` data frame.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.spike_concentration <- function(x, ...) {
  cat("\nSpike concentration\n")
  .print_kv_rows(x, list(
    "Input concentrations" = c("Stock concentration" = "stock_concentration"),
    "Calculated concentrations" = c(
      "Added concentration" = "added_concentration",
      "Final concentration" = "final_concentration"),
    "Volume and fraction" = c(
      "Sample volume" = "sample_volume", "Spike volume" = "spike_volume",
      "Total volume" = "total_volume", "Spike fraction" = "spike_fraction")
  ), row_title = c("Sample concentration" = "sample_concentration"))
  invisible(x)
}

.ep34_numeric_spec <- function(data, value, argument) {
  if (is.character(value)) {
    column <- .ep34_col(data, value, argument, numeric = TRUE)
    out <- as.numeric(data[[column]])
  } else if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
    column <- NULL
    out <- rep(as.numeric(value), nrow(data))
  } else {
    stop("`", argument, "` must be one finite number or a numeric column name.",
         call. = FALSE)
  }
  list(value = out, column = column)
}

#' Calculate spike recovery against a solvent control
#'
#' Calculates the difference between analyte-spiked and equal-volume
#' solvent-control means and expresses that difference as a percentage of the
#' known added concentration. No acceptance criterion is applied.
#'
#' @param data Data frame in long format.
#' @param sample Sample-identification column.
#' @param result Numeric measurement-result column.
#' @param condition Column distinguishing analyte-spiked and solvent-control
#'   preparations.
#' @param analyte_level Value identifying the analyte-spiked preparation.
#' @param solvent_level Value identifying the solvent-control preparation.
#' @param added_concentration A finite scalar or numeric column containing the
#'   theoretical concentration added to the final mixture. When omitted,
#'   `sample_volume`, `stock_concentration`, and `spike_volume` are required.
#' @param sample_volume,stock_concentration,spike_volume Finite scalars or
#'   numeric columns used to calculate the added concentration.
#' @param conf.level Confidence level for the across-sample mean recovery.
#'
#' @return An object of class `spike_recovery`.
#' @references CLSI. Establishing and Verifying an Extended Measuring Interval
#'   Through Specimen Dilution and Spiking. EP34, 1st ed. 2018.
#' @examples
#' spike_data <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 4),
#'   condition = rep(rep(c("spiked", "control"), each = 2), 2),
#'   result = c(20.1, 19.9, 10.0, 10.2, 30.2, 29.8, 20.1, 19.9)
#' )
#' recovery <- spike_recovery(
#'   spike_data, "sample", "result", "condition",
#'   analyte_level = "spiked", solvent_level = "control",
#'   added_concentration = 10
#' )
#' recovery$results
#' @export
spike_recovery <- function(data, sample, result, condition,
                           analyte_level, solvent_level,
                           added_concentration = NULL,
                           sample_volume = NULL, stock_concentration = NULL,
                           spike_volume = NULL, conf.level = 0.95) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  sample <- .ep34_col(data, sample, "sample")
  result <- .ep34_col(data, result, "result", numeric = TRUE)
  condition <- .ep34_col(data, condition, "condition")
  conf.level <- .ep34_conf_level(conf.level)
  if (length(analyte_level) != 1L || is.na(analyte_level) ||
      length(solvent_level) != 1L || is.na(solvent_level) ||
      identical(analyte_level, solvent_level)) {
    stop("`analyte_level` and `solvent_level` must be distinct single values.",
         call. = FALSE)
  }

  if (!is.null(added_concentration)) {
    added <- .ep34_numeric_spec(data, added_concentration,
                                "added_concentration")
    added_source <- added$column %||% "scalar"
  } else {
    if (is.null(sample_volume) || is.null(stock_concentration) ||
        is.null(spike_volume)) {
      stop("Supply `added_concentration`, or all of `sample_volume`, ",
           "`stock_concentration`, and `spike_volume`.", call. = FALSE)
    }
    sv <- .ep34_numeric_spec(data, sample_volume, "sample_volume")
    sc <- .ep34_numeric_spec(data, stock_concentration, "stock_concentration")
    pv <- .ep34_numeric_spec(data, spike_volume, "spike_volume")
    if (any(sv$value <= 0) || any(sc$value < 0) || any(pv$value <= 0)) {
      stop("Volumes must be positive and stock concentration nonnegative.",
           call. = FALSE)
    }
    added <- list(value = sc$value * pv$value / (sv$value + pv$value),
                  column = NULL)
    added_source <- "calculated from concentration and volumes"
  }
  if (any(!is.finite(added$value)) || any(added$value <= 0)) {
    stop("Added concentrations must be positive and finite.", call. = FALSE)
  }

  selected <- data[[condition]] %in% c(analyte_level, solvent_level)
  keep <- selected & !is.na(data[[sample]]) & is.finite(data[[result]]) &
    is.finite(added$value)
  excluded <- which(!keep)
  d <- data[keep, , drop = FALSE]
  if (!nrow(d)) stop("No usable analyte-spike or solvent-control rows remain.",
                     call. = FALSE)
  d$.sample <- as.character(d[[sample]])
  d$.result <- as.numeric(d[[result]])
  d$.condition <- ifelse(d[[condition]] == analyte_level,
                         "analyte", "solvent")
  d$.added <- added$value[keep]

  rows <- lapply(unique(d$.sample), function(id) {
    z <- d[d$.sample == id, , drop = FALSE]
    analyte <- z$.result[z$.condition == "analyte"]
    solvent <- z$.result[z$.condition == "solvent"]
    if (!length(analyte) || !length(solvent)) {
      stop("Sample '", id, "' must contain both analyte-spiked and ",
           "solvent-control results.", call. = FALSE)
    }
    added_value <- .ep34_single_value(z$.added, "added_concentration", id)
    analyte_stats <- .ep34_stats(analyte)
    solvent_stats <- .ep34_stats(solvent)
    difference <- analyte_stats$mean - solvent_stats$mean
    data.frame(
      sample = id,
      analyte_n = analyte_stats$n,
      analyte_mean = analyte_stats$mean,
      analyte_sd = analyte_stats$sd,
      solvent_n = solvent_stats$n,
      solvent_mean = solvent_stats$mean,
      solvent_sd = solvent_stats$sd,
      observed_increase = difference,
      added_concentration = added_value,
      recovery = 100 * difference / added_value,
      row.names = NULL
    )
  })
  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  overall <- .ep34_mean_ci(results$recovery, conf.level)
  names(overall) <- c("n_samples", "mean_recovery", "sd_recovery",
                      "se_recovery", "df", "conf_low", "conf_high")

  structure(list(
    call = match.call(),
    columns = list(sample = sample, result = result, condition = condition),
    levels = list(analyte = analyte_level, solvent = solvent_level),
    added_source = added_source,
    conf.level = conf.level,
    results = results,
    summary = overall,
    data = d,
    excluded = excluded
  ), class = "spike_recovery")
}

#' Print spike-recovery results
#'
#' @param x A `spike_recovery` object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.spike_recovery <- function(x, ...) {
  cat("\nSpike recovery\n")
  cat(sprintf("  Samples: %d\n", nrow(x$results)))
  cat(sprintf("  Confidence level: %.1f%%\n\n", 100 * x$conf.level))
  if (nrow(x$summary) == 1L) {
    .print_kv_sections(x$summary, list(
      "Recovery estimate" = c(
        "Samples" = "n_samples", "Mean recovery" = "mean_recovery",
        "SD of recovery" = "sd_recovery", "SE of recovery" = "se_recovery"),
      "Confidence interval" = c(
        "Degrees of freedom" = "df", "CI lower" = "conf_low",
        "CI upper" = "conf_high")
    ))
  } else {
    .print_kv_rows(x$summary, list(
      "Recovery estimate" = c(
        "Samples" = "n_samples", "Mean recovery" = "mean_recovery",
        "SD recovery" = "sd_recovery", "SE recovery" = "se_recovery"),
      "Confidence interval" = c(
        "Samples" = "n_samples", "Degrees of freedom" = "df",
        "CI lower" = "conf_low", "CI upper" = "conf_high")
    ))
  }
  cat("\n  Sample results\n\n")
  .print_kv_rows(x$results, list(
    "Analyte-spiked results" = c(
      "Analyte replicates" = "analyte_n", "Analyte mean" = "analyte_mean",
      "Analyte SD" = "analyte_sd"),
    "Solvent-control results" = c(
      "Solvent replicates" = "solvent_n", "Solvent mean" = "solvent_mean",
      "Solvent SD" = "solvent_sd"),
    "Recovery result" = c(
      "Observed increase" = "observed_increase",
      "Added concentration" = "added_concentration", "Recovery" = "recovery")
  ), row_title = c("Sample" = "sample"))
  invisible(x)
}

#' Plot spike-recovery results
#'
#' @param x A `spike_recovery` object.
#' @param ... Reserved arguments.
#' @return A ggplot object, invisibly.
#' @examples
#' spike_data <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 4),
#'   condition = rep(rep(c("spiked", "control"), each = 2), 2),
#'   result = c(20.1, 19.9, 10.0, 10.2, 30.2, 29.8, 20.1, 19.9)
#' )
#' recovery <- spike_recovery(
#'   spike_data, "sample", "result", "condition",
#'   analyte_level = "spiked", solvent_level = "control",
#'   added_concentration = 10
#' )
#' plot(recovery)
#' @export
plot.spike_recovery <- function(x, ...) {
  average <- x$summary$mean_recovery
  p <- ggplot2::ggplot(
    x$results,
    ggplot2::aes(x = stats::reorder(.data$sample, .data$recovery),
                 y = .data$recovery)
  ) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = average, color = "steelblue") +
    ggplot2::geom_point(size = 2.7, color = "steelblue") +
    ggplot2::labs(x = "Sample", y = "Spike recovery (%)",
                  title = "Spike recovery by sample") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)
  invisible(p)
}

#' Summarize a high-dose hook-effect experiment
#'
#' Summarizes raw response by known concentration and identifies the first
#' post-peak measured concentration whose mean response is at or below the raw
#' response associated with the ULoQ. The reported `hook_concentration` is the
#' immediately preceding measured concentration; no interpolation is used. No
#' hook-effect or extended-interval acceptability decision is made.
#'
#' @param data Data frame containing known concentrations and raw responses.
#' @param concentration Numeric known-concentration column.
#' @param response Numeric raw-response column.
#' @param series Optional column identifying independent hook experiments.
#' @param uloq Finite concentration representing the upper limit of quantitation.
#' @param uloq_response Optional raw response at the ULoQ. It may be one finite
#'   value for all series or one value per series. When omitted, each series
#'   must contain measurements at `concentration == uloq`.
#' @param conf.level Confidence level for mean raw response at each concentration.
#'
#' @return An object of class `hook_effect` containing concentration-level
#'   summaries and one threshold result per series.
#' @references CLSI. Establishing and Verifying an Extended Measuring Interval
#'   Through Specimen Dilution and Spiking. EP34, 1st ed. 2018.
#' @examples
#' hook_data <- data.frame(
#'   concentration = rep(c(100, 200, 400, 800), each = 2),
#'   response = c(48, 52, 98, 102, 148, 152, 88, 92)
#' )
#' hook <- hook_effect(
#'   hook_data, "concentration", "response", uloq = 200
#' )
#' hook$results
#' @export
hook_effect <- function(data, concentration, response, series = NULL,
                        uloq, uloq_response = NULL, conf.level = 0.95) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  concentration <- .ep34_col(data, concentration, "concentration", numeric = TRUE)
  response <- .ep34_col(data, response, "response", numeric = TRUE)
  series <- .ep34_col(data, series, "series", allow_null = TRUE)
  conf.level <- .ep34_conf_level(conf.level)
  if (!is.numeric(uloq) || length(uloq) != 1L || is.na(uloq) ||
      !is.finite(uloq)) {
    stop("`uloq` must be one finite numeric concentration.", call. = FALSE)
  }

  keep <- is.finite(data[[concentration]]) & is.finite(data[[response]])
  if (!is.null(series)) keep <- keep & !is.na(data[[series]])
  excluded <- which(!keep)
  d <- data[keep, , drop = FALSE]
  if (!nrow(d)) stop("No usable observations remain.", call. = FALSE)
  d$.series <- if (is.null(series)) "all" else as.character(d[[series]])
  d$.concentration <- as.numeric(d[[concentration]])
  d$.response <- as.numeric(d[[response]])
  series_levels <- unique(d$.series)

  groups <- split(seq_len(nrow(d)),
                  interaction(d$.series, d$.concentration,
                              drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(index) {
    z <- d[index, , drop = FALSE]
    descriptive <- .ep34_stats(z$.response)
    ci <- .ep34_mean_ci(z$.response, conf.level)
    data.frame(
      series = z$.series[1L],
      concentration = z$.concentration[1L],
      n = descriptive$n,
      mean = descriptive$mean,
      sd = descriptive$sd,
      cv = descriptive$cv,
      se = ci$se,
      conf_low = ci$conf_low,
      conf_high = ci$conf_high,
      row.names = NULL
    )
  })
  level_summary <- do.call(rbind, rows)
  rownames(level_summary) <- NULL
  level_summary <- level_summary[
    order(match(level_summary$series, series_levels),
          level_summary$concentration), , drop = FALSE]

  if (is.null(uloq_response)) {
    thresholds <- setNames(numeric(length(series_levels)), series_levels)
    tolerance <- .Machine$double.eps^0.5 * max(1, abs(uloq))
    for (id in series_levels) {
      z <- level_summary[level_summary$series == id &
        abs(level_summary$concentration - uloq) <= tolerance, , drop = FALSE]
      if (!nrow(z)) {
        stop("Series '", id, "' has no response measured at `uloq`; supply ",
             "`uloq_response`.", call. = FALSE)
      }
      thresholds[id] <- mean(z$mean)
    }
    threshold_source <- "measured at ULoQ"
  } else {
    if (!is.numeric(uloq_response) || anyNA(uloq_response) ||
        any(!is.finite(uloq_response)) ||
        !length(uloq_response) %in% c(1L, length(series_levels))) {
      stop("`uloq_response` must have length 1 or the number of series and ",
           "contain finite numeric values.", call. = FALSE)
    }
    thresholds <- setNames(rep_len(as.numeric(uloq_response),
                                   length(series_levels)), series_levels)
    threshold_source <- "supplied"
  }

  assessments <- lapply(series_levels, function(id) {
    z <- level_summary[level_summary$series == id, , drop = FALSE]
    z <- z[order(z$concentration), , drop = FALSE]
    threshold <- unname(thresholds[id])
    peak_index <- which.max(z$mean)
    candidate <- which(seq_len(nrow(z)) > peak_index &
                         z$concentration > uloq & z$mean <= threshold)
    highest <- max(z$concentration)
    if (!length(candidate)) {
      return(data.frame(
        series = id, uloq = uloq, uloq_response = threshold,
        hook_concentration = NA_real_, hook_response = NA_real_,
        first_in_range_concentration = NA_real_,
        first_in_range_response = NA_real_,
        highest_tested_concentration = highest,
        hook_lower_bound = highest,
        row.names = NULL
      ))
    }
    crossing <- candidate[1L]
    previous <- crossing - 1L
    data.frame(
      series = id, uloq = uloq, uloq_response = threshold,
      hook_concentration = z$concentration[previous],
      hook_response = z$mean[previous],
      first_in_range_concentration = z$concentration[crossing],
      first_in_range_response = z$mean[crossing],
      highest_tested_concentration = highest,
      hook_lower_bound = NA_real_,
      row.names = NULL
    )
  })
  results <- do.call(rbind, assessments)
  rownames(results) <- NULL
  finite_hook <- results$hook_concentration[is.finite(results$hook_concentration)]

  structure(list(
    call = match.call(),
    columns = list(concentration = concentration, response = response,
                   series = series),
    uloq = as.numeric(uloq),
    threshold_source = threshold_source,
    conf.level = conf.level,
    level_summary = level_summary,
    results = results,
    minimum_hook_concentration = if (length(finite_hook)) min(finite_hook) else
      NA_real_,
    data = d,
    excluded = excluded
  ), class = "hook_effect")
}

#' Print high-dose hook-effect results
#'
#' @param x A `hook_effect` object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.hook_effect <- function(x, ...) {
  cat("\nHigh-dose hook-effect experiment\n")
  cat(sprintf("  ULoQ concentration: %s\n", format(x$uloq)))
  cat(sprintf("  ULoQ response: %s\n", x$threshold_source))
  cat(sprintf("  Series: %d\n\n", nrow(x$results)))
  .print_kv_rows(x$results, list(
    "Hook detection" = c(
      "ULoQ concentration" = "uloq", "ULoQ response" = "uloq_response",
      "Hook concentration" = "hook_concentration", "Hook response" = "hook_response"),
    "First in-range result" = c(
      "First in-range concentration" = "first_in_range_concentration",
      "First in-range response" = "first_in_range_response"),
    "Tested range" = c(
      "Highest tested concentration" = "highest_tested_concentration",
      "Hook lower bound" = "hook_lower_bound")
  ), row_title = c("Series" = "series"))
  cat("\n  Concentration-level results\n\n")
  .print_kv_rows(x$level_summary, list(
    "Response summary" = c(
      "Replicates" = "n", "Mean response" = "mean", "SD response" = "sd",
      "CV response" = "cv"),
    "Confidence interval" = c(
      "SE response" = "se", "CI lower" = "conf_low",
      "CI upper" = "conf_high")
  ), row_title = c("Series" = "series",
                   "Concentration" = "concentration"))
  invisible(x)
}

#' Plot a high-dose hook-effect experiment
#'
#' @param x A `hook_effect` object.
#' @param ... Reserved arguments.
#' @return A ggplot object, invisibly.
#' @examples
#' hook_data <- data.frame(
#'   series = rep("S1", 8),
#'   concentration = rep(c(10, 20, 40, 80, 160, 320, 640, 1280), each = 2),
#'   response = c(10, 11, 20, 21, 39, 40, 70, 71,
#'                95, 94, 80, 79, 55, 54, 30, 29)
#' )
#' hook <- hook_effect(hook_data, "concentration", "response",
#'                     series = "series", uloq = 320)
#' plot(hook)
#' @export
plot.hook_effect <- function(x, ...) {
  thresholds <- x$results[, c("series", "uloq_response"), drop = FALSE]
  hooks <- x$results[is.finite(x$results$hook_concentration),
                     c("series", "hook_concentration"), drop = FALSE]
  p <- ggplot2::ggplot(
    x$level_summary,
    ggplot2::aes(x = .data$concentration, y = .data$mean,
                 group = .data$series, color = .data$series)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
      width = 0, na.rm = TRUE
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = x$uloq, linetype = 2, color = "grey40") +
    ggplot2::geom_hline(
      data = thresholds,
      ggplot2::aes(yintercept = .data$uloq_response, color = .data$series),
      linetype = 3, show.legend = FALSE
    ) +
    ggplot2::labs(x = "Known concentration", y = "Raw response",
                  color = "Series", title = "High-dose hook-effect experiment")
  if (nrow(hooks)) {
    p <- p + ggplot2::geom_vline(
      data = hooks,
      ggplot2::aes(xintercept = .data$hook_concentration,
                   color = .data$series),
      show.legend = FALSE
    )
  }
  if (nrow(thresholds) > 1L) {
    p <- p + ggplot2::facet_wrap(~series, scales = "free_y")
  }
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)
  invisible(p)
}
