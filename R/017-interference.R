# CLSI EP07 interference tools.
#
# The standard paired-difference and replicate-count calculations follow
# CLSI EP07, 3rd edition (2018). The analyte-spike and direct-blank designs
# are explicit extensions for common laboratory preparation schemes.
# Functions report estimates and uncertainty only; they do not make an
# acceptability or interference conclusion.

utils::globalVariables(c(".stratum", ".condition", ".matrix", ".spike",
                         ".result", ".concentration", ".group", ".method"))

.interference_col <- function(data, column, argument, numeric = FALSE,
                              allow_null = FALSE) {
  if (is.null(column) && allow_null) return(NULL)
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !column %in% names(data)) {
    stop("`", argument, "` must be a single column name in `data`.",
         call. = FALSE)
  }
  if (numeric && !is.numeric(data[[column]])) {
    stop("Column '", column, "' supplied to `", argument,
         "` must be numeric.", call. = FALSE)
  }
  column
}

.interference_probability <- function(x, argument) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop("`", argument, "` must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(x)
}

.interference_by <- function(data, by) {
  if (is.null(by)) return(character())
  if (!is.character(by) || anyNA(by) || any(!by %in% names(data))) {
    stop("`by` must contain column names found in `data`.", call. = FALSE)
  }
  unique(by)
}

.interference_strata <- function(data, by) {
  if (!length(by)) return(rep("all", nrow(data)))
  # Build stable, readable keys without relying on factor level ordering.
  do.call(paste, c(lapply(seq_along(by), function(j) {
    paste0(by[j], "=", as.character(data[[by[j]]]))
  }), sep = " | "))
}

.interference_group_values <- function(data, by) {
  if (!length(by)) return(data.frame(.stratum = "all", stringsAsFactors = FALSE))
  out <- data[!duplicated(data$.stratum), c(".stratum", by), drop = FALSE]
  rownames(out) <- NULL
  out
}

.interference_issue <- function(stratum, code, message) {
  data.frame(stratum = as.character(stratum), code = as.character(code),
             message = as.character(message), stringsAsFactors = FALSE)
}

.interference_issues <- function(rows) {
  out <- .interference_bind(rows)
  if (nrow(out)) return(out)
  data.frame(stratum = character(), code = character(), message = character(),
             stringsAsFactors = FALSE)
}

.interference_bind <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.interference_stats <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  data.frame(n = n, mean = if (n) mean(x) else NA_real_,
             sd = if (n > 1L) stats::sd(x) else NA_real_,
             se = if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_)
}

.interference_df <- function(variances, ns) {
  terms <- variances / ns
  denominator <- sum(terms^2 / (ns - 1))
  if (!is.finite(denominator) || denominator <= 0) return(Inf)
  sum(terms)^2 / denominator
}

.interference_ci <- function(estimate, se, df, conf.level) {
  if (!is.finite(se) || se < 0 || (!is.finite(df) && !is.infinite(df))) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  critical <- if (is.infinite(df)) stats::qnorm((1 + conf.level) / 2) else
    stats::qt((1 + conf.level) / 2, df = df)
  c(lower = estimate - critical * se, upper = estimate + critical * se)
}

.interference_add_by <- function(row, lookup, stratum) {
  z <- lookup[lookup$.stratum == stratum, , drop = FALSE]
  if (!nrow(z)) return(row)
  columns <- setdiff(names(z), ".stratum")
  if (!length(columns)) return(row)
  by_values <- z[rep(1L, nrow(row)), columns, drop = FALSE]
  rownames(by_values) <- NULL
  rownames(row) <- NULL
  cbind(by_values, row)
}

#' Calculate replicate counts for an EP07 paired-difference experiment
#'
#' Calculates the number of replicate measurements needed in each of the test
#' and control samples. `sd` and `detectable_difference` must use the same
#' scale: both absolute units, or both relative units (SD/CV and percent).
#' The calculation is vectorized and reports design parameters only.
#'
#' @param sd Assumed repeatability SD, or repeatability CV when a relative
#'   calculation is intended.
#' @param detectable_difference Smallest interference difference the experiment
#'   is designed to detect, in the same units as `sd`.
#' @param alpha Type I error probability.
#' @param power Desired statistical power.
#' @param alternative `"two.sided"` or `"one.sided"`.
#' @param min_replicates Minimum returned replicate count per test/control sample.
#' @param label Optional labels, one per calculation.
#'
#' @return A data frame of class `interference_replicates`.
#' @references CLSI. Interference Testing in Clinical Chemistry. EP07, 3rd ed.
#' @examples
#' interference_replicates(sd = c(5, 5),
#'                         detectable_difference = c(7, 10),
#'                         label = c("7 percent", "10 percent"))
#' @export
interference_replicates <- function(sd, detectable_difference,
                                    alpha = 0.05, power = 0.90,
                                    alternative = c("two.sided", "one.sided"),
                                    min_replicates = 5L, label = NULL) {
  alternative <- match.arg(alternative)
  alpha <- .interference_probability(alpha, "alpha")
  power <- .interference_probability(power, "power")
  if (!is.numeric(sd) || !is.numeric(detectable_difference) ||
      anyNA(sd) || anyNA(detectable_difference) ||
      any(!is.finite(sd)) || any(!is.finite(detectable_difference)) ||
      any(sd <= 0) || any(detectable_difference <= 0)) {
    stop("`sd` and `detectable_difference` must contain positive finite values.",
         call. = FALSE)
  }
  n <- max(length(sd), length(detectable_difference))
  if (any(!c(length(sd), length(detectable_difference)) %in% c(1L, n))) {
    stop("`sd` and `detectable_difference` must have length 1 or a common length.",
         call. = FALSE)
  }
  if (!is.numeric(min_replicates) || length(min_replicates) != 1L ||
      is.na(min_replicates) || min_replicates < 1 ||
      min_replicates != as.integer(min_replicates)) {
    stop("`min_replicates` must be a positive integer.", call. = FALSE)
  }
  sd <- rep_len(as.numeric(sd), n)
  detectable_difference <- rep_len(as.numeric(detectable_difference), n)
  if (is.null(label)) label <- paste0("calculation_", seq_len(n))
  if (length(label) != n || anyNA(label)) {
    stop("`label` must have one value per calculation.", call. = FALSE)
  }
  z_alpha <- stats::qnorm(1 - if (alternative == "two.sided") alpha / 2 else alpha)
  z_power <- stats::qnorm(power)
  raw <- 2 * ((z_alpha + z_power) * sd / detectable_difference)^2
  out <- data.frame(
    label = as.character(label), sd = sd,
    detectable_difference = detectable_difference,
    difference_over_sd = detectable_difference / sd,
    alpha = alpha, power = power, alternative = alternative,
    z_alpha = z_alpha, z_power = z_power,
    calculated_replicates = raw,
    replicates_per_sample = pmax(as.integer(min_replicates), ceiling(raw)),
    min_replicates = as.integer(min_replicates), row.names = NULL
  )
  class(out) <- c("interference_replicates", "data.frame")
  out
}

#' Analyze paired-difference interference data
#'
#' Supports the standard EP07 test/control preparation and two laboratory
#' extensions: analyte spiking in normal/interferent matrices and direct
#' measurement of interferent dissolved in a blank. One row represents one
#' replicate result. Multiple interferents or measurand levels are handled by
#' supplying their columns in `by`.
#'
#' @param data Data frame in long format.
#' @param result Numeric measurement-result column.
#' @param design `"ep07"`, `"analyte_spike"`, or `"blank"`.
#' @param by Optional character vector of stratification columns.
#' @param condition Test/control condition column for `ep07`, or optional
#'   condition column for `blank`.
#' @param test_level Value identifying the test/interferent condition.
#' @param control_level Value identifying the solvent control condition.
#' @param matrix Matrix-type column for `analyte_spike`.
#' @param normal_level Value identifying the normal matrix.
#' @param interferent_level Value identifying the interferent-containing matrix.
#' @param spike Spike-status column for `analyte_spike`.
#' @param unspiked_level Value identifying measurements before analyte spiking.
#' @param spiked_level Value identifying measurements after analyte spiking.
#' @param conf.level Confidence level for intervals.
#'
#' @return An S3 object of class `interference_paired` containing `effects`,
#'   `cell_summary`, `issues`, and excluded row indices.
#' @examples
#' ep07_data <- data.frame(
#'   interferent = rep(c("bilirubin", "hemoglobin"), each = 10),
#'   condition = rep(rep(c("control", "test"), each = 5), 2),
#'   result = c(0.45, 0.56, 0.48, 0.54, 0.47,
#'              0.24, 0.28, 0.35, 0.37, 0.26,
#'              5.02, 4.98, 5.08, 4.95, 5.00,
#'              5.20, 5.12, 5.18, 5.10, 5.15)
#' )
#' interference_paired(ep07_data, "result", design = "ep07",
#'                     condition = "condition", by = "interferent")
#'
#' spike_data <- expand.grid(
#'   matrix = c("normal", "interferent"), spike = c("before", "after"),
#'   replicate = 1:4, KEEP.OUT.ATTRS = FALSE
#' )
#' spike_data$result <- with(spike_data,
#'   ifelse(matrix == "normal", 10, 12) +
#'   ifelse(spike == "after", ifelse(matrix == "normal", 10, 8), 0))
#' interference_paired(spike_data, "result", design = "analyte_spike",
#'                     matrix = "matrix", spike = "spike",
#'                     unspiked_level = "before", spiked_level = "after")
#' @export
interference_paired <- function(data, result,
                                design = c("ep07", "analyte_spike", "blank"),
                                by = NULL, condition = NULL,
                                test_level = "test", control_level = "control",
                                matrix = NULL, normal_level = "normal",
                                interferent_level = "interferent", spike = NULL,
                                unspiked_level = "unspiked",
                                spiked_level = "spiked", conf.level = 0.95) {
  call <- match.call()
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  design <- match.arg(design)
  result <- .interference_col(data, result, "result", numeric = TRUE)
  by <- .interference_by(data, by)
  conf.level <- .interference_probability(conf.level, "conf.level")
  labels <- list(test_level = test_level, control_level = control_level,
                 normal_level = normal_level,
                 interferent_level = interferent_level,
                 unspiked_level = unspiked_level, spiked_level = spiked_level)
  if (any(!vapply(labels, function(z) length(z) == 1L && !is.na(z), logical(1)))) {
    stop("Condition and preparation level arguments must each be one nonmissing value.",
         call. = FALSE)
  }
  if (identical(as.character(test_level), as.character(control_level)) ||
      identical(as.character(normal_level), as.character(interferent_level)) ||
      identical(as.character(unspiked_level), as.character(spiked_level))) {
    stop("The two labels defining each experimental axis must be distinct.",
         call. = FALSE)
  }
  if (design == "ep07") {
    condition <- .interference_col(data, condition, "condition")
  } else if (design == "analyte_spike") {
    matrix <- .interference_col(data, matrix, "matrix")
    spike <- .interference_col(data, spike, "spike")
  } else if (!is.null(condition)) {
    condition <- .interference_col(data, condition, "condition")
  }

  required <- unique(c(result, by, condition, matrix, spike))
  complete <- rep(TRUE, nrow(data))
  for (nm in required) {
    complete <- complete & !is.na(data[[nm]])
    if (is.numeric(data[[nm]])) complete <- complete & is.finite(data[[nm]])
  }
  excluded <- which(!complete)
  d <- data[complete, , drop = FALSE]
  if (!nrow(d)) stop("No complete observations are available.", call. = FALSE)
  d$.result <- as.numeric(d[[result]])
  d$.stratum <- .interference_strata(d, by)
  if (!is.null(condition)) d$.condition <- as.character(d[[condition]])
  if (!is.null(matrix)) d$.matrix <- as.character(d[[matrix]])
  if (!is.null(spike)) d$.spike <- as.character(d[[spike]])
  lookup <- .interference_group_values(d, by)
  split_rows <- split(seq_len(nrow(d)), factor(d$.stratum,
                                               levels = unique(d$.stratum)))
  effects <- cells <- issues <- list()

  add_cells <- function(z, labels, stratum) {
    rows <- lapply(names(labels), function(nm) {
      st <- .interference_stats(labels[[nm]])
      st$cell <- nm
      st$stratum <- stratum
      st
    })
    .interference_bind(rows)
  }

  for (key in names(split_rows)) {
    z <- d[split_rows[[key]], , drop = FALSE]
    effect_row <- NULL
    if (design == "ep07") {
      test <- z$.result[z$.condition == as.character(test_level)]
      control <- z$.result[z$.condition == as.character(control_level)]
      if (length(test) < 2L || length(control) < 2L) {
        issues[[length(issues) + 1L]] <- .interference_issue(
          key, "missing_cell", "EP07 design needs at least two test and two control results.")
        next
      }
      s_t <- stats::sd(test); s_c <- stats::sd(control)
      n_t <- length(test); n_c <- length(control)
      estimate <- mean(test) - mean(control)
      se <- sqrt(s_t^2 / n_t + s_c^2 / n_c)
      df <- n_t + n_c - 2
      ci_abs <- .interference_ci(estimate, se, df, conf.level)
      if (abs(mean(control)) <= .Machine$double.eps^0.5) {
        pct <- pct_se <- NA_real_; ci_pct <- c(lower = NA_real_, upper = NA_real_)
        issues[[length(issues) + 1L]] <- .interference_issue(
          key, "zero_control_mean", "Percent difference is undefined because the control mean is zero.")
      } else {
        pct <- 100 * estimate / mean(control)
        pct_se <- 100 * se / abs(mean(control))
        ci_pct <- sort(100 * ci_abs / mean(control))
        names(ci_pct) <- c("lower", "upper")
      }
      effect_row <- data.frame(
        stratum = key, design = design, n_test = n_t, n_control = n_c,
        test_mean = mean(test), control_mean = mean(control),
        absolute_difference = estimate, absolute_se = se, df = df,
        absolute_ci_lower = ci_abs[1L], absolute_ci_upper = ci_abs[2L],
        percent_difference = pct, percent_se = pct_se,
        percent_ci_lower = ci_pct[1L], percent_ci_upper = ci_pct[2L],
        stringsAsFactors = FALSE
      )
      cells[[length(cells) + 1L]] <- add_cells(
        z, list(test = test, control = control), key)
    } else if (design == "analyte_spike") {
      pick <- function(m, s) z$.result[z$.matrix == as.character(m) &
                                        z$.spike == as.character(s)]
      n_before <- pick(normal_level, unspiked_level)
      n_after <- pick(normal_level, spiked_level)
      i_before <- pick(interferent_level, unspiked_level)
      i_after <- pick(interferent_level, spiked_level)
      values <- list(normal_unspiked = n_before, normal_spiked = n_after,
                     interferent_unspiked = i_before,
                     interferent_spiked = i_after)
      ns <- vapply(values, length, integer(1))
      if (any(ns < 2L)) {
        issues[[length(issues) + 1L]] <- .interference_issue(
          key, "missing_cell",
          "Analyte-spike design needs at least two results in each of its four cells.")
        next
      }
      means <- vapply(values, mean, numeric(1))
      vars <- vapply(values, stats::var, numeric(1))
      normal_increment <- means["normal_spiked"] - means["normal_unspiked"]
      interferent_increment <- means["interferent_spiked"] -
        means["interferent_unspiked"]
      did <- interferent_increment - normal_increment
      abs_terms <- vars / ns
      abs_se <- sqrt(sum(abs_terms))
      abs_df <- .interference_df(vars, ns)
      abs_ci <- .interference_ci(did, abs_se, abs_df, conf.level)
      if (abs(normal_increment) <= .Machine$double.eps^0.5) {
        issues[[length(issues) + 1L]] <- .interference_issue(
          key, "zero_normal_increment",
          "Double-difference percent is undefined because the normal-sample increment is zero.")
        next
      }
      pct <- 100 * (interferent_increment / normal_increment - 1)
      gradients <- c(
        normal_unspiked = unname(100 * interferent_increment / normal_increment^2),
        normal_spiked = unname(-100 * interferent_increment / normal_increment^2),
        interferent_unspiked = unname(-100 / normal_increment),
        interferent_spiked = unname(100 / normal_increment)
      )
      pct_components <- gradients[names(values)]^2 * vars / ns
      pct_se <- sqrt(sum(pct_components))
      pct_df_den <- sum(pct_components^2 / (ns - 1))
      pct_df <- if (pct_df_den > 0) sum(pct_components)^2 / pct_df_den else Inf
      pct_ci <- .interference_ci(pct, pct_se, pct_df, conf.level)
      effect_row <- data.frame(
        stratum = key, design = design,
        normal_increment = unname(normal_increment),
        interferent_increment = unname(interferent_increment),
        absolute_double_difference = unname(did), absolute_se = abs_se,
        absolute_df = abs_df, absolute_ci_lower = abs_ci[1L],
        absolute_ci_upper = abs_ci[2L],
        percent_double_difference = unname(pct), percent_se = pct_se,
        percent_df = pct_df, percent_ci_lower = pct_ci[1L],
        percent_ci_upper = pct_ci[2L],
        stringsAsFactors = FALSE
      )
      cells[[length(cells) + 1L]] <- add_cells(z, values, key)
    } else {
      if (is.null(condition)) {
        test <- z$.result; control <- numeric()
      } else {
        test <- z$.result[z$.condition == as.character(test_level)]
        control <- z$.result[z$.condition == as.character(control_level)]
      }
      if (length(test) < 2L) {
        issues[[length(issues) + 1L]] <- .interference_issue(
          key, "missing_test", "Blank design needs at least two interferent-blank results.")
        next
      }
      if (length(control) >= 2L) {
        estimate <- mean(test) - mean(control)
        se <- sqrt(stats::var(test) / length(test) +
                     stats::var(control) / length(control))
        df <- length(test) + length(control) - 2
        ci_abs <- .interference_ci(estimate, se, df, conf.level)
        comparison <- "solvent_control"
      } else {
        if (length(control) == 1L) {
          issues[[length(issues) + 1L]] <- .interference_issue(
            key, "insufficient_control",
            "Only one solvent-control result was present; the analysis used theoretical zero instead.")
        }
        estimate <- mean(test)
        se <- stats::sd(test) / sqrt(length(test))
        df <- length(test) - 1L
        ci_abs <- .interference_ci(estimate, se, df, conf.level)
        comparison <- "theoretical_zero"
      }
      effect_row <- data.frame(
        stratum = key, design = design, comparison = comparison,
        n_test = length(test), n_control = length(control),
        test_mean = mean(test),
        control_mean = if (length(control)) mean(control) else 0,
        apparent_signal = estimate, se = se, df = df,
        ci_lower = ci_abs[1L], ci_upper = ci_abs[2L],
        stringsAsFactors = FALSE
      )
      cells[[length(cells) + 1L]] <- add_cells(
        z, c(list(interferent_blank = test),
             if (length(control)) list(solvent_control = control) else list()), key)
    }
    effects[[length(effects) + 1L]] <- .interference_add_by(effect_row, lookup, key)
  }

  effects <- .interference_bind(effects)
  if (!nrow(effects)) {
    detail <- .interference_issues(issues)
    stop("No stratum satisfied the data requirements for design '", design,
         "'.", if (nrow(detail)) paste0(" First issue: ", detail$message[1L]) else "",
         call. = FALSE)
  }
  out <- list(call = call, design = design, result = result, by = by,
              conf.level = conf.level, effects = effects,
              cell_summary = .interference_bind(cells),
              issues = .interference_issues(issues), excluded_rows = excluded,
              analysis_data = d)
  class(out) <- "interference_paired"
  out
}

#' Analyze an EP07 interference dose-response experiment
#'
#' Summarizes replicate results at each interferent concentration relative to
#' an observed baseline concentration. The recommended point-to-point method
#' is the default; a straight-line regression can be requested explicitly.
#'
#' @param data Long-format data frame.
#' @param result Numeric measurement-result column.
#' @param concentration Numeric interferent-concentration column.
#' @param by Optional stratification columns, typically interferent name and
#'   measurand level.
#' @param baseline Observed interferent concentration used to estimate M0.
#' @param method `"point_to_point"` or `"linear"`.
#' @param conf.level Confidence level for regression coefficients.
#'
#' @return An object of class `interference_dose_response`.
#' @examples
#' dose_data <- expand.grid(
#'   interferent = c("A", "B"), concentration = c(0, 10, 20, 30, 40),
#'   replicate = 1:5, KEEP.OUT.ATTRS = FALSE
#' )
#' dose_data$result <- 100 + ifelse(dose_data$interferent == "A", 0.5, -0.3) *
#'   dose_data$concentration + stats::rnorm(nrow(dose_data), 0, 0.5)
#' dose_fit <- interference_dose_response(
#'   dose_data, "result", "concentration", by = "interferent")
#' dose_fit$summary
#' @export
interference_dose_response <- function(data, result, concentration, by = NULL,
                                       baseline = 0,
                                       method = c("point_to_point", "linear"),
                                       conf.level = 0.95) {
  call <- match.call()
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  result <- .interference_col(data, result, "result", numeric = TRUE)
  concentration <- .interference_col(data, concentration, "concentration",
                                     numeric = TRUE)
  by <- .interference_by(data, by)
  method <- match.arg(method)
  conf.level <- .interference_probability(conf.level, "conf.level")
  if (!is.numeric(baseline) || length(baseline) != 1L || is.na(baseline) ||
      !is.finite(baseline)) stop("`baseline` must be one finite number.", call. = FALSE)
  required <- unique(c(result, concentration, by))
  complete <- rep(TRUE, nrow(data))
  for (nm in required) {
    complete <- complete & !is.na(data[[nm]])
    if (is.numeric(data[[nm]])) complete <- complete & is.finite(data[[nm]])
  }
  excluded <- which(!complete)
  d <- data[complete, , drop = FALSE]
  if (!nrow(d)) stop("No complete observations are available.", call. = FALSE)
  d$.result <- as.numeric(d[[result]])
  d$.concentration <- as.numeric(d[[concentration]])
  d$.stratum <- .interference_strata(d, by)
  d$stratum <- d$.stratum
  lookup <- .interference_group_values(d, by)
  splits <- split(seq_len(nrow(d)), factor(d$.stratum, levels = unique(d$.stratum)))
  summaries <- diagnostics <- coefficients <- issues <- list()
  models <- list()

  for (key in names(splits)) {
    z <- d[splits[[key]], , drop = FALSE]
    levels <- sort(unique(z$.concentration))
    if (length(levels) < 2L) {
      issues[[length(issues) + 1L]] <- .interference_issue(
        key, "too_few_levels", "At least two interferent concentrations are required.")
      next
    }
    if (!baseline %in% levels) {
      issues[[length(issues) + 1L]] <- .interference_issue(
        key, "baseline_absent", paste0("No observations were found at baseline ", baseline, "."))
      next
    }
    rows <- lapply(levels, function(x) {
      st <- .interference_stats(z$.result[z$.concentration == x])
      st$concentration <- x
      st$cv_percent <- if (is.finite(st$mean) && abs(st$mean) > .Machine$double.eps^0.5)
        100 * st$sd / abs(st$mean) else NA_real_
      st
    })
    sm <- .interference_bind(rows)
    m0 <- sm$mean[sm$concentration == baseline][1L]
    sm$absolute_difference <- sm$mean - m0
    sm$percent_difference <- if (abs(m0) > .Machine$double.eps^0.5)
      100 * sm$absolute_difference / m0 else NA_real_
    sm$stratum <- key
    sm <- .interference_add_by(sm, lookup, key)
    summaries[[length(summaries) + 1L]] <- sm
    diffs <- diff(sm$mean)
    monotonic <- if (all(diffs >= 0)) "increasing" else
      if (all(diffs <= 0)) "decreasing" else "nonmonotonic"
    diag <- data.frame(
      stratum = key, levels = length(levels), total_results = nrow(z),
      baseline = baseline, m0 = m0,
      min_replicates = min(sm$n), max_replicates = max(sm$n),
      monotonicity = monotonic,
      stringsAsFactors = FALSE
    )
    diagnostics[[length(diagnostics) + 1L]] <- .interference_add_by(diag, lookup, key)
    if (method == "linear") {
      fit <- stats::lm(.result ~ .concentration, data = z)
      models[[key]] <- fit
      cf <- summary(fit)$coefficients
      critical <- stats::qt((1 + conf.level) / 2, stats::df.residual(fit))
      cft <- data.frame(
        stratum = key, term = rownames(cf), estimate = cf[, 1L],
        std_error = cf[, 2L], statistic = cf[, 3L], p_value = cf[, 4L],
        ci_lower = cf[, 1L] - critical * cf[, 2L],
        ci_upper = cf[, 1L] + critical * cf[, 2L],
        r_squared = summary(fit)$r.squared,
        adjusted_r_squared = summary(fit)$adj.r.squared,
        residual_sd = summary(fit)$sigma,
        df_residual = stats::df.residual(fit), row.names = NULL
      )
      coefficients[[length(coefficients) + 1L]] <- .interference_add_by(cft, lookup, key)
    }
  }
  summary_table <- .interference_bind(summaries)
  if (!nrow(summary_table)) stop("No stratum could be analyzed.", call. = FALSE)
  out <- list(call = call, result = result, concentration = concentration,
              by = by, baseline = baseline, method = method,
              conf.level = conf.level, summary = summary_table,
              diagnostics = .interference_bind(diagnostics),
              coefficients = .interference_bind(coefficients), models = models,
              issues = .interference_issues(issues), excluded_rows = excluded,
              analysis_data = d)
  class(out) <- "interference_dose_response"
  out
}

#' Predict from an interference dose-response analysis
#'
#' @param object An `interference_dose_response` object.
#' @param newdata Numeric interferent concentrations for forward prediction.
#' @param target Numeric effects for inverse prediction. Supply instead of
#'   `newdata`.
#' @param metric `"result"`, `"absolute"`, or `"percent"`.
#' @param stratum Optional stratum key; by default all strata are processed.
#' @param ... Unused.
#'
#' @return A data frame of predictions or inverse interpolations.
#' @examples
#' d <- expand.grid(conc = c(0, 10, 20, 30, 40), rep = 1:5)
#' d$result <- 100 + 0.5 * d$conc
#' fit <- interference_dose_response(d, "result", "conc")
#' predict(fit, newdata = c(5, 15), metric = "percent")
#' predict(fit, target = 10, metric = "percent")
#' @export
predict.interference_dose_response <- function(object, newdata = NULL,
                                                target = NULL,
                                                metric = c("result", "absolute", "percent"),
                                                stratum = NULL, ...) {
  metric <- match.arg(metric)
  if (!inherits(object, "interference_dose_response"))
    stop("`object` must be an interference_dose_response object.", call. = FALSE)
  if (is.null(newdata) == is.null(target))
    stop("Supply exactly one of `newdata` or `target`.", call. = FALSE)
  values <- if (!is.null(newdata)) newdata else target
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)))
    stop("Prediction values must be finite numeric values.", call. = FALSE)
  keys <- unique(object$summary$stratum)
  if (!is.null(stratum)) {
    if (!all(stratum %in% keys)) stop("Unknown `stratum` value.", call. = FALSE)
    keys <- as.character(stratum)
  }
  rows <- list()
  for (key in keys) {
    sm <- object$summary[object$summary$stratum == key, , drop = FALSE]
    sm <- sm[order(sm$concentration), ]
    y <- switch(metric, result = sm$mean,
                absolute = sm$absolute_difference,
                percent = sm$percent_difference)
    if (any(!is.finite(y))) {
      pred <- rep(NA_real_, length(values))
    } else if (is.null(target)) {
      if (object$method == "linear") {
        fit <- object$models[[key]]
        response <- unname(stats::coef(fit)[1L] + stats::coef(fit)[2L] * values)
        m0 <- sm$mean[sm$concentration == object$baseline][1L]
        pred <- switch(metric, result = response, absolute = response - m0,
                       percent = 100 * (response - m0) / m0)
        pred[values < min(sm$concentration) | values > max(sm$concentration)] <- NA_real_
      } else {
        pred <- stats::approx(sm$concentration, y, xout = values, rule = 1,
                              ties = "ordered")$y
      }
    } else {
      if (object$method == "linear") {
        fit <- object$models[[key]]
        slope <- unname(stats::coef(fit)[2L])
        intercept <- unname(stats::coef(fit)[1L])
        m0 <- sm$mean[sm$concentration == object$baseline][1L]
        response_target <- switch(
          metric, result = values, absolute = m0 + values,
          percent = m0 * (1 + values / 100))
        pred <- if (is.finite(slope) && abs(slope) > .Machine$double.eps^0.5)
          (response_target - intercept) / slope else rep(NA_real_, length(values))
        pred[pred < min(sm$concentration) | pred > max(sm$concentration)] <- NA_real_
      } else if (!(all(diff(y) >= 0) || all(diff(y) <= 0))) {
        pred <- rep(NA_real_, length(values))
      } else {
        ord <- order(y)
        yy <- y[ord]; xx <- sm$concentration[ord]
        keep <- !duplicated(yy)
        pred <- stats::approx(yy[keep], xx[keep], xout = values, rule = 1,
                              ties = "ordered")$y
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      stratum = key, direction = if (is.null(target)) "forward" else "inverse",
      metric = metric, input = values, estimate = pred,
      within_observed_range = is.finite(pred), stringsAsFactors = FALSE)
  }
  .interference_bind(rows)
}

#' Evaluate interference with patient specimens
#'
#' Averages replicates within specimen and method, then calculates the evaluated
#' procedure result minus the comparative procedure result. Potential
#' interferent concentrations can be supplied as multiple numeric columns and
#' are modeled separately.
#'
#' @param data Long-format patient-result data.
#' @param id Specimen identifier column.
#' @param group Patient test/control group column.
#' @param method Measurement-procedure column.
#' @param result Numeric result column.
#' @param test_group Value identifying the selected/test patient group.
#' @param control_group Value identifying the control patient group.
#' @param evaluated_method Value identifying the evaluated procedure.
#' @param comparative_method Value identifying the comparative procedure.
#' @param interferent_cols Optional numeric interferent-concentration columns.
#' @param conf.level Confidence level for summaries and regression coefficients.
#'
#' @return An object of class `interference_patient`.
#' @examples
#' patient_data <- expand.grid(
#'   id = paste0("S", 1:12), method = c("evaluated", "comparative"),
#'   replicate = 1:2, KEEP.OUT.ATTRS = FALSE
#' )
#' patient_data$group <- ifelse(as.integer(sub("S", "", patient_data$id)) <= 6,
#'                              "control", "test")
#' patient_data$bilirubin <- 2 * (as.integer(sub("S", "", patient_data$id)) - 1)
#' truth <- 50 + as.integer(sub("S", "", patient_data$id))
#' patient_data$result <- truth +
#'   ifelse(patient_data$method == "evaluated", 0.08 * patient_data$bilirubin, 0)
#' patient_fit <- interference_patient(
#'   patient_data, "id", "group", "method", "result",
#'   test_group = "test", control_group = "control",
#'   evaluated_method = "evaluated", comparative_method = "comparative",
#'   interferent_cols = "bilirubin")
#' patient_fit$group_summary
#' @export
interference_patient <- function(data, id, group, method, result,
                                 test_group, control_group,
                                 evaluated_method, comparative_method,
                                 interferent_cols = NULL,
                                 conf.level = 0.95) {
  call <- match.call()
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  id <- .interference_col(data, id, "id")
  group <- .interference_col(data, group, "group")
  method <- .interference_col(data, method, "method")
  result <- .interference_col(data, result, "result", numeric = TRUE)
  conf.level <- .interference_probability(conf.level, "conf.level")
  labels <- list(test_group = test_group, control_group = control_group,
                 evaluated_method = evaluated_method,
                 comparative_method = comparative_method)
  if (any(!vapply(labels, function(z) length(z) == 1L && !is.na(z), logical(1)))) {
    stop("Group and method labels must each be one nonmissing value.", call. = FALSE)
  }
  if (identical(as.character(test_group), as.character(control_group)) ||
      identical(as.character(evaluated_method), as.character(comparative_method))) {
    stop("Test/control groups and evaluated/comparative methods must be distinct.",
         call. = FALSE)
  }
  if (!is.null(interferent_cols)) {
    if (!is.character(interferent_cols) || anyNA(interferent_cols) ||
        any(!interferent_cols %in% names(data))) {
      stop("`interferent_cols` must contain column names found in `data`.", call. = FALSE)
    }
    bad <- interferent_cols[!vapply(data[interferent_cols], is.numeric, logical(1))]
    if (length(bad)) stop("Interferent concentration columns must be numeric: ",
                          paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  complete <- !is.na(data[[id]]) & !is.na(data[[group]]) & !is.na(data[[method]]) &
    !is.na(data[[result]]) & is.finite(data[[result]])
  excluded <- which(!complete)
  d <- data[complete, , drop = FALSE]
  d$.id <- as.character(d[[id]])
  d$.group <- as.character(d[[group]])
  d$.method <- as.character(d[[method]])
  d$.result <- as.numeric(d[[result]])
  keep_groups <- c(as.character(test_group), as.character(control_group))
  keep_methods <- c(as.character(evaluated_method), as.character(comparative_method))
  d <- d[d$.group %in% keep_groups & d$.method %in% keep_methods, , drop = FALSE]
  if (!nrow(d)) stop("No rows match the requested groups and methods.", call. = FALSE)
  group_by_id <- split(d$.group, d$.id)
  inconsistent_group <- names(Filter(function(z) length(unique(z)) != 1L, group_by_id))
  if (length(inconsistent_group)) {
    stop("Patient group is not constant within specimen(s): ",
         paste(inconsistent_group, collapse = ", "), ".", call. = FALSE)
  }
  key <- interaction(d$.id, d$.method, drop = TRUE, lex.order = TRUE)
  means <- lapply(split(seq_len(nrow(d)), key), function(i) {
    data.frame(id = d$.id[i[1L]], group = d$.group[i[1L]],
               method = d$.method[i[1L]], n = length(i),
               mean = mean(d$.result[i]),
               sd = if (length(i) > 1L) stats::sd(d$.result[i]) else NA_real_,
               stringsAsFactors = FALSE)
  })
  means <- .interference_bind(means)
  ev <- means[means$method == as.character(evaluated_method), ]
  cp <- means[means$method == as.character(comparative_method), ]
  pairs <- merge(ev, cp, by = c("id", "group"), suffixes = c("_evaluated", "_comparative"))
  all_ids <- unique(d$.id)
  incomplete <- setdiff(all_ids, pairs$id)
  if (!nrow(pairs)) stop("No specimen has results from both requested methods.", call. = FALSE)
  pairs$absolute_difference <- pairs$mean_evaluated - pairs$mean_comparative
  pairs$percent_difference <- ifelse(
    abs(pairs$mean_comparative) > .Machine$double.eps^0.5,
    100 * pairs$absolute_difference / pairs$mean_comparative, NA_real_)

  summaries <- lapply(keep_groups, function(g) {
    z <- pairs$absolute_difference[pairs$group == g]
    zp <- pairs$percent_difference[pairs$group == g]
    n <- length(z)
    if (!n) return(NULL)
    se <- if (n > 1L) stats::sd(z) / sqrt(n) else NA_real_
    ci <- if (n > 1L) .interference_ci(mean(z), se, n - 1L, conf.level) else c(NA, NA)
    se_p <- if (sum(is.finite(zp)) > 1L) stats::sd(zp, na.rm = TRUE) /
      sqrt(sum(is.finite(zp))) else NA_real_
    ci_p <- if (sum(is.finite(zp)) > 1L)
      .interference_ci(mean(zp, na.rm = TRUE), se_p,
                       sum(is.finite(zp)) - 1L, conf.level) else c(NA, NA)
    data.frame(
      group = g, n = n, mean_difference = mean(z),
      sd_difference = if (n > 1L) stats::sd(z) else NA_real_, se_difference = se,
      median_difference = stats::median(z), q25_difference = stats::quantile(z, .25),
      q75_difference = stats::quantile(z, .75),
      ci_lower = ci[1L], ci_upper = ci[2L],
      mean_percent_difference = mean(zp, na.rm = TRUE),
      sd_percent_difference = stats::sd(zp, na.rm = TRUE),
      percent_ci_lower = ci_p[1L], percent_ci_upper = ci_p[2L],
      stringsAsFactors = FALSE, row.names = NULL)
  })
  group_summary <- .interference_bind(summaries)

  interferent_data <- data.frame(id = unique(d$.id), stringsAsFactors = FALSE)
  issues <- list(); regressions <- list(); models <- list()
  if (length(interferent_cols)) {
    id_rows <- split(seq_len(nrow(d)), d$.id)
    concentration_rows <- lapply(names(id_rows), function(identifier) {
      i <- id_rows[[identifier]]
      vals <- lapply(interferent_cols, function(nm) {
        v <- unique(d[[nm]][i][is.finite(d[[nm]][i])])
        if (length(v) > 1L) {
          issues[[length(issues) + 1L]] <<- .interference_issue(
            identifier, "inconsistent_concentration",
            paste0("Column '", nm, "' has multiple concentrations within this specimen; their mean was used."))
        }
        if (length(v)) mean(v) else NA_real_
      })
      out <- data.frame(id = identifier, stringsAsFactors = FALSE)
      for (j in seq_along(interferent_cols)) out[[interferent_cols[j]]] <- vals[[j]]
      out
    })
    interferent_data <- .interference_bind(concentration_rows)
    pairs <- merge(pairs, interferent_data, by = "id", all.x = TRUE)
    for (nm in interferent_cols) {
      ok <- is.finite(pairs[[nm]]) & is.finite(pairs$absolute_difference)
      if (sum(ok) < 3L || length(unique(pairs[[nm]][ok])) < 2L) {
        issues[[length(issues) + 1L]] <- .interference_issue(
          nm, "insufficient_regression_data",
          "At least three complete specimens at two interferent concentrations are required for regression.")
        next
      }
      fit_data <- data.frame(difference = pairs$absolute_difference[ok],
                             concentration = pairs[[nm]][ok])
      fit <- stats::lm(difference ~ concentration, data = fit_data)
      models[[nm]] <- fit
      cf <- summary(fit)$coefficients
      critical <- stats::qt((1 + conf.level) / 2, stats::df.residual(fit))
      regressions[[length(regressions) + 1L]] <- data.frame(
        interferent = nm, term = rownames(cf), estimate = cf[, 1L],
        std_error = cf[, 2L], statistic = cf[, 3L], p_value = cf[, 4L],
        ci_lower = cf[, 1L] - critical * cf[, 2L],
        ci_upper = cf[, 1L] + critical * cf[, 2L],
        n = sum(ok), r_squared = summary(fit)$r.squared,
        adjusted_r_squared = summary(fit)$adj.r.squared,
        residual_sd = summary(fit)$sigma, row.names = NULL)
    }
  }
  out <- list(call = call, id = id, group = group, method = method,
              result = result, test_group = as.character(test_group),
              control_group = as.character(control_group),
              evaluated_method = as.character(evaluated_method),
              comparative_method = as.character(comparative_method),
              interferent_cols = interferent_cols, conf.level = conf.level,
              specimen_results = pairs, method_summary = means,
              group_summary = group_summary,
              regressions = .interference_bind(regressions), models = models,
              issues = .interference_issues(issues), incomplete_ids = incomplete,
              excluded_rows = excluded)
  class(out) <- "interference_patient"
  out
}

#' Print EP07 replicate-count results
#' @param x An `interference_replicates` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(interference_replicates(5, 7))
#' @export
print.interference_replicates <- function(x, ...) {
  cat("EP07 paired-difference replicate design\n")
  .print_kv_rows(x, list(
    "Design inputs" = c(
      "SD" = "sd", "Detectable difference" = "detectable_difference",
      "Difference / SD" = "difference_over_sd"),
    "Significance and power" = c(
      "Alpha" = "alpha", "Power" = "power", "Alternative" = "alternative"),
    "Replicate calculation" = c(
      "Z alpha" = "z_alpha", "Z power" = "z_power",
      "Calculated replicates" = "calculated_replicates"),
    "Replicate recommendation" = c(
      "Replicates per sample" = "replicates_per_sample",
      "Minimum replicates" = "min_replicates")
  ), row_title = c("Label" = "label"))
  invisible(x)
}

#' Print paired interference estimates
#' @param x An `interference_paired` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' d <- data.frame(group = rep(c("control", "test"), each = 5),
#'                 result = c(10, 10.1, 9.9, 10.2, 9.8, 11, 11.1, 10.9, 11.2, 10.8))
#' print(interference_paired(d, "result", "ep07", condition = "group"))
#' @export
print.interference_paired <- function(x, ...) {
  cat("Interference paired-difference analysis\n")
  cat("  Design: ", x$design, "\n", sep = "")
  cat("  Confidence level: ", format(100 * x$conf.level), "%\n", sep = "")
  section <- function(columns) intersect(columns, names(x$effects))
  sections <- switch(
    x$design,
    ep07 = list(
      "Analysis strata" = section(c("interferent", "stratum", "design")),
      "Sample counts" = section(c("interferent", "n_test", "n_control")),
      "Sample means" = section(c("interferent", "test_mean", "control_mean")),
      "Absolute estimate" = section(c("interferent", "absolute_difference",
                                       "absolute_se", "df")),
      "Absolute interval" = section(c("interferent", "absolute_ci_lower",
                                       "absolute_ci_upper")),
      "Percent estimate" = section(c("interferent", "percent_difference",
                                      "percent_se")),
      "Percent interval" = section(c("interferent", "percent_ci_lower",
                                      "percent_ci_upper"))
    ),
    analyte_spike = list(
      "Analysis strata" = section(c("stratum", "design")),
      "Preparation increments" = section(c("normal_increment", "interferent_increment")),
      "Absolute estimate" = section(c("stratum", "absolute_double_difference",
                                       "absolute_se", "absolute_df")),
      "Absolute interval" = section(c("stratum", "absolute_ci_lower",
                                       "absolute_ci_upper")),
      "Percent estimate" = section(c("stratum", "percent_double_difference",
                                      "percent_se", "percent_df")),
      "Percent interval" = section(c("stratum", "percent_ci_lower",
                                      "percent_ci_upper"))
    ),
    blank = list(
      "Analysis strata" = section(c("stratum", "design", "comparison")),
      "Sample counts" = section(c("stratum", "n_test", "n_control")),
      "Sample means" = section(c("stratum", "test_mean", "control_mean")),
      "Apparent estimate" = section(c("stratum", "apparent_signal", "se", "df")),
      "Apparent interval" = section(c("stratum", "ci_lower", "ci_upper"))
    )
  )
  row_title <- switch(
    x$design,
    ep07 = if ("interferent" %in% names(x$effects)) {
      c("Interferent" = "interferent")
    } else if (length(x$by) && x$by[1L] %in% names(x$effects)) {
      stats::setNames(x$by[1L], x$by[1L])
    } else {
      NULL
    },
    analyte_spike = c("Stratum" = "stratum"),
    blank = c("Stratum" = "stratum")
  )
  sections <- lapply(sections, function(fields) {
    fields[!unname(fields) %in% unname(row_title)]
  })
  .print_kv_rows(x$effects, sections, row_title = row_title)
  if (nrow(x$issues)) cat("  Issues recorded: ", nrow(x$issues), "\n", sep = "")
  invisible(x)
}

#' Plot paired interference estimates
#' @param x An `interference_paired` object.
#' @param ... Unused.
#' @return A ggplot object.
#' @examples
#' d <- data.frame(group = rep(c("control", "test"), each = 5),
#'                 result = c(10, 10.1, 9.9, 10.2, 9.8, 11, 11.1, 10.9, 11.2, 10.8))
#' plot(interference_paired(d, "result", "ep07", condition = "group"))
#' @export
plot.interference_paired <- function(x, ...) {
  d <- x$cell_summary
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$cell, y = .data$mean, color = .data$cell)
  ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$mean - .data$se,
                   ymax = .data$mean + .data$se),
      width = 0.15
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$stratum), scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Mean measurement result",
                  title = "Interference preparation cell summaries") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Print interference dose-response results
#' @param x An `interference_dose_response` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' d <- expand.grid(conc = c(0, 10, 20), rep = 1:5)
#' d$result <- 100 + 0.5 * d$conc
#' print(interference_dose_response(d, "result", "conc"))
#' @export
print.interference_dose_response <- function(x, ...) {
  cat("Interference dose-response analysis\n")
  cat("  Method: ", x$method, "\n", sep = "")
  cat("  Baseline concentration: ", x$baseline, "\n", sep = "")
  .print_kv_rows(x$summary, list(
    "Response summary" = c("Observations" = "n", "Mean" = "mean"),
    "Response precision" = c("SD" = "sd", "SE" = "se"),
    "Difference from baseline" = c(
      "Absolute difference" = "absolute_difference",
      "Percent difference" = "percent_difference"),
    "Relative precision" = c("CV percent" = "cv_percent")
  ), row_title = c("Stratum" = "stratum", "Concentration" = "concentration"))
  invisible(x)
}

#' Plot interference dose-response results
#' @param x An `interference_dose_response` object.
#' @param type `"response"`, `"absolute"`, or `"percent"`.
#' @param ... Unused.
#' @return A ggplot object.
#' @examples
#' d <- expand.grid(conc = c(0, 10, 20), rep = 1:5)
#' d$result <- 100 + 0.5 * d$conc
#' plot(interference_dose_response(d, "result", "conc"), type = "percent")
#' @export
plot.interference_dose_response <- function(x,
                                            type = c("response", "absolute", "percent"),
                                            ...) {
  type <- match.arg(type)
  sm <- x$summary
  y_col <- switch(type, response = "mean", absolute = "absolute_difference",
                  percent = "percent_difference")
  y_lab <- switch(type, response = "Mean measurement result",
                  absolute = "Absolute difference from M0",
                  percent = "Percent difference from M0")
  p <- ggplot2::ggplot(
    sm,
    ggplot2::aes(x = .data$concentration,
                 y = .data[[y_col]], group = .data$stratum)
  ) +
    ggplot2::geom_line(color = "#2166AC", linewidth = 0.8) +
    ggplot2::geom_point(color = "#2166AC", size = 2.3) +
    ggplot2::facet_wrap(ggplot2::vars(.data$stratum), scales = "free_y") +
    ggplot2::labs(x = "Interferent concentration", y = y_lab,
                  title = "Interference dose-response") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (type == "response") {
    p <- p + ggplot2::geom_point(
      data = x$analysis_data,
      ggplot2::aes(x = .data$.concentration, y = .data$.result),
      inherit.aes = FALSE, alpha = 0.35, color = "grey35")
  } else {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = "grey45",
                                  linetype = "dotted")
  }
  p
}

#' Print patient-specimen interference results
#' @param x An `interference_patient` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' d <- expand.grid(id = paste0("S", 1:6), method = c("test", "ref"), rep = 1:2)
#' d$group <- ifelse(as.integer(sub("S", "", d$id)) <= 3, "control", "selected")
#' d$result <- 10 + as.integer(sub("S", "", d$id)) + (d$method == "test")
#' print(interference_patient(d, "id", "group", "method", "result",
#'   "selected", "control", "test", "ref"))
#' @export
print.interference_patient <- function(x, ...) {
  cat("Patient-specimen interference evaluation\n")
  cat("  Evaluated - comparative: ", x$evaluated_method, " - ",
      x$comparative_method, "\n", sep = "")
  .print_kv_rows(x$group_summary, list(
    "Difference summary" = c(
      "Observations" = "n", "Mean difference" = "mean_difference",
      "SD difference" = "sd_difference"),
    "Difference uncertainty" = c(
      "SE difference" = "se_difference", "Median difference" = "median_difference"),
    "Difference quartiles" = c(
      "25th percentile" = "q25_difference", "75th percentile" = "q75_difference"),
    "Difference interval" = c("CI lower" = "ci_lower", "CI upper" = "ci_upper"),
    "Percent summary" = c(
      "Mean percent difference" = "mean_percent_difference",
      "SD percent difference" = "sd_percent_difference"),
    "Percent interval" = c(
      "Percent CI lower" = "percent_ci_lower",
      "Percent CI upper" = "percent_ci_upper")
  ), row_title = c("Group" = "group"))
  if (length(x$incomplete_ids))
    cat("  Incomplete method pairs: ", length(x$incomplete_ids), "\n", sep = "")
  invisible(x)
}

#' Plot patient-specimen interference results
#' @param x An `interference_patient` object.
#' @param type `"difference"` or `"interferent"`.
#' @param interferent Interferent concentration column to plot when
#'   `type = "interferent"`.
#' @param ... Unused.
#' @return A ggplot object.
#' @examples
#' d <- expand.grid(id = paste0("S", 1:8), method = c("test", "ref"), rep = 1:2)
#' d$group <- ifelse(as.integer(sub("S", "", d$id)) <= 4, "control", "selected")
#' d$bilirubin <- 2 * (as.integer(sub("S", "", d$id)) - 1)
#' d$result <- 20 + as.integer(sub("S", "", d$id)) +
#'   ifelse(d$method == "test", 0.1 * d$bilirubin, 0)
#' fit <- interference_patient(d, "id", "group", "method", "result",
#'   "selected", "control", "test", "ref", "bilirubin")
#' plot(fit, type = "difference")
#' @export
plot.interference_patient <- function(x, type = c("difference", "interferent"),
                                      interferent = NULL, ...) {
  type <- match.arg(type)
  d <- x$specimen_results
  if (type == "difference") {
    return(ggplot2::ggplot(
      d, ggplot2::aes(x = .data$mean_comparative,
                      y = .data$absolute_difference,
                      color = .data$group)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey45", linetype = "dotted") +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::labs(x = "Comparative procedure mean",
                    y = "Evaluated - comparative",
                    color = "Patient group",
                    title = "Patient-specimen difference plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (is.null(interferent)) {
    if (length(x$interferent_cols) != 1L)
      stop("Supply `interferent` when more than one concentration column is available.",
           call. = FALSE)
    interferent <- x$interferent_cols
  }
  if (!interferent %in% x$interferent_cols)
    stop("Requested interferent was not analyzed.", call. = FALSE)
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data[[interferent]],
                 y = .data$absolute_difference, color = .data$group)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey45", linetype = "dotted") +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::labs(x = interferent, y = "Evaluated - comparative",
                  color = "Patient group",
                  title = "Difference vs interferent concentration") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  fit <- x$models[[interferent]]
  if (!is.null(fit)) {
    grid <- data.frame(concentration = seq(min(d[[interferent]], na.rm = TRUE),
                                           max(d[[interferent]], na.rm = TRUE),
                                           length.out = 100L))
    grid$fit <- stats::predict(fit, newdata = grid)
    p <- p + ggplot2::geom_line(data = grid,
      ggplot2::aes(x = .data$concentration, y = .data$fit),
      inherit.aes = FALSE,
      color = "#B2182B", linewidth = 0.8)
  }
  p
}
