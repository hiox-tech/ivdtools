# EP28 simple nonparametric point estimate and rank metadata. R quantile type 6
# uses the same h = p * (n + 1) interpolation rule.
.ri_rank_quantile <- function(sorted_values, p) {
  n <- length(sorted_values)
  rank <- p * (n + 1)
  if (rank <= 1) {
    return(list(value = sorted_values[1L], rank = rank,
                bracket = c(1L, 1L), fraction = 0))
  }
  if (rank >= n) {
    return(list(value = sorted_values[n], rank = rank,
                bracket = c(n, n), fraction = 0))
  }
  lower_rank <- floor(rank)
  upper_rank <- ceiling(rank)
  fraction <- rank - lower_rank
  value <- sorted_values[lower_rank] +
    fraction * (sorted_values[upper_rank] - sorted_values[lower_rank])
  list(value = unname(value), rank = rank,
       bracket = c(lower_rank, upper_rank), fraction = fraction)
}

# Exact binomial order-statistic CI for a population quantile. For n = 120,
# p = .025, and confidence .90 this gives ranks 1 and 7 as in EP28 Table 8.
.ri_binomial_rank_ci <- function(sorted_values, p, ci) {
  n <- length(sorted_values)
  tail <- (1 - ci) / 2
  raw_lower <- stats::qbinom(tail, n, p)
  raw_upper <- stats::qbinom(1 - tail, n, p) + 1
  lower_rank <- max(1L, as.integer(raw_lower))
  upper_rank <- min(n, as.integer(raw_upper))
  actual <- stats::pbinom(upper_rank - 1L, n, p) -
    stats::pbinom(lower_rank - 1L, n, p)
  list(
    lower = sorted_values[lower_rank],
    upper = sorted_values[upper_rank],
    rank = c(lower_rank, upper_rank),
    raw_rank = c(raw_lower, raw_upper),
    requested_coverage = ci,
    actual_coverage = unname(actual),
    attainable = raw_lower >= 1 && raw_upper <= n
  )
}

.ri_normal_rank_ci <- function(sorted_values, p, ci) {
  n <- length(sorted_values)
  z <- stats::qnorm(1 - (1 - ci) / 2)
  lower_rank <- max(1L, floor(n * p - z * sqrt(n * p * (1 - p))))
  upper_rank <- min(n, ceiling(n * p + z * sqrt(n * p * (1 - p))))
  list(lower = sorted_values[lower_rank], upper = sorted_values[upper_rank],
       rank = c(lower_rank, upper_rank))
}

.ri_nonparametric <- function(values, interval, ci) {
  n <- length(values)
  alpha <- (1 - interval) / 2
  sorted <- sort(values)
  lower_point <- .ri_rank_quantile(sorted, alpha)
  upper_point <- .ri_rank_quantile(sorted, 1 - alpha)
  lower_ci <- .ri_binomial_rank_ci(sorted, alpha, ci)
  upper_ci <- .ri_binomial_rank_ci(sorted, 1 - alpha, ci)
  list(
    method = "Nonparametric",
    method_key = "nonparametric",
    standard_basis = "CLSI EP28-A3c simple nonparametric",
    n = n,
    lower = lower_point$value,
    upper = upper_point$value,
    lower_rank = lower_point$rank,
    upper_rank = upper_point$rank,
    lower_bracket = lower_point$bracket,
    upper_bracket = upper_point$bracket,
    lower_fraction = lower_point$fraction,
    upper_fraction = upper_point$fraction,
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    ci_attainable = lower_ci$attainable && upper_ci$attainable,
    quantile_rule = "h = p * (n + 1), linear interpolation",
    ci_method = "exact binomial order-statistic ranks",
    alpha = alpha
  )
}

.ri_type6 <- function(values, interval, ci) {
  n <- length(values)
  alpha <- (1 - interval) / 2
  sorted <- sort(values)
  lower <- stats::quantile(values, alpha, type = 6)
  upper <- stats::quantile(values, 1 - alpha, type = 6)
  lower_ci <- .ri_normal_rank_ci(sorted, alpha, ci)
  upper_ci <- .ri_normal_rank_ci(sorted, 1 - alpha, ci)
  list(
    method = "Type 6",
    method_key = "type6",
    n = n,
    lower = unname(lower),
    upper = unname(upper),
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    quantile_type = 6L,
    quantile_rule = "R quantile type 6",
    ci_method = "normal approximation to binomial ranks",
    alpha = alpha
  )
}

.ri_parametric <- function(values, interval, ci) {
  n <- length(values)
  alpha <- (1 - interval) / 2
  z_value <- stats::qnorm(1 - alpha)
  sw <- stats::shapiro.test(values)
  transformed <- FALSE
  transform_type <- "none"
  working <- values
  if (sw$p.value < 0.05) {
    if (all(values > 0)) {
      log_values <- log(values)
      sw_log <- tryCatch(stats::shapiro.test(log_values),
                         error = function(e) NULL)
      if (!is.null(sw_log) && sw_log$p.value >= 0.05) {
        working <- log_values
        transformed <- TRUE
        transform_type <- "log"
        sw <- sw_log
      } else {
        warning("Data not normal (Shapiro-Wilk p = ", .format_num(sw$p.value),
                "). Log transformation also failed. Proceeding with original ",
                "values.", call. = FALSE)
      }
    } else {
      warning("Data not normal (Shapiro-Wilk p = ", .format_num(sw$p.value),
              ") and log transformation not possible (non-positive values). ",
              "Proceeding with original values.", call. = FALSE)
    }
  }
  center <- mean(working)
  spread <- stats::sd(working)
  lower <- center - z_value * spread
  upper <- center + z_value * spread
  se_limit <- spread * sqrt(1 / n + z_value^2 / (2 * n))
  critical <- stats::qt(1 - (1 - ci) / 2, df = max(n - 2, 1))
  lower_ci <- list(lower = lower - critical * se_limit,
                   upper = lower + critical * se_limit)
  upper_ci <- list(lower = upper - critical * se_limit,
                   upper = upper + critical * se_limit)
  if (transformed) {
    lower <- exp(lower); upper <- exp(upper)
    lower_ci <- lapply(lower_ci, exp)
    upper_ci <- lapply(upper_ci, exp)
  }
  list(
    method = "Parametric",
    method_key = "parametric",
    n = n,
    lower = lower,
    upper = upper,
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    mean = if (transformed) exp(center) else center,
    sd = if (transformed) NA_real_ else spread,
    shapiro_wilk = list(W = unname(sw$statistic), p = sw$p.value),
    transformed = transformed,
    transform_type = transform_type,
    z_value = z_value,
    alpha = alpha
  )
}

.ri_biweight_scale <- function(values, tuning, center, scale,
                               include_n = TRUE) {
  u <- (values - center) / (tuning * scale)
  keep <- is.finite(u) & abs(u) < 1
  if (!any(keep)) return(NA_real_)
  u <- u[keep]
  denominator_sum <- sum((1 - u^2) * (1 - 5 * u^2))
  denominator <- denominator_sum * max(1, -1 + denominator_sum)
  numerator <- sum(u^2 * (1 - u^2)^4)
  if (include_n) numerator <- length(values) * numerator
  if (!is.finite(denominator) || denominator <= 0 || numerator < 0)
    return(NA_real_)
  tuning * scale * sqrt(numerator / denominator)
}

.ri_biweight_core <- function(values, interval, max_iterations = 1000L,
                              tolerance = 1e-5) {
  n <- length(values)
  alpha <- (1 - interval) / 2
  median_value <- stats::median(values)
  initial_scale <- stats::median(abs(values - median_value)) / 0.6745
  if (!is.finite(initial_scale) || initial_scale <= 0)
    stop("Biweight estimation requires a positive MAD/0.6745 scale.",
         call. = FALSE)
  center <- median_value
  converged <- FALSE
  iterations <- 0L
  for (iteration in seq_len(max_iterations)) {
    u <- (values - center) / (3.7 * initial_scale)
    weights <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
    if (!is.finite(sum(weights)) || sum(weights) <= 0)
      stop("Biweight location weights are degenerate.", call. = FALSE)
    new_center <- sum(weights * values) / sum(weights)
    iterations <- iteration
    relative_change <- abs(new_center - center) /
      max(abs(center), .Machine$double.eps)
    center <- new_center
    if (relative_change < tolerance) {
      converged <- TRUE
      break
    }
  }
  s_bi_205_6 <- .ri_biweight_scale(
    values, tuning = 205.6, center = median_value,
    scale = initial_scale, include_n = TRUE
  )
  s_bi_3_7 <- .ri_biweight_scale(
    values, tuning = 3.7, center = median_value,
    scale = initial_scale, include_n = TRUE
  )
  s_t_3_7 <- .ri_biweight_scale(
    values, tuning = 3.7, center = center,
    scale = s_bi_3_7, include_n = FALSE
  )
  components <- c(s_bi_205_6, s_bi_3_7, s_t_3_7)
  if (any(!is.finite(components)) || any(components <= 0))
    stop("Biweight scale components could not be estimated.", call. = FALSE)
  combined_scale <- sqrt(s_bi_205_6^2 + s_t_3_7^2)
  critical <- stats::qt(1 - alpha, df = n - 1)
  list(
    lower = center - critical * combined_scale,
    upper = center + critical * combined_scale,
    center = center,
    initial_median = median_value,
    initial_scale = initial_scale,
    s_bi_205_6 = s_bi_205_6,
    s_bi_3_7 = s_bi_3_7,
    s_t_3_7 = s_t_3_7,
    combined_scale = combined_scale,
    critical_value = critical,
    iterations = iterations,
    converged = converged
  )
}

.ri_huber_core <- function(values, interval, max_iterations = 100L,
                           tolerance = 1e-8) {
  n <- length(values)
  alpha <- (1 - interval) / 2
  center <- stats::median(values)
  spread <- stats::mad(values, constant = 1 / 0.6745)
  if (!is.finite(spread) || spread <= 0)
    stop("Huber estimation requires a positive MAD scale.", call. = FALSE)
  converged <- FALSE
  iterations <- 0L
  for (iteration in seq_len(max_iterations)) {
    u <- (values - center) / spread
    weights <- ifelse(abs(u) <= 1.5, 1, 1.5 / abs(u))
    new_center <- sum(weights * values) / sum(weights)
    new_spread <- sqrt(sum(weights * (values - new_center)^2) /
                         max(n - 1, 1))
    iterations <- iteration
    if (abs(new_center - center) < tolerance &&
        abs(new_spread - spread) < tolerance) {
      center <- new_center
      spread <- new_spread
      converged <- TRUE
      break
    }
    center <- new_center
    spread <- new_spread
  }
  if (!is.finite(spread) || spread <= 0)
    stop("Huber scale could not be estimated.", call. = FALSE)
  z_value <- stats::qnorm(1 - alpha)
  list(lower = center - z_value * spread,
       upper = center + z_value * spread,
       center = center, spread = spread, z_value = z_value,
       iterations = iterations, converged = converged)
}

.ri_bootstrap_ci <- function(values, interval, ci, replicates, seed,
                             estimator, method_name) {
  n <- length(values)
  draws <- .with_seed_preserved(seed, replicate(replicates, {
    sample_values <- values[sample.int(n, n, replace = TRUE)]
    tryCatch({
      estimate <- estimator(sample_values, interval)
      c(estimate$lower, estimate$upper)
    }, error = function(e) c(NA_real_, NA_real_))
  }))
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = 2L)
  valid <- is.finite(draws[1L, ]) & is.finite(draws[2L, ])
  n_valid <- sum(valid)
  if (n_valid < max(100L, ceiling(0.8 * replicates)))
    warning(method_name, " bootstrap produced only ", n_valid, " valid of ",
            replicates, " replicates.", call. = FALSE)
  if (!n_valid) {
    return(list(ci_lower = list(lower = NA_real_, upper = NA_real_),
                ci_upper = list(lower = NA_real_, upper = NA_real_),
                valid_replicates = 0L))
  }
  probs <- c((1 - ci) / 2, 1 - (1 - ci) / 2)
  lower_ci <- stats::quantile(draws[1L, valid], probs, names = FALSE,
                              type = 7)
  upper_ci <- stats::quantile(draws[2L, valid], probs, names = FALSE,
                              type = 7)
  list(ci_lower = list(lower = lower_ci[1L], upper = lower_ci[2L]),
       ci_upper = list(lower = upper_ci[1L], upper = upper_ci[2L]),
       valid_replicates = n_valid)
}

.ri_biweight <- function(values, interval, ci, bootstrap_replicates, seed) {
  estimate <- .ri_biweight_core(values, interval)
  boot <- .ri_bootstrap_ci(values, interval, ci, bootstrap_replicates, seed,
                           .ri_biweight_core, "Biweight")
  c(list(method = "Biweight", method_key = "biweight",
         standard_basis = "CLSI EP28-A3c Appendix B", n = length(values),
         ci_lower = boot$ci_lower, ci_upper = boot$ci_upper,
         bootstrap_reps = bootstrap_replicates,
         bootstrap_valid = boot$valid_replicates,
         seed = seed, alpha = (1 - interval) / 2), estimate)
}

.ri_huber <- function(values, interval, ci, bootstrap_replicates, seed) {
  estimate <- .ri_huber_core(values, interval)
  boot <- .ri_bootstrap_ci(values, interval, ci, bootstrap_replicates, seed,
                           .ri_huber_core, "Huber")
  list(
    method = "Huber",
    method_key = "huber",
    n = length(values),
    lower = estimate$lower,
    upper = estimate$upper,
    ci_lower = boot$ci_lower,
    ci_upper = boot$ci_upper,
    robust_mean = estimate$center,
    robust_sd = estimate$spread,
    algorithm = paste0("Huber M-estimation (k=1.5, ",
                       estimate$iterations, " iter) + MAD"),
    iterations = estimate$iterations,
    converged = estimate$converged,
    bootstrap_reps = bootstrap_replicates,
    bootstrap_valid = boot$valid_replicates,
    seed = seed,
    z_value = estimate$z_value,
    alpha = (1 - interval) / 2
  )
}

#' Reference interval analysis
#'
#' Calculates reference intervals using the CLSI EP28-A3c simple
#' nonparametric or biweight procedures, a normal/log-normal parametric
#' procedure, or the retained R type-6 and Huber procedures.
#'
#' `method = "nonparametric"` uses ranks `p * (n + 1)`, linear interpolation,
#' and exact binomial order-statistic confidence intervals. `method =
#' "biweight"` implements the EP28 Appendix B biweight location and scale
#' calculation and percentile bootstrap confidence intervals. `method =
#' "type6"` retains the previous R type-6 quantiles with normal-approximation
#' rank intervals. `method = "huber"` retains the previous Huber M-estimator.
#'
#' @param data A data frame.
#' @param col A single character string naming the numeric result column.
#' @param interval Central reference interval coverage in `(0, 1)`; the
#'   default is `0.95`.
#' @param ci Confidence level for each reference-limit interval in `(0, 1)`;
#'   the default is `0.95`.
#' @param method One or more of `"nonparametric"` (default), `"biweight"`,
#'   `"type6"`, `"parametric"`, or `"huber"`; `"all"` selects all five.
#' @param id Optional single character string naming an ID column used only
#'   to report duplicate IDs.
#' @param bootstrap_replicates Number of percentile-bootstrap resamples for
#'   `"biweight"` and `"huber"`. Must be an integer of at least 100.
#' @param seed Optional integer seed for biweight and Huber bootstrap
#'   resampling.
#'   The caller's random-number state is preserved.
#'
#' @return An object of class `reference_interval`. Selected results are stored
#'   in correspondingly named elements: `nonparametric`, `biweight`, `type6`,
#'   `parametric`, and `huber`. Data-quality counts and cleaned `values` are
#'   also retained.
#'
#' @references
#' Clinical and Laboratory Standards Institute. *Defining, Establishing, and
#' Verifying Reference Intervals in the Clinical Laboratory; Approved
#' Guideline—Third Edition*. CLSI document EP28-A3c. Wayne, PA: CLSI; 2010.
#'
#' @examples
#' set.seed(17)
#' d <- data.frame(id = seq_len(150), value = rnorm(150, 100, 15))
#' reference_interval(d, "value", id = "id")
#' reference_interval(d, "value", method = c("nonparametric", "parametric"))
#' reference_interval(d, "value", method = "all",
#'                    bootstrap_replicates = 199, seed = 28)
#'
#' @importFrom stats mad median pbinom qbinom qnorm qt quantile sd shapiro.test
#' @importFrom ggplot2 ggplot aes geom_jitter geom_linerange geom_point
#'   geom_errorbar scale_color_manual labs theme_bw theme element_text margin
#'   position_dodge
#' @export
reference_interval <- function(data, col, interval = 0.95, ci = 0.95,
                               method = "nonparametric", id = NULL,
                               bootstrap_replicates = 1999L, seed = NULL) {
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(col) || length(col) != 1L || is.na(col) ||
      !nzchar(col) || !col %in% names(data))
    stop("'col' must be a single column name in 'data'.", call. = FALSE)
  if (!is.numeric(data[[col]]))
    stop("Column '", col, "' must be numeric.", call. = FALSE)
  if (!is.numeric(interval) || length(interval) != 1L ||
      !is.finite(interval) || interval <= 0 || interval >= 1)
    stop("'interval' must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(ci) || length(ci) != 1L || !is.finite(ci) ||
      ci <= 0 || ci >= 1)
    stop("'ci' must be a single number in (0, 1).", call. = FALSE)
  if (!is.null(id) && (!is.character(id) || length(id) != 1L ||
                       is.na(id) || !nzchar(id) || !id %in% names(data)))
    stop("'id' must be NULL or a single column name in 'data'.", call. = FALSE)

  valid_methods <- c("nonparametric", "biweight", "type6", "parametric",
                     "huber")
  if (!is.character(method) || !length(method) || anyNA(method) ||
      any(!nzchar(method)))
    stop("'method' must contain one or more method names.", call. = FALSE)
  if ("all" %in% method) {
    if (!identical(method, "all"))
      stop("Use 'method = \"all\"' alone.", call. = FALSE)
    methods <- valid_methods
  } else {
    invalid <- setdiff(method, valid_methods)
    if (length(invalid))
      stop("Invalid method(s): ", paste(invalid, collapse = ", "),
           ". Valid: ", paste(c(valid_methods, "all"), collapse = ", "),
           call. = FALSE)
    methods <- unique(method)
  }

  if (any(methods %in% c("biweight", "huber"))) {
    if (!is.numeric(bootstrap_replicates) ||
        length(bootstrap_replicates) != 1L ||
        !is.finite(bootstrap_replicates) ||
        bootstrap_replicates != as.integer(bootstrap_replicates) ||
        bootstrap_replicates < 100L)
      stop("'bootstrap_replicates' must be an integer of at least 100.",
           call. = FALSE)
    bootstrap_replicates <- as.integer(bootstrap_replicates)
    if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                          !is.finite(seed) || seed != as.integer(seed)))
      stop("'seed' must be NULL or one finite integer.", call. = FALSE)
    if (!is.null(seed)) seed <- as.integer(seed)
  }

  n_total <- nrow(data)
  complete <- is.finite(data[[col]])
  n_miss <- sum(!complete)
  values <- data[[col]][complete]
  n_complete <- length(values)
  if (n_complete < 3L)
    stop("Need at least 3 complete observations, but only ", n_complete,
         " available.", call. = FALSE)
  if (length(unique(values)) < 2L)
    stop("Zero variance in data; cannot compute reference interval.",
         call. = FALSE)

  n_duplicates <- 0L
  duplicate_rows <- integer(0)
  id_name <- id %||% "--"
  if (!is.null(id)) {
    complete_ids <- data[[id]][complete]
    duplicated_ids <- duplicated(complete_ids) |
      duplicated(complete_ids, fromLast = TRUE)
    n_duplicates <- sum(duplicated_ids)
    duplicate_rows <- which(complete)[duplicated_ids]
  }

  results <- setNames(vector("list", length(valid_methods)), valid_methods)
  if ("nonparametric" %in% methods) {
    if (n_complete < 120L)
      warning("Only ", n_complete, " complete observations; CLSI EP28-A3c ",
              "recommends at least 120 for the simple nonparametric method.",
              call. = FALSE)
    results$nonparametric <- .ri_nonparametric(values, interval, ci)
    if (!results$nonparametric$ci_attainable)
      warning("The requested nonparametric reference-limit confidence ",
              "interval is not fully attainable at n = ", n_complete,
              "; boundary order statistics were used and actual coverage is ",
              "reported.", call. = FALSE)
  }
  if ("biweight" %in% methods)
    results$biweight <- .ri_biweight(
      values, interval, ci, bootstrap_replicates, seed
    )
  if ("type6" %in% methods)
    results$type6 <- .ri_type6(values, interval, ci)
  if ("parametric" %in% methods)
    results$parametric <- .ri_parametric(values, interval, ci)
  if ("huber" %in% methods)
    results$huber <- .ri_huber(
      values, interval, ci, bootstrap_replicates, seed
    )

  out <- c(list(
    call = match.call(),
    col = col,
    id_name = id_name,
    n_total = n_total,
    n_miss = n_miss,
    n_complete = n_complete,
    n_duplicates = n_duplicates,
    dup_rows = duplicate_rows,
    values = values,
    interval = interval,
    ci = ci,
    methods = methods,
    bootstrap_replicates = bootstrap_replicates,
    seed = seed
  ), results)
  class(out) <- "reference_interval"
  out
}


# S3 print method

#' Print reference interval analysis results
#'
#' @param x    A "reference_interval" object.
#' @param ...  Additional arguments (reserved for generic).
#'
#' @return The unchanged \code{reference_interval} object, invisibly.
#' @export
print.reference_interval <- function(x, ...) {
  alpha <- (1 - x$interval) / 2

  cat("\nReference Interval Analysis\n")
  cat(sprintf("  %-20s %s\n", "Column:", x$col))
  cat(sprintf("  %-20s %s%% (alpha = %.3f)\n", "Reference interval:", x$interval * 100, alpha))
  cat(sprintf("  %-20s %s\n", "Confidence level:", x$ci))
  cat(sprintf("  %-20s %s\n", "ID variable:", x$id_name))
  cat(sprintf("  %-20s %d / %d\n", "Complete cases:", x$n_complete, x$n_total))
  cat("\n")

  # Missing values
  if (x$n_miss > 0) {
    cat(sprintf("  Missing values: %d / %d rows excluded\n", x$n_miss, x$n_total))
  } else {
    cat("  No missing values.\n")
  }

  # Duplicate IDs
  if (x$id_name != "--") {
    if (x$n_duplicates > 0) {
      cat(sprintf("  Duplicate IDs: %d found at row(s): %s\n",
                  x$n_duplicates, paste(x$dup_rows, collapse = ", ")))
    } else {
      cat("  No duplicate IDs found.\n")
    }
  } else {
    cat("  Note: 'id' not provided, skip duplicate detection.\n")
  }

  if ("nonparametric" %in% x$methods && x$n_complete < 120L)
    cat(sprintf("  Note: %d observations (CLSI EP28 recommends >= 120 for the simple nonparametric method)\n",
                x$n_complete))

  # Per-method results
  for (method_key in x$methods) {
    res <- x[[method_key]]
    if (is.null(res)) next

    # Method description
    desc <- res$method
    if (method_key == "nonparametric") {
      desc <- "Nonparametric (CLSI EP28-A3c; exact binomial rank CI)"
    } else if (method_key == "biweight") {
      desc <- sprintf("Biweight (CLSI EP28-A3c Appendix B; %d iterations)",
                      res$iterations)
    } else if (method_key == "type6") {
      desc <- sprintf("Type 6 (R quantile type %d; normal-approximation rank CI)",
                      res$quantile_type)
    } else if (res$method == "Parametric") {
      sw <- res$shapiro_wilk
      norm_note <- if (res$transformed) {
        sprintf("log-transformed, Shapiro-Wilk p = %s", formatC(sw$p, format = "f", digits = 4))
      } else if (sw$p >= 0.05) {
        sprintf("Shapiro-Wilk p = %s, normal", formatC(sw$p, format = "f", digits = 4))
      } else {
        "WARNING: non-normal data"
      }
      desc <- sprintf("Parametric (%s)", norm_note)
    } else if (method_key == "huber") {
      desc <- sprintf("Huber (%s)", res$algorithm)
    }

    cat("\n")
    cat(sprintf("  %s\n", desc))
    cat(sprintf("    Lower = %s  [%s, %s]\n",
                formatC(res$lower, format = "f", digits = 4),
                formatC(res$ci_lower$lower, format = "f", digits = 4),
                formatC(res$ci_lower$upper, format = "f", digits = 4)))
    cat(sprintf("    Upper = %s  [%s, %s]\n",
                formatC(res$upper, format = "f", digits = 4),
                formatC(res$ci_upper$lower, format = "f", digits = 4),
                formatC(res$ci_upper$upper, format = "f", digits = 4)))
    if (method_key == "nonparametric") {
      cat(sprintf("    Point ranks: %.4f and %.4f\n",
                  res$lower_rank, res$upper_rank))
      cat(sprintf("    CI ranks: lower %d-%d; upper %d-%d\n",
                  res$ci_lower$rank[1L], res$ci_lower$rank[2L],
                  res$ci_upper$rank[1L], res$ci_upper$rank[2L]))
      cat(sprintf("    Actual rank coverage: lower %.4f; upper %.4f\n",
                  res$ci_lower$actual_coverage,
                  res$ci_upper$actual_coverage))
      if (!isTRUE(res$ci_attainable))
        cat("    Warning: requested rank-CI coverage is not fully attainable at this sample size.\n")
    }
    if (method_key == "biweight") {
      cat(sprintf("    Tbi = %.6f; sbi[205.6] = %.6f; ST[3.7] = %.6f\n",
                  res$center, res$s_bi_205_6, res$s_t_3_7))
      cat(sprintf("    Bootstrap valid: %d / %d\n",
                  res$bootstrap_valid, res$bootstrap_reps))
    } else if (method_key == "huber") {
      cat(sprintf("    Bootstrap valid: %d / %d\n",
                  res$bootstrap_valid, res$bootstrap_reps))
    }
  }
  cat("\n")

  invisible(x)
}


# S3 plot method

#' Plot reference interval results
#'
#' One panel per method: jitter strip chart overlaid with
#' reference interval limits, endpoints, and CI error bars.
#'
#' @param x    A "reference_interval" object.
#' @param ...  Additional arguments (reserved for generic).
#'
#' @return A \code{ggplot} object, invisibly. The plot is also drawn as a side
#'   effect.
#' @examples
#' set.seed(11)
#' d <- data.frame(value = rnorm(120, 100, 15))
#' ri <- reference_interval(d, "value")
#' plot(ri)
#' @export
plot.reference_interval <- function(x, ...) {

  values <- x$values
  col_name <- x$col
  methods <- x$methods

  method_colors <- c(
    nonparametric = "#2166AC",
    biweight      = "#1B7837",
    type6         = "#67A9CF",
    parametric    = "#B2182B",
    huber         = "#762A83"
  )
  color_labels <- c(
    nonparametric = "Nonparametric",
    biweight      = "Biweight",
    type6         = "Type 6",
    parametric    = "Parametric",
    huber         = "Huber"
  )

  # Build method limits data
  limits_list <- list()
  for (m in methods) {
    res <- x[[m]]
    if (is.null(res)) next
    limits_list[[length(limits_list) + 1L]] <- data.frame(
      method      = factor(unname(color_labels[m]), levels = color_labels[methods]),
      lower       = res$lower,
      upper       = res$upper,
      ci_lo_lower = res$ci_lower$lower,
      ci_lo_upper = res$ci_lower$upper,
      ci_up_lower = res$ci_upper$lower,
      ci_up_upper = res$ci_upper$upper,
      stringsAsFactors = FALSE
    )
  }
  limits_df <- do.call(rbind, limits_list)

  # Build data frame with one row per value per method (for jitter)
  data_list <- list()
  for (m in methods) {
    if (is.null(x[[m]])) next
    data_list[[length(data_list) + 1L]] <- data.frame(
      method = factor(unname(color_labels[m]), levels = color_labels[methods]),
      value  = as.numeric(values),
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, data_list)

  pd <- ggplot2::position_dodge(width = 0.25)

  # Single panel: each column is one method, showing jitter + RI limits
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_jitter(ggplot2::aes(x = .data$method, y = .data$value),
      width = 0.15, size = 1.2, alpha = 0.3, color = "grey40") +
    ggplot2::geom_linerange(data = limits_df,
      ggplot2::aes(x = .data$method, ymin = .data$lower, ymax = .data$upper,
                   color = .data$method),
      linewidth = 1, position = pd) +
    ggplot2::geom_point(data = limits_df,
      ggplot2::aes(x = .data$method, y = .data$lower, color = .data$method),
      size = 2.5, position = pd) +
    ggplot2::geom_point(data = limits_df,
      ggplot2::aes(x = .data$method, y = .data$upper, color = .data$method),
      size = 2.5, position = pd) +
    ggplot2::geom_errorbar(data = limits_df,
      ggplot2::aes(x = .data$method, ymin = .data$ci_lo_lower,
                   ymax = .data$ci_lo_upper, color = .data$method),
      width = 0.1, linewidth = 0.4, position = pd) +
    ggplot2::geom_errorbar(data = limits_df,
      ggplot2::aes(x = .data$method, ymin = .data$ci_up_lower,
                   ymax = .data$ci_up_upper, color = .data$method),
      width = 0.1, linewidth = 0.4, position = pd) +
    ggplot2::scale_color_manual(
      values = setNames(method_colors[methods], color_labels[methods]),
      labels = color_labels[methods]
    ) +
    ggplot2::labs(x = NULL, y = col_name,
                  title = paste("Reference Interval -", col_name)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position   = "none",
      plot.title        = ggplot2::element_text(hjust = 0.5),
      plot.margin       = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )

  print(p)
  invisible(p)
}
