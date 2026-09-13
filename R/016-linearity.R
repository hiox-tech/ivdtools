# CLSI EP06 linearity tools.
#
# The straight-line workflow follows CLSI EP06-Ed2 (2020). The optional
# polynomial workflow follows the superseded EP6-A polynomial procedure and
# is retained as an exploratory analysis. Neither workflow supplies a
# clinical acceptance limit or a regulatory conclusion.

utils::globalVariables(".weight")

.linearity_match_col <- function(data, column, argument) {
  if (!is.character(column) || length(column) != 1L || is.na(column)) {
    stop("`", argument, "` must be a single column name.", call. = FALSE)
  }
  if (!column %in% names(data)) {
    stop("Column '", column, "' supplied to `", argument,
         "` was not found in `data`.", call. = FALSE)
  }
  column
}

.linearity_scalar_probability <- function(x, name, open_lower = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      (open_lower && x <= 0) || (!open_lower && x < 0) || x >= 1) {
    bounds <- if (open_lower) "strictly between 0 and 1" else
      "greater than or equal to 0 and less than 1"
    stop("`", name, "` must be ", bounds, ".", call. = FALSE)
  }
  as.numeric(x)
}

.linearity_recycle <- function(x, n, name, allow_null = TRUE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop("`", name, "` must be supplied.", call. = FALSE)
  }
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` must contain finite numeric values.", call. = FALSE)
  }
  if (length(x) == 1L) return(rep(as.numeric(x), n))
  if (length(x) != n) {
    stop("`", name, "` must have length 1 or the number of sample levels (",
         n, ").", call. = FALSE)
  }
  as.numeric(x)
}

.linearity_model_formula <- function(degree, intercept) {
  rhs <- c("x", if (degree >= 2L) paste0("I(x^", 2:degree, ")"))
  stats::as.formula(paste("y ~", if (!intercept) "0 +" else "",
                          paste(rhs, collapse = " + ")))
}

.linearity_term_number <- function(term) {
  if (identical(term, "(Intercept)")) return(0L)
  if (identical(term, "x")) return(1L)
  out <- suppressWarnings(as.integer(sub("^I\\(x\\^([0-9]+)\\)$", "\\1", term)))
  ifelse(is.na(out), NA_integer_, out)
}

.linearity_fit_metrics <- function(fit, degree, weights) {
  sm <- summary(fit)
  cf <- sm$coefficients
  terms <- rownames(cf)
  powers <- vapply(terms, .linearity_term_number, integer(1))
  labels <- ifelse(powers == 0L, "b0", paste0("b", powers))
  coef_table <- data.frame(
    degree = degree,
    term = labels,
    model_term = terms,
    estimate = unname(cf[, 1L]),
    std_error = unname(cf[, 2L]),
    statistic = unname(cf[, 3L]),
    p_value = unname(cf[, 4L]),
    nonlinear = powers >= 2L,
    row.names = NULL,
    check.names = FALSE
  )
  residual <- stats::residuals(fit)
  df <- stats::df.residual(fit)
  rss <- sum(residual^2)
  wrss <- sum(weights * residual^2)
  data.frame(
    degree = degree,
    n = length(residual),
    parameters = length(stats::coef(fit)),
    df_residual = df,
    rss = rss,
    weighted_rss = wrss,
    syx = sqrt(rss / df),
    weighted_syx = sqrt(wrss / df),
    r_squared = unname(sm$r.squared),
    adjusted_r_squared = unname(sm$adj.r.squared),
    row.names = NULL
  ) -> model_table
  list(coefficients = coef_table, model = model_table)
}

.linearity_pool_variance <- function(summary) {
  group <- as.character(summary$variance_group)
  runs <- rle(group)
  if (anyDuplicated(runs$values)) {
    stop("Each `variance_group` must form one contiguous block after ordering ",
         "sample levels by x.", call. = FALSE)
  }
  group_order <- unique(group)
  group_sizes <- table(factor(group, levels = group_order))
  if (any(group_sizes < 2L)) {
    bad <- names(group_sizes)[group_sizes < 2L]
    stop("Each pooled variance group must contain at least two sample levels. ",
         "Problem group(s): ", paste(bad, collapse = ", "), ".",
         call. = FALSE)
  }
  pooled <- lapply(group_order, function(group_name) {
    z <- summary[group == group_name, , drop = FALSE]
    denominator <- sum(z$n - 1)
    variance <- sum((z$n - 1) * z$sd^2) / denominator
    data.frame(variance_group = group_name,
               x_min = min(z$x), x_max = max(z$x),
               pooled_variance = variance,
               pooled_sd = sqrt(variance),
               levels = nrow(z),
               df = denominator)
  })
  pooled <- do.call(rbind, pooled)
  rownames(pooled) <- NULL
  summary$variance_used <- pooled$pooled_variance[
    match(summary$variance_group, pooled$variance_group)]
  list(summary = summary, pooled = pooled)
}

.linearity_adl <- function(levels, predicted, adl, adl_abs, adl_rel) {
  n <- nrow(levels)
  direct <- .linearity_recycle(adl, n, "adl")
  absolute <- .linearity_recycle(adl_abs, n, "adl_abs")
  relative <- .linearity_recycle(adl_rel, n, "adl_rel")
  for (z in list(direct, absolute, relative)) {
    if (!is.null(z) && any(z < 0)) {
      stop("ADL values must be nonnegative.", call. = FALSE)
    }
  }
  candidates <- list()
  if (!is.null(direct)) candidates[[length(candidates) + 1L]] <- direct
  if (!is.null(absolute)) candidates[[length(candidates) + 1L]] <- absolute
  if (!is.null(relative)) {
    candidates[[length(candidates) + 1L]] <- abs(predicted) * relative / 100
  }
  if (!length(candidates)) return(NULL)
  limit <- do.call(pmax, candidates)
  data.frame(adl_lower = -limit, adl_upper = limit, adl = limit,
             row.names = NULL)
}

#' Calculate replicate counts for a linearity study
#'
#' Calculates the minimum number of replicates from the imprecision and
#' allowable deviation from linearity (ADL), following CLSI EP06-Ed2. Relative
#' inputs use percentage units (for example, `5` means 5%). Absolute inputs
#' must use the same measurement unit for `imprecision`, `adl`, and
#' `true_deviation`.
#'
#' @param imprecision Repeatability CV (percent) or SD (absolute units).
#' @param adl Allowable deviation from linearity.
#' @param type Either `"relative"` or `"absolute"`.
#' @param level_probability Desired central probability for an individual
#'   level. The EP06 default is 0.99.
#' @param levels Optional number of concentration levels. Required when
#'   `family_error` is supplied.
#' @param family_error Optional overall risk that at least one level is outside
#'   the ADL. When supplied, it replaces `level_probability`.
#' @param true_deviation Expected true deviation from linearity. It is
#'   subtracted from `adl` before calculating the replicate count.
#' @param min_replicates Smallest returned replicate count.
#'
#' @return A data frame of class `linearity_replicates`.
#' @examples
#' # Relative ADL: both CV and ADL are expressed in percent.
#' linearity_replicates(imprecision = c(3, 5), adl = 5)
#'
#' # Choose the per-level probability from a 10% family-wise failure risk.
#' linearity_replicates(5, 5, levels = 9, family_error = 0.10)
#' @export
linearity_replicates <- function(imprecision, adl,
                                 type = c("relative", "absolute"),
                                 level_probability = 0.99,
                                 levels = NULL, family_error = NULL,
                                 true_deviation = 0,
                                 min_replicates = 2L) {
  type <- match.arg(type)
  args <- list(imprecision = imprecision, adl = adl,
               true_deviation = true_deviation)
  if (any(!vapply(args, is.numeric, logical(1))) ||
      any(vapply(args, function(x) anyNA(x) || any(!is.finite(x)), logical(1)))) {
    stop("`imprecision`, `adl`, and `true_deviation` must be finite numeric values.",
         call. = FALSE)
  }
  n <- max(length(imprecision), length(adl), length(true_deviation))
  lengths <- vapply(args, length, integer(1))
  if (any(!lengths %in% c(1L, n))) {
    stop("`imprecision`, `adl`, and `true_deviation` must each have length 1 ",
         "or the common maximum length.", call. = FALSE)
  }
  imprecision <- rep_len(imprecision, n)
  adl <- rep_len(adl, n)
  true_deviation <- rep_len(true_deviation, n)
  if (any(imprecision <= 0) || any(adl <= 0) || any(true_deviation < 0)) {
    stop("Imprecision and ADL must be positive; true deviation must be nonnegative.",
         call. = FALSE)
  }
  effective_adl <- adl - true_deviation
  if (any(effective_adl <= 0)) {
    stop("`adl - true_deviation` must be positive.", call. = FALSE)
  }
  if (!is.null(family_error)) {
    family_error <- .linearity_scalar_probability(family_error, "family_error")
    if (!is.numeric(levels) || length(levels) != 1L || is.na(levels) ||
        levels < 2 || levels != as.integer(levels)) {
      stop("`levels` must be an integer of at least 2 when `family_error` is used.",
           call. = FALSE)
    }
    level_probability <- (1 - family_error)^(1 / as.integer(levels))
  } else {
    level_probability <- .linearity_scalar_probability(
      level_probability, "level_probability")
  }
  if (!is.numeric(min_replicates) || length(min_replicates) != 1L ||
      is.na(min_replicates) || min_replicates < 1 ||
      min_replicates != as.integer(min_replicates)) {
    stop("`min_replicates` must be a positive integer.", call. = FALSE)
  }
  z <- stats::qnorm((1 + level_probability) / 2)
  raw_n <- (z * imprecision / effective_adl)^2
  replicates <- pmax(as.integer(min_replicates), ceiling(raw_n))
  out <- data.frame(
    type = type,
    imprecision = imprecision,
    adl = adl,
    true_deviation = true_deviation,
    effective_adl = effective_adl,
    level_probability = level_probability,
    z = z,
    calculated_replicates = raw_n,
    replicates = replicates,
    recommend_four = replicates %in% 2:3,
    row.names = NULL
  )
  class(out) <- c("linearity_replicates", "data.frame")
  out
}

#' Design a HIGH/LOW linearity panel
#'
#' @param proportion_high Proportion of the HIGH material at each level,
#'   between 0 and 1.
#' @param high Assigned or measured value of the HIGH material.
#' @param low Assigned or measured value of the LOW material; use zero for a
#'   HIGH/blank design.
#' @param total_volume Optional total volume prepared at each level.
#' @param sample Optional sample labels.
#'
#' @return A data frame containing proportions, expected values, relative
#'   concentrations, and (when requested) component volumes.
#' @examples
#' linearity_panel(
#'   proportion_high = c(1, 0.75, 0.5, 0.25, 0),
#'   high = 210, low = 10, total_volume = 1
#' )
#' @export
linearity_panel <- function(proportion_high, high, low = 0,
                            total_volume = NULL, sample = NULL) {
  if (!is.numeric(proportion_high) || anyNA(proportion_high) ||
      any(!is.finite(proportion_high)) ||
      any(proportion_high < 0 | proportion_high > 1)) {
    stop("`proportion_high` must contain finite values between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(high) || length(high) != 1L || is.na(high) ||
      !is.finite(high) || high <= 0 || !is.numeric(low) ||
      length(low) != 1L || is.na(low) || !is.finite(low) || low < 0 ||
      low >= high) {
    stop("`high` must be positive and greater than the nonnegative `low` value.",
         call. = FALSE)
  }
  n <- length(proportion_high)
  if (is.null(sample)) sample <- paste0("S", seq_len(n))
  if (length(sample) != n || anyNA(sample)) {
    stop("`sample` must have the same length as `proportion_high`.", call. = FALSE)
  }
  expected <- proportion_high * high + (1 - proportion_high) * low
  out <- data.frame(
    sample = sample,
    proportion_high = proportion_high,
    proportion_low = 1 - proportion_high,
    expected = expected,
    relative_concentration = expected / high,
    row.names = NULL,
    check.names = FALSE
  )
  if (!is.null(total_volume)) {
    total_volume <- .linearity_recycle(total_volume, n, "total_volume",
                                       allow_null = FALSE)
    if (any(total_volume <= 0)) stop("`total_volume` must be positive.", call. = FALSE)
    out$total_volume <- total_volume
    out$high_volume <- total_volume * proportion_high
    out$low_volume <- total_volume * (1 - proportion_high)
  }
  class(out) <- c("linearity_panel", "data.frame")
  out
}

#' Calculate target LOW and HIGH concentrations for verification
#'
#' Implements the equations in Appendix J of CLSI EP06-Ed2 for studies in
#' which results outside the measuring interval cannot be obtained.
#'
#' @param lloq Lower limit of quantitation. Supply with `cv_at_lloq` and
#'   `cv_at_k_lloq` to calculate the LOW target.
#' @param uloq Upper limit of quantitation. Supply with `cv_high` to calculate
#'   the HIGH target.
#' @param cv_high Repeatability CV near the upper limit.
#' @param cv_at_lloq CV at the LLoQ.
#' @param cv_at_k_lloq CV at `k * lloq`.
#' @param k Concentration multiple for `cv_at_k_lloq`; EP06 uses 3.
#' @param cv_scale Whether CV inputs are percentages or proportions.
#'
#' @return A list of class `linearity_endpoints`.
#' @examples
#' linearity_endpoints(
#'   lloq = 10, uloq = 1000,
#'   cv_high = 5, cv_at_lloq = 20, cv_at_k_lloq = 15
#' )
#' @export
linearity_endpoints <- function(lloq = NULL, uloq = NULL, cv_high = NULL,
                                cv_at_lloq = NULL, cv_at_k_lloq = NULL,
                                k = 3, cv_scale = c("percent", "proportion")) {
  cv_scale <- match.arg(cv_scale)
  z <- stats::qnorm(0.95)
  convert_cv <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
      stop("`", name, "` must be one positive finite value.", call. = FALSE)
    }
    if (cv_scale == "percent") x / 100 else x
  }
  high_target <- NULL
  if (xor(is.null(uloq), is.null(cv_high))) {
    stop("Supply both `uloq` and `cv_high`, or neither.", call. = FALSE)
  }
  if (!is.null(uloq)) {
    if (!is.numeric(uloq) || length(uloq) != 1L || is.na(uloq) ||
        !is.finite(uloq) || uloq <= 0) stop("`uloq` must be positive.", call. = FALSE)
    qh <- convert_cv(cv_high, "cv_high")
    high_target <- uloq / (1 + z * qh)
  }
  low_args <- c(is.null(lloq), is.null(cv_at_lloq), is.null(cv_at_k_lloq))
  if (any(low_args) && !all(low_args)) {
    stop("Supply `lloq`, `cv_at_lloq`, and `cv_at_k_lloq` together.", call. = FALSE)
  }
  low_target <- a <- b <- NULL
  if (!all(low_args)) {
    if (!is.numeric(lloq) || length(lloq) != 1L || is.na(lloq) ||
        !is.finite(lloq) || lloq <= 0) stop("`lloq` must be positive.", call. = FALSE)
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) || k <= 1) {
      stop("`k` must be one finite value greater than 1.", call. = FALSE)
    }
    q <- convert_cv(cv_at_lloq, "cv_at_lloq")
    p <- convert_cv(cv_at_k_lloq, "cv_at_k_lloq")
    a <- (p * k - q) / (k - 1)
    b <- lloq * k * (q - p) / (k - 1)
    if (1 - z * a <= 0) {
      stop("The supplied low-end precision profile does not yield a finite LOW target.",
           call. = FALSE)
    }
    low_target <- (lloq + z * b) / (1 - z * a)
  }
  if (is.null(high_target) && is.null(low_target)) {
    stop("Supply inputs for at least one endpoint.", call. = FALSE)
  }
  out <- list(
    low_target = low_target,
    high_target = high_target,
    low_profile = if (is.null(low_target)) NULL else c(intercept = b, slope = a),
    one_sided_probability = 0.95,
    z = z,
    cv_scale = cv_scale,
    call = match.call()
  )
  class(out) <- "linearity_endpoints"
  out
}

#' Evaluate linearity using straight-line or polynomial regression
#'
#' `linearity()` is a parameter-driven analysis engine rather than a separate
#' validation/verification workflow. The straight-line calculations implement
#' the approach in CLSI EP06-Ed2. Polynomial regression reproduces the
#' superseded EP6-A procedure and is reported as exploratory. Straight-line
#' models are fitted to sample-level means; EP6-A polynomial models are fitted
#' to the retained replicate observations.
#'
#' With `weights = "pooled"`, variance groups must be supplied explicitly.
#' Within group `g`, the pooled variance is
#' `sum((n[j] - 1) * s[j]^2) / sum(n[j] - 1)`, and the level-mean regression
#' weight is `n[j] / pooled_variance[g]`. The function does not infer groups
#' from the number of concentration levels because EP06 grouping also depends
#' on the observed repeatability pattern.
#'
#' @param data Data frame in long format, with one row per replicate result.
#' @param sample Sample-level identifier column name.
#' @param result Numeric measurement-result column name.
#' @param x Numeric expected value, assigned value, relative concentration, or
#'   HIGH-proportion column name. Its value must be constant within a sample.
#' @param method `"linear"` or `"polynomial"`.
#' @param intercept Include an intercept in every fitted model?
#' @param weights One of `"equal"`, `"level_variance"`, `"pooled"`, or
#'   `"precision_profile"`; alternatively, a numeric vector with one weight per
#'   observation or one weight per sample level (in increasing x order).
#' @param variance_group When `weights = "pooled"`, a single character string
#'   naming the column that assigns every sample level to an explicitly chosen
#'   contiguous variance group. The value must be constant within a sample,
#'   and every group must contain at least two sample levels. It must be `NULL`
#'   for other weighting methods.
#' @param max_degree Highest polynomial degree, 2 or 3.
#' @param poly_alpha Significance level for EP6-A nonlinear coefficients.
#' @param percent_base Denominator for percentage residuals and polynomial
#'   deviation from linearity: the straight-line prediction, x value, or
#'   observed sample mean. EP06-Ed2 examples use the prediction.
#' @param conf.level Optional overall or pointwise confidence level for
#'   straight-line residuals. Use `NULL` to omit confidence intervals.
#' @param multiplicity `"familywise"` adjusts intervals so their joint
#'   coverage is `conf.level`; `"none"` applies `conf.level` to each level.
#' @param adl Optional direct absolute ADL, scalar or one value per level.
#' @param adl_abs Optional fixed absolute ADL.
#' @param adl_rel Optional relative ADL in percent. If more than one ADL form
#'   is supplied, the largest allowable value is used at each level.
#' @param profile_exclude_low Exclude the lowest-x level from the precision
#'   profile and use its observed SD directly?
#' @param profile_intercept Include an intercept in the SD-on-mean precision
#'   profile? EP06 uses a zero intercept.
#'
#' @return An S3 object of class `linearity`.
#' @examples
#' set.seed(2020)
#' linearity_data <- data.frame(
#'   level = rep(1:6, each = 4),
#'   expected = rep(c(10, 30, 50, 70, 90, 110), each = 4)
#' )
#' linearity_data$result <- 2 + 0.98 * linearity_data$expected +
#'   stats::rnorm(nrow(linearity_data), sd = 1.5)
#'
#' # Parameter-driven straight-line analysis.
#' fit_linear <- linearity(
#'   linearity_data, sample = "level", result = "result", x = "expected",
#'   method = "linear", intercept = TRUE, weights = "level_variance",
#'   adl_rel = 5
#' )
#' fit_linear$model_comparison
#' fit_linear$level_summary
#'
#' # EP6-A polynomial analysis (reported as exploratory under EP06-Ed2).
#' fit_polynomial <- linearity(
#'   linearity_data, sample = "level", result = "result", x = "expected",
#'   method = "polynomial", max_degree = 3
#' )
#' fit_polynomial$coefficients
#' @export
linearity <- function(data, sample, result, x,
                      method = c("linear", "polynomial"),
                      intercept = TRUE,
                      weights = c("equal", "level_variance", "pooled",
                                  "precision_profile"),
                      variance_group = NULL,
                      max_degree = 3L, poly_alpha = 0.05,
                      percent_base = c("predicted", "x", "mean"),
                      conf.level = NULL,
                      multiplicity = c("familywise", "none"),
                      adl = NULL, adl_abs = NULL, adl_rel = NULL,
                      profile_exclude_low = TRUE,
                      profile_intercept = FALSE) {
  cl <- match.call()
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  sample <- .linearity_match_col(data, sample, "sample")
  result <- .linearity_match_col(data, result, "result")
  x <- .linearity_match_col(data, x, "x")
  method <- match.arg(method)
  custom_weights <- is.numeric(weights)
  weight_method <- if (custom_weights) "custom" else match.arg(weights)
  if (weight_method == "pooled") {
    if (is.null(variance_group)) {
      stop("`variance_group` must name a data column when ",
           "`weights = \"pooled\"`; pooled groups are not selected ",
           "automatically.", call. = FALSE)
    }
    variance_group <- .linearity_match_col(
      data, variance_group, "variance_group"
    )
  } else if (!is.null(variance_group)) {
    stop("`variance_group` is used only when `weights = \"pooled\"`.",
         call. = FALSE)
  }
  percent_base <- match.arg(percent_base)
  multiplicity <- match.arg(multiplicity)
  if (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept)) {
    stop("`intercept` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(profile_exclude_low) || length(profile_exclude_low) != 1L ||
      is.na(profile_exclude_low) || !is.logical(profile_intercept) ||
      length(profile_intercept) != 1L || is.na(profile_intercept)) {
    stop("Precision-profile options must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(max_degree) || length(max_degree) != 1L || is.na(max_degree) ||
      !max_degree %in% 2:3) stop("`max_degree` must be 2 or 3.", call. = FALSE)
  max_degree <- as.integer(max_degree)
  poly_alpha <- .linearity_scalar_probability(poly_alpha, "poly_alpha")
  if (!is.null(conf.level)) {
    conf.level <- .linearity_scalar_probability(conf.level, "conf.level")
  }

  raw <- data
  missing_rows <- is.na(raw[[sample]]) | is.na(raw[[result]]) | is.na(raw[[x]])
  if (!is.null(variance_group)) {
    missing_rows <- missing_rows | is.na(raw[[variance_group]]) |
      !nzchar(trimws(as.character(raw[[variance_group]])))
  }
  missing_count <- sum(missing_rows)
  if (missing_count > 0L) {
    stop("Found ", missing_count, " rows with missing required values; ",
         "linearity() does not silently exclude or impute observations.", call. = FALSE)
  }
  if (!is.numeric(raw[[result]]) || !is.numeric(raw[[x]]) ||
      any(!is.finite(raw[[result]])) || any(!is.finite(raw[[x]]))) {
    stop("The result and x columns must contain finite numeric values.", call. = FALSE)
  }
  duplicated_rows <- sum(duplicated(raw[, unique(c(sample, result, x)), drop = FALSE]))
  analysis <- data.frame(
    .sample = as.character(raw[[sample]]),
    .result = as.numeric(raw[[result]]),
    .x = as.numeric(raw[[x]]),
    .variance_group = if (is.null(variance_group)) NA_character_ else
      trimws(as.character(raw[[variance_group]])),
    .row = seq_len(nrow(raw)),
    stringsAsFactors = FALSE
  )
  if (nrow(analysis) < 4L) stop("At least four observations are required.", call. = FALSE)
  x_counts <- tapply(analysis$.x, analysis$.sample,
                     function(z) length(unique(z)))
  if (any(x_counts != 1L)) {
    bad <- names(x_counts)[x_counts != 1L]
    stop("The x value must be constant within each sample. Problem sample(s): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  if (weight_method == "pooled") {
    group_counts <- tapply(analysis$.variance_group, analysis$.sample,
                           function(z) length(unique(z)))
    if (any(group_counts != 1L)) {
      bad <- names(group_counts)[group_counts != 1L]
      stop("The variance group must be constant within each sample. Problem ",
           "sample(s): ", paste(bad, collapse = ", "), ".", call. = FALSE)
    }
  }

  split_result <- split(analysis$.result, analysis$.sample)
  level_summary <- data.frame(
    sample = names(split_result),
    x = vapply(names(split_result), function(s) analysis$.x[
      match(s, analysis$.sample)], numeric(1)),
    mean = vapply(split_result, mean, numeric(1)),
    sd = vapply(split_result, stats::sd, numeric(1)),
    n = vapply(split_result, length, integer(1)),
    row.names = NULL,
    check.names = FALSE
  )
  level_summary$cv <- 100 * level_summary$sd / level_summary$mean
  if (weight_method == "pooled") {
    group_map <- setNames(analysis$.variance_group, analysis$.sample)
    level_summary$variance_group <- unname(group_map[level_summary$sample])
  }
  level_summary <- level_summary[order(level_summary$x), , drop = FALSE]
  rownames(level_summary) <- NULL
  n_levels <- nrow(level_summary)
  if (n_levels < 2L) stop("At least two distinct sample levels are required.", call. = FALSE)
  if (method == "polynomial" && n_levels < 5L) {
    stop("EP6-A polynomial evaluation requires at least five sample levels.",
         call. = FALSE)
  }
  if (method == "polynomial" && any(level_summary$n < 2L)) {
    stop("EP6-A polynomial evaluation requires at least two replicates per level.",
         call. = FALSE)
  }

  custom_observation_weights <- FALSE
  precision <- list(method = weight_method, pooled = NULL, profile = NULL)
  if (custom_weights) {
    if (anyNA(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
      stop("Custom weights must be positive finite numbers.", call. = FALSE)
    }
    if (length(weights) == nrow(analysis)) {
      observation_weight <- as.numeric(weights)
      level_fit_weight <- rep(NA_real_, n_levels)
      custom_observation_weights <- TRUE
      level_summary$variance_used <- NA_real_
      level_summary$sigma_used <- NA_real_
    } else if (length(weights) == n_levels) {
      level_weight <- as.numeric(weights)
      names(level_weight) <- level_summary$sample
      observation_weight <- level_weight[analysis$.sample]
      level_fit_weight <- as.numeric(level_weight[level_summary$sample])
      level_summary$variance_used <- NA_real_
      level_summary$sigma_used <- NA_real_
    } else {
      stop("Custom weights must have one value per observation or sample level.",
           call. = FALSE)
    }
  } else if (weight_method == "equal") {
    observation_weight <- rep(1, nrow(analysis))
    level_fit_weight <- rep(1, n_levels)
    level_summary$variance_used <- level_summary$sd^2
    level_summary$sigma_used <- level_summary$sd
  } else if (weight_method == "level_variance") {
    if (any(level_summary$n < 2L) || any(!is.finite(level_summary$sd)) ||
        any(level_summary$sd <= 0)) {
      stop("`level_variance` requires at least two replicates and a positive SD ",
           "at every level.", call. = FALSE)
    }
    level_summary$variance_used <- level_summary$sd^2
    level_summary$sigma_used <- level_summary$sd
    variance_map <- setNames(level_summary$variance_used, level_summary$sample)
    observation_weight <- 1 / variance_map[analysis$.sample]
    level_fit_weight <- level_summary$n / level_summary$variance_used
  } else if (weight_method == "pooled") {
    if (any(level_summary$n < 2L) || any(!is.finite(level_summary$sd))) {
      stop("`pooled` requires at least two replicates at every level.", call. = FALSE)
    }
    pooled <- .linearity_pool_variance(level_summary)
    level_summary <- pooled$summary
    if (any(!is.finite(level_summary$variance_used)) ||
        any(level_summary$variance_used <= 0)) {
      stop("Pooled precision groups must have positive finite variances.",
           call. = FALSE)
    }
    level_summary$sigma_used <- sqrt(level_summary$variance_used)
    variance_map <- setNames(level_summary$variance_used, level_summary$sample)
    observation_weight <- 1 / variance_map[analysis$.sample]
    level_fit_weight <- level_summary$n / level_summary$variance_used
    precision$pooled <- pooled$pooled
    group_weights <- split(
      level_fit_weight,
      factor(level_summary$variance_group,
             levels = precision$pooled$variance_group)
    )
    precision$pooled$fit_weight <- vapply(group_weights, function(z) {
      if (length(unique(z)) == 1L) z[1L] else NA_real_
    }, numeric(1))
    precision$pooled$fit_weight_min <- vapply(group_weights, min, numeric(1))
    precision$pooled$fit_weight_max <- vapply(group_weights, max, numeric(1))
  } else {
    if (any(level_summary$n < 2L) || any(!is.finite(level_summary$sd)) ||
        any(level_summary$sd < 0)) {
      stop("`precision_profile` requires at least two replicates and finite SDs.",
           call. = FALSE)
    }
    profile_data <- level_summary
    excluded <- rep(FALSE, n_levels)
    if (profile_exclude_low) excluded[which.min(profile_data$x)] <- TRUE
    fit_data <- profile_data[!excluded, , drop = FALSE]
    if (nrow(fit_data) < if (profile_intercept) 3L else 2L) {
      stop("Too few levels remain to fit the requested precision profile.",
           call. = FALSE)
    }
    profile_fit <- stats::lm(sd ~ mean, data = fit_data,
                             x = FALSE, y = FALSE,
                             model = TRUE,
                             contrasts = NULL)
    if (!profile_intercept) {
      profile_fit <- stats::lm(sd ~ 0 + mean, data = fit_data,
                               x = FALSE, y = FALSE,
                               model = TRUE,
                               contrasts = NULL)
    }
    sigma <- as.numeric(stats::predict(profile_fit, newdata = profile_data))
    if (profile_exclude_low) sigma[excluded] <- profile_data$sd[excluded]
    if (any(!is.finite(sigma)) || any(sigma <= 0)) {
      stop("The precision profile produced nonpositive or nonfinite sigma values.",
           call. = FALSE)
    }
    level_summary$variance_used <- sigma^2
    level_summary$sigma_used <- sigma
    level_summary$profile_excluded <- excluded
    precision$profile <- list(fit = profile_fit, data = fit_data,
                              excluded_samples = profile_data$sample[excluded])
    sigma_map <- setNames(sigma, level_summary$sample)
    observation_weight <- 1 / sigma_map[analysis$.sample]^2
    level_fit_weight <- 1 / level_summary$sigma_used^2
  }
  observation_weight <- as.numeric(observation_weight)
  analysis$.weight <- observation_weight
  level_summary$fit_weight <- level_fit_weight

  fit_uses_observations <- method == "polynomial" ||
    custom_observation_weights
  if (fit_uses_observations) {
    fit_data <- data.frame(y = analysis$.result, x = analysis$.x,
                           .weight = observation_weight)
  } else {
    fit_data <- data.frame(y = level_summary$mean, x = level_summary$x,
                           .weight = level_fit_weight)
  }
  degrees <- if (method == "linear") 1L else seq_len(max_degree)
  fits <- lapply(degrees, function(d) {
    stats::lm(.linearity_model_formula(d, intercept), data = fit_data,
              weights = .weight, model = TRUE, x = TRUE, y = TRUE)
  })
  names(fits) <- paste0("degree_", degrees)
  singular <- vapply(fits, function(f) anyNA(stats::coef(f)), logical(1))
  if (any(singular)) {
    stop("At least one requested regression model is rank deficient; check the ",
         "number and spacing of x levels.", call. = FALSE)
  }
  metric_parts <- Map(.linearity_fit_metrics, fits, degrees,
                      MoreArgs = list(weights = fit_data$.weight))
  coefficients <- do.call(rbind, lapply(metric_parts, `[[`, "coefficients"))
  model_comparison <- do.call(rbind, lapply(metric_parts, `[[`, "model"))
  rownames(coefficients) <- rownames(model_comparison) <- NULL
  coefficients$significant <- coefficients$nonlinear &
    coefficients$p_value < poly_alpha

  selected_degree <- 1L
  if (method == "polynomial") {
    significant_degree <- unique(coefficients$degree[coefficients$significant])
    if (length(significant_degree)) {
      candidates <- model_comparison[
        model_comparison$degree %in% significant_degree, , drop = FALSE]
      selected_degree <- candidates$degree[
        order(candidates$syx, candidates$degree)][1L]
    }
  }
  selected_fit <- fits[[match(selected_degree, degrees)]]
  linear_fit <- fits[[1L]]

  level_newdata <- data.frame(x = level_summary$x)
  predictions <- lapply(fits, stats::predict, newdata = level_newdata)
  for (i in seq_along(predictions)) {
    level_summary[[paste0("pred_degree_", degrees[i])]] <- as.numeric(predictions[[i]])
  }
  level_summary$pred_linear <- as.numeric(stats::predict(linear_fit,
                                                          newdata = level_newdata))
  level_summary$residual_linear <- level_summary$mean - level_summary$pred_linear
  level_summary$pred_selected <- as.numeric(stats::predict(selected_fit,
                                                            newdata = level_newdata))
  level_summary$residual_selected <- level_summary$mean - level_summary$pred_selected
  level_summary$dl <- level_summary$pred_selected - level_summary$pred_linear
  denominator <- switch(percent_base,
    predicted = level_summary$pred_linear,
    x = level_summary$x,
    mean = level_summary$mean
  )
  level_summary$percent_residual_linear <- 100 * level_summary$residual_linear / denominator
  level_summary$percent_residual_selected <- 100 * level_summary$residual_selected / denominator
  level_summary$percent_dl <- 100 * level_summary$dl / denominator
  if (any(!is.finite(level_summary$percent_residual_linear)) ||
      any(!is.finite(level_summary$percent_residual_selected)) ||
      any(!is.finite(level_summary$percent_dl))) {
    warnings <- "At least one percentage metric is undefined because its denominator is zero."
  } else {
    warnings <- character()
  }

  observation_metrics <- analysis
  observation_metrics$pred_linear <- as.numeric(stats::predict(linear_fit))
  observation_metrics$residual_linear <- stats::residuals(linear_fit)
  observation_metrics$pred_selected <- as.numeric(stats::predict(selected_fit))
  observation_metrics$residual_selected <- stats::residuals(selected_fit)

  if (duplicated_rows > 0L) {
    warnings <- c(warnings, paste0(duplicated_rows,
      " duplicate sample/result/x row(s) were retained."))
  }
  if (method == "polynomial") {
    warnings <- c(warnings,
      "Polynomial regression follows superseded EP6-A and is exploratory under EP06-Ed2.")
  }
  if (any(!is.finite(level_summary$cv))) {
    warnings <- c(warnings,
      "At least one level CV is undefined because its observed mean is zero.")
  }
  if (any(level_summary$sd == 0, na.rm = TRUE)) {
    warnings <- c(warnings,
      "At least one level has zero observed SD; verify the raw replicate results.")
  }

  if (!is.null(conf.level)) {
    sigma <- level_summary$sigma_used
    if (is.null(sigma) || any(!is.finite(sigma)) || any(sigma <= 0)) {
      stop("Confidence intervals require finite level sigma values; custom ",
           "weights do not define sigma.", call. = FALSE)
    }
    if (multiplicity == "familywise") {
      point_confidence <- conf.level^(1 / n_levels)
    } else {
      point_confidence <- conf.level
    }
    z <- stats::qnorm((1 + point_confidence) / 2)
    half_width <- z * sigma / sqrt(level_summary$n)
    if (method == "polynomial") {
      level_summary$linear_ci_lower <- level_summary$residual_linear - half_width
      level_summary$linear_ci_upper <- level_summary$residual_linear + half_width
      warnings <- c(warnings,
        "Confidence intervals describe observed mean minus the straight-line prediction, not polynomial DL.")
    } else {
      level_summary$ci_lower <- level_summary$residual_linear - half_width
      level_summary$ci_upper <- level_summary$residual_linear + half_width
    }
  } else {
    point_confidence <- z <- NULL
  }

  adl_table <- .linearity_adl(level_summary, level_summary$pred_linear,
                              adl, adl_abs, adl_rel)
  if (!is.null(adl_table)) {
    level_summary <- cbind(level_summary, adl_table)
    metric <- if (method == "polynomial") level_summary$dl else
      level_summary$residual_linear
    level_summary$adl_metric <- metric
    level_summary$adl_margin <- level_summary$adl - abs(metric)
    level_summary$adl_exceeded <- abs(metric) > level_summary$adl
    if (!is.null(conf.level) && method == "linear") {
      level_summary$ci_overlaps_adl <-
        level_summary$ci_lower <= level_summary$adl_upper &
        level_summary$ci_upper >= level_summary$adl_lower
    }
  }

  out <- list(
    call = cl,
    data = raw,
    columns = c(sample = sample, result = result, x = x,
                if (!is.null(variance_group))
                  c(variance_group = variance_group)),
    method = method,
    settings = list(
      intercept = intercept,
      weights = weight_method,
      variance_group = variance_group,
      fit_unit = if (fit_uses_observations) "observation" else "level_mean",
      max_degree = max_degree,
      poly_alpha = poly_alpha,
      percent_base = percent_base,
      conf.level = conf.level,
      multiplicity = multiplicity,
      point_confidence = point_confidence,
      z = z,
      profile_exclude_low = profile_exclude_low,
      profile_intercept = profile_intercept
    ),
    n = nrow(analysis),
    levels = n_levels,
    duplicated_rows = duplicated_rows,
    level_summary = level_summary,
    observation_metrics = observation_metrics,
    precision = precision,
    fits = fits,
    linear_fit = linear_fit,
    selected_fit = selected_fit,
    selected_degree = selected_degree,
    coefficients = coefficients,
    model_comparison = model_comparison,
    warnings = unique(warnings)
  )
  class(out) <- "linearity"
  out
}

#' @rdname linearity_replicates
#' @param x A `linearity_replicates` object.
#' @param ... Reserved arguments.
#' @return The unchanged `x`, invisibly.
#' @export
print.linearity_replicates <- function(x, ...) {
  .print_kv_rows(x, list(
    "Design inputs" = c(
      "Imprecision" = "imprecision", "ADL" = "adl",
      "True deviation" = "true_deviation", "Effective ADL" = "effective_adl",
      "Level probability" = "level_probability"),
    "Replicate calculation" = c(
      "Z" = "z", "Calculated replicates" = "calculated_replicates",
      "Replicates" = "replicates", "Recommend four" = "recommend_four")
  ), row_title = c("Type" = "type"))
  if (any(x$recommend_four)) {
    cat("Note: EP06 recommends at least four replicates when possible.\n")
  }
  invisible(x)
}

#' @rdname linearity_panel
#' @param x A `linearity_panel` object.
#' @param ... Reserved arguments.
#' @return The unchanged `x`, invisibly.
#' @export
print.linearity_panel <- function(x, ...) {
  sections <- list(
    "Concentration design" = c(
      "Proportion high" = "proportion_high", "Proportion low" = "proportion_low",
      "Expected" = "expected", "Relative concentration" = "relative_concentration")
  )
  volume_columns <- intersect(c("total_volume", "high_volume", "low_volume"),
                              names(x))
  if (length(volume_columns) > 1L) {
    sections[["Preparation volumes"]] <- setNames(volume_columns, c(
      "Total volume", "High volume", "Low volume")[match(volume_columns,
                                                            c("total_volume", "high_volume", "low_volume"))])
  }
  .print_kv_rows(x, sections, row_title = c("Sample" = "sample"))
  invisible(x)
}

#' @rdname linearity_endpoints
#' @param x A `linearity_endpoints` object.
#' @param ... Reserved arguments.
#' @return The unchanged `x`, invisibly.
#' @export
print.linearity_endpoints <- function(x, ...) {
  cat("CLSI EP06-Ed2 target concentrations\n")
  if (!is.null(x$low_target)) cat("  LOW target: ", format(x$low_target), "\n", sep = "")
  if (!is.null(x$high_target)) cat("  HIGH target: ", format(x$high_target), "\n", sep = "")
  invisible(x)
}

#' Print a linearity analysis
#'
#' @param x A `linearity` object.
#' @param ... Reserved arguments.
#' @return The unchanged `x`, invisibly.
#' @export
print.linearity <- function(x, ...) {
  print_coefficients <- function(coefficients, degree, selected = FALSE) {
    powers <- vapply(
      coefficients$model_term, .linearity_term_number, integer(1)
    )
    coefficient_names <- ifelse(
      is.na(powers), coefficients$model_term, letters[powers + 1L]
    )
    coefficient_matrix <- as.matrix(coefficients[c(
      "estimate", "std_error", "statistic", "p_value"
    )])
    colnames(coefficient_matrix) <- c(
      "Estimate", "Std. Error", "t value", "Pr(>|t|)"
    )
    rownames(coefficient_matrix) <- coefficient_names

    heading <- if (x$method == "polynomial") {
      paste0(
        "Regression coefficients (degree ", degree,
        if (selected) ", selected" else "", ")"
      )
    } else {
      "Regression coefficients"
    }
    cat("\n", heading, "\n\n", sep = "")
    stats::printCoefmat(
      coefficient_matrix,
      digits = 6,
      signif.stars = TRUE,
      signif.legend = TRUE
    )
  }

  sample_weights <- function() {
    if (identical(x$settings$fit_unit, "level_mean")) {
      return(as.numeric(x$level_summary$fit_weight))
    }

    observation_weights <- split(
      x$observation_metrics$.weight,
      x$observation_metrics$.sample
    )
    vapply(x$level_summary$sample, function(sample_name) {
      values <- unique(observation_weights[[sample_name]])
      values <- values[is.finite(values)]
      if (length(values) == 1L) values else NA_real_
    }, numeric(1))
  }

  print_table <- function(table, digits = 6L) {
    numeric_columns <- vapply(table, is.numeric, logical(1))
    formatted <- lapply(table, function(column) {
      if (is.numeric(column)) {
        format(column, digits = digits, trim = TRUE, scientific = NA)
      } else {
        as.character(column)
      }
    })
    widths <- vapply(seq_along(formatted), function(i) {
      max(nchar(c(names(table)[[i]], formatted[[i]])), na.rm = TRUE)
    }, integer(1))
    format_cell <- function(value, width, right) {
      value[is.na(value)] <- "NA"
      if (right) {
        sprintf("%*s", width, value)
      } else {
        sprintf("%-*s", width, value)
      }
    }
    header <- vapply(seq_along(formatted), function(i) {
      format_cell(names(table)[[i]], widths[[i]], numeric_columns[[i]])
    }, character(1))
    cat(paste(header, collapse = "  "), "\n", sep = "")
    for (row in seq_len(nrow(table))) {
      values <- vapply(seq_along(formatted), function(i) {
        format_cell(formatted[[i]][[row]], widths[[i]], numeric_columns[[i]])
      }, character(1))
      cat(paste(values, collapse = "  "), "\n", sep = "")
    }
    invisible(table)
  }

  print_key_value_row <- function(row, fields, row_label, row_column) {
    labels <- names(fields)
    columns <- unname(fields)
    width <- max(nchar(labels))
    cat(row_label, ": ", .print_value_text(row[[row_column]][[1L]]),
        "\n", sep = "")
    for (i in seq_along(columns)) {
      cat(sprintf("%-*s", width, labels[[i]]), " : ",
          .print_value_text(row[[columns[[i]]]][[1L]]), "\n", sep = "")
    }
    invisible(row)
  }

  cat("Linearity regression analysis\n\n")
  cat("Analysis settings\n\n")
  cat("Method: ", x$method,
      if (x$method == "polynomial") " (EP6-A, exploratory)" else " (EP06-Ed2)",
      "\n", sep = "")
  cat("Observations / levels: ", x$n, " / ", x$levels, "\n", sep = "")
  cat("Intercept: ", if (x$settings$intercept) "estimated" else "fixed at zero",
      "\n", sep = "")
  cat("Weights: ", x$settings$weights, "\n", sep = "")
  if (identical(x$settings$weights, "pooled")) {
    cat("Variance group: ", x$settings$variance_group, "\n", sep = "")
  }
  cat("Fit unit: ", x$settings$fit_unit, "\n", sep = "")
  if (x$method == "polynomial") {
    cat("Selected degree: ", x$selected_degree, "\n", sep = "")
  }

  cat("\nModel comparison\n\n")
  comparison_fields <- c(
    "Observations" = "n",
    "Parameters" = "parameters",
    "Residual degrees of freedom" = "df_residual",
    "RSS" = "rss",
    "Weighted RSS" = "weighted_rss",
    "SY.X" = "syx",
    "Weighted SY.X" = "weighted_syx",
    "R-squared" = "r_squared",
    "Adjusted R-squared" = "adjusted_r_squared"
  )
  for (i in seq_len(nrow(x$model_comparison))) {
    comparison_row <- x$model_comparison[i, , drop = FALSE]
    fields <- comparison_fields
    if (x$method == "polynomial" &&
        identical(as.integer(comparison_row$degree),
                  as.integer(x$selected_degree))) {
      comparison_row$selected <- "Yes"
      fields <- c(fields, "Selected" = "selected")
    }
    print_key_value_row(
      comparison_row, fields,
      row_label = "Degree", row_column = "degree"
    )
    if (i < nrow(x$model_comparison)) cat("\n")
  }

  for (i in seq_along(x$fits)) {
    degree <- x$model_comparison$degree[[i]]
    print_coefficients(
      x$coefficients[x$coefficients$degree == degree, , drop = FALSE],
      degree,
      selected = identical(as.integer(degree), as.integer(x$selected_degree))
    )
  }

  if (!is.null(x$precision$pooled)) {
    cat("\nPooled variance groups\n\n")
    pooled_fields <- c(
      "Levels" = "levels",
      "Degrees of freedom" = "df",
      "Minimum x" = "x_min",
      "Maximum x" = "x_max",
      "Pooled variance" = "pooled_variance",
      "Pooled SD" = "pooled_sd",
      "Level-mean weight" = "fit_weight",
      "Minimum weight" = "fit_weight_min",
      "Maximum weight" = "fit_weight_max"
    )
    for (i in seq_len(nrow(x$precision$pooled))) {
      print_key_value_row(
        x$precision$pooled[i, , drop = FALSE], pooled_fields,
        row_label = "Variance group", row_column = "variance_group"
      )
      if (i < nrow(x$precision$pooled)) cat("\n")
    }
  }

  level_summary <- x$level_summary
  predicted <- if (x$method == "polynomial") {
    level_summary$pred_selected
  } else {
    level_summary$pred_linear
  }
  residual <- level_summary$mean - predicted
  percent_residual <- if (x$method == "polynomial") {
    level_summary$percent_residual_selected
  } else {
    level_summary$percent_residual_linear
  }
  recovery_bias <- level_summary$mean - level_summary$x
  percent_recovery_bias <- ifelse(
    level_summary$x == 0,
    NA_real_,
    100 * recovery_bias / level_summary$x
  )
  adl <- if ("adl" %in% names(level_summary)) {
    level_summary$adl
  } else {
    rep(NA_real_, nrow(level_summary))
  }
  decision <- if ("adl_exceeded" %in% names(level_summary)) {
    ifelse(level_summary$adl_exceeded, "Exceeds ADL", "Within ADL")
  } else {
    rep("Not assessed", nrow(level_summary))
  }

  fit_table <- data.frame(
    Sample = level_summary$sample,
    Expected = level_summary$x,
    Predicted = predicted,
    `Observed mean` = level_summary$mean,
    Weight = sample_weights(),
    `Observed SD` = level_summary$sd,
    n = level_summary$n,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  bias_table <- data.frame(
    Sample = level_summary$sample,
    Residual = residual,
    `Recovery bias` = recovery_bias,
    `Residual (%)` = percent_residual,
    `Recovery bias (%)` = percent_recovery_bias,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (x$method == "polynomial") {
    bias_table$DL <- level_summary$dl
  }
  bias_table$ADL <- adl
  bias_table$Decision <- decision

  cat("\nSample fit information\n\n")
  print_table(fit_table)
  cat("\nSample bias information\n\n")
  print_table(bias_table)

  cat("\nDefinitions\n\n")
  cat("Residual = observed mean - predicted value.\n")
  cat("Recovery bias = observed mean - expected value.\n")
  cat("Recovery bias (%) = 100 * recovery bias / expected value.\n")
  cat("Residual (%) denominator: ", x$settings$percent_base, ".\n", sep = "")
  if (x$method == "polynomial") {
    cat("ADL decision metric = deviation from linearity (DL).\n")
  } else {
    cat("ADL decision metric = residual.\n")
  }
  if (any(level_summary$x == 0)) {
    cat("Recovery bias (%) is NA when the expected value is zero.\n")
  }
  if (identical(x$settings$fit_unit, "observation") &&
      anyNA(fit_table$Weight)) {
    cat("Weight is NA when observation weights vary within a sample.\n")
  }

  if (length(x$warnings)) {
    cat("\nWarnings / interpretation notes\n\n")
    cat(paste(x$warnings, collapse = "\n\n"), "\n")
  }
  invisible(x)
}

#' Plot a linearity analysis
#'
#' @param x A `linearity` object.
#' @param type One of `"fit"`, `"replicates"`, `"precision"`,
#'   `"residual"`, or `"difference"`.
#' @param ... Reserved for future use.
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(6)
#' d <- data.frame(level = rep(1:5, each = 3), x = rep(1:5, each = 3))
#' d$y <- 1 + 2 * d$x + stats::rnorm(nrow(d), sd = 0.2)
#' fit <- linearity(d, "level", "y", "x")
#' plot(fit, type = "fit")
#' plot(fit, type = "residual")
#' @export
plot.linearity <- function(x, type = c("fit", "replicates", "precision",
                                       "residual", "difference"), ...) {
  type <- match.arg(type)
  .require_pkg("ggplot2")
  lev <- x$level_summary
  obs <- x$observation_metrics
  if (type == "replicates") {
    return(ggplot2::ggplot(obs, ggplot2::aes(x = .data$.x, y = .data$.result,
                                             color = .data$.sample)) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::labs(x = x$columns[["x"]], y = x$columns[["result"]],
                    color = x$columns[["sample"]],
                    title = "Linearity replicates") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (type == "precision") {
    if (identical(x$settings$weights, "pooled")) {
      p <- ggplot2::ggplot(
        lev,
        ggplot2::aes(x = .data$mean, y = .data$sd,
                     color = .data$variance_group)
      ) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::geom_line(
          ggplot2::aes(y = .data$sigma_used,
                       group = .data$variance_group),
          linewidth = 0.8, linetype = 2
        ) +
        ggplot2::labs(x = "Observed mean", y = "SD",
                      color = x$settings$variance_group,
                      title = "Observed and pooled precision") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      p <- ggplot2::ggplot(lev, ggplot2::aes(x = .data$mean, y = .data$sd)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::labs(x = "Observed mean", y = "SD",
                      title = "Linearity precision profile") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    if (!is.null(x$precision$profile)) {
      grid <- data.frame(mean = seq(min(lev$mean), max(lev$mean), length.out = 200L))
      grid$sd <- stats::predict(x$precision$profile$fit, newdata = grid)
      p <- p + ggplot2::geom_line(data = grid,
                                  ggplot2::aes(x = .data$mean, y = .data$sd),
                                  color = "#0072B2")
    }
    return(p)
  }
  if (type == "fit") {
    grid <- data.frame(x = seq(min(lev$x), max(lev$x), length.out = 300L))
    degree_values <- as.integer(sub("degree_", "", names(x$fits)))
    curves <- do.call(rbind, lapply(seq_along(x$fits), function(i) {
      data.frame(x = grid$x,
                 fitted = as.numeric(stats::predict(x$fits[[i]], newdata = grid)),
                 model = paste0("Degree ", degree_values[i]))
    }))
    return(ggplot2::ggplot() +
      ggplot2::geom_point(data = obs,
                          ggplot2::aes(x = .data$.x, y = .data$.result),
                          alpha = 0.35) +
      ggplot2::geom_point(data = lev,
                          ggplot2::aes(x = .data$x, y = .data$mean), size = 2.5) +
      ggplot2::geom_line(data = curves,
                         ggplot2::aes(x = .data$x, y = .data$fitted,
                                      color = .data$model),
                         linewidth = 0.8) +
      ggplot2::labs(x = x$columns[["x"]], y = x$columns[["result"]],
                    color = "Model", title = "Linearity fitted models") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (type == "residual") {
    metric <- if (x$method == "polynomial") "residual_selected" else "residual_linear"
    d <- data.frame(x = lev$x, residual = lev[[metric]])
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x,
                                         y = .data$residual)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey50") +
      ggplot2::geom_line() + ggplot2::geom_point(size = 2.5) +
      ggplot2::labs(x = x$columns[["x"]], y = "Observed mean - fitted",
                    title = "Linearity residuals") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    if (all(c("ci_lower", "ci_upper") %in% names(lev))) {
      d$lower <- lev$ci_lower
      d$upper <- lev$ci_upper
      p <- p + ggplot2::geom_errorbar(data = d,
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper), width = 0)
    }
    return(p)
  }
  if (x$method != "polynomial") {
    stop("`type = \"difference\"` is available only for polynomial analyses.",
         call. = FALSE)
  }
  d <- data.frame(x = lev$x, dl = lev$dl)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$dl)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = x$columns[["x"]],
                  y = "Selected polynomial - straight line",
                  title = "Linearity difference plot") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (all(c("adl_lower", "adl_upper") %in% names(lev))) {
    bounds <- data.frame(x = lev$x, lower = lev$adl_lower, upper = lev$adl_upper)
    p <- p + ggplot2::geom_ribbon(data = bounds,
      ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
      inherit.aes = FALSE, alpha = 0.15, fill = "#009E73")
  }
  p
}
