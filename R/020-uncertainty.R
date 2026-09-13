# Measurement-uncertainty tools based on CLSI EP29-A, EP30-A, and EP32-R.
#
# The functions estimate and propagate uncertainty. They do not compare the
# results with analytical or clinical acceptance limits.

utils::globalVariables(c(".name", ".contribution_percent", ".stage",
                         ".cumulative", ".simulated"))

.unc_probability <- function(x, argument = "coverage") {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop("`", argument, "` must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(x)
}

.unc_scalar <- function(x, argument, positive = FALSE, nonnegative = FALSE,
                        infinite = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
    (is.finite(x) || (infinite && is.infinite(x)))
  if (positive) valid <- valid && x > 0
  if (nonnegative) valid <- valid && x >= 0
  if (!valid) {
    qualifier <- if (positive) "positive " else if (nonnegative) "nonnegative " else ""
    stop("`", argument, "` must be one ", qualifier,
         if (infinite) "numeric" else "finite numeric", " value.",
         call. = FALSE)
  }
  as.numeric(x)
}

.unc_name <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be one nonempty character value.", call. = FALSE)
  }
  name
}

.unc_component <- function(name, estimate, standard_uncertainty, df = Inf,
                           evaluation_type, distribution, source,
                           relative = FALSE, parameters = list()) {
  name <- .unc_name(name)
  estimate <- .unc_scalar(estimate, "estimate")
  standard_uncertainty <- .unc_scalar(
    standard_uncertainty, "standard_uncertainty", nonnegative = TRUE)
  df <- .unc_scalar(df, "df", positive = TRUE, infinite = TRUE)
  if (!is.logical(relative) || length(relative) != 1L || is.na(relative)) {
    stop("`relative` must be TRUE or FALSE.", call. = FALSE)
  }
  structure(list(
    name = name,
    estimate = estimate,
    standard_uncertainty = standard_uncertainty,
    relative_uncertainty = if (relative) standard_uncertainty else
      if (abs(estimate) > .Machine$double.eps^0.5)
        standard_uncertainty / abs(estimate) else NA_real_,
    df = df,
    evaluation_type = evaluation_type,
    distribution = distribution,
    source = source,
    relative = relative,
    parameters = parameters
  ), class = "uncertainty_component")
}

.unc_component_row <- function(x) {
  data.frame(
    name = x$name,
    estimate = x$estimate,
    standard_uncertainty = x$standard_uncertainty,
    relative_uncertainty = x$relative_uncertainty,
    df = x$df,
    evaluation_type = x$evaluation_type,
    distribution = x$distribution,
    source = x$source,
    relative = x$relative,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Print an uncertainty component
#'
#' @param x An uncertainty_component object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_component <- function(x, ...) {
  cat("\nUncertainty component\n")
  .print_kv_sections(.unc_component_row(x), list(
    "Estimate" = c(
      "Name" = "name", "Estimate" = "estimate",
      "Standard uncertainty" = "standard_uncertainty",
      "Relative uncertainty" = "relative_uncertainty",
      "Degrees of freedom" = "df"),
    "Component metadata" = c(
      "Evaluation type" = "evaluation_type", "Distribution" = "distribution",
      "Source" = "source", "Relative" = "relative")
  ))
  invisible(x)
}

#' Estimate a Type A uncertainty component
#'
#' Uses the sample standard deviation for a future single result and the
#' standard error for a mean. No normality or fitness-for-purpose decision is
#' made.
#'
#' @param x Numeric repeated measurements.
#' @param name Component name.
#' @param quantity `"single"` or `"mean"`.
#' @param relative Store the standard uncertainty relative to the absolute mean.
#' @return An `uncertainty_component` object.
#' @references CLSI EP29-A.
#' @examples
#' type_a <- uncertainty_type_a(
#'   c(10.1, 9.9, 10.0, 10.2, 9.8),
#'   name = "repeatability", quantity = "mean"
#' )
#' type_a
#' @export
uncertainty_type_a <- function(x, name,
                               quantity = c("single", "mean"),
                               relative = FALSE) {
  quantity <- match.arg(quantity)
  name <- .unc_name(name)
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  x <- x[is.finite(x)]
  if (length(x) < 2L) stop("At least two finite measurements are required.",
                           call. = FALSE)
  estimate <- mean(x)
  uncertainty <- stats::sd(x) / if (quantity == "mean") sqrt(length(x)) else 1
  if (relative) {
    if (abs(estimate) <= .Machine$double.eps^0.5) {
      stop("A relative component requires a nonzero estimate.", call. = FALSE)
    }
    stored <- uncertainty / abs(estimate)
  } else {
    stored <- uncertainty
  }
  .unc_component(
    name, estimate, stored, df = length(x) - 1L,
    evaluation_type = "A", distribution = "normal",
    source = paste0("replicate measurements: ", quantity), relative = relative,
    parameters = list(mean = estimate, sd = uncertainty,
                      n = length(x), values = x, quantity = quantity)
  )
}

#' Estimate a Type B uncertainty component
#'
#' Converts limits, an expanded uncertainty, or an already standardized value
#' to standard uncertainty. For bounded distributions `lower` and `upper`
#' define the interval. For normal and Student t distributions they define a
#' central interval with probability `coverage`.
#'
#' @param name Component name.
#' @param estimate Best estimate of the input quantity.
#' @param lower,upper Optional interval limits.
#' @param uncertainty Optional standard or expanded uncertainty.
#' @param distribution Distribution model.
#' @param coverage Central probability associated with normal or t limits.
#' @param k Coverage factor when `uncertainty` is expanded.
#' @param df Degrees of freedom for Student t or the component.
#' @param relative Store uncertainty relative to the absolute estimate.
#' @return An `uncertainty_component` object.
#' @references CLSI EP29-A.
#' @examples
#' rectangular <- uncertainty_type_b(
#'   name = "calibrator", estimate = 100, lower = 98, upper = 102,
#'   distribution = "rectangular"
#' )
#' rectangular
#'
#' uncertainty_type_b(
#'   "certificate", estimate = 100, uncertainty = 2, k = 2,
#'   distribution = "normal"
#' )
#' @export
uncertainty_type_b <- function(name, estimate, lower = NULL, upper = NULL,
                               uncertainty = NULL,
                               distribution = c("rectangular", "triangular",
                                                "normal", "student_t",
                                                "u_shaped", "standard"),
                               coverage = 0.95, k = NULL, df = Inf,
                               relative = FALSE) {
  distribution <- match.arg(distribution)
  name <- .unc_name(name)
  estimate <- .unc_scalar(estimate, "estimate")
  df <- .unc_scalar(df, "df", positive = TRUE, infinite = TRUE)
  coverage <- .unc_probability(coverage)
  limits_given <- !is.null(lower) || !is.null(upper)
  if (limits_given) {
    if (is.null(lower) || is.null(upper)) {
      stop("Supply both `lower` and `upper`.", call. = FALSE)
    }
    lower <- .unc_scalar(lower, "lower")
    upper <- .unc_scalar(upper, "upper")
    if (lower >= upper) stop("`lower` must be smaller than `upper`.", call. = FALSE)
    if (estimate < lower || estimate > upper) {
      stop("`estimate` must lie between `lower` and `upper`.", call. = FALSE)
    }
    midpoint <- (lower + upper) / 2
    if (abs(estimate - midpoint) >
        .Machine$double.eps^0.5 * max(1, abs(midpoint))) {
      stop("`estimate` must be the midpoint of `lower` and `upper`.",
           call. = FALSE)
    }
    half_width <- (upper - lower) / 2
    divisor <- switch(
      distribution,
      rectangular = sqrt(3),
      triangular = sqrt(6),
      u_shaped = sqrt(2),
      normal = stats::qnorm((1 + coverage) / 2),
      student_t = {
        if (!is.finite(df)) stop("Student t limits require finite `df`.",
                                 call. = FALSE)
        stats::qt((1 + coverage) / 2, df = df)
      },
      standard = stop("Use `uncertainty` with `distribution = 'standard'`.",
                      call. = FALSE)
    )
    standard <- half_width / divisor
    source <- "interval limits"
  } else {
    if (is.null(uncertainty)) {
      stop("Supply interval limits or `uncertainty`.", call. = FALSE)
    }
    uncertainty <- .unc_scalar(uncertainty, "uncertainty", nonnegative = TRUE)
    if (!is.null(k)) {
      k <- .unc_scalar(k, "k", positive = TRUE)
      standard <- uncertainty / k
      source <- "expanded uncertainty"
    } else {
      standard <- uncertainty
      source <- "standard uncertainty"
    }
    lower <- upper <- NA_real_
    half_width <- NA_real_
    if (distribution != "standard" && distribution != "normal") {
      stop("Without limits, `distribution` must be 'standard' or 'normal'.",
           call. = FALSE)
    }
  }
  if (relative) {
    if (abs(estimate) <= .Machine$double.eps^0.5) {
      stop("A relative component requires a nonzero estimate.", call. = FALSE)
    }
    stored <- standard / abs(estimate)
  } else {
    stored <- standard
  }
  .unc_component(
    name, estimate, stored, df = df, evaluation_type = "B",
    distribution = distribution, source = source, relative = relative,
    parameters = list(lower = lower, upper = upper,
                      half_width = half_width, coverage = coverage,
                      k = k, absolute_standard_uncertainty = standard)
  )
}

.unc_col <- function(data, column, argument, numeric = FALSE,
                     allow_null = FALSE) {
  if (is.null(column) && allow_null) return(NULL)
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !column %in% names(data)) {
    stop("`", argument, "` must be one column name in `data`.", call. = FALSE)
  }
  if (numeric && !is.numeric(data[[column]])) {
    stop("Column '", column, "' must be numeric.", call. = FALSE)
  }
  column
}

.unc_anova <- function(value, group) {
  keep <- is.finite(value) & !is.na(group)
  value <- as.numeric(value[keep])
  group <- as.character(group[keep])
  if (length(value) < 2L || length(unique(group)) < 2L) {
    stop("At least two groups and two finite observations are required.",
         call. = FALSE)
  }
  split_values <- split(value, group, drop = TRUE)
  sizes <- vapply(split_values, length, integer(1))
  means <- vapply(split_values, mean, numeric(1))
  m <- length(split_values)
  total <- length(value)
  grand <- mean(value)
  ss_within <- sum(vapply(split_values, function(z) sum((z - mean(z))^2),
                          numeric(1)))
  ss_between <- sum(sizes * (means - grand)^2)
  df_within <- total - m
  df_between <- m - 1L
  if (df_within <= 0L) stop("Within-group replication is required.", call. = FALSE)
  ms_within <- ss_within / df_within
  ms_between <- ss_between / df_between
  n_eff <- (total - sum(sizes^2) / total) / (m - 1L)
  between_variance <- max(0, (ms_between - ms_within) / n_eff)
  list(
    n = total, groups = m, sizes = sizes, means = means, grand_mean = grand,
    ss_within = ss_within, ss_between = ss_between,
    df_within = df_within, df_between = df_between,
    ms_within = ms_within, ms_between = ms_between, n_eff = n_eff,
    within_sd = sqrt(ms_within), between_sd = sqrt(between_variance)
  )
}

.unc_ws <- function(u, df) {
  keep <- is.finite(u) & u > 0 & is.finite(df) & df > 0
  total_variance <- sum(u^2)
  if (!any(keep) || total_variance <= 0) return(Inf)
  denominator <- sum(u[keep]^4 / df[keep])
  if (denominator <= 0) Inf else total_variance^2 / denominator
}

#' Estimate top-down uncertainty from internal quality-control data
#'
#' Uses a one-way random-effects ANOVA when runs contain replicates. The object
#' reports within-run, between-run, single-result, and overall-mean estimates.
#'
#' @param data Data frame.
#' @param result Numeric result column.
#' @param run Optional run column. If omitted, the SD of all results is used.
#' @param quantity Component returned in `$component`: `"single"` or `"mean"`.
#' @param name Component name.
#' @param relative Return the selected component on a relative scale.
#' @return An `uncertainty_iqc` object.
#' @references CLSI EP29-A.
#' @examples
#' iqc_data <- data.frame(
#'   run = rep(1:5, each = 2),
#'   result = c(99, 101, 100, 102, 98, 100, 101, 103, 99, 100)
#' )
#' iqc <- uncertainty_iqc(iqc_data, "result", run = "run")
#' iqc$variance_components
#' @export
uncertainty_iqc <- function(data, result, run = NULL,
                            quantity = c("single", "mean"), name = "IQC",
                            relative = FALSE) {
  quantity <- match.arg(quantity)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  result <- .unc_col(data, result, "result", numeric = TRUE)
  run <- .unc_col(data, run, "run", allow_null = TRUE)
  keep <- is.finite(data[[result]])
  if (!is.null(run)) keep <- keep & !is.na(data[[run]])
  values <- as.numeric(data[[result]][keep])
  if (length(values) < 2L) stop("At least two finite results are required.",
                                call. = FALSE)
  estimate <- mean(values)
  if (is.null(run) || all(table(data[[run]][keep]) == 1L)) {
    single <- stats::sd(values)
    mean_u <- single / sqrt(length(values))
    details <- list(method = "long-term SD", within_sd = NA_real_,
                    between_sd = NA_real_, u_single = single,
                    u_mean = mean_u, n = length(values), runs = length(values),
                    df_single = length(values) - 1L,
                    df_mean = length(values) - 1L)
  } else {
    a <- .unc_anova(values, data[[run]][keep])
    single <- sqrt(a$within_sd^2 + a$between_sd^2)
    mean_u <- if (a$between_sd > 0) sqrt(a$ms_between / a$n) else
      stats::sd(values) / sqrt(a$n)
    details <- c(a, list(method = "one-way ANOVA", runs = a$groups,
                         u_single = single,
                         u_mean = mean_u,
                         df_single = .unc_ws(c(a$within_sd, a$between_sd),
                                             c(a$df_within, a$df_between)),
                         df_mean = if (a$between_sd > 0) a$df_between else
                           a$n - 1L))
  }
  selected <- if (quantity == "single") details$u_single else details$u_mean
  selected_df <- if (quantity == "single") details$df_single else details$df_mean
  stored <- if (relative) selected / abs(estimate) else selected
  component <- .unc_component(
    name, estimate, stored, df = selected_df, evaluation_type = "A",
    distribution = "normal", source = paste0("IQC: ", quantity),
    relative = relative,
    parameters = list(absolute_standard_uncertainty = selected)
  )
  structure(list(call = match.call(), result = result, run = run,
                 quantity = quantity, details = details,
                 component = component, excluded = which(!keep)),
            class = "uncertainty_iqc")
}

#' Estimate the between-unit homogeneity uncertainty of a reference material
#'
#' Uses the EP30 one-way ANOVA calculations, including the effective replicate
#' count for an unbalanced design and the inhomogeneity that can be hidden by
#' measurement repeatability. The reported component is the larger of the
#' observed between-unit SD and the hidden-inhomogeneity estimate.
#'
#' @param data Data frame in long format.
#' @param unit Bottle, vial, or material-unit column.
#' @param result Numeric measurement-result column.
#' @param name Component name.
#' @param relative Return the component on a relative scale.
#' @return An `uncertainty_homogeneity` object.
#' @references CLSI EP30-A.
#' @examples
#' homogeneity_data <- data.frame(
#'   unit = rep(paste0("V", 1:4), each = 2),
#'   result = c(99, 100, 101, 100, 98, 99, 102, 101)
#' )
#' homogeneity <- uncertainty_homogeneity(
#'   homogeneity_data, "unit", "result"
#' )
#' homogeneity$component
#' @export
uncertainty_homogeneity <- function(data, unit, result,
                                    name = "between-unit homogeneity",
                                    relative = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  unit <- .unc_col(data, unit, "unit")
  result <- .unc_col(data, result, "result", numeric = TRUE)
  keep <- !is.na(data[[unit]]) & is.finite(data[[result]])
  a <- .unc_anova(data[[result]][keep], data[[unit]][keep])
  mean_nonzero <- is.finite(a$grand_mean) &&
    abs(a$grand_mean) > .Machine$double.eps^0.5
  if (relative && !mean_nonzero) {
    stop("The grand mean must be nonzero for EP30 relative calculations.",
         call. = FALSE)
  }
  average_replicates <- a$n / a$groups
  u_hidden <- a$within_sd / sqrt(average_replicates) *
    (2 / a$df_within)^0.25
  rsd_method <- if (mean_nonzero) a$within_sd / abs(a$grand_mean) else NA_real_
  u_hidden_relative <- if (mean_nonzero) u_hidden / abs(a$grand_mean) else NA_real_
  u_bb <- max(a$between_sd, u_hidden)
  selected_df <- if (a$between_sd >= u_hidden && a$between_sd > 0)
    a$df_between else a$df_within
  stored <- if (relative) u_bb / abs(a$grand_mean) else u_bb
  component <- .unc_component(
    name, a$grand_mean, stored, df = selected_df,
    evaluation_type = "A", distribution = "normal",
    source = "EP30 homogeneity ANOVA", relative = relative,
    parameters = list(absolute_standard_uncertainty = u_bb)
  )
  details <- data.frame(
    n_units = a$groups, n_results = a$n,
    mean = a$grand_mean, n_eff = a$n_eff,
    average_replicates = average_replicates,
    ms_within = a$ms_within, ms_between = a$ms_between,
    df_within = a$df_within, df_between = a$df_between,
    repeatability_sd = a$within_sd,
    between_unit_sd = a$between_sd,
    hidden_inhomogeneity = u_hidden,
    u_bb = u_bb,
    relative_repeatability = rsd_method,
    relative_between_unit = if (mean_nonzero) a$between_sd / abs(a$grand_mean) else NA_real_,
    relative_hidden_inhomogeneity = u_hidden_relative,
    relative_u_bb = if (mean_nonzero) u_bb / abs(a$grand_mean) else NA_real_,
    row.names = NULL
  )
  structure(list(call = match.call(), unit = unit, result = result,
                 details = details, unit_sizes = a$sizes,
                 unit_means = a$means, component = component,
                 excluded = which(!keep)),
            class = "uncertainty_homogeneity")
}

#' Estimate long-term stability uncertainty
#'
#' Fits result on time and calculates the standard uncertainty at the requested
#' shelf life as the slope standard error multiplied by shelf life. Separate
#' components are returned for each optional storage condition.
#'
#' @param data Data frame.
#' @param time Numeric storage-time column.
#' @param result Numeric result column.
#' @param shelf_life Positive time on the same scale as `time`.
#' @param condition Optional storage-condition column.
#' @param name Base component name.
#' @param relative Return components relative to their fitted intercept-level
#'   mean.
#' @return An `uncertainty_stability` object.
#' @references CLSI EP30-A.
#' @examples
#' stability_data <- data.frame(
#'   time = rep(0:4, each = 2),
#'   result = c(100.1, 99.9, 100.0, 99.8, 99.7, 99.9,
#'              99.6, 99.8, 99.5, 99.7)
#' )
#' stability_u <- uncertainty_stability(
#'   stability_data, "time", "result", shelf_life = 6
#' )
#' stability_u$summary
#' @export
uncertainty_stability <- function(data, time, result, shelf_life,
                                  condition = NULL,
                                  name = "long-term stability",
                                  relative = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  time <- .unc_col(data, time, "time", numeric = TRUE)
  result <- .unc_col(data, result, "result", numeric = TRUE)
  condition <- .unc_col(data, condition, "condition", allow_null = TRUE)
  shelf_life <- .unc_scalar(shelf_life, "shelf_life", positive = TRUE)
  keep <- is.finite(data[[time]]) & is.finite(data[[result]])
  if (!is.null(condition)) keep <- keep & !is.na(data[[condition]])
  d <- data[keep, , drop = FALSE]
  if (!nrow(d)) stop("No usable observations remain.", call. = FALSE)
  d$.condition <- if (is.null(condition)) "all" else as.character(d[[condition]])
  groups <- split(d, d$.condition, drop = TRUE)
  rows <- list()
  components <- list()
  for (id in names(groups)) {
    z <- groups[[id]]
    if (nrow(z) < 3L || length(unique(z[[time]])) < 2L) {
      stop("Each condition requires at least three results at two time points.",
           call. = FALSE)
    }
    fit_data <- data.frame(.time = as.numeric(z[[time]]),
                           .result = as.numeric(z[[result]]))
    fit <- stats::lm(.result ~ .time, data = fit_data)
    coefficients <- summary(fit)$coefficients
    slope <- unname(coefficients[".time", "Estimate"])
    slope_se <- unname(coefficients[".time", "Std. Error"])
    standard <- slope_se * shelf_life
    center <- mean(fit_data$.result)
    if (relative && abs(center) <= .Machine$double.eps^0.5) {
      stop("A relative stability component requires a nonzero result mean.",
           call. = FALSE)
    }
    component_name <- if (length(groups) == 1L) name else paste(name, id, sep = ": ")
    components[[id]] <- .unc_component(
      component_name, center,
      if (relative) standard / abs(center) else standard,
      df = stats::df.residual(fit), evaluation_type = "A",
      distribution = "normal", source = "EP30 stability regression",
      relative = relative,
      parameters = list(absolute_standard_uncertainty = standard)
    )
    rows[[id]] <- data.frame(
      condition = id, n = nrow(fit_data), mean = center,
      time_min = min(fit_data$.time), time_max = max(fit_data$.time),
      shelf_life = shelf_life, intercept = stats::coef(fit)[1L],
      slope = slope, slope_se = slope_se,
      residual_sd = summary(fit)$sigma,
      df = stats::df.residual(fit), u_lts = standard,
      relative_u_lts = if (abs(center) > .Machine$double.eps^0.5)
        standard / abs(center) else NA_real_, row.names = NULL
    )
  }
  details <- do.call(rbind, rows)
  rownames(details) <- NULL
  structure(list(call = match.call(), time = time, result = result,
                 condition = condition, shelf_life = shelf_life,
                 details = details, components = unname(components),
                 component = if (length(components) == 1L) components[[1L]] else NULL,
                 excluded = which(!keep)),
            class = "uncertainty_stability")
}

#' Estimate uncertainty from interlaboratory characterization
#'
#' Calculates the unweighted mean of laboratory means and its standard error.
#' Common calibration or group-correlated components should be supplied
#' separately to [uncertainty_combine()] so they are not counted repeatedly.
#'
#' @param data Data frame.
#' @param laboratory Laboratory or measurement-procedure column.
#' @param result Numeric result column.
#' @param name Component name.
#' @param relative Return the component relative to the assigned value.
#' @return An `uncertainty_characterization` object.
#' @references CLSI EP30-A; CLSI EP32-R.
#' @examples
#' characterization_data <- data.frame(
#'   laboratory = rep(paste0("Lab", 1:4), each = 3),
#'   result = c(99, 100, 101, 100, 101, 100,
#'              98, 99, 100, 101, 102, 101)
#' )
#' characterization <- uncertainty_characterization(
#'   characterization_data, "laboratory", "result"
#' )
#' characterization$laboratory_summary
#' @export
uncertainty_characterization <- function(data, laboratory, result,
                                         name = "characterization",
                                         relative = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  laboratory <- .unc_col(data, laboratory, "laboratory")
  result <- .unc_col(data, result, "result", numeric = TRUE)
  keep <- !is.na(data[[laboratory]]) & is.finite(data[[result]])
  d <- data[keep, , drop = FALSE]
  if (!nrow(d)) stop("No usable observations remain.", call. = FALSE)
  groups <- split(as.numeric(d[[result]]), as.character(d[[laboratory]]),
                  drop = TRUE)
  if (length(groups) < 2L) stop("At least two laboratories are required.",
                                call. = FALSE)
  means <- vapply(groups, mean, numeric(1))
  assigned <- mean(means)
  laboratory_sd <- stats::sd(means)
  standard <- laboratory_sd / sqrt(length(means))
  if (relative && abs(assigned) <= .Machine$double.eps^0.5) {
    stop("A relative component requires a nonzero assigned value.", call. = FALSE)
  }
  component <- .unc_component(
    name, assigned, if (relative) standard / abs(assigned) else standard,
    df = length(means) - 1L, evaluation_type = "A",
    distribution = "normal", source = "laboratory means",
    relative = relative,
    parameters = list(absolute_standard_uncertainty = standard)
  )
  laboratory_summary <- data.frame(
    laboratory = names(groups), n = vapply(groups, length, integer(1)),
    mean = means, sd = vapply(groups, function(z) if (length(z) > 1L)
      stats::sd(z) else NA_real_, numeric(1)), row.names = NULL
  )
  structure(list(call = match.call(), laboratory = laboratory, result = result,
                 assigned_value = assigned, laboratory_sd = laboratory_sd,
                 standard_uncertainty = standard,
                 laboratory_summary = laboratory_summary,
                 component = component, excluded = which(!keep)),
            class = "uncertainty_characterization")
}

.unc_print_analysis <- function(x, title, table) {
  cat("\n", title, "\n", sep = "")
  if (nrow(table) == 1L) {
    .print_kv_sections(table, list("Summary" = names(table)))
  } else {
    .print_kv_rows(table, list("Summary" = names(table)))
  }
  invisible(x)
}

#' Print IQC uncertainty results
#'
#' @param x An uncertainty_iqc object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_iqc <- function(x, ...) {
  d <- x$details
  table <- data.frame(method = d$method, n = d$n, runs = d$runs,
                      within_sd = d$within_sd, between_sd = d$between_sd,
                      u_single = d$u_single, u_mean = d$u_mean)
  cat("\nIQC uncertainty\n")
  .print_kv_sections(table, list(
    "Study size" = c("Method" = "method", "Observations" = "n", "Runs" = "runs"),
    "Uncertainty estimates" = c(
      "Within-run SD" = "within_sd", "Between-run SD" = "between_sd",
      "Single-result uncertainty" = "u_single", "Mean uncertainty" = "u_mean")
  ))
  invisible(x)
}

#' Print between-unit homogeneity uncertainty results
#'
#' @param x An uncertainty_homogeneity object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_homogeneity <- function(x, ...) {
  cat("\nBetween-unit homogeneity uncertainty\n")
  if (nrow(x$details) == 1L) {
    .print_kv_sections(x$details, list(
      "Study summary" = c(
        "Units" = "n_units", "Results" = "n_results", "Mean" = "mean",
        "Effective degrees of freedom" = "n_eff",
        "Average replicates" = "average_replicates"),
      "Variance components" = c(
        "Within-unit mean square" = "ms_within",
        "Between-unit mean square" = "ms_between",
        "Within-unit degrees of freedom" = "df_within",
        "Between-unit degrees of freedom" = "df_between",
        "Repeatability SD" = "repeatability_sd",
        "Between-unit SD" = "between_unit_sd"),
      "Selected uncertainty" = c(
        "Hidden inhomogeneity" = "hidden_inhomogeneity",
        "Between-unit uncertainty" = "u_bb",
        "Relative repeatability" = "relative_repeatability",
        "Relative between-unit" = "relative_between_unit",
        "Relative hidden inhomogeneity" = "relative_hidden_inhomogeneity",
        "Relative between-unit uncertainty" = "relative_u_bb")
    ))
  } else {
    .print_kv_rows(x$details, list(
      "Study summary" = c(
        "Units" = "n_units", "Results" = "n_results", "Mean" = "mean",
        "Effective degrees of freedom" = "n_eff"),
      "Replicate summary" = c("Average replicates" = "average_replicates"),
      "Variance components" = c(
        "Within-unit mean square" = "ms_within",
        "Between-unit mean square" = "ms_between",
        "Within-unit degrees of freedom" = "df_within",
        "Between-unit degrees of freedom" = "df_between"),
      "Standard deviations" = c(
        "Repeatability SD" = "repeatability_sd",
        "Between-unit SD" = "between_unit_sd"),
      "Selected uncertainty" = c(
        "Hidden inhomogeneity" = "hidden_inhomogeneity", "Between-unit uncertainty" = "u_bb",
        "Relative repeatability" = "relative_repeatability",
        "Relative between-unit" = "relative_between_unit"),
      "Relative uncertainty" = c(
        "Relative hidden inhomogeneity" = "relative_hidden_inhomogeneity",
        "Relative between-unit uncertainty" = "relative_u_bb")
    ))
  }
  invisible(x)
}

#' Print stability uncertainty results
#'
#' @param x An uncertainty_stability object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_stability <- function(x, ...) {
  cat("\nStability uncertainty\n")
  if (nrow(x$details) == 1L) {
    .print_kv_sections(x$details, list(
      "Study and model" = c(
        "Condition" = "condition", "Observations" = "n",
        "Minimum time" = "time_min", "Maximum time" = "time_max",
        "Shelf life" = "shelf_life", "Intercept" = "intercept",
        "Slope" = "slope", "Mean" = "mean"),
      "Uncertainty" = c(
        "Slope SE" = "slope_se", "Residual SD" = "residual_sd",
        "Degrees of freedom" = "df", "Long-term uncertainty" = "u_lts",
        "Relative long-term uncertainty" = "relative_u_lts")
    ))
  } else {
    .print_kv_rows(x$details, list(
      "Study and model" = c(
        "Observations" = "n", "Minimum time" = "time_min",
        "Maximum time" = "time_max"),
      "Fit estimates" = c(
        "Shelf life" = "shelf_life", "Intercept" = "intercept", "Slope" = "slope"),
      "Mean" = c("Mean" = "mean"),
      "Uncertainty" = c(
        "Slope SE" = "slope_se", "Residual SD" = "residual_sd",
        "Degrees of freedom" = "df"),
      "Long-term uncertainty" = c(
        "Long-term uncertainty" = "u_lts",
        "Relative long-term uncertainty" = "relative_u_lts")
    ), row_title = c("Condition" = "condition"))
  }
  invisible(x)
}

#' Print characterization uncertainty results
#'
#' @param x An uncertainty_characterization object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_characterization <- function(x, ...) {
  table <- data.frame(laboratories = nrow(x$laboratory_summary),
                      assigned_value = x$assigned_value,
                      laboratory_sd = x$laboratory_sd,
                      standard_uncertainty = x$standard_uncertainty)
  .unc_print_analysis(x, "Characterization uncertainty", table)
}

.unc_flatten <- function(x) {
  if (inherits(x, "uncertainty_component")) return(list(x))
  if (is.list(x) && !is.null(x$component) &&
      inherits(x$component, "uncertainty_component")) {
    return(list(x$component))
  }
  if (is.list(x) && !is.null(x$components) &&
      all(vapply(x$components, inherits, logical(1), "uncertainty_component"))) {
    return(x$components)
  }
  if (is.list(x)) return(unlist(lapply(x, .unc_flatten), recursive = FALSE))
  stop("All inputs must be uncertainty components or objects containing them.",
       call. = FALSE)
}

.unc_validate_correlation <- function(correlation, n) {
  if (is.null(correlation)) return(diag(n))
  if (!is.matrix(correlation) || !is.numeric(correlation) ||
      any(dim(correlation) != n) || anyNA(correlation) ||
      any(!is.finite(correlation))) {
    stop("`correlation` must be a finite numeric square matrix matching the ",
         "number of components.", call. = FALSE)
  }
  if (!isTRUE(all.equal(correlation, t(correlation), tolerance = 1e-10)) ||
      any(abs(diag(correlation) - 1) > 1e-10) ||
      any(abs(correlation) > 1 + 1e-10)) {
    stop("`correlation` must be symmetric with a unit diagonal and entries in ",
         "[-1, 1].", call. = FALSE)
  }
  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < -1e-8) {
    stop("`correlation` must be positive semidefinite.", call. = FALSE)
  }
  correlation
}

.unc_combine_list <- function(components, sensitivity = NULL,
                              correlation = NULL, value = NULL,
                              coverage = 0.95, k = NULL) {
  components <- .unc_flatten(components)
  if (!length(components)) stop("At least one component is required.", call. = FALSE)
  coverage <- .unc_probability(coverage)
  n <- length(components)
  if (is.null(sensitivity)) sensitivity <- rep(1, n)
  if (!is.numeric(sensitivity) || length(sensitivity) != n ||
      anyNA(sensitivity) || any(!is.finite(sensitivity))) {
    stop("`sensitivity` must contain one finite value per component.",
         call. = FALSE)
  }
  if (!is.null(value)) value <- .unc_scalar(value, "value")
  relative <- vapply(components, `[[`, logical(1), "relative")
  stored_u <- vapply(components, `[[`, numeric(1), "standard_uncertainty")
  if (any(relative) && any(!relative) && is.null(value)) {
    stop("Combining absolute and relative components requires `value`.",
         call. = FALSE)
  }
  output_relative <- all(relative) && is.null(value)
  u <- stored_u
  if (!output_relative && any(relative)) u[relative] <- u[relative] * abs(value)
  correlation <- .unc_validate_correlation(correlation, n)
  covariance <- outer(u, u) * correlation
  combined_variance <- drop(t(sensitivity) %*% covariance %*% sensitivity)
  if (combined_variance < -1e-10) {
    stop("The propagated variance is negative; inspect the correlation matrix.",
         call. = FALSE)
  }
  combined_variance <- max(0, combined_variance)
  combined <- sqrt(combined_variance)
  marginal <- sensitivity * drop(covariance %*% sensitivity)
  contribution_percent <- if (combined_variance > 0)
    100 * marginal / combined_variance else rep(NA_real_, n)
  independent <- isTRUE(all.equal(correlation, diag(n), tolerance = 1e-12))
  df <- vapply(components, `[[`, numeric(1), "df")
  effective_df <- if (independent && combined > 0) {
    finite <- is.finite(df) & df > 0 & abs(sensitivity * u) > 0
    denominator <- sum((sensitivity[finite] * u[finite])^4 / df[finite])
    if (denominator > 0) combined^4 / denominator else Inf
  } else if (independent) Inf else NA_real_
  if (is.null(k)) {
    coverage_factor <- if (is.finite(effective_df))
      stats::qt((1 + coverage) / 2, df = effective_df) else
        stats::qnorm((1 + coverage) / 2)
    factor_source <- if (is.finite(effective_df)) "Student t" else "normal"
  } else {
    coverage_factor <- .unc_scalar(k, "k", positive = TRUE)
    factor_source <- "supplied"
  }
  expanded <- coverage_factor * combined
  rows <- lapply(seq_len(n), function(i) {
    z <- .unc_component_row(components[[i]])
    z$sensitivity <- sensitivity[i]
    z$uncertainty_on_output_scale <- u[i]
    z$signed_variance_contribution <- marginal[i]
    z$contribution_percent <- contribution_percent[i]
    z
  })
  budget <- do.call(rbind, rows)
  rownames(budget) <- NULL
  relative_combined <- if (output_relative) combined else
    if (!is.null(value) && abs(value) > .Machine$double.eps^0.5)
      combined / abs(value) else NA_real_
  absolute_combined <- if (output_relative) {
    if (is.null(value)) NA_real_ else combined * abs(value)
  } else combined
  structure(list(
    call = NULL, components = components, budget = budget,
    correlation = correlation, covariance = covariance,
    output_scale = if (output_relative) "relative" else "absolute",
    value = value, combined_standard_uncertainty = combined,
    absolute_standard_uncertainty = absolute_combined,
    relative_standard_uncertainty = relative_combined,
    effective_df = effective_df, coverage = coverage,
    coverage_factor = coverage_factor, factor_source = factor_source,
    expanded_uncertainty = expanded
  ), class = "uncertainty_budget")
}

#' Combine standard uncertainty components
#'
#' Applies the GUM first-order covariance propagation formula `c' Sigma c`.
#' Components from any estimator in this module can be supplied directly.
#'
#' @param ... Uncertainty components, analysis objects containing `$component`
#'   or `$components`, or lists of such objects.
#' @param sensitivity Optional sensitivity coefficient per component.
#' @param correlation Optional correlation matrix.
#' @param value Optional output quantity value, needed to mix relative and
#'   absolute components.
#' @param coverage Coverage probability.
#' @param k Optional fixed coverage factor. Otherwise a normal or effective-df
#'   Student t factor is used.
#' @return An `uncertainty_budget` object.
#' @references CLSI EP29-A.
#' @examples
#' repeatability <- uncertainty_type_a(
#'   c(9.9, 10.0, 10.1, 10.0), "repeatability", quantity = "mean"
#' )
#' calibrator <- uncertainty_type_b(
#'   "calibrator", estimate = 10, uncertainty = 0.2, k = 2,
#'   distribution = "normal"
#' )
#' budget <- uncertainty_combine(repeatability, calibrator, value = 10)
#' budget$components
#' @export
uncertainty_combine <- function(..., sensitivity = NULL, correlation = NULL,
                                value = NULL, coverage = 0.95, k = NULL) {
  out <- .unc_combine_list(list(...), sensitivity = sensitivity,
                           correlation = correlation, value = value,
                           coverage = coverage, k = k)
  out$call <- match.call()
  out
}

#' Print a combined uncertainty budget
#'
#' @param x An uncertainty_budget object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_budget <- function(x, ...) {
  cat("\nCombined uncertainty budget\n")
  cat(sprintf("  Scale: %s\n", x$output_scale))
  cat(sprintf("  Combined standard uncertainty: %.6g\n",
              x$combined_standard_uncertainty))
  cat(sprintf("  Expanded uncertainty: %.6g (k = %.5g, %.1f%%)\n\n",
              x$expanded_uncertainty, x$coverage_factor, 100 * x$coverage))
  .print_kv_rows(x$budget, list(
    "Component inputs" = c(
      "Estimate" = "estimate", "Standard uncertainty" = "standard_uncertainty"),
    "Relative component inputs" = c(
      "Relative uncertainty" = "relative_uncertainty", "Degrees of freedom" = "df"),
    "Component metadata" = c(
      "Evaluation type" = "evaluation_type", "Distribution" = "distribution"),
    "Component source and scale" = c("Source" = "source", "Relative" = "relative"),
    "Contribution" = c(
      "Sensitivity" = "sensitivity",
      "Uncertainty on output scale" = "uncertainty_on_output_scale"),
    "Contribution share" = c(
      "Signed variance contribution" = "signed_variance_contribution",
      "Contribution percent" = "contribution_percent")
  ), row_title = c("Component" = "name"))
  invisible(x)
}

#' Plot an uncertainty budget
#'
#' @param x An `uncertainty_budget` object.
#' @param ... Reserved arguments.
#' @return A ggplot object, invisibly.
#' @examples
#' repeatability <- uncertainty_type_a(
#'   c(9.9, 10.0, 10.1, 10.0), "repeatability", quantity = "mean"
#' )
#' calibrator <- uncertainty_type_b(
#'   "calibrator", estimate = 10, uncertainty = 0.2, k = 2,
#'   distribution = "normal"
#' )
#' budget <- uncertainty_combine(repeatability, calibrator, value = 10)
#' plot(budget)
#' @export
plot.uncertainty_budget <- function(x, ...) {
  d <- x$budget
  d$.name <- factor(d$name, levels = d$name[order(d$contribution_percent)])
  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$.name, y = .data$contribution_percent)
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Variance contribution (%)",
                  title = "Uncertainty budget") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)
  invisible(p)
}

.unc_absolute_u <- function(component) {
  if (component$relative) {
    supplied <- component$parameters$absolute_standard_uncertainty
    if (!is.null(supplied) && length(supplied) == 1L && is.finite(supplied)) {
      supplied
    } else component$standard_uncertainty * abs(component$estimate)
  } else component$standard_uncertainty
}

.unc_evaluate_model <- function(model, values, names) {
  result <- do.call(model, as.list(stats::setNames(values, names)))
  if (!is.numeric(result) || length(result) != 1L || !is.finite(result)) {
    stop("`model` must return one finite numeric value for scalar inputs.",
         call. = FALSE)
  }
  as.numeric(result)
}

.unc_sensitivity <- function(model, components) {
  estimates <- vapply(components, `[[`, numeric(1), "estimate")
  names <- vapply(components, `[[`, character(1), "name")
  vapply(seq_along(estimates), function(i) {
    step <- .Machine$double.eps^(1 / 3) *
      max(1, abs(estimates[i]), .unc_absolute_u(components[[i]]))
    upper <- lower <- estimates
    upper[i] <- upper[i] + step
    lower[i] <- lower[i] - step
    (.unc_evaluate_model(model, upper, names) -
       .unc_evaluate_model(model, lower, names)) / (2 * step)
  }, numeric(1))
}

.unc_draw_component <- function(component, n) {
  estimate <- component$estimate
  u <- .unc_absolute_u(component)
  p <- component$parameters
  distribution <- component$distribution
  if (distribution %in% c("normal", "standard")) {
    return(stats::rnorm(n, mean = estimate, sd = u))
  }
  if (distribution == "rectangular") {
    return(stats::runif(n, min = p$lower, max = p$upper))
  }
  if (distribution == "triangular") {
    lower <- p$lower; upper <- p$upper; mode <- estimate
    probability_mode <- (mode - lower) / (upper - lower)
    probability <- stats::runif(n)
    return(ifelse(
      probability < probability_mode,
      lower + sqrt(probability * (upper - lower) * (mode - lower)),
      upper - sqrt((1 - probability) * (upper - lower) * (upper - mode))
    ))
  }
  if (distribution == "u_shaped") {
    return(estimate + p$half_width * cos(stats::runif(n, 0, pi)))
  }
  if (distribution == "student_t") {
    return(estimate + u * stats::rt(n, df = component$df))
  }
  stop("Monte Carlo sampling is unavailable for distribution '", distribution,
       "'.", call. = FALSE)
}

#' Propagate uncertainty through a general measurement model
#'
#' The delta method numerically estimates sensitivity coefficients and applies
#' GUM first-order propagation. Monte Carlo propagation samples the stored
#' component distributions and returns a central quantile interval.
#'
#' @param model Function whose named arguments match component names.
#' @param components Components or analysis objects containing components.
#' @param method `"delta"` or `"monte_carlo"`.
#' @param correlation Optional correlation matrix. Nonidentity correlation in
#'   Monte Carlo mode currently requires normal or standardized components.
#' @param simulations Number of Monte Carlo draws.
#' @param coverage Coverage probability.
#' @param seed Optional random seed.
#' @return An `uncertainty_propagation` object.
#' @references CLSI EP29-A.
#' @examples
#' numerator <- uncertainty_type_b(
#'   "numerator", estimate = 20, uncertainty = 0.4,
#'   distribution = "standard"
#' )
#' denominator <- uncertainty_type_b(
#'   "denominator", estimate = 10, uncertainty = 0.1,
#'   distribution = "standard"
#' )
#' propagated <- uncertainty_propagate(
#'   function(numerator, denominator) numerator / denominator,
#'   list(numerator, denominator), method = "delta"
#' )
#' propagated$estimate
#' @export
uncertainty_propagate <- function(model, components,
                                  method = c("delta", "monte_carlo"),
                                  correlation = NULL, simulations = 100000L,
                                  coverage = 0.95, seed = NULL) {
  method <- match.arg(method)
  if (!is.function(model)) stop("`model` must be a function.", call. = FALSE)
  components <- .unc_flatten(components)
  if (!length(components)) stop("At least one component is required.", call. = FALSE)
  component_names <- vapply(components, `[[`, character(1), "name")
  if (anyDuplicated(component_names)) stop("Component names must be unique.",
                                           call. = FALSE)
  estimates <- vapply(components, `[[`, numeric(1), "estimate")
  output_estimate <- .unc_evaluate_model(model, estimates, component_names)
  coverage <- .unc_probability(coverage)
  correlation_matrix <- .unc_validate_correlation(correlation,
                                                   length(components))
  sensitivity <- .unc_sensitivity(model, components)
  absolute_components <- lapply(components, function(z) {
    z$standard_uncertainty <- .unc_absolute_u(z)
    z$relative <- FALSE
    z$relative_uncertainty <- if (abs(z$estimate) > .Machine$double.eps^0.5)
      z$standard_uncertainty / abs(z$estimate) else NA_real_
    z
  })
  budget <- .unc_combine_list(
    absolute_components, sensitivity = sensitivity,
    correlation = correlation_matrix, value = output_estimate,
    coverage = coverage
  )
  if (method == "delta") {
    return(structure(list(
      call = match.call(), method = method, estimate = output_estimate,
      standard_uncertainty = budget$combined_standard_uncertainty,
      coverage = coverage,
      interval = output_estimate + c(-1, 1) * budget$expanded_uncertainty,
      sensitivity = stats::setNames(sensitivity, component_names),
      budget = budget, simulations = NULL
    ), class = "uncertainty_propagation"))
  }
  if (!is.numeric(simulations) || length(simulations) != 1L ||
      is.na(simulations) || simulations < 1000 ||
      simulations != as.integer(simulations)) {
    stop("`simulations` must be an integer of at least 1000.", call. = FALSE)
  }
  simulations <- as.integer(simulations)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be one nonmissing numeric value.", call. = FALSE)
    }
    set.seed(seed)
  }
  correlated <- !isTRUE(all.equal(correlation_matrix,
                                   diag(length(components)), tolerance = 1e-12))
  if (correlated) {
    distributions <- vapply(components, `[[`, character(1), "distribution")
    if (any(!distributions %in% c("normal", "standard"))) {
      stop("Correlated Monte Carlo propagation currently requires normal or ",
           "standard components.", call. = FALSE)
    }
    root <- chol(correlation_matrix + diag(1e-12, length(components)))
    z <- matrix(stats::rnorm(simulations * length(components)),
                nrow = simulations) %*% root
    draws <- sweep(z, 2L, vapply(components, .unc_absolute_u, numeric(1)), `*`)
    draws <- sweep(draws, 2L, estimates, `+`)
  } else {
    draws <- do.call(cbind, lapply(components, .unc_draw_component,
                                  n = simulations))
  }
  colnames(draws) <- component_names
  simulated <- tryCatch(
    do.call(model, as.data.frame(draws, check.names = FALSE)),
    error = function(e) NULL
  )
  if (is.null(simulated) || !is.numeric(simulated) ||
      length(simulated) != simulations) {
    simulated <- apply(draws, 1L, function(z)
      .unc_evaluate_model(model, z, component_names))
  }
  if (any(!is.finite(simulated))) {
    stop("Monte Carlo propagation produced nonfinite output values.",
         call. = FALSE)
  }
  alpha <- (1 - coverage) / 2
  interval <- unname(stats::quantile(simulated, c(alpha, 1 - alpha),
                                     names = FALSE, type = 8))
  structure(list(
    call = match.call(), method = method, estimate = output_estimate,
    simulated_mean = mean(simulated), simulated_median = stats::median(simulated),
    standard_uncertainty = stats::sd(simulated), coverage = coverage,
    interval = interval,
    lower_uncertainty = output_estimate - interval[1L],
    upper_uncertainty = interval[2L] - output_estimate,
    sensitivity = stats::setNames(sensitivity, component_names),
    budget = budget, simulations = simulated, draws = draws
  ), class = "uncertainty_propagation")
}

#' Print propagated uncertainty results
#'
#' @param x An uncertainty_propagation object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_propagation <- function(x, ...) {
  cat("\nUncertainty propagation\n")
  cat(sprintf("  Method: %s\n", x$method))
  cat(sprintf("  Model estimate: %.6g\n", x$estimate))
  cat(sprintf("  Standard uncertainty: %.6g\n", x$standard_uncertainty))
  cat(sprintf("  %.1f%% interval: [%.6g, %.6g]\n",
              100 * x$coverage, x$interval[1L], x$interval[2L]))
  invisible(x)
}

#' Plot propagated uncertainty
#'
#' @param x An `uncertainty_propagation` object.
#' @param type `"distribution"` or `"contribution"`.
#' @param ... Reserved arguments.
#' @return A ggplot object, invisibly.
#' @examples
#' numerator <- uncertainty_type_b(
#'   "numerator", estimate = 10, uncertainty = 0.2, k = 2,
#'   distribution = "normal"
#' )
#' denominator <- uncertainty_type_b(
#'   "denominator", estimate = 2, uncertainty = 0.04, k = 2,
#'   distribution = "normal"
#' )
#' propagated <- uncertainty_propagate(
#'   function(numerator, denominator) numerator / denominator,
#'   list(numerator, denominator), method = "delta"
#' )
#' plot(propagated, type = "contribution")
#' @export
plot.uncertainty_propagation <- function(
    x, type = c("distribution", "contribution"), ...) {
  type <- match.arg(type)
  if (type == "contribution") return(plot(x$budget))
  if (is.null(x$simulations)) {
    stop("A distribution plot requires Monte Carlo propagation.", call. = FALSE)
  }
  d <- data.frame(.simulated = x$simulations)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$.simulated)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::geom_vline(xintercept = x$interval, linetype = 2,
                        color = "firebrick") +
    ggplot2::labs(x = "Model output", y = "Count",
                  title = "Propagated output distribution") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)
  invisible(p)
}

#' Accumulate uncertainty along a metrological traceability chain
#'
#' @param stages Named list. Each element contains the new uncertainty
#'   components introduced at that traceability stage.
#' @param value Optional assigned value for absolute/relative conversion.
#' @param coverage Coverage probability.
#' @param k Optional fixed coverage factor.
#' @return An `uncertainty_traceability` object.
#' @references CLSI EP29-A; CLSI EP32-R.
#' @examples
#' calibration <- uncertainty_type_b(
#'   "calibration", estimate = 100, uncertainty = 1, k = 2,
#'   distribution = "normal"
#' )
#' transport <- uncertainty_type_b(
#'   "transport", estimate = 100, lower = 99.5, upper = 100.5,
#'   distribution = "rectangular"
#' )
#' chain <- uncertainty_traceability(
#'   list(calibration = list(calibration), transport = list(transport)),
#'   value = 100
#' )
#' chain$summary
#' @export
uncertainty_traceability <- function(stages, value = NULL,
                                     coverage = 0.95, k = NULL) {
  if (!is.list(stages) || !length(stages)) {
    stop("`stages` must be a nonempty named list.", call. = FALSE)
  }
  if (is.null(names(stages)) || any(!nzchar(names(stages))) ||
      anyDuplicated(names(stages))) {
    stop("Every traceability stage must have a unique nonempty name.",
         call. = FALSE)
  }
  accumulated <- list()
  budgets <- vector("list", length(stages))
  rows <- vector("list", length(stages))
  for (i in seq_along(stages)) {
    introduced <- .unc_flatten(stages[[i]])
    accumulated <- c(accumulated, introduced)
    budget <- .unc_combine_list(accumulated, value = value,
                                coverage = coverage, k = k)
    budgets[[i]] <- budget
    rows[[i]] <- data.frame(
      stage = names(stages)[i], introduced_components = length(introduced),
      cumulative_components = length(accumulated),
      cumulative_standard_uncertainty = budget$combined_standard_uncertainty,
      relative_standard_uncertainty = budget$relative_standard_uncertainty,
      coverage_factor = budget$coverage_factor,
      expanded_uncertainty = budget$expanded_uncertainty,
      row.names = NULL
    )
  }
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  structure(list(call = match.call(), stages = stages, value = value,
                 coverage = coverage, budgets = budgets, summary = summary,
                 final = budgets[[length(budgets)]]),
            class = "uncertainty_traceability")
}

#' Print traceability-chain uncertainty results
#'
#' @param x An uncertainty_traceability object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.uncertainty_traceability <- function(x, ...) {
  cat("\nTraceability-chain uncertainty\n")
  .print_kv_rows(x$summary, list(
    "Traceability stages" = c(
      "Introduced components" = "introduced_components",
      "Cumulative components" = "cumulative_components"),
    "Cumulative standard uncertainty" = c(
      "Cumulative standard uncertainty" = "cumulative_standard_uncertainty"),
    "Relative standard uncertainty" = c(
      "Relative standard uncertainty" = "relative_standard_uncertainty"),
    "Coverage and expanded uncertainty" = c(
      "Coverage factor" = "coverage_factor",
      "Expanded uncertainty" = "expanded_uncertainty")
  ), row_title = c("Stage" = "stage"))
  invisible(x)
}

#' Plot cumulative uncertainty through a traceability chain
#'
#' @param x An `uncertainty_traceability` object.
#' @param ... Reserved arguments.
#' @return A ggplot object, invisibly.
#' @examples
#' calibration <- uncertainty_type_b(
#'   "calibration", estimate = 100, uncertainty = 1, k = 2,
#'   distribution = "normal"
#' )
#' transport <- uncertainty_type_b(
#'   "transport", estimate = 100, lower = 99.5, upper = 100.5,
#'   distribution = "rectangular"
#' )
#' chain <- uncertainty_traceability(
#'   list(calibration = list(calibration), transport = list(transport)),
#'   value = 100
#' )
#' plot(chain)
#' @export
plot.uncertainty_traceability <- function(x, ...) {
  d <- x$summary
  d$.stage <- factor(d$stage, levels = d$stage)
  d$.cumulative <- d$cumulative_standard_uncertainty
  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$.stage, y = .data$.cumulative, group = 1)
  ) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::geom_point(color = "steelblue", size = 2.7) +
    ggplot2::labs(x = "Traceability stage", y = "Cumulative standard uncertainty",
                  title = "Traceability-chain uncertainty") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  print(p)
  invisible(p)
}
