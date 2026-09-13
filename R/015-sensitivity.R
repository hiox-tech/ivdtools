# Analytical sensitivity tools based on CLSI EP17-A2.

.sensitivity_column <- function(data, column, argument, optional = FALSE) {
  if (optional && is.null(column)) return(NULL)
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column)) {
    stop("`", argument, "` must be a single column name.", call. = FALSE)
  }
  if (!column %in% names(data)) {
    stop("Column '", column, "' supplied to `", argument,
         "` was not found in `data`.", call. = FALSE)
  }
  column
}

.sensitivity_probability <- function(x, argument) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop("`", argument, "` must be strictly between 0 and 1.", call. = FALSE)
  }
  as.numeric(x)
}

.sensitivity_number <- function(x, argument, lower = -Inf, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      if (strict) x <= lower else x < lower) {
    relation <- if (strict) "greater than" else "greater than or equal to"
    stop("`", argument, "` must be a finite number ", relation, " ", lower,
         ".", call. = FALSE)
  }
  as.numeric(x)
}

.sensitivity_prepare <- function(data, result, sample, lot = NULL,
                                 extra = character(), exclude = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  result <- .sensitivity_column(data, result, "result")
  sample <- .sensitivity_column(data, sample, "sample")
  lot <- .sensitivity_column(data, lot, "lot", optional = TRUE)
  for (nm in extra) .sensitivity_column(data, nm, nm)
  if (!is.numeric(data[[result]])) {
    stop("`result` must identify a numeric column.", call. = FALSE)
  }
  keep <- rep(TRUE, nrow(data))
  if (!is.null(exclude)) {
    if (is.logical(exclude) && length(exclude) == nrow(data)) {
      keep <- !exclude
    } else if (is.numeric(exclude) && all(is.finite(exclude)) &&
               all(exclude == as.integer(exclude))) {
      idx <- as.integer(exclude)
      if (any(idx < 1L | idx > nrow(data))) {
        stop("Numeric `exclude` indices are outside `data`.", call. = FALSE)
      }
      keep[idx] <- FALSE
    } else {
      stop("`exclude` must be row indices or a logical vector with one value per row.",
           call. = FALSE)
    }
  }
  needed <- unique(c(result, sample, lot, extra))
  complete <- stats::complete.cases(data[, needed, drop = FALSE])
  used <- data[keep & complete, , drop = FALSE]
  if (!nrow(used)) stop("No complete, nonexcluded observations remain.", call. = FALSE)
  list(data = used, excluded = which(!keep), incomplete = which(keep & !complete),
       result = result, sample = sample, lot = lot)
}

.sensitivity_lot_groups <- function(data, lot, rule = c("ep17", "separate", "combined")) {
  rule <- match.arg(rule)
  if (is.null(lot)) {
    return(list(groups = list(combined = data), mode = "combined", lots = 1L))
  }
  lot_values <- unique(as.character(data[[lot]]))
  n_lots <- length(lot_values)
  separate <- rule == "separate" || (rule == "ep17" && n_lots %in% 2:3)
  if (!separate) {
    return(list(groups = list(combined = data), mode = "combined", lots = n_lots))
  }
  groups <- split(data, as.character(data[[lot]]), drop = TRUE)
  list(groups = groups, mode = "separate", lots = n_lots)
}

.sensitivity_extract_limit <- function(x, type, argument) {
  if (inherits(x, paste0("sensitivity_", type))) x <- x$estimate
  .sensitivity_number(x, argument)
}

.sensitivity_cp <- function(risk, df) {
  if (!is.finite(df) || df <= 0) return(NA_real_)
  stats::qnorm(1 - risk) / sqrt(1 - 1 / (4 * df))
}

.sensitivity_new <- function(type, method, ...) {
  structure(c(list(call = match.call(), type = type, method = method), list(...)),
            class = c(paste0("sensitivity_", type), "sensitivity"))
}

.sensitivity_rank_quantile <- function(x, probability) {
  x <- sort(x)
  rank <- 0.5 + length(x) * probability
  lo <- max(1L, min(length(x), floor(rank)))
  hi <- max(1L, min(length(x), ceiling(rank)))
  fraction <- rank - floor(rank)
  value <- x[lo] + fraction * (x[hi] - x[lo])
  list(value = unname(value), rank = rank, lower_rank = lo,
       upper_rank = hi, fraction = fraction)
}

.sensitivity_root <- function(fun, observed, lower = 0) {
  observed <- observed[is.finite(observed)]
  if (!length(observed)) return(NA_real_)
  lo <- max(lower, min(observed) / 10, .Machine$double.eps^0.5)
  hi <- max(observed) * 2
  if (!is.finite(hi) || hi <= lo) hi <- lo * 10
  for (multiplier in c(1, 2, 5, 10, 100)) {
    grid <- exp(seq(log(lo), log(hi * multiplier), length.out = 1000L))
    values <- suppressWarnings(fun(grid))
    ok <- is.finite(values)
    if (sum(ok) < 2L) next
    grid <- grid[ok]
    values <- values[ok]
    exact <- which(values == 0)
    if (length(exact)) return(min(grid[exact]))
    crossing <- which(values[-length(values)] * values[-1L] < 0)
    if (!length(crossing)) next
    roots <- vapply(crossing, function(i) {
      tryCatch(stats::uniroot(fun, c(grid[i], grid[i + 1L]), maxiter = 1000L)$root,
               error = function(e) NA_real_)
    }, numeric(1))
    roots <- roots[is.finite(roots) & roots >= lower]
    if (length(roots)) return(min(roots))
  }
  NA_real_
}

#' Establish a limit of blank
#'
#' Establishes a LoB from replicate-level blank results using the CLSI EP17-A2
#' nonparametric, parametric, or zero-LoB procedure. With two or three reagent
#' lots, lots are analyzed separately and the maximum estimate is reported;
#' four or more lots are combined by default.
#'
#' @param data A data frame containing replicate-level results.
#' @param result,sample Column names for the numeric result and blank sample ID.
#' @param lot Optional reagent-lot column.
#' @param method One of `"nonparametric"`, `"parametric"`, or `"zero"`.
#' @param alpha Type I error risk.
#' @param lot_rule EP17 lot handling, forced separate analysis, or forced pooling.
#' @param exclude Optional row indices or logical exclusion indicator. Exclusions
#'   are recorded in the returned object.
#' @return An object of class `sensitivity_lob`.
#' @references CLSI. Evaluation of Detection Capability for Clinical
#'   Laboratory Measurement Procedures. EP17-A2, 2nd ed. 2012.
#' @examples
#' blank <- data.frame(
#'   sample = rep(c("blank_1", "blank_2"), each = 5),
#'   result = c(0.02, 0.04, 0.03, 0.05, 0.01,
#'              0.03, 0.06, 0.04, 0.02, 0.05)
#' )
#' lob <- sensitivity_lob(blank, result = "result", sample = "sample")
#' lob$estimate
#' @export
sensitivity_lob <- function(data, result, sample, lot = NULL,
                            method = c("nonparametric", "parametric", "zero"),
                            alpha = 0.05,
                            lot_rule = c("ep17", "separate", "combined"),
                            exclude = NULL) {
  method <- match.arg(method)
  lot_rule <- match.arg(lot_rule)
  alpha <- .sensitivity_probability(alpha, "alpha")
  p <- .sensitivity_prepare(data, result, sample, lot, exclude = exclude)
  grouping <- .sensitivity_lot_groups(p$data, p$lot, lot_rule)

  details <- lapply(names(grouping$groups), function(label) {
    d <- grouping$groups[[label]]
    y <- d[[p$result]]
    B <- length(y)
    K <- length(unique(d[[p$sample]]))
    if (method == "nonparametric") {
      q <- .sensitivity_rank_quantile(y, 1 - alpha)
      data.frame(lot = label, lob = q$value, B = B, K = K, df = B - K,
                 mean = mean(y), sd = stats::sd(y), cp = NA_real_,
                 rank = q$rank, lower_rank = q$lower_rank,
                 upper_rank = q$upper_rank, false_positive_rate = NA_real_)
    } else if (method == "parametric") {
      df <- B - K
      if (df < 1L) stop("Blank degrees of freedom (`B - K`) must be positive.",
                        call. = FALSE)
      cp <- .sensitivity_cp(alpha, df)
      data.frame(lot = label, lob = mean(y) + cp * stats::sd(y), B = B, K = K,
                 df = df, mean = mean(y), sd = stats::sd(y), cp = cp,
                 rank = NA_real_, lower_rank = NA_integer_,
                 upper_rank = NA_integer_, false_positive_rate = NA_real_)
    } else {
      data.frame(lot = label, lob = 0, B = B, K = K, df = B - K,
                 mean = mean(y), sd = stats::sd(y), cp = NA_real_,
                 rank = NA_real_, lower_rank = NA_integer_,
                 upper_rank = NA_integer_, false_positive_rate = mean(y > 0))
    }
  })
  by_lot <- do.call(rbind, details)
  rownames(by_lot) <- NULL
  estimate <- max(by_lot$lob, na.rm = TRUE)
  .sensitivity_new("lob", method, estimate = estimate, alpha = alpha,
                   lot_rule = lot_rule, lot_mode = grouping$mode,
                   by_lot = by_lot, data = p$data, result = p$result,
                   sample = p$sample, lot = p$lot, excluded = p$excluded,
                   incomplete = p$incomplete)
}

.sensitivity_summary <- function(data, result, sample) {
  groups <- split(data[[result]], as.character(data[[sample]]), drop = TRUE)
  out <- data.frame(
    sample = names(groups),
    mean = vapply(groups, mean, numeric(1)),
    sd = vapply(groups, stats::sd, numeric(1)),
    n = vapply(groups, length, integer(1)),
    stringsAsFactors = FALSE
  )
  out$df <- out$n - 1L
  out
}

.sensitivity_fit_profile <- function(profile_data, model, model_no = 1:9) {
  x <- profile_data$mean
  y <- profile_data$sd
  if (nrow(profile_data) < if (model == "quadratic") 4L else 3L) {
    stop("Too few sample levels for the requested precision profile.", call. = FALSE)
  }
  if (model == "linear") {
    fit <- stats::lm(sd ~ mean, data = profile_data)
    predict_sd <- function(z) pmax(0, as.numeric(stats::predict(fit,
      newdata = data.frame(mean = z))))
    metrics <- data.frame(model = "linear", AIC = stats::AIC(fit),
                          r_squared = base::summary(fit)$r.squared)
  } else if (model == "quadratic") {
    fit <- stats::lm(sd ~ mean + I(mean^2), data = profile_data)
    predict_sd <- function(z) pmax(0, as.numeric(stats::predict(fit,
      newdata = data.frame(mean = z))))
    metrics <- data.frame(model = "quadratic", AIC = stats::AIC(fit),
                          r_squared = base::summary(fit)$r.squared)
  } else {
    vd <- data.frame(Mean = x, VC = y^2, DF = pmax(1, profile_data$df))
    candidates <- if (model == "sadler") 1:9 else model_no
    fit <- suppressMessages(VFP::fit_vfp(Data = vd, model.no = candidates))
    selected <- if (length(candidates) == 1L) candidates else
      as.integer(sub("Model_", "", names(fit$AIC)[which.min(fit$AIC)]))
    predict_sd <- function(z) {
      pred <- suppressWarnings(stats::predict(fit, model.no = selected,
                                               newdata = as.numeric(z)))
      sqrt(pmax(0, as.numeric(pred$Fitted)))
    }
    best_name <- paste0("Model_", selected)
    metrics <- data.frame(model = paste0(model, ":", best_name),
                          AIC = unname(fit$AIC[best_name]), r_squared = NA_real_)
    attr(fit, "selected_model") <- selected
  }
  list(fit = fit, predict_sd = predict_sd, metrics = metrics)
}

.sensitivity_precision_summary <- function(data, result, sample) {
  if (!inherits(data, "precision")) return(NULL)
  if (is.null(data$results)) {
    stop("A `precision` object must be processed by `variance()` first.",
         call. = FALSE)
  }
  rows <- lapply(data$results, function(z) {
    total_sd <- if (!is.null(z$vc)) sqrt(sum(z$vc, na.rm = TRUE)) else NA_real_
    data.frame(sample = z$sample, mean = z$mean, sd = total_sd,
               n = if (!is.null(z$n)) z$n else NA_integer_,
               df = if (!is.null(z$n)) z$n - 1L else 1L)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Establish a limit of detection
#'
#' Establishes LoD using the classical, nonparametric trial, precision-profile,
#' or probit approach described in CLSI EP17-A2.
#'
#' @param data A replicate-level data frame, or a completed `precision` object
#'   for `method = "profile"`.
#' @param result,sample Column names for result and sample ID. `result` may be
#'   `NULL` only when a `precision` object is supplied.
#' @param lob Numeric LoB or a `sensitivity_lob` object.
#' @param lot Optional reagent-lot column.
#' @param concentration Concentration column for probit analysis.
#' @param detected Optional detection indicator or hit-count column for probit.
#' @param method LoD method.
#' @param beta Type II error risk.
#' @param profile_model Precision-profile model.
#' @param concentration_scale Probit concentration scale.
#' @param total Optional total-count column for summarized probit data.
#' @param strata Optional additional probit stratification columns.
#' @param model_no VFP model numbers when `profile_model = "vfp"`.
#' @param lot_rule EP17, separate, or combined lot handling.
#' @param exclude Optional excluded rows.
#' @return An object of class `sensitivity_lod`.
#' @references CLSI. Evaluation of Detection Capability for Clinical
#'   Laboratory Measurement Procedures. EP17-A2, 2nd ed. 2012.
#' @examples
#' low <- data.frame(
#'   sample = rep(c("low_1", "low_2"), each = 5),
#'   result = c(0.22, 0.26, 0.24, 0.25, 0.23,
#'              0.31, 0.29, 0.33, 0.30, 0.32)
#' )
#' lod <- sensitivity_lod(low, "result", "sample", lob = 0.10,
#'                        method = "classical")
#' lod$estimate
#' @export
sensitivity_lod <- function(data, result = NULL, sample = NULL, lob,
                            lot = NULL, concentration = NULL, detected = NULL,
                            method = c("classical", "nonparametric", "profile", "probit"),
                            beta = 0.05,
                            profile_model = c("linear", "quadratic", "sadler", "vfp"),
                            concentration_scale = c("log10", "identity"),
                            total = NULL, strata = NULL, model_no = 1:9,
                            lot_rule = c("ep17", "separate", "combined"),
                            exclude = NULL) {
  method <- match.arg(method)
  profile_model <- match.arg(profile_model)
  concentration_scale <- match.arg(concentration_scale)
  lot_rule <- match.arg(lot_rule)
  beta <- .sensitivity_probability(beta, "beta")
  lob_value <- .sensitivity_extract_limit(lob, "lob", "lob")

  precision_summary <- .sensitivity_precision_summary(data, result, sample)
  if (!is.null(precision_summary)) {
    if (method != "profile") stop("A `precision` object is valid only for profile LoD.",
                                  call. = FALSE)
    groups <- list(combined = precision_summary)
    p <- list(result = "result", sample = "sample", lot = NULL,
              data = NULL, excluded = integer(), incomplete = integer())
    grouping <- list(groups = groups, mode = "combined", lots = 1L)
  } else {
    extra <- character()
    if (method == "probit") {
      if (is.null(concentration)) stop("`concentration` is required for probit LoD.",
                                       call. = FALSE)
      extra <- c(concentration, detected, total, strata)
    }
    p <- .sensitivity_prepare(data, result, sample, lot,
                              extra = extra, exclude = exclude)
    grouping <- .sensitivity_lot_groups(p$data, p$lot, lot_rule)
  }

  if (method %in% c("classical", "nonparametric")) {
    per_sample <- list()
    rows <- lapply(names(grouping$groups), function(label) {
      d <- grouping$groups[[label]]
      sm <- .sensitivity_summary(d, p$result, p$sample)
      per_sample[[label]] <<- sm
      if (method == "classical") {
        df <- sum(sm$df)
        if (df < 1L || any(!is.finite(sm$sd))) {
          stop("Each low-level sample needs at least two results.", call. = FALSE)
        }
        pooled_sd <- sqrt(sum(sm$df * sm$sd^2) / df)
        cp <- .sensitivity_cp(beta, df)
        cochran <- max(sm$sd^2) / sum(sm$sd^2)
        data.frame(lot = label, lod = lob_value + cp * pooled_sd,
                   lob = lob_value, L = sum(sm$n), J = nrow(sm), df = df,
                   pooled_sd = pooled_sd, cp = cp, cochran = cochran,
                   below_lob = NA_integer_, proportion_below = NA_real_,
                   median = NA_real_)
      } else {
        y <- d[[p$result]]
        below <- sum(y < lob_value)
        proportion <- below / length(y)
        candidate <- stats::median(y)
        data.frame(lot = label,
                   lod = if (proportion <= beta) candidate else NA_real_,
                   lob = lob_value, L = length(y), J = nrow(sm),
                   df = length(y) - nrow(sm), pooled_sd = NA_real_, cp = NA_real_,
                   cochran = NA_real_, below_lob = below,
                   proportion_below = proportion, median = candidate)
      }
    })
    by_lot <- do.call(rbind, rows)
    estimate <- if (any(is.finite(by_lot$lod))) max(by_lot$lod, na.rm = TRUE) else NA_real_
    return(.sensitivity_new("lod", method, estimate = estimate, lob = lob_value,
      beta = beta, lot_rule = lot_rule, lot_mode = grouping$mode,
      by_lot = by_lot, per_sample = per_sample, data = p$data,
      result = p$result, sample = p$sample, lot = p$lot,
      excluded = p$excluded, incomplete = p$incomplete))
  }

  if (method == "profile") {
    fits <- list()
    curves <- list()
    rows <- lapply(names(grouping$groups), function(label) {
      d <- grouping$groups[[label]]
      sm <- if (all(c("mean", "sd", "n", "df") %in% names(d))) d else
        .sensitivity_summary(d, p$result, p$sample)
      sm <- sm[is.finite(sm$mean) & is.finite(sm$sd) & sm$sd >= 0, , drop = FALSE]
      fitted <- .sensitivity_fit_profile(sm, profile_model, model_no)
      df <- sum(sm$df, na.rm = TRUE)
      cp <- .sensitivity_cp(beta, df)
      root <- .sensitivity_root(function(x) x - lob_value - cp * fitted$predict_sd(x),
                                sm$mean, lower = lob_value)
      grid <- seq(max(0, min(sm$mean) / 2), max(sm$mean) * 1.5, length.out = 400L)
      curve <- data.frame(lot = label, concentration = grid,
                          predicted_sd = fitted$predict_sd(grid))
      curve$trial_lod <- lob_value + cp * curve$predicted_sd
      fits[[label]] <<- fitted
      curves[[label]] <<- curve
      data.frame(lot = label, lod = root, lob = lob_value, df = df, cp = cp,
                 model = fitted$metrics$model, AIC = fitted$metrics$AIC,
                 r_squared = fitted$metrics$r_squared,
                 min_observed = min(sm$mean), max_observed = max(sm$mean),
                 extrapolated = is.finite(root) &&
                   (root < min(sm$mean) || root > max(sm$mean)))
    })
    by_lot <- do.call(rbind, rows)
    estimate <- if (any(is.finite(by_lot$lod))) max(by_lot$lod, na.rm = TRUE) else NA_real_
    return(.sensitivity_new("lod", method, estimate = estimate, lob = lob_value,
      beta = beta, profile_model = profile_model, lot_rule = lot_rule,
      lot_mode = grouping$mode, by_lot = by_lot, fits = fits,
      curves = curves, data = if (is.null(precision_summary)) p$data else precision_summary,
      result = p$result, sample = p$sample, lot = p$lot,
      excluded = p$excluded, incomplete = p$incomplete))
  }

  concentration <- .sensitivity_column(p$data, concentration, "concentration")
  detected <- .sensitivity_column(p$data, detected, "detected")
  if (!is.numeric(p$data[[concentration]])) {
    stop("`concentration` must identify a numeric column.", call. = FALSE)
  }
  if (concentration_scale == "log10" && any(p$data[[concentration]] <= 0)) {
    stop("Probit concentrations must be positive on the log10 scale.", call. = FALSE)
  }
  split_data <- list()
  for (base_label in names(grouping$groups)) {
    base_data <- grouping$groups[[base_label]]
    if (length(strata)) {
      subgroup <- interaction(base_data[strata], drop = TRUE)
      pieces <- split(base_data, subgroup, drop = TRUE)
      names(pieces) <- paste(base_label, names(pieces), sep = ":")
      split_data <- c(split_data, pieces)
    } else {
      split_data[[base_label]] <- base_data
    }
  }
  fits <- list()
  rates <- list()
  rows <- lapply(names(split_data), function(label) {
    d <- split_data[[label]]
    if (is.null(total)) {
      z <- split(d[[detected]], d[[concentration]], drop = TRUE)
      rate <- data.frame(concentration = as.numeric(names(z)),
                         hits = vapply(z, function(v) sum(as.numeric(v) != 0), numeric(1)),
                         total = vapply(z, length, integer(1)))
    } else {
      total_col <- .sensitivity_column(d, total, "total")
      rate <- data.frame(concentration = d[[concentration]],
                         hits = d[[detected]], total = d[[total_col]])
    }
    if (any(rate$hits < 0 | rate$total <= 0 | rate$hits > rate$total)) {
      stop("Probit hit and total counts are invalid.", call. = FALSE)
    }
    rate$x <- if (concentration_scale == "log10") log10(rate$concentration) else
      rate$concentration
    fit <- stats::glm(cbind(hits, total - hits) ~ x, data = rate,
                      family = stats::binomial(link = "probit"))
    cf <- stats::coef(fit)
    target_x <- (stats::qnorm(1 - beta) - cf[1L]) / cf[2L]
    estimate <- if (concentration_scale == "log10") 10^target_x else target_x
    grid <- seq(min(rate$x), max(rate$x), length.out = 200L)
    fitted_rate <- stats::pnorm(cf[1L] + cf[2L] * grid)
    rate$observed_proportion <- rate$hits / rate$total
    rate$group <- label
    rates[[label]] <<- list(observed = rate,
      curve = data.frame(group = label,
        concentration = if (concentration_scale == "log10") 10^grid else grid,
        probability = fitted_rate))
    fits[[label]] <<- fit
    pearson <- sum(stats::residuals(fit, type = "pearson")^2)
    data.frame(group = label, lod = estimate, target_probability = 1 - beta,
               intercept = cf[1L], slope = cf[2L], deviance = fit$deviance,
               df_residual = fit$df.residual,
               deviance_p = stats::pchisq(fit$deviance, fit$df.residual,
                                          lower.tail = FALSE),
               pearson = pearson,
               pearson_p = stats::pchisq(pearson, fit$df.residual,
                                         lower.tail = FALSE))
  })
  by_group <- do.call(rbind, rows)
  estimate <- max(by_group$lod, na.rm = TRUE)
  .sensitivity_new("lod", method, estimate = estimate, lob = lob_value,
    beta = beta, concentration_scale = concentration_scale, by_group = by_group,
    fits = fits, rates = rates, data = p$data, result = p$result,
    sample = p$sample, lot = p$lot, concentration = concentration,
    detected = detected, excluded = p$excluded, incomplete = p$incomplete)
}

.sensitivity_assigned <- function(data, assigned) {
  if (is.character(assigned) && length(assigned) == 1L) {
    .sensitivity_column(data, assigned, "assigned")
    value <- data[[assigned]]
    name <- assigned
  } else if (is.numeric(assigned) && length(assigned) %in% c(1L, nrow(data))) {
    value <- rep_len(assigned, nrow(data))
    name <- ".assigned"
  } else {
    stop("`assigned` must be a numeric value/vector or a single column name.",
         call. = FALSE)
  }
  if (anyNA(value) || any(!is.finite(value))) {
    stop("Assigned values must be finite and nonmissing.", call. = FALSE)
  }
  list(value = as.numeric(value), name = name)
}

#' Establish a limit of quantitation
#'
#' Establishes LoQ from Westgard total error, RMS error, a precision goal, or
#' a user-supplied accuracy function. The original candidate and the LoD-
#' constrained reported LoQ are retained separately.
#'
#' @param data Replicate-level data frame.
#' @param result,sample Column names for result and sample ID.
#' @param assigned Assigned-value column or numeric value/vector.
#' @param lot Optional reagent-lot column.
#' @param method Accuracy definition.
#' @param goal Accuracy goal; proportions are used for relative goals.
#' @param goal_scale Relative or absolute scale.
#' @param lod Optional numeric LoD or `sensitivity_lod` object.
#' @param k SD multiplier for the Westgard definition.
#' @param accuracy_function Function used by `method = "custom"`.
#' @param profile Fit an error profile and calculate its goal crossing?
#' @param profile_model Linear or quadratic error profile.
#' @param lot_rule EP17, separate, or combined lot handling.
#' @param exclude Optional excluded rows.
#' @param ... Additional arguments passed to `accuracy_function`.
#' @return An object of class `sensitivity_loq`.
#' @references CLSI. Evaluation of Detection Capability for Clinical
#'   Laboratory Measurement Procedures. EP17-A2, 2nd ed. 2012.
#' @examples
#' loq_data <- data.frame(
#'   sample = rep(paste0("L", 1:3), each = 4),
#'   assigned = rep(c(0.5, 1, 2), each = 4),
#'   result = c(0.42, 0.48, 0.51, 0.46,
#'              0.91, 1.02, 0.96, 1.05,
#'              1.91, 2.04, 1.98, 2.08)
#' )
#' loq <- sensitivity_loq(
#'   loq_data, "result", "sample", assigned = "assigned",
#'   method = "precision", goal = 0.15, goal_scale = "relative"
#' )
#' loq$summary
#' @export
sensitivity_loq <- function(data, result, sample, assigned, lot = NULL,
                            method = c("westgard", "rms", "precision", "custom"),
                            goal, goal_scale = c("relative", "absolute"),
                            lod = NULL, k = 2, accuracy_function = NULL,
                            profile = FALSE,
                            profile_model = c("linear", "quadratic"),
                            lot_rule = c("ep17", "separate", "combined"),
                            exclude = NULL, ...) {
  method <- match.arg(method)
  goal_scale <- match.arg(goal_scale)
  profile_model <- match.arg(profile_model)
  lot_rule <- match.arg(lot_rule)
  goal <- .sensitivity_number(goal, "goal", 0, strict = TRUE)
  k <- .sensitivity_number(k, "k", 0)
  if (method == "custom" && !is.function(accuracy_function)) {
    stop("`accuracy_function` is required for custom LoQ.", call. = FALSE)
  }
  p <- .sensitivity_prepare(data, result, sample, lot, exclude = exclude)
  av <- .sensitivity_assigned(p$data, assigned)
  p$data$.sensitivity_assigned <- av$value
  grouping <- .sensitivity_lot_groups(p$data, p$lot, lot_rule)
  summaries <- list()
  profiles <- list()
  candidates <- lapply(names(grouping$groups), function(label) {
    d <- grouping$groups[[label]]
    key <- interaction(d[[p$sample]], d$.sensitivity_assigned, drop = TRUE)
    pieces <- split(d, key, drop = TRUE)
    sm <- do.call(rbind, lapply(pieces, function(z) {
      y <- z[[p$result]]
      assigned_value <- z$.sensitivity_assigned[1L]
      measured_mean <- mean(y)
      measured_sd <- stats::sd(y)
      bias <- measured_mean - assigned_value
      raw_metric <- switch(method,
        westgard = abs(bias) + k * measured_sd,
        rms = sqrt(measured_sd^2 + bias^2),
        precision = measured_sd,
        custom = accuracy_function(mean = measured_mean, sd = measured_sd,
                                   assigned = assigned_value, n = length(y), ...))
      if (!is.numeric(raw_metric) || length(raw_metric) != 1L ||
          !is.finite(raw_metric)) {
        stop("The accuracy metric must be one finite numeric value per sample.",
             call. = FALSE)
      }
      metric <- if (goal_scale == "relative") raw_metric / abs(assigned_value) else
        raw_metric
      data.frame(lot = label, sample = as.character(z[[p$sample]][1L]),
                 assigned = assigned_value, mean = measured_mean, sd = measured_sd,
                 n = length(y), bias = bias,
                 relative_bias = bias / assigned_value,
                 raw_metric = raw_metric, metric = metric,
                 goal = goal, within_goal = metric <= goal)
    }))
    sm <- sm[order(sm$mean), , drop = FALSE]
    summaries[[label]] <<- sm
    observed_candidate <- if (any(sm$within_goal)) min(sm$mean[sm$within_goal]) else NA_real_
    profile_candidate <- NA_real_
    if (isTRUE(profile)) {
      degree <- if (profile_model == "linear") 1L else 2L
      if (nrow(sm) <= degree + 1L) {
        warning("Too few levels to fit the requested LoQ error profile.", call. = FALSE)
      } else {
        fit <- if (profile_model == "linear") stats::lm(metric ~ assigned, sm) else
          stats::lm(metric ~ assigned + I(assigned^2), sm)
        root <- .sensitivity_root(function(x) as.numeric(stats::predict(
          fit, newdata = data.frame(assigned = x))) - goal, sm$assigned)
        grid <- seq(min(sm$assigned), max(sm$assigned), length.out = 300L)
        profiles[[label]] <<- list(fit = fit,
          curve = data.frame(lot = label, assigned = grid,
            metric = as.numeric(stats::predict(fit,
              newdata = data.frame(assigned = grid)))))
        profile_candidate <- root
      }
    }
    data.frame(lot = label, observed_candidate = observed_candidate,
               profile_candidate = profile_candidate,
               candidate_loq = if (isTRUE(profile) && is.finite(profile_candidate))
                 profile_candidate else observed_candidate)
  })
  by_lot <- do.call(rbind, candidates)
  finite_candidates <- by_lot$candidate_loq[is.finite(by_lot$candidate_loq)]
  candidate_loq <- if (length(finite_candidates)) max(finite_candidates) else NA_real_
  lod_value <- if (is.null(lod)) NA_real_ else
    .sensitivity_extract_limit(lod, "lod", "lod")
  reported_loq <- if (is.finite(candidate_loq) && is.finite(lod_value))
    max(candidate_loq, lod_value) else candidate_loq
  .sensitivity_new("loq", method, estimate = reported_loq,
    candidate_loq = candidate_loq, reported_loq = reported_loq,
    lod = lod_value, goal = goal, goal_scale = goal_scale, k = k,
    lot_rule = lot_rule, lot_mode = grouping$mode, by_lot = by_lot,
    summary = do.call(rbind, summaries), profiles = profiles,
    data = p$data, result = p$result, sample = p$sample, lot = p$lot,
    assigned = av$name, excluded = p$excluded, incomplete = p$incomplete)
}

.sensitivity_wilson <- function(successes, n, alpha) {
  z <- stats::qnorm(1 - alpha)
  p <- successes / n
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  c(lower = max(0, center - half), upper = min(1, center + half))
}

.sensitivity_boundary <- function(n) {
  table <- data.frame(
    N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 400, 500, 1000),
    boundary = c(.85, .87, .88, .88, .90, .90, .90, .91, .91, .92,
                 .92, .92, .93, .93, .93, .94)
  )
  exact <- which(table$N == n)
  if (length(exact)) return(list(value = table$boundary[exact], rows = table[exact, ]))
  lower_candidates <- which(table$N < n)
  upper_candidates <- which(table$N > n)
  lower <- if (length(lower_candidates)) max(lower_candidates) else NA_integer_
  upper <- if (length(upper_candidates)) min(upper_candidates) else NA_integer_
  if (is.na(lower)) lower <- upper
  if (is.na(upper)) upper <- lower
  rows <- unique(table[c(lower, upper), , drop = FALSE])
  list(value = max(rows$boundary), rows = rows)
}

#' Verify an LoB, LoD, or LoQ claim
#'
#' Calculates the proportion of verification results consistent with a stated
#' claim, its one-sided Wilson score limits, and the CLSI EP17-A2 Table 1
#' boundary. No pass/fail field or regulatory conclusion is produced.
#'
#' @param data Replicate-level verification data.
#' @param result,sample Column names for result and sample ID.
#' @param limit One of `"lob"`, `"lod"`, or `"loq"`.
#' @param claim Stated limit being verified.
#' @param lob LoB threshold required for LoD verification.
#' @param assigned Assigned-value column or vector for LoQ verification.
#' @param allowable_error Accuracy window half-width for LoQ verification.
#' @param error_scale Relative or absolute LoQ error window.
#' @param alpha One-sided confidence error probability.
#' @param exclude Optional excluded rows.
#' @return An object of class `sensitivity_verification`.
#' @references CLSI. Evaluation of Detection Capability for Clinical
#'   Laboratory Measurement Procedures. EP17-A2, 2nd ed. 2012.
#' @examples
#' verification <- data.frame(
#'   sample = rep(c("blank_1", "blank_2"), each = 10),
#'   result = c(rep(0.03, 9), 0.07, rep(0.04, 9), 0.08)
#' )
#' verified <- sensitivity_verify(
#'   verification, "result", "sample", limit = "lob", claim = 0.05
#' )
#' verified$result
#' @export
sensitivity_verify <- function(data, result, sample,
                               limit = c("lob", "lod", "loq"), claim,
                               lob = NULL, assigned = NULL,
                               allowable_error = NULL,
                               error_scale = c("relative", "absolute"),
                               alpha = 0.05, exclude = NULL) {
  limit <- match.arg(limit)
  error_scale <- match.arg(error_scale)
  claim <- .sensitivity_number(claim, "claim")
  alpha <- .sensitivity_probability(alpha, "alpha")
  p <- .sensitivity_prepare(data, result, sample, exclude = exclude)
  y <- p$data[[p$result]]
  lower <- upper <- rep(NA_real_, length(y))
  threshold <- NA_real_
  if (limit == "lob") {
    consistent <- y <= claim
    threshold <- claim
  } else if (limit == "lod") {
    if (is.null(lob)) stop("`lob` is required for LoD verification.", call. = FALSE)
    threshold <- .sensitivity_extract_limit(lob, "lob", "lob")
    consistent <- y >= threshold
  } else {
    if (is.null(assigned) || is.null(allowable_error)) {
      stop("`assigned` and `allowable_error` are required for LoQ verification.",
           call. = FALSE)
    }
    error <- .sensitivity_number(allowable_error, "allowable_error", 0,
                                 strict = TRUE)
    av <- .sensitivity_assigned(p$data, assigned)
    if (error_scale == "relative") {
      lower <- av$value * (1 - error)
      upper <- av$value * (1 + error)
      swap <- lower > upper
      if (any(swap)) {
        tmp <- lower[swap]; lower[swap] <- upper[swap]; upper[swap] <- tmp
      }
    } else {
      lower <- av$value - error
      upper <- av$value + error
    }
    consistent <- y >= lower & y <= upper
  }
  n <- length(consistent)
  successes <- sum(consistent)
  ci <- .sensitivity_wilson(successes, n, alpha)
  boundary <- .sensitivity_boundary(n)
  row_data <- data.frame(sample = as.character(p$data[[p$sample]]), result = y,
                         lower = lower, upper = upper, consistent = consistent)
  by_sample <- do.call(rbind, lapply(split(row_data, row_data$sample), function(z) {
    data.frame(sample = z$sample[1L], n_total = nrow(z),
               n_consistent = sum(z$consistent),
               observed_proportion = mean(z$consistent))
  }))
  rownames(by_sample) <- NULL
  result_table <- data.frame(limit = limit, claim = claim, lob = threshold,
    n_total = n, n_consistent = successes,
    observed_proportion = successes / n,
    score_ci_lower = ci["lower"], score_ci_upper = ci["upper"],
    ep17_boundary = boundary$value,
    difference_from_boundary = successes / n - boundary$value,
    row.names = NULL)
  .sensitivity_new("verification", limit, result = result_table,
    limit = limit, claim = claim, lob = threshold, alpha = alpha,
    error_scale = if (limit == "loq") error_scale else NA_character_,
    boundary_rows = boundary$rows, by_sample = by_sample,
    observations = row_data, data = p$data, result_col = p$result,
    sample = p$sample, excluded = p$excluded, incomplete = p$incomplete)
}

#' Summarize an analytical-sensitivity study design
#'
#' Returns minimum CLSI EP17-A2 design requirements and, optionally, counts
#' observed in a supplied data set. Differences are descriptive only; the
#' function does not issue a pass/fail conclusion.
#'
#' @param method One or more of `"classical"`, `"profile"`, `"probit"`,
#'   `"loq"`, `"verification"`, or `"all"`.
#' @param data Optional study data.
#' @param result,sample,lot,day,concentration Optional column names used to
#'   summarize actual study counts.
#' @param limit Detection limit for a verification design.
#' @return An object of class `sensitivity_design`.
#' @references CLSI. Evaluation of Detection Capability for Clinical
#'   Laboratory Measurement Procedures. EP17-A2, 2nd ed. 2012.
#' @examples
#' sensitivity_design(c("classical", "loq"))
#'
#' design_data <- data.frame(
#'   result = 1:8, sample = rep(c("S1", "S2"), each = 4),
#'   lot = rep(c("A", "B"), 4), day = rep(1:2, each = 4)
#' )
#' sensitivity_design("verification", data = design_data,
#'                    result = "result", sample = "sample",
#'                    lot = "lot", day = "day", limit = "lod")
#' @export
sensitivity_design <- function(method = c("classical", "profile", "probit",
                                           "loq", "verification"),
                               data = NULL, result = NULL, sample = NULL,
                               lot = NULL, day = NULL, concentration = NULL,
                               limit = c("lob", "lod", "loq")) {
  allowed <- c("classical", "profile", "probit", "loq", "verification")
  if (identical(method, "all")) method <- allowed
  if (!is.character(method) || !length(method) || any(!method %in% allowed)) {
    stop("Unknown analytical-sensitivity design method.", call. = FALSE)
  }
  method <- unique(method)
  limit <- match.arg(limit)
  requirements <- data.frame(
    method = allowed,
    reagent_lots = c(2, 2, 2, 2, 1),
    instruments = c(1, 1, NA, 1, 1),
    days = c(3, 5, NA, 3, 3),
    blank_samples = c(4, 0, 0, 0, if (limit == "lob") 2 else 0),
    low_samples = c(4, 5, 3, 4, if (limit == "lob") 0 else 2),
    concentration_levels = c(NA, 5, 5, 1, 1),
    replicates_per_lot = c(60, 40, 20, 36, 20),
    replicate_scope = c("blank and low results separately",
      "per sample", "per sample and concentration level", "total low results",
      "total verification results"),
    stringsAsFactors = FALSE
  )
  requirements <- requirements[match(method, requirements$method), , drop = FALSE]
  actual <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
    columns <- list(result = result, sample = sample, lot = lot, day = day,
                    concentration = concentration)
    for (nm in names(columns)) {
      if (!is.null(columns[[nm]])) .sensitivity_column(data, columns[[nm]], nm)
    }
    count_unique <- function(column) if (is.null(column)) NA_integer_ else
      length(unique(data[[column]][!is.na(data[[column]])]))
    actual <- data.frame(
      observations = if (is.null(result)) nrow(data) else sum(!is.na(data[[result]])),
      reagent_lots = count_unique(lot), days = count_unique(day),
      samples = count_unique(sample), concentration_levels = count_unique(concentration),
      stringsAsFactors = FALSE
    )
    if (!is.null(lot)) {
      actual$observations_per_lot <- list(as.data.frame(table(data[[lot]]),
        stringsAsFactors = FALSE))
    }
  }
  structure(list(call = match.call(), methods = method, limit = limit,
                 requirements = requirements, actual = actual),
            class = "sensitivity_design")
}

#' Print an analytical-sensitivity analysis result
#'
#' @param x A `sensitivity_lob`, `sensitivity_lod`,
#'   `sensitivity_loq`, or `sensitivity_verification` object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.sensitivity <- function(x, ...) {
  cat("Analytical sensitivity:", toupper(x$type), "\n")
  cat("  Method:", x$method, "\n")
  if (!is.null(x$estimate)) cat("  Estimate:", format(x$estimate, digits = 6), "\n")
  if (x$type == "loq") {
    cat("  Candidate LoQ:", format(x$candidate_loq, digits = 6), "\n")
    cat("  Reported LoQ :", format(x$reported_loq, digits = 6), "\n")
  }
  if (x$type == "verification") {
    if (nrow(x$result) == 1L) {
      .print_kv_sections(x$result, list(
        "Verification summary" = c(
          "Limit" = "limit", "Claim" = "claim", "LoB" = "lob",
          "Total results" = "n_total", "Consistent results" = "n_consistent",
          "Observed proportion" = "observed_proportion"),
        "Confidence interval and EP17 boundary" = c(
          "Score CI lower" = "score_ci_lower", "Score CI upper" = "score_ci_upper",
          "EP17 boundary" = "ep17_boundary",
          "Difference from boundary" = "difference_from_boundary")
      ))
    } else {
      .print_kv_rows(x$result, list(
        "Verification summary" = c(
          "Claim" = "claim", "LoB" = "lob", "Total results" = "n_total",
          "Consistent results" = "n_consistent",
          "Observed proportion" = "observed_proportion"),
        "Confidence interval and EP17 boundary" = c(
          "Score CI lower" = "score_ci_lower", "Score CI upper" = "score_ci_upper",
          "EP17 boundary" = "ep17_boundary",
          "Difference from boundary" = "difference_from_boundary")
      ), row_title = c("Limit" = "limit"))
    }
  }
  invisible(x)
}

#' Print an analytical-sensitivity study design
#'
#' @param x A `sensitivity_design` object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.sensitivity_design <- function(x, ...) {
  cat("Analytical sensitivity study design\n")
  if (nrow(x$requirements) == 1L) {
    .print_kv_sections(x$requirements, list(
      "Study setup" = c(
        "Method" = "method", "Reagent lots" = "reagent_lots",
        "Instruments" = "instruments", "Days" = "days"),
      "Sample and replicate requirements" = c(
        "Blank samples" = "blank_samples", "Low samples" = "low_samples",
        "Concentration levels" = "concentration_levels",
        "Replicates per lot" = "replicates_per_lot"),
      "Replicate scope" = c("Replicate scope" = "replicate_scope")
    ))
  } else {
    .print_kv_rows(x$requirements, list(
      "Study setup" = c("Reagent lots" = "reagent_lots",
                         "Instruments" = "instruments", "Days" = "days"),
      "Sample counts" = c("Blank samples" = "blank_samples",
                           "Low samples" = "low_samples"),
      "Concentration and replicates" = c(
        "Concentration levels" = "concentration_levels",
        "Replicates per lot" = "replicates_per_lot"),
      "Replicate scope" = c("Replicate scope" = "replicate_scope")
    ), row_title = c("Method" = "method"))
  }
  if (!is.null(x$actual)) {
    cat("\nObserved study counts\n")
    actual <- x$actual
    per_lot <- NULL
    if ("observations_per_lot" %in% names(actual) &&
        length(actual$observations_per_lot)) {
      per_lot <- actual$observations_per_lot[[1L]]
      actual$observations_per_lot <- NULL
    }
    if (nrow(actual) == 1L) {
      .print_kv_sections(actual, list(
        "Observed counts" = c(
          "Observations" = "observations", "Reagent lots" = "reagent_lots",
          "Days" = "days", "Samples" = "samples",
          "Concentration levels" = "concentration_levels")
      ))
    } else {
      .print_kv_rows(actual, list("Observed counts" = names(actual)))
    }
    if (!is.null(per_lot)) {
      .print_kv_rows(
        per_lot,
        list("Observations" = setNames(names(per_lot)[-1L],
                                        names(per_lot)[-1L])),
        row_title = setNames(names(per_lot)[1L], "Lot"))
    }
  }
  invisible(x)
}

#' Plot an analytical-sensitivity analysis
#'
#' @param x A sensitivity-analysis result object.
#' @param ... Reserved for the S3 generic.
#' @return A \code{ggplot} object.
#' @examples
#' d <- data.frame(sample = rep("blank", 10), result = seq(0.1, 1, length.out = 10))
#' fit <- sensitivity_lob(d, "result", "sample")
#' plot(fit)
#' @export
plot.sensitivity <- function(x, ...) {
  if (x$type == "lob") {
    d <- data.frame(value = x$data[[x$result]])
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data$value)) +
      ggplot2::stat_ecdf(geom = "step") +
      ggplot2::geom_vline(xintercept = x$estimate, linetype = "dashed",
                          color = "darkred") +
      ggplot2::labs(x = "Result", y = "Empirical cumulative probability",
                    title = "Limit of blank") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (x$type == "lod" && x$method == "profile") {
    curve <- do.call(rbind, x$curves)
    return(ggplot2::ggplot(curve, ggplot2::aes(x = .data$concentration,
                                               y = .data$trial_lod,
                                               color = .data$lot)) +
      ggplot2::geom_line() + ggplot2::geom_abline(slope = 1, intercept = 0,
        linetype = "dashed") +
      ggplot2::labs(x = "Concentration", y = "Trial LoD",
                    title = "LoD precision profile") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (x$type == "lod" && x$method == "probit") {
    observed <- do.call(rbind, lapply(x$rates, `[[`, "observed"))
    curve <- do.call(rbind, lapply(x$rates, `[[`, "curve"))
    return(ggplot2::ggplot() +
      ggplot2::geom_point(data = observed,
        ggplot2::aes(x = .data$concentration, y = .data$observed_proportion,
                     color = .data$group)) +
      ggplot2::geom_line(data = curve,
        ggplot2::aes(x = .data$concentration, y = .data$probability,
                     color = .data$group)) +
      ggplot2::geom_hline(yintercept = 1 - x$beta, linetype = "dashed") +
      ggplot2::labs(x = "Concentration", y = "Detection probability",
                    title = "Probit LoD") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (x$type == "lod") {
    sm <- do.call(rbind, lapply(names(x$per_sample), function(nm) {
      z <- x$per_sample[[nm]]; z$lot <- nm; z
    }))
    return(ggplot2::ggplot(sm, ggplot2::aes(x = .data$sample, y = .data$sd,
                                             color = .data$lot)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::labs(x = "Low-level sample", y = "SD", title = "Classical LoD") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (x$type == "loq") {
    return(ggplot2::ggplot(x$summary,
      ggplot2::aes(x = .data$assigned, y = .data$metric, color = .data$lot)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_hline(yintercept = x$goal, linetype = "dashed") +
      ggplot2::labs(x = "Assigned value", y = "Accuracy metric",
                    title = "Limit of quantitation") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  if (x$type == "verification") {
    return(ggplot2::ggplot(x$by_sample,
      ggplot2::aes(x = .data$sample, y = .data$observed_proportion)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_hline(yintercept = x$result$ep17_boundary,
                          linetype = "dashed", color = "darkred") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Sample", y = "Observed proportion",
                    title = "Claim verification") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  stop("No plot is available for this sensitivity object.", call. = FALSE)
}
