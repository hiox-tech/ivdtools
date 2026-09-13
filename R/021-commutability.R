# Commutability assessment -------------------------------------------------

.comm_col <- function(data, column, argument, required = TRUE) {
  if (is.null(column) && !required) return(NULL)
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column) || !column %in% names(data)) {
    stop("'", argument, "' must be a single column name in 'data'.",
         call. = FALSE)
  }
  data[[column]]
}

.comm_material_type <- function(type, material_type_map = NULL) {
  raw_type <- as.character(type)
  if (anyNA(raw_type))
    stop("Material type values must be complete.", call. = FALSE)

  if (is.null(material_type_map)) {
    normalized <- tolower(raw_type)
    normalized[normalized %in% c("patient", "native")] <- "clinical"
    normalized[normalized == "reference_material"] <- "rm"
  } else {
    valid_map <- is.character(material_type_map) &&
      length(material_type_map) > 0L && !is.null(names(material_type_map)) &&
      !anyNA(material_type_map) && !anyNA(names(material_type_map)) &&
      all(nzchar(names(material_type_map))) && !anyDuplicated(names(material_type_map))
    if (!valid_map) {
      stop("'material_type_map' must be a nonempty named character vector with unique, nonempty names.",
           call. = FALSE)
    }
    mapped_values <- tolower(unname(material_type_map))
    if (!all(mapped_values %in% c("clinical", "rm"))) {
      stop("All 'material_type_map' values must be 'clinical' or 'rm'.",
           call. = FALSE)
    }
    normalized <- mapped_values[match(raw_type, names(material_type_map))]
  }

  unknown <- unique(raw_type[is.na(normalized) |
    !normalized %in% c("clinical", "rm")])
  if (length(unknown)) {
    quoted <- paste0("'", unknown, "'")
    stop("Unrecognized material type value(s): ", paste(quoted, collapse = ", "),
         ". Supply 'material_type_map' as a named character vector mapping each observed value to 'clinical' or 'rm'.",
         call. = FALSE)
  }
  normalized
}

.comm_transform <- function(x, scale) {
  if (scale == "raw") return(x)
  if (any(x <= 0))
    stop("All results must be positive for logarithmic analysis.", call. = FALSE)
  if (scale == "log10") log10(x) else log(x)
}

.comm_summarise <- function(d) {
  keys <- unique(d[c("sample", "type", "procedure")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    z <- d$result[d$sample == keys$sample[i] &
                  d$type == keys$type[i] &
                  d$procedure == keys$procedure[i]]
    data.frame(sample = keys$sample[i], type = keys$type[i],
               procedure = keys$procedure[i], n = length(z), mean = mean(z),
               sd = if (length(z) > 1L) stats::sd(z) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.comm_pair_data <- function(summary, x_name, y_name, type) {
  x <- summary[summary$type == type & summary$procedure == x_name,
               c("sample", "n", "mean", "sd")]
  y <- summary[summary$type == type & summary$procedure == y_name,
               c("sample", "n", "mean", "sd")]
  names(x)[-1L] <- paste0(c("n", "mean", "sd"), "_x")
  names(y)[-1L] <- paste0(c("n", "mean", "sd"), "_y")
  merge(x, y, by = "sample", all = FALSE, sort = FALSE)
}

.comm_pairs <- function(procedures, procedure_pairs, reference_procedure) {
  if (!is.null(reference_procedure)) {
    if (!is.character(reference_procedure) || length(reference_procedure) != 1L ||
        !reference_procedure %in% procedures) {
      stop("'reference_procedure' must identify one observed procedure.",
           call. = FALSE)
    }
    others <- setdiff(procedures, reference_procedure)
    if (!length(others)) stop("At least two procedures are required.", call. = FALSE)
    return(lapply(others, function(z) c(reference_procedure, z)))
  }
  if (is.null(procedure_pairs)) {
    if (length(procedures) < 2L)
      stop("At least two procedures are required.", call. = FALSE)
    return(utils::combn(procedures, 2L, simplify = FALSE))
  }
  if (is.matrix(procedure_pairs) || is.data.frame(procedure_pairs)) {
    procedure_pairs <- lapply(seq_len(nrow(procedure_pairs)),
                              function(i) as.character(procedure_pairs[i, ]))
  }
  if (!is.list(procedure_pairs) ||
      any(vapply(procedure_pairs, length, integer(1)) != 2L)) {
    stop("'procedure_pairs' must be a list of two-procedure character vectors.",
         call. = FALSE)
  }
  procedure_pairs <- lapply(procedure_pairs, as.character)
  if (any(!unlist(procedure_pairs, use.names = FALSE) %in% procedures))
    stop("Every procedure in 'procedure_pairs' must occur in the data.",
         call. = FALSE)
  procedure_pairs
}

.comm_pooled_variance <- function(pair, side) {
  n <- pair[[paste0("n_", side)]]
  sd <- pair[[paste0("sd_", side)]]
  ok <- n > 1L & is.finite(sd)
  df <- sum(n[ok] - 1L)
  if (df <= 0L) return(c(variance = NA_real_, df = 0))
  c(variance = sum((n[ok] - 1L) * sd[ok]^2) / df, df = df)
}

.comm_range_text <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("not available")
  limits <- range(x)
  if (limits[1L] == limits[2L]) return(as.character(limits[1L]))
  paste(limits, collapse = "-")
}

.comm_design_check <- function(pair, requirement, observed, recommended,
                               met, source) {
  data.frame(
    pair = pair, requirement = requirement,
    observed = as.character(observed), recommended = recommended,
    met = isTRUE(met), source = source, stringsAsFactors = FALSE
  )
}

.comm_bind <- function(x, name) {
  rows <- lapply(x, `[[`, name)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

.comm_deming <- function(x, y, lambda) {
  n <- length(x)
  xb <- mean(x); yb <- mean(y)
  sxx <- sum((x - xb)^2) / n
  syy <- sum((y - yb)^2) / n
  sxy <- sum((x - xb) * (y - yb)) / n
  if (!is.finite(sxy) || sxy == 0)
    stop("Deming regression cannot be fitted because covariance is zero.",
         call. = FALSE)
  disc <- (syy - lambda * sxx)^2 + 4 * lambda * sxy^2
  slope <- (syy - lambda * sxx + sign(sxy) * sqrt(disc)) / (2 * sxy)
  intercept <- yb - slope * xb
  var_slope <- slope^2 * max(0, sxx * syy - sxy^2) / (n * sxy^2)
  list(intercept = intercept, slope = slope, variance_slope = var_slope,
       x_mean = xb, y_mean = yb, sxx = sxx, syy = syy, sxy = sxy, n = n)
}

.comm_ep14_pair <- function(summary, pair, conf.level) {
  x_name <- pair[1L]; y_name <- pair[2L]
  cs <- .comm_pair_data(summary, x_name, y_name, "clinical")
  rm <- .comm_pair_data(summary, x_name, y_name, "rm")
  if (nrow(cs) < 3L)
    stop("At least three complete clinical samples are needed to fit the EP14 regression for this procedure pair.",
         call. = FALSE)
  if (!nrow(rm))
    stop("No complete RM results are available for procedure pair '",
         x_name, "' and '", y_name, "'.", call. = FALSE)
  vx <- .comm_pooled_variance(cs, "x")
  vy <- .comm_pooled_variance(cs, "y")
  if (!is.finite(vx["variance"]) || !is.finite(vy["variance"]) ||
      vx["variance"] <= 0 || vy["variance"] <= 0) {
    stop("The EP14 calculation needs estimable positive pooled repeatability variance for both procedures.",
         call. = FALSE)
  }
  lambda <- unname(vy["variance"] / vx["variance"])
  fit <- .comm_deming(cs$mean_x, cs$mean_y, lambda)
  alpha <- 1 - conf.level
  df <- max(1, floor(min(vx["df"], vy["df"])))
  critical <- stats::qt(1 - alpha / 2, df)
  result <- lapply(seq_len(nrow(rm)), function(i) {
    pred <- fit$intercept + fit$slope * rm$mean_x[i]
    var_pred <- (rm$mean_x[i] - fit$x_mean)^2 * fit$variance_slope +
      (1 + 1 / fit$n) *
      (fit$slope^2 * vx["variance"] / rm$n_x[i] +
       vy["variance"] / rm$n_y[i])
    se <- sqrt(max(0, var_pred))
    lower <- pred - critical * se
    upper <- pred + critical * se
    inside <- rm$mean_y[i] >= lower && rm$mean_y[i] <= upper
    data.frame(pair = paste(x_name, y_name, sep = " vs "), rm = rm$sample[i],
      procedure_x = x_name, procedure_y = y_name, mean_x = rm$mean_x[i],
      mean_y = rm$mean_y[i], predicted_y = pred, residual = rm$mean_y[i] - pred,
      prediction_se = se, lower = lower, upper = upper,
      extrapolated = rm$mean_x[i] < min(cs$mean_x) || rm$mean_x[i] > max(cs$mean_x),
      classification = if (inside) "commutable" else "noncommutable",
      stringsAsFactors = FALSE)
  })
  pair_name <- paste(x_name, y_name, sep = " vs ")
  min_replicates <- min(c(cs$n_x, cs$n_y))
  design_checks <- rbind(
    .comm_design_check(
      pair_name, "Complete clinical samples", nrow(cs), ">= 20",
      nrow(cs) >= 20L, "CLSI EP14-A3"
    ),
    .comm_design_check(
      pair_name, "Replicates per clinical sample",
      .comm_range_text(c(cs$n_x, cs$n_y)), ">= 3",
      min_replicates >= 3L, "CLSI EP14-A3"
    )
  )
  list(
    result = do.call(rbind, result),
    model = data.frame(pair = paste(x_name, y_name, sep = " vs "),
      procedure_x = x_name, procedure_y = y_name, n_clinical = nrow(cs),
      variance_x = unname(vx["variance"]), variance_y = unname(vy["variance"]),
      lambda = lambda, intercept = fit$intercept, slope = fit$slope,
      slope_se = sqrt(fit$variance_slope), df = df, stringsAsFactors = FALSE),
    clinical = transform(cs, pair = paste(x_name, y_name, sep = " vs ")),
    fit = fit, design_checks = design_checks
  )
}

.comm_ep30_pair <- function(summary, pair, conf.level, rm_set_size,
                            relative_residual_limit) {
  x_name <- pair[1L]; y_name <- pair[2L]
  cs <- .comm_pair_data(summary, x_name, y_name, "clinical")
  rm <- .comm_pair_data(summary, x_name, y_name, "rm")
  pair_name <- paste(x_name, y_name, sep = " vs ")
  if (nrow(cs) < 3L)
    stop("At least three complete clinical samples are needed to fit the EP30 regression for this procedure pair.",
         call. = FALSE)
  if (!nrow(rm))
    stop("No complete RM results are available for procedure pair '",
         x_name, "' and '", y_name, "'.", call. = FALSE)
  clinical_replicates <- c(cs$n_x, cs$n_y)
  rm_replicates <- c(rm$n_x, rm$n_y)
  clinical_replicates_equal <- length(unique(clinical_replicates)) == 1L
  rm_replicates_equal <- length(unique(rm_replicates)) == 1L
  vx <- .comm_pooled_variance(cs, "x")
  vy <- .comm_pooled_variance(cs, "y")
  rm_vx <- .comm_pooled_variance(rm, "x")
  rm_vy <- .comm_pooled_variance(rm, "y")
  variances <- c(vx["variance"], vy["variance"],
                 rm_vx["variance"], rm_vy["variance"])
  if (any(!is.finite(variances)) || any(variances <= 0)) {
    stop("The EP30 prediction calculation needs estimable positive pooled repeatability variance for clinical samples and reference materials in both procedures.",
         call. = FALSE)
  }
  lambda <- unname(vy["variance"] / vx["variance"])
  fit <- .comm_deming(cs$mean_x, cs$mean_y, lambda)
  residual_variance <- fit$syy - 2 * fit$slope * fit$sxy +
    fit$slope^2 * fit$sxx
  variance_slope <- fit$variance_slope
  variance_intercept <- residual_variance / fit$n +
    fit$x_mean^2 * variance_slope
  covariance_intercept_slope <- -fit$x_mean * variance_slope
  syx <- sqrt(sum((fit$intercept + fit$slope * cs$mean_x -
                    cs$mean_y)^2) / (fit$n - 1L))
  relative_residual_evaluable <- is.finite(syx) && syx > 0
  critical <- stats::qnorm(
    1 - (1 - conf.level) / (2 * rm_set_size)
  )
  prediction_variance <- function(x, nx, ny) {
    variance_intercept + 2 * x * covariance_intercept_slope +
      x^2 * variance_slope +
      fit$slope^2 * unname(rm_vx["variance"]) / nx +
      unname(rm_vy["variance"]) / ny
  }
  results <- lapply(seq_len(nrow(rm)), function(i) {
    predicted <- fit$intercept + fit$slope * rm$mean_x[i]
    variance <- prediction_variance(rm$mean_x[i], rm$n_x[i], rm$n_y[i])
    if (!is.finite(variance) || variance < -sqrt(.Machine$double.eps)) {
      stop("EP30 C16 produced an invalid prediction variance.",
           call. = FALSE)
    }
    prediction_se <- sqrt(max(0, variance))
    lower <- predicted - critical * prediction_se
    upper <- predicted + critical * prediction_se
    extrapolated <- rm$mean_x[i] < min(cs$mean_x) ||
      rm$mean_x[i] > max(cs$mean_x) ||
      rm$mean_y[i] < min(cs$mean_y) || rm$mean_y[i] > max(cs$mean_y)
    prediction_classification <- if (extrapolated) {
      "not_evaluable"
    } else if (rm$mean_y[i] >= lower && rm$mean_y[i] <= upper) {
      "commutable"
    } else {
      "noncommutable"
    }
    relative_residual <- if (relative_residual_evaluable)
      (predicted - rm$mean_y[i]) / syx else NA_real_
    relative_classification <- if (extrapolated ||
                                        !relative_residual_evaluable) {
      "not_evaluable"
    } else if (abs(relative_residual) <= relative_residual_limit) {
      "commutable"
    } else {
      "noncommutable"
    }
    data.frame(
      pair = pair_name, rm = rm$sample[i], procedure_x = x_name,
      procedure_y = y_name, mean_x = rm$mean_x[i], mean_y = rm$mean_y[i],
      predicted_y = predicted, residual = rm$mean_y[i] - predicted,
      prediction_se = prediction_se, lower = lower, upper = upper,
      critical_value = critical, rm_set_size = rm_set_size,
      extrapolated = extrapolated,
      prediction_classification = prediction_classification,
      syx = syx, relative_residual = relative_residual,
      relative_residual_limit = relative_residual_limit,
      relative_residual_classification = relative_classification,
      classification = prediction_classification,
      stringsAsFactors = FALSE
    )
  })
  grid_x <- seq(min(cs$mean_x), max(cs$mean_x), length.out = 200L)
  grid_nx <- min(rm$n_x)
  grid_ny <- min(rm$n_y)
  grid_variance <- prediction_variance(
    grid_x, rep(grid_nx, length(grid_x)), rep(grid_ny, length(grid_x))
  )
  if (any(!is.finite(grid_variance)) ||
      any(grid_variance < -sqrt(.Machine$double.eps))) {
    stop("EP30 C16 produced an invalid prediction-grid variance.",
         call. = FALSE)
  }
  grid_fitted <- fit$intercept + fit$slope * grid_x
  grid_se <- sqrt(pmax(0, grid_variance))
  design_checks <- rbind(
    .comm_design_check(
      pair_name, "Complete clinical samples", nrow(cs), ">= 20",
      nrow(cs) >= 20L, "CLSI EP30-A Appendix D"
    ),
    .comm_design_check(
      pair_name, "Replicates per clinical sample",
      .comm_range_text(clinical_replicates), ">= 3",
      min(clinical_replicates) >= 3L, "CLSI EP30-A Appendix D"
    ),
    .comm_design_check(
      pair_name, "Equal clinical-sample replicate counts",
      if (clinical_replicates_equal) "yes" else "no", "yes",
      clinical_replicates_equal, "CLSI EP30-A Appendix D"
    ),
    .comm_design_check(
      pair_name, "Equal reference-material replicate counts",
      if (rm_replicates_equal) "yes" else "no", "yes",
      rm_replicates_equal, "CLSI EP30-A Appendix D"
    )
  )
  list(
    result = do.call(rbind, results),
    model = data.frame(
      pair = pair_name, procedure_x = x_name, procedure_y = y_name,
      n_clinical = nrow(cs),
      clinical_replicates = if (clinical_replicates_equal)
        clinical_replicates[1L] else NA_integer_,
      clinical_replicates_min = min(clinical_replicates),
      clinical_replicates_max = max(clinical_replicates),
      clinical_replicate_range = .comm_range_text(clinical_replicates),
      clinical_replicates_equal = clinical_replicates_equal,
      rm_replicates = if (rm_replicates_equal)
        rm_replicates[1L] else NA_integer_,
      rm_replicates_min = min(rm_replicates),
      rm_replicates_max = max(rm_replicates),
      rm_replicate_range = .comm_range_text(rm_replicates),
      rm_replicates_equal = rm_replicates_equal,
      variance_x = unname(vx["variance"]),
      variance_y = unname(vy["variance"]), lambda = lambda,
      intercept = fit$intercept, slope = fit$slope,
      variance_intercept = variance_intercept,
      variance_slope = variance_slope,
      covariance_intercept_slope = covariance_intercept_slope,
      slope_se = sqrt(variance_slope),
      rm_variance_x = unname(rm_vx["variance"]),
      rm_variance_y = unname(rm_vy["variance"]),
      syx = syx, relative_residual_evaluable = relative_residual_evaluable,
      critical_value = critical, rm_set_size = rm_set_size,
      relative_residual_limit = relative_residual_limit,
      stringsAsFactors = FALSE
    ),
    clinical = transform(cs, pair = pair_name),
    prediction_grid = data.frame(
      pair = pair_name, x = grid_x, fitted = grid_fitted,
      prediction_se = grid_se,
      lower = grid_fitted - critical * grid_se,
      upper = grid_fitted + critical * grid_se,
      replicates_x = grid_nx, replicates_y = grid_ny,
      stringsAsFactors = FALSE
    ),
    fit = fit, design_checks = design_checks
  )
}

.comm_position_stats <- function(d, rm_id, procedure) {
  z <- d[d$type == "rm" & d$sample == rm_id & d$procedure == procedure, ]
  if (!nrow(z) || all(is.na(z$position))) {
    return(list(
      p = 0L, total_mean = NA_real_, position_mean_sd = NA_real_,
      within_sd = NA_real_, means = numeric(0), k = integer(0),
      positions = character(0)
    ))
  }
  split_z <- split(z$result, z$position, drop = TRUE)
  split_z <- split_z[lengths(split_z) > 0L]
  means <- vapply(split_z, mean, numeric(1))
  vars <- vapply(split_z, function(v) if (length(v) > 1L) stats::var(v) else NA_real_,
                  numeric(1))
  k <- lengths(split_z)
  pooled <- if (sum(k - 1L) > 0L)
    sum((k - 1L) * vars, na.rm = TRUE) / sum(k - 1L) else NA_real_
  list(p = length(means), total_mean = mean(unlist(split_z, use.names = FALSE)),
       position_mean_sd = if (length(means) > 1L) stats::sd(means) else NA_real_,
       within_sd = sqrt(pooled), means = means, k = k,
       positions = names(split_z))
}

.comm_position_pool <- function(d, rm_ids, procedure) {
  rm_ids <- unique(as.character(rm_ids))
  stats <- lapply(rm_ids, function(rm_id) {
    .comm_position_stats(d, rm_id, procedure)
  })
  names(stats) <- rm_ids
  by_rm <- do.call(rbind, lapply(rm_ids, function(rm_id) {
    z <- stats[[rm_id]]
    data.frame(
      rm = rm_id, procedure = procedure, positions = z$p,
      position_df = z$p - 1L, position_mean_sd = z$position_mean_sd,
      within_sd = z$within_sd, stringsAsFactors = FALSE
    )
  }))
  valid <- by_rm$position_df > 0L & is.finite(by_rm$position_mean_sd)
  pooled_df <- sum(by_rm$position_df[valid])
  pooled_variance <- if (pooled_df > 0L) {
    sum(by_rm$position_df[valid] * by_rm$position_mean_sd[valid]^2) /
      pooled_df
  } else {
    NA_real_
  }
  list(
    stats = stats,
    by_rm = by_rm,
    pooled_position_mean_sd = sqrt(pooled_variance),
    pooled_df = pooled_df
  )
}

.comm_ifcc_class <- function(d, u, criterion) {
  lo <- d - u; hi <- d + u
  if (lo >= -criterion && hi <= criterion) return("commutable")
  if (lo > criterion || hi < -criterion) return("noncommutable")
  "inconclusive"
}

.comm_ifcc_pair <- function(d, summary, pair, criterion, bias_model,
                            local_n, coverage_factor, reference_procedure,
                            rm_position) {
  x_name <- pair[1L]; y_name <- pair[2L]
  cs <- .comm_pair_data(summary, x_name, y_name, "clinical")
  rm <- .comm_pair_data(summary, x_name, y_name, "rm")
  if (nrow(cs) < 2L)
    stop("At least two complete clinical samples are needed to estimate clinical-sample bias variation for this IFCC procedure pair.",
         call. = FALSE)
  if (!nrow(rm))
    stop("No complete RM results are available for procedure pair '",
         x_name, "' and '", y_name, "'.", call. = FALSE)
  cs$coordinate <- if (!is.null(reference_procedure) && x_name == reference_procedure)
    cs$mean_x else (cs$mean_x + cs$mean_y) / 2
  cs$bias <- cs$mean_y - cs$mean_x
  cs <- cs[order(cs$coordinate), ]
  n <- nrow(cs)
  s_mssd <- sqrt(sum(diff(cs$bias)^2) / (2 * (n - 1L)))
  s_b <- stats::sd(cs$bias)
  vx <- .comm_pooled_variance(cs, "x")
  vy <- .comm_pooled_variance(cs, "y")
  k_eff <- 2 / mean(1 / cs$n_x + 1 / cs$n_y)
  random_mean_variance <- if (is.finite(vx["variance"]) &&
                              is.finite(vy["variance"]))
    mean(vx["variance"] / cs$n_x + vy["variance"] / cs$n_y) else NA_real_
  s_d <- if (is.finite(vx["variance"]) && is.finite(vy["variance"]))
    sqrt(max(0, s_mssd^2 - random_mean_variance)) else NA_real_

  position_x <- .comm_position_pool(d, rm$sample, x_name)
  position_y <- .comm_position_pool(d, rm$sample, y_name)

  results <- lapply(seq_len(nrow(rm)), function(i) {
    px <- position_x$stats[[rm$sample[i]]]
    py <- position_y$stats[[rm$sample[i]]]
    rm_coordinate <- if (!is.null(reference_procedure) && x_name == reference_procedure)
      rm$mean_x[i] else (rm$mean_x[i] + rm$mean_y[i]) / 2
    if (rm_position == "pooled") {
      position_sd_x <- position_x$pooled_position_mean_sd
      position_sd_y <- position_y$pooled_position_mean_sd
      position_valid <- px$p >= 1L && py$p >= 1L &&
        is.finite(position_sd_x) && is.finite(position_sd_y)
    } else {
      position_sd_x <- px$position_mean_sd
      position_sd_y <- py$position_mean_sd
      position_valid <- px$p >= 2L && py$p >= 2L &&
        is.finite(position_sd_x) && is.finite(position_sd_y)
    }
    if (!position_valid) {
      return(data.frame(
        pair = paste(x_name, y_name, sep = " vs "), rm = rm$sample[i],
        procedure_x = x_name, procedure_y = y_name,
        rm_coordinate = rm_coordinate,
        rm_bias = rm$mean_y[i] - rm$mean_x[i], clinical_bias = NA_real_,
        bias_difference = NA_real_, rm_standard_uncertainty = NA_real_,
        clinical_standard_uncertainty = NA_real_,
        standard_uncertainty = NA_real_, expanded_uncertainty = NA_real_,
        lower = NA_real_, upper = NA_real_, criterion = criterion,
        n_clinical_local = 0L, position_mean_sd_x = position_sd_x,
        position_mean_sd_y = position_sd_y, positions_x = px$p,
        positions_y = py$p, position_pooling = rm_position,
        classification = "not_evaluable", stringsAsFactors = FALSE
      ))
    }
    u_rm <- sqrt(position_sd_x^2 / px$p + position_sd_y^2 / py$p)
    if (bias_model == "constant_global") {
      selected <- seq_len(n)
      scatter <- s_b
    } else {
      q <- local_n %||% min(12L, n)
      if (!is.numeric(q) || length(q) != 1L || q < 2L || q > n)
        stop("'local_n' must be between 2 and the number of clinical samples.",
             call. = FALSE)
      q <- as.integer(q)
      if (bias_model == "linear_local") {
        if (q %% 2L != 0L)
          stop("'local_n' must be even for 'linear_local'.", call. = FALSE)
        below <- which(cs$coordinate <= rm_coordinate)
        above <- which(cs$coordinate >= rm_coordinate)
        half <- q / 2L
        if (length(below) < half || length(above) < half) {
          selected <- integer(0)
        } else {
          selected <- c(utils::tail(below, half), utils::head(above, half))
        }
      } else {
        selected <- order(abs(cs$coordinate - rm_coordinate))[seq_len(q)]
      }
      if (!length(selected)) {
        return(data.frame(pair = paste(x_name, y_name, sep = " vs "),
          rm = rm$sample[i], procedure_x = x_name, procedure_y = y_name,
          rm_coordinate = rm_coordinate, rm_bias = rm$mean_y[i] - rm$mean_x[i],
          clinical_bias = NA_real_, bias_difference = NA_real_,
          rm_standard_uncertainty = u_rm,
          clinical_standard_uncertainty = NA_real_,
          standard_uncertainty = NA_real_, expanded_uncertainty = NA_real_,
          lower = NA_real_, upper = NA_real_, criterion = criterion,
          n_clinical_local = 0L, position_mean_sd_x = position_sd_x,
          position_mean_sd_y = position_sd_y, positions_x = px$p,
          positions_y = py$p, position_pooling = rm_position,
          classification = "not_evaluable",
          stringsAsFactors = FALSE))
      }
      if (min(cs$coordinate[selected]) > rm_coordinate ||
          max(cs$coordinate[selected]) < rm_coordinate) {
        return(data.frame(pair = paste(x_name, y_name, sep = " vs "),
          rm = rm$sample[i], procedure_x = x_name, procedure_y = y_name,
          rm_coordinate = rm_coordinate, rm_bias = rm$mean_y[i] - rm$mean_x[i],
          clinical_bias = NA_real_, bias_difference = NA_real_,
          rm_standard_uncertainty = u_rm,
          clinical_standard_uncertainty = NA_real_,
          standard_uncertainty = NA_real_, expanded_uncertainty = NA_real_,
          lower = NA_real_, upper = NA_real_, criterion = criterion,
          n_clinical_local = length(selected),
          position_mean_sd_x = position_sd_x,
          position_mean_sd_y = position_sd_y, positions_x = px$p,
          positions_y = py$p, position_pooling = rm_position,
          classification = "not_evaluable",
          stringsAsFactors = FALSE))
      }
      scatter <- if (bias_model == "linear_local") s_mssd else stats::sd(cs$bias[selected])
    }
    b_cs <- mean(cs$bias[selected])
    b_rm <- rm$mean_y[i] - rm$mean_x[i]
    d_rm <- b_rm - b_cs
    u_cs <- scatter / sqrt(length(selected))
    u <- sqrt(u_rm^2 + u_cs^2)
    U <- coverage_factor * u
    data.frame(pair = paste(x_name, y_name, sep = " vs "), rm = rm$sample[i],
      procedure_x = x_name, procedure_y = y_name, rm_coordinate = rm_coordinate,
      rm_bias = b_rm, clinical_bias = b_cs, bias_difference = d_rm,
      rm_standard_uncertainty = u_rm,
      clinical_standard_uncertainty = u_cs,
      standard_uncertainty = u, expanded_uncertainty = U,
      lower = d_rm - U, upper = d_rm + U, criterion = criterion,
      n_clinical_local = length(selected),
      position_mean_sd_x = position_sd_x,
      position_mean_sd_y = position_sd_y, positions_x = px$p,
      positions_y = py$p, position_pooling = rm_position,
      classification = .comm_ifcc_class(d_rm, U, criterion),
      stringsAsFactors = FALSE)
  })
  position_diagnostics <- rbind(
    transform(position_x$by_rm, side = "x",
      pooled_position_mean_sd = position_x$pooled_position_mean_sd,
      pooled_position_df = position_x$pooled_df),
    transform(position_y$by_rm, side = "y",
      pooled_position_mean_sd = position_y$pooled_position_mean_sd,
      pooled_position_df = position_y$pooled_df)
  )
  position_diagnostics$pair <- paste(x_name, y_name, sep = " vs ")
  position_counts <- c(
    vapply(position_x$stats, `[[`, integer(1), "p"),
    vapply(position_y$stats, `[[`, integer(1), "p")
  )
  replicate_counts <- c(cs$n_x, cs$n_y)
  matching_positions <- all(vapply(rm$sample, function(id) {
    setequal(position_x$stats[[id]]$positions,
             position_y$stats[[id]]$positions)
  }, logical(1)))
  pair_name <- paste(x_name, y_name, sep = " vs ")
  design_checks <- rbind(
    .comm_design_check(
      pair_name, "Complete clinical samples", n, ">= 30 usually",
      n >= 30L, "IFCC Part 2"
    ),
    .comm_design_check(
      pair_name, "Replicates per clinical-sample result",
      .comm_range_text(replicate_counts), ">= 2 (3 recommended)",
      min(replicate_counts) >= 2L, "IFCC Part 2"
    ),
    .comm_design_check(
      pair_name, "Reference-material run positions",
      .comm_range_text(position_counts), ">= 5 recommended",
      min(position_counts) >= 5L, "IFCC Part 2"
    ),
    .comm_design_check(
      pair_name, "Matching RM positions between procedures",
      if (matching_positions) "yes" else "no", "yes",
      matching_positions, "IFCC Part 2"
    )
  )
  list(
    result = do.call(rbind, results),
    diagnostics = data.frame(pair = paste(x_name, y_name, sep = " vs "),
      procedure_x = x_name, procedure_y = y_name, n_clinical = n,
      s_mssd = s_mssd, s_b = s_b, variance_x = unname(vx["variance"]),
      variance_y = unname(vy["variance"]), effective_replicates = k_eff,
      sample_specific_sd = s_d,
      pooled_position_mean_sd_x = position_x$pooled_position_mean_sd,
      pooled_position_mean_sd_y = position_y$pooled_position_mean_sd,
      pooled_position_df_x = position_x$pooled_df,
      pooled_position_df_y = position_y$pooled_df,
      position_pooling = rm_position, stringsAsFactors = FALSE),
    clinical = transform(cs, pair = paste(x_name, y_name, sep = " vs ")),
    position_diagnostics = position_diagnostics,
    design_checks = design_checks
  )
}

.comm_calibration <- function(d, criterion, exclude_procedures,
                              w3_reference_procedure,
                              calibration_replicate_summary) {
  stages <- unique(d$stage)
  if (!setequal(stages, c("before", "after")) || length(stages) != 2L) {
    stop("'calibration_stage' must contain exactly 'before' and 'after'.",
         call. = FALSE)
  }
  procedures <- unique(d$procedure)
  if (length(procedures) < 3L) {
    stop("At least three measurement procedures are needed to calculate trimmed-mean targets.",
         call. = FALSE)
  }

  keys <- unique(d[c("sample", "procedure", "stage")])
  cell_rows <- lapply(seq_len(nrow(keys)), function(i) {
    keep <- d$sample == keys$sample[i] &
      d$procedure == keys$procedure[i] & d$stage == keys$stage[i]
    values <- d$result[keep]
    if (length(values) > 1L && calibration_replicate_summary == "error") {
      stop("More than one calibration result was supplied for sample '",
           keys$sample[i], "', procedure '", keys$procedure[i],
           "', and stage '", keys$stage[i], "'.", call. = FALSE)
    }
    value <- if (calibration_replicate_summary == "median")
      stats::median(values) else mean(values)
    data.frame(
      sample = keys$sample[i], type = "clinical",
      procedure = keys$procedure[i], result = value,
      position = NA_character_, stage = keys$stage[i],
      replicates = length(values), stringsAsFactors = FALSE
    )
  })
  d <- do.call(rbind, cell_rows)
  all_samples <- unique(d$sample)
  required_cells <- 2L * length(procedures)
  cells_per_sample <- table(d$sample)
  complete_samples <- names(cells_per_sample)[cells_per_sample == required_cells]
  dropped_samples <- setdiff(all_samples, complete_samples)
  if (!length(complete_samples)) {
    stop("No clinical sample has a complete before-and-after result grid across all procedures.",
         call. = FALSE)
  }
  d <- d[d$sample %in% complete_samples, , drop = FALSE]
  samples <- unique(d$sample[d$stage == "before"])
  if (is.null(exclude_procedures)) exclude_procedures <- character(0)
  if (!is.character(exclude_procedures) || anyNA(exclude_procedures) ||
      any(!exclude_procedures %in% procedures)) {
    stop("'exclude_procedures' must contain observed procedure names.",
         call. = FALSE)
  }
  exclude_procedures <- unique(exclude_procedures)
  if (length(setdiff(procedures, exclude_procedures)) < 2L) {
    stop("At least two procedures must remain after exclusions to calculate a between-procedure range.",
         call. = FALSE)
  }
  if (is.null(w3_reference_procedure)) {
    w3_reference_procedure <- procedures[1L]
  }
  if (!is.character(w3_reference_procedure) ||
      length(w3_reference_procedure) != 1L ||
      is.na(w3_reference_procedure) ||
      !w3_reference_procedure %in% procedures) {
    stop("'w3_reference_procedure' must identify one observed procedure.",
         call. = FALSE)
  }
  targets <- vapply(samples, function(sample) {
    z <- d$result[d$stage == "before" & d$sample == sample]
    z <- sort(z)
    mean(z[2L:(length(z) - 1L)])
  }, numeric(1))
  if (any(!is.finite(targets)) || any(targets <= 0)) {
    stop("Initial trimmed-mean targets must be positive and finite.",
         call. = FALSE)
  }
  target_data <- data.frame(
    sample = names(targets), target = unname(targets),
    stringsAsFactors = FALSE
  )
  sample_results <- merge(d, target_data, by = "sample", all.x = TRUE,
                          sort = FALSE)
  sample_results$bias_percent <- 100 *
    (sample_results$result - sample_results$target) / sample_results$target
  sample_results$included <-
    !sample_results$procedure %in% exclude_procedures
  sample_order <- order(target_data$target)
  sample_rank <- setNames(seq_along(sample_order),
                          target_data$sample[sample_order])
  sample_results$sample_order <- unname(sample_rank[sample_results$sample])

  stage_order <- c("before", "after")
  result_rows <- lapply(stage_order, function(stage) {
    lapply(procedures, function(procedure) {
      z <- sample_results[sample_results$stage == stage &
                            sample_results$procedure == procedure, ]
      ordered_bias <- sort(z$bias_percent)
      w3 <- if (length(ordered_bias) >= 6L)
        ordered_bias[length(ordered_bias) - 2L] - ordered_bias[3L] else NA_real_
      trend_intercept <- trend_slope <- trend_p_value <- NA_real_
      if (nrow(z) >= 3L && length(unique(z$target)) >= 2L) {
        fit <- stats::lm(bias_percent ~ target, data = z)
        fit_summary <- summary(fit)$coefficients
        trend_intercept <- unname(stats::coef(fit)[1L])
        trend_slope <- unname(stats::coef(fit)[2L])
        if (nrow(fit_summary) >= 2L)
          trend_p_value <- unname(fit_summary[2L, 4L])
      }
      data.frame(
        procedure = procedure, calibration_stage = stage,
        n_clinical = nrow(z),
        median_bias_percent = stats::median(z$bias_percent),
        sd_bias_percent = if (nrow(z) >= 2L)
          stats::sd(z$bias_percent) else NA_real_,
        w3_percent = w3, adjusted_w3_percent = NA_real_,
        trend_intercept = trend_intercept, trend_slope = trend_slope,
        trend_p_value = trend_p_value,
        included = !procedure %in% exclude_procedures,
        stringsAsFactors = FALSE
      )
    })
  })
  results <- do.call(rbind, unlist(result_rows, recursive = FALSE))
  for (stage in stage_order) {
    ref <- results$calibration_stage == stage &
      results$procedure == w3_reference_procedure
    ratio <- results$sd_bias_percent[ref] / results$w3_percent[ref]
    if (!length(ratio) || !is.finite(ratio) || ratio <= 0) ratio <- NA_real_
    rows <- results$calibration_stage == stage
    results$adjusted_w3_percent[rows] <- results$w3_percent[rows] * ratio
  }

  range_row <- function(stage, assessment_set, excluded) {
    z <- results[results$calibration_stage == stage &
                   !results$procedure %in% excluded, ]
    minimum <- which.min(z$median_bias_percent)
    maximum <- which.max(z$median_bias_percent)
    impbr <- z$median_bias_percent[maximum] - z$median_bias_percent[minimum]
    classification <- if (stage == "before") {
      "screening"
    } else if (impbr <= criterion) {
      "commutable"
    } else {
      "requires_investigation"
    }
    data.frame(
      calibration_stage = stage, assessment_set = assessment_set,
      n_procedures = nrow(z), impbr_percent = impbr,
      min_bias_percent = z$median_bias_percent[minimum],
      max_bias_percent = z$median_bias_percent[maximum],
      min_procedure = z$procedure[minimum],
      max_procedure = z$procedure[maximum], criterion = criterion,
      excluded_procedures = if (length(excluded))
        paste(excluded, collapse = ", ") else "",
      classification = classification,
      stringsAsFactors = FALSE
    )
  }
  diagnostics <- do.call(rbind, lapply(stage_order, function(stage) {
    rows <- list(range_row(stage, "all", character(0)))
    if (length(exclude_procedures)) {
      rows[[2L]] <- range_row(stage, "specified_exclusions",
                              exclude_procedures)
    }
    do.call(rbind, rows)
  }))
  after <- results[results$calibration_stage == "after", ]
  leave_one_out <- do.call(rbind, lapply(after$procedure, function(omitted) {
    z <- after[after$procedure != omitted, ]
    impbr <- max(z$median_bias_percent) - min(z$median_bias_percent)
    data.frame(
      omitted_procedure = omitted, n_procedures = nrow(z),
      impbr_percent = impbr, criterion = criterion,
      classification = if (impbr <= criterion)
        "commutable" else "requires_investigation",
      stringsAsFactors = FALSE
    )
  }))
  max_replicates <- max(vapply(cell_rows, function(z) z$replicates,
                               integer(1)))
  input_complete <- !length(dropped_samples)
  design_checks <- rbind(
    .comm_design_check(
      "all procedures", "Measurement procedures", length(procedures),
      "3-20", length(procedures) <= 20L, "IFCC Part 3"
    ),
    .comm_design_check(
      "all procedures", "Complete clinical samples", length(samples),
      ">= 32", length(samples) >= 32L, "IFCC Part 3"
    ),
    .comm_design_check(
      "all procedures", "One final result per sample/procedure/stage",
      paste0("maximum ", max_replicates), "1",
      max_replicates == 1L, "IFCC Part 3"
    ),
    .comm_design_check(
      "all procedures", "Complete matched before/after grid",
      if (input_complete) "yes" else paste0(
        "no; ", length(dropped_samples), " sample(s) excluded"),
      "yes", input_complete, "IFCC Part 3"
    )
  )
  list(
    result = results, diagnostics = diagnostics,
    sample_results = sample_results,
    leave_one_out = leave_one_out,
    targets = target_data,
    w3_reference_procedure = w3_reference_procedure,
    design_checks = design_checks, dropped_samples = dropped_samples,
    data = d
  )
}

#' Assess reference-material commutability using four standard approaches
#'
#' Implements the statistical calculations for CLSI EP14-A3 Deming-regression
#' prediction intervals, CLSI EP30-A prediction intervals and relative
#' residuals, the IFCC Part 2 difference-in-bias approach, or the IFCC Part 3
#' calibration-effectiveness approach. The function expects long-format data
#' and does not assess specimen collection, run randomization, reagent lots, or
#' other experimental conduct.
#'
#' Sample-size, replicate-count, position-matching, and method-count targets
#' from CLSI and IFCC publications are reported in `design_checks`; they are not
#' treated as proof of compliance and do not stop a calculation that is
#' statistically defined. Conditions required by the formula itself still
#' produce an error. A local result that cannot be calculated is returned as
#' `not_evaluable`, while an unavailable secondary statistic is returned as
#' `NA`.
#'
#' @param data A data frame in long format.
#' @param sample_id,material_type,procedure,result Column names identifying the
#'   material, material type, measurement procedure, and numeric result.
#'   `material_type` may be `NULL` for `approach = "calibration"`.
#' @param material_type_map Optional named character vector that explicitly maps
#'   observed material-type values to `"clinical"` or `"rm"`, for example
#'   `c("patient sample" = "clinical", "candidate RM" = "rm")`. When supplied,
#'   this mapping is used as given. When `NULL`, the function automatically
#'   recognizes `"clinical"`, `"patient"`, and `"native"` as clinical samples,
#'   and `"rm"` and `"reference_material"` as reference materials. If automatic
#'   recognition fails, the error lists the unrecognized values and requests
#'   this argument.
#' @param position Optional column identifying RM positions within a run;
#'   required for `approach = "ifcc"`.
#' @param approach One of `"ep14"`, `"ep30"`, `"ifcc"` (IFCC Part 2), or
#'   `"calibration"` (IFCC Part 3).
#' @param scale Analysis scale: `"raw"`, `"log10"`, or `"ln"`. EP14 permits
#'   raw or log10; EP30 permits all three scales; IFCC Part 2 permits raw or
#'   natural-log analysis; calibration effectiveness requires raw results.
#' @param reference_procedure Optional reference measurement procedure. When
#'   supplied, it is compared with every other procedure.
#' @param procedure_pairs Optional list of two-element procedure vectors.
#' @param conf.level Prediction-interval level for EP14 and EP30.
#' @param criterion Predefined IFCC Part 2 criterion on the selected analysis
#'   scale, or maximum IMPBR percentage for calibration effectiveness. It is
#'   mandatory for the two IFCC-derived approaches and is never estimated from
#'   data.
#' @param bias_model IFCC clinical-sample bias model: `"constant_global"`,
#'   `"constant_local"`, or `"linear_local"`.
#' @param local_n Number of clinical samples around each RM for a local IFCC
#'   model. Ignored for a global model.
#' @param coverage_factor IFCC uncertainty coverage factor; default 1.9.
#' @param rm_position IFCC position-mean uncertainty mode. `"pooled"` pools
#'   position-mean variances across reference materials within each procedure,
#'   as in the IFCC Part 2 worked example. `"material_specific"` uses each
#'   reference material's own position-mean variance.
#' @param rm_set_size Number of related reference materials used as a set for
#'   the EP30 normal critical-value adjustment. Independent materials use 1.
#' @param relative_residual_limit Positive EP30 relative-residual limit;
#'   default 2.
#' @param calibration_stage Column identifying `"before"` and `"after"`
#'   candidate-RM recalibration results; required for `approach = "calibration"`.
#' @param calibration_replicate_summary How multiple calibration results in the
#'   same sample, procedure, and stage are reduced to one reported result:
#'   `"mean"` (default), `"median"`, or `"error"` to reject duplicates.
#' @param exclude_procedures Optional procedures explicitly excluded after
#'   investigation in a calibration-effectiveness analysis. Full-set results
#'   are retained.
#' @param w3_reference_procedure Procedure whose SD/W(3) ratio scales the robust
#'   W(3) statistics. The first observed procedure is used when `NULL`.
#'
#' @return An object of class `commutability_result` containing material-level
#'   results, model estimates, design checks, IFCC position diagnostics, and
#'   the transformed analysis data. Design checks report whether the supplied
#'   study meets the sample-size and replication recommendations of the named
#'   method; unmet checks do not prevent otherwise computable analyses.
#' @examples
#' # Long-format replicate results for eight clinical samples.
#' clinical <- expand.grid(
#'   sample = paste0("CS", 1:8), procedure = c("MP1", "MP2"),
#'   replicate = 1:3, stringsAsFactors = FALSE
#' )
#' level <- setNames(seq(10, 80, length.out = 8), paste0("CS", 1:8))
#' clinical$result <- level[clinical$sample] +
#'   ifelse(clinical$procedure == "MP2", 1, 0) +
#'   rep(c(-0.2, 0, 0.2), each = 16)
#' clinical$material_type <- "clinical"
#' clinical$position <- NA_integer_
#'
#' # One reference material measured at three positions, in triplicate.
#' rm <- expand.grid(
#'   sample = "RM1", procedure = c("MP1", "MP2"),
#'   position = 1:3, replicate = 1:3, stringsAsFactors = FALSE
#' )
#' rm$result <- 45 + ifelse(rm$procedure == "MP2", 1, 0) +
#'   rep(c(-0.1, 0, 0.1), each = 6)
#' rm$material_type <- "rm"
#' columns <- c("sample", "material_type", "procedure", "result",
#'              "position", "replicate")
#' example_data <- rbind(clinical[columns], rm[columns])
#'
#' ep14 <- commutability(
#'   example_data, "sample", "material_type", "procedure", "result",
#'   approach = "ep14"
#' )
#' ep14$results
#'
#' # Explicitly map nonstandard material-type labels.
#' mapped_data <- example_data
#' mapped_data$material_type <- ifelse(
#'   mapped_data$material_type == "clinical", "patient_sample", "candidate_rm"
#' )
#' ep14_mapped <- commutability(
#'   mapped_data, "sample", "material_type", "procedure", "result",
#'   approach = "ep14",
#'   material_type_map = c(patient_sample = "clinical", candidate_rm = "rm")
#' )
#' ep14_mapped$results
#'
#' ifcc <- commutability(
#'   example_data, "sample", "material_type", "procedure", "result",
#'   position = "position", approach = "ifcc", criterion = 2,
#'   bias_model = "constant_global"
#' )
#' ifcc$results
#' @export
commutability <- function(data, sample_id, material_type = NULL,
                          procedure, result,
                          position = NULL,
                          approach = c("ep14", "ep30", "ifcc", "calibration"),
                          scale = c("raw", "log10", "ln"),
                          reference_procedure = NULL,
                          procedure_pairs = NULL,
                          conf.level = 0.95,
                          criterion = NULL,
                          bias_model = c("constant_global", "constant_local",
                                         "linear_local"),
                          local_n = NULL, coverage_factor = 1.9,
                          rm_position = c("pooled", "material_specific"),
                          rm_set_size = 1L,
                          relative_residual_limit = 2,
                          calibration_stage = NULL,
                          calibration_replicate_summary = c("mean", "median",
                                                            "error"),
                          exclude_procedures = NULL,
                          w3_reference_procedure = NULL,
                          material_type_map = NULL) {
  bias_model_supplied <- !missing(bias_model)
  local_n_supplied <- !missing(local_n)
  coverage_factor_supplied <- !missing(coverage_factor)
  rm_position_supplied <- !missing(rm_position)
  approach <- match.arg(approach)
  scale <- match.arg(scale)
  bias_model <- match.arg(bias_model)
  rm_position <- match.arg(rm_position)
  calibration_replicate_summary <- match.arg(calibration_replicate_summary)
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  id <- .comm_col(data, sample_id, "sample_id")
  type <- .comm_col(data, material_type, "material_type",
                    required = approach != "calibration")
  proc <- .comm_col(data, procedure, "procedure")
  value <- .comm_col(data, result, "result")
  pos <- .comm_col(data, position, "position", required = FALSE)
  stage <- .comm_col(data, calibration_stage, "calibration_stage",
                     required = approach == "calibration")
  if (is.null(type) && !is.null(material_type_map)) {
    stop("'material_type_map' requires a 'material_type' column.", call. = FALSE)
  }
  if (!is.numeric(value) || any(!is.finite(value)))
    stop("'result' must contain complete finite numeric values.", call. = FALSE)
  if (anyNA(id) || anyNA(proc))
    stop("Sample and procedure values must be complete.", call. = FALSE)
  if (is.null(type)) {
    type <- rep("clinical", nrow(data))
  } else {
    type <- .comm_material_type(type, material_type_map)
  }
  if (approach == "calibration") {
    if (any(type != "clinical")) {
      stop("Calibration analysis accepts clinical-sample results only.",
           call. = FALSE)
    }
  } else if (!all(c("clinical", "rm") %in% type)) {
    stop("Both clinical samples and RM materials are required.", call. = FALSE)
  }
  if (approach == "ep14" && scale == "ln")
    stop("EP14 supports 'raw' or 'log10' scale in this implementation.",
         call. = FALSE)
  if (approach == "ifcc" && scale == "log10")
    stop("IFCC supports 'raw' or 'ln' scale in this implementation.",
         call. = FALSE)
  if (approach == "calibration" && scale != "raw")
    stop("Calibration analysis requires the 'raw' scale.", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number in (0, 1).", call. = FALSE)
  if (approach == "ifcc" && is.null(pos))
    stop("'position' is required for IFCC analysis.", call. = FALSE)
  if (approach %in% c("ifcc", "calibration")) {
    if (!is.numeric(criterion) || length(criterion) != 1L ||
        !is.finite(criterion) || criterion <= 0) {
      stop("'criterion' must be a predefined positive number for this analysis.",
           call. = FALSE)
    }
  }
  if (approach == "ifcc") {
    if (!is.numeric(coverage_factor) || length(coverage_factor) != 1L ||
        !is.finite(coverage_factor) || coverage_factor <= 0)
      stop("'coverage_factor' must be a positive number.", call. = FALSE)
  }
  if (approach == "ep30") {
    if (!is.numeric(rm_set_size) || length(rm_set_size) != 1L ||
        !is.finite(rm_set_size) || rm_set_size < 1 ||
        rm_set_size != as.integer(rm_set_size)) {
      stop("'rm_set_size' must be a positive integer.", call. = FALSE)
    }
    rm_set_size <- as.integer(rm_set_size)
    if (!is.numeric(relative_residual_limit) ||
        length(relative_residual_limit) != 1L ||
        !is.finite(relative_residual_limit) ||
        relative_residual_limit <= 0) {
      stop("'relative_residual_limit' must be a positive number.",
           call. = FALSE)
    }
  }
  if (approach == "calibration") {
    if (!is.null(pos) || !is.null(reference_procedure) ||
        !is.null(procedure_pairs) || bias_model_supplied || local_n_supplied ||
        coverage_factor_supplied || rm_position_supplied) {
      stop("Pairwise-analysis arguments are not applicable to calibration analysis.",
           call. = FALSE)
    }
    if (anyNA(stage))
      stop("Calibration-stage values must be complete.", call. = FALSE)
    stage <- tolower(as.character(stage))
  }
  d <- data.frame(sample = as.character(id), type = type,
                  procedure = as.character(proc),
                  result = .comm_transform(value, scale),
                  position = if (is.null(pos)) NA_character_ else as.character(pos),
                  stage = if (is.null(stage)) NA_character_ else stage,
                  stringsAsFactors = FALSE)
  if (approach == "calibration") {
    analysis <- .comm_calibration(
      d, criterion = criterion,
      exclude_procedures = exclude_procedures,
      w3_reference_procedure = w3_reference_procedure,
      calibration_replicate_summary = calibration_replicate_summary
    )
    out <- list(
      call = match.call(), approach = approach, scale = scale,
      settings = list(
        criterion = criterion,
        material_type_map = material_type_map,
        calibration_stage = calibration_stage,
        calibration_replicate_summary = calibration_replicate_summary,
        exclude_procedures = exclude_procedures,
        w3_reference_procedure = analysis$w3_reference_procedure
      ),
      results = analysis$result, models = NULL,
      diagnostics = analysis$diagnostics,
      position_diagnostics = NULL,
      sample_results = analysis$sample_results,
      leave_one_out = analysis$leave_one_out,
      targets = analysis$targets,
      design_checks = analysis$design_checks,
      analysis_errors = NULL,
      dropped_samples = analysis$dropped_samples,
      clinical = analysis$sample_results,
      prediction_grid = NULL,
      data = analysis$data,
      input_data = d
    )
    class(out) <- "commutability_result"
    return(out)
  }
  procedures <- unique(d$procedure)
  pairs <- .comm_pairs(procedures, procedure_pairs, reference_procedure)
  summary <- .comm_summarise(d)
  analyse_pair <- function(pair) {
    tryCatch({
      if (approach == "ep14") {
        .comm_ep14_pair(summary, pair, conf.level)
      } else if (approach == "ep30") {
        .comm_ep30_pair(
          summary, pair, conf.level, rm_set_size, relative_residual_limit)
      } else {
        .comm_ifcc_pair(
          d, summary, pair, criterion, bias_model, local_n, coverage_factor,
          reference_procedure, rm_position)
      }
    }, error = function(e) {
      list(error = data.frame(
        pair = paste(pair, collapse = " vs "), message = conditionMessage(e),
        stringsAsFactors = FALSE
      ))
    })
  }
  attempted <- lapply(pairs, analyse_pair)
  failed <- vapply(attempted, function(z) !is.null(z$error), logical(1))
  analysis_errors <- .comm_bind(attempted[failed], "error")
  analyses <- attempted[!failed]
  if (!length(analyses)) {
    stop("No procedure pair could be analyzed: ",
         paste(analysis_errors$message, collapse = "; "), call. = FALSE)
  }
  out <- list(
    call = match.call(), approach = approach, scale = scale,
    settings = list(reference_procedure = reference_procedure,
      procedure_pairs = pairs, conf.level = conf.level, criterion = criterion,
      material_type_map = material_type_map,
      bias_model = if (approach == "ifcc") bias_model else NULL,
      local_n = local_n,
      coverage_factor = if (approach == "ifcc") coverage_factor else NULL,
      rm_position = if (approach == "ifcc") rm_position else NULL,
      rm_set_size = if (approach == "ep30") rm_set_size else NULL,
      relative_residual_limit = if (approach == "ep30")
        relative_residual_limit else NULL),
    results = do.call(rbind, lapply(analyses, `[[`, "result")),
    models = if (approach %in% c("ep14", "ep30"))
      do.call(rbind, lapply(analyses, `[[`, "model")) else NULL,
    diagnostics = if (approach == "ifcc")
      do.call(rbind, lapply(analyses, `[[`, "diagnostics")) else NULL,
    position_diagnostics = if (approach == "ifcc")
      do.call(rbind, lapply(analyses, `[[`, "position_diagnostics")) else NULL,
    design_checks = .comm_bind(analyses, "design_checks"),
    analysis_errors = analysis_errors,
    dropped_samples = NULL,
    clinical = do.call(rbind, lapply(analyses, `[[`, "clinical")),
    prediction_grid = if (approach == "ep30")
      do.call(rbind, lapply(analyses, `[[`, "prediction_grid")) else NULL,
    data = d, input_data = d
  )
  class(out) <- "commutability_result"
  out
}

#' Print a commutability assessment
#'
#' @param x A commutability_result object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged x, invisibly.
#' @export
print.commutability_result <- function(x, ...) {
  cat("\nCommutability assessment\n")
  if (!is.null(x$design_checks)) {
    unmet <- sum(!x$design_checks$met)
    cat("Study-design checks: ", nrow(x$design_checks) - unmet, " met, ",
        unmet, " unmet. Unmet checks do not invalidate computations; inspect ",
        "'$design_checks'.\n", sep = "")
  }
  if (!is.null(x$analysis_errors) && nrow(x$analysis_errors)) {
    cat("Procedure pairs not analyzed: ", nrow(x$analysis_errors),
        ". Inspect '$analysis_errors'.\n", sep = "")
  }
  if (x$approach == "calibration") {
    header <- data.frame(
      approach = "CALIBRATION (IFCC PART 3)", scale = x$scale,
      procedures = length(unique(x$results$procedure)),
      clinical_samples = nrow(x$targets),
      criterion_percent = x$settings$criterion,
      stringsAsFactors = FALSE
    )
    .print_kv_sections(header, list(
      "Assessment" = c(
        "Approach" = "approach", "Scale" = "scale",
        "Procedures" = "procedures",
        "Clinical samples" = "clinical_samples",
        "IMPBR criterion (%)" = "criterion_percent"
      )
    ))
    .print_kv_rows(x$results, list(
      "Procedure results" = c(
        "Clinical samples" = "n_clinical",
        "Median bias (%)" = "median_bias_percent",
        "Bias SD (%)" = "sd_bias_percent",
        "Adjusted W(3) (%)" = "adjusted_w3_percent",
        "Included after investigation" = "included"
      )
    ), row_title = c(
      "Stage" = "calibration_stage", "Procedure" = "procedure"
    ))
    .print_kv_rows(x$diagnostics, list(
      "IMPBR assessment" = c(
        "Procedures" = "n_procedures", "IMPBR (%)" = "impbr_percent",
        "Minimum procedure" = "min_procedure",
        "Maximum procedure" = "max_procedure",
        "Criterion (%)" = "criterion",
        "Excluded procedures" = "excluded_procedures",
        "Classification" = "classification"
      )
    ), row_title = c(
      "Stage" = "calibration_stage", "Assessment set" = "assessment_set"
    ))
    return(invisible(x))
  }
  header <- data.frame(
    approach = toupper(x$approach), scale = x$scale,
    procedure_pairs = length(unique(x$results$pair)),
    reference_materials = nrow(x$results), stringsAsFactors = FALSE
  )
  .print_kv_sections(header, list(
    "Assessment" = c(
      "Approach" = "approach", "Scale" = "scale",
      "Procedure pairs" = "procedure_pairs",
      "Reference materials" = "reference_materials"
    )
  ))
  if (x$approach == "ep14") {
    .print_kv_rows(x$results, list(
      "Results" = c(
        "Mean X" = "mean_x", "Mean Y" = "mean_y",
        "Predicted Y" = "predicted_y", "Residual" = "residual",
        "Prediction SE" = "prediction_se", "Lower" = "lower",
        "Upper" = "upper", "Extrapolated" = "extrapolated",
        "Classification" = "classification"
      )
    ), row_title = c(
      "Procedure pair" = "pair", "Reference material" = "rm"
    ))
  } else if (x$approach == "ep30") {
    .print_kv_rows(x$results, list(
      "Results" = c(
        "Mean X" = "mean_x", "Mean Y" = "mean_y",
        "Predicted Y" = "predicted_y", "Residual" = "residual",
        "Prediction SE" = "prediction_se", "Lower" = "lower",
        "Upper" = "upper", "Extrapolated" = "extrapolated",
        "Prediction classification" = "prediction_classification",
        "Sy.x" = "syx", "Relative residual" = "relative_residual",
        "Relative-residual limit" = "relative_residual_limit",
        "Relative-residual classification" =
          "relative_residual_classification"
      )
    ), row_title = c(
      "Procedure pair" = "pair", "Reference material" = "rm"
    ))
  } else {
    .print_kv_rows(x$results, list(
      "Results" = c(
        "RM coordinate" = "rm_coordinate", "RM bias" = "rm_bias",
        "Clinical bias" = "clinical_bias",
        "Difference in bias" = "bias_difference",
        "RM standard uncertainty" = "rm_standard_uncertainty",
        "Clinical standard uncertainty" = "clinical_standard_uncertainty",
        "Standard uncertainty" = "standard_uncertainty",
        "Expanded uncertainty" = "expanded_uncertainty",
        "Lower" = "lower", "Upper" = "upper",
        "Criterion" = "criterion", "Clinical samples" = "n_clinical_local",
        "Classification" = "classification"
      )
    ), row_title = c(
      "Procedure pair" = "pair", "Reference material" = "rm"
    ))
  }
  invisible(x)
}

#' Summarize a commutability assessment
#'
#' @param object A commutability_result object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged object, invisibly.
#' @export
summary.commutability_result <- function(object, ...) {
  print(object)
  if (!is.null(object$design_checks) && any(!object$design_checks$met)) {
    cat("\nUnmet study-design checks\n")
    print(object$design_checks[!object$design_checks$met,
      c("pair", "requirement", "observed", "recommended", "source")],
      row.names = FALSE)
  }
  if (!is.null(object$analysis_errors) && nrow(object$analysis_errors)) {
    cat("\nProcedure-pair analysis errors\n")
    print(object$analysis_errors, row.names = FALSE)
  }
  if (!is.null(object$models)) {
    if (object$approach == "ep30") {
      cat("\nEP30 model estimates\n")
      .print_kv_rows(object$models, list(
        "Model" = c(
          "Clinical samples" = "n_clinical",
          "Clinical replicate range" = "clinical_replicate_range",
          "Equal clinical replicates" = "clinical_replicates_equal",
          "RM replicate range" = "rm_replicate_range",
          "Equal RM replicates" = "rm_replicates_equal",
          "Variance X" = "variance_x", "Variance Y" = "variance_y",
          "Lambda" = "lambda", "Intercept" = "intercept",
          "Slope" = "slope",
          "Intercept variance" = "variance_intercept",
          "Slope variance" = "variance_slope",
          "Intercept-slope covariance" = "covariance_intercept_slope",
          "RM variance X" = "rm_variance_x",
          "RM variance Y" = "rm_variance_y",
          "Sy.x" = "syx", "Critical value" = "critical_value",
          "RM set size" = "rm_set_size"
        )
      ), row_title = c("Procedure pair" = "pair"))
    } else {
      cat("\nEP14 model estimates\n")
      .print_kv_rows(object$models, list(
        "Model" = c(
          "Clinical samples" = "n_clinical",
          "Variance X" = "variance_x", "Variance Y" = "variance_y",
          "Lambda" = "lambda", "Intercept" = "intercept",
          "Slope" = "slope", "Slope SE" = "slope_se", "DF" = "df"
        )
      ), row_title = c("Procedure pair" = "pair"))
    }
  }
  if (object$approach == "calibration") {
    cat("\nCalibration trend and robust-variation diagnostics\n")
    .print_kv_rows(object$results, list(
      "Diagnostics" = c(
        "W(3) (%)" = "w3_percent",
        "Adjusted W(3) (%)" = "adjusted_w3_percent",
        "Trend intercept" = "trend_intercept",
        "Trend slope" = "trend_slope",
        "Trend P value" = "trend_p_value"
      )
    ), row_title = c(
      "Stage" = "calibration_stage", "Procedure" = "procedure"
    ))
    cat("\nLeave-one-procedure-out sensitivity analysis\n")
    .print_kv_rows(object$leave_one_out, list(
      "Sensitivity" = c(
        "Remaining procedures" = "n_procedures",
        "IMPBR (%)" = "impbr_percent", "Criterion (%)" = "criterion",
        "Classification" = "classification"
      )
    ), row_title = c("Omitted procedure" = "omitted_procedure"))
  } else if (!is.null(object$diagnostics)) {
    cat("\nIFCC diagnostics\n")
    .print_kv_rows(object$diagnostics, list(
      "Diagnostics" = c(
        "Clinical samples" = "n_clinical", "MSSD SD" = "s_mssd",
        "Bias SD" = "s_b", "Variance X" = "variance_x",
        "Variance Y" = "variance_y",
        "Effective replicates" = "effective_replicates",
        "Sample-specific SD" = "sample_specific_sd",
        "Pooled position-mean SD X" = "pooled_position_mean_sd_x",
        "Pooled position-mean SD Y" = "pooled_position_mean_sd_y",
        "Position pooling" = "position_pooling"
      )
    ), row_title = c("Procedure pair" = "pair"))
  }
  invisible(object)
}

#' Plot a commutability assessment
#'
#' @param x A `commutability_result` object.
#' @param method For EP30, either `"prediction"` or `"relative_residual"`.
#'   Ignored by the other approaches.
#' @param ... Reserved arguments.
#' @return A ggplot object, invisibly.
#' @examples
#' clinical <- expand.grid(
#'   sample = paste0("CS", 1:8), procedure = c("MP1", "MP2"),
#'   replicate = 1:3, stringsAsFactors = FALSE
#' )
#' level <- setNames(seq(10, 80, length.out = 8), paste0("CS", 1:8))
#' clinical$result <- level[clinical$sample] +
#'   ifelse(clinical$procedure == "MP2", 1, 0) +
#'   rep(c(-0.2, 0, 0.2), each = 16)
#' clinical$material_type <- "clinical"
#' clinical$position <- NA_integer_
#' rm <- expand.grid(
#'   sample = "RM1", procedure = c("MP1", "MP2"),
#'   position = 1:3, replicate = 1:3, stringsAsFactors = FALSE
#' )
#' rm$result <- 45 + ifelse(rm$procedure == "MP2", 1, 0) +
#'   rep(c(-0.1, 0, 0.1), each = 6)
#' rm$material_type <- "rm"
#' columns <- c("sample", "material_type", "procedure", "result",
#'              "position", "replicate")
#' example_data <- rbind(clinical[columns], rm[columns])
#' assessment <- commutability(
#'   example_data, "sample", "material_type", "procedure", "result"
#' )
#' plot(assessment)
#' @export
plot.commutability_result <- function(x, method = NULL, ...) {
  if (!inherits(x, "commutability_result"))
    stop("'x' must be a commutability_result object.", call. = FALSE)
  .require_pkg("ggplot2")
  if (x$approach == "ep14") {
    lines <- do.call(rbind, lapply(seq_len(nrow(x$models)), function(i) {
      model <- x$models[i, ]
      cs <- x$clinical[x$clinical$pair == model$pair, ]
      grid <- seq(min(cs$mean_x), max(cs$mean_x), length.out = 200L)
      data.frame(pair = model$pair, x = grid,
                 fitted = model$intercept + model$slope * grid)
    }))
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(data = lines,
        ggplot2::aes(x = .data$x, y = .data$fitted), color = "steelblue") +
      ggplot2::geom_point(data = x$clinical,
        ggplot2::aes(x = .data$mean_x, y = .data$mean_y), color = "grey40") +
      ggplot2::geom_errorbar(data = x$results,
        ggplot2::aes(x = .data$mean_x, ymin = .data$lower,
                     ymax = .data$upper), width = 0) +
      ggplot2::geom_point(data = x$results,
        ggplot2::aes(x = .data$mean_x, y = .data$mean_y,
                     color = .data$classification), size = 3) +
      ggplot2::facet_wrap(~pair, scales = "free") +
      ggplot2::labs(x = "Procedure X", y = "Procedure Y",
                    title = "EP14 commutability", color = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else if (x$approach == "ep30") {
    if (is.null(method)) method <- "prediction"
    method <- match.arg(method, c("prediction", "relative_residual"))
    if (method == "relative_residual") {
      p <- ggplot2::ggplot(x$results,
        ggplot2::aes(x = .data$rm, y = .data$relative_residual,
                     color = .data$relative_residual_classification)) +
        ggplot2::geom_hline(
          yintercept = c(-x$settings$relative_residual_limit,
                          x$settings$relative_residual_limit),
          linetype = "dashed", color = "firebrick"
        ) +
        ggplot2::geom_point(size = 3) +
        ggplot2::facet_wrap(~pair, scales = "free_x") +
        ggplot2::labs(
          x = "Reference material", y = "Relative residual",
          title = "EP30 relative-residual assessment", color = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = x$prediction_grid,
          ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
          fill = "steelblue", alpha = 0.15
        ) +
        ggplot2::geom_line(
          data = x$prediction_grid,
          ggplot2::aes(x = .data$x, y = .data$fitted), color = "steelblue"
        ) +
        ggplot2::geom_point(
          data = x$clinical,
          ggplot2::aes(x = .data$mean_x, y = .data$mean_y), color = "grey40"
        ) +
        ggplot2::geom_point(
          data = x$results,
          ggplot2::aes(x = .data$mean_x, y = .data$mean_y,
                       color = .data$prediction_classification), size = 3
        ) +
        ggplot2::facet_wrap(~pair, scales = "free") +
        ggplot2::labs(
          x = "Procedure X", y = "Procedure Y",
          title = "EP30 commutability", color = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
  } else if (x$approach == "calibration") {
    p <- ggplot2::ggplot(
      x$sample_results,
      ggplot2::aes(x = .data$sample_order, y = .data$bias_percent,
                   color = .data$procedure, group = .data$procedure)
    ) +
      ggplot2::geom_hline(yintercept = 0, color = "grey50") +
      ggplot2::geom_line(alpha = 0.65) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::geom_hline(
        data = x$results,
        ggplot2::aes(yintercept = .data$median_bias_percent,
                     color = .data$procedure),
        linetype = "dashed"
      ) +
      ggplot2::facet_wrap(~stage, ncol = 1, scales = "free_y") +
      ggplot2::labs(
        x = "Clinical sample (ordered by initial target)",
        y = "Difference from initial target (%)",
        title = "Calibration effectiveness (IFCC Part 3)", color = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    p <- ggplot2::ggplot(x$results,
      ggplot2::aes(x = .data$rm, y = .data$bias_difference,
                   color = .data$classification)) +
      ggplot2::geom_hline(yintercept = c(-x$settings$criterion,
                                         x$settings$criterion),
                          linetype = "dashed", color = "firebrick") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                          ymax = .data$upper), width = 0.15) +
      ggplot2::geom_point(size = 3) +
      ggplot2::facet_wrap(~pair, scales = "free_x") +
      ggplot2::labs(x = "Reference material", y = "Difference in bias",
                    title = "IFCC commutability", color = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  print(p)
  invisible(p)
}
