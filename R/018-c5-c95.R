# C5 and C95 estimation for qualitative examinations with an internal
# continuous response. The calculation uses a normal response distribution
# and a VFP precision profile for SD as a function of the expected response.

.c5c95_col <- function(data, column, argument) {
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

.c5c95_profile_from_data <- function(data, result, sample, model.no, K) {
  result <- .c5c95_col(data, result, "result")
  sample <- .c5c95_col(data, sample, "sample")
  keep <- !is.na(data[[sample]]) & !is.na(data[[result]])
  if (is.numeric(data[[result]])) keep <- keep & is.finite(data[[result]])
  excluded <- which(!keep)
  d <- data[keep, , drop = FALSE]
  if (!is.numeric(d[[result]])) {
    stop("The `result` column must be numeric.", call. = FALSE)
  }
  if (!nrow(d)) stop("No complete finite observations are available.", call. = FALSE)

  values <- split(as.numeric(d[[result]]), as.character(d[[sample]]), drop = TRUE)
  summary <- data.frame(
    sample = names(values),
    n = vapply(values, length, integer(1)),
    mean = vapply(values, mean, numeric(1)),
    SD = vapply(values, stats::sd, numeric(1)),
    stringsAsFactors = FALSE
  )
  summary$DF <- summary$n - 1
  summary$CV <- ifelse(abs(summary$mean) > .Machine$double.eps^0.5,
                       100 * summary$SD / abs(summary$mean), NA_real_)
  summary$component <- "repeatability"
  valid <- summary$n >= 2L & is.finite(summary$mean) &
    is.finite(summary$SD) & summary$SD > 0
  if (sum(valid) < 3L) {
    stop("At least three sample levels with two or more replicates and a ",
         "positive finite SD are required.", call. = FALSE)
  }
  omitted_levels <- summary$sample[!valid]
  fit <- NULL
  utils::capture.output(
    fit <- .sadler(summary[valid, , drop = FALSE], type = "sd",
                   model.no = model.no, K = K, quiet = TRUE)
  )
  if (!length(fit) || is.null(fit$vfp) || is.null(fit$best_no)) {
    stop("The VFP precision profile could not be fitted.", call. = FALSE)
  }
  list(fit = fit, summary = summary, excluded = excluded,
       omitted_levels = omitted_levels, result = result, sample = sample,
       source = "raw data", component = "repeatability")
}

.c5c95_profile_from_precision <- function(object, component) {
  if (is.null(object$profile) || !length(object$profile$fits)) {
    stop("The `precision` object has no fitted precision profile; run ",
         "variance() and profile() first.", call. = FALSE)
  }
  fits <- object$profile$fits
  if (is.null(component)) component <- if ("total" %in% names(fits)) "total" else names(fits)[1L]
  if (!is.character(component) || length(component) != 1L ||
      !component %in% names(fits)) {
    stop("`component` must be one of: ", paste(names(fits), collapse = ", "),
         ".", call. = FALSE)
  }
  fit <- fits[[component]]
  if (is.null(fit$vfp) || is.null(fit$best_no)) {
    stop("The selected precision-profile component is incomplete.", call. = FALSE)
  }
  list(fit = fit, summary = fit$data, excluded = integer(),
       omitted_levels = character(), result = object$response,
       sample = object$by, source = "precision object", component = component)
}

.c5c95_predict_sd <- function(fit, mean_response) {
  out <- tryCatch(
    suppressWarnings(VFP::predict.VFP(
      fit$vfp, model.no = fit$best_no,
      newdata = mean_response, type = "sd")$Fitted),
    error = function(e) rep(NA_real_, length(mean_response))
  )
  out <- as.numeric(out)
  out[!is.finite(out) | out <= 0] <- NA_real_
  out
}

.c5c95_probability <- function(mean_response, cutoff, direction, fit) {
  sd <- .c5c95_predict_sd(fit, mean_response)
  if (direction == "geq") {
    stats::pnorm(cutoff, mean = mean_response, sd = sd, lower.tail = FALSE)
  } else {
    stats::pnorm(cutoff, mean = mean_response, sd = sd, lower.tail = TRUE)
  }
}

.c5c95_find_root <- function(target, cutoff, direction, fit, data_range,
                             search_range, tol) {
  fixed_range <- !is.null(search_range)
  if (fixed_range) {
    bounds <- sort(as.numeric(search_range))
  } else {
    bounds <- range(c(data_range, cutoff), finite = TRUE)
    span <- diff(bounds)
    if (!is.finite(span) || span <= 0) {
      span <- max(1, abs(cutoff) * 0.1)
      bounds <- bounds + c(-1, 1) * span
    }
  }
  roots <- numeric()
  iterations <- if (fixed_range) 1L else 16L
  for (iteration in seq_len(iterations)) {
    grid <- seq(bounds[1L], bounds[2L], length.out = 1001L)
    value <- .c5c95_probability(grid, cutoff, direction, fit) - target
    ok <- is.finite(value)
    exact <- which(ok & abs(value) <= tol)
    if (length(exact)) roots <- c(roots, grid[exact])
    pairs <- which(ok[-length(ok)] & ok[-1L] &
                     value[-length(value)] * value[-1L] < 0)
    if (length(pairs)) {
      for (i in pairs) {
        root <- tryCatch(stats::uniroot(
          function(z) .c5c95_probability(z, cutoff, direction, fit) - target,
          interval = grid[c(i, i + 1L)], tol = tol)$root,
          error = function(e) NA_real_)
        roots <- c(roots, root)
      }
    }
    roots <- unique(roots[is.finite(roots)])
    if (length(roots) || fixed_range) break
    width <- diff(bounds)
    bounds <- bounds + c(-1, 1) * width
  }
  if (!length(roots)) {
    return(list(root = NA_real_, sd = NA_real_, actual = NA_real_,
                roots = 0L, search_lower = bounds[1L],
                search_upper = bounds[2L]))
  }
  root <- roots[which.min(abs(roots - cutoff))]
  sd <- .c5c95_predict_sd(fit, root)
  actual <- .c5c95_probability(root, cutoff, direction, fit)
  list(root = root, sd = sd, actual = actual, roots = length(roots),
       search_lower = bounds[1L], search_upper = bounds[2L])
}

#' Estimate C5 and C95 from an internal continuous response
#'
#' Estimates response values associated with selected probabilities of a
#' qualitative result. The response SD is modeled as a function of the expected
#' response using a VFP precision profile. For `direction = "geq"`, the event
#' is `result >= cutoff`; for `direction = "leq"`, the event is
#' `result <= cutoff`. No acceptance decision is made.
#'
#' @param data A data frame of replicate continuous responses, or a
#'   `precision` object containing a fitted precision profile.
#' @param result Name of the numeric response column. Required for a data frame.
#' @param sample Name of the sample or level column. Required for a data frame.
#' @param cutoff One or more finite decision cutoffs.
#' @param direction Direction defining the condition of interest: `"geq"`
#'   means response greater than or equal to the cutoff; `"leq"` means
#'   response less than or equal to the cutoff.
#' @param probabilities Probabilities to estimate, default `c(0.05, 0.95)`.
#' @param component Precision-profile component used when `data` is a
#'   `precision` object. Defaults to `"total"` when available.
#' @param model.no Candidate VFP model numbers for raw data.
#' @param K Fixed exponent supplied to applicable VFP models.
#' @param search_range Optional two-element finite numeric range for root
#'   finding. By default, the observed mean range is expanded automatically.
#' @param tol Positive numeric root-finding tolerance.
#'
#' @return An object of class `c5_c95`. Its `estimates` element contains
#'   the estimated response, fitted SD, achieved probability, extrapolation flag,
#'   and root-search diagnostics for each cutoff and probability.
#' @export
#' @examples
#' set.seed(1203)
#' cx_data <- data.frame(
#'   level = rep(paste0("L", 1:7), each = 20),
#'   response = unlist(lapply(seq(34, 46, length.out = 7), function(mu)
#'     stats::rnorm(20, mean = mu, sd = 0.7 + 0.01 * mu)))
#' )
#' cx <- c5_c95(cx_data, result = "response", sample = "level",
#'              cutoff = 40, model.no = 1)
#' cx$estimates
#' plot(cx)
c5_c95 <- function(data, result = NULL, sample = NULL, cutoff,
                   direction = c("geq", "leq"),
                   probabilities = c(0.05, 0.95), component = NULL,
                   model.no = 1:10, K = 2, search_range = NULL,
                   tol = 1e-7) {
  call <- match.call()
  direction <- match.arg(direction)
  if (!is.numeric(cutoff) || !length(cutoff) || anyNA(cutoff) ||
      any(!is.finite(cutoff))) {
    stop("`cutoff` must contain finite numeric values.", call. = FALSE)
  }
  if (!is.numeric(probabilities) || !length(probabilities) ||
      anyNA(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1)) {
    stop("`probabilities` must contain values strictly between 0 and 1.",
         call. = FALSE)
  }
  probabilities <- unique(as.numeric(probabilities))
  cutoff <- unique(as.numeric(cutoff))
  if (!is.numeric(model.no) || !length(model.no) || anyNA(model.no) ||
      any(!model.no %in% 1:10)) {
    stop("`model.no` must contain integers from 1 to 10.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1L || is.na(K) || !is.finite(K)) {
    stop("`K` must be one finite numeric value.", call. = FALSE)
  }
  if (!is.null(search_range) &&
      (!is.numeric(search_range) || length(search_range) != 2L ||
       anyNA(search_range) || any(!is.finite(search_range)) ||
       diff(range(search_range)) <= 0)) {
    stop("`search_range` must contain two distinct finite numeric values.",
         call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("`tol` must be one positive finite number.", call. = FALSE)
  }

  info <- if (inherits(data, "precision")) {
    .c5c95_profile_from_precision(data, component)
  } else if (is.data.frame(data)) {
    .c5c95_profile_from_data(data, result, sample, model.no, K)
  } else {
    stop("`data` must be a data frame or a `precision` object.", call. = FALSE)
  }
  fit <- info$fit
  data_range <- range(fit$data$mean, finite = TRUE)
  combinations <- expand.grid(cutoff = cutoff, probability = probabilities,
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  solved <- lapply(seq_len(nrow(combinations)), function(i) {
    z <- combinations[i, ]
    .c5c95_find_root(z$probability, z$cutoff, direction, fit,
                     data_range, search_range, tol)
  })
  estimates <- do.call(rbind, lapply(seq_along(solved), function(i) {
    z <- solved[[i]]
    p <- combinations$probability[i]
    data.frame(
      cutoff = combinations$cutoff[i], probability = p,
      label = paste0("C", format(100 * p, trim = TRUE, scientific = FALSE)),
      mean_response = z$root, SD = z$sd, actual_probability = z$actual,
      extrapolated = is.finite(z$root) &&
        (z$root < data_range[1L] || z$root > data_range[2L]),
      roots_found = z$roots, search_lower = z$search_lower,
      search_upper = z$search_upper, stringsAsFactors = FALSE
    )
  }))
  issues <- character()
  if (any(!is.finite(estimates$mean_response))) {
    issues <- c(issues, "One or more requested probabilities had no root in the search range.")
  }
  if (any(estimates$roots_found > 1L)) {
    issues <- c(issues, "Multiple roots were found; the root nearest the cutoff was selected.")
  }
  if (any(estimates$extrapolated, na.rm = TRUE)) {
    issues <- c(issues, "One or more estimates extrapolate beyond the observed mean-response range.")
  }
  if (length(info$omitted_levels)) {
    issues <- c(issues, paste0("Levels omitted from profile fitting: ",
                               paste(info$omitted_levels, collapse = ", "), "."))
  }

  curve_range <- range(c(data_range, cutoff, estimates$mean_response), finite = TRUE)
  span <- diff(curve_range)
  if (!is.finite(span) || span <= 0) span <- max(1, abs(curve_range[1L]) * 0.1)
  curve_x <- seq(curve_range[1L] - 0.05 * span,
                 curve_range[2L] + 0.05 * span, length.out = 401L)
  probability_curve <- do.call(rbind, lapply(cutoff, function(cut) {
    data.frame(mean_response = curve_x, cutoff = cut,
      SD = .c5c95_predict_sd(fit, curve_x),
      probability = .c5c95_probability(curve_x, cut, direction, fit),
      stringsAsFactors = FALSE)
  }))
  out <- list(call = call, estimates = estimates,
              profile_data = info$summary,
              probability_curve = probability_curve,
              profile = fit, source = info$source,
              component = info$component, direction = direction,
              data_range = data_range, excluded_rows = info$excluded,
              issues = unique(issues))
  class(out) <- "c5_c95"
  out
}

#' Print C5 and C95 estimates
#'
#' @param x A `c5_c95` object.
#' @param digits Number of digits shown.
#' @param ... Reserved arguments.
#' @return Invisibly returns `x`.
#' @export
#' @examples
#' set.seed(12)
#' d <- data.frame(level = rep(1:5, each = 10),
#'                 y = stats::rnorm(50, rep(seq(37, 43, length.out = 5), each = 10), 1))
#' print(c5_c95(d, "y", "level", cutoff = 40, model.no = 1))
print.c5_c95 <- function(x, digits = 4L, ...) {
  cat("C5/C95 estimation from an internal continuous response\n")
  cat(sprintf("  Source: %s\n", x$source))
  cat(sprintf("  Precision component: %s\n", x$component))
  cat(sprintf("  Event: response %s cutoff\n",
              if (x$direction == "geq") ">=" else "<="))
  shown <- x$estimates[, c("cutoff", "label", "probability", "mean_response",
                           "SD", "actual_probability", "extrapolated"), drop = FALSE]
  numeric <- vapply(shown, is.numeric, logical(1))
  shown[numeric] <- lapply(shown[numeric], round, digits = digits)
  .print_df(shown)
  if (length(x$excluded_rows)) cat("  Excluded incomplete rows:", length(x$excluded_rows), "\n")
  if (length(x$issues)) {
    cat("Notes\n")
    for (note in x$issues) cat("  - ", note, "\n", sep = "")
  }
  invisible(x)
}

#' Plot C5 and C95 results
#'
#' @param x A `c5_c95` object.
#' @param type Plot type: `"probability"` or `"profile"`.
#' @param ... Reserved arguments.
#' @return A ggplot object.
#' @export
#' @examples
#' set.seed(95)
#' d <- data.frame(level = rep(1:6, each = 12),
#'                 y = stats::rnorm(72, rep(seq(36, 44, length.out = 6), each = 12), 1))
#' z <- c5_c95(d, "y", "level", cutoff = 40, model.no = 1)
#' plot(z, type = "probability")
plot.c5_c95 <- function(x, type = c("probability", "profile"), ...) {
  type <- match.arg(type)
  estimates <- x$estimates[is.finite(x$estimates$mean_response), , drop = FALSE]
  if (type == "probability") {
    p <- ggplot2::ggplot(x$probability_curve,
      ggplot2::aes(x = .data$mean_response, y = .data$probability)) +
      ggplot2::geom_line(color = "#2166AC", linewidth = 0.9, na.rm = TRUE) +
      ggplot2::geom_hline(data = estimates,
        ggplot2::aes(yintercept = .data$probability), color = "grey55",
        linetype = "dotted") +
      ggplot2::geom_vline(data = estimates,
        ggplot2::aes(xintercept = .data$mean_response, color = .data$label),
        linewidth = 0.7, linetype = "dashed") +
      ggplot2::geom_vline(data = unique(x$probability_curve["cutoff"]),
        ggplot2::aes(xintercept = .data$cutoff), color = "black", linewidth = 0.6) +
      ggplot2::scale_y_continuous(limits = c(0, 1),
        breaks = c(0, 0.05, 0.5, 0.95, 1)) +
      ggplot2::labs(x = "Expected internal continuous response",
        y = "Probability of condition of interest", color = NULL,
        title = "C5 and C95 probability profile") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    if (length(unique(x$probability_curve$cutoff)) > 1L) {
      p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$cutoff),
                                   scales = "free_x")
    }
    return(p)
  }

  observed <- x$profile$data
  grid <- seq(min(c(x$data_range, estimates$mean_response), na.rm = TRUE),
              max(c(x$data_range, estimates$mean_response), na.rm = TRUE),
              length.out = 401L)
  fitted <- data.frame(mean = grid, SD = .c5c95_predict_sd(x$profile, grid))
  ggplot2::ggplot(observed, ggplot2::aes(x = .data$mean, y = .data$SD)) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_line(data = fitted, color = "#2166AC", linewidth = 0.9) +
    ggplot2::geom_vline(data = estimates,
      ggplot2::aes(xintercept = .data$mean_response, color = .data$label),
      linewidth = 0.7, linetype = "dashed") +
    ggplot2::labs(x = "Expected internal continuous response", y = "SD",
      color = NULL, title = "Precision profile and C5/C95") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}
