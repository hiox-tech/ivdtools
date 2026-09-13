# Internal validation helpers -------------------------------------------------

.ss_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop("'", name, "' must be a single finite value in (0, 1), got ",
         deparse(x), call. = FALSE)
  x
}

.ss_proportion <- function(x, name, endpoints = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (endpoints) {
    valid <- valid && x >= 0 && x <= 1
    range_text <- "[0, 1]"
  } else {
    valid <- valid && x > 0 && x < 1
    range_text <- "(0, 1)"
  }
  if (!valid)
    stop("'", name, "' must be a single finite value in ", range_text,
         ", got ", deparse(x), call. = FALSE)
  x
}

.ss_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop("'", name, "' must be a single finite positive value, got ",
         deparse(x), call. = FALSE)
  x
}

.ss_integer <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != floor(x))
    stop("'", name, "' must be a single integer >= ", minimum,
         ", got ", deparse(x), call. = FALSE)
  as.integer(x)
}


# Bland-Altman sample size ----------------------------------------------------

.ss_bland_altman_power <- function(n, mu, sd, delta, conf.level, agree.level) {
  z_gamma <- qnorm(1 - (1 - agree.level) / 2)
  df <- n - 1
  t_crit <- qt(1 - (1 - conf.level) / 2, df)
  se <- sd * sqrt(1 / n + z_gamma^2 / (2 * (n - 1)))
  lambda_1 <- (delta - mu - z_gamma * sd) / se
  lambda_2 <- (delta + mu - z_gamma * sd) / se

  value <- pt(t_crit, df, ncp = lambda_1, lower.tail = FALSE) +
    pt(-t_crit, df, ncp = lambda_2, lower.tail = TRUE)
  min(max(value, 0), 1)
}

#' Sample size for Bland-Altman agreement assessment
#'
#' Calculates the minimum number of measurement pairs required to attain a
#' target power using the noncentral-t approximation of Lu et al. (2016).
#'
#' @param power Target power in `(0, 1)`.
#' @param mu Anticipated mean difference between the two measurement methods.
#' @param sd Anticipated standard deviation of the differences.
#' @param delta Positive symmetric clinical agreement limit.
#' @param conf.level Confidence level for the confidence intervals around the
#'   limits of agreement. Default 0.95.
#' @param agree.level Central proportion defining the limits of agreement.
#'   Default 0.95.
#' @param method Calculation method. Currently only `"lu"` is available.
#' @param n.min Minimum sample size considered. Default 3.
#' @param n.max Maximum sample size considered. Default 100000.
#'
#' @return An object of class `sample_size_bland_altman`. The required
#'   sample size is available as `x$n`; the requested and attained powers
#'   are available as `x$target_power` and `x$achieved_power`.
#'
#' @references
#' Lu MJ, Zhong WH, Liu YX, Miao HZ, Li YC, Ji MH (2016).
#' Sample size for assessing agreement between two methods of measurement
#' by Bland-Altman method. *International Journal of Biostatistics*, 12(2).
#' \doi{10.1515/ijb-2015-0039}
#'
#' @examples
#' sample_size_bland_altman(
#'   power = 0.80, mu = 0.2, sd = 1, delta = 2.5
#' )
#'
#' @export
sample_size_bland_altman <- function(
    power = 0.80,
    mu,
    sd,
    delta,
    conf.level = 0.95,
    agree.level = 0.95,
    method = "lu",
    n.min = 3L,
    n.max = 1e5) {

  power <- .ss_probability(power, "power")
  .ss_positive(sd, "sd")
  .ss_positive(delta, "delta")
  .ss_probability(conf.level, "conf.level")
  .ss_probability(agree.level, "agree.level")
  n.min <- .ss_integer(n.min, "n.min", 3L)
  n.max <- .ss_integer(n.max, "n.max", n.min)

  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu))
    stop("'mu' must be a single finite numeric value, got ", deparse(mu),
         call. = FALSE)
  if (!identical(method, "lu"))
    stop("'method' must be \"lu\".", call. = FALSE)

  z_gamma <- qnorm(1 - (1 - agree.level) / 2)
  population_margin <- delta - abs(mu) - z_gamma * sd
  if (population_margin <= 0) {
    stop(sprintf(
      paste0("Agreement cannot be demonstrated under the assumed parameters: ",
             "delta (%.6g) must be greater than abs(mu) + z * sd (%.6g)."),
      delta, abs(mu) + z_gamma * sd), call. = FALSE)
  }

  power_at <- function(n) {
    .ss_bland_altman_power(n, mu, sd, delta, conf.level, agree.level)
  }

  if (power_at(n.max) < power) {
    stop(sprintf(
      "Target power %.4f is not reached at n.max = %d.", power, n.max),
      call. = FALSE)
  }

  lo <- n.min
  if (power_at(lo) < power) {
    hi <- n.max
    while (hi - lo > 1L) {
      mid <- (lo + hi) %/% 2L
      if (power_at(mid) >= power) hi <- mid else lo <- mid
    }
    n_result <- hi
  } else {
    n_result <- lo
  }

  structure(list(
    n = as.integer(n_result),
    target_power = power,
    achieved_power = power_at(n_result),
    mu = mu,
    sd = sd,
    delta = delta,
    conf.level = conf.level,
    agree.level = agree.level,
    method = method,
    population_margin = population_margin
  ), class = "sample_size_bland_altman")
}

#' Print Bland-Altman sample-size results
#'
#' @param x A `sample_size_bland_altman` result.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.sample_size_bland_altman <- function(x, ...) {
  cat("Sample size for Bland-Altman agreement assessment\n")
  cat("  Method            : Lu et al. (2016), noncentral-t approximation\n")
  cat(sprintf("  Mean difference   : %.4f\n", x$mu))
  cat(sprintf("  SD of differences : %.4f\n", x$sd))
  cat(sprintf("  Clinical limit    : %.4f\n", x$delta))
  cat(sprintf("  Confidence level  : %.4f\n", x$conf.level))
  cat(sprintf("  Agreement level   : %.4f\n", x$agree.level))
  cat(sprintf("  Sample size (n)   : %d\n", x$n))
  cat(sprintf("  Target power      : %.6f\n", x$target_power))
  cat(sprintf("  Achieved power    : %.6f\n", x$achieved_power))
  invisible(x)
}


# Proportion confidence-interval sample size --------------------------------

.ss_proportion_ci_bounds <- function(x, n, conf.level, interval, method) {
  alpha <- 1 - conf.level
  tail_alpha <- if (interval == "two.sided") alpha / 2 else alpha
  z <- qnorm(1 - tail_alpha)
  p_hat <- x / n

  bounds <- switch(method,
    wald = {
      se <- sqrt(p_hat * (1 - p_hat) / n)
      list(lower = p_hat - z * se, upper = p_hat + z * se)
    },
    `wald-cc` = {
      se <- sqrt(p_hat * (1 - p_hat) / n)
      list(lower = p_hat - z * se - 1 / (2 * n),
           upper = p_hat + z * se + 1 / (2 * n))
    },
    wilson = {
      denominator <- 1 + z^2 / n
      center <- (p_hat + z^2 / (2 * n)) / denominator
      radius <- z / denominator *
        sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
      list(lower = center - radius, upper = center + radius)
    },
    `wilson-cc` = {
      # Newcombe's continuity-corrected score interval.  The correction is
      # inside the score inversion; adding 1/(2n) to an uncorrected Wilson
      # interval is not algebraically equivalent.
      denominator <- 2 * (n + z^2)
      lower_radicand <- z^2 - (2 + 1 / n) +
        4 * p_hat * (n * (1 - p_hat) + 1)
      upper_radicand <- z^2 + (2 - 1 / n) +
        4 * p_hat * (n * (1 - p_hat) - 1)
      lower <- (2 * n * p_hat + z^2 - 1 -
        z * sqrt(pmax(lower_radicand, 0))) / denominator
      upper <- (2 * n * p_hat + z^2 + 1 +
        z * sqrt(pmax(upper_radicand, 0))) / denominator
      lower[lower_radicand < 0] <- 0
      upper[upper_radicand < 0] <- 1
      list(lower = lower, upper = upper)
    },
    `agresti-coull` = {
      n_tilde <- n + z^2
      p_tilde <- (x + z^2 / 2) / n_tilde
      radius <- z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
      list(lower = p_tilde - radius, upper = p_tilde + radius)
    },
    jeffreys = {
      lower <- qbeta(tail_alpha, x + 0.5, n - x + 0.5)
      upper <- qbeta(1 - tail_alpha, x + 0.5, n - x + 0.5)
      lower[x == 0] <- 0
      upper[x == n] <- 1
      list(lower = lower, upper = upper)
    },
    `clopper-pearson` = {
      lower <- numeric(length(x))
      upper <- numeric(length(x))
      lower[x == 0] <- 0
      lower[x > 0] <- qbeta(
        tail_alpha, x[x > 0], n - x[x > 0] + 1
      )
      upper[x == n] <- 1
      upper[x < n] <- qbeta(
        1 - tail_alpha, x[x < n] + 1, n - x[x < n]
      )
      list(lower = lower, upper = upper)
    }
  )

  bounds$lower <- pmax(0, bounds$lower)
  bounds$upper <- pmin(1, bounds$upper)
  if (interval == "lower") bounds$upper <- rep(1, length(x))
  if (interval == "upper") bounds$lower <- rep(0, length(x))
  bounds
}

.ss_proportion_ci_precision <- function(n, p, conf.level, interval, method) {
  # Sample-size planning evaluates the interval at the anticipated proportion.
  # Keeping n * p continuous avoids degenerate p-hat = 0/1 widths caused by
  # rounding at small n. Exact coverage is evaluated separately over integer x.
  x <- n * p
  p_hat <- p
  bounds <- .ss_proportion_ci_bounds(x, n, conf.level, interval, method)

  value <- switch(interval,
    two.sided = (bounds$upper - bounds$lower) / 2,
    lower = p_hat - bounds$lower,
    upper = bounds$upper - p_hat
  )
  list(value = value, x = x, p_hat = p_hat,
       lower = bounds$lower, upper = bounds$upper)
}

.ss_proportion_ci_coverage <- function(n, p, conf.level, interval, method) {
  x <- 0:n
  bounds <- .ss_proportion_ci_bounds(x, n, conf.level, interval, method)
  covered <- bounds$lower <= p & p <= bounds$upper
  sum(dbinom(x[covered], n, p))
}

#' Sample size for a single-proportion confidence interval
#'
#' Calculates the minimum sample size whose two-sided confidence interval has
#' total width no greater than twice `half_width`, or whose one-sided limit is
#' no farther than `half_width` from the anticipated sample proportion.
#'
#' @param p Anticipated proportion in `[0, 1]`. Boundary values are useful for
#'   one-sided exact intervals.
#' @param half_width For a two-sided interval, one half of the total interval
#'   width. For a one-sided interval, the maximum distance from the anticipated
#'   proportion to the confidence limit.
#' @param conf.level Confidence level. Default 0.95.
#' @param interval One of `"two.sided"`, `"lower"`, or `"upper"`.
#' @param method Confidence interval method. One of `"wald"`,
#'   `"wald-cc"`, `"wilson"`, `"wilson-cc"`, `"agresti-coull"`,
#'   `"jeffreys"`, or `"clopper-pearson"`.
#' @param n.max Maximum sample size considered. Default 1000000.
#'
#' @return An object of class `sample_size_proportion_ci`. The required
#'   sample size is available as `x$n`. The returned object also contains
#'   the achieved half-width and the exact binomial coverage at the assumed
#'   value of `p`.
#'
#' @examples
#' sample_size_proportion_ci(p = 0.3, half_width = 0.05)
#' sample_size_proportion_ci(
#'   p = 0.3, half_width = 0.03, method = "clopper-pearson"
#' )
#'
#' @export
sample_size_proportion_ci <- function(
    p,
    half_width,
    conf.level = 0.95,
    interval = c("two.sided", "lower", "upper"),
    method = c("wilson", "wald", "wald-cc", "agresti-coull",
               "jeffreys", "wilson-cc", "clopper-pearson"),
    n.max = 1e6) {

  interval <- match.arg(interval)
  method <- match.arg(method)
  .ss_proportion(p, "p", endpoints = TRUE)
  .ss_probability(conf.level, "conf.level")
  .ss_positive(half_width, "half_width")
  if (half_width > 1)
    stop("'half_width' must not exceed 1.", call. = FALSE)
  n.max <- .ss_integer(n.max, "n.max", 2L)

  n_result <- NULL
  precision <- NULL
  for (n in seq.int(2L, n.max)) {
    candidate <- .ss_proportion_ci_precision(
      n, p, conf.level, interval, method
    )
    if (candidate$value <= half_width) {
      n_result <- n
      precision <- candidate
      break
    }
  }
  if (is.null(n_result)) {
    stop(sprintf(
      "Target half-width %.6g is not reached at n.max = %d.",
      half_width, n.max), call. = FALSE)
  }

  actual_coverage <- .ss_proportion_ci_coverage(
    n_result, p, conf.level, interval, method
  )

  structure(list(
    n = as.integer(n_result),
    p = p,
    anticipated_successes = precision$x,
    anticipated_p_hat = precision$p_hat,
    interval = interval,
    method = method,
    conf.level = conf.level,
    target_half_width = half_width,
    achieved_half_width = precision$value,
    target_width = if (interval == "two.sided") 2 * half_width else NA_real_,
    achieved_width = if (interval == "two.sided")
      precision$upper - precision$lower else NA_real_,
    lower = precision$lower,
    upper = precision$upper,
    actual_coverage = actual_coverage,
    actual_alpha = 1 - actual_coverage
  ), class = "sample_size_proportion_ci")
}

#' Print single-proportion confidence-interval sample-size results
#'
#' @param x A `sample_size_proportion_ci` result.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.sample_size_proportion_ci <- function(x, ...) {
  cat("Sample size for a single-proportion confidence interval\n")
  cat(sprintf("  Method                : %s\n", x$method))
  cat(sprintf("  Interval              : %s\n", x$interval))
  cat(sprintf("  Anticipated proportion: %.4f\n", x$p))
  cat(sprintf("  Confidence level      : %.4f\n", x$conf.level))
  cat(sprintf("  Sample size (n)       : %d\n", x$n))
  cat(sprintf("  Target half-width     : %.6f\n", x$target_half_width))
  cat(sprintf("  Achieved half-width   : %.6f\n", x$achieved_half_width))
  if (x$interval == "two.sided") {
    cat(sprintf("  Target total width    : %.6f\n", x$target_width))
    cat(sprintf("  Achieved total width  : %.6f\n", x$achieved_width))
  }
  cat(sprintf("  Planning interval     : [%.6f, %.6f]\n", x$lower, x$upper))
  cat(sprintf("  Actual coverage       : %.6f\n", x$actual_coverage))
  cat(sprintf("  Actual alpha          : %.6f\n", x$actual_alpha))
  invisible(x)
}


# Single-proportion target-value test ----------------------------------------

.ss_test_lower_bound <- function(k, n, alpha, test) {
  p_hat <- k / n
  z <- qnorm(1 - alpha)
  switch(test,
    wald = p_hat - z * sqrt(p_hat * (1 - p_hat) / n),
    `wald-cc` = p_hat - z * sqrt(p_hat * (1 - p_hat) / n) - 1 / (2 * n),
    score = {
      denominator <- 1 + z^2 / n
      (p_hat + z^2 / (2 * n) -
         z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))) /
        denominator
    },
    `score-cc` = .ss_proportion_ci_bounds(
      k, n, 1 - alpha, "lower", "wilson-cc"
    )$lower,
    `agresti-coull` = {
      n_tilde <- n + z^2
      p_tilde <- (k + z^2 / 2) / n_tilde
      p_tilde - z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
    },
    jeffreys = if (k == 0L) 0 else qbeta(alpha, k + 0.5, n - k + 0.5),
    exact = if (k == 0L) 0 else qbeta(alpha, k, n - k + 1)
  )
}

.ss_test_upper_bound <- function(k, n, alpha, test) {
  p_hat <- k / n
  z <- qnorm(1 - alpha)
  switch(test,
    wald = p_hat + z * sqrt(p_hat * (1 - p_hat) / n),
    `wald-cc` = p_hat + z * sqrt(p_hat * (1 - p_hat) / n) + 1 / (2 * n),
    score = {
      denominator <- 1 + z^2 / n
      (p_hat + z^2 / (2 * n) +
         z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))) /
        denominator
    },
    `score-cc` = .ss_proportion_ci_bounds(
      k, n, 1 - alpha, "upper", "wilson-cc"
    )$upper,
    `agresti-coull` = {
      n_tilde <- n + z^2
      p_tilde <- (k + z^2 / 2) / n_tilde
      p_tilde + z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
    },
    jeffreys = if (k == n) 1 else qbeta(1 - alpha, k + 0.5, n - k + 0.5),
    exact = if (k == n) 1 else qbeta(1 - alpha, k + 1, n - k)
  )
}

.ss_test_critical <- function(n, p0, alpha, alternative, test) {
  find_upper <- function(tail_alpha) {
    if (.ss_test_lower_bound(n, n, tail_alpha, test) <= p0)
      return(n + 1L)
    lo <- 0L
    hi <- n
    while (lo < hi) {
      mid <- (lo + hi) %/% 2L
      if (.ss_test_lower_bound(mid, n, tail_alpha, test) > p0)
        hi <- mid
      else
        lo <- mid + 1L
    }
    lo
  }

  find_lower <- function(tail_alpha) {
    if (.ss_test_upper_bound(0L, n, tail_alpha, test) >= p0)
      return(-1L)
    lo <- 0L
    hi <- n
    while (lo < hi) {
      mid <- (lo + hi + 1L) %/% 2L
      if (.ss_test_upper_bound(mid, n, tail_alpha, test) < p0)
        lo <- mid
      else
        hi <- mid - 1L
    }
    lo
  }

  switch(alternative,
    greater = c(lower = NA_integer_, upper = find_upper(alpha)),
    less = c(lower = find_lower(alpha), upper = NA_integer_),
    two.sided = c(lower = find_lower(alpha / 2),
                  upper = find_upper(alpha / 2))
  )
}

.ss_proportion_oc <- function(n, p0, p1, alpha, alternative, test) {
  critical <- .ss_test_critical(n, p0, alpha, alternative, test)

  probabilities <- function(p) {
    lower <- if (!is.na(critical[["lower"]]) && critical[["lower"]] >= 0L)
      pbinom(critical[["lower"]], n, p) else 0
    upper <- if (!is.na(critical[["upper"]]) && critical[["upper"]] <= n)
      pbinom(critical[["upper"]] - 1L, n, p, lower.tail = FALSE) else 0
    min(lower + upper, 1)
  }

  list(
    critical = critical,
    actual_alpha = probabilities(p0),
    actual_power = probabilities(p1)
  )
}

.ss_proportion_normal_power <- function(
    n, p0, p1, alpha, alternative, test) {
  normal_tests <- c("exact", "score", "score-cc", "wald", "wald-cc")
  if (!test %in% normal_tests)
    stop("Normal power calculation is unavailable for this test.",
         call. = FALSE)

  z <- qnorm(if (alternative == "two.sided") 1 - alpha / 2 else 1 - alpha)
  variance_critical <- if (test %in% c("exact", "score", "score-cc"))
    p0 * (1 - p0) else p1 * (1 - p1)
  variance_alternative <- p1 * (1 - p1)
  denominator <- sqrt(variance_alternative)
  correction <- if (
      test %in% c("score-cc", "wald-cc") &&
      abs(p1 - p0) > 1 / (2 * n)
    ) 1 / (2 * sqrt(n)) else 0

  switch(alternative,
    greater = 1 - pnorm(
      (sqrt(n) * (p0 - p1) + z * sqrt(variance_critical) + correction) /
        denominator
    ),
    less = pnorm(
      (sqrt(n) * (p0 - p1) - z * sqrt(variance_critical) - correction) /
        denominator
    ),
    two.sided = {
      left <- pnorm(
        (sqrt(n) * (p0 - p1) - z * sqrt(variance_critical) - correction) /
          denominator
      )
      right <- 1 - pnorm(
        (sqrt(n) * (p0 - p1) + z * sqrt(variance_critical) + correction) /
          denominator
      )
      min(left + right, 1)
    }
  )
}

#' Sample size for a single-proportion target-value test
#'
#' Calculates the minimum sample size required to test a single proportion
#' against a target value. The test statistic and the power-calculation method
#' are selected independently. Regardless of the planning method, the returned
#' result includes the exact binomial Type I error and power of the resulting
#' discrete rejection region.
#'
#' @param p0 Null or target proportion in `(0, 1)`.
#' @param p1 Anticipated proportion under the alternative in `(0, 1)`.
#' @param alpha Target significance level in `(0, 1)`.
#' @param power Target power in `(0, 1)`.
#' @param alternative One of `"greater"`, `"less"`, or `"two.sided"`.
#' @param test Rejection-region method. One of `"exact"`, `"score"`,
#'   `"score-cc"`, `"wald"`, `"wald-cc"`, `"agresti-coull"`, or
#'   `"jeffreys"`.
#' @param power_method Power calculation used to select the sample size:
#'   `"binomial"` enumerates the discrete binomial distribution and
#'   requires both actual alpha and actual power to meet their targets;
#'   `"normal"` uses the continuous normal approximation and is available
#'   for exact, Score, and Wald tests, with or without continuity correction.
#' @param n.max Maximum sample size considered. Default 1000000.
#'
#' @return An object of class `sample_size_proportion`. The required
#'   sample size is available as `x$n`. The object also contains target
#'   and planning power, exact achieved power, actual alpha, and critical
#'   success counts.
#'
#' @details
#' `test = "score", power_method = "normal"` uses the null variance for
#' the alpha term and the alternative variance for the beta term. This is the
#' normal-approximation procedure used by the CMDE target-value formula.
#' With an infinite population, `test = "exact"` and `test = "score"` have
#' the same normal-approximation formula, matching PASS. The test choice still
#' determines the discrete rejection region used for the reported actual alpha
#' and achieved power.
#'
#' @examples
#' # Exact binomial design
#' sample_size_proportion(
#'   p0 = 0.2, p1 = 0.35, alpha = 0.05, power = 0.8,
#'   alternative = "greater"
#' )
#'
#' # CMDE continuous normal approximation
#' sample_size_proportion(
#'   p0 = 0.85, p1 = 0.90, alpha = 0.05, power = 0.8,
#'   alternative = "two.sided", test = "score", power_method = "normal"
#' )
#'
#' @export
sample_size_proportion <- function(
    p0,
    p1,
    alpha = 0.05,
    power = 0.80,
    alternative = c("greater", "less", "two.sided"),
    test = c("exact", "score", "score-cc", "wald", "wald-cc",
             "agresti-coull", "jeffreys"),
    power_method = c("binomial", "normal"),
    n.max = 1e6) {

  alternative <- match.arg(alternative)
  test <- match.arg(test)
  power_method <- match.arg(power_method)
  .ss_probability(p0, "p0")
  .ss_probability(p1, "p1")
  alpha <- .ss_probability(alpha, "alpha")
  power <- .ss_probability(power, "power")
  n.max <- .ss_integer(n.max, "n.max", 2L)

  if (alternative == "greater" && p1 <= p0)
    stop("For alternative = 'greater', p1 must be greater than p0.", call. = FALSE)
  if (alternative == "less" && p1 >= p0)
    stop("For alternative = 'less', p1 must be less than p0.", call. = FALSE)
  if (alternative == "two.sided" && p1 == p0)
    stop("For alternative = 'two.sided', p1 must differ from p0.", call. = FALSE)
  if (power_method == "normal" &&
      !test %in% c("exact", "score", "score-cc", "wald", "wald-cc"))
    stop("power_method = 'normal' is unavailable for this test.",
         call. = FALSE)

  n_result <- NULL
  planning_power <- NULL
  if (power_method == "normal") {
    for (n in seq.int(2L, n.max)) {
      candidate_power <- .ss_proportion_normal_power(
        n, p0, p1, alpha, alternative, test
      )
      if (candidate_power >= power) {
        n_result <- n
        planning_power <- candidate_power
        break
      }
    }
  } else {
    for (n in seq.int(2L, n.max)) {
      candidate <- .ss_proportion_oc(
        n, p0, p1, alpha, alternative, test
      )
      if (candidate$actual_power >= power &&
          candidate$actual_alpha <= alpha + 1e-12) {
        n_result <- n
        planning_power <- candidate$actual_power
        break
      }
    }
  }
  if (is.null(n_result)) {
    stop(sprintf("Target power %.4f is not reached at n.max = %d.",
                 power, n.max), call. = FALSE)
  }

  operating <- .ss_proportion_oc(
    n_result, p0, p1, alpha, alternative, test
  )
  critical <- operating$critical

  structure(list(
    n = as.integer(n_result),
    p0 = p0,
    p1 = p1,
    alternative = alternative,
    test = test,
    power_method = power_method,
    target_alpha = alpha,
    actual_alpha = operating$actual_alpha,
    target_power = power,
    planning_power = planning_power,
    achieved_power = operating$actual_power,
    alpha_met = operating$actual_alpha <= alpha + 1e-12,
    critical_lower = unname(critical[["lower"]]),
    critical_upper = unname(critical[["upper"]])
  ), class = "sample_size_proportion")
}

.ss_test_label <- function(test) {
  switch(test,
    exact = "Exact binomial test",
    score = "Score Z test using S(P0)",
    `score-cc` = "Score Z test using S(P0), continuity corrected",
    wald = "Wald Z test using S(Phat)",
    `wald-cc` = "Wald Z test using S(Phat), continuity corrected",
    `agresti-coull` = "Agresti-Coull interval inversion",
    jeffreys = "Jeffreys interval inversion"
  )
}

#' Print single-proportion target-value sample-size results
#'
#' @param x A `sample_size_proportion` result.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.sample_size_proportion <- function(x, ...) {
  relation_0 <- switch(x$alternative,
    greater = "<=", less = ">=", two.sided = "=")
  relation_1 <- switch(x$alternative,
    greater = ">", less = "<", two.sided = "!=")

  cat("Sample size for a single-proportion target-value test\n")
  cat(sprintf("  Test                  : %s\n", .ss_test_label(x$test)))
  cat(sprintf("  Power calculation     : %s\n",
              if (x$power_method == "normal") "Normal approximation"
              else "Binomial enumeration"))
  cat(sprintf("  H0: p %s %.4f\n", relation_0, x$p0))
  cat(sprintf("  H1: p %s %.4f\n", relation_1, x$p0))
  cat(sprintf("  Expected p            : %.4f\n", x$p1))
  cat(sprintf("  Sample size (n)       : %d\n", x$n))
  cat(sprintf("  Target alpha          : %.6f\n", x$target_alpha))
  cat(sprintf("  Actual alpha          : %.6f\n", x$actual_alpha))
  cat(sprintf("  Alpha criterion       : %s\n",
              if (x$alpha_met) "met" else "not met"))
  cat(sprintf("  Target power          : %.6f\n", x$target_power))
  if (x$power_method == "normal")
    cat(sprintf("  Planning power        : %.6f (normal approximation)\n",
                x$planning_power))
  cat(sprintf("  Achieved power        : %.6f (binomial evaluation)\n",
              x$achieved_power))

  if (x$alternative == "greater") {
    cat(sprintf("  Critical successes    : >= %d\n", x$critical_upper))
  } else if (x$alternative == "less") {
    cat(sprintf("  Critical successes    : <= %d\n", x$critical_lower))
  } else {
    cat(sprintf("  Critical successes    : <= %d or >= %d\n",
                x$critical_lower, x$critical_upper))
  }
  invisible(x)
}
