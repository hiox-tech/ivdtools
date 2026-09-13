linearity_test_data <- function(replicates = 3L) {
  x <- rep(seq(10, 60, by = 10), each = replicates)
  level <- rep(paste0("L", 1:6), each = replicates)
  residual <- rep(seq(-0.3, 0.3, length.out = replicates), 6)
  data.frame(level = level, expected = x, result = 2 + 1.5 * x + residual)
}

test_that("linearity study-design helpers follow their defining equations", {
  z <- qnorm(0.995)
  counts <- linearity_replicates(4, 5, level_probability = 0.99,
                                 min_replicates = 2)
  expect_s3_class(counts, "linearity_replicates")
  expect_equal(counts$calculated_replicates, (z * 4 / 5)^2)
  expect_equal(counts$replicates, max(2, ceiling((z * 4 / 5)^2)))

  family <- linearity_replicates(4, 5, levels = 6, family_error = 0.1)
  expect_equal(family$level_probability, 0.9^(1 / 6))

  panel <- linearity_panel(c(1, 0.5, 0), high = 110, low = 10,
                           total_volume = 2, sample = c("H", "M", "L"))
  expect_s3_class(panel, "linearity_panel")
  expect_equal(panel$expected, c(110, 60, 10))
  expect_equal(panel$high_volume + panel$low_volume, rep(2, 3))

  endpoints <- linearity_endpoints(
    lloq = 10, uloq = 1000, cv_high = 5,
    cv_at_lloq = 20, cv_at_k_lloq = 15
  )
  one_sided_z <- qnorm(0.95)
  expect_equal(endpoints$high_target, 1000 / (1 + one_sided_z * 0.05))
  expect_true(is.finite(endpoints$low_target))
  expect_gt(endpoints$low_target, 10)
})

test_that("straight-line analysis agrees with an independent level-mean fit", {
  d <- linearity_test_data()
  fit <- suppressWarnings(linearity(
    d, "level", "result", "expected", method = "linear",
    weights = "equal", conf.level = 0.95,
    multiplicity = "none", adl_abs = 1
  ))
  expected <- lm(mean ~ x, data = fit$level_summary)

  expect_s3_class(fit, "linearity")
  expect_equal(unname(coef(fit$linear_fit)), unname(coef(expected)))
  expect_equal(fit$level_summary$residual_linear, rep(0, 6), tolerance = 1e-10)
  expect_equal(fit$settings$point_confidence, 0.95)
  expect_false(any(fit$level_summary$adl_exceeded))
  expect_true(all(c("ci_lower", "ci_upper") %in% names(fit$level_summary)))
})

test_that("linearity supports variance, pooled, profile, and custom weighting", {
  d <- linearity_test_data()
  d$result <- d$result + rep(c(-0.1, 0, 0.1), 6) * rep(1:6, each = 3)

  level <- suppressWarnings(linearity(d, "level", "result", "expected",
                                      weights = "level_variance"))
  pooled_data <- data.frame(
    level = rep(paste0("P", 1:9), each = 2),
    expected = rep(seq(10, 90, by = 10), each = 2),
    result = 2 + 1.5 * rep(seq(10, 90, by = 10), each = 2) +
      rep(c(-0.2, 0.2), 9),
    variance_group = rep(rep(c("low", "middle", "high"), each = 3), each = 2)
  )
  pooled <- suppressWarnings(linearity(
    pooled_data, "level", "result", "expected", weights = "pooled",
    variance_group = "variance_group"
  ))
  profile <- suppressWarnings(linearity(
    d, "level", "result", "expected", weights = "precision_profile"
  ))
  custom <- suppressWarnings(linearity(
    d, "level", "result", "expected", weights = seq_len(6)
  ))

  expect_equal(level$settings$weights, "level_variance")
  expect_true(all(level$level_summary$variance_used > 0))
  expect_equal(nrow(pooled$precision$pooled), 3L)
  expect_equal(pooled$level_summary$variance_group,
               rep(c("low", "middle", "high"), each = 3))
  expect_equal(pooled$settings$variance_group, "variance_group")
  expect_true(all(pooled$level_summary$fit_weight > 0))
  expect_s3_class(profile$precision$profile$fit, "lm")
  expect_equal(custom$settings$fit_unit, "level_mean")
})

test_that("explicit pooled groups reproduce EP06 Appendix G", {
  result <- c(
    525.0, 533.0, 483.3, 510.0, 453.8, 482.6,
    383.5, 399.4, 346.7, 353.6, 287.6, 296.5, 229.2, 237.1,
    171.7, 172.2, 110.0, 112.6, 51.0, 51.5
  )
  d <- data.frame(
    level = rep(1:10, each = 2),
    expected = rep(529 * seq(1, 0.1, by = -0.1), each = 2),
    result = result,
    variance_group = rep(c(rep("high", 3), rep("middle", 4), rep("low", 3)),
                         each = 2)
  )
  fit <- linearity(
    d, "level", "result", "expected", method = "linear",
    intercept = FALSE, weights = "pooled",
    variance_group = "variance_group"
  )

  pooled <- fit$precision$pooled
  pooled <- pooled[match(c("high", "middle", "low"),
                         pooled$variance_group), ]
  expect_equal(pooled$pooled_variance,
               c(267.721666666667, 55.255, 1.21), tolerance = 1e-10)
  expect_equal(pooled$fit_weight,
               2 / c(267.721666666667, 55.255, 1.21), tolerance = 1e-10)
  expect_equal(fit$coefficients$estimate, 1.0677, tolerance = 5e-5)
  expect_equal(max(abs(fit$level_summary$percent_residual_linear)),
               9.3, tolerance = 0.05)

  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "Variance group: variance_group")
  expect_match(printed, "Pooled variance groups")
  expect_centered_plot(plot(fit, type = "precision"))
})

test_that("pooled weighting rejects absent or invalid variance groups", {
  d <- data.frame(
    level = rep(1:6, each = 2),
    expected = rep(1:6, each = 2),
    result = rep(1:6, each = 2) + rep(c(-0.1, 0.1), 6),
    variance_group = rep(rep(c("a", "b", "c"), each = 2), each = 2)
  )
  expect_error(
    linearity(d, "level", "result", "expected", weights = "pooled"),
    "must name a data column"
  )
  expect_error(
    linearity(d, "level", "result", "expected", weights = "equal",
              variance_group = "variance_group"),
    "used only"
  )

  inconsistent <- d
  inconsistent$variance_group[2] <- "b"
  expect_error(
    linearity(inconsistent, "level", "result", "expected",
              weights = "pooled", variance_group = "variance_group"),
    "constant within each sample"
  )

  noncontiguous <- d
  noncontiguous$variance_group <- rep(c("a", "b", "a", "c", "c", "c"),
                                      each = 2)
  expect_error(
    linearity(noncontiguous, "level", "result", "expected",
              weights = "pooled", variance_group = "variance_group"),
    "contiguous block"
  )

  singleton <- d
  singleton$variance_group <- rep(c("a", "b", "b", "c", "c", "c"), each = 2)
  expect_error(
    linearity(singleton, "level", "result", "expected",
              weights = "pooled", variance_group = "variance_group"),
    "at least two sample levels"
  )
})

test_that("polynomial analysis detects curvature and exposes all plot types", {
  d <- linearity_test_data(4)
  d$result <- 2 + 0.9 * d$expected + 0.02 * d$expected^2 +
    rep(c(-0.15, -0.05, 0.05, 0.15), 6)
  fit <- linearity(d, "level", "result", "expected",
                   method = "polynomial", max_degree = 3, adl_rel = 5)

  expect_equal(fit$selected_degree, 2)
  expect_true(any(fit$coefficients$nonlinear & fit$coefficients$significant))
  expect_true(any(grepl("exploratory", fit$warnings)))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (type in c("fit", "replicates", "precision", "residual",
                 "difference")) {
    expect_centered_plot(plot(fit, type = type))
  }
})

test_that("linearity rejects ambiguous and inadequate study data", {
  d <- linearity_test_data()
  expect_error(linearity(transform(d, result = replace(result, 1, NA)),
                         "level", "result", "expected"), "does not silently")
  expect_error(linearity(transform(d, expected = replace(expected, 1, 11)),
                         "level", "result", "expected"), "constant within")
  expect_error(linearity(d[d$level %in% c("L1", "L2", "L3", "L4"), ],
                         "level", "result", "expected", method = "polynomial"),
               "at least five")
  expect_error(linearity(d, "level", "result", "expected",
                         weights = c(1, 2)), "one value per observation")
  expect_error(linearity_panel(c(0, 1.2), 100), "between 0 and 1")
  expect_error(linearity_endpoints(uloq = 100), "both `uloq` and `cv_high`")
  expect_error(linearity_replicates(5, 4, true_deviation = 4),
               "must be positive")
})
