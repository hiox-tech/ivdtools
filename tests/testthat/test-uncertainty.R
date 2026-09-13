test_that("Type A components distinguish a single result from a mean", {
  x <- c(8, 9, 10, 11, 12)
  single <- uncertainty_type_a(x, "repeatability", "single")
  average <- uncertainty_type_a(x, "repeatability mean", "mean")

  expect_s3_class(single, "uncertainty_component")
  expect_equal(single$estimate, mean(x))
  expect_equal(single$standard_uncertainty, sd(x))
  expect_equal(average$standard_uncertainty, sd(x) / sqrt(length(x)))
  expect_equal(single$df, length(x) - 1)

  relative <- uncertainty_type_a(x, "relative", "single", relative = TRUE)
  expect_equal(relative$standard_uncertainty, sd(x) / mean(x))
  expect_equal(relative$relative_uncertainty, sd(x) / mean(x))
})

test_that("Type B distributions are transformed to standard uncertainty", {
  rectangular <- uncertainty_type_b(
    "rectangular", 10, lower = 8, upper = 12,
    distribution = "rectangular"
  )
  triangular <- uncertainty_type_b(
    "triangular", 10, lower = 8, upper = 12,
    distribution = "triangular"
  )
  normal <- uncertainty_type_b(
    "normal", 10, lower = 8, upper = 12,
    distribution = "normal", coverage = 0.95
  )
  student <- uncertainty_type_b(
    "student", 10, lower = 8, upper = 12,
    distribution = "student_t", coverage = 0.95, df = 9
  )
  u_shaped <- uncertainty_type_b(
    "u", 10, lower = 8, upper = 12,
    distribution = "u_shaped"
  )
  expanded <- uncertainty_type_b(
    "certificate", 10, uncertainty = 0.8,
    distribution = "standard", k = 2
  )

  expect_equal(rectangular$standard_uncertainty, 2 / sqrt(3))
  expect_equal(triangular$standard_uncertainty, 2 / sqrt(6))
  expect_equal(normal$standard_uncertainty, 2 / qnorm(0.975))
  expect_equal(student$standard_uncertainty, 2 / qt(0.975, 9))
  expect_equal(u_shaped$standard_uncertainty, 2 / sqrt(2))
  expect_equal(expanded$standard_uncertainty, 0.4)
  expect_error(
    uncertainty_type_b("off center", 9, lower = 8, upper = 12),
    "midpoint"
  )
})

test_that("IQC one-way ANOVA reproduces the EP29 worked example", {
  d <- data.frame(
    run = rep(1:5, each = 5),
    result = c(140, 140, 140, 141, 140,
               138, 139, 138, 137, 139,
               143, 144, 144, 145, 143,
               143, 143, 142, 143, 142,
               142, 143, 141, 142, 143)
  )
  out <- uncertainty_iqc(d, "result", "run")
  expect_s3_class(out, "uncertainty_iqc")
  expect_equal(out$details$ms_within, 0.52, tolerance = 1e-12)
  expect_equal(out$details$ms_between, 24.4, tolerance = 1e-12)
  expect_equal(out$details$within_sd, sqrt(0.52), tolerance = 1e-12)
  expect_equal(out$details$between_sd, sqrt((24.4 - 0.52) / 5),
               tolerance = 1e-12)
  expect_equal(out$details$u_single,
               sqrt(0.52 + (24.4 - 0.52) / 5), tolerance = 1e-12)
  expect_equal(out$details$u_mean, sqrt(24.4 / 25), tolerance = 1e-12)

  long_term <- uncertainty_iqc(data.frame(result = 1:5), "result")
  expect_equal(long_term$details$u_single, sd(1:5))
  expect_equal(long_term$details$u_mean, sd(1:5) / sqrt(5))
})

test_that("homogeneity analysis reports observed and hidden components", {
  d <- data.frame(
    unit = rep(LETTERS[1:5], each = 3),
    result = c(10, 10.1, 9.9, 10.2, 10.1, 10.2,
               9.9, 9.8, 10, 10.1, 10, 10.2, 9.8, 9.9, 9.9)
  )
  out <- uncertainty_homogeneity(d, "unit", "result")
  z <- out$details
  groups <- split(d$result, d$unit)
  n_total <- sum(lengths(groups))
  df_within <- n_total - length(groups)
  ss_within <- sum(vapply(groups, function(values) {
    sum((values - mean(values))^2)
  }, numeric(1)))
  within_sd <- sqrt(ss_within / df_within)
  expected_hidden <- within_sd / sqrt(n_total / length(groups)) *
    (2 / df_within)^0.25
  expect_s3_class(out, "uncertainty_homogeneity")
  expect_equal(z$hidden_inhomogeneity, expected_hidden, tolerance = 1e-12)
  expect_equal(z$u_bb, max(z$between_unit_sd, z$hidden_inhomogeneity))
  expect_equal(out$component$standard_uncertainty, z$u_bb)
})

test_that("homogeneity and stability handle zero-centered absolute results", {
  homogeneity_data <- data.frame(
    unit = rep(LETTERS[1:4], each = 2),
    result = rep(c(-1, 1), 4)
  )
  homogeneity <- uncertainty_homogeneity(
    homogeneity_data, "unit", "result", relative = FALSE
  )
  expect_true(is.finite(homogeneity$details$u_bb))
  expect_true(all(is.na(homogeneity$details[
    c("relative_repeatability", "relative_between_unit",
      "relative_hidden_inhomogeneity", "relative_u_bb")
  ])))
  expect_error(
    uncertainty_homogeneity(homogeneity_data, "unit", "result", relative = TRUE),
    "grand mean must be nonzero"
  )

  stability_data <- data.frame(
    time = rep(0:2, each = 2),
    result = rep(c(-1, 1), 3)
  )
  stability <- uncertainty_stability(
    stability_data, "time", "result", shelf_life = 1, relative = FALSE
  )
  expect_true(is.finite(stability$details$u_lts))
  expect_true(is.na(stability$details$relative_u_lts))
})

test_that("stability and characterization return reusable components", {
  stability_data <- data.frame(
    time = rep(c(0, 3, 6, 9), each = 3),
    result = c(100, 101, 99, 100, 99, 100,
               98, 99, 99, 97, 98, 98)
  )
  stability <- uncertainty_stability(
    stability_data, "time", "result", shelf_life = 12
  )
  expect_equal(stability$details$u_lts,
               stability$details$slope_se * 12)
  expect_s3_class(stability$component, "uncertainty_component")

  characterization_data <- data.frame(
    lab = rep(LETTERS[1:4], each = 3),
    result = c(10, 10.1, 9.9, 10.2, 10.1, 10.3,
               9.8, 9.9, 9.7, 10.1, 10, 10.2)
  )
  characterization <- uncertainty_characterization(
    characterization_data, "lab", "result"
  )
  laboratory_means <- tapply(characterization_data$result,
                             characterization_data$lab, mean)
  expect_equal(characterization$assigned_value, mean(laboratory_means))
  expect_equal(characterization$standard_uncertainty,
               sd(laboratory_means) / sqrt(4))
})

test_that("analytical combination handles sensitivity and correlation", {
  a <- uncertainty_type_b("a", 10, uncertainty = 3,
                          distribution = "standard")
  b <- uncertainty_type_b("b", 20, uncertainty = 4,
                          distribution = "standard")
  independent <- uncertainty_combine(a, b, k = 2)
  expect_s3_class(independent, "uncertainty_budget")
  expect_equal(independent$combined_standard_uncertainty, 5)
  expect_equal(independent$expanded_uncertainty, 10)
  expect_equal(sum(independent$budget$contribution_percent), 100)

  correlated <- uncertainty_combine(
    a, b, correlation = matrix(c(1, 0.5, 0.5, 1), 2), k = 2
  )
  expect_equal(correlated$combined_standard_uncertainty, sqrt(37))
  expect_true(is.na(correlated$effective_df))

  transformed <- uncertainty_combine(a, b, sensitivity = c(2, -1), k = 2)
  expect_equal(transformed$combined_standard_uncertainty, sqrt(52))
  expect_error(
    uncertainty_combine(a, b, correlation = matrix(c(1, 2, 2, 1), 2)),
    "entries"
  )
})

test_that("relative components combine on relative and absolute scales", {
  a <- uncertainty_type_b("a", 100, uncertainty = 5,
                          distribution = "standard", relative = TRUE)
  b <- uncertainty_type_b("b", 100, uncertainty = 2,
                          distribution = "standard", relative = TRUE)
  relative <- uncertainty_combine(a, b, k = 2)
  expect_identical(relative$output_scale, "relative")
  expect_equal(relative$combined_standard_uncertainty,
               sqrt(0.05^2 + 0.02^2))

  absolute <- uncertainty_combine(a, b, value = 250, k = 2)
  expect_identical(absolute$output_scale, "absolute")
  expect_equal(absolute$combined_standard_uncertainty,
               250 * sqrt(0.05^2 + 0.02^2))

  c <- uncertainty_type_b("c", 100, uncertainty = 1,
                          distribution = "standard")
  expect_error(uncertainty_combine(a, c), "requires `value`")
})

test_that("delta propagation follows derivatives of the measurement model", {
  mass <- uncertainty_type_b("mass", 20, uncertainty = 0.2,
                             distribution = "standard")
  volume <- uncertainty_type_b("volume", 4, uncertainty = 0.1,
                               distribution = "standard")
  out <- uncertainty_propagate(
    function(mass, volume) mass / volume,
    list(mass, volume), method = "delta"
  )
  expected <- sqrt((1 / 4 * 0.2)^2 + (-20 / 4^2 * 0.1)^2)
  expect_equal(out$estimate, 5)
  expect_equal(unname(out$sensitivity), c(1 / 4, -20 / 4^2),
               tolerance = 1e-7)
  expect_equal(out$standard_uncertainty, expected, tolerance = 1e-8)
})

test_that("Monte Carlo propagation retains bounded component distributions", {
  rectangular <- uncertainty_type_b(
    "rectangular", 10, lower = 8, upper = 12,
    distribution = "rectangular"
  )
  triangular <- uncertainty_type_b(
    "triangular", 5, lower = 3, upper = 7,
    distribution = "triangular"
  )
  out <- uncertainty_propagate(
    function(rectangular, triangular) rectangular + triangular,
    list(rectangular, triangular), method = "monte_carlo",
    simulations = 20000, seed = 2026
  )
  expected <- sqrt((2 / sqrt(3))^2 + (2 / sqrt(6))^2)
  expect_s3_class(out, "uncertainty_propagation")
  expect_equal(out$simulated_mean, 15, tolerance = 0.04)
  expect_equal(out$standard_uncertainty, expected, tolerance = 0.04)
  expect_true(all(out$simulations >= 11 & out$simulations <= 19))
  expect_true(out$lower_uncertainty > 0)
  expect_true(out$upper_uncertainty > 0)
})

test_that("traceability chain accumulates newly introduced components", {
  reference <- uncertainty_type_b("reference", 1, uncertainty = 0.01,
                                  distribution = "standard", relative = TRUE)
  preparation <- uncertainty_type_b("preparation", 1, uncertainty = 0.02,
                                    distribution = "standard", relative = TRUE)
  transfer <- uncertainty_type_b("transfer", 1, uncertainty = 0.03,
                                 distribution = "standard", relative = TRUE)
  out <- uncertainty_traceability(
    list(reference = reference, preparation = preparation, transfer = transfer),
    k = 2
  )
  expect_s3_class(out, "uncertainty_traceability")
  expect_equal(out$summary$cumulative_standard_uncertainty,
               sqrt(cumsum(c(0.01, 0.02, 0.03)^2)))
  expect_equal(out$final$combined_standard_uncertainty, sqrt(0.0014))
})

test_that("uncertainty objects do not make acceptance decisions", {
  component <- uncertainty_type_b("specification", 10, lower = 9, upper = 11)
  budget <- uncertainty_combine(component)
  expect_false(any(grepl("accept|pass|fail", names(component), ignore.case = TRUE)))
  expect_false(any(grepl("accept|pass|fail", names(budget), ignore.case = TRUE)))
  expect_false(any(grepl("accept|pass|fail", names(budget$budget),
                         ignore.case = TRUE)))
})

test_that("uncertainty plot methods return ggplot objects", {
  a <- uncertainty_type_b("a", 10, uncertainty = 1,
                          distribution = "standard")
  b <- uncertainty_type_b("b", 20, lower = 18, upper = 22,
                          distribution = "rectangular")
  budget <- uncertainty_combine(a, b)
  propagated <- uncertainty_propagate(
    function(a, b) a + b, list(a, b), method = "monte_carlo",
    simulations = 2000, seed = 1
  )
  chain <- uncertainty_traceability(list(first = a, second = b))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_centered_plot(plot(budget))
  expect_centered_plot(plot(propagated, "distribution"))
  expect_centered_plot(plot(propagated, "contribution"))
  expect_centered_plot(plot(chain))
})

test_that("uncertainty methods validate dimensions and print invisibly", {
  a <- uncertainty_type_b("a", 10, uncertainty = 1,
                          distribution = "standard")
  b <- uncertainty_type_b("b", 20, uncertainty = 2,
                          distribution = "standard")
  budget <- uncertainty_combine(a, b)
  propagated <- uncertainty_propagate(
    function(a, b) a + b, list(a, b), method = "delta"
  )
  chain <- uncertainty_traceability(list(a = a, b = b))

  for (object in list(a, budget, propagated, chain)) {
    capture.output(printed <- withVisible(print(object)))
    expect_false(printed$visible)
    expect_identical(printed$value, object)
  }

  expect_error(uncertainty_combine(a, b, sensitivity = 1),
               "finite value per component")
  expect_error(uncertainty_combine(a, b, correlation = diag(3)),
               "matching the number")
  expect_error(uncertainty_propagate(function(a) a, list(a, b)),
               "unused argument")
  expect_error(uncertainty_propagate(function(a, b) a + b, list(a, b),
                                     method = "monte_carlo", simulations = 10),
               "at least")
})
