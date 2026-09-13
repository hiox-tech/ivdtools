test_that("nonparametric reference interval uses EP28 ranks and exact rank CIs", {
  values <- seq_len(120)
  ri <- reference_interval(data.frame(value = values), "value", ci = 0.90)

  expect_identical(ri$methods, "nonparametric")
  expect_equal(ri$nonparametric$lower_rank, 3.025, tolerance = 1e-12)
  expect_equal(ri$nonparametric$upper_rank, 117.975, tolerance = 1e-12)
  expect_equal(ri$nonparametric$lower, 3.025, tolerance = 1e-12)
  expect_equal(ri$nonparametric$upper, 117.975, tolerance = 1e-12)
  expect_equal(ri$nonparametric$ci_lower$rank, c(1L, 7L))
  expect_equal(ri$nonparametric$ci_upper$rank, c(114L, 120L))
  expect_equal(ri$nonparametric$ci_lower$actual_coverage,
               0.920466554592726, tolerance = 1e-12)
  expect_true(ri$nonparametric$ci_attainable)

  values_240 <- seq_len(240)
  ri_240 <- reference_interval(data.frame(value = values_240), "value",
                               ci = 0.90)
  expect_equal(c(ri_240$nonparametric$lower_rank,
                 ri_240$nonparametric$upper_rank), c(6.025, 234.975),
               tolerance = 1e-12)
  expect_equal(ri_240$nonparametric$ci_lower$rank, c(2L, 11L))
  expect_equal(ri_240$nonparametric$ci_upper$rank, c(230L, 239L))
})

test_that("type6 retains the prior quantile and normal-rank procedure", {
  values <- seq_len(120)
  ri <- reference_interval(data.frame(value = values), "value",
                           method = "type6")
  expect_equal(ri$type6$lower,
               unname(stats::quantile(values, 0.025, type = 6)))
  expect_equal(ri$type6$upper,
               unname(stats::quantile(values, 0.975, type = 6)))
  expect_equal(ri$type6$ci_lower$rank, c(1, 7))
  expect_equal(ri$type6$ci_upper$rank, c(113, 120))
})

test_that("parametric reference interval agrees with mean plus or minus z SD", {
  values <- c(-2, -1, 0, 1, 2, -1.5, 1.5, -0.5, 0.5)
  ri <- suppressWarnings(reference_interval(
    data.frame(value = values), "value", method = "parametric"
  ))
  z <- stats::qnorm(0.975)
  expected <- mean(values) + c(-1, 1) * z * stats::sd(values)
  expect_num_equal(c(ri$parametric$lower, ri$parametric$upper), expected, 1e-12)
})

test_that("biweight reproduces EP28 Appendix B", {
  values <- c(
    8.9, 9.2, 9.4, 9.4, 9.5, 9.5, 9.5, 9.6, 9.6, 9.6,
    9.6, 9.7, 9.7, 9.7, 9.7, 9.7, 9.8, 9.9, 9.9, 10.2
  )
  ri <- suppressWarnings(reference_interval(
    data.frame(value = values), "value", method = "biweight", ci = 0.90,
    bootstrap_replicates = 199, seed = 28
  ))

  expect_equal(ri$biweight$center, 9.6244, tolerance = 5e-5)
  expect_equal(ri$biweight$s_bi_205_6, 0.27043, tolerance = 5e-5)
  expect_lt(abs(ri$biweight$s_t_3_7 - 0.04816), 5e-5)
  expect_equal(ri$biweight$lower, 9.05, tolerance = 5e-3)
  expect_equal(ri$biweight$upper, 10.20, tolerance = 5e-3)
  expect_true(ri$biweight$converged)
  expect_equal(ri$biweight$iterations, 6L)
  expect_equal(ri$biweight$bootstrap_valid, 199L)
})

test_that("huber retains the prior robust procedure under its true name", {
  set.seed(20260830)
  values <- c(seq(9, 11, length.out = 119), 50)
  ri <- suppressWarnings(reference_interval(
    data.frame(value = values), "value", method = "huber",
    bootstrap_replicates = 199, seed = 17
  ))

  expect_s3_class(ri, "reference_interval")
  expect_true(all(is.finite(c(ri$huber$lower, ri$huber$upper))))
  expect_true(ri$huber$lower <= ri$huber$upper)
  expect_match(ri$huber$algorithm, "Huber M-estimation")
})

test_that("reference interval method names are explicit and old names fail", {
  d <- data.frame(value = seq_len(120))
  expect_error(reference_interval(d, "value", method = "percentile"),
               "Invalid method")
  expect_error(reference_interval(d, "value", method = "robust"),
               "Invalid method")
  expect_error(reference_interval(d, "value", method = c("all", "type6")),
               "alone")
  expect_error(reference_interval(d, "value", method = "huber",
                                  bootstrap_replicates = 99),
               "at least 100")
  expect_error(reference_interval(d, "value", method = "biweight",
                                  seed = 1.5),
               "finite integer")
})

test_that("nonparametric method reports unattainable requested coverage", {
  d <- data.frame(value = seq_len(120))
  expect_warning(
    ri <- reference_interval(d, "value", method = "nonparametric", ci = 0.95),
    "not fully attainable"
  )
  expect_false(ri$nonparametric$ci_attainable)
  expect_lt(ri$nonparametric$ci_lower$actual_coverage, 0.95)
  expect_equal(ri$nonparametric$ci_lower$rank[1L], 1L)
})

test_that("biweight fails explicitly when the MAD scale is zero", {
  d <- data.frame(value = c(rep(10, 19), 11))
  expect_error(
    reference_interval(d, "value", method = "biweight",
                       bootstrap_replicates = 100),
    "positive MAD"
  )
})

test_that("all reference interval methods print and plot", {
  set.seed(1701)
  d <- data.frame(value = stats::rnorm(120, 100, 8))
  ri <- suppressWarnings(reference_interval(
    d, "value", method = "all", ci = 0.90,
    bootstrap_replicates = 100, seed = 17
  ))
  expect_identical(
    ri$methods,
    c("nonparametric", "biweight", "type6", "parametric", "huber")
  )
  expect_true(all(vapply(ri[ri$methods], Negate(is.null), logical(1))))
  printed <- paste(capture.output(visible <- withVisible(print(ri))),
                   collapse = "\n")
  expect_false(visible$visible)
  expect_match(printed, "Nonparametric.*exact binomial")
  expect_match(printed, "Biweight.*Appendix B")
  expect_match(printed, "Type 6")
  expect_match(printed, "Huber")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- quietly_value(plot(ri))
  expect_s3_class(plotted, "ggplot")
})

test_that("robust bootstrap seed is reproducible and preserves RNG state", {
  d <- data.frame(value = seq(80, 120, length.out = 80))
  set.seed(909)
  before <- .Random.seed
  first <- reference_interval(d, "value", method = "biweight",
                              bootstrap_replicates = 100, seed = 44)
  expect_identical(.Random.seed, before)
  second <- reference_interval(d, "value", method = "biweight",
                               bootstrap_replicates = 100, seed = 44)
  expect_equal(first$biweight$ci_lower, second$biweight$ci_lower)
  expect_equal(first$biweight$ci_upper, second$biweight$ci_upper)
})

test_that("reference intervals report missing and duplicate IDs", {
  d <- data.frame(id = c(1, 1, 2, 2), value = c(1, NA, 3, 4))
  ri <- suppressWarnings(reference_interval(d, "value", id = "id"))
  expect_equal(ri$n_miss, 1)
  expect_equal(ri$n_complete, 3)
  expect_equal(ri$n_duplicates, 2)
  expect_equal(ri$dup_rows, c(3L, 4L))
})

test_that("outlier procedures identify a clear extreme observation", {
  d <- data.frame(value = c(rep(10, 8), 100))
  for (method in c("grubbs", "esd", "dixon", "iqr")) {
    result <- outliers_test(d, "value", method = method)
    expect_s3_class(result, "outliers_test")
    expect_true(9L %in% result$indices, info = method)
  }
})

test_that("normality statistics agree with their source packages", {
  set.seed(20260726)
  d <- data.frame(value = stats::rnorm(50))

  shapiro <- normal_test(d, "value", method = "shapiro")
  direct_shapiro <- stats::shapiro.test(d$value)
  expect_equal(shapiro$statistic, unname(direct_shapiro$statistic), tolerance = 1e-12)
  expect_equal(shapiro$p_value, direct_shapiro$p.value, tolerance = 1e-12)

  ad <- normal_test(d, "value", method = "ad")
  direct_ad <- nortest::ad.test(d$value)
  expect_equal(ad$statistic, unname(direct_ad$statistic), tolerance = 1e-12)
  expect_equal(ad$p_value, direct_ad$p.value, tolerance = 1e-12)
})

test_that("precision normality keeps the automatic method and statistic name", {
  set.seed(20260907)
  d <- data.frame(
    sample = rep("sample-1", 75),
    day = factor(rep(seq_len(5), each = 15)),
    run = factor(rep(rep(seq_len(3), each = 5), times = 5)),
    y = stats::rnorm(75)
  )
  p <- precision(d, y ~ day/run, by = "sample")
  auto <- quietly_value(normal(p, method = "auto"))
  result <- auto$normal[["sample-1"]]

  expect_identical(result$test_method, "ad")
  expect_identical(result$statistic_name, "A")
  expect_true(is.finite(result$statistic))
  expect_true(is.na(result$W))
  expect_match(
    paste(capture.output(print(auto$normal)), collapse = "\n"),
    "Anderson-Darling"
  )

  shapiro <- quietly_value(normal(p, method = "shapiro"))
  shapiro_result <- shapiro$normal[["sample-1"]]
  expect_identical(shapiro_result$test_method, "shapiro")
  expect_identical(shapiro_result$statistic_name, "W")
  expect_true(is.finite(shapiro_result$W) && shapiro_result$W > 0 && shapiro_result$W <= 1)
})

test_that("reference and distribution tools reject degenerate inputs", {
  expect_error(reference_interval(data.frame(x = letters[1:4]), "x"), "numeric")
  expect_equal(
    outliers_test(data.frame(x = 1:2), "x", method = "grubbs")$n_outliers,
    0L
  )
  constant <- normal_test(data.frame(x = rep(1, 5)), "x")
  expect_true(is.na(constant$test_method))
  expect_true(is.na(constant$statistic))
  expect_true(is.na(constant$p_value))
})
