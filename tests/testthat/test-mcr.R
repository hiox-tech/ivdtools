test_that("MCR constructor filters only incomplete pairs and validates weights", {
  d <- data.frame(
    id = 1:5,
    x = c(1, 2, NA, 4, 5),
    y = c(1.1, 2.1, 3.1, Inf, 5.1),
    w = c(1, 2, 3, 4, 5)
  )
  m <- mcr(d, "id", "x", "y", weights = "w")
  expect_equal(m$complete_idx, c(1L, 2L, 5L))
  expect_equal(m$weights, c(1, 2, 5))
  expect_equal(m$n, 3)

  d$w[2] <- -1
  expect_error(mcr(d, "id", "x", "y", weights = "w"), "non-negative")
  expect_error(mcr(list(), "id", "x", "y"), "data frame")
  expect_error(mcr(data.frame(id = 1:2, x = 1:2), "id", "x", "y"), "column name")
})

test_that("MCR requires one summarized row per sample", {
  raw <- data.frame(
    id = rep(1:4, each = 2),
    candidate = 1:8,
    reference = 2:9
  )
  expect_error(
    mcr(raw, "id", "candidate", "reference"),
    "one summarized row|replicate_to_mean"
  )
})

test_that("replicate SD columns estimate lambda on the mean scale", {
  raw <- data.frame(
    id = rep(1:8, each = 3),
    reference = rep(seq(10, 80, 10), each = 3) + rep(c(-1, 0, 1), 8),
    candidate = rep(2 + 1.05 * seq(10, 80, 10), each = 3) +
      rep(c(-0.5, 0, 0.5), 8)
  )
  summarized <- replicate_to_mean(
    raw, "id", c("candidate", "reference")
  )
  object <- mcr(
    summarized, "x", "candidate_mean", "reference_mean",
    candidate_sd = "candidate_sd", reference_sd = "reference_sd",
    candidate_n = "candidate_n", reference_n = "reference_n"
  )
  fit <- quietly_value(regression(
    object, method = "deming", lambda = "replicates"
  ))$regression
  expected_lambda <- (0.5^2 / 3) / (1^2 / 3)
  direct <- quietly_value(regression(
    object, method = "deming", lambda = expected_lambda
  ))$regression
  expect_equal(fit$lambda, expected_lambda, tolerance = 1e-12)
  expect_equal(c(fit$intercept, fit$slope),
               c(direct$intercept, direct$slope), tolerance = 1e-12)
  expect_equal(fit$lambda_metadata$source, "summarized SD data")

  without_n <- mcr(
    summarized, "x", "candidate_mean", "reference_mean",
    candidate_sd = "candidate_sd", reference_sd = "reference_sd"
  )
  fit_without_n <- quietly_value(regression(
    without_n, "deming", lambda = "replicates"
  ))$regression
  expect_equal(fit_without_n$lambda, 0.25, tolerance = 1e-12)
  expect_error(
    regression(mcr(summarized, "x", "candidate_mean", "reference_mean"),
               "deming", lambda = "replicates"),
    "SD|replicate_to_mean"
  )

  varying <- summarized
  varying$candidate_n[1] <- 4L
  expect_error(
    regression(mcr(
      varying, "x", "candidate_mean", "reference_mean",
      candidate_sd = "candidate_sd", reference_sd = "reference_sd",
      candidate_n = "candidate_n", reference_n = "reference_n"
    ), "deming", lambda = "replicates"),
    "constant across samples|precision_profile"
  )
})

test_that("split_by_cutpoints creates stable left-closed data-frame lists", {
  d <- data.frame(id = 1:9, reference = 1:9,
                  candidate = 2 + 1.1 * (1:9))
  parts <- split_by_cutpoints(d, "reference", c(4, 7))

  expect_s3_class(parts, "cutpoint_split")
  expect_named(parts, paste0("segment_", 1:3))
  expect_true(all(vapply(parts, is.data.frame, logical(1))))
  expect_equal(lapply(parts, `[[`, "reference"),
               list(segment_1 = 1:3, segment_2 = 4:6, segment_3 = 7:9))
  info <- attr(parts, "segment_info")
  expect_equal(info$interval, c("[-Inf, 4)", "[4, 7)", "[7, Inf)"))
  expect_equal(info$n, c(3L, 3L, 3L))
  expect_equal(attr(parts, "cutpoints"), c(4, 7))
  expect_match(paste(capture.output(print(parts)), collapse = "\n"),
               "segment_2.*\\[4, 7\\).*n = 3")

  labelled <- split_by_cutpoints(
    d, "reference", c(4, 7), labels = c("low", "middle", "high")
  )
  expect_named(labelled, c("low", "middle", "high"))
  expect_error(
    split_by_cutpoints(d, "reference", c(7, 4)),
    "strictly increasing"
  )
  expect_error(
    split_by_cutpoints(transform(d, reference = as.character(reference)),
                       "reference", 4),
    "numeric column"
  )
})

test_that("split_by_cutpoints handles empty and non-finite groups explicitly", {
  d <- data.frame(id = 1:4, reference = c(1, 2, NA, 4), candidate = 1:4)
  expect_error(split_by_cutpoints(d, "reference", c(0, 3, 10)),
               "non-finite")

  parts <- split_by_cutpoints(
    d, "reference", c(0, 3, 10), na_action = "drop"
  )
  expect_equal(attr(parts, "n_excluded"), 1L)
  expect_equal(attr(parts, "segment_info")$n, c(0L, 2L, 1L, 0L))
  expect_length(parts, 4L)

  nonempty <- split_by_cutpoints(
    d, "reference", c(0, 3, 10), na_action = "drop", drop_empty = TRUE
  )
  expect_named(nonempty, c("segment_2", "segment_3"))
})

test_that("each pre-split data frame supports a complete MCR workflow", {
  d <- data.frame(
    id = 1:12,
    reference = 1:12,
    candidate = 1 + 1.2 * (1:12) + rep(c(-0.2, 0.2), 6)
  )
  parts <- split_by_cutpoints(d, "reference", c(5, 9))
  objects <- lapply(parts, function(part) {
    mcr(part, "id", "candidate", "reference")
  })
  expect_named(objects, names(parts))

  analysed <- quietly_value(describe(objects$segment_2))
  analysed <- quietly_value(correlation(analysed))
  analysed <- quietly_value(regression(analysed, "ols"))
  analysed <- quietly_value(bland_altman(analysed))
  expected <- stats::lm(candidate ~ reference, data = d[5:8, ])

  expect_equal(analysed$n, 4L)
  expect_equal(c(analysed$regression$intercept,
                 analysed$regression$slope),
               unname(coef(expected)), tolerance = 1e-12)
  expect_equal(analysed$regression$complete_idx, 1:4)

  ba <- analysed$bland_altman
  expect_equal(ba$n, 4)
  expect_equal(ba$y_vals, d$candidate[5:8] - d$reference[5:8])
})

test_that("MCR correlations agree with stats::cor.test", {
  m <- mcr(ivd_mcr_example, "sid", "test", "ref")
  for (method in c("pearson", "spearman", "kendall")) {
    observed <- quietly_value(correlation(m, method = method))$correlation
    expected <- suppressWarnings(stats::cor.test(
      ivd_mcr_example$test,
      ivd_mcr_example$ref,
      method = method
    ))
    expect_num_equal(observed$estimate, unname(expected$estimate), 1e-12)
    expect_num_equal(observed$p.value, expected$p.value, 1e-10)
  }
})

test_that("OLS and WLS regression agree with stats::lm", {
  d <- transform(ivd_mcr_example, w = seq_len(nrow(ivd_mcr_example)))
  m <- mcr(d, "sid", "test", "ref", weights = "w")

  observed_ols <- quietly_value(regression(m, method = "ols"))$regression
  expected_ols <- stats::lm(test ~ ref, data = d)
  expect_num_equal(
    c(observed_ols$intercept, observed_ols$slope),
    stats::coef(expected_ols),
    1e-12
  )
  expect_num_equal(observed_ols$residuals, stats::residuals(expected_ols), 1e-12)

  observed_wls <- quietly_value(
    regression(m, method = "wls", weights = "w")
  )$regression
  expected_wls <- stats::lm(test ~ ref, data = d, weights = w)
  expect_num_equal(
    c(observed_wls$intercept, observed_wls$slope),
    stats::coef(expected_wls),
    1e-12
  )
})

test_that("Deming regression recovers an exact linear relationship", {
  d <- data.frame(id = 1:10, reference = 1:10, candidate = 1 + 2 * (1:10))
  m <- mcr(d, "id", "candidate", "reference")
  fit <- quietly_value(regression(m, method = "deming"))$regression
  expect_equal(fit$intercept, 1, tolerance = 1e-10)
  expect_equal(fit$slope, 2, tolerance = 1e-10)
})

test_that("Bland-Altman estimates follow independent formulas", {
  m <- mcr(ivd_mcr_example, "sid", "test", "ref")
  result <- quietly_value(bland_altman(m))$bland_altman
  differences <- ivd_mcr_example$test - ivd_mcr_example$ref
  z <- stats::qnorm(0.975)

  expect_equal(result$mean_diff, mean(differences), tolerance = 1e-14)
  expect_equal(result$sd_diff, stats::sd(differences), tolerance = 1e-14)
  expect_equal(
    unname(result$loa),
    mean(differences) + c(-1, 1) * z * stats::sd(differences),
    tolerance = 1e-12
  )
  t_value <- stats::qt(0.975, length(differences) - 1)
  expect_equal(
    unname(result$mean_diff_ci),
    mean(differences) + c(-1, 1) * t_value *
      stats::sd(differences) / sqrt(length(differences)),
    tolerance = 1e-12
  )
  expect_equal(result$center_ci, result$mean_diff_ci)
})

test_that("Bland-Altman percentage denominator is independent of x-axis", {
  d <- data.frame(
    id = 1:6,
    reference = c(10, 20, 30, 40, 50, 60),
    candidate = c(11, 18, 33, 38, 55, 63)
  )
  m <- mcr(d, "id", "candidate", "reference")

  by_reference <- quietly_value(bland_altman(
    m, type = "percent", x_axis = "mean",
    percent_denominator = "reference"
  ))$bland_altman
  by_mean <- quietly_value(bland_altman(
    m, type = "percent", x_axis = "mean",
    percent_denominator = "mean"
  ))$bland_altman
  by_candidate <- quietly_value(bland_altman(
    m, type = "percent", x_axis = "mean",
    percent_denominator = "candidate"
  ))$bland_altman

  expect_equal(by_reference$x_vals, by_mean$x_vals)
  expect_equal(by_reference$y_vals,
               100 * (d$candidate - d$reference) / d$reference)
  expect_equal(by_mean$y_vals,
               100 * (d$candidate - d$reference) /
                 ((d$candidate + d$reference) / 2))
  expect_equal(by_candidate$y_vals,
               100 * (d$candidate - d$reference) / d$candidate)
  expect_equal(by_reference$percent_denominator, "reference")
  expect_equal(by_mean$percent_denominator, "mean")
  expect_match(by_mean$y_label, "mean", fixed = TRUE)

  default <- quietly_value(
    bland_altman(m, type = "percent", x_axis = "mean")
  )$bland_altman
  expect_equal(default$y_vals, by_reference$y_vals)
  expect_error(
    bland_altman(m, type = "percent", percent_denominator = "other"),
    "arg"
  )
})

test_that("Bland-Altman nonparametric method uses empirical LoA and reproducible bootstrap", {
  m <- mcr(ivd_mcr_example, "sid", "test", "ref")
  a <- quietly_value(bland_altman(m, method = "nonparametric", bootstrap_replicates = 200, seed = 42))$bland_altman
  b <- quietly_value(bland_altman(m, method = "nonparametric", bootstrap_replicates = 200, seed = 42))$bland_altman
  d <- ivd_mcr_example$test - ivd_mcr_example$ref
  expect_equal(a$method, "nonparametric")
  expect_equal(unname(a$loa), as.numeric(quantile(d, c(.025, .975))))
  expect_equal(a$loa_ci, b$loa_ci)
  ranks <- seq_len(floor((length(d) + 1) / 2))
  coverage <- 1 - 2 * stats::pbinom(ranks - 1, length(d), 0.5)
  lower_rank <- max(ranks[coverage >= 0.95])
  expected_ci <- sort(d)[c(lower_rank, length(d) - lower_rank + 1)]
  expect_equal(unname(a$median_diff_ci), expected_ci)
  expect_equal(a$center_ci, a$median_diff_ci)
})

test_that("Passing-Bablok uses Algorithm I vertical slopes and K adjustment", {
  d <- data.frame(
    id = 1:6,
    reference = c(1, 1, 2, 3, 4, 5),
    candidate = c(4, 1, 2, 5, 5, 7)
  )
  fit <- quietly_value(
    regression(mcr(d, "id", "candidate", "reference"), method = "pb")
  )$regression

  expect_equal(fit$k_adjustment, 2)
  expect_equal(fit$n_slopes, 15)
  expect_equal(fit$slope, 1.5, tolerance = 1e-12)
  expect_equal(fit$intercept, -0.5, tolerance = 1e-12)
})

test_that("residual WLS reproduces four EP09-style weight updates", {
  d <- data.frame(
    id = 1:20,
    reference = 1:20,
    candidate = 3 + 1.05 * (1:20) +
      c(-.2, .1, -.4, .3, -.5, .6, -.7, .8, -.9, 1,
        -1.1, 1.2, -1.3, 1.4, -1.5, 1.6, -1.7, 1.8, -1.9, 2)
  )
  observed <- quietly_value(regression(
    mcr(d, "id", "candidate", "reference"), method = "wls",
    weights = "residual", weight_iterations = 4
  ))$regression

  expected <- stats::lm(candidate ~ reference, data = d)
  for (i in 1:4) {
    sd_model <- stats::lm(abs(stats::residuals(expected)) ~ d$reference)
    expected_weights <- 1 / stats::fitted(sd_model)^2
    expected <- stats::lm(candidate ~ reference, data = d,
                          weights = expected_weights)
  }
  expect_equal(c(observed$intercept, observed$slope),
               unname(coef(expected)), tolerance = 1e-12)
  expect_equal(observed$weights, unname(expected_weights), tolerance = 1e-12)
  expect_equal(nrow(observed$residual_weight_history), 5)
  expect_equal(observed$weight_iterations, 4L)
  expect_error(
    regression(mcr(d, "id", "candidate", "reference"), "ols",
               weights = "residual"),
    "only for method"
  )
})

test_that("MCR estimates are invariant to row order and transform predictably", {
  d <- ivd_mcr_example
  base <- quietly_value(regression(mcr(d, "sid", "test", "ref"), "ols"))$regression
  perm <- quietly_value(regression(
    mcr(d[c(12:1), ], "sid", "test", "ref"),
    "ols"
  ))$regression
  expect_equal(base$slope, perm$slope, tolerance = 1e-12)
  expect_equal(base$intercept, perm$intercept, tolerance = 1e-12)

  scaled <- transform(d, test = test * 10, ref = ref * 5)
  scaled_fit <- quietly_value(regression(
    mcr(scaled, "sid", "test", "ref"),
    "ols"
  ))$regression
  expect_equal(scaled_fit$slope, base$slope * 2, tolerance = 1e-12)
  expect_equal(scaled_fit$intercept, base$intercept * 10, tolerance = 1e-12)
})

test_that("reference is X and candidate is Y throughout regression and prediction", {
  d <- data.frame(
    id = 1:10,
    reference = 1:10,
    candidate = 3 + 2 * (1:10) +
      c(0.1, -0.1, 0.2, -0.2, 0.1, 0, -0.1, 0.2, -0.2, 0)
  )
  m <- mcr(d, "id", candidate = "candidate", reference = "reference")
  fit <- quietly_value(regression(m, method = "ols"))
  expected <- stats::lm(candidate ~ reference, data = d)

  expect_equal(
    c(fit$regression$intercept, fit$regression$slope),
    unname(coef(expected)), tolerance = 1e-12
  )

  forward <- predict(fit, reference = c(4, 8))
  expect_named(forward, c("reference", "candidate_pred"))
  expect_equal(
    forward$candidate_pred,
    unname(stats::predict(expected, data.frame(reference = c(4, 8)))),
    tolerance = 1e-12
  )

  inverse <- predict(fit, candidate = forward$candidate_pred, inverse = TRUE)
  expect_named(inverse, c("candidate", "reference_pred"))
  expect_equal(inverse$reference_pred, c(4, 8), tolerance = 1e-12)

  biased <- quietly_value(bias(fit, mdl = 5))
  expect_equal(
    biased$bias$bias,
    unname(stats::predict(expected, data.frame(reference = 5))) - 5,
    tolerance = 1e-12
  )
})

test_that("built-in x and y weights follow reference X and candidate Y", {
  d <- data.frame(
    id = 1:8,
    reference = 1:8,
    candidate = c(4, 5, 8, 10, 12, 15, 16, 20)
  )
  m <- mcr(d, "id", candidate = "candidate", reference = "reference")

  by_x <- quietly_value(regression(m, method = "wls", weights = "1/x^2"))$regression
  expected_x <- stats::lm(candidate ~ reference, data = d,
                          weights = 1 / reference^2)
  expect_equal(c(by_x$intercept, by_x$slope), unname(coef(expected_x)),
               tolerance = 1e-12)

  by_y <- quietly_value(regression(m, method = "wls", weights = "1/y^2"))$regression
  expected_y <- stats::lm(candidate ~ reference, data = d,
                          weights = 1 / candidate^2)
  expect_equal(c(by_y$intercept, by_y$slope), unname(coef(expected_y)),
               tolerance = 1e-12)
})

test_that("scatter, regression, and bias plots use reference X and candidate Y", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(
    id = 1:10,
    reference = 1:10,
    candidate = 3 + 2 * (1:10) +
      c(0.1, -0.1, 0.2, -0.2, 0.1, 0, -0.1, 0.2, -0.2, 0)
  )
  m <- mcr(d, "id", candidate = "candidate", reference = "reference")
  fit <- quietly_value(regression(m, method = "ols"))

  capture.output(scatter <- plot(m, type = "scatter"))
  capture.output(regression_plot <- plot(fit, type = "regression"))
  biased <- quietly_value(bias(fit, mdl = 5))
  capture.output(bias_plot <- plot(biased, type = "bias"))

  expect_equal(scatter$labels$x, "reference")
  expect_equal(scatter$labels$y, "candidate")
  expect_equal(regression_plot$labels$x, "reference")
  expect_equal(regression_plot$labels$y, "candidate")
  expect_equal(bias_plot$labels$x, "reference")
  expect_match(bias_plot$labels$y, "candidate - reference", fixed = TRUE)

  line_data <- regression_plot$layers[[3]]$data
  expect_equal(
    line_data$y,
    fit$regression$intercept + fit$regression$slope * line_data$x,
    tolerance = 1e-12
  )
})

test_that("Bland-Altman plot includes the center confidence ribbon", {
  skip_if_not_installed("ggplot2")
  m <- mcr(ivd_mcr_example, "sid", "test", "ref")
  analysed <- quietly_value(bland_altman(m))
  capture.output(p <- plot(analysed, type = "bland_altman"))

  ribbon_data <- p$layers[[2]]$data
  expect_equal(unique(ribbon_data$center_lo),
               unname(analysed$bland_altman$center_ci[1]))
  expect_equal(unique(ribbon_data$center_hi),
               unname(analysed$bland_altman$center_ci[2]))
})

test_that("constant-CV WLS and Deming use their EP09c algorithms", {
  skip_if_not_installed("ppwdeming", minimum_version = "3.0.2")
  d <- data.frame(
    id = 1:12,
    reference = seq(2, 24, 2),
    candidate = 0.4 + 1.08 * seq(2, 24, 2) *
      c(1.02, .98, 1.01, .97, 1.03, .99, 1.01, .98, 1.02, .99, 1.01, 1)
  )
  object <- mcr(d, "id", "candidate", "reference")

  wls <- quietly_value(regression(
    object, "wls", weights = "constant_cv"
  ))$regression
  expected_wls <- stats::lm(candidate ~ reference, d,
                            weights = 1 / reference^2)
  expect_equal(c(wls$intercept, wls$slope), unname(coef(expected_wls)),
               tolerance = 1e-12)

  fit <- quietly_value(regression(
    object, "wdeming", weights = "constant_cv", lambda = 2
  ))$regression
  direct <- ppwdeming::WD_Linnet(
    d$reference, d$candidate, lambda = 1 / 2, getCI = FALSE
  )
  expect_equal(unname(fit$intercept), direct$alpha, tolerance = 1e-10)
  expect_equal(unname(fit$slope), direct$beta, tolerance = 1e-10)
  expect_true(fit$convergence$converged)
  expect_length(fit$latent_concentration, nrow(d))
  expect_true(all(is.finite(fit$weights) & fit$weights > 0))
  expect_error(
    regression(mcr(transform(d, reference = reference - 2), "id",
                   "candidate", "reference"), "wdeming",
               weights = "constant_cv"),
    "strictly positive"
  )
})

test_that("precision-profile Deming matches ppwdeming for variance functions", {
  skip_if_not_installed("ppwdeming", minimum_version = "3.0.2")
  d <- data.frame(
    id = 1:10, reference = seq(5, 45, length.out = 10),
    candidate = 1 + 1.04 * seq(5, 45, length.out = 10) +
      c(-.4, .2, -.3, .5, -.2, .3, -.1, .4, -.3, .2)
  )
  gv <- function(z) 0.04 + 0.002 * z^2
  hv <- function(z) 0.09 + 0.003 * z^2
  fit <- quietly_value(regression(
    mcr(d, "id", "candidate", "reference"), "wdeming",
    weights = "precision_profile",
    variance_models = list(reference = gv, candidate = hv)
  ))$regression
  direct <- ppwdeming::PWD_known(
    d$reference, d$candidate,
    gfun = function(z, p) gv(z), hfun = function(z, p) hv(z),
    gparms = numeric(), hparms = numeric(), getCI = FALSE
  )
  expect_equal(unname(fit$intercept), direct$alpha, tolerance = 1e-8)
  expect_equal(unname(fit$slope), direct$beta, tolerance = 1e-8)
  expect_length(fit$latent_concentration, nrow(d))
  expect_true(all(fit$reference_variance > 0))
  expect_error(
    regression(mcr(d, "id", "candidate", "reference"), "wdeming",
      weights = "precision_profile",
      variance_models = list(reference = function(z) rep(0, length(z)),
                             candidate = hv)),
    "positive finite|non-positive"
  )

  old_wd <- setwd(tempdir())
  on.exit(setwd(old_wd), add = TRUE)
  file.create("stdout.log")
  profile_data <- data.frame(
    mean = c(5, 10, 20, 30, 40, 50),
    SD = c(.24, .53, .96, 1.60, 1.90, 2.70)
  )
  profile_data$CV <- 100 * profile_data$SD / profile_data$mean
  profile_fit <- suppressMessages(
    getFromNamespace(".sadler", "ivdtools")(
      profile_data, type = "sd", model.no = 2, quiet = TRUE
    )
  )
  precision_object <- structure(
    list(profile = list(fits = list(total = profile_fit))),
    class = "precision"
  )
  beta1 <- unname(stats::coef(profile_fit$models$model2)[1L])
  from_object <- quietly_value(regression(
    mcr(d, "id", "candidate", "reference"), "wdeming",
    weights = "precision_profile",
    variance_models = list(reference = precision_object,
                           candidate = precision_object)
  ))$regression
  from_function <- quietly_value(regression(
    mcr(d, "id", "candidate", "reference"), "wdeming",
    weights = "precision_profile",
    variance_models = list(
      reference = function(z) beta1 * z^2,
      candidate = function(z) beta1 * z^2
    )
  ))$regression
  expect_equal(c(from_object$intercept, from_object$slope),
               c(from_function$intercept, from_function$slope),
               tolerance = 1e-6)
  expect_equal(unname(from_object$profile_sources),
               rep("precision:total", 2))
})

test_that("Bland-Altman supports interpolated median and HL pseudomedian", {
  d <- data.frame(
    id = 1:20, reference = seq(10, 200, 10),
    candidate = seq(10, 200, 10) *
      (1 + c(-4, -3, -2.5, -2, -1.5, -1, -.8, -.5, -.2, 0,
             .1, .3, .5, .8, 1, 1.5, 2, 2.5, 3, 4) / 100)
  )
  object <- mcr(d, "id", "candidate", "reference")
  ordinary <- quietly_value(bland_altman(
    object, type = "percent", method = "nonparametric",
    median_ci_method = "interpolated", bootstrap_replicates = 100, seed = 1
  ))$bland_altman
  sorted <- sort(100 * (d$candidate - d$reference) / d$reference)
  ranks <- c((21 - qnorm(.975) * sqrt(20)) / 2,
             (21 + qnorm(.975) * sqrt(20)) / 2)
  interpolate <- function(r) sorted[floor(r)] + (r - floor(r)) *
    (sorted[ceiling(r)] - sorted[floor(r)])
  expect_equal(unname(ordinary$center_ci), vapply(ranks, interpolate, 0),
               tolerance = 1e-12)
  expect_equal(ordinary$center_ci_ranks, setNames(ranks, c("lower", "upper")))

  hl <- quietly_value(bland_altman(
    object, type = "percent", method = "hodges_lehmann",
    bootstrap_replicates = 100, seed = 1
  ))$bland_altman
  expected <- suppressWarnings(stats::wilcox.test(
    ordinary$y_vals, conf.int = TRUE, exact = FALSE, correct = FALSE
  ))
  expect_equal(hl$hl_estimate, unname(expected$estimate), tolerance = 1e-10)
  expect_equal(unname(hl$hl_ci), as.numeric(expected$conf.int),
               tolerance = 1e-10)
  expect_equal(hl$median_diff, stats::median(ordinary$y_vals))
  expect_equal(hl$center_estimand, "Hodges-Lehmann pseudomedian")
})

test_that("regression bootstrap methods are paired and reproducible", {
  d <- data.frame(
    id = 1:10, reference = 1:10,
    candidate = 1 + 1.1 * (1:10) +
      c(-.2, .3, -.1, .2, -.3, .1, -.2, .2, -.1, .1),
    weight = seq(1, 2, length.out = 10)
  )
  object <- mcr(d, "id", "candidate", "reference", weights = "weight")
  percentile_a <- quietly_value(regression(
    object, "ols", ci_method = "bootstrap_percentile",
    bootstrap_replicates = 5000, seed = 42, mdl = c(3, 7)
  ))$regression
  percentile_b <- quietly_value(regression(
    object, "ols", ci_method = "bootstrap_percentile",
    bootstrap_replicates = 5000, seed = 42, mdl = c(3, 7)
  ))$regression
  expect_equal(percentile_a$parameter_samples,
               percentile_b$parameter_samples)
  expect_equal(percentile_a$bootstrap$effective, 5000)
  expect_equal(dim(percentile_a$bootstrap$bias_samples), c(5000, 2))

  standard <- quietly_value(regression(
    object, "wls", weights = "weight", ci_method = "bootstrap_standard",
    bootstrap_replicates = 5000, seed = 43
  ))$regression
  expect_equal(standard$ci_method, "bootstrap_standard")
  expect_true(all(is.finite(c(standard$intercept_ci, standard$slope_ci))))
  expect_error(
    regression(object, "ols", ci_method = "bootstrap_percentile",
               bootstrap_replicates = 4999),
    "at least 5000"
  )

  for (regression_method in c("deming", "wdeming", "pb")) {
    extra <- if (regression_method == "wdeming")
      list(weights = "weight") else list()
    boot_fit <- quietly_value(do.call(regression, c(list(
      x = object, method = regression_method,
      ci_method = "bootstrap_percentile",
      bootstrap_replicates = 5000, seed = 44
    ), extra)))$regression
    expect_equal(boot_fit$bootstrap$effective, 5000,
                 info = regression_method)
    expect_true(all(is.finite(c(boot_fit$intercept_ci,
                                boot_fit$slope_ci))),
                info = regression_method)
  }
})

test_that("Passing-Bablok bias defaults to percentile bootstrap", {
  d <- data.frame(
    id = 1:9, reference = 1:9,
    candidate = c(2.1, 3.0, 4.2, 5.1, 6.4, 7.2, 8.0, 9.3, 10.1)
  )
  fit <- quietly_value(regression(
    mcr(d, "id", "candidate", "reference"), "pb"
  ))
  biased <- quietly_value(bias(
    fit, mdl = c(3, 7), interval = "confidence",
    bootstrap_replicates = 5000, seed = 99
  ))
  expect_equal(attr(biased$bias, "ci_method"), "bootstrap_percentile")
  expect_true(all(is.finite(biased$bias$ci_lower)))
  predicted <- predict(
    fit, reference = c(3, 7), interval = "confidence",
    bootstrap_replicates = 5000, seed = 99
  )
  expect_equal(attr(predicted, "ci_method"), "bootstrap_percentile")
  expect_true(all(predicted$ci_lower <= predicted$candidate_pred &
                  predicted$candidate_pred <= predicted$ci_upper))
})

test_that("outlier percent denominator matches Bland-Altman definition", {
  d <- data.frame(id = 1:8, reference = 10:17,
                  candidate = c(11:17, 30))
  object <- mcr(d, "id", "candidate", "reference")
  for (base in c("reference", "mean", "candidate")) {
    out <- quietly_value(outlier(
      object, method = "iqr", type = "percent",
      percent_denominator = base
    ))$outlier
    ba <- quietly_value(bland_altman(
      object, type = "percent", percent_denominator = base
    ))$bland_altman
    expect_equal(out$values, ba$y_vals[match(out$indices, object$complete_idx)])
    expect_equal(out$percent_denominator, base)
  }
})
