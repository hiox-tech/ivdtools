test_that("ROC AUC agrees with the independent rank definition", {
  marker <- c(0.1, 0.4, 0.4, 0.8, 0.9, 1.2, 1.2, 1.5)
  truth <- c(0, 0, 1, 0, 1, 1, 0, 1)
  d <- data.frame(marker = marker, truth = truth)
  r0 <- roc(d, "marker", "truth")
  r <- quietly_value(auc(r0))

  expect_identical(r0$col_direction$marker, "geq")
  expect_equal(
    r$auc_list$marker$auc,
    manual_auc(marker, truth),
    tolerance = 1e-12
  )
})

test_that("ROC handles perfect, reversed, and tied markers", {
  truth <- rep(c(0, 1), each = 5)
  perfect <- data.frame(marker = c(1:5, 11:15), truth = truth)
  reverse <- data.frame(marker = -perfect$marker, truth = truth)
  tied <- data.frame(marker = rep(1, 10), truth = truth)

  a1 <- quietly_value(auc(roc(perfect, "marker", "truth")))$auc_list$marker$auc
  a2 <- quietly_value(auc(roc(reverse, "marker", "truth")))$auc_list$marker$auc
  a3 <- quietly_value(auc(roc(tied, "marker", "truth")))$auc_list$marker$auc

  expect_equal(a1, 1)
  expect_equal(a2, 1)
  expect_equal(a3, 0.5)
})

test_that("paired AUC comparison reproduces EP24-A2 Appendix D", {
  d <- data.frame(
    diagnosis = c(rep(0L, 22), rep(1L, 28)),
    oxldl = c(
      37,44,42,62,42,61,77,51,52,60,74,73,70,64,54,66,63,54,66,48,59,22,
      83,86,57,76,96,77,72,71,41,95,116,60,77,66,76,60,143,88,64,73,78,
      53,60,78,82,66,76,45
    ),
    ldl = c(
      2.1,2.35,3.91,5.4,3.31,3.9,4.38,2.85,3.67,1.48,2.6,3.25,3.76,
      3.5,2.66,4.45,5.27,3.57,3.74,2.78,3.15,3.01,5.88,4.05,3.75,3.21,
      4.11,4.15,2.31,2.57,2.6,4.22,7.55,2.74,4.57,3.51,3.08,2.95,5.71,
      3.92,3.38,1.29,3.71,3.22,3.4,3.4,4.03,3.09,3.47,3.57
    )
  )
  object <- quietly_value(auc(roc(d, c("oxldl", "ldl"), "diagnosis")))
  object <- quietly_value(auc_compare(
    object, c("oxldl", "ldl"), method = "ep24"
  ))
  ep24 <- object$auc_comparisons$oxldl_vs_ldl_ep24

  expect_s3_class(ep24, "roc_auc_comparison")
  expect_equal(ep24$difference, 0.237824675324675, tolerance = 1e-9)
  expect_equal(ep24$se, 0.075261470714225, tolerance = 1e-9)
  expect_equal(ep24$conf.int,
               c(lower = 0.0903149033012784, upper = 0.385334447348072),
               tolerance = 1e-12)
  expect_equal(ep24$statistic, 3.15997911106094, tolerance = 1e-9)
  expect_equal(ep24$p.value, 0.00157780450507665, tolerance = 1e-9)
  expect_equal(ep24$rating_correlation[["average"]], 0.5118177,
               tolerance = 1e-7)
  expect_equal(ep24$area_correlation, 0.4818177, tolerance = 1e-7)
  expect_true(ep24$lookup$clamped)
  expect_equal(ep24$lookup$average_auc, 0.7)

  object <- quietly_value(auc_compare(
    object, c("oxldl", "ldl"), method = "delong"
  ))
  delong <- object$auc_comparisons$oxldl_vs_ldl_delong
  expect_equal(delong$difference, 0.237824675324675, tolerance = 1e-9)
  expect_equal(delong$se, 0.0790442503119987, tolerance = 1e-9)
  expect_equal(delong$conf.int,
               c(lower = 0.082900791528189, upper = 0.392748559121162),
               tolerance = 1e-12)
  summary_text <- paste(capture.output(summary(object)), collapse = "\n")
  expect_match(summary_text, "EP24-A2 Hanley-McNeil")
  expect_match(summary_text, "DeLong")
})

test_that("paired AUC comparison accepts fitted MLR scores and validates curves", {
  d <- data.frame(
    x1 = seq(-2, 2, length.out = 20),
    x2 = rep(c(-1, 1), 10),
    truth = rep(c(0, 0, 1, 0, 1), 4)
  )
  raw <- roc(d, c("x1", "x2"), "truth")
  expect_error(auc_compare(raw, c("x1", "x2")), "Run auc")
  object <- quietly_value(auc(raw))
  object <- quietly_value(mlr(object, name = "combined"))
  object <- quietly_value(auc_compare(
    object, c("combined", "x1"), method = "ep24"
  ))
  object <- quietly_value(auc_compare(
    object, c("x1", "combined"), method = "delong"
  ))
  comparison <- object$auc_comparisons$x1_vs_combined_delong

  expect_true(
    object$auc_comparisons$combined_vs_x1_ep24$apparent_model_comparison
  )
  expect_identical(comparison$sources,
                   c(x1 = "marker", combined = "model"))
  expect_true(comparison$apparent_model_comparison)
  expect_match(paste(capture.output(print(comparison)), collapse = "\n"),
               "model optimism")
  expect_error(auc_compare(object, c("x1", "x1")), "distinct")
  expect_error(auc_compare(object, c("x1", "absent")), "not found")
  expect_error(auc_compare(object, c("x1", "combined"), method = "delong"),
               "already exists")
  expect_error(mlr(object, name = "x1"), "conflicts")
})

test_that("ROC cutoff rows contain independently computed confusion metrics", {
  d <- data.frame(marker = 1:8, truth = c(0, 0, 1, 0, 1, 1, 0, 1))
  r <- quietly_value(auc(roc(d, "marker", "truth")))
  r <- quietly_value(cutoff(r))
  row <- r$cutoff_table[r$cutoff_table$column == "marker" &
                          r$cutoff_table$cutoff == 5, ]
  pred <- as.integer(d$marker >= 5)

  expect_equal(row$tp, sum(pred == 1 & d$truth == 1))
  expect_equal(row$fp, sum(pred == 1 & d$truth == 0))
  expect_equal(row$tn, sum(pred == 0 & d$truth == 0))
  expect_equal(row$fn, sum(pred == 0 & d$truth == 1))
})

test_that("multivariate ROC logistic fit and prediction agree with glm", {
  d <- data.frame(
    x1 = seq(-2, 2, length.out = 20),
    x2 = rep(c(-1, 1), 10),
    truth = rep(c(0, 0, 1, 0, 1), 4)
  )
  r <- quietly_value(auc(roc(d, c("x1", "x2"), "truth")))
  r <- quietly_value(mlr(r, name = "model"))
  expected <- stats::glm(truth ~ x1 + x2, data = d, family = stats::binomial())

  expect_equal(
    stats::coef(r$mlr_results$model$fit),
    stats::coef(expected),
    tolerance = 1e-10
  )

  newdata <- data.frame(x1 = c(-0.5, 0.5), x2 = c(1, -1))
  observed <- predict(r, newdata = newdata, mlr = "model")
  expected_prob <- stats::predict(expected, newdata, type = "response")
  expect_num_equal(observed$pred_prob, expected_prob, 1e-10)
})

test_that("ROC validates reference and analysis ordering", {
  expect_error(roc(data.frame(x = 1:4, ref = 1:4), "x", "ref"), "quantitative")
  expect_error(roc(data.frame(x = letters[1:4], ref = c(0, 1, 0, 1)), "x", "ref"), "numeric")
  expect_error(
    roc(data.frame(x = 1:4, ref = rep(0, 4)), "x", "ref"),
    "at least one positive and one negative"
  )
  r <- roc(data.frame(x = 1:4, ref = c(0, 1, 0, 1)), "x", "ref")
  expect_error(cutoff(r), "auc")
  expect_error(mlr(r), "auc")

  converted <- continuous_to_binary(data.frame(x = c(1, 2, 3)), "x", 2)
  expect_identical(converted$x_binary, c(0L, 1L, 1L))
  expect_error(
    continuous_to_binary(data.frame(x = 1:3), "x", NA_real_),
    "finite"
  )
  expect_error(
    continuous_to_binary(data.frame(x = 1:3), "x", Inf),
    "finite"
  )
})

test_that("ROC plots cutoff-dependent metrics and stored optimal points", {
  obj <- roc(ivd_roc_example, cols = c("x1", "x2"), reference = "ref")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  cutoff_plot <- quietly_value(plot(obj, curves = "x1", type = "cutoff"))
  expect_s3_class(cutoff_plot, "ggplot")
  expect_error(
    suppressMessages(plot(obj, curves = "x1", type = "cutoff",
                          cutpoint = "Youden")),
    "Run cutoff"
  )

  obj <- quietly_value(auc(obj))
  obj <- quietly_value(cutoff(obj))
  for (method in c("Youden", "Corner", "all")) {
    expect_warning(
      p <- quietly_value(plot(obj, curves = "x1", type = "cutoff",
                              cutpoint = method)),
      NA
    )
    expect_s3_class(p, "ggplot")
  }

  expect_s3_class(
    quietly_value(plot(obj, curves = c("x1", "x2"), type = "roc",
                       cutpoint = "all")),
    "ggplot"
  )
})

test_that("ROC plots use one curve namespace for markers and MLR models", {
  obj <- roc(ivd_roc_example, cols = c("x1", "x2"), reference = "ref")
  obj <- quietly_value(auc(obj))
  obj <- quietly_value(mlr(obj, cols = c("x1", "x2"), name = "combined"))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(
    quietly_value(plot(obj, curves = c("x1", "combined"))),
    "ggplot"
  )
  expect_s3_class(quietly_value(plot(obj, curves = "all")), "ggplot")
  expect_s3_class(quietly_value(plot(obj, curves = "combined")), "ggplot")
  expect_error(
    suppressMessages(plot(obj, type = "cutoff")),
    "exactly one curve"
  )
  expect_error(
    suppressMessages(plot(obj, curves = "combined", type = "cutoff")),
    "only original marker curves"
  )
  expect_error(
    suppressMessages(plot(obj, curves = "absent", type = "cutoff")),
    "was not found"
  )
  expect_error(plot(obj, curves = c("x1", "x1")), "duplicate")
  expect_error(plot(obj, curves = c("all", "x1")), "alone")
  expect_error(plot(obj, cols = "x1"), "Unused argument.*cols")
  expect_error(plot(obj, mlr = "all"), "Unused argument.*mlr")
})
