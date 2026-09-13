test_that("interference replicate counts follow the two-sample normal formula", {
  out <- interference_replicates(5, c(7, 10), alpha = 0.05, power = 0.9,
                                 min_replicates = 2,
                                 label = c("seven", "ten"))
  expected <- 2 * ((qnorm(0.975) + qnorm(0.9)) * 5 / c(7, 10))^2

  expect_s3_class(out, "interference_replicates")
  expect_equal(out$calculated_replicates, expected)
  expect_equal(out$replicates_per_sample, pmax(2, ceiling(expected)))

  one_sided <- interference_replicates(5, 7, alternative = "one.sided")
  expect_equal(one_sided$z_alpha, qnorm(0.95))
})

test_that("EP07 paired analysis reproduces independent differences and intervals", {
  d <- data.frame(
    analyte = rep(c("A", "B"), each = 8),
    condition = rep(rep(c("control", "test"), each = 4), 2),
    result = c(9.8, 10.0, 10.1, 10.1, 11.0, 11.2, 11.1, 11.3,
               20.0, 20.2, 19.8, 20.0, 19.0, 19.2, 18.8, 19.0)
  )
  out <- interference_paired(d, "result", condition = "condition",
                             by = "analyte")
  a <- subset(d, analyte == "A")
  test <- a$result[a$condition == "test"]
  control <- a$result[a$condition == "control"]
  se <- sqrt(var(test) / length(test) + var(control) / length(control))
  estimate <- mean(test) - mean(control)
  critical <- qt(0.975, length(test) + length(control) - 2)

  expect_s3_class(out, "interference_paired")
  expect_equal(out$effects$absolute_difference[1], estimate)
  expect_equal(out$effects$absolute_se[1], se)
  expect_equal(out$effects$absolute_ci_lower[1], estimate - critical * se)
  expect_equal(out$effects$percent_difference[1], 100 * estimate / mean(control))
  expect_equal(nrow(out$cell_summary), 4)
})

test_that("spike and blank designs preserve their defining contrasts", {
  spike <- expand.grid(matrix = c("normal", "interferent"),
                       spike = c("before", "after"), replicate = 1:3,
                       KEEP.OUT.ATTRS = FALSE)
  spike$result <- with(spike,
    ifelse(matrix == "normal", 10, 12) +
      ifelse(spike == "after", ifelse(matrix == "normal", 10, 8), 0) +
      rep(c(-0.1, 0, 0.1), each = 4))
  spiked <- interference_paired(
    spike, "result", design = "analyte_spike", matrix = "matrix",
    spike = "spike", unspiked_level = "before", spiked_level = "after"
  )
  expect_equal(spiked$effects$normal_increment, 10, tolerance = 0.2)
  expect_equal(spiked$effects$interferent_increment, 8, tolerance = 0.2)
  expect_equal(spiked$effects$percent_double_difference, -20, tolerance = 3)

  blank_data <- data.frame(condition = rep(c("control", "test"), each = 3),
                           result = c(-0.1, 0, 0.1, 0.9, 1, 1.1))
  blank <- interference_paired(blank_data, "result", design = "blank",
                               condition = "condition")
  expect_identical(blank$effects$comparison, "solvent_control")
  expect_equal(blank$effects$apparent_signal, 1)

  theoretical <- interference_paired(
    subset(blank_data, condition == "test"), "result", design = "blank"
  )
  expect_identical(theoretical$effects$comparison, "theoretical_zero")
  expect_equal(theoretical$effects$apparent_signal, 1)
})

test_that("dose response supports point interpolation and linear inversion", {
  d <- expand.grid(interferent = c("A", "B"), concentration = c(0, 10, 20, 30),
                   replicate = 1:3, KEEP.OUT.ATTRS = FALSE)
  d$result <- 100 + ifelse(d$interferent == "A", 0.5, -0.25) * d$concentration +
    rep(c(-0.1, 0, 0.1), length.out = nrow(d))

  point <- interference_dose_response(d, "result", "concentration",
                                      by = "interferent")
  expect_s3_class(point, "interference_dose_response")
  a <- point$summary$stratum[point$summary$interferent == "A"][1]
  forward <- predict(point, newdata = 15, metric = "absolute", stratum = a)
  inverse <- predict(point, target = 7.5, metric = "absolute", stratum = a)
  expect_equal(forward$estimate, 7.5, tolerance = 0.2)
  expect_equal(inverse$estimate, 15, tolerance = 0.4)

  linear <- interference_dose_response(d, "result", "concentration",
                                       by = "interferent", method = "linear")
  slope <- subset(linear$coefficients,
                  interferent == "A" & term == ".concentration")$estimate
  expect_equal(slope, 0.5, tolerance = 0.02)
  expect_equal(predict(linear, newdata = 10, metric = "percent",
                       stratum = a)$estimate, 5, tolerance = 0.3)
  expect_true(is.na(predict(linear, newdata = 100, stratum = a)$estimate))
})

test_that("patient-specimen analysis averages replicates and models interference", {
  d <- expand.grid(id = paste0("S", 1:12),
                   method = c("evaluated", "comparative"), replicate = 1:2,
                   KEEP.OUT.ATTRS = FALSE)
  index <- as.integer(sub("S", "", d$id))
  d$group <- ifelse(index <= 6, "control", "test")
  d$bilirubin <- 2 * (index - 1)
  truth <- 50 + index
  d$result <- truth + ifelse(d$method == "evaluated", 0.08 * d$bilirubin, 0) +
    ifelse(d$replicate == 1, -0.05, 0.05)
  out <- interference_patient(
    d, "id", "group", "method", "result",
    test_group = "test", control_group = "control",
    evaluated_method = "evaluated", comparative_method = "comparative",
    interferent_cols = "bilirubin"
  )

  expect_s3_class(out, "interference_patient")
  expect_equal(out$specimen_results$absolute_difference,
               0.08 * out$specimen_results$bilirubin, tolerance = 1e-10)
  slope <- subset(out$regressions,
                  interferent == "bilirubin" & term == "concentration")$estimate
  expect_equal(slope, 0.08, tolerance = 1e-10)
  expect_equal(out$group_summary$n, c(6, 6))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_centered_plot(plot(out, type = "difference"))
  expect_centered_plot(plot(out, type = "interferent", interferent = "bilirubin"))
})

test_that("interference methods plot and reject incomplete designs", {
  d <- data.frame(condition = rep(c("control", "test"), each = 3),
                  result = c(1, 1.1, 0.9, 1.4, 1.5, 1.6))
  paired <- interference_paired(d, "result", condition = "condition")
  dose <- interference_dose_response(
    data.frame(concentration = rep(c(0, 10, 20), each = 3),
               result = rep(c(100, 105, 110), each = 3) + rep(c(-0.1, 0, 0.1), 3)),
    "result", "concentration"
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  capture.output(printed <- withVisible(print(paired)))
  expect_false(printed$visible)
  expect_identical(printed$value, paired)
  expect_centered_plot(plot(paired))
  for (type in c("response", "absolute", "percent"))
    expect_centered_plot(plot(dose, type = type))

  expect_error(interference_paired(d[1:3, ], "result", condition = "condition"),
               "No stratum satisfied")
  expect_error(interference_dose_response(
    transform(dose$analysis_data, concentration = 1), "result", "concentration"),
    "No stratum could")
  expect_error(predict(dose), "exactly one")
  expect_error(predict(dose, newdata = 1, target = 1), "exactly one")
  expect_error(interference_replicates(0, 1), "positive finite")
})
