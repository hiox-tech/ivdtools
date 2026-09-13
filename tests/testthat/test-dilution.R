test_that("dilution recovery follows the EP34 equations", {
  d <- data.frame(
    sample = rep(c("A", "B"), each = 6),
    factor = rep(rep(c(1, 2), each = 3), 2),
    result = c(99, 100, 101, 49, 50, 51,
               198, 200, 202, 99, 100, 101)
  )
  out <- dilution_recovery(d, "sample", "result", "factor")

  expect_s3_class(out, "dilution_recovery")
  expect_equal(out$results$recovery, rep(100, 4), tolerance = 1e-12)
  expect_equal(out$results$corrected_concentration,
               c(100, 200, 100, 200), tolerance = 1e-12)
  expect_equal(out$summary$mean_recovery, c(100, 100), tolerance = 1e-12)
  expect_false(any(grepl("acceptable", names(out$results), ignore.case = TRUE)))
  expect_false(any(grepl("acceptable", names(out$summary), ignore.case = TRUE)))
})

test_that("assigned values and nonzero diluent are handled by mass balance", {
  d <- data.frame(sample = "A", factor = rep(2, 3),
                  assigned = 100, result = c(54, 55, 56))
  out <- dilution_recovery(d, "sample", "result", "factor",
                           assigned = "assigned", diluent_value = 10)
  expect_equal(out$results$target, 55)
  expect_equal(out$results$corrected_concentration, 100)
  expect_equal(out$results$recovery, 100)

  expect_error(
    dilution_recovery(transform(d, assigned = c(100, 101, 100)),
                      "sample", "result", "factor", assigned = "assigned"),
    "exactly one"
  )
  expect_error(
    dilution_recovery(subset(d, factor > 1), "sample", "result", "factor"),
    "no undiluted"
  )
})

test_that("spike concentration uses the EP34 mixture equation", {
  out <- spike_concentration(sample_concentration = 10, sample_volume = 1,
                             stock_concentration = 40, spike_volume = 0.1)
  expect_equal(out$added_concentration, 40 * 0.1 / 1.1)
  expect_equal(out$final_concentration, (10 * 1 + 40 * 0.1) / 1.1)
  expect_error(spike_concentration(10, 0, 40, 0.1), "positive")
})

test_that("spike recovery compares analyte spike with solvent control", {
  d <- data.frame(
    sample = rep(c("S1", "S2", "S3"), each = 4),
    condition = rep(rep(c("solvent", "analyte"), each = 2), 3),
    result = c(9.7, 9.9, 12.9, 13.3,
               9.5, 9.5, 12.8, 12.6,
               8.7, 9.1, 12.2, 12.4)
  )
  added <- 40 * 0.1 / 1.1
  out <- spike_recovery(d, "sample", "result", "condition",
                        analyte_level = "analyte", solvent_level = "solvent",
                        added_concentration = added)
  expect_s3_class(out, "spike_recovery")
  expect_equal(out$results$observed_increase, c(3.3, 3.2, 3.4),
               tolerance = 1e-12)
  expect_equal(out$results$recovery, 100 * c(3.3, 3.2, 3.4) / added,
               tolerance = 1e-12)
  expect_equal(out$summary$mean_recovery, mean(out$results$recovery))
  expect_false(any(grepl("acceptable", names(out$results), ignore.case = TRUE)))

  calculated <- spike_recovery(
    d, "sample", "result", "condition", "analyte", "solvent",
    sample_volume = 1, stock_concentration = 40, spike_volume = 0.1
  )
  expect_equal(calculated$results$recovery, out$results$recovery)
})

test_that("hook concentration is the preceding measured level without interpolation", {
  concentrations <- c(500, 1000, 3000, 5000, 6000, 7000, 8000)
  d <- rbind(
    data.frame(series = "A", concentration = concentrations,
               response = c(1, 1.5, 3, 2, 1.5, 0.9, 0.5)),
    data.frame(series = "B", concentration = concentrations,
               response = c(1, 1.6, 3.2, 2.5, 1.8, 1.2, 0.8))
  )
  out <- hook_effect(d, "concentration", "response", series = "series",
                     uloq = 500)

  expect_s3_class(out, "hook_effect")
  expect_equal(out$results$hook_concentration, c(6000, 7000))
  expect_equal(out$results$first_in_range_concentration, c(7000, 8000))
  expect_equal(out$results$hook_response, c(1.5, 1.2))
  expect_equal(out$results$first_in_range_response, c(0.9, 0.8))
  expect_equal(out$minimum_hook_concentration, 6000)
  expect_false("downturn_concentration" %in% names(out$results))
})

test_that("hook experiments without a threshold crossing report a lower bound", {
  d <- data.frame(concentration = c(500, 1000, 2000, 4000),
                  response = c(1, 2, 3, 4))
  out <- hook_effect(d, "concentration", "response", uloq = 500)
  expect_true(is.na(out$results$hook_concentration))
  expect_true(is.na(out$results$first_in_range_concentration))
  expect_equal(out$results$hook_lower_bound, 4000)
  expect_equal(out$results$highest_tested_concentration, 4000)

  supplied <- hook_effect(
    subset(d, concentration != 500), "concentration", "response",
    uloq = 500, uloq_response = 1
  )
  expect_equal(supplied$results$hook_lower_bound, 4000)
  expect_error(
    hook_effect(subset(d, concentration != 500), "concentration", "response",
                uloq = 500),
    "no response measured"
  )
})

test_that("EP34 plot methods return ggplot objects", {
  dilution_data <- data.frame(
    sample = rep(c("A", "B"), each = 4),
    factor = rep(rep(c(1, 2), each = 2), 2),
    result = c(99, 101, 49, 51, 198, 202, 99, 101)
  )
  dilution <- dilution_recovery(dilution_data, "sample", "result", "factor")

  spike_data <- data.frame(
    sample = rep(c("A", "B"), each = 4),
    condition = rep(rep(c("solvent", "analyte"), each = 2), 2),
    result = c(10, 10.2, 13, 13.2, 20, 20.2, 23, 23.2)
  )
  spike <- spike_recovery(spike_data, "sample", "result", "condition",
                          "analyte", "solvent", added_concentration = 3)

  hook_data <- data.frame(concentration = c(500, 1000, 2000, 3000),
                          response = c(1, 3, 2, 0.8))
  hook <- hook_effect(hook_data, "concentration", "response", uloq = 500)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_centered_plot(plot(dilution, "recovery"))
  expect_centered_plot(plot(dilution, "summary"))
  expect_centered_plot(plot(dilution, "observed"))
  expect_centered_plot(plot(spike))
  expect_centered_plot(plot(hook))
})

test_that("EP34 objects print invisibly and validate experimental axes", {
  dilution_data <- data.frame(sample = rep("A", 4),
                              factor = rep(c(1, 2), each = 2),
                              result = c(99, 101, 49, 51))
  dilution <- dilution_recovery(dilution_data, "sample", "result", "factor")
  spike_data <- data.frame(sample = rep("A", 4),
                           condition = rep(c("solvent", "analyte"), each = 2),
                           result = c(10, 10.2, 13, 13.2))
  spike <- spike_recovery(spike_data, "sample", "result", "condition",
                          "analyte", "solvent", added_concentration = 3)
  hook <- hook_effect(
    data.frame(concentration = c(1, 2, 3), response = c(1, 2, 0.5)),
    "concentration", "response", uloq = 1
  )

  for (object in list(dilution, spike, hook)) {
    capture.output(printed <- withVisible(print(object)))
    expect_false(printed$visible)
    expect_identical(printed$value, object)
  }

  dilution_output <- paste(capture.output(print(dilution)), collapse = "\n")
  expect_match(dilution_output, "Sample results", fixed = TRUE)
  expect_match(dilution_output, "Sample: A | Dilution factor: 1", fixed = TRUE)
  expect_match(dilution_output, "Corrected concentration", fixed = TRUE)

  spike_output <- paste(capture.output(print(spike)), collapse = "\n")
  expect_match(spike_output, "Sample results", fixed = TRUE)
  expect_match(spike_output, "Sample: A", fixed = TRUE)
  expect_match(spike_output, "Observed increase", fixed = TRUE)

  hook_output <- paste(capture.output(print(hook)), collapse = "\n")
  expect_match(hook_output, "Concentration-level results", fixed = TRUE)
  expect_match(hook_output, "Series: all | Concentration: 1", fixed = TRUE)
  expect_match(hook_output, "Mean response", fixed = TRUE)

  filtered_dilution <- dilution_recovery(
    transform(dilution_data, factor = c(1, 1, 0, 2)),
    "sample", "result", "factor"
  )
  expect_identical(filtered_dilution$excluded, 3L)
  expect_error(spike_recovery(spike_data, "sample", "result", "condition",
                              "analyte", "analyte", added_concentration = 3),
               "distinct")
  filtered_hook <- hook_effect(
    transform(data.frame(concentration = c(1, 2), response = c(1, 2)),
              concentration = c(1, NA)),
    "concentration", "response", uloq = 1
  )
  expect_identical(filtered_hook$excluded, 2L)
})
