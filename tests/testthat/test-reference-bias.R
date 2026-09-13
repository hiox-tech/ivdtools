test_that("reference-material bias follows the YY/T uncertainty equations", {
  d <- data.frame(level = rep(c("low", "high"), each = 4),
                  result = c(9.8, 10.0, 10.2, NA, 19.8, 20.0, 20.2, 20.4))
  out <- reference_bias(d, "result", reference = c(10, 20),
                        reference_uncertainty = c(0.2, 0.3),
                        level = "level", k = 2)
  low <- d$result[d$level == "low" & is.finite(d$result)]
  expected_u <- sd(low) / sqrt(length(low))
  expected_half_width <- sqrt((2 * expected_u)^2 + 0.2^2)

  expect_s3_class(out, "reference_bias")
  expect_equal(out$result$mean[1], mean(low))
  expect_equal(out$result$standard_uncertainty[1], expected_u)
  expect_equal(out$result$interval_half_width[1], expected_half_width)
  expect_equal(out$result$lower[1], out$result$bias[1] - expected_half_width)
  expect_equal(out$result$upper[1], out$result$bias[1] + expected_half_width)
  expect_equal(out$result$excluded, c(1, 0))
})

test_that("reference inputs can be scalars, row vectors, or data columns", {
  d <- data.frame(level = rep(c("zero", "positive"), each = 3),
                  result = c(-0.1, 0, 0.1, 9.9, 10, 10.1),
                  assigned = rep(c(0, 10), each = 3),
                  assigned_u = rep(c(0.05, 0.2), each = 3))
  from_columns <- reference_bias(d, "result", assigned, assigned_u,
                                 level = "level")
  from_names <- reference_bias(d, "result", "assigned", "assigned_u",
                               level = "level")
  overall <- reference_bias(subset(d, level == "positive"), "result", 10, 0.2)

  expect_equal(from_columns$result, from_names$result)
  expect_true(is.na(from_columns$result$relative_bias[1]))
  expect_identical(overall$result$level, "Overall")
})

test_that("reference bias validates level constants and reusable printing", {
  d <- data.frame(level = rep(c("A", "B"), each = 3),
                  result = c(1, 2, 3, 4, 5, 6))
  expect_error(reference_bias(d, "result", c(1, 2, 1, 4, 4, 4), 0.1,
                              level = "level"), "constant within level")
  expect_error(reference_bias(d, "result", c(1, 2, 3), 0.1,
                              level = "level"), "length 1")
  expect_error(reference_bias(d, "result", c(2, 5), c(-0.1, 0.2),
                              level = "level"), "must not contain negative")
  expect_error(reference_bias(transform(d, level = replace(level, 1, NA)),
                              "result", c(2, 5), 0.1, level = "level"),
               "must not contain missing")
  expect_error(reference_bias(data.frame(result = c(1, NA)), "result", 1, 0.1),
               "at least two finite")

  out <- reference_bias(d, "result", c(2, 5), 0.1, level = "level")
  capture.output(printed <- withVisible(print(out)))
  expect_false(printed$visible)
  expect_identical(printed$value, out)
  expect_error(ivdtools:::print.reference_bias(list()), "reference_bias object")
})
