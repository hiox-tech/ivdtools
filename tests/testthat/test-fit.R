test_that("replicate summaries match base calculations and weighting formulas", {
  d <- data.frame(
    x = rep(1:3, each = 3),
    y = c(1, 2, 3, 2, 3, 4, 4, 5, 6)
  )
  observed <- replicate_to_mean(d, "x", "y", weights = "1/sd^2")
  expected_mean <- vapply(split(d$y, d$x), mean, numeric(1))
  expected_sd <- vapply(split(d$y, d$x), stats::sd, numeric(1))

  expect_num_equal(observed$y_mean, expected_mean)
  expect_num_equal(observed$y_sd, expected_sd)
  expect_num_equal(observed$weight, 1 / expected_sd^2)
  expect_equal(observed$y_n, rep(3L, 3))
})

test_that("replicate summaries support multiple method columns", {
  d <- data.frame(
    sample = rep(c("a", "b", "c"), each = 3),
    candidate = c(10, 11, 12, 20, 22, 24, 30, 33, 36),
    reference = c(9, 10, 11, 18, 20, 22, 27, 30, 33)
  )
  observed <- replicate_to_mean(
    d, "sample", c("candidate", "reference")
  )
  expect_named(observed, c(
    "x", "candidate_mean", "candidate_sd", "candidate_n",
    "reference_mean", "reference_sd", "reference_n"
  ))
  expect_equal(observed$candidate_mean,
               unname(vapply(split(d$candidate, d$sample), mean, numeric(1))))
  expect_equal(observed$reference_sd,
               unname(vapply(split(d$reference, d$sample), stats::sd, numeric(1))))
  expect_equal(observed$candidate_n, rep(3L, 3))
  expect_error(
    replicate_to_mean(d, "sample", c("candidate", "reference"),
                      weights = "n"),
    "single response"
  )
})

test_that("linear equation fitting agrees with stats::lm", {
  d <- data.frame(x = 1:8, y = c(3.1, 5, 7.2, 8.9, 11.1, 13, 15.2, 16.9))
  observed <- fit_equation("E01", d)
  expected <- stats::lm(y ~ x, data = d)

  expect_s3_class(observed, "fit_equation")
  expect_num_equal(coef(observed), stats::coef(expected), 1e-12)
  expect_num_equal(residuals(observed), stats::residuals(expected), 1e-12)
  expect_num_equal(
    quietly_value(predict(
      observed,
      newdata = data.frame(x = c(2.5, 4.5))
    ))$y,
    stats::predict(expected, data.frame(x = c(2.5, 4.5))),
    1e-12
  )
})

test_that("exact nonlinear data recover known parameters", {
  d <- data.frame(x = 0:6)
  d$y <- 10 * exp(-0.4 * d$x)
  fit <- fit_equation("E04", d, start = list(a = 9, b = 0.3))
  expect_equal(unname(coef(fit)["a"]), 10, tolerance = 1e-6)
  expect_equal(unname(coef(fit)["b"]), 0.4, tolerance = 1e-6)
})

test_that("box and equality constraints are satisfied", {
  d <- data.frame(x = 0:6)
  d$y <- 10 * exp(-0.4 * d$x)

  bounded <- fit_equation(
    "E04", d,
    start = list(a = 8, b = 0.2),
    lower = list(a = 0, b = 0),
    upper = list(a = 20, b = 1)
  )
  expect_true(all(coef(bounded) >= c(a = 0, b = 0)))
  expect_true(all(coef(bounded) <= c(a = 20, b = 1)))

  equal <- fit_equation(
    "E04", d,
    start = list(a = 8, b = 0.4),
    constraints = list(eq = function(p) p["b"] - 0.4)
  )
  expect_equal(unname(coef(equal)["b"]), 0.4, tolerance = 1e-5)
})

test_that("fit comparison ranks an exact generating model first", {
  d <- data.frame(x = 1:10)
  d$y <- 2 + 3 * d$x + c(0.10, -0.05, 0.02, -0.08, 0.04,
                          0.03, -0.06, 0.08, -0.02, -0.06)
  cmp <- compare_equation(d, eqs = c("E01", "E02", "E06"))
  expect_s3_class(cmp, "compare_equation")
  expect_equal(cmp$table$Eq[1], "E01")
  expect_true(all(diff(cmp$table$AIC) >= 0))
  printed <- paste(capture.output(print(cmp)), collapse = "\n")
  expect_true(grepl("Formula", printed, fixed = TRUE))
  expect_true(grepl(cmp$table$Formula[1], printed, fixed = TRUE))
})

test_that("fit functions validate columns, methods, and weights", {
  d <- data.frame(x = 1:4, y = 1:4)
  expect_error(replicate_to_mean(d, "missing", "y"), "[Cc]olumn")
  expect_error(fit_equation("unknown", d), "not found|No equation")
  expect_error(fit_equation("E01", d, weights = c(1, -1, 1, 1)), "weight")
})

test_that("list_equation does not change the console width", {
  old_width <- getOption("width")
  on.exit(options(width = old_width), add = TRUE)
  options(width = 173L)

  invisible(capture.output(list_equation()))
  expect_identical(getOption("width"), 173L)
})

test_that("list_equation does not change the console width after a printing error", {
  old_width <- getOption("width")
  on.exit(options(width = old_width), add = TRUE)
  options(width = 173L)
  testthat::local_mocked_bindings(
    .print_df = function(...) stop("simulated printing failure"),
    .package = "ivdtools"
  )

  expect_error(list_equation(), "simulated printing failure")
  expect_identical(getOption("width"), 173L)
})

test_that("Sadler variance models are removed from the fit API", {
  d <- data.frame(x = c(1, 2, 5, 10, 20, 50),
                  var = c(0.1, 0.3, 0.8, 3.2, 12.5, 78.4))
  # Sadler is no longer part of fit_equation() / list_equation()
  expect_error(fit_equation("sadler", d, "x", "var"), "not found|No equation")
  expect_error(fit_equation("V03", d, "x", "var"), "not found|No equation")
  tbl <- eq_registry()
  expect_false(any(tbl$category == "Sadler"))
  expect_false(any(grepl("^V", tbl$id)))
})
