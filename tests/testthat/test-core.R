test_that("package exports load and generics dispatch", {
  expected <- c(
    "mcr", "roc", "raw_to_table", "precision", "reference_interval",
    "fit_equation", "outliers_test", "normal_test", "qc_chart",
    "sample_size_proportion_ci", "stability_regression"
  )
  expect_true(all(vapply(expected, exists, logical(1), mode = "function")))
})

test_that("method comparison performs core analyses", {
  d <- data.frame(
    id = 1:12,
    candidate = c(10.1, 11.8, 15.2, 18.9, 21.1, 25.2,
                  29.8, 31.9, 35.4, 40.2, 44.8, 49.7),
    reference = c(10, 12, 15, 19, 21, 25, 30, 32, 35, 40, 45, 50)
  )
  x <- mcr(d, "id", "candidate", "reference")
  expect_s3_class(x, "mcr")
  expect_equal(x$n, 12)

  x <- correlation(x)
  expect_true(is.list(x$correlation))
  x <- regression(x, method = "ols")
  expect_true(is.list(x$regression))
  x <- bland_altman(x, plot = FALSE)
  expect_true(is.list(x$bland_altman))
})

test_that("qualitative and ROC workflows return S3 objects", {
  d <- data.frame(
    candidate = rep(c("positive", "negative"), 10),
    reference = rep(c("positive", "positive", "negative", "negative"), 5)
  )
  q <- raw_to_table(d, "candidate", "reference", positive = "positive")
  expect_s3_class(q, "fourfold_table")

  rdat <- data.frame(
    marker = c(1, 2, 3, 4, 5, 6, 8, 9),
    truth = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
  r <- roc(rdat, "marker", "truth")
  expect_s3_class(r, "roc")
  r <- auc(r)
  expect_s3_class(r$auc_table, "roc_auc_table")
})

test_that("standalone tools validate and calculate", {
  expect_error(reference_interval(data.frame(x = letters[1:5]), "x"))
  ri <- reference_interval(data.frame(x = seq_len(30)), "x")
  expect_s3_class(ri, "reference_interval")

  f <- fit_equation(
    "E01",
    data.frame(x = 1:5, y = c(2.1, 4, 6.2, 7.9, 10.1))
  )
  expect_s3_class(f, "fit_equation")
  expect_length(coef(f), 2)
  expect_length(residuals(f), 5)
})
