test_that("Westgard rules flag exactly the intended observations", {
  z <- c(0, 2.1, 2.2, -2.1, rep(0.5, 10))
  q <- qc_chart(data.frame(z), "z", mean = 0, sd = 1)

  expect_equal(q$violations[["1-2s"]], 2:4)
  expect_length(q$violations[["1-3s"]], 0)
  expect_equal(q$violations[["2-2s"]], 2:3)
  expect_equal(q$violations[["R-4s"]], 3:4)
  expect_equal(q$violations[["10x"]], 5:14)
})

test_that("QC cleaning preserves run labels and reports missing values", {
  d <- data.frame(run = 11:15, lot = c("a", "a", "b", "b", "b"),
                  value = c(0, NA, 3.5, 0, 0))
  q <- qc_chart(d, "value", rules = "1-3s", mean = 0, sd = 1,
                group = "lot", run = "run")
  expect_equal(q$n_miss, 1L)
  expect_equal(q$plot_data$run, c(11L, 13L, 14L, 15L))
  expect_equal(q$violations[["1-3s"]], 2L)
  expect_equal(q$n_oc, 1L)
})

test_that("Youden statistics and two-SD classification are independent", {
  d <- data.frame(a = c(0, 1, 2, NA), b = c(0, 2, 5, 1))
  y <- youden_plot(d, "a", "b", mean1 = 1, mean2 = 1, sd1 = 1, sd2 = 1)
  expect_equal(y$correlation, stats::cor(c(0, 1, 2), c(0, 2, 5)))
  expect_equal(y$outside_rows, 3L)
  expect_equal(y$n_miss, 1L)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(quietly_value(plot(y)), y)
})

test_that("QC functions reject invalid rule and scale inputs", {
  expect_error(qc_chart(data.frame(x = 1:3), "x", rules = "bad"), "Unknown")
  expect_error(qc_chart(data.frame(x = 1:3), "x", mean = 2, sd = 0), "positive")
  expect_error(youden_plot(data.frame(x = 1:3), "x", "missing"), "column")
})
