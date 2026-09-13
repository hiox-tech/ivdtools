c5c95_test_data <- function() {
  means <- seq(36, 44, length.out = 7)
  residual <- seq(-1, 1, length.out = 12)
  data.frame(
    level = rep(paste0("L", 1:7), each = 12),
    response = unlist(lapply(means, function(mu)
      mu + (0.5 + 0.01 * mu) * residual))
  )
}

test_that("C5 and C95 roots achieve requested probabilities", {
  out <- c5_c95(c5c95_test_data(), "response", "level",
                cutoff = 40, model.no = 1)

  expect_s3_class(out, "c5_c95")
  expect_equal(out$estimates$label, c("C5", "C95"))
  expect_equal(out$estimates$actual_probability,
               out$estimates$probability, tolerance = 1e-5)
  expect_lt(out$estimates$mean_response[1], 40)
  expect_gt(out$estimates$mean_response[2], 40)
  expect_true(all(out$estimates$SD > 0))
})

test_that("C5/C95 supports directions, cutoffs, probabilities, and precision input", {
  d <- c5c95_test_data()
  geq <- c5_c95(d, "response", "level", cutoff = c(39, 41),
                probabilities = c(0.1, 0.9), model.no = 1)
  leq <- c5_c95(d, "response", "level", cutoff = 40,
                direction = "leq", model.no = 1)
  expect_equal(nrow(geq$estimates), 4)
  expect_gt(leq$estimates$mean_response[1], 40)
  expect_lt(leq$estimates$mean_response[2], 40)

  precision_object <- structure(
    list(profile = list(fits = list(total = geq$profile)),
         response = "response", by = "level"),
    class = "precision"
  )
  reused <- c5_c95(precision_object, cutoff = 40,
                   component = "total", probabilities = 0.5)
  expect_identical(reused$source, "precision object")
  expect_equal(reused$estimates$mean_response, 40, tolerance = 0.1)
})

test_that("C5/C95 reports bounded searches without a root", {
  out <- c5_c95(c5c95_test_data(), "response", "level", cutoff = 40,
                model.no = 1, search_range = c(39.9, 40.1),
                probabilities = 0.999)
  expect_true(is.na(out$estimates$mean_response))
  expect_equal(out$estimates$roots_found, 0)
  expect_true(any(grepl("no root", out$issues)))
})

test_that("C5/C95 validates profiles and returns reusable plots and printing", {
  d <- c5c95_test_data()
  out <- c5_c95(d, "response", "level", cutoff = 40, model.no = 1)
  capture.output(printed <- withVisible(print(out)))
  expect_false(printed$visible)
  expect_identical(printed$value, out)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_centered_plot(plot(out, type = "probability"))
  expect_centered_plot(plot(out, type = "profile"))

  expect_error(c5_c95(d, "response", "level", cutoff = Inf), "finite numeric")
  expect_error(c5_c95(d, "response", "level", cutoff = 40,
                      probabilities = c(0, 0.95)), "strictly between")
  expect_error(c5_c95(d, "response", "level", cutoff = 40,
                      model.no = 11), "integers from 1 to 10")
  expect_error(c5_c95(d, "response", "level", cutoff = 40,
                      search_range = c(1, 1)), "distinct")
  expect_error(c5_c95(structure(list(), class = "precision"), cutoff = 40),
               "no fitted precision profile")
})
