test_that("stability bias uses condition-specific baselines", {
  d <- data.frame(
    condition = rep(c("ref", "test"), each = 3),
    time = rep(0:2, 2),
    value = c(100, 98, 96, 110, 105, 99)
  )
  s <- stability_bias(d, "condition", "time", "value")
  expect_equal(s$summary$relative_bias, c(0, -2, -4, 0, -100 * 5 / 110,
                                          -100 * 11 / 110))
  expect_equal(s$between_bias$absolute_bias, c(10, 7, 3))
})

test_that("self stability regression agrees with lm on node means", {
  d <- expand.grid(time = 0:4, rep = 1:3)
  d$value <- 100 - 2 * d$time + c(-1, 0, 1)
  s <- stability_regression(d, "time", "value", mode = "self",
                            bias_type = "absolute", direction = "decrease",
                            limit = 8)
  direct <- stats::lm(value ~ time, aggregate(value ~ time, d, mean))
  expect_equal(unname(stats::coef(s$fit)), unname(stats::coef(direct)),
               tolerance = 1e-12)
  st <- stability_time(s, max_time = 20)
  expect_equal(st$result$time[st$result$target == "estimate"], 4,
               tolerance = 1e-6)
})

test_that("MKT agrees with direct Boltzmann averaging", {
  d <- data.frame(p1 = c(20, 25, 30), p2 = c(22, 27, 32))
  observed <- mkt(d, c("p1", "p2"), ea = 80)
  r <- 8.3145e-3
  direct <- function(x) -80 / (r * log(mean(exp(-80 / (r * (x + 273.15))))))
  expect_equal(observed$probes$mkt_K, c(direct(d$p1), direct(d$p2)),
               tolerance = 1e-10)
  expect_equal(observed$combined$mkt_K,
               -80 / (r * log(mean(exp(-80 / (r * (unlist(d) + 273.15)))))),
               tolerance = 1e-10)
})

test_that("Arrhenius recovers a known first-order activation model", {
  r <- 8.3145e-3
  ea <- 60
  temps <- c(25, 35, 45)
  rates <- 0.02 * exp(-ea / r * (1 / (temps + 273.15) - 1 / 298.15))
  d <- do.call(rbind, lapply(seq_along(temps), function(i) {
    time <- 0:4
    data.frame(temp = temps[i], time = time, value = 100 * exp(-rates[i] * time))
  }))
  observed <- suppressWarnings(arrhenius(
    d, "temp", "time", "value", target_temp = 25,
    order = "first", direction = "decrease", limit = 10
  ))
  expect_equal(unname(observed$activation_energy), ea, tolerance = 1e-8)
  expect_equal(observed$target_rate, 0.02, tolerance = 1e-8)
  expect_equal(observed$duration, -log(0.9) / 0.02, tolerance = 1e-8)
})

test_that("stability tools reject insufficient or invalid inputs", {
  expect_error(stability_regression(
    data.frame(t = 0:2, y = 3:1), "t", "y", mode = "self"
  ), "At least 4")
  expect_error(mkt(data.frame(t = c(1, 1), x = c(20, 21)), "x", time = "t"),
               "strictly increasing")
  expect_error(arrhenius(
    data.frame(temp = rep(c(25, 35), each = 3), time = rep(0:2, 2), value = 6:1),
    "temp", "time", "value", target_temp = 5, limit = 10
  ), "At least 3")
})
