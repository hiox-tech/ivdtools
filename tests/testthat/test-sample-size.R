test_that("Bland-Altman function only returns a sample-size design", {
  result <- sample_size_bland_altman(
    power = 0.8, mu = 0.2, sd = 1, delta = 2.5
  )

  expect_s3_class(result, "sample_size_bland_altman")
  expect_identical(result$n, 201L)
  expect_equal(result$target_power, 0.8)
  expect_equal(result$achieved_power, 0.800300, tolerance = 1e-6)
  expect_lt(
    .ss_bland_altman_power(200, 0.2, 1, 2.5, 0.95, 0.95),
    0.8
  )
})

test_that("Bland-Altman validates population agreement feasibility", {
  expect_error(
    sample_size_bland_altman(power = 0.8, mu = 0, sd = 1, delta = 1.9),
    "cannot be demonstrated"
  )
  expect_error(
    sample_size_bland_altman(power = 0.8, mu = 0, sd = 1, delta = 2,
                             method = "other"),
    "method"
  )
})

test_that("PASS Bland-Altman sample-size matrix is transparently reproduced", {
  scenarios <- expand.grid(
    sd = c(2.5, 2.6, 2.7),
    power = c(0.8, 0.9)
  )
  package_n <- vapply(
    seq_len(nrow(scenarios)),
    function(i) sample_size_bland_altman(
      power = scenarios$power[i], mu = 0.5, sd = scenarios$sd[i],
      delta = 7, conf.level = 0.95, agree.level = 0.95
    )$n,
    integer(1)
  )

  # PASS Example 1 reports 60 rather than 59 in the first scenario. The
  # remaining five scenarios agree exactly with the Lu formula search.
  expect_identical(package_n, c(59L, 82L, 118L, 78L, 108L, 156L))
  expect_identical(c(60L, 82L, 118L, 78L, 108L, 156L) - package_n,
                   c(1L, 0L, 0L, 0L, 0L, 0L))
})

test_that("proportion CI function returns the first half-width design", {
  result <- sample_size_proportion_ci(
    p = 0.3, half_width = 0.05, method = "wilson"
  )

  expect_s3_class(result, "sample_size_proportion_ci")
  expect_identical(result$n, 320L)
  expect_lte(result$achieved_half_width, 0.05)
  previous <- .ss_proportion_ci_precision(
    result$n - 1L, 0.3, 0.95, "two.sided", "wilson"
  )
  expect_gt(previous$value, 0.05)
  expect_equal(result$actual_alpha, 1 - result$actual_coverage)
  expect_equal(result$achieved_width, 2 * result$achieved_half_width)
})

test_that("PASS exact proportion-CI sample-size matrix is reproduced", {
  scenarios <- expand.grid(
    p = seq(0.1, 0.5, by = 0.1),
    width = c(0.04, 0.06, 0.10)
  )
  pass_n <- c(
    914L, 1585L, 2065L, 2353L, 2449L,
    417L, 715L, 928L, 1056L, 1098L,
    158L, 264L, 341L, 387L, 402L
  )
  package_n <- vapply(
    seq_len(nrow(scenarios)),
    function(i) sample_size_proportion_ci(
      p = scenarios$p[i], half_width = scenarios$width[i] / 2,
      conf.level = 0.95, interval = "two.sided",
      method = "clopper-pearson"
    )$n,
    integer(1)
  )

  expect_identical(package_n, pass_n)
})

test_that("PASS one-sided, continuity-corrected, and boundary CI cases match", {
  pass_2 <- sample_size_proportion_ci(
    p = 0.92, half_width = 0.15104, conf.level = 0.95,
    interval = "lower", method = "clopper-pearson"
  )
  pass_3 <- sample_size_proportion_ci(
    p = 0.034483, half_width = 0.1945 / 2, conf.level = 0.95,
    interval = "two.sided", method = "wilson-cc"
  )
  pass_4 <- sample_size_proportion_ci(
    p = 0, half_width = 0.01, conf.level = 0.95,
    interval = "upper", method = "clopper-pearson"
  )

  expect_identical(c(pass_2$n, pass_3$n, pass_4$n), c(25L, 29L, 299L))
  expect_equal(pass_3$achieved_width, 0.19448, tolerance = 1e-5)
})

test_that("proportion CI methods and interval directions are supported", {
  methods <- c(
    "wald", "wald-cc", "wilson", "wilson-cc",
    "agresti-coull", "jeffreys", "clopper-pearson"
  )
  for (method in methods) {
    result <- sample_size_proportion_ci(
      p = 0.3, half_width = 0.08, method = method
    )
    expect_true(result$achieved_half_width <= 0.08)
    expect_true(result$actual_coverage >= 0 && result$actual_coverage <= 1)
  }

  lower <- sample_size_proportion_ci(
    p = 0.3, half_width = 0.05, interval = "lower"
  )
  upper <- sample_size_proportion_ci(
    p = 0.3, half_width = 0.05, interval = "upper"
  )
  expect_equal(lower$upper, 1)
  expect_equal(upper$lower, 0)
})

test_that("CMDE normal score calculations match PASS results", {
  positive <- sample_size_proportion(
    p0 = 0.85, p1 = 0.90, alpha = 0.05, power = 0.80,
    alternative = "two.sided", test = "score", power_method = "normal"
  )
  negative <- sample_size_proportion(
    p0 = 0.90, p1 = 0.94, alpha = 0.05, power = 0.80,
    alternative = "two.sided", test = "score", power_method = "normal"
  )

  expect_identical(positive$n, 363L)
  expect_identical(negative$n, 388L)
  expect_equal(positive$planning_power, 0.8002747, tolerance = 1e-7)
  expect_equal(negative$planning_power, 0.8000536, tolerance = 1e-7)
  expect_false(positive$alpha_met)
})

test_that("PASS normal-approximation validation examples are reproduced", {
  score <- sample_size_proportion(
    p0 = 0.5, p1 = 0.6, alpha = 0.05, power = 0.8,
    alternative = "greater", test = "score", power_method = "normal"
  )
  wald <- sample_size_proportion(
    p0 = 0.3, p1 = 0.5, alpha = 0.05, power = 0.8,
    alternative = "two.sided", test = "wald", power_method = "normal"
  )

  expect_identical(score$n, 153L)
  expect_equal(score$planning_power, 0.80125, tolerance = 1e-5)
  expect_identical(wald$n, 50L)
  expect_equal(wald$planning_power, 0.80743, tolerance = 1e-5)

  sas_exact <- sample_size_proportion(
    p0 = 0.3, p1 = 0.2, alpha = 0.05, power = 0.8,
    alternative = "less", test = "exact", power_method = "normal"
  )
  sas_cc <- sample_size_proportion(
    p0 = 0.3, p1 = 0.2, alpha = 0.05, power = 0.8,
    alternative = "less", test = "score-cc", power_method = "normal"
  )
  expect_identical(sas_exact$n, 119L)
  expect_identical(sas_cc$n, 129L)
  expect_equal(sas_cc$planning_power, 0.801, tolerance = 5e-4)
})

test_that("binomial designs control actual alpha and attain target power", {
  exact <- sample_size_proportion(
    p0 = 0.85, p1 = 0.90, alpha = 0.05, power = 0.80,
    alternative = "two.sided", test = "exact", power_method = "binomial"
  )
  wald <- sample_size_proportion(
    p0 = 0.85, p1 = 0.90, alpha = 0.05, power = 0.80,
    alternative = "two.sided", test = "wald", power_method = "binomial"
  )

  expect_identical(exact$n, 356L)
  expect_equal(round(exact$actual_alpha, 8), 0.04516485)
  expect_equal(exact$achieved_power, 0.8083882, tolerance = 1e-7)
  expect_identical(wald$n, 327L)
  expect_lte(wald$actual_alpha, wald$target_alpha + 1e-12)
  expect_gte(wald$achieved_power, wald$target_power)
  expect_true(wald$alpha_met)
})

test_that("PASS binomial-enumeration operating characteristics are reproduced", {
  ns <- c(10L, 11L, 12L, 25L, 50L, 70L)
  tests <- c("exact", "score", "score-cc", "wald", "wald-cc")
  pass_power <- rbind(
    c(0.04804, 0.03097, 0.08625, 0.15476, 0.23706, 0.36009),
    c(0.04804, 0.12484, 0.08625, 0.15476, 0.33613, 0.36009),
    c(0.04804, 0.03097, 0.08625, 0.15476, 0.23706, 0.36009),
    c(0.17958, 0.12484, 0.24060, 0.15476, 0.33613, 0.45495),
    c(0.17958, 0.12484, 0.08625, 0.15476, 0.23706, 0.36009)
  )
  pass_alpha <- rbind(
    c(0.0215, 0.0117, 0.0386, 0.0433, 0.0328, 0.0414),
    c(0.0215, 0.0654, 0.0386, 0.0433, 0.0649, 0.0414),
    c(0.0215, 0.0117, 0.0386, 0.0433, 0.0328, 0.0414),
    c(0.1094, 0.0654, 0.1460, 0.0433, 0.0649, 0.0722),
    c(0.1094, 0.0654, 0.0386, 0.0433, 0.0328, 0.0414)
  )

  for (i in seq_along(tests)) {
    observed <- lapply(
      ns, .ss_proportion_oc,
      p0 = 0.5, p1 = 0.6, alpha = 0.05,
      alternative = "two.sided", test = tests[i]
    )
    actual_power <- vapply(observed, `[[`, numeric(1), "actual_power")
    actual_alpha <- vapply(observed, `[[`, numeric(1), "actual_alpha")
    expect_equal(round(actual_power, 5), pass_power[i, ])
    expect_equal(round(actual_alpha, 4), pass_alpha[i, ])
  }
})

test_that("unsupported normal combinations and invalid directions fail", {
  expect_error(
    sample_size_proportion(
      p0 = 0.2, p1 = 0.35, test = "agresti-coull", power_method = "normal"
    ),
    "unavailable"
  )
  expect_error(
    sample_size_proportion(p0 = 0.4, p1 = 0.3, alternative = "greater"),
    "greater than p0"
  )
})

test_that("sample-size print methods expose achieved design quantities", {
  ba <- sample_size_bland_altman(power = 0.8, mu = 0.2, sd = 1, delta = 2.5)
  ci <- sample_size_proportion_ci(p = 0.3, half_width = 0.08)
  target <- sample_size_proportion(
    p0 = 0.2, p1 = 0.35, alternative = "greater"
  )

  expect_match(paste(capture.output(print(ba)), collapse = "\n"),
               "Achieved power")
  expect_match(paste(capture.output(print(ci)), collapse = "\n"),
               "Actual coverage")
  expect_match(paste(capture.output(print(target)), collapse = "\n"),
               "Actual alpha")
})
