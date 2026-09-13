test_that("raw and count constructors produce the same 2x2 table", {
  raw <- data.frame(
    candidate = c(rep("positive", 12), rep("negative", 23)),
    reference = c(rep("positive", 10), rep("negative", 2),
                  rep("positive", 3), rep("negative", 20))
  )
  from_raw <- raw_to_table(raw, "candidate", "reference", positive = "positive")
  from_default <- raw_to_table(raw, "candidate", "reference")
  from_counts <- counts_to_table(10, 2, 20, 3)
  expect_equal(unname(from_raw$table), unname(from_counts$table))
  expect_equal(unname(from_default$table), unname(from_counts$table))
  default_diagnostics <- quietly_value(diagnostics(from_default))$diagnostics
  expect_equal(
    unname(c(default_diagnostics$tp, default_diagnostics$fp,
             default_diagnostics$fn, default_diagnostics$tn)),
    c(10, 2, 3, 20)
  )
})

test_that("diagnostic point estimates agree with hand calculations", {
  q <- counts_to_table(tp = 45, fp = 3, tn = 97, fn = 5)
  d <- quietly_value(diagnostics(q))$diagnostics

  expect_equal(unname(d$sensitivity[1]), 45 / 50)
  expect_equal(unname(d$specificity[1]), 97 / 100)
  expect_equal(unname(d$ppv[1]), 45 / 48)
  expect_equal(unname(d$npv[1]), 97 / 102)
  expect_equal(unname(d$accuracy[1]), 142 / 150)
  expect_equal(unname(d$prevalence[1]), 50 / 150)
  expect_equal(unname(d$lr_positive[1]), (45 / 50) / (3 / 100))
  expect_equal(unname(d$lr_negative[1]), (5 / 50) / (97 / 100))
  expect_equal(unname(d$odds_ratio[1]), (45 * 97) / (3 * 5))
})

test_that("Cohen kappa and McNemar agree with independent formulas", {
  q <- counts_to_table(tp = 45, fp = 3, tn = 97, fn = 5)
  k <- quietly_value(kappa(q))$kappa
  observed_agreement <- (45 + 97) / 150
  expected_agreement <- ((48 * 50) + (102 * 100)) / 150^2
  expected_kappa <- (observed_agreement - expected_agreement) /
    (1 - expected_agreement)
  expect_equal(k$kappa, expected_kappa, tolerance = 1e-12)

  m <- quietly_value(mcnemar(q))$mcnemar
  expected_p <- stats::binom.test(min(3, 5), 8, p = 0.5)$p.value
  expect_equal(m$p_value, expected_p, tolerance = 1e-12)
})

test_that("diagnostic confidence interval methods remain bounded", {
  q <- counts_to_table(1, 0, 9, 0)
  methods <- c(
    "wald", "wald-cc", "wilson", "wilson-cc",
    "agresti-coull", "jeffreys", "clopper-pearson"
  )
  for (method in methods) {
    d <- quietly_value(diagnostics(q, ci.method = method))$diagnostics
    for (metric in c("sensitivity", "specificity", "accuracy", "prevalence")) {
      expect_true(all(d[[metric]][2:3] >= 0 & d[[metric]][2:3] <= 1),
                  info = paste(method, metric))
    }
  }
})

test_that("qualitative analysis rejects impossible counts and levels", {
  expect_error(counts_to_table(-1, 0, 1, 0), "non-negative")
  expect_error(counts_to_table(1.2, 0, 1, 0), "integer")
  bad <- data.frame(a = c("a", "b", "c"), b = c("a", "b", "c"))
  expect_error(raw_to_table(bad, "a", "b"), "exactly 2")
})

test_that("qualitative print methods return their objects invisibly", {
  q <- counts_to_table(tp = 45, fp = 3, tn = 97, fn = 5)
  q <- quietly_value(diagnostics(q))
  q <- quietly_value(kappa(q))
  q <- quietly_value(mcnemar(q))

  for (object in list(q, q$diagnostics, q$kappa, q$mcnemar)) {
    captured <- capture.output(visible <- withVisible(print(object)))
    expect_false(visible$visible)
    expect_identical(visible$value, object)
    expect_true(length(captured) > 0)
  }
})
