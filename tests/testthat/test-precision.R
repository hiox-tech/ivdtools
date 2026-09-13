test_that("precision variance components agree with direct VCA fits", {
  d <- balanced_precision_data()
  p <- precision(d, y ~ day/run, by = "sample")
  p <- quietly_value(variance(p))

  for (sample_name in c("low", "high")) {
    sample_data <- d[d$sample == sample_name, ]
    direct <- suppressMessages(VCA::anovaVCA(y ~ day/run, Data = sample_data))
    direct_vc <- direct$aov.tab[, "VC"]
    names(direct_vc) <- rownames(direct$aov.tab)
    direct_vc <- direct_vc[names(p$results[[sample_name]]$vc)]
    expect_equal(
      p$results[[sample_name]]$vc,
      direct_vc,
      tolerance = 1e-10,
      info = sample_name
    )
  }
})

test_that("precision SD, CV, and percentages are internally coherent", {
  d <- balanced_precision_data()
  p <- quietly_value(variance(precision(d, y ~ day/run, by = "sample")))

  for (result in p$results) {
    expect_equal(result$sd_comp, sqrt(pmax(result$vc, 0)), tolerance = 1e-12)
    expect_equal(
      result$cv_comp,
      100 * result$sd_comp / result$mean,
      tolerance = 1e-12
    )
    expect_equal(
      sum(result$vc / sum(result$vc) * 100),
      100,
      tolerance = 1e-10
    )
  }
})

test_that("variance alias and confidence intervals preserve precision objects", {
  d <- balanced_precision_data()
  p <- precision(d, y ~ day/run, by = "sample")
  by_variance <- quietly_value(variance(p))
  by_alias <- quietly_value(vc(p))
  expect_equal(by_alias$results, by_variance$results)

  with_ci <- quietly_value(ci(by_variance))
  expect_s3_class(with_ci, "precision")
  expect_s3_class(with_ci$ci, "precision_ci")
  expect_true(all(as.numeric(with_ci$ci$SD$lower) >= 0))
  expect_true(all(
    as.numeric(with_ci$ci$SD$lower) <= as.numeric(with_ci$ci$SD$upper)
  ))
})

test_that("precision is invariant to row order", {
  d <- balanced_precision_data()
  p1 <- quietly_value(variance(precision(d, y ~ day/run, by = "sample")))
  p2 <- quietly_value(variance(precision(
    d[rev(seq_len(nrow(d))), ],
    y ~ day/run,
    by = "sample"
  )))
  expect_equal(p1$results$low$vc, p2$results$low$vc, tolerance = 1e-10)
  expect_equal(p1$results$high$vc, p2$results$high$vc, tolerance = 1e-10)
})

test_that("precision rejects invalid designs and incomplete prerequisites", {
  d <- balanced_precision_data()
  expect_error(precision(list(), y ~ day/run), "is.data.frame")
  expect_error(precision(d, y ~ missing/run), "not found|undefined")
  p <- precision(d, y ~ day/run, by = "sample")
  expect_output(profile(p, model.no = 1), "Run variance")
})

test_that("variance method argument selects ANOVA vs REML", {
  d <- balanced_precision_data()   # balanced design

  # auto (balanced by VCA::isBalanced) -> anova
  p <- quietly_value(variance(precision(d, y ~ day/run, by = "sample")))
  expect_equal(unique(vapply(p$results, `[[`, "", "method")), "anova")
  expect_true("Method" %in% names(p$vc))

  # missing cell (unbalanced nested) -> auto picks reml for that group only
  dm <- d[!(d$sample == "low" & d$day == 2 & d$run == 2), , drop = FALSE]
  pm <- quietly_value(variance(precision(dm, y ~ day/run, by = "sample")))
  expect_equal(pm$results$low$method,  "reml")
  expect_equal(pm$results$high$method, "anova")

  # explicit reml on balanced data
  pr <- quietly_value(variance(precision(d, y ~ run/rep, by = "sample"), method = "reml"))
  expect_equal(unique(vapply(pr$results, `[[`, "", "method")), "reml")
  expect_true(all(is.finite(as.numeric(pr$vc$VC))))

  # explicitly anova
  pa <- quietly_value(variance(precision(d, y ~ day/run/rep, by = "sample"), method = "anova"))
  expect_equal(unique(vapply(pa$results, `[[`, "", "method")), "anova")

  # auto (unequal cell sizes) -> reml
  ub <- data.frame(
    sample = "s1",
    day    = c("D1","D1","D2","D2"),
    run    = c("R1","R2","R1","R2"),
    stringsAsFactors = FALSE
  )
  ub <- ub[rep(seq_len(4), c(5, 1, 3, 1)), , drop = FALSE]
  set.seed(20260830)
  ub$y <- rnorm(nrow(ub), 50, 3)
  pu <- quietly_value(variance(precision(ub, y ~ day/run, by = "sample")))
  expect_equal(unique(vapply(pu$results, `[[`, "", "method")), "reml")
  expect_true(all(is.finite(as.numeric(pu$vc$VC))))
})

test_that("precision report requires variance and supports simple, detailed and CI tables", {
  d <- balanced_precision_data()
  p0 <- precision(d, y ~ day/run, by = "sample")
  expect_error(report(p0, day = "day", run = "run"), "Run variance")
  p <- quietly_value(variance(p0))
  simple <- report(p, day = "day", run = "run", table = "simple")
  detailed <- report(p, day = "day", run = "run", table = "detailed")
  expect_s3_class(simple, "data.frame")
  expect_true(nrow(simple) == 2L)
  expect_true(all(c("Repeatability SD", "Within-Laboratory Precision SD") %in% names(simple)))
  expect_true(all(c("Between-Run SD", "Between-Day SD") %in% names(detailed)))
  expect_error(report(p, day = "day", run = "run", table = "CI"), "Run ci")

  with_ci <- quietly_value(ci(p))
  ci_report <- report(
    with_ci, day = "day", run = "run", table = "CI"
  )
  expected_ci_columns <- c(
    "Repeatability SD Lower", "Repeatability SD Upper",
    "Repeatability CV Lower", "Repeatability CV Upper",
    "Within-Laboratory Precision SD Lower",
    "Within-Laboratory Precision SD Upper",
    "Within-Laboratory Precision CV Lower",
    "Within-Laboratory Precision CV Upper"
  )
  expect_true(all(expected_ci_columns %in% names(ci_report)))
  expect_equal(
    ci_report$`Within-Laboratory Precision SD Lower`,
    ci_report$`Reproducibility SD Lower`
  )
})

test_that("precision report supports EP05-A3 multisite nesting and combined CIs", {
  data("CA19_9", package = "VCA", envir = environment())
  p <- precision(CA19_9, result ~ site/day, by = "sample")
  p <- quietly_value(variance(p, NegVC = FALSE))
  p <- quietly_value(ci(p))

  expect_identical(
    .precision_term_key(c("site:day", "day:site")),
    c("day:site", "day:site")
  )
  detailed <- report(p, day = "day", site = "site", table = "detailed")
  ci_report <- report(p, day = "day", site = "site", table = "CI")

  expect_equal(nrow(detailed), 6L)
  expect_equal(
    detailed$`Within-Laboratory Precision SD`,
    c(0.8382, 1.3259, 7.7558, 1.4433, 3.1111, 8.7738),
    tolerance = 1e-4
  )
  expect_equal(
    detailed$`Reproducibility SD`,
    c(1.0425, 1.8376, 9.2228, 2.2929, 6.3050, 15.5271),
    tolerance = 1e-4
  )
  expect_true(all(c(
    "Within-Laboratory Precision SD Lower",
    "Within-Laboratory Precision SD Upper",
    "Within-Laboratory Precision CV Lower",
    "Within-Laboratory Precision CV Upper",
    "Reproducibility SD Lower", "Reproducibility SD Upper"
  ) %in% names(ci_report)))
  expect_equal(
    ci_report$`Within-Laboratory Precision SD Lower`,
    c(0.7029, 1.1358, 6.6510, 1.2107, 2.6317, 7.5297),
    tolerance = 1e-4
  )
  expect_equal(
    ci_report$`Within-Laboratory Precision SD Upper`,
    c(1.0385, 1.5931, 9.3043, 1.7873, 3.8057, 10.5144),
    tolerance = 1e-4
  )
  expect_false(isTRUE(all.equal(
    ci_report$`Within-Laboratory Precision SD Lower`,
    ci_report$`Reproducibility SD Lower`
  )))

  first <- p$results[["P1"]]
  total_ci <- .precision_combined_interval(
    first, names(first$vc), alpha = 0.05
  )
  inference <- suppressMessages(VCA::VCAinference(
    first$vca, VarVC = TRUE, alpha = 0.05
  ))
  vca_total <- inference$ConfInt$SD$TwoSided["total", c("LCL", "UCL")]
  expect_equal(unname(total_ci[c("lower", "upper")]),
               unname(as.numeric(vca_total)),
               tolerance = 1e-10)
})

test_that("precision report validates original RHS variables", {
  d <- balanced_precision_data(); d$lot <- "L1"
  p <- precision(d, y ~ lot/day/run, by = "sample")
  p$results <- list()
  expect_error(report(p, day = "day", run = "run"), "Unmapped design variables.*lot")
})
