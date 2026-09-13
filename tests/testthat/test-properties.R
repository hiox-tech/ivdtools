test_that("ROC AUC equals the rank definition across deterministic random cases", {
  for (seed in 1:12) {
    set.seed(seed)
    truth <- rep(c(0L, 1L), each = 15)
    marker <- round(stats::rnorm(30, mean = truth * 0.8), 1)
    d <- data.frame(marker, truth)
    object <- roc(d, "marker", "truth")
    signed_marker <- if (object$col_direction[["marker"]] == "leq") -marker else marker
    expect_equal(
      quietly_value(auc(object))$auc_list$marker$auc,
      manual_auc(signed_marker, truth),
      tolerance = 1e-12,
      info = paste("seed", seed)
    )
  }
})

test_that("qualitative raw and count interfaces are equivalent across tables", {
  set.seed(20260726)
  for (i in seq_len(10)) {
    counts <- sample(1:20, 4, replace = TRUE)
    raw <- data.frame(
      reference = rep(c("1", "1", "0", "0"), counts),
      test = rep(c("1", "0", "1", "0"), counts)
    )
    by_raw <- raw_to_table(raw, "test", "reference", positive = "1")
    by_count <- counts_to_table(
      tp = counts[1], fn = counts[2], fp = counts[3], tn = counts[4]
    )
    expect_equal(unname(by_raw$table), unname(by_count$table),
                 info = paste("iteration", i))
    expect_equal(
      quietly_value(diagnostics(by_raw))$diagnostics,
      quietly_value(diagnostics(by_count))$diagnostics,
      info = paste("iteration", i)
    )
  }
})

test_that("method comparison location shifts alter only the intercept", {
  d <- data.frame(id = 1:20, reference = 1:20,
                  candidate = 3 + 1.2 * (1:20))
  base <- suppressWarnings(quietly_value(
    regression(mcr(d, "id", "candidate", "reference"), "ols")
  ))
  shifted <- transform(d, candidate = candidate + 10)
  moved <- suppressWarnings(quietly_value(
    regression(mcr(shifted, "id", "candidate", "reference"), "ols")
  ))
  expect_equal(moved$regression$slope, base$regression$slope, tolerance = 1e-12)
  expect_equal(moved$regression$intercept - base$regression$intercept,
               10, tolerance = 1e-12)
})
