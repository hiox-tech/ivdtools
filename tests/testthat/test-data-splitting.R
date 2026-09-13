test_that("split_by_factors preserves factor order and data frames", {
  d <- data.frame(
    batch = factor(c("Batch-B", "Batch-A", "Batch-B", "Batch-C"),
                   levels = c("Batch-A", "Batch-B", "Batch-C", "Unused")),
    result = c(10.2, 9.9, 10.4, 10.1),
    run = 1:4
  )
  parts <- split_by_factors(d, "batch")

  expect_s3_class(parts, "factor_split")
  expect_named(parts, c("Batch-A", "Batch-B", "Batch-C"))
  expect_true(all(vapply(parts, is.data.frame, logical(1))))
  expect_equal(parts[["Batch-B"]]$run, c(1L, 3L))
  expect_true("batch" %in% names(parts[["Batch-A"]]))
  expect_equal(attr(parts, "variable"), "batch")
  expect_equal(attr(parts, "group_info")$n, c(1L, 2L, 1L))
  expect_equal(attr(parts, "n_excluded"), 0L)

  printed <- paste(capture.output(print(parts)), collapse = "\n")
  expect_match(printed, "Data split by factors")
  expect_match(printed, "Variable: batch")
  expect_match(printed, "Groups: 3")
  expect_match(printed, "Batch-B.*n = 2")
})

test_that("split_by_factors handles unused levels and missing groups", {
  d <- data.frame(
    batch = factor(c("A", NA, "B", "A"), levels = c("A", "B", "C")),
    result = 1:4
  )

  parts <- split_by_factors(d, "batch")
  expect_named(parts, c("A", "B"))
  expect_equal(attr(parts, "n_excluded"), 1L)
  expect_equal(sum(vapply(parts, nrow, integer(1))), 3L)
  expect_match(paste(capture.output(print(parts)), collapse = "\n"),
               "Missing/excluded: 1")

  with_empty <- split_by_factors(d, "batch", drop_empty = FALSE)
  expect_named(with_empty, c("A", "B", "C"))
  expect_equal(nrow(with_empty[["C"]]), 0L)
  expect_equal(attr(with_empty, "group_info")$n, c(2L, 1L, 0L))

  expect_error(split_by_factors(d, "batch", na_action = "error"),
               "missing value")
  expect_error(split_by_factors(d, "batch", drop_empty = NA),
               "TRUE or FALSE")
  expect_error(
    split_by_factors(data.frame(batch = c(NA_character_, NA_character_)),
                     "batch"),
    "No non-missing groups"
  )
})

test_that("split_by_factors uses first appearance for character and logical data", {
  character_data <- data.frame(
    sample = c("Sample B", "Sample A", "Sample B", "Sample C"),
    result = 1:4
  )
  character_parts <- split_by_factors(character_data, "sample")
  expect_named(character_parts, c("Sample B", "Sample A", "Sample C"))
  expect_equal(character_parts[["Sample B"]]$result, c(1L, 3L))

  logical_data <- data.frame(flag = c(TRUE, FALSE, TRUE), result = 1:3)
  logical_parts <- split_by_factors(logical_data, "flag")
  expect_named(logical_parts, c("TRUE", "FALSE"))
})

test_that("split_by_factors validates the grouping variable", {
  d <- data.frame(batch = c("A", "B"), result = 1:2)
  expect_error(split_by_factors(as.list(d), "batch"), "data frame")
  expect_error(split_by_factors(d, "absent"), "name one column")
  expect_error(split_by_factors(d, "result"),
               "factor, character, or logical")
})

test_that("factor-split data frames support identical downstream workflows", {
  set.seed(20260904)
  d <- data.frame(
    analyte = rep(c("A", "B"), each = 8),
    result = c(stats::rnorm(8, 10, 1), stats::rnorm(8, 20, 2))
  )
  parts <- split_by_factors(d, "analyte")
  results <- lapply(parts, function(part) {
    normal_test(part, "result", method = "shapiro")
  })

  expect_named(results, c("A", "B"))
  expect_true(all(vapply(results, inherits, logical(1), "normal_test")))
  expect_equal(vapply(results, `[[`, integer(1), "n"), c(A = 8L, B = 8L))
})
