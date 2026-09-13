quietly_value <- function(expr) {
  out <- NULL
  capture.output(out <- force(expr))
  out
}

expect_num_equal <- function(object, expected, tolerance = 1e-8) {
  testthat::expect_equal(
    as.numeric(object),
    as.numeric(expected),
    tolerance = tolerance
  )
}

expect_centered_plot <- function(object) {
  testthat::expect_s3_class(object, "ggplot")
  testthat::expect_true(
    is.character(object$labels$title) && length(object$labels$title) == 1L &&
      nzchar(object$labels$title)
  )
  testthat::expect_equal(object$theme$plot.title$hjust, 0.5)
  invisible(object)
}

manual_auc <- function(marker, truth) {
  pos <- marker[truth == 1]
  neg <- marker[truth == 0]
  comparisons <- outer(pos, neg, "-")
  (sum(comparisons > 0) + 0.5 * sum(comparisons == 0)) /
    (length(pos) * length(neg))
}

balanced_precision_data <- function() {
  design <- expand.grid(
    sample = c("low", "high"),
    day = factor(seq_len(3)),
    run = factor(seq_len(2)),
    rep = seq_len(2)
  )
  sample_mean <- ifelse(design$sample == "low", 10, 30)
  day_effect <- c(-0.3, 0, 0.3)[as.integer(design$day)]
  run_effect <- c(-0.1, 0.1)[as.integer(design$run)]
  rep_effect <- ifelse(design$rep == 1, -0.05, 0.05)
  design$y <- sample_mean + day_effect + run_effect + rep_effect
  design
}
