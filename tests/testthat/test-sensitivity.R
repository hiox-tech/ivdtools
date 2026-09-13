ep17_long <- function(values, lot, prefix) {
  values <- as.matrix(values)
  data.frame(
    lot = lot,
    sample = rep(paste(prefix, seq_len(ncol(values))), each = nrow(values)),
    result = as.vector(values),
    stringsAsFactors = FALSE
  )
}

ep17_appendix_a <- function() {
  blank_1 <- matrix(c(
     2.6, 1.0,-4.4, 1.5, 1.2, -0.8, 2.9,-3.4,-1.9,-0.7,
     5.5, 4.9, 7.0, 5.1, 6.1,  6.0, 8.0, 6.9, 5.7, 5.1,
     4.5, 6.9, 4.3, 4.1, 4.8,  0.6, 5.0, 3.2, 4.5, 3.3,
    -2.3, 3.4,-1.4,-0.6,-2.8,  3.4, 1.2,-4.2, 0.5,-1.4,
     5.9, 6.5, 5.9, 5.4, 8.7,  7.6, 5.6, 7.6, 7.6, 3.6,
     4.1,-2.2, 3.8, 4.4, 5.1, -1.4, 2.3, 5.8, 6.6, 3.5
  ), ncol = 5, byrow = TRUE)
  blank_2 <- matrix(c(
     4.6, 9.2, 6.1, 4.0, 4.0,  4.1, 8.3, 3.2,11.5, 6.2,
     1.6, 4.8, 3.9, 4.5,-0.2,  3.7, 5.4, 1.4, 3.6, 2.3,
     2.2, 4.8, 3.1, 4.4, 1.6,  0.7, 6.3, 4.1, 6.8, 2.6,
     4.6, 5.4, 1.0, 7.1, 6.4,  2.6, 9.6, 3.4, 4.2, 5.7,
     1.1, 7.7, 0.1, 3.7, 4.2, -4.4, 3.1, 0.4, 3.7, 3.7,
     0.9, 6.1, 2.9, 5.3, 1.4,  0.7,10.0,-1.6, 4.5, 1.5
  ), ncol = 5, byrow = TRUE)
  low_1 <- matrix(c(
    21.0,13.3,12.8,17.3,19.2, 22.8,12.6,12.9,19.2,22.7,
    28.2,18.2,17.4,21.5,28.3, 25.9,14.7,16.0,22.2,26.2,
    26.4,17.8,15.9,24.1,25.1, 28.3,14.0,14.1,25.8,30.3,
    20.7,14.1,11.3,16.0,23.4, 21.9,12.5, 9.4,16.4,19.2,
    24.7,11.3,10.6,24.9,26.3, 22.5,12.2,13.6,23.8,23.1,
    28.5,16.2,17.6,22.1,27.5, 29.2,13.9,14.9,26.1,30.1
  ), ncol = 5, byrow = TRUE)
  low_2 <- matrix(c(
    22.0,15.6,13.0,18.8,32.9, 22.5,21.2,15.9,17.6,30.4,
    21.8,14.8, 9.0,14.1,29.4, 22.1,14.9, 7.0,14.9,27.6,
    20.3,16.0,13.4,19.2,27.7, 21.0,15.8, 8.5,15.8,30.6,
    25.3,21.6,16.3,19.8,31.4, 26.0,22.8,18.1,21.4,30.4,
    27.2,15.3,12.4,18.0,32.5, 25.1,18.7,11.1,18.0,28.9,
    25.3,18.3,11.3,19.6,29.8, 25.3,19.5,10.1,23.1,35.1
  ), ncol = 5, byrow = TRUE)
  list(
    blank = rbind(ep17_long(blank_1, "1", "Blank"),
                  ep17_long(blank_2, "2", "Blank")),
    low = rbind(ep17_long(low_1, "1", "Low"),
                ep17_long(low_2, "2", "Low"))
  )
}

test_that("Appendix A LoB and classical LoD are reproduced", {
  a <- ep17_appendix_a()
  lob_np <- sensitivity_lob(a$blank, "result", "sample", "lot",
                            method = "nonparametric")
  expect_s3_class(lob_np, "sensitivity_lob")
  expect_equal(lob_np$by_lot$lob, c(7.6, 9.4), tolerance = 1e-12)
  expect_equal(lob_np$estimate, 9.4, tolerance = 1e-12)

  lod_np <- sensitivity_lod(a$low, "result", "sample", lob_np, "lot",
                            method = "classical")
  expect_s3_class(lod_np, "sensitivity_lod")
  expect_equal(round(lod_np$by_lot$lod, 1), c(14.5, 13.8))
  expect_equal(round(lod_np$estimate, 1), 14.5)

  lob_p <- sensitivity_lob(a$blank, "result", "sample", "lot",
                           method = "parametric")
  expect_equal(round(lob_p$by_lot$lob, 1), c(8.8, 8.6))
  expect_equal(round(lob_p$estimate, 1), 8.8)
  lod_p <- sensitivity_lod(a$low, "result", "sample", lob_p, "lot",
                           method = "classical")
  expect_equal(round(lod_p$estimate, 1), 13.9)
})

test_that("nonparametric trial LoD reports the median only at the beta target", {
  d <- data.frame(sample = rep(c("A", "B"), each = 10),
                  result = c(4, rep(8, 19)))
  ok <- sensitivity_lod(d, "result", "sample", lob = 5,
                        method = "nonparametric", beta = 0.051)
  expect_equal(ok$by_lot$proportion_below, 0.05)
  expect_equal(ok$estimate, 8)
  not_ok <- sensitivity_lod(d, "result", "sample", lob = 5,
                            method = "nonparametric", beta = 0.049)
  expect_true(is.na(not_ok$estimate))
})

ep17_profile_data <- function() {
  means <- rbind(
    data.frame(lot = "1", sample = LETTERS[1:6],
               mean = c(.69, 1.42, 2.65, 4.08, 6.08, 10.36),
               sd = c(.39, .39, .46, .55, .64, 1.12)),
    data.frame(lot = "2", sample = LETTERS[1:6],
               mean = c(.78, 1.73, 2.89, 3.82, 6.33, 10.92),
               sd = c(.29, .54, .55, .63, .82, 1.38))
  )
  z <- as.numeric(scale(seq_len(80), center = TRUE, scale = TRUE))
  do.call(rbind, lapply(seq_len(nrow(means)), function(i) {
    data.frame(lot = means$lot[i], sample = means$sample[i],
               result = means$mean[i] + means$sd[i] * z)
  }))
}

test_that("Appendix B quadratic precision-profile LoD is reproduced", {
  res <- sensitivity_lod(ep17_profile_data(), "result", "sample", lob = .51,
                         lot = "lot", method = "profile",
                         profile_model = "quadratic")
  expect_equal(res$estimate, 1.17, tolerance = 0.03)
  expect_true(all(abs(res$by_lot$lod - 1.17) < 0.03))
  expect_false(any(res$by_lot$extrapolated))
})

ep17_appendix_d <- function() {
  assigned <- c(38.2, 47.1, 44.7, 36.5, 42.8)
  lot_a <- matrix(c(
    36.7,49.9,46.1,33.3,42.9, 37.9,50.0,43.1,34.2,41.8,
    38.3,48.1,39.4,34.5,43.8, 36.8,47.8,47.3,43.1,46.3,
    33.5,43.9,45.8,34.0,43.3, 39.2,45.6,44.8,37.1,46.0,
    41.3,45.4,44.6,35.3,42.6, 37.9,51.5,47.3,32.4,41.4,
    34.9,45.8,38.9,36.0,42.8), ncol = 5, byrow = TRUE)
  lot_b <- matrix(c(
    38.5,45.8,46.7,35.5,42.0, 41.0,47.8,43.6,40.0,44.1,
    43.2,46.6,42.4,34.0,43.2, 36.8,46.9,46.5,32.9,46.6,
    42.1,51.3,47.9,33.1,45.5, 35.8,50.5,42.7,38.6,43.5,
    36.8,44.3,42.1,36.2,41.4, 44.1,47.5,43.4,41.4,48.2,
    39.5,52.4,44.7,33.0,45.7), ncol = 5, byrow = TRUE)
  make <- function(x, lot) {
    d <- ep17_long(x, lot, "Pool")
    d$assigned <- rep(assigned, each = nrow(x))
    d
  }
  rbind(make(lot_a, "A"), make(lot_b, "B"))
}

test_that("Appendix D Westgard total-error LoQ is reproduced", {
  res <- sensitivity_loq(ep17_appendix_d(), "result", "sample", "assigned",
                         lot = "lot", method = "westgard", goal = .216,
                         goal_scale = "relative")
  expect_s3_class(res, "sensitivity_loq")
  expect_equal(round(res$by_lot$candidate_loq, 1), c(35.5, 36.1))
  expect_equal(round(res$candidate_loq, 1), 36.1)
  constrained <- sensitivity_loq(ep17_appendix_d(), "result", "sample",
    "assigned", "lot", "westgard", .216, "relative", lod = 40)
  expect_equal(round(constrained$candidate_loq, 1), 36.1)
  expect_equal(constrained$reported_loq, 40)
})

test_that("RMS, precision, and custom LoQ definitions are available", {
  d <- ep17_appendix_d()
  rms <- sensitivity_loq(d, "result", "sample", "assigned", "lot",
                         "rms", .20, "relative")
  precision <- sensitivity_loq(d, "result", "sample", "assigned", "lot",
                               "precision", .10, "relative")
  custom <- sensitivity_loq(d, "result", "sample", "assigned", "lot",
    "custom", .25, "relative",
    accuracy_function = function(mean, sd, assigned, n) abs(mean-assigned)+sd)
  expect_s3_class(rms, "sensitivity_loq")
  expect_s3_class(precision, "sensitivity_loq")
  expect_s3_class(custom, "sensitivity_loq")
})

test_that("Appendix E LoB and LoD claim verification is reproduced", {
  blanks <- c(rep(0, 8), 1.08,1.92,2.38,2.98,3.80,4.78,7.30,8.81,
              10.31,11.29,13.48,14.39,16.97,17.40,18.01,22.65)
  positives <- c(18.80,19.02,26.63,26.91,31.08,33.99,35.11,35.90,
    41.67,43.90,46.32,47.77,47.99,48.83,54.67,57.30,59.10,61.17,
    61.96,62.97,66.44,73.44,73.80,75.71)
  lob <- sensitivity_verify(data.frame(sample = rep(c("B1", "B2"), each = 12),
    result = blanks), "result", "sample", "lob", claim = 20)
  lod <- sensitivity_verify(data.frame(sample = rep(c("L1", "L2"), each = 12),
    result = positives), "result", "sample", "lod", claim = 45, lob = 20)
  expect_equal(lob$result$n_consistent, 23)
  expect_equal(lod$result$n_consistent, 22)
  expect_equal(lob$result$ep17_boundary, .87)
  expect_equal(lod$result$ep17_boundary, .87)
  expect_false(any(c("pass", "verified") %in% names(lob$result)))
  expect_equal(lob$boundary_rows$N, c(20, 30))
})

test_that("Appendix F LoQ claim verification is reproduced", {
  target <- c(46.4,45.8,49.1,46.3,49.7)
  values <- matrix(c(
    47.8,47.3,49.7,50.4,53.7, 44.6,48.8,51.2,49.5,52.7,
    47.1,47.6,57.3,44.0,55.9, 50.8,54.7,54.6,51.5,55.1,
    48.2,50.7,49.3,51.4,55.5, 52.5,50.8,53.3,49.8,57.3,
    49.4,52.5,58.0,46.1,51.8, 52.0,50.4,49.5,45.7,48.8,
    46.3,49.6,52.2,50.9,51.7), ncol = 5, byrow = TRUE)
  d <- ep17_long(values, "1", "Sample")
  d$target <- rep(target, each = nrow(values))
  res <- sensitivity_verify(d, "result", "sample", "loq", claim = 50,
                            assigned = "target", allowable_error = .15)
  expect_equal(res$result$n_consistent, 41)
  expect_equal(res$result$n_total, 45)
  expect_equal(res$result$ep17_boundary, .88)
  expect_equal(res$boundary_rows$N, c(40, 50))
})

test_that("probit LoD supports summarized hit counts", {
  d <- data.frame(sample = paste0("L", 1:5), concentration = c(1,2,4,8,16),
                  hits = c(0,2,10,18,20), total = rep(20, 5), result = 1:5)
  res <- sensitivity_lod(d, "result", "sample", lob = 0,
    concentration = "concentration", detected = "hits", method = "probit",
    total = "total", concentration_scale = "log10")
  expect_true(is.finite(res$estimate))
  expect_true(res$estimate > 4)
  expect_true(res$estimate < 20)
  expect_true(all(c("deviance_p", "pearson_p") %in% names(res$by_group)))
})

test_that("study design reports requirements and observed counts without conclusions", {
  d <- data.frame(lot = rep(c("A", "B"), each = 12), day = rep(1:3, 8),
                  sample = rep(letters[1:4], 6), result = seq_len(24))
  design <- sensitivity_design(c("classical", "verification"), d,
    "result", "sample", "lot", "day", limit = "lob")
  expect_s3_class(design, "sensitivity_design")
  expect_equal(design$requirements$reagent_lots, c(2, 1))
  expect_equal(design$actual$observations, 24)
  expect_false(any(c("pass", "met", "verified") %in% names(design$actual)))
})

test_that("input validation and plotting methods work", {
  a <- ep17_appendix_a()
  expect_error(sensitivity_lob(a$blank, "missing", "sample"), "not found")
  expect_error(sensitivity_lob(a$blank, "result", "sample", alpha = 1),
               "strictly between")
  expect_error(sensitivity_verify(a$low, "result", "sample", "lod", 10),
               "lob")
  expect_error(sensitivity_verify(a$low, "result", "sample", "loq", 10),
               "assigned")
  lob <- sensitivity_lob(a$blank, "result", "sample", "lot")
  expect_centered_plot(plot(lob))
  expect_output(printed_lob <- withVisible(print(lob)), "LOB")
  expect_false(printed_lob$visible)
  expect_identical(printed_lob$value, lob)
  design <- sensitivity_design("all")
  expect_output(printed_design <- withVisible(print(design)), "study design")
  expect_false(printed_design$visible)
  expect_identical(printed_design$value, design)
})
