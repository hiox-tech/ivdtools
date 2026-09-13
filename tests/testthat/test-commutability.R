commutability_test_data <- function() {
  clinical <- expand.grid(
    sample_id = sprintf("CS%02d", 1:20),
    procedure = c("MP1", "MP2"),
    replicate = 1:3,
    stringsAsFactors = FALSE
  )
  concentration <- setNames(seq(10, 100, length.out = 20),
                            sprintf("CS%02d", 1:20))
  error <- rep(c(-0.2, 0, 0.2), each = 40)
  clinical$result <- concentration[clinical$sample_id] +
    ifelse(clinical$procedure == "MP2", 1 + 0.02 *
             concentration[clinical$sample_id], 0) + error
  clinical$material_type <- "clinical"
  clinical$position <- NA_integer_

  rm <- expand.grid(
    sample_id = c("RM1", "RM2"), procedure = c("MP1", "MP2"),
    position = 1:5, replicate = 1:3, stringsAsFactors = FALSE
  )
  rm_level <- ifelse(rm$sample_id == "RM1", 40, 70)
  rm$result <- rm_level + ifelse(rm$procedure == "MP2",
    1 + 0.02 * rm_level, 0) +
    rep(c(-0.1, 0, 0.1), each = 20)
  rm$material_type <- "rm"
  rbind(clinical, rm)
}

test_that("EP14 uses replicate precision and Deming prediction limits", {
  d <- commutability_test_data()
  x <- commutability(d, "sample_id", "material_type", "procedure", "result",
                     approach = "ep14")

  expect_s3_class(x, "commutability_result")
  expect_identical(x$approach, "ep14")
  expect_equal(nrow(x$results), 2)
  expect_equal(x$models$n_clinical, 20)
  expect_true(is.finite(x$models$lambda))
  expect_true(x$models$lambda > 0)
  expect_true(all(x$results$lower < x$results$predicted_y))
  expect_true(all(x$results$upper > x$results$predicted_y))
  expect_true(all(x$results$classification == "commutable"))
})

test_that("EP14 marks a matrix-shifted material noncommutable", {
  d <- commutability_test_data()
  shifted <- d$sample_id == "RM2" & d$procedure == "MP2"
  d$result[shifted] <- d$result[shifted] + 10
  x <- commutability(d, "sample_id", "material_type", "procedure", "result",
                     approach = "ep14")
  expect_identical(x$results$classification,
                   c("commutable", "noncommutable"))
})

test_that("IFCC MSSD and global difference in bias follow direct formulas", {
  d <- commutability_test_data()
  x <- commutability(d, "sample_id", "material_type", "procedure", "result",
    position = "position", approach = "ifcc", criterion = 1,
    bias_model = "constant_global", coverage_factor = 1.9)

  cs <- x$clinical
  expected_mssd <- sqrt(sum(diff(cs$bias)^2) / (2 * (nrow(cs) - 1)))
  expect_equal(x$diagnostics$s_mssd, expected_mssd, tolerance = 1e-14)
  expect_equal(x$results$clinical_bias,
               rep(mean(cs$bias), nrow(x$results)), tolerance = 1e-14)
  expect_equal(x$results$bias_difference,
               x$results$rm_bias - mean(cs$bias), tolerance = 1e-14)
  expect_equal(x$results$expanded_uncertainty,
               1.9 * x$results$standard_uncertainty, tolerance = 1e-14)
  expect_true(all(x$results$classification == "commutable"))
})

test_that("IFCC pooled position variance follows the official Part 2 workbook", {
  # Nilsson et al., Clin Chem 2018;64:455-464, official supplemental workbook.
  # The fixture omits MPx replicate 3 for CS6 and all results for CS21, as
  # directed by the supplemental explanation of calculations.
  d <- read.csv(testthat::test_path("fixtures", "ifcc-part2-example.csv"),
                stringsAsFactors = FALSE)
  x <- commutability(
    d, "sample_id", "material_type", "procedure", "result",
    position = "position", approach = "ifcc", scale = "ln",
    criterion = 0.12, bias_model = "constant_global",
    coverage_factor = 1.9, rm_position = "pooled"
  )

  expect_equal(x$diagnostics$n_clinical, 49)
  expect_equal(x$results$clinical_bias,
               rep(0.1809523809523809, 5), tolerance = 1e-12)
  expect_equal(x$diagnostics$s_b, 0.1086154920022030, tolerance = 1e-12)
  expect_equal(x$diagnostics$pooled_position_mean_sd_x,
               0.0673550468187914, tolerance = 1e-12)
  expect_equal(x$diagnostics$pooled_position_mean_sd_y,
               0.0504183866999046, tolerance = 1e-12)
  expect_equal(x$results$rm_standard_uncertainty,
               rep(0.0376263632294770, 5), tolerance = 1e-12)
  expect_equal(x$results$clinical_standard_uncertainty,
               rep(0.0155164988574576, 5), tolerance = 1e-12)
  expect_equal(x$results$standard_uncertainty,
               rep(0.0407001836196106, 5), tolerance = 1e-12)
  expect_equal(x$results$expanded_uncertainty,
               rep(0.0773303488772602, 5), tolerance = 1e-12)
  expect_equal(x$results$bias_difference, c(
    0.1583068783068778, -0.0519153439153432, -0.0209523809523808,
    -0.1076190476190476, 0.6777142857142863
  ), tolerance = 1e-12)
  expect_identical(x$results$classification, c(
    "inconclusive", "inconclusive", "commutable", "inconclusive",
    "noncommutable"
  ))
  expect_true(all(x$results$position_pooling == "pooled"))
  expect_equal(nrow(x$position_diagnostics), 10)
})

test_that("IFCC material-specific position mode preserves the former formula", {
  d <- commutability_test_data()
  x <- commutability(
    d, "sample_id", "material_type", "procedure", "result",
    position = "position", approach = "ifcc", criterion = 1,
    bias_model = "constant_global", rm_position = "material_specific"
  )
  expected_rm_u <- sqrt(
    x$results$position_mean_sd_x^2 / x$results$positions_x +
      x$results$position_mean_sd_y^2 / x$results$positions_y
  )
  expect_equal(x$results$rm_standard_uncertainty, expected_rm_u,
               tolerance = 1e-14)
  expect_true(all(x$results$position_pooling == "material_specific"))
})

test_that("IFCC pooled analysis isolates materials without usable positions", {
  d <- commutability_test_data()
  missing <- d$sample_id == "RM1" & d$procedure == "MP2"
  d$position[missing] <- NA_integer_
  x <- commutability(
    d, "sample_id", "material_type", "procedure", "result",
    position = "position", approach = "ifcc", criterion = 1,
    bias_model = "constant_global", rm_position = "pooled"
  )
  expect_identical(x$results$classification[x$results$rm == "RM1"],
                   "not_evaluable")
  expect_true(x$results$classification[x$results$rm == "RM2"] !=
                "not_evaluable")
})

test_that("IFCC pooled position variance uses position degrees of freedom", {
  d <- rbind(
    data.frame(sample = "RM1", type = "rm", procedure = "MP1",
               result = rep(c(0, 2), each = 2),
               position = rep(1:2, each = 2)),
    data.frame(sample = "RM2", type = "rm", procedure = "MP1",
               result = rep(c(0, 1, 2), each = 2),
               position = rep(1:3, each = 2))
  )
  pooled <- ivdtools:::.comm_position_pool(d, c("RM1", "RM2"), "MP1")
  expect_equal(pooled$pooled_position_mean_sd, sqrt(4 / 3),
               tolerance = 1e-14)
  expect_equal(pooled$pooled_df, 3)
})

test_that("IFCC classification includes noncommutable and inconclusive states", {
  expect_identical(ivdtools:::.comm_ifcc_class(0, 0.2, 1), "commutable")
  expect_identical(ivdtools:::.comm_ifcc_class(2, 0.2, 1), "noncommutable")
  expect_identical(ivdtools:::.comm_ifcc_class(0.9, 0.2, 1), "inconclusive")
  expect_identical(ivdtools:::.comm_ifcc_class(1, 0, 1), "commutable")
})

test_that("local IFCC analysis requires bracketing clinical samples", {
  d <- commutability_test_data()
  d$result[d$sample_id == "RM2"] <- 200
  x <- commutability(d, "sample_id", "material_type", "procedure", "result",
    position = "position", approach = "ifcc", criterion = 1,
    bias_model = "constant_local", local_n = 8)
  expect_true("not_evaluable" %in% x$results$classification)
})

test_that("commutability validates scales, criteria, and replicate structure", {
  d <- commutability_test_data()
  expect_error(commutability(d, "sample_id", "material_type", "procedure",
                             "result", approach = "ifcc"), "position")
  expect_error(commutability(d, "sample_id", "material_type", "procedure",
    "result", position = "position", approach = "ifcc"), "criterion")
  expect_error(commutability(d, "sample_id", "material_type", "procedure",
                             "result", approach = "ep14", scale = "ln"),
               "EP14 supports")
  d$result[1] <- 0
  expect_error(commutability(d, "sample_id", "material_type", "procedure",
                             "result", approach = "ep14", scale = "log10"),
               "positive")

  d <- commutability_test_data()
  expect_error(commutability(
    d, "sample_id", "material_type", "procedure", "result",
    position = "position", approach = "ifcc", criterion = 1,
    rm_position = "invalid"
  ), "arg")
  incomplete <- d$sample_id == "RM1" & d$procedure == "MP2" &
    d$position == 5
  relaxed <- commutability(
    d[!incomplete, ], "sample_id", "material_type", "procedure", "result",
    position = "position", approach = "ifcc", criterion = 1
  )
  matching <- relaxed$design_checks$requirement ==
    "Matching RM positions between procedures"
  expect_false(relaxed$design_checks$met[matching])
  expect_true(all(relaxed$results$classification != "not_evaluable"))
})

test_that("material types support automatic recognition and explicit mapping", {
  d <- commutability_test_data()
  d$type_alias <- ifelse(d$material_type == "clinical", "patient",
                         "reference_material")
  automatic <- commutability(
    d, "sample_id", "type_alias", "procedure", "result",
    approach = "ep14"
  )
  expect_setequal(unique(automatic$input_data$type), c("clinical", "rm"))

  d$type_custom <- ifelse(d$material_type == "clinical", "patient_sample",
                          "candidate_rm")
  mapping <- c(patient_sample = "clinical", candidate_rm = "rm")
  explicit <- commutability(
    d, "sample_id", "type_custom", "procedure", "result",
    approach = "ep14", material_type_map = mapping
  )
  expect_setequal(unique(explicit$input_data$type), c("clinical", "rm"))
  expect_identical(explicit$settings$material_type_map, mapping)

  expect_error(
    commutability(d, "sample_id", "type_custom", "procedure", "result",
                  approach = "ep14"),
    "Unrecognized material type.*material_type_map"
  )
  expect_error(
    commutability(d, "sample_id", "type_alias", "procedure", "result",
                  approach = "ep14",
                  material_type_map = c(patient = "clinical")),
    "reference_material.*material_type_map"
  )
  expect_error(
    commutability(d, "sample_id", "type_custom", "procedure", "result",
                  approach = "ep14",
                  material_type_map = c(patient_sample = "other",
                                        candidate_rm = "rm")),
    "values must be 'clinical' or 'rm'"
  )
})

test_that("commutability print, summary, and plot methods preserve contracts", {
  d <- commutability_test_data()
  x <- commutability(d, "sample_id", "material_type", "procedure", "result",
                     approach = "ep14")
  capture.output(p <- withVisible(print(x)))
  capture.output(s <- withVisible(summary(x)))
  expect_false(p$visible)
  expect_identical(p$value, x)
  expect_false(s$visible)
  expect_identical(s$value, x)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_centered_plot(suppressMessages(plot(x)))
})

test_that("EP30 reproduces Appendix B regression and both classifications", {
  d <- read.csv(testthat::test_path("fixtures", "ep30-appendix-b.csv"),
                stringsAsFactors = FALSE)
  d$material_type <- ifelse(d$type == "patient", "clinical", "rm")
  x <- commutability(
    d, "specimen", "material_type", "method", "result",
    approach = "ep30", scale = "raw",
    procedure_pairs = list(c("MA", "MB"), c("MA", "MD"))
  )

  expect_equal(x$models$intercept, c(3.660, 0.396), tolerance = 0.02)
  expect_equal(x$models$slope, c(0.961, 0.834), tolerance = 0.002)
  expect_equal(x$models$syx, c(1.188, 1.033), tolerance = 0.01)
  expect_equal(x$models$critical_value, rep(stats::qnorm(0.975), 2),
               tolerance = 1e-14)
  expect_true(all(is.finite(x$models$variance_intercept)))
  expect_true(all(is.finite(x$models$variance_slope)))
  expect_true(all(is.finite(x$models$covariance_intercept_slope)))

  official_mb <- c(1.088, -0.117, 0.601, 2.164, -0.523, 3.241, -0.923)
  official_md <- c(-0.603, 0.518, 0.326, 0.997, -0.021, -6.758, -0.066)
  expect_equal(x$results$relative_residual[x$results$pair == "MA vs MB"],
               official_mb, tolerance = 0.07)
  expect_equal(x$results$relative_residual[x$results$pair == "MA vs MD"],
               official_md, tolerance = 0.08)
  expect_identical(
    x$results$prediction_classification,
    rep(c(rep("commutable", 5), "noncommutable", "commutable"), 2)
  )
  expect_identical(
    x$results$relative_residual_classification,
    c("commutable", "commutable", "commutable", "noncommutable",
      "commutable", "noncommutable", "commutable",
      rep("commutable", 5), "noncommutable", "commutable")
  )
  expect_equal(nrow(x$prediction_grid), 400)
})

test_that("EP30 validates study structure, adjustment, and extrapolation", {
  d <- read.csv(testthat::test_path("fixtures", "ep30-appendix-b.csv"),
                stringsAsFactors = FALSE)
  d$material_type <- ifelse(d$type == "patient", "clinical", "rm")
  x <- commutability(
    d, "specimen", "material_type", "method", "result",
    approach = "ep30", procedure_pairs = list(c("MA", "MB")),
    rm_set_size = 7
  )
  expect_equal(x$models$critical_value,
               stats::qnorm(1 - 0.05 / (2 * 7)), tolerance = 1e-14)

  d2 <- d
  outside <- d2$specimen == "RMA" & d2$type == "reference_material"
  d2$result[outside] <- d2$result[outside] + 1000
  x2 <- commutability(
    d2, "specimen", "material_type", "method", "result",
    approach = "ep30", procedure_pairs = list(c("MA", "MB"))
  )
  expect_identical(x2$results$classification[x2$results$rm == "RMA"],
                   "not_evaluable")
  expect_error(commutability(
    d, "specimen", "material_type", "method", "result",
    approach = "ep30", procedure_pairs = list(c("MA", "MB")),
    rm_set_size = 0
  ), "positive integer")
  expect_error(commutability(
    d, "specimen", "material_type", "method", "result",
    approach = "ep30", procedure_pairs = list(c("MA", "MB")),
    relative_residual_limit = 0
  ), "positive number")
})

test_that("calibration reproduces the IFCC Part 3 official example", {
  d <- read.csv(testthat::test_path("fixtures", "ifcc-part3-example.csv"),
                stringsAsFactors = FALSE)
  x <- commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    criterion = 6, w3_reference_procedure = "MP1"
  )

  expect_equal(x$targets$target[c(1, 40)], c(29.4, 114.2),
               tolerance = 1e-14)
  expected_medians <- c(
    -0.4, 9.5, 20.7, 46.6, -30.0, -19.4, -11.4,
    -0.8, -0.3, -0.2, -0.3, 0.4, -21.3, -1.1
  )
  expect_true(all(abs(x$results$median_bias_percent - expected_medians) < 0.11))
  expect_equal(x$diagnostics$impbr_percent, c(76.6, 21.6),
               tolerance = 0.05)
  expect_identical(x$diagnostics$classification,
                   c("screening", "requires_investigation"))

  expected_sd <- c(4.01, 5.72, 9.54, 6.73, 4.59, 5.30, 5.65,
                   4.05, 4.99, 7.83, 3.80, 6.34, 4.89, 6.45)
  expected_w3 <- c(4.01, 5.57, 5.90, 5.89, 3.49, 3.96, 4.63,
                   4.05, 4.98, 4.95, 3.45, 4.39, 3.57, 5.06)
  expect_equal(x$results$sd_bias_percent, expected_sd, tolerance = 0.02)
  expect_equal(x$results$adjusted_w3_percent, expected_w3, tolerance = 0.10)
  expect_lt(x$results$trend_p_value[
    x$results$calibration_stage == "before" & x$results$procedure == "MP4"
  ], 0.0001)
  expect_lt(abs(x$results$trend_p_value[
    x$results$calibration_stage == "before" & x$results$procedure == "MP7"
  ] - 0.0497), 0.0001)
})

test_that("calibration retains the full result when MP6 is explicitly excluded", {
  d <- read.csv(testthat::test_path("fixtures", "ifcc-part3-example.csv"),
                stringsAsFactors = FALSE)
  x <- commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    criterion = 6, exclude_procedures = "MP6",
    w3_reference_procedure = "MP1"
  )
  after_all <- x$diagnostics$calibration_stage == "after" &
    x$diagnostics$assessment_set == "all"
  after_excluded <- x$diagnostics$calibration_stage == "after" &
    x$diagnostics$assessment_set == "specified_exclusions"
  expect_equal(x$diagnostics$impbr_percent[after_all], 21.6, tolerance = 0.05)
  expect_equal(x$diagnostics$impbr_percent[after_excluded], 1.5,
               tolerance = 0.11)
  expect_identical(x$diagnostics$classification[after_excluded], "commutable")
  expect_false(x$results$included[x$results$procedure == "MP6"][1])
  expect_equal(x$leave_one_out$impbr_percent[
    x$leave_one_out$omitted_procedure == "MP6"
  ], 1.5, tolerance = 0.11)

  capture.output(p <- withVisible(print(x)))
  capture.output(s <- withVisible(summary(x)))
  expect_false(p$visible)
  expect_false(s$visible)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(suppressWarnings(plot(x)), "ggplot")
})

test_that("calibration validates stages, grids, and pairwise arguments", {
  d <- read.csv(testthat::test_path("fixtures", "ifcc-part3-example.csv"),
                stringsAsFactors = FALSE)
  expect_error(commutability(
    d[d$calibration_stage == "before", ], sample_id = "sample_id",
    procedure = "procedure", result = "result", approach = "calibration",
    calibration_stage = "calibration_stage", criterion = 6
  ), "before.*after")
  relaxed <- commutability(
    d[-1, ], sample_id = "sample_id", procedure = "procedure",
    result = "result", approach = "calibration",
    calibration_stage = "calibration_stage", criterion = 6
  )
  expect_length(relaxed$dropped_samples, 1)
  expect_equal(nrow(relaxed$targets), 39)
  grid_check <- relaxed$design_checks$requirement ==
    "Complete matched before/after grid"
  expect_false(relaxed$design_checks$met[grid_check])
  expect_error(commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    criterion = 6, procedure_pairs = list(c("MP1", "MP2"))
  ), "not applicable")
  expect_error(commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    criterion = 6, scale = "ln"
  ), "raw")
  expect_error(commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    criterion = 6, exclude_procedures = "missing"
  ), "observed procedure")
})

test_that("EP30 computes below design targets and with unequal replicates", {
  d <- read.csv(testthat::test_path("fixtures", "ep30-appendix-b.csv"),
                stringsAsFactors = FALSE)
  d$material_type <- ifelse(d$type == "patient", "clinical", "rm")
  patients <- unique(d$specimen[d$type == "patient"])[1:8]
  keep <- d$type == "reference_material" | d$specimen %in% patients
  d <- d[keep, ]
  d <- d[-which(d$type == "patient" & d$method == "MB")[1], ]
  d <- d[-which(d$type == "reference_material" & d$method == "MB")[1], ]

  x <- commutability(
    d, "specimen", "material_type", "method", "result",
    approach = "ep30", procedure_pairs = list(c("MA", "MB"))
  )
  expect_equal(x$models$n_clinical, 8)
  expect_false(x$models$clinical_replicates_equal)
  expect_false(x$models$rm_replicates_equal)
  expect_true(any(!x$design_checks$met))
  expect_true(all(is.finite(x$results$predicted_y)))
})

test_that("calibration reports small studies and safely omits unavailable W3", {
  d <- read.csv(testthat::test_path("fixtures", "ifcc-part3-example.csv"),
                stringsAsFactors = FALSE)
  samples <- unique(d$sample_id)[1:5]
  d <- d[d$sample_id %in% samples, ]
  duplicate <- d[1, ]
  d <- rbind(d, duplicate)

  x <- commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    criterion = 6
  )
  expect_equal(nrow(x$targets), 5)
  expect_true(all(is.na(x$results$w3_percent)))
  expect_true(all(is.na(x$results$adjusted_w3_percent)))
  expect_true(all(is.finite(x$diagnostics$impbr_percent)))
  sample_check <- x$design_checks$requirement == "Complete clinical samples"
  final_check <- x$design_checks$requirement ==
    "One final result per sample/procedure/stage"
  expect_false(x$design_checks$met[sample_check])
  expect_false(x$design_checks$met[final_check])
  expect_error(commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    calibration_replicate_summary = "error", criterion = 6
  ), "More than one calibration result")
})

test_that("calibration computes with more than 20 procedures", {
  d <- expand.grid(
    sample_id = paste0("CS", 1:6),
    procedure = paste0("MP", 1:21),
    calibration_stage = c("before", "after"),
    stringsAsFactors = FALSE
  )
  sample_level <- setNames(seq(50, 100, length.out = 6), paste0("CS", 1:6))
  procedure_bias <- setNames(seq(-2, 2, length.out = 21), paste0("MP", 1:21))
  sample_index <- as.integer(sub("CS", "", d$sample_id))
  procedure_index <- as.integer(sub("MP", "", d$procedure))
  d$result <- sample_level[d$sample_id] *
    (1 + procedure_bias[d$procedure] / 100) *
    ifelse(d$calibration_stage == "after", 0.999, 1) +
    0.05 * sin(sample_index * procedure_index)
  x <- commutability(
    d, sample_id = "sample_id", procedure = "procedure", result = "result",
    approach = "calibration", calibration_stage = "calibration_stage",
    criterion = 6
  )
  expect_equal(length(unique(x$results$procedure)), 21)
  procedure_check <- x$design_checks$requirement == "Measurement procedures"
  expect_false(x$design_checks$met[procedure_check])
})
