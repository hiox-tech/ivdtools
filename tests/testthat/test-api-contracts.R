test_that("all namespace exports are available with expected object kinds", {
  exports <- getNamespaceExports("ivdtools")
  expect_setequal(
    grep("^sensitivity_", exports, value = TRUE),
    c("sensitivity_lob", "sensitivity_lod", "sensitivity_loq",
      "sensitivity_verify", "sensitivity_design")
  )
  expect_false("lob_lod_loq" %in% exports)
  expect_false("uncertainty_bias" %in% exports)

  objects <- mget(exports, envir = asNamespace("ivdtools"), inherits = FALSE)
  data_exports <- grep("^ivd_.*_example$", exports, value = TRUE)
  function_exports <- setdiff(exports, data_exports)

  expect_true(all(vapply(objects[function_exports], is.function, logical(1))))
  expect_true(all(vapply(objects[data_exports], is.data.frame, logical(1))))
})

test_that("documented S3 methods are registered and dispatchable", {
  registered <- c(
    "auc.roc", "auc_compare.roc", "bias.mcr", "bland_altman.mcr", "ci.precision",
    "correlation.mcr", "cutoff.roc", "describe.fourfold_table",
    "describe.mcr", "describe.roc", "diagnostics.fourfold_table",
    "kappa.fourfold_table", "mcnemar.fourfold_table", "mlr.roc",
    "normal.precision", "outlier.mcr", "outlier.precision",
    "profile.precision", "regression.mcr", "variance.precision",
    "vc.precision", "print.sensitivity", "plot.sensitivity",
    "print.sensitivity_design", "print.commutability_result",
    "summary.commutability_result", "plot.commutability_result",
    "print.dilution_recovery", "plot.dilution_recovery",
    "print.spike_recovery", "plot.spike_recovery",
    "print.hook_effect", "plot.hook_effect",
    "print.uncertainty_component", "print.uncertainty_iqc",
    "print.uncertainty_homogeneity",
    "print.uncertainty_stability", "print.uncertainty_characterization",
    "print.uncertainty_budget", "plot.uncertainty_budget",
    "print.uncertainty_propagation", "plot.uncertainty_propagation",
    "print.uncertainty_traceability", "plot.uncertainty_traceability",
    "print.reference_bias",
    "print.linearity_replicates", "print.linearity_panel",
    "print.linearity_endpoints", "print.linearity", "plot.linearity",
    "print.interference_replicates", "print.interference_paired",
    "plot.interference_paired", "print.interference_dose_response",
    "plot.interference_dose_response", "predict.interference_dose_response",
    "print.interference_patient", "plot.interference_patient",
    "print.c5_c95", "plot.c5_c95", "print.cutpoint_split",
    "print.factor_split",
    "print.roc_auc_comparison"
  )
  for (method in registered) {
    parts <- strsplit(method, ".", fixed = TRUE)[[1]]
    generic <- parts[1]
    class <- paste(parts[-1], collapse = ".")
    expect_true(
      is.function(getS3method(
        generic, class, envir = asNamespace("ivdtools"), optional = TRUE
      )),
      info = method
    )
  }
})

test_that("constructors preserve stable S3 object contracts", {
  m <- mcr(ivd_mcr_example, "sid", "test", "ref")
  expect_s3_class(m, "mcr")
  expect_named(
    m,
    c(
      "call", "data", "id", "id_name", "candidate", "reference",
      "candidate_name", "reference_name", "weights", "weights_name", "n",
      "replicate_columns", "candidate_sd", "reference_sd",
      "candidate_n", "reference_n", "complete_idx", "print",
      "describe", "correlation", "regression",
      "outlier", "bland_altman", "bias"
    ),
    ignore.order = TRUE
  )

  r <- roc(ivd_roc_example, c("x1", "x2"), "ref")
  expect_s3_class(r, "roc")
  expect_equal(r$n, nrow(ivd_roc_example))
  expect_equal(r$n_pos + r$n_neg, r$n)
  expect_identical(r$auc_comparisons, list())

  q <- counts_to_table(10, 2, 20, 3)
  expect_s3_class(q, "fourfold_table")
  expect_equal(unname(q$table[3, 3]), 35)

  expect_s3_class(linearity_replicates(5, 5), "linearity_replicates")
  expect_s3_class(interference_replicates(5, 10), "interference_replicates")
  expect_s3_class(uncertainty_type_b("calibrator", 10, lower = 9, upper = 11),
                   "uncertainty_component")
})

test_that("print and summary methods return objects invisibly", {
  m <- mcr(ivd_mcr_example, "sid", "test", "ref")
  capture.output(printed <- withVisible(print(m)))
  expect_false(printed$visible)
  expect_identical(printed$value, m)

  capture.output(summarised <- withVisible(summary(m)))
  expect_false(summarised$visible)
  expect_identical(summarised$value, m)
})

test_that("invalid classes fail through generics", {
  expect_error(correlation(list()), "no applicable method")
  expect_error(auc(list()), "no applicable method")
  expect_error(auc_compare(list()), "no applicable method")
  expect_error(diagnostics(list()), "no applicable method")
  expect_error(variance(list()), "no applicable method")
})
