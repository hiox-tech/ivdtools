# ivdtools 0.2.5

## Method comparison

* Added Linnet constant-CV Deming and precision-profile Deming.
* Added replicate summaries and `lambda = "replicates"` for Deming regression.
* Added analytic, jackknife, rank-based, and Bootstrap confidence intervals.
* Added segmented analysis with `breakpoints`.
* Added median and Hodges–Lehmann Bland–Altman intervals.
* Added selectable percent-difference denominators for Bland–Altman and
  outlier analysis.

## Other analyses

* Added paired ROC AUC comparisons using EP24 or DeLong methods.
* Added EP28-style nonparametric and biweight reference intervals, with
  parametric, type-6, and Huber methods retained.
* Added `split_by_cutpoints()` and `split_by_factors()` for preparing analyses.
* Added pooled-variance weighting to linearity analysis.
* Added sample-size calculations for Bland–Altman studies and proportion
  confidence intervals and target-value tests.

## Precision and robustness

* Added EP05 precision reports, including confidence intervals and multisite
  designs.
* Expanded normality-test results and improved reference-outlier handling.
* Improved validation, printing, plotting, documentation, and regression tests.

# ivdtools 0.2.4

## Method comparison

* Added reference-scale breakpoints for segmented regression and Bland–Altman
  analysis.
* Added residual-SD weighted least squares and the corrected Passing–Bablok
  Algorithm I.
* Added confidence intervals for mean and median differences.
* Separated the Bland–Altman plotting axis from its percent-difference base.

# ivdtools 0.2.3

## Method comparison

* Standardized `mcr()` orientation: reference is X, candidate is Y, and bias is
  candidate minus reference.
* Added parametric and nonparametric Bland–Altman limits of agreement.

## Commutability

* Added `material_type_map` and the EP30 and IFCC Part 3 workflows.
* Added unequal-replicate handling, position diagnostics, calibration checks,
  and non-blocking study-design checks.
* Added component-level IFCC uncertainty output and expanded validation,
  summaries, plots, and documentation.

## Precision and maintenance

* Added `report.precision()` for EP05 summary tables.
* Improved validation, plotting, printing, documentation, and tests across the
  package.

# ivdtools 0.2.2

## Measurement uncertainty

* Removed `uncertainty_bias()`. Use `uncertainty_type_a()`,
  `uncertainty_type_b()`, and `uncertainty_combine()` instead.

# ivdtools 0.2.1

## Commutability

* Replaced the former new-point method with a replicate-data workflow.
* Added EP14 and IFCC Part 2 calculations, prediction intervals, bias
  comparisons, and expanded uncertainty.

# ivdtools 0.2.0

## New workflows

* Added linearity, interference, C5/C95, dilution, spike-recovery, and hook-
  effect analyses.
* Added measurement-uncertainty workflows and `reference_bias()`.

## Enhancements

* Expanded ROC graphics and analytical-sensitivity workflows.
* Expanded validation, S3 methods, plotting, and tests.

# ivdtools 0.1.3

## Analytical sensitivity

* Added explicit precision-profile model selection to the former
  `lob_lod_loq()` workflow.

# ivdtools 0.1.2

## Initial release

* Initial release.
