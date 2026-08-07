# Stability, QC, Fitting, Data Checks, and Sample Size

## Stability

Use `stability_regression()` for single-condition drift (`mode="self"`) or paired test/reference-condition bias (`mode="compare"`). Require explicit condition mapping, time unit, bias scale, direction when scientifically known, and allowable limit when making an acceptability conclusion.

```r
s <- ivdtools::stability_regression(data, time = "day", value = "result",
  condition = "condition", mode = "compare", test_condition = "25C",
  reference_condition = "2-8C", bias_type = "relative", direction = "decrease",
  conf.level = 0.95, limit = 5, time_unit = "d")
t <- ivdtools::stability_time(s, limit = 5)
plot(s)
```

Other functions include `stability_bias()`, `stability_plan()`, `mkt()`, and `arrhenius()`. For MKT/Arrhenius, confirm temperature units, exposure durations, activation energy assumptions, reaction order, target temperature, and bias limit. Export aggregated time nodes, model coefficients, confidence limits, estimated stability time, assumptions, and plots. Do not extrapolate past observed time without prominent qualification.

## Quality control

```r
q <- ivdtools::qc_chart(data, value = "result",
  rules = "1-2s,1-3s,2-2s,R-4s,4-1s,10x", mean = target_mean,
  sd = target_sd, group = "lot", run = "run")
print(q)
plot(q)
```

Use `list_westgard()` to verify rule names. Prefer protocol target mean/SD; if estimated from the same observations, label the chart as retrospective/exploratory. Preserve run order and grouping. Export violations with rule and run identifiers plus the Levey-Jennings plot.

## Equation fitting

Use `list_equation()` and `list_sadler()` before selecting a built-in model. Fit with `fit_equation(eq, data, x, y, start, lower, upper, constraints, weights)`. Use `compare_equation()` only across scientifically plausible candidate models. Record starting values, bounds, weighting, convergence, coefficients, fit metrics, residual diagnostics, and prediction range. Do not use predictions outside the calibration range without warning.

For replicated responses, `replicate_to_mean()` can construct means and weights; retain replicate counts and state the weighting rule. Available weight strings include `equal`, `1/y`, `1/y^2`, `1/x`, and `1/x^2`.

## Analytical sensitivity (LoB / LoD / LoQ)

Use `lob_lod_loq()` for CLSI EP17-A2 analytical sensitivity limits. Input is a **summarized** data.frame with one row per sample: a sample-name column, a mean column (the concentration axis), an SD column, and a replicate-count (n) column. Column names are free-form; pass them via `sample_col`, `mean_col`, `sd_col`, `n_col`. `blank` is a vector of blank sample names in `sample_col`. The precision profile is fit internally with the Sadler variance models (best by AIC).

```r
s <- ivdtools::lob_lod_loq(summ,
  sample_col = "sample", mean_col = "mean",
  sd_col = "sd", n_col = "n",
  blank = c("blank_A", "blank_B"), target_cv = 0.20)
print(s)
plot(s)
```

- **LoB** (5.3.3.1, parametric): pooled blanks, `LoB = M_B + cp*SD_B` with the EP17 small-sample correction. Blank samples are **excluded** from the profile fit.
- **LoD** (5.4.3): fixed-point solution `X = LoB + cp*SD(X)` on the fitted precision profile.
- **LoQ** (Appendix D1, functional sensitivity): concentration where `CV(X) = target_cv`. This LoQ reflects **precision only, not bias**. If LoQ < LoD, it is reported as `max(LoQ, LoD)` with a warning.
- When only blanks are given, only LoB is computed; when no blanks are given, only LoQ is computed.
- Report: LoB/LoD/LoQ, blank counts (B, K), the best Sadler model equation, its parameter estimates, and fit metrics (AIC, RSS, GoF p). State the target CV and that the LoQ is precision-only.

## Outliers and normality

Use `outliers_test(data, col, method=...)` with `grubbs`, `esd`, `dixon`, or `iqr`; use `normal_test(data, col, method=..., level=...)` with `shapiro`, `ad`, `lillie`, or `cvm`. Confirm method assumptions and sample-size suitability. Treat flags as review candidates, never automatic deletion instructions.

## Sample size

Available functions are:

- `sample_size_bland_altman(n=NULL, power=NULL, mu, sd, delta, conf.level=0.95, agree.level=0.95)`
- `sample_size_proportion(p0, p1, n=NULL, alpha=NULL, power=NULL, alternative=..., method=...)`
- `sample_size_proportion_ci(p, n=NULL, d=NULL, conf.level=0.95, method=...)`
- `stability_plan(allowable_drift, variability, replicates=1, expected_drift=NULL, power=c(0.8,0.9), ...)`

State every assumption, effect size, precision target, confidence level, power, sidedness, interval method, expected attrition, and whether the returned number is evaluable or enrolled subjects. Do not invent inputs from desired conclusions.
