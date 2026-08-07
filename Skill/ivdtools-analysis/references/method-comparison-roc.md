# Method Comparison and ROC

## Method comparison

Require one row per paired candidate/reference observation, numeric candidate and reference columns, and preferably a stable sample ID.

```r
x <- ivdtools::mcr(data, id = "sample_id", candidate = "candidate", reference = "reference")
x <- ivdtools::describe(x)
x <- ivdtools::correlation(x, method = "pearson")
x <- ivdtools::regression(x, method = "deming", conf.level = 0.95, lambda = 1)
x <- ivdtools::bland_altman(x, type = "difference", x_axis = "mean",
                            agree.level = 0.95, conf.level = 0.95)
summary(x)
```

Available regression methods are `ols`, `wls`, `deming`, `wdeming`, and `pb`. Do not default to Deming merely because both axes are measured: confirm the protocol and the variance ratio `lambda`. Use `bias(x, mdl=...)` only when the medical decision level is supplied. Available Bland-Altman scales are `difference`, `ratio`, and `percent`; ratio methods require strictly positive values.

Check pairing, duplicate IDs, range coverage, replicate handling, missing pairs, influential observations, heteroscedasticity, and regression convergence. Do not delete an outlier based only on `outlier()`; report sensitivity analyses with and without exclusions when exclusion is justified independently.

Export regression coefficients and confidence intervals, correlation, Bland-Altman mean bias and limits of agreement, bias at supplied decision levels, exclusions, and the comparison/Bland-Altman plots.

## ROC analysis

Require numeric marker columns and exactly one binary reference column. Explicitly set `positive` for character/factor outcomes. If the reference is continuous, use `continuous_to_binary()` only with a user- or protocol-supplied cutoff.

```r
r <- ivdtools::roc(data, cols = c("marker1", "marker2"),
                   reference = "truth", id = "sample_id", positive = "positive")
r <- ivdtools::describe(r)
r <- ivdtools::auc(r)
r <- ivdtools::cutoff(r)
summary(r)
plot(r)
```

Use `mlr(r, cols=..., name=...)` for a prespecified multivariable logistic marker and `predict()` with an explicit `column_map` for new data. Avoid reporting a data-selected cutoff as externally validated. State the selection rule, positive class, inferred marker direction, class counts, missing exclusions, and whether evaluation reused the fitting sample.

Export AUC results, cutoff/sensitivity/specificity results, ROC coordinates where available, class counts, the ROC plot, and any prediction table. Interpret AUC separately from performance at a clinically chosen cutoff.
