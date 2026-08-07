# Qualitative Agreement, Precision, and Reference Analyses

## Qualitative agreement and diagnostic accuracy

Use `raw_to_table()` for row-level paired categorical results or `counts_to_table()` for known TP/FP/TN/FN counts. Explicitly identify the positive level.

```r
tab <- ivdtools::raw_to_table(data, candidate = "candidate", reference = "reference",
                              id = "sample_id", positive = "positive")
tab <- ivdtools::describe(tab)
tab <- ivdtools::diagnostics(tab, conf.level = 0.95, ci.method = "wilson")
tab <- ivdtools::kappa(tab, conf.level = 0.95)
tab <- ivdtools::mcnemar(tab, conf.level = 0.95)
summary(tab)
```

`diagnostics()` can accept an external prevalence for predictive values; label it clearly. Report the complete 2x2 table, sensitivity, specificity, predictive values, agreement, kappa, McNemar result, confidence method, positive level, missing pairs, and duplicate IDs. Do not call agreement diagnostic accuracy unless the reference is an appropriate truth standard.

## Precision and variance components

Require an explicit design formula matching the study structure and factor columns, for example `y ~ day/run/rep`. Do not infer nesting from column names alone.

```r
p <- ivdtools::precision(data, form = result ~ day/run, by = "sample", NegVC = FALSE)
ivdtools::variance(p)
ivdtools::ci(p)
ivdtools::profile(p)
ivdtools::normal(p)
ivdtools::outlier(p)
summary(p)
plot(p)
```

`by` may contain multiple grouping columns. State whether negative variance components were allowed, show observations and design levels per group, and flag failed or singular fits. Export variance components, SD/CV estimates and confidence intervals, profile results when requested, and plots. Never reinterpret a malformed or unbalanced design without user confirmation.

## Reference intervals

```r
ri <- ivdtools::reference_interval(data, col = "result", interval = 0.95,
                                   ci = 0.95, method = "all", id = "sample_id")
print(ri)
plot(ri)
```

Methods are `percentile`, `parametric`, `robust`, or `all`. Do not select the parametric method without checking distributional assumptions, and do not substitute automated normality testing for protocol judgment. Report sample size, missing and duplicate counts, coverage, confidence level, method-specific limits and their confidence intervals, and the plot. Discuss whether partitioning or an independent eligibility review was performed when relevant.

## Bottle ANOVA and Youden plots

Use `bottle_anova(data, formula, conf.level=...)`, followed by `tukey()` when the study asks for pairwise bottle comparisons. Use `youden_plot()` for paired sample measurements with supplied or estimated target mean/SD. Report the formula, group sizes, ANOVA/Tukey results, multiplicity context, and plot.
