---
name: ivdtools-analysis
description: Analyze in-vitro diagnostic (IVD) evaluation data with the R package ivdtools (>= 0.1.1, from CRAN) and produce reproducible Chinese HTML reports, R code, result tables, and figures. Use when inspecting CSV, TSV, Excel, or RDS data and performing method comparison, Bland-Altman, ROC, qualitative agreement, diagnostic accuracy, precision or variance components, analytical sensitivity (LoB / LoD / LoQ), reference intervals, stability, quality control, response-curve fitting, outlier or normality assessment, bottle ANOVA, or sample-size calculations.
---

# IVDTools Analysis

Use `ivdtools` as the statistical engine while preserving every analysis decision and output. Treat package results as statistical evidence, not clinical validation or regulatory approval.

## Workflow

1. Resolve this skill directory from the location of this `SKILL.md`. Create a new analysis output directory outside the skill. Never edit the source data.
2. Run `Rscript <skill-dir>/scripts/check_environment.R`. If it fails, explain the missing components. Run `install_ivdtools.R --yes` only after the user explicitly authorizes package installation or environment changes.
3. Run `inspect_data.R` on each input file before writing analysis code. For Excel, select a sheet explicitly when the workbook has more than one relevant sheet.
4. Select the analysis family and read only the matching reference:
   - Method comparison, Bland-Altman, bias, ROC: `references/method-comparison-roc.md`
   - Qualitative agreement, precision, reference intervals, bottle ANOVA: `references/qualitative-precision-reference.md`
   - Stability, QC, equation fitting, analytical sensitivity, data checks, sample size: `references/stability-qc-fitting-sample-size.md`
5. Confirm only decisions that cannot be derived safely: column mapping, positive class, pairing or replicate structure, test/reference condition, statistical method, clinically allowable limits, target mean/SD, confidence level, or desired power. Never invent a clinical cutoff, allowable bias, stability limit, or acceptance criterion.
6. Copy `assets/report-template.Rmd` into the output directory. Create `analysis.R` for reusable analysis functions and complete the report with explicit parameters and `ivdtools::` calls. Use namespace-qualified calls because `kappa()` and `profile()` can mask functions from base or stats.
7. Execute the analysis, export machine-readable tables, save plots, and render with `render_report.R`. Capture warnings and errors in the report instead of silently suppressing them.
8. Validate row counts, excluded observations, direction/level choices, convergence status, and agreement between reported values and exported tables. Record `sessionInfo()` and package versions.

## Input Rules

- Accept `.csv`, `.tsv`, `.txt`, `.xls`, `.xlsx`, and `.rds`; require an RDS object to be a data frame.
- Preserve original column names with `check.names = FALSE`.
- Report missingness, duplicate rows, duplicate IDs, nonnumeric analysis values, factor levels, and rows excluded by each operation.
- Do not impute, remove outliers, average replicates, convert a continuous reference to binary, or reverse outcome direction without stating the operation and its rationale.
- Ask when multiple mappings are plausible or a binary positive class is not explicit. Otherwise proceed after validation.

## Output Contract

Create a timestamped directory named `<input-stem>-ivdtools-analysis-YYYYMMDD-HHMMSS` unless the user supplies a destination. Do not overwrite an existing directory. Include:

```text
analysis.R
report.Rmd
report.html
analysis-manifest.json
session-info.txt
results/*.csv
figures/*.png
```

Write the report in Chinese unless the user requests another language. Include purpose, data provenance, data-quality findings, methods and parameters, results, plots, exclusions/warnings, interpretation, limitations, and reproducibility information. Keep raw observations out of the narrative unless needed; result tables may contain only fields required by the analysis.

## Environment and Version Policy

- Require R >= 4.1.0 and `ivdtools` >= 0.1.1 (the version that introduced `lob_lod_loq()` for analytical sensitivity) for this skill revision.
- Install `ivdtools` and its dependencies from CRAN. `check_environment.R` verifies the minimum version; `install_ivdtools.R` installs from CRAN.
- Treat dependency installation, upgrades, and writes to an R library as environment changes requiring explicit user authorization.
- Do not replace a newer installed `ivdtools`; if the installed version is below 0.1.1, install into a user-selected library or stop with instructions.

## Interpretation Boundaries

- Distinguish statistical significance from clinical acceptability.
- Tie conclusions to the supplied study protocol and acceptance limits.
- Label exploratory analyses and post-hoc choices.
- State that synthetic/package examples verify software operation only.
- Recommend independent review for clinical, regulatory, or release decisions.
