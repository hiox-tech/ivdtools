<p align="center">

<img src="logo.jpeg" alt="ivdtools logo" width="200"/>

</p>

<h1 align="center">

ivdtools

</h1>

<p align="center">

<em>Statistical Tools for Evaluation of in Vitro Diagnostic Reagents - R Package</em>

</p>

<p align="center">

<em>体外诊断试剂性能评估统计工具 — R 语言工具集</em>

</p>

<p align="center">

<a href="#"><img src="https://img.shields.io/badge/R-%3E%3D%204.0-blue" alt="R &gt;= 4.0"/></a> <a href="#"><img src="https://img.shields.io/badge/CLSI-EP05%20%7C%20EP06%20%7C%20EP09%20%7C%20EP12%20%7C%20EP15%20%7C%20EP17%20%7C%20EP25%20%7C%20EP28-blueviolet" alt="CLSI standards"/></a> <a href="#"><img src="https://img.shields.io/badge/license-MIT-green" alt="License MIT"/></a>

</p>

------------------------------------------------------------------------

## Overview · 概述

**ivdtools** is an R toolset for **in vitro diagnostic (IVD) reagent performance verification and method comparison**, covering the full workflow from data exploration to statistical analysis and visualization. It follows the **S3 generic + slots** design pattern, providing a unified, consistent API and slot-based result storage.

**ivdtools** 是一套面向**体外诊断（IVD）试剂性能验证与方法学比较**的 R 工具集，覆盖从数据探索、统计分析到可视化呈现的完整工作流。采用 **S3 泛型 + 槽位存储**设计模式，接口风格统一，分析结果自动存入对象槽位。

------------------------------------------------------------------------

## Modules · 模块一览

| Module · 模块 | CLSI Standard | S3 Class | Key Functions · 核心函数 | Description · 功能描述 |
|----|----|----|----|----|
| **MCR** | EP09 | `mcr` | `mcr()`, `describe()`, `correlation()`, `regression()`, `bland_altman()`, `bias()`, `outlier()`, `plot()`, `summary()` | Method comparison regression (OLS / Deming / Passing-Bablok), Bland-Altman analysis, bias analysis · 方法学比较回归分析，Bland-Altman 分析，偏差分析 |
| **Precision** | EP05 | `precision` | `precision()`, `variance()`, `ci()`, `profile()`, `outlier()`, `normal()`, `plot()`, `summary()` | Variance component analysis (ANOVA-based VCA), Sadler precision profile fitting · 精密度方差分量分析，Sadler 精密度剖面拟合 |
| **Qualitative** | EP12 | `fourfold_table` | `raw_to_table()`, `counts_to_table()`, `describe()`, `diagnostics()`, `kappa()`, `mcnemar()`, `summary()` | 2×2 contingency table, diagnostic parameters (sensitivity, specificity, etc.), Kappa agreement, McNemar test · 2×2 四格表、诊断参数、Kappa 一致性、McNemar 检验 |
| **ROC** | General · 通用 | `roc` | `roc()`, `describe()`, `auc()`, `cutoff()`, `mlr()`, `plot()`, `summary()` | ROC curve, AUC (95% CI), optimal cutoff (Youden index / closest-to-(0,1)), multivariate logistic regression · ROC 曲线、AUC、最佳截断值、多变量 Logistic 回归 |
| **Reference Interval** | EP28 | `reference_interval` | `reference_interval()`, `print()`, `plot()` | Non-parametric (percentile), parametric (normal), robust (Horn & Pesce) methods · 非参数（百分位数）、参数法（正态）、稳健法 |
| **Stability** | EP25 | — | `stability_bias()`, `stability_regression()`, `stability_time()`, `mkt()`, `arrhenius()`, `stability_plan()` | Bias analysis, stability regression, time-to-out-of-spec prediction, mean kinetic temperature (MKT), Arrhenius accelerated model · 偏差分析、稳定性回归、超期预测、平均动力学温度、Arrhenius 加速模型 |
| **Fit/Regression** | EP06 | `fit_equation` | `fit_equation()`, `list_equation()`, `coef()`, `predict()`, `residuals()`, `plot()`, `compare_equation()` | Response curve fitting (linear, exponential, 4PLC, 5PLC, etc.), Sadler variance function, weighted fitting · 响应曲线拟合、权重拟合、参数约束 |
| **Bottle/Batch ANOVA** | EP15 | `bottle_anova` | `bottle_anova()`, `tukey()`, `print()` | One-way & nested ANOVA, post-hoc Tukey HSD, compact letter display · 单向和嵌套 ANOVA、事后 Tukey HSD、紧凑字母标记 |
| **QC** | Westgard | — | `list_westgard()`, `qc_chart()`, `youden_plot()` | Westgard rules, Levey-Jennings chart, Youden plot · Westgard 规则、Levey-Jennings 图、Youden 图 |
| **Analytical Sensitivity** | EP17 | `sensitivity` | `lob_lod_loq()` | LoB, LoD, LoQ 空白限、检出限、定量限 |
| **Sample Size** | General · 通用 | — | `sample_size_bland_altman()`, `sample_size_proportion_ci()`, `sample_size_proportion()` | Bland-Altman agreement sample size (Lu et al.), single proportion CI, one-arm target value test · Bland-Altman 一致性样本量、单比例置信区间、单臂目标值检验 |
| **Outliers & Normality** | General · 通用 | — | `outliers_test()`, `normal_test()` | Grubbs / ESD / Dixon Q / IQR outlier tests; Shapiro-Wilk / AD / Lilliefors / Cramer–von Mises normality tests with QQ plots · 离群值检测和正态性检验（含 QQ 图） |

------------------------------------------------------------------------

## Design Philosophy · 设计理念

### S3 Generics + Slot Storage · S3 泛型 + 槽位存储

Each analysis module centers on an **S3 class** object. Results from every analysis step are automatically stored in **named slots** within the object, allowing:

- **Step-by-step pipelining**: run analyses incrementally, results accumulate in the object
- **Unified plotting**: `plot()` provides a single entry point for all visualizations
- **Summary output**: `summary()` prints all stored results at once

每个分析模块围绕一个 **S3 类对象**构建。每次分析的中间与最终结果自动存入对象的命名**槽位**中，支持：

- **逐步流水线式分析**：分析结果不断追加到对象中
- **统一绘图入口**：所有可视化通过 `plot()` 统一调用
- **一键汇总输出**：`summary()` 打印全部已存储的分析结果

### Module Structure · 模块结构

```         
Constructor        → S3 object with print slot auto-filled
  ↓
Generics           → describe / regression / correlation / … (each fills one slot)
  ↓
summary()          → prints all filled slots (does not fill any itself)
plot()             → unified plotting entry (requires prior analysis for some plot types)
```

------------------------------------------------------------------------

## Getting Started · 快速开始

``` r
install.packages("ivdtools")
library("ivdtools")
```

### Example: MCR (Method Comparison) · 方法学比较示例

``` r
# Load data and create mcr object · 载入数据并创建对象
data <- read.csv("comparison_data.csv")
obj <- mcr(data, id = "sample_id",
           candidate = "new_method",
           reference = "standard_method")

# Step-by-step analysis · 逐步分析
obj <- describe(obj)           # data description · 数据描述
obj <- correlation(obj)        # correlation analysis · 相关分析
obj <- regression(obj)         # regression (OLS/Deming/PB) · 回归分析
obj <- bland_altman(obj)       # Bland-Altman analysis · BA 分析
obj <- outlier(obj)            # outlier detection · 离群值检测
obj <- bias(obj)               # bias analysis · 偏差分析

# Unified outputs · 统一输出
summary(obj)                   # print all results · 打印全部结果
plot(obj)
plot(obj, type = "regression") # specific plot · 指定绘图类型
plot(obj, type = "bland_altman")
```

### Example: Precision (VCA) · 精密度分析示例

``` r
obj <- precision(data, form = y ~ day/run, by = "sample")

obj <- outlier(obj)            # outlier detection · 离群值检测
obj <- normal(obj)             # normality test · 正态性检验
obj <- variance(obj)           # VCA · 方差分量分析
obj <- ci(obj)                 # confidence intervals · 置信区间
obj <- profile(obj)            # Sadler precision profile · 精密度剖面

summary(obj)
plot(obj, type = "dot")        # run-order scatter · 运行顺序散点图
plot(obj, type = "var")    # VCA bar chart · VCA 条形图
plot(obj, type = "profile")    # precision profile · 精密度剖面图
```

### Example: ROC · ROC 曲线示例

``` r
obj <- roc(data, cols = c("marker1", "marker2"),
           reference = "gold_standard", ref_cutoff = 10)

obj <- describe(obj)           # data description · 数据描述
obj <- auc(obj)                # AUC · AUC 计算
obj <- cutoff(obj)             # optimal cutoff · 最佳截断值
obj <- mlr(obj)                # multivariate logistic regression · 多变量 Logistic 回归

summary(obj)
plot(obj)
```

### Example: Qualitative Analysis · 定性分析示例

``` r
# From raw data · 从原始数据构建
obj <- raw_to_table(data,
                    candidate = "method_new",
                    reference = "method_standard",
                    positive = "positive")

# Or from counts · 或从计数值构建  
obj <- counts_to_table(tp = 45, fp = 3, tn = 97, fn = 5,
                       positive = "positive")

obj <- describe(obj)           # data summary · 数据描述
obj <- diagnostics(obj)        # diagnostic parameters · 诊断参数
obj <- kappa(obj)              # Kappa agreement · Kappa 一致性
obj <- mcnemar(obj)            # McNemar test · McNemar 检验

summary(obj)
```

---
## AI Skills (Codex, Claude) · AI自动分析

The `ivdtools` R package provides statistical functions, while the `ivdtools-analysis` Skill provides a constrained automated analysis workflow for AI. Once the user supplies the data file, study design, and analysis requirements, the Skill uses `ivdtools` as its statistical engine to organize data pre-check, parameter confirmation, statistical analysis, result verification, and Chinese report generation.

`ivdtools` R 包提供统计函数，`ivdtools-analysis` Skill 则面向 AI 提供一套受约束的 自动分析工作流。用户给出数据文件、研究设计和分析要求后，Skill 以 `ivdtools` 为统计 引擎，组织数据预检、参数确认、统计分析、结果核验和中文报告生成。

The value of using the Skill in a conversation is not to "replace statistical judgment" but to standardize repetitive technical steps, and to save the data source, analysis parameters, excluded records, warnings, result tables, figures, and software environment together, reducing the risk of omissions and manual copy errors.

在对话中使用 Skill 的价值不是“替代统计判断”，而是把重复的技术步骤标准化，并将数据 来源、分析参数、排除记录、警告、结果表、图形和软件环境一起保存，降低遗漏和手工复制 错误的风险。
---

## Dependencies · 依赖

| Package | Module | Purpose |
|----|----|----|
| `ggplot2` | QC, Fit, Stability, Reference Interval, MCR, ROC | Visualization · 可视化 |
| `VCA` | Precision | Variance component analysis · 方差分量分析 |
| `VFP` | Precision, Fit | Sadler precision profile fitting · Sadler 精密度剖面拟合 |
| `minpack.lm` | Fit | Nonlinear least squares (Levenberg-Marquardt) · 非线性最小二乘 |
| `nloptr` | Fit | Constrained optimization · 约束优化 |
| `nls2` | Fit | Robust starting values for NLS · NLS 稳健初值 |

------------------------------------------------------------------------

## License · 许可

MIT License

## Author · 作者

[hiox-tech](https://github.com/hiox-tech)

------------------------------------------------------------------------

<p align="center">

<sub>Built for IVD reagent verification — designed for reproducibility and clarity.</sub> <br> <sub>为 IVD 试剂验证构建 — 追求可复现与清晰。</sub>

</p>
