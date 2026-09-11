![ivdtools logo](logo.jpeg)

# R Package "ivdtools" 

## What is ivdtools · 功能介绍

`ivdtools` is an R package providing statistical workflows used in the evaluation of in vitro diagnostic (IVD) reagents, covering data inspection, statistical analysis, result summarization, and visualization. The current version requires R >= 4.1.0 and is distributed under the MIT license.

`ivdtools` 是一个用于体外诊断（IVD）试剂评价的 R 软件包，提供数据检查、统计分析、结果汇总和可视化等统计工作流程。当前版本要求 R >= 4.1.0，并采用 MIT 许可证发布。

The package provides the following analysis modules:

该软件包提供以下分析模块：

| Module · 模块 | Main entry points · 主要入口函数 | Related CLSI guideline · 相关 CLSI 指南 |
| :---: | :---: | :---: |
| Method comparison · 方法学比较 | `mcr()` | EP09 |
| Reference-material commutability · 参考物质互换性 | `commutability()` | EP14, EP30, and IFCC Part 2 |
| Precision · 精密度 | `precision()` | EP05 and EP15 |
| Qualitative agreement · 定性一致性 | `raw_to_table()` and others | EP12 |
| C5/C95 estimation · C5/C95 估计 | `c5_c95()` | EP12 |
| Reference interval · 参考区间 | `reference_interval()` | EP28 |
| Stability · 稳定性 | `stability_*()`, `arrhenius()`, `mkt()` | EP25 |
| Linearity · 线性 | `linearity()` and `linearity_*()` | EP06 |
| Curve fitting · 曲线拟合 | `fit_equation()` family | General · 通用 |
| Interference · 干扰 | `interference_*()` | EP07 |
| Dilution, spiking, and hook effect · 稀释、加标和钩状效应 | `dilution_recovery()`, `spike_recovery()`, `hook_effect()` | EP34 |
| Measurement uncertainty · 测量不确定度 | `uncertainty_*()` | EP29 |
| Reference-material bias · 参考物质偏倚 | `reference_bias()` | YY/T 1789.2 |
| LoB/LoD/LoQ · 空白限、检出限、定量限| `sensitivity_lob()` / `sensitivity_lod()` / `sensitivity_loq()` | EP17 |
| Bottle-to-bottle ANOVA · 瓶间差方差分析 | `bottle_anova()` | EP15 |
| ROC analysis · ROC 分析 | `roc()` family | EP24 |
| Quality control · 质量控制 | `qc_chart()`, `youden_plot()` | Westgard rules · Westgard 规则 |
| Outliers and normality tests · 异常值与正态性检验 | `outliers_test()`, `normal_test()` | General · 通用 |
| Sample size · 样本量 | `sample_size_*()` | General · 通用 |

Statistical results should be interpreted in the context of a pre-specified study protocol, applicable standards, and clinical or analytical acceptance limits; they are not a substitute for professional judgment.

统计结果应结合预先规定的研究方案、适用标准以及临床或分析性能接受限进行解释；统计结果不能替代专业判断。

## Online documentation · 在线文档

Detailed tutorials covering data import, analysis environment setup, and complete worked examples for every analysis module are available:

这里提供涵盖数据导入、分析环境配置以及各分析模块完整操作示例的详细教程：

English handbook · 英文手册：<https://hiox-tech.github.io/ivdtools/>

Chinese handbook · 中文手册：<https://ivdtools.hiox-tech.cn/>

Manual · 软件包手册：<https://cran.r-project.org/web/packages/ivdtools/refman/ivdtools.html>

## Homepage · 项目主页

For PDF manuals, cheatsheets and more resources, visit:

如需获取 PDF 手册、速查表及其他资源，请访问：

GitHub · 英文：<https://github.com/hiox-tech/ivdtools>

Gitee · 中文：<https://gitee.com/hiox-tech/ivdtools>

CRAN：<https://CRAN.R-project.org/package=ivdtools>

## Feedback · 反馈

Bug reports, feature requests, and feedback:

如需报告错误、提出功能请求或反馈意见，请访问：

<https://github.com/hiox-tech/ivdtools/issues>

## Author · 作者与联系方式

**hiox-tech** <GeorgeBinDragon@outlook.com>

---

*Built for IVD reagent evaluation — designed for reproducibility and clarity.*

*为 IVD 试剂评估构建——追求清晰与可复现。*
