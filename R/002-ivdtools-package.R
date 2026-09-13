#' ivdtools: Statistical Tools for In Vitro Diagnostic Method Evaluation
#'
#' The package provides a collection of statistical workflows for evaluating
#' in vitro diagnostic methods and reagents. Most multi-step workflows use S3
#' objects: construct an analysis object, add results with analysis generics,
#' and inspect the accumulated results with `summary()` or `plot()`.
#'
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom stats addmargins binomial cor dbinom density df.residual fitted pbeta
#'   pbinom plogis pnorm predict.glm pt qbeta qbinom sigma time uniroot
#' @importFrom utils tail
"_PACKAGE"

## Silence R CMD check notes for ggplot2 non-standard evaluation symbols used
## throughout the plotting methods.
utils::globalVariables(c(
  ".data", "..density..", "x", "y", "sample", "component", "VC", "SD", "CV",
  "mean", "lower", "upper", "fit", "residual", "method", "value", "variable",
  "group", "rule", "status", "level", "type", ".condition", ".time",
  "ci_limit", "density", "estimate", "inv_temp1000", "loa_lower_hi",
  "loa_lower_lo", "loa_upper_hi", "loa_upper_lo", "log_rate", "observed",
  "observed_bias", "relative_bias", "time", "ymax", "ymin"
))
