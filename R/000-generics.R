#' Describe an analysis object
#'
#' @param x An S3 analysis object.
#' @param ... Additional arguments passed to a method.
#' @return The value returned by the dispatched method.
#' @examples
#' obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#' describe(obj)
#' @export
describe <- function(x, ...) UseMethod("describe")

#' Detect outliers in an analysis object
#'
#' @inheritParams describe
#' @return The value returned by the dispatched method.
#' @examples
#'   obj <- mcr(ivd_mcr_example, "sid", "test", "ref")
#'   outlier(obj, method = "grubbs")
#' @export
outlier <- function(x, ...) UseMethod("outlier")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
correlation <- function(x, ...) UseMethod("correlation")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
regression <- function(x, ...) UseMethod("regression")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
bias <- function(x, ...) UseMethod("bias")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
bland_altman <- function(x, ...) UseMethod("bland_altman")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
report <- function(x, ...) UseMethod("report")

#' @param object An S3 analysis object.
#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
profile <- function(object, ...) UseMethod("profile")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
normal <- function(x, ...) UseMethod("normal")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
ci <- function(x, ...) UseMethod("ci")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
variance <- function(x, ...) UseMethod("variance")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
vc <- function(x, ...) UseMethod("vc")

#' @param object An S3 analysis object.
#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
diagnostics <- function(object, ...) UseMethod("diagnostics")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
kappa <- function(x, ...) UseMethod("kappa")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
mcnemar <- function(x, ...) UseMethod("mcnemar")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
cutoff <- function(x, ...) UseMethod("cutoff")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
mlr <- function(x, ...) UseMethod("mlr")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
auc <- function(x, ...) UseMethod("auc")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
auc_compare <- function(x, ...) UseMethod("auc_compare")

#' @rdname describe
#' @return The value returned by the dispatched method.
#' @export
tukey <- function(x, ...) UseMethod("tukey")
