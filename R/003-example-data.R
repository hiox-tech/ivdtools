#' Small deterministic method-comparison example
#'
#' @format A data frame with 12 rows and columns `sid`, `test`, and `ref`.
#' @export
ivd_mcr_example <- data.frame(
  sid = seq_len(12),
  test = seq(10, 32, 2) + c(-0.1, 0.1),
  ref = seq(10, 32, 2)
)

#' Small deterministic ROC example
#'
#' @format A data frame with 12 rows, two markers, and a binary reference.
#' @export
ivd_roc_example <- data.frame(
  x1 = seq_len(12),
  x2 = c(1:6, 8:13),
  ref = rep(c(0, 1), 6)
)

#' Small deterministic qualitative-method example
#'
#' @format A data frame with paired binary results and sample identifiers.
#' @export
ivd_qualitative_example <- data.frame(
  id = seq_len(8),
  new = rep(c("positive", "negative"), 4),
  gold = rep(c("positive", "positive", "negative", "negative"), 2)
)

#' Small deterministic one-way ANOVA example
#'
#' @format A data frame with 12 observations from three bottles.
#' @export
ivd_bottle_example <- data.frame(
  bottle = rep(c("A", "B", "C"), each = 4),
  value = c(10.1, 9.9, 10.2, 10, 11.1, 10.9, 11.2, 11,
            12.1, 11.9, 12.2, 12)
)
