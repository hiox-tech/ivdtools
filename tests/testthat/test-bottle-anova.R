test_that("one-way bottle ANOVA agrees with base R", {
  d <- data.frame(
    bottle = rep(letters[1:3], each = 4),
    value = c(1, 2, 1, 2, 4, 5, 4, 5, 7, 8, 7, 8)
  )
  observed <- bottle_anova(d, value ~ bottle)
  direct <- summary(stats::aov(value ~ bottle, data = d))[[1]]

  expect_equal(observed$aov_table$Df, direct$Df)
  expect_equal(observed$aov_table$Sum_Sq, direct$`Sum Sq`, tolerance = 1e-12)
  expect_equal(observed$aov_table$Mean_Sq, direct$`Mean Sq`, tolerance = 1e-12)
  expect_equal(observed$desc_stats$Mean, c(1.5, 4.5, 7.5))
})

test_that("Tukey pairwise estimates agree with TukeyHSD", {
  d <- data.frame(
    bottle = rep(letters[1:3], each = 4),
    value = c(1, 2, 1, 2, 4, 5, 4, 5, 7, 8, 7, 8)
  )
  observed <- tukey(bottle_anova(d, value ~ bottle))
  direct <- stats::TukeyHSD(stats::aov(value ~ bottle, d))$bottle
  actual <- observed$tukey_list$bottle$matrix
  expect_equal(unname(actual[, "diff"]), unname(direct[, "diff"]),
               tolerance = 1e-12)
  expect_equal(unname(actual[, "p adj"]), unname(direct[, "p adj"]),
               tolerance = 1e-12)
})

test_that("nested designs and missing rows are recorded", {
  d <- expand.grid(lot = c("A", "B"), vial = c("1", "2"), rep = 1:3)
  d$value <- seq_len(nrow(d))
  d$value[2] <- NA
  observed <- bottle_anova(d, value ~ lot/vial)
  expect_equal(observed$design_type, "nested")
  expect_equal(observed$missing_rows, 2L)
  expect_equal(observed$n_complete, nrow(d) - 1L)
})

test_that("bottle ANOVA validates design variables and sample size", {
  expect_error(bottle_anova(data.frame(y = 1:3), y ~ group), "not found")
  expect_error(
    bottle_anova(data.frame(y = 1:2, group = c("a", "b")), y ~ group),
    "at least 3"
  )
})
