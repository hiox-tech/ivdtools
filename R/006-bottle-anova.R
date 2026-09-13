
#' bottle_anova constructor
#'
#' Create a bottle/batch ANOVA analysis object.
#' Supports one-way and nested designs, with automatic assumption checks,
#' missing value reporting, and ANOVA table output.
#'
#' @param data       A data frame
#' @param formula    A formula, e.g. \code{value ~ batch} or \code{value ~ batch/vial}
#' @param conf.level Confidence level (default 0.95)
#'
#' @return An S3 object of class \code{"bottle_anova"} with components:
#'   \item{call}{Function call}
#'   \item{data}{Original data}
#'   \item{formula}{Input formula}
#'   \item{conf.level}{Confidence level}
#'   \item{response_name}{Response variable name}
#'   \item{design_type}{Design type (\code{"one_way"} or \code{"nested"})}
#'   \item{rhs_vars}{Right-hand side grouping variable names}
#'   \item{term_labels}{Model term labels}
#'   \item{n_total}{Total sample size}
#'   \item{n_complete}{Number of complete cases}
#'   \item{missing_rows}{Indices of excluded missing rows}
#'   \item{cc_data}{Complete-case data frame}
#'   \item{aov_fit}{\code{aov} fitted object}
#'   \item{aov_table}{ANOVA table}
#'   \item{desc_stats}{Descriptive statistics table}
#'   \item{desc_group_col}{Grouping column label}
#'   \item{assumptions}{Assumption test results (Bartlett + Shapiro-Wilk)}
#'   \item{group_means}{Group means}
#'   \item{tukey_results}{Tukey HSD results (NULL initially, filled by \code{tukey()})}
#'
#' @importFrom stats aov as.formula bartlett.test lm na.omit qt
#'             residuals sd shapiro.test terms
#'             TukeyHSD
#' @export
#'
#' @examples
#'   # One-way design
#'   df1 <- data.frame(
#'     bottle = rep(c("A", "B", "C"), each = 5),
#'     value  = c(rnorm(5, 10, 1), rnorm(5, 12, 1), rnorm(5, 11, 1))
#'   )
#'   obj <- bottle_anova(df1, value ~ bottle)
#'
#'   # Nested design
#'   df2 <- data.frame(
#'     lot   = rep(c("L1", "L2"), each = 10),
#'     vial  = rep(c("V1", "V2"), each = 5, times = 2),
#'     value = c(rnorm(5, 10, 1), rnorm(5, 10.5, 1),
#'               rnorm(5, 12, 1), rnorm(5, 12.5, 1))
#'   )
#'   obj2 <- bottle_anova(df2, value ~ lot/vial)
bottle_anova <- function(data, formula, conf.level = 0.95) {

  # Internal helpers

  # Descriptive statistics
  .desc_stats <- function(x, conf.level = 0.95) {
    x <- x[is.finite(x)]
    n <- length(x)

    if (n < 2L) {
      return(list(n = n, mean = if (n > 0) mean(x) else NA_real_,
                  sd = NA_real_, sem = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_))
    }

    m <- mean(x)
    s <- sd(x)
    sem <- s / sqrt(n)
    alpha <- 1 - conf.level
    t_val <- qt(1 - alpha / 2, df = n - 1)

    list(n = n, mean = m, sd = s, sem = sem,
         ci_lower = m - t_val * sem,
         ci_upper = m + t_val * sem)
  }

  # Assumption tests: Bartlett (variance homogeneity) + Shapiro-Wilk (residual normality)
  .check_assumptions <- function(aov_fit, resp_name, rhs_vars, data) {
    shapiro_res <- tryCatch({
      r <- residuals(aov_fit)
      sw <- shapiro.test(r)
      list(statistic = sw$statistic, p.value = sw$p.value)
    }, error = function(e) NULL)

    bartlett_res <- tryCatch({
      if (length(rhs_vars) == 1L) {
        bt <- bartlett.test(data[[resp_name]], data[[rhs_vars[1L]]])
      } else {
        interact <- interaction(data[, rhs_vars, drop = FALSE], drop = TRUE)
        bt <- bartlett.test(data[[resp_name]], interact)
      }
      list(statistic = bt$statistic, parameter = bt$parameter, p.value = bt$p.value)
    }, error = function(e) NULL)

    list(bartlett = bartlett_res, shapiro = shapiro_res)
  }

  # Input validation

  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)

  if (!inherits(formula, "formula"))
    stop("'formula' must be a formula object.", call. = FALSE)

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single numeric value in (0, 1).", call. = FALSE)

  # Parse formula

  resp_name <- all.vars(formula[[2L]])
  if (length(resp_name) != 1L)
    stop("Formula must have exactly one response variable.", call. = FALSE)

  if (!resp_name %in% names(data))
    stop("Response variable '", resp_name, "' not found in data.", call. = FALSE)

  if (!is.numeric(data[[resp_name]]))
    stop("Response variable '", resp_name, "' must be numeric.", call. = FALSE)

  trms <- terms(formula)
  rhs_vars <- all.vars(formula[[3L]])

  bad_vars <- setdiff(rhs_vars, names(data))
  if (length(bad_vars) > 0)
    stop("Variable(s) on RHS not found in data: ",
         paste(bad_vars, collapse = ", "), call. = FALSE)

  # Detect nested design
  term_labels <- attr(trms, "term.labels")
  has_interaction <- any(grepl(":", term_labels))
  design_type <- if (has_interaction) "nested" else "one_way"

  model_formula <- formula

  # Extract data & handle missing values

  resp_vec <- data[[resp_name]]
  n_total <- nrow(data)

  miss_resp <- is.na(resp_vec) | !is.finite(resp_vec)
  miss_group <- rep(FALSE, n_total)
  for (v in rhs_vars) {
    miss_group <- miss_group | is.na(data[[v]])
  }
  miss_any <- miss_resp | miss_group
  all_ok <- !miss_any

  n_complete <- sum(all_ok)
  n_miss <- n_total - n_complete

  if (n_complete < 3L)
    stop("Need at least 3 complete cases for ANOVA, got ", n_complete, ".", call. = FALSE)

  cc_data <- data[all_ok, , drop = FALSE]

  for (v in rhs_vars) {
    cc_data[[v]] <- as.factor(cc_data[[v]])
  }

  # Check grouping factor levels
  for (v in rhs_vars) {
    n_levels <- nlevels(cc_data[[v]])
    if (n_levels < 2L)
      stop("Grouping variable '", v, "' must have at least 2 levels, got ",
           n_levels, ".", call. = FALSE)
  }

  if (design_type == "nested") {
    interact <- interaction(cc_data[, rhs_vars, drop = FALSE], drop = TRUE)
    if (nlevels(interact) < 2L)
      stop("Nested design: interaction of grouping variables must have at least 2 levels.",
           call. = FALSE)
  }

  # Fit ANOVA

  aov_fit <- aov(model_formula, data = cc_data)
  aov_summary <- summary(aov_fit)

  aov_tab <- aov_summary[[1L]]
  aov_terms <- trimws(rownames(aov_tab))

  # Descriptive statistics

  if (design_type == "one_way") {
    group_factor <- cc_data[[rhs_vars[1L]]]
    group_label <- rhs_vars[1L]
  } else {
    group_factor <- interaction(cc_data[, rhs_vars, drop = FALSE], drop = TRUE,
                                sep = ":")
    group_label <- paste(rhs_vars, collapse = ":")
  }

  group_levels <- levels(group_factor)
  desc_list <- lapply(group_levels, function(gl) {
    vals <- cc_data[[resp_name]][group_factor == gl]
    ds <- .desc_stats(vals, conf.level)
    data.frame(Group = gl, n = ds$n,
               Mean = ds$mean, SD = ds$sd,
               SEM = ds$sem,
               CI_lower = ds$ci_lower, CI_upper = ds$ci_upper,
               stringsAsFactors = FALSE)
  })
  desc_df <- do.call(rbind, desc_list)

  # Assumption tests

  assumptions <- .check_assumptions(aov_fit, resp_name, rhs_vars, cc_data)

  # ANOVA table with significance stars

  .sig_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    if (p < 0.1)   return(".")
    return(" ")
  }

  anova_df <- data.frame(
    Term     = aov_terms,
    Df       = aov_tab$Df,
    Sum_Sq   = aov_tab$`Sum Sq`,
    Mean_Sq  = aov_tab$`Mean Sq`,
    F_value  = if (!is.null(aov_tab$`F value`)) aov_tab$`F value` else NA_real_,
    p_value  = if (!is.null(aov_tab$`Pr(>F)`)) aov_tab$`Pr(>F)` else NA_real_,
    Sig      = sapply(
               if (!is.null(aov_tab$`Pr(>F)`)) aov_tab$`Pr(>F)` else NA_real_,
               .sig_stars),
    stringsAsFactors = FALSE
  )

  # Group means (for post-hoc CLD)

  group_means <- tapply(cc_data[[resp_name]], group_factor, mean, na.rm = TRUE)

  # Construct object

  out <- list(
    call            = match.call(),
    data            = data,
    formula         = formula,
    conf.level      = conf.level,

    response_name   = resp_name,
    design_type     = design_type,
    rhs_vars        = rhs_vars,
    term_labels     = term_labels,

    n_total         = n_total,
    n_complete      = n_complete,
    missing_rows    = if (n_miss > 0) which(miss_any) else integer(0),

    cc_data         = cc_data,
    aov_fit         = aov_fit,
    aov_table       = anova_df,

    desc_stats      = desc_df,
    desc_group_col  = group_label,

    assumptions     = assumptions,
    group_means     = group_means,

    tukey_results   = NULL
  )
  class(out) <- "bottle_anova"

  out
}

#' S3 print method: bottle_anova
#'
#' Print bottle ANOVA analysis results including descriptive statistics,
#' ANOVA table, and assumption checks.
#'
#' @param x   A \code{bottle_anova} object
#' @param ... Other arguments (unused, for S3 generic signature)
#'
#' @return Invisibly returns \code{x}.
#' @export
print.bottle_anova <- function(x, ...) {

  .format_num <- function(x, digits = 4L) {
    formatC(x, format = "f", digits = digits)
  }

  sep <- paste(rep("-", 60), collapse = "")

  # Header
  cat("\nBottle ANOVA Analysis\n")
  cat(sprintf("  Response:        %s\n", x$response_name))
  if (x$design_type == "nested") {
    cat(sprintf("  Group structure: %s (nested)\n",
                paste(x$rhs_vars, collapse = " / ")))
  } else {
    cat(sprintf("  Group variable:  %s (%d groups)\n",
                x$rhs_vars[1L], nlevels(x$cc_data[[x$rhs_vars[1L]]])))
  }
  cat(sprintf("  Complete cases:  %d / %d\n", x$n_complete, x$n_total))
  cat(sprintf("  Confidence level: %.2f\n", x$conf.level))

  if (length(x$missing_rows) > 0) {
    n_miss <- x$n_total - x$n_complete
    cat(sprintf("  Missing:         %d row(s) excluded\n", n_miss))
    cat(sprintf("    Row(s): %s\n", paste(x$missing_rows, collapse = ", ")))
  }

  # Descriptive statistics table
  cat("\n  Descriptive Statistics:\n")
  cat(" ", sep, "\n", sep = "")
  cat(sprintf("  %-15s %5s %9s %9s\n", x$desc_group_col, "n", "Mean", "SD"))
  cat(" ", sep, "\n", sep = "")
  for (i in seq_len(nrow(x$desc_stats))) {
    cat(sprintf("  %-15s %5d %9s %9s\n",
                x$desc_stats$Group[i], x$desc_stats$n[i],
                .format_num(x$desc_stats$Mean[i]),
                .format_num(x$desc_stats$SD[i])))
  }

  # ANOVA table
  cat("\n  ANOVA Table:\n")
  cat(" ", sep, "\n", sep = "")
  cat(sprintf("  %-15s %4s %9s %10s %9s %10s  %s\n",
              "Term", "Df", "Sum Sq", "Mean Sq", "F-value", "p-value", "Sig"))
  cat(" ", sep, "\n", sep = "")

  at <- x$aov_table
  for (i in seq_len(nrow(at))) {
    f_str <- if (!is.na(at$F_value[i])) .format_num(at$F_value[i]) else ""
    p_str <- if (!is.na(at$p_value[i])) .format_num(at$p_value[i]) else ""
    cat(sprintf("  %-15s %4d %9s %10s %9s %10s  %s\n",
                at$Term[i], at$Df[i],
                .format_num(at$Sum_Sq[i]),
                .format_num(at$Mean_Sq[i]),
                f_str, p_str, at$Sig[i]))
  }
  cat(" ", sep, "\n", sep = "")
  cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  # Assumption checks
  cat("\n  Assumption Checks:\n")
  ass <- x$assumptions
  if (!is.null(ass$bartlett)) {
    b <- ass$bartlett
    b_ok <- b$p.value >= 0.05
    cat(sprintf("    Bartlett test: K-squared = %.4f, df = %d, p-value = %.4f %s\n",
                b$statistic, b$parameter, b$p.value,
                if (b_ok) "(OK)" else "(WARNING: variances may be unequal)"))
  } else {
    cat("    Bartlett test: could not compute\n")
  }
  if (!is.null(ass$shapiro)) {
    s <- ass$shapiro
    s_ok <- s$p.value >= 0.05
    cat(sprintf("    Shapiro-Wilk:  W = %.4f, p-value = %.4f %s\n",
                s$statistic, s$p.value,
                if (s_ok) "(OK)" else "(WARNING: residuals may not be normal)"))
  } else {
    cat("    Shapiro-Wilk:  could not compute\n")
  }

  cat("\n  Post-hoc: use tukey() for pairwise comparisons\n\n")

  invisible(x)
}

#' S3 print method: tukey
#'
#' Print Tukey HSD post-hoc test results including pairwise comparison
#' table, significance markers, and compact letter display.
#'
#' @param x   A \code{tukey} object
#' @param ... Other arguments (unused, for S3 generic signature)
#'
#' @return Invisibly returns \code{x}.
#' @export
print.tukey <- function(x, ...) {

  .format_num <- function(x, digits = 4L) {
    formatC(x, format = "f", digits = digits)
  }

  .sig_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    if (p < 0.1)   return(".")
    return(" ")
  }

  sep <- paste(rep("-", 60), collapse = "")

  cat("\nTukey HSD Post-hoc Test\n")
  cat(sprintf("  Response:        %s\n", x$response_name))
  cat(sprintf("  Complete cases:  %d / %d\n", x$n_complete, x$n_total))
  cat(sprintf("  Confidence level: %.2f\n", x$conf.level))

  for (tm in names(x$tukey_list)) {
    entry <- x$tukey_list[[tm]]

    cat("\n  Term:", tm, "\n")
    cat(" ", sep, "\n", sep = "")

    cat(sprintf("  %-20s %9s %9s %9s %9s  %s\n",
                "Comparison", "Diff", "lwr", "upr", "p-value", "Sig"))
    cat(" ", sep, "\n", sep = "")
    tm_mat <- entry$matrix
    for (i in seq_len(nrow(tm_mat))) {
      cmp <- rownames(tm_mat)[i]
      p_val <- tm_mat[i, "p adj"]
      cat(sprintf("  %-20s %9s %9s %9s %9s  %s\n",
                  cmp,
                  .format_num(tm_mat[i, "diff"]),
                  .format_num(tm_mat[i, "lwr"]),
                  .format_num(tm_mat[i, "upr"]),
                  .format_num(p_val),
                  if (!is.na(p_val) && p_val < 1 - entry$conf.level) .sig_stars(p_val) else ""))
    }
    cat(" ", sep, "\n", sep = "")
    cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

    cat("\n  Compact Letter Display:\n")
    cat(" ", sep, "\n", sep = "")
    sorted_groups <- names(sort(entry$means, decreasing = TRUE))
    for (g in sorted_groups) {
      cat(sprintf("    %s: %s\n", g, entry$cld[g]))
    }
    cat("\n")
  }

  invisible(x)
}

#' Tukey HSD post-hoc test (S3 method)
#'
#' Perform Tukey HSD pairwise comparisons on an ANOVA result,
#' with compact letter display.
#'
#' @param x         A \code{bottle_anova} object
#' @param term      Effect term to analyse (character vector); \code{NULL} for all non-residual terms
#' @param conf.level Confidence level; \code{NULL} uses the value stored in the object
#' @param ...       Other arguments (unused)
#'
#' @return An S3 object of class \code{"tukey"} with components:
#'   \item{response_name}{Response variable name}
#'   \item{conf.level}{Confidence level}
#'   \item{tukey_list}{List, one element per term, each with matrix, means, cld}
#'
#' Note: The original \code{bottle_anova} object's \code{$tukey_results} is also updated.
#' @export
#'
#' @examples
#'   obj <- bottle_anova(ivd_bottle_example, value ~ bottle)
#'   res <- tukey(obj)
#'   res
tukey.bottle_anova <- function(x, term = NULL, conf.level = NULL, ...) {

  if (!inherits(x, "bottle_anova"))
    stop("Input must be a 'bottle_anova' object.", call. = FALSE)

  # Internal helpers

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Compact Letter Display
  # Generate CLD from TukeyHSD results without depending on multcompView.
  .compact_letters <- function(tukey_mat, group_means, alpha = 0.05) {
    comps <- rownames(tukey_mat)
    all_groups <- unique(unlist(strsplit(comps, "-")))

    groups <- names(sort(group_means[all_groups], decreasing = TRUE))

    p_mat <- matrix(NA, nrow = length(groups), ncol = length(groups),
                    dimnames = list(groups, groups))
    diag(p_mat) <- 1

    for (cmp in comps) {
      gs <- strsplit(cmp, "-")[[1]]
      p_val <- tukey_mat[cmp, "p adj"]
      p_mat[gs[1], gs[2]] <- p_mat[gs[2], gs[1]] <- p_val
    }

    # Round 1: initial assignment
    group_letters <- setNames(character(length(groups)), groups)
    used_letters <- character(0)

    for (i in seq_along(groups)) {
      gi <- groups[i]
      qualifies <- character(0)

      for (lid in seq_along(used_letters)) {
        l <- used_letters[lid]
        groups_with_l <- names(group_letters[grepl(l, group_letters)])

        all_ok <- TRUE
        for (g_with in groups_with_l) {
          if (!is.na(p_mat[gi, g_with]) && p_mat[gi, g_with] < alpha) {
            all_ok <- FALSE
            break
          }
        }
        if (all_ok) qualifies <- c(qualifies, l)
      }

      if (length(qualifies) > 0) {
        group_letters[gi] <- paste(qualifies, collapse = "")
      } else {
        new_letter <- letters[length(used_letters) + 1L]
        used_letters <- c(used_letters, new_letter)
        group_letters[gi] <- new_letter
      }
    }

    # Round 2: propagate/merge
    changed <- TRUE
    while (changed) {
      changed <- FALSE
      for (i in seq_along(groups)) {
        gi <- groups[i]
        for (j in seq_along(groups)) {
          if (i == j) next
          gj <- groups[j]
          if (is.na(p_mat[gi, gj]) || p_mat[gi, gj] >= alpha) {
            for (k in seq_len(nchar(group_letters[gi]))) {
              letter <- substr(group_letters[gi], k, k)
              if (!grepl(letter, group_letters[gj], fixed = TRUE)) {
                groups_with_letter <- names(group_letters[grepl(letter, group_letters)])
                has_conflict <- FALSE
                for (gw in groups_with_letter) {
                  if (gw != gj && !is.na(p_mat[gj, gw]) && p_mat[gj, gw] < alpha) {
                    has_conflict <- TRUE
                    break
                  }
                }
                if (!has_conflict) {
                  group_letters[gj] <- paste0(group_letters[gj], letter)
                  changed <- TRUE
                }
              }
            }
          }
        }
      }
    }

    # Sort and deduplicate letters within each group
    for (g in groups) {
      chars <- unique(strsplit(group_letters[g], "")[[1]])
      group_letters[g] <- paste(sort(chars), collapse = "")
    }

    group_letters
  }

  # Parameter handling

  cl <- conf.level %||% x$conf.level

  if (!is.numeric(cl) || length(cl) != 1L || cl <= 0 || cl >= 1)
    stop("'conf.level' must be a single numeric value in (0, 1).", call. = FALSE)

  all_terms <- setdiff(x$term_labels, "Residuals")

  if (is.null(term)) {
    terms_to_do <- all_terms
  } else {
    bad <- setdiff(term, all_terms)
    if (length(bad) > 0)
      stop("Term(s) not found in ANOVA model: ", paste(bad, collapse = ", "),
           call. = FALSE)
    terms_to_do <- term
  }

  if (length(terms_to_do) < 1L)
    stop("No terms available for Tukey HSD.", call. = FALSE)

  # Run TukeyHSD for each term

  tukey_list <- list()

  for (tm in terms_to_do) {
    thsd <- TukeyHSD(x$aov_fit, which = tm, conf.level = cl)
    tukey_mat <- thsd[[tm]]

    if (tm %in% x$rhs_vars) {
      gm <- tapply(x$cc_data[[x$response_name]],
                   x$cc_data[[tm]], mean, na.rm = TRUE)
    } else {
      vars_in <- all.vars(as.formula(paste("~", tm)))
      gm <- tapply(x$cc_data[[x$response_name]],
                   interaction(x$cc_data[, vars_in, drop = FALSE],
                               drop = TRUE, sep = ":"),
                   mean, na.rm = TRUE)
    }

    cld <- .compact_letters(tukey_mat, gm, alpha = 1 - cl)

    tukey_list[[tm]] <- list(
      term       = tm,
      matrix     = tukey_mat,
      means      = gm,
      cld        = cld,
      conf.level = cl
    )
  }

  # Sync back to original object
  x$tukey_results <- tukey_list

  # Build standalone tukey object
  result <- list(
    call          = match.call(),
    response_name = x$response_name,
    conf.level    = cl,
    formula       = x$formula,
    n_total       = x$n_total,
    n_complete    = x$n_complete,
    tukey_list    = tukey_list
  )
  class(result) <- "tukey"
  result
}
