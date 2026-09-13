# fit.R — General equation fitting tool
# Supported: common response curves (Linear/Exponential/4PLC/5PLC etc.),
#     user-defined arbitrary formulas.
#     Automatic or manual initial parameters.
#     Weighted fitting, parameter constraints (box/inequality/equality).
#
# Dependencies: minpack.lm, nloptr, nls2, ggplot2

# Print a data.frame as a plain-text block (cat + .print_df combined into one block)
# so knitr/qmd does not split table from surrounding text.
#' @noRd
.print_block <- function(title, tbl) {
  txt <- paste0(title, "\n")
  df_lines <- utils::capture.output(print(tbl, row.names = FALSE))
  if (length(df_lines) > 0) {
    nc <- nchar(df_lines[1L])
    if (is.na(nc)) nc <- 0L
    txt <- paste0(txt, df_lines[1L], "\n")
    if (nc > 0L)
      txt <- paste0(txt, paste(rep("-", nc), collapse = ""), "\n")
    txt <- paste0(txt, paste(df_lines[-1L], collapse = "\n"), "\n")
  }
  cat(txt)
  invisible(tbl)
}

# 0. Replicate summary

#' Summarize replicate measurements
#'
#' Groups data by x, computes mean, SD, and sample size for each group.
#' When `weights` is specified, returns an additional `weight` column.
#' Summary metadata (call, weighting method, column names) is stored as
#' attributes for downstream reference.
#'
#' @param data     a data frame
#' @param x        x column name (grouping variable)
#' @param y        One or more y column names (response variables). A single
#'   response retains the historical `y_mean`, `y_sd`, and `y_n` column names;
#'   multiple responses use `<name>_mean`, `<name>_sd`, and `<name>_n`.
#' @param weights  weighting method: NULL / "n" / "1/sd" / "1/sd^2"
#'                   NULL     -- summary only, no weight column
#'                   "n"      -- w = number of replicates
#'                   "1/sd"   -- w = 1/sd(y)
#'                   "1/sd^2" -- w = 1/var(y) (inverse-variance weighting)
#' @param na.rm    remove NA? (default TRUE)
#'
#' @return A data.frame with columns:
#'   \item{x}{grouping variable}
#'   \item{y_mean, y_sd, y_n}{Historical names used for one response.}
#'   \item{<name>_mean, <name>_sd, <name>_n}{Names used for each response when
#'   two or more response columns are supplied.}
#'   \item{weight}{(only when weights is specified) weight values}
#'   Attributes: call, weights (method name), x_col, y_col
#' @examples
#' # Basic summary
#' df <- data.frame(
#'   conc = rep(c(1, 2, 5, 10), each = 3),
#'   resp = c(3.1, 3.2, 2.9, 5.8, 6.1, 5.9,
#'            14.2, 14.5, 14.0, 28.1, 27.9, 28.3)
#' )
#' rep <- replicate_to_mean(df, "conc", "resp")
#' rep
#'
#' # Inverse-variance weighting
#' rep_w <- replicate_to_mean(df, "conc", "resp", weights = "1/sd^2")
#' rep_w
#' @export
replicate_to_mean <- function(data, x, y,
                     weights = NULL,
                     na.rm = TRUE) {
  cl <- match.call()

  if (!is.character(x) || length(x) != 1L || !x %in% names(data))
    stop("Column '", x, "' not found in data.")
  if (!is.character(y) || !length(y) || any(!y %in% names(data)))
    stop("Response column(s) not found in data: ",
         paste(setdiff(y, names(data)), collapse = ", "), ".")
  if (anyDuplicated(y))
    stop("'y' must contain unique response column names.", call. = FALSE)
  if (length(y) > 1L && !is.null(weights))
    stop("'weights' is available only when a single response column is supplied.",
         call. = FALSE)

  xv <- data[[x]]
  keep_x <- if (na.rm) !is.na(xv) else rep(TRUE, length(xv))
  x_groups <- split(seq_along(xv)[keep_x], xv[keep_x])
  x_vals <- if (is.numeric(xv)) as.numeric(names(x_groups)) else names(x_groups)

  summarize_response <- function(response) {
    yv <- data[[response]]
    if (!is.numeric(yv))
      stop("Response column '", response, "' must be numeric.", call. = FALSE)
    stats <- lapply(x_groups, function(ii) {
      values <- yv[ii]
      if (na.rm) values <- values[!is.na(values)]
      c(mean = mean(values, na.rm = na.rm),
        sd = stats::sd(values, na.rm = na.rm), n = length(values))
    })
    stats <- do.call(rbind, stats)
    y_mean <- stats[, "mean"]
    y_sd <- stats[, "sd"]
    y_n <- as.integer(stats[, "n"])

    # Preserve the historical treatment of singleton/zero-SD groups, which is
    # needed by the optional weighting modes. The original SDs are otherwise
    # untouched and therefore remain suitable for replicate-lambda estimates.
    replace_sd <- is.na(y_sd) | y_sd == 0
    if (length(y) == 1L && any(replace_sd)) {
      sd_ok <- !replace_sd & y_n > 1L
      pooled_sd <- if (any(sd_ok)) {
        sqrt(sum((y_n[sd_ok] - 1) * y_sd[sd_ok]^2) /
             sum(y_n[sd_ok] - 1))
      } else {
        values <- yv[keep_x & (if (na.rm) !is.na(yv) else TRUE)]
        spread <- diff(range(values, na.rm = TRUE))
        if (!is.finite(spread) || spread == 0) NA_real_
        else spread / (2 * sqrt(2 * length(values)))
      }
      y_sd[replace_sd] <- pooled_sd
    }
    list(mean = y_mean, sd = y_sd, n = y_n)
  }

  summaries <- lapply(y, summarize_response)
  result <- data.frame(x = x_vals, row.names = NULL, check.names = FALSE)
  if (length(y) == 1L) {
    result$y_mean <- summaries[[1L]]$mean
    result$y_sd <- summaries[[1L]]$sd
    result$y_n <- summaries[[1L]]$n
  } else {
    for (i in seq_along(y)) {
      result[[paste0(y[i], "_mean")]] <- summaries[[i]]$mean
      result[[paste0(y[i], "_sd")]] <- summaries[[i]]$sd
      result[[paste0(y[i], "_n")]] <- summaries[[i]]$n
    }
  }

  # weights
  wt <- NULL
  if (!is.null(weights)) {
    weights <- match.arg(weights, c("1/sd^2", "1/sd", "n"))
    # if all y_sd are still NA/0 (all data identical), fall back to equal weights
    y_sd <- summaries[[1L]]$sd
    y_n <- summaries[[1L]]$n
    if (all(is.na(y_sd) | y_sd == 0)) {
      warning("All groups have identical y values; cannot compute sd-based weights. Using equal weights.")
      wt <- rep(1, length(y_n))
    } else {
      wt <- switch(weights,
        "n"      = y_n,
        "1/sd"   = 1 / y_sd,
        "1/sd^2" = 1 / y_sd^2
      )
    }
  }

  if (!is.null(weights)) {
    result$weight <- wt
  }
  attr(result, "call")     <- cl
  attr(result, "weights")  <- weights
  attr(result, "x_col")    <- x
  attr(result, "y_col")    <- y
  result
}

# 1. Equation registry

#' Return the complete equation registry (data.frame)
#' @noRd
eq_registry <- function() {
  data.frame(
    id = c(
      # Response curves
      "E01", "E02", "E03", "E04", "E05", "E06",
      "E07", "E08", "E09", "E10", "E11", "E12"
    ),
    category = c(
      rep("Response", 12)
    ),
    name = c(
      "Linear", "Quadratic", "ExpoGrowth", "ExpoDecay", "Power", "Log",
      "4PLC", "5PLC", "Gaussian", "Michaelis", "Sinusoidal", "Cubic"
    ),
    formula_str = c(
      "a + b*x",
      "a + b*x + c*x^2",
      "a * exp(b * x)",
      "a * exp(-b * x)",
      "a * x^b",
      "a + b * log(x)",
      "A + (B - A) / (1 + (x / C)^D)",
      "A + (B - A) / (1 + (x / C)^D)^E",
      "a * exp(-(x - b)^2 / (2 * c^2))",
      "Vmax * x / (Km + x)",
      "a * sin(b * x + c) + d",
      "a + b*x + c*x^2 + d*x^3"
    ),
    params = c(
      "a, b", "a, b, c", "a, b", "a, b", "a, b", "a, b",
      "A, B, C, D", "A, B, C, D, E", "a, b, c", "Vmax, Km", "a, b, c, d",
      "a, b, c, d"
    ),
    n_params = c(2L, 3L, 2L, 2L, 2L, 2L,
                 4L, 5L, 3L, 2L, 4L, 4L),
    engine = c(
      "lm", "lm", "nls", "nls", "nls", "lm",
      "nls", "nls", "nls", "nls", "nls", "lm"
    ),
    description = c(
      "y = a + bx",
      "y = a + bx + cx^2",
      "y = a * exp(bx)",
      "y = a * exp(-bx)",
      "y = a * x^b",
      "y = a + b * ln(x)",
      "y = A + (B-A)/(1+(x/C)^D)  [4PLC]",
      "y = A + (B-A)/(1+(x/C)^D)^E  [5PLC]",
      "y = a * exp(-(x-b)^2/(2c^2))  [Gaussian]",
      "y = Vmax * x/(Km+x)  [Michaelis-Menten]",
      "y = a * sin(bx+c)+d",
      "y = a + bx + cx^2 + dx^3"
    )
  )
}

# 2. Print equation registry

#' Print the equation registry
#'
#' @param category filter category: "Response", or NULL (all)
#' @param engine   filter engine: "lm", "nls", or NULL (all)
#' @return invisible data.frame (matching rows), prints a formatted table to the console
#' @examples
#' # All equations
#' list_equation()
#'
#' # Response only
#' list_equation(category = "Response")
#'
#' # lm engine only
#' list_equation(engine = "lm")
#' @export
list_equation <- function(category = NULL, engine = NULL) {
  reg <- eq_registry()

  if (!is.null(category))
    reg <- reg[reg$category == category, ]
  if (!is.null(engine))
    reg <- reg[reg$engine == engine, ]

  if (nrow(reg) == 0) {
    cat("No matching equations.\n")
    return(invisible(reg))
  }

  # format table (use formula instead of long description, drop category for one-line-per-entry)
  tbl <- reg[, c("id", "name", "formula_str", "params", "n_params", "engine")]
  names(tbl) <- c("ID", "Name", "Formula", "Params", "nP", "Engine")

  # Print using the user's current console width.
  .print_df(tbl)

  invisible(reg)
}

# 3. Automatic start values

#' Estimate initial parameter values automatically
#'
#' @param eq_id  equation ID (e.g. "E07")
#' @param x,y   numeric vectors
#' @return named list (start list)
#' @noRd
autostart <- function(eq_id, x, y) {
  switch(eq_id,
    # Linear
    "E01" = {
      fit0 <- lm(y ~ x)
      list(a = unname(coef(fit0)[1]), b = unname(coef(fit0)[2]))
    },
    # Quadratic
    "E02" = {
      fit0 <- lm(y ~ x + I(x^2))
      list(a = unname(coef(fit0)[1]),
           b = unname(coef(fit0)[2]),
           c = unname(coef(fit0)[3]))
    },
    # Cubic
    "E12" = {
      fit0 <- lm(y ~ x + I(x^2) + I(x^3))
      list(a = unname(coef(fit0)[1]),
           b = unname(coef(fit0)[2]),
           c = unname(coef(fit0)[3]),
           d = unname(coef(fit0)[4]))
    },
    # Exponential Growth
    "E03" = {
      # log(y) = log(a) + b*x
      ok <- y > 0
      if (sum(ok) < 3) return(list(a = 1, b = 1))
      fit0 <- lm(log(y[ok]) ~ x[ok])
      list(a = exp(unname(coef(fit0)[1])),
           b = unname(coef(fit0)[2]))
    },
    # Exponential Decay
    # formula: y = a * exp(-b * x) ⇒ log(y) = log(a) - b*x
    # so b = -lm(log(y)~x)_slope
    "E04" = {
      y_safe <- pmax(y, min(y[y > 0], na.rm = TRUE) * 0.1, 1e-6)
      fit0 <- lm(log(y_safe) ~ x)
      b_est <- -unname(coef(fit0)[2])
      if (b_est <= 0 || is.na(b_est)) b_est <- 0.1
      list(a = exp(unname(coef(fit0)[1])),
           b = b_est)
    },
    # Power Law
    "E05" = {
      ok <- y > 0 & x > 0
      if (sum(ok) < 3) return(list(a = 1, b = 1))
      fit0 <- lm(log(y[ok]) ~ log(x[ok]))
      list(a = exp(unname(coef(fit0)[1])),
           b = unname(coef(fit0)[2]))
    },
    # Logarithmic
    "E06" = {
      ok <- x > 0
      if (sum(ok) < 3) return(list(a = mean(y), b = 0))
      fit0 <- lm(y[ok] ~ log(x[ok]))
      list(a = unname(coef(fit0)[1]),
           b = unname(coef(fit0)[2]))
    },
    # 4PLC
    "E07" = {
      A <- min(y, na.rm = TRUE)
      B <- max(y, na.rm = TRUE)
      half <- (A + B) / 2
      # C is x at half response (linear interpolation)
      idx <- which.min(abs(y - half))
      C <- x[idx]
      # D sign: D>0 descending (B→A), D<0 ascending (A→B)
      D_sign <- if (cor(x, y, use = "complete.obs") < 0) 1 else -1
      list(A = A, B = B, C = C, D = D_sign)
    },
    # 5PLC
    "E08" = {
      st <- autostart("E07", x, y)
      c(st, list(E = 1))
    },
    # Gaussian
    "E09" = {
      a <- max(y, na.rm = TRUE)
      b <- x[which.max(y)]
      c <- stats::sd(x, na.rm = TRUE)
      list(a = a, b = b, c = c)
    },
    # Michaelis-Menten
    "E10" = {
      list(Vmax = max(y, na.rm = TRUE),
           Km   = median(x, na.rm = TRUE))
    },
    # Sinusoidal
    "E11" = {
      d <- mean(y, na.rm = TRUE)
      y_detrend <- y - d
      a <- (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)) / 2
      # estimate period via zero-crossings
      sign_changes <- which(diff(sign(y_detrend)) != 0)
      period_guess <- if (length(sign_changes) >= 4) {
        # at least 2 full periods: mean zero-crossing spacing × 2
        2 * mean(diff(x[sign_changes]))
      } else if (length(sign_changes) >= 2) {
        mean(diff(x[sign_changes]))
      } else {
        diff(range(x)) / 2
      }
      b <- 2 * pi / max(abs(period_guess), 1e-6)
      c <- 0
      list(a = a, b = b, c = c, d = d)
    },
    # fallback
    list(a = 1, b = 1)
  )
}

# 4. Response curve fitting (internal) — core nls / lm

#' Fit a response curve model
#' @noRd
fit_response <- function(eq_id, reg_row, data, x, y,
                         start = NULL, engine = NULL,
                         lower = NULL, upper = NULL,
                         constraints = NULL,
                         weights = NULL, ...) {
  if (is.null(engine)) engine <- reg_row$engine

  # rename columns so the formula can use "x" and "y"
  d <- data
  names(d)[names(d) == x] <- "x"
  names(d)[names(d) == y] <- "y"

  # parameter constraints → constrained engine (priority over lm/nls)
  if (!is.null(lower) || !is.null(upper) || !is.null(constraints)) {
    # build formula: constrained fitting needs nls style (explicit parameter names)
    frm <- if (engine == "lm") {
      # map lm equations to explicit parameter formulas, matching autostart return value
      switch(eq_id,
        "E01" = stats::as.formula("y ~ a + b*x"),
        "E02" = stats::as.formula("y ~ a + b*x + c*I(x^2)"),
        "E06" = stats::as.formula("y ~ a + b*log(x)"),
        "E12" = stats::as.formula("y ~ a + b*x + c*I(x^2) + d*I(x^3)"),
        stats::as.formula("y ~ x")
      )
    } else {
      rhs <- reg_row$formula_str
      stats::as.formula(paste("y ~", rhs))
    }
    if (is.null(start)) {
      start <- autostart(eq_id, d[["x"]], d[["y"]])
    }
    return(fit_constrained(frm = frm, data = d, start = start,
                            lower = lower, upper = upper,
                            constraints = constraints,
                            weights = weights, ...))
  }

  # lm engine
  if (engine == "lm") {
    frm <- switch(eq_id,
      "E01" = stats::as.formula("y ~ x"),
      "E02" = stats::as.formula("y ~ x + I(x^2)"),
      "E06" = stats::as.formula("y ~ log(x)"),
      "E12" = stats::as.formula("y ~ x + I(x^2) + I(x^3)"),
      stats::as.formula("y ~ x")
    )
    fit <- stats::lm(frm, data = d, weights = weights)
    # rename lm coefficients to match registry parameter names
    lm_rename <- switch(eq_id,
      "E01" = stats::setNames(c("a", "b"), names(stats::coef(fit))),
      "E02" = stats::setNames(c("a", "b", "c"), names(stats::coef(fit))),
      "E06" = stats::setNames(c("a", "b"), names(stats::coef(fit))),
      "E12" = stats::setNames(c("a", "b", "c", "d"), names(stats::coef(fit))),
      NULL
    )
    if (!is.null(lm_rename)) {
      names(fit$coefficients) <- lm_rename
    }
    class(fit) <- c("lm_eq", class(fit))
    return(fit)
  }

  # nls engine
  rhs <- reg_row$formula_str
  frm <- stats::as.formula(paste("y ~", rhs))

  # explicitly remove NA (nls does not auto-handle), sync weights
  ok <- stats::complete.cases(d[["x"]], d[["y"]])
  d <- d[ok, , drop = FALSE]
  if (!is.null(weights)) weights <- weights[ok]

  if (is.null(start)) {
    start <- autostart(eq_id, d[["x"]], d[["y"]])
  }

  if (nrow(d) < length(start) + 1) {
    stop("Too few complete observations (", nrow(d),
         ") for model with ", length(start), " parameters.")
  }

  # unconstrained → nlsLM
  .nls_args <- c(list(formula = frm, data = d, start = start),
                 list(...))
  if (!is.null(weights)) .nls_args$weights <- weights
  fit <- tryCatch(
    do.call(minpack.lm::nlsLM, .nls_args),
    error = function(e) {
      stop("nlsLM fit failed: ", e$message, "\n",
           "  Suggestions: (1) Check the formula  (2) Provide better start values  (3) Add lower/upper bounds")
    }
  )
  class(fit) <- c("nls_eq", class(fit))
  fit
}

# 6. Custom equation fitting

#' Fit a user-defined custom equation
#' @noRd
fit_custom <- function(formula_expr, data, x, y,
                       start = NULL,
                       lower = NULL, upper = NULL,
                       constraints = NULL,
                       weights = NULL, ...) {
  # rename columns so the formula can use "x" and "y" (consistent with fit_response)
  d <- data
  names(d)[names(d) == x] <- "x"
  names(d)[names(d) == y] <- "y"

  # formula_expr can be string or formula
  # unify: extract formula text, replace original column names with "x"/"y" (to match renamed data)
  if (is.character(formula_expr)) {
    # string: if no "~", prepend LHS
    if (grepl("~", formula_expr)) {
      frm_text <- formula_expr
    } else {
      frm_text <- paste("y ~", formula_expr)
    }
  } else if (inherits(formula_expr, "formula")) {
    frm_text <- deparse(formula_expr)
  } else {
    stop("formula_expr must be a character string or a formula object.")
  }
  # replace original column names with "x"/"y" in formula
  if (y != "y") frm_text <- gsub(paste0("\\b", y, "\\b"), "y", frm_text)
  if (x != "x") frm_text <- gsub(paste0("\\b", x, "\\b"), "x", frm_text)
  frm <- stats::as.formula(frm_text)

  # extract parameter names (symbols not x or y on RHS)
  all_vars <- all.vars(frm)
  free_vars <- setdiff(all_vars, c("x", "y"))

  if (length(free_vars) == 0) {
    return(stats::lm(frm, data = d, weights = weights))
  }

  # explicitly remove NA (nls does not auto-handle), sync weights
  ok <- stats::complete.cases(d[["x"]], d[["y"]])
  d <- d[ok, , drop = FALSE]
  if (!is.null(weights)) weights <- weights[ok]

  if (nrow(d) < length(free_vars) + 1) {
    stop("Too few complete observations (", nrow(d),
         ") for custom model with ", length(free_vars), " parameters.")
  }

  # attempt automatic start values (inline)
  if (is.null(start)) {
    .autostart_custom_nls2 <- function(frm, data, free_vars) {
      lower <- list(); upper <- list()
      y <- data[["y"]]
      for (v in free_vars) {
        yrange <- diff(range(y, na.rm = TRUE))
        if (yrange == 0) yrange <- 1
        lower[[v]] <- -abs(yrange) * 10
        upper[[v]] <-  abs(yrange) * 10
      }
      start_grid <- expand.grid(lapply(free_vars, function(v) c(lower[[v]], 0, upper[[v]])))
      names(start_grid) <- free_vars
      grid <- nls2::nls2(frm, data = data, start = start_grid,
                         algorithm = "brute-force",
                         control = stats::nls.control(warnOnly = TRUE))
      as.list(stats::coef(grid))
    }
    start <- .autostart_custom_nls2(frm, d, free_vars)
  }

  # parameter constraints → constrained engine
  if (!is.null(lower) || !is.null(upper) || !is.null(constraints)) {
    return(fit_constrained(frm, d, start, lower, upper, constraints,
                           weights = weights, ...))
  }

  # unconstrained → nlsLM
  .nls_args <- c(list(formula = frm, data = d, start = start),
                 list(...))
  if (!is.null(weights)) .nls_args$weights <- weights
  fit <- tryCatch(
    do.call(minpack.lm::nlsLM, .nls_args),
    error = function(e) {
      stop("nlsLM fit failed: ", e$message, "\n",
           "  Suggestions: (1) Check the formula  (2) Provide better start values")
    }
  )
  class(fit) <- c("nls_eq", class(fit))
  fit
}



# 6b. Constrained fitting engine

#' General-purpose constrained parameter fitting (via nloptr or nlsLM)
#'
#' @param frm        formula (e.g. y ~ a * exp(-b * x))
#' @param data       data frame (must contain x, y columns)
#' @param start      named list of starting values
#' @param lower      named list/vector of lower bounds (NULL = none)
#' @param upper      named list/vector of upper bounds (NULL = none)
#' @param constraints list of constraint functions: list(ineq = fun, eq = fun)
#'                    ineq(par) returns a vector; all values must be non-positive
#'                    eq(par) returns a vector; must satisfy all(eq = 0)
#' @param ...        additional arguments passed to the solver
#' @return list with coefficients, fitted.values, residuals, df.residual, rss, conv
#' @noRd
fit_constrained <- function(frm, data, start, lower = NULL,
                             upper = NULL, constraints = NULL,
                             weights = NULL, ...) {
  dots <- list(...)
  rhs <- frm[[3]]
  lhs_name <- all.vars(frm[[2]])
  par_names <- names(start)
  w <- weights  # already preprocessed in fit_equation

  # build weighted RSS objective
  rss_fun <- function(par) {
    names(par) <- par_names
    ee <- list2env(as.list(par), parent = baseenv())
    ee[[lhs_name]] <- data[[lhs_name]]
    for (nm in names(data)) {
      if (!nm %in% par_names) ee[[nm]] <- data[[nm]]
    }
    y_pred <- tryCatch(eval(rhs, envir = ee), error = function(e) {
      rep(NA_real_, nrow(data))
    })
    y_obs <- data[[lhs_name]]
    resid2 <- (y_obs - y_pred)^2
    # handle NA element-wise: single NA contributes a large penalty instead of column-wide Inf
    resid2[is.na(resid2) | is.infinite(resid2)] <- .Machine$double.xmax / 1e10
    if (!is.null(w)) sum(w * resid2) else sum(resid2)
  }

  # check if constraints are needed
  has_bounds <- !is.null(lower) || !is.null(upper)
  has_cons <- !is.null(constraints)

  if (!has_bounds && !has_cons) {
    stop("fit_constrained called without constraints; use fit_response or nlsLM directly.")
  }

  # coerce lower/upper to named vectors
  .to_named_vec <- function(x, default) {
    if (is.null(x)) return(setNames(rep(default, length(start)), par_names))
    if (is.list(x)) x <- unlist(x)
    if (is.null(names(x))) {
      # unnamed scalar → broadcast to all parameters
      if (length(x) == 1L) {
        x <- setNames(rep(x, length(par_names)), par_names)
      } else if (length(x) != length(par_names)) {
        stop("Unnamed bound vector has length ", length(x),
             " but there are ", length(par_names), " parameters (",
             paste(par_names, collapse = ", "), ").",
             "\n  Provide a named list/vector (e.g. list(a = 0, b = 1))",
             " or a single value to apply to all parameters.")
      } else {
        names(x) <- par_names
      }
    }
    # check for unknown parameter names
    unknown <- setdiff(names(x), par_names)
    if (length(unknown) > 0) {
      stop("Unknown parameter name(s) in bounds: ",
           paste(unknown, collapse = ", "),
           ". Valid parameters: ", paste(par_names, collapse = ", "))
    }
    # fill in missing parameter names
    for (nm in par_names) {
      if (!nm %in% names(x)) x[[nm]] <- default
    }
    x[par_names]
  }
  lo <- .to_named_vec(lower, -Inf)
  hi <- .to_named_vec(upper, Inf)

  # 1) box constraints only → nlsLM (preferred)
  if (!has_cons) {
    fit <- tryCatch(
      minpack.lm::nlsLM(frm, data = data, start = start,
                         lower = lo, upper = hi, weights = w, ...),
      error = function(e) NULL
    )
    if (!is.null(fit)) return(fit)
  }

  # 2) constraints present or nlsLM failed → nloptr
    # build constraint functions
    eval_g_ineq <- NULL
    eval_g_eq <- NULL

    if (!is.null(constraints$ineq)) {
      eval_g_ineq <- function(x) {
        names(x) <- par_names
        constraints$ineq(x)
      }
    }
    if (!is.null(constraints$eq)) {
      eval_g_eq <- function(x) {
        names(x) <- par_names
        constraints$eq(x)
      }
    }

    start_vec <- unlist(start)
    par_names <- names(start_vec)

    # attempt nloptr solve with multi-level retry
    .try_nloptr <- function(xtol, maxeval, constr_tol, local_xtol, local_maxeval) {
      # select algorithm
      algo <- if (has_cons) {
        if (!is.null(constraints$ineq) && is.null(constraints$eq)) {
          "NLOPT_LN_COBYLA"
        } else if (!is.null(constraints$eq)) {
          "NLOPT_LN_AUGLAG"
        } else {
          "NLOPT_LN_COBYLA"
        }
      } else {
        "NLOPT_LN_BOBYQA"
      }

      opts <- list(
        algorithm = algo,
        xtol_rel = xtol,
        maxeval = maxeval,
        constr_tol_abs = constr_tol,
        constr_tol_rel = constr_tol,
        print_level = 0
      )
      # user can override opts via ...
      for (nm in names(dots)) {
        if (nm %in% names(opts)) opts[[nm]] <- dots[[nm]]
      }

      # AUGLAG requires a local solver
      if (algo == "NLOPT_LN_AUGLAG") {
        local_algo <- if (is.null(constraints$ineq)) "NLOPT_LN_BOBYQA" else "NLOPT_LN_COBYLA"
        opts$local_opts <- list(
          algorithm = local_algo,
          xtol_rel = local_xtol,
          maxeval = local_maxeval
        )
      }

      tryCatch(
        nloptr::nloptr(
          x0 = start_vec,
          eval_f = rss_fun,
          lb = lo,
          ub = hi,
          eval_g_ineq = eval_g_ineq,
          eval_g_eq = eval_g_eq,
          opts = opts
        ),
        error = function(e) NULL
      )
    }

    # round 1: tight tolerance
    fit <- .try_nloptr(xtol = 1e-7, maxeval = 10000,
                       constr_tol = 1e-6,
                       local_xtol = 1e-7, local_maxeval = 5000)

    # if failed, retry: relax tolerance + increase iterations
    if (is.null(fit) || fit$status < 0) {
      fit <- .try_nloptr(xtol = 1e-6, maxeval = 20000,
                         constr_tol = 1e-4,
                         local_xtol = 1e-6, local_maxeval = 10000)
    }

    # if still failed, relax further
    if (is.null(fit) || fit$status < 0) {
      fit <- .try_nloptr(xtol = 1e-5, maxeval = 50000,
                         constr_tol = 1e-3,
                         local_xtol = 1e-5, local_maxeval = 20000)
    }

    # if nloptr still fails, try penalty method as final fallback
    .try_penalty <- function(start_vec, rss_fun, lo, hi,
                             constraints, par_names) {
      penalty_weight <- 1e8
      .penalized_fun <- function(x) {
        obj <- rss_fun(x)
        pen <- 0
        if (!is.null(constraints$ineq)) {
          ineq_vals <- constraints$ineq(x)
          pen <- pen + sum(pmax(ineq_vals, 0)^2)
        }
        if (!is.null(constraints$eq)) {
          eq_vals <- constraints$eq(x)
          pen <- pen + sum(eq_vals^2)
        }
        obj + penalty_weight * pen
      }
      nloptr::nloptr(
        x0 = start_vec,
        eval_f = .penalized_fun,
        lb = lo,
        ub = hi,
        opts = list(
          algorithm = "NLOPT_LN_BOBYQA",
          xtol_rel = 1e-5,
          maxeval = 50000,
          print_level = 0
        )
      )
    }
    if (is.null(fit) || fit$status < 0) {
      fit <- .try_penalty(start_vec, rss_fun, lo, hi,
                          constraints, par_names)
    }

    if (!is.null(fit) && fit$status >= 0) {
      # build compatible result object
      par_est <- fit$solution
      names(par_est) <- par_names

      # compute fitted values and residuals
      ee <- list2env(as.list(par_est), parent = baseenv())
      ee[[lhs_name]] <- data[[lhs_name]]
      for (nm in names(data)) {
        if (!nm %in% par_names) ee[[nm]] <- data[[nm]]
      }
      y_pred <- eval(rhs, envir = ee)
      y_obs <- data[[lhs_name]]
      res <- y_obs - y_pred

      # approximate vcov (numeric Jacobian + NLS formula)
      .eval_r_at <- function(par_vec) {
        ee2 <- list2env(as.list(par_vec), parent = baseenv())
        ee2[[lhs_name]] <- data[[lhs_name]]
        for (nm in names(data)) {
          if (!nm %in% par_names) ee2[[nm]] <- data[[nm]]
        }
        y_obs - eval(rhs, envir = ee2)
      }
      h2 <- pmax(abs(par_est) * 1e-6, 1e-8)
      J <- matrix(0, nrow = nrow(data), ncol = length(par_est))
      for (i in seq_along(par_est)) {
        p_plus <- par_est; p_plus[i] <- par_est[i] + h2[i]
        p_minus <- par_est; p_minus[i] <- par_est[i] - h2[i]
        J[, i] <- (.eval_r_at(p_plus) - .eval_r_at(p_minus)) / (2 * h2[i])
      }
      s2 <- fit$objective / max(nrow(data) - length(par_est), 1)
      JJ <- tryCatch(solve(crossprod(J)), error = function(e) NULL)
      vcov_mat <- if (!is.null(JJ)) JJ * s2 else NULL

      result <- list(
        coefficients = par_est,
        fitted.values = y_pred,
        residuals = res,
        df.residual = nrow(data) - length(par_est) -
          if (!is.null(eval_g_eq)) length(eval_g_eq(par_est)) else 0,
        rss = fit$objective,
        convergence = fit$status,
        message = fit$message,
        formula = frm,
        data = data,
        vcov = vcov_mat,
        nloptr_result = fit,
        call = match.call()
      )
      class(result) <- "fit_constrained"
      return(result)
    }

    # nloptr failed
    stop("Constrained fit did not converge (nloptr status=", fit$status, "): ",
         fit$message, "\n  Suggestions: (1) Check constraints  (2) Adjust start values  (3) Relax bounds")
  }

# S3 methods: fit_constrained
#' @noRd
coef.fit_constrained <- function(object, ...) object$coefficients

#' @noRd
predict.fit_constrained <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  # reconstruct prediction from fit result
  if (is.null(object$formula) || is.null(object$data)) {
    warning("predict.fit_constrained: formula and data needed for newdata; returning fitted values")
    return(object$fitted.values)
  }
  par <- object$coefficients
  par_names <- names(par)
  frm <- object$formula
  rhs <- frm[[3]]

  ee <- list2env(as.list(par), parent = baseenv())
  for (nm in names(newdata)) {
    if (!nm %in% par_names) ee[[nm]] <- newdata[[nm]]
  }
  eval(rhs, envir = ee)
}

#' @noRd
residuals.fit_constrained <- function(object, ...) object$residuals

#' @noRd
print.fit_constrained <- function(x, ...) {
  rse <- sqrt(x$rss / max(x$df.residual, 1))
  tbl <- data.frame(term = names(x$coefficients), estimate = round(x$coefficients, 6),
                    row.names = NULL)
  title <- "Constrained fit results\n  Coefficients"
  .print_block(title, tbl)
  cat(sprintf("  Residual standard error: %g on %d degrees of freedom\n", rse, x$df.residual))
  cat(sprintf("  Convergence: %d%s\n", x$convergence,
              if (!is.null(x$message)) paste0(" (", x$message, ")") else ""))
  invisible(x)
}

# 7. Main function fit_equation

#' General-purpose equation fitting
#'
#' @param eq       equation identifier: ID("E07"), name("4PLC"),
#'                  or a custom formula string ("y ~ a*exp(-b*x)")
#' @param data     a data frame
#' @param x        x variable column name (default "x")
#' @param y        y variable column name (default "y")
#' @param start    named list of starting values (NULL=auto)
#' @param lower    named list/vector of lower bounds (e.g. list(Vmax=0, Km=0))
#' @param upper    named list/vector of upper bounds
#' @param constraints list of constraint functions: list(ineq = f, eq = g)
#'                    f(par) returns a vector; all values must be non-positive
#'                    g(par) returns a vector; must satisfy all(g(par) = 0)
#' @param weights  weights: NULL(equal) / numeric vector / method name string
#'                 methods: "equal", "1/y", "1/y^2", "1/x", "1/x^2"
#' @param ...      additional arguments passed to nlsLM() / nloptr()
#'
#' @importFrom minpack.lm nlsLM
#' @importFrom nloptr nloptr
#' @importFrom nls2 nls2
#'
#' @return an S3 object of class "fit_equation"
#' @examples
#' # ==== Response curve fitting ====
#'
#' # 1) Linear
#' df <- data.frame(x = 1:5, y = c(2.1, 4.0, 6.2, 7.9, 10.1))
#' f1 <- fit_equation("E01", df, "x", "y")
#' print(f1)
#'
#' # 2) Logarithmic
#' f2 <- fit_equation("E06", df, "x", "y")
#' print(f2)
#'
#' # 3) 4PLC (dose-response, descending)
#' df4plc <- data.frame(
#'   conc = c(0.1, 0.3, 1, 3, 10, 30, 100),
#'   resp = c(99, 96, 88, 72, 48, 25, 12)
#' )
#' f3 <- fit_equation("E07", df4plc, "conc", "resp")
#' print(f3)
#'
#' # 4) Custom formula
#' f4 <- fit_equation("a + b * log(x)", df4plc, "conc", "resp")
#' print(f4)
#'
#' # 5) Box constraints (parameter bounds)
#' f5 <- fit_equation("E07", df4plc, "conc", "resp",
#'   lower = list(D = 0), upper = list(D = 5))
#' print(f5)
#'
#' # 6) Inequality constraints (parameter relationship)
#' \donttest{
#' # require D > 1 (Hill slope cannot be too shallow)
#' f6 <- fit_equation("E07", df4plc, "conc", "resp",
#'   constraints = list(ineq = function(p) 1 - p["D"]))
#' print(f6)
#' }
#'
#' # 7) Equality constraints (fixed parameter relationship)
#' \donttest{
#' # fix b = 0.5 in exponential decay y = a*exp(-b*x)
#' df_exp <- data.frame(x = 0:5, y = c(10, 6.1, 3.7, 2.2, 1.4, 0.8))
#' f7 <- fit_equation("E04", df_exp, "x", "y",
#'   constraints = list(eq = function(p) p["b"] - 0.5))
#' print(f7)
#' }
#'
#' # 8) Weighted fitting via replicate_to_mean
#' \donttest{
#' # data with replicate measurements (8 dose levels, 3 replicates each)
#' df_rep <- data.frame(
#'   dose = rep(c(0.3, 0.6, 1.5, 4, 10, 25, 60, 150), each = 3),
#'   resp = c(4.1, 3.8, 3.9, 8.5, 8.9, 8.6,
#'            18.2, 17.9, 18.5,
#'            34.6, 35.1, 34.8,
#'            54.3, 53.9, 54.7,
#'            73.0, 73.6, 73.2,
#'            87.5, 87.1, 87.8,
#'            97.9, 98.3, 98.0)
#' )
#' # first summarize, generate inverse-variance weights
#' rep <- replicate_to_mean(df_rep, "dose", "resp", weights = "1/sd^2")
#'
#' # fit 4PLC with summarized means and explicit start values
#' f8 <- fit_equation("E07", rep, x = "x", y = "y_mean",
#'   weights = rep$weight,
#'   start = list(A = 2, B = 100, C = 8, D = -1.5))
#' print(f8)
#' }
#' @export
fit_equation <- function(eq, data, x = "x", y = "y",
                   start = NULL,
                   lower = NULL, upper = NULL,
                   constraints = NULL,
                   weights = NULL, ...) {
  cl <- match.call()
  reg <- eq_registry()

  # parse eq argument
  # handle formula object: convert to length-1 string
  if (inherits(eq, "formula")) {
    eq_str <- deparse(eq)
  } else {
    eq_str <- as.character(eq)
  }

  # case 2: ID or Name → look up registry (case-insensitive)
  row_id <- if (length(eq_str) == 1L) {
    which(
      tolower(reg$id) == tolower(eq_str) |
      tolower(reg$name) == tolower(eq_str)
    )
  } else integer(0)

  if (length(row_id) > 0) {
    eq_id   <- reg$id[row_id[1]]
    eq_name <- reg$name[row_id[1]]
    eq_cat  <- reg$category[row_id[1]]
    engine  <- reg$engine[row_id[1]]
    reg_row <- reg[row_id[1], ]
  } else {
    # case 3: custom formula — check for R formula operators
    .is_formula_like <- function(s) {
      grepl("[~+*/^()]", s)
    }
    if (length(eq_str) != 1L || !.is_formula_like(eq_str)) {
      stop("Equation '", paste(eq_str, collapse = " "), "' not found.",
           " Use list_equation() to see available equations,",
           " or supply a formula string like 'a*x + b'.")
    }
    eq_id   <- "CUSTOM"
    eq_name <- "Custom"
    eq_cat  <- "Custom"
    engine  <- "custom"
    reg_row <- NULL
  }

  # check data columns
  if (!x %in% names(data)) stop("Column '", x, "' not found in data.")
  if (!y %in% names(data)) stop("Column '", y, "' not found in data.")
  if (x == y) stop("x and y must refer to different columns, both are '", x, "'.")

  # check numeric (Response requires numeric x, y)
  if (eq_cat == "Response") {
    if (!is.numeric(data[[x]]))
      stop("x column '", x, "' must be numeric for ", eq_cat, " models.")
    if (!is.numeric(data[[y]]))
      stop("y column '", y, "' must be numeric for ", eq_cat, " models.")
  }

  # resolve weights (inline)
  .resolve_weights <- function(weights, data, x, y) {
    if (is.null(weights)) return(NULL)
    if (is.numeric(weights)) {
      if (length(weights) == 1) weights <- rep(weights, nrow(data))
      if (length(weights) != nrow(data))
        stop("weights length (", length(weights), ") does not match number of rows (", nrow(data), ")")
      if (any(weights < 0, na.rm = TRUE))
        stop("weights must be non-negative")
      if (anyNA(weights))
        stop("weights must not contain NA or NaN values")
      if (!all(is.finite(weights[!is.na(weights)])))
        stop("weights must be finite (no Inf or -Inf)")
      if (all(weights == 0))
        stop("all weights are zero; fitting is impossible")
      return(weights)
    }
    if (is.character(weights) && length(weights) == 1) {
      .xv <- function(nm) data[[nm]] %||% rep(1, nrow(data))
      `%||%` <- function(a, b) if (is.null(a)) b else a
      w <- switch(weights,
        "equal"    = rep(1, nrow(data)),
        "1/y"      = 1 / pmax(abs(.xv(y)), 1e-10),
        "1/y^2"    = 1 / pmax(.xv(y)^2, 1e-10),
        "1/x"      = 1 / pmax(abs(.xv(x)), 1e-10),
        "1/x^2"    = 1 / pmax(.xv(x)^2, 1e-10),
        "inverse"  = 1 / pmax(abs(.xv(y)), 1e-10),
        "inverse2" = 1 / pmax(.xv(y)^2, 1e-10),
        stop("Unknown weights method: '", weights, "'.",
             " Available: equal, 1/y, 1/y^2, 1/x, 1/x^2, inverse, inverse2")
      )
      ok_finite <- is.finite(w)
      if (any(ok_finite) && all(w[ok_finite] == 0))
        stop("The computed weights from method '", weights, "' are all zero; fitting is impossible.")
      if (any(w[ok_finite] < 0))
        stop("The computed weights from method '", weights, "' contain negative values.")
      return(w)
    }
    stop("weights must be NULL, a numeric vector, or a method name string.")
  }
  wt <- .resolve_weights(weights, data, x, y)

  # dispatch by category
  result <- switch(eq_cat,
    "Response" = fit_response(eq_id, reg_row, data, x, y,
                               start = start, engine = engine,
                               lower = lower, upper = upper,
                               constraints = constraints,
                               weights = wt, ...),
    "Custom"   = fit_custom(eq, data, x, y, start = start,
                             lower = lower, upper = upper,
                             constraints = constraints,
                             weights = wt, ...),
    stop("Unknown equation category: ", eq_cat)
  )

  # wrap S3 object
  formula_str <- if (!is.null(reg_row)) reg_row$formula_str[1L] else eq_str
  .make_fit_curve(result, cl, eq_id, eq_name, eq_cat, data, x, y,
                  formula_str = formula_str)
}

# Internal helper to construct S3 object
.make_fit_curve <- function(fit, call, eq_id, eq_name,
                          eq_cat, data, x, y,
                          formula_str = NULL) {
  # attempt to extract coefficients
  coefs <- tryCatch(stats::coef(fit), error = function(e) NULL)

  # attempt to extract AIC
  aic_val <- tryCatch({
    if (inherits(fit, "fit_constrained")) {
      n <- nrow(data)
      k <- length(coefs)
      n * log(fit$rss / n) + 2 * k
    } else {
      stats::AIC(fit)
    }
  }, error = function(e) NA)

  # attempt to extract R² (response curves only, not variance functions)
  r2 <- tryCatch({
    y_obs <- data[[y]]
    y_pred <- tryCatch({
      stats::predict(fit)
    }, error = function(e) NULL)
    if (!is.null(y_pred) && length(y_pred) == length(y_obs)) {
      # align predicted and observed values (handle NA)
      ok <- !is.na(y_obs) & !is.na(y_pred)
      if (sum(ok) < 2) {
        NA
      } else {
        tss <- sum((y_obs[ok] - mean(y_obs[ok]))^2)
        if (tss < .Machine$double.eps * sum(ok)) { NA_real_ }
        else { 1 - sum((y_obs[ok] - y_pred[ok])^2) / tss }
      }
    } else if (!is.null(y_pred) && length(y_pred) < length(y_obs)) {
      # with na.omit, predict.lm only returns non-NA rows
      ok <- !is.na(y_obs)
      if (sum(ok) == length(y_pred) && sum(ok) >= 2) {
        tss <- sum((y_obs[ok] - mean(y_obs[ok]))^2)
        if (tss < .Machine$double.eps * sum(ok)) { NA_real_ }
        else { 1 - sum((y_obs[ok] - y_pred)^2) / tss }
      } else {
        NA
      }
    } else {
      NA
    }
  }, error = function(e) NA)

  structure(list(
    call      = call,
    fit       = fit,
    eq_id     = eq_id,
    eq_name   = eq_name,
    category  = eq_cat,
    formula   = formula_str,
    data      = data,
    x_col     = x,
    y_col     = y,
    coefs     = coefs,
    aic       = aic_val,
    r_squared = r2
  ), class = "fit_equation")
}

# 8. S3 methods: coef / print / predict / residuals / plot

# 8a. print
#' Extract coefficients from an equation fit
#'
#' @param object A `fit_equation` object.
#' @param ... Additional arguments (currently unused).
#' @return A named numeric vector of model coefficients.
#' @export
#' @examples
#' df <- data.frame(x = 1:5, y = c(2.1, 4.0, 6.2, 7.9, 10.1))
#' f <- fit_equation("E01", df, "x", "y")
#' coef(f)
coef.fit_equation <- function(object, ...) {
  # Default: lm / nls / fit_constrained
  object$coefs
}

#' Print an equation fit
#'
#' @param x A `fit_equation` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
#' @export
print.fit_equation <- function(x, ...) {
  s <- x$fit

  # Section helpers

  .title <- function() {
    cat(sprintf("%s (%s)\n", x$eq_id, x$eq_name))
    if (!is.null(x$formula) && !is.na(x$formula) && nchar(x$formula) > 0)
      cat(sprintf("Formula: %s\n", x$formula))
    cat("\n")
  }

  .call <- function(call_text) {
    if (is.null(call_text)) return()
    cat("Call:\n  ", paste(deparse(call_text), collapse = "\n  "), "\n\n", sep = "")
  }

  .residuals <- function(resid_values) {
    if (is.null(resid_values) || length(resid_values) == 0) return()
    qq <- stats::quantile(resid_values, na.rm = TRUE)
    cat("Residuals:\n")
    cat(sprintf("    %5s  %7s  %7s  %7s  %7s\n",
        "Min", "1Q", "Median", "3Q", "Max"))
    cat(sprintf("    %5.3f  %7.3f  %7.3f  %7.3f  %7.3f\n\n",
        qq[1], qq[2], qq[3], qq[4], qq[5]))
  }

  .coefficients <- function(co, se, t_val, p_val, col1_width = 12) {
    cat("Coefficients:\n")
    cat(sprintf("  %-*s %10s %10s %8s %s\n", col1_width, "", "Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    stars <- rep("", length(p_val))
    stars[!is.na(p_val) & p_val < 0.1]   <- "."
    stars[!is.na(p_val) & p_val < 0.05]  <- "*"
    stars[!is.na(p_val) & p_val < 0.01]  <- "**"
    stars[!is.na(p_val) & p_val < 0.001] <- "***"
    for (i in seq_along(co)) {
      p_str <- if (is.na(p_val[i])) "      " else sprintf("%-7s", format.pval(p_val[i], digits = 3))
      cat(sprintf("  %-*s %10.6f %10.6f %8.3f %s\n",
          col1_width, names(co)[i], co[i], se[i] %||% NA_real_, t_val[i] %||% NA_real_,
          sprintf("%s %s", p_str, stars[i])))
    }
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  }

  .footer <- function(rse, df_r, extra) {
    parts <- character()
    if (is.finite(rse %||% NA_integer_) && is.finite(df_r %||% NA_integer_))
      parts <- c(parts, sprintf("Residual standard error: %g on %d degrees of freedom", rse, df_r))
    if (!is.null(extra)) parts <- c(parts, extra)
    cat(paste(parts, collapse = "\n"), "\n", sep = "")
  }

  # lm
  if (inherits(s, "lm")) {
    .title()
    .call(s$call)
    ssumm <- summary(s)
    .residuals(stats::residuals(s))
    # coefficients table with significance stars
    cat("Coefficients:\n")
    stats::printCoefmat(stats::coef(ssumm), digits = 6, signif.stars = TRUE, cs.ind = 1L)
    cat("\n")
    # footer
    extra <- sprintf("Multiple R-squared: %.4f,  Adjusted R-squared: %.4f",
                     ssumm$r.squared, ssumm$adj.r.squared)
    f <- ssumm$fstatistic
    if (!is.null(f) && is.vector(f) && length(f) >= 3) {
      p_f <- stats::pf(f[1L], f[2L], f[3L], lower.tail = FALSE)
      extra <- paste0(extra, sprintf("\nF-statistic: %g on %d and %d DF,  p-value: %s",
                      f[1L], f[2L], f[3L], format.pval(p_f, digits = 3)))
    }
    .footer(ssumm$sigma, ssumm$df[2L], extra)
    return(invisible(x))
  }

  # nls
  if (inherits(s, "nls")) {
    .title()
    .call(x$call)
    .residuals(stats::residuals(s))
    ssumm <- summary(s)
    cat("Coefficients:\n")
    stats::printCoefmat(ssumm$coefficients, digits = 6, signif.stars = TRUE, cs.ind = 1L)
    cat("\n")
    extra <- sprintf("AIC: %.2f", x$aic %||% NA_real_)
    n_iter <- tryCatch(s$convInfo$finIter, error = function(e) NULL)
    if (!is.null(n_iter)) extra <- paste0(extra, sprintf("  Iterations: %d", n_iter))
    .footer(ssumm$sigma, ssumm$df[2L], extra)
    return(invisible(x))
  }


  # fit_constrained
  if (inherits(s, "fit_constrained")) {
    .title()
    .call(x$call)
    rval <- s$residuals
    if (!is.null(rval) && length(rval) > 0) {
      .residuals(rval)
    }
    co <- s$coefficients
    if (!is.null(s$vcov)) {
      se <- sqrt(diag(s$vcov))
      t_val <- co / se
      p_val <- 2 * stats::pt(-abs(t_val), max(s$df.residual, 3))
    } else {
      se <- rep(NA_real_, length(co)); t_val <- rep(NA_real_, length(co)); p_val <- rep(NA_real_, length(co))
    }
    .coefficients(co, se, t_val, p_val)
    rse <- sqrt(s$rss / max(s$df.residual, 1))
    extra <- sprintf("AIC: %.2f", x$aic %||% NA_real_)
    if (is.finite(x$r_squared %||% NA_integer_))
      extra <- paste0(extra, sprintf(", R^2: %.4f", x$r_squared))
    extra <- paste0(extra, sprintf("\nConvergence: %d%s", s$convergence,
                    if (!is.null(s$message)) paste0(" (", s$message, ")") else ""))
    .footer(rse, s$df.residual, extra)
    return(invisible(x))
  }

  # generic fallback
  .title()

  if (!is.null(x$coefs)) {
    nm <- names(x$coefs)
    se <- rep(NA_real_, length(nm)); tv <- rep(NA_real_, length(nm)); pv <- rep(NA_real_, length(nm))
    .coefficients(x$coefs, se, tv, pv)
  }

  rse <- tryCatch({
    if (inherits(s, "lm") || inherits(s, "nls")) sigma(s) else NA_real_
  }, error = function(e) NA_real_)
  df_r <- tryCatch(stats::df.residual(s), error = function(e) NA_integer_)
  extra <- character()
  if (is.finite(x$aic %||% NA_integer_))
    extra <- c(extra, sprintf("AIC: %.2f", x$aic))
  if (is.finite(x$r_squared %||% NA_integer_))
    extra <- c(extra, sprintf("R^2: %.4f", x$r_squared))
  .footer(rse, df_r, paste(extra, collapse = ", "))

  invisible(x)
}

# 8b. predict
#' Predict fitted values or inverse-predict (from y to x)
#'
#' @param object   a fit_equation object
#' @param newdata  new data as a data.frame
#' @param x        x column name; NULL uses the original x column from the fit
#' @param y        y column name (inverse only); NULL uses the original y column
#' @param inverse  logical; if TRUE, solve for x given y values in newdata
#' @param interval "confidence" or "prediction"
#' @param level    confidence level (default 0.95)
#' @param ...      additional arguments passed to the underlying predict()
#' @return data.frame with x, y, and (if requested) confidence/prediction interval columns
#' @examples
#' df4plc <- data.frame(
#'   conc = c(0.1, 0.3, 1, 3, 10, 30, 100),
#'   resp = c(12, 25, 48, 72, 88, 96, 99)
#' )
#' f <- fit_equation("E07", df4plc, "conc", "resp")
#'
#' # Fitted values
#' predict(f)
#'
#' # Confidence interval
#' predict(f, interval = "confidence")
#'
#' # Prediction interval
#' predict(f, interval = "prediction")
#'
#' # Inverse prediction: estimate x for y=50
#' predict(f, newdata = data.frame(resp = 50), inverse = TRUE)
#' @export
predict.fit_equation <- function(object, newdata = NULL,
                              x = NULL, y = NULL,
                              inverse = FALSE,
                              interval = NULL, level = 0.95, ...) {
  fit <- object$fit
  x_col <- if (!is.null(x)) x else object$x_col
  y_col <- if (!is.null(y)) y else object$y_col

  # Build result data.frame
  .make_tbl <- function(xv, yv, ci = NULL, pi = NULL) {
    tbl <- data.frame(x = xv, y = yv, row.names = NULL)
    if (!is.null(ci)) { tbl$y_ci_lwr <- ci$y_lwr; tbl$y_ci_upr <- ci$y_upr }
    if (!is.null(pi)) { tbl$y_pi_lwr <- pi$y_lwr; tbl$y_pi_upr <- pi$y_upr }
    tbl
  }

  # helper: print data.frame and return invisibly
  .show_pred <- function(tbl) {
    print(tbl, row.names = FALSE)
    invisible(tbl)
  }

  # inline inverse prediction function
  .inverse_predict <- function(object, y_target) {
    eq_id <- object$eq_id
    co <- object$coefs
    if (is.null(co)) stop("No coefficient estimates available.")
    x_vals <- object$data[[object$x_col]]
    x_range <- range(x_vals, na.rm = TRUE)
    x_ext <- diff(x_range) * 3
    x_est <- tryCatch({
      # auto-detect coefficient names: prefer mapped names, fallback to lm defaults
      .c <- function(nm, def = NULL) {
        v <- co[nm]
        if (!is.na(v)) v else co[def]
      }
      switch(eq_id,
      "E01" = { (y_target - .c("a", "(Intercept)")) / .c("b", "x") },
      "E02" = {
        a <- .c("a", "(Intercept)"); b <- .c("b", "x"); c <- .c("c", "I(x^2)")
        D <- b^2 - 4*c*(a - y_target)
        if (D < 0) stop("No real solution.")
        x1 <- (-b + sqrt(D)) / (2*c)
        x2 <- (-b - sqrt(D)) / (2*c)
        if (abs(x1 - mean(x_vals)) < abs(x2 - mean(x_vals))) x1 else x2
      },
      "E03" = { if (y_target / co["a"] <= 0) stop("y/a must be > 0.")
                log(y_target / co["a"]) / co["b"] },
      "E04" = { if (y_target / co["a"] <= 0) stop("y/a must be > 0.")
                -log(y_target / co["a"]) / co["b"] },
      "E05" = { if (y_target / co["a"] <= 0) stop("y/a must be > 0.")
                (y_target / co["a"])^(1 / co["b"]) },
      "E06" = { exp((y_target - .c("a", "(Intercept)")) / .c("b", "log(x)")) },
      "E07" = {
        A <- co["A"]; B <- co["B"]; C <- co["C"]; D <- co["D"]
        if (y_target <= min(A, B) || y_target >= max(A, B))
          stop("y must be between A (", round(A, 3), ") and B (", round(B, 3), ").")
        ratio <- (B - A) / (y_target - A) - 1
        if (ratio <= 0) stop("Computed ratio <= 0; y may be outside valid range.")
        C * ratio^(1 / D)
      },
      "E08" = {
        A <- co["A"]; B <- co["B"]; C <- co["C"]; D <- co["D"]; E <- co["E"]
        if (y_target <= min(A, B) || y_target >= max(A, B))
          stop("y must be between A (", round(A, 3), ") and B (", round(B, 3), ").")
        ratio <- (B - A) / (y_target - A) - 1
        if (ratio <= 0) stop("Computed ratio <= 0; y may be outside valid range.")
        C * ratio^(1 / (D * E))
      },
      "E09" = { if (y_target / co["a"] <= 0) return(NA_real_)
                t <- -2 * co["c"]^2 * log(y_target / co["a"])
                if (t < 0) return(NA_real_)
                co["b"] - sqrt(t) },
      "E10" = { Vmax <- co["Vmax"]; Km <- co["Km"]
                if (y_target >= Vmax) stop("y must be < Vmax = ", round(Vmax, 3))
                Km * y_target / (Vmax - y_target) },
      "E11" = { a <- co["a"]; b <- co["b"]; c <- co["c"]; d <- co["d"]
                (asin((y_target - d) / a) - c) / b },
      "E12" = {
        a <- .c("a", "(Intercept)"); b <- .c("b", "x"); c <- .c("c", "I(x^2)"); d <- .c("d", "I(x^3)")
        polyroot(c(a - y_target, b, c, d)) |> Re() |> (\(x) x[which.min(abs(x - mean(x_vals)))])()
      },
      NULL
    )}, error = function(e) NULL)
    if (!is.null(x_est) && length(x_est) == 1L && is.finite(x_est)) return(x_est)
    .forward <- function(x) {
      nd <- data.frame(x = x); names(nd) <- object$x_col
      utils::capture.output(pred <- predict(object, newdata = nd))
      as.numeric(pred$y[1L])
    }
    lo <- min(x_range[1] - x_ext, if (x_range[1] > 0) 0 else x_range[1] - x_ext)
    hi <- x_range[2] + x_ext
    for (mult in c(1, 2, 5, 10, 50)) {
      r_lo <- lo * mult; r_hi <- hi * mult
      n_grid <- 200; grid_x <- seq(r_lo, r_hi, length.out = n_grid)
      grid_f <- vapply(grid_x, function(x) suppressWarnings(.forward(x) - y_target), numeric(1))
      finite_ok <- is.finite(grid_f)
      if (sum(finite_ok) < 2) next
      signs <- sign(grid_f[finite_ok])
      idx <- which(diff(signs) != 0)
      solutions <- c()
      for (i in idx) {
        x1 <- grid_x[finite_ok][i]; x2 <- grid_x[finite_ok][i + 1]
        sol <- tryCatch(stats::uniroot(function(x) suppressWarnings(.forward(x) - y_target),
                                       c(x1, x2), maxiter = 1000),
                        error = function(e) NULL)
        if (!is.null(sol)) solutions <- c(solutions, sol$root)
      }
      if (length(solutions) >= 1) {
        if (length(solutions) == 1) return(solutions[1])
        center <- mean(x_range)
        return(solutions[which.min(abs(solutions - center))])
      }
    }
    stop("Inverse prediction failed: no solution found for y = ", y_target,
         " in the search range.")
  }

  # inverse prediction (y → x)
  if (isTRUE(inverse)) {
    if (is.null(newdata)) stop("'newdata' is required for inverse prediction.")
    y_vals <- if (is.data.frame(newdata)) {
      if (y_col %in% names(newdata)) newdata[[y_col]] else newdata[[1L]]
    } else {
      newdata
    }
    if (is.null(y_vals)) stop("No y values found in newdata.")

    x_vals <- vapply(y_vals, function(y_target) {
      .inverse_predict(object, y_target)
    }, numeric(1))
    tbl <- .make_tbl(xv = round(x_vals, 6), yv = y_vals)
    return(.show_pred(tbl))
  }

  # get x values and predictions
  nd <- if (is.null(newdata)) object$data else newdata
  x_vals <- if (is.data.frame(nd)) nd[[x_col]] else seq_len(nrow(nd))

  # nls / lm: rename columns and predict
  nd_pred <- nd
  if (is.data.frame(nd_pred) && x_col %in% names(nd_pred)) {
    names(nd_pred)[names(nd_pred) == x_col] <- "x"
  }
  if (is.data.frame(nd_pred) && y_col %in% names(nd_pred)) {
    names(nd_pred)[names(nd_pred) == y_col] <- "y"
  }
  y_pred <- as.numeric(stats::predict(fit, newdata = nd_pred, ...))

  # confidence / prediction intervals
  ci <- NULL
  pi <- NULL
  if (!is.null(interval)) {
    interval <- match.arg(interval, c("confidence", "prediction"), several.ok = TRUE)
    if ("confidence" %in% interval) {
      ci <- .predict_interval(object, nd, level = level, interval = "confidence")
    }
    if ("prediction" %in% interval) {
      pi <- .predict_interval(object, nd, level = level, interval = "prediction")
    }
  }

  tbl <- .make_tbl(xv = x_vals, yv = y_pred, ci = ci, pi = pi)
  .show_pred(tbl)
}


#' Compute fitted confidence/prediction intervals
#'
#' @param object  a fit_equation object
#' @param newdata new data as a data.frame
#' @param level   confidence level (default 0.95)
#' @param interval "confidence" or "prediction"
#' @noRd
.predict_interval <- function(object, newdata, level = 0.95,
                         interval = c("confidence", "prediction")) {
  interval <- match.arg(interval)
  fit <- object$fit
  is_pred <- interval == "prediction"

  # rename columns
  nd <- newdata
  if (is.data.frame(nd) && object$x_col %in% names(nd)) {
    names(nd)[names(nd) == object$x_col] <- "x"
  }

  alpha <- 1 - level

  # lm — native se.fit + interval support
  if (inherits(fit, "lm")) {
    pred <- stats::predict(fit, newdata = nd, se.fit = TRUE,
                           interval = interval, level = level)
    return(data.frame(
      y     = as.numeric(pred$fit[, "fit"]),
      y_lwr = as.numeric(pred$fit[, "lwr"]),
      y_upr = as.numeric(pred$fit[, "upr"])
    ))
  }


  # generic: numeric Jacobian + Delta method (nls, fit_constrained)
  # extract parameters, formula, vcov
  if (inherits(fit, "fit_constrained")) {
    V <- fit$vcov
    if (is.null(V)) return(NULL)
    coefs <- fit$coefficients
    frm <- fit$formula
    df_resid <- fit$df.residual
  } else if (inherits(fit, "nls")) {
    V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    if (is.null(V)) return(NULL)
    coefs <- stats::coef(fit)
    frm <- stats::formula(fit)
    df_resid <- stats::df.residual(fit)
  } else {
    return(NULL)
  }

  par_names <- names(coefs)
  n_par <- length(coefs)
  n_obs <- nrow(nd)
  rhs <- frm[[3]]

  .eval_f <- function(par_vec) {
    names(par_vec) <- par_names
    ee <- list2env(as.list(par_vec), parent = baseenv())
    for (nm in names(nd)) {
      if (is.numeric(nd[[nm]])) ee[[nm]] <- nd[[nm]]
    }
    tryCatch(eval(rhs, envir = ee), error = function(e) rep(NA_real_, n_obs))
  }

  h <- pmax(abs(coefs) * 1e-6, 1e-8)
  grad <- matrix(0, nrow = n_obs, ncol = n_par)
  for (i in seq_len(n_par)) {
    p_plus <- coefs; p_plus[i] <- coefs[i] + h[i]
    p_minus <- coefs; p_minus[i] <- coefs[i] - h[i]
    grad[, i] <- (.eval_f(p_plus) - .eval_f(p_minus)) / (2 * h[i])
  }

  y_pred <- .eval_f(coefs)

  # confidence interval variance: Var(ŷ) = G'·V·G
  se <- numeric(n_obs)
  for (j in seq_len(n_obs)) {
    g <- grad[j, , drop = FALSE]
    se[j] <- sqrt(max(g %*% V %*% t(g), 0))
  }

  # prediction interval: Var(ŷ - y_new) = G'·V·G + σ²
  if (is_pred) {
    rss <- if (inherits(fit, "fit_constrained")) fit$rss else
             sum(stats::residuals(fit)^2, na.rm = TRUE)
    df_r <- max(df_resid, 1)
    sigma2 <- rss / df_r
    se <- sqrt(se^2 + sigma2)
  }

  df <- max(df_resid, 3)
  t_val <- stats::qt(1 - alpha / 2, df)

  data.frame(
    y     = as.numeric(y_pred),
    y_lwr = as.numeric(y_pred) - t_val * se,
    y_upr = as.numeric(y_pred) + t_val * se
  )
}

# 8c. residuals
#' Extract residuals from an equation fit
#'
#' @param object A `fit_equation` object.
#' @param ... Additional arguments passed to the fitted model.
#' @return A numeric vector of residuals.
#' @examples
#' df <- data.frame(x = 1:5, y = c(2.1, 4.0, 6.2, 7.9, 10.1))
#' f <- fit_equation("E01", df, "x", "y")
#' residuals(f)
#' @export
residuals.fit_equation <- function(object, ...) {
  fit <- object$fit
  raw <- if (inherits(fit, "fit_constrained")) {
    fit$residuals
  } else {
    stats::residuals(fit, ...)
  }
  # Strip attributes to return a plain vector
  if (!is.numeric(raw)) raw <- as.numeric(raw)
  attributes(raw) <- NULL
  raw
}

# 8d. plot
#' Plot an equation fit
#'
#' @param x A `fit_equation` object.
#' @param interval Interval type, either `"confidence"` or `"prediction"`.
#' @param level Confidence level.
#' @param frame Whether to draw a frame.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
#' @import ggplot2
#' @export
#' @examples
#' df4plc <- data.frame(
#'   conc = c(0.1, 0.3, 1, 3, 10, 30, 100),
#'   resp = c(12, 25, 48, 72, 88, 96, 99)
#' )
#' f <- fit_equation("E07", df4plc, "conc", "resp")
#' plot(f)
#' plot(f, interval = "prediction", frame = TRUE)
plot.fit_equation <- function(x, interval = "confidence", level = 0.95,
                            frame = FALSE, ...) {
  fit  <- x$fit
  data <- x$data


  # extract data
  y_obs <- data[[x$y_col]]
  x_val <- data[[x$x_col]]
  y_pred <- tryCatch(unname(predict(fit)), error = function(e) NULL)

  # use ggplot2 for plotting
  .plot_ggplot(x, x_val, y_obs, y_pred,
               interval = interval, level = level, frame = frame)
  invisible(x)
}

#' ggplot2 version of plot (with confidence/prediction intervals)
#' @noRd
.plot_ggplot <- function(fit_curve, x_val, y_obs, y_pred,
                         interval = "confidence", level = 0.95, frame = FALSE) {
  # inline: extract formula string for annotation
  .get_formula_text <- function(fit_curve) {
    eq_id <- fit_curve$eq_id
    if (eq_id != "CUSTOM") {
      reg <- eq_registry()
      row <- reg[tolower(reg$id) == tolower(eq_id), ]
      if (nrow(row) > 0) return(paste0("y = ", row$formula_str[1L]))
    }
    fit <- fit_curve$fit
    if (inherits(fit, "lm") || inherits(fit, "nls")) {
      return(deparse(stats::formula(fit)))
    }
    NULL
  }

  df <- data.frame(x = x_val, y = y_obs)
  if (!is.null(y_pred)) df$y_pred <- y_pred

  # ordered x for curve drawing
  x_sort <- seq(min(x_val, na.rm = TRUE),
                max(x_val, na.rm = TRUE), length.out = 200)
  nd <- data.frame(x = x_sort)
  df_line <- data.frame(x = x_sort)
  fit <- fit_curve$fit
  df_line$y <- tryCatch(
    stats::predict(fit, newdata = nd),
    error = function(e) NULL
  )
  nd_pred <- nd
  names(nd_pred)[1] <- fit_curve$x_col  # for .predict_interval (uses object$x_col for rename)

  # confidence/prediction interval
  df_ci <- NULL
  if (is.character(interval) && length(interval) == 1L) {
    df_ci <- tryCatch(
      .predict_interval(fit_curve, nd_pred, level = level, interval = interval),
      error = function(e) NULL
    )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(size = 2.5, color = "steelblue", alpha = 0.7, na.rm = TRUE) +
    ggplot2::labs(title = paste0(fit_curve$eq_id, " -- ", fit_curve$eq_name),
                  subtitle = paste0("AIC = ", round(fit_curve$aic, 2),
                                    "  R^2 = ", round(fit_curve$r_squared, 4)),
                  x = fit_curve$x_col, y = fit_curve$y_col) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # ribbon
  if (!is.null(df_ci)) {
    p <- p + ggplot2::geom_ribbon(
      data = data.frame(x = x_sort, df_ci),
      ggplot2::aes(x = .data$x, ymin = .data$y_lwr, ymax = .data$y_upr),
      fill = "darkorange", alpha = 0.15, na.rm = TRUE
    )
  }

  if (!is.null(df_line$y)) {
    p <- p + ggplot2::geom_line(data = df_line,
                                 ggplot2::aes(x = .data$x, y = .data$y),
                                 color = "darkorange", linewidth = 1, na.rm = TRUE)
  }

  # formula annotation
  if (isTRUE(frame)) {
    eq_str <- .get_formula_text(fit_curve)
    if (!is.null(eq_str)) {
      p <- p + ggplot2::geom_text(
        data = data.frame(x = -Inf, y = Inf, label = eq_str),
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        hjust = -0.05, vjust = 1.5, size = 3.5, inherit.aes = FALSE
      )
    }
  }

  print(p)
}

# 9. Multi-equation comparison

#' Fit and compare multiple equations
#'
#' @param data a data frame
#' @param x    x column name
#' @param y    y column name
#' @param eqs  equation name/ID vector; NULL = all Response equations
#' @param ...  additional arguments passed to fit_equation
#' @return An object of class "compare_equation" with components:
#'   \item{table}{data.frame sorted by AIC: Eq, Name, Formula, nPar, AIC, deltaAIC, R2, RSE}
#'   \item{fits}{list of fit_equation objects, named by eq_id}
#'   \item{failed}{character vector of eq IDs that failed to fit}
#'   \item{data, x, y}{input data reference}
#' @examples
#' df <- data.frame(
#'   x = c(0.1, 0.3, 1, 3, 10, 30, 100),
#'   y = c(12, 25, 48, 72, 88, 96, 99)
#' )
#' compare_equation(df, "x", "y", eqs = c("E01", "E06", "E07"))
#' @export
compare_equation <- function(data, x = "x", y = "y",
                        eqs = NULL, ...) {
  cl <- match.call()
  reg <- eq_registry()
  response_reg <- reg[reg$category == "Response", ]

  if (is.null(eqs)) {
    eqs <- response_reg$name
    eq_ids <- response_reg$id
  } else {
    # parse eqs (case-insensitive)
    eq_ids <- sapply(eqs, function(e, reg) {
      r <- reg[tolower(reg$id) == tolower(e) |
               tolower(reg$name) == tolower(e), ]
      if (nrow(r) == 0) return(NA)
      r$id[1]
    }, reg = reg, USE.NAMES = FALSE)
    eqs <- sapply(eqs, function(e, reg) {
      r <- reg[tolower(reg$id) == tolower(e) |
               tolower(reg$name) == tolower(e), ]
      if (nrow(r) == 0) return(NA)
      r$name[1]
    }, reg = reg, USE.NAMES = FALSE)
    ok <- !is.na(eq_ids)
    eqs <- eqs[ok]
    eq_ids <- eq_ids[ok]
  }

  results <- list()
  failed <- character()
  n_total <- length(eq_ids)
  for (i in seq_along(eq_ids)) {
    f <- tryCatch(
      fit_equation(eq_ids[i], data, x = x, y = y, ...),
      error = function(e) NULL
    )
    if (!is.null(f)) {
      results[[eq_ids[i]]] <- f
    } else {
      failed <- c(failed, eq_ids[i])
    }
  }

  if (length(results) == 0) {
    cat("All equations failed to fit.\n")
    return(invisible(NULL))
  }

  # extract formula strings from fit objects
  formula_str <- sapply(results, function(r) r$formula %||% "")

  # residual standard error
  rse <- sapply(results, function(r) {
    tryCatch({
      fit <- r$fit
      if (inherits(fit, "lm") || inherits(fit, "nls")) {
        sigma(fit)
      } else if (inherits(fit, "fit_constrained")) {
        sqrt(fit$rss / max(fit$df.residual, 1))
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
  })

  aic_vals <- sapply(results, `[[`, "aic")
  r2_vals  <- sapply(results, `[[`, "r_squared")
  n_par    <- sapply(results, function(r) length(r$coefs))

  # summary table — truncate formula to keep table compact
  .trunc <- function(s, w = 24) {
    ifelse(nchar(s) > w, paste0(substr(s, 1, w - 3), "..."), s)
  }
  summ <- data.frame(
    Eq     = sapply(results, `[[`, "eq_id"),
    Name   = sapply(results, `[[`, "eq_name"),
    nPar   = n_par,
    AIC    = round(aic_vals, 2),
    deltaAIC = 0,
    R2     = round(r2_vals, 4),
    RSE    = round(rse, 6),
    Formula = .trunc(unname(formula_str)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  # sort by AIC ascending, then fill deltaAIC
  ord <- order(aic_vals, na.last = TRUE)
  summ <- summ[ord, ]
  best_aic <- min(aic_vals, na.rm = TRUE)
  summ$deltaAIC <- round(summ$AIC - best_aic, 2)

  result <- list(
    table  = summ,
    fits   = results[ord],
    failed = failed,
    data   = data,
    x      = x,
    y      = y,
    call   = cl,
    n_total = n_total
  )
  class(result) <- "compare_equation"
  result
}

#' Print a multiple-equation comparison
#'
#' @param x A `compare_equation` object.
#' @param ... Reserved for the S3 generic.
#' @return The unchanged `x`, invisibly.
#' @export
print.compare_equation <- function(x, ...) {
  n_ok   <- nrow(x$table)
  n_fail <- length(x$failed)
  best   <- x$table[1, ]

  # build header
  txt <- "Equation Comparison\n"
  txt <- paste0(txt, sprintf("  Best model           : %s (%s)\n", best$Name, best$Eq))
  txt <- paste0(txt, sprintf("  Best AIC             : %.2f\n", best$AIC))
  if (n_ok > 1)
    txt <- paste0(txt, sprintf("  AIC range            : %.2f\n", diff(range(x$table$AIC, na.rm = TRUE))))
  txt <- paste0(txt, sprintf("  Models fitted        : %d / %d\n", n_ok, x$n_total))
  if (n_fail > 0)
    txt <- paste0(txt, sprintf("  Failed models        : %s\n", paste(x$failed, collapse = ", ")))
  txt <- paste0(txt, "\n")

  # build comparison table as formatted text
  tbl <- x$table
  eq_col   <- format(tbl$Eq, justify = "left")
  name_col <- format(tbl$Name, justify = "left")
  npar_col <- sprintf("%d", tbl$nPar)
  aic_col  <- sprintf("%.2f", tbl$AIC)
  daic_col <- sprintf("%.2f", tbl$deltaAIC)
  r2_col   <- sprintf("%.4f", tbl$R2)
  rse_col  <- sprintf("%.6f", tbl$RSE)
  frm_col  <- format(tbl$Formula, justify = "left")

  hdr <- sprintf("    %s %s %s %s %s %s %s %s %s",
    "*", format("Eq", width = max(nchar(eq_col))),
    format("Name", width = max(nchar(name_col))),
    format("nPar", width = max(nchar(npar_col)), justify = "right"),
    format("AIC", width = max(nchar(aic_col)), justify = "right"),
    format("deltaAIC", width = max(nchar(daic_col)), justify = "right"),
    format("R2", width = max(nchar(r2_col)), justify = "right"),
    format("RSE", width = max(nchar(rse_col)), justify = "right"),
    format("Formula", width = max(nchar(frm_col)), justify = "left"))

  sep <- paste(rep("-", nchar(hdr)), collapse = "")

  txt <- paste0(txt, hdr, "\n")
  txt <- paste0(txt, sep, "\n")
  for (i in seq_len(n_ok)) {
    flag <- if (i == 1) "*" else " "
    txt <- paste0(txt, sprintf("    %s %s %s %s %s %s %s %s %s",
      flag,
      format(eq_col[i], width = max(nchar(eq_col))),
      format(name_col[i], width = max(nchar(name_col))),
      format(npar_col[i], width = max(nchar(npar_col)), justify = "right"),
      format(aic_col[i], width = max(nchar(aic_col)), justify = "right"),
      format(daic_col[i], width = max(nchar(daic_col)), justify = "right"),
      format(r2_col[i], width = max(nchar(r2_col)), justify = "right"),
      format(rse_col[i], width = max(nchar(rse_col)), justify = "right"),
      format(frm_col[i], width = max(nchar(frm_col)), justify = "left")), "\n")
  }

  txt <- paste0(txt, "\n  * = best model by AIC\n")
  cat(txt)
  invisible(x)
}
