#' ROC Curve Analysis —— Receiver Operating Characteristic
#'
#' Functions: construct ROC analysis object, AUC calculation, optimal cutoff analysis,
#' multivariate logistic regression, ROC curve visualization.
#'
#' @importFrom stats sd median quantile qnorm glm predict
#'   coef as.formula setNames
#' @importFrom utils capture.output
#' @importFrom ggplot2 ggplot aes geom_abline geom_line geom_step geom_vline
#'   geom_point geom_ribbon geom_text labs coord_equal theme_bw theme
#'   element_text margin
#' @noRd
NULL

# Compute diagnostic metrics for a single cutoff value
#
# @param values   numeric vector of evaluation column values
# @param truth    gold standard 0/1 vector
# @param cutoff   cutoff value
# @param direction "geq": values >= cutoff predicted positive; "leq": values <= cutoff
# @return list, containing cutoff, tp, fp, tn, fn, sensitivity, specificity, youden,
#         precision, npv, accuracy
.roc_curve <- function(values, truth, direction = c("geq", "leq")) {
  direction <- match.arg(direction)

  # Include both sentinel cutoffs in the same ordering as finite candidates.
  # In ROC coordinates (FPR, sensitivity), this traverses from (0, 0) to
  # (1, 1): Inf -> finite -> -Inf for "geq", and the reverse for "leq".
  cutoff_candidates <- sort(unique(c(-Inf, values, Inf)),
                            decreasing = (direction == "geq"))

  # Compute diagnostic metrics for each cutoff (inlined from .roc_at_cutoff)
  results <- lapply(cutoff_candidates, function(c) {
    pred <- if (direction == "geq") as.integer(values >= c) else as.integer(values <= c)
    tp <- sum(pred == 1 & truth == 1, na.rm = TRUE)
    fp <- sum(pred == 1 & truth == 0, na.rm = TRUE)
    tn <- sum(pred == 0 & truth == 0, na.rm = TRUE)
    fn <- sum(pred == 0 & truth == 1, na.rm = TRUE)
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    list(
      cutoff = c, tp = tp, fp = fp, tn = tn, fn = fn,
      sensitivity = sensitivity, specificity = specificity,
      youden = sensitivity + specificity - 1,
      precision = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
      npv = if ((tn + fn) > 0) tn / (tn + fn) else NA_real_,
      accuracy = (tp + tn) / (tp + fp + tn + fn)
    )
  })

  do.call(rbind, lapply(results, function(r) {
    data.frame(
      cutoff      = r$cutoff,
      sensitivity = r$sensitivity,
      specificity = r$specificity,
      tp          = r$tp,
      fp          = r$fp,
      tn          = r$tn,
      fn          = r$fn,
      youden      = r$youden,
      precision   = r$precision,
      npv         = r$npv,
      accuracy    = r$accuracy
    )
  }))
}

# Select finite Youden/corner optima from a cutoff table.
# Sentinel cutoffs are retained for complete ROC curves but are not useful
# reported decision thresholds.
.roc_optimal_indices <- function(curve_df) {
  finite <- which(is.finite(curve_df$cutoff))
  if (!length(finite)) return(list(youden = NA_integer_, corner = NA_integer_))
  finite_youden <- finite[is.finite(curve_df$youden[finite])]
  youden <- if (length(finite_youden))
    finite_youden[which.max(curve_df$youden[finite_youden])] else NA_integer_
  distance <- sqrt((1 - curve_df$sensitivity)^2 +
                     (1 - curve_df$specificity)^2)
  finite_corner <- finite[is.finite(distance[finite])]
  corner <- if (length(finite_corner))
    finite_corner[which.min(distance[finite_corner])] else NA_integer_
  list(youden = youden, corner = corner)
}

# Return one row per distinct selected cutoff and combine method labels when
# Youden and corner select the same threshold.
.roc_cutpoint_rows <- function(cutoff_table, column,
                               cutpoint = c("Youden", "Corner", "all")) {
  cutpoint <- match.arg(cutpoint)
  ct <- cutoff_table[cutoff_table$column == column, , drop = FALSE]
  if (!nrow(ct)) return(data.frame())
  idx <- .roc_optimal_indices(ct)
  methods <- switch(cutpoint, Youden = "Youden", Corner = "Corner",
                    all = c("Youden", "Corner"))
  rows <- lapply(methods, function(method) {
    i <- if (method == "Youden") idx$youden else idx$corner
    if (is.na(i)) return(NULL)
    out <- ct[i, , drop = FALSE]
    out$method <- method
    out
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  if (anyDuplicated(out$cutoff)) {
    split_rows <- split(seq_len(nrow(out)), out$cutoff)
    out <- do.call(rbind, lapply(split_rows, function(i) {
      z <- out[i[1L], , drop = FALSE]
      z$method <- paste(unique(out$method[i]), collapse = " / ")
      z
    }))
  }
  rownames(out) <- NULL
  out
}

# Trapezoidal AUC
#
# @param sensitivity  sensitivity vector
# @param specificity  specificity vector
# @return AUC value (in [0, 1] range)
.auc <- function(sensitivity, specificity) {
  fpr <- 1 - specificity
  sens <- sensitivity
  # At identical false-positive rates, traverse the vertical ROC segment from
  # lower to higher sensitivity. Ordering only by FPR can leave tied points in
  # cutoff order and use the wrong endpoint for the next horizontal segment,
  # underestimating AUC (including for a perfectly separating marker).
  ord <- order(fpr, sens)
  fpr <- fpr[ord]
  sens <- sens[ord]

  ok <- is.finite(fpr) & is.finite(sens)
  fpr <- fpr[ok]
  sens <- sens[ok]

  if (length(fpr) < 2L) return(NA_real_)

  # Trapezoidal integration
  sum(diff(fpr) * (sens[-1L] + sens[-length(sens)]) / 2)
}

# Compute AUC standard error (Hanley & McNeil normal approximation)
#
# Reference: Hanley JA, McNeil BJ. Radiology 1982; 143(1):29-36
#
# @param auc       AUC value
# @param n_pos     number of positive samples
# @param n_neg     number of negative samples
# @return AUC standard error
.auc_se <- function(auc, n_pos, n_neg) {
  if (n_pos < 1L || n_neg < 1L) return(NA_real_)
  if (is.na(auc) || auc <= 0 || auc >= 1) return(NA_real_)

  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)

  se <- sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) /
               (n_pos * n_neg))
  se
}

# EP24-A2 Table 5 / Hanley-McNeil (1983) conversion from the average
# within-group rating correlation and average AUC to the correlation between
# two paired AUC estimates. Rows cover 0.02 to 0.90; columns cover 0.700 to
# 0.975. Linear interpolation avoids discontinuities from manual rounding.
.ep24_auc_correlation_table <- matrix(c(
  .02,.02,.02,.02,.02,.02,.02,.02,.01,.01,.01,.01,
  .04,.04,.03,.03,.03,.03,.03,.03,.03,.02,.02,.02,
  .05,.05,.05,.05,.05,.05,.05,.04,.04,.04,.03,.02,
  .07,.07,.07,.07,.07,.06,.06,.06,.06,.05,.04,.03,
  .09,.09,.09,.09,.08,.08,.08,.07,.07,.06,.06,.04,
  .11,.11,.11,.10,.10,.10,.09,.09,.08,.08,.07,.05,
  .13,.12,.12,.12,.12,.11,.11,.11,.10,.09,.08,.06,
  .14,.14,.14,.14,.13,.13,.13,.12,.11,.11,.09,.07,
  .16,.16,.16,.16,.15,.15,.14,.14,.13,.12,.11,.09,
  .18,.18,.18,.17,.17,.17,.16,.15,.15,.14,.12,.10,
  .20,.20,.19,.19,.19,.18,.18,.17,.16,.15,.14,.11,
  .22,.22,.21,.21,.21,.20,.19,.19,.18,.17,.15,.12,
  .24,.23,.23,.23,.22,.22,.21,.20,.19,.18,.16,.13,
  .26,.25,.25,.25,.24,.24,.23,.22,.21,.20,.18,.15,
  .27,.27,.27,.26,.26,.25,.25,.24,.23,.21,.19,.16,
  .29,.29,.29,.28,.28,.27,.26,.26,.24,.23,.21,.18,
  .31,.31,.31,.30,.30,.29,.28,.27,.26,.25,.23,.19,
  .33,.33,.32,.32,.31,.31,.30,.29,.28,.26,.24,.21,
  .35,.35,.34,.34,.33,.33,.32,.31,.30,.28,.26,.22,
  .37,.37,.36,.36,.35,.35,.34,.33,.32,.30,.28,.24,
  .39,.39,.38,.38,.37,.36,.36,.35,.33,.32,.29,.25,
  .41,.40,.40,.40,.39,.38,.38,.37,.35,.34,.31,.27,
  .43,.42,.42,.42,.41,.40,.39,.38,.37,.35,.33,.29,
  .45,.44,.44,.43,.43,.42,.41,.40,.39,.37,.35,.30,
  .47,.46,.46,.45,.45,.44,.43,.42,.41,.39,.37,.32,
  .49,.48,.48,.47,.47,.46,.45,.44,.43,.41,.39,.34,
  .51,.50,.50,.49,.49,.48,.47,.46,.45,.43,.41,.36,
  .53,.52,.52,.51,.51,.50,.49,.48,.47,.45,.43,.38,
  .55,.54,.54,.53,.53,.52,.51,.50,.49,.47,.45,.40,
  .57,.56,.56,.55,.55,.54,.53,.52,.51,.49,.47,.42,
  .59,.58,.58,.57,.57,.56,.55,.54,.53,.51,.49,.45,
  .61,.60,.60,.59,.59,.58,.58,.57,.55,.54,.51,.47,
  .63,.62,.62,.62,.61,.60,.60,.59,.57,.56,.53,.49,
  .65,.64,.64,.64,.63,.62,.62,.61,.60,.58,.56,.51,
  .67,.66,.66,.66,.65,.65,.64,.63,.62,.60,.58,.54,
  .69,.69,.68,.68,.67,.67,.66,.65,.64,.63,.60,.56,
  .71,.71,.70,.70,.69,.69,.68,.67,.66,.65,.63,.59,
  .73,.73,.72,.72,.72,.71,.71,.70,.69,.67,.65,.61,
  .75,.75,.75,.74,.74,.73,.73,.72,.71,.70,.68,.64,
  .77,.77,.77,.76,.76,.76,.75,.74,.73,.72,.70,.67,
  .79,.79,.79,.79,.78,.78,.77,.77,.76,.75,.73,.70,
  .82,.81,.81,.81,.81,.80,.80,.79,.78,.77,.76,.73,
  .84,.84,.83,.83,.83,.82,.82,.81,.81,.80,.78,.75,
  .86,.86,.86,.85,.85,.85,.84,.84,.83,.82,.81,.79,
  .88,.88,.88,.88,.87,.87,.87,.86,.86,.85,.84,.82
), nrow = 45L, byrow = TRUE,
dimnames = list(sprintf("%.2f", seq(.02, .90, .02)),
                sprintf("%.3f", seq(.700, .975, .025))))

.ep24_auc_correlation <- function(rating_correlation, average_auc) {
  row_grid <- seq(.02, .90, .02)
  col_grid <- seq(.700, .975, .025)
  used_rating <- min(max(rating_correlation, min(row_grid)), max(row_grid))
  used_auc <- min(max(average_auc, min(col_grid)), max(col_grid))
  by_row <- apply(.ep24_auc_correlation_table, 1L, function(values) {
    stats::approx(col_grid, values, xout = used_auc, rule = 2)$y
  })
  value <- stats::approx(row_grid, by_row, xout = used_rating, rule = 2)$y
  list(
    value = unname(value),
    rating_correlation = used_rating,
    average_auc = used_auc,
    clamped = !isTRUE(all.equal(used_rating, rating_correlation)) ||
      !isTRUE(all.equal(used_auc, average_auc))
  )
}

.roc_curve_source <- function(x, name) {
  marker_match <- name %in% x$cols
  model_match <- name %in% names(x$mlr_results)
  if (!marker_match && !model_match)
    stop("ROC curve '", name, "' was not found.", call. = FALSE)
  if (marker_match && model_match)
    stop("ROC curve name '", name, "' is ambiguous between a marker and an ",
         "MLR model. Rename the MLR model.", call. = FALSE)
  if (marker_match) "marker" else "model"
}

.roc_comparison_entry <- function(x, name) {
  source <- .roc_curve_source(x, name)
  if (source == "marker") {
    entry <- x$auc_list[[name]]
    if (is.null(entry))
      stop("ROC marker curve '", name, "' has no AUC result. Run auc() for ",
           "this marker before comparison.", call. = FALSE)
    direction <- entry$direction %||% "geq"
    score <- x$cols_values[[name]]
    if (direction == "leq") score <- -score
    return(list(name = name, source = "marker", score = score,
                auc = entry$auc, se = entry$auc_se,
                direction = direction))
  }
  entry <- x$mlr_results[[name]]
  list(name = name, source = "model", score = entry$predicted,
       auc = entry$auc, se = entry$auc_se, direction = "geq")
}

.delong_components <- function(scores, truth) {
  positive <- truth == 1L
  m <- sum(positive); n <- sum(!positive)
  v10 <- vapply(scores, function(score) {
    comparison <- outer(score[positive], score[!positive], function(a, b) {
      (a > b) + 0.5 * (a == b)
    })
    rowMeans(comparison)
  }, numeric(m))
  v01 <- vapply(scores, function(score) {
    comparison <- outer(score[positive], score[!positive], function(a, b) {
      (a > b) + 0.5 * (a == b)
    })
    colMeans(comparison)
  }, numeric(n))
  if (is.null(dim(v10))) v10 <- matrix(v10, ncol = length(scores))
  if (is.null(dim(v01))) v01 <- matrix(v01, ncol = length(scores))
  list(
    auc = colMeans(v10),
    covariance = stats::cov(v10) / m + stats::cov(v01) / n,
    v10 = v10,
    v01 = v01
  )
}

# Convert reference column to binary 0/1
#
# @param ref       original reference values
# @param positive  specify positive level (character or numeric)
# @return list, containing binary (0/1 vector) and message (diagnostic message)
.ref_to_binary <- function(ref, positive = NULL) {
  msg <- NULL

  # Case 1: already numeric 0/1
  if (is.numeric(ref) && all(na.omit(ref) %in% c(0, 1))) {
    if (is.null(positive)) {
      bin <- as.integer(ref)
    } else {
      pos_val <- if (is.character(positive)) as.integer(positive) else positive
      if (!pos_val %in% c(0, 1))
        stop("'positive' must be 0 or 1 when reference is already 0/1.", call. = FALSE)
      bin <- as.integer(ref == pos_val)
    }
    return(list(binary = bin, message = "binary (0/1)"))
  }

  # Case 2: factor
  if (is.factor(ref)) {
    lvls <- levels(ref)
    if (length(lvls) != 2L)
      stop("Factor reference must have exactly 2 levels, got ", length(lvls), ".",
           call. = FALSE)
    if (is.null(positive)) {
      # Default: last level as positive
      bin <- as.integer(ref == lvls[2L])
      msg <- sprintf("'%s' treated as positive (default; use 'positive' to change)", lvls[2L])
    } else {
      if (!positive %in% lvls)
        stop("'positive = ", positive, "' is not a level of 'reference'. Levels: ",
             paste(lvls, collapse = ", "), call. = FALSE)
      bin <- as.integer(ref == positive)
      msg <- sprintf("'%s' treated as positive", positive)
    }
    return(list(binary = bin, message = msg))
  }

  # Case 3: character
  if (is.character(ref)) {
    uniq <- unique(na.omit(ref))
    if (length(uniq) != 2L)
      stop("Character reference must have exactly 2 unique values, got ", length(uniq), ".",
           call. = FALSE)
    if (is.null(positive)) {
      # Last level after sorting is positive
      lvls <- sort(uniq)
      bin <- as.integer(ref == lvls[2L])
      msg <- sprintf("'%s' treated as positive (default; use 'positive' to change)", lvls[2L])
    } else {
      if (!positive %in% uniq)
        stop("'positive = ", positive, "' not found in reference values. Values: ",
             paste(uniq, collapse = ", "), call. = FALSE)
      bin <- as.integer(ref == positive)
      msg <- sprintf("'%s' treated as positive", positive)
    }
    return(list(binary = bin, message = msg))
  }

  # Case 4: quantitative numeric (must be pre-converted with continuous_to_binary)
  if (is.numeric(ref))
    stop("'reference' column appears to be quantitative (non-binary numeric).\n",
         "  Use continuous_to_binary() first to convert it to a 0/1 binary column,\n",
         "  then pass the new binary column as 'reference'.", call. = FALSE)

  stop("Unsupported 'reference' column type.", call. = FALSE)
}

#' Convert quantitative columns to binary based on cutoff values
#'
#' Creates new binary columns (0/1) from existing quantitative columns using
#' specified cutoff values, appending them to the original data frame.
#' New columns are named \code{<col>_binary}.
#'
#' @param data   data frame
#' @param cols   character vector of column names to binarize
#' @param cutoff numeric vector of cutoff values; recycled to \code{length(cols)}
#'               if a single value is given
#'
#' @return data frame with original columns plus new \code{<col>_binary} columns
#'
#' @examples
#'   df <- data.frame(x = 1:10, y = rnorm(10))
#'   continuous_to_binary(df, "x", 5)
#'   continuous_to_binary(df, c("x", "y"), c(5, 0))
#' @export
continuous_to_binary <- function(data, cols, cutoff) {
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(cols) || length(cols) < 1L)
    stop("'cols' must be a character vector of column names.", call. = FALSE)
  bad <- setdiff(cols, names(data))
  if (length(bad) > 0)
    stop("Column(s) not found: ", paste(bad, collapse = ", "), call. = FALSE)
  if (!is.numeric(cutoff))
    stop("'cutoff' must be numeric.", call. = FALSE)
  if (length(cutoff) == 1L)
    cutoff <- rep(cutoff, length(cols))
  if (length(cutoff) != length(cols))
    stop("Length of 'cutoff' must be 1 or equal to length of 'cols'.", call. = FALSE)
  if (anyNA(cutoff) || any(!is.finite(cutoff)))
    stop("'cutoff' must contain only finite values.", call. = FALSE)

  for (i in seq_along(cols)) {
    col <- cols[i]
    if (!is.numeric(data[[col]]))
      stop("Column '", col, "' must be numeric.", call. = FALSE)
    new_col <- paste0(col, "_binary")
    data[[new_col]] <- as.integer(data[[col]] >= cutoff[i])
  }

  data
}

#' roc constructor
#'
#' Create an ROC analysis object.
#'
#' Supports multiple evaluation columns (quantitative) and a single reference column
#' (binary 0/1, factor, or character with two levels). Automatically detects direction,
#' reports missing values and ID duplicates.
#'
#' @param data       data frame
#' @param cols       evaluation column names (character vector, quantitative, can be multiple)
#' @param reference  reference column name (character, single column). Must be binary 0/1,
#'                   a factor with 2 levels, or a character vector with 2 unique values.
#'                   If the original reference is quantitative, use \code{continuous_to_binary()}
#'                   first to convert it to 0/1 and pass the new column name.
#' @param id         ID column name (optional, for duplicate detection)
#' @param positive   specify the positive level of reference (for factor/character)
#'
#' @return Returns an S3 object of class "roc"
#'
#' @examples
#'   df <- data.frame(
#'     sid = 1:100,
#'     x1  = c(rnorm(50, 10, 2), rnorm(50, 12, 2)),
#'     ref = rep(c(0, 1), each = 50),
#'     age = 1:100,
#'     sex = rep(c("F", "M"), each = 50)
#'   )
#'   obj <- roc(df, cols = "x1", reference = "ref", id = "sid")
#'   print(obj)
#' @export
roc <- function(data, cols, reference, id = NULL, positive = NULL) {

  # Input validation
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(cols) || length(cols) < 1L)
    stop("'cols' must be a character vector of column names.", call. = FALSE)
  bad_cols <- setdiff(cols, names(data))
  if (length(bad_cols) > 0)
    stop("Column(s) in 'cols' not found: ", paste(bad_cols, collapse = ", "), call. = FALSE)
  if (!is.character(reference) || length(reference) != 1L || !reference %in% names(data))
    stop("'reference' must be a single column name in 'data'.", call. = FALSE)
  if (!is.null(id) && (!is.character(id) || length(id) != 1L || !id %in% names(data)))
    stop("'id' must be NULL or a single column name in 'data'.", call. = FALSE)
  if (reference %in% cols)
    stop("'reference' column cannot also be in 'cols'.", call. = FALSE)

  # Extract data
  id_vec   <- if (!is.null(id)) data[[id]] else seq_len(nrow(data))
  ref_orig <- data[[reference]]

  # Verify all cols are quantitative
  for (col in cols) {
    if (!is.numeric(data[[col]]))
      stop("Column '", col, "' in 'cols' must be numeric.", call. = FALSE)
  }

  # Convert reference to binary
  ref_result <- .ref_to_binary(ref_orig, positive)
  ref_binary <- ref_result$binary
  ref_msg <- ref_result$message

  # Missing value detection
  n_total <- nrow(data)
  miss_info <- list()
  all_ok <- rep(TRUE, n_total)

  for (col in cols) {
    miss <- is.na(data[[col]]) | !is.finite(data[[col]])
    miss_info[[col]] <- list(n_miss = sum(miss), miss_rows = which(miss))
    all_ok <- all_ok & !miss
  }

  miss_ref <- is.na(ref_binary) | !is.finite(ref_binary)
  miss_info[[reference]] <- list(n_miss = sum(miss_ref), miss_rows = which(miss_ref))
  all_ok <- all_ok & !miss_ref

  n_complete <- sum(all_ok)
  n_miss_total <- n_total - n_complete

  if (n_complete < 3L)
    stop("Need at least 3 complete cases for ROC analysis, got ", n_complete, ".", call. = FALSE)

  # ID duplicate detection
  n_dup <- 0L
  dup_ids <- integer(0)
  if (!is.null(id)) {
    dups <- duplicated(id_vec) | duplicated(id_vec, fromLast = TRUE)
    n_dup <- sum(dups)
    if (n_dup > 0)
      dup_ids <- which(dups)
  }

  # Extract complete case data
  id_vec_cc    <- id_vec[all_ok]
  ref_binary_cc <- ref_binary[all_ok]

  # Count positive and negative samples
  n_pos <- sum(ref_binary_cc == 1)
  n_neg <- sum(ref_binary_cc == 0)
  if (n_pos < 1L || n_neg < 1L)
    stop("Reference must contain at least one positive and one negative complete case.",
         call. = FALSE)

  cols_values <- list()
  col_direction <- list()

  # Compute direction for each col (direction only, no AUC)
  for (col in cols) {
    v <- data[[col]][all_ok]
    cols_values[[col]] <- v

    # Auto-detect direction
    mean_pos <- mean(v[ref_binary_cc == 1], na.rm = TRUE)
    mean_neg <- mean(v[ref_binary_cc == 0], na.rm = TRUE)
    dir <- "geq"
    if (is.finite(mean_pos) && is.finite(mean_neg) && mean_pos < mean_neg) {
      dir <- "leq"
    }
    col_direction[[col]] <- dir
  }

  # print slot: data overview
  ref_type <- if (is.factor(ref_orig) || is.character(ref_orig) ||
                  (is.numeric(ref_orig) && all(na.omit(ref_orig) %in% c(0,1)))) {
    if (is.numeric(ref_orig) && all(na.omit(ref_orig) %in% c(0,1))) "binary"
    else if (is.factor(ref_orig)) "factor"
    else "character"
  } else {
    "quantitative"
  }

  miss_ids_cols <- list()
  for (col in cols) {
    mi <- miss_info[[col]]
    miss_ids_cols[[col]] <- mi$miss_rows
  }

  print_slot <- list(
    data_class  = class(data)[1L],
    n           = n_complete,
    n_total     = n_total,
    n_miss      = n_miss_total,
    n_pos       = n_pos,
    n_neg       = n_neg,
    n_dup       = n_dup,
    var_names   = list(cols = cols, reference = reference),
    ref_type    = ref_type,
    reference_message = ref_msg,
    dup_ids     = dup_ids,
    miss_ids    = list(cols = miss_ids_cols, reference = miss_info[[reference]]$miss_rows)
  )

  # Construct object
  out <- list(
    call               = match.call(),
    data               = data,
    id                 = id_vec_cc,
    id_name            = id %||% "row_number",
    cols               = cols,
    cols_values        = cols_values,
    col_direction      = col_direction,
    reference_binary   = ref_binary_cc,
    reference_name     = reference,
    reference_original = ref_orig,
    n                  = n_complete,
    n_pos              = n_pos,
    n_neg              = n_neg,
    n_total            = n_total,
    print              = print_slot,
    describe           = NULL,
    auc_list           = list(),
    auc_table          = NULL,
    auc_comparisons    = list(),
    cutoff_table       = NULL,
    mlr_results        = list()
  )
  class(out) <- "roc"

  invisible(out)
}

#' Print a ROC analysis object
#'
#' @param x A `roc` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
#' @export
print.roc <- function(x, ...) {
  if (!inherits(x, "roc"))
    stop("Input must be a 'roc' object.", call. = FALSE)

  p <- x$print %||% list(
    data_class = class(x$data)[1L], n = x$n, n_total = x$n_total,
    n_miss = x$n_total - x$n, n_pos = x$n_pos, n_neg = x$n_neg,
    n_dup = 0L, var_names = list(cols = x$cols, reference = x$reference_name),
    ref_type = "binary", reference_message = NULL, dup_ids = integer(0),
    miss_ids = list(cols = list(), reference = integer(0))
  )

  cat("\nROC Analysis\n")
  cat(sprintf("  Data class:        %s\n", p$data_class))
  cat(sprintf("  Total rows:        %d\n", p$n_total))
  cat(sprintf("  Complete cases:    %d\n", p$n))
  if (p$n_miss > 0)
    cat(sprintf("  Missing/excluded:  %d\n", p$n_miss))
  cat(sprintf("  Reference:         %s (%s)\n", p$var_names$reference, p$ref_type))
  if (!is.null(p$reference_message))
    cat(sprintf("  Reference coding:  %s\n", p$reference_message))
  cat(sprintf("  Positive / Neg:    %d / %d\n", p$n_pos, p$n_neg))
  cat(sprintf("  Evaluation cols:   %d  (%s)\n",
              length(p$var_names$cols), paste(p$var_names$cols, collapse = ", ")))
  if (p$n_dup > 0)
    cat(sprintf("  Duplicate IDs:     %d (rows: %s)\n",
                p$n_dup, paste(p$dup_ids, collapse = ", ")))
  else
    cat("  Duplicate IDs:     none\n")

  # Missing detail
  if (p$n_miss > 0) {
    cat("  Missing value detail:\n")
    for (col in names(p$miss_ids$cols)) {
      mi_rows <- p$miss_ids$cols[[col]]
      if (length(mi_rows) > 0)
        cat(sprintf("    %s: %d row(s): %s\n", col, length(mi_rows),
                    paste(mi_rows, collapse = ", ")))
    }
    if (length(p$miss_ids$reference) > 0)
      cat(sprintf("    %s: %d row(s): %s\n", p$var_names$reference,
                  length(p$miss_ids$reference),
                  paste(p$miss_ids$reference, collapse = ", ")))
  }

  # Analysis slot status
  cat("\nAnalysis slots:\n")
  slots_desc <- list(
    describe  = "Describe",
    auc       = "AUC",
    auc_compare = "Paired AUC comparison",
    cutoff    = "Cutoff",
    mlr       = "MLR"
  )
  for (nm in names(slots_desc)) {
    if (nm == "auc") {
      val <- length(x$auc_list) > 0
    } else if (nm == "auc_compare") {
      val <- length(x$auc_comparisons) > 0
    } else if (nm == "mlr") {
      val <- length(x$mlr_results) > 0
    } else {
      val <- !is.null(x[[nm]])
    }
    cat(sprintf("    %s %s\n", if (val) "[x]" else "[ ]", slots_desc[[nm]]))
  }

  cat("\n")
  invisible(x)
}


#' Print method for roc AUC table
#' @param x roc_auc_table object (data.frame)
#' @param ... additional arguments
#' @return The unchanged \code{roc_auc_table} data frame, invisibly.
#' @export
print.roc_auc_table <- function(x, ...) {
  cat("\nAUC Table\n")
  cat(sprintf("  %-15s %8s  %14s  %8s\n", "Column", "AUC", "95% CI", "Direction"))
  cat(" ", paste(rep("-", 60), collapse = ""), "\n")
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  %-15s %8s  [%s, %s]  %8s\n",
                x$Column[i], x$AUC[i], x$CI_lower[i], x$CI_upper[i], x$Direction[i]))
  }
  cat("\n")
  invisible(x)
}

#' Print a paired AUC comparison
#'
#' @param x A `roc_auc_comparison` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object, invisibly.
#' @export
print.roc_auc_comparison <- function(x, ...) {
  cat("\nPaired AUC Comparison\n")
  cat(sprintf("  Curves: %s - %s\n", x$curve1, x$curve2))
  cat(sprintf("  Method: %s\n", x$method_label))
  cat(sprintf("  Subjects: %d (positive = %d, negative = %d)\n",
              x$n, x$n_pos, x$n_neg))
  cat(sprintf("  AUC %s: %.6f\n", x$curve1, x$auc1))
  cat(sprintf("  AUC %s: %.6f\n", x$curve2, x$auc2))
  cat(sprintf("  Difference: %.6f\n", x$difference))
  cat(sprintf("  SE difference: %.6f\n", x$se))
  cat(sprintf("  %.1f%% CI: [%.6f, %.6f]\n",
              100 * x$conf.level, x$conf.int[1L], x$conf.int[2L]))
  cat(sprintf("  z = %.6f\n", x$statistic))
  cat(sprintf("  p-value: %s\n", format.pval(x$p.value, digits = 6)))
  if (identical(x$method, "ep24")) {
    cat(sprintf("  Mean within-group rating correlation: %.6f\n",
                x$rating_correlation[["average"]]))
    cat(sprintf("  EP24 area correlation: %.6f\n",
                x$area_correlation))
    if (isTRUE(x$lookup$clamped))
      cat(sprintf("  Table lookup used bounded coordinates: correlation %.4f, AUC %.4f\n",
                  x$lookup$rating_correlation, x$lookup$average_auc))
  }
  if (isTRUE(x$apparent_model_comparison))
    cat("  Note: fitted training probabilities are treated as fixed scores;\n",
        "        model optimism and fitting uncertainty are not corrected.\n",
        sep = "")
  cat("\n")
  invisible(x)
}

#' Print method for roc cutoff table
#' @param x roc_cutoff_table object (data.frame)
#' @param ... additional arguments
#' @return The unchanged \code{roc_cutoff_table} data frame, invisibly.
#' @export
print.roc_cutoff_table <- function(x, ...) {
  all_tables <- attr(x, "all_tables")
  cols <- if (!is.null(all_tables)) names(all_tables) else unique(x$column)

  cat("\nOptimal Cutoff Analysis\n")

  for (col in cols) {
    ct <- if (!is.null(all_tables)) all_tables[[col]] else x[x$column == col, ]
    optimal_rows <- ct[ct$optimal != "", , drop = FALSE]
    cat(sprintf("  Column: %s\n", col))
    for (i in seq_len(nrow(optimal_rows))) {
      or <- optimal_rows[i, ]
      cat(sprintf("    %s: cutoff = %s  (sens = %.4f, spec = %.4f, Youden = %.4f)\n",
                  or$optimal, .format_num(or$cutoff),
                  or$sensitivity, or$specificity, or$youden))
    }

    cat("\n  Full cutoff table:\n")
    cat(sprintf("  %12s  %8s  %8s  %5s %5s %5s %5s  %8s\n",
                "Cutoff", "Sens", "Spec", "TP", "FP", "TN", "FN", "Youden"))
    cat(" ", paste(rep("-", 70), collapse = ""), "\n")
    for (j in seq_len(nrow(ct))) {
      cat(sprintf("  %12s  %8s  %8s  %5d %5d %5d %5d  %8s\n",
                  if (is.infinite(ct$cutoff[j])) as.character(ct$cutoff[j])
                  else .format_num(ct$cutoff[j], 4L),
                  .format_num(ct$sensitivity[j], 4L),
                  .format_num(ct$specificity[j], 4L),
                  ct$tp[j], ct$fp[j], ct$tn[j], ct$fn[j],
                  .format_num(ct$youden[j], 4L)))
    }
  }

  cat("\n")
  invisible(x)
}

#' Print a multivariate ROC model
#'
#' @param x A `roc_mlr` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
#' @export
print.roc_mlr <- function(x, ...) {
  cat("\nMultivariate Logistic Regression (", x$name, ")\n", sep = "")
  cat(" ", paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("  Formula: %s\n", deparse(x$formula)))
  cat(sprintf("  n = %d  (positive = %d, negative = %d)\n",
              length(x$predicted), sum(x$fit$y == 1), sum(x$fit$y == 0)))
  ci_l <- x$ci_lower %||% NA_real_
  ci_u <- x$ci_upper %||% NA_real_
  if (!is.na(ci_l))
    cat(sprintf("  AUC = %.4f  [%.4f, %.4f]\n", x$auc, ci_l, ci_u))
  else
    cat(sprintf("  AUC = %.4f\n", x$auc))
  cat("\n")

  s <- x$fit_summary
  cat("  Coefficients:\n")
  coef_tab <- coef(s)
  txt <- capture.output(print(coef_tab))
  for (line in txt) cat("   ", line, "\n")

  cat("\n  * Note: AUC is evaluated on the same data used for fitting;\n")
  cat("    the estimate may be optimistic.\n\n")
  invisible(x)
}

#' Descriptive Statistics (S3 method)
#'
#' Print distribution summaries for evaluation columns and the reference column (if numeric).
#' When the reference is binary, evaluation column distributions are shown grouped by positive/negative.
#'
#' @param x       roc object
#' @param cols    columns to describe; defaults to all evaluation columns + reference column (if numeric)
#' @param digits  number of decimal places, default 4
#' @param ...     reserved arguments
#'
#' @return Updated roc object with results stored in $describe
#'
#' @examples
#'   obj <- roc(ivd_roc_example, cols = "x1", reference = "ref")
#'   describe(obj)
#' @export
describe.roc <- function(x, cols = NULL, digits = 4L, ...) {

  if (!inherits(x, "roc"))
    stop("Input must be a 'roc' object.")

  data <- x$data
  all_names <- names(data)

  # Determine which columns to analyze
  if (is.null(cols)) {
    exclude <- if (!is.null(x$id_name) && x$id_name != "row_number") x$id_name else NULL
    cols <- setdiff(all_names, exclude)
  } else {
    bad <- setdiff(cols, all_names)
    if (length(bad) > 0)
      stop("Column(s) not found: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  n_total <- nrow(data)
  n_col <- length(cols)

  # Per-column summary
  columns <- list()

  for (col in cols) {
    vec     <- data[[col]]
    n_na    <- sum(is.na(vec))
    pct_na  <- round(n_na / n_total * 100, 1)

    # Force reference column as categorical
    if (col == x$reference_name) {
      tbl   <- table(vec, useNA = "ifany")
      n_ok  <- sum(!is.na(vec))
      freqs <- as.data.frame(tbl, responseName = "Freq")
      names(freqs)[1] <- "Level"
      freqs$Pct <- round(freqs$Freq / n_ok * 100, 1)
      na_row <- which(is.na(freqs$Level))
      if (length(na_row))
        freqs$Pct[na_row] <- round(freqs$Freq[na_row] / n_total * 100, 1)

      columns[[col]] <- list(
        type     = "categorical",
        n_na     = n_na,
        pct_na   = pct_na,
        n_levels = nrow(freqs),
        table    = freqs
      )

    } else if (is.factor(vec) || is.character(vec)) {
      tbl   <- table(vec, useNA = "ifany")
      n_ok  <- sum(!is.na(vec))
      freqs <- as.data.frame(tbl, responseName = "Freq")
      names(freqs)[1] <- "Level"
      freqs$Pct <- round(freqs$Freq / n_ok * 100, 1)
      na_row <- which(is.na(freqs$Level))
      if (length(na_row))
        freqs$Pct[na_row] <- round(freqs$Freq[na_row] / n_total * 100, 1)

      columns[[col]] <- list(
        type     = "categorical",
        n_na     = n_na,
        pct_na   = pct_na,
        n_levels = nrow(freqs),
        table    = freqs
      )

    } else if (is.numeric(vec)) {
      ok  <- vec[!is.na(vec)]
      if (length(ok) == 0) {
        columns[[col]] <- list(
          type = "numeric", n_na = n_na, pct_na = pct_na,
          n = 0L, mean = NA_real_, sd = NA_real_,
          min = NA_real_, q25 = NA_real_, median = NA_real_,
          q75 = NA_real_, max = NA_real_
        )
      } else {
        columns[[col]] <- list(
          type   = "numeric",
          n_na   = n_na,
          pct_na = pct_na,
          n      = length(ok),
          mean   = mean(ok),
          sd     = sd(ok),
          min    = min(ok),
          q25    = quantile(ok, 0.25, names = FALSE),
          median = median(ok),
          q75    = quantile(ok, 0.75, names = FALSE),
          max    = max(ok)
        )
      }
    } else {
      columns[[col]] <- list(type = "other", n_na = n_na, pct_na = pct_na)
    }
  }

  result <- structure(list(
    n_complete     = x$n,
    n_total        = n_total,
    n_col          = n_col,
    columns        = columns,
    data           = data,
    digits         = digits
  ), class = "roc_describe")
  x$describe <- result
  print(result)
  invisible(x)
}

#' Print method for roc describe results
#' @param x roc_describe object
#' @param digits number of decimal places, default 4
#' @param ... additional arguments
#' @return The unchanged \code{roc_describe} object, invisibly.
#' @export
print.roc_describe <- function(x, digits = NULL, ...) {
  d <- x
  columns <- d$columns
  data    <- d$data
  if (is.null(digits)) digits <- d$digits %||% 4L

  cat("\nDescriptive Statistics\n")
  sep75 <- paste(rep("-", 75), collapse = "")
  cat(sprintf("  %-18s %5s %8s %8s %8s %8s %8s %8s %8s  %s\n",
              "Variable", "n", "Mean", "SD", "Min", "Q1", "Median", "Q3", "Max", "NA"))
  cat(" ", sep75, "\n")

  for (col in names(columns)) {
    ci <- columns[[col]]
    if (!is.null(ci) && ci$type == "numeric") {
      vec <- data[[col]]
      n_na <- ci$n_na
      na_str <- if (n_na > 0) sprintf("%d", n_na) else ""
      if (length(vec) > 0) {
        ok <- vec[!is.na(vec)]
        cat(sprintf("  %-18s %5d %8s %8s %8s %8s %8s %8s %8s  %s\n",
                    col, length(ok),
                    .format_num(mean(ok), digits),
                    .format_num(sd(ok), digits),
                    .format_num(min(ok), digits),
                    .format_num(quantile(ok, 0.25, names = FALSE), digits),
                    .format_num(median(ok), digits),
                    .format_num(quantile(ok, 0.75, names = FALSE), digits),
                    .format_num(max(ok), digits),
                    na_str))
      }
    }
  }

  # Categorical columns
  for (col in names(columns)) {
    ci <- columns[[col]]
    if (!is.null(ci) && ci$type == "categorical") {
      na_str <- if (ci$n_na > 0)
        sprintf(" (%d NA, %.1f%%)", ci$n_na, ci$pct_na) else ""
      cat(sprintf("\n  %s%s\n", col, na_str))
      cat(sprintf("  %-20s %5s  %7s\n", "Level", "Count", "Percent"))
      cat("  ", paste(rep("-", 38), collapse = ""), "\n")
      for (r in seq_len(nrow(ci$table))) {
        lv <- as.character(ci$table$Level[r])
        if (is.na(lv)) lv <- "<NA>"
        cat(sprintf("  %-20s %5d  %6.1f%%\n",
                    lv, ci$table$Freq[r], ci$table$Pct[r]))
      }
    }
  }

  cat("\n")
  invisible(x)
}

#'
#' Compute ROC curves and AUC for each evaluation column.
#'
#' @param x    roc object
#' @param cols column names to analyze; default is all
#' @param ...  reserved arguments
#'
#' @return Updated roc object with results stored in $auc_list and $auc_table
#'
#' @examples
#'   obj <- roc(ivd_roc_example, cols = c("x1", "x2"), reference = "ref")
#'   obj <- auc(obj)
#'   obj <- auc(obj, cols = "x1")
#' @export
auc.roc <- function(x, cols = NULL, ...) {

  if (!inherits(x, "roc"))
    stop("Input must be a 'roc' object.")

  if (is.null(cols)) {
    cols <- x$cols
  } else {
    bad <- setdiff(cols, x$cols)
    if (length(bad) > 0)
      stop("Column(s) not in roc object: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  n_pos <- x$n_pos
  n_neg <- x$n_neg
  auc_list <- x$auc_list  # preserve existing results

  for (col in cols) {
    v <- x$cols_values[[col]]
    dir <- x$col_direction[[col]] %||% "geq"

    # Compute ROC curve and AUC
    curve <- .roc_curve(v, x$reference_binary, direction = dir)
    auc_val <- .auc(curve$sensitivity, curve$specificity)
    auc_se  <- .auc_se(auc_val, n_pos, n_neg)

    # 95% CI
    ci_lower <- if (!is.na(auc_se)) max(0, auc_val - 1.96 * auc_se) else NA_real_
    ci_upper <- if (!is.na(auc_se)) min(1, auc_val + 1.96 * auc_se) else NA_real_

    auc_list[[col]] <- list(
      col       = col,
      auc       = auc_val,
      auc_se    = auc_se,
      ci_lower  = ci_lower,
      ci_upper  = ci_upper,
      direction = dir
    )
  }

  # Build AUC table
  auc_df <- do.call(rbind, lapply(auc_list[intersect(names(auc_list), x$cols)], function(a) {
    data.frame(
      Column    = a$col,
      AUC       = .format_num(a$auc, 4L),
      CI_lower  = .format_num(a$ci_lower, 4L),
      CI_upper  = .format_num(a$ci_upper, 4L),
      Direction = a$direction,
      stringsAsFactors = FALSE
    )
  }))
  rownames(auc_df) <- NULL

  x$auc_list <- auc_list
  x$auc_table <- structure(auc_df, class = c("roc_auc_table", "data.frame"))

  # Print
  print(x$auc_table)

  # Note on AUC < 0.5
  for (col in cols) {
    a <- auc_list[[col]]
    if (!is.na(a$auc) && a$auc < 0.5) {
      cat(sprintf("  Note: AUC for '%s' = %s (< 0.5). Consider checking direction.\n",
                  col, .format_num(a$auc)))
    }
  }

  cat("\n")
  invisible(x)
}

#' Compare two paired ROC AUCs
#'
#' Compare two ROC curves evaluated on the common complete cases stored in a
#' single [roc()] object. A curve can be a marker previously processed by
#' [auc()] or a named multivariate model produced by [mlr()].
#'
#' `method = "ep24"` implements the paired Hanley--McNeil approximation in
#' CLSI EP24-A2. It combines the two single-curve Hanley--McNeil standard
#' errors with the correlation between the paired AUCs. The latter is obtained
#' by linear interpolation of EP24-A2 Table 5 from the average within-group
#' score correlation and average AUC. Coordinates outside the published table
#' are bounded to its range and recorded in the result.
#'
#' `method = "delong"` uses the nonparametric DeLong covariance estimator and
#' handles ties with half credit. Both methods test the two-sided null
#' hypothesis that the AUC difference is zero.
#'
#' @param x A `roc` object.
#' @param curves Character vector naming exactly two marker or MLR curves. The
#'   reported difference is `curves[1] - curves[2]`.
#' @param method Paired comparison method: `"ep24"` (default) or `"delong"`.
#' @param conf.level Confidence level for the normal-approximation interval.
#' @param rating_method Correlation method for EP24 paired scores: `"pearson"`
#'   for interval-scale laboratory results (default) or `"kendall"` for
#'   ordinal ratings. Ignored for DeLong.
#' @param name Optional unique name for storing the comparison. By default it
#'   is generated from the two curve names and method.
#' @param ... Additional arguments (currently unused).
#'
#' @return The updated `roc` object, invisibly. The comparison is stored in
#'   `x$auc_comparisons[[name]]` as a `roc_auc_comparison` object.
#'
#' @examples
#' obj <- roc(ivd_roc_example, c("x1", "x2"), "ref")
#' obj <- auc(obj)
#' obj <- auc_compare(obj, c("x1", "x2"), method = "ep24")
#' obj$auc_comparisons[[1L]]
#'
#' obj <- mlr(obj, name = "combined")
#' obj <- auc_compare(obj, c("x1", "combined"), method = "delong")
#' @export
auc_compare.roc <- function(x, curves, method = c("ep24", "delong"),
                            conf.level = 0.95,
                            rating_method = c("pearson", "kendall"),
                            name = NULL, ...) {
  if (!inherits(x, "roc"))
    stop("Input must be a 'roc' object.", call. = FALSE)
  if (!is.character(curves) || length(curves) != 2L || anyNA(curves) ||
      any(!nzchar(curves)) || curves[1L] == curves[2L])
    stop("'curves' must contain two distinct non-empty curve names.",
         call. = FALSE)
  method <- match.arg(method)
  rating_method <- match.arg(rating_method)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number in (0, 1).", call. = FALSE)
  if (x$n_pos < 2L || x$n_neg < 2L)
    stop("Paired AUC comparison requires at least two positive and two ",
         "negative subjects.", call. = FALSE)

  entries <- lapply(curves, function(curve) .roc_comparison_entry(x, curve))
  names(entries) <- curves
  truth <- x$reference_binary
  scores <- lapply(entries, `[[`, "score")
  if (any(vapply(scores, length, integer(1)) != length(truth)))
    stop("Compared curves do not contain the same subjects.", call. = FALSE)

  rating <- lookup <- covariance <- NULL
  if (method == "ep24") {
    auc_values <- vapply(entries, `[[`, numeric(1), "auc")
    auc_se <- vapply(entries, `[[`, numeric(1), "se")
    if (any(!is.finite(auc_se)))
      stop("EP24 comparison requires finite single-curve Hanley-McNeil ",
           "standard errors; perfect AUCs are not supported.", call. = FALSE)
    group_correlation <- function(group) {
      suppressWarnings(stats::cor(scores[[1L]][group], scores[[2L]][group],
                                  method = rating_method))
    }
    r_neg <- group_correlation(truth == 0L)
    r_pos <- group_correlation(truth == 1L)
    if (!is.finite(r_neg) || !is.finite(r_pos))
      stop("Within-group score correlations could not be estimated.",
           call. = FALSE)
    rating <- c(negative = r_neg, positive = r_pos,
                average = mean(c(r_neg, r_pos)))
    lookup <- .ep24_auc_correlation(rating[["average"]], mean(auc_values))
    area_correlation <- lookup$value
    variance <- auc_se[1L]^2 + auc_se[2L]^2 -
      2 * area_correlation * auc_se[1L] * auc_se[2L]
    method_label <- "EP24-A2 Hanley-McNeil paired approximation"
  } else {
    delong <- .delong_components(scores, truth)
    auc_values <- delong$auc
    covariance <- delong$covariance
    auc_se <- sqrt(diag(covariance))
    contrast <- c(1, -1)
    variance <- drop(t(contrast) %*% covariance %*% contrast)
    area_correlation <- covariance[1L, 2L] /
      sqrt(covariance[1L, 1L] * covariance[2L, 2L])
    method_label <- "DeLong nonparametric paired covariance"
  }
  if (!is.finite(variance) || variance <= 0)
    stop("The estimated variance of the paired AUC difference is not ",
         "positive.", call. = FALSE)

  difference <- unname(auc_values[1L] - auc_values[2L])
  se <- unname(sqrt(variance))
  statistic <- difference / se
  p_value <- 2 * stats::pnorm(-abs(statistic))
  critical <- stats::qnorm(1 - (1 - conf.level) / 2)
  conf_int <- difference + c(-1, 1) * critical * se

  if (is.null(name)) {
    name <- paste(curves[1L], "vs", curves[2L], method, sep = "_")
  } else if (!is.character(name) || length(name) != 1L || is.na(name) ||
             !nzchar(name)) {
    stop("'name' must be NULL or one non-empty character string.",
         call. = FALSE)
  }
  comparisons <- x$auc_comparisons %||% list()
  if (name %in% names(comparisons))
    stop("AUC comparison '", name, "' already exists; use a different name.",
         call. = FALSE)

  result <- structure(list(
    name = name,
    curves = curves,
    curve1 = curves[1L],
    curve2 = curves[2L],
    sources = vapply(entries, `[[`, character(1), "source"),
    directions = vapply(entries, `[[`, character(1), "direction"),
    method = method,
    method_label = method_label,
    rating_method = if (method == "ep24") rating_method else NULL,
    n = x$n,
    n_pos = x$n_pos,
    n_neg = x$n_neg,
    auc1 = unname(auc_values[1L]),
    auc2 = unname(auc_values[2L]),
    auc_se = setNames(unname(auc_se), curves),
    difference = difference,
    se = se,
    conf.level = conf.level,
    conf.int = setNames(conf_int, c("lower", "upper")),
    statistic = unname(statistic),
    p.value = unname(p_value),
    rating_correlation = rating,
    area_correlation = unname(area_correlation),
    lookup = lookup,
    covariance = covariance,
    apparent_model_comparison = any(vapply(entries, `[[`, character(1),
                                             "source") == "model")
  ), class = "roc_auc_comparison")

  comparisons[[name]] <- result
  x$auc_comparisons <- comparisons
  print(result)
  invisible(x)
}

#' Optimal Cutoff Analysis (S3 method)
#'
#' For each evaluation column, compute diagnostic metrics at all cutoff values,
#' and identify the optimal cutoff (maximum Youden index, closest to (0,1) corner).
#'
#' @param x    roc object
#' @param cols column names to analyze; default is all
#' @param ...  reserved arguments
#'
#' @return Updated roc object with results stored in $cutoff_table
#'
#' @examples
#'   obj <- roc(ivd_roc_example, cols = c("x1", "x2"), reference = "ref")
#'   obj <- auc(obj)
#'   cutoff(obj)
#'   cutoff(obj, cols = "x1")
#' @export
cutoff.roc <- function(x, cols = NULL, ...) {
  if (!inherits(x, "roc"))
    stop("Input must be a 'roc' object.", call. = FALSE)

  if (length(x$auc_list) == 0L)
    stop("Run auc() first before cutoff().", call. = FALSE)

  if (is.null(cols)) {
    cols <- x$cols
  } else {
    bad <- setdiff(cols, x$cols)
    if (length(bad) > 0)
      stop("Column(s) not in roc object: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # Compute cutoff table for each col
  all_tables <- list()

  for (col in cols) {
    values    <- x$cols_values[[col]]
    truth     <- x$reference_binary
    direction <- x$col_direction[[col]] %||% "geq"

    curve_df <- .roc_curve(values, truth, direction = direction)
    curve_df$column <- col

    # Select clinically reportable finite cutoffs. Infinite sentinels remain in
    # the table solely to complete the ROC curve endpoints.
    optimal_idx <- .roc_optimal_indices(curve_df)
    youden_max_idx <- optimal_idx$youden
    corner_min_idx <- optimal_idx$corner

    # Mark optimal rows
    curve_df$optimal <- ""
    if (!is.na(youden_max_idx)) curve_df$optimal[youden_max_idx] <- "Youden"
    if (!is.na(corner_min_idx) && corner_min_idx == youden_max_idx) {
      curve_df$optimal[youden_max_idx] <- "Youden / Corner"
    } else if (!is.na(corner_min_idx)) {
      curve_df$optimal[corner_min_idx] <- "Corner"
    }

    # Reorder columns
    curve_df <- curve_df[, c("column", "cutoff", "sensitivity", "specificity",
                              "youden", "tp", "fp", "tn", "fn",
                              "precision", "npv", "accuracy", "optimal")]

    all_tables[[col]] <- curve_df
  }

  cutoff_table <- do.call(rbind, all_tables)
  rownames(cutoff_table) <- NULL

  x$cutoff_table <- structure(cutoff_table, class = c("roc_cutoff_table", "data.frame"),
                               all_tables = all_tables)
  print(x$cutoff_table)


  invisible(x)
}

#' Multivariate Logistic Regression (S3 method)
#'
#' Fit logistic regression using selected evaluation columns,
#' compute AUC of predicted probabilities.
#' Names must be unique; auto-generated as MLR01, MLR02 ... when not specified.
#'
#' @param x    roc object
#' @param cols column names to fit; default is all
#' @param name fit name (must be unique); auto-generated by default
#' @param ...  additional arguments passed to glm
#'
#' @return Updated roc object with results stored in $mlr_results[[name]]
#'
#' @examples
#'   obj <- roc(ivd_roc_example, cols = c("x1", "x2"), reference = "ref")
#'   obj <- auc(obj)
#'   obj <- mlr(obj)                     # -> MLR01
#'   obj <- mlr(obj)                     # -> MLR02
#'   obj <- mlr(obj, name = "my_model")  # -> my_model
#' @export
mlr.roc <- function(x, cols = NULL, name = NULL, ...) {
  if (!inherits(x, "roc"))
    stop("Input must be a 'roc' object.", call. = FALSE)

  if (length(x$auc_list) == 0L)
    stop("Run auc() first before mlr().", call. = FALSE)

  if (is.null(cols)) {
    cols <- x$cols
  } else {
    bad <- setdiff(cols, x$cols)
    if (length(bad) > 0) stop("Column(s) not in roc object: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  if (length(cols) < 1L)
    stop("Need at least one column for logistic regression.", call. = FALSE)

  # Generate or validate name
  existing_names <- names(x$mlr_results)

  if (is.null(name)) {
    if (length(existing_names) == 0) {
      name <- "MLR01"
    } else {
      nums <- as.integer(gsub("^MLR", "", existing_names))
      nums <- nums[!is.na(nums)]
      next_num <- if (length(nums) > 0) max(nums) + 1L else 1L
      name <- sprintf("MLR%02d", next_num)
      # Avoid collision (if user manually created MLR99 and MLR01 already exists)
      while (name %in% existing_names) {
        next_num <- next_num + 1L
        name <- sprintf("MLR%02d", next_num)
      }
    }
  } else {
    if (!is.character(name) || length(name) != 1L || nchar(name) == 0)
      stop("'name' must be a non-empty character string.", call. = FALSE)
    if (name %in% x$cols)
      stop("MLR result name '", name, "' conflicts with an evaluation column. ",
           "Please use a different name.", call. = FALSE)
    if (name %in% existing_names)
      stop("MLR result '", name, "' already exists. Please use a different name.",
           call. = FALSE)
  }

  # Fit Logistic Regression
  formula <- as.formula(paste("reference_binary ~", paste(cols, collapse = " + ")))
  model_data <- as.data.frame(x$cols_values[cols])
  model_data$reference_binary <- x$reference_binary

  fit <- glm(formula, data = model_data, family = binomial, ...)

  # Compute ROC and AUC of predicted probabilities
  predicted_probs <- predict(fit, type = "response")
  mlr_curve <- .roc_curve(predicted_probs, x$reference_binary, direction = "geq")
  mlr_auc <- .auc(mlr_curve$sensitivity, mlr_curve$specificity)
  mlr_se <- .auc_se(mlr_auc, x$n_pos, x$n_neg)
  mlr_ci_lower <- if (!is.na(mlr_se)) max(0, mlr_auc - 1.96 * mlr_se) else NA_real_
  mlr_ci_upper <- if (!is.na(mlr_se)) min(1, mlr_auc + 1.96 * mlr_se) else NA_real_

  # Store results
  mlr_entry <- list(
    name        = name,
    formula     = formula,
    cols        = cols,
    fit         = fit,
    fit_summary = summary(fit),
    predicted   = predicted_probs,
    auc         = mlr_auc,
    auc_se      = mlr_se,
    ci_lower    = mlr_ci_lower,
    ci_upper    = mlr_ci_upper,
    roc_curve   = mlr_curve
  )

  x$mlr_results[[name]] <- structure(mlr_entry, class = "roc_mlr")
  print(x$mlr_results[[name]])
  invisible(x)
}

#' Predict (S3 method) -- predict from multivariate logistic regression
#'
#' Predict positive-class probability for new data using a stored MLR model.
#' Returns the input predictor columns together with predicted probability,
#' standard error, confidence interval, and 0/1 classification.
#'
#' @param object     roc object (must run mlr() first)
#' @param mlr        MLR model name (character). If only one MLR exists, defaults to it;
#'                   if multiple exist, must specify.
#' @param newdata    data.frame of new observations. If NULL (default), uses the
#'                   training data (complete cases) from the roc object.
#' @param column_map Named character vector mapping newdata column names to training
#'                   column names, e.g. \code{c(x1_new = "x1", x2_new = "x2")}.
#'                   If NULL, tries to match by name automatically.
#' @param ...        Additional arguments passed to stats::predict.glm
#'
#' @return data.frame containing the predictor columns used in the model, plus:
#'   \item{pred_prob}{predicted probability of the positive class}
#'   \item{pred_se}{standard error of the predicted probability}
#'   \item{ci_lower}{lower bound of the confidence interval (probability scale)}
#'   \item{ci_upper}{upper bound of the confidence interval (probability scale)}
#'   \item{pred_class}{binary prediction (1 if pred_prob >= 0.5, 0 otherwise)}
#'
#' @examples
#'   df <- data.frame(sid = 1:100,
#'                    x1  = c(rnorm(50, 10, 2), rnorm(50, 12, 2)),
#'                    x2  = rnorm(100, 5, 1),
#'                    ref = rep(c(0, 1), each = 50))
#'   obj <- roc(df, cols = c("x1", "x2"), reference = "ref")
#'   obj <- auc(obj)
#'   obj <- mlr(obj)
#'   predict(obj)
#'   predict(obj, newdata = df[1:5, ])
#'   ## column_map example: rename x1 to new_x1 in newdata
#'   new_df <- df[1:5, c("x2", "ref")]
#'   new_df$new_x1 <- df$x1[1:5]
#'   predict(obj, newdata = new_df,
#'           column_map = c(new_x1 = "x1", x2 = "x2"))
#' @export
predict.roc <- function(object, mlr = NULL, newdata = NULL,
                         column_map = NULL, ...) {

  if (!inherits(object, "roc"))
    stop("Input must be a 'roc' object.", call. = FALSE)

  if (length(object$mlr_results) == 0L)
    stop("No MLR results found. Run mlr() first.", call. = FALSE)

  # Resolve mlr name
  mlr_names <- names(object$mlr_results)
  if (is.null(mlr)) {
    if (length(mlr_names) == 1L) {
      mlr <- mlr_names[1L]
    } else {
      stop("'mlr' must be specified when multiple MLR results exist. Available: ",
           paste(mlr_names, collapse = ", "), call. = FALSE)
    }
  } else {
    if (!is.character(mlr) || length(mlr) != 1L || !mlr %in% mlr_names)
      stop("'mlr' must be one of: ", paste(mlr_names, collapse = ", "), call. = FALSE)
  }

  entry <- object$mlr_results[[mlr]]
  fit <- entry$fit
  train_cols <- entry$cols  # training column names used in the model

  # Prepare prediction data
  if (is.null(newdata)) {
    # Use original training data (already complete cases)
    pred_data <- as.data.frame(object$cols_values[train_cols])
  } else {
    if (!is.data.frame(newdata))
      stop("'newdata' must be a data.frame.", call. = FALSE)

    # Build column mapping: newdata names -> training names
    col_map <- list()
    remaining_train <- train_cols

    # Apply explicit column_map first
    if (!is.null(column_map)) {
      if (is.null(names(column_map)) || any(names(column_map) == ""))
        stop("'column_map' must be a named character vector.", call. = FALSE)
      bad <- setdiff(column_map, train_cols)
      if (length(bad) > 0)
        stop("In 'column_map', training column(s) not found: ",
             paste(bad, collapse = ", "), call. = FALSE)
      for (nm in names(column_map)) {
        if (!nm %in% names(newdata))
          stop("Column '", nm, "' in 'column_map' not found in newdata.", call. = FALSE)
        col_map[[column_map[nm]]] <- nm
      }
      remaining_train <- setdiff(remaining_train, column_map)
    }

    # Auto-match remaining by name
    for (col in remaining_train) {
      if (col %in% names(newdata)) {
        col_map[[col]] <- col
      }
    }

    # Check for missing columns
    found <- names(col_map)
    missing <- setdiff(train_cols, found)
    if (length(missing) > 0)
      stop("Required column(s) not found in newdata or column_map: ",
           paste(missing, collapse = ", "), call. = FALSE)

    # Build prediction data in training column order
    pred_data <- as.data.frame(lapply(train_cols, function(col) newdata[[col_map[[col]]]]))
    names(pred_data) <- train_cols
  }

  # Predict on linear predictor scale for proper CI computation
  pred <- predict.glm(fit, newdata = pred_data, type = "link", se.fit = TRUE, ...)
  linear_pred <- pred$fit
  se_linear   <- pred$se.fit

  # CI on linear predictor scale, then transform to probability scale
  z_val <- qnorm(0.975)
  prob      <- plogis(linear_pred)
  ci_lower  <- plogis(linear_pred - z_val * se_linear)
  ci_upper  <- plogis(linear_pred + z_val * se_linear)

  # SE on probability scale via Delta method
  prob_se <- se_linear * prob * (1 - prob)

  # Binary class (threshold 0.5)
  binary_class <- ifelse(prob >= 0.5, 1L, 0L)

  # Assemble result
  result <- cbind(
    pred_data,
    data.frame(
      pred_prob  = prob,
      pred_se    = prob_se,
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      pred_class = binary_class
    )
  )

  result
}

#' Plot ROC or cutoff-dependent diagnostic metrics
#'
#' Draw ROC curves for evaluation columns and/or MLR model predictions, or draw
#' sensitivity and specificity as functions of the cutoff for one evaluation
#' column. Optionally mark optimal cutoff points (requires \code{cutoff.roc()}
#' first).
#'
#' @param x      roc object
#' @param curves ROC curve names to plot. Names can refer to original marker
#'               columns or stored MLR models. \code{NULL} (default) plots all
#'               original marker curves; \code{"all"} plots all original and
#'               MLR curves.
#' @param cutpoint mark optimal cutoff points: \code{NULL} (default, none),
#'               \code{"Youden"}, \code{"Corner"}, or \code{"all"}.
#' @param type \code{"roc"} (default) draws the conventional ROC curve;
#'               \code{"cutoff"} draws sensitivity and specificity against
#'               finite cutoff values for exactly one evaluation column.
#' @param ...    unused; supplying additional arguments is an error.
#'
#' @return invisible ggplot object
#'
#' @examples
#'   obj <- roc(ivd_roc_example, cols = c("x1", "x2"), reference = "ref")
#'   obj <- auc(obj)
#'   plot(obj)
#'   plot(obj, curves = "x1")
#'   # AUC calculation is not needed for the cutoff-dependent metric plot.
#'   plot(roc(ivd_roc_example, cols = "x1", reference = "ref"),
#'        curves = "x1", type = "cutoff")
#'   obj <- auc(obj); obj <- cutoff(obj); obj <- mlr(obj)
#'   plot(obj, curves = "all", cutpoint = "Youden")
#'   plot(obj, curves = "x1", type = "cutoff", cutpoint = "Youden")
#' @export
plot.roc <- function(x,
                     curves = NULL, cutpoint = NULL,
                     type = c("roc", "cutoff"), ...) {
  if (!inherits(x, "roc"))
    stop("Input must be a 'roc' object.", call. = FALSE)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required but not installed.", call. = FALSE)

  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    dot_names[is.na(dot_names) | !nzchar(dot_names)] <- "<unnamed>"
    stop("Unused argument(s): ", paste(dot_names, collapse = ", "),
         call. = FALSE)
  }

  type <- match.arg(type)

  # A single curve namespace is shared by original markers and stored models.
  if (is.null(curves)) {
    curves <- x$cols
  } else {
    if (!is.character(curves) || anyNA(curves) || any(!nzchar(curves)))
      stop("'curves' must be NULL or a character vector of curve names.",
           call. = FALSE)
    if ("all" %in% curves) {
      if (!identical(curves, "all"))
        stop("Use 'curves = \"all\"' alone.", call. = FALSE)
      curves <- c(x$cols, names(x$mlr_results))
    }
  }
  if (!length(curves))
    stop("Nothing to plot. Provide at least one curve in 'curves'.",
         call. = FALSE)
  if (anyDuplicated(curves))
    stop("'curves' must not contain duplicate names.", call. = FALSE)

  curve_sources <- vapply(curves, function(curve) {
    .roc_curve_source(x, curve)
  }, character(1))
  marker_curves <- curves[curve_sources == "marker"]
  model_curves <- curves[curve_sources == "model"]

  # Validate cutoff parameter
  if (!is.null(cutpoint))
    cutpoint <- match.arg(cutpoint, c("Youden", "Corner", "all"))

  # Sensitivity/specificity across finite cutoffs for one raw marker.
  if (type == "cutoff") {
    if (length(curves) != 1L)
      stop("'type = \"cutoff\"' requires exactly one curve.",
           call. = FALSE)
    if (length(model_curves))
      stop("'type = \"cutoff\"' currently supports only original marker curves, ",
           "not MLR model curves.",
           call. = FALSE)

    col <- marker_curves[[1L]]
    direction <- x$col_direction[[col]] %||% "geq"
    curve <- .roc_curve(x$cols_values[[col]], x$reference_binary,
                        direction = direction)
    curve <- curve[is.finite(curve$cutoff), , drop = FALSE]
    curve <- curve[order(curve$cutoff), , drop = FALSE]
    if (!nrow(curve))
      stop("No finite cutoff values are available for plotting.", call. = FALSE)

    cutoff_plot_data <- rbind(
      data.frame(cutoff = curve$cutoff, metric = "Sensitivity",
                 value = curve$sensitivity, stringsAsFactors = FALSE),
      data.frame(cutoff = curve$cutoff, metric = "Specificity",
                 value = curve$specificity, stringsAsFactors = FALSE)
    )
    direction_text <- if (direction == "geq")
      sprintf("%s: positive when result >= cutoff", col) else
      sprintf("%s: positive when result <= cutoff", col)

    p <- ggplot2::ggplot(
      cutoff_plot_data,
      ggplot2::aes(x = .data$cutoff, y = .data$value,
                   color = .data$metric, group = .data$metric)) +
      ggplot2::geom_step(linewidth = 0.85, direction = "hv") +
      ggplot2::labs(
        x = "Cutoff", y = "Diagnostic metric",
        title = "Sensitivity and Specificity vs Cutoff",
        subtitle = direction_text, color = NULL
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5)
      )

    if (!is.null(cutpoint)) {
      if (is.null(x$cutoff_table))
        stop("Run cutoff() first before using 'cutpoint' in plot().",
             call. = FALSE)
      marks <- .roc_cutpoint_rows(x$cutoff_table, col, cutpoint)
      if (!nrow(marks))
        stop("No stored finite cutoff result is available for column '", col,
             "'. Run cutoff() for this column first.", call. = FALSE)
      mark_points <- rbind(
        data.frame(cutoff = marks$cutoff, metric = "Sensitivity",
                   value = marks$sensitivity),
        data.frame(cutoff = marks$cutoff, metric = "Specificity",
                   value = marks$specificity)
      )
      mark_labels <- data.frame(
        cutoff = marks$cutoff,
        value = pmin(0.96, pmax(marks$sensitivity, marks$specificity) + 0.04),
        label = sprintf(
          "%s\ncutoff = %s\nSensitivity = %.3f, Specificity = %.3f",
          marks$method, .format_num(marks$cutoff),
          marks$sensitivity, marks$specificity),
        stringsAsFactors = FALSE
      )
      p <- p +
        ggplot2::geom_vline(
          data = marks, ggplot2::aes(xintercept = .data$cutoff),
          color = "grey35", linetype = "dashed") +
        ggplot2::geom_point(
          data = mark_points,
          ggplot2::aes(x = .data$cutoff, y = .data$value),
          inherit.aes = FALSE, color = "black", size = 2.3)
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = mark_labels,
          ggplot2::aes(x = .data$cutoff, y = .data$value,
                       label = .data$label),
          inherit.aes = FALSE, color = "black", size = 3,
          min.segment.length = 0.2)
      } else {
        p <- p + ggplot2::geom_text(
          data = mark_labels,
          ggplot2::aes(x = .data$cutoff, y = .data$value,
                       label = .data$label),
          inherit.aes = FALSE, color = "black", size = 3, vjust = -0.2)
      }
    }

    print(p)
    return(invisible(p))
  }

  # Original marker ROC plots retain the existing requirement that auc() has
  # been run for every selected marker.
  missing_auc <- setdiff(marker_curves, names(x$auc_list))
  if (length(missing_auc))
    stop("Run auc() first for marker curve(s): ",
         paste(missing_auc, collapse = ", "), ".", call. = FALSE)

  plot_lines <- list()

  if (length(marker_curves)) {
    for (col in marker_curves) {
      a <- x$auc_list[[col]]
      if (is.null(a)) next
      values    <- x$cols_values[[col]]
      direction <- x$col_direction[[col]] %||% "geq"
      curve <- .roc_curve(values, x$reference_binary, direction = direction)
      plot_lines[[length(plot_lines) + 1L]] <- data.frame(
        fpr         = 1 - curve$specificity,
        sensitivity = curve$sensitivity,
        group       = col,
        stringsAsFactors = FALSE
      )
    }
  }

  # Build plot lines: MLR models
  if (length(model_curves)) {
    for (nm in model_curves) {
      entry <- x$mlr_results[[nm]]
      curve <- entry$roc_curve
      label <- sprintf("%s (AUC=%.4f)", nm, entry$auc)
      plot_lines[[length(plot_lines) + 1L]] <- data.frame(
        fpr         = 1 - curve$specificity,
        sensitivity = curve$sensitivity,
        group       = label,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(plot_lines) == 0L)
    stop("No data to plot.", call. = FALSE)

  plot_data <- do.call(rbind, plot_lines)

  # Determine title
  if (length(model_curves) > 0 && !length(marker_curves))
    title <- "MLR ROC Curve"
  else
    title <- "ROC Curve"

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$fpr, y = .data$sensitivity,
                                                color = .data$group)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dotted", color = "gray60", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(
      x        = "1 - Specificity (False Positive Rate)",
      y        = "Sensitivity (True Positive Rate)",
      title    = title,
      color    = NULL
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.text       = ggplot2::element_text(size = 10),
      plot.margin       = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
      plot.title        = ggplot2::element_text(hjust = 0.5)
    )

  # Cutoff annotations
  if (!is.null(cutpoint) && length(marker_curves)) {
    if (is.null(x$cutoff_table))
      stop("Run cutoff() first before using 'cutpoint' in plot().", call. = FALSE)

    plot_cutoff_points <- list()

    for (col in marker_curves) {
      marks <- .roc_cutpoint_rows(x$cutoff_table, col, cutpoint)
      if (nrow(marks)) {
        plot_cutoff_points[[length(plot_cutoff_points) + 1L]] <- data.frame(
          fpr = 1 - marks$specificity,
          sensitivity = marks$sensitivity,
          group = col,
          cutoff_lab = sprintf("%s cutoff = %s (%s)", col,
                               .format_num(marks$cutoff), marks$method),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(plot_cutoff_points) > 0) {
      cutoff_points <- do.call(rbind, plot_cutoff_points)
      p <- p +
        ggplot2::geom_point(data = cutoff_points,
                            ggplot2::aes(x = .data$fpr, y = .data$sensitivity),
                            color = "black", size = 2.5, inherit.aes = FALSE)
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(data = cutoff_points,
                   ggplot2::aes(x = .data$fpr, y = .data$sensitivity,
                                label = .data$cutoff_lab),
                   size = 3, color = "black", inherit.aes = FALSE,
                   min.segment.length = 0.2)
      }
    }
  } else if (!is.null(cutpoint)) {
    stop("'cutpoint' requires at least one original marker curve.",
         call. = FALSE)
  }

  print(p)
  invisible(p)
}

#' Summary (S3 method)
#'
#' Output key information from all analyses performed on the roc object.
#'
#' @param object roc object
#' @param ...    reserved arguments
#'
#' @return invisible roc object
#'
#' @examples
#'   obj <- roc(ivd_roc_example, cols = c("x1", "x2"), reference = "ref")
#'   summary(obj)
#'   obj <- auc(obj)
#'   obj <- cutoff(obj)
#'   obj <- mlr(obj)
#'   summary(obj)
#' @export
summary.roc <- function(object, ...) {
  if (!inherits(object, "roc"))
    stop("Input must be a 'roc' object.", call. = FALSE)

  cat("\nROC Analysis -- Summary\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")

  # Basic info (from $print slot)
  print.roc(object)

  # Describe details
  if (!is.null(object$describe)) print(object$describe)

  # AUC details
  if (length(object$auc_list) > 0 && !is.null(object$auc_table)) print(object$auc_table)

  # Paired AUC comparisons
  if (length(object$auc_comparisons) > 0) {
    for (nm in names(object$auc_comparisons))
      print(object$auc_comparisons[[nm]])
  }

  # Cutoff analysis
  if (!is.null(object$cutoff_table)) print(object$cutoff_table)

  # MLR results
  if (length(object$mlr_results) > 0) {
    cat("MLR Results\n")
    for (nm in names(object$mlr_results))
      print(object$mlr_results[[nm]])
    cat("\n")
  }

  invisible(object)
}
