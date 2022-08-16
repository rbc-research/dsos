#' @noRd
#' @keywords  internal
classifier_args <- function(data,
                            num_trees = 500L,
                            response_name = "label") {

  # Check inputs
  stopifnot(all(response_name %in% names(data)))

  # Fix hyperparameters
  # Optimized defaults for the Brier score taken from:
  #   https://philipppro.shinyapps.io/tunability/
  num_cols <- ncol(data)
  num_features <- num_cols - 1L
  num_obs <- nrow(data)

  # Return list of function arguments
  list(
    data = data,
    dependent.variable.name = response_name,
    num.trees = num_trees, # 198L value for 'optimal' default
    mtry = ceiling(0.666 * num_features),
    min.node.size = max(ceiling(num_obs^0.110), 10),
    sample.fraction = 0.667,
    importance = "none",
    probability = TRUE,
    write.forest = TRUE,
    replace = FALSE,
    save.memory = FALSE,
    respect.unordered.factors = TRUE,
    verbose = FALSE,
    keep.inbag = TRUE
  )
}

#' @noRd
#' @keywords  internal
residual_args <- function(data, num_trees = 500L, response_name = "label") {

  # Check inputs
  stopifnot(all(response_name %in% names(data)))

  # Fix hyperparameters
  # Optimized defaults for the Brier score taken from:
  #   https://philipppro.shinyapps.io/tunability/
  num_cols <- ncol(data)
  num_features <- num_cols - 1L
  num_obs <- nrow(data)

  # Return list of function arguments
  list(
    data = data,
    dependent.variable.name = response_name,
    num.trees = num_trees, # 198L value for 'optimal' default
    mtry = ceiling(0.666 * num_features),
    min.node.size = max(ceiling(num_obs^0.110), 10),
    sample.fraction = 0.667,
    importance = "none",
    write.forest = TRUE,
    replace = FALSE,
    save.memory = FALSE,
    respect.unordered.factors = TRUE,
    verbose = FALSE,
    keep.inbag = TRUE
  )
}

#' @noRd
#' @keywords  internal
#' @importFrom stats predict
predict_oob <- function(rf_obj, rf_data, is_training = TRUE) {
  if (is_training) {
    preds <- rf_obj$predictions
  } else {
    preds <- predict(rf_obj, rf_data)$predictions
  }
  stopifnot(nrow(preds) == nrow(rf_data))
  return(preds)
}

#' @noRd
#' @keywords  internal
score_predicted <- function(rf_obj, rf_data, is_training = TRUE) {
  preds <- predict_oob(rf_obj, rf_data, is_training = is_training)
  col_idx <- apply(preds, 1, which.max)
  col_names <- colnames(preds)[col_idx]
  class_score <- data.table::data.table(
    row = 1:length(col_idx),
    col = col_names,
    score = apply(preds, 1, max)
  )
  return(class_score)
}

#' @noRd
#' @keywords  internal
residuals_from_rf <- function(rf_obj,
                              rf_data,
                              is_training = TRUE,
                              response_name = "y") {
  preds <- predict_oob(rf_obj, rf_data, is_training = is_training)
  labels <- 1:ncol(preds)
  names(labels) <- colnames(preds)
  col_key <- as.character(rf_data[[response_name]])
  col_idx <- unname(labels[col_key])
  row_idx <- 1:nrow(preds)
  class_prob <- preds[cbind(row_idx, col_idx)]
  residuals <- 1.0 - class_prob
  return(residuals)
}

#' @noRd
#' @keywords  internal
#' @import data.table
oob_indices <- function(rf_obj) {
  # Identify OOB cases
  oob_idx <- ifelse(simplify2array(rf_obj$inbag.counts) == 0, TRUE, FALSE)
  oob_idx <- data.table::as.data.table(oob_idx)
  mapping <- seq_len(ncol(oob_idx))
  names(mapping) <- colnames(oob_idx)
  oob_idx[, row := .I]
  oob_idx <- data.table::melt.data.table(
    data = oob_idx,
    id.vars = "row",
    variable.name = "tree",
    variable.factor = FALSE,
    value.name = "is_oob"
  )
  # To silence NOTEs from R CMD check
  tree <- is_oob <- NULL
  oob_idx <- oob_idx[is_oob == TRUE]
  oob_idx[, is_oob := NULL]
  oob_idx[, tree := mapping[tree]]
  return(oob_idx)
}

#' @noRd
#' @keywords  internal
#' @importFrom stats predict
predict_all_trees <- function(rf_obj, rf_data, is_training = TRUE) {
  predicted <- tree <- row <- . <- NULL
  predicted <- predict(rf_obj, rf_data, predict.all = TRUE)$predictions
  from_trees <- data.table::as.data.table(predicted)
  data.table::setnames(
    from_trees,
    new = c("row", "col", "tree", "predicted")
  )
  if (is_training) {
    oob_idx <- oob_indices(rf_obj)
    from_trees <- from_trees[oob_idx, on = c("row", "tree")]
  }
  selected_df <- score_predicted(rf_obj, rf_data, is_training = is_training)
  from_trees <- from_trees[selected_df, on = c("row", "col")]
  return(from_trees[, .(row, predicted)])
}

#' @noRd
#' @keywords  internal
#' @importFrom stats var
mean_se <- function(x) sqrt(var(x) / length(x))

#' @noRd
#' @keywords  internal
boot_se <- function(x, oob_prop, n_iter = 200) {
  n <- length(x)
  size <- floor(n * oob_prop)
  se_sub <- future.apply::future_replicate(
    n = n_iter,
    expr = {
      x_sub <- sample(x, size = size)
      mean_se(x_sub)
    }
  )
  se_avg <- mean(se_sub)
  return(se_avg)
}

#' @noRd
#' @keywords  internal
dsos_se <- function(x, oob_prop, n_iter = 30, is_training = TRUE) {
  if (is_training) {
    return(mean_se(x))
  }
  return(boot_se(x, oob_prop, n_iter))
}

#' @noRd
#' @keywords  internal
#' @importFrom stats median
oob_ratio <- function(rf_obj) {
  oob_idx <- ifelse(simplify2array(rf_obj$inbag.counts) == 0, TRUE, FALSE)
  oob_per_instance <- rowMeans(oob_idx)
  oob_prop <- median(oob_per_instance)
  return(oob_prop)
}

#' @noRd
#' @keywords  internal
se_predicted <- function(rf_obj, rf_data, is_training = TRUE) {
  se <- predicted <- row <- . <- NULL
  trees_df <- predict_all_trees(rf_obj, rf_data, is_training = is_training)
  oob_prop <- oob_ratio(rf_obj)
  se_df <- trees_df[
    ,
    .(se = dsos_se(predicted, oob_prop = oob_prop, is_training = is_training)),
    by = c("row")
  ]
  return(se_df[order(row), se])
}
