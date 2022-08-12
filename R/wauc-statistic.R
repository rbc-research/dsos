#' @noRd
#' @keywords  internal
weights_from_contamined <- function(contamination_rate) {
  weight <- (contamination_rate - 1)^2
  return(weight)
}

#' @noRd
#' @keywords  internal
make_threshold_fn <- function(os_train, weight = rep(1, length(os_train))) {
  # Create empirical cumulative distribution function
  fn <- ewcdf(os_train, weight)

  # Create function closure
  weight_fn <- function(os) {
    contamined <- 1 - fn(os)
    weights <- weights_from_contamined(contamined)
    return(weights)
  }

  return(weight_fn)
}

#' @noRd
#' @keywords  internal
weighted_roc <- function(score, label, weight = rep(1, length(label))) {
  # Copied this function from WeightedROC package. See:
  # https://github.com/tdhock/WeightedROC/blob/master/R/ROC.R
  if (is.factor(label)) {
    label <- as.integer(label)
  }
  stopifnot(is.numeric(label))
  label_tab <- table(label)
  if (length(label_tab) == 1) {
    print(label_tab)
    stop("only one label value")
  }
  if (all(label %in% c(0, 1))) {
    label[label == 0] <- -1
  }
  if (all(label %in% c(1, 2))) {
    label[label == 1] <- -1
    label[label == 2] <- 1
  }
  stopifnot(label %in% c(-1, 1))
  stopifnot(is.numeric(score))
  stopifnot(length(label) == length(score))
  if (any(is.na(score))) {
    stop("ROC curve undefined for NA score")
  }
  stopifnot(is.numeric(weight))
  stopifnot(length(label) == length(weight))
  stopifnot(weight > 0)
  ord <- order(score)
  y <- label[ord]
  w <- weight[ord]
  y_hat <- score[ord]
  is_positive <- y == 1
  is_negative <- y == -1
  w_positive <- w_negative <- w
  w_positive[is_negative] <- 0
  w_negative[is_positive] <- 0
  cum_positive <- cumsum(w_positive)
  cum_negative <- cumsum(w_negative)
  is_end <- c(diff(y_hat) != 0, TRUE)
  n <- length(y)
  threshold <- c(y_hat[is_end], Inf)
  total_positive <- cum_positive[n]
  total_negative <- cum_negative[n]
  fn <- c(0, cum_positive[is_end])
  fnr <- fn / total_positive
  tpr <- 1 - fnr
  tn <- c(0, cum_negative[is_end])
  fp <- total_negative - tn
  fpr <- fp / total_negative
  d <- data.frame(tpr, fpr, threshold, fn, fp)
  return(d)
}

#' @noRd
#' @keywords  internal
wauc <- function(label,
                 score,
                 threshold_fn,
                 weight = rep(1, length(label))) {
  # Handle edge cases: same score is assigned to all observations.
  n_unique <- length(unique(score))
  if (n_unique < 2L) {
    return(1 / 12)
  }

  # Get tpr and fpr
  tpr_fpr <- weighted_roc(score, label, weight = weight)

  # Get threshold weights
  tpr_fpr$width_factor <- threshold_fn(tpr_fpr$threshold)

  # Calculate widths
  right <- tpr_fpr[-nrow(tpr_fpr), ]
  left <- tpr_fpr[-1, ]
  width <- right$fpr - left$fpr
  width_factor <- right$width_factor
  adj_width <- width_factor * width

  # Calculate areas
  rect_area <- left$tpr * adj_width
  triangle_height <- right$tpr - left$tpr
  triangle_area <- triangle_height * adj_width / 2

  # Sum areas
  final_stat <- sum(rect_area, triangle_area)
  return(final_stat)
}

#' @title
#' Weighted AUC from Outlier Scores
#'
#' @description
#' Computes the weighted AUC with the weighting scheme described in
#' Kamulete, V. M. (2021). This assumes that the training set is the reference
#' distribution and specifies a particular functional form to derive weights
#' from threshold scores.
#'
#' @param os_train Outlier scores in training set.
#' @param os_test Outlier scores in test set.
#' @param weight Numeric vector of weights of length
#' \code{length(os_train) + length(os_test)}. The first \code{length(os_train)}
#' weights belongs to the training set, the rest is for the test set. If
#' \code{NULL}, the default, all weights are set to 1.
#'
#' @return
#' The value (scalar) of the weighted AUC given the weighting scheme.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' os_train <- runif(n = 100)
#' os_test <- runif(n = 100)
#' test_stat <- wauc_from_os(os_train, os_test)
#' }
#'
#' @family statistic
#'
#' @references Kamulete, V. M. (2022).
#' \emph{Test for non-negligible adverse shifts}.
#' In The 38th Conference on Uncertainty in Artificial Intelligence. PMLR.
#'
#' @export
wauc_from_os <- function(os_train, os_test, weight = NULL) {
  # Pool instance labels
  n_train <- length(os_train)
  n_test <- length(os_test)
  label <- c(rep(0L, n_train), rep(1L, n_test))

  # Create weights if missing
  if (is.null(weight)) weight <- rep(1, n_train + n_test)
  stopifnot(
    length(weight) == n_train + n_test,
    all(weight >= 0)
  )

  # Pool outlier scores
  score <- c(os_train, os_test)

  # Create function to weigh thresholds
  w_train <- weight[1:n_train]
  fn <- make_threshold_fn(os_train, weight = w_train)

  # Calculate WAUC
  wauc_stat <- wauc(label, score, threshold_fn = fn, weight = weight)
  return(wauc_stat)
}

#' @noRd
#' @keywords  internal
wauc_and_os <- function(os_train, os_test, weight = NULL) {

  # Run test on outlier scores
  wauc_stat <- wauc_from_os(
    os_train = os_train,
    os_test = os_test,
    weight = weight
  )
  os_list <- list(train = os_train, test = os_test)
  result_list <- list(wauc = wauc_stat, outlier_scores = os_list)
  return(result_list)
}

#' @noRd
#' @keywords  internal
wauc_from_data <- function(x_train, x_test, scorer, weight = NULL) {
  # Get list of outlier scores
  os_list <- scorer(x_train, x_test)

  # Calculate stats for two-sample test
  result_list <- wauc_and_os(
    os_list[["train"]],
    os_list[["test"]],
    weight = weight
  )
  return(result_list)
}

#' @noRd
#' @keywords  internal
wauc_helper <- function(x_train, x_test, scorer, weight = NULL) {

  # Get observed wauc and outlier scores
  observed <- wauc_from_data(x_train, x_test, scorer, weight)
  test_stat <- observed$wauc
  os_list <- observed$outlier_scores

  # Create function to permute test statistic
  n_test <- nrow(x_test)
  pooled_os <- c(os_list$train, os_list$test)
  pooled_data <- data.table::rbindlist(list(x_train, x_test))
  os_fn <- function() permute_from_os(pooled_os, n_test)
  data_fn <- function() permute_from_data(pooled_data, n_test, scorer)$wauc

  # Collect info in list
  results <- list()
  results[["permute_os_fn"]] <- os_fn
  results[["permute_data_fn"]] <- data_fn
  results[["test_stat"]] <- test_stat
  results[["os_list"]] <- os_list
  return(results)
}
