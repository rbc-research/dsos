#' @noRd
#' @keywords  internal
weights_from_contamination_rates <- function(contamination_rate) {
  weight <- (contamination_rate - 1)^2
  return(weight)
}

#' @noRd
#' @keywords  internal
weight_function <- function(os_train) {
  # Create empirical cumulative distribution function
  fn <- stats::ecdf(os_train)

  # Create function closure
  weighter <- function(os) {
    contamination_rate <- 1 - fn(os)
    weights <- weights_from_contamination_rates(contamination_rate)
    return(weights)
  }

  return(weighter)
}

#' @noRd
#' @keywords  internal
wauc <- function(y, prob, weighter) {
  # Handle edge cases
  n_prob <- length(unique(prob))
  if (n_prob < 2L) {
    warning("The same outlier score is assigned to all observations.")
    return(1 / 12)
  }

  # Get TPR and FPR
  tpr_fpr <- WeightedROC::WeightedROC(guess = prob, label = y)

  # Get threshold weights
  tpr_fpr$width_factor <- weighter(tpr_fpr$threshold)

  # Calculate widths
  right <- tpr_fpr[-nrow(tpr_fpr), ]
  left <- tpr_fpr[-1, ]
  width <- right$FPR - left$FPR
  width_factor <- right$width_factor
  adj_width <- width_factor * width

  # Calculate areas
  rect_area <- left$TPR * adj_width
  triangle_height <- right$TPR - left$TPR
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
#' In Uncertainty in Artificial Intelligence. PMLR.
#'
#' @export
wauc_from_os <- function(os_train, os_test) {
  # Pool instance labels and weights
  actual_pooled <- c(
    rep(0L, length(os_train)),
    rep(1L, length(os_test))
  )

  # Pool outlier scores
  predicted_pooled <- c(os_train, os_test)

  # Create function to weigh thresholds
  weighter <- weight_function(os_train)

  # Calculate WAUC
  wauc_stat <- wauc(actual_pooled, predicted_pooled, weighter)
  return(wauc_stat)
}

#' @noRd
#' @keywords  internal
wauc_and_os <- function(os_train, os_test) {

  # Run test on outlier scores
  wauc_stat <- wauc_from_os(
    os_train = os_train,
    os_test = os_test
  )
  os_list <- list(train = os_train, test = os_test)
  result_list <- list(wauc = wauc_stat, outlier_scores = os_list)
  return(result_list)
}

#' @noRd
#' @keywords  internal
wauc_from_data <- function(x_train, x_test, scorer) {
  # Get list of outlier scores
  os_list <- scorer(
    x_train = x_train,
    x_test = x_test
  )

  # Calculate stats for two-sample test
  result_list <- wauc_and_os(os_list[["train"]], os_list[["test"]])
  return(result_list)
}

#' @noRd
#' @keywords  internal
wauc_helper <- function(x_train, x_test, scorer) {

  # Get observed wauc and outlier scores
  observed <- wauc_from_data(x_train, x_test, scorer)
  test_stat <- observed$wauc
  os_list <- observed$outlier_scores

  # Create function to permute test statistic
  n_test <- nrow(x_test)
  pooled_os <- c(os_list$train, os_list$test)
  pooled_data <- data.table::rbindlist(list(x_train, x_test))
  permute_os_fn <- function() permute_wauc_from_os(pooled_os, n_test)
  permute_data_fn <- function() permute_wauc_from_data(pooled_data, n_test, scorer)$wauc

  # Collect info in list
  results <- list()
  results[["permute_os_fn"]] <- permute_os_fn
  results[["permute_data_fn"]] <- permute_data_fn
  results[["test_stat"]] <- test_stat
  results[["os_list"]] <- os_list
  return(results)
}
