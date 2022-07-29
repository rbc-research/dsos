#' @title
#' Predict Isolation Scores
#'
#' @description
#' Predict isolation scores using (extended) isolation forest with the
#' \pkg{isotree} package. The prefix \emph{od} stands for outlier detection,
#' the relevant notion of outlyingness. This function is useful to test for
#' dataset shift via density-based scores: isolation scores are inversely
#' related, if not quite proportional, to densities.
#'
#' @param x_train Training (reference) sample.
#' @param x_test Test sample.
#' @param n_trees The number of trees in isolation forest.
#'
#' @return
#' A named list or object of class \code{outlier.test} containing:
#' \itemize{
#'    \item \code{train}: vector of scores in training set
#'    \item \code{test}: vector of scores in test set
#' }
#'
#' @section Notes:
#' Isolation forest detects \emph{isolated} points that are typically
#' out-of-distribution relative to the high-density regions of the data
#' distribution. Any performant method for density-based out-of-distribution
#' detection can replace isolation forest, the default in this implementation.
#'
#' @references Liu, F. T., Ting, K. M., & Zhou, Z. H. (2008, December).
#' \emph{Isolation forest}.
#' In 2008 Eighth IEEE International Conference on Data Mining (pp. 413-422).
#' IEEE.
#'
#' @references Liu, F. T., Ting, K. M., & Zhou, Z. H. (2008, December).
#' \emph{Isolation forest}.
#' In 2008 Eighth IEEE International Conference on Data Mining (pp. 413-422).
#' IEEE.
#'
#' #' @details
#' \code{od_scorer} first fits to the training data and then predict in-sample
#' for this reference sample. Then it predicts out-of-sample for the test set.
#' As a result, estimating p-value via permutations require refitting the
#' algorithm for every permutation.
#'
#' @importFrom stats predict
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' setosa <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' versicolor <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' od_scorer(setosa, versicolor)
#' }
#' @family scoring
#'
#' @export
od_scorer <- function(x_train, x_test, n_trees = 500) {

  # First fit models
  iso_fit <- isotree::isolation.forest(
    data = x_train,
    ntrees = n_trees # 100 in original iForest paper
  )

  # Then predict
  os_train <- predict(iso_fit, newdata = x_train)
  os_test <- predict(iso_fit, newdata = x_test)
  return(list(test = os_test, train = os_train))
}

#' @title
#' Predict Out-of-bag Errors (Residuals)
#'
#' @description
#' Predict out-of-bag errors (residuals) using random forest with the
#' \pkg{ranger} package. The prefix \emph{rd} stands for residual diagnostic,
#' the relevant notion of outlyingness. This function is useful to test for
#' dataset shift via prediction errors from the underlying supervised
#' algorithm, namely random forest in this implementation.
#'
#' @param x_train Training (reference) sample.
#' @param x_test Test sample.
#' @param response_name The column name of the categorical outcome to predict.
#' @param n_trees The number of trees in random forest.
#'
#' @inherit od_scorer return
#'
#' @section Notes:
#' Residuals traditionally underpin diagnostics (misspecification) tests in
#' supervised learning. For a contemporaneous example of this approach using
#' machine learning, see Janková et al. (2020) and references therein.
#'
#' @references Janková, J., Shah, R. D., Bühlmann, P., & Samworth, R. J. (2020).
#' \emph{Goodness-of-fit testing in high dimensional generalized linear models}.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(3), 773-795.
#'
#' @references Hediger, S., Michel, L., & Näf, J. (2022).
#' \emph{On the use of random forest for two-sample testing}.
#' Computational Statistics & Data Analysis, 170, 107435.
#'
#' @details
#' \code{rd_scorer} first fits to the training data and uses out-of-bag
#' predictions to estimate errors (residuals) for this reference sample. Then it
#' leverages out-of-sample predictions to calculate errors for the test set.
#' As a result, estimating p-value via permutations does not require refitting
#' the algorithm for every permutation. See Hediger et al. (2022) for details
#' on this approach.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' idx <- sample(nrow(iris), 2 / 3 * nrow(iris))
#' xy_train <- iris[idx, ]
#' xy_test <- iris[-idx, ]
#' rd_scorer(xy_train, xy_test, n_trees = 500L, response_name = "Species")
#' }
#'
#' @family scoring
#'
#' @export
rd_scorer <- function(x_train, x_test, n_trees = 500L, response_name = "label") {

  # First fit models
  rf_args <- classifier_args(
    data = x_train,
    num_trees = n_trees,
    response_name = response_name
  )
  rf_model <- do.call(what = ranger::ranger, args = rf_args)

  # Then get residuals
  os_train <- residuals_from_rf(
    rf_model,
    x_train,
    response_name = response_name,
    is_training = TRUE
  )
  os_test <- residuals_from_rf(
    rf_model,
    x_test,
    response_name = response_name,
    is_training = FALSE
  )
  return(list(test = os_test, train = os_train))
}

#' @title
#' Estimate Resampling Uncertainty (Prediction Confidence)
#'
#' @description
#' Estimate prediction uncertainty using random forest with the \pkg{ranger}
#' package. The prefix \emph{rue} stands for resampling uncertainty estimation,
#' the relevant notion of outlyingness. This function is useful to test for
#' dataset shift via prediction uncertainty from the underlying supervised
#' algorithm, namely random forest in this implementation.
#'
#' @inherit rd_scorer return param
#'
#' @section Notes:
#' For prediction uncertainty, we essentially implement the approach in
#' Schulam & Saria (2019) with random forest. The standard errors of the mean
#' predictions are the specified outlier scores. Any performant method for
#' confidence-based out-of-distribution detection can replace random forest,
#' the baseline in this implementation. Other methods for confidence-based
#' out-of-distribution detection are discussed in Berger et al. (2021).
#'
#' @references Schulam, P., & Saria, S. (2019, April).
#' \emph{Can you trust this prediction? Auditing pointwise reliability after learning}.
#' In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 1022-1031).
#' PMLR.
#'
#' @references Berger, C., Paschali, M., Glocker, B., & Kamnitsas, K. (2021).
#' \emph{Confidence-based out-of-distribution detection: a comparative study and analysis}.
#' In Uncertainty for Safe Utilization of Machine Learning in Medical Imaging, and Perinatal Imaging, Placental and Preterm Image Analysis (pp. 122-132).
#' Springer, Cham.
#'
#' @details
#' \code{rue_scorer} first fits to the training data and uses out-of-bag
#' predictions to estimate prediction uncertainty for this reference
#' sample. Then it leverages out-of-sample predictions to do the same for
#' the test set. As a result, estimating p-value via permutations does not
#' require refitting the algorithm for every permutation.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' idx <- sample(nrow(iris), 2 / 3 * nrow(iris))
#' xy_train <- iris[idx, ]
#' xy_test <- iris[-idx, ]
#' rue_scorer(xy_train, xy_test, n_trees = 500L, response_name = "Species")
#' }
#'
#' @family scoring
#'
#' @export
rue_scorer <- function(x_train, x_test, n_trees = 500L, response_name = "label") {

  # First fit models
  rue_args <- classifier_args(
    data = x_train,
    num_trees = n_trees,
    response_name = response_name
  )
  rue_model <- do.call(what = ranger::ranger, args = rue_args)

  # Then get prediction uncertainty (standard errors for mean predictions)
  rue_train <- se_predicted(rue_model, x_train, is_training = TRUE)
  rue_test <- se_predicted(rue_model, x_test, is_training = FALSE)
  return(list(test = rue_test, train = rue_train))
}