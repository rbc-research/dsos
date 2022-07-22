#' @noRd
#' @keywords  internal
merge_for_c2st <- function(x_train, x_test, response_name = "label") {

  # Get number of observations in splits
  n_train <- nrow(x_train)
  n_test <- nrow(x_test)

  # Create labels for each dataset
  y_train <- as.factor(rep(0, n_train))
  y_test <- as.factor(rep(1, n_test))
  data.table::setDT(x_train)
  data.table::setDT(x_test)
  data.table::set(x_train, j = response_name, value = y_train)
  data.table::set(x_test, j = response_name, value = y_test)

  # Merge datasets
  x_merged <- data.table::rbindlist(list(x_train, x_test))
  return(x_merged)
}

#' @noRd
#' @keywords  internal
memberships_no_split <- function(x_train,
                                 x_test,
                                 num_trees = 198,
                                 response_name = "label") {

  # First fit models
  rf_data <- merge_for_c2st(x_train, x_test, response_name = response_name)
  rf_args <- classifier_args(
    data = rf_data,
    num_trees = num_trees,
    response_name = response_name
  )

  # Then predict
  predictions <- do.call(what = ranger::ranger, args = rf_args)$predictions[, 2]

  # Split scores
  n_train <- nrow(x_train)
  n_test <- nrow(x_test)
  n_all <- n_train + n_test
  os_train <- predictions[1:n_train]
  os_test <- predictions[(n_train + 1L):n_all]
  return(list(train = os_train, test = os_test))
}

#' @noRd
#' @keywords  internal
#' @importFrom stats predict
memberships_after_split <- function(x_train,
                                    x_test,
                                    num_trees = 198,
                                    sub_ratio = 1 / 2,
                                    response_name = "label") {

  # Split datasets
  train_splits <- split_data(x_train, sub_ratio = sub_ratio)
  test_splits <- split_data(x_test, sub_ratio = sub_ratio)

  # First fit models
  rf_data <- merge_for_c2st(
    train_splits$first,
    test_splits$first,
    response_name = response_name
  )
  rf_args <- classifier_args(
    data = rf_data,
    num_trees = num_trees,
    response_name = response_name
  )
  rf <- do.call(what = ranger::ranger, args = rf_args)

  # Get prediction scores
  os_train <- predict(rf, train_splits$second)$predictions[, 2]
  os_test <- predict(rf, test_splits$second)$predictions[, 2]
  return(list(train = os_train, test = os_test))
}

#' @title
#' Dataset Shift via Class Probabilities
#'
#' @description
#' Test for no adverse shift via class probabilities for two-sample comparison.
#' The scores are out-of-bag predictions from random forests with the package
#' \pkg{ranger}. The prefix \emph{cp} stands for class probability, whether
#' the instance belongs to the training or test set. The probability of
#' belonging to the test set is the relevant notion of outlyingness.
#'
#' @param x_train Training sample.
#' @param x_test Test sample.
#' @param R The number of permutations. May be ignored.
#' @param num_trees The number of trees in random forests.
#' @param sub_ratio Subsampling ratio for sample splitting. May be ignored.
#'
#' @return
#' A named list or object of class \code{outlier.test} containing:
#' \itemize{
#'    \item \code{statistic}: observed WAUC statistic
#'    \item \code{seq_mct}: sequential Monte Carlo test, if applicable
#'    \item \code{p_value}: p-value
#'    \item \code{outlier_scores}: outlier scores from training and test set
#' }
#'
#' @details
#' The suffix \emph{at} refers to the asymptotic test statistic. This variant
#' uses the asymptotic null distribution for the weighted AUC (WAUC), the test
#' statistic. Li & Fine (2010) derives its null distribution. This approximation
#' is reliable in large sample; otherwise, prefer permutations for inference.
#' The example below uses datasets with small samples, which is generally not
#' advisable, is for illustration only.
#'
#' @section Notes:
#' Please see references for the classifier two-sample test, the
#' inspiration behind this approach. Note that Ciemencon et al. (2009) uses
#' both sample splitting for inference and the AUC, rather than the WAUC. Most
#' supervised method for binary classification can replace random forests, the
#' default in this implementation.
#'
#' @references Kamulete, V. M. (2022). 
#' \emph{Test for non-negligible adverse shifts}.
#' In Uncertainty in Artificial Intelligence. PMLR.
#'
#' @references Ciemencon, S., Depecker, M., & Vayatis, N. (2009, December).
#' \emph{AUC optimization and the two-sample problem}.
#' In Proceedings of the 22nd International Conference on Neural Information Processing Systems (pp. 360-368).
#'
#' @references Lopez-Paz, D., & Oquab, M. (2016).
#' \emph{Revisiting classifier two-sample tests}.
#' arXiv preprint arXiv:1610.06545.
#'
#' @references Friedman, J. (2004).
#' \emph{On multivariate goodness-of-fit and two-sample testing}.
#'
#' @references Gandy, A. (2009).
#' \emph{Sequential implementation of Monte Carlo tests with uniformly bounded resampling risk}.
#' Journal of the American Statistical Association, 104(488), 1504-1511.
#'
#' @references Li, J., & Fine, J. P. (2010).
#' \emph{Weighted area under the receiver operating characteristic curve and its application to gene selection}.
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 59(4), 673-692.
#'
#' @references Rinaldo, A., Wasserman, L., & G'Sell, M. (2019).
#' \emph{Bootstrapping and sample splitting for high-dimensional, assumption-lean inference}.
#' Annals of Statistics, 47(6), 3438-3469.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' x_train <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' x_test <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' iris_test <- cp_at(x_train, x_test) # Can also use: cp_ss and cp_pt
#' str(iris_test)
#' }
#' @family classifiers
#'
#' @seealso
#' [cp_ss()] for asymptotic p-value via sample splitting.
#' [cp_pt()] for p-value approximation via permutations.
#'
#' @export
cp_at <- function(x_train,
                  x_test,
                  R = 1e3,
                  num_trees = 500,
                  sub_ratio = 1 / 2) {

  # Create scorer function
  scorer <- function(x_train, x_test) {
    memberships_no_split(
      x_train = x_train,
      x_test = x_test,
      num_trees = num_trees
    )
  }

  # Calculate stats for two-sample test
  result <- asymptotic_null(x_train, x_test, scorer)
  return(result)
}

#' @inherit cp_at title return examples params references description
#' @inheritSection cp_at Notes
#'
#' @details
#' The empirical null distribution uses \code{R} permutations to estimate
#' the p-value. For speed, this is implemented as a sequential Monte Carlo test
#' with the \pkg{simctest} package. See Gandy (2009) for details. The suffix
#' \emph{pt} refers to permutation test. It does not use the asymptotic
#' (theoretical) null distribution for the weighted AUC (WAUC), the test
#' statistic. This is the recommended approach for small samples.
#'
#' @family classifiers
#' @export
cp_pt <- function(x_train,
                  x_test,
                  R = 1e3,
                  num_trees = 500,
                  sub_ratio = 1 / 2) {

  # Create scorer function
  scorer <- function(x_train, x_test) {
    memberships_no_split(
      x_train = x_train,
      x_test = x_test,
      num_trees = num_trees
    )
  }

  # Calculate stats for two-sample test
  result <- exchangeable_null(x_train, x_test, scorer = scorer, R = R, is_oob = TRUE)
  return(result)
}

#' @inherit cp_at title return examples params references description
#' @inheritSection cp_at Notes
#'
#' @details
#' This approach uses sample splitting to compute the p-value for inference.
#' \code{sub_ratio} splits each sample into two parts: one half for estimation
#' (calibration) and the half for inference, if \code{sub_ratio} is 1/2. In other
#' words, it sacrifices some predictive accuracy for inferential robustness as
#' in Rinaldo et al. (2019). The suffix \emph{ss} refers to sample splitting.
#' Sample splitting relies on the asymptotic null distribution for the weighted
#' AUC (WAUC), the test statistic. Li & Fine (2010) derives its null
#' distribution.
#'
#' @family classifiers
#' @export
cp_ss <- function(x_train,
                  x_test,
                  sub_ratio = 1 / 2,
                  R = 1e3,
                  num_trees = 500) {

  # Create scorer function
  scorer <- function(x_train, x_test) {
    memberships_after_split(
      x_train = x_train,
      x_test = x_test,
      num_trees = num_trees,
      sub_ratio = sub_ratio
    )
  }

  # Calculate stats for two-sample test
  result <- asymptotic_null(x_train, x_test, scorer)
  return(result)
}
