#' @noRd
#' @keywords  internal
rue_no_split <- function(x_train,
                         x_test,
                         num_trees = 500L,
                         response_name = "label") {

  # First fit models
  rue_args <- classifier_args(
    data = x_train,
    num_trees = num_trees,
    response_name = response_name
  )
  rue_model <- do.call(what = ranger::ranger, args = rue_args)

  # Then get prediction variances
  rue_train <- se_predicted(rue_model, x_train, is_training = TRUE)
  rue_test <- se_predicted(rue_model, x_test, is_training = FALSE)
  return(list(test = rue_test, train = rue_train))
}

#' @title
#' Dataset Shift via Resampling (Prediction) Uncertainty
#'
#' @description
#' Test for no adverse shift via prediction uncertainty for two-sample comparison.
#' The scores are out-of-bag predictions from random forests with the package
#' \pkg{ranger}. The prefix \emph{rue} stands for resampling uncertainty, the
#' relevant notion of outlier. This uncertainty is the standard error of the
#' mean predictions. This assumes that both training and test sets are labeled.
#'
#' @inherit rd_pt return params details
#'
#' @section Notes:
#' For resampling uncertainty, we essentially implement the approach in
#' Schulam & Saria (2019) with random forests. The standard errors of the mean
#' predictions are the underlying scores. Any performant method for 
#' confidence-based out-of-distribution detection can replace random forests,
#' the default in this implementation.
#'
#' @references Kamulete, V. M. (2022). 
#' \emph{Test for non-negligible adverse shifts}.
#' In Uncertainty in Artificial Intelligence. PMLR.
#'
#' @references Schulam, P., & Saria, S. (2019, April).
#' Can you trust this prediction? Auditing pointwise reliability after learning.
#' In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 1022-1031). PMLR.
#'
#' @references Berger, C., Paschali, M., Glocker, B., & Kamnitsas, K. (2021). 
#' Confidence-based Out-of-Distribution Detection: A Comparative Study and Analysis. 
#' arXiv preprint arXiv:2107.02568.
#'
#' @references Gandy, A. (2009).
#' \emph{Sequential implementation of Monte Carlo tests with uniformly bounded resampling risk}.
#' Journal of the American Statistical Association, 104(488), 1504-1511.
#'
#' @references Li, J., & Fine, J. P. (2010).
#' \emph{Weighted area under the receiver operating characteristic curve and its application to gene selection}.
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 59(4), 673-692.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' idx <- sample(nrow(iris), 2 / 3 * nrow(iris))
#' xy_train <- iris[idx, ]
#' xy_test <- iris[-idx, ]
#' iris_test <- rue_pt(xy_train, xy_test, response_name = "Species")
#' str(iris_test)
#' }
#' @family uncertainty
#'
#' @export
rue_pt <- function(x_train,
                   x_test,
                   R = 1e3,
                   sub_ratio = 1 / 2,
                   num_trees = 500L,
                   response_name = "label") {

  # Create scorer function
  scorer <- function(x_train, x_test) {
    rue_no_split(
      x_train = x_train,
      x_test = x_test,
      response_name = response_name,
      num_trees = num_trees
    )
  }

  # Calculate stats for two-sample test
  result <- exchangeable_null(x_train, x_test, scorer, R = R, is_oob = TRUE)
  return(result)
}
