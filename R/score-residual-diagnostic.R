#' @noRd
#' @keywords  internal
residuals_no_split <- function(x_train,
                               x_test,
                               num_trees = 500L,
                               response_name = "label") {

  # First fit models
  rf_args <- classifier_args(
    data = x_train,
    num_trees = num_trees,
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
#' Dataset Shift via Residuals
#'
#' @description
#' Test for no adverse shift via residuals for multivariate two-sample comparison.
#' The scores are obtained using out-of-bag predictions from random forest with
#' the package \pkg{ranger} to get the residuals. The prefix \emph{rd} stands for
#' residual diagnostics, the relevant notion of outlier. This test assumes that
#' both training and test sets are labeled.
#'
#' @inherit cp_pt details return
#'
#' @param x_train Training sample.
#' @param x_test Test sample.
#' @param R The number of permutations. May be ignored.
#' @param response_name The column name of the categorical outcome to predict.
#' @param num_trees The number of trees in random forests.
#' @param sub_ratio Subsampling ratio for sample splitting. May be ignored.
#'
#' @section Notes:
#' Residuals traditionally underpin diagnostics (misspecification) tests in
#' supervised learning. For a contemporaneous example of this approach also using
#' machine learning, see see Janková et al. (2020) and references therein.
#'
#' @references Kamulete, V. M. (2022). 
#' \emph{Test for non-negligible adverse shifts}.
#' In Uncertainty in Artificial Intelligence. PMLR.
#'
#' @references Janková, J., Shah, R. D., Bühlmann, P., & Samworth, R. J. (2020).
#' \emph{Goodness-of-fit testing in high dimensional generalized linear models}.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(3), 773-795.
#'
#' @references Li, J., & Fine, J. P. (2010).
#' \emph{Weighted area under the receiver operating characteristic curve and its application to gene selection}.
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 59(4), 673-692.
#'
#' @references Gandy, A. (2009).
#' \emph{Sequential implementation of Monte Carlo tests with uniformly bounded resampling risk}.
#' Journal of the American Statistical Association, 104(488), 1504-1511.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' idx <- sample(nrow(iris), 2 / 3 * nrow(iris))
#' xy_train <- iris[idx, ]
#' xy_test <- iris[-idx, ]
#' iris_test <- rd_pt(xy_train, xy_test, response_name = "Species")
#' str(iris_test)
#' }
#' @family residuals
#'
#' @export
rd_pt <- function(x_train,
                  x_test,
                  R = 1e3,
                  num_trees = 500L,
                  sub_ratio = 1 / 2,
                  response_name = "label") {

  # Create scorer function
  scorer <- function(x_train, x_test) {
    residuals_no_split(
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
