#' @noRd
#' @keywords  internal
#' @importFrom stats predict
outliers_no_split <- function(x_train, x_test, num_trees = 500) {

  # First fit models
  iso_fit <- isotree::isolation.forest(
    data = x_train,
    ntrees = num_trees # 100 in original iForest paper
  )

  # Then predict
  os_train <- predict(iso_fit, newdata = x_train)
  os_test <- predict(iso_fit, newdata = x_test)
  return(list(test = os_test, train = os_train))
}

#' @title
#' Dataset Shift via Isolation Scores
#'
#' @description
#' Test for no adverse shift via isolation scores for two-sample comparison.
#' The scores are predictions from extended isolation forest with the package
#' \pkg{isotree}. The prefix \emph{od} stands for outlier detection, the
#' relevant notion of outlyingness.
#'
#' @inherit cp_pt return params details
#'
#' @section Notes:
#' Isolation forest detects \emph{isolated} points, instances that are typically
#' out-of-distribution relative to the high-density regions of the data distribution.
#' Any performant method for density-based out-of-distribution detection can
#' replace isolation forest, the default in this implementation.
#'
#' @references Kamulete, V. M. (2022). 
#' \emph{Test for non-negligible adverse shifts}.
#' In Uncertainty in Artificial Intelligence. PMLR.
#'
#' @references Liu, F. T., Ting, K. M., & Zhou, Z. H. (2008, December).
#' \emph{Isolation forest}.
#' In 2008 Eighth IEEE International Conference on Data Mining (pp. 413-422).
#' IEEE.
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
#' x_train <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' x_test <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' iris_test <- od_pt(x_train, x_test) # Can also use: od_ss
#' str(iris_test)
#' }
#' @family anomalies
#'
#' @export
od_pt <- function(x_train,
                  x_test,
                  R = 1e3,
                  num_trees = 500,
                  sub_ratio = 1 / 2) {

  # Create scorer function
  scorer <- function(x_train, x_test) {
    outliers_no_split(
      x_train = x_train,
      x_test = x_test,
      num_trees = num_trees
    )
  }

  # Calculate stats for two-sample test
  result <- exchangeable_null(x_train, x_test, scorer, R = R, is_oob = FALSE)
  return(result)
}

