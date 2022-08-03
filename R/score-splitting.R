#' @title
#' Split Samples And Predict Class Probability (Sample Membership)
#'
#' @inherit score_cp description return
#' @inheritParams score_cp
#'
#' @section Notes:
#' See the docs for \code{score_cp} for more information. \code{split_cp} uses
#' sample splitting (half samples), rather than out-of-bag predictions as in
#' \code{score_cp}, for inference. Rinaldo et al. (2019) discusses how sample
#' splitting can be used for valid inference (p-value estimation).
#'
#' @references Rinaldo, A., Wasserman, L., & G’Sell, M. (2019).
#' \emph{Bootstrapping and sample splitting for high-dimensional, assumption-lean inference}.
#' The Annals of Statistics, 47(6), 3438-3469.
#'
#' @details
#' \code{split_cp} fits a classifier to discriminate between training and test
#' sets. It splits training and test sets into half samples so that the first
#' halves are for model fitting and the second halves, for out-of-sample
#' predictions. As a result, estimating the p-value can take advantage of the
#' asymptotic null distribution.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' setosa <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' versicolor <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' outlier_scores <- split_cp(setosa, versicolor, response_name = "label")
#' str(outlier_scores)
#' }
#'
#' @importFrom stats predict
#'
#' @family splitting
#'
#' @seealso
#' [score_cp()] for the out-of-bag variant, rather than sample splitting.
#'
#' @export
split_cp <- function(x_train, x_test, n_trees = 500L, response_name = "label") {

  # Split datasets
  train_splits <- split_data(x_train, sub_ratio = 1 / 2)
  test_splits <- split_data(x_test, sub_ratio = 1 / 2)

  # First fit models
  rf_data <- stack_data(
    train_splits$first,
    test_splits$first,
    response_name = response_name
  )
  rf_args <- classifier_args(
    data = rf_data,
    num_trees = n_trees,
    response_name = response_name
  )
  rf <- do.call(what = ranger::ranger, args = rf_args)

  # Get out-of-bag predictions
  os_train <- predict(rf, train_splits$second)$predictions[, 2]
  os_test <- predict(rf, test_splits$second)$predictions[, 2]
  return(list(train = os_train, test = os_test))
}
