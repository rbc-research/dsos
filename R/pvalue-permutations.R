#' @noRd
#' @keywords  internal
smc_pvalue <- function(test_stat, permute_fn, n_pt = 2e3) {

  # Generate null distribution via SMC
  gen_smc <- function() {
    gen_stat <- permute_fn()
    return(gen_stat >= test_stat)
  }
  res <- simctest::simctest(gen_smc, maxsteps = n_pt)
  p_value <- res@pos / res@steps
  dsos_test <- list()
  dsos_test[["seq_mct"]] <- res
  dsos_test[["statistic"]] <- test_stat
  dsos_test[["p_value"]] <- p_value
  dsos_test[["outlier_scores"]] <- NULL
  class(dsos_test) <- "outlier.test"
  return(dsos_test)
}

#' @noRd
#' @keywords  internal
exchangeable_null <- function(x_train,
                              x_test,
                              scorer,
                              is_oob = TRUE,
                              n_pt = 2e3) {
  # Set as data.tables
  x_train <- data.table::as.data.table(x_train)
  x_test <- data.table::as.data.table(x_test)

  # Get observed wauc and helper functions
  helper <- wauc_helper(x_train, x_test, scorer, weight = NULL)
  test_stat <- helper$test_stat
  if (is_oob) {
    permute_fn <- helper$permute_os_fn
  } else {
    permute_fn <- helper$permute_data_fn
  }

  # Gather test info
  dsos_test <- smc_pvalue(test_stat, permute_fn, n_pt)
  dsos_test[["outlier_scores"]] <- helper$os_list
  return(dsos_test)
}

#' @title
#' Permutation Test With Out-Of-Bag Scores
#'
#' @description
#' Test for no adverse shift with outlier scores. Like goodness-of-fit testing,
#' this two-sample comparison takes the training set, \code{x_train} or
#' \code{os_train}, as the reference. The method checks whether the test set,
#' \code{x_test} or \code{os_test}, is worse off relative to this reference
#' set.
#'
#' @param x_train Training (reference/validation) sample.
#' @param x_test Test sample.
#' @param scorer Function which returns a named list with outlier scores from
#' the training and test sample. The first argument to \code{scorer} must be
#' \code{x_train}; the second, \code{x_test}. The returned named list contains
#' two elements: \emph{train} and \emph{test}, each of which is a vector of
#' (outlier) scores. See notes below for more information.
#' @param n_pt The number of permutations.
#'
#' @return
#' A named list of class \code{outlier.test} containing:
#' \itemize{
#'    \item \code{statistic}: observed WAUC statistic
#'    \item \code{seq_mct}: sequential Monte Carlo test, when applicable
#'    \item \code{p_value}: p-value
#'    \item \code{outlier_scores}: outlier scores from training and test set
#' }
#'
#' @references Kamulete, V. M. (2022).
#' \emph{Test for non-negligible adverse shifts}.
#' In The 38th Conference on Uncertainty in Artificial Intelligence. PMLR.
#'
#' @references Gandy, A. (2009).
#' \emph{Sequential implementation of Monte Carlo tests with uniformly bounded resampling risk}.
#' Journal of the American Statistical Association, 104(488), 1504-1511.
#'
#' @references Li, J., & Fine, J. P. (2010).
#' \emph{Weighted area under the receiver operating characteristic curve and its application to gene selection}.
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 59(4), 673-692.
#'
#' @details
#' The empirical null distribution uses \code{n_pt} permutations to estimate
#' the p-value. For speed, this is implemented as a sequential Monte Carlo test
#' with the \pkg{simctest} package. See Gandy (2009) for details. The prefix
#' \emph{pt} refers to permutation test. This approach does not use the
#' asymptotic null distribution for the weighted AUC (WAUC), the test
#' statistic. This is the recommended approach for small samples.
#'
#' @section Notes:
#' The scoring function, \code{scorer}, predicts out-of-bag scores to mimic
#' out-of-sample behaviour. The suffix \emph{oob} stands for out-of-bag to
#' highlight this point. This out-of-bag variant avoids refitting the
#' underlying algorithm from \code{scorer} at every permutation. It can, as a
#' result, be computationally appealing.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' idx <- sample(nrow(iris), 2 / 3 * nrow(iris))
#' xy_train <- iris[idx, ]
#' xy_test <- iris[-idx, ]
#'
#' # First example: residual diagnostics
#' scorer <- function(tr, te) list(train=runif(nrow(tr)), test=runif(nrow(te)))
#' pt_test <- pt_oob(xy_train, xy_test, scorer = scorer)
#' pt_test
#' }
#'
#' @family permutation-test
#'
#' @seealso
#' [pt_refit()] for (slower) p-value approximation via refitting.
#' [at_oob()] for p-value approximation from asymptotic null distribution.
#'
#' @export
pt_oob <- function(x_train, x_test, scorer, n_pt = 2e3) {
  result <- exchangeable_null(
    x_train,
    x_test,
    scorer = scorer,
    n_pt = n_pt,
    is_oob = TRUE
  )
  return(result)
}

#' @title
#' Permutation Test By Refitting
#'
#' @inherit pt_oob description return references details
#' @inheritParams pt_oob
#'
#' @section Notes:
#' The scoring function, \code{scorer}, predicts out-of-sample scores by
#' refitting the underlying algorithm from \code{scorer} at every permutation
#' The suffix \emph{refit} emphasizes this point. This is in contrast to the
#' out-of-bag variant, \code{pt_oob}, which only fits once. This method can be
#' be computationally expensive.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' setosa <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' versicolor <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' scorer <- function(tr, te) list(train=runif(nrow(tr)), test=runif(nrow(te)))
#' pt_test <- pt_refit(setosa, versicolor, scorer = scorer)
#' pt_test
#' }
#'
#' @family permutation-test
#'
#' @seealso
#' [pt_oob()] for (faster) p-value approximation via out-of-bag predictions.
#' [at_oob()] for p-value approximation from asymptotic null distribution.
#'
#' @export
pt_refit <- function(x_train, x_test, scorer, n_pt = 2e3) {
  result <- exchangeable_null(
    x_train,
    x_test,
    scorer = scorer,
    n_pt = n_pt,
    is_oob = FALSE
  )
  return(result)
}

#' @title
#' Permutation Test from Outlier Scores
#'
#' @param os_train Outlier scores in training (reference) set.
#' @param os_test Outlier scores in test set.
#' @param n_pt The number of permutations.
#'
#' @inherit at_from_os description return references
#' @inherit pt_oob details
#' @inheritSection at_from_os Notes
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' os_train <- rnorm(n = 100)
#' os_test <- rnorm(n = 100)
#' null_test <- pt_from_os(os_train, os_test)
#' null_test
#' }
#'
#' @family permutation-test
#'
#' @seealso
#' [pt_oob()] for variant requiring a scoring function.
#' [at_from_os()] for asymptotic test with the outlier scores.
#'
#' @export
pt_from_os <- function(os_train, os_test, n_pt = 2e3) {
  wauc_stat <- wauc_from_os(os_train, os_test)
  # Create function to permute test statistic
  n_test <- length(os_test)
  pooled_os <- c(os_train, os_test)
  permute_fn <- function() permute_from_os(pooled_os, n_test)
  dsos_test <- smc_pvalue(wauc_stat, permute_fn, n_pt)
  dsos_test[["outlier_scores"]] <- list(train = os_train, test = os_test)
  return(dsos_test)
}