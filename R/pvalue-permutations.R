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

  # Wrap result in a list
  test_list <- list()
  test_list[["seq_mct"]] <- res
  test_list[["statistic"]] <- test_stat
  test_list[["p_value"]] <- p_value
  test_list[["outlier_scores"]] <- NULL
  class(test_list) <- "outlier.test"
  return(test_list)
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
  helper <- wauc_helper(x_train, x_test, scorer)
  test_stat <- helper$test_stat
  if (is_oob) {
    permute_fn <- helper$permute_os_fn
  } else {
    permute_fn <- helper$permute_data_fn
  }

  # Gather test info
  test_list <- smc_pvalue(test_stat, permute_fn, n_pt)
  test_list[["outlier_scores"]] <- helper$os_list
  return(test_list)
}

#' @title
#' Permutation Test With Out-Of-Bag Scores
#'
#' @description
#' Test for no adverse shift with outlier scores. Like goodness-of-fit testing,
#' this two-sample comparison takes the training set, \code{x_train}, as the
#' the reference. The method checks whether the test set, \code{x_test}, is
#' worse off relative to this reference distribution.
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
#' A named list or object of class \code{outlier.test} containing:
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
#' result, be computionally appealing.
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
#' scorer_1 <- function(x_train, x_test) score_rd(x_train, x_test, response_name = "Species")
#' rd_test <- pt_oob(x_train, x_test, scorer = scorer_1)
#' str(rd_test)
#'
#' # Second example: prediction uncertainty
#' scorer_2 <- function(x_train, x_test) score_rue(x_train, x_test, response_name = "Species")
#' rue_test <- pt_oob(x_train, x_test, scorer = scorer_2)
#' str(rue_test)
#'
#' # Third example: sample memberships (class probabilities)
#' setosa <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' versicolor <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' scorer_3 <- function(x_train, x_test) score_cp(x_train, x_test)
#' cp_test <- pt_oob(setosa, versicolor, scorer = scorer_3)
#' str(cp_test)
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
#' be computionally expensive.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' setosa <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' versicolor <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' scorer <- function(x_train, x_test) score_od(x_train, x_test)
#' iris_test <- pt_refit(setosa, versicolor, scorer = scorer)
#' str(iris_test)
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
