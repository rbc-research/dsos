#' @noRd
#' @keywords  internal
asymptotic_sd <- function(n_train = 100, n_test = 100) {
  # Calculate each expression from analytical formula
  # See formula (9) in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4129959/
  # N.B: Use wolfram alpha for constants
  lambda <- n_test / n_train
  mean_null <- 1 / 12
  second <- mean_null^2
  fourth <- (1 / 4)^2
  first <- 1 / 63
  third <- 2 * (1 / 28)
  var_null <- lambda * (first - second) + (third - fourth)
  sd_null <- sqrt(var_null / n_test)
  return(sd_null)
}

#' @noRd
#' @keywords  internal
asymptotic_pvalue <- function(n_train = 100, n_test = 100, obs = 1 / 12) {

  # Get asymptotically normal null
  mean_null <- 1 / 12
  sd_null <- asymptotic_sd(n_train, n_test)

  # Get p-value from normal approximation
  pvalue <- 1 - stats::pnorm(obs, mean = mean_null, sd = sd_null)
  return(pvalue)
}

#' @noRd
#' @keywords  internal
asymptotic_os <- function(os_train, os_test) {

  # Calculate stats for two-sample test
  result_list <- wauc_and_os(os_train, os_test)

  # Compute p-value from null distribution
  p_value <- asymptotic_pvalue(
    n_train = length(os_train),
    n_test = length(os_test),
    obs = result_list$wauc
  )
  dsos_test <- list()
  dsos_test[["seq_mct"]] <- NULL
  dsos_test[["statistic"]] <- result_list$wauc
  dsos_test[["p_value"]] <- p_value
  dsos_test[["outlier_scores"]] <- result_list$outlier_scores
  class(dsos_test) <- "outlier.test"
  return(dsos_test)
}

#' @title
#' Asymptotic Test With Out-Of-Bag Scores
#'
#' @inherit pt_oob description return references
#' @inheritSection pt_oob Notes
#'
#' @param x_train Training (reference/validation) sample.
#' @param x_test Test sample.
#' @param scorer Function which returns a named list with outlier scores from
#' the training and test sample. The first argument to \code{scorer} must be
#' \code{x_train}; the second, \code{x_test}. The returned named list contains
#' two elements: \emph{train} and \emph{test}, each of which is a vector of
#' (outlier) scores. See notes for more information.
#'
#' @details
#' Li and Fine (2010) derives the asymptotic null distribution for the weighted
#' AUC (WAUC), the test statistic. This approach does not use permutations
#' and can, as a result, be much faster because it sidesteps the need to refit
#' the scoring function \code{scorer}. This works well for large samples. The
#' prefix \emph{at} stands for asymptotic test to tell it apart from the
#' prefix \emph{pt}, the permutation test.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' data(iris)
#' setosa <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' versicolor <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#'
#' # Using fake scoring function
#' scorer <- function(tr, te) list(train=runif(nrow(tr)), test=runif(nrow(te)))
#' oob_test <- at_oob(setosa, versicolor, scorer = scorer)
#' oob_test
#'
#' }
#'
#' @family asymptotic-test
#'
#' @seealso
#' [pt_oob()] for (faster) p-value approximation via out-of-bag predictions.
#' [pt_refit()] for (slower) p-value approximation via refitting.
#'
#' @export
at_oob <- function(x_train, x_test, scorer) {

  # Get list of outlier scores
  data.table::setDT(x_train)
  data.table::setDT(x_test)
  os_list <- scorer(x_train, x_test)

  # Gather test info
  result <- at_from_os(
    os_train = os_list[["train"]],
    os_test = os_list[["test"]]
  )
  return(result)
}

#' @title
#' Asymptotic Test from Outlier Scores
#'
#' @param os_train Outlier scores in training (reference) set.
#' @param os_test Outlier scores in test set.
#'
#' @inherit pt_oob description return references
#' @inherit at_oob details
#'
#' @section Notes:
#' The outlier scores should all mimic out-of-sample behaviour. Mind that the
#' training scores are not in-sample and thus, biased (overfitted) while the
#' test scores are out-of-sample. The mismatch -- in-sample versus out-of-sample
#' scores -- voids the test validity. A simple fix for this is to get the
#' training scores from an indepedent (fresh) validation set; this follows
#' the train/validation/test sample splitting convention and the validation set
#' is effectively the reference set or distribution in this case.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' os_train <- rnorm(n = 100)
#' os_test <- rnorm(n = 100)
#' test_result <- at_from_os(os_train, os_test)
#' test_result
#' }
#'
#' @family asymptotic-test
#'
#' @seealso
#' [at_oob()] for variant requiring a scoring function.
#' [pt_from_os()] for permutation test with the outlier scores.
#'
#' @export
at_from_os <- function(os_train, os_test) {
  result <- asymptotic_os(os_train, os_test)
  return(result)
}