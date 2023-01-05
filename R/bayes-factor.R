#' @noRd
#' @keywords  internal
rudirichlet <- function(n) {
  # Function adapted (inspired) from bayesboot package
  weights <- stats::rexp(n, rate = 1)
  weights <- weights / sum(weights)
  return(weights * n)
}

#' @noRd
#' @keywords  internal
draw_bb <- function(os_train, os_test) {
  # _bb suffix stands for bayesian bootstrap
  n_train <- length(os_train)
  n_test <- length(os_test)
  w_train <- rudirichlet(n_train)
  w_test <- rudirichlet(n_test)
  weight <- c(w_train, w_test)
  stat_bb <- wauc_from_os(os_train, os_test, weight = weight)
  return(stat_bb)
}

#' @noRd
#' @keywords  internal
repeat_fn <- function(os_train, os_test, fn, n_pt = 4e3) {
  wauc_dist <- future.apply::future_replicate(
    n_pt,
    fn(os_train, os_test)
  )
  return(wauc_dist)
}

#' @noRd
#' @keywords  internal
wauc_bb <- function(os_train, os_test, n_pt = 4e3) {
  wauc_dist <- repeat_fn(
    os_train,
    os_test,
    fn = draw_bb,
    n_pt = n_pt
  )
  return(wauc_dist)
}

#' @noRd
#' @keywords  internal
draw_bb_and_perm <- function(os_train, os_test) {
  n_train <- length(os_train)
  n_test <- length(os_test)
  w_train <- rudirichlet(n_train)
  w_test <- rudirichlet(n_test)
  weight <- c(w_train, w_test)
  wauc_bb <- wauc_from_os(os_train, os_test, weight = weight)
  shuffled <- shuffle_os(c(os_train, os_test), n_test)
  wauc_perm <- wauc_from_os(shuffled$train, shuffled$test, weight = NULL)
  return(c(permuted = wauc_perm, posterior = wauc_bb))
}

#' @noRd
#' @keywords  internal
wauc_samples <- function(os_train, os_test, n_pt = 4e3) {
  wauc_dist <- repeat_fn(
    os_train,
    os_test,
    fn = draw_bb_and_perm,
    n_pt = n_pt
  )
  return(list(permuted = wauc_dist[1, ], posterior = wauc_dist[2, ]))
}

#' @noRd
#' @keywords  internal
bf_and_pvalue <- function(os_train,
                          os_test,
                          adverse_threshold = 1 / 12,
                          n_pt = 4e3) {
  draws <- wauc_samples(os_train, os_test, n_pt = n_pt)
  # Get p-value
  permuted <- draws$permuted
  test_stat <- wauc_from_os(os_train, os_test)
  p_value <- 1 - stats::ecdf(permuted)(test_stat)
  # Get bayes factor from asymptotic threshold
  posterior <- draws$posterior
  cdf_fn <- stats::ecdf(posterior)
  asym_prob <- 1 - cdf_fn(adverse_threshold)
  asym_bf <- asym_prob / (1 - asym_prob)
  # Get bayes factor from exchangeable threshold
  pt_threshold <- stats::quantile(permuted, probs = 0.5)
  pt_prob <- 1 - cdf_fn(pt_threshold)
  pt_bf <- pt_prob / (1 - pt_prob)
  bayes_factor <- list(exchangeable = pt_bf, asymptotic = asym_bf)
  return(list(bayes_factor = bayes_factor, p_value = p_value))
}

#' @title
#' Bayesian Test from Outlier Scores
#'
#' @description
#' Test for no adverse shift with outlier scores. Like goodness-of-fit testing,
#' this two-sample comparison takes the training (outlier) scores,
#' \code{os_train}, as the reference. The method checks whether the test
#' scores, \code{os_test}, are worse off relative to the training set.
#'
#' @param os_train Outlier scores in training (reference) set.
#' @param os_test Outlier scores in test set.
#' @param n_pt The number of permutations.
#' @param adverse_threshold WAUC threshold for adverse shift.
#'
#' @return
#' A named list of class \code{outlier.bayes} containing:
#' \itemize{
#'    \item \code{posterior}: Posterior distribution of WAUC test statistic
#'    \item \code{adverse_threshold}: WAUC threshold for adverse shift
#'    \item \code{adverse_probability}: probability of adverse shift
#'    \item \code{bayes_factor}: Bayes factor
#'    \item \code{outlier_scores}: outlier scores from training and test set
#' }
#'
#' @references Kamulete, V. M. (2023).
#' \emph{Are you OK? A Bayesian sequential test for adverse shift}.
#' Manuscript in preparation.
#'
#' @references Johnson, V. E. (2005).
#' \emph{Bayes factors based on test statistics}.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 67(5), 689-701.
#'
#' @references Gu, J., Ghosal, S., & Roy, A. (2008).
#' \emph{Bayesian bootstrap estimation of ROC curve}.
#' Statistics in medicine, 27(26), 5407-5420.
#'
#' @details
#' The posterior distribution of the WAUC test statistic is based on \code{n_pt}
#' (boostrap) permutations. The method relies on the Bayesian bootstrap as
#' a resampling procedure to smooth the instance weights as in Gu et al (2008).
#' This framework borrows from Johnson (2005), which shows how to turn a test
#' statistic, which may very well be frequentist, into its Bayesian counterpart,
#' the bayes factor.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' os_train <- rnorm(n = 100)
#' os_test <- rnorm(n = 100)
#' bayes_test <- bf_from_os(os_train, os_test)
#' bayes_test
#' # Run in parallel on local cluster
#' library(future)
#' future::plan(future::multisession)
#' parallel_test <- bf_from_os(os_train, os_test)
#' parallel_test
#' }
#'
#' @family bayesian-test
#'
#' @export
bf_from_os <- function(os_train,
                       os_test,
                       n_pt = 4e3,
                       adverse_threshold = 1 / 12) {
  posterior <- wauc_bb(os_train, os_test, n_pt = n_pt)
  adverse_prob <- 1 - stats::ecdf(posterior)(adverse_threshold)
  bayes_factor <- adverse_prob / (1 - adverse_prob)
  result <- list(
    posterior = posterior,
    adverse_threshold = adverse_threshold,
    adverse_probability = adverse_prob,
    bayes_factor = bayes_factor,
    outlier_scores = list(train = os_train, test = os_test)
  )
  class(result) <- "outlier.bayes"
  return(result)
}

#' @title
#' Convert P-value to Bayes Factor
#'
#' @param pvalue P-value.
#'
#' @return Bayes Factor (scalar value).
#'
#' @references Marsman, M., & Wagenmakers, E. J. (2017).
#' \emph{Three insights from a Bayesian interpretation of the one-sided P value}.
#' Educational and Psychological Measurement, 77(3), 529-539.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' bf_from_pvalue <- as_bf(pvalue = 0.5)
#' bf_from_pvalue
#' }
#'
#' @family bayesian-test
#'
#' @seealso
#' [as_pvalue()] to convert Bayes factor to p-value.
#'
#' @export
as_bf <- function(pvalue) {
  bf <- exp(stats::qlogis(pvalue))
  inv_bf <- 1. / bf
  return(inv_bf)
}

#' @title
#' Convert Bayes Factor to P-value
#'
#' @param bf Bayes factor.
#'
#' @return p-value (scalar value).
#'
#' @references Marsman, M., & Wagenmakers, E. J. (2017).
#' \emph{Three insights from a Bayesian interpretation of the one-sided P value}.
#' Educational and Psychological Measurement, 77(3), 529-539.
#' 
#' @examples
#' \donttest{
#' library(dsos)
#' pvalue_from_bf <- as_pvalue(bf = 1)
#' pvalue_from_bf
#' }
#'
#' @family bayesian-test
#'
#' @seealso
#' [as_bf()] to convert p-value to Bayes factor.
#'
#' @export
as_pvalue <- function(bf) {
  inv_bf <- 1. / bf
  pvalue <- stats::plogis(log(inv_bf))
  return(pvalue)
}
