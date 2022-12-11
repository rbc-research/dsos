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

#' @title
#' Bayesian Test from Outlier Scores
#'
#' @param os_train Outlier scores in training (reference) set.
#' @param os_test Outlier scores in test set.
#' @param n_pt The number of permutations.
#' @param percentile_cutoff The number of permutations.
#' @param use_asymptotic Numeric vector of weights of length
#' @param asymptotic_cutoff Numeric vector of weights of length
#'
#' @return
#' A named list of class \code{outlier.bayes} containing:
#' \itemize{
#'    \item \code{permuted}: WAUC from permutations, if applicable
#'    \item \code{posterior}: WAUC from posterior distribution
#'    \item \code{adverse_threshold}: WAUC threshold for adverse shift
#'    \item \code{adverse_probability}: probability of adverse shift
#'    \item \code{bayes_factor}: Bayes factor
#'    \item \code{outlier_scores}: outlier scores from training and test set
#' }
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' os_train <- rnorm(n = 100)
#' os_test <- rnorm(n = 100)
#' bayes_pt <- bf_from_os(os_train, os_test)
#' bayes_pt
#' # Use the (faster) asymptotic cutoff
#' bayes_at <- bf_from_os(os_train, os_test, use_asymptotic = TRUE)
#' bayes_at
#' # Run in parallel on local cluster
#' library(future)
#' future::plan(future::multisession)
#' bayes_parallel <- bf_from_os(os_train, os_test)
#' bayes_parallel
#' }
#'
#' @family bayesian-test
#'
#' @export
bf_from_os <- function(os_train,
                       os_test,
                       n_pt = 4e3,
                       percentile_cutoff = 0.5,
                       use_asymptotic = FALSE,
                       asymptotic_cutoff = 1 / 12) {
  if (use_asymptotic) {
    permuted <- NULL
    posterior <- wauc_bb(os_train, os_test, n_pt = n_pt)
    adverse_threshold <- asymptotic_cutoff
  } else {
    draws <- wauc_samples(os_train, os_test, n_pt = n_pt)
    permuted <- draws$permuted
    posterior <- draws$posterior
    adverse_threshold <- stats::quantile(permuted, probs = percentile_cutoff)
  }
  adverse_prob <- 1 - stats::ecdf(posterior)(adverse_threshold)
  bayes_factor <- adverse_prob / (1 - adverse_prob)
  result <- list(
    permuted = permuted,
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

#' @title
#' Convert Bayesian to Permutation Test.
#'
#' @param bayes_test A \code{outlier.bayes} object from a Bayesian D-SOS test.
#'
#' @inherit pt_oob return
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' os_train <- rnorm(n = 3e2)
#' os_test <- rnorm(n = 3e2)
#' bayes_test <- bf_from_os(os_train, os_test, use_asymptotic = FALSE)
#' freq_test <- pt_from_bayes(bayes_test)
#' freq_test
#' }
#'
#' @family bayesian-test
#'
#' @export
pt_from_bayes <- function(bayes_test) {
  stopifnot(inherits(bayes_test, what = "outlier.bayes"))
  wauc_samples <- bayes_test[["permuted"]]
  if (is.null(wauc_samples)) {
    stop(
      paste0(
        "Bayesian test should use permutations. Set \"use_asymptotic\" ",
        "argument to \"bf_from_os\" to FALSE."
      ),
      call. = FALSE
    )
  }
  os <- bayes_test[["outlier_scores"]]
  test_stat <- wauc_from_os(
    os_train = os[["train"]],
    os_test = os[["test"]]
  )
  p_value <- 1 - stats::ecdf(wauc_samples)(test_stat)
  dsos_test <- list()
  dsos_test[["seq_mct"]] <- NULL
  dsos_test[["statistic"]] <- test_stat
  dsos_test[["p_value"]] <- p_value
  dsos_test[["outlier_scores"]] <- os
  class(dsos_test) <- "outlier.test"
  return(dsos_test)
}