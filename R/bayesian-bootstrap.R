# Functions adapted (inspired) from bayesboot package

#' @noRd
#' @keywords  internal
rudirichlet <- function(n) {
  weights <- stats::rexp(n, rate = 1)
  weights <- weights / sum(weights)
  return(weights * n)
}

#' @noRd
#' @keywords  internal
iter_posterior <- function(os_train, os_test) {
  # Create instance weights
  n_train <- length(os_train)
  n_test <- length(os_test)
  w_train <- rudirichlet(n_train)
  w_test <- rudirichlet(n_test)
  weight <- c(w_train, w_test)
  # calculate WAUC with instance weights
  wauc_stat <- wauc_from_os(os_train, os_test, weight = weight)
  return(wauc_stat)
}

#' @noRd
#' @keywords  internal
iter_prior <- function(os_train, os_test) {
  # Permute scores (assumes exchangeable)
  n_test <- length(os_test)
  shuffled <- shuffle_os(c(os_train, os_test), n_test)
  wauc_stat <- iter_posterior(shuffled$train, shuffled$test)
  return(wauc_stat)
}

#' @noRd
#' @keywords  internal
repeat_fn <- function(os_train, os_test, fn, n_pt = 4e3) {
  wauc_dist <- replicate(n_pt, fn(os_train, os_test))
  return(wauc_dist)
}

#' @noRd
#' @keywords  internal
bb_posterior <- function(os_train, os_test, n_pt = 4e3) {
  wauc_dist <- repeat_fn(
    os_train,
    os_test,
    fn = iter_posterior,
    n_pt = n_pt
  )
  return(wauc_dist)
}

#' @noRd
#' @keywords  internal
bb_prior <- function(os_train, os_test, n_pt = 4e3) {
  wauc_dist <- repeat_fn(
    os_train,
    os_test,
    fn = iter_prior,
    n_pt = n_pt
  )
  return(wauc_dist)
}

#' @noRd
#' @keywords  internal
iter_both <- function(os_train, os_test) {
  # Create instance weights
  n_train <- length(os_train)
  n_test <- length(os_test)
  w_train <- rudirichlet(n_train)
  w_test <- rudirichlet(n_test)
  weight <- c(w_train, w_test)
  # calculate WAUC with instance weights
  posterior <- wauc_from_os(os_train, os_test, weight = weight)
  shuffled <- shuffle_os(c(os_train, os_test), n_test)
  prior <- wauc_from_os(shuffled$train, shuffled$test, weight = weight)
  return(c(prior = prior, posterior = posterior))
}

#' @noRd
#' @keywords  internal
bb_prior_and_posterior <- function(os_train, os_test, n_pt = 4e3) {
  wauc_dist <- repeat_fn(
    os_train,
    os_test,
    fn = iter_both,
    n_pt = n_pt
  )
  return(list(prior = wauc_dist[1, ], posterior = wauc_dist[2, ]))
}

#' @noRd
#' @keywords  internal
bf_empirical_null <- function(os_train,
                              os_test,
                              n_pt = 4e3,
                              alpha = 0.5) {
  sampler <- bb_prior_and_posterior(os_train, os_test, n_pt = n_pt)
  critical_value <- quantile(sampler$prior, probs = 1 - alpha, names = FALSE)
  tail_prob <- 1 - ecdf(sampler$posterior)(critical_value)
  bayes_factor <- tail_prob / (1 - tail_prob)
  result <- list(
    prior = sampler$prior,
    posterior = sampler$posterior,
    alpha = alpha,
    critical_value = critical_value,
    tail_prob = tail_prob,
    bayes_factor = bayes_factor
  )
  return(result)
}

#' @noRd
#' @keywords  internal
bf_asymptotic_null <- function(os_train,
                               os_test,
                               n_pt = 4e3,
                               critical_value = 1 / 12) {
  posterior <- bb_posterior(os_train, os_test, n_pt = n_pt)
  tail_prob <- 1 - ecdf(posterior)(critical_value)
  bayes_factor <- tail_prob / (1 - tail_prob)
  result <- list(
    posterior = posterior,
    critical_value = critical_value,
    bayes_factor = bayes_factor
  )
  return(result)
}

#' @noRd
#' @keywords  internal
bf_compare <- function(os_train,
                       os_test,
                       n_pt = 4e3,
                       alpha = 0.5,
                       cutoff_asymptotic = 1 / 12) {
  sampler <- bb_prior_and_posterior(os_train, os_test, n_pt = n_pt)
  cutoff_empirical <- quantile(sampler$prior, probs = 1 - alpha, names = FALSE)
  tail_empirical <- 1 - ecdf(sampler$posterior)(cutoff_empirical)
  tail_asymptotic <- 1 - ecdf(sampler$posterior)(cutoff_asymptotic)
  bf_empirical <- tail_empirical / (1 - tail_empirical)
  bf_asymptotic <- tail_asymptotic / (1 - tail_asymptotic)
  result <- list(
    prior = sampler$prior,
    posterior = sampler$posterior,
    alpha = alpha,
    cutoff_empirical = cutoff_empirical,
    cutoff_asymptotic = cutoff_asymptotic,
    bf_empirical = bf_empirical,
    bf_asymptotic = bf_asymptotic
  )
  return(result)
}
