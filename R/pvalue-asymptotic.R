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

  # Wrap result in a list
  test_list <- list()
  test_list[["seq_mct"]] <- NULL
  test_list[["statistic"]] <- result_list$wauc
  test_list[["p_value"]] <- p_value
  test_list[["outlier_scores"]] <- result_list$outlier_scores
  class(test_list) <- "outlier.test"
  return(test_list)
}

#' @noRd
#' @keywords  internal
asymptotic_null <- function(x_train, x_test, scorer) {

  # Get list of outlier scores
  data.table::setDT(x_train)
  data.table::setDT(x_test)
  os_list <- scorer(x_train, x_test)

  # Gather test info
  result <- asymptotic_os(
    os_train = os_list[["train"]],
    os_test = os_list[["test"]]
  )
  return(result)
}
