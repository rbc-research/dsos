#' @noRd
#' @keywords  internal
smc_pvalue <- function(test_stat, permute_fn, R = 1e3) {

  # Generate null distribution via SMC
  gen_smc <- function() {
    gen_stat <- permute_fn()
    return(gen_stat >= test_stat)
  }
  res <- simctest::simctest(gen_smc, maxsteps = R)
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
                              R = 1e3) {
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
  test_list <- smc_pvalue(test_stat, permute_fn, R)
  test_list[["outlier_scores"]] <- helper$os_list
  return(test_list)
}

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
