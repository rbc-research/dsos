suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})
set.seed(123456)

diff_pct <- 0.02
n_obs <- 1e3
n_rep <- 1e4
mean_asym <- 1 / 12

test_that("n_test = n_train = 1e3", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- replicate(
    n = n_rep,
    expr = {
      os_train <- runif(n = n_obs)
      os_test <- runif(n = n_obs)
      wauc_from_os(os_train, os_test)
    }
  )

  mean_emp <- mean(wauc_null)
  sd_emp <- sd(wauc_null)
  sd_asymp <- dsos:::asymptotic_sd(n_obs, n_obs)
  expect_equal(mean_emp, mean_asym, tolerance = diff_pct * mean_emp)
  expect_equal(sd_emp, sd_asymp, tolerance = diff_pct * sd_emp)
})

test_that("n_test = 2 * n_train; n_train = 1e3", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- replicate(
    n = n_rep,
    expr = {
      os_train <- runif(n = n_obs)
      os_test <- runif(n = 2 * n_obs)
      wauc_from_os(os_train, os_test)
    }
  )

  mean_emp <- mean(wauc_null)
  sd_emp <- sd(wauc_null)
  sd_asymp <- dsos:::asymptotic_sd(n_obs, 2 * n_obs)
  expect_equal(mean_emp, mean_asym, tolerance = diff_pct * mean_emp)
  expect_equal(sd_emp, sd_asymp, tolerance = diff_pct * sd_emp)
})

test_that("n_test = 0.1 * n_train; n_train = 1e3 (runif)", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- replicate(
    n = n_rep,
    expr = {
      os_train <- runif(n = n_obs)
      os_test <- runif(n = 0.1 * n_obs)
      wauc_from_os(os_train, os_test)
    }
  )

  mean_emp <- mean(wauc_null)
  sd_emp <- sd(wauc_null)
  sd_asymp <- dsos:::asymptotic_sd(n_obs, 0.1 * n_obs)
  expect_equal(mean_emp, mean_asym, tolerance = diff_pct * mean_emp)
  expect_equal(sd_emp, sd_asymp, tolerance = diff_pct * sd_emp)
})

test_that("n_test = 0.1 * n_train; n_train = 1e3 (rlnorm)", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- replicate(
    n = n_rep,
    expr = {
      os_train <- rlnorm(n = n_obs)
      os_test <- rlnorm(n = 0.1 * n_obs)
      wauc_from_os(os_train, os_test)
    }
  )

  mean_emp <- mean(wauc_null)
  sd_emp <- sd(wauc_null)
  sd_asymp <- dsos:::asymptotic_sd(n_obs, 0.1 * n_obs)
  expect_equal(mean_emp, mean_asym, tolerance = diff_pct * mean_emp)
  expect_equal(sd_emp, sd_asymp, tolerance = diff_pct * sd_emp)
})
