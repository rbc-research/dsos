suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})
set.seed(123456)

diff_pct <- 0.02
n_obs <- 1e4
n_rep <- 1e4
mean_asym <- 1 / 12

test_that("Equal sample size", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- future.apply::future_replicate(
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

test_that("Upsample test set by a factor of 2", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- future.apply::future_replicate(
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

test_that("Downsample test set by a factor of 10 (runif)", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- future.apply::future_replicate(
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

test_that("Downsample test set by a factor of 10 (rlnorm)", {
  suppressWarnings(testthat::skip_on_cran())
  wauc_null <- future.apply::future_replicate(
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
