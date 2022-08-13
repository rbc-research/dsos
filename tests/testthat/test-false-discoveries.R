suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})

set.seed(123456)
null_proportion <- 0.9
n_obs <- 5e2
n_reps <- 2e3

test_that("High null proportion", {
  suppressWarnings(skip_on_cran())

  pvalues <- replicate(
    n = n_reps,
    expr = {
      os_train <- rnorm(n_obs)
      os_test <- rnorm(n_obs)
      test_at <- at_from_os(os_train, os_test)
      test_pt <- pt_from_os(os_train, os_test)
      c(test_at$p_value, test_pt$p_value)
    }
  )
  fdr_at <- fdrtool::fdrtool(
    pvalues[1, ],
    statistic = "pvalue",
    plot = FALSE,
    verbose = FALSE
  )
  expect_gt(fdr_at$param[, "eta0"], null_proportion)
  fdr_pt <- fdrtool::fdrtool(
    pvalues[2, ],
    statistic = "pvalue",
    plot = FALSE,
    verbose = FALSE
  )
  expect_gt(fdr_pt$param[, "eta0"], null_proportion)
})