suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})

set.seed(123456)
null_proportion <- 0.75
n_obs <- 4e2
n_reps <- 2e3

test_that("High null proportion", {
  suppressWarnings(skip_on_cran())

  pvalues <- replicate(
    n = n_reps,
    expr = {
      x1 <- data.frame(x = rnorm(n_obs))
      x2 <- data.frame(x = rnorm(n_obs))
      cp_split <- cp_ss(x1, x2)
      c(cp_split$p_value)
    }
  )
  fdr_cp <- fdrtool::fdrtool(pvalues, statistic = "pvalue", plot = FALSE, verbose = FALSE)
  expect_gt(fdr_cp$param[, "eta0"], null_proportion)
  expect_gt(1.0, null_proportion)
})
