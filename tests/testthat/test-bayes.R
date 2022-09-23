suppressPackageStartupMessages({
  library(testthat)
  library(dsos)
})
set.seed(12345)

check_conversion <- function(pvalue) {
    bf_from_pvalue <- as_bf(pvalue)
    pvalue_from_bf <- as_pvalue(bf_from_pvalue)
    is_same <- all.equal(pvalue_from_bf, pvalue)
    return(is_same)
}

test_that("convert p-value to Bayes factor", {
  suppressWarnings(testthat::skip_on_cran())
  n_val <- 100
  pvalues <- seq(0, 1, length.out = n_val)[c(-1, -n_val)]
  expect_true(check_conversion(pvalues))
})