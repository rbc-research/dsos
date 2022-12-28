suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})

set.seed(123456)
n <- 1e4
alpha <- 0.05
test_that("Inliers", {
  tr <- rlnorm(n, sdlog = 1)
  te <- rlnorm(n, sdlog = 1 / 2)
  expect_gt(at_from_os(tr, te)$p_value, alpha)
})

test_that("Outliers", {
  tr <- rlnorm(n, sdlog = 1)
  te <- rlnorm(n, sdlog = 3 / 2)
  expect_lt(at_from_os(tr, te)$p_value, alpha)
})

test_that("Same distribution", {
  tr <- rlnorm(n)
  te <- rlnorm(n)
  expect_gt(at_from_os(tr, te)$p_value, alpha)
})
