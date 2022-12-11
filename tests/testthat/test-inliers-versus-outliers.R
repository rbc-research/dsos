suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})

set.seed(123456)
n <- 1e4
alpha <- 0.05
test_that("Inliers", {
  tr <- rnorm(n, sd = 1)
  te <- rnorm(n, sd = 1 / 2)
  expect_lt(at_from_os(tr, te)$p_value, alpha)
})

test_that("Outliers", {
  tr <- rnorm(n, sd = 1)
  te <- rnorm(n, sd = 3 / 2)
  expect_lt(at_from_os(tr, te)$p_value, alpha)
})

test_that("Same distribution", {
  tr <- rnorm(n)
  te <- rnorm(n)
  expect_gt(at_from_os(tr, te)$p_value, alpha)
})
