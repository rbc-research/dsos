suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})

set.seed(123456)
n <- 1e4
alpha <- 0.05
test_that("Inliers", {
  x1 <- data.frame(x = rnorm(n, sd = 1))
  x2 <- data.frame(x = rnorm(n, sd = 1 / 2))
  expect_lt(cp_ss(x1, x2)$p_value, alpha)
})

test_that("Outliers", {
  x1 <- data.frame(x = rnorm(n, sd = 1))
  x2 <- data.frame(x = rnorm(n, sd = 3 / 2))
  expect_lt(cp_ss(x1, x2)$p_value, alpha)
})

test_that("Same distribution", {
  x1 <- data.frame(x = rnorm(n))
  x2 <- data.frame(x = rnorm(n))
  expect_gt(cp_ss(x1, x2)$p_value, alpha)
})
