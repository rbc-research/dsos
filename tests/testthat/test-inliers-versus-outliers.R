suppressPackageStartupMessages({
  library(dsos)
  library(testthat)
})

set.seed(123456)
n <- 1e4
alpha <- 0.05
test_that("Inliers", {
  skip_if_not_installed("ranger")

  x1 <- data.frame(x = rnorm(n, sd = 1))
  x2 <- data.frame(x = rnorm(n, sd = 1 / 2))
  expect_lt(at_oob(x1, x2, scorer = split_cp)$p_value, alpha)
})

test_that("Outliers", {
  skip_if_not_installed("ranger")

  x1 <- data.frame(x = rnorm(n, sd = 1))
  x2 <- data.frame(x = rnorm(n, sd = 3 / 2))
  expect_lt(at_oob(x1, x2, scorer = split_cp)$p_value, alpha)
})

test_that("Same distribution", {
  skip_if_not_installed("ranger")

  x1 <- data.frame(x = rnorm(n))
  x2 <- data.frame(x = rnorm(n))
  expect_gt(at_oob(x1, x2, scorer = split_cp)$p_value, alpha)
})
