suppressPackageStartupMessages({
  library(testthat)
  library(dsos)
})
set.seed(12345)

# Split samples by Species
data(iris)
setosa <- iris[1:50, 1:4] # Species == 'setosa'
versicolor <- iris[51:100, 1:4] # Species == 'versicolor'
virginica <- iris[101:150, 1:4] # Species == 'virginica'

# Random split
idx <- sample(nrow(iris), 2 / 3 * nrow(iris))
iris_train <- iris[idx, ]
iris_test <- iris[-idx, ]

# Stratified split
splits <- split(iris, list(iris$Species))
subsample_n <- function(x) ceiling(nrow(x) * 2 / 3)
indices <- lapply(splits, function(x) sample.int(nrow(x), subsample_n(x)))
include <- function(x, idx) x[idx, ]
exclude <- function(x, idx) x[-idx, ]
splits_train <- mapply(include, splits, indices, SIMPLIFY = FALSE)
splits_test <- mapply(exclude, splits, indices, SIMPLIFY = FALSE)
stratefied_train <- do.call(rbind, splits_train)
stratefied_test <- do.call(rbind, splits_test)

# Fix significance level and region of practical equivalence for s-values
alpha <- 0.1
s_rope <- 1.00

# Helper functions for s-values
abs_diff <- function(p1, p2, eps = 1e-3) {
  if (p1 <= eps) p1 <- eps
  if (p2 <= eps) p2 <- eps
  abs(-log(p1, 2) + log(p2, 2))
}

diff_from_samples <- function(x_train, x_test) {
  permutations <- cp_pt(x_train, x_test)
  asymptotic <- cp_at(x_train, x_test)
  abs_diff(permutations$p_value, asymptotic$p_value)
}

test_that("Separate species", {
  expect_lt(cp_at(setosa, virginica)$p_value, alpha)
  expect_lt(cp_at(setosa, versicolor)$p_value, alpha)
  expect_lt(cp_at(versicolor, virginica)$p_value, alpha)
  expect_lt(od_pt(setosa, virginica)$p_value, alpha)
  expect_lt(od_pt(setosa, versicolor)$p_value, alpha)
  expect_lt(od_pt(versicolor, virginica)$p_value, alpha)
})

test_that("Same species", {
  expect_gt(cp_at(setosa, setosa)$p_value, alpha)
  expect_gt(cp_at(versicolor, versicolor)$p_value, alpha)
  expect_gt(cp_at(virginica, virginica)$p_value, alpha)
  expect_gt(od_pt(setosa, setosa)$p_value, alpha)
  expect_gt(od_pt(versicolor, versicolor)$p_value, alpha)
  expect_gt(od_pt(virginica, virginica)$p_value, alpha)
})

test_that("Permutations versus asymptotic p-values (CP)", {
  # Iris: random splits
  expect_lte(diff_from_samples(iris_train, iris_test), s_rope)
  expect_lte(diff_from_samples(stratefied_train, stratefied_test), s_rope)

  # Iris: different species
  expect_lte(diff_from_samples(setosa, virginica), s_rope)
  expect_lte(diff_from_samples(setosa, versicolor), s_rope)
  expect_lte(diff_from_samples(versicolor, virginica), s_rope)
})

test_that("Random splits with class probabilities", {
  expect_lte(
    abs_diff(
      cp_pt(iris_train, iris_test)$p_value,
      cp_pt(stratefied_train, stratefied_test)$p_value
    ),
    s_rope
  )
})

test_that("Random splits with outlier scores", {
  expect_lte(
    abs_diff(
      od_pt(iris_train, iris_test)$p_value,
      od_pt(stratefied_train, stratefied_test)$p_value
    ),
    s_rope
  )
})

test_that("Random splits with confidence intervals", {
  expect_gte(
    abs_diff(
      rue_pt(iris_train, iris_test, response_name = "Species")$p_value,
      rue_pt(stratefied_train, stratefied_test, response_name = "Species")$p_value
    ),
    s_rope
  )
})


test_that("Random splits with out-of-bag residuals", {
  expect_gte(
    abs_diff(
      rd_pt(iris_train, iris_test, response_name = "Species")$p_value,
      rd_pt(stratefied_train, stratefied_test, response_name = "Species")$p_value
    ),
    s_rope
  )
})
