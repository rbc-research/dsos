suppressPackageStartupMessages({
  library(testthat)
  library(dsos)
})
set.seed(12345)

# Split samples by Species
data(iris)
setosa <- iris[1:50, 1:4]
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]

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
  permutations <- pt_oob(x_train, x_test, scorer = score_cp)
  asymptotic <- at_oob(x_train, x_test, scorer = score_cp)
  abs_diff(permutations$p_value, asymptotic$p_value)
}

test_that("Separate species (ranger)", {
  skip_if_not_installed("ranger")

  expect_lt(at_oob(setosa, virginica, scorer = score_cp)$p_value, alpha)
  expect_lt(at_oob(setosa, versicolor, scorer = score_cp)$p_value, alpha)
  expect_lt(at_oob(versicolor, virginica, scorer = score_cp)$p_value, alpha)
})

test_that("Separate species (isotree)", {
  skip_if_not_installed("isotree")

  expect_lt(pt_refit(setosa, virginica, scorer = score_od)$p_value, alpha)
  expect_lt(pt_refit(setosa, versicolor, scorer = score_od)$p_value, alpha)
  expect_lt(pt_refit(versicolor, virginica, scorer = score_od)$p_value, alpha)
})

test_that("Same species (ranger)", {
  skip_if_not_installed("ranger")

  expect_gt(at_oob(setosa, setosa, scorer = score_cp)$p_value, alpha)
  expect_gt(at_oob(versicolor, versicolor, scorer = score_cp)$p_value, alpha)
  expect_gt(at_oob(virginica, virginica, scorer = score_cp)$p_value, alpha)
})

test_that("Same species (isotree)", {
  skip_if_not_installed("isotree")

  expect_gt(pt_refit(setosa, setosa, scorer = score_od)$p_value, alpha)
  expect_gt(pt_refit(versicolor, versicolor, scorer = score_od)$p_value, alpha)
  expect_gt(pt_refit(virginica, virginica, scorer = score_od)$p_value, alpha)
})

test_that("Permutation versus asymptotic", {
  skip_if_not_installed("ranger")

  # Iris: random splits
  expect_lte(diff_from_samples(iris_train, iris_test), s_rope)
  expect_lte(diff_from_samples(stratefied_train, stratefied_test), s_rope)

  # Iris: different species
  expect_lte(diff_from_samples(setosa, virginica), s_rope)
  expect_lte(diff_from_samples(setosa, versicolor), s_rope)
  expect_lte(diff_from_samples(versicolor, virginica), s_rope)
})

test_that("Random splits with class probabilities", {
  skip_if_not_installed("ranger")

  expect_gte(
    pt_oob(iris_train, iris_test, scorer = score_cp)$p_value,
    alpha
  )
})

test_that("Stratefied (balanced) splits with class probabilities", {
  skip_if_not_installed("ranger")

  expect_gte(
    pt_oob(stratefied_train, stratefied_test, scorer = score_cp)$p_value,
    alpha
  )
})

test_that("Random splits with outlier scores", {
  skip_if_not_installed("isotree")

  expect_gte(
    pt_refit(iris_train, iris_test, scorer = score_od)$p_value,
    alpha
  )
})

test_that("Stratefied (balanced) splits with outlier scores", {
  skip_if_not_installed("isotree")

  expect_gte(
    pt_refit(stratefied_train, stratefied_test, scorer = score_od)$p_value,
    alpha
  )
})

test_that("Splits with prediction uncertainty", {
  skip_if_not_installed("ranger")

  scorer <- function(tr, te) score_rue(tr, te, response_name = "Species")
  expect_gte(
    abs_diff(
      pt_oob(iris_train, iris_test, scorer = scorer)$p_value,
      pt_oob(stratefied_train, stratefied_test, scorer = scorer)$p_value
    ),
    s_rope
  )
})


test_that("Random splits with out-of-bag residuals", {
  skip_if_not_installed("ranger")

  scorer <- function(tr, te) score_rd(tr, te, response_name = "Species")
  expect_gte(
    abs_diff(
      pt_oob(iris_train, iris_test, scorer = scorer)$p_value,
      pt_oob(stratefied_train, stratefied_test, scorer = scorer)$p_value
    ),
    s_rope
  )
})
