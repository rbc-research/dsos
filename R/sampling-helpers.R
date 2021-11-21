#' @noRd
#' @keywords  internal
shuffle_os <- function(os_pooled, n_shuffle) {
  n <- length(os_pooled)
  idx_test <- sample(n, n_shuffle, replace = FALSE)
  os_test <- os_pooled[idx_test]
  os_train <- os_pooled[-idx_test]
  return(list(train = os_train, test = os_test))
}

#' @noRd
#' @keywords  internal
shuffle_data <- function(data, n_shuffle) {
  n <- nrow(data)
  idx_test <- sample(n, n_shuffle, replace = FALSE)
  permuted_test <- as.data.frame(data)[idx_test, , drop = FALSE]
  permuted_train <- as.data.frame(data)[-idx_test, , drop = FALSE]
  data.table::setDT(permuted_train)
  data.table::setDT(permuted_test)
  return(list(train = permuted_train, test = permuted_test))
}

#' @noRd
#' @keywords  internal
split_data <- function(data, sub_ratio = 1 / 2) {
  n <- nrow(data)
  n_split <- floor(sub_ratio * n)
  idx_split <- sample(n, n_split, replace = FALSE)
  x_fit <- as.data.frame(data)[idx_split, , drop = FALSE]
  x_train <- as.data.frame(data)[-idx_split, , drop = FALSE]
  data.table::setDT(x_fit)
  data.table::setDT(x_train)
  return(list(first = x_fit, second = x_train))
}

#' @noRd
#' @keywords  internal
permute_wauc_from_data <- function(data, n_shuffle, scorer) {
  shuffled <- shuffle_data(data, n_shuffle)
  wauc_stat <- wauc_from_data(
    x_train = shuffled[["train"]],
    x_test = shuffled[["test"]],
    scorer = scorer
  )
  return(wauc_stat)
}

#' @noRd
#' @keywords  internal
permute_wauc_from_os <- function(data, n_shuffle) {
  shuffled <- shuffle_os(data, n_shuffle)
  wauc_stat <- wauc_from_os(
    os_train = shuffled[["train"]],
    os_test = shuffled[["test"]]
  )
  return(wauc_stat)
}
