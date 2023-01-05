#' @noRd
#' @keywords  internal
create_x_ticks <- function(os_range) {
  x_ticks <- round(seq(os_range[1], os_range[2], length.out = 4), 2)
  return(sort(x_ticks))
}

#' @noRd
#' @keywords  internal
plot_outliers <- function(data, os_range) {
  mapping <- ggplot2::aes_string(
    x = "score",
    color = "source"
  )
  g <- ggplot2::ggplot(data = data, mapping = mapping)
  g <- g + ggplot2::stat_ecdf(geom = "step", pad = TRUE, size = 1.2)
  g <- g + ggplot2::scale_y_continuous(
    name = "Percentile",
    trans = "exp",
    labels = scales::percent
  )
  g <- g + ggplot2::scale_color_brewer(name = NULL, type = "qual")
  x_ticks <- create_x_ticks(os_range)
  g <- g + ggplot2::scale_x_continuous(
    name = "Outlier scores",
    breaks = x_ticks
  )

  # Fix legend
  os_legend <- set_legend_manual()
  g <- g + os_legend[["theme"]]
  return(g)
}

#' @noRd
#' @keywords  internal
set_legend_manual <- function() {
  dsos_theme <- ggplot2::theme_light() + ggplot2::theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.key = ggplot2::element_blank()
  )
  return(list(theme = dsos_theme))
}

#' @noRd
#' @keywords  internal
annotate_pvalue <- function(p_value) {
  pvalue_layer <- ggplot2::annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = scales::pvalue(p_value, add_p = TRUE),
    hjust = 0,
    vjust = 1,
    size = 9
  )
  return(pvalue_layer)
}

#' @noRd
#' @keywords  internal
annotate_bf <- function(bf) {
  bf_layer <- ggplot2::annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = paste0("BF = ", round(bf, 2)),
    hjust = 0,
    vjust = 1,
    size = 9
  )
  return(bf_layer)
}

#' @noRd
#' @keywords  internal
prep_os_data <- function(x) {
  stopifnot(inherits(x, what = c("outlier.test", "outlier.bayes")))
  os <- x[["outlier_scores"]]
  os_range <- range(os)
  os_train <- os[["train"]]
  os_test <- os[["test"]]
  score <- c(os_train, os_test)
  sc_train <- rep("Train", length(os_train))
  sc_test <- rep("Test", length(os_test))
  source <- as.factor(c(sc_train, sc_test))
  data <- data.frame(score = score, source = source)
  return(list(data = data, os_range = os_range))
}


#' Plot frequentist test for no adverse shift.
#'
#' @param x A \code{outlier.test} result from test of no adverse shift.
#' @param ... Placeholder to be compatible with S3 method \code{plot}.
#'
#' @return A \pkg{ggplot2} plot with outlier scores and p-value.
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' os_train <- rnorm(n = 3e2)
#' os_test <- rnorm(n = 3e2)
#' test_to_plot <- at_from_os(os_train, os_test)
#' # Also: pt_from_os(os_train, os_test) for permutation test
#' plot(test_to_plot)
#' }
#'
#' @family s3-method
#'
#' @export
plot.outlier.test <- function(x, ...) {

  # Prep data from D-SOS test
  os_list <- prep_os_data(x)
  os_range <- os_list[["os_range"]]
  os_df <- os_list[["data"]]
  p_value <- x[["p_value"]]

  # Plot outliers
  os_plot <- plot_outliers(os_df, os_range)

  # Add p-value
  os_plot <- os_plot + annotate_pvalue(p_value)
  return(os_plot)
}

#' Plot Bayesian test for no adverse shift.
#'
#' @param x A \code{outlier.bayes} result from test of no adverse shift.
#' @param ... Placeholder to be compatible with S3 method \code{plot}.
#'
#' @return A \pkg{ggplot2} plot with outlier scores and p-value.
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' os_train <- rnorm(n = 3e2)
#' os_test <- rnorm(n = 3e2)
#' test_to_plot <- bf_from_os(os_train, os_test)
#' plot(test_to_plot)
#' }
#'
#' @family s3-method
#'
#' @export
plot.outlier.bayes <- function(x, ...) {

  # Prep data from D-SOS test
  os_list <- prep_os_data(x)
  os_range <- os_list[["os_range"]]
  os_df <- os_list[["data"]]
  bf <- x[["bayes_factor"]]

  # Plot outliers
  os_plot <- plot_outliers(os_df, os_range)

  # Add p-value
  os_plot <- os_plot + annotate_bf(bf)
  return(os_plot)
}

#' Print frequentist test for no adverse shift.
#'
#' @param x A \code{outlier.test} object from a D-SOS test.
#' @param ... Placeholder to be compatible with S3 method \code{plot}.
#'
#' @return Print to screen: display p-value and other information.
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' os_train <- rnorm(n = 3e2)
#' os_test <- rnorm(n = 3e2)
#' test_to_print <- at_from_os(os_train, os_test)
#' # Also: pt_from_os(os_train, os_test) for permutation test
#' test_to_print
#' }
#'
#' @family s3-method
#'
#' @export
print.outlier.test <- function(x, ...) {
  cat(strwrap("Frequentist test for no adverse shift", prefix = "\t"), "\n")
  cat("\n")
  out <- character()
  if (!is.null(x[["p_value"]])) {
    p_str <- format.pval(x[["p_value"]])
    p_fmt <- if (startsWith(p_str, "<")) p_str else paste("=", p_str)
    out <- c(out, paste("p-value", p_fmt))
  }
  if (!is.null(x[["statistic"]])) {
    wauc_str <- round(x[["statistic"]], 4)
    wauc_fmt <- paste("=", wauc_str)
    out <- c(out, paste("test statistic (weighted AUC/WAUC)", wauc_fmt))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("\n")
  alt_hypothesis <- paste0(
    "Alternative hypothesis: Pr(WAUC >= ", wauc_str, ")"
  )
  cat(strwrap(alt_hypothesis), sep = "\n")
  cat(strwrap("=> the test set is worse off than training."), sep = "\n")
  n_train <- length(x[["outlier_scores"]][["train"]])
  n_test <- length(x[["outlier_scores"]][["test"]])
  samples <- paste0(
    "Sample sizes: ",
    n_train,
    " in training and ",
    n_test,
    " in test set."
  )
  cat(strwrap(paste(samples, collapse = ", ")), sep = "\n")
  invisible(x)
}

#' Print Bayesian test for no adverse shift.
#'
#' @param x A \code{outlier.test} object from a D-SOS test.
#' @param ... Placeholder to be compatible with S3 method \code{plot}.
#'
#' @return Print to screen: display Bayes factor and other information.
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' os_train <- rnorm(n = 3e2)
#' os_test <- rnorm(n = 3e2)
#' test_to_print <- bf_from_os(os_train, os_test)
#' test_to_print
#' }
#'
#' @family s3-method
#'
#' @export
print.outlier.bayes <- function(x, ...) {
  cat(strwrap("Bayesian test for no adverse shift", prefix = "\t"), "\n")
  cat("\n")
  out <- character()
  if (!is.null(x[["bayes_factor"]])) {
    bf_str <- round(x[["bayes_factor"]], 2)
    bf_fmt <- paste("=", bf_str)
    out <- c(out, paste("Bayes factor (BF)", bf_fmt))
  }
  if (!is.null(x[["adverse_threshold"]])) {
    wauc_str <- round(x[["adverse_threshold"]], 4)
    wauc_fmt <- paste("=", wauc_str)
    out <- c(out, paste("cutoff (weighted AUC/WAUC)", wauc_fmt))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("\n")
  model_str <- sprintf(
    "Model: bayesian bootstrap with %d replicates (simulations)",
    length(x[["posterior"]])
  )
  cat(strwrap(model_str), "\n")
  nmtr_str <- paste0("BF's numerator: Pr(WAUC >= ", wauc_str, ")")
  dnmt_str <- paste0("BF's denominator: Pr(WAUC < ", wauc_str, ")")
  cat(strwrap(nmtr_str), "\n")
  cat(strwrap(dnmt_str), "\n")
  cat(
    strwrap(
      paste0(
        "=> BF > 3 favors view that ",
        "the test set is worse off than training.")
        ),
    sep = "\n"
  )
  n_train <- length(x[["outlier_scores"]][["train"]])
  n_test <- length(x[["outlier_scores"]][["test"]])
  samples <- paste0(
    "Sample sizes: ",
    n_train,
    " in training and ",
    n_test,
    " in test set."
  )
  cat(strwrap(paste(samples, collapse = ", ")), sep = "\n")
  invisible(x)
}
