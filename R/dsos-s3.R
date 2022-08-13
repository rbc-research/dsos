#' @noRd
#' @keywords  internal
create_x_ticks <- function(os_range) {
  x_ticks <- round(seq(os_range[1], os_range[2], length.out = 4), 2)
  return(sort(x_ticks))
}

#' @noRd
#' @keywords  internal
plot_outliers <- function(data, os_range) {
  dens_mapping <- ggplot2::aes_string(
    x = "score",
    color = "source"
  )
  g <- ggplot2::ggplot(data = data, mapping = dens_mapping)
  g <- g + ggplot2::stat_ecdf(geom = "step", pad = TRUE, size = 1.2)
  g <- g + ggplot2::scale_y_continuous(
    name = "ECDF",
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
prep_os_data <- function(x) {
  stopifnot(inherits(x, what = "outlier.test"))

  # Extract outlier scores
  os <- x[["outlier_scores"]]
  os_range <- range(os)
  os_train <- os[["train"]]
  os_test <- os[["test"]]
  score <- c(os_train, os_test)

  # Mark sample (origign) source
  sc_train <- rep("Train", length(os_train))
  sc_test <- rep("Test", length(os_test))
  source <- as.factor(c(sc_train, sc_test))

  # Create data.frame for plotting
  data <- data.frame(score = score, source = source)
  return(list(data = data, os_range = os_range))
}


#' Plot result of test for no adverse shift.
#'
#' @param x A \code{outlier.test} object from a D-SOS test.
#' @param ... Placeholder to be comptatible with S3 `plot` generic.
#'
#' @return A \pkg{ggplot2} plot with outlier scores and p-value.
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' os_train <- rnorm(n=3e2)
#' os_test <- rnorm(n=3e2)
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


#' Print result of test for no adverse shift.
#'
#' @param x A \code{outlier.test} object from a D-SOS test.
#' @param n The number of outlier scores to print for each sample.
#' @param ... Placeholder to be comptatible with S3 `print` generic.
#'
#' @return A \pkg{ggplot2} plot with outlier scores and p-value.
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' os_train <- rnorm(n=3e2)
#' os_test <- rnorm(n=3e2)
#' test_to_print <- at_from_os(os_train, os_test)
#' # Also: pt_from_os(os_train, os_test) for permutation test
#' test_to_print
#' }
#'
#' @family s3-method
#'
#' @export
print.outlier.test <- function(x, n = 5, ...) {
  cat("\n")
  cat(strwrap("Test for no adverse shift", prefix = "\t"), sep = "\n")
  cat("\n")

  out <- character()
  if (!is.null(x[["p_value"]])) {
    fp <- format.pval(x[["p_value"]])
    fp <- if (startsWith(fp, "<")) fp else paste("=", fp)
    out <- c(out, paste("p-value", fp))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("\n")
  cat(
    paste0("The first ",
            n,
            " outlier scores (out of ",
            length(x[["outlier_scores"]][["train"]]),
            ") from the training (reference) sample:\n"
          )
    )
  print(utils::head(unname(x[["outlier_scores"]][["train"]]), n))
  cat("\n")
  cat(
    paste0("The first ",
            n,
            " outlier scores (out of ",
            length(x[["outlier_scores"]][["test"]]),
            ") from the test sample:\n"
          )
    )
  print(utils::head(unname(x[["outlier_scores"]][["test"]]), n))
  cat("\n")
  invisible(x)
}