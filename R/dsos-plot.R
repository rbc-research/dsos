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
    label = scales::pvalue(p_value, add_p = T),
    hjust = 0,
    vjust = 1
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


#' Plot the result of the D-SOS test.
#'
#' @param x A \code{outlier.test} object from a D-SOS test.
#' @param ... Placeholder to be comptatible with S3 `plot` generic.
#'
#' @return A \pkg{ggplot2} plot with outlier scores and p-value.
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' data(iris)
#' x_train <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
#' x_test <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
#' iris_test <- pt_refit(x_train, x_test, scorer = score_od)
#' plot(iris_test)
#' }
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
