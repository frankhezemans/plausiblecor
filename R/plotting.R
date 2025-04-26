#' Plot sample correlation coefficients from plausible value analysis
#'
#' @description
#' Visualises the distribution of sample Pearson correlation coefficients
#' obtained from `run_plausible_cor()`, using either a dotplot or histogram
#' style. This function wraps [ggdist::stat_dotsinterval()] or
#' [ggdist::stat_histinterval()], providing a quick way to plot uncertainty
#' about the sample correlations.
#'
#' @param .data A data frame returned by [run_plausible_cor()], containing at
#'        least `.draw` and `r` columns.
#' @param n_draws Integer specifying the maximum number of MCMC draws to plot.
#'        Default is `Inf` (use all draws).
#' @param style Character, either `"dots"` (default) for a dotplot or `"hist"`
#'        for a histogram.
#' @param binwidth,breaks Passed to [ggdist::stat_dotsinterval()] or
#'        [ggdist::stat_histinterval()]. Control bin width or histogram breaks.
#' @param colour,fill Colours for the plotted intervals/dots. Default greys.
#' @param alpha Alpha transparency level for intervals/dots. Default is `0.75`.
#' @param zero_refline Logical. If `TRUE` (default), adds a vertical dashed line
#'        at 0 to help reference null correlations.
#' @param point_method Character, `"mean"` (default) or `"median"`, indicating
#'        how to summarize the sample draws.
#' @param interval_width Numeric vector giving widths of uncertainty intervals
#'        to display. Default is `c(0.5, 0.8, 0.95)`.
#' @param interval_method Character, `"hdci"` (default) for highest-density
#'        continuous intervals or `"qi"` for quantile intervals.
#' @param x_title Label for the x-axis. Defaults to waiver (no label).
#' @param x_axis_limits Optional numeric vector of length 2 to control x-axis
#'        limits using [ggplot2::coord_cartesian()].
#' @param plot_text_scaling Numeric scaling factor for axis text and title
#'        size. Default is `1`.
#' @param rng_seed Optional integer seed for reproducible sampling of draws.
#' @param ... Additional arguments passed to [ggdist::stat_dotsinterval()] or
#'        [ggdist::stat_histinterval()].
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @details
#' Designed for visualizing the output of [run_plausible_cor()] — a plausible
#' value analysis for correlations between model parameters and covariates.
#'
#' If `n_draws` is less than the number of available MCMC samples, a random
#' subset is selected (optionally with reproducible seeding via `rng_seed`).
#'
#' @seealso [run_plausible_cor()], [ggdist::stat_dotsinterval()], [ggdist::stat_histinterval()]
#'
#' @export
plot_sample_cor <- function(
    .data,
    n_draws = Inf,
    style = c("dots", "hist"),
    binwidth = grid::unit(c(1, Inf), "mm"),
    breaks = ggdist::waiver(),
    colour = "#4C4C4C",
    fill = "#7F7F7F",
    alpha = 0.75,
    zero_refline = TRUE,
    point_method = c("mean", "median"),
    interval_width = c(0.5, 0.8, 0.95),
    interval_method = c("hdci", "qi"),
    x_title = ggplot2::waiver(),
    x_axis_limits = NULL,
    plot_text_scaling = 1,
    rng_seed = NULL,
    ...
) {

  style <- rlang::arg_match(arg = style)
  point_method <- rlang::arg_match(arg = point_method)
  interval_method <- rlang::arg_match(arg = interval_method)

  validate_column_inputs(
    col_names = c(".draw", "r"),
    data_frame = .data,
    data_name = ".data"
  )

  draw_rows <- function(data, n = n_draws) {
    dplyr::slice_sample(
      .data = data,
      n = n,
      replace = FALSE
    )
  }

  plot_data <- .data %>%
    dplyr::distinct(
      dplyr::all_of(c(".draw", "r"))
    )

  if (n_draws < nrow(plot_data)) {
    if (!is.null(rng_seed)) {
      plot_data <- withr::with_seed(
        seed = rng_seed,
        code = draw_rows(data = plot_data)
      )
    } else {
      plot_data <- draw_rows(data = plot_data)
    }
  }

  base_plot <- ggplot2::ggplot(
    data = plot_data,
    mapping = ggplot2::aes(
      x = .data[["r"]]
    )
  )

  if (zero_refline) {
    base_plot <- base_plot +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 1.25,
        colour = "black"
      )
  }

  if (style == "dots") {

    result <- base_plot +
      ggdist::stat_dotsinterval(
        point_interval = paste0(point_method, "_", interval_method),
        .width = interval_width,
        interval_size_range = c(1, 3),
        shape = 21,
        stroke = 1.75,
        point_size = 5,
        point_fill = "white",
        point_colour = "black",
        point_alpha = 1,
        interval_colour = "black",
        interval_alpha = 1,
        binwidth = binwidth,
        overflow = "compress",
        alpha = alpha,
        colour = colour,
        fill = fill,
        ...
      )

  } else {

    result <- base_plot +
      ggdist::stat_histinterval(
        point_interval = paste0(point_method, "_", interval_method),
        .width = interval_width,
        interval_size_range = c(1, 3),
        shape = 21,
        stroke = 1.75,
        point_size = 5,
        point_fill = "white",
        point_colour = "black",
        point_alpha = 1,
        interval_colour = "black",
        interval_alpha = 1,
        breaks = breaks,
        alpha = alpha,
        colour = colour,
        fill = fill,
        ...
      )

  }

  result <- result +
    ggplot2::labs(
      x = x_title
    )

  result <- add_theme_elements(
    x = result,
    plot_text_scaling = plot_text_scaling
  )

  if (!is.null(x_axis_limits)) {
    result <- result +
      ggplot2::coord_cartesian(xlim = x_axis_limits)
  }

  return(result)

}

#' @noRd
add_theme_elements <- function(x, plot_text_scaling) {
  result <- x +
    ggplot2::scale_x_continuous(
      # https://doi.org/10.1016/j.paid.2016.06.069
      breaks = c(-1, -0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5, 1),
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "#CCCCCC", linewidth = 0.25
      ),
      legend.position = "null",
      axis.line.x = ggplot2::element_line(colour = "#CCCCCC", linewidth = 1),
      axis.ticks.x = ggplot2::element_line(colour = "#CCCCCC", linewidth = 1),
      axis.title.x = ggplot2::element_text(
        size = 30 * plot_text_scaling, face = "bold.italic", vjust = -0.75
      ),
      axis.text.x = ggplot2::element_text(
        size = 18 * plot_text_scaling, colour = "black"
      ),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
    )
  return(result)
}
