#' Plot method for plausible correlation objects
#'
#' @param x An object of class `"plausible_cor"` (output from [run_plausible_cor()])
#' @param type Character string specifying the plot type. Options:
#'        - `"population"` (default): Shows posterior distribution of population-level correlation
#'        - `"sample"`: Shows distribution of sample-level correlations
#' @param ... Additional arguments passed to the underlying plotting function:
#'        [plot_population_cor()] for `type == "population"` or
#'        [plot_sample_cor()] for `type == "sample"`. See those functions for details.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#' @method plot plausible_cor
plot.plausible_cor <- function(
    x,
    type = c("population", "sample"),
    ...
) {

  type <- match.arg(type)
  result <- switch(
    type,
    "population" = {
      plot_population_cor(
        .data = x,
        ...
      )
    },
    "sample" = {
      plot_sample_cor(
        .data = x,
        ...
      )
    }
  )

  return(result)
}

#' Plot correlation posterior distributions from plausible value analysis
#'
#' @description
#' Visualises density traces of posterior distributions of Pearson correlation
#' coefficients obtained from [run_plausible_cor()], along with the mean
#' posterior density across all MCMC samples.
#'
#' @param .data A data frame returned by [run_plausible_cor()], containing at
#'        least `.draw`, `r`, and `posterior_updf` columns.
#' @param n_draws Integer specifying the maximum number of MCMC draws to plot.
#'        Defaults to `500` to avoid overplotting of density traces.
#' @param trace_aes A named list of aesthetics for the individual posterior
#'        density traces (e.g., `colour`, `alpha`, `linewidth`), passed to
#'        [ggplot2::geom_line()]. Defaults to
#'        `list(colour = "#4C4C4C", alpha = 0.1, linewidth = 0.1)`. Set to
#'        `FALSE` to omit these lines entirely.
#' @param mean_aes A named list of aesthetics for the mean posterior density
#'        trace (e.g., `colour`, `alpha`, `linewidth`), passed to
#'        [ggplot2::geom_line()]. Defaults to
#'        `list(colour = "#377EB8", linewidth = 1.5)`. Set to
#'        `FALSE` to omit the mean trace entirely.
#' @param mean_interval_args A named list that specifies how an interval for the
#'        mean posterior density is computed. Should include the elements `width`,
#'        which is a single numeric value of the desired interval width, and `method`,
#'        which is a character specifying which interval type to compute.
#'        Defaults to `list(width = 0.95, method = "hdci")`.
#' @param mean_interval_aes A named list of aesthetics for a light shaded
#'        background rectangle that represents the requested interval of the
#'        mean posterior density trace. Arguments are passed to [ggplot2::annotate()].
#'        Defaults to `list(fill = "#BDD7E7", alpha = 0.4)`.
#'        Set to `FALSE` to omit the interval entirely.
#' @param sample_rug_aes A named list of aesthetics for "rug" plot of the sample
#'        correlation coefficients on which the individual posterior density
#'        traces are based. Arguments are passed to [ggplot2::geom_rug()].
#'        Defaults to `list(alpha = 0.1, linewidth = 0.1, length = grid::unit(0.03, "npc"))`.
#'        Set to `FALSE` to omit the rug plot entirely.
#' @param zero_refline_aes A named list of aesthetics for the vertical reference
#'        line at zero (e.g., `linetype`, `linewidth`, `colour`), passed to
#'        [ggplot2::geom_vline()]. Defaults to
#'        `list(linetype = "dashed", linewidth = 0.75)`. Set to `FALSE` to omit
#'        the reference line entirely.
#' @param x_title Label for the x-axis. Defaults to `"population-level correlation"`
#' @param x_axis_limits Numeric vector of length 2 that defines x-axis limits
#'        to "zoom" into, using [ggplot2::coord_cartesian()]. Alternatively,
#'        `NULL` (default) retains the full scale.
#' @param plot_text_scaling Numeric scaling factor for axis text and title
#'        size. Default is `1`.
#' @param rng_seed Optional integer seed to ensure reproducible sampling of
#'        posterior draws.
#'
#' @details
#' This function plots the plausible posterior distributions of the Pearson
#' correlation coefficient obtained via [run_plausible_cor()]. It displays:
#'
#' - Density traces for a subset (`n_draws`) of individual posterior draws,
#'   capturing both parameter uncertainty and population variability.
#' - The mean posterior density across all draws, offering a summary of the
#'   plausible population-level correlation.
#' - An optional vertical reference line at 0 to facilitate visual
#'   interpretation.
#'
#' Subsampling with `n_draws` is used to limit visual clutter when many MCMC
#' draws are available. However, the mean posterior density is always computed
#' from the full set of draws.
#'
#' Each plotted layer can be customized using a named list of aesthetic values,
#' or disabled entirely by passing `FALSE`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @seealso [run_plausible_cor()]
#'
#' @export
plot_population_cor <- function(
    .data,
    n_draws = 500,
    trace_aes = NULL,
    mean_aes = NULL,
    mean_interval_args = NULL,
    mean_interval_aes = NULL,
    sample_rug_aes = NULL,
    zero_refline_aes = NULL,
    x_title = "population-level correlation",
    x_axis_limits = NULL,
    plot_text_scaling = 1,
    rng_seed = NULL
) {

  input_data <- .data

  validate_column_inputs(
    col_names = c(".draw", "r", "posterior_updf"),
    data_frame = input_data,
    data_name = ".data"
  )

  trace_data <- input_data %>%
    get_posterior_rho_densities()

  mean_data <- trace_data %>%
    get_mean_posterior_rho()

  plot_data <- input_data %>%
    dplyr::select(
      dplyr::all_of(".draw")
    ) %>%
    draw_rows(
      n_draws = n_draws,
      rng_seed = rng_seed
    ) %>%
    dplyr::left_join(
      y = trace_data,
      by = ".draw"
    )

  result <- ggplot2::ggplot(
    data = plot_data,
    mapping = ggplot2::aes(
      x = .data[["x"]]
    )
  )

  if (!isFALSE(mean_interval_aes)) {
    mean_interval_args <- mean_interval_args %||% list(
      width = 0.95, method = "hdci"
    )
    validate_geom_args(
      aes_list = mean_interval_args,
      allowed_names = c("width", "method"),
      arg_name = "mean_interval_args"
    )
    mean_interval <- get_interval(
      val = mean_data[["x"]],
      dens = mean_data[["mean_density"]],
      width = mean_interval_args[["width"]],
      method = mean_interval_args[["method"]]
    )

    mean_interval_aes <- mean_interval_aes %||% list(
      fill = "#BDD7E7", alpha = 0.4
    )
    validate_geom_args(
      aes_list = mean_interval_aes,
      disallowed_names = c(
        "geom", "x", "y", "xmin", "xmax", "ymin", "ymax", "xend", "yend"
      ),
      arg_name = "mean_interval_aes"
    )

    result <- result +
      rlang::exec(
        .fn = ggplot2::annotate,
        geom = "rect",
        xmin = min(mean_interval[["lower"]]),
        xmax = max(mean_interval[["upper"]]),
        ymin = -Inf,
        ymax = Inf,
        !!!mean_interval_aes
      )
  }

  if (!isFALSE(zero_refline_aes)) {
    zero_refline_aes <- zero_refline_aes %||% list(
      linetype = "dashed", linewidth = 0.75
    )
    validate_geom_args(
      aes_list = zero_refline_aes,
      disallowed_names = "xintercept",
      arg_name = "zero_refline_aes"
    )
    result <- result +
      rlang::exec(
        .fn = ggplot2::geom_vline,
        xintercept = 0,
        !!!zero_refline_aes
      )
  }

  if (!isFALSE(sample_rug_aes)) {
    sample_rug_aes <- sample_rug_aes %||% list(
      alpha = 0.1, linewidth = 0.1, length = grid::unit(0.03, "npc")
    )
    validate_geom_args(
      aes_list = sample_rug_aes,
      disallowed_names = c(
        "mapping", "data", "stat", "position", "outside", "sides"
      ),
      arg_name = "sample_rug_aes"
    )
    result <- result +
      rlang::exec(
        .fn = ggplot2::geom_rug,
        data = input_data,
        mapping = ggplot2::aes(
          x = .data[["r"]]
        ),
        sides = "b",
        !!!sample_rug_aes
      )
  }

  if (!isFALSE(trace_aes)) {
    trace_aes <- trace_aes %||% list(
      colour = "#4C4C4C", alpha = 0.1, linewidth = 0.1
    )
    validate_geom_args(
      aes_list = trace_aes,
      disallowed_names = c("mapping", "data", "stat", "position"),
      arg_name = "trace_aes"
    )
    result <- result +
      rlang::exec(
        .fn = ggplot2::geom_line,
        mapping = ggplot2::aes(
          y = .data[["density"]],
          group = .data[["r"]]
        ),
        !!!trace_aes
      )
  }

  if (!isFALSE(mean_aes)) {
    mean_aes <- mean_aes %||% list(
      colour = "#377EB8", linewidth = 1.5
    )
    validate_geom_args(
      aes_list = mean_aes,
      disallowed_names = c("mapping", "data", "stat", "position"),
      arg_name = "mean_aes"
    )
    result <- result +
      rlang::exec(
        .fn = ggplot2::geom_line,
        data = mean_data,
        mapping = ggplot2::aes(
          x = .data[["x"]],
          y = .data[["mean_density"]]
        ),
        !!!mean_aes
      )
  }

  result <- result +
    ggplot2::labs(
      x = x_title
    ) +
    ggplot2::coord_cartesian(
      xlim = x_axis_limits,
      expand = FALSE
    ) +
    add_theme_elements(plot_text_scaling)

  return(result)

}

#' Plot sample correlation coefficients from plausible value analysis
#'
#' @description
#' Visualises the distribution of sample Pearson correlation coefficients
#' obtained from [run_plausible_cor()], using either a dotplot or histogram
#' style. This function wraps [ggdist::stat_dotsinterval()] or
#' [ggdist::stat_histinterval()], providing a quick way to plot uncertainty
#' about the sample correlations.
#'
#' @param .data A data frame returned by [run_plausible_cor()], containing at
#'        least `.draw` and `r` columns.
#' @param n_draws Integer specifying the maximum number of MCMC draws to plot.
#'        Default is `Inf` (use all draws).
#' @param style Character, either `"dots"` (default) for a dotplot using
#'        [ggdist::stat_dotsinterval()] or `"hist"` for a histogram using
#'        [ggdist::stat_histinterval()].
#' @param plot_aes A named list of aesthetics, passed to the underlying plotting
#'        function. Defaults to
#'        `list(binwidth = grid::unit(c(1, Inf), "mm"), overflow = "compress", colour = "#377EB8", fill = "#7F7F7F", alpha = 0.75, point_alpha = 1, point_size = 3, interval_alpha = 1, interval_size_range = c(1.25, 2.5))`
#'        if `style == "dots"`, or to
#'        `list(breaks = ggdist::waiver(), colour = "#377EB8", fill = "#7F7F7F", alpha = 0.75, point_alpha = 1, point_size = 3, interval_alpha = 1, interval_size_range = c(1.25, 2.5))`
#'        if `style == "hist"`.
#' @param zero_refline_aes A named list of aesthetics for the vertical reference
#'        line at zero (e.g., `linetype`, `linewidth`, `colour`), passed to
#'        [ggplot2::geom_vline()]. Defaults to
#'        `list(linetype = "dashed", linewidth = 0.75)`. Set to `FALSE` to omit
#'        the reference line entirely.
#' @param point_interval_args A named list that specifies how the point +
#'        multiple-interval plot should be created. Should include the elements
#'        `point_method`, which is a string that specifies the point summary
#'        (e.g., `"mean"`), `interval_method`, which is a string that specifies
#'        the interval type (e.g., `"hdci"`), and `interval_width`, which is a
#'        numeric vector specifying the desired interval width(s). Defaults to
#'        `list(point_method = "mean", interval_method = "hdci", width = c(0.5, 0.8, 0.95))`.
#' @param x_title Label for the x-axis. Defaults to `"sample-level correlation"`.
#' @param x_axis_limits Numeric vector of length 2 that defines x-axis limits
#'        to "zoom" into, using [ggplot2::coord_cartesian()]. Alternatively,
#'        `NULL` (default) retains the full scale.
#' @param plot_text_scaling Numeric scaling factor for axis text and title
#'        size. Default is `1`.
#' @param rng_seed Optional integer seed for reproducible sampling of draws.
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
    plot_aes = NULL,
    zero_refline_aes = NULL,
    point_interval_args = NULL,
    x_title = "sample-level correlation",
    x_axis_limits = NULL,
    plot_text_scaling = 1,
    rng_seed = NULL
) {

  style <- rlang::arg_match(arg = style)

  validate_column_inputs(
    col_names = c(".draw", "r"),
    data_frame = .data,
    data_name = ".data"
  )

  plot_data <- .data %>%
    dplyr::select(
      dplyr::all_of(c(".draw", "r"))
    ) %>%
    draw_rows(
      n_draws = n_draws,
      rng_seed = rng_seed
    )

  result <- ggplot2::ggplot(
    data = plot_data,
    mapping = ggplot2::aes(
      x = .data[["r"]]
    )
  )

  if (!isFALSE(zero_refline_aes)) {
    zero_refline_aes <- zero_refline_aes %||% list(
      linetype = "dashed", linewidth = 0.75, colour = "black"
    )
    validate_geom_args(
      aes_list = zero_refline_aes,
      disallowed_names = "xintercept",
      arg_name = "zero_refline_aes"
    )
    result <- result +
      rlang::exec(
        .fn = ggplot2::geom_vline,
        xintercept = 0,
        !!!zero_refline_aes
      )
  }

  point_interval_args <- point_interval_args %||% list(
    point_method = "mean",
    interval_method = "hdci", interval_width = c(0.5, 0.8, 0.95)
  )
  point_interval_args <- validate_point_interval_args(point_interval_args)

  if (style == "dots") {
    plot_aes <- plot_aes %||% list(
      binwidth = grid::unit(c(1, Inf), "mm"), overflow = "compress",
      colour = "#377EB8", fill = "#7F7F7F", alpha = 0.75,
      point_alpha = 1, point_size = 3,
      interval_alpha = 1, interval_size_range = c(1.25, 2.5)
    )
    validate_geom_args(
      aes_list = plot_aes,
      disallowed_names = c(
        "mapping", "data", "geom", "position", "point_interval", ".width"
      ),
      arg_name = "plot_aes"
    )
    result <- result +
      rlang::exec(
        .fn = ggdist::stat_dotsinterval,
        point_interval = paste0(
          point_interval_args[["point_method"]], "_",
          point_interval_args[["interval_method"]]
        ),
        .width = point_interval_args[["interval_width"]],
        !!!plot_aes
      )
  } else {
    plot_aes <- plot_aes %||% list(
      breaks = ggdist::waiver(),
      colour = "#377EB8", fill = "#7F7F7F", alpha = 0.75,
      point_alpha = 1, point_size = 3,
      interval_alpha = 1, interval_size_range = c(1.25, 2.5)
    )
    validate_geom_args(
      aes_list = plot_aes,
      disallowed_names = c(
        "mapping", "data", "geom", "position", "point_interval", ".width"
      ),
      arg_name = "plot_aes"
    )
    result <- result +
      rlang::exec(
        .fn = ggdist::stat_histinterval,
        point_interval = paste0(
          point_interval_args[["point_method"]], "_",
          point_interval_args[["interval_method"]]
        ),
        .width = point_interval_args[["interval_width"]],
        !!!plot_aes
      )
  }

  result <- result +
    ggplot2::labs(
      x = x_title
    ) +
    ggplot2::coord_cartesian(
      xlim = x_axis_limits
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = 0
      )
    ) +
    add_theme_elements(plot_text_scaling)

  return(result)

}

#' @noRd
draw_rows <- function(.data, n_draws, rng_seed) {
  if (n_draws >= nrow(.data)) {
    return(.data)
  }
  if (!is.null(rng_seed) && is.finite(rng_seed)) {
    result <- withr::with_seed(
      seed = rng_seed,
      code = .data %>%
        dplyr::slice_sample(
          n = n_draws,
          replace = FALSE
        )
    )
  } else {
    result <- .data %>%
      dplyr::slice_sample(
        n = n_draws,
        replace = FALSE
      )
  }
  return(result)
}

#' @noRd
add_theme_elements <- function(plot_text_scaling) {
  result <- ggplot2::theme_minimal() +
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
        size = 30 * plot_text_scaling, vjust = -0.75
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
