% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_population_cor}
\alias{plot_population_cor}
\title{Plot correlation posterior distributions from plausible value analysis}
\usage{
plot_population_cor(
  .data,
  n_draws = 500,
  trace_colour = "#4C4C4C",
  trace_alpha = 0.1,
  trace_linewidth = 0.1,
  mean_colour = "#377EB8",
  mean_alpha = 1,
  mean_linewidth = 1.25,
  zero_refline = TRUE,
  zero_refline_linetype = "dashed",
  zero_refline_linewidth = 1.25,
  zero_refline_colour = "black",
  x_title = ggplot2::waiver(),
  x_axis_limits = NULL,
  plot_text_scaling = 1,
  rng_seed = NULL
)
}
\arguments{
\item{.data}{A data frame returned by [run_plausible_cor()], containing
`.draw`, `r`, and `posterior_updf` columns.}

\item{n_draws}{Integer specifying the maximum number of MCMC draws to plot.
Default is `500` to avoid overplotting of density traces.}

\item{trace_colour}{Colour for individual trace lines. Defaults to dark gray
(`"#4C4C4C"`).}

\item{trace_alpha}{Alpha (transparency) for individual trace lines. Defaults
to `0.1`.}

\item{trace_linewidth}{Line width for individual trace lines. Defaults to
`0.1`.}

\item{mean_colour}{Colour for the mean density trace. Defaults to blue
(`"#377EB8"`).}

\item{mean_alpha}{Alpha (transparency) for the mean density trace. Defaults
to `1`.}

\item{mean_linewidth}{Line width for the mean density trace. Defaults to `3`.}

\item{zero_refline}{Logical. If `TRUE` (default), adds a vertical dashed line
at 0 to help reference null correlations.}

\item{zero_refline_linetype, zero_refline_linewidth, zero_refline_colour}{Characters used to control the appearance of the reference line at
zero (if applicable).}

\item{x_title}{Label for the x-axis. Defaults to [`ggplot2::waiver()`waiver]
(no label).}

\item{x_axis_limits}{Optional numeric vector of length 2 to control x-axis
limits using [ggplot2::coord_cartesian()].}

\item{plot_text_scaling}{Numeric scaling factor for axis text and title
size. Default is `1`.}

\item{rng_seed}{Optional integer seed for reproducible sampling of draws.}
}
\value{
A [ggplot2::ggplot] object.
}
\description{
Visualises density traces of posterior distributions of Pearson correlation
coefficients obtained from [`run_plausible_cor()`], along with the mean
posterior density across all MCMC samples.
}
\details{
This function plots the plausible posterior distributions of the Pearson correlation
coefficient obtained via [run_plausible_cor()]. It shows:

- Density traces for a random subset (`n_draws`) of individual posterior distributions,
  reflecting uncertainty both in individual parameter estimates and in generalising to the population.
- The mean posterior density across all draws, providing a population-level summary of the
  plausible correlation.

Subsampling with `n_draws` is used to reduce overplotting when many MCMC draws are available.
However, the mean posterior density is always based on the full data.
}
\seealso{
[run_plausible_cor()].
}
