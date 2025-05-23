% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plausible_cor.R
\name{get_mean_posterior_rho}
\alias{get_mean_posterior_rho}
\title{Compute mean posterior correlation density}
\usage{
get_mean_posterior_rho(.data)
}
\arguments{
\item{.data}{Data frame output from [get_posterior_rho_densities()].}
}
\value{
A tibble with columns `x` (grid points) and `density`
        (mean posterior density).
}
\description{
Helper function that computes the mean posterior density across MCMC samples,
based on the evaluated densities on a shared grid. Typical users would not
have to call this function, instead relying on [summarise_plausible_cor()],
[compare_plausible_cors()], and [plot_population_cor()].
}
