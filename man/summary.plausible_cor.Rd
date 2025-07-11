% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plausible_cor.R
\name{summary.plausible_cor}
\alias{summary.plausible_cor}
\title{Summary method for plausible correlation objects}
\usage{
\method{summary}{plausible_cor}(object, point_interval_args = NULL, rope_range = NULL, ...)
}
\arguments{
\item{object}{An object of class "plausible_cor" (output from \code{\link[=run_plausible_cor]{run_plausible_cor()}}).}

\item{point_interval_args}{A named list specifying how the central tendency
and credible interval(s) should be computed. See \code{\link[=summarise_plausible_cor]{summarise_plausible_cor()}}
for details.}

\item{rope_range}{Optional numeric vector of length 2 specifying the ROPE bounds.
See \code{\link[=summarise_plausible_cor]{summarise_plausible_cor()}} for details.}

\item{...}{Additional arguments (currently unused)}
}
\value{
A summary object with class "summary.plausible_cor"
}
\description{
Summary method for plausible correlation objects
}
