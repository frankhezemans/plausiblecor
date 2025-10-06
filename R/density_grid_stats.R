#' @title Compute Point Estimate from Discretised Density
#'
#' @description Internal helper function to compute a point estimate (mean,
#' median, or arbitrary quantile) from a numeric vector of values and a
#' corresponding density over a discretised grid.
#'
#' @param val Numeric vector representing the grid of values over which the
#'        density is defined.
#' @param dens Numeric vector of the same length as `val`, representing the
#'        density at each grid point.
#' @param method Character string specifying which point estimate to compute.
#'        One of `"mean"`, `"median"`, or `"quantile"`.
#' @param prob Numeric value between 0 and 1, only required when
#'        `method = "quantile"`. Specifies the cumulative probability at which
#'        to evaluate the quantile.
#' @param dx Grid spacing (step size). If `NULL` (default), it is estimated
#'        from the unique values of `val` (grid points) using [diff()].
#'
#' @return A numeric value corresponding to the requested point estimate.
#'
#' @keywords internal
get_point_estimate <- function(
    val,
    dens,
    method = c("mean", "median", "quantile"),
    prob = NULL,
    dx = NULL
) {

  if (length(val) != length(dens)) {
    rlang::abort(
      message = "Input arguments 'val' and 'dens' must be of equal length."
    )
  }
  method <- rlang::arg_match(arg = method)
  dx <- dx %||% get_dx(val)

  if (method == "mean") {
    result <- sum(val * dens) * dx
  } else {
    cdf <- cumsum(dens) * dx
    if (method == "median" && is.null(prob)) {
      prob <- 0.5
    }
    if (is.null(prob)) {
      rlang::abort(message = "No value for 'prob' provided.")
    }
    result <- stats::approx(
      x = cdf,
      y = val,
      xout = prob,
      ties = "ordered"
    )
    result <- result[["y"]]
  }

  return(result)

}

#' @title Compute Posterior Interval from Discretised Density
#'
#' @description Internal helper function to compute credible intervals
#' (quantile interval or highest-density continuous interval) from a numeric
#' vector of values and a corresponding density over a discretised grid.
#'
#' @param val Numeric vector representing the grid of values over which the
#'        density is defined.
#' @param dens Numeric vector of the same length as `val`, representing the
#'        density at each grid point.
#' @param width Numeric vector of desired interval widths
#'        (e.g., `0.95` for a 95% interval).
#' @param method Character string specifying which interval type to compute.
#'        One of `"qi"` (quantile interval) or `"hdci"`
#'        (highest-density continuous interval, default).
#' @param dx Grid spacing (step size). If `NULL` (default), it is estimated
#'        from the unique values of `val` (grid points) using [diff()].
#'
#' @return A data frame with columns `"lower"`, `"upper"`, and `"width"` for
#'         each requested interval.
#'
#' @keywords internal
get_interval <- function(
    val,
    dens,
    width,
    method = c("hdci", "qi"),
    dx = NULL
) {

  if (length(val) != length(dens)) {
    rlang::abort(
      message = "Input arguments 'val' and 'dens' must be of equal length."
    )
  }
  method <- rlang::arg_match(arg = method)
  dx <- dx %||% get_dx(val)

  cdf <- cumsum(dens) * dx

  compute_interval <- function(w) {
    if (method == "qi") {
      qi_approx <- stats::approx(
        x = cdf,
        y = val,
        xout = c((1 - w) / 2, (1 + w) / 2),
        ties = "ordered"
      )
      out <- qi_approx[["y"]]
    } else {
      df <- data.frame(val = val, dens = dens)
      df <- df[order(-df$dens), ]
      df$cum_dens <- cumsum(df$dens * dx)
      included_vals <- df$cum_dens <= w
      out <- range(df$val[included_vals])
    }
    out <- c(out, w)
    names(out) <- c("lower", "upper", "width")
    return(out)
  }

  result <- purrr::map(
    .x = width,
    .f = compute_interval
  ) %>%
    dplyr::bind_rows()

  return(result)

}

#' @noRd
get_dx <- function(x, tol = 1e-12) {
  dxs <- diff(sort(unique(x)))
  if (max(dxs) - min(dxs) > tol) {
    rlang::warn(
      message = "Grid appears uneven; using first diff() as dx."
    )
  }
  dx <- dxs[1]
  return(dx)
}
