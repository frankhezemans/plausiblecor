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

  assert_val_dens_lengths(val, dens)
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

  assert_val_dens_lengths(val, dens)
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

#' @title Compute Probability of Direction from Discretised Density
#'
#' @description Internal helper function to compute the probability of direction
#' from a numeric vector of values and a corresponding density over a
#' discretised grid.
#'
#' @param val Numeric vector representing the grid of values over which the
#'        density is defined.
#' @param dens Numeric vector of the same length as `val`, representing the
#'        density at each grid point.
#' @param alternative Character string specifying the alternative hypothesis.
#'        One of `"two.sided"` (default), `"greater"`, `"less"`.
#' @param dx Grid spacing (step size). If `NULL` (default), it is estimated
#'        from the unique values of `val` (grid points) using [diff()].
#'
#' @details
#' The probability of direction (PD) represents the proportion of the
#' distribution that lies on one side of zero. For example, a PD of `0.92` means
#' that 92% of the distribution's mass lies one one side (strictly positive or
#' strictly negative values, whichever is greatest). In the context of
#' hypothesis testing, it can be taken as an index of the _existence_ of an
#' effect, but does _not_ necessarily address the size, importance, or
#' _"significance"_ of an effect. For example, a vanishingly small effect whose
#' distribution is precisely concentrated within the
#' \eqn{\left[0.0001, 0.0002\right]} range would still have a PD of
#' approximately `1`.
#'
#' Note that the notion of the probability of direction is undefined for one-sided
#' hypotheses. If `alternative != "two.sided`, the returned value is `NA_real_`.
#'
#' @return A numeric value corresponding to the probability of direction.
#'
#' @keywords internal
get_p_dir <- function(
    val,
    dens,
    alternative = c("two.sided", "greater", "less"),
    dx = NULL
) {

  assert_val_dens_lengths(val, dens)
  alternative <- rlang::arg_match(arg = alternative)
  if (alternative != "two.sided") {
    return(NA_real_)
  }
  dx <- dx %||% get_dx(val)

  result <- max(
    sum(dens[val > 0]) * dx,
    sum(dens[val < 0]) * dx
  )
  return(result)

}

#' @title Perform ROPE-based Equivalence Test from Discretised Density
#'
#' @description Internal helper function to perform an equivalence test, based
#' on a numeric vector of values and a corresponding density over a discretised
#' grid. In equivalence testing, the null hypothesis is that a parameter value
#' falls within a specified region of practical equivalence (ROPE). The
#' proportion of the distribution (or a credible interval thereof) contained
#' within the ROPE can be used as a decision rule regarding the "significance"
#' of a parameter.
#'
#' @param val Numeric vector representing the grid of values over which the
#'        density is defined.
#' @param dens Numeric vector of the same length as `val`, representing the
#'        density at each grid point.
#' @param rope_range Numeric vector of length 2 specifying the lower and upper
#'        bounds of the ROPE.
#' @param ci_range Optional numeric vector of length 2 specifying the lower and
#'        bounds of a credible interval (CI) of the distribution. If `NULL`
#'        (default), the returned value represents the proportion of the *full distribution*
#'        contained within the ROPE. If a CI of the distribution is provided,
#'        the returned value represents the proportion of the *distribution's CI*
#'        contained within the ROPE.
#' @param alternative Character string specifying the alternative hypothesis.
#'        One of `"two.sided"` (default), `"greater"`, `"less"`.
#' @param dx Grid spacing (step size). If `NULL` (default), it is estimated
#'        from the unique values of `val` (grid points) using [diff()].
#'
#' @details
#' The ROPE-based p-value allows for a pragmatic interpretation of the
#' "significance" of an effect, indicating if a parameter is effectively zero
#' for all practical purposes. It is commonly computed as the proportion of the
#' distribution's 95% CI contained within the ROPE (Kruschke & Liddell, 2018).
#' If `p_rope` is approximately equal to 1, (almost) all of the distribution's
#' most-credible values are inside the ROPE, suggesting that the null hypothesis
#' can be accepted. If `p_rope` is approximately equal to 0, (almost) all of
#' the distribution's most-credible values are outside the ROPE, suggesting
#' that the null hypothesis can be rejected. Otherwise, it remains unclear if
#' the null hypothesis should be accepted or rejected.
#' However, some have argued that computing the proportion of the *full distribution*
#' (as opposed to the *distribution's CI*) contained within the ROPE provides a
#' more sensitive test (Makowski et al., 2019). In this case, `p_rope` values
#' less than 0.025 and greater than 0.975 are used as thresholds to reject or
#' accept the null hypothesis, respectively.
#'
#' @return A numeric value corresponding to the proportion of the distribution
#'  (or credible interval thereof) contained within the ROPE.
#'
#' @references
#' Makowski, D., Ben-Shachar, M. S., Chen, S. H. A., & Lüdecke, D. (2019).
#' Indices of effect existence and significance in the Bayesian framework.
#' *Frontiers in Psychology*, *10*. \doi{10.3389/fpsyg.2019.02767}
#'
#' Kruschke, J. K., Liddell, T. M. (2018). The Bayesian new statistics:
#' Hypothesis testing, estimation, meta-analysis, and power analysis from a
#' Bayesian perspective. *Psychonomic Bulletin & Review*, *25*(1), 178-206.
#' \doi{10.3758/s13423-016-1221-4}
#'
#' @keywords internal
get_p_rope <- function(
    val,
    dens,
    rope_range,
    ci_range = NULL,
    alternative = c("two.sided", "greater", "less"),
    dx = NULL
) {

  assert_val_dens_lengths(val, dens)
  alternative <- rlang::arg_match(arg = alternative)
  rope_range <- validate_rope_range(x = rope_range, alternative = alternative)
  if (is.null(rope_range)) {
    return(NULL)
  }
  dx <- dx %||% get_dx(val)

  dens <- dens / sum(dens * dx)

  if (!is.null(ci_range)) {
    ci_mask <- val >= ci_range[1] & val <= ci_range[2]
  } else {
    ci_mask <- rep(TRUE, times = length(val))
  }
  val_ci <- val[ci_mask]
  dens_ci <- dens[ci_mask]
  dens_ci <- dens_ci / sum(dens_ci * dx)

  result <- get_area_in_range(
    x = val_ci,
    density = dens_ci,
    lower = rope_range[1],
    upper = rope_range[2],
    dx = dx
  )
  return(result)

}

#' @importFrom utils tail
#' @noRd
get_area_in_range <- function(x, density, lower, upper, dx) {
  in_range <- x >= lower & x <= upper
  x_sub <- x[in_range]
  d_sub <- density[in_range]
  if (length(x_sub) < 2) {
    return(0)
  }
  result <- dx * sum(
    (utils::head(d_sub, n = -1) + utils::tail(d_sub, n = -1)) / 2
  )
  return(result)
}


#' @noRd
assert_val_dens_lengths <- function(val, dens) {
  if (length(val) != length(dens)) {
    rlang::abort(
      message = "Input arguments 'val' and 'dens' must be of equal length."
    )
  }
  return(invisible(x = NULL))
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
