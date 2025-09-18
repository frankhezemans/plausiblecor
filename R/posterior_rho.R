#' Calculate Unnormalised Posterior Density for Pearson's Correlation
#' Coefficient
#'
#' @description
#' This function calculates the unnormalised posterior density for Pearson's
#' correlation coefficient based on methods introduced by Ly et al. (2018).
#' It returns an approximation function that can be evaluated at any correlation
#' value between -1 and 1.
#'
#' @param r Numeric value. The observed sample correlation coefficient, must be
#'        between `-1` and `1`.
#' @param n Numeric value. The sample size, must be at least `3`.
#' @param kappa Numeric value. Parameter controlling the "concentration" of the
#'        stretched beta prior on the correlation coefficient (see details
#'        below). Default is `1`, resulting in a uniform prior.
#' @param alternative Character string specifying the alternative hypothesis:
#'   \describe{
#'     \item{"two.sided"}{ (default) Prior is supported on \eqn{\left[-1, 1\right]}}
#'     \item{"greater"}{ Prior truncated to \eqn{\left[0, 1\right]} and renormalised}
#'     \item{"less"}{ Prior truncated to \eqn{\left[-1, 0\right]} and renormalised}
#'   }
#' @param n_bins Integer. Number of grid points for the approximation, default
#'        is `1000`.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions that are necessary to compute
#'        the posterior density. Default is `1e7`. See details.
#' @param ... Additional arguments passed to [stats::approxfun()].
#'
#' @return A function that evaluates the unnormalized posterior density at any
#'         correlation value between `-1` and `1`.
#'
#' @details
#' The unnormalised probability density of the posterior distribution of the
#' correlation coefficient is computed for a discretised grid of values. Then,
#' linear interpolation between these density values is used to construct an
#' approximation to the unnormalised probability density *function* of the
#' posterior, using [stats::approxfun()].
#'
#' The discrete grid of correlation values is set adaptively, so that there are
#' many estimates of the density near the posterior mode (i.e., the observed
#' correlation) and near the endpoints of +1 and -1.
#'
#' Computing the posterior density involves several evaluations of the
#' generalised hypergeometric function, which is handled internally by
#' [hypergeo::genhypergeo()]. In some cases, especially when the observed
#' correlation approaches +1 or -1, [hypergeo::genhypergeo()] may need several
#' thousand attempts to find the solution. Here, we cautiously set the default
#' `max_iter` to `1e7`, but note that the original default in
#' [hypergeo::genhypergeo()] is much lower (`2e3`).
#'
#' Any undefined or non-finite estimate of the posterior density is treated as
#' `NA` and ignored by [stats::approxfun()] when constructing the function.
#'
#' The prior is a stretched beta distribution with shape parameters
#' \eqn{\alpha = \beta = \frac{1}{\kappa}}, scaled to the interval (-1, 1). This
#' creates a symmetric distribution centered at zero and its domain stretched to
#' cover the full range of the correlation coefficient. The prior can optionally
#' also be truncated (and renormalised) to support a directional (strictly positive or negative)
#' alternative hypothesis.
#'
#' This function was adapted from code previously released with the Dynamic
#' Models of Choice toolbox (Heathcote et al., 2019).
#'
#' @references
#' Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Analytic posteriors for
#' Pearson's correlation coefficient. *Statistica Neerlandica*, 72, 4–13.
#' https://doi.org/10.1111/stan.12111
#'
#' Heathcote, A., Lin, Y.S., Reynolds, A., Strickland, L., Gretton, M., &
#' Matzke, D. (2019). Dynamic models of choice. *Behavior Research Methods*,
#' 51, 961-985. https://doi.org/10.3758/s13428-018-1067-y
#'
#' @export
posterior_rho_updf <- function(
    r, n, kappa = 1,
    alternative = c("two.sided", "greater", "less"),
    n_bins = 1e3, max_iter = 1e7, ...
) {

  alternative <- rlang::arg_match(alternative)

  params <- validate_posterior_rho_updf_input(r, n, kappa, n_bins, max_iter)
  r <- params[["r"]]
  n <- params[["n"]]
  kappa <- params[["kappa"]]
  n_bins <- params[["n_bins"]]
  max_iter <- params[["max_iter"]]

  null_approxfun <- function(x) rep(NA_real_, length(x))

  if (abs(r) == 1) {
    rlang::warn(
      message = paste0(
        "Perfect correlation (|r| = 1) detected. ",
        "Returning degenerate function that always returns NA."
      ),
      .frequency = "once",
      .frequency_id = "perfect_correlation"
    )
    return(null_approxfun)
  }

  rho_grid <- create_rho_grid(r, n, n_bins)

  if (abs(r) < sqrt(.Machine$double.eps)) {
    d <- posterior_rho_jeffreys(r, n, rho_grid, kappa, alternative, max_iter)
  } else {
    d <- posterior_rho_exact(r, n, rho_grid, kappa, alternative, max_iter)
  }

  if (any(!is.finite(d))) {
    if (all(!is.finite(d))) {
      rlang::warn(
        message = paste0(
          "Computation of posterior density failed. ",
          "Returning degenerate function that always returns NA."
        )
      )
      return(null_approxfun)
    }
    d[!is.finite(d)] <- NA_real_
    rlang::warn(
      message = paste0(
        "Some posterior densities were non-finite; will be ignored for ",
        "constructing the UPDF."
      )
    )
  }

  result <- stats::approxfun(x = rho_grid, y = d, na.rm = TRUE, ...)

  return(result)

}

#' Calculate exact posterior density for Pearson's correlation
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param rho_grid Numeric vector. Grid of correlation values to evaluate.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param alternative Character string specifying the alternative hypothesis.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric vector of unnormalized posterior density values.
#' @noRd
posterior_rho_exact <- function(r, n, rho_grid, kappa, alternative, max_iter) {
  result <- bf_rho_exact(r, n, kappa, max_iter) *
    likelihood_rho_exact(r, n, rho_grid, max_iter) *
    prior_rho(rho_grid, kappa, alternative)
  return(result)
}

#' Calculate Jeffreys approximation of posterior density for Pearson's correlation
#'
#' @description
#' Used when the observed correlation coefficient is effectively zero.
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param rho_grid Numeric vector. Grid of correlation values to evaluate.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param alternative Character string specifying the alternative hypothesis.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric vector of unnormalized posterior density values.
#' @noRd
posterior_rho_jeffreys <- function(r, n, rho_grid, kappa, alternative, max_iter) {
  result <- bf_rho_jeffreys(r, n, kappa, max_iter) *
    likelihood_rho_jeffreys(r, n, rho_grid) *
    prior_rho(rho_grid, kappa, alternative)
  return(result)
}

#' Calculate Bayes factor for exact posterior correlation
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric value representing the Bayes factor.
#' @noRd
bf_rho_exact <- function(r, n, kappa, max_iter) {

  hyper_term <- hypergeo::genhypergeo(
    U = c((n - 1) / 2, (n - 1) / 2),
    L = (n + 2 / kappa) / 2,
    z = r^2,
    maxiter = max_iter
  )

  log_result <- log(2^(1 - 2 / kappa)) +
    0.5 * log(pi) -
    lbeta(a = 1 / kappa, b = 1 / kappa) +
    lgamma((n + 2 / kappa - 1) / 2) -
    lgamma((n + 2 / kappa) / 2) +
    log(hyper_term)

  result <- exp(Re(log_result))

  return(result)

}

#' Calculate Bayes factor for Jeffreys approximation of correlation
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric value representing the Bayes factor.
#' @noRd
bf_rho_jeffreys <- function(r, n, kappa, max_iter) {

  hyper_term <- hypergeo::genhypergeo(
    U = c((2 * n - 3) / 4, (2 * n - 1) / 4),
    L = (n + 2 / kappa) / 2,
    z = r^2,
    maxiter = max_iter
  )

  log_term <- lgamma((n + 2 / kappa - 1) / 2) -
    lgamma((n + 2 / kappa) / 2) -
    lbeta(a = 1 / kappa, b = 1 / kappa)

  result <- sqrt(pi) *
    2^(1 - 2 / kappa) *
    exp(log_term) *
    Re(hyper_term)

  return(result)

}

#' Calculate exact likelihood for Pearson's correlation
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param rho_grid Numeric vector. Grid of correlation values to evaluate.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric vector of likelihood values.
#' @noRd
likelihood_rho_exact <- function(r, n, rho_grid, max_iter) {

  hyper_term_even <- hypergeo::genhypergeo(
    U = c((n - 1) / 2, (n - 1) / 2),
    L = 1 / 2,
    z = (r * rho_grid)^2,
    maxiter = max_iter
  )

  likelihood_rho_even <- (1 - rho_grid^2)^((n - 1) / 2) *
    Re(hyper_term_even)

  hyper_term_odd <- hypergeo::genhypergeo(
    U = c(n / 2, n / 2),
    L = 3 / 2,
    z = (r * rho_grid)^2,
    maxiter = max_iter
  )

  log_term_odd <- 2 * (lgamma(n / 2) - lgamma((n - 1) / 2)) +
    ((n - 1) / 2) *
    log(1 - rho_grid^2)

  likelihood_rho_odd <- 2 * r * rho_grid *
    exp(log_term_odd) *
    Re(hyper_term_odd)

  result <- likelihood_rho_even + likelihood_rho_odd

  return(result)

}

#' Calculate Jeffreys approximation of likelihood for Pearson's correlation
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param rho_grid Numeric vector. Grid of correlation values to evaluate.
#'
#' @return Numeric vector of likelihood values.
#' @noRd
likelihood_rho_jeffreys <- function(r, n, rho_grid) {
  numerator <- (1 - rho_grid^2)^((n - 1) / 2)
  denominator <- (1 - rho_grid * r)^(n - 1 - (1 / 2))
  result <- numerator / denominator
  return(result)
}

#' Calculate prior density for correlation coefficient
#'
#' @description
#' Calculates the stretched beta prior density for the correlation coefficient,
#' with optional truncation for directional alternatives.
#'
#' @param rho_grid Numeric vector. Grid of correlation values to evaluate.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param alternative Character string specifying the alternative hypothesis.
#'
#' @return Numeric vector of prior density values.
#' @noRd
prior_rho <- function(rho_grid, kappa, alternative) {

  base_density <- (1 / 2) *
    stats::dbeta(
      x = (rho_grid + 1) / 2,
      shape1 = 1 / kappa,
      shape2 = 1 / kappa
    )

  if (alternative == "two.sided") {
    result <- base_density
  } else {
    result <- numeric(length(rho_grid)) # initialise with zero density
    if (alternative == "greater") {
      idx <- rho_grid >= 0
    } else {
      idx <- rho_grid <= 0
    }
    # overwrite valid rho values with base density, normalised to account for
    # truncation by multiplying by 2 (NB with base prior being symmetric around
    # zero, normalisation constant for truncating at 0 is always exactly 1/2).
    result[idx] <- 2 * base_density[idx]
  }

  return(result)
}

#' Create grid of correlation coefficient values
#'
#' @description
#' Creates a grid of correlation coefficient values for numerical approximation
#' of posterior densities. The grid begins with evenly spaced coverage from
#' -0.999 to 0.999, and then more points are added around the observed
#' correlation for greater precision around the peak of the distribution.
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#'   Must be between -1 and 1.
#' @param n Numeric value. The sample size. Used to determine the spread of
#'   points around the observed correlation via Fisher z-transformation.
#' @param n_bins Integer. Total number of grid points desired. Will be allocated
#'   proportionally: 70% for base coverage and 30% around observed correlation.
#'   Minimum value of 100 is enforced.
#'
#' @return Numeric vector of correlation values between -0.999 and 0.999,
#'   sorted in ascending order. The actual number of returned points may
#'   be slightly less than `n_bins` due to de-duplication.
#'
#' @noRd
create_rho_grid <- function(r, n, n_bins) {

  n_bins <- max(n_bins, 100)

  n_base <- round(n_bins * 0.7)
  base_grid <- seq(from = -0.999, to = 0.999, length.out = n_base)

  n_peak <- n_bins - n_base
  unit_grid <- stats::qlogis(
    p = seq(from = 0, to = 1, length.out = n_peak + 2)
  )
  unit_grid <- unit_grid[-c(1, length(unit_grid))]
  peak_grid <- tanh(
    atanh(r) + unit_grid / sqrt(n)
  )
  peak_grid <- pmax(-0.999, pmin(0.999, peak_grid))

  result <- sort(unique(c(
    base_grid, peak_grid
  )))

  return(result)

}
