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
#' @param n_bins Integer. Number of grid points for the approximation, default
#'        is `1000`.
#' @param ... Additional arguments passed to [stats::approxfun()].
#'
#' @return A function that evaluates the unnormalized posterior density at any
#'         correlation value between `-1` and `1`.
#'
#' @details
#' The function implements the analytical posterior for Pearson's correlation
#' coefficient as described in Ly et al. (2018). For r = 0, Jeffreys'
#' approxmation is used; otherwise the exact solution given by Ly et al. (2018)
#' is used.
#'
#' The prior on the correlation coefficient is a so-called stretched beta
#' distribution, which is a beta distribution with its shape parameters alpha
#' and beta set to alpha = beta = 1/kappa, and its minimum and maximum
#' stretched to -1 and 1, respectively. This creates a symmetric distribution
#' centered at zero and its domain stretched to cover the full range of the
#' correlation coefficient.
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
posterior_rho_updf <- function(r, n, kappa = 1, n_bins = 1e3, ...) {

  validate_posterior_rho_updf_input(r, n, kappa, n_bins)

  input_warn_rounded <- function(param) {
    rlang::warn(
      message = paste0(
        "Input '", param, "' has been rounded to the nearest integer."
      )
    )
  }

  if (n != round(n)) {
    input_warn_rounded(param = "n")
    n <- round(n)
  }

  if (n_bins != round(n_bins)) {
    input_warn_rounded(param = "n_bins")
    n_bins <- round(n_bins)
  }

  rho_grid <- create_rho_grid(r, n, n_bins)

  if (abs(r) < sqrt(.Machine$double.eps)) {
    d <- posterior_rho_jeffreys(r, n, rho_grid, kappa)
  } else {
    d <- posterior_rho_exact(r, n, rho_grid, kappa)
  }

  if (all(!is.finite(d))) {
    rlang::abort(message = "Computation of posterior density failed.")
  }
  if (any(!is.finite(d))) {
    rlang::warn(message = "Some posterior density values were non-finite.")
  }

  result <- stats::approxfun(x = rho_grid, y = d, ...)

  return(result)

}

#' Validate inputs for posterior_rho_updf function
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param kappa Numeric value. Parameter controlling the concentration of the
#'        prior.
#' @param n_bins Integer. Number of grid points for the approximation.
#'
#' @return NULL invisibly if validation passes, otherwise aborts with an error
#'         message.
#' @noRd
validate_posterior_rho_updf_input <- function(r, n, kappa, n_bins) {

  input_abort <- function(param, conditions) {
    rlang::abort(
      message = paste0(
        "Input '", param, "' must be a single finite value ", conditions, "."
      )
    )
  }

  if (length(r) != 1 || !is.numeric(r) || !is.finite(r) || abs(r) > 1) {
    input_abort(param = "r", conditions = "between -1 and 1 (inclusive)")
  }
  if (length(n) != 1 || !is.numeric(n) || !is.finite(n) || n < 3) {
    input_abort(param = "n", conditions = "greater than or equal to 3")
  }
  if (
    length(kappa) != 1 || !is.numeric(kappa) || !is.finite(kappa) || kappa <= 0
  ) {
    input_abort(param = "kappa", conditions = "that is strictly positive")
  }
  if (
    length(n_bins) != 1 || !is.numeric(n_bins) || !is.finite(n_bins) ||
    n_bins < 10
  ) {
    input_abort(param = "n_bins", conditions = "greater than or equal to 10")
  }

  return(invisible(x = NULL))

}

#' Create grid of correlation coefficient values
#'
#' @description
#' Creates a grid of correlation coefficient values for numerical approximation,
#' with higher density of points around the observed correlation coefficient
#' and around the endpoints -1 and +1.
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param n_bins Integer. Number of grid points desired.
#'
#' @return Numeric vector of correlation values between -1 and 1.
#' @noRd
create_rho_grid <- function(r, n, n_bins) {

  unit_grid <- stats::qlogis(
    p = seq(from = 0, to = 1, length.out = n_bins + 2)
  )
  unit_grid <- unit_grid[-c(1, length(unit_grid))]

  result <- unique(sort(c(
    -1,
    seq(-0.999, -0.9, by = 0.01),
    tanh(atanh(r) + unit_grid / sqrt(n)),
    seq(0.9, 0.999, by = 0.01),
    1
  )))

  return(result)

}

#' Calculate exact posterior density for Pearson's correlation
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param rho_grid Numeric vector. Grid of correlation values to evaluate.
#' @param kappa Numeric value. Prior concentration parameter.
#'
#' @return Numeric vector of unnormalized posterior density values.
#' @noRd
posterior_rho_exact <- function(r, n, rho_grid, kappa) {
  result <- bf_rho_exact(r, n, kappa) *
    likelihood_rho_exact(r, n, rho_grid) *
    prior_rho(rho_grid, kappa)
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
#'
#' @return Numeric vector of unnormalized posterior density values.
#' @noRd
posterior_rho_jeffreys <- function(r, n, rho_grid, kappa) {
  result <- bf_rho_jeffreys(r, n, kappa) *
    likelihood_rho_jeffreys(r, n, rho_grid) *
    prior_rho(rho_grid, kappa)
  return(result)
}

#' Calculate Bayes factor for exact posterior correlation
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param kappa Numeric value. Prior concentration parameter.
#'
#' @return Numeric value representing the Bayes factor.
#' @noRd
bf_rho_exact <- function(r, n, kappa) {

  if (n <= 2) {
    return(1)
  }

  if (kappa >= 1 && abs(r) >= 1) {
    return(Inf)
  }

  hyper_term <- hypergeo::genhypergeo(
    U = c((n - 1) / 2, (n - 1) / 2),
    L = (n + 2 / kappa) / 2,
    z = r^2
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
#'
#' @return Numeric value representing the Bayes factor.
#' @noRd
bf_rho_jeffreys <- function(r, n, kappa) {

  if (n <= 2) {
    return(1)
  }

  if (abs(r) >= 1) {
    return(Inf)
  }

  hyper_term <- hypergeo::genhypergeo(
    U = c((2 * n - 3) / 4, (2 * n - 1) / 4),
    L = (n + 2 / kappa) / 2,
    z = r^2
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
#'
#' @return Numeric vector of likelihood values.
#' @noRd
likelihood_rho_exact <- function(r, n, rho_grid) {

  hyper_term_even <- hypergeo::genhypergeo(
    U = c((n - 1) / 2, (n - 1) / 2),
    L = 1 / 2,
    z = (r * rho_grid)^2
  )

  likelihood_rho_even <- (1 - rho_grid^2)^((n - 1) / 2) *
    Re(hyper_term_even)

  hyper_term_odd <- hypergeo::genhypergeo(
    U = c(n / 2, n / 2),
    L = 3 / 2,
    z = (r * rho_grid)^2
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
#' Calculates the stretched beta prior density for the correlation coefficient.
#'
#' @param rho_grid Numeric vector. Grid of correlation values to evaluate.
#' @param kappa Numeric value. Prior concentration parameter.
#'
#' @return Numeric vector of prior density values.
#' @noRd
prior_rho <- function(rho_grid, kappa) {
  result <- (1 / 2) *
    stats::dbeta(
      x = (rho_grid + 1) / 2,
      shape1 = 1 / kappa,
      shape2 = 1 / kappa
    )
  return(result)
}
