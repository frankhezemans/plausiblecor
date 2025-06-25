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
#' The prior is a stretched beta distribution with shape parameters alpha =
#' beta = 1/kappa, scaled to the interval [-1, 1]. This creates a symmetric
#' distribution centered at zero and its domain stretched to cover the full
#' range of the correlation coefficient.
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
    r, n, kappa = 1, n_bins = 1e3, max_iter = 1e7, ...
) {

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
    d <- posterior_rho_jeffreys(r, n, rho_grid, kappa, max_iter)
  } else {
    d <- posterior_rho_exact(r, n, rho_grid, kappa, max_iter)
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
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric vector of unnormalized posterior density values.
#' @noRd
posterior_rho_exact <- function(r, n, rho_grid, kappa, max_iter) {
  result <- bf_rho_exact(r, n, kappa, max_iter) *
    likelihood_rho_exact(r, n, rho_grid, max_iter) *
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
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric vector of unnormalized posterior density values.
#' @noRd
posterior_rho_jeffreys <- function(r, n, rho_grid, kappa, max_iter) {
  result <- bf_rho_jeffreys(r, n, kappa, max_iter) *
    likelihood_rho_jeffreys(r, n, rho_grid) *
    prior_rho(rho_grid, kappa)
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

#' Create grid of correlation coefficient values
#'
#' @description
#' Creates a grid of correlation coefficient values for numerical approximation
#' of posterior densities. The grid uses adaptive allocation to place points
#' strategically across three regions: around the observed correlation, near the
#' extreme values (-1 and 1), and a base coverage of the middle range.
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#'   Must be between -1 and 1.
#' @param n Numeric value. The sample size. Used to determine the spread of
#'   points around the observed correlation via Fisher z-transformation.
#' @param n_bins Integer. Total number of grid points desired. Will be allocated
#'   proportionally: ~30% for extremes, ~40% around observed correlation,
#'   remainder for base coverage. Minimum value of 10 is enforced.
#'
#' @return Numeric vector of correlation values between -0.999 and 0.999,
#'   sorted in ascending order. The function aims to return approximately
#'   n_bins points, but the actual number may be somewhat fewer due to
#'   de-duplication.
#'
#' @details
#' The function uses Fisher's z-transformation to create normally distributed
#' points around the observed correlation in z-space, then transforms back to
#' correlation space. The spread in z-space is 1 standard error (1/sqrt(n)),
#' providing conservative coverage of plausible correlation values.
#'
#' For extreme regions, the function uses the beta distribution to create
#' higher point density near the actual boundaries (±1). The boundary
#' locations adapt based on the observed correlation to provide better
#' coverage on the relevant side.
#'
#' @noRd
create_rho_grid <- function(r, n, n_bins) {

  n_bins <- max(n_bins, 10)

  extreme_range <- 0.15 * (1 - abs(r))
  if (r >= 0) {
    lower_extreme <- -0.8 - extreme_range
    upper_extreme <- 0.999
  } else {
    lower_extreme <- -0.999
    upper_extreme <- 0.8 + extreme_range
  }

  n_extremes <- round(n_bins * 0.3)
  n_lower_extreme <- ceiling(n_extremes / 2)
  n_upper_extreme <- n_extremes - n_lower_extreme
  n_peak <- round(n_bins * 0.4)
  n_base <- n_bins - n_extremes - n_peak

  peak_grid <- tanh(seq(
    from = atanh(r) - 1 / sqrt(n),
    to = atanh(r) + 1 / sqrt(n),
    length.out = n_peak
  ))
  peak_grid <- pmax(-0.999, pmin(0.999, peak_grid))

  base_start <- lower_extreme + 0.1
  base_end <- upper_extreme - 0.1
  if (base_end - base_start < 0.3) {
    base_start <- -0.7
    base_end <- 0.7
  }
  base_grid <- seq(
    from = base_start,
    to = base_end,
    length.out = n_base
  )

  if (n_lower_extreme == 1) {
    lower_extreme_grid <- (lower_extreme + base_start - 0.1) / 2
  } else {
    lower_props <- stats::pbeta(
      q = seq(from = 0.05, to = 0.95, length.out = n_lower_extreme),
      shape1 = 0.5,
      shape2 = 2
    )
    lower_extreme_grid <- lower_extreme +
      (base_start - 0.1 - lower_extreme) * lower_props
  }

  if (n_upper_extreme == 1) {
    upper_extreme_grid <- (base_end + 0.1 + upper_extreme) / 2
  } else {
    upper_props <- stats::pbeta(
      q = seq(from = 0.05, to = 0.95, length.out = n_upper_extreme),
      shape1 = 2,
      shape2 = 0.5
    )
    upper_extreme_grid <- (base_end + 0.1) +
      (upper_extreme - base_end - 0.1) * upper_props
  }

  result <- sort(unique(c(
    lower_extreme_grid, base_grid, peak_grid, upper_extreme_grid
  )))
  result <- result[result >= -0.999 & result <= 0.999]

  return(result)

}

#' Validate inputs for posterior_rho_updf function
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Integer The sample size.
#' @param kappa Numeric value. Parameter controlling the concentration of the
#'        prior.
#' @param n_bins Integer. Number of grid points for the approximation.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Named list containing the validated inputs.
#'
#' @noRd
validate_posterior_rho_updf_input <- function(r, n, kappa, n_bins, max_iter) {

  rules <- list(
    r = list(
      fun = function(x) abs(x) <= 1,
      fail_msg = "between -1 and 1 (inclusive)",
      round = FALSE
    ),
    n = list(
      fun = function(x) x > 2,
      fail_msg = "greater than 2",
      round = TRUE
    ),
    kappa = list(
      fun = function(x) x > 0,
      fail_msg = "that is strictly positive",
      round = FALSE
    ),
    n_bins = list(
      fun = function(x) x >= 100,
      fail_msg = "greater than or equal to 100",
      round = TRUE
    ),
    max_iter = list(
      fun = function(x) x >= 1,
      fail_msg = "greater than or equal to 1",
      round = TRUE
    )
  )

  process_param <- function(value, name, rule_list = rules) {
    rule <- rule_list[[name]]
    if (is.null(rule)) {
      return(value)
    }
    if (
      length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      !rule[["fun"]](value)
    ) {
      rlang::abort(
        message = paste0(
          "Input '", name, "' must be a single finite value ",
          rule[["fail_msg"]], "."
        )
      )
    }
    if (rule[["round"]] && value != round(value)) {
      value <- round(value)
      rlang::warn(
        message = paste0(
          "Input '", name, "' has been rounded to the nearest integer."
        )
      )
    }
    return(value)
  }

  result <- purrr::imap(
    .x = list(
      r = r, n = n, kappa = kappa, n_bins = n_bins, max_iter = max_iter
    ),
    .f = process_param
  )

  return(result)

}
