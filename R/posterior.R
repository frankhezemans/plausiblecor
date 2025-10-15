#' Approximation of Posterior Density Function for Correlation Coefficient
#'
#' @description
#' This function calculates an approximation to the posterior density of
#' Pearson's correlation coefficient (Ly et al., 2018) or Kendall's rank
#' correlation coefficient (Van Doorn et al., 2018). The result is returned as a
#' function that is proportional to the posterior density (i.e., an
#' *unnormalised*) posterior density function), constructed via linear
#' interpolation between grid points.
#'
#' @param r Numeric value. The observed sample correlation coefficient, must be
#'        between `-1` and `1`.
#' @param n Numeric value. The sample size, must be at least `3`.
#' @param kappa Numeric value. Parameter controlling the "concentration" of the
#'        stretched beta prior on the correlation coefficient (see details
#'        below). If `NULL` (default), it is set to a value that induces a
#'        uniform prior; namely `1` for `method = "pearson"` and `2` for
#'        `method = "kendall"`.
#' @param alternative Character string specifying the alternative hypothesis:
#'   \describe{
#'     \item{"two.sided"}{ (default) Prior is supported on \eqn{\left[-1, 1\right]}}
#'     \item{"greater"}{ Prior truncated to \eqn{\left[0, 1\right]} and re-normalised}
#'     \item{"less"}{ Prior truncated to \eqn{\left[-1, 0\right]} and re-normalised}
#'   }
#' @param method Character string specifying which correlation coefficient is to
#'        be used for the test. One of `"pearson"` (default) or `"kendall"`.
#' @param n_bins Integer. Number of grid points for the approximation, default
#'        is `1e3`.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions that are necessary to compute
#'        the posterior density for `method = "pearson"`. Default is `1e7`. See
#'        details.
#'
#' @return A function that evaluates the unnormalized posterior density at any
#'         correlation value between `-1` and `1`. The function is proportional
#'         to the true posterior density but may not integrate to 1 due to
#'         interpolation.
#'
#' @details
#' The probability density of the posterior distribution of the correlation
#' coefficient is computed for a discretised grid of values. Then,
#' linear interpolation between these density values is used to construct an
#' approximation to the probability density *function* of the posterior, using
#' [stats::approxfun()].
#'
#' The discrete grid of correlation values is set adaptively, so that there are
#' many estimates of the density near the posterior mode (i.e., the observed
#' correlation) and near the endpoints of +1 and -1.
#'
#' For `method = "pearson"`, computing the posterior density involves several
#' evaluations of the generalised hypergeometric function, which is handled
#' internally by [hypergeo::genhypergeo()]. In some cases, especially when the
#' observed correlation approaches +1 or -1, [hypergeo::genhypergeo()] may need
#' several thousand attempts to find the solution. Here, we cautiously set the
#' default `max_iter` to `1e7`, but note that the original default in
#' [hypergeo::genhypergeo()] is much lower (`2e3`).
#'
#' For `method = "kendall"`, the posterior density is based on the asymptotic
#' normal approximation to the sampling distribution of Kendall's rank
#' correlation coefficient, tau (see Van Doorn, 2018, and references therein
#' for details). In practice, the quality of this normal approximation depends
#' on both the value of the observed tau `r` and the sample size `n`: the
#' approximation improves exponentially as `n` increases, but larger values of
#' `r` require larger `n` to achieve the same accuracy. Users should be cautious
#' when applying this method to small samples and/or strong correlations.
#' A warning is issued if the given values of `r` and `n` fail to meet
#' rule-of-thumb criteria based on simulation results from Van Doorn et al. (2018).
#'
#' Any undefined or non-finite estimate of the posterior density is treated as
#' `NA` and ignored by [stats::approxfun()] when constructing the function.
#'
#' The prior is a stretched beta distribution with shape parameters
#' \eqn{\alpha = \beta = \frac{1}{\kappa}}, scaled to the interval (-1, 1). This
#' creates a symmetric distribution centered at zero and its domain stretched to
#' cover the full range of the correlation coefficient. The prior can optionally
#' also be truncated (and re-normalised) to support a directional
#' (i.e., non-negative or non-positive) alternative hypothesis.
#'
#' This function was primarily based on code that was originally written by
#' Alexander Ly, adapted by Dora Matzke, and then released with the Dynamic
#' Models of Choice toolbox (Heathcote et al., 2019). Furthermore, code for
#' Kendall's rank correlation was adapted from code written by Alexander Ly and
#' Johnny van Doorn that was released with the `bstats` R package (which is used
#' under the hood by the stand-alone statistics programme JASP).
#'
#' @references
#' Heathcote, A., Lin, Y.S., Reynolds, A., Strickland, L., Gretton, M., &
#' Matzke, D. (2019). Dynamic models of choice. *Behavior Research Methods*,
#' *51*, 961-985. \doi{10.3758/s13428-018-1067-y}
#'
#' Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Analytic posteriors for
#' Pearson's correlation coefficient. *Statistica Neerlandica*, *72*, 4–13.
#' \doi{10.1111/stan.12111}
#'
#' Van Doorn, J., Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Bayesian
#' inference for Kendall's rank correlation coefficient. *The American Statistician*,
#' *72*, 303-308. \doi{10.1080/00031305.2016.1264998}
#'
#' @examples
#' # Posterior for Pearson's correlation coefficient --------------------------
#'
#' # example: correlation between cars' horsepower and quarter mile race time
#' # (will be strongly negative: cars with more HP tend to need less time)
#' r <- cor(mtcars$hp, mtcars$qsec)
#' n <- nrow(mtcars)
#' post_fun <- posterior_cor_updf(r, n)
#' # this function can be used to obtain the (unnormalised) posterior density
#' # of a correlation value
#' post_fun(r)
#' post_fun(c(-0.25, 0.25))
#' # for visualisation, obtain the posterior density for an equally-spaced grid
#' # of correlation values
#' grid_res <- 1e-3
#' grid <- seq(from = -1, to = 1, by = grid_res)
#' post_dens <- post_fun(grid)
#' plot(
#'   x = grid, y = post_dens, type = "l",
#'   xlab = "Pearson correlation", ylab = "",
#'   bty = "n", yaxt = "n"
#' )
#'
#' # for inference, we first need to normalise the grid of posterior densities
#' normalise <- function(dens, dx = grid_res) {
#'   return(dens / sum(dens * dx))
#' }
#' post_dens <- normalise(post_dens)
#' post_mode <- grid[which.max(post_dens)]
#' post_median <- approx(
#'   x = cumsum(post_dens) * grid_res,
#'   y = grid,
#'   xout = 0.5,
#'   ties = "ordered"
#' )[["y"]]
#' post_mean <- sum(grid * post_dens) * grid_res
#'
#' # Directional hypotheses ---------------------------------------------------
#'
#' # the alternative hypothesis can be two-sided or one-sided (strictly non-
#' # negative or strictly non-positive).
#' # example: more ambiguous correlation (cars' weight and quarter mile time)
#' r_weak <- cor(mtcars$wt, mtcars$qsec)
#' post_fun <- posterior_cor_updf(r_weak, n)
#' # specify hypothesis that correlation is strictly non-positive
#' post_fun_lower <- posterior_cor_updf(r_weak, n, alternative = "less")
#' # compare posterior distributions for two-sided and one-sided hypotheses
#' plot(
#'   x = grid, y = normalise(post_fun_lower(grid)), type = "l", col = "red",
#'   xlab = "Pearson correlation", ylab = "",
#'   bty = "n", yaxt = "n"
#' )
#' lines(x = grid, y = normalise(post_fun(grid)))
#'
#' # Kendall's rank correlation coefficient -----------------------------------
#'
#' # Kendall's correlation is invariant to monotonic (rank-preserving)
#' # transformations of data, as it is based on ranks rather than values
#' # example: cars' HP is positively skewed (most cars with low-to-moderate HP,
#' # some cars with high HP); use square root to transform
#' r_pearson <- cor(mtcars$hp, mtcars$drat)
#' r_pearson_trans <- cor(sqrt(mtcars$hp), mtcars$drat)
#' identical(r_pearson, r_pearson_trans)
#' r_kendall <- cor(mtcars$hp, mtcars$drat, method = "kendall")
#' r_kendall_trans <- cor(sqrt(mtcars$hp), mtcars$drat, method = "kendall")
#' identical(r_kendall, r_kendall_trans)
#'
#' post_fun_pearson <- posterior_cor_updf(r_pearson, n)
#' post_fun_pearson_trans <- posterior_cor_updf(r_pearson_trans, n)
#' post_fun_kendall <- posterior_cor_updf(r_kendall, n, method = "kendall")
#' post_fun_kendall_trans <- posterior_cor_updf(r_kendall_trans, n, method = "kendall")
#'
#' # the square-root transformation affects inference for Pearson's correlation
#' # but not for Kendall's correlation
#' original_par <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' plot(
#'   x = grid, y = normalise(post_fun_pearson_trans(grid)),
#'   type = "l", col = "red",
#'   xlab = "Pearson correlation", ylab = "",
#'   bty = "n", yaxt = "n"
#' )
#' lines(x = grid, y = normalise(post_fun_pearson(grid)))
#' plot(
#'   x = grid, y = normalise(post_fun_kendall(grid)), type = "l",
#'   xlab = "Kendall correlation", ylab = "",
#'   bty = "n", yaxt = "n"
#' )
#' lines(x = grid, y = normalise(post_fun_kendall_trans(grid)), col = "red")
#' par(original_par)
#'
#' @export
posterior_cor_updf <- function(
    r,
    n,
    kappa = NULL,
    alternative = c("two.sided", "greater", "less"),
    method = c("pearson", "kendall"),
    n_bins = 1e3,
    max_iter = 1e7
) {

  method <- rlang::arg_match(method)
  alternative <- rlang::arg_match(alternative)
  kappa <- kappa %||% ifelse(
    test = method == "pearson",
    yes = 1, no = 2
  )

  params <- validate_posterior_cor_updf_input(r, n, kappa, n_bins, max_iter)
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

  cor_grid <- create_cor_grid(r, n, n_bins)

  if (method == "kendall") {
    check_kendall_normal_approx(r, n)
    d <- posterior_tau(r, n, cor_grid, kappa, alternative)
  } else {
    d <- posterior_rho(r, n, cor_grid, kappa, alternative, max_iter)
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

  base_fun <- stats::approxfun(
    x = cor_grid, y = d, method = "linear", yleft = 0, yright = 0, na.rm = TRUE
  )

  if (alternative == "two.sided") {
    result <- base_fun
  } else {
    if (alternative == "greater") {
      result <- function(x) {
        y <- base_fun(x)
        y[x < 0] <- 0
        return(y)
      }
    } else {
      result <- function(x) {
        y <- base_fun(x)
        y[x > 0] <- 0
        return(y)
      }
    }
  }

  return(result)

}


# PEARSON ---------------------------------------------------------------------

#' Posterior density for Pearson's correlation
#'
#' @description
#' Computes the normalised posterior density of the population correlation
#' coefficient \eqn{\rho}. The function combines the exact sampling
#' distribution of the sample correlation, the stretched beta prior on \eqn{\rho},
#' and a closed-form normalizing constant. If the sample correlation is
#' (practically) equal to zero, Jeffrey's approximation is used.
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param cor_grid Numeric vector. Grid of correlation values to evaluate.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param alternative Character string specifying the alternative hypothesis.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric vector of posterior density values.
#' @keywords internal
posterior_rho <- function(r, n, cor_grid, kappa, alternative, max_iter) {
  if (abs(r) < sqrt(.Machine$double.eps)) {
    result <- posterior_rho_jeffreys(r, n, cor_grid, kappa, alternative, max_iter)
  } else {
    result <- posterior_rho_exact(r, n, cor_grid, kappa, alternative, max_iter)
  }
  return(result)
}

#' @noRd
posterior_rho_exact <- function(r, n, cor_grid, kappa, alternative, max_iter) {
  result <- bf_rho_exact(r, n, kappa, max_iter) *
    likelihood_rho_exact(r, n, cor_grid, max_iter) *
    prior_cor(cor_grid, kappa, alternative, method = "pearson")
  return(result)
}

#' @noRd
posterior_rho_jeffreys <- function(r, n, cor_grid, kappa, alternative, max_iter) {
  result <- bf_rho_jeffreys(r, n, kappa, max_iter) *
    likelihood_rho_jeffreys(r, n, cor_grid) *
    prior_cor(cor_grid, kappa, alternative, method = "pearson")
  return(result)
}

#' @name bf_rho_exact
#' @title Normalizing constant for the exact posterior density of Pearson's
#' correlation
#'
#' @description
#' Computes the analytic normalising constant for the posterior density of the
#' correlation coefficient. This expression coincides with the Bayes factor
#' reported by Ly et al. (2018), because the marginal likelihood under the null
#' hypothesis equals one. Instead of the exact solution, Jeffrey's approximation
#' can be used, which is accurate when the sample correlation is (practically)
#' equal to zero.
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric scalar giving the normalising constant.
#' @keywords internal
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

#' @rdname bf_rho_exact
#' @keywords internal
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

#' @name likelihood_rho_exact
#' @title Likelihood of the sample Pearson's correlation
#'
#' @description
#' Evaluates the exact sampling distribution of the sample correlation
#' coefficient \eqn{r} given the population correlation \eqn{\rho}.
#' This is expressed in terms of generalized hypergeometric functions. Instead
#' of the exact solution, Jeffrey's approximation can be used, which is accurate
#' when the sample correlation is (practically) equal to zero.
#'
#' @param r Numeric value. The observed sample correlation coefficient.
#' @param n Numeric value. The sample size.
#' @param cor_grid Numeric vector. Grid of correlation values to evaluate.
#' @param max_iter Integer. Maximum number of iterations (attempts) to solve
#'        generalised hypergeometric functions.
#'
#' @return Numeric vector of likelihood values.
#' @keywords internal
likelihood_rho_exact <- function(r, n, cor_grid, max_iter) {

  hyper_term_even <- hypergeo::genhypergeo(
    U = c((n - 1) / 2, (n - 1) / 2),
    L = 1 / 2,
    z = (r * cor_grid)^2,
    maxiter = max_iter
  )

  likelihood_rho_even <- (1 - cor_grid^2)^((n - 1) / 2) *
    Re(hyper_term_even)

  hyper_term_odd <- hypergeo::genhypergeo(
    U = c(n / 2, n / 2),
    L = 3 / 2,
    z = (r * cor_grid)^2,
    maxiter = max_iter
  )

  log_term_odd <- 2 * (lgamma(n / 2) - lgamma((n - 1) / 2)) +
    ((n - 1) / 2) *
    log(1 - cor_grid^2)

  likelihood_rho_odd <- 2 * r * cor_grid *
    exp(log_term_odd) *
    Re(hyper_term_odd)

  result <- likelihood_rho_even + likelihood_rho_odd

  return(result)

}

#' @rdname likelihood_rho_exact
#' @keywords internal
likelihood_rho_jeffreys <- function(r, n, cor_grid) {
  numerator <- (1 - cor_grid^2)^((n - 1) / 2)
  denominator <- (1 - cor_grid * r)^(n - 1 - (1 / 2))
  result <- numerator / denominator
  return(result)
}


# KENDALL ---------------------------------------------------------------------

#' Compute posterior density for Kendall's tau
#'
#' Computes the posterior density for Kendall's rank correlation coefficient.
#' The method relies on the asymptotic normal approximation to the sampling
#' distribution of Kendall's tau. Numerical integration (trapezoidal rule) is
#' used to normalise the posterior density over a discretised grid.
#'
#' @param r Observed Kendall's tau correlation coefficient (scalar).
#' @param n Sample size.
#' @param cor_grid Numeric vector of grid values for tau over which the
#'   posterior density is evaluated.
#' @param kappa Concentration parameter of the prior distribution.
#' @param alternative Character string specifying the alternative hypothesis.
#' @return Numeric vector of posterior density values over `cor_grid`.
#' @keywords internal
posterior_tau <- function(r, n, cor_grid, kappa, alternative) {
  density <- posterior_tau_integrand(r, n, cor_grid, kappa, alternative)
  dx <- diff(cor_grid)
  normalising_constant <- sum(
    (density[-1] + density[-length(density)]) / 2 * dx
  )
  result <- density / normalising_constant
  return(result)
}

#' Posterior integrand for Kendall's tau
#'
#' Computes the unnormalised posterior density, i.e. the product of the
#' likelihood and prior for each grid value of tau.
#'
#' @inheritParams posterior_tau
#' @return Numeric vector of unnormalised posterior densities.
#' @noRd
posterior_tau_integrand <- function(r, n, cor_grid, kappa, alternative) {
  result <- likelihood_tau(r, n, cor_grid) *
    prior_cor(cor_grid, kappa, alternative, method = "kendall")
  return(result)
}

#' Likelihood function for Kendall's tau
#'
#' Evaluates the normal approximation likelihood of observing Kendall's
#' tau, expressed in terms of the standardised test statistic \eqn{T^*}.
#'
#' @param r Observed Kendall's tau correlation coefficient (scalar).
#' @param n Sample size.
#' @param cor_grid Numeric vector of candidate population tau values.
#' @param log_d Logical; if TRUE, return log-likelihood values.
#' @return Numeric vector of (log-)likelihood values.
#' @keywords internal
likelihood_tau <- function(r, n, cor_grid, log_d = FALSE) {
  t_star <- standardise_tau(r, n)
  mu_vec <- (3 / 2) * cor_grid * sqrt(n)
  result <- stats::dnorm(
    x = t_star,
    mean = mu_vec,
    sd = 1,
    log = log_d
  )
  return(result)
}

#' Standardised test statistic for Kendall's tau
#'
#' Computes the test statistic \eqn{T^*} used in the normal approximation
#' to the sampling distribution of Kendall's tau.
#'
#' @param r Observed Kendall's tau correlation coefficient (scalar).
#' @param n Sample size.
#' @return Numeric scalar, the standardised statistic \eqn{T^*}.
#' @noRd
standardise_tau <- function(r, n) {
  numerator <- r * ((n * (n - 1)) / 2)
  denominator <- sqrt(n * (n - 1) * (2 * n + 5) / 18)
  result <- numerator / denominator
  return(result)
}

#' Check adequacy of the normal approximation for Kendall's tau
#'
#' The posterior density is based on the asymptotic normal approximation to the
#' sampling distribution of Kendall's tau. The adequacy of this approximation
#' depends on both the sample size `n` and the magnitude of the observed
#' correlation coefficient `r`. Van Doorn et al. (2018) provide simulation
#' results indicating the smallest `n` required for different ranges of `abs(r)`
#' such that the Kolmogorov–Smirnov statistic does not exceed 0.038. This
#' function checks whether the observed combination of `r` and `n` meets these
#' criteria, and issues a warning if not.
#'
#' @param r Numeric scalar. Observed Kendall's tau.
#' @param n Integer scalar. Sample size.
#'
#' @details
#' The simulation data from Van Doorn et al. (2018) are available at:
#' \url{https://osf.io/download/7u7p2/}.
#'
#' @return Invisibly returns `NULL`. Called for its side effect of issuing a
#' warning if the (r, n) combination falls below the recommended thresholds.
#'
#' @keywords internal
check_kendall_normal_approx <- function(r, n) {
  tau_abs <- abs(r)
  min_n <- dplyr::case_when(
    tau_abs < 0.1 ~ 10L,
    tau_abs < 0.3 ~ 15L,
    tau_abs < 0.5 ~ 20L,
    tau_abs < 0.8 ~ 30L,
    TRUE          ~ 100L
  )
  if (n < min_n) {
    rlang::warn(
      message = paste0(
        "Normal approximation for Kendall's tau may be inaccurate: ",
        "observed |tau| = ", round(tau_abs, 3), ", n = ", n,
        " (requires n >= ", min_n, ")."
      )
    )
  }
  invisible(NULL)
}


# SHARED ----------------------------------------------------------------------

#' Calculate prior density for correlation coefficient
#'
#' @description
#' Calculates the stretched beta prior density for the correlation coefficient,
#' with optional truncation for directional alternatives.
#'
#' @param cor_grid Numeric vector. Grid of correlation values to evaluate.
#' @param kappa Numeric value. Prior concentration parameter.
#' @param alternative Character string specifying the alternative hypothesis.
#'
#' @return Numeric vector of prior density values.
#' @keywords internal
prior_cor <- function(cor_grid, kappa, alternative, method) {
  if (method == "kendall") {
    beta_term <- (2^(-2 / kappa)) / beta(a = 1 / kappa, b = 1 / kappa)
    cos_term <- cos((pi * cor_grid) / 2)^(2 / kappa - 1)
    base_density <- pi * beta_term * cos_term
  } else {
    base_density <- (1 / 2) *
      stats::dbeta(
        x = (cor_grid + 1) / 2,
        shape1 = 1 / kappa,
        shape2 = 1 / kappa
      )
  }
  result <- truncate_prior_density(
    cor_grid = cor_grid,
    base_density = base_density,
    alternative = alternative
  )
  return(result)
}

#' @noRd
truncate_prior_density <- function(cor_grid, base_density, alternative) {
  if (alternative == "two.sided") {
    return(base_density)
  }
  result <- numeric(length(cor_grid)) # initialise with zero density
  if (alternative == "greater") {
    idx <- cor_grid >= 0
  } else {
    idx <- cor_grid <= 0
  }
  # with prior being symmetric around zero, normalisation constant for
  # truncating at 0 is always exactly 1/2, hence multiply by 2.
  result[idx] <- 2 * base_density[idx]
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
#' @keywords internal
create_cor_grid <- function(r, n, n_bins) {

  n_bins <- max(n_bins, 100)

  max_abs_r <- 0.999
  grid_spacing <- 2 * max_abs_r / (round(n_bins * 0.7) - 1)
  n_steps <- floor(max_abs_r / grid_spacing)
  base_grid <- seq(
    from = -n_steps * grid_spacing,
    to = n_steps * grid_spacing,
    by = grid_spacing
  )

  n_peak <- n_bins - length(base_grid)
  if (n_peak > 0) {
    unit_grid <- stats::qlogis(
      p = seq(from = 0, to = 1, length.out = n_peak + 2)
    )
    unit_grid <- unit_grid[-c(1, length(unit_grid))]
    peak_grid <- tanh(
      atanh(r) + unit_grid / sqrt(n)
    )
    peak_grid <- pmax(-max_abs_r, pmin(max_abs_r, peak_grid))
  } else {
    peak_grid <- numeric(0)
  }

  result <- sort(unique(c(
    base_grid, peak_grid
  )))

  return(result)

}
