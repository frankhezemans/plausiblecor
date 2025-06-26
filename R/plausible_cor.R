#' Calculate plausible correlation values between model parameter and covariate
#'
#' @description
#' Implements a Bayesian approach for exploring individual differences in a
#' model parameter that was estimated without reference to the observed
#' covariate of interest. This function computes the Pearson correlation
#' coefficient between subject-wise parameter values (estimated using MCMC
#' sampling methods) and an observed covariate (e.g., questionnaire) for each
#' MCMC sample, resulting in a distribution of "plausible values" for the
#' correlation coefficient. It additionally provides the posterior density
#' function for each correlation coefficient, to enable inference about the
#' latent correlation in the population (as opposed to the given sample of
#' subjects).
#'
#' @param mcmc_data Data frame of MCMC samples in long format, where each row
#'        represents a unique combination of MCMC draw, subject, and
#'        parameter value.
#' @param covariate_data Optional data frame containing covariate data (if not
#'        in `mcmc_data`). Should contain one row per participant.
#' @param draw_id String denoting the column with MCMC sample IDs.
#' @param subject_id String denoting the column with subject IDs.
#' @param parameter String denoting the column with the estimated model
#'        parameter.
#' @param covariate String denoting the column with the observed covariate.
#' @param cor_method String indicating the method for computing the sample
#'        correlation coefficient. Currently ignored (fixed to "pearson").
#' @param posterior_args Named list of arguments passed to
#'        [posterior_rho_updf()], which controls how the unnormalised posterior
#'        density function is computed for each Pearson correlation coefficient.
#'        The following arguments can be included:
#'        * `kappa` Numeric value controlling the "concentration" of the
#'        stretched beta prior on the correlation coefficient. Default is `1`,
#'        resulting in a uniform prior.
#'        * `n_bins` Integer denoting the number of bins ("resolution") used
#'        for the discretised grid of correlation values, for which the
#'        unnormalised posterior density values are computed. Default is `1e3`.
#'        * `max_iter` Integer denoting the maximum number of iterations
#'        (attempts) to obtain a posterior density function for a given MCMC
#'        sample. Default is `1e7`.
#'        * `...` additional arguments passed forward to [stats::approxfun()].
#'
#' @return A [tibble::tbl_df-class] with one row per MCMC sample containing:
#'   \item{.draw}{The MCMC sample ID}
#'   \item{r}{Pearson correlation coefficient}
#'   \item{n}{Number of complete observations used for correlation calculation}
#'   \item{posterior_updf}{Function object for evaluating the unnormalised
#'         posterior density at any correlation value between -1 and 1}
#'
#' @details
#' This function implements the approach described by Ly et al. (2017) for
#' exploring individual differences in cognitive-model-based neuroscience:
#'
#' 1. First, it computes the Pearson correlation between an estimated model
#'    parameters and observed covariate separately for each MCMC sample,
#'    resulting in a set of plausible correlation values for the given sample
#'    of subjects.
#'    This handles uncertainty in the estimation of individual subjects'
#'    model parameters, but does not account for uncertainty in generalising
#'    from the sample to the population.
#'
#' 2. Second, it calculates the posterior density function for each correlation
#'    using the analytical solution derived by Ly et al. (2018), resulting in
#'    a set of plausible posterior distributions of the correlation.
#'    This handles uncertainty in generalizing from the sample to the
#'    population, and is therefore suitable for inferences about the
#'    correlation coefficient in the population.
#'
#' The function requires at least 3 valid observations per MCMC sample and
#' filters out any correlation values that are not finite or outside the valid
#' range of [-1, 1].
#'
#' @references
#' Ly, A., Boehm, U., Heathcote, A., Turner, B. M., Forstmann, B., Marsman, M.,
#' & Matzke, D. (2017). A flexible and efficient hierarchical bayesian approach
#' to the exploration of individual differences in cognitive‐model‐based
#' neuroscience. *Computational models of brain and behavior*, 467-479.
#' https://doi.org/10.1002/9781119159193.ch34
#'
#' Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Analytic posteriors for
#' Pearson's correlation coefficient. *Statistica Neerlandica*, 72, 4–13.
#' https://doi.org/10.1111/stan.12111
#'
#' @export
run_plausible_cor <- function(
    mcmc_data,
    covariate_data = NULL,
    draw_id,
    subject_id,
    parameter,
    covariate,
    cor_method = "pearson",
    posterior_args = NULL
) {

  if (cor_method != "pearson") {
    rlang::abort(message = "'cor_method' must be 'pearson'")
  }

  data <- prep_plausible_cor_data(
    mcmc_data = mcmc_data,
    covariate_data = covariate_data,
    draw_id = draw_id,
    subject_id = subject_id,
    parameter = parameter,
    covariate = covariate
  )

  posterior_args <- validate_posterior_args(posterior_args)

  result <- run_plausible_cor_pearson(
    data = data,
    parameter = parameter,
    covariate = covariate,
    posterior_args = posterior_args
  )

  return(result)

}

#' @noRd
run_plausible_cor_pearson <- function(
    data,
    parameter,
    covariate,
    posterior_args
) {

  result <- data %>%
    dplyr::group_by(.data[[".draw"]]) %>%
    dplyr::summarise(
      r = stats::cor(
        x = .data[[parameter]],
        y = .data[[covariate]],
        use = "everything",
        method = "pearson"
      ),
      n = sum(!is.na(.data[[parameter]]) & !is.na(.data[[covariate]])),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      .data[["n"]] >= 3 & is.finite(.data[["r"]]) & abs(.data[["r"]]) <= 1
    )

  if (nrow(result) == 0) {
    rlang::abort(message = "No MCMC samples with valid correlation values.")
  }

  result <- result %>%
    dplyr::mutate(
      posterior_updf = purrr::map2(
        .x = .data[["r"]],
        .y = .data[["n"]],
        .f = function(r_val, n_val) {
          rlang::exec(
            .fn = posterior_rho_updf,
            r = r_val,
            n = n_val,
            !!!posterior_args
          )
        }
      )
    )

  return(result)

}

#' Summarise plausible correlation estimates
#'
#' @description
#' Takes the output of [run_plausible_cor()] and summarises the distribution of
#' plausible correlation values at both the sample level (across MCMC samples)
#' and the population level (via the mean posterior density across MCMC
#' samples).
#'
#' @param .data A data frame output from [run_plausible_cor()].
#' @param point_method Method used to compute the central tendency of the
#'        correlation estimate. One of `"mean"` or `"median"`. Default is
#'        `"mean"`.
#' @param interval_width Numeric vector of desired credible interval widths
#'        (e.g., `c(0.5, 0.8, 0.95)`). Must be values between 0 and 1.
#' @param interval_method Method used for computing intervals. One of `"hdci"`
#'        (highest-density continuous interval, a.k.a. shortest probability
#'        interval) or `"qi"` (quantile interval). Default is `"hdci"`.
#'
#' @return A tibble with one row per summary type ("sample" and "posterior"),
#'         and credible interval width, containing:
#'   \item{type}{Whether the row corresponds to the sample-based or
#'         posterior-based summary}
#'   \item{mean / median}{Point estimate of the correlation coefficient}
#'   \item{lower / upper}{Credible interval bounds for each `interval_width`}
#'   \item{p_dir}{Directional probability (proportion of mass > 0 or < 0)}
#'
#' @export
summarise_plausible_cor <- function(
    .data,
    point_method = c("mean", "median"),
    interval_width = c(0.5, 0.8, 0.95),
    interval_method = c("hdci", "qi")
) {

  validate_column_inputs(
    col_names = c(".draw", "r", "posterior_updf"),
    data_frame = .data,
    data_name = ".data"
  )

  point_method <- rlang::arg_match(arg = point_method)
  interval_method <- rlang::arg_match(arg = interval_method)

  sample_summary <- .data %>%
    dplyr::select(dplyr::all_of("r")) %>%
    summarise_samples(
      point_method = point_method,
      interval_width = interval_width,
      interval_method = interval_method
    ) %>%
    dplyr::mutate(
      type = "sample"
    ) %>%
    dplyr::relocate(dplyr::all_of("type"))

  mean_posterior_rho <- .data %>%
    get_posterior_rho_densities() %>%
    get_mean_posterior_rho()

  rho_grid <- mean_posterior_rho[["x"]]
  dx <- diff(rho_grid)[1]
  mean_density <- mean_posterior_rho[["mean_density"]]

  population_point <- get_point_estimate(
    val = rho_grid,
    dens = mean_density,
    method = point_method
  )

  population_interval <- get_interval(
    val = rho_grid,
    dens = mean_density,
    width = interval_width,
    method = interval_method
  )

  population_summary <- population_interval %>%
    dplyr::mutate(
      type = "population",
      !!point_method := population_point,
      p_dir = max(
        sum(mean_density[rho_grid > 0]) * dx,
        sum(mean_density[rho_grid < 0]) * dx
      )
    ) %>%
    dplyr::relocate(dplyr::all_of("type"))

  result <- dplyr::bind_rows(
    sample_summary, population_summary
  ) %>%
    dplyr::mutate(
      type = factor(.data[["type"]])
    )

  return(result)

}

#' Compare two distributions of plausible correlation values
#'
#' @description
#' Compares the posterior distributions of plausible correlation values, as
#' returned by two separate calls to [run_plausible_cors()]. Uses two methods
#' methods for defining the distribution of the difference: a **sample-based**
#' comparison, which subtracts the correlation values from matching MCMC
#' samples, and a **population-based** comparison, which draws quantiles from
#' each posterior distribution (accounting for population uncertainty) before
#' computing the difference.
#'
#' @param x,y Data frames returned by [run_plausible_cors()], each containing
#'   the columns `.draw`, `r`, and `posterior_updf`.
#' @param point_method Character string specifying the summary point estimate
#'   of the delta distribution. Either `"mean"` (default) or `"median"`.
#' @param interval_width Numeric vector of interval widths (between 0 and 1)
#'   for the uncertainty intervals to be returned. Default is
#'   `c(0.5, 0.8, 0.95)`.
#' @param interval_method Character string specifying the method used to
#'   compute intervals. Either `"hdci"` (highest-density continuous interval;
#'   default) or `"qi"` (quantile interval).
#' @param n_samples Integer indicating how many quantiles to draw per posterior
#'   distribution when `comparison_method = "population"`. Default is 1.
#' @param rng_seed Optional numeric vector of length 2 to control the random
#'   seed for quantile sampling of `x` and `y`, respectively. If `NULL`,
#'   random seeds will not be set.
#'
#' @return A [tibble::tibble()] containing summary statistics of the difference
#'  between correlation coefficients.
#'
#' @details
#' For the sample-based comparison, the function computes the difference in the
#' Pearson correlation coefficient for each `.draw` shared between the two data
#' frames.
#'
#' For the population-based comparison, the function draws quantile value(s)
#' from the posterior distribution of each Pearson correlation (as defined by
#' `posterior_updf`) and computes the difference. This allows for comparison at
#' the population level rather than the specific sample of subjects.
#'
#' This function was adapted from code previously released with the Dynamic
#' Models of Choice toolbox (Heathcote et al., 2019).
#'
#' @references
#' Ly, A., Boehm, U., Heathcote, A., Turner, B. M., Forstmann, B., Marsman, M.,
#' & Matzke, D. (2017). A flexible and efficient hierarchical bayesian approach
#' to the exploration of individual differences in cognitive‐model‐based
#' neuroscience. *Computational models of brain and behavior*, 467-479.
#' https://doi.org/10.1002/9781119159193.ch34
#'
#' Heathcote, A., Lin, Y.S., Reynolds, A., Strickland, L., Gretton, M., &
#' Matzke, D. (2019). Dynamic models of choice. *Behavior Research Methods*,
#' 51, 961-985. https://doi.org/10.3758/s13428-018-1067-y
#'
#' @export
compare_plausible_cors <- function(
    x,
    y,
    point_method = c("mean", "median"),
    interval_width = c(0.5, 0.8, 0.95),
    interval_method = c("hdci", "qi"),
    n_samples = 1,
    rng_seed = NULL
) {

  validate_column_inputs(
    col_names = c(".draw", "r", "posterior_updf"),
    data_frame = x,
    data_name = "x"
  )
  validate_column_inputs(
    col_names = c(".draw", "r", "posterior_updf"),
    data_frame = y,
    data_name = "y"
  )

  if (nrow(x) != nrow(y)) {
    rlang::abort(
      message = "Data frames 'x' and 'y' should have the same number of rows."
    )
  }

  if (is.null(rng_seed) || any(!is.finite(rng_seed))) {
    rng_seed <- c(NA, NA)
  }
  if (length(rng_seed) != 2) {
    rlang::abort(
      message = "'rng_seed' should contain two elements (one per data frame)."
    )
  }

  point_method <- rlang::arg_match(arg = point_method)
  interval_method <- rlang::arg_match(arg = interval_method)

  sample_delta_summary <- dplyr::inner_join(
    x = dplyr::select(
      .data = x,
      dplyr::all_of(c(".draw", "r"))
    ),
    y = dplyr::select(
      .data = y,
      dplyr::all_of(c(".draw", "r"))
    ),
    by = ".draw",
    suffix = c("_x", "_y")
  ) %>%
    dplyr::mutate(
      delta = .data[["r_x"]] - .data[["r_y"]],
      .keep = "none"
    ) %>%
    summarise_samples(
      point_method = point_method,
      interval_width = interval_width,
      interval_method = interval_method
    ) %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::relocate(dplyr::all_of("type"))

  population_delta_summary <- purrr::map2(
    .x = list(x = x, y = y),
    .y = stats::setNames(object = rng_seed, nm = c("x", "y")),
    .f = function(data, seed, n = n_samples) {
      get_sampled_quantiles(
        .data = data,
        n_samples = n,
        starter_seed = seed
      )
    }
  ) %>%
    dplyr::bind_rows(
      .id = "dataset"
    ) %>%
    dplyr::select(
      dplyr::all_of(c("dataset", ".draw", "quantile"))
    ) %>%
    tidyr::pivot_wider(
      id_cols = tidyr::all_of(".draw"),
      names_from = "dataset",
      values_from = "quantile"
    ) %>%
    dplyr::mutate(
      delta = .data[["x"]] - .data[["y"]],
      .keep = "none"
    ) %>%
    summarise_samples(
      point_method = point_method,
      interval_width = interval_width,
      interval_method = interval_method
    ) %>%
    dplyr::mutate(type = "population") %>%
    dplyr::relocate(dplyr::all_of("type"))

  result <- dplyr::bind_rows(
    sample_delta_summary, population_delta_summary
  ) %>%
    dplyr::mutate(
      type = factor(.data[["type"]])
    )

  return(result)

}

#' Evaluate posterior density functions over a grid
#'
#' @description
#' Helper function that evaluates each posterior density function in the output
#' of [run_plausible_cor()] over a shared grid from -1 to 1. Typical users
#' would not have to call this function, instead relying on
#' [summarise_plausible_cor()], [compare_plausible_cors()], and
#' [plot_population_cor()].
#'
#' @param .data Data frame output from [run_plausible_cor()].
#' @param grid_spacing The step size for the grid of correlation values to
#'        evaluate.
#'
#' @return A data frame where each row corresponds to one density value for a
#'         given posterior density function.
#'
#' @export
get_posterior_rho_densities <- function(.data, grid_spacing = 1e-3) {

  validate_column_inputs(
    col_names = c(".draw", "r", "posterior_updf"),
    data_frame = .data,
    data_name = ".data"
  )

  rho_grid <- seq(from = -0.999, to = 0.999, by = grid_spacing)

  single_density_grid <- function(
    draw_id, r_val, updf, grid = rho_grid, dx = grid_spacing
  ) {
    if (!rlang::is_function(updf)) {
      return(NULL)
    }
    dens <- updf(grid)
    if (any(!is.finite(dens))) {
      return(NULL)
    }
    dens_norm <- dens / sum(dens * dx)
    result <- tibble::tibble(
      .draw = draw_id,
      r = r_val,
      x = grid,
      density = dens_norm
    )
    return(result)
  }

  result <- purrr::pmap(
    .l = .data,
    .f = function(.draw, r, posterior_updf, ...) {
      single_density_grid(
        draw_id = .draw,
        r_val = r,
        updf = posterior_updf
      )
    }
  ) %>%
    dplyr::bind_rows()

  return(result)

}

#' Compute mean posterior correlation density
#'
#' @description
#' Helper function that computes the mean posterior density across MCMC samples,
#' based on the evaluated densities on a shared grid. Typical users would not
#' have to call this function, instead relying on [summarise_plausible_cor()],
#' [compare_plausible_cors()], and [plot_population_cor()].
#'
#' @param .data Data frame output from [get_posterior_rho_densities()].
#'
#' @return A tibble with columns `x` (grid points) and `density`
#'         (mean posterior density).
#'
#' @export
get_mean_posterior_rho <- function(.data) {

  validate_column_inputs(
    col_names = c(".draw", "r", "x", "density"),
    data_frame = .data,
    data_name = ".data"
  )

  result <- .data %>%
    dplyr::group_by(.data[["x"]]) %>%
    dplyr::summarise(
      mean_density = mean(.data[["density"]]),
      .groups = "drop"
    )
  dx <- diff(result[["x"]])[1]

  result <- result %>%
    dplyr::mutate(
      mean_density = .data[["mean_density"]] /
        sum(.data[["mean_density"]] * dx)
    )

  return(result)

}

#' @title Compute Point Estimate from Discretised Density
#'
#' @description Internal helper function to compute a point estimate (mean,
#' median, or arbitrary quantile) from a numeric vector of values and a
#' corresponding density over a discretised grid.
#'
#' @param val Numeric vector representing the grid of values over which the
#'        density is defined. Assumed to be equally spaced.
#' @param dens Numeric vector of the same length as `val`, representing the
#'        density at each grid point.
#' @param method Character string specifying which point estimate to compute.
#'        One of `"mean"`, `"median"`, or `"quantile"`.
#' @param prob Numeric value between 0 and 1, only required when
#'        `method = "quantile"`. Specifies the cumulative probability at which
#'        to evaluate the quantile.
#'
#' @return A numeric value corresponding to the requested point estimate.
#'
#' @noRd
get_point_estimate <- function(
    val,
    dens,
    method = c("mean", "median", "quantile"),
    prob = NULL
) {

  if (length(val) != length(dens)) {
    rlang::abort(
      message = "Input arguments 'val' and 'dens' must be of equal length."
    )
  }
  method <- rlang::arg_match(arg = method)
  dx <- diff(val)[1]

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
#'        density is defined. Assumed to be equally spaced.
#' @param dens Numeric vector of the same length as `val`, representing the
#'        density at each grid point.
#' @param width Numeric vector of desired interval widths
#'        (e.g., `0.95` for a 95% interval).
#' @param method Character string specifying which interval type to compute.
#'        One of `"qi"` (quantile interval) or `"hdci"`
#'        (highest-density continuous interval, default).
#'
#' @return A data frame with columns `"lower"`, `"upper"`, and `"width"` for
#'         each requested interval.
#'
#' @noRd
get_interval <- function(
    val,
    dens,
    width,
    method = c("hdci", "qi")
) {

  if (length(val) != length(dens)) {
    rlang::abort(
      message = "Input arguments 'val' and 'dens' must be of equal length."
    )
  }
  method <- rlang::arg_match(arg = method)
  dx <- diff(val)[1]
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

#' Sample quantiles from posterior correlation densities
#'
#' @description
#' Generates random quantiles from the posterior density function of a
#' correlation coefficient for each MCMC draw. This enables approximate
#' sampling from the full posterior distribution over correlations for
#' population-level inference.
#'
#' @param .data A data frame output by [run_plausible_cors()], containing one
#'   row per MCMC draw with columns `.draw`, `r`, and `posterior_updf`.
#' @param n_samples Integer specifying the number of quantiles to sample per
#'   MCMC draw. Defaults to `1`.
#' @param starter_seed Optional scalar numeric used to seed the random number
#'   generator for reproducibility. If `NA` (default), no seed is set.
#'
#' @return A [tibble::tibble] with one row per sampled quantile and columns:
#'   \item{.draw}{The MCMC draw ID}
#'   \item{r}{The original correlation estimate for the draw}
#'   \item{quantile}{A randomly sampled quantile from the posterior density}
#'
#' @details
#' The posterior density function is first approximated using
#' [get_posterior_rho_densities()], and quantiles are then sampled using
#' inverse transform sampling. If a seed is provided in `starter_seed`, it is
#' used to generate different seeds for each MCMC draw, to ensure independence
#' in the quantiles drawn for each posterior distribution, while retaining
#' reproducibility.
#'
#' @noRd
get_sampled_quantiles <- function(
    .data,
    n_samples = 1,
    starter_seed = NA
) {

  validate_column_inputs(
    col_names = c(".draw", "r", "posterior_updf"),
    data_frame = .data,
    data_name = ".data"
  )

  draw_quantiles <- function(density_grid, n_quantiles, seed = NA) {
    if (!is.na(seed)) {
      probs <- withr::with_seed(
        seed = seed,
        code = stats::runif(n = n_quantiles, min = 0, max = 1)
      )
    } else {
      probs <- stats::runif(n = n_quantiles, min = 0, max = 1)
    }
    quantiles <- get_point_estimate(
      val = density_grid[["x"]],
      dens = density_grid[["density"]],
      method = "quantile",
      prob = probs
    )
    result <- tibble::tibble(quantile = quantiles)
    return(result)
  }

  seeds <- NA
  if (!is.na(starter_seed)) {
    seeds <- withr::with_seed(
      seed = starter_seed,
      code = sample.int(
        n = nrow(.data) * 1e3,
        size = nrow(.data),
        replace = FALSE
      )
    )
  }

  result <- .data %>%
    get_posterior_rho_densities() %>%
    tidyr::nest(
      density_grid = tidyr::all_of(c("x", "density")),
      .by = tidyr::all_of(c(".draw", "r"))
    ) %>%
    dplyr::mutate(
      seed = seeds
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      samples = list(
        draw_quantiles(
          density_grid = .data[["density_grid"]],
          n_quantiles = n_samples,
          seed = .data[["seed"]]
        )
      ),
      .keep = "unused"
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(tidyr::all_of("samples"))

  return(result)

}

#' Summarise posterior samples
#'
#' @description
#' Summarises a single numeric column from a data frame containing posterior
#' samples (or any numeric data), providing a point estimate (mean or median),
#' credible intervals (e.g., highest-density continous intervals or quantiles),
#' and the directionality of the distribution (positive or negative).
#'
#' @param .data A data frame containing a single numeric column representing
#'   posterior samples or any distribution of interest.
#' @param ... Additional arguments passed to [ggdist::point_interval()].
#' @param point_method A string indicating the point estimate to compute. Can
#'   either be `"mean"` (default) or `"median"`.
#' @param interval_width Numeric vector of values between 0 and 1 specifying the
#'   widths of the credible intervals to calculate. Default is
#'   `c(0.5, 0.8, 0.95)`.
#' @param interval_method A string indicating the method for computing intervals.
#'   Can either be `"hdci"` (highest-density continous intervals; default) or
#'   `"qi"` (quantile intervals).
#'
#' @return A [tibble::tibble] summarizing the posterior samples.
#'
#' @noRd
summarise_samples <- function(
    .data,
    ...,
    point_method = c("mean", "median"),
    interval_width = c(0.5, 0.8, 0.95),
    interval_method = c("hdci", "qi")
) {

  data <- .data %>%
    dplyr::select(dplyr::where(is.numeric))

  varname <- colnames(data)
  if (length(varname) != 1) {
    rlang::abort(
      message = "Input data frame should only contain a single numeric column."
    )
  }

  p_dir <- max(
    mean(data[[varname]] > 0),
    mean(data[[varname]] < 0)
  )

  point_fun <- switch(
    point_method,
    "mean" = mean,
    "median" = stats::median
  )
  interval_fun <- switch(
    interval_method,
    "hdci" = ggdist::hdci,
    "qi" = ggdist::qi
  )

  result <- data %>%
    ggdist::point_interval(
      .width = interval_width,
      .point = point_fun,
      .interval = interval_fun,
      ...
    ) %>%
    dplyr::select(
      -dplyr::all_of(c(".point", ".interval"))
    ) %>%
    dplyr::rename_with(
      .fn = ~ gsub("^\\.", "", .x),
      .cols = dplyr::starts_with(".")
    ) %>%
    dplyr::rename(
      !!point_method := !!rlang::sym(varname)
    ) %>%
    dplyr::mutate(
      p_dir = p_dir
    )

  return(result)

}

#' Validate column name inputs for functions working with data frames
#'
#' @description
#' Internal helper function that checks whether user-supplied column names are
#' valid. Optionally verifies that the named columns exist in a provided data
#' frame.
#'
#' @param col_names Named character vector of column names to validate.
#' @param data_frame Optional data frame to check column existence for.
#' @param data_name Optional name of the data frame for error messages.
#'
#' @return Invisibly returns NULL if validation passes, otherwise raises an
#'         error.
#' @noRd
validate_column_inputs <- function(
    col_names,
    data_frame = NULL,
    data_name = NULL
) {
  invalid_cols <- purrr::map_lgl(
    .x = col_names,
    .f = function(col_name) {
      length(col_name) != 1 || !is.character(col_name)
    }
  )
  if (any(invalid_cols)) {
    bad_names <- names(col_names)[invalid_cols]
    rlang::abort(
      message = paste0(
        "Input '", paste(bad_names, collapse = "', '"),
        "' must be a single string denoting a column."
      )
    )
  }

  if (!is.null(data_frame)) {
    if (is.null(data_name)) {
      data_name <- rlang::as_label(rlang::enquo(data_frame))
    }
    missing_cols <- setdiff(unname(col_names), names(data_frame))
    if (length(missing_cols) > 0) {
      rlang::abort(
        message = paste0(
          "Missing required columns in '", data_name, "': ",
          paste(missing_cols, collapse = ", "), "."
        )
      )
    }
  }

  return(invisible(x = NULL))

}


#' Validate posterior_args parameter
#'
#' @description
#' Internal helper function to validate the structure and contents of the
#' posterior_args parameter passed to functions that use [posterior_rho_updf()].
#' Allowed arguments are: kappa, n_bins, max_iter, and [stats::approxfun()]
#' arguments (method, yleft, yright, rule, f, ties).
#'
#' @param posterior_args A named list of arguments or NULL
#'
#' @return The validated (and potentially cleaned) posterior_args list
#'
#' @noRd
validate_posterior_args <- function(posterior_args) {

  if (is.null(posterior_args)) {
    return(NULL)
  }
  if (!is.list(posterior_args)) {
    rlang::abort(message = "'posterior_args' must be a named list or NULL")
  }
  if (is.null(names(posterior_args)) && length(posterior_args) > 0) {
    rlang::abort(message = "'posterior_args' must be a named list")
  }

  arg_names <- names(posterior_args)
  allowed_args <- c(
    # posterior_rho_updf args
    "kappa", "n_bins", "max_iter",
    # stats::approxfun args
    "method", "yleft", "yright", "rule", "f", "ties"
  )
  disallowed_args <- setdiff(arg_names, allowed_args)
  if ("na.rm" %in% arg_names) {
    rlang::warn(message = "'na.rm' argument is ignored; forced to be TRUE")
    posterior_args[["na.rm"]] <- NULL
    disallowed_args <- setdiff(disallowed_args, "na.rm")
  }
  if (length(disallowed_args) > 0) {
    rlang::warn(
      message = paste0(
        "Removing unknown argument(s) from 'posterior_args': ",
        paste(disallowed_args, collapse = ", "),
        ".\nAllowed arguments are: ", paste(allowed_args, collapse = ", ")
      )
    )
    for (arg in disallowed_args) {
      posterior_args[[arg]] <- NULL
    }
    arg_names <- names(posterior_args)
  }

  return(posterior_args)

}

#' @noRd
prep_plausible_cor_data <- function(
    mcmc_data,
    covariate_data,
    draw_id,
    subject_id,
    parameter,
    covariate
) {

  validate_column_inputs(
    col_names = c(
      draw_id = draw_id,
      subject_id = subject_id,
      parameter = parameter,
      covariate = covariate
    )
  )

  validate_column_inputs(
    col_names = c(draw_id, subject_id, parameter),
    data_frame = mcmc_data,
    data_name = "mcmc_data"
  )

  if (covariate %in% names(mcmc_data)) {

    data <- mcmc_data

  } else if (!is.null(covariate_data)) {

    validate_column_inputs(
      col_names = c(subject_id, covariate),
      data_frame = covariate_data,
      data_name = "covariate_data"
    )

    covariate_data_clean <- covariate_data %>%
      dplyr::select(dplyr::all_of(c(subject_id, covariate))) %>%
      dplyr::distinct()

    n_subs <- length(unique(covariate_data_clean[[subject_id]]))
    if (nrow(covariate_data_clean) != n_subs) {
      rlang::abort(
        message = "Covariate data should have one unique value per subject."
      )
    }

    data <- dplyr::left_join(
      x = mcmc_data,
      y = covariate_data_clean,
      by = subject_id
    )

  } else {
    rlang::abort(
      message = paste0(
        "Covariate column not found in 'mcmc_data', ",
        "and no 'covariate_data' provided."
      )
    )
  }

  data <- data %>%
    dplyr::rename(dplyr::all_of(c(.draw = draw_id))) %>%
    dplyr::select(dplyr::all_of(c(".draw", parameter, covariate)))

  if (anyNA(data)) {
    rlang::warn(message = paste0(
      "Removing any rows with missing values in ",
      paste(draw_id, parameter, covariate, collapse = ", ")
    ))
    data <- tidyr::drop_na(data)
  }

  data <- data %>%
    dplyr::arrange(.data[[".draw"]])

  return(data)

}
