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
#' @param cor_use String indicating how to handle missing values in the
#'        Pearson correlation analysis, passed to the `use` argument of
#'        [stats::cor()]. Default is `"complete.obs"`.
#' @param kappa Numeric value controlling the "concentration" of the stretched
#'        beta prior on the correlation coefficient. Default is `1`, resulting
#'        in a uniform prior. See [posterior_rho_updf()] for details.
#' @param ... Additional arguments passed to [posterior_rho_updf()].
#'
#' @return A [tibble::tbl_df-class] with one row per MCMC sample containing:
#'   \item{r}{Pearson correlation coefficient for each MCMC sample}
#'   \item{n}{Number of complete observations used for correlation calculation}
#'   \item{posterior_density}{Function object for evaluating the unnormalized
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
    cor_use = "complete.obs",
    kappa = 1,
    ...
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

  result <- data %>%
    dplyr::group_by(.data[[draw_id]]) %>%
    dplyr::summarise(
      r = stats::cor(
        x = .data[[parameter]],
        y = .data[[covariate]],
        use = cor_use,
        method = "pearson"
      ),
      n = sum(!is.na(.data[[parameter]]) & !is.na(.data[[covariate]])),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      .data[["n"]] >= 3 & is.finite(.data[["r"]]) & abs(.data[["r"]]) <= 1
    )

  if (nrow(result) == 0) {
    rlang::abort(
      message = "No MCMC samples with valid correlation values."
    )
  }

  result <- result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      posterior_updf = list(
        posterior_rho_updf(
          r = .data[["r"]],
          n = .data[["n"]],
          kappa = kappa,
          ...
        )
      )
    ) %>%
    dplyr::ungroup()

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
  invalid_cols <- vapply(
    X = col_names,
    FUN = function(col_name) {
      length(col_name) != 1 || !is.character(col_name)
    },
    FUN.VALUE = logical(1)
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

#' Summarise plausible correlation estimates
#'
#' @description
#' Takes the output of [run_plausible_cor()] and summarises the distribution of
#' plausible correlation values at both the sample level (across MCMC samples)
#' and the population level (via the mean posterior density across MCMC
#' samples).
#'
#' @param x A data frame output from [run_plausible_cor()].
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
    x,
    point_method = c("mean", "median"),
    interval_width = c(0.5, 0.8, 0.95),
    interval_method = c("hdci", "qi")
) {

  point_method <- rlang::arg_match(arg = point_method)
  interval_method <- rlang::arg_match(arg = interval_method)

  sample_summary <- x %>%
    ggdist::point_interval(
      .data[["r"]],
      .width = interval_width,
      .point = eval(rlang::sym(point_method)),
      .interval = eval(rlang::sym(interval_method))
    ) %>%
    dplyr::select(
      -dplyr::all_of(c(".point", ".interval"))
    ) %>%
    dplyr::rename_with(
      .fn = ~ gsub("^\\.", "", .x),
      .cols = dplyr::starts_with(".")
    ) %>%
    dplyr::rename_with(
      .fn = ~ point_method,
      .cols = "r"
    ) %>%
    dplyr::mutate(
      type = "sample",
      p_dir = max(
        x[["r"]] > 0,
        x[["r"]] < 0
      )
    ) %>%
    dplyr::relocate(dplyr::all_of("type"))

  mean_posterior_rho <- x %>%
    get_posterior_rho_densities() %>%
    get_mean_posterior_rho()

  rho_grid <- mean_posterior_rho[["x"]]
  grid_spacing <- diff(rho_grid)[1]
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
        sum(mean_density[rho_grid > 0]) * grid_spacing,
        sum(mean_density[rho_grid < 0]) * grid_spacing
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

#' Evaluate posterior density functions over a grid
#'
#' @description
#' Internal helper that evaluates each posterior density function in the output
#' of [run_plausible_cor()] over a shared grid from -1 to 1.
#'
#' @param x Data frame with columns `posterior_updf` of function objects and
#'        `r` containing the sample correlation coefficients.
#' @param grid_spacing The step size for the grid of correlation values to
#'        evaluate.
#'
#' @return A data frame where each row corresponds to one density value for a
#'         given posterior density function.
#' @noRd
get_posterior_rho_densities <- function(
    x,
    grid_spacing = 1e-3
) {

  rho_grid <- seq(from = -1, to = 1, by = grid_spacing)

  posterior_rho_densities_list <- Map(
    f = function(updf, r_val, grid = rho_grid) {
      if (!rlang::is_function(updf)) {
        return(NULL)
      }
      dens <- updf(grid)
      dens_norm <- dens / sum(dens * diff(grid)[1])
      result <- tibble::tibble(
        !!rlang::sym(as.character(r_val)) := dens_norm
      )
      return(result)
    },
    x[["posterior_updf"]],
    x[["r"]]
  )

  result <- posterior_rho_densities_list %>%
    dplyr::bind_cols() %>%
    dplyr::mutate(x = rho_grid) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of("x"),
      names_to = "r",
      values_to = "density"
    ) %>%
    dplyr::mutate(
      r = as.numeric(.data[["r"]])
    ) %>%
    dplyr::select(
      dplyr::all_of(c("r", "x", "density"))
    )

  return(result)

}

#' Compute mean posterior correlation density
#'
#' @description
#' Internal helper that computes the mean posterior density across MCMC samples,
#' based on the evaluated densities on a shared grid.
#'
#' @param posterior_rho_densities Data frame output from
#'        [get_posterior_rho_densities()].
#'
#' @return A tibble with columns `x` (grid points) and `density`
#'         (mean posterior density).
#' @noRd
get_mean_posterior_rho <- function(
    posterior_rho_densities
) {
  result <- posterior_rho_densities %>%
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
#' @description Internal helper function to compute a point estimate (mean or
#' median) from a numeric vector of values and a corresponding density over a
#' discretised grid.
#'
#' @param val Numeric vector representing the grid of values over which the
#'        density is defined. Assumed to be equally spaced.
#' @param dens Numeric vector of the same length as `val`, representing the
#'        density at each grid point.
#' @param method Character string specifying which point estimate to compute.
#'        One of `"mean"` (default) or `"median"`.
#'
#' @return A numeric value corresponding to the requested point estimate.
#'
#' @noRd
get_point_estimate <- function(
    val,
    dens,
    method = c("mean", "median")
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
    median_approx <- stats::approx(
      x = cdf,
      y = val,
      xout = 0.5,
      ties = "ordered"
    )
    result <- median_approx[["y"]]
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

  intervals <- lapply(
    X = width,
    FUN = function(w) {
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
  )

  result <- dplyr::bind_rows(intervals)
  return(result)

}
