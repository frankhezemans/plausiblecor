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
#'        Alternatively, an object of class `"emc"` can be provided, which corresponds
#'        to a model fit using the `EMC2` package ([EMC2::EMC2()]). In this case,
#'        [EMC2::parameters()] is used to extract subject-wise MCMC samples of the
#'        model parameters.
#' @param covariate_data Optional data frame containing covariate data (if not
#'        already in `mcmc_data`). Should contain one row per participant.
#' @param draw_id String denoting the column with MCMC sample IDs in `mcmc_data`.
#'        If `mcmc_data` is an object of class `"emc"`, the provided `draw_id` is
#'        used to create an MCMC sample ID column, since the output of
#'        [EMC2::parameters()] does not include such a column by default.
#' @param subject_id String denoting the column with subject IDs. If both `mcmc_data`
#'        and `covariate_data` are provided, this column name should be consistent
#'        across data frames.
#'        If `mcmc_data` is an object of class `"emc"`, `subject_id` must be
#'        equal to `"subjects"`.
#' @param parameter String denoting the column with the estimated model
#'        parameter.
#' @param covariate String denoting the column with the observed covariate of interest.
#' @param confounders Optional character vector denoting the column(s) in `covariate_data`
#'        that are confounding or "controlling" variables, that is, variables that are
#'        associated with both `parameter` and `covariate`, and therefore contaminate
#'        the direct association between `parameter` and `covariate`.
#' @param alternative Character string specifying the alternative hypothesis,
#'        which influences how the unnormalised posterior density function is
#'        computed for each Pearson correlation coefficient:
#'   \describe{
#'     \item{"two.sided"}{ (default) Prior is supported on \eqn{\left[-1, 1\right]}}
#'     \item{"greater"}{ Prior truncated to \eqn{\left[0, 1\right]} and renormalised}
#'     \item{"less"}{ Prior truncated to \eqn{\left[-1, 0\right]} and renormalised}
#'   }
#' @param posterior_args Named list of arguments passed to
#'        [posterior_rho_updf()], which influence how the unnormalised posterior
#'        density function is computed for each Pearson correlation coefficient.
#'        The following arguments can be included:
#'  * `kappa` Numeric value controlling the "concentration" of the
#'        stretched beta prior on the correlation coefficient. Default is `1`,
#'        resulting in a uniform prior.
#'  * `n_bins` Integer denoting the number of bins ("resolution") used
#'        for the discretised grid of correlation values, for which the
#'        unnormalised posterior density values are computed. Default is `1e3`.
#'  * `max_iter` Integer denoting the maximum number of iterations
#'        (attempts) to obtain a posterior density function for a given MCMC
#'        sample. Default is `1e7`.
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
#'    If confounding variables are specified, the Pearson *partial* correlation
#'    is computed, that is, the Pearson correlation between the model parameter
#'    and covariate while accounting for the effects of confounding variables.
#'    This step handles uncertainty in the estimation of individual subjects'
#'    model parameters, and is suitable for inferences about the correlation
#'    coefficient in the given sample of subjects. However, it does not account
#'    for uncertainty in generalising from the sample to the population.
#'
#' 2. Second, it calculates the posterior density function for each (partial)
#'    correlation using the analytical solution derived by Ly et al. (2018),
#'    resulting in a set of plausible posterior distributions of the correlation.
#'    This handles uncertainty in generalizing from the sample to the
#'    population, and is therefore suitable for inferences about the
#'    correlation coefficient in the population.
#'
#' The function requires at least 3 valid observations per MCMC sample and
#' filters out any correlation values that are not finite or outside the valid
#' range of \eqn{\left[-1, 1\right]}.
#'
#' @references
#' Ly, A., Boehm, U., Heathcote, A., Turner, B. M., Forstmann, B., Marsman, M.,
#' & Matzke, D. (2017). A flexible and efficient hierarchical bayesian approach
#' to the exploration of individual differences in cognitive‐model‐based
#' neuroscience. *Computational models of brain and behavior*, 467-479. DOI:
#' [10.1002/9781119159193.ch34](https://doi.org/10.1002/9781119159193.ch34)
#'
#' Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Analytic posteriors for
#' Pearson's correlation coefficient. *Statistica Neerlandica*, *72*, 4–13. DOI:
#' [10.1111/stan.12111](https://doi.org/10.1111/stan.12111)
#'
#' Kucharský, S., Wagenmakers, E.-J., Van den Bergh, D., & Ly, A. (2023).
#' Analytic posterior distributions and Bayes Factor for Pearson Partial
#' Correlation. *PsyArXiv Preprints*. DOI:
#' [10.31234/osf.io/6muwy](https://doi.org/10.31234/osf.io/6muwy)
#'
#' @export
run_plausible_cor <- function(
    mcmc_data,
    covariate_data = NULL,
    draw_id,
    subject_id,
    parameter,
    covariate,
    confounders = NULL,
    alternative = c("two.sided", "greater", "less"),
    posterior_args = NULL
) {

  posterior_args <- posterior_args %||% list()
  posterior_args[["alternative"]] <- rlang::arg_match(alternative)
  posterior_args <- assert_posterior_args(posterior_args)

  column_names <- list(
    draw_id = draw_id,
    subject_id = subject_id,
    parameter = parameter,
    covariate = covariate,
    confounders = confounders
  )

  if (inherits(mcmc_data, "emc")) {
    mcmc_data <- prep_emc_data(
      mcmc_data = mcmc_data,
      column_names = column_names
    )
    if (column_names[["subject_id"]] != "subjects") {
      rlang::warn(
        message = paste0(
          "Overriding 'subject_id' input. ",
          "For 'emc' objects this is always 'subjects'."
        )
      )
      column_names[["subject_id"]] <- "subjects"
    }
  }

  data <- prep_plausible_cor_data(
    mcmc_data = mcmc_data,
    covariate_data = covariate_data,
    column_names = column_names
  )

  result <- run_plausible_cor_engine(
    data = data,
    column_names = column_names,
    posterior_args = posterior_args
  )

  class(result) <- c("plausible_cor", class(result))

  return(result)

}

#' @noRd
run_plausible_cor_engine <- function(
    data,
    column_names,
    posterior_args
) {

  result <- data %>%
    dplyr::group_by(.data[[".draw"]]) %>%
    dplyr::summarise(
      cor_result = list(compute_cor(
        data = dplyr::pick(dplyr::everything()),
        column_names = column_names
      )),
      .groups = "drop"
    ) %>%
    tidyr::unnest_wider(
      dplyr::all_of("cor_result")
    ) %>%
    dplyr::filter(
      .data[["n"]] >= 3 & .data[["n"]] > (.data[["k"]] + 2) &
        is.finite(.data[["r"]]) & abs(.data[["r"]]) <= 1
    )

  checkmate::assert_data_frame(
    x = result,
    all.missing = FALSE,
    min.rows = 1
  )

  result <- result %>%
    dplyr::mutate(
      posterior_updf = purrr::pmap(
        .l = list(
          .data[["r"]],
          .data[["n"]],
          .data[["k"]]
        ),
        .f = function(r_val, n_val, k_val) {
          rlang::exec(
            .fn = posterior_rho_updf,
            r = r_val,
            n = n_val - k_val,
            !!!posterior_args
          )
        }
      )
    )

  return(result)

}


#' @noRd
compute_cor <- function(
    data,
    column_names
) {

  parameter_values <- data[[column_names[["parameter"]]]]
  covariate_values <- data[[column_names[["covariate"]]]]
  confounders <- column_names[["confounders"]]
  is_partial_cor <- !is.null(confounders)

  if (!is_partial_cor) {
    r <- stats::cor(
      x = parameter_values,
      y = covariate_values,
      method = "pearson"
    )
    k <- 0
  } else {
    cor_mat <- stats::cor(
      x = cbind(
        parameter_values,
        covariate_values,
        data[confounders]
      ),
      method = "pearson"
    )
    pcor_mat <- corpcor::cor2pcor(cor_mat)
    r <- pcor_mat[1, 2, drop = TRUE]
    k <- length(confounders)
  }

  result <- list(
    r = r,
    n = nrow(data), # NB upstream check guarantees complete cases
    k = k
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
#' @param point_interval_args A named list specifying how the central tendency
#'        and credible interval(s) should be computed, or `NULL` to use defaults.
#'        Valid elements are: `point_method` (character string, either `"mean"`
#'        or `"median"` for the point summary), `interval_method` (character
#'        string, either `"hdci"` for highest density credible interval or `"qi"`
#'        for quantile interval), and `interval_width` (numeric vector with
#'        values between 0 and 1 specifying the desired interval width(s)).
#'        Any invalid arguments will be ignored with a warning. Partial
#'        specification is supported; unspecified elements will use defaults.
#'        Defaults to
#'        `list(point_method = "mean", interval_method = "hdci", interval_width = 0.95)`.
#' @param rope_range Optional numeric vector of length 2 specifying the lower
#'        and upper bounds of the region of practical equivalence (ROPE), which
#'        is used to compute the proportion of the distribution contained within
#'        the ROPE. Defaults to `NULL` in which case the ROPE is ignored. See
#'        Gignac & Szodorai (2016) for guidance on what constitutes a correlation
#'        coefficient practically equivalent to null.
#' @param posterior_grid_spacing For summarising the population-level plausible
#'        correlation: The step size used to discretise each MCMC sample's
#'        posterior density function into a grid of posterior density values.
#'
#' @return A [tibble::tbl_df-class] with one row per summary type
#'         ("sample" and "posterior") and credible interval width, containing:
#'   \item{type}{Whether the row corresponds to the sample-based or
#'         population-based summary.}
#'   \item{mean / median}{Point estimate of the correlation coefficient.}
#'   \item{lower / upper}{Credible interval bounds for each `interval_width`.}
#'   \item{width}{If the length of `interval_width` is greater than 1: The width of the credible interval.}
#'   \item{p_dir}{Directional probability (proportion of mass > 0 or < 0, whichever is greater).}
#'   \item{p_rope}{If specified: The proportion contained within the ROPE.}
#'
#' @references
#' Gignac, G. E., & Szodorai, E. T. (2016). Effect size guidelines for
#' individual differences researchers. *Personality and individual differences*,
#' *102*, 74-78. DOI: [10.1016/j.paid.2016.06.069](https://doi.org/10.1016/j.paid.2016.06.069)
#'
#' @export
summarise_plausible_cor <- function(
    .data,
    point_interval_args = NULL,
    rope_range = NULL,
    posterior_grid_spacing = 1e-3
) {

  validate_column_inputs(
    col_names = c(".draw", "r", "posterior_updf"),
    data_frame = .data,
    data_name = ".data"
  )

  point_interval_args <- validate_point_interval_args(point_interval_args)

  sample_summary <- .data %>%
    dplyr::select(dplyr::all_of("r")) %>%
    summarise_samples(
      point_method = point_interval_args[["point_method"]],
      interval_width = point_interval_args[["interval_width"]],
      interval_method = point_interval_args[["interval_method"]],
      rope_range = rope_range
    ) %>%
    dplyr::mutate(
      type = "sample"
    ) %>%
    dplyr::relocate(dplyr::all_of("type"))

  mean_posterior_rho <- .data %>%
    get_posterior_rho_densities(grid_spacing = posterior_grid_spacing) %>%
    get_mean_posterior_rho()

  rho_grid <- mean_posterior_rho[["x"]]
  mean_density <- mean_posterior_rho[["mean_density"]]

  population_point <- get_point_estimate(
    val = rho_grid,
    dens = mean_density,
    method = point_interval_args[["point_method"]]
  )

  population_interval <- get_interval(
    val = rho_grid,
    dens = mean_density,
    width = point_interval_args[["interval_width"]],
    method = point_interval_args[["interval_method"]]
  )

  population_summary <- population_interval %>%
    dplyr::mutate(
      type = "population",
      !!point_interval_args[["point_method"]] := population_point,
      p_dir = max(
        sum(mean_density[rho_grid > 0]) * posterior_grid_spacing,
        sum(mean_density[rho_grid < 0]) * posterior_grid_spacing
      )
    ) %>%
    dplyr::relocate(dplyr::all_of("type"))

  if (test_rope_range(rope_range)) {
    population_summary <- population_summary %>%
      dplyr::mutate(
        p_rope = sum(
          mean_density[rho_grid >= rope_range[1] & rho_grid <= rope_range[2]]
        ) * posterior_grid_spacing
      )
  }

  result <- dplyr::bind_rows(
    sample_summary, population_summary
  ) %>%
    dplyr::mutate(
      type = factor(.data[["type"]])
    )

  if (length(point_interval_args[["interval_width"]]) == 1) {
    result <- result %>%
      dplyr::select(-dplyr::all_of("width"))
  }

  return(result)

}

#' Compare two distributions of plausible correlation values
#'
#' @description
#' Compares the posterior distributions of plausible correlation values, as
#' returned by two separate calls to [run_plausible_cor()]. Uses two methods
#' methods for defining the distribution of the difference: a **sample-based**
#' comparison, which subtracts the correlation values from matching MCMC
#' samples, and a **population-based** comparison, which draws quantiles from
#' each posterior distribution (accounting for population uncertainty) before
#' computing the difference.
#'
#' @param x,y Data frames returned by [run_plausible_cor()], each containing
#'   the columns `.draw`, `r`, and `posterior_updf`.
#' @param point_interval_args A named list specifying how the central tendency
#'        and credible interval(s) of the difference should be computed, or
#'        `NULL` to use defaults. Valid elements are: `point_method`
#'        (character string, either `"mean"` or `"median"` for the point summary),
#'        `interval_method` (character string, either `"hdci"` for highest
#'        density credible interval or `"qi"` for quantile interval), and
#'        `interval_width` (numeric vector with values between 0 and 1
#'        specifying the desired interval width(s)).
#'        Any invalid arguments will be ignored with a warning. Partial
#'        specification is supported; unspecified elements will use defaults.
#'        Defaults to
#'        `list(point_method = "mean", interval_method = "hdci", interval_width = 0.95)`.
#' @param rope_range Optional numeric vector of length 2 specifying the lower
#'        and upper bounds of the region of practical equivalence (ROPE), which
#'        is used to compute the proportion of the distribution contained within
#'        the ROPE. Defaults to `NULL` in which case the ROPE is ignored.
#' @param n_samples Integer indicating how many quantiles to draw per posterior
#'   distribution when `comparison_method = "population"`. Default is 1.
#' @param rng_seed Optional numeric vector of length 2 to control the random
#'   seed for quantile sampling of `x` and `y`, respectively. If `NULL`,
#'   random seeds will not be set.
#'
#' @return A [tibble::tibble()] containing summary statistics of the difference
#'  between correlation coefficients, with the same format as [summarise_plausible_cor()].
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
    point_interval_args = NULL,
    rope_range = NULL,
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

  if (!test_rng_seed(rng_seed)) {
    rng_seed <- c(NA, NA)
  }

  point_interval_args <- validate_point_interval_args(point_interval_args)

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
      point_method = point_interval_args[["point_method"]],
      interval_width = point_interval_args[["interval_width"]],
      interval_method = point_interval_args[["interval_method"]],
      rope_range = rope_range
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
      point_method = point_interval_args[["point_method"]],
      interval_width = point_interval_args[["interval_width"]],
      interval_method = point_interval_args[["interval_method"]],
      rope_range = rope_range
    ) %>%
    dplyr::mutate(type = "population") %>%
    dplyr::relocate(dplyr::all_of("type"))

  result <- dplyr::bind_rows(
    sample_delta_summary, population_delta_summary
  ) %>%
    dplyr::mutate(
      type = factor(.data[["type"]])
    )

  if (length(point_interval_args[["interval_width"]]) == 1) {
    result <- result %>%
      dplyr::select(-dplyr::all_of("width"))
  }

  return(result)

}

#' Evaluate posterior density functions over a grid
#'
#' @description
#' Helper function that evaluates each posterior density function in the output
#' of [run_plausible_cor()] over a shared grid in the range \eqn{\left[-1, 1\right]}.
#' Typical users would not have to call this function, instead relying on
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

  max_abs_r <- 0.999
  n_steps <- floor(max_abs_r / grid_spacing)
  rho_grid <- seq(
    from = -n_steps * grid_spacing,
    to = n_steps * grid_spacing,
    by = grid_spacing
  )

  single_density_grid <- function(
    draw_id, r_val, updf, grid = rho_grid, dx = grid_spacing
  ) {
    if (!checkmate::test_function(updf)) {
      return(NULL)
    }
    dens <- updf(grid)
    if (!test_densities(dens, test_len = length(grid))) {
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
#' Accepts an optional `dx` argument to avoid relying on floating-point
#' differences.
#'
#' @param .data Data frame output from [get_posterior_rho_densities()].
#' @param dx Grid spacing (step size). If `NULL` (default), it is estimated
#'        from the unique values of `x` (grid points) using [diff()].
#'
#' @return A tibble with columns `x` (grid points) and `density`
#'         (mean posterior density).
#'
#' @export
get_mean_posterior_rho <- function(.data, dx = NULL) {

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

  if (is.null(dx)) {
    dxs <- diff(sort(unique(result[["x"]])))
    dxs_rounded <- round(dxs, digits = 10)
    if (length(unique(dxs_rounded)) > 1) {
      rlang::warn(
        message = "Grid appears irregular; using first diff() as dx."
      )
    }
    dx <- dxs_rounded[1]
  }

  result <- result %>%
    dplyr::mutate(
      mean_density = .data[["mean_density"]] /
        sum(.data[["mean_density"]] * dx)
    )

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
#' @param .data A data frame output by [run_plausible_cor()], containing one
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
#' credible interval(s) (e.g., highest-density continous intervals or quantiles),
#' the directionality of the distribution (positive or negative), and optionally
#' the proportion contained in a region of practical equivalence (ROPE).
#'
#' @param .data A data frame containing a single numeric column representing
#'   posterior samples or any distribution of interest.
#' @param ... Additional arguments passed to [ggdist::point_interval()].
#' @param point_method A string indicating the point estimate to compute. Can
#'   either be `"mean"` (default) or `"median"`.
#' @param interval_width Numeric vector of values between 0 and 1 specifying the
#'   width(s) of the credible interval(s) to calculate. Default is `0.95`.
#' @param interval_method A string indicating the method for computing intervals.
#'   Can either be `"hdci"` (highest-density continous intervals; default) or
#'   `"qi"` (quantile intervals).
#' @param rope_range Optional numeric vector of length 2 specifying the lower
#'        and upper bounds of the region of practical equivalence (ROPE), which
#'        is used to compute the proportion of the distribution contained within
#'        the ROPE. Defaults to `NULL` in which case the ROPE is ignored.
#'
#' @return A [tibble::tibble] summarizing the posterior samples.
#'
#' @noRd
summarise_samples <- function(
    .data,
    ...,
    point_method = c("mean", "median"),
    interval_width = 0.95,
    interval_method = c("hdci", "qi"),
    rope_range = NULL
) {

  data <- .data %>%
    dplyr::select(dplyr::where(is.numeric))

  checkmate::assert_data_frame(
    x = data,
    any.missing = FALSE,
    ncols = 1
  )
  varname <- colnames(data)

  p_dir <- max(
    mean(data[[varname]] > 0),
    mean(data[[varname]] < 0)
  )

  p_rope <- NULL
  if (test_rope_range(rope_range)) {
    p_rope <- mean(
      data[[varname]] >= rope_range[1] & data[[varname]] <= rope_range[2]
    )
  }

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

  if (!is.null(p_rope)) {
    result <- result %>%
      dplyr::mutate(
        p_rope = p_rope
      )
  }

  return(result)

}


#' @noRd
prep_plausible_cor_data <- function(
    mcmc_data,
    covariate_data,
    column_names
) {

  validate_column_inputs(
    col_names = column_names[c("draw_id", "subject_id", "parameter", "covariate")]
  )

  validate_column_inputs(
    col_names = column_names[c("draw_id", "subject_id", "parameter")],
    data_frame = mcmc_data,
    data_name = "mcmc_data"
  )

  covariate <- column_names[["covariate"]]
  subject_id <- column_names[["subject_id"]]
  confounders <- column_names[["confounders"]]

  if (covariate %in% names(mcmc_data)) {

    data <- mcmc_data

  } else if (!is.null(covariate_data)) {

    validate_column_inputs(
      col_names = c(
        subject_id = subject_id,
        covariate = covariate,
        if (!is.null(confounders)) {
          stats::setNames(confounders, paste0("confounder_", confounders))
        }
      ),
      data_frame = covariate_data,
      data_name = "covariate_data"
    )

    covariate_data_clean <- covariate_data %>%
      dplyr::select(dplyr::all_of(c(subject_id, covariate, confounders))) %>%
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
        "Covariate column '", covariate, "' not found in 'mcmc_data', ",
        "and no 'covariate_data' provided."
      )
    )
  }

  keep_cols <- c(
    ".draw",
    column_names[["parameter"]],
    covariate,
    confounders
  )

  data <- data %>%
    dplyr::rename(dplyr::all_of(c(.draw = column_names[["draw_id"]]))) %>%
    dplyr::select(dplyr::all_of(keep_cols))

  if (anyNA(data)) {
    rlang::warn(message = paste0(
      "Removing any rows with missing values in ",
      paste(keep_cols, collapse = ", ")
    ))
    data <- tidyr::drop_na(data)
  }

  data <- data %>%
    dplyr::arrange(.data[[".draw"]])

  return(data)

}

#' @noRd
prep_emc_data <- function(
    mcmc_data,
    column_names
) {

  rlang::check_installed(
    pkg = "EMC2",
    reason = "'EMC2' must be installed if input `mcmc_data` is of class 'emc'."
  )

  mcmc_data <- EMC2::parameters(mcmc_data, selection = "alpha")

  validate_column_inputs(
    col_names = c("subjects", column_names[["parameter"]]),
    data_frame = mcmc_data,
    data_name = "mcmc_data"
  )

  n_draws <- nrow(mcmc_data) / dplyr::n_distinct(mcmc_data[["subjects"]])

  mcmc_data <- mcmc_data %>%
    dplyr::group_by(
      .data[["subjects"]]
    ) %>%
    dplyr::mutate(
      .draw = seq_len(n_draws)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(
      .data[[".draw"]],
      .before = 1
    ) %>%
    dplyr::rename(
      dplyr::all_of(c(.draw = column_names[["draw_id"]]))
    )

  return(mcmc_data)

}
