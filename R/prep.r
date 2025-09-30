#' @noRd
prep_plausible_cor_data <- function(
    mcmc_data,
    covariate_data,
    column_names,
    n_draws,
    rng_seed
) {

  validate_column_inputs(
    col_names = column_names[c("draw_id", "subject_id", "parameter", "covariate")]
  )
  validate_column_inputs(
    col_names = column_names[["draw_id"]],
    data_frame = mcmc_data,
    data_name = "mcmc_data"
  )
  mcmc_data <- mcmc_data %>%
    dplyr::rename(
      dplyr::all_of(c(.draw = column_names[["draw_id"]]))
    )

  parameter_data <- get_parameter_data(
    mcmc_data = mcmc_data,
    subject_id = column_names[["subject_id"]],
    parameter = column_names[["parameter"]]
  )

  covariate_data <- get_covariate_data(
    mcmc_data = mcmc_data,
    covariate_data = covariate_data,
    subject_id = column_names[["subject_id"]],
    covariate = column_names[["covariate"]]
  )

  join_by <- column_names[["subject_id"]]
  if (".draw" %in% names(covariate_data)) {
    join_by <- c(".draw", join_by)
  }
  result <- dplyr::left_join(
    x = parameter_data,
    y = covariate_data,
    by = join_by
  )

  if (!is.null(column_names[["confounders"]])) {
    validate_column_inputs(col_names = column_names[["confounders"]])
    confounder_data <- get_confounder_data(
      mcmc_data = mcmc_data,
      covariate_data = covariate_data,
      subject_id = column_names[["subject_id"]],
      confounders = column_names[["confounders"]]
    )
    join_by <- column_names[["subject_id"]]
    if (".draw" %in% names(confounder_data)) {
      join_by <- c(".draw", join_by)
    }
    result <- dplyr::left_join(
      x = result,
      y = confounder_data,
      by = join_by
    )
  }

  if (anyNA(result)) {
    n_removed <- sum(!stats::complete.cases(result))
    rlang::warn(
      message = paste0(
        "Removing ", n_removed, " row(s) with missing values after joining."
      )
    )
    result <- tidyr::drop_na(result)
  }

  result <- result %>%
    dplyr::arrange(
      .data[[".draw"]]
    ) %>%
    draw_rows(
      n_draws = n_draws,
      rng_seed = rng_seed
    )

  return(result)

}

#' @noRd
prep_emc_data <- function(x) {
  rlang::check_installed(
    pkg = "EMC2",
    reason = "'EMC2' must be installed if input `mcmc_data` is of class 'emc'."
  )
  result <- EMC2::parameters(x, selection = "alpha")
  n_draws <- nrow(result) / dplyr::n_distinct(result[["subjects"]])
  result <- result %>%
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
    )
  return(result)
}

#' @noRd
get_parameter_data <- function(
    mcmc_data,
    subject_id,
    parameter
) {
  validate_column_inputs(
    col_names = c(".draw", subject_id, parameter),
    data_frame = mcmc_data,
    data_name = "mcmc_data"
  )
  result <- mcmc_data %>%
    dplyr::select(
      dplyr::all_of(c(".draw", subject_id, parameter))
    )
  return(result)
}

#' @noRd
get_covariate_data <- function(
    mcmc_data,
    covariate_data,
    subject_id,
    covariate
) {
  # starting assumption: covariate is in mcmc_data
  if (covariate %in% names(mcmc_data)) {
    validate_column_inputs(
      col_names = c(".draw", subject_id, covariate),
      data_frame = mcmc_data,
      data_name = "mcmc_data"
    )
    result <- mcmc_data %>%
      dplyr::select(
        dplyr::all_of(c(".draw", subject_id, covariate))
      )
    return(result)
  }
  # if we're not done yet, then covariate must be in covariate_data
  if (!is.null(covariate_data)) {
    if (covariate %in% names(covariate_data)) {
      validate_column_inputs(
        col_names = c(subject_id, covariate),
        data_frame = covariate_data,
        data_name = "covariate_data"
      )
      result <- covariate_data %>%
        dplyr::select(
          dplyr::all_of(c(subject_id, covariate))
        )
      n_subjects <- dplyr::n_distinct(
        result[[subject_id]]
      )
      if (nrow(result) != n_subjects) {
        rlang::abort(
          message = "'covariate_data' should have one unique value per subject."
        )
      }
      return(result)
    } else {
      rlang::abort(
        message = paste0(
          "Covariate '", covariate, "' not found in 'mcmc_data' ",
          "or 'covariate_data'."
        )
      )
    }
  } else {
    rlang::abort(
      message = paste0(
        "Covariate '", covariate, "' not found in 'mcmc_data' ",
        "and no 'covariate_data' provided."
      )
    )
  }
  return(invisible(x = NULL))
}

#' @noRd
get_confounder_data <- function(
  mcmc_data,
  covariate_data,
  subject_id,
  confounders
) {
  # initialize holders to avoid "object not found" ambiguity
  confounders_mcmc_data <- NULL
  confounders_cov_data  <- NULL
  confs_in_mcmc <- character(0)
  confs_in_cov  <- character(0)

  # first check for any confounders in mcmc_data
  confs_in_mcmc <- intersect(confounders, names(mcmc_data))
  if (length(confs_in_mcmc) > 0) {
    validate_column_inputs(
      col_names = c(".draw", subject_id, confs_in_mcmc),
      data_frame = mcmc_data,
      data_name = "mcmc_data"
    )
    confounders_mcmc_data <- mcmc_data %>%
      dplyr::select(
        dplyr::all_of(c(".draw", subject_id, confs_in_mcmc))
      )
    confounders_data_complete <- checkmate::test_names(
      x = confs_in_mcmc,
      permutation.of = confounders
    )
    if (confounders_data_complete) {
      return(confounders_mcmc_data)
    }
  }

  # if we're not done yet, there must be at least some confounders in
  # covariate_data
  if (!is.null(covariate_data)) {
    confs_in_cov <- intersect(confounders, names(covariate_data))
    confs_in_cov <- setdiff(confs_in_cov, confs_in_mcmc)
    if (length(confs_in_cov) > 0) {
      validate_column_inputs(
        col_names = c(subject_id, confs_in_cov),
        data_frame = covariate_data,
        data_name = "covariate_data"
      )
      confounders_cov_data <- covariate_data %>%
        dplyr::select(
          dplyr::all_of(c(subject_id, confs_in_cov))
        )
      confounders_data_complete <- checkmate::test_names(
        x = confs_in_cov,
        permutation.of = confounders
      )
      if (confounders_data_complete) {
        return(confounders_cov_data)
      }
    }
  } else if (length(confs_in_mcmc) == 0) {
    rlang::abort(
      message = paste0(
        "Confounder(s) ", paste(confounders, collapse = ", "),
        " not found in 'mcmc_data', and no 'covariate_data' provided."
      )
    )
  }

  # merge confounders from mcmc_data and covariate_data (if applicable)
  if (length(confs_in_mcmc) > 0 && length(confs_in_cov) > 0) {
    result <- dplyr::left_join(
      x = confounders_mcmc_data,
      y = confounders_cov_data,
      by = subject_id
    )
  } else if (length(confs_in_mcmc) > 0) {
    result <- confounders_mcmc_data
  } else if (length(confs_in_cov) > 0) {
    result <- confounders_cov_data
  } else {
    rlang::abort(
      message = paste0(
        "None of the confounder(s) were found in either 'mcmc_data' ",
        "or 'covariate_data'."
      )
    )
  }

  # final completeness check
  confs_missing <- setdiff(confounders, c(confs_in_mcmc, confs_in_cov))
  if (length(confs_missing) > 0) {
    rlang::abort(
      message = paste0(
        "Confounder(s) ", paste(confs_missing, collapse = ", "),
        " not found in either 'mcmc_data' or 'covariate_data'."
      )
    )
  }

  return(result)

}


#' @noRd
draw_rows <- function(.data, n_draws, rng_seed) {
  if (is.null(n_draws) || !is.finite(n_draws)) {
    return(.data)
  }
  n_draws <- checkmate::assert_count(
    x = n_draws,
    coerce = TRUE
  )
  if (n_draws >= nrow(.data)) {
    return(.data)
  }
  if (!is.null(rng_seed)) {
    rng_seed <- checkmate::assert_integerish(
      x = rng_seed,
      any.missing = FALSE,
      len = 1,
      coerce = TRUE
    )
    result <- withr::with_seed(
      seed = rng_seed,
      code = .data %>%
        dplyr::slice_sample(
          n = n_draws,
          replace = FALSE
        )
    )
  } else {
    result <- .data %>%
      dplyr::slice_sample(
        n = n_draws,
        replace = FALSE
      )
  }
  return(result)
}
