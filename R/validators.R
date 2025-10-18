# PLAUSIBLE_COR ---------------------------------------------------------------

#' Check posterior_args parameter
#'
#' @description
#' Internal helper function to validate the structure and contents of the
#' posterior_args parameter passed to functions that use [posterior_cor_updf()].
#' Allowed arguments are: kappa, alternative, n_bins, and max_iter.
#'
#' @param posterior_args A named list of arguments or NULL
#'
#' @return The checked posterior_args list
#'
#' @noRd
assert_posterior_args <- function(posterior_args) {

  if (is.null(posterior_args)) {
    return(NULL)
  }

  checkmate::assert_list(
    x = posterior_args,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_names(
    x = names(posterior_args),
    type = "strict",
    subset.of = c(
      "alternative", "method", "kappa", "n_bins", "max_iter"
    )
  )

  return(posterior_args)

}

#' @noRd
validate_column_names <- function(
    parameter,
    covariate,
    draw_id = NULL,
    subject_id = NULL,
    confounders = NULL
) {

  checkmate::assert_character(
    x = parameter, any.missing = FALSE, min.len = 1
  )
  checkmate::assert_character(
    x = covariate, any.missing = FALSE, min.len = 1
  )
  if (!is.null(confounders)) {
    checkmate::assert_character(
      x = confounders, any.missing = FALSE, min.len = 1, unique = TRUE
    )
  }

  if (length(parameter) == 1L && length(covariate) > 1L) {
    parameter <- rep(parameter, length(covariate))
  } else if (length(covariate) == 1L && length(parameter) > 1L) {
    covariate <- rep(covariate, length(parameter))
  }
  if (length(parameter) != length(covariate)) {
    rlang::abort(
      message = paste0(
        "`parameter` and `covariate` must have equal lengths ",
        "or one must have length 1."
      )
    )
  }

  overlap_idx <- which(parameter == covariate)
  if (length(overlap_idx) > 0L) {
    rlang::abort(
      paste0(
        "Parameter and covariate names overlap in the following set(s): ",
        paste(overlap_idx, collapse = ", ")
      )
    )
  }

  if (!is.null(confounders)) {
    if (any(parameter %in% confounders) || any(covariate %in% confounders)) {
      rlang::abort(
        "A variable cannot be both a parameter/covariate and a confounder."
      )
    }
  }

  result <- purrr::map2(
    .x = parameter,
    .y = covariate,
    .f = function(x, y) {
      return(
        list(
          draw_id = draw_id %||% ".draw",
          subject_id = subject_id %||% "subjects",
          parameter = x,
          covariate = y,
          confounders = confounders
        )
      )
    }
  )
  names(result) <- paste0(parameter, " ~ ", covariate)

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
    col_names_vec <- as.character(col_names)
    checkmate::assert_data_frame(
      x = data_frame,
      min.rows = 1,
      min.cols = length(col_names_vec),
      .var.name = data_name
    )
    checkmate::assert_names(
      x = names(data_frame),
      type = "strict",
      must.include = col_names_vec,
      .var.name = data_name
    )
  }

  return(invisible(x = NULL))
}

#' @noRd
validate_point_interval_args <- function(x) {

  defaults <- list(
    point_method = "mean",
    interval_method = "hdci",
    interval_width = 0.95
  )

  if (is.null(x)) {
    return(defaults)
  }

  checkmate::assert_list(
    x = x,
    any.missing = FALSE,
    min.len = 1,
    max.len = 3
  )
  checkmate::assert_names(
    x = names(x),
    type = "strict",
    subset.of = c("point_method", "interval_method", "interval_width")
  )

  result <- defaults
  result[names(x)] <- x

  checkmate::assert_string(
    x = result[["point_method"]]
  )
  checkmate::assert_names(
    x = result[["point_method"]],
    subset.of = c("mean", "median")
  )

  checkmate::assert_string(
    x = result[["interval_method"]]
  )
  checkmate::assert_names(
    x = result[["interval_method"]],
    subset.of = c("hdci", "qi")
  )

  checkmate::assert_double(
    x = result[["interval_width"]],
    lower = 0 + sqrt(.Machine$double.eps),
    upper = 1 - sqrt(.Machine$double.eps),
    any.missing = FALSE
  )

  return(result)
}


#' @noRd
test_rope_range <- function(x) {
  return(
    checkmate::test_numeric(
      x = x,
      lower = -1,
      upper = 1,
      finite = TRUE,
      any.missing = FALSE,
      sorted = TRUE,
      len = 2,
      null.ok = FALSE
    )
  )
}

#' @noRd
test_rng_seed <- function(x) {
  return(
    checkmate::test_integerish(
      x = x,
      any.missing = FALSE,
      len = 2
    )
  )
}

#' @noRd
test_densities <- function(x, test_len = NULL) {
  return(
    checkmate::test_double(
      x = x,
      finite = TRUE,
      any.missing = FALSE,
      len = test_len
    )
  )
}

# POSTERIOR -------------------------------------------------------------------

#' Validate inputs for posterior_cor_updf function
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
validate_posterior_cor_updf_input <- function(r, n, kappa, n_bins, max_iter) {

  assert_double_wrap <- function(x, l = -Inf, u = Inf) {
    return(
      checkmate::assert_double(
        x = x, lower = l, upper = u, any.missing = FALSE, len = 1
      )
    )
  }
  assert_integerish_wrap <- function(x, l = -Inf, u = Inf) {
    return(
      checkmate::assert_integerish(
        x = x, lower = l, upper = u, any.missing = FALSE, len = 1, coerce = TRUE
      )
    )
  }

  result <- list(
    r = assert_double_wrap(x = r, l = -1, u = 1),
    n = assert_integerish_wrap(x = n, l = 3),
    kappa = assert_double_wrap(x = kappa, l = 0 + sqrt(.Machine$double.eps)),
    n_bins = assert_integerish_wrap(x = n_bins, l = 100),
    max_iter = assert_integerish_wrap(x = max_iter, l = 1)
  )

  return(result)
}

# PLOTTING --------------------------------------------------------------------

#' @noRd
validate_geom_args <- function(
    aes_list,
    allowed_names = character(),
    disallowed_names = character(),
    arg_name = deparse(substitute(aes_list))
) {
  if (!is.list(aes_list)) {
    rlang::abort(
      message = sprintf("`%s` must be a list or FALSE.", arg_name)
    )
  }
  if (isFALSE(aes_list) || is.null(aes_list)) {
    return(invisible(x = NULL))
  }

  actual_names <- names(aes_list)
  disallowed_found <- intersect(actual_names, disallowed_names)
  if (length(disallowed_found) > 0) {
    rlang::abort(
      message = sprintf(
        "In `%s`, the following names are not allowed: %s",
        arg_name, paste(disallowed_found, collapse = ", ")
      )
    )
  }

  if (length(allowed_names) >= 1) {
    unknown_names <- setdiff(actual_names, allowed_names)
    if (length(unknown_names) > 0) {
      rlang::abort(
        message = sprintf(
          "In `%s`, unknown aesthetic name(s): %s\nAllowed: %s",
          arg_name,
          paste(unknown_names, collapse = ", "),
          paste(allowed_names, collapse = ", ")
        )
      )
    }
  }

  return(invisible(x = NULL))
}
