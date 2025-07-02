# PLAUSIBLE_COR ---------------------------------------------------------------

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

  valid_point_methods <- c("mean", "median")
  valid_interval_methods <- c("hdci", "qi")
  valid_arg_names <- c("point_method", "interval_method", "interval_width")

  if (!is.list(x)) {
    rlang::abort(message = "point_interval_args must be a list or NULL")
  }

  provided_names <- names(x)
  if (!is.null(provided_names)) {
    invalid_names <- setdiff(provided_names, valid_arg_names)
    if (length(invalid_names) > 0) {
      rlang::warn(
        message = paste(
          "Ignoring invalid arguments in point_interval_args:",
          paste(invalid_names, collapse = ", ")
        )
      )
      x <- x[valid_arg_names]
    }
  }

  result <- defaults
  result[names(x)] <- x

  if (!is.null(result[["point_method"]])) {
    if (
      !is.character(result[["point_method"]]) ||
      length(result[["point_method"]]) != 1
    ) {
      rlang::abort(message = "point_method must be a single character string")
    }
    if (!result[["point_method"]] %in% valid_point_methods) {
      rlang::abort(
        message = paste(
          "point_method must be one of:",
          paste(valid_point_methods, collapse = ", ")
        )
      )
    }
  }

  if (!is.null(result[["interval_method"]])) {
    if (
      !is.character(result[["interval_method"]]) ||
      length(result[["interval_method"]]) != 1
    ) {
      rlang::abort(message = "interval_method must be a single character string")
    }
    if (!result[["interval_method"]] %in% valid_interval_methods) {
      rlang::abort(
        message = paste(
          "interval_method must be one of:",
          paste(valid_interval_methods, collapse = ", ")
        )
      )
    }
  }

  if (!is.null(result[["interval_width"]])) {
    if (!is.numeric(result[["interval_width"]])) {
      rlang::abort(message = "interval_width must be numeric")
    }
    if (
      any(result[["interval_width"]] <= 0 | result[["interval_width"]] >= 1)
    ) {
      rlang::abort(
        message = "interval_width values must be between 0 and 1 (exclusive)"
      )
    }
  }

  return(result)
}


# POSTERIOR_RHO ---------------------------------------------------------------

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
