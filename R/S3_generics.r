#' Print method for plausible correlation list objects
#'
#' @param x An object of class "plausible_cor_list" (a list of plausible_cor objects).
#' @param ... Optional additional arguments (ignored).
#'
#' @export
#' @method print plausible_cor_list
print.plausible_cor_list <- function(x, ...) {
  cat("<plausible_cor_list object>\n")
  cat(" Number of plausible correlation analyses: ", length(x), "\n", sep = "")
  cat(" Shared settings:", "\n", sep = "")
  purrr::walk(
    .x = c(
      "method", "alternative", "confounders"
    ),
    .f = function(name) {
      print_plausible_cor_attr(x = x[[1]], name = name, n_space = 2L)
    }
  )
  cat("\n", sep = "")

  print_out <- purrr::map(
    .x = seq_along(x),
    .f = function(idx) {
      obj <- x[[idx]]
      return(
        tibble::tibble_row(
          parameter = attr(obj, "parameter", exact = TRUE),
          covariate = attr(obj, "covariate", exact = TRUE),
          n_draws = dplyr::n_distinct(obj[[".draw"]]),
          n_valid_post = n_valid_updf(obj)
        )
      )
    }
  ) %>%
    purrr::list_rbind()

  cat(" Analysis-specific details:", "\n", sep = "")
  print(print_out, ...)
  return(invisible(x))
}

#' Summary method for plausible correlation list objects
#'
#' @description
#' Summarises a list of plausible correlation objects returned by
#' [run_plausible_cor()]. Each element of the list is summarised
#' via [summary.plausible_cor()], and the results are combined into a single
#' tibble for convenient inspection or downstream processing.
#'
#' @param object An object of class "plausible_cor_list": the output from
#'    [run_plausible_cor()] when multiple parameter–covariate pairs are analysed.
#' @param point_interval_args Optional named list specifying how the central
#'    tendency and credible interval(s) should be computed. See
#'    [summarise_plausible_cor()] for details.
#' @param rope_range Optional numeric vector of length 2 specifying the ROPE
#'    bounds. See [summarise_plausible_cor()] for details.
#' @param ... Additional arguments passed to [summary.plausible_cor()].
#'
#' @return A tibble with class `"summary.plausible_cor_list"`.
#' @export
#' @method summary plausible_cor_list
summary.plausible_cor_list <- function(
    object,
    point_interval_args = NULL,
    rope_range = NULL,
    ...
) {
  result <- purrr::map(
    .x = object,
    summary,
    point_interval_args = point_interval_args,
    rope_range = rope_range,
    ...
  ) %>%
    purrr::list_rbind(
      names_to = "analysis"
    ) %>%
    tidyr::separate_wider_delim(
      cols = tidyr::all_of("analysis"),
      delim = " ~ ",
      names = c("parameter", "covariate")
    )

  class(result) <- c("summary.plausible_cor_list", class(result))
  attr(result, "n_analyses") <- length(object)
  attr(result, "method") <- attr(object[[1]], "method", exact = TRUE)
  attr(result, "alternative") <- attr(object[[1]], "alternative", exact = TRUE)
  attr(result, "confounders") <- attr(object[[1]], "confounders", exact = TRUE)
  attr(result, "rope_range") <- rope_range
  return(result)
}

#' Print method for summary.plausible_cor_list objects
#'
#' @param x An object of class "summary.plausible_cor_list".
#' @param ... Additional arguments passed to tibble's print method.
#'
#' @export
#' @method print summary.plausible_cor_list
print.summary.plausible_cor_list <- function(x, ...) {
  cat("<summary.plausible_cor_list object>\n")
  cat(
    " Number of plausible correlation analyses: ",
    attr(x, "n_analyses", exact = TRUE), "\n", sep = ""
  )
  cat(" Shared settings:", "\n", sep = "")
  purrr::walk(
    .x = c(
      "method", "alternative", "confounders", "rope_range"
    ),
    .f = function(name) {
      print_plausible_cor_attr(x = x, name = name, n_space = 2L)
    }
  )
  cat("\nSummary of plausible correlation estimates:\n")
  NextMethod("print", x, ...)
  return(invisible(x))
}


#' Print method for plausible correlation objects
#'
#' @param x An object of class "plausible_cor" (output from [run_plausible_cor()]).
#' @param ... Optional additional arguments, passed to tibble's print method.
#'
#' @importFrom utils head
#' @export
#' @method print plausible_cor
print.plausible_cor <- function(x, ...) {
  cat("<plausible_cor object>\n")
  purrr::walk(
    .x = c(
      "parameter", "covariate", "confounders", "method", "alternative"
    ),
    .f = function(name) {print_plausible_cor_attr(x = x, name = name)}
  )
  cat(" Number of MCMC samples: ", dplyr::n_distinct(x[[".draw"]]), "\n", sep = "")
  cat(" Number of valid posterior density functions: ", n_valid_updf(x), "\n", sep = "")

  r_preview <- utils::head(x[["r"]], n = 5)
  cat(
    " Example r values: ",
    paste0(round(r_preview, digits = 3), collapse = ", ")
  )
  if (nrow(x) > 5) {
    cat(", ...")
  }
  cat("\n\n")

  NextMethod("print", x, ...)
  return(invisible(x))
}


#' Summary method for plausible correlation objects
#'
#' @param object An object of class "plausible_cor" (output from [run_plausible_cor()]).
#' @param point_interval_args A named list specifying how the central tendency
#'        and credible interval(s) should be computed. See [summarise_plausible_cor()]
#'        for details.
#' @param rope_range Optional numeric vector of length 2 specifying the ROPE bounds.
#'        See [summarise_plausible_cor()] for details.
#' @param ... Additional arguments (currently unused)
#'
#' @return A summary object with class "summary.plausible_cor"
#' @export
#' @method summary plausible_cor
summary.plausible_cor <- function(
    object,
    point_interval_args = NULL,
    rope_range = NULL,
    ...
) {
  result <- summarise_plausible_cor(
    .data = object,
    point_interval_args = point_interval_args,
    rope_range = rope_range
  )
  class(result) <- c("summary.plausible_cor", class(result))
  attr(result, "parameter") <- attr(object, "parameter", exact = TRUE)
  attr(result, "covariate") <- attr(object, "covariate", exact = TRUE)
  attr(result, "confounders") <- attr(object, "confounders", exact = TRUE)
  attr(result, "method") <- attr(object, "method", exact = TRUE)
  attr(result, "alternative") <- attr(object, "alternative", exact = TRUE)
  attr(result, "rope_range") <- rope_range
  return(result)
}

#' Print method for summary.plausible_cor objects
#'
#' @param x An object of class "summary.plausible_cor".
#' @param ... Optional additional arguments, passed to tibble's print method.
#'
#' @export
#' @method print summary.plausible_cor
print.summary.plausible_cor <- function(x, ...) {
  cat("<summary.plausible_cor object>\n")
  purrr::walk(
    .x = c(
      "parameter", "covariate", "confounders", "method", "alternative", "rope_range"
    ),
    .f = function(name) {print_plausible_cor_attr(x = x, name = name)}
  )
  cat("\nSummary of plausible correlation estimates:\n")
  NextMethod("print", x, ...)
  return(invisible(x))
}


#' Plot method for plausible correlation objects
#'
#' @param x An object of class `"plausible_cor"` (output from [run_plausible_cor()])
#' @param type Character string specifying the plot type. Options:
#'        - `"population"` (default): Shows posterior distribution of population-level correlation
#'        - `"sample"`: Shows distribution of sample-level correlations
#' @param ... Additional arguments passed to the underlying plotting function:
#'        [plot_population_cor()] for `type == "population"` or
#'        [plot_sample_cor()] for `type == "sample"`. See those functions for details.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#' @method plot plausible_cor
plot.plausible_cor <- function(
    x,
    type = c("population", "sample"),
    ...
) {
  type <- match.arg(type)
  result <- switch(
    type,
    "population" = {
      plot_population_cor(
        .data = x,
        ...
      )
    },
    "sample" = {
      plot_sample_cor(
        .data = x,
        ...
      )
    }
  )
  return(result)
}


#' @noRd
print_plausible_cor_attr <- function(x, name, n_space = 1L) {
  if (name == "confounders") {
    return(print_confounders(x))
  }
  attribute <- attr(x, name, exact = TRUE)
  if (is.null(attribute)) {
    return(invisible(x = NULL))
  }
  name_label <- dplyr::case_when(
    name == "parameter" ~ "Parameter",
    name == "covariate" ~ "Covariate",
    name == "method" ~ "Correlation method",
    name == "alternative" ~ "Alternative hypothesis",
    name == "rope_range" ~ "ROPE bounds",
    TRUE ~ name
  )
  if (name == "rope_range") {
    attribute <- paste0(attribute, collapse = ", ")
  }
  cat(
    paste0(rep(" ", times = n_space), collapse = ""),
    name_label, ": ", attribute, "\n",
    sep = ""
  )
  return(invisible(x = NULL))
}

#' @noRd
print_confounders <- function(x) {
  confounders <- attr(x, "confounders", exact = TRUE)
  if (!is.null(confounders)) {
    if (length(confounders) > 5L) {
      confounders_print <- paste0(
        paste(confounders[1:4], collapse = ", "),
        ", ..., ",
        confounders[length(confounders)]
      )
    } else {
      confounders_print <- paste(confounders, collapse = ", ")
    }
    cat(" Confounders: ", confounders_print, "\n", sep = "")
  }
  return(invisible(x = NULL))
}

#' @noRd
n_valid_updf <- function(x) {
  is_valid <- x[["posterior_updf"]] %>%
    purrr::map_lgl(
      .f = function(updf) {
        out <- tryCatch(
          updf(0),
          error = function(e) {NA_real_}
        )
        result <- is.finite(out) && !is.na(out)
        return(result)
      }
    )
  result <- sum(is_valid)
  return(result)
}
