#' Print method for plausible correlation objects
#'
#' @param x An object of class "plausible_cor" (output from [run_plausible_cor()]).
#' @param ... Optional additional arguments, passed to tibble's print method.
#'
#' @importFrom utils head
#' @export
#' @method print plausible_cor
print.plausible_cor <- function(x, ...) {
  # Check posterior_updf validity (cheap probe)
  if ("posterior_updf" %in% names(x)) {
    n_valid <- sum(
      purrr::map_lgl(
        .x = x[["posterior_updf"]],
        .f = function(updf) {
          out <- tryCatch(
            updf(0),
            error = function(e) {NA_real_}
          )
          result <- is.finite(out) && !is.na(out)
          return(result)
        }
      )
    )
  } else {
    n_valid <- nrow(x)
  }

  cat("<plausible_cor object>\n")
  cat(" Parameter: ", attr(x, "parameter", exact = TRUE), "\n", sep = "")
  cat(" Covariate: ", attr(x, "covariate", exact = TRUE), "\n", sep = "")
  if (!is.null(attr(x, "confounders", exact = TRUE))) {
    cat(
      " Confounders: ",
      paste(attr(x, "confounders", exact = TRUE), collapse = ", "),
      "\n", sep = ""
    )
  }
  cat(" Method: ", attr(x, "method", exact = TRUE), "\n", sep = "")
  cat(" Alternative hypothesis: ", attr(x, "alternative", exact = TRUE), "\n", sep = "")
  cat(" Number of MCMC samples: ", dplyr::n_distinct(x[[".draw"]]), "\n", sep = "")
  cat(" Number of valid posterior density functions: ", n_valid, "\n", sep = "")

  if ("r" %in% names(x)) {
    r_preview <- utils::head(x[["r"]], n = 5)
    cat(
      " Example r values: ",
      paste0(round(r_preview, digits = 3), collapse = ", ")
    )
    if (nrow(x) > 5) {
      cat(", ...")
    }
    cat("\n\n")
  }

  # delegate further input arguments to tibble's print method
  NextMethod("print", x, ...)
  invisible(x)
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
  cat(" Parameter: ", attr(x, "parameter", exact = TRUE), "\n", sep = "")
  cat(" Covariate: ", attr(x, "covariate", exact = TRUE), "\n", sep = "")
  if (!is.null(attr(x, "confounders", exact = TRUE))) {
    cat(
      " Confounders: ",
      paste(attr(x, "confounders", exact = TRUE), collapse = ", "),
      "\n", sep = ""
    )
  }
  cat(" Method: ", attr(x, "method", exact = TRUE), "\n", sep = "")
  cat(" Alternative hypothesis: ", attr(x, "alternative", exact = TRUE), "\n", sep = "")
  if (!is.null(attr(x, "rope_range", exact = TRUE))) {
    cat(
      " ROPE: [",
      paste0(attr(x, "rope_range", exact = TRUE), collapse = ", "),
      "]\n", sep = ""
    )
  }
  cat("\nSummary of plausible correlation estimates:\n")
  # delegate to tibble's print method
  NextMethod("print", x, ...)
  invisible(x)
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
