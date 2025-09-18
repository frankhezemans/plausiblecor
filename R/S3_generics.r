#' Print method for plausible correlation objects
#'
#' @param x An object of class "plausible_cor" (output from [run_plausible_cor()]).
#' @param ... Additional arguments (currently unused).
#'
#' @export
#' @method print plausible_cor
print.plausible_cor <- function(x, ...) {

  n_draws <- dplyr::n_distinct(x[[".draw"]])
  n_valid <- nrow(x)

  cat("<plausible_cor object>\n")
  cat(" Number of posterior draws: ", n_draws, "\n", sep = "")
  cat(" Valid correlations: ", n_valid, "\n", sep = "")

  if ("r" %in% names(x)) {
    r_preview <- utils::head(x[["r"]], n = 5)
    cat(" Example r values: ", paste0(round(r_preview, digits = 3), collapse = ", "))
    if (n_valid > 5) cat(", ...")
    cat("\n")
  }

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
  return(result)
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
