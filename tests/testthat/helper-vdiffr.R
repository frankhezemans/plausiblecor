# Following code adapted from `ggplot2` package:
# (tidyverse/ggplot2/tests/testthat/helper-vdiffr.R)
# -----------------------------------------------------------------------------
if (
  rlang::is_installed("vdiffr (>= 1.0.0)") &&
  rlang::is_installed("testthat (>= 3.0.0)")
) {
  expect_doppelganger <- vdiffr::expect_doppelganger
} else {
  if (identical(Sys.getenv("VDIFFR_RUN_TESTS"), "true")) {
    rlang::abort("'vdiffr' is not installed")
  }
  expect_doppelganger <- function(...) {
    testthat::skip("'vdiffr' is not installed")
  }
}
