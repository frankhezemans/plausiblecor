testthat::test_that(
  desc = "appropriately validates input arguments",
  code = {
    testthat::expect_error(
      object = plausiblecor::run_plausible_cor(
        parameter = "Lorum_ipsum",
        covariate = "striatum",
        mcmc_data = plausiblecor::Forstmann_LBA,
        covariate_data = plausiblecor::Forstmann_fMRI
      ),
      regexp = "^Assertion on 'mcmc_data' failed:"
    )
    testthat::expect_error(
      object = plausiblecor::run_plausible_cor(
        parameter = "caution_effect_speed",
        covariate = "Lorum_ipsum",
        mcmc_data = plausiblecor::Forstmann_LBA,
        covariate_data = plausiblecor::Forstmann_fMRI
      ),
      regexp = "^Covariate 'Lorum_ipsum' not found"
    )
    testthat::expect_error(
      object = plausiblecor::run_plausible_cor(
        parameter = "caution_effect_speed",
        covariate = "striatum",
        mcmc_data = plausiblecor::Forstmann_LBA
      ),
      regexp = "^Covariate 'striatum' not found in 'mcmc_data'"
    )
    testthat::expect_error(
      object = plausiblecor::run_plausible_cor(
        parameter = "caution_effect_speed",
        covariate = "striatum",
        mcmc_data = as.matrix(plausiblecor::Forstmann_LBA),
        covariate_data = plausiblecor::Forstmann_fMRI
      ),
      regexp = "^Assertion on 'mcmc_data' failed:"
    )
    testthat::expect_error(
      object = plausiblecor::run_plausible_cor(
        parameter = "caution_effect_speed",
        covariate = "striatum",
        mcmc_data = plausiblecor::Forstmann_LBA,
        covariate_data = plausiblecor::Forstmann_fMRI,
        confounders = "pre_sma",
        method = "kendall"
      ),
      regexp = "^Partial correlation is not implemented for Kendall's"
    )
  }
)

testthat::test_that(
  desc = "returns 'plausible_cor' object that can be summarised",
  code = {
    out <- plausiblecor::run_plausible_cor(
      parameter = "caution_effect_speed",
      covariate = "striatum",
      mcmc_data = plausiblecor::Forstmann_LBA,
      covariate_data = plausiblecor::Forstmann_fMRI,
      n_draws = 5L,
      rng_seed = 123L
    )
    testthat::expect_s3_class(
      object = out,
      class = c("plausible_cor", "tbl_df", "tbl", "data.frame"),
      exact = TRUE
    )
    testthat::expect_named(
      object = out,
      expected = c(".draw", "r", "n", "k", "posterior_updf")
    )
    testthat::expect_identical(
      object = nrow(out),
      expected = 5L
    )
    testthat::expect_identical(
      object = out[[".draw"]],
      expected = c(526L, 2227L, 2463L, 2511L, 4291L)
    )
    testthat::expect_identical(
      object = round(out[["r"]], digits = 3),
      expected = c(-0.377, -0.125, -0.488, -0.335, -0.326)
    )
    testthat::expect_identical(
      object = unique(out[["n"]]),
      expected = 19L
    )
    testthat::expect_identical(
      object = unique(out[["k"]]),
      expected = 0L
    )
    summary_out <- summarise_plausible_cor(
      out,
      point_interval_args = list(interval_width = c(0.8, 0.95)),
      rope_range = c(-0.1, 0.1)
    )
    testthat::expect_s3_class(
      object = summary_out,
      class = c("tbl_df", "tbl", "data.frame"),
      exact = TRUE
    )
    testthat::expect_named(
      object = summary_out,
      expected = c("type", "mean", "lower", "upper", "width", "p_dir", "p_rope")
    )
    testthat::expect_identical(
      object = nrow(summary_out),
      expected = 4L
    )
    testthat::expect_identical(
      object = unique(summary_out[["type"]]),
      expected = factor(c("sample", "population"))
    )
    testthat::expect_identical(
      object = round(unique(summary_out[["mean"]]), digits = 3),
      expected = c(-0.330, -0.281)
    )
    testthat::expect_identical(
      object = round(summary_out[["lower"]], digits = 3),
      expected = c(-0.488, -0.488, -0.595, -0.695)
    )
    testthat::expect_identical(
      object = round(summary_out[["upper"]], digits = 3),
      expected = c(-0.225, -0.125, -0.017, 0.166)
    )
    testthat::expect_identical(
      object = unique(summary_out[["width"]]),
      expected = c(0.8, 0.95)
    )
    testthat::expect_identical(
      object = round(unique(summary_out[["p_dir"]]), digits = 3),
      expected = c(1, 0.882)
    )
    testthat::expect_identical(
      object = round(unique(summary_out[["p_rope"]]), digits = 3),
      expected = c(0, 0.153)
    )
  }
)

testthat::test_that(
  desc = "returns full output with expected S3 methods",
  code = {
    testthat::skip_on_cran()
    out <- plausiblecor::run_plausible_cor(
      parameter = "caution_effect_speed",
      covariate = "striatum",
      mcmc_data = plausiblecor::Forstmann_LBA,
      covariate_data = plausiblecor::Forstmann_fMRI
    )
    summary_out <- summary(out)
    population_plot <- plot(out, type = "population")
    sample_plot <- plot(out, type = "sample")
    sample_plot_hist <- plot(out, type = "sample", style = "hist")
    testthat::expect_snapshot(out)
    testthat::expect_snapshot(print(out))
    testthat::expect_snapshot(summary_out)
    testthat::expect_snapshot(print(summary_out))
    expect_doppelganger(
      title = "population-level distribution",
      fig = population_plot
    )
    expect_doppelganger(
      title = "sample-level dotplot",
      fig = sample_plot
    )
    expect_doppelganger(
      title = "sample-level histogram",
      fig = sample_plot_hist
    )
  }
)

testthat::test_that(
  desc = "returns 'plausible_cor' object for plausible partial correlation",
  code = {
    out <- plausiblecor::run_plausible_cor(
      parameter = "caution_effect_speed",
      covariate = "striatum",
      mcmc_data = plausiblecor::Forstmann_LBA,
      covariate_data = plausiblecor::Forstmann_fMRI,
      confounders = "pre_sma",
      n_draws = 5L,
      rng_seed = 123L
    )
    testthat::expect_identical(
      object = unique(out[["k"]]),
      expected = 1L
    )
  }
)

testthat::test_that(
  desc = "performs comparison of 'plausible_cor' objects",
  code = {
    out_a <- plausiblecor::run_plausible_cor(
      parameter = "caution_effect_speed",
      covariate = "striatum",
      mcmc_data = plausiblecor::Forstmann_LBA,
      covariate_data = plausiblecor::Forstmann_fMRI,
      n_draws = 5L,
      rng_seed = 123L
    )
    out_b <- plausiblecor::run_plausible_cor(
      parameter = "caution_effect_speed",
      covariate = "pre_sma",
      mcmc_data = plausiblecor::Forstmann_LBA,
      covariate_data = plausiblecor::Forstmann_fMRI,
      n_draws = 5L,
      rng_seed = 123L
    )
    out <- plausiblecor::compare_plausible_cors(
      x = out_a,
      y = out_b,
      rng_seed = c(123L, 456L)
    )
    testthat::expect_s3_class(
      object = out,
      class = c("tbl_df", "tbl", "data.frame"),
      exact = TRUE
    )
    testthat::expect_named(
      object = out,
      expected = c("type", "mean", "lower", "upper", "p_dir")
    )
    testthat::expect_identical(
      object = out[["type"]],
      expected = factor(c("sample", "population"))
    )
    testthat::expect_identical(
      object = round(out[["mean"]], digits = 3),
      expected = c(0.077, 0.200)
    )
    testthat::expect_identical(
      object = round(out[["lower"]], digits = 3),
      expected = c(0.055, 0.069)
    )
    testthat::expect_identical(
      object = round(out[["upper"]], digits = 3),
      expected = c(0.095, 0.310)
    )
    testthat::expect_identical(
      object = round(out[["p_dir"]], digits = 3),
      expected = c(1, 1)
    )
  }
)

testthat::test_that(
  desc = "non-default methods",
  code = {
    out_a <- plausiblecor::run_plausible_cor(
      parameter = "caution_effect_speed",
      covariate = "striatum",
      mcmc_data = plausiblecor::Forstmann_LBA,
      covariate_data = plausiblecor::Forstmann_fMRI,
      alternative = "less",
      method = "kendall",
      n_draws = 5L,
      rng_seed = 123L
    )
    out <- summary(
      out_a,
      point_interval_args = list(
        point_method = "median",
        interval_method = "qi",
        interval_width = c(0.8, 0.95)
      ),
      rope_range = c(-0.1, 0.1)
    )
    testthat::expect_named(
      object = out,
      expected = c("type", "median", "lower", "upper", "width", "p_dir", "p_rope")
    )
    testthat::expect_identical(
      object = round(unique(out[["median"]]), digits = 3),
      expected = c(-0.196, -0.213)
    )
    testthat::expect_identical(
      object = round(out[["lower"]], digits = 3),
      expected = c(-0.262, -0.277, -0.404, -0.509)
    )
    testthat::expect_identical(
      object = round(out[["upper"]], digits = 3),
      expected = c(-0.171, -0.171, -0.059, -0.017)
    )
    testthat::expect_identical(
      object = unique(out[["p_dir"]]),
      expected = NA_real_
    )
    testthat::expect_identical(
      object = round(unique(out[["p_rope"]]), digits = 3),
      expected = c(0, 0.193)
    )
  }
)
