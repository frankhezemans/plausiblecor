testthat::test_that(
  desc = "appropriately validates input arguments",
  code = {
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0L,
        n = 30
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = -1.5,
        n = 30
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 27.3
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 1L
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 30,
        kappa = -1
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 30,
        kappa = 1L
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 30,
        n_bins = 25
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 30,
        n_bins = 112.3
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 30,
        max_iter = 0L
      )
    )
    testthat::expect_error(
      object = plausiblecor::posterior_cor_updf(
        r = 0.3,
        n = 30,
        max_iter = 112.3
      )
    )
  }
)

testthat::test_that(
  desc = "returns a probability density function",
  code = {
    posterior_pdf <- plausiblecor::posterior_cor_updf(
      r = 0.5,
      n = 30
    )
    testthat::expect_true(
      object = inherits(posterior_pdf, "function")
    )
  }
)

testthat::test_that(
  desc = "returned PDF returns numeric vector of finite values",
  code = {
    posterior_pdf <- plausiblecor::posterior_cor_updf(
      r = 0.5,
      n = 30
    )
    posterior_pdf_approx <- plausiblecor::posterior_cor_updf(
      r = 0,
      n = 30
    )
    cor_vals <- seq(from = -0.9, to = 0.9, by = 0.1)
    testthat::expect_vector(
      object = posterior_pdf(cor_vals),
      ptype = numeric(),
      size = length(cor_vals)
    )
    testthat::expect_true(
      object = all(is.finite(posterior_pdf(cor_vals)))
    )
    testthat::expect_vector(
      object = posterior_pdf_approx(cor_vals),
      ptype = numeric(),
      size = length(cor_vals)
    )
    testthat::expect_true(
      object = all(is.finite(posterior_pdf_approx(cor_vals)))
    )
  }
)

testthat::test_that(
  desc = "returned PDF peaks at true mode and returns zero outside [-1, 1]",
  code = {
    posterior_pdf <- plausiblecor::posterior_cor_updf(
      r = 0.236,
      n = 30
    )
    cor_vals <- seq(from = -0.9, to = 0.9, by = 0.1)
    cor_vals <- c(0.236, cor_vals)
    densities <- posterior_pdf(cor_vals)
    testthat::expect_true(
      object = which.max(densities) == 1L
    )
    testthat::expect_identical(
      object = posterior_pdf(c(-1.1, 1.1)),
      expected = numeric(length = 2L)
    )
    posterior_pdf_approx <- plausiblecor::posterior_cor_updf(
      r = 0,
      n = 30
    )
    cor_vals <- c(0, cor_vals)
    densities_approx <- posterior_pdf_approx(cor_vals)
    testthat::expect_true(
      object = which.max(densities_approx) == 1L
    )
    testthat::expect_identical(
      object = posterior_pdf_approx(c(-1.1, 1.1)),
      expected = numeric(length = 2L)
    )
  }
)

testthat::test_that(
  desc = "returned PDF is truncated for directional hypotheses",
  code = {
    pos_pdf <- plausiblecor::posterior_cor_updf(
      r = 0.3,
      n = 30,
      alternative = "greater"
    )
    neg_pdf <- plausiblecor::posterior_cor_updf(
      r = -0.3,
      n = 30,
      alternative = "less"
    )
    testthat::expect_identical(
      object = pos_pdf(-0.3),
      expected = numeric(length = 1L)
    )
    testthat::expect_identical(
      object = neg_pdf(0.3),
      expected = numeric(length = 1L)
    )
  }
)

testthat::test_that(
  desc = "warns and returns degenerate function for perfect sample correlation",
  code = {
    testthat::expect_warning(
      object = plausiblecor::posterior_cor_updf(
        r = 1,
        n = 30
      ),
      regexp = "^Perfect correlation"
    )
    testthat::expect_identical(
      object = plausiblecor::posterior_cor_updf(
        r = 1,
        n = 30
      )(c(-0.5, 0, 0.5)),
      expected = rep(NA_real_, times = 3L)
    )
  }
)

testthat::test_that(
  desc = "warns when Kendall's tau method may be inappropriate",
  code = {
    testthat::expect_warning(
      object = plausiblecor::posterior_cor_updf(
        r = 0.95,
        n = 10,
        method = "kendall"
      ),
      regexp = "^Normal approximation"
    )
  }
)
