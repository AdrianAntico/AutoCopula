testthat::test_that("preferred API preserves fitted family identity", {
  set.seed(42)
  z <- matrix(stats::rnorm(800), ncol = 2L)
  data <- data.table::data.table(a = z[, 1L],
    b = 0.6 * z[, 1L] + sqrt(1 - 0.6^2) * z[, 2L])
  fit <- copula_fit(data, families = c("Gaussian", "tCopula"))
  testthat::expect_s3_class(fit, "autocopula_fit")
  testthat::expect_setequal(fit$families_fitted, c("Gaussian", "tCopula"))

  draws <- copula_simulate(fit, "Gaussian", n = 100L, seed = 3L)
  testthat::expect_s3_class(draws, "data.table")
  testthat::expect_equal(nrow(draws), 100L)
  testthat::expect_named(draws, names(data))

  conditional <- copula_conditional(fit, "Gaussian",
    known_ranges = list(a = c(-0.5, 0.5)), n = 40L, seed = 4L)
  testthat::expect_s3_class(conditional, "data.table")
  testthat::expect_true(nrow(conditional) > 0L)

  diagnosed <- copula_diagnose(fit)
  testthat::expect_true(all(c("metrics", "families_fitted", "fitted_models") %in%
    names(diagnosed)))
})
