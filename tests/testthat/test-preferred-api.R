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

testthat::test_that("expert API exposes vines and hides R6 generators", {
  catalog <- copula_families()
  testthat::expect_s3_class(catalog, "data.table")
  testthat::expect_true("RVine" %in% catalog$family)
  testthat::expect_true(all(c("vine", "CVine", "DVine") %in% catalog$family))
  fit_formals <- names(formals(copula_fit))
  testthat::expect_equal(fit_formals[1:2], c("data", "families"))
  testthat::expect_true(all(c("family_set", "trunc_lvl", "treecrit", "vine_type",
    "family_definitions") %in% fit_formals))
  curated <- c("copula_fit", "copula_simulate", "copula_conditional",
    "copula_diagnose", "copula_families")
  exports <- getNamespaceExports("AutoCopula")
  testthat::expect_true(all(curated %in% exports))
  if (length(exports) <= 10L) {
    testthat::expect_false("ModelFitter" %in% exports)
    testthat::expect_false("EDA" %in% exports)
    testthat::expect_false("qa_autocopula_conditional_elliptical" %in% exports)
    testthat::expect_setequal(exports, curated)
  }
  ns_path <- system.file("NAMESPACE", package = "AutoCopula")
  ns_exports <- sub("^export\\((.*)\\)$", "\\1",
    grep("^export\\(", readLines(ns_path), value = TRUE))
  testthat::expect_setequal(ns_exports, curated)
  testthat::expect_false("ModelFitter" %in% ns_exports)
})

testthat::test_that("vine families fit, simulate, condition, and diagnose", {
  testthat::skip_if_not_installed("VineCopula")
  set.seed(7)
  z <- matrix(stats::rnorm(900), ncol = 3L)
  data <- data.table::data.table(
    a = z[, 1L],
    b = 0.5 * z[, 1L] + sqrt(1 - 0.5^2) * z[, 2L],
    c = 0.3 * z[, 1L] + 0.3 * z[, 2L] + sqrt(1 - 0.3^2 - 0.3^2) * z[, 3L]
  )
  fit <- copula_fit(data, families = "RVine", family_set = c(1L, 5L),
    trunc_lvl = 1L, treecrit = "tau", vine_type = "rvine")
  testthat::expect_s3_class(fit, "autocopula_fit")
  testthat::expect_true("RVine" %in% fit$families_fitted)
  testthat::expect_true(inherits(fit$models$RVine, "autocopula_vine"))
  testthat::expect_true(inherits(fit$models$RVine$rvine, "RVineMatrix"))

  draws <- copula_simulate(fit, "RVine", n = 40L, seed = 8L)
  testthat::expect_s3_class(draws, "data.table")
  testthat::expect_equal(nrow(draws), 40L)
  testthat::expect_named(draws, names(data))

  conditional <- copula_conditional(fit, "RVine",
    known_ranges = list(a = median(data$a)), n = 25L, seed = 9L)
  testthat::expect_s3_class(conditional, "data.table")
  testthat::expect_true(nrow(conditional) > 0L)
  testthat::expect_true(!is.null(attr(conditional, "conditional_provenance")))

  diagnosed <- copula_diagnose(fit)
  testthat::expect_true("RVine" %in% names(diagnosed$vines))
  testthat::expect_true(is.data.frame(diagnosed$vines$RVine$pair_copulas) ||
    data.table::is.data.table(diagnosed$vines$RVine$pair_copulas))
  testthat::expect_true(grepl("nonexchangeable", diagnosed$vines$RVine$exchangeability,
    ignore.case = TRUE))
})

