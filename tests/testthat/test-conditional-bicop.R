testthat::test_that("bivariate conditional laws invert authoritative h-functions", {
  out <- AutoCopula:::qa_autocopula_conditional_bicop()
  if (identical(attr(out, "status"), "dependency_missing")) {
    testthat::skip("VineCopula is not installed in this qualification runtime.")
  }
  testthat::expect_true(isTRUE(attr(out, "passed")))
  testthat::expect_true(all(out$passed))
})
