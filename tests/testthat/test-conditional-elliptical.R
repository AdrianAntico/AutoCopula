testthat::test_that("elliptical conditionals remain on copula scale", {
  out <- qa_autocopula_conditional_elliptical()
  testthat::expect_true(isTRUE(attr(out, "passed")),
    info = paste(out[!passed, detail], collapse = "; "))
})
