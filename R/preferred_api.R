#' Fit one or more copula families
#'
#' @description
#' Preferred expert fitting entry point. Family identity is explicit and the
#' returned object retains fitted models, empirical marginal transforms, source
#' data, and the underlying fitter for advanced inspection.
#'
#' @param data A `data.table` of continuous variables.
#' @param families Character vector of explicit copula families. Use
#'   `ModelFitter$new(data)$list_models()` to inspect the complete built-in set.
#' @param family_definitions Optional named model definitions, each containing a
#'   description and `fit_function(data)`. These extend or override the built-in
#'   library and expose engine-native fitting controls without changing the
#'   compact public API.
#' @return An `autocopula_fit` with fitted family models, marginal transforms,
#'   source data, fit failures, and the underlying `ModelFitter`.
#' @export
copula_fit <- function(data, families = c("Gaussian", "tCopula"),
    family_definitions = list()) {
  data <- data.table::as.data.table(data)
  if (ncol(data) < 2L || nrow(data) < 3L) {
    stop("Copula fitting requires at least two variables and three rows.", call. = FALSE)
  }
  if (any(!vapply(data, is.numeric, logical(1L)))) {
    stop("Copula fitting requires numeric variables.", call. = FALSE)
  }
  fitter <- ModelFitter$new(data.table::copy(data))
  if (length(family_definitions)) {
    if (is.null(names(family_definitions)) || any(!nzchar(names(family_definitions)))) {
      stop("family_definitions must be a named list.", call. = FALSE)
    }
    fitter$copula_library[names(family_definitions)] <- family_definitions
  }
  families <- unique(as.character(families))
  unknown <- setdiff(families, names(fitter$copula_library))
  if (length(unknown)) stop("Unknown copula families: ",
    paste(unknown, collapse = ", "), call. = FALSE)
  fitter$fit_models(families)
  models <- fitter$fit_results[families]
  failures <- families[vapply(models, is.null, logical(1L))]
  out <- list(
    models = models[!vapply(models, is.null, logical(1L))],
    marginals = fitter$fit_results$marginals,
    families_requested = families,
    families_fitted = setdiff(families, failures),
    failures = failures,
    data = data.table::copy(data),
    fitter = fitter,
    call = match.call()
  )
  class(out) <- c("autocopula_fit", "list")
  out
}

autocopula_require_fit <- function(fit, family) {
  if (!inherits(fit, "autocopula_fit")) stop("fit must be returned by copula_fit().", call. = FALSE)
  family <- as.character(family)[1L]
  if (!family %in% names(fit$models)) stop("Family was not fitted: ", family, call. = FALSE)
  family
}

#' Simulate from a fitted copula
#' @param fit An `autocopula_fit`.
#' @param family Fitted family to use.
#' @param n Number of draws.
#' @param seed Reproducibility seed.
#' @param batches Number of batches; values above one use the batch engine.
#' @param parallel Whether batches run in parallel.
#' @param threads Optional worker count.
#' @return A `data.table` on the original marginal scale.
#' @export
copula_simulate <- function(fit, family, n = 1000L, seed = 1L,
    batches = 1L, parallel = FALSE, threads = NULL) {
  family <- autocopula_require_fit(fit, family)
  n <- as.integer(n)[1L]
  batches <- as.integer(batches)[1L]
  if (n < 1L || batches < 1L) stop("n and batches must be positive.", call. = FALSE)
  set.seed(as.integer(seed)[1L])
  scorer <- ModelScorer$new(c(fit$models, list(marginals = fit$marginals)), fit$data)
  if (batches == 1L) return(scorer$batch_prediction(family, n = n))
  scorer$large_scale_simulation(family, batches = batches,
    batch_size = n, parallel = parallel, threads = threads)
}

#' Simulate a fitted conditional copula law
#' @param fit An `autocopula_fit`.
#' @param family Fitted family to use.
#' @param known_ranges Named list of exact or ranged conditioning values on the
#'   original marginal scale.
#' @param n Draws per conditioning combination.
#' @param seed Reproducibility seed.
#' @param parallel Whether combinations run in parallel.
#' @param threads Optional worker count.
#' @return A `data.table` of conditional draws with batch identities.
#' @export
copula_conditional <- function(fit, family, known_ranges, n = 1000L,
    seed = 1L, parallel = FALSE, threads = NULL) {
  family <- autocopula_require_fit(fit, family)
  set.seed(as.integer(seed)[1L])
  scorer <- ModelScorer$new(c(fit$models, list(marginals = fit$marginals)), fit$data)
  scorer$conditional_range_prediction(family, known_ranges = known_ranges,
    n = as.integer(n), parallel = parallel, threads = threads)
}

#' Diagnose fitted copula families
#' @param fit An `autocopula_fit`.
#' @param families Optional fitted families to diagnose.
#' @return A list containing fit coverage, AIC/BIC/log-likelihood and dependence
#'   metrics where supported, and retained fit failures.
#' @export
copula_diagnose <- function(fit, families = names(fit$models)) {
  if (!inherits(fit, "autocopula_fit")) stop("fit must be returned by copula_fit().", call. = FALSE)
  families <- intersect(as.character(families), names(fit$models))
  evaluator <- ModelEvaluation$new(fit$models[families], fit$data)
  list(
    metrics = evaluator$generate_metrics(),
    families_requested = fit$families_requested,
    families_fitted = fit$families_fitted,
    failures = fit$failures,
    fitted_models = fit$models[families]
  )
}
