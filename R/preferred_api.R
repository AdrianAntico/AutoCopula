#' Fit one or more copula families
#'
#' @description
#' Preferred expert fitting entry point. Family identity is explicit and the
#' returned object retains fitted models, empirical marginal transforms, source
#' data, and the underlying fitter for advanced inspection.
#'
#' Vines (`"RVine"`, `"vine"`, `"CVine"`, `"DVine"`) are first-class families.
#' They are the multivariate nonexchangeable path: each pair (and conditional
#' pair) may have its own copula family and tail behavior. Clayton, Gumbel,
#' Frank, and Joe fits from the copula package remain exchangeable for `d > 2`;
#' experts who need pair-specific tails should request a vine.
#'
#' @param data A `data.table` of continuous variables.
#' @param families Character vector of explicit copula families. Use
#'   [copula_families()] to inspect the complete built-in set, including vines.
#' @param family_set VineCopula pair-copula family codes or names used when a
#'   vine family is requested. The default is a serious set: Gaussian, t,
#'   Clayton, Gumbel, Frank, Joe, BB1, BB7, and their survival (180-degree)
#'   rotations. It is not Gaussian+t only.
#' @param trunc_lvl Vine truncation level passed to VineCopula. `NA` fits all
#'   trees.
#' @param treecrit Tree-selection criterion for R-/C-vine structure selection
#'   (`"tau"` by default).
#' @param vine_type Structure used when `families` contains the generic
#'   `"vine"` alias: `"rvine"` (default), `"cvine"`, or `"dvine"`. `"RVine"`,
#'   `"CVine"`, and `"DVine"` family names override this.
#' @param family_definitions Optional named model definitions, each containing a
#'   description and `fit_function(data)`. These extend or override the built-in
#'   library and expose engine-native fitting controls without changing the
#'   compact public API.
#' @return An `autocopula_fit` with fitted family models, marginal transforms,
#'   source data, fit failures, and the underlying `ModelFitter`. Vine fits
#'   store the `RVineMatrix` on the model object.
#' @export
copula_fit <- function(data, families = c("Gaussian", "tCopula"),
    family_set = NULL, trunc_lvl = NA_integer_, treecrit = "tau",
    vine_type = "rvine", family_definitions = list()) {
  data <- data.table::as.data.table(data)
  if (ncol(data) < 2L || nrow(data) < 3L) {
    stop("Copula fitting requires at least two variables and three rows.", call. = FALSE)
  }
  if (any(!vapply(data, is.numeric, logical(1L)))) {
    stop("Copula fitting requires numeric variables.", call. = FALSE)
  }
  fitter <- ModelFitter$new(data.table::copy(data))
  families <- unique(as.character(families))
  if (any(vapply(families, autocopula_is_vine_family_name, logical(1L)))) {
    autocopula_require_vinecopula("vine copula fitting")
  }
  autocopula_register_vine_families(fitter, families, vine_type, family_set,
    trunc_lvl, treecrit)
  if (length(family_definitions)) {
    if (is.null(names(family_definitions)) || any(!nzchar(names(family_definitions)))) {
      stop("family_definitions must be a named list.", call. = FALSE)
    }
    fitter$copula_library[names(family_definitions)] <- family_definitions
  }
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
    vine_controls = list(family_set = autocopula_resolve_family_set(family_set),
      trunc_lvl = autocopula_trunclevel(trunc_lvl),
      treecrit = as.character(treecrit)[1L],
      vine_type = autocopula_resolve_vine_type(vine_type)),
    call = match.call()
  )
  class(out) <- c("autocopula_fit", "list")
  out
}

#' Built-in copula family catalog
#'
#' Specialized expert helper. Returns the package family library, including
#' first-class vines, without constructing `ModelFitter`.
#'
#' Archimedean families Clayton, Gumbel, Frank, and Joe are exchangeable for
#' `d > 2` in the copula package. Vines (`RVine`, `CVine`, `DVine`) are the
#' nonexchangeable multivariate path with pair-specific tails.
#'
#' @return A `data.table` with family names, descriptions, model class,
#'   dimension, exchangeability, and expert notes.
#' @export
copula_families <- function() {
  autocopula_builtin_family_rows()
}

autocopula_require_fit <- function(fit, family) {
  if (!inherits(fit, "autocopula_fit")) stop("fit must be returned by copula_fit().", call. = FALSE)
  family <- as.character(family)[1L]
  if (!family %in% names(fit$models)) stop("Family was not fitted: ", family, call. = FALSE)
  family
}

autocopula_vine_simulate_dt <- function(fit, family, n) {
  model <- fit$models[[family]]
  u <- autocopula_vine_sim_u(model, n)
  autocopula_back_transform_u(u, fit$marginals)
}

#' Simulate from a fitted copula
#' @param fit An `autocopula_fit`.
#' @param family Fitted family to use.
#' @param n Number of draws.
#' @param seed Reproducibility seed.
#' @param batches Number of batches; values above one use the batch engine.
#' @param parallel Whether batches run in parallel.
#' @param threads Optional worker count.
#' @return A `data.table` on the original marginal scale. Vine families are
#'   simulated with `VineCopula::RVineSim` from the stored `RVineMatrix`.
#' @export
copula_simulate <- function(fit, family, n = 1000L, seed = 1L,
    batches = 1L, parallel = FALSE, threads = NULL) {
  family <- autocopula_require_fit(fit, family)
  n <- as.integer(n)[1L]
  batches <- as.integer(batches)[1L]
  if (n < 1L || batches < 1L) stop("n and batches must be positive.", call. = FALSE)
  set.seed(as.integer(seed)[1L])
  if (autocopula_is_vine_model(fit$models[[family]])) {
    if (batches == 1L) return(autocopula_vine_simulate_dt(fit, family, n))
    batch_fun <- function(batch_id) {
      draws <- autocopula_vine_simulate_dt(fit, family, n)
      draws[, BatchID := batch_id]
      draws
    }
    results <- if (isTRUE(parallel)) {
      num_threads <- if (is.null(threads)) future::availableCores() - 1L else
        min(as.integer(threads), future::availableCores())
      future::plan(future::multisession, workers = max(1L, num_threads))
      on.exit(future::plan("sequential"), add = TRUE)
      future.apply::future_lapply(seq_len(batches), batch_fun, future.seed = TRUE)
    } else {
      lapply(seq_len(batches), batch_fun)
    }
    return(data.table::rbindlist(results, use.names = TRUE, fill = TRUE))
  }
  scorer <- ModelScorer$new(c(fit$models, list(marginals = fit$marginals)), fit$data)
  if (batches == 1L) return(scorer$batch_prediction(family, n = n))
  scorer$large_scale_simulation(family, batches = batches,
    batch_size = n, parallel = parallel, threads = threads)
}

autocopula_vine_conditional_dt <- function(fit, family, known_ranges, n) {
  if (!is.list(known_ranges) || !length(known_ranges) ||
      is.null(names(known_ranges)) || any(!nzchar(names(known_ranges)))) {
    stop("known_ranges must be a non-empty named list.", call. = FALSE)
  }
  if (!all(names(known_ranges) %in% names(fit$data))) {
    stop("known_ranges names must match fitted variables.", call. = FALSE)
  }
  combinations <- do.call(expand.grid, known_ranges)
  model <- fit$models[[family]]
  var_names <- names(fit$data)
  batches <- lapply(seq_len(nrow(combinations)), function(i) {
    known_values <- as.list(combinations[i, , drop = FALSE])
    known_vars <- names(known_values)
    known_indices <- match(known_vars, var_names)
    known_u <- vapply(known_vars, function(col) {
      autocopula_x_to_u(fit$data[[col]], known_values[[col]])
    }, numeric(1L))
    sampled <- autocopula_vine_conditional_u(model, known_indices, known_u,
      as.integer(n)[1L])
    dt <- autocopula_back_transform_u(sampled, fit$marginals)
    dt[, batch_id := i]
    list(dt = dt, provenance = attr(sampled, "conditional_provenance"))
  })
  out <- data.table::rbindlist(lapply(batches, `[[`, "dt"), use.names = TRUE,
    fill = TRUE)
  attr(out, "conditional_provenance") <- lapply(batches, `[[`, "provenance")
  out
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
#' @return A `data.table` of conditional draws with batch identities. Vine
#'   conditionals use inverse h-functions or `RVineCondSim` when they apply to
#'   the fitted structure; otherwise Monte Carlo kernel SIR from
#'   `VineCopula::RVineSim` with recorded provenance. Analytic conditionals are
#'   never faked.
#' @export
copula_conditional <- function(fit, family, known_ranges, n = 1000L,
    seed = 1L, parallel = FALSE, threads = NULL) {
  family <- autocopula_require_fit(fit, family)
  set.seed(as.integer(seed)[1L])
  if (autocopula_is_vine_model(fit$models[[family]])) {
    return(autocopula_vine_conditional_dt(fit, family, known_ranges, n))
  }
  scorer <- ModelScorer$new(c(fit$models, list(marginals = fit$marginals)), fit$data)
  scorer$conditional_range_prediction(family, known_ranges = known_ranges,
    n = as.integer(n), parallel = parallel, threads = threads)
}

#' Diagnose fitted copula families
#' @param fit An `autocopula_fit`.
#' @param families Optional fitted families to diagnose.
#' @return A list containing fit coverage, AIC/BIC/log-likelihood and dependence
#'   metrics where supported, vine pair-family tables and truncation notes, and
#'   retained fit failures.
#' @export
copula_diagnose <- function(fit, families = names(fit$models)) {
  if (!inherits(fit, "autocopula_fit")) stop("fit must be returned by copula_fit().", call. = FALSE)
  families <- intersect(as.character(families), names(fit$models))
  vine_names <- families[vapply(fit$models[families], autocopula_is_vine_model,
    logical(1L))]
  other_names <- setdiff(families, vine_names)
  metrics <- if (length(other_names)) {
    ModelEvaluation$new(fit$models[other_names], fit$data)$generate_metrics()
  } else {
    data.table::data.table()
  }
  vines <- lapply(vine_names, function(nm) autocopula_vine_diagnose(fit, nm))
  names(vines) <- vine_names
  list(
    metrics = metrics,
    vines = vines,
    families_requested = fit$families_requested,
    families_fitted = fit$families_fitted,
    failures = fit$failures,
    fitted_models = fit$models[families]
  )
}
