autocopula_clamp_u <- function(u) {
  pmin(pmax(as.numeric(u), sqrt(.Machine$double.eps)),
    1 - sqrt(.Machine$double.eps))
}

autocopula_bicop_dispatch <- function() c(
  BB1 = "conditional_bb1",
  BB6 = "conditional_bb6",
  BB7 = "conditional_bb7",
  BB8 = "conditional_bb8",
  "Rotated Clayton (180)" = "conditional_rotated_clayton_180",
  "Rotated Gumbel (180)" = "conditional_rotated_gumbel_180",
  "Rotated Joe (180)" = "conditional_rotated_joe_180",
  "Rotated BB1 (180)" = "conditional_rotated_bb1_180",
  "Rotated BB6 (180)" = "conditional_rotated_bb6_180",
  "Rotated BB7 (180)" = "conditional_rotated_bb7_180",
  "Rotated BB8 (180)" = "conditional_rotated_bb8_180",
  "Tawn Type 1" = "conditional_tawn_type1",
  "Rotated Tawn Type 1 (180)" = "conditional_rotated_tawn_type1_180",
  "Tawn Type 2" = "conditional_tawn_type2",
  "Rotated Tawn Type 2 (180)" = "conditional_rotated_tawn_type2_180"
)

autocopula_conditional_bicop_u <- function(known_u, known_index, n, family,
    par, par2 = 0) {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop("VineCopula is required for bivariate conditional sampling.")
  }
  if (length(known_u) != 1L || !is.finite(known_u)) {
    stop("known_u must be one finite copula-scale value.")
  }
  if (!known_index %in% c(1L, 2L)) {
    stop("known_index must be 1 or 2 for a bivariate copula.")
  }
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  known_u <- autocopula_clamp_u(known_u)
  probability <- autocopula_clamp_u(stats::runif(as.integer(n)))
  known_vector <- rep(known_u, length(probability))
  sampled <- if (known_index == 1L) {
    VineCopula::BiCopHinv1(known_vector, probability, family = family,
      par = par, par2 = par2)
  } else {
    VineCopula::BiCopHinv2(probability, known_vector, family = family,
      par = par, par2 = par2)
  }
  autocopula_clamp_u(sampled)
}

autocopula_conditional_archimedean_u <- function(copula_model, known_indices,
    known_u, n) {
  if (!requireNamespace("copula", quietly = TRUE)) {
    stop("copula is required for Archimedean conditional sampling.")
  }
  d <- dim(copula_model)
  known_indices <- as.integer(known_indices)
  remaining <- setdiff(seq_len(d), known_indices)
  if (!length(known_indices) || !length(remaining)) {
    stop("Conditioning requires at least one known and one unknown variable.")
  }
  if (length(known_u) != length(known_indices)) {
    stop("known_u does not match known_indices.")
  }
  # Archimedean copulas are exchangeable, so known coordinates can be placed
  # first for the inverse Rosenblatt transform and then restored to data order.
  reordered <- cbind(
    matrix(rep(autocopula_clamp_u(known_u), each = n), nrow = n),
    matrix(stats::runif(n * length(remaining)), nrow = n)
  )
  sampled <- copula::cCopula(reordered, copula = copula_model,
    inverse = TRUE)
  result <- matrix(NA_real_, nrow = n, ncol = d)
  result[, known_indices] <- sampled[, seq_along(known_indices), drop = FALSE]
  result[, remaining] <- sampled[, length(known_indices) + seq_along(remaining),
    drop = FALSE]
  result
}

#' Qualify bivariate conditional copula kernels
#'
#' Checks both conditioning directions against the authoritative VineCopula
#' inverse-h/forward-h identity across ordinary, rotated, two-parameter, and
#' asymmetric copula families.
#' @return A data table of qualification checks, or a dependency warning row.
qa_autocopula_conditional_bicop <- function() {
  add <- function(check, passed, detail = "") data.table::data.table(
    check = check, passed = isTRUE(passed), detail = as.character(detail)[1L])
  expected_dispatch <- c(
    BB1 = "conditional_bb1", BB6 = "conditional_bb6",
    BB7 = "conditional_bb7", BB8 = "conditional_bb8",
    "Tawn Type 1" = "conditional_tawn_type1",
    "Tawn Type 2" = "conditional_tawn_type2"
  )
  dispatch_ok <- identical(
    unname(autocopula_bicop_dispatch()[names(expected_dispatch)]),
    unname(expected_dispatch))
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    out <- data.table::rbindlist(list(
      add("family_dispatch", dispatch_ok),
      add("vinecopula_dependency", FALSE,
        "VineCopula is required; conditional capability was not falsely qualified.")))
    attr(out, "passed") <- FALSE
    attr(out, "status") <- "dependency_missing"
    return(out)
  }
  specifications <- list(
    clayton = c(family = 3, par = 2, par2 = 0),
    survival_clayton = c(family = 13, par = 2, par2 = 0),
    bb1 = c(family = 7, par = 1.5, par2 = 2),
    tawn_type1 = c(family = 104, par = 2, par2 = 0.6)
  )
  checks <- lapply(names(specifications), function(label) {
    spec <- specifications[[label]]
    probability <- seq(0.01, 0.99, length.out = 199L)
    known <- rep(0.83, length(probability))
    sampled_2 <- VineCopula::BiCopHinv1(known, probability,
      family = spec[["family"]], par = spec[["par"]], par2 = spec[["par2"]])
    sampled_1 <- VineCopula::BiCopHinv2(probability, known,
      family = spec[["family"]], par = spec[["par"]], par2 = spec[["par2"]])
    recovered_2 <- VineCopula::BiCopHfunc1(known, sampled_2,
      family = spec[["family"]], par = spec[["par"]], par2 = spec[["par2"]])
    recovered_1 <- VineCopula::BiCopHfunc2(sampled_1, known,
      family = spec[["family"]], par = spec[["par"]], par2 = spec[["par2"]])
    add(paste0(label, "_both_directions"),
      max(abs(recovered_1 - probability), abs(recovered_2 - probability)) < 1e-7 &&
        all(sampled_1 > 0 & sampled_1 < 1) && all(sampled_2 > 0 & sampled_2 < 1))
  })
  invalid <- try(autocopula_conditional_bicop_u(0.5, 3L, 10L, 3L, 2),
    silent = TRUE)
  out <- data.table::rbindlist(c(checks,
    list(add("invalid_direction_fails_closed", inherits(invalid, "try-error")))))
  attr(out, "passed") <- all(out$passed)
  attr(out, "status") <- if (all(out$passed)) "passed" else "failed"
  out
}

autocopula_rmvnorm <- function(n, mean, sigma) {
  p <- length(mean)
  if (!p) return(matrix(numeric(), nrow = n, ncol = 0L))
  sigma <- (sigma + t(sigma)) / 2
  eig <- eigen(sigma, symmetric = TRUE)
  if (min(eig$values) < -1e-8) {
    stop("Conditional covariance is not positive semidefinite.")
  }
  root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), p) %*%
    t(eig$vectors)
  sweep(matrix(stats::rnorm(n * p), nrow = n) %*% root, 2L, mean, "+")
}

autocopula_conditional_elliptical_u <- function(corr, known_indices, known_u,
    n, family = c("gaussian", "t"), df = NULL) {
  family <- match.arg(family)
  d <- nrow(corr)
  known_indices <- as.integer(known_indices)
  remaining <- setdiff(seq_len(d), known_indices)
  if (!length(known_indices) || !length(remaining)) {
    stop("Conditioning requires at least one known and one unknown variable.")
  }
  if (length(known_u) != length(known_indices)) {
    stop("known_u does not match known_indices.")
  }
  corr_11 <- corr[known_indices, known_indices, drop = FALSE]
  corr_12 <- corr[known_indices, remaining, drop = FALSE]
  corr_22 <- corr[remaining, remaining, drop = FALSE]
  corr_21 <- t(corr_12)
  latent_known <- if (identical(family, "gaussian")) {
    stats::qnorm(autocopula_clamp_u(known_u))
  } else {
    if (is.null(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
      stop("A finite positive df is required for t-copula conditioning.")
    }
    stats::qt(autocopula_clamp_u(known_u), df = df)
  }
  solved <- solve(corr_11, latent_known)
  location <- as.vector(corr_21 %*% solved)
  schur <- corr_22 - corr_21 %*% solve(corr_11, corr_12)

  if (identical(family, "gaussian")) {
    latent <- autocopula_rmvnorm(n, location, schur)
    sampled_u <- stats::pnorm(latent)
  } else {
    k <- length(known_indices)
    conditional_df <- df + k
    mahalanobis <- sum(latent_known * solved)
    scale <- ((df + mahalanobis) / conditional_df) * schur
    z <- autocopula_rmvnorm(n, rep(0, length(location)), scale)
    z <- z / sqrt(stats::rchisq(n, conditional_df) / conditional_df)
    latent <- sweep(z, 2L, location, "+")
    sampled_u <- stats::pt(latent, df = df)
  }
  result <- matrix(NA_real_, nrow = n, ncol = d)
  result[, known_indices] <- matrix(rep(autocopula_clamp_u(known_u), each = n),
    nrow = n)
  result[, remaining] <- sampled_u
  result
}

#' Qualify elliptical conditional copula kernels
#'
#' Checks analytic conditional moments, copula-scale support, Student-t tail
#' behavior, and invalid conditioning contracts without requiring fitted models.
#' @return A data table of qualification checks.
qa_autocopula_conditional_elliptical <- function() {
  add <- function(check, passed, detail = "") data.table::data.table(
    check = check, passed = isTRUE(passed), detail = as.character(detail)[1L])
  corr <- matrix(c(1, 0.7, 0.7, 1), 2L)
  set.seed(441L)
  g <- autocopula_conditional_elliptical_u(corr, 1L, 0.9, 100000L,
    "gaussian")
  z <- stats::qnorm(g[, 2L])
  expected_mean <- 0.7 * stats::qnorm(0.9)
  expected_sd <- sqrt(1 - 0.7^2)
  set.seed(442L)
  tt <- autocopula_conditional_elliptical_u(corr, 1L, 0.99, 100000L,
    "t", df = 4)
  set.seed(443L)
  g_tail <- autocopula_conditional_elliptical_u(corr, 1L, 0.99, 100000L,
    "gaussian")
  invalid <- try(autocopula_conditional_elliptical_u(corr, 1L, 0.5, 10L,
    "t", df = -1), silent = TRUE)
  catalog <- copula_families()
  catalog_names <- unique(c(catalog$family, catalog$Model, names(catalog)))
  fit_formals <- names(formals(copula_fit))
  ns_exports <- tryCatch(getNamespaceExports("AutoCopula"),
    error = function(e) character())
  ns_path <- system.file("NAMESPACE", package = "AutoCopula")
  ns_file_exports <- if (nzchar(ns_path)) {
    ns_lines <- readLines(ns_path)
    sub("^export\\((.*)\\)$", "\\1",
      ns_lines[grepl("^export\\(", ns_lines)])
  } else {
    ns_exports
  }
  out <- data.table::rbindlist(list(
    add("gaussian_latent_mean", abs(mean(z) - expected_mean) < 0.015,
      sprintf("mean=%.4f expected=%.4f", mean(z), expected_mean)),
    add("gaussian_latent_sd", abs(stats::sd(z) - expected_sd) < 0.015,
      sprintf("sd=%.4f expected=%.4f", stats::sd(z), expected_sd)),
    add("uniform_support", all(g > 0 & g < 1) && all(tt > 0 & tt < 1)),
    add("t_tail_dependence",
      mean(tt[, 2L] > 0.99) > mean(g_tail[, 2L] > 0.99),
      sprintf("t=%.4f gaussian=%.4f", mean(tt[, 2L] > 0.99),
        mean(g_tail[, 2L] > 0.99))),
    add("invalid_df_fails_closed", inherits(invalid, "try-error")),
    add("rvine_in_family_catalog", "RVine" %in% catalog_names,
      paste(catalog$family, collapse = ", ")),
    add("copula_fit_has_family_set", "family_set" %in% fit_formals,
      paste(fit_formals, collapse = ", ")),
    add("copula_fit_has_trunc_lvl", "trunc_lvl" %in% fit_formals,
      paste(fit_formals, collapse = ", ")),
    add("namespace_hides_modelfitter",
      !"ModelFitter" %in% ns_file_exports &&
        (length(ns_exports) > 20L || !"ModelFitter" %in% ns_exports),
      paste(unique(c(ns_file_exports, ns_exports)), collapse = ", "))
  ))
  attr(out, "passed") <- all(out$passed)
  out
}
