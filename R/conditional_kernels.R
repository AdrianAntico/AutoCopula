autocopula_clamp_u <- function(u) {
  pmin(pmax(as.numeric(u), sqrt(.Machine$double.eps)),
    1 - sqrt(.Machine$double.eps))
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
#' @export
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
    add("invalid_df_fails_closed", inherits(invalid, "try-error"))
  ))
  attr(out, "passed") <- all(out$passed)
  out
}
