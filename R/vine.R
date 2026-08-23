# Vine copulas are the multivariate nonexchangeable path. Clayton/Gumbel/Frank/Joe
# from the copula package remain exchangeable for d > 2.

autocopula_vine_aliases <- function() {
  c(rvine = "RVine", vine = "vine", cvine = "CVine", dvine = "DVine")
}

autocopula_is_vine_family_name <- function(name) {
  key <- tolower(gsub("[^A-Za-z]", "", as.character(name)[1L]))
  key %in% names(autocopula_vine_aliases())
}

autocopula_is_vine_model <- function(model) {
  inherits(model, "autocopula_vine") || inherits(model, "RVineMatrix")
}

autocopula_as_rvine <- function(model) {
  if (inherits(model, "RVineMatrix")) return(model)
  if (inherits(model, "autocopula_vine") && inherits(model$rvine, "RVineMatrix")) {
    return(model$rvine)
  }
  stop("Fitted vine model does not contain an RVineMatrix.", call. = FALSE)
}

autocopula_require_vinecopula <- function(action = "vine copula operations") {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop("VineCopula is required for ", action, ".", call. = FALSE)
  }
  invisible(TRUE)
}

autocopula_default_vine_family_set <- function() {
  # Gaussian, t, Clayton, Gumbel, Frank, Joe, BB1, BB7 and survival (180) rotations.
  c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 9L, 13L, 14L, 16L, 17L, 19L)
}

autocopula_pair_family_name_codes <- function() {
  c(
    gaussian = 1L, normal = 1L, n = 1L,
    t = 2L, tcopula = 2L, student = 2L, studentt = 2L,
    clayton = 3L, c = 3L,
    gumbel = 4L, g = 4L,
    frank = 5L, f = 5L,
    joe = 6L, j = 6L,
    bb1 = 7L, bb6 = 8L, bb7 = 9L, bb8 = 10L,
    `survival clayton` = 13L, `rotated clayton (180)` = 13L, sc = 13L,
    `survival gumbel` = 14L, `rotated gumbel (180)` = 14L, sg = 14L,
    `survival joe` = 16L, `rotated joe (180)` = 16L, sj = 16L,
    `survival bb1` = 17L, `rotated bb1 (180)` = 17L, sbb1 = 17L,
    `survival bb6` = 18L, sbb6 = 18L,
    `survival bb7` = 19L, `rotated bb7 (180)` = 19L, sbb7 = 19L,
    `survival bb8` = 20L, sbb8 = 20L
  )
}

autocopula_resolve_vine_type <- function(vine_type) {
  key <- tolower(gsub("[^A-Za-z]", "", as.character(vine_type)[1L]))
  if (!nzchar(key)) key <- "rvine"
  if (!key %in% c("rvine", "cvine", "dvine", "vine")) {
    stop("vine_type must be 'rvine', 'cvine', or 'dvine'.", call. = FALSE)
  }
  if (identical(key, "vine")) "rvine" else key
}

autocopula_vine_type_for_family <- function(family, vine_type = "rvine") {
  key <- tolower(gsub("[^A-Za-z]", "", as.character(family)[1L]))
  if (identical(key, "cvine")) return("cvine")
  if (identical(key, "dvine")) return("dvine")
  if (identical(key, "rvine")) return("rvine")
  if (identical(key, "vine")) return(autocopula_resolve_vine_type(vine_type))
  stop("Not a vine family name: ", family, call. = FALSE)
}

autocopula_vine_structure_code <- function(vine_type) {
  switch(autocopula_resolve_vine_type(vine_type),
    rvine = 0L,
    cvine = 1L,
    dvine = NA_integer_,
    stop("Unsupported vine_type: ", vine_type, call. = FALSE))
}

autocopula_resolve_family_set <- function(family_set) {
  if (is.null(family_set) || (length(family_set) == 1L && is.na(family_set)[1L])) {
    return(autocopula_default_vine_family_set())
  }
  if (is.numeric(family_set)) {
    codes <- as.integer(family_set)
    if (any(!is.finite(codes))) {
      stop("family_set numeric codes must be finite integers.", call. = FALSE)
    }
    return(unique(codes))
  }
  labels <- trimws(as.character(family_set))
  labels <- labels[nzchar(labels)]
  if (!length(labels)) {
    return(autocopula_default_vine_family_set())
  }
  aliases <- autocopula_pair_family_name_codes()
  vapply(labels, function(label) {
    key <- tolower(gsub("[_-]", " ", label))
    if (key %in% names(aliases)) return(unname(aliases[[key]]))
    compact <- tolower(gsub("[^A-Za-z0-9]", "", label))
    compact_names <- tolower(gsub("[^A-Za-z0-9]", "", names(aliases)))
    hit <- match(compact, compact_names)
    if (!is.na(hit)) return(unname(aliases[[hit]]))
    if (requireNamespace("VineCopula", quietly = TRUE)) {
      converted <- try(VineCopula::BiCopName(label), silent = TRUE)
      if (!inherits(converted, "try-error")) {
        code <- suppressWarnings(as.integer(unname(converted)[1L]))
        if (length(code) && is.finite(code)) return(code)
      }
    }
    stop("Unknown vine pair-copula family in family_set: ", label, call. = FALSE)
  }, integer(1L), USE.NAMES = FALSE)
}

autocopula_pair_family_name <- function(code) {
  code <- as.integer(code)[1L]
  if (!is.finite(code)) return(NA_character_)
  if (requireNamespace("VineCopula", quietly = TRUE)) {
    nm <- try(VineCopula::BiCopName(code, short = FALSE), silent = TRUE)
    if (!inherits(nm, "try-error") && length(nm)) return(as.character(nm)[1L])
  }
  as.character(code)
}

autocopula_vine_description <- function(vine_type) {
  switch(autocopula_resolve_vine_type(vine_type),
    rvine = paste(
      "Regular vine copula (nonexchangeable pair-copula construction).",
      "Use for pair-specific tails in dimension d > 2."),
    cvine = paste(
      "C-vine copula (star trees; nonexchangeable pair-copula construction).",
      "Use for pair-specific tails in dimension d > 2."),
    dvine = paste(
      "D-vine copula (path trees; nonexchangeable pair-copula construction).",
      "Use for pair-specific tails in dimension d > 2."),
    "Vine copula")
}

autocopula_dvine_order <- function(u, treecrit = "tau") {
  d <- ncol(u)
  if (d < 2L) stop("D-vine fitting requires at least two variables.", call. = FALSE)
  if (d == 2L) return(c(1L, 2L))
  method <- if (identical(as.character(treecrit)[1L], "rho")) "spearman" else
    "kendall"
  tau <- abs(stats::cor(u, method = method, use = "pairwise.complete.obs"))
  tau[is.na(tau)] <- 0
  diag(tau) <- 0
  best <- NULL
  best_score <- -Inf
  for (start in seq_len(d)) {
    used <- logical(d)
    path <- integer(d)
    path[1L] <- start
    used[start] <- TRUE
    score <- 0
    for (i in 2:d) {
      prev <- path[i - 1L]
      cand <- which(!used)
      nxt <- cand[which.max(tau[prev, cand])]
      score <- score + tau[prev, nxt]
      path[i] <- nxt
      used[nxt] <- TRUE
    }
    if (score > best_score) {
      best_score <- score
      best <- path
    }
  }
  as.integer(best)
}

autocopula_trunclevel <- function(trunc_lvl) {
  if (is.null(trunc_lvl) || (length(trunc_lvl) == 1L && is.na(trunc_lvl))) {
    return(NA_integer_)
  }
  lvl <- as.integer(trunc_lvl)[1L]
  if (!is.finite(lvl) || lvl < 1L) {
    stop("trunc_lvl must be a positive integer or NA.", call. = FALSE)
  }
  lvl
}

autocopula_pack_vine <- function(rvm, vine_type, family_set, trunc_lvl,
    treecrit, var_names, pobs = NULL) {
  if (!is.null(var_names) && (is.null(rvm$names) || !length(rvm$names))) {
    rvm$names <- as.character(var_names)
  }
  d <- nrow(rvm$Matrix)
  trees_possible <- max(d - 1L, 0L)
  trunc_used <- autocopula_trunclevel(trunc_lvl)
  if (is.na(trunc_used)) trunc_used <- trees_possible
  structure(list(
    rvine = rvm,
    vine_type = autocopula_resolve_vine_type(vine_type),
    vine_type_fitted = if (!is.null(rvm$type)) rvm$type else vine_type,
    family_set = as.integer(family_set),
    trunc_lvl = autocopula_trunclevel(trunc_lvl),
    trunc_lvl_used = min(as.integer(trunc_used), trees_possible),
    treecrit = as.character(treecrit)[1L],
    var_names = as.character(var_names),
    nobs = if (!is.null(rvm$nobs)) as.integer(rvm$nobs) else
      if (!is.null(pobs)) nrow(pobs) else NA_integer_,
    pobs = pobs
  ), class = c("autocopula_vine", "list"))
}

autocopula_fit_vine <- function(u, vine_type = "rvine", family_set = NULL,
    trunc_lvl = NA_integer_, treecrit = "tau", var_names = NULL) {
  autocopula_require_vinecopula("vine copula fitting")
  u <- as.matrix(u)
  storage.mode(u) <- "double"
  if (ncol(u) < 2L) {
    stop("Vine copula fitting requires at least two variables.", call. = FALSE)
  }
  if (is.null(var_names)) {
    var_names <- colnames(u)
    if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(u)))
  }
  colnames(u) <- var_names
  family_set <- autocopula_resolve_family_set(family_set)
  trunclevel <- autocopula_trunclevel(trunc_lvl)
  vine_type <- autocopula_resolve_vine_type(vine_type)
  treecrit <- as.character(treecrit)[1L]
  if (!nzchar(treecrit)) treecrit <- "tau"

  rvm <- if (identical(vine_type, "dvine")) {
    d <- ncol(u)
    order <- autocopula_dvine_order(u, treecrit = treecrit)
    n_pairs <- as.integer(d * (d - 1L) / 2L)
    skeleton <- VineCopula::D2RVine(order, family = rep(0, n_pairs),
      par = rep(0, n_pairs))
    VineCopula::RVineCopSelect(u, familyset = family_set,
      Matrix = skeleton$Matrix, trunclevel = trunclevel, rotations = TRUE)
  } else {
    VineCopula::RVineStructureSelect(u, familyset = family_set,
      type = autocopula_vine_structure_code(vine_type),
      trunclevel = trunclevel, treecrit = treecrit, rotations = TRUE,
      progress = FALSE)
  }
  if (!inherits(rvm, "RVineMatrix")) {
    stop("Vine fitting did not return an RVineMatrix.", call. = FALSE)
  }
  rvm$names <- var_names
  autocopula_pack_vine(rvm, vine_type, family_set, trunc_lvl, treecrit,
    var_names, pobs = u)
}

autocopula_make_vine_fit_function <- function(vine_type, family_set, trunc_lvl,
    treecrit, var_names) {
  force(vine_type)
  force(family_set)
  force(trunc_lvl)
  force(treecrit)
  force(var_names)
  function(data) {
    autocopula_fit_vine(data, vine_type = vine_type, family_set = family_set,
      trunc_lvl = trunc_lvl, treecrit = treecrit, var_names = var_names)
  }
}

autocopula_register_vine_families <- function(fitter, families, vine_type,
    family_set, trunc_lvl, treecrit) {
  vine_names <- families[vapply(families, autocopula_is_vine_family_name,
    logical(1L))]
  if (!length(vine_names)) return(invisible(fitter))
  var_names <- names(fitter$data)
  for (nm in unique(vine_names)) {
    this_type <- autocopula_vine_type_for_family(nm, vine_type)
    fitter$copula_library[[nm]] <- list(
      description = autocopula_vine_description(this_type),
      fit_function = autocopula_make_vine_fit_function(this_type, family_set,
        trunc_lvl, treecrit, var_names)
    )
  }
  invisible(fitter)
}

autocopula_back_transform_u <- function(u, marginals) {
  u <- data.table::as.data.table(u)
  for (col_name in names(marginals)) {
    if (col_name %in% names(u)) {
      u[, (col_name) := marginals[[col_name]](autocopula_clamp_u(get(col_name)))]
    }
  }
  u
}

autocopula_x_to_u <- function(x, value) {
  x <- as.numeric(x)
  n <- length(x)
  u <- (sum(x < value, na.rm = TRUE) + 0.5 * sum(x == value, na.rm = TRUE)) /
    (n + 1)
  autocopula_clamp_u(u)
}

autocopula_vine_sim_u <- function(model, n) {
  autocopula_require_vinecopula("vine copula simulation")
  rvm <- autocopula_as_rvine(model)
  n <- as.integer(n)[1L]
  if (!is.finite(n) || n < 1L) stop("n must be a positive integer.", call. = FALSE)
  out <- VineCopula::RVineSim(n, rvm)
  colnames(out) <- if (!is.null(rvm$names)) rvm$names else
    if (inherits(model, "autocopula_vine")) model$var_names else
      paste0("V", seq_len(ncol(out)))
  out
}

autocopula_rosenblatt_first_variable <- function(rvm) {
  as.integer(diag(rvm$Matrix)[nrow(rvm$Matrix)])
}

autocopula_vine_conditional_bicop_u <- function(rvm, known_indices, known_u, n) {
  fam <- as.integer(rvm$family[2L, 1L])
  par <- as.numeric(rvm$par[2L, 1L])
  par2 <- as.numeric(rvm$par2[2L, 1L])
  if (!is.finite(par2)) par2 <- 0
  sampled <- autocopula_conditional_bicop_u(known_u, known_indices, n, fam, par,
    par2)
  result <- matrix(NA_real_, nrow = n, ncol = 2L)
  result[, known_indices] <- known_u
  result[, setdiff(c(1L, 2L), known_indices)] <- sampled
  result
}

autocopula_try_rvine_cond_sim <- function(rvm, known_indices, known_u, n) {
  if (!exists("RVineCondSim", envir = asNamespace("VineCopula"),
    inherits = FALSE)) {
    return(NULL)
  }
  fn <- get("RVineCondSim", envir = asNamespace("VineCopula"), inherits = FALSE)
  fml <- names(formals(fn))
  d <- nrow(rvm$Matrix)
  remaining <- setdiff(seq_len(d), known_indices)
  n <- as.integer(n)[1L]
  try_call <- function(...) {
    out <- try(fn(...), silent = TRUE)
    if (inherits(out, "try-error")) return(NULL)
    out <- as.matrix(out)
    if (nrow(out) != n) return(NULL)
    if (ncol(out) == d) return(out)
    if (ncol(out) == length(remaining)) {
      result <- matrix(NA_real_, nrow = n, ncol = d)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining] <- out
      return(result)
    }
    NULL
  }
  cond_u <- matrix(rep(autocopula_clamp_u(known_u), each = n), nrow = n)
  args <- list()
  if ("N" %in% fml) args$N <- n
  if ("n" %in% fml && !"N" %in% fml) args$n <- n
  if ("RVM" %in% fml) args$RVM <- rvm
  if ("rvm" %in% fml) args$rvm <- rvm
  if ("cond.vars" %in% fml) args[["cond.vars"]] <- as.integer(known_indices)
  if ("condVars" %in% fml) args$condVars <- as.integer(known_indices)
  if ("cond.idx" %in% fml) args[["cond.idx"]] <- as.integer(known_indices)
  if ("cond.u" %in% fml) args[["cond.u"]] <- cond_u
  if ("U" %in% fml) {
    simdata <- matrix(stats::runif(n * length(remaining)), nrow = n)
    args$U <- simdata
  }
  if ("simdata" %in% fml) {
    args$simdata <- matrix(stats::runif(n * length(remaining)), nrow = n)
  }
  do.call(try_call, args)
}

autocopula_vine_conditional_prefix_u <- function(rvm, known_indices, known_u, n) {
  first <- autocopula_rosenblatt_first_variable(rvm)
  if (length(known_indices) != 1L || known_indices[1L] != first) return(NULL)
  d <- nrow(rvm$Matrix)
  u <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  colnames(u) <- if (!is.null(rvm$names)) rvm$names else paste0("V", seq_len(d))
  u[, first] <- autocopula_clamp_u(known_u[1L])
  sampled <- VineCopula::RVineSim(n, rvm, U = u)
  # Inverse Rosenblatt only preserves a coordinate that is first in sampling
  # order. If VineCopula transformed it, this is not an analytic conditional.
  if (max(abs(sampled[, first] - u[, first])) > 1e-6) return(NULL)
  sampled
}

autocopula_clamp_u_matrix <- function(u) {
  u <- as.matrix(u)
  storage.mode(u) <- "double"
  u[] <- autocopula_clamp_u(u)
  u
}

autocopula_vine_conditional_mc_u <- function(rvm, known_indices, known_u, n) {
  n <- as.integer(n)[1L]
  d <- nrow(rvm$Matrix)
  remaining <- setdiff(seq_len(d), known_indices)
  n_pool <- max(n * 80L, 4000L)
  pool <- VineCopula::RVineSim(n_pool, rvm)
  z_known <- stats::qnorm(autocopula_clamp_u(known_u))
  z_pool <- stats::qnorm(autocopula_clamp_u_matrix(
    pool[, known_indices, drop = FALSE]))
  sds <- apply(z_pool, 2L, stats::sd)
  bandwidth <- max(0.08, 1.06 * mean(sds, na.rm = TRUE) *
    n_pool^(-1 / (4 + length(known_indices))))
  delta <- sweep(z_pool, 2L, z_known, "-")
  log_w <- -0.5 * rowSums((delta / bandwidth)^2)
  log_w <- log_w - max(log_w)
  w <- exp(log_w)
  w <- w / sum(w)
  ess <- 1 / sum(w^2)
  idx <- sample.int(n_pool, n, replace = TRUE, prob = w)
  sampled <- pool[idx, , drop = FALSE]
  sampled[, known_indices] <- matrix(rep(autocopula_clamp_u(known_u), each = n),
    nrow = n)
  sampled[, remaining] <- autocopula_clamp_u_matrix(
    sampled[, remaining, drop = FALSE])
  list(u = sampled, provenance = list(
    method = "monte_carlo_kernel_sir",
    n_pool = n_pool,
    bandwidth = bandwidth,
    ess = ess,
    note = paste("Analytic vine h-functions were not available for this",
      "conditioning set; samples are kernel-weighted draws from the fitted vine.")
  ))
}

autocopula_vine_conditional_u <- function(model, known_indices, known_u, n) {
  autocopula_require_vinecopula("vine copula conditional simulation")
  rvm <- autocopula_as_rvine(model)
  d <- nrow(rvm$Matrix)
  known_indices <- as.integer(known_indices)
  remaining <- setdiff(seq_len(d), known_indices)
  if (!length(known_indices) || !length(remaining)) {
    stop("Conditioning requires at least one known and one unknown variable.",
      call. = FALSE)
  }
  if (length(known_u) != length(known_indices)) {
    stop("known_u does not match known_indices.", call. = FALSE)
  }
  n <- as.integer(n)[1L]
  known_u <- autocopula_clamp_u(known_u)
  provenance <- list(method = "analytic", n_pool = NA_integer_,
    bandwidth = NA_real_, ess = NA_real_, note = "")

  sampled <- NULL
  if (d == 2L) {
    sampled <- autocopula_vine_conditional_bicop_u(rvm, known_indices, known_u, n)
    provenance$method <- "hfunction_bicop"
    provenance$note <- "VineCopula inverse h-functions on the fitted pair-copula."
  }
  if (is.null(sampled)) {
    cond_sim <- autocopula_try_rvine_cond_sim(rvm, known_indices, known_u, n)
    if (!is.null(cond_sim)) {
      sampled <- cond_sim
      provenance$method <- "RVineCondSim"
      provenance$note <- "VineCopula::RVineCondSim on the fitted RVineMatrix."
    }
  }
  if (is.null(sampled)) {
    prefix <- autocopula_vine_conditional_prefix_u(rvm, known_indices, known_u, n)
    if (!is.null(prefix)) {
      sampled <- prefix
      provenance$method <- "inverse_rosenblatt_prefix"
      provenance$note <- paste("Inverse Rosenblatt with the conditioned variable",
        "first in the fitted vine sampling order.")
    }
  }
  if (is.null(sampled)) {
    mc <- autocopula_vine_conditional_mc_u(rvm, known_indices, known_u, n)
    sampled <- mc$u
    provenance <- mc$provenance
  }
  colnames(sampled) <- if (!is.null(rvm$names)) rvm$names else
    if (inherits(model, "autocopula_vine")) model$var_names else
      paste0("V", seq_len(d))
  attr(sampled, "conditional_provenance") <- provenance
  sampled
}

autocopula_vine_pair_table <- function(rvm) {
  tab <- NULL
  printed <- try(utils::capture.output(tab <- summary(rvm)), silent = TRUE)
  if (!inherits(printed, "try-error") && is.data.frame(tab) && nrow(tab)) {
    out <- data.table::as.data.table(tab)
    if ("utd" %in% names(out)) data.table::setnames(out, "utd", "lambda_upper")
    if ("ltd" %in% names(out)) data.table::setnames(out, "ltd", "lambda_lower")
    if ("cop" %in% names(out) && !"family_name" %in% names(out)) {
      data.table::setnames(out, "cop", "family_name")
    }
    if ("family" %in% names(out) && !"family_code" %in% names(out)) {
      data.table::setnames(out, "family", "family_code")
    }
    return(out)
  }
  d <- nrow(rvm$Matrix)
  rows <- list()
  for (tree in seq_len(d - 1L)) {
    i <- d + 1L - tree
    for (edge in seq_len(d - tree)) {
      fam <- as.integer(rvm$family[i, edge])
      rows[[length(rows) + 1L]] <- data.table::data.table(
        tree = tree,
        family_code = fam,
        family_name = autocopula_pair_family_name(fam),
        par = as.numeric(rvm$par[i, edge]),
        par2 = as.numeric(rvm$par2[i, edge]),
        tau = as.numeric(rvm$tau[i, edge]),
        lambda_upper = as.numeric(rvm$taildep$upper[i, edge]),
        lambda_lower = as.numeric(rvm$taildep$lower[i, edge])
      )
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

autocopula_vine_tail_notes <- function(pairs) {
  if (!nrow(pairs)) return("No pair-copulas recorded.")
  upper <- if ("lambda_upper" %in% names(pairs)) pairs$lambda_upper else
    if ("utd" %in% names(pairs)) pairs$utd else rep(NA_real_, nrow(pairs))
  lower <- if ("lambda_lower" %in% names(pairs)) pairs$lambda_lower else
    if ("ltd" %in% names(pairs)) pairs$ltd else rep(NA_real_, nrow(pairs))
  names <- if ("family_name" %in% names(pairs)) pairs$family_name else
    as.character(pairs$family_code)
  upper_pairs <- names[is.finite(upper) & upper > 1e-8]
  lower_pairs <- names[is.finite(lower) & lower > 1e-8]
  paste0(
    "Pair-specific tails are allowed; this is the nonexchangeable path for d > 2. ",
    "Upper tail dependence in: ",
    if (length(upper_pairs)) paste(unique(upper_pairs), collapse = ", ") else "none",
    ". Lower tail dependence in: ",
    if (length(lower_pairs)) paste(unique(lower_pairs), collapse = ", ") else "none",
    ".")
}

autocopula_vine_aic_bic <- function(model) {
  rvm <- autocopula_as_rvine(model)
  aic <- if (!is.null(rvm$AIC)) rvm$AIC else NA_real_
  bic <- if (!is.null(rvm$BIC)) rvm$BIC else NA_real_
  loglik <- if (!is.null(rvm$logLik)) rvm$logLik else NA_real_
  pobs <- if (inherits(model, "autocopula_vine")) model$pobs else NULL
  missing_ic <- !is.finite(aic) || !is.finite(bic)
  if (missing_ic && !is.null(pobs) && requireNamespace("VineCopula", quietly = TRUE)) {
    if (!is.finite(aic)) {
      got <- try(VineCopula::RVineAIC(pobs, rvm), silent = TRUE)
      if (!inherits(got, "try-error")) aic <- unname(got$AIC)
    }
    if (!is.finite(bic)) {
      got <- try(VineCopula::RVineBIC(pobs, rvm), silent = TRUE)
      if (!inherits(got, "try-error")) bic <- unname(got$BIC)
    }
    if (!is.finite(loglik)) {
      got <- try(VineCopula::RVineLogLik(pobs, rvm), silent = TRUE)
      if (!inherits(got, "try-error")) {
        loglik <- if (is.list(got) && !is.null(got$loglik)) got$loglik else
          unname(as.numeric(got)[1L])
      }
    }
  }
  list(logLik = loglik, AIC = aic, BIC = bic)
}

autocopula_vine_diagnose <- function(fit, family) {
  model <- fit$models[[family]]
  rvm <- autocopula_as_rvine(model)
  pairs <- autocopula_vine_pair_table(rvm)
  ic <- autocopula_vine_aic_bic(model)
  trunc_lvl <- if (inherits(model, "autocopula_vine")) model$trunc_lvl else
    NA_integer_
  trunc_used <- if (inherits(model, "autocopula_vine")) model$trunc_lvl_used else
    nrow(rvm$Matrix) - 1L
  list(
    family = family,
    vine_type = if (inherits(model, "autocopula_vine")) model$vine_type else
      rvm$type,
    vine_type_fitted = rvm$type,
    pair_copulas = pairs,
    logLik = ic$logLik,
    AIC = ic$AIC,
    BIC = ic$BIC,
    trunc_lvl = trunc_lvl,
    trunc_lvl_used = trunc_used,
    treecrit = if (inherits(model, "autocopula_vine")) model$treecrit else NA_character_,
    family_set = if (inherits(model, "autocopula_vine")) model$family_set else
      NA_integer_,
    tail_dependence = autocopula_vine_tail_notes(pairs),
    exchangeability = paste(
      "Vines are nonexchangeable. d>2 Clayton/Gumbel/Frank/Joe fits from the",
      "copula package remain exchangeable; use RVine/CVine/DVine for pair-specific tails.")
  )
}

autocopula_builtin_family_rows <- function() {
  data.table::data.table(
    family = c(
      "Gaussian", "tCopula", "Clayton", "Gumbel", "Frank", "Joe",
      "Galambos", "HuslerReiss", "tEV", "Plackett", "FGM",
      "BB1", "BB6", "BB7", "BB8",
      "Rotated Clayton (180)", "Rotated Gumbel (180)", "Rotated Joe (180)",
      "Rotated BB1 (180)", "Rotated BB6 (180)", "Rotated BB7 (180)",
      "Rotated BB8 (180)",
      "Tawn Type 1", "Rotated Tawn Type 1 (180)",
      "Tawn Type 2", "Rotated Tawn Type 2 (180)",
      "RVine", "vine", "CVine", "DVine"
    ),
    description = c(
      "Gaussian copula with a correlation matrix.",
      "t-copula with degrees of freedom and a correlation matrix.",
      "Clayton copula for lower tail dependence. Exchangeable for d > 2.",
      "Gumbel copula for upper tail dependence. Exchangeable for d > 2.",
      "Frank copula with no tail dependence. Exchangeable for d > 2.",
      "Joe copula with upper tail dependence. Exchangeable for d > 2.",
      "Galambos copula for modeling extreme values.",
      "Husler-Reiss copula for extreme value modeling.",
      "t-EV copula for extreme value dependence.",
      "Plackett copula for symmetric dependence.",
      "FGM (Farlie-Gumbel-Morgenstern) copula for weak dependence.",
      "BB1 copula for flexible upper and lower tail dependence.",
      "BB6 copula for flexible tail dependence.",
      "BB7 copula for asymmetrical dependence.",
      "BB8 copula for modeling both tail dependencies.",
      "180-degree rotated Clayton copula (survival Clayton).",
      "180-degree rotated Gumbel copula (survival Gumbel).",
      "180-degree rotated Joe copula (survival Joe).",
      "180-degree rotated BB1 copula.",
      "180-degree rotated BB6 copula.",
      "180-degree rotated BB7 copula.",
      "180-degree rotated BB8 copula.",
      "Tawn Type 1 copula for asymmetrical lower tail dependence.",
      "180-degree rotated Tawn Type 1 copula.",
      "Tawn Type 2 copula for asymmetrical upper tail dependence.",
      "180-degree rotated Tawn Type 2 copula.",
      "Regular vine: nonexchangeable pair-copula construction for d >= 2.",
      "Alias for a vine copula; vine_type chooses R-/C-/D-vine.",
      "C-vine: star trees, nonexchangeable pair-copula construction.",
      "D-vine: path trees, nonexchangeable pair-copula construction."
    ),
    class = c(
      "elliptical", "elliptical", "archimedean", "archimedean", "archimedean",
      "archimedean", "extreme_value", "extreme_value", "extreme_value",
      "other", "other", "bb", "bb", "bb", "bb",
      "rotated", "rotated", "rotated", "rotated", "rotated", "rotated", "rotated",
      "tawn", "tawn", "tawn", "tawn",
      "vine", "vine", "vine", "vine"
    ),
    dimension = c(
      "d>=2", "d>=2", "d>=2", "d>=2", "d>=2", "d>=2",
      "d=2", "d=2", "d=2", "d=2", "d=2",
      "d=2", "d=2", "d=2", "d=2",
      "d=2", "d=2", "d=2", "d=2", "d=2", "d=2", "d=2",
      "d=2", "d=2", "d=2", "d=2",
      "d>=2", "d>=2", "d>=2", "d>=2"
    ),
    exchangeable = c(
      FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
      NA, NA, NA, NA, NA,
      NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA,
      FALSE, FALSE, FALSE, FALSE
    ),
    notes = c(
      "Correlation matrix is nonexchangeable when pair correlations differ.",
      "Correlation matrix plus df; nonexchangeable when pair correlations differ.",
      "d>2 Clayton from copula is exchangeable; use RVine for pair-specific tails.",
      "d>2 Gumbel from copula is exchangeable; use RVine for pair-specific tails.",
      "d>2 Frank from copula is exchangeable; use RVine for pair-specific tails.",
      "d>2 Joe from copula is exchangeable; use RVine for pair-specific tails.",
      "Bivariate extreme-value copula.",
      "Bivariate extreme-value copula.",
      "Bivariate extreme-value copula.",
      "Bivariate.",
      "Bivariate; weak dependence only.",
      "Bivariate two-parameter family.",
      "Bivariate two-parameter family.",
      "Bivariate two-parameter family.",
      "Bivariate two-parameter family.",
      "Bivariate survival rotation.",
      "Bivariate survival rotation.",
      "Bivariate survival rotation.",
      "Bivariate survival rotation.",
      "Bivariate survival rotation.",
      "Bivariate survival rotation.",
      "Bivariate survival rotation.",
      "Bivariate asymmetric family.",
      "Bivariate asymmetric family.",
      "Bivariate asymmetric family.",
      "Bivariate asymmetric family.",
      "Preferred nonexchangeable multivariate model; pair-family set via family_set.",
      "Generic vine alias resolved through vine_type.",
      "Star vine; still pair-specific and nonexchangeable.",
      "Path vine; still pair-specific and nonexchangeable."
    )
  )
}
