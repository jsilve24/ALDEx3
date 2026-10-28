# Adaptive covariance refinement. The forward checks and eight-step budget
# retain the forward-inference thresholds with a general bounded search.
# The context is private to one validated design; direct feature calls construct
# and validate their own context instead of trusting a public skip-check flag.
.blmm_context_key <- "ALDEx3 validated inference context v1"

blmm_inference_context <- function(X, contrasts, reTrms) {
  context <- new.env(parent = emptyenv())
  contrasts <- .aldex_validate_contrasts(contrasts, X)
  context$key <- .blmm_context_key
  context$cnms <- reTrms$cnms
  context$X <- X
  context$contrasts <- contrasts
  context$qr <- qr(X)
  context$random <- blmm_random_effect_metadata(reTrms)
  context$vectors <- if (is.null(contrasts)) diag(ncol(X)) else
    cbind(diag(ncol(X)), t(contrasts))
  context
}

blmm_valid_design_context <- function(context, X, contrasts, reTrms) {
  is.environment(context) && identical(context$key, .blmm_context_key) &&
    identical(context$X, X) && identical(context$contrasts, contrasts) &&
    identical(context$cnms, reTrms$cnms)
}

# Evaluate the shared covariance solve in natural Cholesky coordinates without
# passing theta through the log-coordinate conversion used by default BLMM.
blmm_adaptive_common_pieces <- function(theta, basis, lower) {
  if (length(theta) != length(lower) || any(!is.finite(theta)) || any(theta < lower))
    stop("blmm: proposed natural covariance parameters are non-finite or infeasible")
  T <- basis$Lambdat
  T@x <- theta[basis$Lind]
  A <- Matrix::forceSymmetric(T %*% basis$G %*% Matrix::t(T) +
                              Matrix::Diagonal(nrow(T)))
  if (basis$backend == "sparse") {
    factor <- Matrix::update(basis$symbolic, A)
    ldL2 <- as.numeric(Matrix::determinant(factor, logarithm = TRUE, sqrt = FALSE)$modulus)
  } else {
    factor <- chol(as.matrix(A))
    ldL2 <- 2 * sum(log(diag(factor)))
  }
  BX <- as.matrix(T %*% basis$ZX)
  solved_X <- blmm_factor_solve(factor, BX, basis$backend)
  M <- basis$XX - crossprod(BX, solved_X)
  LX <- t(chol((M + t(M)) / 2))
  dimnames(LX) <- NULL
  list(factor = factor, backend = basis$backend, BX = BX, LX = LX,
       T = T, design = basis, ldL2 = ldL2)
}

# Raw-response exact rescue is deliberately outside the approximation handler.
# An exact fitting error propagates without retrying accepted draws.
blmm_fit_feature_adaptive <- function(d, Y_d, X, basis, is_log, phi_init,
                                      lower, reTrms, formula, data, contrasts,
                                      strict, context) {
  N <- nrow(X)
  p <- ncol(X)
  S <- ncol(Y_d)
  vectors <- context$vectors
  fields <- c("estimate", "std.error", "df", "p.lower", "p.upper")
  random_rows <- length(context$random$names)
  fit <- stats::setNames(lapply(fields, function(field)
    array(NA_real_, c(p, 1L, S))), fields)
  fit$random.eff <- array(NA_real_, c(random_rows, 1L, S))
  if (!is.null(contrasts)) for (field in fields) {
    fit[[paste0("contrast.", field)]] <- array(NA_real_, c(nrow(contrasts), 1L, S),
      dimnames = list(rownames(contrasts), NULL, NULL))
  }
  diagnostics <- list(draw = seq_len(S), rescued = rep(FALSE, S),
    reason = rep("converged", S), status = rep("unresolved", S),
    anchor_status = rep("unresolved", S), iterations = integer(S), evaluations = integer(S),
    hessians = integer(S), checks = integer(S), objective_calls = integer(S),
    gradient_calls = integer(S), inference_calls = integer(S))
  for (field in c("decrement", "movement", "stability", "tail_change",
                 "statistic_change", "probe_objective_change"))
    diagnostics[[field]] <- rep(NA_real_, S)
  rescue_seconds <- 0
  initial <- tryCatch({
    if (context$qr$rank < p)
      stop("blmm: adaptive centering requires a full-rank fixed-effect design")
    offsets <- qr.coef(context$qr, Y_d)
    residual <- qr.resid(context$qr, Y_d)
    if (any(!is.finite(offsets)) || any(!is.finite(residual)))
      stop("blmm: adaptive centering produced non-finite values")
    response <- blmm_response_cache(basis, X, residual)
    natural_init <- reTrms$theta
    natural_log <- integer(length(natural_init))
    batch <- blmm_make_adfun(X, residual, basis, natural_log, natural_init, response)
    anchor <- blmm_bounded_optimize(batch, natural_init, lower, max_steps = 500L)
    if (anchor$status == "unresolved") stop("blmm adaptive anchor failed: ", anchor$reason)
    list(anchor = anchor, residual = residual, offsets = offsets, response = response)
  }, error = identity)
  if (!inherits(initial, "error"))
    diagnostics$anchor_status[] <- initial$anchor$status
  if (inherits(initial, "error")) {
    if (strict) stop(initial)
    diagnostics$rescued[] <- TRUE
    diagnostics$reason[] <- conditionMessage(initial)
  } else for (s in seq_len(S)) {
    iterations <- evaluations <- hessians <- checks <- objective_calls <-
      gradient_calls <- inference_calls <- 0L
    decrement <- movement <- stability <- tail_change <- statistic_change <-
      probe_objective_change <- NA_real_
    attempted <- tryCatch({
      response_s <- list(ZY = initial$response$ZY[, s, drop = FALSE],
        XY = initial$response$XY[, s, drop = FALSE], yy = initial$response$yy[s])
      residual_s <- initial$residual[, s, drop = FALSE]
      objective <- blmm_make_adfun(X, residual_s, basis, integer(length(phi_init)),
                                   reTrms$theta, response_s)
      last_value <- last_gradient <- last_raw_H <- last_information <- last_inference <- NULL
      value_at <- function(location) {
        if (!is.null(last_value) && identical(location, last_value$location))
          return(last_value$value)
        objective_calls <<- objective_calls + 1L
        value <- objective$fn(location)
        last_value <<- list(location = location, value = value)
        value
      }
      raw_hessian_at <- function(location) {
        if (!is.null(last_raw_H) && identical(location, last_raw_H$location))
          return(last_raw_H$H)
        hessians <<- hessians + 1L
        H <- objective$he(location)
        last_raw_H <<- list(location = location, H = H)
        H
      }
      information_at <- function(location) {
        if (!is.null(last_information) && identical(location, last_information$location))
          return(last_information$information)
        H <- blmm_validate_H(raw_hessian_at(location))
        information <- list(H = H, chol = chol(H))
        last_information <<- list(location = location, information = information)
        information
      }
      gradient_at <- function(location) {
        if (!is.null(last_gradient) && identical(location, last_gradient$location))
          return(last_gradient$gradient)
        gradient_calls <<- gradient_calls + 1L
        gradient <- as.vector(objective$gr(location))
        last_gradient <<- list(location = location, gradient = gradient)
        gradient
      }
      inference_at <- function(location, information) {
        if (!is.null(last_inference) && identical(location, last_inference$location) &&
            identical(information$H, last_inference$H)) return(last_inference$inference)
        inference_calls <<- inference_calls + 1L
        report <- objective$report(location)
        sensitivity <- blmm_scores_batch(objective, location, report$PWRSS_s,
          N, p, contrasts = vectors, inference = TRUE)
        dfs <- blmm_anchor_df(sensitivity$u[, 1, drop = FALSE], sensitivity$b,
          information$chol, N - p)
        pieces <- blmm_adaptive_common_pieces(location, basis, lower)
        fixed <- blmm_fixed_effects_draw(pieces, X, residual_s, response_s)
        estimates <- as.vector(crossprod(vectors, fixed$beta + initial$offsets[, s]))
        standard_errors <- sqrt(fixed$sigma2 * colSums(forwardsolve(pieces$LX, vectors)^2))
        values <- list(estimate = estimates, std.error = standard_errors, df = dfs,
          p.lower = stats::pt(estimates / standard_errors, dfs))
        if (any(!is.finite(unlist(values))) || any(standard_errors <= 0))
          stop("nonfinite inference")
        values$p.upper <- stats::pt(estimates / standard_errors, dfs, lower.tail = FALSE)
        random <- blmm_random_effect_vector(location,
          fixed$sigma2, reTrms, context$random)
        if (any(!is.finite(random))) stop("nonfinite random-effect inference")
        inference <- list(values = values, random = random)
        last_inference <<- list(location = location, H = information$H, inference = inference)
        inference
      }
      validate_interior <- function(location, H, value, gradient) {
        checks <<- checks + 1L
        information <- list(H = H, chol = chol(H))
        correction <- as.vector(backsolve(information$chol,
          forwardsolve(t(information$chol), gradient)))
        decrement <<- sqrt(sum(forwardsolve(t(information$chol), gradient)^2))
        movement <<- max(abs(correction))
        current <- inference_at(location, information)
        probe_location <- location - correction
        if (any(probe_location <= lower & is.finite(lower)))
          stop("infeasible forward inference probe")
        probe_objective_change <<- value_at(probe_location) - value
        if (!is.finite(probe_objective_change) || probe_objective_change >
            sqrt(.Machine$double.eps) * (1 + abs(value)))
          stop("probe objective validation failed")
        probe <- inference_at(probe_location, information_at(probe_location))
        stability <<- max(abs(log(probe$values$std.error / current$values$std.error)),
          abs(log(probe$values$df / current$values$df)))
        tail_change <<- max(abs(probe$values$p.lower - current$values$p.lower))
        statistic_change <<- max(abs(probe$values$estimate / probe$values$std.error -
          current$values$estimate / current$values$std.error))
        list(ok = stability <= 1e-3 && tail_change <= 1e-5 && statistic_change <= 1e-4,
          inference = current)
      }
      optimized <- blmm_bounded_optimize(
        list(fn = value_at, gr = gradient_at, he = raw_hessian_at),
        initial$anchor$par, lower, max_steps = 8L,
        validate_interior = validate_interior)
      iterations <- optimized$iterations
      evaluations <- optimized$diagnostics$evaluations
      diagnostics$status[s] <- optimized$status
      if (optimized$status == "boundary" && all(optimized$par == 0)) {
        blmm_validate_H(optimized$H)
        sigma2 <- sum(residual_s^2) / (N - p)
        R <- qr.R(context$qr)
        unpivoted <- matrix(0, p, ncol(vectors))
        unpivoted[seq_len(p), ] <- vectors[context$qr$pivot, , drop = FALSE]
        standard_errors <- sqrt(sigma2 * colSums(forwardsolve(t(R), unpivoted)^2))
        estimates <- as.vector(crossprod(vectors, initial$offsets[, s]))
        dfs <- rep(N - p, length(estimates))
        values <- list(estimate = estimates, std.error = standard_errors, df = dfs,
          p.lower = stats::pt(estimates / standard_errors, dfs),
          p.upper = stats::pt(estimates / standard_errors, dfs, lower.tail = FALSE))
        if (any(!is.finite(unlist(values))) || any(standard_errors <= 0))
          stop("nonfinite boundary inference")
        current <- list(values = values, random = blmm_random_effect_vector(
          optimized$par, sigma2, reTrms, context$random))
        diagnostics$status[s] <- "boundary_zero"
        diagnostics$reason[s] <- "stationary zero covariance"
      } else {
        if (optimized$status != "interior") stop(optimized$reason)
        current <- optimized$inference
      }
      current
    }, error = identity)
    for (field in c("iterations", "evaluations", "hessians", "checks", "objective_calls",
                    "gradient_calls", "inference_calls", "decrement", "movement", "stability",
                    "tail_change", "statistic_change", "probe_objective_change"))
      diagnostics[[field]][s] <- get(field)
    if (inherits(attempted, "error")) {
      if (strict) stop(attempted)
      diagnostics$rescued[s] <- TRUE
      diagnostics$reason[s] <- conditionMessage(attempted)
    } else {
      for (field in fields) {
        fit[[field]][, 1, s] <- attempted$values[[field]][seq_len(p)]
        if (!is.null(contrasts)) fit[[paste0("contrast.", field)]][, 1, s] <-
          attempted$values[[field]][p + seq_len(nrow(contrasts))]
      }
      fit$random.eff[, 1, s] <- attempted$random
    }
  }
  rejected <- which(diagnostics$rescued)
  if (length(rejected)) {
    started <- proc.time()[[3L]]
    exact <- tryCatch(blmm_exact_feature(Y_d[, rejected, drop = FALSE], formula,
      data, contrasts = contrasts), error = function(e) stop(
        "blmm: exact refit failed for feature ", d, " (draws ",
        paste(rejected, collapse = ", "), "): ", conditionMessage(e), call. = FALSE))
    rescue_seconds <- proc.time()[[3L]] - started
    for (field in names(exact)) fit[[field]][, , rejected] <- exact[[field]]
  }
  examples <- utils::head(rejected, 3L)
  message <- paste(sprintf("draw %d: %s", examples, diagnostics$reason[examples]), collapse = "; ")
  if (length(rejected) > length(examples)) message <- paste0(message, " ...")
  list(fit = fit, fallback = length(rejected) > 0L,
    fallback_message = if (inherits(initial, "error")) conditionMessage(initial) else message,
    fallback_draws = length(rejected), rejected_draws = rejected,
    full_fallback = inherits(initial, "error"), rescue_seconds = rescue_seconds,
    adaptive_diagnostics = if (isTRUE(getOption("ALDEx3.blmm.adaptive.diagnostics")) ||
      is.function(getOption("ALDEx3.blmm.adaptive.diagnostics")))
      list(diagnostics = as.data.frame(diagnostics), rescue_seconds = rescue_seconds) else NULL)
}

# Cache the random-effect blocks and lme4-compatible packaging order once.
blmm_random_effect_metadata <- function(reTrms) {
  blocks <- vector("list", length(reTrms$cnms))
  names_out <- character(0)
  theta_offset <- output_offset <- 0L
  for (g in seq_along(reTrms$cnms)) {
    terms <- reTrms$cnms[[g]]
    size <- length(terms)
    group <- names(reTrms$cnms)[g] %||% paste0("grp", g)
    template <- matrix(seq_len(size * size), size, size)
    lower <- which(lower.tri(template, diag = TRUE))
    covariance <- c(diag(template), template[lower.tri(template)])
    block_names <- paste(group, terms, sep = ".")
    if (size > 1L) for (j in seq_len(size - 1L)) for (i in seq.int(j + 1L, size))
      block_names <- c(block_names, paste(group, terms[j], terms[i], sep = "."))
    blocks[[g]] <- list(size = size, lower = lower,
      theta = theta_offset + seq_along(lower), covariance = covariance,
      output = output_offset + seq_along(covariance))
    names_out <- c(names_out, block_names)
    theta_offset <- theta_offset + length(lower)
    output_offset <- output_offset + length(covariance)
  }
  list(blocks = blocks, names = c(names_out, "Residual"))
}
