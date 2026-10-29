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

# Raw-response exact rescue is deliberately outside the approximation handler.
# An exact fitting error propagates without retrying accepted draws.
blmm_fit_feature_adaptive <- function(d, Y_d, X, basis, lower, reTrms, formula,
                                      data, contrasts, strict, context) {
  N <- nrow(X)
  p <- ncol(X)
  S <- ncol(Y_d)
  q <- length(reTrms$theta)
  vectors <- context$vectors
  fit <- blmm_output_arrays(p, length(context$random$names), 1L, S, contrasts)
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
    batch <- blmm_make_adfun(X, residual, basis, integer(q), reTrms$theta, response)
    anchor <- blmm_bounded_optimize(batch, reTrms$theta, lower, max_steps = 500L)
    # The shared fit is only a starting point: every draw is validated on its
    # own below, so an unresolved but finite anchor is still a usable start.
    if (isTRUE(anchor$failed) || length(anchor$par) != q || any(!is.finite(anchor$par)))
      stop("blmm adaptive anchor failed: ", anchor$reason)
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
    stale_report_error <- NULL
    attempted <- tryCatch({
      response_s <- list(ZY = initial$response$ZY[, s, drop = FALSE],
        XY = initial$response$XY[, s, drop = FALSE], yy = initial$response$yy[s])
      residual_s <- initial$residual[, s, drop = FALSE]
      objective <- blmm_make_adfun(X, residual_s, basis, integer(q),
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
      partial_information_at <- function(location) {
        blmm_partial_information(location, lower, context$random,
          raw_hessian_at(location), gradient_at(location))
      }
      inference_at <- function(location, information) {
        indices <- if (is.null(information$indices)) seq_along(location) else information$indices
        if (!is.null(last_inference) && identical(location, last_inference$location) &&
            identical(information$H, last_inference$H) &&
            identical(indices, last_inference$indices)) return(last_inference$inference)
        inference_calls <<- inference_calls + 1L
        report <- objective$report(location)
        if (!is.matrix(report$LX) || !identical(dim(report$LX), c(p, p)) ||
            any(!is.finite(report$LX)) || length(report$PWRSS_s) != 1L ||
            !is.finite(report$PWRSS_s) || report$PWRSS_s <= 0)
          stop("blmm: invalid native inference report")
        if (is.null(report$beta)) {
          error <- simpleError("blmm: native report lacks beta; reinstall the package to update its compiled kernel")
          class(error) <- c("blmm_missing_report_beta", class(error))
          stale_report_error <<- error
          stop(error)
        }
        if (!is.matrix(report$beta) || !identical(dim(report$beta), c(p, 1L)) ||
            any(!is.finite(report$beta)))
          stop("blmm: invalid native inference report")
        sensitivity <- blmm_scores_batch(objective, location, report$PWRSS_s,
          N, p, contrasts = vectors, inference = TRUE)
        dropped <- setdiff(seq_along(location), indices)
        if (length(dropped) && (any(sensitivity$u[dropped, , drop = FALSE] != 0) ||
                                any(sensitivity$b[dropped, , drop = FALSE] != 0)))
          stop("partial boundary sensitivity does not decouple")
        dfs <- blmm_anchor_df(sensitivity$u[indices, 1, drop = FALSE],
          sensitivity$b[indices, , drop = FALSE],
          information$chol, N - p)
        beta <- as.vector(report$beta[, 1])
        sigma2 <- as.numeric(report$PWRSS_s) / (N - p)
        estimates <- as.vector(crossprod(vectors, beta + initial$offsets[, s]))
        standard_errors <- sqrt(sigma2 * colSums(forwardsolve(report$LX, vectors)^2))
        values <- list(estimate = estimates, std.error = standard_errors, df = dfs,
          p.lower = stats::pt(estimates / standard_errors, dfs))
        if (any(!is.finite(unlist(values))) || any(standard_errors <= 0))
          stop("nonfinite inference")
        values$p.upper <- stats::pt(estimates / standard_errors, dfs, lower.tail = FALSE)
        random <- blmm_random_effect_vector(location,
          sigma2, reTrms, context$random)
        if (any(!is.finite(random))) stop("nonfinite random-effect inference")
        inference <- list(values = values, random = random)
        last_inference <<- list(location = location, H = information$H,
          indices = indices, inference = inference)
        inference
      }
      # Forward check shared by interior and partial-boundary fits: one Newton
      # probe in the retained coordinates must neither worsen the objective nor
      # materially move SEs, df, tails, or statistics. Partial boundaries must
      # also be stationary in the retained coordinates.
      validate_probe <- function(location, value, information, gradient,
                                 probe_information_at, partial) {
        checks <<- checks + 1L
        indices <- if (partial) information$indices else seq_along(location)
        scaled_gradient <- forwardsolve(t(information$chol), gradient)
        correction <- as.vector(backsolve(information$chol, scaled_gradient))
        decrement <<- sqrt(sum(scaled_gradient^2))
        movement <<- max(abs(correction))
        if (partial && (decrement > 1e-4 || movement > 1e-3))
          stop("partial boundary stationarity validation failed")
        label <- if (partial) "partial boundary " else ""
        current <- inference_at(location, information)
        probe_location <- location
        probe_location[indices] <- location[indices] - correction
        if (any(probe_location[indices] <= lower[indices] & is.finite(lower[indices])))
          stop("infeasible ", if (partial) label else "forward ", "inference probe")
        probe_objective_change <<- value_at(probe_location) - value
        if (!is.finite(probe_objective_change) || probe_objective_change >
            sqrt(.Machine$double.eps) * (1 + abs(value)))
          stop(label, "probe objective validation failed")
        probe_information <- probe_information_at(probe_location)
        if (partial && !identical(probe_information$indices, information$indices))
          stop("partial boundary probe changed zero columns")
        probe <- inference_at(probe_location, probe_information)
        stability <<- max(abs(log(probe$values$std.error / current$values$std.error)),
          abs(log(probe$values$df / current$values$df)))
        tail_change <<- max(abs(probe$values$p.lower - current$values$p.lower))
        statistic_change <<- max(abs(probe$values$estimate / probe$values$std.error -
          current$values$estimate / current$values$std.error))
        list(ok = stability <= 1e-3 && tail_change <= 1e-5 && statistic_change <= 1e-4,
          inference = current)
      }
      validate_interior <- function(location, H, value, gradient) {
        validate_probe(location, value, list(H = H, chol = chol(H)), gradient,
                       information_at, partial = FALSE)
      }
      validate_boundary <- function(location, H, value, gradient) {
        if (all(location == 0)) {
          blmm_validate_H(H)
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
            location, sigma2, reTrms, context$random))
          return(list(ok = TRUE, inference = current))
        }
        information <- partial_information_at(location)
        validate_probe(location, value, information,
                       gradient_at(location)[information$indices],
                       partial_information_at, partial = TRUE)
      }
      optimized <- blmm_bounded_optimize(
        list(fn = value_at, gr = gradient_at, he = raw_hessian_at),
        initial$anchor$par, lower, max_steps = 8L,
        validate_interior = validate_interior, validate_boundary = validate_boundary)
      iterations <- optimized$iterations
      evaluations <- optimized$diagnostics$evaluations
      diagnostics$status[s] <- optimized$status
      if (optimized$status == "boundary") {
        if (all(optimized$par == 0)) {
          diagnostics$status[s] <- "boundary_zero"
          diagnostics$reason[s] <- "stationary zero covariance"
        } else {
          diagnostics$status[s] <- "boundary_partial"
          diagnostics$reason[s] <- "stationary partial zero covariance"
        }
      } else if (optimized$status != "interior") stop(optimized$reason)
      current <- optimized$inference
      if (is.null(current)) stop("blmm: optimizer returned no validated inference")
      current
    }, error = identity)
    if (!is.null(stale_report_error)) stop(stale_report_error)
    for (field in c("iterations", "evaluations", "hessians", "checks", "objective_calls",
                    "gradient_calls", "inference_calls", "decrement", "movement", "stability",
                    "tail_change", "statistic_change", "probe_objective_change"))
      diagnostics[[field]][s] <- get(field)
    if (inherits(attempted, "error")) {
      if (strict) stop(attempted)
      diagnostics$rescued[s] <- TRUE
      diagnostics$reason[s] <- conditionMessage(attempted)
    } else {
      for (field in .blmm_fields) {
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
