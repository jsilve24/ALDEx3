# Adaptive covariance refinement. The forward checks and eight-step budget
# reproduce diagnostics/deniz-runtime/curvature-candidate.R (refresh strategy).
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
    reason = rep("converged", S), iterations = integer(S), evaluations = integer(S),
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
    batch <- blmm_make_adfun(X, residual, basis, is_log, phi_init, response)
    anchor <- blmm_fit_anchor(batch, phi_init)
    scores <- blmm_scores_batch(batch, anchor$phi_bar,
      batch$report(anchor$phi_bar)$PWRSS_s, N, p)
    cholH <- chol(anchor$H_d)
    updates <- blmm_phi_updates(anchor$phi_bar, anchor$H_d, scores, cholH)
    list(anchor = anchor, cholH = cholH, updates = updates,
         residual = residual, offsets = offsets, response = response)
  }, error = identity)
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
      objective <- blmm_make_adfun(X, residual_s, basis, is_log, phi_init, response_s)
      last_value <- last_gradient <- last_information <- last_inference <- NULL
      value_at <- function(location) {
        if (!is.null(last_value) && identical(location, last_value$location))
          return(last_value$value)
        objective_calls <<- objective_calls + 1L
        value <- objective$fn(location)
        last_value <<- list(location = location, value = value)
        value
      }
      information_at <- function(location) {
        if (!is.null(last_information) && identical(location, last_information$location))
          return(last_information$information)
        hessians <<- hessians + 1L
        H <- blmm_validate_H(objective$he(location))
        information <- list(H = H, chol = chol(H))
        last_information <<- list(location = location, information = information)
        information
      }
      correction_at <- function(information, gradient)
        as.vector(backsolve(information$chol, forwardsolve(t(information$chol), gradient)))
      scaled_gradient <- function(information, gradient)
        sqrt(sum(forwardsolve(t(information$chol), gradient)^2))
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
        pieces <- blmm_common_pieces(location, X, basis, lower)
        fixed <- blmm_fixed_effects_draw(pieces, X, residual_s, response_s)
        estimates <- as.vector(crossprod(vectors, fixed$beta + initial$offsets[, s]))
        standard_errors <- sqrt(fixed$sigma2 * colSums(forwardsolve(pieces$LX, vectors)^2))
        values <- list(estimate = estimates, std.error = standard_errors, df = dfs,
          p.lower = stats::pt(estimates / standard_errors, dfs))
        if (any(!is.finite(unlist(values))) || any(standard_errors <= 0))
          stop("nonfinite inference")
        values$p.upper <- stats::pt(estimates / standard_errors, dfs, lower.tail = FALSE)
        random <- blmm_random_effect_vector(phi_to_theta_blmm(location, lower),
          fixed$sigma2, reTrms, context$random)
        if (any(!is.finite(random))) stop("nonfinite random-effect inference")
        inference <- list(values = values, random = random)
        last_inference <<- list(location = location, H = information$H, inference = inference)
        inference
      }
      location <- initial$anchor$phi_bar
      working <- list(H = initial$anchor$H_d, chol = initial$cholH)
      gradient <- gradient_at(location)
      value <- value_at(location)
      converged <- FALSE
      for (iteration in seq_len(8L)) {
        iterations <- iteration
        checked_information <- NULL
        direction <- if (iteration == 1L) initial$updates[, s] - location else
          -correction_at(working, gradient)
        accepted <- FALSE
        for (halving in 0:12) {
          candidate <- location + direction * 2^(-halving)
          proposed_value <- value_at(candidate)
          evaluations <- evaluations + 1L
          if (is.finite(proposed_value) && proposed_value <= value) {
            accepted <- TRUE
            break
          }
        }
        if (accepted) { location <- candidate; value <- proposed_value }
        gradient <- gradient_at(location)
        working_correction <- correction_at(working, gradient)
        working_decrement <- scaled_gradient(working, gradient)
        if (working_decrement <= 1e-4 && max(abs(working_correction)) <= 1e-3) {
          checks <- checks + 1L
          information <- information_at(location)
          checked_information <- information
          correction <- correction_at(information, gradient)
          decrement <- scaled_gradient(information, gradient)
          movement <- max(abs(correction))
          if (decrement <= 1e-4 && movement <= 1e-3) {
            current <- inference_at(location, information)
            probe_location <- location - correction
            probe_objective_change <- value_at(probe_location) - value
            if (!is.finite(probe_objective_change) || probe_objective_change >
                sqrt(.Machine$double.eps) * (1 + abs(value)))
              stop("probe objective validation failed")
            probe <- inference_at(probe_location, information_at(probe_location))
            stability <- max(abs(log(probe$values$std.error / current$values$std.error)),
              abs(log(probe$values$df / current$values$df)))
            tail_change <- max(abs(probe$values$p.lower - current$values$p.lower))
            statistic_change <- max(abs(probe$values$estimate / probe$values$std.error -
              current$values$estimate / current$values$std.error))
            converged <- stability <= 1e-3 && tail_change <= 1e-5 && statistic_change <= 1e-4
          }
        }
        if (converged) break
        if (!accepted) stop("line search stalled")
        if (iteration < 8L) working <- if (is.null(checked_information))
          information_at(location) else checked_information
      }
      if (!converged) stop("iteration budget exhausted")
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
