# Natural relative Cholesky coordinates: lower bounds apply to the diagonal
# entries. Search curvature may be damped; returned/inference curvature is raw.
blmm_bounded_optimize <- function(obj, start, lower, max_steps = 8L,
                                  validate_interior = NULL) {
  location <- as.numeric(start)
  lower <- as.numeric(lower)
  if (length(lower) != length(location) || !length(location) ||
      any(!is.finite(location)) || anyNA(lower) || any(lower == Inf) ||
      length(max_steps) != 1L || !is.finite(max_steps) || max_steps < 0)
    stop("blmm: invalid bounded optimisation inputs")
  location <- pmax(location, lower)
  max_steps <- as.integer(max_steps)
  iterations <- evaluations <- hessians <- gradient_calls <- checks <- 0L
  value <- NA_real_
  H <- matrix(NA_real_, length(location), length(location))
  diagnostics <- list()
  finish <- function(status, reason, inference = NULL) {
    list(par = location, objective = value, H = H, status = status,
      iterations = iterations, reason = reason, inference = inference,
      diagnostics = c(list(evaluations = evaluations, hessians = hessians,
        gradient_calls = gradient_calls, checks = checks), diagnostics))
  }
  value_at <- function(x) {
    evaluations <<- evaluations + 1L
    as.numeric(obj$fn(x))
  }
  run <- tryCatch({
    value <- value_at(location)
    if (length(value) != 1L || !is.finite(value)) stop("nonfinite objective")
    for (attempt in 0:max_steps) {
      iterations <- attempt
      gradient_calls <- gradient_calls + 1L
      gradient <- as.vector(obj$gr(location))
      hessians <- hessians + 1L
      H <- as.matrix(obj$he(location))
      if (!identical(dim(H), c(length(location), length(location))) ||
          any(!is.finite(H)) || length(gradient) != length(location) ||
          any(!is.finite(gradient))) stop("nonfinite derivatives")
      H <- (H + t(H)) / 2
      active <- is.finite(lower) & location <= lower +
        1e-10 * (1 + abs(replace(lower, !is.finite(lower), 0)))
      projected_gradient <- gradient
      projected_gradient[active & gradient >= 0] <- 0
      eigenH <- eigen(H, symmetric = TRUE)
      escape <- NULL
      if (any(active) && max(abs(projected_gradient)) <= 1e-6) {
        checks <- checks + 1L
        # Check inward coordinates and both signs of raw curvature directions.
        # Projection is needed for mixed free/bounded covariance blocks.
        directions <- cbind(diag(length(location))[, active, drop = FALSE],
          eigenH$vectors, -eigenH$vectors)
        probe_scale <- 1e-3 * (1 + max(abs(location)))
        best <- value
        for (j in seq_len(ncol(directions))) {
          candidate <- pmax(lower, location + probe_scale * directions[, j])
          if (!any(candidate != location)) next
          proposed <- value_at(candidate)
          if (length(proposed) != 1L || !is.finite(proposed))
            return(finish("unresolved", "nonfinite boundary probe"))
          if (proposed < best -
              .Machine$double.eps * (1 + abs(value))) {
            best <- proposed
            escape <- candidate - location
          }
        }
        if (is.null(escape)) {
          # Negative raw curvature can be infeasible (e.g. an outward cone
          # direction). Conservatively leave such ambiguous cases to rescue.
          if (min(eigenH$values) < -1e-8 * max(1, max(abs(eigenH$values))))
            return(finish("unresolved", "boundary negative curvature"))
          return(finish("boundary", "stationary boundary"))
        }
      }
      information <- tryCatch(blmm_validate_H(H), error = identity)
      if (!any(active) && !inherits(information, "error")) {
        factor <- chol(information)
        correction <- as.vector(backsolve(factor,
          forwardsolve(t(factor), gradient)))
        decrement <- sqrt(sum(forwardsolve(t(factor), gradient)^2))
        movement <- max(abs(correction))
        diagnostics$decrement <- decrement
        diagnostics$movement <- movement
        if (decrement <= 1e-4 && movement <= 1e-3) {
          checks <- checks + 1L
          # A small Newton correction alone cannot distinguish a nearby
          # boundary optimum from an interior optimum. Probe each bound before
          # allowing interior inference, without changing the other coordinates.
          bounded <- which(is.finite(lower))
          best_boundary <- NULL
          best_value <- Inf
          indistinguishable <- FALSE
          tolerance <- .Machine$double.eps * (1 + abs(value))
          for (j in bounded) {
            candidate <- location
            candidate[j] <- lower[j]
            proposed <- value_at(candidate)
            if (length(proposed) != 1L || !is.finite(proposed))
              return(finish("unresolved", "nonfinite boundary objective"))
            if (proposed < best_value) {
              best_value <- proposed
              best_boundary <- candidate
            }
            if (abs(proposed - value) <= tolerance) indistinguishable <- TRUE
          }
          if (!is.null(best_boundary) && best_value < value - tolerance) {
            if (attempt == max_steps)
              return(finish("unresolved", "iteration budget exhausted"))
            location <- best_boundary
            value <- best_value
            next
          }
          if (indistinguishable)
            return(finish("unresolved", "boundary objective indistinguishable"))
          validation <- if (is.null(validate_interior)) TRUE else
            validate_interior(location, information, value, gradient)
          ok <- if (is.list(validation)) isTRUE(validation$ok) else isTRUE(validation)
          if (is.list(validation) && !is.null(validation$diagnostics))
            diagnostics <- utils::modifyList(diagnostics, validation$diagnostics)
          if (ok) return(finish("interior", "converged",
            if (is.list(validation)) validation$inference else NULL))
        }
      }
      if (attempt == max_steps)
        return(finish("unresolved", "iteration budget exhausted"))
      if (!is.null(escape)) {
        direction <- escape
      } else {
        # Strict positive definiteness is only required for inference. A
        # shifted Hessian supplies a descent search direction away from it.
        floor <- 1e-4 * max(1, max(abs(eigenH$values)))
        curvature <- pmax(eigenH$values, floor)
        direction <- -as.vector(eigenH$vectors %*%
          (crossprod(eigenH$vectors, gradient) / curvature))
        radius <- 1 + max(abs(location))
        if (max(abs(direction)) > radius)
          direction <- direction * radius / max(abs(direction))
        direction <- pmax(lower, location + direction) - location
        if (sum(gradient * direction) >= 0)
          direction <- pmax(lower, location - projected_gradient) - location
      }
      accepted <- FALSE
      for (halving in 0:12) {
        candidate <- pmax(lower, location + direction * 2^(-halving))
        if (!any(candidate != location)) next
        proposed <- value_at(candidate)
        if (is.finite(proposed) && proposed <= value) {
          accepted <- TRUE
          break
        }
      }
      if (!accepted) return(finish("unresolved", "line search stalled"))
      location <- candidate
      value <- proposed
    }
  }, error = function(e) finish("unresolved", conditionMessage(e)))
  run
}
