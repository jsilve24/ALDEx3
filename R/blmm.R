## BLMM: Batched Linear Mixed-Effects Model engine for ALDEx3.
##
## Internal implementation. Not exported. One batched profiled REML anchor fit
## per feature, draw-specific local covariance updates around that anchor, and
## exact GLS fixed-effect solves conditional on the updated covariance
## parameters. The approximation is only in the variance-component step.

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------------------------------------------------------------------------
# Parameterisation helpers
# ---------------------------------------------------------------------------

#' Map unconstrained phi to lme4 variance-component parameters theta.
#'
#' @details
#' lme4 constrains variance-component parameters theta >= 0 (lower > -Inf).
#' To allow unconstrained optimisation those entries are log-transformed:
#'   phi[j] = log(theta[j])  for constrained parameters
#'   phi[j] = theta[j]       for unconstrained parameters (e.g. correlations)
#' phi is the optimisation variable throughout BLMM; theta is only used when
#' building LambdaZt or reporting random-effect variances.
#'
#' @param phi numeric vector of unconstrained optimisation parameters
#' @param lower numeric vector of lower bounds from \code{reTrms$lower};
#'   entries \code{> -Inf} indicate log-transformed (constrained) parameters
#' @return numeric vector of lme4-scale variance-component parameters theta
#' @noRd
phi_to_theta_blmm <- function(phi, lower) {
  theta <- phi
  bounded <- lower > -Inf
  theta[bounded] <- exp(phi[bounded])
  theta
}

#' Map lme4 theta back to unconstrained phi.
#'
#' @details Inverse of \code{phi_to_theta_blmm}.
#'
#' @param theta numeric vector of lme4-scale parameters
#' @param lower numeric vector of lower bounds from \code{reTrms$lower}
#' @return numeric vector of unconstrained optimisation parameters phi
#' @noRd
theta_to_phi_blmm <- function(theta, lower) {
  phi <- theta
  bounded <- lower > -Inf
  phi[bounded] <- log(pmax(theta[bounded], .Machine$double.eps))
  phi
}

# ---------------------------------------------------------------------------
# Basis tensor for LambdaZt
# ---------------------------------------------------------------------------

#' Precompute LambdaZt basis tensors for the TMB kernel.
#'
#' @details
#' LambdaZt = Lambda(theta) \%*\% Zt is linear in theta because
#' \code{Lambdat@x = theta[Lind]}. We precompute one basis slice per theta_j
#' by evaluating Lambda(e_j) \%*\% Zt, where e_j is the j-th unit vector.
#' Then at any theta: LambdaZt = sum_j theta_j * basis[j,,].
#' This factorisation lets the TMB kernel reconstruct LambdaZt cheaply without
#' re-parsing the lme4 formula structures in C++.
#'
#' @param reTrms random-effects terms list from \code{lme4::lFormula()}
#' @return array of dimension \code{(n_theta, q, n)}: basis[j,,] is the
#'   contribution of theta_j to LambdaZt
#' @noRd
compute_basis_blmm <- function(reTrms) {
  n_theta <- length(reTrms$theta)
  q <- nrow(reTrms$Zt)
  n <- ncol(reTrms$Zt)

  basis <- array(0, dim = c(n_theta, q, n))
  for (j in seq_len(n_theta)) {
    theta_ej <- rep(0, n_theta)
    theta_ej[j] <- 1
    Lam <- reTrms$Lambdat
    Lam@x <- theta_ej[reTrms$Lind]
    basis[j, , ] <- as.matrix(Lam %*% reTrms$Zt)
  }
  basis
}

# ---------------------------------------------------------------------------
# Random-effects covariance packaging
# ---------------------------------------------------------------------------

#' Convert lme4 Cholesky entries and sigma2 to named variance components.
#'
#' @details
#' lme4 stores random-effect covariance as a scaled Cholesky factor: the
#' covariance block for group g is sigma^2 * L_g \%*\% t(L_g), where L_g is a
#' lower-triangular matrix whose entries are the corresponding theta values.
#' This function converts those entries and the profiled sigma^2 back into
#' interpretable variance and covariance terms, using the same naming
#' convention as sr.mem/lme4: group.var for variances,
#' group.var1.var2 for covariances, and Residual.
#'
#' @param theta numeric vector of lme4 Cholesky factor entries
#' @param sigma2 profiled residual variance estimate
#' @param reTrms random-effects terms list from \code{lme4::lFormula()}
#' @return named numeric vector of variance and covariance components
#' @noRd
blmm_random_effect_vector <- function(theta, sigma2, reTrms) {
  cnms <- reTrms$cnms
  if (length(cnms) == 0L) {
    return(stats::setNames(sigma2, "Residual"))
  }

  block_sizes <- lengths(cnms, use.names = FALSE)
  block_npar <- block_sizes * (block_sizes + 1L) / 2L
  theta_split <- split(theta, rep(seq_along(block_sizes), block_npar))

  values <- numeric(0)
  names_out <- character(0)

  for (g in seq_along(cnms)) {
    grp_name <- names(cnms)[g] %||% paste0("grp", g)
    term_names <- cnms[[g]]
    block_size <- block_sizes[g]
    block_theta <- theta_split[[g]]

    L <- matrix(0, nrow = block_size, ncol = block_size)
    L[lower.tri(L, diag = TRUE)] <- block_theta
    Sigma <- sigma2 * tcrossprod(L)

    for (i in seq_len(block_size)) {
      values <- c(values, Sigma[i, i])
      names_out <- c(names_out, paste(grp_name, term_names[i], sep = "."))
    }

    if (block_size > 1L) {
      for (j in seq_len(block_size - 1L)) {
        for (i in seq.int(j + 1L, block_size)) {
          values <- c(values, Sigma[i, j])
          names_out <- c(
            names_out,
            paste(grp_name, term_names[j], term_names[i], sep = ".")
          )
        }
      }
    }
  }

  stats::setNames(c(values, sigma2), c(names_out, "Residual"))
}

# ---------------------------------------------------------------------------
# Formula parsing
# ---------------------------------------------------------------------------

#' Parse a mixed-effects formula and extract all objects needed by BLMM.
#'
#' @details
#' Uses a dummy zero response to drive \code{lme4::lFormula()} purely for its
#' side-effect of constructing X, Zt, and Lambdat. The actual response (logW)
#' is passed separately. lme4 validity checks are suppressed because the dummy
#' data dimensions would otherwise be rejected.
#'
#' @param formula lme4 mixed-effects formula
#' @param data data.frame used for formula evaluation
#' @return list with X, reTrms, basis, phi_init, is_log, lower, random_template
#' @noRd
blmm_parse_formula <- function(formula, data) {
  data$.blmm_y <- 0
  parsed <- lme4::lFormula(
    stats::update(formula, .blmm_y ~ .),
    data = data,
    control = lme4::lmerControl(
      check.nobs.vs.nRE = "ignore",
      check.nobs.vs.rankZ = "ignore"
    )
  )

  reTrms <- parsed$reTrms
  lower <- reTrms$lower

  list(
    X = parsed$X,
    reTrms = reTrms,
    basis = compute_basis_blmm(reTrms),
    phi_init = theta_to_phi_blmm(reTrms$theta, lower),
    is_log = as.integer(lower > -Inf),
    lower = lower,
    random_template = blmm_random_effect_vector(reTrms$theta, 1, reTrms)
  )
}

# ---------------------------------------------------------------------------
# TMB object construction and profiled objective helpers
# ---------------------------------------------------------------------------

#' Construct a TMB autodiff object for the batched profiled REML objective.
#'
#' @param X fixed-effects design matrix (n x p)
#' @param Y_d response matrix for one feature (n x S), draws stacked as columns
#' @param basis precomputed basis tensor from \code{compute_basis_blmm()}
#' @param is_log integer vector flagging which phi entries are log-transformed
#' @param phi_init starting values for the optimisation
#' @return TMB autodiff object with \code{$fn}, \code{$gr}, \code{$he},
#'   and \code{$report()} methods
#' @noRd
blmm_make_adfun <- function(X, Y_d, basis, is_log, phi_init) {
  TMB::MakeADFun(
    data = list(
      X = X,
      Y = Y_d,
      basis = basis,
      is_log = is_log
    ),
    parameters = list(phi = phi_init),
    DLL = "ALDEx3",
    silent = TRUE
  )
}

#' Per-draw profiled REML objective values.
#'
#' @details
#' The shared log-determinant terms (ldL2 + ldRX2) are constant across draws;
#' only the PWRSS term differs per draw since each draw has its own response.
#'
#' @param report list returned by \code{obj$report()} after evaluating the TMB
#'   objective at a given phi
#' @param n number of observations
#' @param p number of fixed-effect columns
#' @return numeric vector of length S: per-draw profiled REML objective values
#' @noRd
blmm_profiled_objectives <- function(report, n, p) {
  np <- n - p
  common <- 0.5 * (report$ldL2 + report$ldRX2)
  common + np / 2 * log(report$PWRSS_s / np)
}

#' Average profiled REML objective over all MC draws.
#'
#' @details
#' This is the quantity the TMB kernel minimises and what
#' \code{blmm_fit_anchor()} optimises.
#'
#' @inheritParams blmm_profiled_objectives
#' @return scalar: mean of per-draw profiled REML objectives
#' @noRd
blmm_anchor_objective_from_report <- function(report, n, p) {
  mean(blmm_profiled_objectives(report, n, p))
}

# ---------------------------------------------------------------------------
# Anchor optimisation and curvature handling
# ---------------------------------------------------------------------------

#' Symmetrise and regularise the observed Hessian at the anchor.
#'
#' @details
#' Adds a diagonal ridge just large enough to make the smallest eigenvalue
#' exceed \code{tol}. Small negative eigenvalues from numerical noise are
#' handled conservatively and explicitly rather than silently ignored.
#' Sets \code{attr(H, "ridge_delta")} when regularisation is applied so
#' callers can detect and diagnose near-flat curvature.
#'
#' @param H square numeric matrix: raw observed Hessian from \code{obj$he()}
#' @param tol minimum acceptable eigenvalue; values at or below trigger ridge
#' @return symmetric positive-definite matrix of the same dimension as H
#' @noRd
blmm_regularise_H <- function(H, tol = 1e-8) {
  H <- (H + t(H)) / 2
  if (any(!is.finite(H))) {
    stop("blmm: observed Hessian contains non-finite values")
  }

  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min_ev <- min(ev)
  if (min_ev <= tol) {
    delta <- tol - min_ev + max(abs(diag(H)), 1) * 1e-8
    H <- H + diag(delta, nrow(H))
    attr(H, "ridge_delta") <- delta
  }

  H
}

#' Optimise the batched profiled REML objective and return the anchor.
#'
#' @details
#' Uses \code{stats::nlminb()} with TMB-computed gradients. Errors and
#' non-convergence are promoted to hard stops so they propagate cleanly to the
#' per-feature fallback handler in \code{blmm_fit_feature()}.
#'
#' @param obj TMB autodiff object from \code{blmm_make_adfun()}
#' @param phi_init starting values for the optimiser
#' @return list with \code{phi_bar} (optimal phi) and \code{H_d} (regularised
#'   observed Hessian at the optimum)
#' @noRd
blmm_fit_anchor <- function(obj, phi_init) {
  fit <- tryCatch(
    stats::nlminb(
      start = phi_init,
      objective = obj$fn,
      gradient = obj$gr,
      control = list(eval.max = 1000L, iter.max = 500L, rel.tol = 1e-10)
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    stop("blmm anchor optimisation failed: ", conditionMessage(fit))
  }
  if (!identical(fit$convergence, 0L)) {
    stop("blmm anchor optimisation did not converge (code ", fit$convergence,
         "): ", fit$message %||% "no message")
  }

  H_raw <- tryCatch(obj$he(fit$par), error = function(e) e)
  if (inherits(H_raw, "error")) {
    stop("blmm Hessian evaluation failed at the anchor: ",
         conditionMessage(H_raw))
  }

  list(phi_bar = fit$par, H_d = blmm_regularise_H(H_raw))
}

# ---------------------------------------------------------------------------
# Per-draw score corrections around the anchor
# ---------------------------------------------------------------------------

#' Compute per-draw score corrections at the anchor.
#'
#' @details
#' The anchor \code{phi_bar} minimises the \emph{average} profiled objective
#' over all S draws. For draw s, the draw-specific objective f_s(phi) differs
#' from the average only through its PWRSS_s term. The gradient of f_s at
#' phi_bar is the "score" for draw s: the direction and magnitude by which that
#' draw wants to move away from the shared anchor.
#'
#' Scores are computed via centred finite differences of log(PWRSS_s). The
#' ldL2/ldRX2 log-determinant terms are identical for all draws at any fixed
#' phi, so their derivative cancels out of the per-draw correction:
#' \deqn{score_s[j] = \frac{n-p}{2}
#'   \left(\frac{\partial \log PWRSS_s}{\partial \phi_j}
#'         - \overline{\frac{\partial \log PWRSS}{\partial \phi_j}}\right)}
#'
#' @param obj TMB autodiff object, with \code{obj$fn} already evaluated at
#'   \code{phi_bar} so that \code{obj$report()} reflects the anchor state
#' @param phi_bar numeric vector: the anchor (optimal phi for the average objective)
#' @param PWRSS_s numeric vector of length S: per-draw penalised residual SS at
#'   the anchor, from \code{obj$report()$PWRSS_s}
#' @param n number of observations
#' @param p number of fixed-effect columns
#' @param eps relative step size for centred finite differences
#' @return matrix of dimension \code{(length(phi_bar), S)}: column s is the
#'   score vector for draw s
#' @noRd
blmm_scores_batch <- function(obj, phi_bar, PWRSS_s, n, p, eps = 1e-4) {
  if (any(!is.finite(PWRSS_s)) || any(PWRSS_s <= 0)) {
    stop("blmm: invalid profiled residual sum of squares at the anchor")
  }

  S <- length(PWRSS_s)
  np <- n - p
  q_phi <- length(phi_bar)
  scores <- matrix(0, nrow = q_phi, ncol = S)

  for (j in seq_len(q_phi)) {
    step <- eps * max(1, abs(phi_bar[j]))
    phi_p <- phi_bar
    phi_m <- phi_bar
    phi_p[j] <- phi_p[j] + step
    phi_m[j] <- phi_m[j] - step

    obj$fn(phi_p)
    rp <- obj$report()
    obj$fn(phi_m)
    rm <- obj$report()

    PWRSS_p <- rp$PWRSS_s
    PWRSS_m <- rm$PWRSS_s
    if (any(!is.finite(PWRSS_p)) || any(!is.finite(PWRSS_m)) ||
        any(PWRSS_p <= 0) || any(PWRSS_m <= 0)) {
      stop("blmm: score perturbation produced an invalid profiled residual sum of squares")
    }

    # Relative step in the profiled parameter scale; centred finite differences
    # are stable here because the variance-component dimension is small.
    dlogPWRSS_s <- (log(PWRSS_p) - log(PWRSS_m)) / (2 * step)
    scores[j, ] <- np / 2 * (mean(dlogPWRSS_s) - dlogPWRSS_s)
  }

  scores
}

# ---------------------------------------------------------------------------
# Newton update around the anchor
# ---------------------------------------------------------------------------

#' Apply a Newton step to get per-draw covariance parameters.
#'
#' @details
#' Each draw s has a score vector \code{scores[, s]} at the shared anchor.
#' A single Newton step using the anchor's Hessian gives the draw-specific
#' optimal covariance parameters:
#' \deqn{\tilde{\phi}_s = \bar{\phi} + H_d^{-1} \cdot \text{scores}_s}
#' This is the core approximation: one shared curvature is reused rather than
#' optimising separately for each of the S draws.
#'
#' @param phi_bar numeric vector: anchor (shared optimal phi)
#' @param H_d positive-definite matrix: regularised Hessian at the anchor
#' @param scores matrix (length(phi_bar) x S): per-draw score vectors from
#'   \code{blmm_scores_batch()}
#' @return matrix of the same dimension as \code{scores}: column s is the
#'   draw-specific updated phi
#' @noRd
blmm_phi_updates <- function(phi_bar, H_d, scores) {
  cholH <- tryCatch(chol(H_d), error = function(e) e)
  if (inherits(cholH, "error")) {
    stop("blmm: curvature matrix is not numerically positive definite")
  }

  H_inv_scores <- backsolve(cholH, forwardsolve(t(cholH), scores))
  phi_bar + H_inv_scores
}

# ---------------------------------------------------------------------------
# Exact conditional fixed-effect solves
# ---------------------------------------------------------------------------

#' Exact GLS fixed-effect solve for one MC draw.
#'
#' @details
#' Given the precomputed Cholesky pieces from \code{blmm_common_pieces()},
#' this is an exact GLS solve — no approximation. Standard errors come from
#' the inverse of the Schur complement M = X'X - RZX'RZX (equivalent to
#' the lme4 RX factor). Residual df is n - p; Satterthwaite approximation
#' is not used.
#'
#' @param pieces list from \code{blmm_common_pieces()}: L, RZX, LX, LambdaZt
#' @param X fixed-effects design matrix (n x p)
#' @param y numeric vector of length n: response for this draw
#' @return list with beta, se, df, sigma2
#' @noRd
blmm_fixed_effects_draw <- function(pieces, X, y) {
  n <- nrow(X)
  p <- ncol(X)

  L <- pieces$L
  RZX <- pieces$RZX
  LX <- pieces$LX
  LambdaZt <- pieces$LambdaZt

  cu <- forwardsolve(L, LambdaZt %*% y)
  Cu_y <- crossprod(X, y) - crossprod(RZX, cu)
  Cbeta <- forwardsolve(LX, Cu_y)
  beta <- backsolve(t(LX), Cbeta)

  PWRSS <- as.numeric(sum(y^2) - sum(cu^2) - sum(Cbeta^2))
  if (!is.finite(PWRSS) || PWRSS <= 0) {
    stop("blmm: conditional fixed-effect solve produced a non-positive profiled residual sum of squares")
  }

  sigma2 <- PWRSS / (n - p)
  if (!is.finite(sigma2) || sigma2 <= 0) {
    stop("blmm: conditional residual variance is not positive")
  }

  Minv <- chol2inv(t(LX))
  diag_m <- diag(Minv)
  if (any(diag_m < -1e-10)) {
    stop("blmm: conditional fixed-effect covariance matrix is not positive semidefinite")
  }

  se <- sqrt(pmax(sigma2 * diag_m, 0))

  list(
    beta = as.vector(beta),
    se = as.vector(se),
    df = n - p,
    sigma2 = sigma2
  )
}

#' Recompute Cholesky pieces for a per-draw GLS solve.
#'
#' @details
#' Once per-draw covariance parameters phi_tilde[, s] are available, we need
#' L, RZX, LX, and LambdaZt to run the exact conditional solve. These depend
#' only on phi, X, and the basis tensors — not on Y — so we recompute them
#' directly in R rather than re-evaluating the full batched TMB objective,
#' which would perform unnecessary Y-block work for all S draws at once.
#'
#' @param phi numeric vector: draw-specific unconstrained covariance parameters
#' @param X fixed-effects design matrix (n x p)
#' @param basis precomputed basis tensor from \code{compute_basis_blmm()}
#' @param lower numeric vector of lower bounds from \code{reTrms$lower}
#' @return list with L (lower Cholesky of A = I + LambdaZt LambdaZt'),
#'   RZX (L^{-1} LambdaZt X), LX (lower Cholesky of Schur complement M),
#'   and LambdaZt
#' @noRd
blmm_common_pieces <- function(phi, X, basis, lower) {
  theta <- phi_to_theta_blmm(phi, lower)
  q <- dim(basis)[2]
  n <- nrow(X)
  LambdaZt <- matrix(0, nrow = q, ncol = n)

  for (j in seq_along(theta)) {
    basis_j <- matrix(basis[j, , ], nrow = q, ncol = n)
    LambdaZt <- LambdaZt + theta[j] * basis_j
  }

  A <- LambdaZt %*% t(LambdaZt)
  diag(A) <- diag(A) + 1
  L_raw <- tryCatch(chol(A), error = function(e) e)
  if (inherits(L_raw, "error")) {
    stop("blmm: local covariance factorisation failed")
  }
  L <- t(L_raw)

  RZX <- forwardsolve(L, LambdaZt %*% X)
  M <- crossprod(X) - crossprod(RZX)
  LX_raw <- tryCatch(chol(M), error = function(e) e)
  if (inherits(LX_raw, "error")) {
    stop("blmm: local fixed-effect factorisation failed")
  }

  dimnames(L) <- NULL
  dimnames(RZX) <- NULL
  dimnames(LX_raw) <- NULL
  dimnames(LambdaZt) <- NULL

  list(
    L = L,
    RZX = RZX,
    LX = t(LX_raw),
    LambdaZt = LambdaZt
  )
}

# ---------------------------------------------------------------------------
# Exact lme4 fallback
# ---------------------------------------------------------------------------

#' Exact lme4 fallback fit for one feature (all S draws).
#'
#' @details
#' Called when the approximate BLMM path fails for a feature. Wraps
#' \code{sr.mem()} with a singleton feature dimension so the output layout
#' matches what \code{blmm_fit_feature()} expects.
#'
#' @param Y_d numeric matrix (n x S): log-composition draws for this feature
#' @param formula lme4 mixed-effects formula
#' @param data data.frame for formula evaluation
#' @return list with estimate, std.error, df, p.lower, p.upper, random.eff
#'   arrays — each retaining a singleton feature dimension (D = 1)
#' @noRd
blmm_exact_feature <- function(Y_d, formula, data) {
  exact <- sr.mem(
    logW = array(Y_d, dim = c(nrow(Y_d), 1L, ncol(Y_d))),
    formula = formula,
    data = data,
    n.cores = 1L,
    method = "lme4",
    mem.args = list()
  )

  list(
    estimate = exact$estimate[, 1, , drop = FALSE],
    std.error = exact$std.error[, 1, , drop = FALSE],
    df = exact$df[, 1, , drop = FALSE],
    p.lower = exact$p.lower[, 1, , drop = FALSE],
    p.upper = exact$p.upper[, 1, , drop = FALSE],
    random.eff = exact$random.eff[, 1, , drop = FALSE]
  )
}

#' Approximate BLMM fit for one feature across all S MC draws.
#'
#' @details
#' Runs the full approximate pipeline: (1) optimise the average-over-draws
#' profiled REML objective to get the shared anchor and its Hessian;
#' (2) compute per-draw score corrections via centred finite differences and
#' apply a Newton step to get draw-specific covariance parameters phi_tilde;
#' (3) for each draw, recompute the Cholesky pieces from phi_tilde[, s] and
#' run an exact GLS solve for fixed effects.
#'
#' On any numerical failure the entire feature falls back silently to exact
#' lme4 via \code{blmm_exact_feature()}, and \code{fallback = TRUE} is
#' returned so the caller can warn the user.
#'
#' @param d integer feature index (column of logW)
#' @param logW numeric array (N x D x S)
#' @param X fixed-effects design matrix
#' @param basis precomputed basis tensor
#' @param is_log integer vector for phi/theta transformation
#' @param phi_init starting values
#' @param lower lower bounds for theta
#' @param reTrms random-effects terms list
#' @param formula lme4 formula
#' @param data data.frame
#' @param forced_fallback integer vector of feature indices to force onto the
#'   exact lme4 path; a testing hook, not used in normal operation
#' @return list with \code{fit} (arrays), \code{fallback} (logical), and
#'   \code{fallback_message} (character)
#' @noRd
blmm_fit_feature <- function(d, logW, X, basis, is_log, phi_init, lower,
                             reTrms, formula, data, forced_fallback = integer(0)) {
  N <- nrow(X)
  p <- ncol(X)
  S <- dim(logW)[3]
  Y_d <- logW[, d, , drop = FALSE]
  dim(Y_d) <- c(N, S)

  feature_random_rows <- length(names(blmm_random_effect_vector(reTrms$theta, 1, reTrms)))
  if (length(forced_fallback) > 0L && d %in% forced_fallback) {
    return(list(
      fit = blmm_exact_feature(Y_d, formula, data),
      fallback = TRUE,
      fallback_message = "forced approximate failure for testing"
    ))
  }

  result <- tryCatch({
    # Step 1: anchor — optimise the average-over-draws objective.
    obj_full <- blmm_make_adfun(X, Y_d, basis, is_log, phi_init)
    anchor <- blmm_fit_anchor(obj_full, phi_init)

    # Step 2: score corrections — differentiate per-draw PWRSS w.r.t. phi at
    # the anchor, then Newton-step to draw-specific covariance parameters.
    obj_full$fn(anchor$phi_bar)
    anchor_report <- obj_full$report()
    g_mat <- blmm_scores_batch(obj_full, anchor$phi_bar,
                               anchor_report$PWRSS_s, N, p)
    phi_tilde <- blmm_phi_updates(anchor$phi_bar, anchor$H_d, g_mat)

    feature_estimate <- matrix(NA_real_, nrow = p, ncol = S)
    feature_std.error <- matrix(NA_real_, nrow = p, ncol = S)
    feature_df <- matrix(NA_real_, nrow = p, ncol = S)
    feature_p.lower <- matrix(NA_real_, nrow = p, ncol = S)
    feature_p.upper <- matrix(NA_real_, nrow = p, ncol = S)
    feature_random <- matrix(NA_real_, nrow = feature_random_rows, ncol = S)

    # Step 3: exact conditional GLS solve per draw. p-values use a
    # t-distribution with residual df (n - p), not Satterthwaite.
    for (s in seq_len(S)) {
      pieces <- blmm_common_pieces(phi_tilde[, s], X, basis, lower)
      fe <- blmm_fixed_effects_draw(pieces, X, Y_d[, s])

      feature_estimate[, s] <- fe$beta
      feature_std.error[, s] <- fe$se
      feature_df[, s] <- fe$df

      t_stat <- fe$beta / fe$se
      feature_p.lower[, s] <- pt(t_stat, df = fe$df, lower.tail = TRUE)
      feature_p.upper[, s] <- 1 - feature_p.lower[, s]

      theta_draw <- phi_to_theta_blmm(phi_tilde[, s], lower)
      feature_random[, s] <- blmm_random_effect_vector(
        theta = theta_draw,
        sigma2 = fe$sigma2,
        reTrms = reTrms
      )
    }

    list(
      fit = list(
        estimate = array(feature_estimate, c(p, 1L, S)),
        std.error = array(feature_std.error, c(p, 1L, S)),
        df = array(feature_df, c(p, 1L, S)),
        p.lower = array(feature_p.lower, c(p, 1L, S)),
        p.upper = array(feature_p.upper, c(p, 1L, S)),
        random.eff = array(feature_random, c(feature_random_rows, 1L, S))
      ),
      fallback = FALSE,
      fallback_message = ""
    )
  }, error = function(e) {
    list(
      fit = blmm_exact_feature(Y_d, formula, data),
      fallback = TRUE,
      fallback_message = conditionMessage(e)
    )
  })

  result
}

# ---------------------------------------------------------------------------
# Main engine
# ---------------------------------------------------------------------------

#' Fast approximate mixed-effects engine for ALDEx3.
#'
#' @details
#' Internal engine called by \code{\link{aldex}} when \code{method = "blmm"}.
#'
#' Uses an \code{lme4}-style profiled Gaussian LMM formulation. One batched
#' anchor REML problem is fit per feature across all S MC draws; draw-specific
#' local covariance updates are obtained via a single Newton step; exact
#' conditional GLS fixed-effect solves are then run per draw. The approximation
#' is only in the variance-component step — fixed effects and standard errors
#' are exact conditional on the updated covariance parameters.
#'
#' Features for which the approximate path fails are silently re-fit with the
#' exact \code{lme4} engine and a consolidated warning is issued naming the
#' affected features and their error messages.
#'
#' @param logW numeric array (N x D x S)
#' @param formula lme4 mixed-effects formula
#' @param data data.frame for formula evaluation
#' @param n.cores number of worker processes for feature-level parallelism.
#'   If \code{n.cores > 1}, features are distributed across workers; the
#'   per-feature draw loop remains serial.
#' @return named list of arrays (P x D x S): estimate, std.error, df,
#'   p.lower, p.upper, and (Pr x D x S) random.eff — same layout as
#'   \code{sr.mem}.
#' @importFrom lme4 lFormula lmerControl
#' @importFrom TMB MakeADFun
#' @useDynLib ALDEx3
blmm <- function(logW, formula, data, n.cores = 1L) {
  N <- dim(logW)[1]
  D <- dim(logW)[2]
  S <- dim(logW)[3]

  parsed <- blmm_parse_formula(formula, data)
  X <- parsed$X
  basis <- parsed$basis
  is_log <- parsed$is_log
  phi_init <- parsed$phi_init
  lower <- parsed$lower
  reTrms <- parsed$reTrms

  p <- ncol(X)
  fnames <- colnames(X)
  rnames <- names(parsed$random_template)
  Pr <- length(rnames)

  # Testing hook: set attr(logW, "blmm.force_fail_features") to an integer
  # vector to force specific features onto the exact lme4 path, exercising
  # the fallback warning without needing pathological data.
  forced_fallback <- attr(logW, "blmm.force_fail_features")
  if (is.null(forced_fallback)) {
    forced_fallback <- integer(0)
  } else {
    forced_fallback <- unique(as.integer(forced_fallback))
    forced_fallback <- forced_fallback[forced_fallback >= 1L & forced_fallback <= D]
  }

  estimate <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  std.error <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  df_arr <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  p.lower <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  p.upper <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  random.eff <- array(NA_real_, c(Pr, D, S), dimnames = list(rnames, NULL, NULL))

  # Closure over the shared parsed objects so parallel dispatch only passes
  # the feature index d. On Windows (PSOCK cluster) this avoids serialising
  # the large logW array per task; on POSIX mclapply uses forking instead.
  feature_worker <- function(d) {
    blmm_fit_feature(
      d = d,
      logW = logW,
      X = X,
      basis = basis,
      is_log = is_log,
      phi_init = phi_init,
      lower = lower,
      reTrms = reTrms,
      formula = formula,
      data = data,
      forced_fallback = forced_fallback
    )
  }

  feature_indices <- seq_len(D)
  if (n.cores > 1L && D > 1L) {
    n_workers <- min(as.integer(n.cores), D)
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(n_workers)
      on.exit(stopCluster(cl), add = TRUE)
      clusterEvalQ(cl, {
        suppressPackageStartupMessages(library(ALDEx3))
        NULL
      })
      feature_results <- parLapply(cl, feature_indices, feature_worker)
    } else {
      feature_results <- parallel::mclapply(
        feature_indices,
        feature_worker,
        mc.cores = n_workers
      )
    }
  } else {
    feature_results <- lapply(feature_indices, feature_worker)
  }

  # Any feature that threw an error was already re-fit by blmm_exact_feature()
  # inside the tryCatch in blmm_fit_feature(). Emit one consolidated warning
  # with up to 3 example error messages so the caller knows what failed.
  fallback_idx <- which(vapply(feature_results, `[[`, logical(1), "fallback"))
  if (length(fallback_idx) > 0L) {
    n_examples <- min(3L, length(fallback_idx))
    example_msgs <- vapply(
      feature_results[fallback_idx[seq_len(n_examples)]],
      `[[`,
      "",
      "fallback_message"
    )
    example_text <- paste(
      sprintf("feature %d: %s", fallback_idx[seq_len(n_examples)], example_msgs),
      collapse = "; "
    )
    suffix <- if (length(fallback_idx) > n_examples) " ..." else ""
    warning(
      sprintf(
        "blmm: %d feature(s) fell back to exact lme4 because the approximate path failed. Examples: %s%s",
        length(fallback_idx), example_text, suffix
      ),
      call. = FALSE
    )
  }

  # Each blmm_fit_feature() result has a singleton feature dimension (D = 1);
  # drop it when writing into the D-wide output arrays.
  for (d in seq_len(D)) {
    feature_fit <- feature_results[[d]]$fit
    estimate[, d, ] <- feature_fit$estimate[, 1, ]
    std.error[, d, ] <- feature_fit$std.error[, 1, ]
    df_arr[, d, ] <- feature_fit$df[, 1, ]
    p.lower[, d, ] <- feature_fit$p.lower[, 1, ]
    p.upper[, d, ] <- feature_fit$p.upper[, 1, ]
    random.eff[, d, ] <- feature_fit$random.eff[, 1, ]
  }

  list(
    estimate = estimate,
    std.error = std.error,
    df = df_arr,
    p.lower = p.lower,
    p.upper = p.upper,
    random.eff = random.eff
  )
}
