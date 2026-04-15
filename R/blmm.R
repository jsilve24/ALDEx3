## BLMM: Batched Linear Mixed-Effects Model engine for ALDEx3.
## Internal — not exported.
##
## Supports all lme4 model structures (random intercepts, random slopes,
## correlated random effects, nested groupings) by using lme4's reTrms
## parameterisation via basis matrices dLambdaZt/dtheta_j.
##
## Speed strategy:
##   - TMB (C++ with automatic differentiation) for the anchor objective,
##     exact gradients, and exact Hessian.
##   - Batched q x q Cholesky (not n x n V): one factorisation per REML eval.
##   - Per-draw scores via the identity:
##       g_{ds,j} = (n-p)/2 * [d log(PWRSS_bar)/dphi_j - d log(PWRSS_s)/dphi_j]
##     evaluated by 2*len(phi) batched TMB function calls (no per-draw loop).
##   - Per-draw fixed effects use the same q x q structure.

# ---------------------------------------------------------------------------
# TMB compilation / loading (done once at package load)
# ---------------------------------------------------------------------------

.blmm_tmb_loaded <- FALSE

blmm_load_tmb <- function() {
  if (isTRUE(.blmm_tmb_loaded)) return(invisible(NULL))

  # Try installed package path first, then development (devtools::load_all)
  src <- system.file("src", "blmm.cpp", package = "ALDEx3")
  if (!nzchar(src) || !file.exists(src)) {
    pkg_path <- tryCatch(find.package("ALDEx3"), error = function(e) NULL)
    if (!is.null(pkg_path))
      src <- file.path(pkg_path, "src", "blmm.cpp")
  }
  if (!nzchar(src) || !file.exists(src))
    stop("blmm: cannot locate src/blmm.cpp. ",
         "Ensure the package is properly installed or loaded via devtools.")

  dll_name <- sub("\\.cpp$", "", basename(src))
  dll_dir  <- dirname(src)
  dll_path <- file.path(dll_dir, paste0(dll_name, .Platform$dynlib.ext))

  if (!file.exists(dll_path)) {
    message("blmm: compiling TMB template (first use only)...")
    TMB::compile(src, "-O2")
  }

  if (!(dll_name %in% names(getLoadedDLLs())))
    dyn.load(TMB::dynlib(file.path(dll_dir, dll_name)))

  utils::assignInMyNamespace(".blmm_tmb_loaded", TRUE)
  invisible(NULL)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------------------------------------------------------------------------
# Parameterisation helpers
# ---------------------------------------------------------------------------

phi_to_theta_blmm <- function(phi, lower) {
  theta <- phi
  bounded <- lower > -Inf
  theta[bounded] <- exp(phi[bounded])
  theta
}

theta_to_phi_blmm <- function(theta, lower) {
  phi <- theta
  bounded <- lower > -Inf
  phi[bounded] <- log(pmax(theta[bounded], .Machine$double.eps))
  phi
}

# ---------------------------------------------------------------------------
# Compute basis matrices: dLambdaZt/dtheta_j  (exact, linear map)
# ---------------------------------------------------------------------------
# Because Lambdat@x = theta[Lind] is a linear map, d(LambdaZt)/dtheta_j
# is the q x n matrix obtained by setting theta = e_j (unit vector) in the
# construction Lambdat -> Lambdat %*% Zt.
#
# Returns a 3D array [n_theta, q, n].

compute_basis_blmm <- function(reTrms) {
  n_theta <- length(reTrms$theta)
  q <- nrow(reTrms$Zt)
  n <- ncol(reTrms$Zt)

  basis <- array(0, dim = c(n_theta, q, n))
  for (j in seq_len(n_theta)) {
    theta_ej <- rep(0, n_theta); theta_ej[j] <- 1
    Lam <- reTrms$Lambdat
    Lam@x <- theta_ej[reTrms$Lind]
    basis[j, , ] <- as.matrix(Lam %*% reTrms$Zt)
  }
  basis
}

# ---------------------------------------------------------------------------
# Formula parsing: returns X, reTrms, basis, phi_init, is_log, re_varnames
# ---------------------------------------------------------------------------

blmm_parse_formula <- function(formula, data) {
  data$.blmm_y <- rnorm(nrow(data))
  parsed <- lme4::lFormula(
    stats::update(formula, .blmm_y ~ .),
    data    = data,
    control = lme4::lmerControl(
      check.nobs.vs.nRE   = "ignore",
      check.nobs.vs.rankZ = "ignore"
    )
  )

  reTrms  <- parsed$reTrms
  lower   <- reTrms$lower
  phi_init <- theta_to_phi_blmm(reTrms$theta, lower)
  is_log   <- as.integer(lower > -Inf)

  # Build basis matrices once
  basis <- compute_basis_blmm(reTrms)

  # Human-readable names for diagonal variance components
  re_varnames <- character(0)
  for (grp in names(reTrms$cnms))
    for (col in reTrms$cnms[[grp]])
      re_varnames <- c(re_varnames, paste0(grp, ".", col))
  diag_idx <- which(lower > -Inf)
  re_diag_names <- if (length(re_varnames) >= length(diag_idx))
                     re_varnames[seq_along(diag_idx)]
                   else
                     paste0("re_comp_", diag_idx)

  list(
    X            = parsed$X,
    reTrms       = reTrms,
    basis        = basis,
    phi_init     = phi_init,
    is_log       = is_log,
    re_diag_idx  = diag_idx,
    re_varnames  = re_diag_names
  )
}

# ---------------------------------------------------------------------------
# Build (or rebuild) a TMB ADFun object for one feature's response block.
# ---------------------------------------------------------------------------

blmm_make_adfun <- function(X, Y_d, basis, is_log, phi_init) {
  blmm_load_tmb()
  data_list <- list(
    X      = X,
    Y      = Y_d,
    basis  = basis,
    is_log = is_log
  )
  param_list <- list(phi = phi_init)
  TMB::MakeADFun(
    data       = data_list,
    parameters = param_list,
    DLL        = "blmm",
    silent     = TRUE
  )
}

# ---------------------------------------------------------------------------
# Anchor optimisation
# ---------------------------------------------------------------------------

blmm_regularise_H <- function(H) {
  H <- (H + t(H)) / 2
  ev <- tryCatch(
    eigen(H, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(-1, nrow(H))
  )
  if (min(ev) <= 0) {
    delta <- abs(min(ev)) + max(abs(diag(H))) * 1e-6 + 1e-8
    H <- H + delta * diag(nrow(H))
  }
  H
}

blmm_fit_anchor <- function(obj, phi_init) {
  fit <- tryCatch(
    stats::nlminb(
      start     = phi_init,
      objective = obj$fn,
      gradient  = obj$gr,
      control   = list(eval.max = 1000L, iter.max = 500L, rel.tol = 1e-10)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || fit$convergence > 1L) {
    warning("blmm: anchor optimisation did not converge; ",
            "results for this feature may be unreliable. ",
            "Convergence: ", if (is.null(fit)) "NULL" else fit$convergence)
  }

  phi_bar   <- if (is.null(fit)) phi_init else fit$par
  H_raw     <- tryCatch(obj$he(phi_bar), error = function(e) diag(length(phi_bar)))
  H_d       <- blmm_regularise_H(H_raw)
  converged <- !is.null(fit) && fit$convergence == 0L

  list(phi_bar = phi_bar, H_d = H_d, converged = converged)
}

# ---------------------------------------------------------------------------
# Per-draw scores via analytical formula + batched TMB function calls
#
# Score identity (derived from zero average gradient at anchor):
#   g_{ds,j} = (n-p)/2 * [d log(PWRSS_bar)/dphi_j - d log(PWRSS_s)/dphi_j]
#            = (n-p)/2 * [dPWRSS_bar_j/PWRSS_bar - dPWRSS_s_j/PWRSS_s]
#
# Needs 2*len(phi) batched TMB fn() calls (each covers all S draws).
# ---------------------------------------------------------------------------

blmm_scores_batch <- function(obj, phi_bar, PWRSS_s, n, p, eps = 1e-4) {
  S     <- length(PWRSS_s)
  np    <- n - p
  q_phi <- length(phi_bar)
  PWRSS_bar <- mean(PWRSS_s)
  scores    <- matrix(0, nrow = q_phi, ncol = S)

  for (j in seq_len(q_phi)) {
    phi_p <- phi_bar; phi_p[j] <- phi_p[j] + eps
    phi_m <- phi_bar; phi_m[j] <- phi_m[j] - eps

    obj$fn(phi_p); rp <- obj$report()
    obj$fn(phi_m); rm <- obj$report()

    PWRSS_p <- pmax(rp$PWRSS_s, .Machine$double.eps)
    PWRSS_m <- pmax(rm$PWRSS_s, .Machine$double.eps)

    dPWRSS_s   <- (PWRSS_p - PWRSS_m) / (2 * eps)
    dPWRSS_bar <- mean(dPWRSS_s)

    scores[j, ] <- np / 2 * (dPWRSS_bar / PWRSS_bar - dPWRSS_s / PWRSS_s)
  }
  scores
}

# ---------------------------------------------------------------------------
# Newton step: phi_tilde_ds = phi_bar + H_d^{-1} g_ds
# ---------------------------------------------------------------------------

blmm_phi_updates <- function(phi_bar, H_d, scores) {
  cholH        <- chol(H_d)
  H_inv_scores <- backsolve(cholH, forwardsolve(t(cholH), scores))
  phi_bar + H_inv_scores
}

# ---------------------------------------------------------------------------
# Exact conditional fixed effects for one draw given phi_draw
# Uses the same q x q structure as the TMB anchor; reuses obj$report()
# ---------------------------------------------------------------------------

blmm_fixed_effects_draw <- function(obj, phi_draw, X, y) {
  n <- nrow(X); p <- ncol(X)
  phi_draw <- pmax(pmin(phi_draw, 15), -15)

  pieces <- tryCatch({
    obj$fn(phi_draw)
    obj$report()
  }, error = function(e) NULL)

  if (is.null(pieces)) {
    # OLS fallback
    warning("blmm: TMB evaluation failed for a draw; falling back to OLS")
    A     <- crossprod(X)
    cholA <- chol(A)
    beta  <- backsolve(cholA, forwardsolve(t(cholA), crossprod(X, y)))
    resid <- y - X %*% beta
    s2    <- sum(resid^2) / (n - p)
    se    <- sqrt(s2 * diag(chol2inv(cholA)))
    return(list(beta = as.vector(beta), se = as.vector(se),
                df = n - p, sigma2 = s2))
  }

  L   <- pieces$L    # q x q lower Cholesky
  RZX <- pieces$RZX  # q x p
  LX  <- pieces$LX   # p x p lower Cholesky

  # Since Y_d was set to a single column (y), use RZY and Cbeta from report
  # Or recompute for this single y (obj was called with Y = matrix(y, ncol=1))
  # We need to call obj with Y = matrix(y, ncol=1) first.
  # (This is handled by the caller: see blmm() inner loop.)

  LambdaZt <- pieces$LambdaZt   # q x n

  cu    <- forwardsolve(L, LambdaZt %*% y)             # q-vector
  Cu_y  <- crossprod(X, y) - crossprod(RZX, cu)         # p-vector
  Cbeta <- forwardsolve(LX, Cu_y)                       # p-vector
  beta  <- backsolve(t(LX), Cbeta)                      # p-vector

  PWRSS <- max(sum(y^2) - sum(cu^2) - sum(Cbeta^2), 0)
  s2    <- PWRSS / (n - p)
  Minv  <- chol2inv(t(LX))
  se    <- sqrt(pmax(s2 * diag(Minv), 0))

  list(beta = as.vector(beta), se = as.vector(se),
       df = n - p, sigma2 = s2)
}

# ---------------------------------------------------------------------------
# Main engine function
# ---------------------------------------------------------------------------

##' Fast Approximate Mixed-Effects Engine (BLMM)
##'
##' Internal engine called by \code{\link{aldex}} when \code{method = "blmm"}.
##'
##' Replaces S separate \code{lmer()} optimisations per feature with:
##' \enumerate{
##'   \item ONE batched TMB anchor optimisation per feature (exact gradients
##'         via automatic differentiation, q x q Cholesky, not n x n V).
##'   \item S Newton-step corrections (no per-draw optimisation).
##'   \item S exact GLS fixed-effect solves conditional on approximated
##'         covariance.
##' }
##'
##' Supports all lme4 model structures (random intercepts, slopes, nested
##' groupings, correlated random effects) via lme4's reTrms.
##'
##' \strong{The approximation is ONLY in the variance-component step.}
##' Fixed-effect estimates are exact GLS conditional on the approximate
##' covariance.
##'
##' \strong{Known limitations:} degrees of freedom use the conservative
##' N-P approximation; feature-level parallelism not yet implemented.
##'
##' @param logW    numeric array (N x D x S)
##' @param formula lme4 mixed-effects formula
##' @param data    data.frame for formula evaluation
##' @param n.cores currently ignored (reserved for future parallelism)
##' @return same list as \code{sr.mem}: named arrays (P x D x S) for
##'   estimate, std.error, df, p.lower, p.upper; (Pr x D x S) for
##'   random.eff
##' @importFrom lme4 lFormula lmerControl
##' @importFrom TMB MakeADFun compile dynlib
blmm <- function(logW, formula, data, n.cores = 1L) {
  N <- dim(logW)[1]
  D <- dim(logW)[2]
  S <- dim(logW)[3]

  ## Parse formula once (expensive: lFormula + basis computation)
  parsed  <- blmm_parse_formula(formula, data)
  X       <- parsed$X
  basis   <- parsed$basis
  is_log  <- parsed$is_log
  phi_init <- parsed$phi_init
  reTrms  <- parsed$reTrms
  lower   <- reTrms$lower
  re_diag <- parsed$re_diag_idx
  p       <- ncol(X)
  n_obs   <- nrow(X)
  fnames  <- colnames(X)
  Pr      <- length(re_diag) + 1L
  rnames  <- c(parsed$re_varnames, "Residual")

  ## Output arrays
  estimate   <- array(NA_real_, c(p,  D, S), dimnames = list(fnames, NULL, NULL))
  std.error  <- array(NA_real_, c(p,  D, S), dimnames = list(fnames, NULL, NULL))
  df_arr     <- array(NA_real_, c(p,  D, S), dimnames = list(fnames, NULL, NULL))
  p.lower    <- array(NA_real_, c(p,  D, S), dimnames = list(fnames, NULL, NULL))
  p.upper    <- array(NA_real_, c(p,  D, S), dimnames = list(fnames, NULL, NULL))
  random.eff <- array(NA_real_, c(Pr, D, S), dimnames = list(rnames, NULL, NULL))

  ## Feature loop
  for (d in seq_len(D)) {
    Y_d <- logW[, d, , drop = FALSE]
    dim(Y_d) <- c(N, S)

    ## 1. Build TMB object for this feature's response block
    obj_full <- blmm_make_adfun(X, Y_d, basis, is_log, phi_init)

    ## 2. Anchor optimisation (TMB provides exact gradients via AD)
    anchor  <- blmm_fit_anchor(obj_full, phi_init)
    phi_bar <- anchor$phi_bar
    H_d     <- anchor$H_d

    ## 3. Get PWRSS_s at anchor (needed for scores)
    obj_full$fn(phi_bar)
    PWRSS_s_bar <- obj_full$report()$PWRSS_s

    ## 4. Per-draw scores via batched formula (2*len(phi) TMB fn() calls)
    g_mat     <- blmm_scores_batch(obj_full, phi_bar, PWRSS_s_bar, N, p)

    ## 5. Newton-step phi updates
    phi_tilde <- blmm_phi_updates(phi_bar, H_d, g_mat)  # len(phi) x S

    ## 6. Exact conditional fixed effects per draw
    ##    (rebuild obj with single-column Y for per-draw fixed effects)
    for (s in seq_len(S)) {
      y_s    <- Y_d[, s]
      obj_1  <- blmm_make_adfun(X, matrix(y_s, ncol = 1L), basis, is_log,
                                 phi_tilde[, s])
      fe <- blmm_fixed_effects_draw(obj_1, phi_tilde[, s], X, y_s)

      estimate [, d, s] <- fe$beta
      std.error[, d, s] <- fe$se
      df_arr   [, d, s] <- fe$df

      t_stat <- fe$beta / fe$se
      p.lower[, d, s] <- pt(t_stat, df = fe$df, lower.tail = TRUE)
      p.upper[, d, s] <- 1 - p.lower[, d, s]

      theta_draw <- phi_to_theta_blmm(phi_tilde[, s], lower)
      random.eff[seq_along(re_diag), d, s] <-
        fe$sigma2 * theta_draw[re_diag]^2
      random.eff[Pr, d, s] <- fe$sigma2
    }
  }

  list(
    estimate   = estimate,
    std.error  = std.error,
    df         = df_arr,
    p.lower    = p.lower,
    p.upper    = p.upper,
    random.eff = random.eff
  )
}
