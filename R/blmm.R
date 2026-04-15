## BLMM: Batched Linear Mixed-Effects Model engine for ALDEx3.
## Internal, not exported.
##
## This implementation keeps the existing lme4-style profiled formulation:
## one batched profiled REML anchor fit per feature, draw-specific local
## covariance updates around that anchor, and exact GLS fixed-effect solves
## conditional on the updated covariance parameters. The approximation is
## only in the variance-component optimisation step.

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
# Compute basis matrices: dLambdaZt/dtheta_j
# ---------------------------------------------------------------------------
# Lambdat@x = theta[Lind] is linear in theta, so each basis slice is the
# exact map obtained by setting theta = e_j, then forming Lambdat %*% Zt.

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
# Return values follow the same naming convention as sr.mem/lme4:
# group.var for variances, group.var1.var2 for covariance terms, and Residual.

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
# TMB object construction
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Anchor optimisation and curvature handling
# ---------------------------------------------------------------------------
# The local update uses the profiled objective's observed Hessian. Small
# negative eigenvalues from numerical noise are regularised conservatively and
# explicitly, rather than changing the model silently.

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
# Per-draw scores around the anchor
# ---------------------------------------------------------------------------
# This keeps the batched anchor trick intact: each perturbation reuses the same
# shared factorisation logic and evaluates all draws at once.

blmm_scores_batch <- function(obj, phi_bar, PWRSS_s, n, p, eps = 1e-4) {
  if (any(!is.finite(PWRSS_s)) || any(PWRSS_s <= 0)) {
    stop("blmm: invalid profiled residual sum of squares at the anchor")
  }

  S <- length(PWRSS_s)
  np <- n - p
  q_phi <- length(phi_bar)
  PWRSS_bar <- mean(PWRSS_s)
  scores <- matrix(0, nrow = q_phi, ncol = S)

  if (!is.finite(PWRSS_bar) || PWRSS_bar <= 0) {
    stop("blmm: invalid batched profiled residual sum of squares at the anchor")
  }

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

    dPWRSS_s <- (PWRSS_p - PWRSS_m) / (2 * step)
    dPWRSS_bar <- mean(dPWRSS_s)
    scores[j, ] <- np / 2 * (dPWRSS_bar / PWRSS_bar - dPWRSS_s / PWRSS_s)
  }

  scores
}

# ---------------------------------------------------------------------------
# Newton update around the anchor
# ---------------------------------------------------------------------------

blmm_phi_updates <- function(phi_bar, H_d, scores) {
  cholH <- tryCatch(chol(H_d), error = function(e) e)
  if (inherits(cholH, "error")) {
    stop("blmm: curvature matrix is not numerically positive definite")
  }

  H_inv_scores <- backsolve(cholH, forwardsolve(t(cholH), scores))
  phi_bar + H_inv_scores
}

# ---------------------------------------------------------------------------
# Exact conditional fixed effects given covariance parameters
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Exact lme4 fallback for features the approximation cannot handle cleanly
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Main engine
# ---------------------------------------------------------------------------

##' Fast Approximate Mixed-Effects Engine (BLMM)
##'
##' Internal engine called by \code{\link{aldex}} when \code{method = "blmm"}.
##'
##' The implementation uses an \code{lme4}-style profiled Gaussian LMM
##' formulation for efficiency and compatibility with \code{lFormula()}.
##' Conceptually, ALDEx3 fits one batched anchor REML problem per feature,
##' takes draw-specific local covariance updates around that anchor, then
##' performs exact conditional GLS fixed-effect solves for each draw.
##'
##' \strong{The approximation is ONLY in the variance-component step.}
##' Fixed-effect estimates and standard errors remain exact conditional on the
##' updated covariance parameters.
##'
##' If the approximate path cannot be evaluated cleanly for a feature, BLMM
##' warns explicitly and falls back to the exact \code{lme4} engine for that
##' feature rather than silently changing the model.
##'
##' @param logW numeric array (N x D x S)
##' @param formula lme4 mixed-effects formula
##' @param data data.frame for formula evaluation
##' @param n.cores currently ignored; reserved for future feature-level
##'   parallelism
##' @return The same list layout as \code{sr.mem}: named arrays (P x D x S) for
##'   estimate, std.error, df, p.lower, p.upper, and (Pr x D x S) for
##'   random.eff.
##' @importFrom lme4 lFormula lmerControl
##' @importFrom TMB MakeADFun
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

  estimate <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  std.error <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  df_arr <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  p.lower <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  p.upper <- array(NA_real_, c(p, D, S), dimnames = list(fnames, NULL, NULL))
  random.eff <- array(NA_real_, c(Pr, D, S), dimnames = list(rnames, NULL, NULL))

  for (d in seq_len(D)) {
    Y_d <- logW[, d, , drop = FALSE]
    dim(Y_d) <- c(N, S)

    feature_fit <- tryCatch({
      obj_full <- blmm_make_adfun(X, Y_d, basis, is_log, phi_init)
      anchor <- blmm_fit_anchor(obj_full, phi_init)

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
      feature_random <- matrix(NA_real_, nrow = Pr, ncol = S)

      for (s in seq_len(S)) {
        obj_full$fn(phi_tilde[, s])
        pieces <- obj_full$report()
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
        estimate = array(feature_estimate, c(p, 1L, S)),
        std.error = array(feature_std.error, c(p, 1L, S)),
        df = array(feature_df, c(p, 1L, S)),
        p.lower = array(feature_p.lower, c(p, 1L, S)),
        p.upper = array(feature_p.upper, c(p, 1L, S)),
        random.eff = array(feature_random, c(Pr, 1L, S))
      )
    }, error = function(e) {
      warning(
        sprintf(
          paste0(
            "blmm: feature %d fell back to exact lme4 because the approximate ",
            "path failed: %s"
          ),
          d, conditionMessage(e)
        ),
        call. = FALSE
      )
      blmm_exact_feature(Y_d, formula, data)
    })

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
