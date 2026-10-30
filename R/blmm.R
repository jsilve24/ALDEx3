## BLMM: Batched Linear Mixed-Effects Model engine for ALDEx3.
##
## Internal implementation. Not exported. One batched profiled REML anchor fit
## per feature, draw-specific local covariance updates around that anchor, and
## exact GLS fixed-effect solves conditional on the updated covariance
## parameters. Conditional GLS does not imply exact finite-sample inference.

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------------------------------------------------------------------------
# Parameterisation helpers
# ---------------------------------------------------------------------------

#' Map unconstrained phi to lme4 variance-component parameters theta.
#'
#' @details
#' lme4 constrains variance-component parameters theta >= 0 (lower > -Inf).
#' To allow unconstrained optimisation those entries are log-transformed:
#'   \code{phi[j] = log(theta[j])}  for constrained parameters
#'   \code{phi[j] = theta[j]}       for unconstrained parameters (e.g. correlations)
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
blmm_random_effect_vector <- function(theta, sigma2, reTrms,
                                      metadata = blmm_random_effect_metadata(reTrms)) {
  values <- numeric(length(metadata$names))
  for (block in metadata$blocks) {
    L <- matrix(0, block$size, block$size)
    L[block$lower] <- theta[block$theta]
    Sigma <- sigma2 * tcrossprod(L)
    values[block$output] <- Sigma[block$covariance]
  }
  values[length(values)] <- sigma2
  stats::setNames(values, metadata$names)
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

# ---------------------------------------------------------------------------
# Output arrays
# ---------------------------------------------------------------------------

# Per-coefficient inference fields. Contrast arrays repeat them with a
# "contrast." prefix; random-effect variances are stored in "random.eff".
.blmm_fields <- c("estimate", "std.error", "df", "p.lower", "p.upper")

blmm_output_fields <- function(contrasts = NULL) {
  c(.blmm_fields, "random.eff",
    if (!is.null(contrasts)) paste0("contrast.", .blmm_fields))
}

# NA-filled (rows x D x S) output arrays in blmm_output_fields() order. Row
# names are attached only when supplied; contrast rows use rownames(contrasts).
blmm_output_arrays <- function(p, random_rows, D, S, contrasts = NULL,
                               names = NULL, random_names = NULL) {
  allocate <- function(rows, row_names) {
    array(NA_real_, c(rows, D, S),
          dimnames = if (!is.null(row_names)) list(row_names, NULL, NULL))
  }
  out <- stats::setNames(lapply(.blmm_fields, function(field)
    allocate(p, names)), .blmm_fields)
  out$random.eff <- allocate(random_rows, random_names)
  if (!is.null(contrasts)) for (field in .blmm_fields) {
    out[[paste0("contrast.", field)]] <- allocate(nrow(contrasts),
                                                  rownames(contrasts))
  }
  out
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
#' @return list with X, reTrms, cached design (basis), phi_init, is_log, lower, random_template
#' @noRd
blmm_parse_formula <- function(formula, data,
                               backend = getOption("ALDEx3.blmm.backend", "auto"),
                               draws = 1L, workers = 1L, features = 1L) {
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
    basis = blmm_design(reTrms, parsed$X, backend, draws, workers, features),
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
#' @param basis cached design from \code{blmm_design()}
#' @param is_log integer vector flagging which phi entries are log-transformed
#' @param phi_init starting values for the optimisation
#' @param reusable whether response cross-products can change without retaping
#' @param object existing reusable objective for the same design and dimensions
#' @return TMB autodiff object with \code{$fn}, \code{$gr}, \code{$he},
#'   and \code{$report()} methods
#' @noRd
blmm_make_adfun <- function(X, Y_d, basis, is_log, phi_init,
                            response = blmm_response_cache(basis, X, Y_d),
                            reusable = FALSE, object = NULL) {
  if (!is.null(object)) {
    blmm_update_response(object, response)
    return(object)
  }
  TMB::MakeADFun(
    data = list(
      Lambdat = basis$Lambdat, Lind = as.integer(basis$Lind - 1L),
      G = basis$G, ZX = basis$ZX, XX = basis$XX,
      ZY = response$ZY, XY = response$XY, yy = response$yy,
      update_response = as.integer(reusable),
      n = nrow(X), is_log = is_log,
      sparse_backend = as.integer(basis$backend == "sparse")
    ),
    parameters = list(phi = phi_init),
    DLL = "ALDEx3",
    silent = TRUE
  )
}

# Dynamic response data are local to one feature/worker. Validate the complete
# replacement before assigning it: TMB requires unchanged dimensions and types.
blmm_update_response <- function(objective, response) {
  data <- objective$env$data
  if (!isTRUE(data$update_response == 1L))
    stop("blmm: objective does not support response updates")
  for (field in c("ZY", "XY", "yy")) {
    value <- response[[field]]
    if (is.null(value) || !identical(dim(value), dim(data[[field]])) ||
        length(value) != length(data[[field]]) ||
        !identical(typeof(value), typeof(data[[field]])) ||
        any(!is.finite(value)))
      stop("blmm: response update must preserve dimensions, storage and finite values")
    data[[field]] <- value
  }
  objective$env$data <- data
  invisible(objective)
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

# ---------------------------------------------------------------------------
# Anchor optimisation and curvature handling
# ---------------------------------------------------------------------------

#' Symmetrise and validate the observed Hessian at the anchor.
#'
#' @details
#' A non-finite, non-positive-definite, or numerically ill-conditioned Hessian
#' is not safe to invert for the shared Newton update. Such matrices are
#' rejected rather than regularised so the caller can fall back to exact
#' \code{lme4}. The relative eigenvalue criterion is scale invariant; for a
#' one-dimensional positive Hessian its value is necessarily one.
#'
#' @param H square numeric matrix: raw observed Hessian from \code{obj$he()}
#' @param tol minimum acceptable ratio of smallest to largest eigenvalue
#' @return the symmetrised Hessian, unchanged apart from symmetrisation
#' @noRd
blmm_validate_H <- function(H, tol = 1e-6) {
  if (!is.matrix(H) || nrow(H) == 0L || nrow(H) != ncol(H)) {
    stop("blmm: observed Hessian must be a non-empty square matrix")
  }

  H <- (H + t(H)) / 2
  if (any(!is.finite(H))) {
    stop("blmm: observed Hessian contains non-finite values")
  }

  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min_ev <- min(ev)
  max_ev <- max(ev)
  if (min_ev <= 0) {
    stop("blmm: observed Hessian is not positive definite")
  }

  relative_min <- min_ev / max_ev
  if (!is.finite(relative_min) || relative_min < tol) {
    stop(
      "blmm: observed Hessian is numerically ill-conditioned ",
      sprintf("(minimum/maximum eigenvalue ratio %.3g < %.3g)",
              relative_min, tol)
    )
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
#' @return list with \code{phi_bar} (optimal phi) and \code{H_d} (validated
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

  list(phi_bar = fit$par, H_d = blmm_validate_H(H_raw))
}

# ---------------------------------------------------------------------------
# Per-draw score corrections around the anchor
# ---------------------------------------------------------------------------

#' Compute per-draw score corrections at the anchor.
#'
#' @details
#' The anchor \code{phi_bar} minimises the \emph{average} profiled objective
#' over all S draws. For draw s, the draw-specific objective f_s(phi) differs
#' from the average only through its PWRSS_s term. The negative gradient of
#' f_s at phi_bar equals the correction used below (up to
#' the anchor optimisation tolerance): the direction in which that draw wants
#' to move away from the shared anchor.
#'
#' Scores use native analytic derivatives of log(PWRSS_s). The
#' ldL2/ldRX2 log-determinant terms are identical for all draws at any fixed
#' phi, so their derivative cancels out of the per-draw correction:
#' \deqn{score_s[j] = \frac{n-p}{2}
#'   \left(\overline{\frac{\partial \log PWRSS}{\partial \phi_j}}
#'         - \frac{\partial \log PWRSS_s}{\partial \phi_j}\right)}
#'
#' @param obj TMB autodiff object; reports are evaluated with explicit parameters
#' @param phi_bar numeric vector: the anchor (optimal phi for the average objective)
#' @param PWRSS_s numeric vector of length S: per-draw penalised residual SS at
#'   the anchor, from \code{obj$report()$PWRSS_s}
#' @param n number of observations
#' @param p number of fixed-effect columns
#' @param contrasts Matrix with contrast vectors in columns; defaults to identity.
#' @param inference Retain uncentered residual and contrast-variance derivatives.
#' @param report Native report at phi_bar; reuse an existing report when available.
#' @return matrix of dimension \code{(length(phi_bar), S)}: column s is the
#'   score vector for draw s; with inference, a list of scores, u, and b.
#' @noRd
blmm_scores_batch <- function(obj, phi_bar, PWRSS_s, n, p,
                              contrasts = NULL, inference = FALSE,
                              report = obj$report(phi_bar)) {
  if (any(!is.finite(PWRSS_s)) || any(PWRSS_s <= 0)) {
    stop("blmm: invalid profiled residual sum of squares at the anchor")
  }

  S <- length(PWRSS_s)
  np <- n - p
  q_phi <- length(phi_bar)
  u <- report$dlogPWRSS
  if (!is.matrix(u) || !identical(dim(u), c(q_phi, S)) || any(!is.finite(u)))
    stop("blmm: invalid native residual-variance sensitivity")
  scores <- np / 2 * (rowMeans(u) - u)
  if (!inference) return(scores)

  if (is.null(contrasts)) contrasts <- diag(p)
  if (!is.matrix(report$dMinv) || !identical(dim(report$dMinv), c(q_phi, p * p)) ||
      any(!is.finite(report$dMinv)))
    stop("blmm: invalid native contrast-variance sensitivity")
  variance <- colSums(forwardsolve(report$LX, contrasts)^2)
  if (any(!is.finite(variance)) || any(variance <= 0))
    stop("blmm: invalid anchor contrast variance")
  # Contract d M^-1 with c c' for each contrast, retaining singleton dimensions.
  index <- seq_len(p)
  weights <- contrasts[rep(index, p), , drop = FALSE] *
    contrasts[rep(index, each = p), , drop = FALSE]
  b <- sweep(report$dMinv %*% weights, 2L, variance, "/")
  list(scores = scores, u = u, b = b)
}

# Shared-anchor Satterthwaite df; see the derivation at the end of
# vignettes/ALDEx3-mixed-effects.Rmd.
blmm_anchor_df <- function(u, b, cholH, d) {
  h <- b + as.vector(u)
  if (!is.finite(d) || d <= 0 || any(!is.finite(h)))
    stop("blmm: invalid anchor inference sensitivity")
  whitened <- forwardsolve(t(cholH), h)
  quadratic <- colSums(whitened^2)
  df <- 2 / (2 / d + quadratic)
  if (any(!is.finite(quadratic)) || any(!is.finite(df)) || any(df <= 0))
    stop("blmm: invalid anchor inference degrees of freedom")
  df
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
#' @param H_d positive-definite matrix: validated Hessian at the anchor
#' @param scores matrix (length(phi_bar) x S): per-draw score vectors from
#'   \code{blmm_scores_batch()}
#' @param cholH Reusable upper Cholesky factor of H_d.
#' @return matrix of the same dimension as \code{scores}: column s is the
#'   draw-specific updated phi
#' @noRd
blmm_phi_updates <- function(phi_bar, H_d, scores, cholH = chol(H_d)) {
  cholH <- tryCatch(cholH, error = function(e) e)
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
#' @param pieces factor/solve interface from \code{blmm_common_pieces()}
#' @param X fixed-effects design matrix (n x p)
#' @param y numeric vector of length n: response for this draw
#' @param response optional cached cross-products for this draw
#' @return list with beta, se, df, sigma2
#' @noRd
blmm_fixed_effects_draw <- function(pieces, X, y, response = NULL) {
  n <- nrow(X)
  p <- ncol(X)

  LX <- pieces$LX
  if (is.null(response)) response <- blmm_response_cache(pieces$design, X, y)
  by <- as.matrix(pieces$T %*% response$ZY)
  solved_y <- blmm_factor_solve(pieces$factor, by, pieces$backend)
  Cu_y <- response$XY - crossprod(pieces$BX, solved_y)
  Cbeta <- forwardsolve(LX, Cu_y)
  beta <- backsolve(t(LX), Cbeta)

  PWRSS <- as.numeric(response$yy - sum(by * solved_y) - sum(Cbeta^2))
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

#' Profiled REML objective from one draw's conditional solve.
#'
#' @details
#' Reconstructs the same profiled objective used by the TMB anchor from the
#' Cholesky factors and residual variance already computed for a proposed
#' draw-specific covariance update. This avoids a second TMB evaluation.
#'
#' @param pieces list from \code{blmm_common_pieces()}
#' @param sigma2 positive profiled residual variance from
#'   \code{blmm_fixed_effects_draw()}
#' @param n number of observations
#' @param p number of fixed-effect columns
#' @return scalar profiled REML objective, up to the same additive constant as
#'   \code{blmm_profiled_objectives()}
#' @noRd
blmm_profiled_objective_from_pieces <- function(pieces, sigma2, n, p) {
  diag_LX <- diag(pieces$LX)
  if (!is.finite(sigma2) || sigma2 <= 0 ||
      !is.finite(pieces$ldL2) ||
      any(!is.finite(diag_LX)) || any(diag_LX <= 0)) {
    stop("blmm: invalid quantities in the proposed profiled REML objective")
  }

  pieces$ldL2 / 2 + sum(log(diag_LX)) + (n - p) / 2 * log(sigma2)
}

#' Recompute the random-effect factor and dense fixed-effect Schur complement.
#'
#' @details The cached design is independent of covariance parameters and Y.
#' Sparse factors use the structural symbolic factor with fresh numeric values.
#' Permutations are handled by the factor solve, never by callers.
#' @param phi draw-specific unconstrained covariance parameters
#' @param X fixed-effects design
#' @param basis cached design from \code{blmm_design()}
#' @param lower bounds for the lme4 covariance parameters
#' @return factor/solve interface, BX, LX, and log determinant ldL2
#' @noRd
blmm_common_pieces <- function(phi, X, basis, lower) {
  theta <- phi_to_theta_blmm(phi, lower)
  if (any(!is.finite(theta)))
    stop("blmm: proposed covariance parameters are non-finite")
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

# ---------------------------------------------------------------------------
# Exact lme4 fallback
# ---------------------------------------------------------------------------

#' Exact lme4 fallback fit for one feature and one or more draws.
#'
#' @details
#' Called when the approximate BLMM path fails for a feature or selected draws.
#' Wraps
#' \code{sr.mem()} with a singleton feature dimension so the output layout
#' matches what \code{blmm_fit_feature()} expects.
#'
#' @param Y_d numeric matrix (n x S): log-composition draws for this feature
#' @param formula lme4 mixed-effects formula
#' @param data data.frame for formula evaluation
#' @param contrasts optional validated contrast matrix with one row per
#'   contrast and columns matching the fixed-effect design.
#' @return list with estimate, std.error, df, p.lower, p.upper, random.eff
#'   arrays — each retaining a singleton feature dimension (D = 1), plus
#'   contrast arrays when contrasts are supplied.
#' @noRd
blmm_exact_feature <- function(Y_d, formula, data, contrasts = NULL) {
  exact <- sr.mem(
    logW = array(Y_d, dim = c(nrow(Y_d), 1L, ncol(Y_d))),
    formula = formula,
    data = data,
    n.cores = 1L,
    method = "lme4",
    mem.args = list(),
    contrasts = contrasts
  )

  fields <- blmm_output_fields(contrasts)
  stats::setNames(lapply(fields, function(field)
    exact[[field]][, 1, , drop = FALSE]), fields)
}

#' Approximate BLMM fit for one feature across all S MC draws.
#'
#' @details
#' Runs the full approximate pipeline: (1) optimise the average-over-draws
#' profiled REML objective to get the shared anchor and its Hessian;
#' (2) compute per-draw score corrections via analytic derivatives and
#' apply a Newton step to get draw-specific covariance parameters phi_tilde;
#' (3) for each draw, recompute the Cholesky pieces from \code{phi_tilde[, s]} and
#' run an exact GLS solve for fixed effects.
#'
#' Failure of the shared anchor optimisation, or an unusable anchor Hessian,
#' causes the entire feature to be refit with \code{lmerTest}. Once an anchor
#' is available, each proposed draw-specific update must leave that draw's
#' profiled REML objective unchanged within numerical tolerance or improve it.
#' Rejected draws are re-fit exactly and spliced into the feature result;
#' accepted draws retain the fast BLMM approximation.
#'
#' @param d integer feature index (column of logW)
#' @param logW numeric array (N x D x S)
#' @param X fixed-effects design matrix
#' @param basis cached design
#' @param is_log integer vector for phi/theta transformation
#' @param phi_init starting values
#' @param lower lower bounds for theta
#' @param reTrms random-effects terms list
#' @param formula lme4 formula
#' @param data data.frame
#' @param forced_fallback integer vector of feature indices to force onto the
#'   lmerTest refit path; a testing hook, not used in normal operation
#' @param inference internal selector: \code{"adaptive"} (default),
#'   \code{"anchor"}, or deprecated \code{"residual"}.
#' @param contrasts optional validated contrast matrix with one row per
#'   contrast and columns matching \code{X}.
#' @return list with \code{fit} (arrays), fallback status and counts, and a
#'   diagnostic \code{fallback_message}
#' @noRd
blmm_fit_feature <- function(d, logW, X, basis, is_log, phi_init, lower,
                             reTrms, formula, data, forced_fallback = integer(0),
                             inference = c("adaptive", "anchor", "residual"),
                             contrasts = NULL, strict = blmm_strict(),
                             .design_context = NULL) {
  inference <- match.arg(inference)
  N <- nrow(X)
  p <- ncol(X)
  S <- dim(logW)[3]
  if (!blmm_valid_design_context(.design_context, X, contrasts, reTrms)) {
    .design_context <- blmm_inference_context(X, contrasts, reTrms)
    contrasts <- .design_context$contrasts
  }
  if (!is.null(contrasts) && .design_context$qr$rank < p) {
    stop("blmm: contrasts require a full-rank fixed-effect design")
  }
  Y_d <- logW[, d, , drop = FALSE]
  dim(Y_d) <- c(N, S)

  feature_random_rows <- length(.design_context$random$names)
  if (length(forced_fallback) > 0L && d %in% forced_fallback) {
    started <- proc.time()[[3L]]
    exact <- blmm_exact_feature(Y_d, formula, data, contrasts = contrasts)
    return(list(
      fit = exact,
      rescue_seconds = proc.time()[[3L]] - started,
      fallback = TRUE,
      fallback_message = "forced approximate failure for testing",
      fallback_draws = S,
      rejected_draws = seq_len(S),
      full_fallback = TRUE
    ))
  }

  if (inference == "adaptive") {
    return(blmm_fit_feature_adaptive(
      d, Y_d, X, basis, lower, reTrms, formula, data,
      contrasts, strict, .design_context
    ))
  }

  rescue_seconds <- 0
  result <- tryCatch({
    # Step 1: anchor — optimise the average-over-draws objective.
    response <- blmm_response_cache(basis, X, Y_d)
    obj_full <- blmm_make_adfun(X, Y_d, basis, is_log, phi_init, response)
    anchor <- blmm_fit_anchor(obj_full, phi_init)

    # Step 2: score corrections — differentiate per-draw PWRSS w.r.t. phi at
    # the anchor, then Newton-step to draw-specific covariance parameters.
    anchor_report <- obj_full$report(anchor$phi_bar)
    anchor_objectives <- blmm_profiled_objectives(anchor_report, N, p)
    if (any(!is.finite(anchor_objectives))) {
      stop("blmm: non-finite profiled REML objective at the anchor")
    }
    inference_vectors <- .design_context$vectors
    g_mat <- blmm_scores_batch(
      obj_full, anchor$phi_bar, anchor_report$PWRSS_s, N, p,
      contrasts = inference_vectors, inference = inference == "anchor",
      report = anchor_report
    )
    cholH <- chol(anchor$H_d)
    sensitivity <- if (inference == "anchor") g_mat else NULL
    if (inference == "anchor") g_mat <- g_mat$scores
    phi_tilde <- blmm_phi_updates(anchor$phi_bar, anchor$H_d, g_mat, cholH)

    fit <- blmm_output_arrays(p, feature_random_rows, 1L, S, contrasts)
    K <- if (is.null(contrasts)) 0L else nrow(contrasts)
    rejected <- rep(FALSE, S)
    rejection_messages <- rep("", S)

    # Step 3: exact conditional GLS solve per draw. Accept the local covariance
    # update only if it does not materially worsen that draw's profiled REML
    # objective. Anchor inference is validated inside this same draw-specific
    # rescue boundary; residual df remains available for reproducible scripts.
    for (s in seq_len(S)) {
      draw_result <- tryCatch({
        pieces <- blmm_common_pieces(phi_tilde[, s], X, basis, lower)
        response_s <- list(ZY = response$ZY[, s, drop = FALSE],
                           XY = response$XY[, s, drop = FALSE], yy = response$yy[s])
        fe <- blmm_fixed_effects_draw(pieces, X, Y_d[, s], response_s)
        proposed_objective <- blmm_profiled_objective_from_pieces(
          pieces, fe$sigma2, N, p
        )
        objective_delta <- proposed_objective - anchor_objectives[s]
        objective_tol <- sqrt(.Machine$double.eps) *
          (1 + abs(anchor_objectives[s]))
        if (!is.finite(objective_delta) || objective_delta > objective_tol) {
          stop(
            sprintf(
              paste0(
                "proposed covariance update worsened the profiled REML ",
                "objective by %.6g (tolerance %.6g)"
              ),
              objective_delta,
              objective_tol
            )
          )
        }

        if (inference == "anchor") {
          coefficient_df <- blmm_anchor_df(
            sensitivity$u[, s], sensitivity$b[, seq_len(p), drop = FALSE],
            cholH, N - p
          )
          fe$df <- coefficient_df
        }

        contrast_fit <- NULL
        if (!is.null(contrasts)) {
          # LX is lower triangular with LX %*% t(LX) = M. The full
          # covariance is therefore sigma2 * M^-1, including off-diagonals.
          contrast_variance <- fe$sigma2 * colSums(
            forwardsolve(pieces$LX, t(contrasts))^2
          )
          if (any(!is.finite(contrast_variance)) || any(contrast_variance <= 0)) {
            stop("proposed covariance update produced invalid contrast variance")
          }
          contrast_df <- rep(N - p, nrow(contrasts))
          if (inference == "anchor") {
            contrast_df <- blmm_anchor_df(
              sensitivity$u[, s], sensitivity$b[, p + seq_len(nrow(contrasts)), drop = FALSE],
              cholH, N - p
            )
          }
          contrast_fit <- list(
            estimate = as.vector(contrasts %*% fe$beta),
            std.error = sqrt(contrast_variance),
            df = contrast_df
          )
          if (any(!is.finite(contrast_fit$estimate)) ||
              any(!is.finite(contrast_fit$std.error)) ||
              any(!is.finite(contrast_fit$df)) ||
              any(contrast_fit$std.error <= 0) || any(contrast_fit$df <= 0)) {
            stop("proposed covariance update produced invalid contrast inference")
          }
        }

        theta_draw <- phi_to_theta_blmm(phi_tilde[, s], lower)
        random_draw <- blmm_random_effect_vector(
          theta = theta_draw,
          sigma2 = fe$sigma2,
          reTrms = reTrms,
          metadata = .design_context$random
        )
        if (any(!is.finite(fe$beta)) || any(!is.finite(fe$se)) || any(fe$se <= 0) ||
            any(!is.finite(random_draw))) {
          stop("proposed covariance update produced non-finite results")
        }

        list(fe = fe, contrast = contrast_fit, random = random_draw)
      }, error = function(e) e)

      if (inherits(draw_result, "error")) {
        if (strict) stop(draw_result)
        rejected[s] <- TRUE
        rejection_messages[s] <- conditionMessage(draw_result)
        next
      }

      fe <- draw_result$fe
      contrast_fit <- draw_result$contrast
      values <- list(
        estimate = c(fe$beta, contrast_fit$estimate),
        std.error = c(fe$se, contrast_fit$std.error),
        df = c(rep_len(fe$df, p), contrast_fit$df)
      )
      t_stat <- values$estimate / values$std.error
      values$p.lower <- pt(t_stat, df = values$df, lower.tail = TRUE)
      values$p.upper <- pt(t_stat, df = values$df, lower.tail = FALSE)
      for (field in .blmm_fields) {
        fit[[field]][, 1, s] <- values[[field]][seq_len(p)]
        if (K) fit[[paste0("contrast.", field)]][, 1, s] <-
          values[[field]][p + seq_len(K)]
      }
      fit$random.eff[, 1, s] <- draw_result$random
    }

    rejected_idx <- which(rejected)

    fallback_message <- ""
    if (length(rejected_idx) > 0L) {
      example_idx <- rejected_idx[seq_len(min(3L, length(rejected_idx)))]
      fallback_message <- paste(
        sprintf("draw %d: %s", example_idx, rejection_messages[example_idx]),
        collapse = "; "
      )
      if (length(rejected_idx) > length(example_idx)) {
        fallback_message <- paste0(fallback_message, " ...")
      }
    }

    list(
      fit = fit,
      fallback = length(rejected_idx) > 0L,
      fallback_message = fallback_message,
      fallback_draws = length(rejected_idx),
      rejected_draws = rejected_idx,
      full_fallback = FALSE
    )
  }, error = function(e) {
    if (strict) stop(e)
    started <- proc.time()[[3L]]
    exact <- blmm_exact_feature(Y_d, formula, data, contrasts = contrasts)
    rescue_seconds <<- proc.time()[[3L]] - started
    list(
      fit = exact,
      fallback = TRUE,
      fallback_message = conditionMessage(e),
      fallback_draws = S,
      rejected_draws = seq_len(S),
      full_fallback = TRUE
    )
  })

  # Exact rescue is outside the approximation error handler: an exact-fit
  # failure must not discard accepted draws and retry the entire feature.
  if (!result$full_fallback && length(result$rejected_draws) > 0L) {
    rejected_idx <- result$rejected_draws
    started <- proc.time()[[3L]]
    exact_rejected <- tryCatch(
      blmm_exact_feature(Y_d[, rejected_idx, drop = FALSE], formula, data,
                         contrasts = contrasts),
      error = function(e) stop(
        "blmm: exact refit failed for feature ", d, " (draws ",
        paste(rejected_idx, collapse = ", "), "): ", conditionMessage(e),
        call. = FALSE
      )
    )
    rescue_seconds <- proc.time()[[3L]] - started
    for (field in names(exact_rejected)) {
      result$fit[[field]][, , rejected_idx] <- exact_rejected[[field]]
    }
  }
  result$rescue_seconds <- rescue_seconds
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
#' Uses an \code{lme4}-style profiled Gaussian LMM formulation. Public
#' \code{aldex(method = "blmm")} selects adaptive Satterthwaite inference by
#' default. Responses are centered in the fixed-effect design, one covariance
#' anchor is fitted per feature and batch of draws, and each draw is refined
#' from that anchor using at most eight safeguarded steps. Removed coefficient
#' shifts are restored after mixed-model fitting. The model and its random
#' effects are unchanged by centering.
#'
#' Adaptive inference uses the accepted draw's own observed curvature and
#' sensitivities. For a coefficient or contrast vector c, the first-order df are
#' \deqn{2 / \{2/d + (u_s+b_c)^T H_s^{-1}(u_s+b_c)\}.}
#' Here \code{d = N - p}, with N observations and p fixed-effect coefficients;
#' H_s is the Hessian of the draw's profiled REML negative log likelihood at
#' its accepted covariance parameters. At that same point, u_s is the gradient
#' of log PWRSS and b_c is the gradient of log(c' M^-1 c), where
#' M = X' R^-1 X and R is marginal covariance relative to residual variance.
#' Gradients use the covariance-shape parameters. Supported singular covariance
#' fits have additional first-order boundary checks, described in
#' \code{vignette("ALDEx3-mixed-effects", package = "ALDEx3")}.
#'
#' Conditional GLS solves give each draw's coefficients and SEs. Their numerical
#' validation and Satterthwaite approximation do not establish finite-sample
#' calibration. Draws with unresolved fitting or inference checks are fitted
#' fully with \code{lmerTest}; for the shared-anchor selector, an unavailable
#' anchor causes all draws for that feature to be fitted fully. A warning reports affected features/draws and
#' representative reasons. Replacement-fit warnings remain visible; invalid
#' numerical contrast inference retains valid fitted quantities with a warning
#' and \code{NA} contrast p-values. Fitting errors stop the analysis, and failed
#' draw refits are not retried as whole-feature fallbacks.
#'
#' This unexported function defaults to \code{inference = "adaptive"}, matching
#' \code{aldex()}, which always passes its resolved selector explicitly.
#' Shared-anchor inference uses one covariance proposal per draw and evaluates
#' the df formula with anchor curvature and sensitivities. It remains available
#' publicly as \code{test = "t.satterthwaite"}. Deprecated public
#' \code{test = "t.residual"} uses N-p df for accepted draws. Both alternatives
#' still use \code{lmerTest} Satterthwaite inference for rescued draws. Native
#' \code{sr.mem()} fits use their own lmerTest or nlme inference directly.
#'
#' @param logW numeric array (N x D x S)
#' @param formula lme4 mixed-effects formula
#' @param data data.frame for formula evaluation
#' @param inference Internal inference selector. Public \code{aldex()} defaults
#'   to adaptive inference: it centers each draw, refines covariance parameters
#'   from the shared anchor with up to eight safeguarded steps, and uses local
#'   Satterthwaite inference for coefficients and contrasts. Unresolved draws
#'   are refitted with \code{lmerTest}. \code{"anchor"} selects shared-anchor
#'   Satterthwaite inference (public \code{test = "t.satterthwaite"}). Public deprecated
#'   \code{test = "t.residual"} uses residual degrees of freedom for accepted
#'   draws and Satterthwaite degrees of freedom for refitted draws.
#' @param contrasts optional validated contrast matrix with one row per
#'   contrast and columns matching the fixed-effect design.
#' @param n.cores number of worker processes for feature-level parallelism.
#'   If \code{n.cores > 1}, features are distributed across workers; the
#'   per-feature draw loop remains serial.
#' @return named list of arrays (P x D x S): estimate, std.error, df,
#'   p.lower, p.upper, and (Pr x D x S) random.eff; when contrasts are
#'   supplied, flat K x D x S contrast arrays use the same layout as
#'   \code{sr.mem}.
#' @importFrom lme4 lFormula lmerControl
#' @importFrom TMB MakeADFun
#' @useDynLib ALDEx3, .registration = TRUE
#' @keywords internal
blmm <- function(logW, formula, data, n.cores = 1L,
                 inference = c("adaptive", "anchor", "residual"), contrasts = NULL) {
  inference <- match.arg(inference)
  N <- dim(logW)[1]
  D <- dim(logW)[2]
  S <- dim(logW)[3]

  parsed <- blmm_parse_formula(formula, data, draws = S,
                               workers = max(1L, min(as.integer(n.cores), D)),
                               features = D)
  diagnostic <- getOption("ALDEx3.blmm.diagnostics")
  if (is.function(diagnostic)) diagnostic(parsed$basis$selection)
  strict <- blmm_strict()
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

  .design_context <- blmm_inference_context(X, contrasts, reTrms)
  contrasts <- .design_context$contrasts
  if (!is.null(contrasts) && .design_context$qr$rank < p) {
    stop("blmm: contrasts require a full-rank fixed-effect design")
  }

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
      forced_fallback = forced_fallback,
      inference = inference,
      contrasts = contrasts,
      strict = strict,
      .design_context = .design_context
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
      # Socket workers start fresh R sessions and do not inherit options.
      # Collect diagnostics there, but keep any callback in the parent.
      diagnostic_hook <- getOption("ALDEx3.blmm.adaptive.diagnostics")
      parallel::clusterCall(cl, options, list(
        ALDEx3.blmm.adaptive.diagnostics =
          isTRUE(diagnostic_hook) || is.function(diagnostic_hook)
      ))
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

  # Fork workers return try-error values instead of throwing in the parent.
  # Preserve the original failure before result aggregation can obscure it.
  failed_workers <- which(vapply(feature_results, inherits, logical(1), "try-error"))
  if (length(failed_workers)) {
    failed <- feature_results[[failed_workers[1L]]]
    condition <- attr(failed, "condition")
    if (inherits(condition, "condition")) stop(condition)
    stop(as.character(failed), call. = FALSE)
  }
  if (any(vapply(feature_results, is.null, logical(1))))
    stop("blmm: a feature worker returned no result")

  # Failed local updates and whole-feature failures have already been re-fit by
  # blmm_exact_feature(). Emit one consolidated warning with up to 3 example
  # messages so the caller knows what used the exact path.
  fallback_idx <- which(vapply(feature_results, `[[`, logical(1), "fallback"))
  if (length(fallback_idx) > 0L) {
    full_idx <- which(vapply(
      feature_results, `[[`, logical(1), "full_fallback"
    ))
    partial_idx <- setdiff(fallback_idx, full_idx)
    partial_draws <- sum(vapply(
      feature_results[partial_idx], `[[`, integer(1), "fallback_draws"
    ))

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
    warning_text <- if (length(partial_idx) == 0L) {
      sprintf(
        paste0(
          "blmm: %d feature(s) were refit with lmerTest for all draws because ",
          "the approximate path failed."
        ),
        length(full_idx)
      )
    } else if (length(full_idx) == 0L) {
      sprintf(
        paste0(
          "blmm: %d draw(s) across %d feature(s) were refit with lmerTest ",
          "because their local covariance updates or inference failed validation."
        ),
        partial_draws,
        length(partial_idx)
      )
    } else {
      sprintf(
        paste0(
          "blmm: %d feature(s) were refit with lmerTest for all draws, and ",
          "%d draw(s) across %d additional feature(s) were refit ",
          "because their local covariance updates or inference failed validation."
        ),
        length(full_idx),
        partial_draws,
        length(partial_idx)
      )
    }
    warning(
      paste0(warning_text, " Examples: ", example_text, suffix),
      call. = FALSE
    )
  }

  # Each blmm_fit_feature() result has a singleton feature dimension (D = 1);
  # drop it when writing into the D-wide output arrays.
  output <- blmm_output_arrays(p, Pr, D, S, contrasts,
                               names = fnames, random_names = rnames)
  for (d in seq_len(D)) {
    feature_fit <- feature_results[[d]]$fit
    for (field in names(output))
      output[[field]][, d, ] <- feature_fit[[field]][, 1, ]
  }
  hook <- getOption("ALDEx3.blmm.adaptive.diagnostics")
  if (isTRUE(hook) || is.function(hook)) {
    diagnostics <- lapply(feature_results, function(result) {
      result$adaptive_diagnostics %||% list(
        diagnostics = data.frame(draw = seq_len(S),
          rescued = seq_len(S) %in% result$rejected_draws,
          reason = ifelse(seq_len(S) %in% result$rejected_draws, result$fallback_message, "accepted"),
          iterations = integer(S)),
        rescue_seconds = result$rescue_seconds %||% 0
      )
    })
    if (is.function(hook)) hook(diagnostics)
    if (isTRUE(hook)) attr(output, "blmm.adaptive.diagnostics") <- diagnostics
  }
  output
}
