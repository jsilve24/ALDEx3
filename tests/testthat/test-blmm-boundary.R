# Bounded covariance optimisation and singular (boundary) covariance fits used
# by adaptive inference. Sections:
#   1. blmm_bounded_optimize() on synthetic objectives
#   2. snapping near-bound coordinates exactly onto the bound
#   3. natural Cholesky AD tapes and coordinate invariance
#   4. partial (zero-column) boundaries
#   5. full zero-covariance boundaries

# Quadratics give known constrained optima and exact raw derivatives.
boundary_quadratic <- function(center, curvature) {
  list(
    fn = function(x) {
      delta <- x - center
      as.numeric(crossprod(delta, curvature %*% delta) / 2)
    },
    gr = function(x) as.vector(curvature %*% (x - center)),
    he = function(x) curvature
  )
}

# Natural Cholesky coordinates must retain zero variances and free correlations.
blmm_bounded_ad_objects <- function(z, backend, start = NULL, log = FALSE) {
  p <- ALDEx3:::blmm_parse_formula(z$formula, z$data, backend = backend,
                                    draws = ncol(z$Y))
  if (is.null(start)) start <- p$reTrms$theta
  is_log <- if (log) as.integer(is.finite(p$lower)) else integer(length(start))
  if (log) start[is_log == 1L] <- base::log(start[is_log == 1L])
  obj <- ALDEx3:::blmm_make_adfun(p$X, z$Y, p$basis, is_log, start)
  list(parsed = p, obj = obj)
}

# ---- 1. Bounded optimisation on synthetic objectives ----------------------

test_that("bounded optimization distinguishes zero-gradient maxima and minima", {
  maximum <- list(fn = function(x) (x[1]^2 - 1)^2,
                  gr = function(x) 4 * x[1] * (x[1]^2 - 1),
                  he = function(x) matrix(12 * x[1]^2 - 4))
  fit <- ALDEx3:::blmm_bounded_optimize(maximum, 0, 0)
  expect_true(fit$status == "unresolved" ||
                (fit$status == "interior" && fit$objective < maximum$fn(0)))
  expect_false(fit$status == "boundary" && fit$par[1] == 0)

  minimum <- list(fn = function(x) x[1]^2,
                  gr = function(x) 2 * x[1],
                  he = function(x) matrix(2))
  fit <- ALDEx3:::blmm_bounded_optimize(minimum, 0, 0)
  expect_identical(fit$status, "boundary")
  expect_equal(as.numeric(fit$par), 0)
  expect_equal(fit$H, matrix(2))
})

test_that("bounded optimization applies KKT signs only to constrained coordinates", {
  objective <- list(fn = function(x) (x[1] + 1)^2 + (x[2] + 2)^2,
                    gr = function(x) 2 * (x + c(1, 2)),
                    he = function(x) diag(2, 2))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 0), c(0, -Inf))
  expect_identical(fit$status, "boundary")
  expect_equal(as.numeric(fit$par), c(0, -2), tolerance = 1e-7)
  expect_equal(as.numeric(objective$gr(fit$par)), c(2, 0), tolerance = 1e-7)

  interior <- list(fn = function(x) (x[1] - 1)^2 + (x[2] + 2)^2,
                   gr = function(x) 2 * (x - c(1, -2)),
                   he = function(x) diag(2, 2))
  fit <- ALDEx3:::blmm_bounded_optimize(interior, c(0, 0), c(0, -Inf))
  expect_identical(fit$status, "interior")
  expect_equal(as.numeric(fit$par), c(1, -2), tolerance = 1e-7)
  expect_equal(fit$H, interior$he(fit$par))
})

test_that("interior convergence checks nearby boundary objectives", {
  quadratic <- list(fn = function(x) x[1]^2,
                    gr = function(x) 2 * x[1],
                    he = function(x) matrix(2))
  fit <- ALDEx3:::blmm_bounded_optimize(quadratic, 1e-6, 0)
  expect_true(fit$status %in% c("boundary", "unresolved"))
  expect_false(ALDEx3:::blmm_bounded_optimize(quadratic, 1e-6, 0,
                                               max_steps = 0)$status == "interior")

  quartic <- list(fn = function(x) x[1]^4,
                  gr = function(x) 4 * x[1]^3,
                  he = function(x) matrix(12 * x[1]^2))
  fit <- ALDEx3:::blmm_bounded_optimize(quartic, 1e-6, 0)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "boundary objective indistinguishable")

  interior <- list(fn = function(x) (x[1] - 1e-3)^2,
                   gr = function(x) 2 * (x[1] - 1e-3),
                   he = function(x) matrix(2))
  fit <- ALDEx3:::blmm_bounded_optimize(interior, 1e-3, 0)
  expect_identical(fit$status, "interior")
})

test_that("nonfinite inward boundary probes remain unresolved", {
  objective <- list(fn = function(x) if (x[1] == 0) 0 else NaN,
                    gr = function(x) 0, he = function(x) matrix(2))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, 0, 0)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "nonfinite boundary probe")
})

test_that("small scaled tangent correction permits boundary validation", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 1000)))
  start <- c(0, 1 + 1e-8)
  fit <- ALDEx3:::blmm_bounded_optimize(objective, start, c(0, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_gt(abs(objective$gr(start)[2]), 1e-6)
  expect_identical(fit$status, "boundary")
  expect_identical(fit$reason, "stationary boundary")
  expect_equal(fit$par, start)
})

test_that("large tangent correction does not qualify as stationary", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 1)))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1.01), c(0, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})

test_that("negative raw curvature prevents boundary acceptance", {
  curvature <- matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1000), 3, 3)
  objective <- boundary_quadratic(c(0, 0, 1), curvature)
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 0, 1 + 1e-8),
    c(0, 0, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "boundary negative curvature")
})

test_that("ill-conditioned tangent information cannot certify stationarity", {
  objective <- boundary_quadratic(c(-1, 1, 1), diag(c(1, 1, 1e-8)))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-5, 1),
    c(0, -Inf, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})

test_that("a small high-curvature inward optimum is not accepted at the bound", {
  objective <- boundary_quadratic(1e-8, matrix(1e8, 1, 1))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, 0, 0, max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_equal(objective$gr(0), -1)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})

test_that("the scaled route requires an explicit boundary validator", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 1000)))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-8),
    c(0, -Inf), max_steps = 0L)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})

test_that("a rejected boundary validation continues within the iteration budget", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 10)))
  calls <- 0L
  validator <- function(location, ...) {
    calls <<- calls + 1L
    list(ok = abs(location[2] - 1) < 1e-12, inference = "accepted")
  }
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-6),
    c(0, -Inf), max_steps = 1L, validate_boundary = validator)
  expect_identical(calls, 2L)
  expect_identical(fit$status, "boundary")
  expect_identical(fit$inference, "accepted")
  expect_identical(fit$iterations, 1L)
  expect_equal(fit$par, c(0, 1), tolerance = 1e-10)
  exhausted <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-6),
    c(0, -Inf), max_steps = 0L, validate_boundary = validator)
  expect_identical(exhausted$status, "unresolved")
  expect_identical(exhausted$iterations, 0L)
})

# ---- 2. Snapping near-bound coordinates -----------------------------------

test_that("coordinates within the bound tolerance are snapped exactly", {
  expect_identical(ALDEx3:::blmm_snap_to_bound(c(1e-9, 5e-11, -1, 1 + 1e-11),
                                               c(0, 0, -Inf, 1)),
                   c(1e-9, 0, -1, 1))
  expect_identical(ALDEx3:::blmm_near_bound(c(1e-9, 5e-11, -1), c(0, 0, -Inf)),
                   c(FALSE, TRUE, FALSE))
  # A near-zero start is moved onto the bound before its objective is used.
  square <- list(fn = function(x) x[1]^2, gr = function(x) 2 * x[1],
                 he = function(x) matrix(2))
  fit <- ALDEx3:::blmm_bounded_optimize(square, 1e-11, 0, max_steps = 0L)
  expect_identical(fit$status, "boundary")
  expect_identical(fit$par, 0)
  expect_identical(fit$objective, 0)
  # Line-search candidates are snapped too, with the objective re-evaluated.
  shifted <- list(fn = function(x) (x[1] - 1e-12)^2,
                  gr = function(x) 2 * (x[1] - 1e-12),
                  he = function(x) matrix(2))
  fit <- ALDEx3:::blmm_bounded_optimize(shifted, 1, 0)
  expect_identical(fit$status, "boundary")
  expect_identical(fit$par, 0)
  expect_identical(fit$objective, shifted$fn(0))
  expect_false(fit$failed)
})

test_that("a near-zero scalar random-intercept fit is accepted, not rescued", {
  old <- options(ALDEx3.blmm.adaptive.diagnostics = TRUE)
  on.exit(options(old), add = TRUE)
  z <- inference_fixture(seed = 41, groups = 12, size = 6, design = "within", S = 1)
  residual <- qr.resid(qr(cbind(1, z$data$x[seq_len(6)])), c(1, -1, 0.5, -0.5, 1, -1))
  z$Y[, 1] <- 2 + 0.7 * z$data$x + rep(residual, 12)
  z$logW <- array(z$Y, dim(z$logW))
  exact <- ALDEx3:::blmm_exact_feature(z$Y, z$formula, z$data)
  original <- ALDEx3:::blmm_bounded_optimize
  # Start each draw a hair above zero, inside the active-bound tolerance.
  testthat::local_mocked_bindings(blmm_bounded_optimize = function(obj, start, lower, ...) {
    if (!identical(list(...)$max_steps, 500L)) start[] <- 1e-11
    original(obj, start, lower, ...)
  }, .package = "ALDEx3")
  for (backend in c("dense", "sparse")) {
    options(ALDEx3.blmm.backend = backend)
    fit <- feature_fit(z, contrasts = NULL)
    expect_identical(fit$fallback_draws, 0L, info = backend)
    expect_identical(fit$adaptive_diagnostics$diagnostics$status, "boundary_zero",
                     info = backend)
    expect_equal(unname(fit$fit$df[, 1, 1]), rep(nrow(z$Y) - 2, 2), info = backend)
    for (field in names(exact))
      expect_equal(unname(fit$fit[[field]]), unname(exact[[field]]),
                   tolerance = 1e-5, info = paste(backend, field))
  }
})

# ---- 3. Natural Cholesky AD ------------------------------------------------

test_that("natural Cholesky AD tapes retain zero and free off-diagonal entries", {
  for (design in c("slope", "nested")) {
    z <- inference_fixture(seed = 914, groups = 12, design = design, S = 2)
    p <- ALDEx3:::blmm_parse_formula(z$formula, z$data)
    zero <- numeric(length(p$lower))
    dense <- blmm_bounded_ad_objects(z, "dense", zero)
    sparse <- blmm_bounded_ad_objects(z, "sparse", zero)
    bounded <- is.finite(p$lower)
    interior <- p$reTrms$theta
    interior[bounded] <- seq(0.4, 0.8, length.out = sum(bounded))
    interior[!bounded] <- -0.25
    partial <- interior
    partial[which(bounded)[1L]] <- 0

    # Recording at zero and subsequently evaluating nonzero parameters catches
    # lost symbolic entries in sparse factorisation and matrix products.
    for (theta in list(zero, partial, interior)) {
      expect_true(all(is.finite(sparse$obj$gr(theta))), info = design)
      expect_true(all(is.finite(sparse$obj$he(theta))), info = design)
      expect_equal(sparse$obj$fn(theta), dense$obj$fn(theta),
                   tolerance = 1e-9, info = design)
      expect_equal(sparse$obj$gr(theta), dense$obj$gr(theta),
                   tolerance = 1e-7, info = design)
      expect_equal(sparse$obj$he(theta), dense$obj$he(theta),
                   tolerance = 1e-7, info = design)

      step <- 1e-4
      gradient_fd <- numeric(length(theta))
      hessian_fd <- matrix(0, length(theta), length(theta))
      for (j in seq_along(theta)) {
        plus <- minus <- theta
        plus[j] <- plus[j] + step
        minus[j] <- minus[j] - step
        gradient_fd[j] <- (sparse$obj$fn(plus) - sparse$obj$fn(minus)) / (2 * step)
        hessian_fd[, j] <- as.numeric(sparse$obj$gr(plus) -
                                       sparse$obj$gr(minus)) / (2 * step)
      }
      expect_equal(as.numeric(sparse$obj$gr(theta)), gradient_fd,
                   tolerance = 1e-5, info = design)
      expect_equal(sparse$obj$he(theta), hessian_fd,
                   tolerance = 1e-5, info = design)
    }
  }
})

test_that("natural and log AD obey the nonstationary Hessian chain rule", {
  z <- inference_fixture(seed = 915, groups = 12, design = "slope", S = 2)
  theta <- c(0.7, -0.2, 0.4)
  natural <- blmm_bounded_ad_objects(z, "sparse", theta)
  logged <- blmm_bounded_ad_objects(z, "sparse", theta, log = TRUE)
  bounded <- is.finite(natural$parsed$lower)
  phi <- theta
  phi[bounded] <- log(theta[bounded])
  scale <- ifelse(bounded, theta, 1)
  D <- diag(scale)
  gradient <- as.numeric(natural$obj$gr(theta))
  extra <- diag(ifelse(bounded, theta * gradient, 0))
  expect_gt(max(abs(extra)), 1e-3)
  expect_equal(logged$obj$fn(phi), natural$obj$fn(theta), tolerance = 1e-9)
  expect_equal(as.numeric(logged$obj$gr(phi)), scale * gradient,
               tolerance = 1e-8)
  expect_equal(logged$obj$he(phi), D %*% natural$obj$he(theta) %*% D + extra,
               tolerance = 1e-7)
  expect_gt(max(abs(logged$obj$he(phi) -
                      D %*% natural$obj$he(theta) %*% D)), 1e-3)
})

test_that("natural scores retain centered draw gradients away from stationarity", {
  z <- inference_fixture(seed = 917, groups = 12, design = "slope", S = 3)
  theta <- c(0.7, -0.2, 0.4)
  batch <- blmm_bounded_ad_objects(z, "sparse", theta)
  p <- batch$parsed
  report <- batch$obj$report(theta)
  sensitivity <- ALDEx3:::blmm_scores_batch(
    batch$obj, theta, report$PWRSS_s, nrow(p$X), ncol(p$X),
    inference = TRUE)
  expect_gt(max(abs(batch$obj$gr(theta))), 1e-3)
  for (s in seq_len(ncol(z$Y))) {
    draw <- z
    draw$Y <- z$Y[, s, drop = FALSE]
    individual <- blmm_bounded_ad_objects(draw, "sparse", theta)
    expect_equal(sensitivity$scores[, s], -as.numeric(individual$obj$gr(theta)) +
                   as.numeric(batch$obj$gr(theta)),
                 tolerance = 1e-5, info = paste("draw", s))
  }
  expect_equal(rowMeans(sensitivity$scores), rep(0, length(theta)),
               tolerance = 1e-7)
})

test_that("interior slope and nested inference is consistent across coordinates", {
  for (design in c("slope", "nested")) {
    z <- inference_fixture(seed = 916, groups = 16, design = design, S = 2)
    natural <- blmm_bounded_ad_objects(z, "sparse")
    p <- natural$parsed
    fit <- stats::nlminb(p$reTrms$theta, natural$obj$fn, natural$obj$gr,
                         lower = p$lower,
                         control = list(rel.tol = 1e-10,
                                        eval.max = 1000, iter.max = 500))
    expect_identical(fit$convergence, 0L, info = design)
    theta <- fit$par
    bounded <- is.finite(p$lower)
    expect_true(all(theta[bounded] > 1e-3), info = design)
    expect_lt(max(abs(natural$obj$gr(theta))), 1e-4)
    logged <- blmm_bounded_ad_objects(z, "sparse", theta, log = TRUE)
    phi <- theta
    phi[bounded] <- log(theta[bounded])
    n <- nrow(p$X)
    q <- ncol(p$X)
    report <- natural$obj$report(theta)
    C <- diag(q)
    sens_theta <- ALDEx3:::blmm_scores_batch(
      natural$obj, theta, report$PWRSS_s, n, q, contrasts = C, inference = TRUE)
    sens_phi <- ALDEx3:::blmm_scores_batch(
      logged$obj, phi, report$PWRSS_s, n, q, contrasts = C, inference = TRUE)
    H_theta <- ALDEx3:::blmm_validate_H(natural$obj$he(theta))
    H_phi <- ALDEx3:::blmm_validate_H(logged$obj$he(phi))
    for (s in seq_len(ncol(z$Y))) {
      df_theta <- ALDEx3:::blmm_anchor_df(sens_theta$u[, s], sens_theta$b,
                                        chol(H_theta), n - q)
      df_phi <- ALDEx3:::blmm_anchor_df(sens_phi$u[, s], sens_phi$b,
                                      chol(H_phi), n - q)
      expect_equal(df_theta, df_phi, tolerance = 1e-4, info = design)
    }
  }
})

# ---- 4. Partial boundaries -------------------------------------------------

test_that("zero columns include all entries and respect block boundaries", {
  metadata <- ALDEx3:::blmm_random_effect_metadata(list(cnms = list(
    group = c("intercept", "slope", "other"), site = "intercept")))
  # The first block is ordered by columns: (11, 21, 31, 22, 32, 33).
  theta <- c(0, 0, 0, 2, 3, 4, 0)
  columns <- ALDEx3:::blmm_zero_columns(theta, metadata)
  expect_identical(columns$drop, c(1L, 2L, 3L, 7L))
  expect_identical(columns$keep, c(4L, 5L, 6L))
  theta[2] <- 1
  columns <- ALDEx3:::blmm_zero_columns(theta, metadata)
  expect_identical(columns$drop, 7L)
  lower <- c(0, -Inf, -Inf, 0, -Inf, 0, 0)
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    diag(7), rep(0, 7)), "unsupported active")
})

test_that("partial information checks raw curvature and decoupling", {
  metadata <- ALDEx3:::blmm_random_effect_metadata(list(cnms = list(
    group = c("intercept", "slope", "other"), site = "intercept")))
  theta <- c(0, 0, 0, 2, 3, 4, 0)
  lower <- c(0, -Inf, -Inf, 0, -Inf, 0, 0)
  H <- diag(c(0, 1, 2, 3, 4, 5, 0))
  gradient <- c(0, 0, 0, 0.01, -0.02, 0.03, 0)
  information <- ALDEx3:::blmm_partial_information(theta, lower, metadata, H, gradient)
  expect_identical(information$indices, c(4L, 5L, 6L))
  expect_identical(information$dropped, c(1L, 2L, 3L, 7L))
  expect_equal(crossprod(information$chol), H[4:6, 4:6])
  coupled <- H
  coupled[1, 4] <- coupled[4, 1] <- 0.01
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    coupled, gradient), "decouple")
  gradient[1] <- 1e-12
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    H, gradient), "decouple")
  gradient[1] <- 0
  negative <- H
  negative[1, 1] <- -1e-12
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    negative, gradient), "negative dropped curvature")
  bad <- H
  bad[4, 4] <- -1
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    bad, gradient), "not positive definite")
  bad[4, 4] <- 1e-7
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    bad, gradient), "ill-conditioned")
  bad <- H
  bad[4, 5] <- 0.1
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    bad, gradient), "not symmetric")
  bad <- H
  bad[1, 1] <- NA_real_
  expect_error(ALDEx3:::blmm_partial_information(theta, lower, metadata,
    bad, gradient), "invalid partial")
})

test_that("natural AD is symmetric at a partial zero slope column", {
  z <- inference_fixture(seed = 914, groups = 12, design = "slope", S = 1)
  C <- cbind(diag(2), c(1, -0.7))
  theta <- c(0.8, -0.25, 0)
  for (backend in c("dense", "sparse")) {
    parsed <- ALDEx3:::blmm_parse_formula(z$formula, z$data, backend = backend)
    obj <- ALDEx3:::blmm_make_adfun(parsed$X, z$Y, parsed$basis,
      integer(3), theta)
    flipped <- theta * c(-1, -1, 1)
    expect_equal(obj$fn(theta), obj$fn(flipped), tolerance = 1e-10, info = backend)
    plus <- minus <- theta
    plus[3] <- 0.05
    minus[3] <- -0.05
    expect_equal(obj$fn(plus), obj$fn(minus), tolerance = 1e-10, info = backend)
    H <- obj$he(theta)
    g <- as.vector(obj$gr(theta))
    expect_identical(g[3], 0, info = backend)
    expect_true(all(H[3, 1:2] == 0), info = backend)
    report <- obj$report(theta)
    sensitivity <- ALDEx3:::blmm_scores_batch(obj, theta, report$PWRSS_s,
      nrow(parsed$X), ncol(parsed$X), contrasts = C, inference = TRUE)
    expect_true(all(sensitivity$u[3, ] == 0), info = backend)
    expect_true(all(sensitivity$b[3, ] == 0), info = backend)
    h <- sensitivity$b + as.vector(sensitivity$u[, 1])
    full <- colSums(h * solve(H, h))
    reduced <- colSums(h[1:2, , drop = FALSE] *
      solve(H[1:2, 1:2], h[1:2, , drop = FALSE]))
    expect_equal(full, reduced, tolerance = 1e-10, info = backend)
    # The first column has two entries, both of which must flip together.
    theta_first_zero <- c(0, 0, 0.7)
    plus <- c(0.05, -0.03, 0.7)
    minus <- c(-0.05, 0.03, 0.7)
    expect_equal(obj$fn(plus), obj$fn(minus), tolerance = 1e-10, info = backend)
    expect_true(all(obj$gr(theta_first_zero)[1:2] == 0), info = backend)
    expect_true(all(obj$he(theta_first_zero)[1:2, 3] == 0), info = backend)
  }
})

test_that("accepted nested partial boundary agrees with exact inference", {
  old <- options(ALDEx3.blmm.adaptive.diagnostics = TRUE,
                 ALDEx3.blmm.backend = "auto")
  on.exit(options(old), add = TRUE)
  z <- inference_fixture(seed = 1, groups = 12, design = "nested", S = 1, noise = 0)
  C <- matrix(c(1, -0.7), 1, dimnames = list("combination", c("(Intercept)", "x")))
  exact <- ALDEx3:::blmm_exact_feature(z$Y, z$formula, z$data, contrasts = C)
  for (backend in c("dense", "sparse")) {
    options(ALDEx3.blmm.backend = backend)
    fit <- feature_fit(z, inference = "adaptive", contrasts = C)
    expect_identical(fit$adaptive_diagnostics$diagnostics$status, "boundary_partial", info = backend)
    expect_identical(fit$fallback_draws, 0L, info = backend)
    for (field in c("estimate", "std.error", "df", "p.lower", "p.upper",
                    "contrast.estimate", "contrast.std.error", "contrast.df",
                    "contrast.p.lower", "contrast.p.upper"))
      expect_equal(unname(fit$fit[[field]]), unname(exact[[field]]),
        tolerance = 0.002, info = paste(backend, field))
  }
})

test_that("partial boundary inference preserves shifts and worker ordering", {
  old <- options(ALDEx3.blmm.adaptive.diagnostics = TRUE)
  on.exit(options(old), add = TRUE)
  z <- inference_fixture(seed = 1, groups = 12, design = "nested", S = 1, noise = 0)
  X <- model.matrix(~ x, z$data)
  shift <- c(100, -20)
  responses <- array(c(z$Y, z$Y + as.vector(X %*% shift)), c(nrow(X), 2L, 1L))
  serial <- ALDEx3:::blmm(responses, z$formula, z$data, inference = "adaptive")
  workers <- ALDEx3:::blmm(responses, z$formula, z$data, inference = "adaptive", n.cores = 2L)
  expect_equal(workers, serial)
  expect_equal(serial$estimate[, 2L, 1L], serial$estimate[, 1L, 1L] + shift,
               tolerance = 1e-6)
  for (field in c("std.error", "df", "random.eff"))
    expect_equal(serial[[field]][, 2L, 1L], serial[[field]][, 1L, 1L],
                 tolerance = 1e-5, info = field)
})

# ---- 5. Full zero boundaries ----------------------------------------------

test_that("all-zero covariance matches exact inference on both backends", {
  old <- options(ALDEx3.blmm.adaptive.diagnostics = TRUE,
                 ALDEx3.blmm.backend = "auto")
  on.exit(options(old), add = TRUE)
  original <- ALDEx3:::blmm_bounded_optimize
  # Exercise the real draw search at an exact zero anchor. The unconstrained
  # off-diagonal can otherwise retain roundoff after the batch search.
  testthat::local_mocked_bindings(blmm_bounded_optimize = function(obj, start, lower, ...) {
    args <- list(...)
    if (identical(args$max_steps, 500L)) start[] <- 0
    original(obj, start, lower, ...)
  }, .package = "ALDEx3")
  for (design in c("within", "slope")) {
    z <- inference_fixture(seed = 41, groups = 12, size = 6, design = design, S = 1)
    x <- z$data$x[seq_len(6)]
    residual <- qr.resid(qr(cbind(1, x)), c(1, -1, 0.5, -0.5, 1, -1))
    z$Y[, 1] <- 2 + 0.7 * z$data$x + rep(residual, 12)
    z$logW <- array(z$Y, dim(z$logW))
    C <- matrix(c(1, -0.5), 1, dimnames = list("combination", c("(Intercept)", "x")))
    reference <- suppressMessages(lmerTest::lmer(
      update(z$formula, y ~ .), data = transform(z$data, y = z$Y[, 1])))
    # theta is a Cholesky factor relative to the residual SD. Use lme4's
    # default boundary tolerance: optimizer roundoff need not be exactly zero.
    # For these fixtures, this bounds relative variances by at most 2e-10.
    expect_true(all(abs(lme4::getME(reference, "theta")) < 1e-5), info = design)
    exact <- ALDEx3:::blmm_exact_feature(z$Y, z$formula, z$data, contrasts = C)
    for (backend in c("dense", "sparse")) {
      options(ALDEx3.blmm.backend = backend)
      fit <- feature_fit(z, inference = "adaptive", contrasts = C)
      expect_identical(fit$adaptive_diagnostics$diagnostics$status, "boundary_zero", info = design)
      expect_identical(fit$fallback_draws, 0L, info = design)
      for (field in names(exact))
        expect_equal(unname(fit$fit[[field]]), unname(exact[[field]]),
                     tolerance = 1e-5, info = paste(design, backend, field))
    }
  }
})

test_that("only established full zero boundaries bypass rescue", {
  z <- inference_fixture(seed = 41, groups = 12, size = 6, design = "slope", S = 1)
  original <- ALDEx3:::blmm_bounded_optimize
  for (scenario in c("unresolved", "partial", "tiny")) {
    calls <- 0L
    testthat::local_mocked_bindings(blmm_bounded_optimize = function(obj, start, lower, ...) {
      calls <<- calls + 1L
      if (calls == 1L) return(original(obj, start, lower, ...))
      par <- rep(0, length(start))
      if (scenario == "partial") par[2] <- 0.1
      if (scenario == "tiny") par[1] <- 1e-12
      validation <- if (scenario == "unresolved") NULL else
        list(...)$validate_boundary(par, obj$he(par), obj$fn(par), obj$gr(par))
      list(status = if (scenario == "unresolved") "unresolved" else "boundary",
           reason = "injected", par = par, iterations = 0L,
           inference = validation$inference,
           diagnostics = list(evaluations = 0L))
    }, .package = "ALDEx3")
    fit <- feature_fit(z, inference = "adaptive")
    expect_identical(fit$fallback_draws, 1L, info = scenario)
  }
})

test_that("zero boundary with invalid raw curvature is rescued", {
  z <- inference_fixture(seed = 41, groups = 12, size = 6, design = "slope", S = 1)
  original <- ALDEx3:::blmm_bounded_optimize
  for (curvature in c(0, -1)) {
    calls <- 0L
    testthat::local_mocked_bindings(blmm_bounded_optimize = function(obj, start, lower, ...) {
      calls <<- calls + 1L
      if (calls == 1L) return(original(obj, start, lower, ...))
      par <- rep(0, length(start))
      H <- diag(curvature, length(start))
      validation <- list(...)$validate_boundary(par, H, obj$fn(par), obj$gr(par))
      list(status = "boundary", reason = "injected", par = rep(0, length(start)),
           H = H, iterations = 0L, inference = validation$inference,
           diagnostics = list(evaluations = 0L))
    }, .package = "ALDEx3")
    fit <- feature_fit(z, inference = "adaptive")
    expect_identical(fit$fallback_draws, 1L, info = paste("curvature", curvature))
  }
})
