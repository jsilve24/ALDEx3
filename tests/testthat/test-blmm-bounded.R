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

test_that("nonfinite inward boundary probes remain unresolved", {
  objective <- list(fn = function(x) if (x[1] == 0) 0 else NaN,
                    gr = function(x) 0, he = function(x) matrix(2))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, 0, 0)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "nonfinite boundary probe")
})
