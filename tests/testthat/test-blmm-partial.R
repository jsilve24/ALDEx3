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
