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
    expect_true(all(lme4::getME(reference, "theta") == 0), info = design)
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
