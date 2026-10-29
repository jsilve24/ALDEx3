test_that("explicit TMB reports are independent of the previous parameter", {
  z <- inference_fixture(S = 3L)
  for (backend in c("dense", "sparse")) {
    p <- ALDEx3:::blmm_parse_formula(z$formula, z$data, backend = backend)
    obj <- ALDEx3:::blmm_make_adfun(p$X, z$Y, p$basis, p$is_log, p$phi_init)
    for (phi in list(p$phi_init, p$phi_init + .1)) {
      obj$fn(phi)
      expected <- obj$report()
      obj$fn(phi + .3)
      expect_equal(obj$report(phi), expected, tolerance = 0)
    }
  }
})

test_that("adaptive interior inference agrees with exact fitting and restores shifts", {
  z <- inference_fixture(groups = 12, S = 3L, design = "within")
  X <- model.matrix(~ x, z$data)
  C <- matrix(c(1, -.5), 1, dimnames = list("combination", colnames(X)))
  a <- feature_fit(z, inference = "adaptive", contrasts = C)
  expect_identical(a$fallback_draws, 0L)
  e <- ALDEx3:::blmm_exact_feature(z$Y, z$formula, z$data, contrasts = C)
  for (field in names(e))
    expect_equal(unname(a$fit[[field]]), unname(e[[field]]), tolerance = .002,
                 info = field)
  shift <- rbind(c(100, -50, 25), c(-30, 20, 10))
  zz <- z
  zz$Y <- z$Y + X %*% shift
  zz$logW <- array(zz$Y, dim(z$logW))
  b <- feature_fit(zz, inference = "adaptive", contrasts = C)
  expect_identical(b$rejected_draws, a$rejected_draws)
  expect_equal(b$fit$estimate, a$fit$estimate + array(shift, dim(a$fit$estimate)),
               tolerance = 1e-5)
  expect_equal(b$fit$contrast.estimate,
               a$fit$contrast.estimate + array(C %*% shift, dim(a$fit$contrast.estimate)),
               tolerance = 1e-5)
  for (field in c("std.error", "df", "random.eff", "contrast.std.error", "contrast.df"))
    expect_equal(b$fit[[field]], a$fit[[field]], tolerance = 1e-5, info = field)
  expect_equal(b$fit$p.upper,
               pt(b$fit$estimate / b$fit$std.error, b$fit$df, lower.tail = FALSE))
})

test_that("adaptive dense, sparse, and worker paths preserve outputs", {
  z <- inference_fixture(groups = 12, design = "slope", S = 3L)
  old <- options(ALDEx3.blmm.backend = "dense")
  on.exit(options(old), add = TRUE)
  logW <- abind::abind(z$logW, z$logW + .2, along = 2)
  a <- ALDEx3:::blmm(logW, z$formula, z$data, inference = "adaptive")
  options(ALDEx3.blmm.backend = "sparse")
  b <- ALDEx3:::blmm(logW, z$formula, z$data, inference = "adaptive")
  expect_equal(b, a, tolerance = 1e-4)
  c <- ALDEx3:::blmm(logW, z$formula, z$data, inference = "adaptive", n.cores = 2L)
  expect_equal(c, b)
})

test_that("adaptive forced rescue retains exact fields and warning aggregation", {
  z <- inference_fixture(S = 3L)
  C <- matrix(c(0, 1), 1, dimnames = list("slope", c("(Intercept)", "x")))
  attr(z$logW, "blmm.force_fail_features") <- 1L
  expect_warning(a <- ALDEx3:::blmm(z$logW, z$formula, z$data,
                                   inference = "adaptive", contrasts = C), "fell back")
  e <- ALDEx3:::blmm_exact_feature(z$Y, z$formula, z$data, contrasts = C)
  for (field in names(e))
    expect_equal(unname(a[[field]]), unname(e[[field]]), info = field)
})

test_that("adaptive partial rescues are batched and exact failures are not retried", {
  z <- inference_fixture(groups = 12, S = 3L, design = "within")
  original_object <- ALDEx3:::blmm_make_adfun
  original_exact <- ALDEx3:::blmm_exact_feature
  draw_calls <- exact_calls <- 0L
  fail_exact <- FALSE
  testthat::local_mocked_bindings(
    blmm_make_adfun = function(X, Y, ...) {
      if (ncol(Y) == 1L) {
        draw_calls <<- draw_calls + 1L
        if (draw_calls %in% c(1L, 3L)) stop("injected draw failure")
      }
      original_object(X, Y, ...)
    },
    blmm_exact_feature = function(Y_d, ...) {
      exact_calls <<- exact_calls + 1L
      if (fail_exact) stop("injected exact failure")
      expect_equal(ncol(Y_d), 2L)
      expect_equal(Y_d, z$Y[, c(1L, 3L)])
      original_exact(Y_d, ...)
    }, .package = "ALDEx3")
  a <- feature_fit(z, inference = "adaptive")
  expect_identical(a$rejected_draws, c(1L, 3L))
  expect_false(a$full_fallback)
  expect_identical(exact_calls, 1L)
  e <- original_exact(z$Y[, c(1L, 3L)], z$formula, z$data)
  for (field in names(e))
    expect_equal(unname(a$fit[[field]][, , c(1L, 3L), drop = FALSE]), unname(e[[field]]))
  draw_calls <- exact_calls <- 0L
  fail_exact <- TRUE
  expect_error(feature_fit(z, inference = "adaptive"),
               "exact refit failed for feature 1.*draws 1, 3.*injected exact failure")
  expect_identical(exact_calls, 1L)
})

test_that("adaptive default routes through streaming on both backends", {
  z <- inference_fixture(S = 4L, noise = .002)
  N <- nrow(z$Y)
  saved <- aperm(z$logW, c(2, 1, 3))
  cursor <- 0L
  testthat::local_mocked_bindings(aldex.sampler = function(Y, X, nsample, ...) {
    take <- cursor + seq_len(nsample)
    cursor <<- cursor + nsample
    list(logW = saved[, , take, drop = FALSE])
  }, .package = "ALDEx3")
  C <- matrix(c(1, -.5), 1, dimnames = list("combination", c("(Intercept)", "x")))
  run <- function(test = NULL, backend = "dense") {
    old <- options(ALDEx3.blmm.backend = backend)
    on.exit(options(old))
    cursor <<- 0L
    ALDEx3::aldex(matrix(10, 1, N), z$formula, z$data, method = "blmm",
                  test = test, n.cores = 1L, nsample = 4L, scale = matrix(0, N, 4),
                  contrasts = C, streamsize = N * 8 * 2 / 100000,
                  return.pars = c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  }
  for (backend in c("dense", "sparse")) {
    default <- run(backend = backend)
    expect_identical(default$test, "t.satterthwaite.adaptive")
    expect_equal(default, run("t.satterthwaite.adaptive", backend))
  }
  a <- run("t.satterthwaite.adaptive")
  expect_identical(a$test, "t.satterthwaite.adaptive")
  expect_true(a$streaming)
  expect_equal(cursor, 4L)
  expect_equal(a, run("t.satterthwaite.adaptive", "sparse"), tolerance = 1e-4)
  deprecations <- character()
  residual <- withCallingHandlers(run("t.residual"), warning = function(w) {
    if (grepl("t.residual.*deprecated", conditionMessage(w))) {
      deprecations <<- c(deprecations, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  })
  expect_length(deprecations, 1L)
  expect_true(residual$streaming)
  expect_identical(residual$test, "t.residual")
  expect_error(ALDEx3::aldex(matrix(1, 1, N), ~ x, z$data,
                            method = "lm", test = "t.satterthwaite.adaptive"), "not supported")
})
