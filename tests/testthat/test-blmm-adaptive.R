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
                                   inference = "adaptive", contrasts = C), "refit with lmerTest")
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

test_that("an unresolved but finite shared start does not rescue the feature", {
  old <- options(ALDEx3.blmm.adaptive.diagnostics = TRUE)
  on.exit(options(old), add = TRUE)
  z <- inference_fixture(groups = 12, S = 3L, design = "within")
  baseline <- feature_fit(z)
  original <- ALDEx3:::blmm_bounded_optimize
  mode <- "unresolved"
  testthat::local_mocked_bindings(blmm_bounded_optimize = function(obj, start, lower, ...) {
    fit <- original(obj, start, lower, ...)
    if (identical(list(...)$max_steps, 500L)) {
      fit$status <- "unresolved"
      fit$reason <- "boundary objective indistinguishable"
      fit$failed <- mode == "failed"
    }
    fit
  }, .package = "ALDEx3")
  a <- feature_fit(z)
  expect_false(a$full_fallback)
  expect_identical(a$fallback_draws, 0L)
  expect_identical(unique(a$adaptive_diagnostics$diagnostics$anchor_status), "unresolved")
  for (field in names(baseline$fit))
    expect_equal(a$fit[[field]], baseline$fit[[field]], info = field)
  # A start that ended in an error still sends the whole feature to rescue.
  mode <- "failed"
  b <- feature_fit(z)
  expect_true(b$full_fallback)
  expect_match(b$fallback_message, "adaptive anchor failed")
})

test_that("adaptive inference calculates small upper tails directly", {
  z <- inference_fixture(S=1,noise=0,design="within",effect=30)
  a <- feature_fit(z)
  e <- ALDEx3:::blmm_exact_feature(z$Y,z$formula,z$data)
  expect_identical(a$fallback_draws,0L)
  expect_gt(a$fit$p.upper[2],0)
  expect_lt(a$fit$p.upper[2],.Machine$double.eps)
  expect_gt(e$p.upper[2],0)
  expect_lt(abs(log(a$fit$p.upper[2])-log(e$p.upper[2])),0.002)
})

# ---- Agreement with per-draw lmerTest refits --------------------------------
# lmerTest's Satterthwaite df use a numerical Hessian, so agreement is checked
# to 2e-3 relative (estimates relative to max(|estimate|, SE)). Draws where lmer
# reports a convergence problem are not an accuracy oracle and are skipped.

adaptive_lmer_reference <- function(formula, data, y, C = NULL) {
  data$y <- y
  converged <- TRUE
  fit <- withCallingHandlers(
    suppressMessages(lmerTest::lmer(update(formula, y ~ .), data = data, REML = TRUE)),
    warning = function(w) {
      if (!grepl("singular", conditionMessage(w))) converged <<- FALSE
      invokeRestart("muffleWarning")
    })
  if (any(!grepl("singular", fit@optinfo$conv$lme4$messages))) converged <- FALSE
  coefs <- summary(fit)$coefficients
  out <- list(converged = converged, singular = lme4::isSingular(fit),
              estimate = coefs[, "Estimate"], std.error = coefs[, "Std. Error"],
              df = coefs[, "df"])
  if (!is.null(C)) {
    contrast <- lmerTest::contest(fit, C, joint = FALSE)
    out$contrast.estimate <- contrast$Estimate
    out$contrast.std.error <- contrast[["Std. Error"]]
    out$contrast.df <- contrast$df
  }
  out
}

# Fits the default BLMM on each backend and compares every converged draw with
# lmerTest. Returns the number of compared and singular reference draws.
expect_adaptive_matches_lmer <- function(logW, formula, data, C = NULL,
                                         singular_df = FALSE, tolerance = 2e-3) {
  D <- dim(logW)[2]
  S <- dim(logW)[3]
  residual_df <- nrow(data) - ncol(model.matrix(reformulas::nobars(formula), data))
  references <- lapply(seq_len(D), function(d) lapply(seq_len(S), function(s)
    adaptive_lmer_reference(formula, data, logW[, d, s], C)))
  fields <- c("estimate", "std.error", "df",
              if (!is.null(C)) c("contrast.estimate", "contrast.std.error", "contrast.df"))
  old <- options(ALDEx3.blmm.backend = getOption("ALDEx3.blmm.backend"))
  on.exit(options(old), add = TRUE)
  counts <- c(compared = 0, singular = 0)
  for (backend in c("dense", "sparse")) {
    options(ALDEx3.blmm.backend = backend)
    fit <- suppressWarnings(ALDEx3:::blmm(logW, formula, data, contrasts = C))
    worst <- stats::setNames(numeric(length(fields)), fields)
    counts[] <- 0
    for (d in seq_len(D)) for (s in seq_len(S)) {
      reference <- references[[d]][[s]]
      if (!reference$converged) next
      counts["compared"] <- counts["compared"] + 1
      for (field in fields) {
        actual <- fit[[field]][, d, s]
        expected <- reference[[field]]
        scale <- abs(expected)
        if (field == "estimate") scale <- pmax(scale, reference$std.error)
        if (field == "contrast.estimate") scale <- pmax(scale, reference$contrast.std.error)
        worst[field] <- max(worst[field], abs(actual - expected) / scale)
      }
      if (singular_df && reference$singular) {
        counts["singular"] <- counts["singular"] + 1
        expect_equal(unname(fit$df[, d, s]), rep(residual_df, length(fit$df[, d, s])),
                     tolerance = 1e-8, info = paste(backend, d, s))
        if (!is.null(C))
          expect_equal(unname(fit$contrast.df[, d, s]), rep(residual_df, nrow(C)),
                       tolerance = 1e-8, info = paste(backend, d, s))
      }
    }
    for (field in fields)
      expect_lt(worst[[field]], tolerance,
                label = paste(deparse(formula), backend, field, "relative error"))
  }
  expect_gt(counts[["compared"]], 0)
  counts
}

test_that("adaptive random intercepts and contrasts match lmerTest, including singular draws", {
  skip_on_cran()
  set.seed(1)
  N <- 40
  data <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(N),
                     z = rbinom(N, 1, 0.5))
  # Zero and near-zero group SDs give singular draws; the others are interior.
  sds <- c(0, 0.05, 0.3, 1, 0, 0.5)
  logW <- array(NA_real_, c(N, length(sds), 5L))
  for (d in seq_along(sds)) for (s in 1:5)
    logW[, d, s] <- 0.5 * data$x + rnorm(8, sd = sds[d])[data$g] + rnorm(N)
  C <- matrix(c(0, 1, -1), 1, dimnames = list("x_minus_z", c("(Intercept)", "x", "z")))
  counts <- expect_adaptive_matches_lmer(logW, ~ x + z + (1 | g), data, C,
                                         singular_df = TRUE)
  expect_gt(counts[["singular"]], 0)
})

test_that("adaptive slopes, crossed, and uncorrelated designs match lmerTest", {
  skip_on_cran()
  set.seed(2)
  N <- 60
  data <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(N),
                     h = factor(rep(1:6, 10)))
  sdi <- c(0, 0.5, 1, 0.3, 0, 1, 0.1, 0.8)
  sds <- c(0, 0, 0.5, 0.5, 1, 1, 0.05, 0)
  logW <- array(NA_real_, c(N, length(sdi), 4L))
  for (d in seq_along(sdi)) for (s in 1:4)
    logW[, d, s] <- 0.5 * data$x + rnorm(10, sd = sdi[d])[data$g] +
      rnorm(10, sd = sds[d])[data$g] * data$x + rnorm(6, sd = 0.3)[data$h] + rnorm(N)
  C <- matrix(c(1, -0.5), 1, dimnames = list("combination", c("(Intercept)", "x")))
  for (formula in list(~ x + (1 + x | g), ~ x + (1 | g) + (1 | h),
                       ~ x + (1 | g) + (0 + x | g)))
    expect_adaptive_matches_lmer(logW, formula, data, C)
})

aldex_blmm_fixture <- function() {
  set.seed(11)
  data <- data.frame(g = factor(rep(1:8, each = 4)), x = rep(c(-0.5, 0.5), 16))
  u <- matrix(rnorm(8 * 3, sd = 0.6), 8, 3)
  lambda <- exp(4 + outer(c(0.8, -0.4, 0), data$x) + t(u[data$g, ]))
  Y <- matrix(rpois(length(lambda), lambda), nrow(lambda),
              dimnames = list(paste0("taxon", 1:3), paste0("s", seq_len(nrow(data)))))
  list(Y = Y, data = data)
}

test_that("aldex blmm defaults to adaptive inference that matches lmerTest", {
  z <- aldex_blmm_fixture()
  set.seed(12)
  fit <- suppressWarnings(ALDEx3::aldex(
    z$Y, ~ x + (1 | g), z$data, method = "blmm", nsample = 4, n.cores = 1,
    scale = ALDEx3::clr.sm, gamma = 0.5,
    return.pars = c("estimate", "std.error", "df", "logComp", "logScale")))
  expect_identical(fit$test, "t.satterthwaite.adaptive")
  log_scale <- matrix(fit$logScale, nrow(z$data))
  worst <- c(estimate = 0, std.error = 0, df = 0)
  compared <- 0
  for (d in seq_len(nrow(z$Y))) for (s in seq_len(4)) {
    y <- fit$logComp[d, , s] + log_scale[, s]
    reference <- adaptive_lmer_reference(~ x + (1 | g), z$data, y)
    if (!reference$converged) next
    compared <- compared + 1
    scale <- list(estimate = pmax(abs(reference$estimate), reference$std.error),
                  std.error = reference$std.error, df = reference$df)
    for (field in names(worst))
      worst[field] <- max(worst[field],
        abs(fit[[field]][, d, s] - reference[[field]]) / scale[[field]])
  }
  expect_gt(compared, 0)
  for (field in names(worst)) expect_lt(worst[[field]], 2e-3, label = field)
})

test_that("aldex blmm with forked workers is identical to a serial run", {
  skip_on_cran()
  skip_on_os("windows")
  z <- aldex_blmm_fixture()
  C <- matrix(c(1, -0.5), 1, dimnames = list("combination", c("(Intercept)", "x")))
  formula <- ~ x + (1 | g)
  run <- function(cores) {
    set.seed(13)
    suppressWarnings(ALDEx3::aldex(
      z$Y, formula, z$data, method = "blmm", nsample = 6, n.cores = cores,
      scale = ALDEx3::clr.sm, gamma = 0.5, contrasts = C, return.pars = c("X", "estimate", "std.error", "df", "p.val",
                                     "p.val.adj", "random.effects")))
  }
  serial <- run(1L)
  forked <- run(2L)
  expect_identical(forked, serial)
})
