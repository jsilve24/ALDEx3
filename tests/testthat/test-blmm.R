# Tests for the approximate BLMM mixed-effects engine.
#
# Tolerance rationale (set prospectively):
#   A. With nsample = 1, BLMM should closely match exact lme4 because there is
#      no posterior averaging across draws. Empirically the gap is around
#      1e-6 or smaller on these problems, so 1e-5 is a safe but strict bound.
#   B. Identical draws should produce numerically identical outputs because the
#      draw-specific score corrections are zero.
#   C. Across repeated draws, posterior mean error should stay around 1e-4 to
#      1e-3 and posterior SD relative error around 1e-3 to 1e-2 on these
#      shared-draw comparisons, so the approximation remains materially smaller
#      than the posterior uncertainty being reported.
#   D. Random-slope and correlated-random-effect output should remain close
#      enough to exact lme4 to preserve interpretation of the covariance terms.
#   E. If the approximation fails, BLMM must warn and fall back to exact lme4
#      rather than silently changing the model.

blmm_true_scale <- function(sim, nsample) {
  matrix(rep(log2(colSums(sim$W)), nsample), nrow = ncol(sim$Y))
}

blmm_sample_logW <- function(Y, formula, data, nsample, scale, scale.args = list()) {
  X <- t(model.matrix(reformulas::nobars(formula), data))
  draws <- ALDEx3:::aldex.sampler(
    Y, X, nsample, scale, scale.args, return.pars = "logW"
  )
  aperm(draws$logW, c(2, 1, 3))
}

blmm_fit_pair <- function(logW, formula, data) {
  list(
    blmm = ALDEx3:::blmm(logW, formula, data, n.cores = 1L),
    lme4 = ALDEx3:::sr.mem(logW, formula, data, n.cores = 1L,
                           method = "lme4", mem.args = list())
  )
}

blmm_vignette_subset <- function() {
  data(oral_mouthwash_data, package = "ALDEx3")

  Y0 <- oral_mouthwash_data$counts
  keep_names <- row.names(Y0[((rowSums(Y0 == 0)) / ncol(Y0)) <= 0.75, ])
  other <- colSums(Y0[((rowSums(Y0 == 0)) / ncol(Y0)) > 0.75, ])
  Y <- rbind(Y0[keep_names, ], other = other)
  Y <- Y[c(seq_len(min(8, nrow(Y) - 1)), nrow(Y)), ]

  list(
    Y = Y,
    meta = oral_mouthwash_data$metadata
  )
}

test_that("blmm: expected input errors", {
  set.seed(4985)
  nsample <- 4
  sim <- ALDEx3:::aldex.mem.sim(5, 4, 4, 10000, 1, FALSE, 0, 0, 0.1)
  logScale <- blmm_true_scale(sim, nsample)

  expect_error(
    ALDEx3::aldex(sim$Y, sim$meta, data = sim$meta, method = "blmm",
          nsample = nsample, scale = logScale, n.cores = 1),
    "X should be a mixed effects formula"
  )
  expect_error(
    ALDEx3::aldex(sim$Y, ~treatment + (1|subject_ids), data = NULL,
          method = "blmm", nsample = nsample, scale = logScale, n.cores = 1),
    "data should not be null if method is a mixed-effects engine"
  )
})

test_that("blmm A: single-draw agreement with exact lme4", {
  set.seed(42)
  nsample <- 1
  formula <- ~treatment + (1|subject_ids)
  sim <- ALDEx3:::aldex.mem.sim(D = 5, days = 5, subjects = 6,
                       depth = 100000, sd_resid = 0.05)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)
  fits <- blmm_fit_pair(logW, formula, sim$meta)

  expect_equal(fits$blmm$estimate, fits$lme4$estimate, tolerance = 1e-5)
  expect_equal(fits$blmm$std.error, fits$lme4$std.error, tolerance = 1e-5)
})

test_that("blmm A2: anchor objective is the average profiled REML criterion", {
  set.seed(2026)
  sim <- ALDEx3:::aldex.mem.sim(D = 1, days = 4, subjects = 6,
                       depth = 100000, sd_resid = 0.1)
  parsed <- ALDEx3:::blmm_parse_formula(~treatment + (1|subject_ids), sim$meta)

  base_y <- log2(sim$Y[1, ] + 0.5)
  Y_d <- vapply(
    seq_len(6),
    function(s) base_y + stats::rnorm(length(base_y), sd = 0.05 * s),
    numeric(length(base_y))
  )

  obj <- ALDEx3:::blmm_make_adfun(
    parsed$X, Y_d, parsed$basis, parsed$is_log, parsed$phi_init
  )
  anchor <- ALDEx3:::blmm_fit_anchor(obj, parsed$phi_init)
  value <- obj$fn(anchor$phi_bar)
  report <- obj$report()

  np <- nrow(parsed$X) - ncol(parsed$X)
  surrogate <- 0.5 * (report$ldL2 + report$ldRX2 +
                        np * log(mean(report$PWRSS_s) / np))
  exact_average <- ALDEx3:::blmm_anchor_objective_from_report(
    report, nrow(parsed$X), ncol(parsed$X)
  )
  scores <- ALDEx3:::blmm_scores_batch(
    obj, anchor$phi_bar, report$PWRSS_s, nrow(parsed$X), ncol(parsed$X)
  )

  expect_equal(value, exact_average, tolerance = 1e-8)
  expect_gt(abs(value - surrogate), 1e-5)
  expect_true(all(abs(rowMeans(scores)) < 1e-6))
})

test_that("blmm B: identical-draw consistency", {
  set.seed(123)
  nsample <- 20
  sim <- ALDEx3:::aldex.mem.sim(D = 5, days = 4, subjects = 5,
                       depth = 100000, sd_resid = 0.05)
  true.S <- log2(colSums(sim$W))

  logp <- log2(sweep(sim$Y + 0.5, 2, colSums(sim$Y + 0.5), "/"))
  logW1 <- logp + matrix(true.S, nrow(sim$Y), ncol(sim$Y), byrow = TRUE)
  logW <- aperm(array(logW1, c(nrow(sim$Y), ncol(sim$Y), nsample)), c(2, 1, 3))

  res_blmm <- ALDEx3:::blmm(logW, ~treatment + (1|subject_ids),
                            sim$meta, n.cores = 1L)

  expect_true(all(apply(res_blmm$estimate, c(1, 2), var) < 1e-8))
  expect_true(all(apply(res_blmm$std.error, c(1, 2), var) < 1e-8))
  expect_true(all(apply(res_blmm$random.eff, c(1, 2), var) < 1e-8))
})

test_that("blmm C: posterior mean and SD remain close to exact lme4", {
  set.seed(999)
  nsample <- 40
  formula <- ~treatment + (1|subject_ids)
  sim <- ALDEx3:::aldex.mem.sim(D = 8, days = 4, subjects = 8,
                       depth = 100000, sd_resid = 0.1)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)
  fits <- blmm_fit_pair(logW, formula, sim$meta)

  blmm_mean <- apply(fits$blmm$estimate[2, , , drop = FALSE], c(1, 2), mean)
  lme4_mean <- apply(fits$lme4$estimate[2, , , drop = FALSE], c(1, 2), mean)
  expect_equal(as.vector(blmm_mean), as.vector(lme4_mean), tolerance = 0.001)

  blmm_sd <- apply(fits$blmm$estimate[2, , , drop = FALSE], c(1, 2), sd)
  lme4_sd <- apply(fits$lme4$estimate[2, , , drop = FALSE], c(1, 2), sd)
  rel_err <- abs(blmm_sd - lme4_sd) / (abs(lme4_sd) + 0.01)
  expect_true(all(as.vector(rel_err) < 0.005))
})

test_that("blmm C2: small repeated-measures problem stays close to exact lme4", {
  set.seed(314)
  nsample <- 16
  formula <- ~treatment + (1|subject_ids)
  sim <- ALDEx3:::aldex.mem.sim(D = 6, days = 5, subjects = 8,
                       depth = 100000, sd_resid = 0.08)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)
  fits <- blmm_fit_pair(logW, formula, sim$meta)

  est_blmm <- apply(fits$blmm$estimate, c(1, 2), mean)
  est_lme4 <- apply(fits$lme4$estimate, c(1, 2), mean)
  se_blmm <- apply(fits$blmm$std.error, c(1, 2), mean)
  se_lme4 <- apply(fits$lme4$std.error, c(1, 2), mean)

  expect_lt(max(abs(est_lme4 - est_blmm)), 0.04)
  expect_lt(max(abs(se_lme4 - se_blmm) / (abs(se_lme4) + 0.01)), 0.02)
})

test_that("blmm C3: vignette subset discrepancies are materially reduced", {
  set.seed(42)
  dat <- blmm_vignette_subset()
  formula <- ~ treat * timec + (1 | participant_id)
  logW <- blmm_sample_logW(
    dat$Y, formula, dat$meta, 16, ALDEx3::clr.sm, list(gamma = 0)
  )
  fits <- blmm_fit_pair(logW, formula, dat$meta)

  est_blmm <- apply(fits$blmm$estimate, c(1, 2), mean)
  est_lme4 <- apply(fits$lme4$estimate, c(1, 2), mean)
  se_blmm <- apply(fits$blmm$std.error, c(1, 2), mean)
  se_lme4 <- apply(fits$lme4$std.error, c(1, 2), mean)

  expect_lt(max(abs(est_lme4 - est_blmm)), 0.04)
  expect_lt(max(abs(se_lme4 - se_blmm) / (abs(se_lme4) + 0.01)), 0.02)
})

test_that("blmm D: feature-level parallelization preserves results and ordering", {
  set.seed(2468)
  nsample <- 6
  formula <- ~treatment + (1|subject_ids)
  sim <- ALDEx3:::aldex.mem.sim(D = 4, days = 5, subjects = 6,
                       depth = 100000, sd_resid = 0.08)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)

  res_serial <- ALDEx3:::blmm(logW, formula, sim$meta, n.cores = 1L)
  res_parallel <- ALDEx3:::blmm(logW, formula, sim$meta, n.cores = 2L)

  expect_equal(res_serial$estimate, res_parallel$estimate, tolerance = 1e-10)
  expect_equal(res_serial$std.error, res_parallel$std.error, tolerance = 1e-10)
  expect_equal(res_serial$df, res_parallel$df, tolerance = 1e-10)
  expect_equal(res_serial$p.lower, res_parallel$p.lower, tolerance = 1e-10)
  expect_equal(res_serial$p.upper, res_parallel$p.upper, tolerance = 1e-10)
  expect_equal(res_serial$random.eff, res_parallel$random.eff, tolerance = 1e-10)
})

test_that("blmm E: parallel fallback aggregates warnings and remains exact", {
  set.seed(97531)
  nsample <- 4
  formula <- ~treatment + (1|subject_ids)
  sim <- ALDEx3:::aldex.mem.sim(D = 4, days = 4, subjects = 6,
                       depth = 100000, sd_resid = 0.08)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)
  attr(logW, "blmm.force_fail_features") <- c(1L, 2L)

  warns <- character()
  res_parallel <- withCallingHandlers(
    ALDEx3:::blmm(logW, formula, sim$meta, n.cores = 2L),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  res_exact <- ALDEx3:::sr.mem(logW, formula, sim$meta, n.cores = 1L,
                               method = "lme4", mem.args = list())

  expect_length(warns, 1)
  expect_match(warns[[1]], "2 feature\\(s\\) fell back to exact lme4")
  expect_equal(res_parallel$estimate[, 1:2, ], res_exact$estimate[, 1:2, ],
               tolerance = 1e-10)
  expect_equal(res_parallel$std.error[, 1:2, ], res_exact$std.error[, 1:2, ],
               tolerance = 1e-10)
  expect_equal(res_parallel$df[, 1:2, ], res_exact$df[, 1:2, ],
               tolerance = 1e-10)
  expect_equal(res_parallel$random.eff[, 1:2, ], res_exact$random.eff[, 1:2, ],
               tolerance = 1e-10)
})

test_that("blmm D/E/F: correlated random-slope output remains compatible", {
  set.seed(4321)
  nsample <- 1
  sim <- ALDEx3:::aldex.mem.sim(D = 4, days = 6, subjects = 25,
                       depth = 100000, random_slope = TRUE,
                       corr = 0.4, sd_resid = 0.05)
  logScale <- blmm_true_scale(sim, nsample)

  formula <- ~treatment + (1 + treatment | subject_ids)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)
  fits <- blmm_fit_pair(logW, formula, sim$meta)

  expect_equal(
    rownames(fits$blmm$random.eff),
    c("subject_ids.(Intercept)", "subject_ids.treatment",
      "subject_ids.(Intercept).treatment", "Residual")
  )
  expect_equal(rownames(fits$blmm$random.eff), rownames(fits$lme4$random.eff))
  expect_equal(fits$blmm$estimate, fits$lme4$estimate, tolerance = 1e-5)
  expect_equal(fits$blmm$random.eff, fits$lme4$random.eff, tolerance = 1e-3)
})

test_that("blmm G: approximation failures warn and fall back to exact lme4", {
  set.seed(77)
  nsample <- 3
  sim <- ALDEx3:::aldex.mem.sim(D = 2, days = 4, subjects = 6,
                       depth = 10000, sd_resid = 0.1)
  logScale <- blmm_true_scale(sim, nsample)

  local_mocked_bindings(
    blmm_make_adfun = function(...) stop("forced approximate failure"),
    .package = "ALDEx3"
  )

  expect_warning(
    res_blmm <- ALDEx3::aldex(sim$Y, ~treatment + (1|subject_ids),
                      data = sim$meta, method = "blmm",
                      nsample = nsample, scale = logScale, n.cores = 1,
                      return.pars = c("X", "estimate", "std.error",
                                      "p.val", "p.val.adj", "random.effects")),
    "fell back to exact lme4"
  )

  expect_s3_class(res_blmm, "aldex")
  expect_true(all(is.finite(res_blmm$estimate)))
  expect_true(all(is.finite(res_blmm$std.error)))
  expect_true(all(is.finite(res_blmm$random.effects)))
})
