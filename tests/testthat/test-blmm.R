# Engine-level BLMM tests under the default adaptive inference. Selector-
# specific checks live in test-blmm-adaptive.R and test-blmm-anchor.R.
#
# Tolerance rationale (set prospectively):
#   B. Identical draws should produce numerically identical outputs.
#   C. On real data, posterior means and SDs should stay well inside the
#      posterior uncertainty of exact lme4.
#   D. Random-slope and correlated-random-effect output should remain close
#      enough to exact lme4 to preserve interpretation of the covariance terms.
#   E. If the approximation fails, BLMM must warn and fall back to exact lme4
#      rather than silently changing the model.

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
  sim <- aldex.mem.sim(5, 4, 4, 10000, 1, FALSE, 0, 0, 0.1)
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

test_that("blmm: raw anchor Hessians are validated without regularisation", {
  H_raw <- matrix(c(2, 0.9, 1.1, 2), nrow = 2)
  H_expected <- matrix(c(2, 1, 1, 2), nrow = 2)

  expect_equal(ALDEx3:::blmm_validate_H(H_raw), H_expected)
  expect_null(attr(ALDEx3:::blmm_validate_H(H_raw), "ridge_delta"))
  expect_equal(
    ALDEx3:::blmm_validate_H(matrix(1e-12, nrow = 1)),
    matrix(1e-12, nrow = 1)
  )
  expect_error(
    ALDEx3:::blmm_validate_H(matrix(c(1, NA, NA, 1), nrow = 2)),
    "non-finite"
  )
  expect_error(
    ALDEx3:::blmm_validate_H(diag(c(1, -1e-8))),
    "not positive definite"
  )
  expect_error(
    ALDEx3:::blmm_validate_H(diag(c(1, 1e-7))),
    "ill-conditioned"
  )
})

test_that("blmm B: identical-draw consistency", {
  set.seed(123)
  nsample <- 20
  sim <- aldex.mem.sim(D = 5, days = 4, subjects = 5,
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

test_that("blmm C3: vignette subset discrepancies are materially reduced", {
  skip_on_cran()
  set.seed(42)
  dat <- blmm_vignette_subset()
  formula <- ~ treat * timec + (1 | participant_id)
  logW <- blmm_sample_logW(
    dat$Y, formula, dat$meta, 16, ALDEx3::clr.sm, list(gamma = 0)
  )
  # An occasional draw is rescued with an exact refit; that warning is expected.
  fits <- suppressWarnings(blmm_fit_pair(logW, formula, dat$meta))

  est_blmm <- apply(fits$blmm$estimate, c(1, 2), mean)
  est_lme4 <- apply(fits$lme4$estimate, c(1, 2), mean)
  se_blmm <- apply(fits$blmm$std.error, c(1, 2), mean)
  se_lme4 <- apply(fits$lme4$std.error, c(1, 2), mean)

  expect_lt(max(abs(est_lme4 - est_blmm)), 0.04)
  expect_lt(max(abs(se_lme4 - se_blmm) / (abs(se_lme4) + 0.01)), 0.02)
})

test_that("blmm E: parallel fallback aggregates warnings and remains exact", {
  set.seed(97531)
  nsample <- 4
  formula <- ~treatment + (1|subject_ids)
  sim <- aldex.mem.sim(D = 4, days = 4, subjects = 6,
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
  expect_match(warns[[1]], "2 feature\\(s\\) were refit with lmerTest")
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
  sim <- aldex.mem.sim(D = 4, days = 6, subjects = 25,
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
  sim <- aldex.mem.sim(D = 2, days = 4, subjects = 6,
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
    "refit with lmerTest"
  )

  expect_s3_class(res_blmm, "aldex")
  expect_true(all(is.finite(res_blmm$estimate)))
  expect_true(all(is.finite(res_blmm$std.error)))
  expect_true(all(is.finite(res_blmm$random.effects)))
})
