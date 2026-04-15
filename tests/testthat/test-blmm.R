# Tests for the approximate BLMM mixed-effects engine.
#
# Tolerance rationale (set prospectively):
#   A. With nsample = 1, BLMM should closely match exact lme4 because there is
#      no posterior averaging across draws. A 0.05 log2-scale tolerance is well
#      below the posterior SD routinely used for inference.
#   B. Identical draws should produce numerically identical outputs because the
#      draw-specific score corrections are zero.
#   C. Across repeated draws, posterior mean error should remain below 0.05 and
#      posterior SD relative error below 0.25, so the approximation remains
#      materially smaller than the posterior uncertainty being reported.
#   D. Random-slope and correlated-random-effect output should remain close
#      enough to exact lme4 to preserve interpretation of the covariance terms.
#   E. If the approximation fails, BLMM must warn and fall back to exact lme4
#      rather than silently changing the model.

blmm_true_scale <- function(sim, nsample) {
  matrix(rep(log2(colSums(sim$W)), nsample), nrow = ncol(sim$Y))
}

test_that("blmm: expected input errors", {
  set.seed(4985)
  nsample <- 4
  sim <- aldex.mem.sim(5, 4, 4, 10000, 1, FALSE, 0, 0, 0.1)
  logScale <- blmm_true_scale(sim, nsample)

  expect_error(
    aldex(sim$Y, sim$meta, data = sim$meta, method = "blmm",
          nsample = nsample, scale = logScale, n.cores = 1),
    "X should be a mixed effects formula"
  )
  expect_error(
    aldex(sim$Y, ~treatment + (1|subject_ids), data = NULL,
          method = "blmm", nsample = nsample, scale = logScale, n.cores = 1),
    "data should not be null if method is a mixed-effects engine"
  )
})

test_that("blmm A: single-draw agreement with exact lme4", {
  set.seed(42)
  nsample <- 1
  sim <- aldex.mem.sim(D = 5, days = 5, subjects = 6,
                       depth = 100000, sd_resid = 0.05)
  logScale <- blmm_true_scale(sim, nsample)

  res_blmm <- aldex(sim$Y, ~treatment + (1|subject_ids),
                    data = sim$meta, method = "blmm",
                    nsample = nsample, scale = logScale, n.cores = 1)
  res_lme4 <- aldex(sim$Y, ~treatment + (1|subject_ids),
                    data = sim$meta, method = "lme4",
                    nsample = nsample, scale = logScale, n.cores = 1)

  expect_equal(res_blmm$estimate, res_lme4$estimate, tolerance = 0.05)
  expect_equal(res_blmm$std.error, res_lme4$std.error, tolerance = 0.05)
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

test_that("blmm C: posterior mean and SD remain close to exact lme4", {
  set.seed(999)
  nsample <- 40
  sim <- aldex.mem.sim(D = 8, days = 4, subjects = 8,
                       depth = 100000, sd_resid = 0.1)
  logScale <- blmm_true_scale(sim, nsample)

  res_blmm <- aldex(sim$Y, ~treatment + (1|subject_ids),
                    data = sim$meta, method = "blmm",
                    nsample = nsample, scale = logScale, n.cores = 1,
                    return.pars = c("X", "estimate", "std.error",
                                    "p.val", "p.val.adj"))
  res_lme4 <- aldex(sim$Y, ~treatment + (1|subject_ids),
                    data = sim$meta, method = "lme4",
                    nsample = nsample, scale = logScale, n.cores = 1,
                    return.pars = c("X", "estimate", "std.error",
                                    "p.val", "p.val.adj"))

  blmm_mean <- apply(res_blmm$estimate[2, , , drop = FALSE], c(1, 2), mean)
  lme4_mean <- apply(res_lme4$estimate[2, , , drop = FALSE], c(1, 2), mean)
  expect_equal(as.vector(blmm_mean), as.vector(lme4_mean), tolerance = 0.05)

  blmm_sd <- apply(res_blmm$estimate[2, , , drop = FALSE], c(1, 2), sd)
  lme4_sd <- apply(res_lme4$estimate[2, , , drop = FALSE], c(1, 2), sd)
  rel_err <- abs(blmm_sd - lme4_sd) / (abs(lme4_sd) + 0.01)
  expect_true(all(as.vector(rel_err) < 0.25))
})

test_that("blmm D/E/F: correlated random-slope output remains compatible", {
  set.seed(4321)
  nsample <- 1
  sim <- aldex.mem.sim(D = 4, days = 6, subjects = 25,
                       depth = 100000, random_slope = TRUE,
                       corr = 0.4, sd_resid = 0.05)
  logScale <- blmm_true_scale(sim, nsample)

  formula <- ~treatment + (1 + treatment | subject_ids)
  res_blmm <- aldex(sim$Y, formula, data = sim$meta, method = "blmm",
                    nsample = nsample, scale = logScale, n.cores = 1,
                    return.pars = c("X", "estimate", "std.error",
                                    "p.val", "p.val.adj", "random.effects"))
  res_lme4 <- aldex(sim$Y, formula, data = sim$meta, method = "lme4",
                    nsample = nsample, scale = logScale, n.cores = 1,
                    return.pars = c("X", "estimate", "std.error",
                                    "p.val", "p.val.adj", "random.effects"))

  expect_equal(
    rownames(res_blmm$random.effects),
    c("subject_ids.(Intercept)", "subject_ids.treatment",
      "subject_ids.(Intercept).treatment", "Residual")
  )
  expect_equal(rownames(res_blmm$random.effects), rownames(res_lme4$random.effects))
  expect_equal(res_blmm$estimate, res_lme4$estimate, tolerance = 0.1)
  expect_equal(res_blmm$random.effects, res_lme4$random.effects, tolerance = 0.2)

  s <- summary(res_blmm)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("estimate", "std.error", "p.val.adj") %in% names(s)))
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
    res_blmm <- aldex(sim$Y, ~treatment + (1|subject_ids),
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
