# Tests for BLMM engine
#
# Tolerance rationale (set prospectively, must NOT be relaxed after observing results):
#   A: With nsample=1, anchor == draw optimum => lme4 agreement within 0.05 (log2).
#      0.05 << typical posterior SD of 0.1-0.5; does not affect inference.
#   B: Identical draws => zero variance across draws (< 1e-8).
#   C: Posterior mean error < 0.05; posterior SD relative error < 0.25.
#      Rankings across features are preserved; inference is not changed.
#   D: Output structure must be identical to lme4 output.
#   E: Near-singular must not crash.
#
# Do not move these thresholds after seeing failures. Surface the failure.

test_that("blmm: expected input errors", {
  set.seed(4985)
  nsample <- 4
  sim <- aldex.mem.sim(5, 4, 4, 10000, 1, FALSE, 0, 0, 0.1)
  true.S   <- log2(colSums(sim$W))
  logScale <- matrix(rep(true.S, nsample), nrow = length(true.S))

  expect_error(
    aldex(sim$Y, sim$meta, data = sim$meta, method = "blmm",
          nsample = nsample, scale = logScale, n.cores = 1),
    "X should be a mixed effects formula"
  )
  expect_error(
    aldex(sim$Y, ~treatment + (1|subject_ids), data = NULL,
          method = "blmm", nsample = nsample, scale = logScale, n.cores = 1),
    "data should not be null"
  )
})

test_that("blmm A: single-draw agreement with lme4 (tolerance 0.05)", {
  set.seed(42)
  nsample <- 1
  sim <- aldex.mem.sim(D = 5, days = 5, subjects = 6,
                       depth = 100000, sd_resid = 0.05)
  true.S   <- log2(colSums(sim$W))
  logScale <- matrix(true.S, nrow = length(true.S), ncol = nsample)

  res_blmm <- aldex(sim$Y, ~treatment + (1|subject_ids),
                    data = sim$meta, method = "blmm",
                    nsample = nsample, scale = logScale, n.cores = 1)
  res_lme4 <- aldex(sim$Y, ~treatment + (1|subject_ids),
                    data = sim$meta, method = "lme4",
                    nsample = nsample, scale = logScale, n.cores = 1)

  expect_equal(
    as.vector(res_blmm$estimate[2, , ]),
    as.vector(res_lme4$estimate[2, , ]),
    tolerance = 0.05,
    label = "treatment estimates agree within 0.05 for single draw"
  )
})

test_that("blmm B: identical-draw consistency (var < 1e-8)", {
  # Property: blmm() is a deterministic function of its logW input.
  # When all S logW slices are identical, scores = 0, phi_tilde = phi_bar
  # for all s, and all draws must produce exactly the same estimates.
  # We bypass aldex.sampler (which adds Dirichlet noise) and call blmm()
  # directly with a replicated single-draw logW array.
  set.seed(123)
  nsample <- 20
  sim <- aldex.mem.sim(D = 5, days = 4, subjects = 5,
                       depth = 100000, sd_resid = 0.05)
  true.S <- log2(colSums(sim$W))

  # Build one deterministic logW slice (log-proportion + log-scale)
  logp   <- log2(sweep(sim$Y + 0.5, 2, colSums(sim$Y + 0.5), "/"))  # D x N
  logW1  <- logp + matrix(true.S, nrow(sim$Y), ncol(sim$Y),
                           byrow = TRUE)  # D x N
  # Replicate to get N x D x nsample with identical slices
  logW <- aperm(array(logW1, c(nrow(sim$Y), ncol(sim$Y), nsample)),
                c(2, 1, 3))  # N x D x nsample

  res_blmm <- ALDEx3:::blmm(logW, ~treatment + (1|subject_ids),
                             sim$meta, n.cores = 1L)

  est_var <- apply(res_blmm$estimate, c(1, 2), var)
  expect_true(all(est_var < 1e-8),
              label = "variance across identical draws < 1e-8")
})

test_that("blmm C: posterior mean/SD close to lme4 (mean tol 0.05, SD tol 0.25)", {
  set.seed(999)
  nsample <- 50
  sim <- aldex.mem.sim(D = 8, days = 4, subjects = 8,
                       depth = 100000, sd_resid = 0.1)
  true.S   <- log2(colSums(sim$W))
  logScale <- matrix(rep(true.S, nsample), nrow = length(true.S))

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
  expect_equal(as.vector(blmm_mean), as.vector(lme4_mean), tolerance = 0.05,
               label = "posterior mean treatment within 0.05")

  blmm_sd <- apply(res_blmm$estimate[2, , , drop = FALSE], c(1, 2), sd)
  lme4_sd  <- apply(res_lme4$estimate[2, , , drop = FALSE], c(1, 2), sd)
  rel_err  <- abs(blmm_sd - lme4_sd) / (abs(lme4_sd) + 0.01)
  expect_true(all(as.vector(rel_err) < 0.25),
              label = "posterior SD relative error < 0.25")
})

test_that("blmm D: output structure compatible with summary() and accessors", {
  set.seed(77)
  nsample <- 10
  sim <- aldex.mem.sim(D = 5, days = 3, subjects = 4, depth = 10000)
  true.S   <- log2(colSums(sim$W))
  logScale <- matrix(rep(true.S, nsample), nrow = length(true.S))

  res <- aldex(sim$Y, ~treatment + (1|subject_ids),
               data = sim$meta, method = "blmm",
               nsample = nsample, scale = logScale, n.cores = 1,
               return.pars = c("X", "estimate", "std.error",
                               "p.val", "p.val.adj", "logComp", "logScale"))

  expect_s3_class(res, "aldex")
  expect_equal(dim(res$estimate), c(2L, 5L, nsample))
  expect_equal(colnames(res$estimate), rownames(sim$Y))
  expect_equal(rownames(res$estimate), c("(Intercept)", "treatment"))

  s <- summary(res)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("estimate", "std.error", "p.val.adj") %in% names(s)))
})

test_that("blmm E: near-singular / tiny sample does not crash", {
  set.seed(888)
  nsample <- 5
  sim <- aldex.mem.sim(D = 3, days = 2, subjects = 2, depth = 1000)
  true.S   <- log2(colSums(sim$W))
  logScale <- matrix(rep(true.S, nsample), nrow = length(true.S))

  expect_no_error({
    res <- aldex(sim$Y, ~treatment + (1|subject_ids),
                 data = sim$meta, method = "blmm",
                 nsample = nsample, scale = logScale, n.cores = 1)
  })
  expect_s3_class(res, "aldex")
  expect_true(is.array(res$estimate))
})
