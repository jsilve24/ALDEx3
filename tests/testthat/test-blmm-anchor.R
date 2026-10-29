# Shared-anchor inference (public test = "t.satterthwaite"): agreement with
# exact lme4/lmerTest, the df derivation, and fallback warnings.
# Tolerances account for lmerTest's numerical Hessian and optimiser stopping:
# 0.2% df/SE, 0.002 absolute log-tail error in regular single-draw examples.
test_that("anchor inference agrees on genuinely accepted single draws", {
  for (design in c("between", "within", "unbalanced", "slope")) {
    z <- inference_fixture(seed = 731, groups = 12, design = design, S = 1, noise = 0)
    a <- feature_fit(z, inference = "anchor")
    e <- ALDEx3:::blmm_exact_feature(z$Y, z$formula, z$data)
    expect_identical(a$fallback_draws, 0L, info = design)
    for (field in c("estimate", "std.error", "df"))
      expect_equal(unname(a$fit[[field]]), unname(e[[field]]), tolerance = 0.002, info = paste(design, field))
    for (tail in c("p.lower", "p.upper"))
      expect_lt(max(abs(log(a$fit[[tail]]) - log(e[[tail]]))), 0.002)
  }
})

test_that("anchor inference matches exact lme4 on a single simulated draw", {
  set.seed(42)
  nsample <- 1
  formula <- ~treatment + (1|subject_ids)
  sim <- aldex.mem.sim(D = 5, days = 5, subjects = 6,
                       depth = 100000, sd_resid = 0.05)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)
  fits <- blmm_fit_pair(logW, formula, sim$meta, inference = "anchor")

  expect_equal(fits$blmm$estimate, fits$lme4$estimate, tolerance = 1e-5)
  expect_equal(fits$blmm$std.error, fits$lme4$std.error, tolerance = 1e-5)
})

test_that("profiled information agrees with full observed information", {
  skip_if_not_installed("numDeriv")
  # Unequal sizes make M nondiagonal, exercising the LX orientation.
  z <- inference_fixture(S = 1, noise = 0, design = "unbalanced")
  p <- ALDEx3:::blmm_parse_formula(z$formula, z$data)
  obj <- ALDEx3:::blmm_make_adfun(p$X,z$Y,p$basis,p$is_log,p$phi_init)
  anchor <- ALDEx3:::blmm_fit_anchor(obj,p$phi_init)
  obj$fn(anchor$phi_bar)
  report <- obj$report()
  d <- nrow(p$X)-ncol(p$X)
  C <- cbind(diag(ncol(p$X)), c(1, -0.7))
  sens <- ALDEx3:::blmm_scores_batch(obj,anchor$phi_bar,report$PWRSS_s,
    nrow(p$X),ncol(p$X), contrasts=C,inference=TRUE)
  df <- ALDEx3:::blmm_anchor_df(sens$u[,1],sens$b,chol(anchor$H_d),d)
  full <- function(par) {
    obj$fn(par[-length(par)])
    r <- obj$report()
    eta <- tail(par,1)
    0.5*(r$ldL2+r$ldRX2+d*eta+r$PWRSS_s*exp(-eta))
  }
  par <- c(anchor$phi_bar, log(report$PWRSS_s/d))
  H <- numDeriv::hessian(full,par)
  for (j in seq_len(ncol(C))) {
    logv <- function(par) {
      obj$fn(par[-length(par)])
      tail(par,1)+log(sum(forwardsolve(obj$report()$LX,C[,j])^2))
    }
    gradient <- numDeriv::grad(logv,par)
    full_df <- 2 / as.numeric(crossprod(gradient,solve(H,gradient)))
    expect_equal(df[j],full_df,tolerance=1e-5)
  }
  expect_equal(ALDEx3:::blmm_anchor_df(0,matrix(0,1,3),matrix(1),d),rep(d,3))
  expect_error(ALDEx3:::blmm_anchor_df(NA,matrix(0),matrix(1),d),"sensitivity")
})

test_that("an indefinite anchor Hessian causes exact feature fallback", {
  set.seed(8642)
  nsample <- 3
  formula <- ~treatment + (1|subject_ids)
  sim <- aldex.mem.sim(D = 1, days = 4, subjects = 6,
                       depth = 100000, sd_resid = 0.08)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)
  exact <- ALDEx3:::sr.mem(logW, formula, sim$meta, n.cores = 1L,
                           method = "lme4", mem.args = list())

  make_adfun <- ALDEx3:::blmm_make_adfun
  local_mocked_bindings(
    blmm_make_adfun = function(...) {
      obj <- make_adfun(...)
      obj$he <- function(...) matrix(-1, nrow = 1L, ncol = 1L)
      obj
    },
    .package = "ALDEx3"
  )

  expect_warning(
    fit <- ALDEx3:::blmm(logW, formula, sim$meta, n.cores = 1L,
                          inference = "anchor"),
    "1 feature\\(s\\) were refit with lmerTest"
  )
  expect_equal(fit$estimate, exact$estimate, tolerance = 1e-10)
  expect_equal(fit$std.error, exact$std.error, tolerance = 1e-10)
  for (field in c("df", "p.lower", "p.upper"))
    expect_equal(fit[[field]], exact[[field]], tolerance = 1e-10)
  expect_equal(fit$random.eff, exact$random.eff, tolerance = 1e-10)
})

test_that("a materially worse Newton update falls back by draw", {
  set.seed(7531)
  nsample <- 4
  formula <- ~treatment + (1|subject_ids)
  sim <- aldex.mem.sim(D = 1, days = 5, subjects = 8,
                       depth = 100000, sd_resid = 0.08)
  logScale <- blmm_true_scale(sim, nsample)
  logW <- blmm_sample_logW(sim$Y, formula, sim$meta, nsample, logScale)

  baseline <- ALDEx3:::blmm(logW, formula, sim$meta, n.cores = 1L,
                          inference = "anchor")
  exact <- ALDEx3:::sr.mem(logW, formula, sim$meta, n.cores = 1L,
                           method = "lme4", mem.args = list())
  phi_updates <- ALDEx3:::blmm_phi_updates
  local_mocked_bindings(
    blmm_phi_updates = function(phi_bar, H_d, scores, ...) {
      updates <- phi_updates(phi_bar, H_d, scores, ...)
      updates[, 1] <- phi_bar + 10
      updates
    },
    .package = "ALDEx3"
  )

  expect_warning(
    fit <- ALDEx3:::blmm(logW, formula, sim$meta, n.cores = 1L,
                          inference = "anchor"),
    "1 draw\\(s\\) across 1 feature\\(s\\) were refit with lmerTest"
  )
  expect_equal(fit$estimate[, , 1], exact$estimate[, , 1], tolerance = 1e-10)
  expect_equal(fit$std.error[, , 1], exact$std.error[, , 1], tolerance = 1e-10)
  for (field in c("df", "p.lower", "p.upper")) {
    expect_equal(fit[[field]][,,1], exact[[field]][,,1], tolerance = 1e-10)
    expect_equal(fit[[field]][,,-1], baseline[[field]][,,-1], tolerance = 1e-10)
  }
  expect_equal(fit$random.eff[, , 1], exact$random.eff[, , 1], tolerance = 1e-10)
  expect_equal(fit$estimate[, , -1], baseline$estimate[, , -1], tolerance = 1e-10)
  expect_equal(fit$std.error[, , -1], baseline$std.error[, , -1], tolerance = 1e-10)
  expect_equal(fit$random.eff[, , -1], baseline$random.eff[, , -1], tolerance = 1e-10)
})
