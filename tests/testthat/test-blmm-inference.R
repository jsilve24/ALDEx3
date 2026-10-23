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

test_that("few-group inference fixes the residual-df discrepancy", {
  z <- inference_fixture(S = 1, noise = 0)
  old <- feature_fit(z, inference = "residual")
  a <- feature_fit(z, inference = "anchor")
  e <- ALDEx3:::blmm_exact_feature(z$Y, z$formula, z$data)
  expect_identical(a$fallback_draws, 0L)
  expect_equal(a$fit$estimate, old$fit$estimate)
  expect_equal(a$fit$std.error, old$fit$std.error)
  expect_gt(max(abs(old$fit$p.lower - e$p.lower)), 0.01)
  expect_equal(unname(a$fit$df), unname(e$df), tolerance = 0.001)
  expect_equal(unname(a$fit$p.lower), unname(e$p.lower), tolerance = 1e-4)
  expect_equal(unname(a$fit$p.upper), unname(e$p.upper), tolerance = 1e-4)
})

test_that("draw duplication and fixed-effect shifts preserve inference", {
  for (noise in c(0, 0.02)) {
    z <- inference_fixture(noise = noise)
    a <- feature_fit(z, inference = "anchor")
    expect_identical(a$fallback_draws, 0L)
    zz <- z
    zz$Y <- cbind(z$Y, z$Y)
    zz$logW <- array(zz$Y, c(nrow(zz$Y), 1, ncol(zz$Y)))
    dup <- feature_fit(zz, inference = "anchor")
    for (field in names(a$fit))
      expect_equal(dup$fit[[field]][,,seq_len(ncol(z$Y)),drop=FALSE], a$fit[[field]], tolerance = 1e-7)
    if (noise == 0) {
      zz$Y <- z$Y[,1,drop=FALSE]
      zz$logW <- array(zz$Y, c(nrow(zz$Y),1,1))
      single <- feature_fit(zz, inference = "anchor")
      for (field in names(a$fit))
        expect_equal(single$fit[[field]][,,1], a$fit[[field]][,,1], tolerance = 1e-7)
    }
    X <- model.matrix(~x, z$data)
    delta <- rbind(seq_len(ncol(z$Y))/10, -seq_len(ncol(z$Y))/5)
    zz <- z
    zz$Y <- z$Y + X %*% delta
    zz$logW <- array(zz$Y, dim(z$logW))
    shifted <- feature_fit(zz, inference = "anchor")
    expect_identical(shifted$fallback_draws, 0L)
    expect_equal(shifted$fit$estimate, a$fit$estimate + array(delta, dim(a$fit$estimate)), tolerance = 1e-7)
    for (field in c("std.error", "df", "random.eff"))
      expect_equal(shifted$fit[[field]], a$fit[[field]], tolerance = 1e-7)
  }
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

test_that("inference-only failure rescues a draw consistently", {
  z <- inference_fixture()
  baseline <- feature_fit(z,inference="anchor")
  original <- ALDEx3:::blmm_anchor_df
  calls <- 0L
  testthat::local_mocked_bindings(blmm_anchor_df = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("invalid anchor inference sensitivity")
    original(...)
  }, .package = "ALDEx3")
  a <- feature_fit(z,inference="anchor")
  e <- ALDEx3:::blmm_exact_feature(z$Y,z$formula,z$data)
  expect_identical(a$fallback_draws,1L)
  expect_false(a$full_fallback)
  expect_match(a$fallback_message,"inference sensitivity")
  for (field in names(e)) expect_equal(unname(a$fit[[field]][,,2]),unname(e[[field]][,,2]))
  for (field in names(a$fit))
    expect_equal(a$fit[[field]][,,-2],baseline$fit[[field]][,,-2])
})

test_that("shared-anchor inference preserves parallel ordering and aggregation", {
  z <- inference_fixture()
  logW <- abind::abind(z$logW,z$logW+0.2,along=2)
  a <- ALDEx3:::blmm(logW,z$formula,z$data,inference="anchor")
  b <- ALDEx3:::blmm(logW,z$formula,z$data,n.cores=2,inference="anchor")
  expect_equal(a,b)
  e <- ALDEx3:::sr.mem(logW,z$formula,z$data,n.cores=1,method="lme4",mem.args=list())
  agg <- function(x) ALDEx3:::aldex.pvals(x$p.lower,x$p.upper,"BH")
  expect_equal(agg(a),agg(e),tolerance=0.002)
  attr(logW,"blmm.force_fail_features") <- 1:2
  expect_warning(b <- ALDEx3:::blmm(logW,z$formula,z$data,inference="anchor"),"fell back")
  expect_equal(b,e)
})

test_that("both wrappers calculate small upper tails directly", {
  z <- inference_fixture(S=1,noise=0,design="within",effect=30)
  a <- feature_fit(z,inference="anchor")
  e <- ALDEx3:::blmm_exact_feature(z$Y,z$formula,z$data)
  expect_identical(a$fallback_draws,0L)
  expect_gt(a$fit$p.upper[2],0)
  expect_lt(a$fit$p.upper[2],.Machine$double.eps)
  expect_gt(e$p.upper[2],0)
  expect_lt(abs(log(a$fit$p.upper[2])-log(e$p.upper[2])),0.002)
})

test_that("natural slope update failures preserve selective exact rescue", {
  z <- inference_fixture(seed=105002,groups=12,design="slope",noise=.02,S=8,effect=0)
  a <- feature_fit(z,inference="anchor")
  expect_false(a$full_fallback)
  expect_gt(a$fallback_draws,0)
  expect_lt(a$fallback_draws,ncol(z$Y))
  expect_match(a$fallback_message,"worsened the profiled REML")
  # Compare rescue only: a reference fit for one retained draw in this
  # stress fixture has a convergence warning, so it is not an accuracy oracle.
  e <- ALDEx3:::blmm_exact_feature(z$Y[,a$rejected_draws,drop=FALSE],z$formula,z$data)
  for (field in names(e))
    expect_equal(unname(a$fit[[field]][,,a$rejected_draws,drop=FALSE]),
                 unname(e[[field]]))
  retained <- setdiff(seq_len(ncol(z$Y)),a$rejected_draws)
  expect_true(all(is.finite(a$fit$df[,,retained])))
  expect_true(all(a$fit$df[,,retained]>0))
  expect_lt(max(a$fit$df[,,retained]),nrow(z$Y)-2)
})
