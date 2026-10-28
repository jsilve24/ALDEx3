test_that("fflm contrast estimates and covariance projections are correct", {
  set.seed(401)
  n <- 15
  X <- cbind(`(Intercept)`=1, age=rnorm(n), genotype=rep(c(0, 1, 0), length.out=n))
  Y <- array(rnorm(n * 2 * 2), c(n, 2, 2))
  C <- rbind(
    age = c(0, 1, 0),
    genotype = c(0, 0, 1),
    age_by_genotype = c(0, 1, -1),
    dependent = c(0, 2, -2)
  )
  colnames(C) <- colnames(X)

  for (selector in c("t", "t.HC0", "t.HC3")) {
    fit <- fflm(Y, X, test=selector, contrasts=C)
    model <- lm(Y[, 1, 1] ~ X - 1)
    V <- if (selector == "t") vcov(model) else
      sandwich::vcovHC(model, type=sub("t\\.", "", selector))
    expect_equal(fit$contrast.estimate[, 1, 1], drop(C %*% coef(model)))
    expected.se <- sqrt(diag(C %*% V %*% t(C)))
    expect_equal(fit$contrast.std.error[, 1, 1], expected.se, tolerance=1e-12)
    expect_equal(unname(fit$contrast.df[, 1, 1]), rep(n - ncol(X), nrow(C)))
  }
})

test_that("aldex resolves method selectors before sampling", {
  Y <- matrix(1:32, 4, 8)
  X <- rbind(`(Intercept)`=1, group=rep(c(0, 1), each=4))
  testthat::local_mocked_bindings(aldex.sampler=function(...) stop("sampled"),
                                  .package="ALDEx3")
  expect_error(aldex(Y, X, method="LM", scale=clr.sm), "method must")
  expect_error(aldex(Y, X, method="unknown", scale=clr.sm), "method must")
  for (engine in c("lm", "lme4", "nlme")) {
    for (selector in c("t.satterthwaite", "t.residual")) {
      expect_error(aldex(Y, X, method=engine, test=selector, scale=clr.sm),
                   "allowed selectors")
    }
  }
  expect_error(aldex(Y, X, method="lm", test="T", scale=clr.sm), "allowed selectors")
  expect_error(aldex(Y, X, method="blmm", test="t.HC3", scale=clr.sm), "allowed selectors")
})

test_that("BLMM default selector uses anchor degrees of freedom", {
  data <- data.frame(g=factor(rep(1:8, each=5)),
                     x=rep(rep(c(-.5, .5), each=5), 4))
  N <- nrow(data)
  set.seed(44)
  subject <- rnorm(8, sd=.35)[data$g]
  response <- .4 + .7 * data$x + subject + rnorm(N, sd=.08)
  response2 <- .2 + .5 * data$x + subject + rnorm(N, sd=.07)
  saved <- array(c(response, response2), c(2, N, 1))
  Y <- matrix(10, 2, N)
  scale <- matrix(0, N, 1)
  C <- matrix(c(0, 1), 1, dimnames=list("x", c("(Intercept)", "x")))
  testthat::local_mocked_bindings(aldex.sampler=function(...) list(logW=saved),
                                  .package="ALDEx3")
  anchor <- aldex(Y, ~x + (1 | g), data, method="blmm", n.cores=1,
                  nsample=1, scale=scale, contrasts=C,
                  return.pars=c("estimate", "df", "p.val"))
  anchor.explicit <- aldex(Y, ~x + (1 | g), data, method="blmm",
                           test="t.satterthwaite", n.cores=1,
                           nsample=1, scale=scale, contrasts=C,
                           return.pars=c("estimate", "df", "p.val"))
  residual <- aldex(Y, ~x + (1 | g), data, method="blmm", test="t.residual",
                    n.cores=1, nsample=1, scale=scale, contrasts=C,
                    return.pars=c("estimate", "df", "p.val"))
  expect_identical(anchor$test, "t.satterthwaite")
  expect_identical(anchor.explicit$test, "t.satterthwaite")
  expect_identical(residual$test, "t.residual")
  expect_equal(anchor$estimate, anchor.explicit$estimate)
  expect_equal(anchor$df, anchor.explicit$df)
  expect_true(all(vapply(c("estimate", "df", "p.val"),
                         function(field) !is.null(anchor[[field]]) &&
                           all(is.finite(anchor[[field]])), logical(1))))
  expect_lt(max(anchor$contrasts$df), 9)
  expect_equal(unname(residual$contrasts$df), array(38, c(1, 2, 1)))
  parsed <- ALDEx3:::blmm_parse_formula(~x + (1 | g), data)
  feature <- ALDEx3:::blmm_fit_feature(
    d=1, logW=aperm(saved, c(2, 1, 3)), X=parsed$X, basis=parsed$basis,
    is_log=parsed$is_log, phi_init=parsed$phi_init, lower=parsed$lower,
    reTrms=parsed$reTrms, formula=~x + (1 | g), data=data,
    inference="anchor", contrasts=C)
  expect_false(feature$fallback)
})

test_that("all method defaults and explicit native selectors are recorded", {
  data <- data.frame(g=factor(rep(1:5, each=2)), x=rep(c(-.5, .5), 5))
  Y <- rbind(a=10:19, b=20:29)
  C <- matrix(c(0, 1), 1, dimnames=list("x", c("(Intercept)", "x")))
  C.lm <- matrix(c(0, 1), 1, dimnames=list("x", c("intercept", "x")))
  scale <- matrix(0, 10, 1)
  lm.default <- aldex(Y, rbind(intercept=1, x=data$x), nsample=1, scale=scale,
                      contrasts=C.lm, return.pars="df")
  lm.explicit <- aldex(Y, rbind(intercept=1, x=data$x), nsample=1, scale=scale,
                       test="t.HC3", contrasts=C.lm, return.pars="df")
  lme4.default <- aldex(Y, ~x + (1 | g), data, method="lme4", n.cores=1,
                         nsample=1, scale=scale, contrasts=C, return.pars="df")
  lme4.explicit <- aldex(Y, ~x + (1 | g), data, method="lme4", test="t",
                          n.cores=1, nsample=1, scale=scale, contrasts=C,
                          return.pars="df")
  nlme.default <- aldex(Y, ~x, data, method="nlme", random=~1 | g, n.cores=1,
                         nsample=1, scale=scale, contrasts=C, return.pars="df")
  nlme.explicit <- aldex(Y, ~x, data, method="nlme", test="t", random=~1 | g,
                          n.cores=1, nsample=1, scale=scale, contrasts=C,
                          return.pars="df")
  expect_identical(c(lm.default$test, lm.explicit$test), c("t.HC3", "t.HC3"))
  expect_identical(c(lme4.default$test, lme4.explicit$test), c("t", "t"))
  expect_identical(c(nlme.default$test, nlme.explicit$test), c("t", "t"))
})

test_that("aldex returns named contrast results under return controls", {
  set.seed(402)
  Y <- matrix(sample.int(20, 32, replace=TRUE), 4, 8)
  X <- rbind(`(Intercept)`=1, age=seq_len(8), genotype=rep(c(0, 1), each=4))
  C <- rbind(age=c(0, 1, 0), genotype=c(0, 0, 1),
             age_minus_genotype=c(0, 1, -1))
  colnames(C) <- rownames(X)
  fit <- aldex(Y, X, nsample=3, scale=clr.sm, contrasts=C,
               return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  expect_identical(fit$method, "lm")
  expect_identical(fit$test, "t.HC3")
  expect_identical(dim(fit$contrasts$estimate), c(3L, 4L, 3L))
  expect_identical(rownames(fit$contrasts$estimate), rownames(C))
  expect_equal(names(fit$contrasts), c("matrix", "estimate", "std.error", "df", "p.val", "p.val.adj"))
  expect_equal(summary(fit, component="contrasts")$parameter,
               rep(rownames(C), times=nrow(Y)))

  matrix.only <- aldex(Y, X, nsample=2, scale=clr.sm, contrasts=C,
                       return.pars="X")
  expect_equal(names(matrix.only$contrasts), "matrix")
})

test_that("contrast validation reports malformed and nonestimable contrasts", {
  X <- cbind(`(Intercept)`=1, x=1:5, duplicate=1:5)
  estimable <- rbind(shared_slope=c(0, 0.5, 0.5))
  colnames(estimable) <- colnames(X)
  expect_equal(ALDEx3:::.aldex_validate_contrasts(estimable, X), estimable)

  bad <- rbind(alias=c(0, 1, 0))
  colnames(bad) <- colnames(X)
  expect_error(ALDEx3:::.aldex_validate_contrasts(bad, X), "nonestimable")

  X2 <- X[, 1:2, drop=FALSE]
  C <- matrix(c(0, 0), 1, dimnames=list("zero", colnames(X2)))
  expect_error(ALDEx3:::.aldex_validate_contrasts(C, X2), "nonzero")
  C <- matrix(c(0, 1), 1, dimnames=list("ok", rev(colnames(X2))))
  expect_equal(colnames(ALDEx3:::.aldex_validate_contrasts(C, X2)), colnames(X2))
})

test_that("contrast validation handles tall full-rank designs", {
  observation <- seq_len(5000)
  X <- cbind(intercept=1, trend=observation, group=rep(c(-1, 1), 2500))
  contrasts <- rbind(trend=c(0, 1, 0), group=c(0, 0, 1))
  colnames(contrasts) <- rev(colnames(X))
  expected <- contrasts[, colnames(X), drop=FALSE]

  expect_equal(qr(X)$rank, ncol(X))
  validate <- ALDEx3:::.aldex_validate_contrasts
  validation_environment <- new.env(parent=environment(validate))
  validation_environment$qr <- function(design, ...) {
    expect_lte(ncol(design), nrow(design))
    base::qr(design, ...)
  }
  environment(validate) <- validation_environment
  expect_equal(validate(contrasts, X), expected)
})

test_that("contrast validation covers dimensions, names, finite values, and dependence", {
  X <- cbind(intercept=1, x=1:6, z=rep(c(0, 1), 3))
  named <- function(value=c(0, 1, 0)) {
    matrix(value, 1, dimnames=list("c", colnames(X)))
  }
  expect_error(ALDEx3:::.aldex_validate_contrasts(c(0, 1, 0), X), "numeric matrix")
  expect_error(ALDEx3:::.aldex_validate_contrasts(matrix(numeric(), 0, 3), X),
               "at least one")
  unnamed.columns <- matrix(1, 1, 2, dimnames=list("c", NULL))
  expect_error(ALDEx3:::.aldex_validate_contrasts(unnamed.columns, X),
               "column names")
  C <- named(); colnames(C)[3] <- "missing"
  expect_error(ALDEx3:::.aldex_validate_contrasts(C, X), "complete coefficient")
  C <- named(); rownames(C) <- ""
  expect_error(ALDEx3:::.aldex_validate_contrasts(C, X), "row names")
  C <- rbind(age=c(0, 1, 0), age=c(0, 0, 1)); colnames(C) <- colnames(X)
  expect_error(ALDEx3:::.aldex_validate_contrasts(C, X), "row names")
  C <- named(); colnames(C)[2] <- colnames(C)[1]
  expect_error(ALDEx3:::.aldex_validate_contrasts(C, X), "column names")
  C <- named(); C[1] <- Inf
  expect_error(ALDEx3:::.aldex_validate_contrasts(C, X), "finite")
  dependent <- rbind(age=c(0, 1, 0), twice_age=c(0, 2, 0))
  colnames(dependent) <- colnames(X)
  expect_equal(ALDEx3:::.aldex_validate_contrasts(dependent, X), dependent)
})

test_that("rank deficient and nonpositive-residual designs fail before sampling", {
  data <- data.frame(x=c(-1, -1, 1, 1), g=factor(c(1, 1, 2, 2)))
  C <- matrix(c(0, 1, 2), 1,
              dimnames=list("x", c("(Intercept)", "x", "I(2 * x)")))
  testthat::local_mocked_bindings(aldex.sampler=function(...) stop("sampled"),
                                  .package="ALDEx3")
  for (method in c("lme4", "nlme", "blmm")) {
    extra <- if (method == "nlme") list(random=~1 | g) else list()
    expect_error(do.call(aldex, c(list(Y=matrix(10, 1, 4), X=~x + I(2 * x),
                                   data=data, method=method, scale=clr.sm,
                                   contrasts=C), extra)),
                 "full-rank")
  }
  X <- rbind(intercept=1, x=c(-1, 1))
  C2 <- matrix(c(0, 1), 1, dimnames=list("x", rownames(X)))
  expect_error(aldex(matrix(10, 1, 2), X, nsample=1, scale=matrix(0, 2, 1),
                     contrasts=C2), "positive residual")
})

test_that("LM identity contrasts preserve every coefficient field for all selectors", {
  set.seed(403)
  X <- cbind(intercept=1, age=rnorm(12), genotype=rep(c(0, 1), each=6))
  Y <- array(rnorm(12 * 2 * 3), c(12, 2, 3))
  C <- diag(ncol(X)); rownames(C) <- colnames(X); colnames(C) <- colnames(X)
  for (selector in c("t", "t.HC0", "t.HC3")) {
    fit <- fflm(Y, X, selector, C)
    for (field in c("estimate", "std.error", "df", "p.lower", "p.upper")) {
      expect_equal(unname(fit[[paste0("contrast.", field)]]), fit[[field]],
                   info=selector)
    }
  }
})

test_that("three age by genotype hypotheses survive an independent reparameterization", {
  set.seed(404)
  n <- 18
  age <- seq(-1, 1, length.out=n)
  genotype <- rep(c(0, 1), each=n / 2)
  Y <- array(rnorm(n * 2 * 2), c(n, 2, 2))
  X1 <- cbind(intercept=1, age=age, genotype=genotype, age_genotype=age*genotype)
  C1 <- rbind(age_g0=c(0, 1, 0, 0), genotype_at_zero=c(0, 0, 1, 0),
              age_difference=c(0, 0, 0, 1))
  colnames(C1) <- colnames(X1)
  X2 <- cbind(intercept=1, age_g0=age*(1-genotype), age_g1=age*genotype,
              genotype=genotype)
  C2 <- rbind(age_g0=c(0, 1, 0, 0), genotype_at_zero=c(0, 0, 0, 1),
              age_difference=c(0, -1, 1, 0))
  colnames(C2) <- colnames(X2)
  a <- fflm(Y, X1, "t.HC3", C1)
  b <- fflm(Y, X2, "t.HC3", C2)
  for (field in c("estimate", "std.error", "df", "p.lower", "p.upper")) {
    expect_equal(a[[paste0("contrast.", field)]],
                 b[[paste0("contrast.", field)]], tolerance=1e-11)
  }
})

test_that("P14/P56 WT-KO contrasts match an independent cell-means design", {
  age <- factor(rep(c("P14", "P56"), each=8), levels=c("P14", "P56"))
  genotype <- factor(rep(rep(c("WT", "KO"), each=4), 2), levels=c("WT", "KO"))
  cell <- interaction(age, genotype, sep="_")
  set.seed(4041)
  Y <- array(rnorm(16 * 2 * 2), c(16, 2, 2))
  X <- model.matrix(~age * genotype)
  C <- rbind(P14_KO_minus_WT=c(0, 0, 1, 0),
             P56_KO_minus_WT=c(0, 0, 1, 1),
             change_in_KO_effect=c(0, 0, 0, 1))
  colnames(C) <- colnames(X)
  X.cell <- cbind(P14_WT=as.numeric(cell == "P14_WT"),
                  P14_KO=as.numeric(cell == "P14_KO"),
                  P56_WT=as.numeric(cell == "P56_WT"),
                  P56_KO=as.numeric(cell == "P56_KO"))
  C.cell <- rbind(P14_KO_minus_WT=c(-1, 1, 0, 0),
                  P56_KO_minus_WT=c(0, 0, -1, 1),
                  change_in_KO_effect=c(1, -1, -1, 1))
  colnames(C.cell) <- colnames(X.cell)
  factorial <- fflm(Y, X, "t.HC3", C)
  cells <- fflm(Y, X.cell, "t.HC3", C.cell)
  expect_equal(factorial$contrast.estimate, cells$contrast.estimate)
  expect_equal(factorial$contrast.std.error, cells$contrast.std.error)
  expect_equal(factorial$contrast.p.lower, cells$contrast.p.lower)
  means <- c(P14_WT=mean(Y[age == "P14" & genotype == "WT", 1, 1]),
             P14_KO=mean(Y[age == "P14" & genotype == "KO", 1, 1]),
             P56_WT=mean(Y[age == "P56" & genotype == "WT", 1, 1]),
             P56_KO=mean(Y[age == "P56" & genotype == "KO", 1, 1]))
  expect_equal(unname(factorial$contrast.estimate[, 1, 1]),
               unname(c(means[2] - means[1], means[4] - means[3],
                        (means[4] - means[3]) - (means[2] - means[1]))))
})

test_that("contrast sign and rescaling transform estimates, errors, and tails", {
  set.seed(405)
  X <- cbind(intercept=1, x=rnorm(14))
  Y <- array(rnorm(14 * 2 * 2), c(14, 2, 2))
  C <- rbind(base=c(0, 1), doubled=c(0, 2), reversed=c(0, -1))
  colnames(C) <- colnames(X)
  fit <- fflm(Y, X, "t.HC0", C)
  expect_equal(fit$contrast.estimate["doubled", , ], 2 * fit$contrast.estimate["base", , ])
  expect_equal(fit$contrast.std.error["doubled", , ], 2 * fit$contrast.std.error["base", , ])
  expect_equal(fit$contrast.df["doubled", , ], fit$contrast.df["base", , ])
  expect_equal(fit$contrast.estimate["reversed", , ], -fit$contrast.estimate["base", , ])
  expect_equal(fit$contrast.std.error["reversed", , ], fit$contrast.std.error["base", , ])
  expect_equal(fit$contrast.p.lower["reversed", , ], fit$contrast.p.upper["base", , ])
})

test_that("singleton contrast output supports one-sided inference and summary", {
  Y <- matrix(c(10, 12, 15, 17, 21, 23), 1,
              dimnames=list("gene", paste0("s", 1:6)))
  X <- rbind(intercept=1, x=c(-1, -1, -1, 1, 1, 1))
  C <- matrix(c(0, 1), 1, dimnames=list("x_contrast", rownames(X)))
  saved <- array(c(-.1, .1, 0, .8, 1.1, .9), c(1, 6, 1))
  testthat::local_mocked_bindings(aldex.sampler=function(...) list(logW=saved),
                                  .package="ALDEx3")
  fit <- aldex(Y, X, nsample=1, scale=matrix(0, 6, 1), contrasts=C,
               onesided="lower",
               return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  expect_identical(dim(fit$contrasts$estimate), c(1L, 1L, 1L))
  expect_identical(dim(fit$contrasts$p.val), c(1L, 1L))
  expect_equal(fit$contrasts$p.val, fit$contrasts$p.val.adj)
  direct <- fflm(aperm(saved, c(2, 1, 3)), t(X), "t.HC3", C)
  expect_equal(unname(fit$contrasts$p.val),
               matrix(unname(direct$contrast.p.lower[, , 1]), 1, 1))
  for (field in c("estimate", "std.error", "df", "p.val", "p.val.adj")) {
    expect_true(all(is.finite(fit[[field]])), info=paste("coefficient", field))
    expect_true(all(is.finite(fit$contrasts[[field]])), info=field)
  }
  expect_equal(nrow(summary(fit, component="contrasts")), 1L)
  expect_error(aldex.effect(fit, "x_contrast"), "coefficient effects only")
})

test_that("adding contrasts preserves saved samples and coefficient families", {
  set.seed(406)
  Y <- matrix(sample(10:80, 24, replace=TRUE), 3, 8)
  X <- rbind(intercept=1, x=rep(c(-.5, .5), each=4))
  C <- matrix(c(0, 1), 1, dimnames=list("x", rownames(X)))
  common <- list(Y=Y, X=X, nsample=3, scale=clr.sm,
                 return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj",
                               "logComp", "logScale"))
  set.seed(407); ordinary <- do.call(aldex, common)
  set.seed(407); contrasted <- do.call(aldex, c(common, list(contrasts=C)))
  for (field in c("estimate", "std.error", "df", "p.val", "p.val.adj",
                  "logComp", "logScale")) {
    expect_identical(contrasted[[field]], ordinary[[field]], info=field)
  }
})

test_that("BLMM contrasts leave coefficient results unchanged on saved draws", {
  data <- data.frame(g=factor(rep(1:8, each=5)),
                     x=rep(rep(c(-.5, .5), each=5), 4))
  N <- nrow(data)
  set.seed(445)
  subject <- rnorm(8, sd=.35)[data$g]
  saved <- array(c(.3 + .6 * data$x + subject + rnorm(N, sd=.08),
                   .5 + .5 * data$x + subject + rnorm(N, sd=.07)), c(2, N, 1))
  C <- matrix(c(0, 1), 1, dimnames=list("x", c("(Intercept)", "x")))
  common <- list(Y=matrix(10, 2, N), X=~x + (1 | g), data=data,
                 method="blmm", n.cores=1, nsample=1, scale=matrix(0, N, 1),
                 return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  testthat::local_mocked_bindings(aldex.sampler=function(...) list(logW=saved),
                                  .package="ALDEx3")
  baseline <- do.call(aldex, common)
  contrasted <- do.call(aldex, c(common, list(contrasts=C)))
  for (field in c("estimate", "std.error", "df", "p.val", "p.val.adj")) {
    expect_identical(contrasted[[field]], baseline[[field]], info=field)
  }
  expect_true(all(is.finite(contrasted$contrasts$estimate)))
})

test_that("contrast p-value aggregation uses a separate explicit family oracle", {
  lo <- array(c(.02, .4, .1, .8, .3, .6, .5, .05, .9, .2, .7, .15), c(2, 3, 2))
  up <- 1 - lo / 2
  oracle <- function(lower, upper, side=FALSE, method="BH") {
    dims <- dim(lower)
    if (identical(side, FALSE)) {
      lower <- array(pmin(1, 2 * lower), dims)
      upper <- array(pmin(1, 2 * upper), dims)
    }
    K <- dim(lower)[1]; D <- dim(lower)[2]; S <- dim(lower)[3]
    adjusted.lower <- adjusted.upper <- array(NA_real_, c(K, D, S))
    for (k in seq_len(K)) for (s in seq_len(S)) {
      adjusted.lower[k, , s] <- p.adjust(lower[k, , s], method)
      adjusted.upper[k, , s] <- p.adjust(upper[k, , s], method)
    }
    raw <- adjusted <- matrix(NA_real_, K, D)
    for (k in seq_len(K)) for (d in seq_len(D)) {
      if (identical(side, FALSE)) {
        raw[k, d] <- min(mean(lower[k, d, ]), mean(upper[k, d, ]))
        adjusted[k, d] <- min(mean(adjusted.lower[k, d, ]), mean(adjusted.upper[k, d, ]))
      } else if (side == "lower") {
        raw[k, d] <- mean(lower[k, d, ]); adjusted[k, d] <- mean(adjusted.lower[k, d, ])
      } else {
        raw[k, d] <- mean(upper[k, d, ]); adjusted[k, d] <- mean(adjusted.upper[k, d, ])
      }
    }
    list(p.res=raw, p.adj.res=adjusted)
  }
  for (side in list(FALSE, "lower", "upper")) {
    expect_equal(aldex.pvals(lo, up, "BH", side), oracle(lo, up, side))
  }
  coefficient.lo <- array(rev(lo), dim(lo))
  coefficient.up <- array(rev(up), dim(up))
  fake.fit <- list(
    estimate=array(0, c(2, 3, 2)), std.error=array(1, c(2, 3, 2)),
    df=array(8, c(2, 3, 2)), p.lower=coefficient.lo, p.upper=coefficient.up,
    contrast.estimate=array(0, c(2, 3, 2)),
    contrast.std.error=array(1, c(2, 3, 2)),
    contrast.df=array(8, c(2, 3, 2)),
    contrast.p.lower=lo, contrast.p.upper=up
  )
  X <- rbind(intercept=1, x=rep(c(-.5, .5), each=3))
  C <- rbind(first=c(0, 1), second=c(1, -1))
  colnames(C) <- rownames(X)
  testthat::local_mocked_bindings(
    aldex.sampler=function(...) list(logW=array(0, c(3, 6, 2))),
    fflm=function(...) fake.fit, .package="ALDEx3")
  for (adjustment in c("BH", "holm")) for (side in list(FALSE, "lower", "upper")) {
    public <- aldex(matrix(10, 3, 6), X, nsample=2, scale=matrix(0, 6, 2),
                    contrasts=C, p.adjust.method=adjustment, onesided=side,
                    return.pars=c("p.val", "p.val.adj"))
    coefficient.oracle <- oracle(coefficient.lo, coefficient.up, side, adjustment)
    contrast.oracle <- oracle(lo, up, side, adjustment)
    expect_equal(unname(public$p.val), coefficient.oracle$p.res)
    expect_equal(unname(public$contrasts$p.val), contrast.oracle$p.res)
    expect_equal(unname(public$contrasts$p.val.adj), contrast.oracle$p.adj.res)
    if (identical(adjustment, "BH") && identical(side, FALSE)) {
      expect_false(isTRUE(all.equal(unname(public$p.val.adj),
                                    unname(public$contrasts$p.val.adj))))
    }
  }
})

test_that("streaming combines saved sampler blocks through real LM fits", {
  set.seed(408)
  N <- 6L; D <- 2L
  X <- rbind(intercept=1, x=rep(c(-.5, .5), each=3))
  C <- matrix(c(0, 1), 1, dimnames=list("x", rownames(X)))
  blocks <- lapply(seq_len(3L), function(i) array(rnorm(D * N), c(D, N, 1L)))
  calls <- 0L
  testthat::local_mocked_bindings(aldex.sampler=function(...) {
    calls <<- calls + 1L
    list(logW=blocks[[calls]])
  }, .package="ALDEx3")
  streamed <- aldex(matrix(10, D, N), X, nsample=3, scale=matrix(0, N, 3),
                    streamsize=.00096, contrasts=C,
                    return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  per.chunk <- lapply(blocks, function(logW) fflm(aperm(logW, c(2, 1, 3)),
                                                    t(X), "t.HC3", C))
  fields <- c("estimate", "std.error", "df", "p.lower", "p.upper",
              "contrast.estimate", "contrast.std.error", "contrast.df",
              "contrast.p.lower", "contrast.p.upper")
  combined <- lapply(fields, function(field) {
    abind::abind(lapply(per.chunk, function(z) z[[field]]), along=3)
  })
  names(combined) <- fields
  expect_true(streamed$streaming)
  expect_equal(unname(streamed$estimate), unname(combined$estimate))
  expect_equal(unname(streamed$std.error), unname(combined$std.error))
  expect_equal(unname(streamed$df), unname(combined$df))
  expected <- aldex.pvals(combined$p.lower, combined$p.upper, "BH")
  expected.contrast <- aldex.pvals(combined$contrast.p.lower,
                                   combined$contrast.p.upper, "BH")
  expect_equal(unname(streamed$p.val), expected$p.res)
  expect_equal(unname(streamed$contrasts$estimate), unname(combined$contrast.estimate))
  expect_equal(unname(streamed$contrasts$p.val.adj), expected.contrast$p.adj.res)
})

test_that("combine.streams preserves singleton contrast dimensions", {
  chunks <- list(
    list(contrast=array(1, c(1, 2, 1), dimnames=list("c", NULL, NULL))),
    list(contrast=array(2, c(1, 2, 1), dimnames=list("c", NULL, NULL)))
  )
  combined <- ALDEx3:::combine.streams(chunks)$contrast
  expect_identical(dim(combined), c(1L, 2L, 2L))
  expect_equal(combined[1, , ], rbind(c(1, 2), c(1, 2)))
})

test_that("streamed BLMM is checked against its explicit chunks, not one anchor", {
  set.seed(409)
  data <- data.frame(g=factor(rep(1:8, each=5)),
                     x=rep(rep(c(-.5, .5), each=5), 4))
  N <- nrow(data); D <- 2L
  C <- matrix(c(0, 1), 1, dimnames=list("x", c("(Intercept)", "x")))
  subject <- rnorm(8, sd=.35)[data$g]
  base <- .4 + .7 * data$x + subject
  blocks <- list(
    array(rep(base, D * 2L), c(D, N, 2L)) + array(rnorm(D * N * 2L, sd=.04), c(D, N, 2L)),
    array(rep(base, D), c(D, N, 1L)) + array(rnorm(D * N, sd=.04), c(D, N, 1L))
  )
  calls <- 0L
  testthat::local_mocked_bindings(aldex.sampler=function(...) {
    calls <<- calls + 1L
    list(logW=blocks[[calls]])
  }, .package="ALDEx3")
  streamed <- aldex(matrix(10, D, N), ~x + (1 | g), data, method="blmm",
                    n.cores=1, nsample=3, scale=matrix(0, N, 3),
                    streamsize=.013, contrasts=C,
                    return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  chunks <- lapply(blocks, function(logW) {
    blmm(aperm(logW, c(2, 1, 3)), ~x + (1 | g), data, n.cores=1,
         inference="anchor", contrasts=C)
  })
  fields <- c("contrast.estimate", "contrast.std.error", "contrast.df",
              "contrast.p.lower", "contrast.p.upper")
  expected <- stats::setNames(lapply(fields, function(field) {
    abind::abind(lapply(chunks, function(z) z[[field]]), along=3)
  }), fields)
  expect_true(streamed$streaming)
  expect_equal(unname(streamed$contrasts$estimate), unname(expected$contrast.estimate))
  expect_equal(unname(streamed$contrasts$std.error), unname(expected$contrast.std.error))
  expect_equal(unname(streamed$contrasts$df), unname(expected$contrast.df))
  expected.p <- aldex.pvals(expected$contrast.p.lower, expected$contrast.p.upper, "BH")
  expect_equal(unname(streamed$contrasts$p.val), expected.p$p.res)
  expect_equal(unname(streamed$contrasts$p.val.adj), expected.p$p.adj.res)
  expect_true(all(is.finite(streamed$contrasts$estimate)))
  # BLMM anchors are fit per streamed block, so whole-run bitwise equivalence
  # is neither expected nor asserted here.
})
