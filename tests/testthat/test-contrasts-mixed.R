mixed_contrast_matrix <- function(fixed_names) {
  C <- rbind(
    identity_intercept = c(1, 0),
    identity_x = c(0, 1),
    correlated = c(1, -0.5),
    rescaled = c(2, -1),
    reversed = c(-1, 0.5)
  )
  colnames(C) <- fixed_names
  C
}

mixed_contrast_fixture <- function(seed = 812, D = 2L, S = 2L) {
  set.seed(seed)
  g <- factor(rep(seq_len(12L), each = 2L))
  x <- rep(c(0, 1), 12L)
  data <- data.frame(g = g, x = x)
  N <- nrow(data)
  subject_effect <- rnorm(nlevels(g), sd = 0.35)[g]
  logW <- array(
    rnorm(N * D * S, sd = 0.08) +
      rep(0.4 + 0.7 * x + subject_effect, D * S),
    dim = c(N, D, S)
  )
  list(data = data, logW = logW,
       lme4_formula = ~ x + (1 | g), nlme_formula = ~ x)
}

feature_fit_with_contrasts <- function(z, C, inference = "anchor") {
  formula <- if (is.null(z$lme4_formula)) z$formula else z$lme4_formula
  p <- ALDEx3:::blmm_parse_formula(formula, z$data)
  do.call(ALDEx3:::blmm_fit_feature, c(list(
    d = 1L, logW = z$logW, X = p$X, basis = p$basis,
    is_log = p$is_log, phi_init = p$phi_init, lower = p$lower,
    reTrms = p$reTrms, formula = formula, data = z$data,
    inference = inference, contrasts = C
  )))
}

test_that("BLMM contrasts use the full GLS covariance and preserve linear rules", {
  z <- mixed_contrast_fixture(S = 1L)
  C <- mixed_contrast_matrix(colnames(model.matrix(~ x, z$data)))
  anchor <- feature_fit_with_contrasts(z, C, inference = "anchor")
  residual <- feature_fit_with_contrasts(z, C, inference = "residual")
  expect_identical(anchor$fallback_draws, 0L)
  expect_identical(residual$fallback_draws, 0L)

  # Direct lmerTest reference is independent of BLMM's conditional GLS code.
  direct_data <- z$data
  direct_data$y <- z$logW[, 1L, 1L]
  direct_fit <- lmerTest::lmer(y ~ x + (1 | g), data = direct_data)
  direct <- lmerTest::contest1D(direct_fit, L = C["correlated", ],
                                 ddf = "Satterthwaite")
  direct_t <- direct[1L, "Estimate"] / direct[1L, "Std. Error"]
  for (field in c("estimate", "std.error", "df")) {
    reference <- switch(field,
      estimate = direct[1L, "Estimate"],
      std.error = direct[1L, "Std. Error"],
      df = direct[1L, "df"]
    )
    expect_equal(unname(anchor$fit[[paste0("contrast.", field)]]["correlated", 1L, 1L]),
                 reference, tolerance = 0.002, info = field)
  }
  for (tail in c("p.lower", "p.upper")) {
    reference <- pt(direct_t, df = direct[1L, "df"],
                    lower.tail = identical(tail, "p.lower"))
    expect_lt(abs(log(anchor$fit[[paste0("contrast.", tail)]]["correlated", 1L, 1L]) -
                  log(reference)), 0.002)
  }

  for (fit in list(anchor$fit, residual$fit)) {
    for (coefficient_idx in seq_len(2L)) {
      coefficient <- colnames(C)[coefficient_idx]
      contrast <- paste0("identity_", if (coefficient == "(Intercept)") "intercept" else coefficient)
      expect_equal(unname(fit$contrast.estimate[contrast, , ]), fit$estimate[coefficient_idx, , ])
      expect_equal(unname(fit$contrast.std.error[contrast, , ]), fit$std.error[coefficient_idx, , ])
      expect_equal(unname(fit$contrast.df[contrast, , ]), fit$df[coefficient_idx, , ])
      expect_equal(unname(fit$contrast.p.lower[contrast, , ]), fit$p.lower[coefficient_idx, , ])
      expect_equal(unname(fit$contrast.p.upper[contrast, , ]), fit$p.upper[coefficient_idx, , ])
    }
    expect_equal(unname(fit$contrast.estimate["rescaled", , ]),
                 unname(2 * fit$contrast.estimate["correlated", , ]))
    expect_equal(unname(fit$contrast.std.error["rescaled", , ]),
                 unname(2 * fit$contrast.std.error["correlated", , ]))
    expect_equal(unname(fit$contrast.df["rescaled", , ]), unname(fit$contrast.df["correlated", , ]))
    expect_equal(unname(fit$contrast.estimate["reversed", , ]),
                 unname(-fit$contrast.estimate["correlated", , ]))
    expect_equal(unname(fit$contrast.std.error["reversed", , ]),
                 unname(fit$contrast.std.error["correlated", , ]))
    expect_equal(unname(fit$contrast.p.lower["reversed", , ]),
                 unname(fit$contrast.p.upper["correlated", , ]))
  }
  expect_equal(unname(residual$fit$contrast.df["correlated", 1L, 1L]),
               nrow(z$data) - ncol(model.matrix(~ x, z$data)))
  residual_t <- residual$fit$contrast.estimate["correlated", 1L, 1L] /
    residual$fit$contrast.std.error["correlated", 1L, 1L]
  expect_equal(residual$fit$contrast.p.upper["correlated", 1L, 1L],
               pt(residual_t, residual$fit$contrast.df["correlated", 1L, 1L],
                  lower.tail = FALSE))

  fit <- ALDEx3:::blmm(z$logW, z$lme4_formula, z$data,
                        contrasts = C, n.cores = 1L)
  fit_parallel <- ALDEx3:::blmm(z$logW, z$lme4_formula, z$data,
                                 contrasts = C, n.cores = 2L)
  expect_equal(fit_parallel, fit)
})

test_that("BLMM contrast rescue replaces every contrast field consistently", {
  z <- mixed_contrast_fixture(S = 2L)
  C <- mixed_contrast_matrix(colnames(model.matrix(~ x, z$data)))
  exact <- ALDEx3:::blmm_exact_feature(z$logW[, 1, ], z$lme4_formula, z$data,
                                       contrasts = C)
  baseline <- feature_fit_with_contrasts(z, C)
  expect_identical(baseline$fallback_draws, 0L)

  whole <- z$logW
  attr(whole, "blmm.force_fail_features") <- 1L
  expect_warning(rescued <- ALDEx3:::blmm(whole, z$lme4_formula, z$data,
                                          contrasts = C), "fell back")
  for (nm in names(exact)) {
    expect_equal(unname(rescued[[nm]][, 1, , drop = FALSE]), unname(exact[[nm]]),
                 info = nm)
  }

  original <- ALDEx3:::blmm_anchor_df
  calls <- 0L
  testthat::local_mocked_bindings(blmm_anchor_df = function(...) {
    calls <<- calls + 1L
    if (calls == 4L) stop("forced contrast inference failure")
    original(...)
  }, .package = "ALDEx3")
  partial <- feature_fit_with_contrasts(z, C)
  expect_false(partial$full_fallback)
  expect_identical(partial$rejected_draws, 2L)
  for (nm in names(exact)) {
    expect_equal(unname(partial$fit[[nm]][, , 2, drop = FALSE]),
                 unname(exact[[nm]][, , 2, drop = FALSE]), info = nm)
    expect_equal(partial$fit[[nm]][, , 1, drop = FALSE],
                 baseline$fit[[nm]][, , 1, drop = FALSE], info = paste("retained", nm))
  }
})

test_that("natural BLMM partial rescue preserves contrast and coefficient fields", {
  z <- inference_fixture(seed = 105002, groups = 12, design = "slope",
                         noise = .02, S = 8, effect = 0)
  C <- mixed_contrast_matrix(colnames(model.matrix(~ x, z$data)))
  baseline <- feature_fit_with_contrasts(z, C)
  fit <- feature_fit_with_contrasts(z, C)
  expect_false(fit$full_fallback)
  expect_gt(fit$fallback_draws, 0L)
  expect_lt(fit$fallback_draws, ncol(z$Y))
  expect_equal(fit$rejected_draws, baseline$rejected_draws)
  exact <- ALDEx3:::blmm_exact_feature(
    z$Y[, fit$rejected_draws, drop = FALSE], z$formula, z$data, contrasts = C
  )
  retained <- setdiff(seq_len(ncol(z$Y)), fit$rejected_draws)
  for (nm in names(fit$fit)) {
    expect_equal(unname(fit$fit[[nm]][, , fit$rejected_draws, drop = FALSE]),
                 unname(exact[[nm]]), info = paste("rescued", nm))
    expect_equal(fit$fit[[nm]][, , retained, drop = FALSE],
                 baseline$fit[[nm]][, , retained, drop = FALSE],
                 info = paste("retained", nm))
  }
})

test_that("failed exact draw rescue is not retried as a whole-feature fallback", {
  z <- mixed_contrast_fixture(D = 1L, S = 2L)
  C <- mixed_contrast_matrix(colnames(model.matrix(~ x, z$data)))
  expect_identical(feature_fit_with_contrasts(z, C)$fallback_draws, 0L)
  original <- ALDEx3:::blmm_anchor_df
  calls <- 0L
  refit_draws <- integer()
  local_mocked_bindings(
    blmm_anchor_df = function(...) {
      calls <<- calls + 1L
      if (calls == 4L) stop("forced contrast inference failure")
      original(...)
    },
    blmm_exact_feature = function(Y_d, ...) {
      refit_draws <<- c(refit_draws, ncol(Y_d))
      stop("forced exact fit failure")
    },
    .package = "ALDEx3"
  )
  expect_error(feature_fit_with_contrasts(z, C),
               "feature 1.*draw.*2.*forced exact fit failure")
  expect_identical(refit_draws, 1L)
})

test_that("invalid lmerTest contrast df preserves fits and propagates missing inference", {
  z <- mixed_contrast_fixture(D = 1L, S = 2L)
  C <- mixed_contrast_matrix(colnames(model.matrix(~ x, z$data)))
  fit <- function() ALDEx3:::sr.mem(z$logW, z$lme4_formula, z$data,
    n.cores = 1L, method = "lme4", mem.args = list(), contrasts = C)
  baseline <- fit()
  contest1D <- lmerTest::contest1D
  invalidated <- FALSE
  local_mocked_bindings(contest1D = function(model, L, ...) {
    result <- contest1D(model, L, ...)
    if (!invalidated && identical(as.numeric(L), as.numeric(C["correlated", ]))) {
      result[1L, "df"] <- NaN
      invalidated <<- TRUE
    }
    result
  }, .package = "lmerTest")
  expect_warning(actual <- fit(), "invalid contrast inference.*1.*NA")
  for (nm in c("estimate", "std.error", "df", "p.lower", "p.upper", "random.eff",
               "contrast.estimate", "contrast.std.error")) {
    expect_equal(actual[[nm]], baseline[[nm]], info = nm)
  }
  for (nm in c("contrast.df", "contrast.p.lower", "contrast.p.upper")) {
    expected <- baseline[[nm]]
    expected["correlated", 1L, 1L] <- NA_real_
    expect_equal(actual[[nm]], expected, info = nm)
  }
  aggregated <- ALDEx3:::aldex.pvals(actual$contrast.p.lower,
                                    actual$contrast.p.upper, "BH")
  expect_true(is.na(aggregated$p.res[3L, 1L]))
  expect_true(is.na(aggregated$p.adj.res[3L, 1L]))
  expect_true(all(is.finite(aggregated$p.res[-3L, , drop = FALSE])))
})

test_that("lmerTest and nlme calculate named contrast arrays with their native df", {
  z <- mixed_contrast_fixture()
  C <- mixed_contrast_matrix(colnames(model.matrix(~ x, z$data)))

  lmer_serial <- ALDEx3:::sr.mem(z$logW, z$lme4_formula, z$data,
    n.cores = 1L, method = "lme4", mem.args = list(), contrasts = C)
  lmer_parallel <- ALDEx3:::sr.mem(z$logW, z$lme4_formula, z$data,
    n.cores = 2L, method = "lme4", mem.args = list(), contrasts = C)
  expect_equal(lmer_parallel, lmer_serial)
  expect_equal(lmer_serial$contrast.estimate["identity_x", , ],
               lmer_serial$estimate["x", , ])
  expect_equal(lmer_serial$contrast.std.error["identity_x", , ],
               lmer_serial$std.error["x", , ])
  expect_equal(lmer_serial$contrast.df["identity_x", , ],
               lmer_serial$df["x", , ])

  nlme_fit <- ALDEx3:::sr.mem(z$logW, z$nlme_formula, z$data,
    n.cores = 1L, method = "nlme", mem.args = list(random = ~ 1 | g), contrasts = C)
  expect_equal(nlme_fit$contrast.estimate["identity_x", , ],
               nlme_fit$estimate["x", , ])
  expect_equal(nlme_fit$contrast.std.error["identity_x", , ],
               nlme_fit$std.error["x", , ])
  expect_equal(nlme_fit$contrast.df["identity_x", , ], nlme_fit$df["x", , ])
  expect_equal(nlme_fit$contrast.p.lower["reversed", , ],
               nlme_fit$contrast.p.upper["correlated", , ])

  direct_data <- z$data
  direct_data$y <- z$logW[, 1L, 1L]
  direct_fit <- nlme::lme(y ~ x, data = direct_data, random = ~ 1 | g)
  L <- matrix(C["correlated", ], nrow = 1L,
              dimnames = list(NULL, names(nlme::fixef(direct_fit))))
  direct_anova <- stats::anova(direct_fit, L = L)
  direct_estimate <- as.numeric(L %*% nlme::fixef(direct_fit))
  direct_se <- sqrt(as.numeric(L %*% stats::vcov(direct_fit) %*% t(L)))
  direct_t <- direct_estimate / direct_se
  expect_equal(unname(nlme_fit$contrast.estimate["correlated", 1L, 1L]), direct_estimate)
  expect_equal(unname(nlme_fit$contrast.std.error["correlated", 1L, 1L]), direct_se)
  expect_equal(unname(nlme_fit$contrast.df["correlated", 1L, 1L]), direct_anova[1L, "denDF"])
  expect_equal(unname(nlme_fit$contrast.p.upper["correlated", 1L, 1L]),
               pt(direct_t, direct_anova[1L, "denDF"], lower.tail = FALSE))
})

test_that("nlme rejects contrasts with mixed native denominator df", {
  set.seed(34)
  g <- factor(rep(seq_len(12L), each = 3L))
  x <- rep(-1:1, 12L)
  between <- rep(rnorm(12L), each = 3L)
  y <- 1 + 0.5 * x + 0.3 * between + rnorm(length(x)) + rnorm(12L, sd = 0.3)[g]
  data <- data.frame(g = g, x = x, between = between)
  C <- matrix(c(0, 1, 1), nrow = 1L,
              dimnames = list("mixed_df", c("(Intercept)", "x", "between")))
  expect_error(
    ALDEx3:::sr.mem(array(y, c(length(y), 1L, 1L)), ~ x + between, data,
      n.cores = 1L, method = "nlme", mem.args = list(random = ~ 1 | g),
      contrasts = C),
    "does not support.*different denominator"
  )
})

test_that("sr.mem retains a singleton fixed-effect coefficient dimension", {
  set.seed(512)
  g <- factor(rep(seq_len(6L), each = 2L))
  data <- data.frame(g = g)
  y <- rnorm(length(g)) + rnorm(6L, sd = 0.25)[g]
  C <- matrix(1, nrow = 1L,
              dimnames = list("intercept", "(Intercept)"))
  fit <- ALDEx3:::sr.mem(array(y, c(length(y), 1L, 1L)), ~ 1, data,
    n.cores = 1L, method = "nlme", mem.args = list(random = ~ 1 | g),
    contrasts = C)
  expect_identical(dim(fit$estimate), c(1L, 1L, 1L))
  expect_identical(dim(fit$contrast.estimate), c(1L, 1L, 1L))
  expect_equal(unname(fit$contrast.estimate), unname(fit$estimate))
  expect_equal(unname(fit$contrast.std.error), unname(fit$std.error))
  expect_equal(unname(fit$contrast.df), unname(fit$df))
})
