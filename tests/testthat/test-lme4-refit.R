refit_fixture <- function(slope = FALSE, unbalanced = FALSE, boundary_first = FALSE) {
  set.seed(8731 + slope)
  groups <- 16L
  sizes <- if (unbalanced) rep(c(4L, 6L, 8L), length.out = groups) else rep(6L, groups)
  data <- data.frame(g = factor(rep(seq_len(groups), sizes)),
    x = if (unbalanced) runif(sum(sizes), -1, 1) else
      rep(seq(-1, 1, length.out = 6L), groups),
    between = rep(rep(c(-0.5, 0.5), groups / 2L), sizes))
  u <- rnorm(groups)
  v <- 0.3 * u + rnorm(groups, sd = 0.5)
  residual <- rnorm(nrow(data), sd = 0.4)
  # This response has exactly zero group means after fixed-effect adjustment,
  # unlike the other draws. Its covariance and df cannot come from the template.
  centered <- residual - ave(residual, data$g)
  if (slope) centered <- ave(seq_len(nrow(data)), data$g, FUN = function(idx)
    qr.resid(qr(cbind(1, data$x[idx])), residual[idx]))
  base <- 2 + 0.4 * data$x - 0.3 * data$between
  random <- u[data$g] + if (slope) v[data$g] * data$x else 0
  Y <- cbind(base + random + residual, base + centered,
             base + 2 * random + 1.5 * residual)
  if (boundary_first) Y <- Y[, c(2L, 1L, 3L)]
  logW <- array(NA_real_, c(nrow(data), 2L, ncol(Y)))
  logW[, 1L, ] <- Y
  logW[, 2L, ] <- -0.7 * Y + 0.2 * data$x
  formula <- if (slope) ~ x + between + (1 + x | g) else
    ~ x + between + (1 | g)
  C <- rbind(identity_x = c(0, 1, 0), combination = c(1, -0.5, 0.25))
  colnames(C) <- colnames(model.matrix(~ x + between, data))
  list(logW = logW, data = data, formula = formula, C = C)
}

test_that("exact refits recompute covariance and Satterthwaite inference for each response", {
  for (z in list(refit_fixture(),
                 refit_fixture(slope = TRUE, unbalanced = TRUE),
                 refit_fixture(slope = TRUE, unbalanced = TRUE, boundary_first = TRUE))) {
    original_data <- z$data
    fit <- suppressMessages(ALDEx3:::sr.mem(z$logW, z$formula, z$data,
      n.cores = 1L, method = "lme4", mem.args = list(), contrasts = z$C))
    expect_identical(z$data, original_data)
    expect_identical(dim(fit$estimate), c(3L, 2L, 3L))
    expect_identical(dim(fit$contrast.df), c(2L, 2L, 3L))
    expect_gt(abs(fit$df["(Intercept)", 1L, 1L] -
                    fit$df["(Intercept)", 1L, 2L]), 20)

    for (s in seq_len(dim(z$logW)[3L])) for (d in seq_len(dim(z$logW)[2L])) {
      data <- z$data
      data$y <- z$logW[, d, s]
      reference <- suppressMessages(lmerTest::lmer(update(z$formula, y ~ .), data))
      coefs <- coef(summary(reference))
      expect_equal(fit$estimate[, d, s], coefs[, "Estimate"], tolerance = 1e-5)
      expect_equal(fit$std.error[, d, s], coefs[, "Std. Error"], tolerance = 2e-4)
      expect_equal(fit$df[, d, s], coefs[, "df"], tolerance = 2e-3)
      expect_equal(fit$p.upper[, d, s], pt(coefs[, "t value"], coefs[, "df"],
        lower.tail = FALSE), tolerance = 2e-5)
      expect_equal(unname(fit$random.eff[, d, s]),
        as.data.frame(lme4::VarCorr(reference))$vcov, tolerance = 2e-4)
      for (k in seq_len(nrow(z$C))) {
        contrast <- lmerTest::contest1D(reference, z$C[k, ], ddf = "Satterthwaite")
        expect_equal(unname(fit$contrast.estimate[k, d, s]),
          contrast[1L, "Estimate"], tolerance = 1e-5)
        expect_equal(unname(fit$contrast.std.error[k, d, s]),
          contrast[1L, "Std. Error"], tolerance = 2e-4)
        expect_equal(unname(fit$contrast.df[k, d, s]),
          contrast[1L, "df"], tolerance = 2e-3)
      }
    }
  }
})

test_that("exact refits retain transformed formula data and offsets", {
  z <- refit_fixture()
  data <- z$data
  data$exposure <- exp(seq(-0.2, 0.2, length.out = nrow(data)))
  formula <- ~ I(x^2) + between + offset(log(exposure)) + (1 | g)
  actual <- suppressMessages(ALDEx3:::sr.mem(z$logW, formula, data,
    n.cores = 1L, method = "lme4", mem.args = list()))
  for (s in 1:3) {
    data$y <- z$logW[, 2L, s]
    reference <- suppressMessages(lmerTest::lmer(update(formula, y ~ .), data))
    coefs <- coef(summary(reference))
    expect_equal(actual$estimate[, 2L, s], coefs[, "Estimate"], tolerance = 1e-5)
    expect_equal(actual$std.error[, 2L, s], coefs[, "Std. Error"], tolerance = 2e-4)
    expect_equal(actual$df[, 2L, s], coefs[, "df"], tolerance = 2e-3)
  }
})

test_that("exact refits preserve serial and socket-worker results", {
  z <- refit_fixture(slope = TRUE)
  fit <- function(cores) suppressWarnings(suppressMessages(ALDEx3:::sr.mem(
    z$logW, z$formula, z$data, n.cores = cores,
    method = "lme4", mem.args = list(), contrasts = z$C)))
  expect_equal(fit(2L), fit(1L), tolerance = 1e-10)
})

test_that("refit warnings and errors retry fresh fits without hiding fresh failures", {
  z <- refit_fixture()
  data <- z$data
  data$y <- z$logW[, 2L, 1L]
  expected <- coef(summary(lmerTest::lmer(update(z$formula, y ~ .), data)))
  for (signal in c("warning", "error")) local({
    testthat::local_mocked_bindings(refit = function(...) {
      if (signal == "warning") warning("forced refit warning") else
        stop("forced refit error")
    }, .package = "lme4")
    expect_no_warning(actual <- ALDEx3:::sr.mem(z$logW, z$formula, z$data,
      n.cores = 1L, method = "lme4", mem.args = list()))
    expect_equal(actual$estimate[, 2L, 1L], expected[, "Estimate"])
    expect_equal(actual$std.error[, 2L, 1L], expected[, "Std. Error"])
    expect_equal(actual$df[, 2L, 1L], expected[, "df"])
  })
  original <- lmerTest::lmer
  calls <- 0L
  testthat::local_mocked_bindings(lmer = function(...) {
    calls <<- calls + 1L
    if (calls > 1L) stop("forced fresh-fit failure")
    original(...)
  }, .package = "lmerTest")
  testthat::local_mocked_bindings(refit = function(...) stop("forced refit error"),
    .package = "lme4")
  expect_error(ALDEx3:::sr.mem(z$logW, z$formula, z$data,
    n.cores = 1L, method = "lme4", mem.args = list()), "forced fresh-fit failure")
})
