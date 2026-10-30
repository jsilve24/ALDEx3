# Numerical references use observation-space GLS, independently of the native
# random-effect cross-products and factorization used by BLMM.
blmm_sensitivity_fixture <- function(design = c("slope", "nested", "crossed")) {
  design <- match.arg(design)
  set.seed(3091)
  sizes <- c(3L, 6L, 4L, 7L, 5L, 8L, 4L, 6L)
  g <- rep(seq_along(sizes), sizes)
  n <- length(g)
  data <- data.frame(g = factor(g), block = factor(ceiling(g / 2)),
    h = factor(rep(seq_len(5L), length.out = n)), x = runif(n, -1, 1))
  formula <- switch(design,
    slope = ~ x + (1 + x | g),
    nested = ~ x + (1 | block) + (1 | g),
    crossed = ~ x + (1 | g) + (1 | h))
  mean_y <- 1 + 0.4 * data$x + rnorm(length(sizes), sd = 0.6)[g]
  Y <- vapply(seq_len(3L), function(s)
    mean_y + rnorm(n, sd = 0.3 + 0.15 * s), numeric(n))
  list(data = data, formula = formula, Y = Y)
}

blmm_observation_sensitivities <- function(parsed, Y, phi, is_log, contrasts) {
  quantities <- function(location) {
    theta <- location
    theta[is_log == 1L] <- exp(theta[is_log == 1L])
    Lambdat <- parsed$reTrms$Lambdat
    Lambdat@x <- theta[parsed$reTrms$Lind]
    ZLambda <- as.matrix(Matrix::t(parsed$reTrms$Zt) %*% Matrix::t(Lambdat))
    V <- diag(nrow(Y)) + tcrossprod(ZLambda)
    X <- parsed$X
    ViX <- solve(V, X)
    M <- crossprod(X, ViX)
    beta <- solve(M, crossprod(X, solve(V, Y)))
    residual <- Y - X %*% beta
    list(residual = log(colSums(residual * solve(V, residual))),
      variance = log(colSums(contrasts * solve(M, contrasts))))
  }
  u <- matrix(0, length(phi), ncol(Y))
  b <- matrix(0, length(phi), ncol(contrasts))
  for (j in seq_along(phi)) {
    step <- 1e-5 * max(1, abs(phi[j]))
    plus <- minus <- phi
    plus[j] <- plus[j] + step
    minus[j] <- minus[j] - step
    qp <- quantities(plus)
    qm <- quantities(minus)
    u[j, ] <- (qp$residual - qm$residual) / (2 * step)
    b[j, ] <- (qp$variance - qm$variance) / (2 * step)
  }
  list(u = u, b = b,
    scores = (nrow(Y) - ncol(parsed$X)) / 2 * (rowMeans(u) - u))
}

blmm_sensitivity_at <- function(obj, phi, parsed, Y, contrasts = NULL) {
  report <- obj$report(phi)
  ALDEx3:::blmm_scores_batch(obj, phi, report$PWRSS_s,
    nrow(Y), ncol(parsed$X), contrasts = contrasts, inference = TRUE)
}

test_that("analytic sensitivities agree with observation-space GLS derivatives", {
  for (design in c("slope", "nested", "crossed")) {
    z <- blmm_sensitivity_fixture(design)
    for (backend in c("dense", "sparse")) {
      parsed <- ALDEx3:::blmm_parse_formula(z$formula, z$data,
        backend = backend, draws = ncol(z$Y))
      q <- length(parsed$lower)
      bounded <- is.finite(parsed$lower)
      theta <- rep(0.7, q)
      theta[!bounded] <- -0.25
      p <- ncol(parsed$X)
      contrasts <- cbind(diag(p), c(1, -0.7), 3 * c(1, -0.7))
      coordinate_results <- list()
      for (logged in c(FALSE, TRUE)) {
        is_log <- as.integer(logged & bounded)
        phi <- theta
        phi[is_log == 1L] <- log(theta[is_log == 1L])
        obj <- ALDEx3:::blmm_make_adfun(parsed$X, z$Y, parsed$basis,
          is_log, phi)
        actual <- blmm_sensitivity_at(obj, phi, parsed, z$Y, contrasts)
        reference <- blmm_observation_sensitivities(parsed, z$Y, phi,
          is_log, contrasts)
        label <- paste(design, backend, if (logged) "log" else "natural")
        expect_identical(dim(actual$u), c(q, ncol(z$Y)), info = label)
        expect_identical(dim(actual$b), c(q, ncol(contrasts)), info = label)
        for (field in c("u", "b", "scores"))
          expect_equal(actual[[field]], reference[[field]], tolerance = 2e-6,
            info = paste(label, field))
        expect_equal(rowSums(actual$scores), rep(0, q), tolerance = 1e-12,
          info = label)
        expect_equal(actual$b[, p + 1L], actual$b[, p + 2L], tolerance = 1e-12,
          info = paste(label, "contrast rescaling"))
        default <- blmm_sensitivity_at(obj, phi, parsed, z$Y)
        expect_equal(default$b, actual$b[, seq_len(p), drop = FALSE],
          tolerance = 1e-12, info = paste(label, "identity contrasts"))
        scores <- ALDEx3:::blmm_scores_batch(obj, phi,
          obj$report(phi)$PWRSS_s, nrow(z$Y), p)
        expect_equal(scores, actual$scores, tolerance = 1e-12, info = label)
        coordinate_results[[if (logged) "log" else "natural"]] <- actual
      }
      scale <- ifelse(bounded, theta, 1)
      for (field in c("u", "b", "scores"))
        expect_equal(coordinate_results$log[[field]],
          scale * coordinate_results$natural[[field]], tolerance = 1e-10,
          info = paste(design, backend, "chain rule", field))
    }
  }
})

test_that("natural sensitivities retain structurally nonzero boundary derivatives", {
  z <- blmm_sensitivity_fixture("slope")
  contrasts <- cbind(c(1, 0), c(0, 1), c(1, -0.7))
  for (backend in c("dense", "sparse")) {
    parsed <- ALDEx3:::blmm_parse_formula(z$formula, z$data, backend = backend,
      draws = ncol(z$Y))
    expect_identical(length(parsed$lower), 3L)
    # Record at zero to detect lost sparse symbolic entries on later evaluations.
    obj <- ALDEx3:::blmm_make_adfun(parsed$X, z$Y, parsed$basis,
      integer(3L), numeric(3L))
    for (theta in list(c(0, 0, 0), c(0, -0.25, 0.7), c(0.7, -0.25, 0))) {
      actual <- blmm_sensitivity_at(obj, theta, parsed, z$Y, contrasts)
      reference <- blmm_observation_sensitivities(parsed, z$Y, theta,
        integer(3L), contrasts)
      label <- paste(backend, paste(theta, collapse = ","))
      for (field in c("u", "b", "scores"))
        expect_equal(actual[[field]], reference[[field]], tolerance = 2e-6,
          info = paste(label, field))
      if (all(theta == 0)) {
        expect_equal(actual$u, matrix(0, 3L, ncol(z$Y)), tolerance = 1e-12)
        expect_equal(actual$b, matrix(0, 3L, ncol(contrasts)), tolerance = 1e-12)
      } else if (theta[1L] == 0) {
        expect_gt(max(abs(c(actual$u[1L, ], actual$b[1L, ]))), 1e-3,
          label = paste(backend, "zero first diagonal with nonzero off-diagonal"))
      } else {
        expect_equal(actual$u[3L, ], rep(0, ncol(z$Y)), tolerance = 1e-12)
        expect_equal(actual$b[3L, ], rep(0, ncol(contrasts)), tolerance = 1e-12)
      }
    }
  }
})

test_that("analytic sensitivities preserve singleton dimensions and response scaling", {
  z <- blmm_sensitivity_fixture("crossed")
  for (formula in list(~ x + (1 | g), z$formula)) {
    parsed <- ALDEx3:::blmm_parse_formula(formula, z$data, backend = "sparse")
    theta <- rep(0.6, length(parsed$lower))
    C <- matrix(c(1, -0.7), ncol = 1L)
    results <- lapply(c(1, 5), function(scale) {
      Y <- scale * z$Y[, 1L, drop = FALSE]
      obj <- ALDEx3:::blmm_make_adfun(parsed$X, Y, parsed$basis,
        integer(length(theta)), theta)
      blmm_sensitivity_at(obj, theta, parsed, Y, C)
    })
    expect_identical(dim(results[[1L]]$u), c(length(theta), 1L))
    expect_identical(dim(results[[1L]]$b), c(length(theta), 1L))
    expect_identical(dim(results[[1L]]$scores), c(length(theta), 1L))
    expect_equal(results[[1L]], results[[2L]], tolerance = 1e-10)
    expect_equal(results[[1L]]$scores, matrix(0, length(theta), 1L),
      tolerance = 1e-12)
  }
  parsed <- ALDEx3:::blmm_parse_formula(~ 1 + (1 | g), z$data, backend = "sparse")
  obj <- ALDEx3:::blmm_make_adfun(parsed$X, z$Y, parsed$basis, 0L, 0.6)
  result <- blmm_sensitivity_at(obj, 0.6, parsed, z$Y,
    matrix(c(1, 3), nrow = 1L))
  expect_identical(dim(result$b), c(1L, 2L))
  expect_equal(result$b[, 1L], result$b[, 2L], tolerance = 1e-12)
})
