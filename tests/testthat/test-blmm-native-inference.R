test_that("native inference reports agree with independent conditional solves", {
  set.seed(1294)
  data <- data.frame(g = factor(rep(seq_len(12), each = 4)),
    site = factor(rep(seq_len(6), each = 8)), x = rep(seq(-1, 1, length.out = 4), 12))
  Y <- vapply(seq_len(3), function(i)
    0.5 + 0.8 * data$x + rnorm(nrow(data)), numeric(nrow(data)))
  for (formula in list(~ x + (1 + x | g), ~ x + (1 | g) + (1 | site))) {
    for (backend in c("dense", "sparse")) {
      parsed <- ALDEx3:::blmm_parse_formula(formula, data, backend = backend)
      X <- parsed$X
      p <- ncol(X)
      vectors <- cbind(diag(p), c(1, -0.7))
      offsets <- qr.coef(qr(X), Y)
      residual <- qr.resid(qr(X), Y)
      metadata <- ALDEx3:::blmm_random_effect_metadata(parsed$reTrms)
      q <- length(parsed$reTrms$theta)
      interior <- rep(0.7, q)
      interior[!is.finite(parsed$lower)] <- -0.25
      partial <- interior
      first_block <- metadata$blocks[[1]]
      column <- (first_block$lower - 1L) %/% first_block$size + 1L
      partial[first_block$theta[column == first_block$size]] <- 0
      obj <- ALDEx3:::blmm_make_adfun(X, residual, parsed$basis,
        integer(q), interior)
      for (theta in list(interior, partial, rep(0, q))) {
        report <- obj$report(theta)
        expect_identical(dim(report$beta), c(p, ncol(Y)))
        phi <- theta
        bounded <- is.finite(parsed$lower)
        phi[bounded] <- log(theta[bounded])
        pieces <- ALDEx3:::blmm_common_pieces(phi, X, parsed$basis, parsed$lower)
        for (s in seq_len(ncol(Y))) {
          reference <- ALDEx3:::blmm_fixed_effects_draw(pieces, X, residual[, s])
          sigma2 <- report$PWRSS_s[s] / (nrow(X) - p)
          label <- paste(backend, deparse(formula), s)
          expect_equal(as.vector(report$beta[, s]), reference$beta, tolerance = 1e-9, info = label)
          expect_equal(sigma2, reference$sigma2, tolerance = 1e-9, info = label)
          expect_equal(sqrt(sigma2 * colSums(forwardsolve(report$LX, vectors)^2)),
            sqrt(reference$sigma2 * colSums(forwardsolve(pieces$LX, vectors)^2)),
            tolerance = 1e-9, info = label)
          raw <- ALDEx3:::blmm_fixed_effects_draw(pieces, X, Y[, s])
          expect_equal(as.vector(report$beta[, s] + offsets[, s]), raw$beta,
            tolerance = 1e-9, info = label)
        }
      }
    }
  }
})

test_that("invalid reports rescue while a stale native kernel reports its error", {
  old <- options(ALDEx3.blmm.adaptive.diagnostics = TRUE)
  on.exit(options(old), add = TRUE)
  z <- inference_fixture(groups = 12, S = 3L, design = "within")
  original_object <- ALDEx3:::blmm_make_adfun
  original_exact <- ALDEx3:::blmm_exact_feature
  mode <- "invalid"
  report_calls <- exact_calls <- 0L
  testthat::local_mocked_bindings(
    blmm_make_adfun = function(X, Y, ...) {
      obj <- original_object(X, Y, ...)
      if (ncol(Y) == 1L) {
        report <- obj$report
        obj$report <- function(...) {
          report_calls <<- report_calls + 1L
          if (mode == "invalid") return(list())
          value <- report(...)
          value[[mode]] <- NULL
          value
        }
      }
      obj
    },
    blmm_exact_feature = function(...) {
      exact_calls <<- exact_calls + 1L
      original_exact(...)
    }, .package = "ALDEx3")
  fit <- feature_fit(z, inference = "adaptive")
  expect_gt(report_calls, 0L)
  expect_identical(fit$rejected_draws, seq_len(ncol(z$Y)))
  expect_identical(exact_calls, 1L)
  expect_match(fit$adaptive_diagnostics$diagnostics$reason[1], "invalid native inference report")
  for (field in c("beta", "dlogPWRSS", "dMinv")) {
    mode <- field
    report_calls <- exact_calls <- 0L
    expect_error(feature_fit(z, inference = "adaptive"), "reinstall.*compiled kernel")
    expect_gt(report_calls, 0L)
    expect_identical(exact_calls, 0L)
  }
})
