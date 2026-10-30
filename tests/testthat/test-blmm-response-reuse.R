# Changing responses must refresh both recorded derivative tapes as well as the
# double-valued inference report, including at an unchanged covariance vector.
blmm_response_reuse_fixture <- function() {
  set.seed(10307)
  sizes <- c(3L, 6L, 4L, 7L, 5L, 8L)
  g <- rep(seq_along(sizes), sizes)
  data <- data.frame(g = factor(g),
    h = factor(rep(seq_len(4L), length.out = length(g))),
    x = runif(length(g), -1, 1))
  Y <- cbind(0.5 + data$x + rnorm(length(g)),
    7 * (rnorm(length(g)) + rnorm(length(sizes))[g]))
  list(data = data, Y = Y)
}

test_that("response reuse preserves objective, derivatives and inference reports", {
  z <- blmm_response_reuse_fixture()
  for (formula in list(~ x + (1 + x | g), ~ x + (1 | g) + (1 | h))) {
    for (backend in c("dense", "sparse")) {
      parsed <- ALDEx3:::blmm_parse_formula(formula, z$data, backend = backend)
      q <- length(parsed$lower)
      bounded <- is.finite(parsed$lower)
      interior <- rep(0.7, q)
      interior[!bounded] <- -0.2
      partial <- interior
      partial[q] <- 0
      for (logged in c(FALSE, TRUE)) {
        is_log <- as.integer(logged & bounded)
        start <- if (logged) interior else numeric(q)
        if (logged) start[bounded] <- log(start[bounded])
        obj <- ALDEx3:::blmm_make_adfun(parsed$X, z$Y[, 1L, drop = FALSE],
          parsed$basis, is_log, start, reusable = TRUE)
        # Initialize both the objective and gradient tape before changing data.
        obj$fn(start)
        obj$gr(start)
        obj$he(start)
        locations <- if (logged) list(start) else
          list(numeric(q), partial, interior)
        for (s in c(1L, 2L, 1L)) {
          Y <- z$Y[, s, drop = FALSE]
          response <- ALDEx3:::blmm_response_cache(parsed$basis, parsed$X, Y)
          ALDEx3:::blmm_update_response(obj, response)
          reference <- ALDEx3:::blmm_make_adfun(parsed$X, Y, parsed$basis,
            is_log, start)
          for (phi in locations) {
            label <- paste(deparse(formula), backend, logged, s,
              paste(phi, collapse = ","))
            # Evaluate the existing Hessian tape first after a response change.
            expect_equal(obj$he(phi), reference$he(phi), tolerance = 1e-10,
              info = paste(label, "Hessian"))
            expect_equal(obj$gr(phi), reference$gr(phi), tolerance = 1e-10,
              info = paste(label, "gradient"))
            expect_equal(obj$fn(phi), reference$fn(phi), tolerance = 1e-10,
              info = paste(label, "objective"))
            expect_equal(obj$report(phi), reference$report(phi), tolerance = 1e-10,
              info = paste(label, "native report"))
          }
        }
      }
    }
  }
})

test_that("invalid response updates leave all cached data unchanged", {
  z <- blmm_response_reuse_fixture()
  parsed <- ALDEx3:::blmm_parse_formula(~ x + (1 + x | g), z$data)
  Y <- z$Y[, 1L, drop = FALSE]
  phi <- c(0.7, -0.2, 0.7)
  obj <- ALDEx3:::blmm_make_adfun(parsed$X, Y, parsed$basis,
    integer(3L), phi, reusable = TRUE)
  response <- ALDEx3:::blmm_response_cache(parsed$basis, parsed$X, Y)
  original_data <- obj$env$data
  original_value <- obj$fn(phi)
  invalid <- list()
  for (field in c("ZY", "XY", "yy")) {
    for (value in c(NA_real_, Inf, -Inf)) {
      candidate <- response
      candidate[[field]][1L] <- value
      invalid[[length(invalid) + 1L]] <- candidate
    }
    candidate <- response
    storage.mode(candidate[[field]]) <- "integer"
    invalid[[length(invalid) + 1L]] <- candidate
  }
  candidate <- response
  candidate$XY <- cbind(candidate$XY, candidate$XY)
  invalid[[length(invalid) + 1L]] <- candidate
  candidate <- response
  candidate$ZY <- as.vector(candidate$ZY)
  invalid[[length(invalid) + 1L]] <- candidate
  candidate <- response
  candidate$yy <- numeric()
  invalid[[length(invalid) + 1L]] <- candidate
  candidate <- response
  candidate$yy <- NULL
  invalid[[length(invalid) + 1L]] <- candidate
  for (candidate in invalid) {
    expect_error(ALDEx3:::blmm_update_response(obj, candidate))
    expect_identical(obj$env$data, original_data)
    expect_equal(obj$fn(phi), original_value, tolerance = 0)
  }
})

test_that("adaptive draws share one recording without changing accepted inference", {
  old <- options(ALDEx3.blmm.adaptive.diagnostics = TRUE)
  on.exit(options(old), add = TRUE)
  z <- inference_fixture(groups = 12L, S = 4L, design = "within")
  original_make <- TMB::MakeADFun
  original_object <- ALDEx3:::blmm_make_adfun
  constructions <- 0L
  testthat::local_mocked_bindings(
    MakeADFun = function(...) {
      constructions <<- constructions + 1L
      original_make(...)
    }, .package = "TMB")
  reused <- feature_fit(z, inference = "adaptive")
  expect_identical(reused$rejected_draws, integer())
  expect_identical(constructions, 2L) # One anchor and one draw-sized recording.

  # Reproduce the former adaptive implementation with one fresh draw objective
  # per response, retaining the same optimizer and inference checks.
  testthat::local_mocked_bindings(
    blmm_make_adfun = function(..., reusable = FALSE, object = NULL)
      original_object(..., reusable = FALSE, object = NULL),
    .package = "ALDEx3")
  constructions <- 0L
  fresh <- feature_fit(z, inference = "adaptive")
  expect_identical(constructions, 1L + ncol(z$Y))
  expect_equal(reused$fit, fresh$fit, tolerance = 1e-9)
  for (field in c("fallback", "fallback_message", "fallback_draws",
                  "rejected_draws", "full_fallback"))
    expect_identical(reused[[field]], fresh[[field]], info = field)
  expect_equal(reused$adaptive_diagnostics$diagnostics,
    fresh$adaptive_diagnostics$diagnostics, tolerance = 1e-9)
})

test_that("a failed draw discards its recording before the next response", {
  z <- inference_fixture(groups = 12L, S = 4L, design = "within")
  reference <- feature_fit(z, inference = "adaptive")
  original_object <- ALDEx3:::blmm_make_adfun
  draw_calls <- 0L
  testthat::local_mocked_bindings(
    blmm_make_adfun = function(X, Y_d, ...) {
      obj <- original_object(X, Y_d, ...)
      if (ncol(Y_d) == 1L) {
        draw_calls <<- draw_calls + 1L
        if (draw_calls == 2L) obj$report <- function(...) list()
      }
      obj
    }, .package = "ALDEx3")
  fit <- feature_fit(z, inference = "adaptive")
  expect_identical(fit$rejected_draws, 2L)
  expect_identical(fit$full_fallback, FALSE)
  for (field in names(reference$fit))
    expect_equal(fit$fit[[field]][, , c(1L, 3L, 4L), drop = FALSE],
      reference$fit[[field]][, , c(1L, 3L, 4L), drop = FALSE],
      tolerance = 1e-9, info = field)
})
