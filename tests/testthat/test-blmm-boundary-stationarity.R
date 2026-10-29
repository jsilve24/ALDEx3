# Quadratics give known constrained optima and exact raw derivatives.
boundary_quadratic <- function(center, curvature) {
  list(
    fn = function(x) {
      delta <- x - center
      as.numeric(crossprod(delta, curvature %*% delta) / 2)
    },
    gr = function(x) as.vector(curvature %*% (x - center)),
    he = function(x) curvature
  )
}

test_that("small scaled tangent correction permits boundary validation", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 1000)))
  start <- c(0, 1 + 1e-8)
  fit <- ALDEx3:::blmm_bounded_optimize(objective, start, c(0, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_gt(abs(objective$gr(start)[2]), 1e-6)
  expect_identical(fit$status, "boundary")
  expect_identical(fit$reason, "stationary boundary")
  expect_equal(fit$par, start)
})

test_that("large tangent correction does not qualify as stationary", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 1)))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1.01), c(0, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})

test_that("negative raw curvature prevents boundary acceptance", {
  curvature <- matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1000), 3, 3)
  objective <- boundary_quadratic(c(0, 0, 1), curvature)
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 0, 1 + 1e-8),
    c(0, 0, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "boundary negative curvature")
})

test_that("an inward descent direction can leave an active bound", {
  objective <- boundary_quadratic(c(1, 1), diag(2))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1), c(0, -Inf))
  expect_gt(fit$par[1], 0)
  expect_identical(fit$status, "interior")
})

test_that("ill-conditioned tangent information cannot certify stationarity", {
  objective <- boundary_quadratic(c(-1, 1, 1), diag(c(1, 1, 1e-8)))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-5, 1),
    c(0, -Inf, -Inf), max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})

test_that("a small high-curvature inward optimum is not accepted at the bound", {
  objective <- boundary_quadratic(1e-8, matrix(1e8, 1, 1))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, 0, 0, max_steps = 0L, validate_boundary = function(...) TRUE)
  expect_equal(objective$gr(0), -1)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})


test_that("the scaled route requires an explicit boundary validator", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 1000)))
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-8),
    c(0, -Inf), max_steps = 0L)
  expect_identical(fit$status, "unresolved")
  expect_identical(fit$reason, "iteration budget exhausted")
})

test_that("a rejected boundary validation continues within the iteration budget", {
  objective <- boundary_quadratic(c(-1, 1), diag(c(1, 10)))
  calls <- 0L
  validator <- function(location, ...) {
    calls <<- calls + 1L
    list(ok = abs(location[2] - 1) < 1e-12, inference = "accepted")
  }
  fit <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-6),
    c(0, -Inf), max_steps = 1L, validate_boundary = validator)
  expect_identical(calls, 2L)
  expect_identical(fit$status, "boundary")
  expect_identical(fit$inference, "accepted")
  expect_identical(fit$iterations, 1L)
  expect_equal(fit$par, c(0, 1), tolerance = 1e-10)
  exhausted <- ALDEx3:::blmm_bounded_optimize(objective, c(0, 1 + 1e-6),
    c(0, -Inf), max_steps = 0L, validate_boundary = validator)
  expect_identical(exhausted$status, "unresolved")
  expect_identical(exhausted$iterations, 0L)
})
