require(lmtest, warn.conflicts=FALSE)
require(sandwich, warn.conflicts=FALSE)

test_that("fflm correctness", {
  N <- 10
  D <- 10
  S <- 200
  P <- 3

  X <- matrix(runif(N*P), N, P)
  Y <- array(rnorm(N*D*S), c(N, D, S))

  res <- fflm(Y, X, test="t")

  tmp <- coefficients(summary(lm(Y[,3,7]~X-1)))
  ## test coefficients
  expect_equal(res$estimate[,3,7], unname(tmp[,1]))
  ## test std.errors
  expect_equal(res$std.error[,3,7], unname(tmp[,2]))
  ## test p-values
  p <- pmin(res$p.lower, res$p.upper)*2
  expect_equal(p[,3,7], unname(tmp[,4]))
})


test_that("fflm, robust correctness HC0", {
  N <- 10
  D <- 10
  S <- 200
  P <- 3

  X <- matrix(runif(N*P), N, P)
  Y <- array(rnorm(N*D*S), c(N, D, S))

  res <- fflm(Y, X, test="t.HC0")
  m <- lm(Y[,3,7]~X-1)
  tmp <- coefficients(summary(m))
  ## test coefficients
  expect_equal(res$estimate[,3,7], unname(tmp[,1]))
  ## test std.errors
  robust <- coeftest(m, vcov=vcovHC(m, type="HC0"))
  expect_equal(res$std.error[,3,7], unname(robust[,2]))
  ## test p-values
  p <- pmin(res$p.lower, res$p.upper)*2
  expect_equal(p[,3,7], unname(robust[,4]))
})



test_that("fflm, robust correctness HC3", {
  N <- 10
  D <- 10
  S <- 200
  P <- 3

  X <- matrix(runif(N*P), N, P)
  Y <- array(rnorm(N*D*S), c(N, D, S))

  res <- fflm(Y, X, test="t.HC3")
  m <- lm(Y[,3,7]~X-1)
  tmp <- coefficients(summary(m))
  ## test coefficients
  expect_equal(res$estimate[,3,7], unname(tmp[,1]))
  ## test std.errors
  robust <- coeftest(m, vcov=vcovHC(m, type="HC3"))
  expect_equal(res$std.error[,3,7], unname(robust[,2]))
  ## test p-values
  p <- pmin(res$p.lower, res$p.upper)*2
  expect_equal(p[,3,7], unname(robust[,4]))
})
