test_that("fflm correctness", {
  N <- 10
  D <- 10
  S <- 200
  P <- 3

  X <- matrix(runif(N*P), N, P)
  Y <- array(rnorm(N*D*S), c(N, D, S))

  res <- fflm(Y, X)

  tmp <- coefficients(summary(lm(Y[,3,7]~X-1)))
  ## test coefficients
  expect_equal(res$estimate[,3,7], unname(tmp[,1]))
  ## test std.errors
  expect_equal(res$std.error[,3,7], unname(tmp[,2]))
  ## test p-values
  p <- pmin(res$p.lower, res$p.upper)*2
  expect_equal(p[,3,7], unname(tmp[,4]))
})
