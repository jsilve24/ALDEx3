test_that("cohensd runs", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  X <- formula(~condition)
  data <- data.frame(condition=condition)
  nsample <- 2000
  foo <- aldex(Y, X, data, nsample=nsample, scale=clr)
  a <- cohensd(foo, condition)
  b <- cohensd(foo, 2)
  expect_equal(a, b)
  expect_true(!any(is.na(a)))
})


test_that("cohensd correct", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  X <- formula(~condition)
  data <- data.frame(condition=condition)
  nsample <- 2000
  foo <- aldex(Y, X, data, nsample=nsample, scale=clr)
  a <- cohensd(foo, condition)
  ## check a[2,5] for brevity
  s <- 5
  d <- 2
  logW <- sweep(foo$logComp, c(2,3), foo$logScale, FUN=`+`)[2,,5]
  x0 <- logW[condition==0]
  x1 <- logW[condition==1]
  mean0 <- mean(x0)
  mean1 <- mean(x1)
  n0 <- sum(condition==0)
  n1 <- sum(condition==1)
  denom <- sqrt(((n0-1)*var(x0) + (n1-1)*var(x1))/(n0+n1-2))
  expect_equal(unname(a[2,5]), (mean1-mean0)/denom)
})
