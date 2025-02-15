test_that("aldex.lm runs", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  X <- formula(~condition)
  data <- data.frame(condition=condition)
  nsample <- 2000
  aldex.lm(Y, X, data, nsample=nsample, gamma=default)
  expect_true(TRUE)
})

test_that("test-aldex_lm.R",{
  set.seed(4985)
  sim <- aldex.lm.sim.clr(N=1000, depth=10000)
  res <- aldex.lm(Y, X, gamma=default)
  mean.estimate <- apply(res$estimate, c(1,2), FUN=`mean`)
  expect_equal(mean.estimate, Lambda, tolerance=0.05)
})
