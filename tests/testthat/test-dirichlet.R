require(rBeta2009)

test_that("dirichlet matches moment of rDirichlet from rBeta2009", {
  set.seed(64841)
  S <- 100000
  alpha <- c(0.1, 2, 3)
  x <- rDirichlet(S, alpha <- c(0.1, 2, 3))
  y <- t(rBeta2009::rdirichlet(S, alpha <- c(0.1, 2, 3)))
  expect_equal(rowMeans(x), rowMeans(y), tolerance=0.01)
  x <- apply(x, 1, var)
  y <- apply(y, 1, var)
  expect_equal(x, y, tolerance=0.05)
})
