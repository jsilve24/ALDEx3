test_that("miniclo correct", {
  x <- matrix(1:6, ncol=2)
  y <- miniclo(x)
  expect_equal(x[,1]/sum(x[,1]), y[,1])
  expect_equal(x[,2]/sum(x[,2]), y[,2])
})

test_that("center correct", {
  set.seed(165548)
  x <- array(rnorm(24), c(3, 4, 5))
  expect_equal(center(x)[,2,3], x[,2,3]-mean(x[,2,3]))
})

