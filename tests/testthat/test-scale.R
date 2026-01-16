test_that("sample.sm expected error thrown", {
  X <- matrix(0, 5, 10)
  logComp <- array(0, c(5, 10, 2000))
  s.mu <- 1:10
  s.var <- rep(0.1, 10)
  s.var.wrong <- rep(0.1, 5)
  s.mu.wrong <- 1:20
  s.cor <- diag(0.5, 10)
  s.cor.wrong <- matrix(0.1, 11, 10)

  expect_error(sample.sm(X, logComp), "s.mu cannot be NULL")
  expect_error(sample.sm(X, logComp, s.mu=s.mu.wrong,
                         s.var=s.var),
               "s.mu should have same length as")
  expect_error(sample.sm(X, logComp, s.mu=s.mu),
               "One of s.var and s.cor should be NULL")
  expect_error(sample.sm(X, logComp, s.mu=s.mu,
                         s.var=s.var, s.cor=s.cor),
               "One of s.var and s.cor should be NULL")
  expect_error(sample.sm(X, logComp, s.mu=s.mu,
                         s.var=s.var.wrong),
               "s.var should have same length as")
  expect_error(sample.sm(X, logComp, s.mu=s.mu,
                         s.cor=s.cor.wrong),
               "s.cor should be an NxN matrix where N")
})

test_that("sample.sm expected output", {
  set.seed(43643)
  X <- matrix(0, 5, 10)
  logComp <- array(0, c(5, 10, 5000))
  s.mu <- 1:10
  s.var <- rep(0.1, 10)
  s.cor <- diag(0.5, 10)

  res_1 <- sample.sm(X, logComp, s.mu=s.mu, s.var=s.var)
  expect_setequal(dim(res_1), c(10, 5000))
  expect_equal(mean(apply(res_1, 1, var)), 0.1,
               tolerance=0.01)

  res_2 <- sample.sm(X, logComp, s.mu=s.mu, s.cor=s.cor)
  expect_setequal(dim(res_2), c(10, 5000))
  expect_equal(mean(apply(res_2, 1, var)), 0.5,
               tolerance=0.01)
})

test_that("coefficient.sm expected error thrown", {
  X <- rbind(c(1, 1, 1, 1, 1, 1),
             c(0, 0, 0, 1, 1, 1),
             c(0.2, 0.1, 1, 0.8, 1.2, 1.4))
  logComp <- array(0, c(5, 6, 2000))
  c.mu <- c(0, 4, 2)
  c.mu.wrong <- 1:6
  c.cor <- diag(0.25, 3)
  c.cor.wrong <- diag(0.25, 6)

  expect_error(coefficient.sm(X, logComp), "c.mu cannot be NULL")
  expect_error(coefficient.sm(X, logComp, c.mu=c.mu),
               "c.cor cannot be NULL")
  expect_error(coefficient.sm(X, logComp, c.mu=c.mu.wrong, c.cor=c.cor),
                       "c.mu should have length of P")
  expect_error(coefficient.sm(X, logComp, c.mu=c.mu, c.cor=c.cor.wrong),
                       "c.cor should be a PxP")
})

test_that("sample.sm expected output", {
  set.seed(43643)
  X <- rbind(c(1, 1, 1, 1, 1, 1),
             c(0, 0, 0, 1, 1, 1),
             c(0.9, 1.2, 1, 0.8, 1.2, 1.3))
  logComp <- array(0, c(5, 6, 5000))
  c.mu <- c(0, 4, 2)
  c.mu.wrong <- 1:6
  c.cor <- diag(0.25, 3)
  c.cor.wrong <- diag(0.25, 6)

  res_1 <- coefficient.sm(X, logComp,
                   c.mu=c.mu, c.cor=c.cor)
  expect_setequal(dim(res_1), c(6, 5000))
  expect_equal(mean(apply(res_1, 2, function(col) {
    mean(col[4:6])-mean(col[1:3])
  })), 4, tolerance=0.25)
})
