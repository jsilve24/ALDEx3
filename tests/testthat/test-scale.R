test_that("sample.sm validates scale moments", {
  X <- matrix(0, 5, 10)
  logComp <- array(0, c(5, 10, 20))
  s.mu <- seq_len(10)

  expect_error(sample.sm(X, logComp, s.sd=rep(0.5, 10)),
               "s.mu cannot be NULL")
  expect_error(sample.sm(X, logComp, s.mu=1:5, s.sd=rep(0.5, 10)),
               "s.mu should be a numeric vector of length 10")
  expect_error(sample.sm(X, logComp, s.mu=s.mu),
               "Exactly one of s.sd")
  expect_error(sample.sm(X, logComp, s.mu=s.mu,
                         s.sd=rep(0.5, 10), s.cov=diag(0.25, 10)),
               "Exactly one of s.sd")
  expect_error(sample.sm(X, logComp, s.mu=s.mu, s.sd=rep(0.5, 5)),
               "s.sd should be a numeric vector of length 10")
  expect_error(sample.sm(X, logComp, s.mu=s.mu, s.sd=c(-1, rep(0, 9))),
               "non-negative")
  expect_error(sample.sm(X, logComp, s.mu=s.mu,
                         s.sd=c(Inf, rep(0, 9))),
               "finite")
  expect_error(sample.sm(X, logComp, s.mu=s.mu, s.cov=diag(0.25, 9)),
               "10x10 covariance matrix")

  asymmetric <- diag(0.25, 10)
  asymmetric[1, 2] <- 0.1
  expect_error(sample.sm(X, logComp, s.mu=s.mu, s.cov=asymmetric),
               "symmetric")

  non_psd <- diag(0.25, 10)
  non_psd[1, 1] <- -0.25
  expect_error(sample.sm(X, logComp, s.mu=s.mu, s.cov=non_psd),
               "positive semidefinite")
})

test_that("sample.sm uses SDs and preserves deprecated argument semantics", {
  X <- matrix(0, 5, 10)
  logComp <- array(0, c(5, 10, 5000))
  s.mu <- seq_len(10)

  set.seed(43643)
  by_sd <- sample.sm(X, logComp, s.mu=s.mu, s.sd=rep(0.5, 10))
  expect_equal(dim(by_sd), c(10, 5000))
  expect_equal(mean(apply(by_sd, 1, var)), 0.25, tolerance=0.01)

  set.seed(43643)
  expect_warning(
    by_variance <- sample.sm(
      X, logComp, s.mu=s.mu, s.var=rep(0.25, 10)
    ),
    "s.var is deprecated"
  )
  expect_identical(by_variance, by_sd)

  covariance <- diag(0.25, 10)
  set.seed(918)
  by_covariance <- sample.sm(X, logComp, s.mu=s.mu, s.cov=covariance)
  set.seed(918)
  expect_warning(
    by_legacy_covariance <- sample.sm(
      X, logComp, s.mu=s.mu, s.cor=covariance
    ),
    "s.cor is deprecated"
  )
  expect_identical(by_legacy_covariance, by_covariance)

  one_draw <- sample.sm(X, array(0, c(5, 10, 1)), s.mu=s.mu,
                        s.cov=covariance)
  expect_equal(dim(one_draw), c(10, 1))
  one_independent_draw <- sample.sm(
    X, array(0, c(5, 10, 1)), s.mu=s.mu, s.sd=rep(0.5, 10)
  )
  expect_equal(dim(one_independent_draw), c(10, 1))
})

test_that("coefficient.sm validates scale moments", {
  X <- rbind(c(1, 1, 1, 1, 1, 1),
             c(0, 0, 0, 1, 1, 1),
             c(0.2, 0.1, 1, 0.8, 1.2, 1.4))
  logComp <- array(0, c(5, 6, 20))
  c.mu <- c(0, 4, 2)

  expect_error(coefficient.sm(X, logComp, c.sd=rep(0.5, 3)),
               "c.mu cannot be NULL")
  expect_error(coefficient.sm(X, logComp, c.mu=1:6, c.sd=rep(0.5, 3)),
               "c.mu should be a numeric vector of length 3")
  expect_error(coefficient.sm(X, logComp, c.mu=c.mu),
               "Exactly one of c.sd")
  expect_error(coefficient.sm(X, logComp, c.mu=c.mu,
                              c.sd=rep(0.5, 3), c.cov=diag(0.25, 3)),
               "Exactly one of c.sd")
  expect_error(coefficient.sm(X, logComp, c.mu=c.mu, c.sd=rep(0.5, 2)),
               "c.sd should be a numeric vector of length 3")
  expect_error(coefficient.sm(X, logComp, c.mu=c.mu,
                              c.cov=matrix(0.1, 3, 2)),
               "3x3 covariance matrix")
})

test_that("coefficient.sm expresses group offsets on the model contrast", {
  X <- rbind("(Intercept)"=rep(1, 6),
             "condsS"=c(0, 0, 0, 1, 1, 1))
  logComp <- array(0, c(5, 6, 5000))

  deterministic <- coefficient.sm(
    X, logComp, c.mu=c(0, 8), c.sd=c(0, 0)
  )
  expect_equal(deterministic[1:3, ], matrix(0, 3, 5000))
  expect_equal(deterministic[4:6, ], matrix(8, 3, 5000))

  set.seed(43643)
  uncertain <- coefficient.sm(
    X, logComp, c.mu=c(0, 8), c.sd=c(0, 0.5)
  )
  contrast <- uncertain[4, ] - uncertain[1, ]
  expect_equal(mean(contrast), 8, tolerance=0.02)
  expect_equal(sd(contrast), 0.5, tolerance=0.02)
  expect_identical(uncertain[1, ], uncertain[2, ])
  expect_identical(uncertain[4, ], uncertain[5, ])

  covariance <- diag(c(0, 0.25))
  set.seed(918)
  by_covariance <- coefficient.sm(
    X, logComp, c.mu=c(0, 8), c.cov=covariance
  )
  set.seed(918)
  expect_warning(
    by_legacy_covariance <- coefficient.sm(
      X, logComp, c.mu=c(0, 8), c.cor=covariance
    ),
    "c.cor is deprecated"
  )
  expect_identical(by_legacy_covariance, by_covariance)
})

test_that("sample- and coefficient-level models retain distinct dependence", {
  X <- rbind(rep(1, 6), c(0, 0, 0, 1, 1, 1))
  logComp <- array(0, c(5, 6, 100))

  set.seed(123)
  sample_scale <- sample.sm(
    X, logComp, s.mu=rep(0, 6), s.sd=rep(0.5, 6)
  )
  set.seed(123)
  coefficient_scale <- tss.sm(X, logComp, gamma=0.5)

  expect_false(isTRUE(all.equal(sample_scale[1, ], sample_scale[2, ])))
  expect_identical(coefficient_scale[1, ], coefficient_scale[2, ])
  expect_identical(coefficient_scale[4, ], coefficient_scale[5, ])
})

test_that("gamma is a validated standard deviation", {
  X <- matrix(1, 1, 4)
  logComp <- array(0, c(5, 4, 5000))

  set.seed(99)
  tss_scale <- tss.sm(X, logComp, gamma=0.5)
  expect_equal(sd(tss_scale[1, ]), 0.5, tolerance=0.02)

  set.seed(99)
  clr_scale <- clr.sm(X, logComp, gamma=0.5)
  expect_equal(sd(clr_scale[1, ]), 0.5, tolerance=0.02)

  expect_error(tss.sm(X, logComp, gamma=-0.5), "non-negative")
  expect_error(clr.sm(X, logComp, gamma=c(0.5, 1)),
               "numeric vector of length 1")
  expect_error(tss.sm(X, logComp, gamma=Inf), "finite")
})

test_that("aldex warns once for deprecated scale arguments when streaming", {
  Y <- matrix(c(10, 12, 11, 13,
                20, 18, 22, 19), nrow=2, byrow=TRUE)
  condition <- factor(c("N", "N", "S", "S"), levels=c("N", "S"))
  warnings <- character()

  withCallingHandlers(
    aldex(Y, ~condition, data=data.frame(condition), nsample=20,
          scale=sample.sm, s.mu=rep(0, 4), s.var=rep(0.25, 4),
          streamsize=0.003),
    warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(grep("s.var is deprecated", warnings), 1)
})
