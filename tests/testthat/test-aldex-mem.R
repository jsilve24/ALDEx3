test_that("aldex mem runs", {
  Y <- matrix(1:120, 10, 12)
  formula <- ~ condition + (1|subjects)
  condition <- c(rep(0, 6), rep(1, 6))
  subjects <- rep(1:6, each=2)
  data <- data.frame(condition=condition,
                     subjects=subjects)
  nsample <- 10
  foo <- aldex(Y, formula, data=data, method="lme4",
          nsample=nsample, gamma=0.5, scale=clr.sm,
          n.cores=1)
  expect_true(TRUE)
  expect_equal(class(foo), "aldex")
  expect_equal(rownames(foo$logComp), paste0("entity_", 1:10))
})

test_that("aldex mem expected errors", {
  set.seed(4985)
  sim <- aldex.lm.sim.clr(N=1000, depth=10000)
  expect_error(aldex(sim$Y, sim$X, data=sim$X, method="lme4",
          nsample=nsample, gamma=0.5, scale=clr.sm,
          n.cores=1), "X should be a mixed effects formula")
  expect_error(aldex(sim$Y, sim$X, data=NULL, method="lme4",
          nsample=nsample, gamma=0.5, scale=clr.sm,
          n.cores=1), "data should not be null")
})

test_that("aldex mem correct naming", {
  Y <- matrix(1:120, 10, 12)
  Y <- provideDimnames(Y)
  formula <- ~ condition + (1|subjects)
  condition <- c(rep(0, 6), rep(1, 6))
  subjects <- rep(1:6, each=2)
  data <- data.frame(condition=condition,
                     subjects=subjects)
  nsample <- 10
  foo <- aldex(Y, formula, data=data, method="lme4",
               nsample=nsample, gamma=0.5, scale=clr.sm,
               n.cores=2,
               return.pars=c("X", "estimate", "std.error", "p.val",
                             "p.val.adj", "logComp", "logScale",
                             "random.effects"))
  expect_equal(colnames(foo$logComp), colnames(Y))
  expect_equal(row.names(foo$random.effects), c("subjects:(Intercept)",
                                                "Residual:(Intercept)"))
})

test_that("aldex mem correct mean estimate",{
  set.seed(4985)
  sim <- aldex.lm.sim.clr(N=500, depth=10000, subjects=50)
  formula <- ~ 0 + A + B + (1|subjects)
  data <- data.frame(t(sim$X))
  subjects <- as.character(unname(
    apply(sim$Z, 2, function(i) which(i==1))))
  data <- cbind(data, subjects)
  colnames(data) <- c("A", "B", "subjects")
  res <- aldex(sim$Y, formula, data=data, scale=clr.sm, method="lme4",
               nsample=1000)
  mean.estimate <- apply(res$estimate, c(1,2), FUN=`mean`)
  expect_equal(unname(mean.estimate), unname(sim$Lambda), tolerance=0.1)
})


