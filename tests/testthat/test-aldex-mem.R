test_that("sr.mem expected errors", {
  set.seed(4985)
  nsample <- 10
  sim <- aldex.mem.sim(10, 10, 20, 100000, 1, F, 0, 0, 0.1)
  expect_error(aldex(sim$Y, sim$meta, data=sim$meta, method="lme4",
          nsample=nsample, gamma=0.5, scale=clr.sm,
          n.cores=1), "X should be a mixed effects formula")
  expect_error(aldex(sim$Y, sim$meta, data=NULL, method="lme4",
          nsample=nsample, gamma=0.5, scale=clr.sm,
          n.cores=1), "data should not be null")
  expect_error(aldex(sim$Y, ~treatment, data=sim$meta,
                     method="nlme", scale=clr.sm, gamma=0.5,
                     nsample=nsample, n.cores=1),
               "If method=\"nlme\", arg \"random\" must be provided")
})

test_that("aldex mem runs parallel", {
  set.seed(56841)
  Y <- matrix(1:20, 10, 12)
  formula <- ~ condition + (1|subjects)
  condition <- c(rep(0, 6), rep(1, 6))
  subjects <- rep(1:6, each=2)
  data <- data.frame(condition=condition,
                     subjects=subjects)
  nsample <- 10
  foo <- aldex(Y, formula, data=data, method="lme4",
          nsample=nsample, gamma=0.5, scale=clr.sm,
          n.cores=2)
  expect_true(TRUE)
  expect_equal(class(foo), "aldex")
  expect_equal(rownames(foo$logComp), paste0("entity_", 1:10))
})

test_that("aldex mem lme4/nlme correct naming", {
  set.seed(6841)
  sim <- aldex.mem.sim(10, 10, 10, 100000, 1, FALSE, 0, 0, 0.25)
  Y <- sim$Y
  Y <- provideDimnames(Y)
  meta <- sim$meta
  nsample <- 10
  foo <- aldex(Y, ~ treatment + (1|subject_ids), data=meta,
               method="lme4", nsample=nsample, gamma=0.5,
               scale=clr.sm, n.cores=1,
               return.pars=c("X", "estimate", "std.error", "p.val",
                             "p.val.adj", "logComp", "logScale",
                             "random.effects"))
  expect_equal(colnames(foo$estimate), row.names(Y))
  expect_equal(row.names(foo$random.effects), c("subject_ids.(Intercept)",
                                                "Residual"))
  foo <- aldex(Y, ~treatment, data=meta, method="nlme",
               nsample=nsample, gamma=0.5, scale=clr.sm,
               n.cores=1, random=~1|subject_ids,
               return.pars=c("X", "estimate", "std.error", "p.val",
                             "p.val.adj", "logComp", "logScale",
                             "random.effects"))
  expect_equal(colnames(foo$estimate), row.names(Y))
  expect_equal(row.names(foo$random.effects), c("(Intercept)",
                                                "Residual"))
})


test_that("aldex mem nlme correct results", {
  set.seed(43262)
  nsample <- 2
  sim <- aldex.mem.sim(20, 24, 200, 1000000, 1, F, 0,
                       0.5, 0.25)
  true.S <- log2(colSums(sim$W))
  custom.logScale <- replicate(nsample, true.S)
  aldex.res <- aldex(sim$Y, ~treatment, random=~1|subject_ids,
                     data=sim$meta, method="nlme", n.cores=1,
                     nsample=nsample, scale=custom.logScale,
                     correlation=corAR1(form=~time|subject_ids),
                     return.pars=c("X", "estimate", "std.error", "p.val",
                                   "p.val.adj", "logComp", "logScale",
                                   "random.effects"))
  mean.estimate <- apply(aldex.res$estimate, c(1,2), FUN=`mean`)
  expect_equal(mean.estimate[2,], sim$fixed_effs[,2], tolerance=0.1)
  random_effs <- rowMeans(
    apply(aldex.res$random.effects, c(1,2), FUN=`mean`))
  expect_equal(unname(random_effs), c(1.5^2, 0.1^2),
               tolerance=0.15)
})

test_that("aldex mem lme4 correct results", {
  set.seed(43262)
  nsample <- 2
  sim <- aldex.mem.sim(20, 24, 200, 1000000, 50, F, 0,
                       0, 0.1)
  true.S <- log2(colSums(sim$W))
  custom.logScale <- replicate(nsample, true.S)
  aldex.res <- aldex(sim$Y,
                     ~treatment+(1|subject_ids)+(1|location_ids),
                     data=sim$meta, method="lme4", n.cores=1,
                     nsample=nsample, scale=custom.logScale,
                     return.pars=c("X", "estimate", "std.error", "p.val",
                                   "p.val.adj", "logComp", "logScale",
                                   "random.effects"))
  mean.estimate <- apply(aldex.res$estimate, c(1,2), FUN=`mean`)
  expect_equal(mean.estimate[2,], sim$fixed_effs[,2], tolerance=0.1)
  random_effs <- rowMeans(
    apply(aldex.res$random.effects, c(1,2), FUN=`mean`))
  expect_equal(names(random_effs), c("subject_ids.(Intercept)",
                 "location_ids.(Intercept)", "Residual"))
  expect_equal(unname(random_effs), c(1.5^2, 0.75^2, 0.1^2),
               tolerance=0.15)
})
