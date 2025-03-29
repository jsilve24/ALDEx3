test_that("aldex runs", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  X <- formula(~condition)
  data <- data.frame(condition=condition)
  nsample <- 2000
  foo <- aldex(Y, X, data, nsample=nsample, GAMMA=default)
  expect_true(TRUE)
})

test_that("test-aldex.R correct mean estimate",{
  set.seed(4985)
  sim <- aldex.lm.sim.clr(N=1000, depth=10000)
  res <- aldex(sim$Y, sim$X, GAMMA=default)
  mean.estimate <- apply(res$estimate, c(1,2), FUN=`mean`)
  expect_equal(unname(mean.estimate), unname(sim$Lambda), tolerance=0.05)
})


test_that("test-aldex.R errors with too small streaming",{
  set.seed(4985)
  sim <- aldex.lm.sim.clr(N=1000, depth=10000)
  expect_error(aldex(sim$Y, sim$X, GAMMA=default, streamsize=0.001))
})


test_that("test-aldex.R correct mean estimate when streaming",{
  set.seed(4985)
  sim <- aldex.lm.sim.clr(N=1000, depth=10000)
  res <- aldex(sim$Y, sim$X, GAMMA=default, streamsize=10)
  mean.estimate <- apply(res$estimate, c(1,2), FUN=`mean`)
  expect_equal(unname(mean.estimate), unname(sim$Lambda), tolerance=0.05)
})

test_that("test-aldex.R gives similar results to ALDEx2's aldex.glm", {
  set.seed(4985)
  ## Sim params
  mc.samples <- 6000
  N <- 80
  D <- 40
  DE <- 16
  disease <- c(rep(0, N/2), rep(1, N/2))
  metadata <- cbind(disease)
  ## Simulation
  sim_A <- matrix(0, nrow=D, ncol=N)
  DE_taxa <- sample(1:D, DE)
  lfcs <- rep(0, D)
  lfcs[DE_taxa] <- c(rnorm(DE-3, 1.5, 1), rnorm(3, -1.5, 1))
  for(d in 1:D) {
    for(n in 1:N) {
      sim_A[d,n] <- rnorm(1,5,2)+rnorm(1, metadata[n,1]*lfcs[d], 1)
    }
  }
  sim_A <- 2^sim_A
  sim_Y <- apply(sim_A, 2, function(col) rmultinom(
                                           1, 10000, col/sum(col)))
  gamma_func <- function(X, Y, logWpara) {
    z <- replicate(mc.samples, {
      rnorm(1, 0.5, 0.5)*X[2,]
    })
    return(z)
  }
  aldex3.res <- aldex(sim_Y, t(cbind(1, metadata)),
                          nsample=mc.samples,
                          GAMMA=gamma_func)
  ## Generated with this code
  ##  gamma <- replicate(mc.samples, {
  ##    rnorm(1, 0.5, 0.5)*metadata[,1]
  ##  })
  ##  gamma <- 2^gamma
  ##  glm_meta <- cbind(disease)
  ##  aldex.obj <- aldex.clr(sim_Y, glm_meta, mc.samples=mc.samples,
  ##                         gamma=gamma)
  ##  aldex.res <- aldex.glm(aldex.obj, fdr.method="BH")
  aldex2.adj.pvals <- c(0.8327, 0.5044, 0.0000, 0.0545, 0.0302, 0.5443,
                        0.0279, 0.2994, 0.5235, 0.6318, 0.6858, 0.5477,
                        0.5605, 0.5571, 0.4444, 0.3122, 0.4868, 0.4553,
                        0.2407, 0.6137, 0.6312, 0.1409, 0.8608, 0.3610,
                        0.3114, 0.5347, 0.5078, 0.0000, 0.5888, 0.2307,
                        0.6040, 0.8256, 0.7460, 0.0726, 0.8712, 0.6427,
                        0.3661, 0.8184, 0.0761, 0.7657)
  aldex2.pvals <- c(0.7450, 0.2391, 0.0000, 0.0026, 0.0037, 0.3304, 0.0032,
                    0.0814, 0.3056, 0.4665, 0.5694, 0.3265, 0.3492, 0.3351,
                    0.1776, 0.0796, 0.2401, 0.1963, 0.0448, 0.4219, 0.4711,
                    0.0288, 0.7173, 0.1124, 0.0856, 0.3076, 0.1840, 0.0000,
                    0.4145, 0.0315, 0.4228, 0.7550, 0.6594, 0.0053, 0.6207,
                    0.4911, 0.1164, 0.7438, 0.0123, 0.662)
   expect_true(all(abs(round(aldex3.res$p.val[2,], 4)-aldex2.pvals)<0.015))
   expect_true(all(abs(round(aldex3.res$p.val.adj[2,], 4)-aldex2.adj.pvals)<0.015))
})


test_that("aldex retrns posterior samples when it should", {
  sim <- aldex.lm.sim.clr(N=10, depth=100)
  res <- aldex(sim$Y, sim$X, GAMMA=default,
               return.pars=c("X", "estimate", "std.error", "p.val",
                                "p.val.adj", "logWpara"))
  expect_true("logWpara" %in% names(res))
  expect_true("estimate" %in% names(res))
  expect_true("std.error" %in% names(res))
  expect_true("p.val" %in% names(res))
  expect_true("p.val.adj" %in% names(res))
  expect_false("logWperp" %in% names(res))
})
