# Smoke check only: small Monte Carlo example, not a coverage/FDR guarantee.
library(ALDEx3)
set.seed(307)
results <- list()
timing <- system.time(for (replicate in 1:10) {
  N <- 24; D <- 40
  treatment <- rep(0:1,each=N/2)
  truth <- c(rep(1,4),rep(0,D-4))
  biological <- matrix(rnorm(D*N,sd=0.3),D,N)+truth%o%treatment
  composition <- 2^biological
  composition <- sweep(composition,2,colSums(composition),`/`)
  Y <- round(composition*1e7)
  fit <- aldex(Y,~treatment,data.frame(treatment),nsample=100,
    test='ssrv.bootstrap',scale=compose.sm(gaussian.sm('(Intercept)',sd=0),
      sparse.sm('treatment','bootstrap',bandwidth=0.2)),return.pars='posterior')
  p <- fit$posterior
  results[[replicate]] <- c(coverage=mean(p$lower[2,]<=truth & p$upper[2,]>=truth),
    width=mean(p$upper[2,]-p$lower[2,]),
    null_exclusion=mean(p$lower[2,-(1:4)]>0 | p$upper[2,-(1:4)]<0),
    nonnull_exclusion=mean(p$lower[2,1:4]>0 | p$upper[2,1:4]<0))
})
print(timing)
print(colMeans(do.call(rbind,results)))
cat('10 independent simulations; D=40 N=24 S=100; high depth, bootstrap scale, fixed h=0.2\n')
