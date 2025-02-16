##' Simulation function soley for testing and exploring ALDEx3, Truth is in CLR
##' Coordinates.
##'
##' Not designed to create realistic data. Does not add any noise to linear
##' regression! True W is in CLR Corrdinates
##'
##' @param D number of taxa/genes
##' @param N number of samples
##' @param P number of covariates
##' @param depth sum of counts for each multinomial draw
##' @return a list with elements Y, X, W, and Lambda
##' @author Justin Silverman
aldex.lm.sim.clr <- function(D=10, N=11, P=2, depth=10000) {
  ## D=10; N=11; P=2; depth=10000
  X <- matrix(runif(P*N, -2, 2), P, N)
  Lambda <- matrix(rnorm(P*D),P, D) # TODO fix so it has correct dimensions (D x
                                    # P, which seems more sensical)
  Lambda <- sweep(Lambda, 1, rowMeans(Lambda), FUN=`-`)

  W <- t(Lambda) %*% X
  pi <- miniclo(2^W)  

  Y <- matrix(NA, D, N)
  for (n in 1:N) {
    Y[,n] <- rmultinom(1, 1000, pi[,n])
  }

  return(list(Y=Y, X=X, W=W, Lambda=Lambda))
}
