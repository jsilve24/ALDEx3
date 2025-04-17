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
##' @param subjects the number of subjects if simulating random
##'   effects 0 means no random effects, must go into N
##' @return a list with elements Y, X, W, and Lambda
##' @author Justin Silverman, Kyle McGovern
aldex.lm.sim.clr <- function(D=10, N=11, P=2, depth=10000, subjects=0) {
  ## D=10; N=11; P=2; depth=10000
  X <- matrix(runif(P*N, -2, 2), P, N)
  Lambda <- matrix(rnorm(P*D),P, D) # TODO fix so it has correct dimensions (D x
                                    # P, which seems more sensical)
  Lambda <- sweep(Lambda, 1, rowMeans(Lambda), FUN=`-`)
  if(subjects != 0) {
    if((N%%subjects)!=0) stop("N must be divisible by # subjects to")
    sample_id <- rep(1:subjects, each=N/subjects)
    Z <- as.matrix(model.matrix(~ 0 + factor(sample_id)))
    colnames(Z) <- NULL
    u <- replicate(D, rnorm(subjects, 0, 1))
    W <- t(Lambda) %*% X + t(Z %*% u)
    Z <- t(Z)
  } else {
    Z <- NULL
    u <- NULL
    W <- t(Lambda) %*% X
  }

  pi <- miniclo(2^W)

  Y <- matrix(NA, D, N)
  for (n in 1:N) {
    Y[,n] <- rmultinom(1, 1000, pi[,n])
  }

  return(list(Y=Y, X=X, W=W, Lambda=Lambda, Z=Z, u=u))
}
