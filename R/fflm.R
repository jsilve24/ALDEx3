## Freaking Fast Linear Models for ALDEx2
## Author: Justin D. Silverman 
## Email: justinsilverman@psu.edu

##' Freaking Fask Linear Models
##'
##' Tailored for ALDEx2 where covariates are shared between massive numbers of
##' linear regressions where only Y is changing.
##'
##' @param Y a numeric array (N x D X S) where D is the number of taxa/genes, N
##'   is the number of samples, and S is the number of posterior samples
##' @param X a numeric matrix (N x P) where P is number of covariates
##' @return A list of (P x D x S)-arrays with the OLS point estimates, the
##'   standard errors, and the two-sided p-values for each coefficient (P), of
##'   each model fit to each taxa (D) and each posterior sample (S)
##' @author Justin Silverman
fflm <- function(Y, X) {
   N <- dim(Y)[1]
   D <- dim(Y)[2]
   S <- dim(Y)[3]
   P <- dim(X)[2]
   stopifnot(dim(X)[1] == N)
  tmp <- 
   stopifnot (P <= N)

   ## from 3D array to 2D as a hack to make this much faster
   ## Will make Y, N x PS 
   Y <- array(Y, c(N, D*S))
   ## Y <- Y[,,1]
   ## Now compute OLS solution 
   A <- t(X) %*% X
   b <- t(X) %*% Y
   Theta <- solve(A, b)
   sqerr <- (Y - X %*% Theta)^2
   sigmaSq <- colSums(sqerr)/(N-P)
   dvcov <- diag( chol2inv(chol(t(X) %*% X)))
   stderr <- sqrt(dvcov %*% t(sigmaSq))
   t <- Theta / stderr
   dof <- N-P
   p.lower <- pt(t, dof)
   p.upper <- pt(-t, dof)
   return(list(estimate = array(Theta, c(P, D, S)),
               std.error = array(stderr, c(P, D, S)),
               p.lower = array(p.lower, c(P, D, S)),
               p.upper = array(p.upper, c(P, D, S))))
}
