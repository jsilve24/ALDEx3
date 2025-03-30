## Freaking Fast Linear Models for ALDEx3
## Author: Justin D. Silverman 
## Email: justinsilverman@psu.edu

##' Freaking Fask Linear Models
##'
##' Tailored for ALDEx3 where covariates are shared between massive numbers of
##' linear regressions where only Y is changing.
##'
##' @param Y a numeric array (N x D X S) where D is the number of taxa/genes, N
##'   is the number of samples, and S is the number of posterior samples
##' @param X a numeric matrix (N x P) where P is number of covariates
##' @return A list of (P x D x S)-arrays with the OLS point estimates, the
##'   standard errors, and the two-sided p-values (using White's robust standard
##'   errors) for each coefficient (P), of each model fit to each taxa (D) and
##'   each posterior sample (S)
##' @author Justin Silverman
fflm <- function(Y, X) {
  N <- dim(Y)[1]
  D <- dim(Y)[2]
  S <- dim(Y)[3]
  P <- dim(X)[2]
  stopifnot(dim(X)[1] == N)
  stopifnot (P <= N)

  Y <- array(Y, c(N, D * S))

  XtX_inv <- chol2inv(chol(crossprod(X)))
  Theta <- XtX_inv %*% crossprod(X, Y)

  residuals <- Y - X %*% Theta

  hat_matrix <- X %*% XtX_inv
  robust_var <- matrix(0, nrow = P, ncol = D*S)

  for (i in 1:(D*S)) {
    X_resid <- hat_matrix * residuals[, i]
    robust_var[, i] <- colSums(X_resid^2)
  }

  stderr <- sqrt(robust_var)
  t_values <- Theta / stderr
  
  dof <- N - P
  p.upper <- 2 * pt(-abs(t_values), df = dof)
  p.lower <- 1 - p.upper

  return(list(
    estimate = array(Theta, c(P, D, S)),
    std.error = array(stderr, c(P, D, S)),
    p.lower = array(p.lower, c(P, D, S)),
    p.upper = array(p.upper, c(P, D, S))
  ))
}
