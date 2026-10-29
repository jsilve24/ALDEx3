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
##' @param test (default t.HC3), "t", t test is performed for each covariate
##'   (fast); "t.HC0" Heteroskedasticsity-Robust Standard Errors used (HC0;
##'   White's; slower); "t.HC3" (default) Heteroskedasticsity-Robust Standard
##'   Errors used (HC3; unlike HC0, this includes a leverage adjustment and is
##'   better for small sample sizes or when there are data with high leverage;
##'   slowest). To learn more about these, loko at Long and Ervin (2000) Using
##'   Heteroscedasticity Consistent Standard Errors in the Linear Regression
##'   Model, The American Statistician. 
##' @param contrasts optional numeric matrix with named coefficient columns and
##'   one contrast per row.
##' @return A list of (P x D x S)-arrays with the OLS point estimates, standard
##'   errors, degrees of freedom, and tail probabilities. If requested,
##'   contrast arrays have dimensions K x D x S.
##' @author Justin Silverman
##' @keywords internal
fflm <- function(Y, X, test="t.HC3", contrasts=NULL) {
   N <- dim(Y)[1]
   D <- dim(Y)[2]
   S <- dim(Y)[3]
   P <- dim(X)[2]
   stopifnot(dim(X)[1] == N)
   if (!is.character(test) || length(test) != 1L || is.na(test) ||
       !(test %in% c("t", "t.HC0", "t.HC3"))) {
     stop("test must be one of 't', 't.HC0', or 't.HC3'")
   }
   if (P > N || qr(X)$rank < P) {
     stop("rank-deficient designs are not supported by fflm")
   }
   contrasts <- .aldex_validate_contrasts(contrasts, X)

   ## from 3D array to 2D as a hack to make this much faster
   ## Will make Y, N x PS 
   Y <- array(Y, c(N, D*S))
   ## Y <- Y[,,1]
   ## Now compute OLS solution 
   A <- t(X) %*% X
   b <- t(X) %*% Y
   Theta <- solve(A, b)
   sqerr <- (Y - X %*% Theta)^2
   inverseXtX <- chol2inv(chol(A))

   if (test=="t") {
     sigmaSq <- colSums(sqerr)/(N-P)
     dvcov <- diag(inverseXtX)
     stderr <- sqrt(dvcov %*% t(sigmaSq))
   } else if (test=="t.HC0") {
     ## use white's robust standard errors instead of the more common standard
     ## errors implemented above.
     ## older -- slower, easier to read
     ## XtXinvXt2 <- (chol2inv(chol(t(X) %*% X)) %*% t(X))^2
     ## stderr <- matrix(0, P, D*S)
     ## for (i in 1:(D*S)) {
     ##   stderr[,i] <- sqrt(rowSums(sweep(XtXinvXt2, 2, sqerr[,i], FUN=`*`)))
     ## }

     ## optimizes / vectorized
     # Precompute XtXinvXt squared
     XtXinvXt2 <- (inverseXtX %*% t(X))^2  # elementwise square
     # Matrix multiplication instead of looping
     stderr <- sqrt(XtXinvXt2 %*% sqerr)
   } else if (test=="t.HC3") {
     ## old version 
     ## XtXinvXt <- chol2inv(chol(t(X) %*% X)) %*% t(X)
     ## leverage.correction <- (1-diag(X %*% XtXinvXt))^2
     ## XtXinvXt <- XtXinvXt^2
     ## stderr <- matrix(0, P, D*S)
     ## for (i in 1:(D*S)) {
     ##   stderr[,i] <- sqrt(rowSums(sweep(XtXinvXt, 2,
     ##                                    sqerr[,i]/leverage.correction,FUN=`*`)))
     ## }

     # vectorized version
     XtXinvXt <- inverseXtX %*% t(X)
     leverage.correction <- (1 - diag(X %*% XtXinvXt))^2
     XtXinvXt_sq <- XtXinvXt^2  # elementwise square, reused

     # Vectorized computation
     # Divide sqerr by leverage correction (broadcasting over rows)
     adjusted_sqerr <- sweep(sqerr, 1, leverage.correction, `/`)

     # Multiply XtXinvXt_sq %*% adjusted_sqerr and take sqrt
     stderr <- sqrt(XtXinvXt_sq %*% adjusted_sqerr)
   }
   t <- Theta / stderr
   dof <- N-P
   p.upper <- pt(t, dof, lower.tail=FALSE)
   p.lower <- pt(t, dof, lower.tail=TRUE)
   out <- list(estimate = array(Theta, c(P, D, S)),
               std.error = array(stderr, c(P, D, S)),
               df = array(dof, c(P, D, S)),
               p.lower = array(p.lower, c(P, D, S)),
               p.upper = array(p.upper, c(P, D, S)))

   if (!is.null(contrasts)) {
     K <- nrow(contrasts)
     contrast.estimate <- contrasts %*% Theta
     if (test == "t") {
       contrast.variance <- rowSums((contrasts %*% inverseXtX) * contrasts)
       contrast.stderr <- sqrt(contrast.variance %*% t(sigmaSq))
     } else {
       ## This projects the full sandwich covariance to C beta.  Combining
       ## marginal coefficient SEs would discard the covariance terms.
       projection <- contrasts %*% inverseXtX %*% t(X)
       if (test == "t.HC3") {
         contrast.stderr <- sqrt(projection^2 %*%
                                   sweep(sqerr, 1, leverage.correction, `/`))
       } else {
         contrast.stderr <- sqrt(projection^2 %*% sqerr)
       }
     }
     contrast.t <- contrast.estimate / contrast.stderr
     out$contrast.estimate <- array(contrast.estimate, c(K, D, S),
                                    dimnames = list(rownames(contrasts), NULL, NULL))
     out$contrast.std.error <- array(contrast.stderr, c(K, D, S),
                                     dimnames = list(rownames(contrasts), NULL, NULL))
     out$contrast.df <- array(dof, c(K, D, S),
                              dimnames = list(rownames(contrasts), NULL, NULL))
     out$contrast.p.lower <- array(pt(contrast.t, dof, lower.tail=TRUE), c(K, D, S),
                                   dimnames = list(rownames(contrasts), NULL, NULL))
     out$contrast.p.upper <- array(pt(contrast.t, dof, lower.tail=FALSE), c(K, D, S),
                                   dimnames = list(rownames(contrasts), NULL, NULL))
   }
   out
}
