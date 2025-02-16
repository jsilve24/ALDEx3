##' ALDEx3 Linear Models
##'
##' 
##' @title ALDEx3 Linear Modles
##' @param counts an (D x N) matrix of sequence count, N is number of samples, D
##'   is number of taxa or genes
##' @param X either a formula (in which case DATA must be non-null) or a model
##'   matrix of dimension P x N (P is number of linear model covariates). If a
##'   formula is passed it should not include the target Y, e.g., should simply
##'   be "~condition-1" (note the lack of the left hand side).
##' @param data a data frame for use with formula, must have N rows
##' @param nsample number of monte carlo replicates
##' @param gamma the scale model, can be a function or an N x nsample matrix. If
##'   a function is passed, it should take one argument (pi) which is a N x D x
##'   nsample array of log-transformed relative abundances. That function must
##'   in turn output an N x nsample matrix of scale samples on the log2 scale.
##' @return List with elements estimate (P x D x nsample array), std.error (P x
##'   D x nsample array), and p.val (P x D matrix) summarizing over the
##'   posterior. TODO p.value calcluation may be slightly different than in
##'   current ALDEx3 -- need to check.
##' @export
##' @author Justin Silverman
aldex.lm <- function(Y, X, data=NULL, nsample=2000,  gamma=NULL) {
  N <- ncol(Y)
  D <- nrow(Y)

  ## dirichlet sample
  logWpara <- array(NA, c(D, N, nsample))
  for (i in 1:N) {
    logWpara[,i,] <- log2(rDirichlet(nsample, Y[,i]+0.5))
  }

  ## compute model matrix
  if (inherits(X, "formula")) {
    if (is.null(data)) stop("data should not be null if X is a formula")
    X <- t(model.matrix(X, data))
  } 

  ## sample from scale model
  if (is.null(gamma)){
    stop("You probably want a scale model :)")
  } else if (is.function(gamma)) {
    logWperp <- gamma(X, Y, logWpara)
  } 

  ## compute scaled abundances (W)
  logW <- sweep(logWpara, c(2,3), logWperp, FUN=`+`)

  ## fit linear model
  out <- fflm(aperm(logW, c(2,1,3)),t(X)) # TODO change fflm so it gives correct
                                          # output dimensions and takes correct
                                          # inputs dimensions

  ## summarise output
  p.lower <- apply(out$p.lower,c(1,2), FUN=`mean`)
  p.upper <- apply(out$p.upper,c(1,2), FUN=`mean`)
  p <- 2*pmin(p.lower, p.upper) # TODO this is different than in Beyond
                                # Normalizations
  return(list(estimate=out$estimate, std.error=out$std.error, p.val=p))
}
