##' Implementation of SR-MEM: scale reliant mixed effects models.
##'
##' @param logW a numeric array (N x D X S) where D is the number of taxa,
##'   N is the number of samples, and S is the number of posterior samples
##' @param formula an lme4 formula with fixed and random effects for lmer
##' @param data A data.frame with the random and fixed effects items in formula
##' @param n.cores The number of cores for parallelization. If n.cores=1, no
##'   parallelization is used
##' @return A list of (P x D x S)-arrays with the fixed effect point estimates,
##'   the standard errors, degrees of freedom and the lower and upper p-values
##'   for each coefficient (P), of each model fit to each taxa (D) and each
##'   posterior sample (S)
##' @importFrom lmerTest lmer
##' @importFrom lme4 VarCorr
##' @importFrom parallel parLapply makeCluster stopCluster clusterExport
##' @author Kyle McGovern
aldex.mem <- function(logW, formula, data, n.cores) {
  N <- dim(logW)[1]
  D <- dim(logW)[2]
  S <- dim(logW)[3]

  logWm <- array(logW, c(N, D*S))

  ## Only start cluster if n.cores > 1
  if(n.cores>1) {
    cl <- makeCluster(n.cores)
    on.exit(stopCluster(cl), add=TRUE)

    ## Trick to avoid code repeat
    lapply_func <- function(X, FUN) parLapply(cl, X, FUN)
  } else {
    lapply_func <- lapply
  }

  ## Process Cols (i.e., per-taxa)
  results <- lapply_func(1:ncol(logWm), function(j) {
    y <- logWm[,j]
    data_temp <- data
    data_temp$y <- y
    fit <- suppressMessages(lmer(update(formula, y~.), data=data_temp))

    ## Fixed Effects
    coefs <- coef(summary(fit))[,1:3]
    coefs <- cbind(coefs, pt(coefs[,1]/coefs[,2], coefs[,3], lower.tail=T))
    coefs <- cbind(coefs, 1-coefs[,4])
    colnames(coefs) <- c("estimate", "std.error", "df", "p.lower", "p.upper")

    ## Random Effects
    re_df <- as.data.frame(lme4::VarCorr(fit))
    re_df$var1[is.na(re_df$var1)] <- "(Intercept)"
    re_df$var2[is.na(re_df$var2)] <- "(Intercept)"
    re_df <- re_df[re_df$var1 == re_df$var2, ]
    names_vec <- paste(re_df$grp, re_df$var1, sep = ":")
    re_df <- setNames(re_df$vcov, names_vec)
    list(fixed_coefs=coefs,
         random_coefs=t(re_df))
  })

  ## To array, get names, build return list
  P <- nrow(results[[1]]$fixed_coefs)
  Pr <- ncol(results[[1]]$random_coefs)
  fixed_arr <- do.call(rbind, lapply(results, function(x) x$fixed_coefs))
  random_arr <- do.call(rbind, lapply(results, function(x) x$random_coefs))
  fnames <- row.names(results[[1]]$fixed_coefs)
  rnames <- colnames(results[[1]]$random_coefs)
  fixed_dn <- list(fnames, NULL, NULL)
  random_dn <- list(rnames, NULL, NULL)
  return(list(
    estimate = array(fixed_arr[,1], c(P, D, S), dimnames=fixed_dn),
    std.error = array(fixed_arr[,2], c(P, D, S), dimnames=fixed_dn),
    df = array(fixed_arr[,3], c(P, D, S), dimnames=fixed_dn),
    p.lower = array(fixed_arr[,4], c(P, D, S), dimnames=fixed_dn),
    p.upper = array(fixed_arr[,5], c(P, D, S), dimnames=fixed_dn),
    random.eff = array(random_arr, c(Pr, D, S), dimnames=random_dn)
  ))
}
