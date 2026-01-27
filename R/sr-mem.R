##' Implementation of SR-MEM: scale-reliant mixed effects models.
##'
##' @param logW a numeric array (N x D X S) where D is the number of taxa,
##'   N is the number of samples, and S is the number of posterior samples
##' @param formula an lme4 formula with fixed and random effects for lmer
##' @param data A data.frame with the random and fixed effects items in formula
##' @param n.cores The number of cores for parallelization. If n.cores=1, no
##'   parallelization is used
##' @param method The mem method to use: either nlme or lme4
##' @param mem.args Additional arguments including random or correlation
##' @return A list of (P x D x S)-arrays with the fixed effect point estimates,
##'   the standard errors, degrees of freedom and the lower and upper p-values
##'   for each coefficient (P), of each model fit to each taxa (D) and each
##'   posterior sample (S) and a (Pr x D x S)-array with (Pr) random
##'   effects.
##' @import nlme
##' @importFrom lmerTest lmer
##' @importFrom lme4 VarCorr
##' @importFrom parallel parLapply makeCluster stopCluster clusterExport clusterEvalQ
##' @author Kyle McGovern
sr.mem <- function(logW, formula, data, n.cores, method, mem.args) {
  N <- dim(logW)[1]
  D <- dim(logW)[2]
  S <- dim(logW)[3]
  logWm <- array(logW, c(N, D*S))

  ## Checks
  if(!(method%in%c("nlme", "lme4"))) {
    stop("method should be nlme or lme4")
  }

  if(method=="nlme") {
    if(!("random"%in%names(mem.args))) {
      stop("If method=\"nlme\", arg \"random\" must be provided")
    } else {
      random <- mem.args[["random"]]
    }
  }

  ## Get args
  if("correlation"%in%names(mem.args)) {
    correlation <- mem.args[["correlation"]]
  } else {
    correlation <- NULL
  }

  ## Only start cluster if n.cores > 1
  if(n.cores>1) {
    cl <- makeCluster(n.cores)
    on.exit(stopCluster(cl), add=TRUE)
    clusterEvalQ(cl, {
      library(lmerTest)
      library(lme4)
      library(nlme)
    })

    ## Trick to avoid code repeat
    lapply_func <- function(X, FUN) parLapply(cl, X, FUN)
  } else {
    lapply_func <- lapply
  }

  ## Process Cols (i.e., per-taxa)
  logWm_list <- split(logWm, col(logWm))
  results <- lapply_func(logWm_list, function(y) {
    data_temp <- data
    data_temp$y <- y
    if(method=="lme4") {
      fit <- suppressMessages(
        lmerTest::lmer(update(formula, y~.), data=data_temp)
      )
    } else if(method=="nlme") {
       fit <- suppressMessages(
        nlme::lme(update(formula, y~.), data=data_temp,
                  random=random, correlation=correlation)
      )
    }

    ## Fixed Effects for both lme4 & nlme
    coefs <- coef(summary(fit))[,1:3]
    coefs <- cbind(coefs, pt(coefs[,1]/coefs[,2], coefs[,3], lower.tail=TRUE))
    coefs <- cbind(coefs, 1-coefs[,4])
    colnames(coefs) <- c("estimate", "std.error", "df",
                         "p.lower", "p.upper")

    ## Random Effects for lme4/nlme
    if(method=="lme4") {
      re_df <- as.data.frame(VarCorr(fit))
      re_df_cnames <- apply(re_df[,c("grp", "var1", "var2")], 1,
                            function(x) paste(na.omit(x),
                                              collapse = "."))
      re_df <- data.frame(variance=re_df$vcov,
                          row.names=re_df_cnames)
    } else if (method=="nlme") {
      vc <- VarCorr(fit)
      re_df <- data.frame(
        row.names = rownames(vc),
        variance = as.numeric(vc[ , "Variance"]),
        stringsAsFactors = FALSE
      )
    }

    list(fixed_coefs=coefs,
         random_coefs=t(re_df))
  })

  ## To array, get names, build return list
  P <- nrow(results[[1]]$fixed_coefs)
  Pr <- ncol(results[[1]]$random_coefs)

  n_rows <- length(results)*P
  n_cols <- ncol(results[[1]]$fixed_coefs)
  fixed_arr <- matrix(NA, nrow = n_rows, ncol = n_cols)
  colnames(fixed_arr) <- colnames(results[[1]]$fixed_coefs)

  n_rows <- length(results)*Pr
  n_cols <- ncol(results[[1]]$random_coefs)
  random_arr <- matrix(NA, nrow = n_rows, ncol = n_cols)
  colnames(random_arr) <- colnames(results[[1]]$random_coefs)

  current_row_f <- 1
  current_row_r <- 1
  for (i in seq_along(results)) {
    fixed_arr[current_row_f:(current_row_f + P - 1), ] <-
      results[[i]]$fixed_coefs
    random_arr[current_row_r:(current_row_r + Pr - 1), ] <-
      as.matrix(results[[i]]$random_coefs)
    current_row_f <- current_row_f + P
    current_row_r <- current_row_r + Pr
  }

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
