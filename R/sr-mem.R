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
##' @param contrasts optional validated contrast matrix with one row per
##'   contrast and columns matching the fixed-effect design.
##' @return A list of (P x D x S)-arrays with the fixed effect point estimates,
##'   the standard errors, degrees of freedom and the lower and upper p-values
##'   for each coefficient (P), of each model fit to each taxa (D) and each
##'   posterior sample (S) and a (Pr x D x S)-array with (Pr) random
##'   effects, plus flat K x D x S contrast arrays when contrasts are supplied.
##' @import nlme
##' @importFrom lmerTest lmer
##' @importFrom lme4 VarCorr
##' @importFrom parallel parLapply makeCluster stopCluster clusterExport clusterEvalQ
##' @author Kyle McGovern
sr.mem <- function(logW, formula, data, n.cores, method, mem.args,
                   contrasts = NULL) {
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

  X_fixed <- stats::model.matrix(reformulas::nobars(formula), data)
  contrasts <- .aldex_validate_contrasts(contrasts, X_fixed)
  if (!is.null(contrasts) && qr(X_fixed)$rank < ncol(X_fixed)) {
    stop("contrasts require a full-rank fixed-effect design")
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
    coefs <- coef(summary(fit))[, 1:3, drop = FALSE]
    coefs <- cbind(coefs, pt(coefs[,1]/coefs[,2], coefs[,3], lower.tail=TRUE))
    coefs <- cbind(coefs, pt(coefs[,1]/coefs[,2], coefs[,3], lower.tail=FALSE))
    colnames(coefs) <- c("estimate", "std.error", "df",
                         "p.lower", "p.upper")

    contrast_coefs <- NULL
    if (!is.null(contrasts)) {
      fixed_names <- rownames(coefs)
      if (!identical(fixed_names, colnames(contrasts))) {
        stop("fitted fixed-effect coefficients changed after contrast validation")
      }
      C <- contrasts
      K <- nrow(C)
      contrast_coefs <- matrix(NA_real_, nrow = K, ncol = 5L,
        dimnames = list(rownames(C), c("estimate", "std.error", "df", "p.lower", "p.upper")))

      if (method == "lme4") {
        for (k in seq_len(K)) {
          contest <- lmerTest::contest1D(
            fit, L = as.numeric(C[k, ]), ddf = "Satterthwaite"
          )
          estimate <- contest[1L, "Estimate"]
          std.error <- contest[1L, "Std. Error"]
          df <- contest[1L, "df"]
          if (!all(is.finite(c(estimate, std.error, df))) || std.error <= 0 || df <= 0) {
            stop("lmerTest produced invalid contrast inference")
          }
          t_stat <- estimate / std.error
          contrast_coefs[k, ] <- c(
            estimate, std.error, df,
            pt(t_stat, df = df, lower.tail = TRUE),
            pt(t_stat, df = df, lower.tail = FALSE)
          )
        }
      } else {
        V <- as.matrix(stats::vcov(fit))
        beta <- nlme::fixef(fit)
        for (k in seq_len(K)) {
          L <- matrix(C[k, ], nrow = 1L, dimnames = list(NULL, fixed_names))
          anova_L <- tryCatch(stats::anova(fit, L = L), error = function(e) e)
          if (inherits(anova_L, "error")) {
            stop(
              "nlme does not support this contrast because its nonzero fixed effects have different denominator degrees of freedom: ",
              conditionMessage(anova_L)
            )
          }
          df <- anova_L[1L, "denDF"]
          estimate <- as.numeric(L %*% beta)
          variance <- as.numeric(L %*% V %*% t(L))
          if (!all(is.finite(c(estimate, variance, df))) || variance <= 0 || df <= 0) {
            stop("nlme produced invalid contrast inference")
          }
          std.error <- sqrt(variance)
          t_stat <- estimate / std.error
          contrast_coefs[k, ] <- c(
            estimate, std.error, df,
            pt(t_stat, df = df, lower.tail = TRUE),
            pt(t_stat, df = df, lower.tail = FALSE)
          )
        }
      }
    }

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
         contrast_coefs=contrast_coefs,
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

  if (!is.null(contrasts)) {
    K <- nrow(results[[1]]$contrast_coefs)
    contrast_arr <- matrix(NA_real_, nrow = length(results) * K, ncol = 5L)
    colnames(contrast_arr) <- colnames(results[[1]]$contrast_coefs)
  }

  current_row_f <- 1
  current_row_r <- 1
  for (i in seq_along(results)) {
    fixed_arr[current_row_f:(current_row_f + P - 1), ] <-
      results[[i]]$fixed_coefs
    random_arr[current_row_r:(current_row_r + Pr - 1), ] <-
      as.matrix(results[[i]]$random_coefs)
    if (!is.null(contrasts)) {
      contrast_arr[((i - 1L) * K + 1L):(i * K), ] <- results[[i]]$contrast_coefs
    }
    current_row_f <- current_row_f + P
    current_row_r <- current_row_r + Pr
  }

  fnames <- row.names(results[[1]]$fixed_coefs)
  rnames <- colnames(results[[1]]$random_coefs)
  fixed_dn <- list(fnames, NULL, NULL)
  random_dn <- list(rnames, NULL, NULL)
  out <- list(
    estimate = array(fixed_arr[,1], c(P, D, S), dimnames=fixed_dn),
    std.error = array(fixed_arr[,2], c(P, D, S), dimnames=fixed_dn),
    df = array(fixed_arr[,3], c(P, D, S), dimnames=fixed_dn),
    p.lower = array(fixed_arr[,4], c(P, D, S), dimnames=fixed_dn),
    p.upper = array(fixed_arr[,5], c(P, D, S), dimnames=fixed_dn),
    random.eff = array(random_arr, c(Pr, D, S), dimnames=random_dn)
  )
  if (!is.null(contrasts)) {
    contrast_dn <- list(rownames(results[[1]]$contrast_coefs), NULL, NULL)
    out$contrast.estimate <- array(contrast_arr[, 1], c(K, D, S), dimnames = contrast_dn)
    out$contrast.std.error <- array(contrast_arr[, 2], c(K, D, S), dimnames = contrast_dn)
    out$contrast.df <- array(contrast_arr[, 3], c(K, D, S), dimnames = contrast_dn)
    out$contrast.p.lower <- array(contrast_arr[, 4], c(K, D, S), dimnames = contrast_dn)
    out$contrast.p.upper <- array(contrast_arr[, 5], c(K, D, S), dimnames = contrast_dn)
  }
  return(out)
}
