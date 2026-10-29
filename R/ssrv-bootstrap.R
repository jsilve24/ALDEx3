.ssrv_subject_indices <- function(groups) {
  unlist(lapply(groups, function(g)
    g[sample.int(length(g), length(g), replace=TRUE)]), use.names=FALSE)
}

.aldex_guard_ssrv <- function(object, helper) {
  if (is.list(object) && identical(object$test, "ssrv.bootstrap"))
    stop(helper, "() is not yet available for test = 'ssrv.bootstrap'.", call.=FALSE)
}

.ssrv_posterior <- function(draws) {
  dims <- dim(draws)[1:2]
  dn <- dimnames(draws)[1:2]
  mat <- function(x) matrix(x, dims[1], dims[2], dimnames=dn)
  list(mean=mat(apply(draws, 1:2, function(x) if(all(is.finite(x))) mean(x) else NA_real_)),
       lower=mat(apply(draws, 1:2, function(x) if(all(is.finite(x))) quantile(x,0.025,names=FALSE) else NA_real_)),
       upper=mat(apply(draws, 1:2, function(x) if(all(is.finite(x))) quantile(x,0.975,names=FALSE) else NA_real_)),
       level=0.95, interval="equal-tailed")
}

.aldex_ssrv_summary <- function(object, ignore.intercept, component) {
  extract <- function(values, coefficient) {
    p <- values$posterior
    if (is.null(p) && !is.null(values$estimate)) p <- .ssrv_posterior(values$estimate)
    if (is.null(p)) stop("posterior summaries or retained estimate draws are required")
    out <- array_to_df(p$mean)
    names(out) <- c("parameter", "entity", "estimate")
    out$credible.lower <- as.vector(p$lower)
    out$credible.upper <- as.vector(p$upper)
    if (length(object$diagnostics$unavailable.coefficients))
      out$status <- ifelse(is.finite(out$estimate),"available","unavailable")
    if (coefficient && ignore.intercept) out <- out[out$parameter != "(Intercept)", , drop=FALSE]
    rownames(out) <- NULL
    attr(out, "level") <- p$level
    attr(out, "interval") <- p$interval
    out
  }
  if (component == "coefficients") return(extract(object, TRUE))
  if (component == "contrasts") {
    if (is.null(object$contrasts)) stop("object does not contain requested contrasts")
    return(extract(object$contrasts, FALSE))
  }
  out <- extract(object, TRUE)
  if (!is.null(object$contrasts)) {
    con <- extract(object$contrasts, FALSE)
    out$type <- rep("coefficient", nrow(out)); con$type <- rep("contrast", nrow(con))
    out <- rbind(out, con)
  }
  attr(out, "level") <- 0.95
  attr(out, "interval") <- "equal-tailed"
  out
}

.aldex_ssrv_bootstrap <- function(Y, X, data, nsample, scale, streamsize,
                                 return.pars, contrasts, strata, max.attempts, args) {
  integer.scalar <- function(x) is.numeric(x) && length(x)==1L && is.finite(x) && x==floor(x)
  if (!integer.scalar(nsample) || nsample < 2L) stop("ssrv.bootstrap requires nsample >= 2")
  if (!integer.scalar(max.attempts) || max.attempts < 1L)
    stop("bootstrap.max.attempts must be a positive integer")
  if (!is.matrix(Y) || !is.numeric(Y) || any(!is.finite(Y)) || any(Y < 0) || any(dim(Y) < 1))
    stop("Y must be a nonempty finite nonnegative count matrix")
  if (!is.numeric(streamsize) || length(streamsize)!=1L || is.na(streamsize) || streamsize <= 0)
    stop("streamsize must be positive")
  if (!is.character(return.pars) || !length(return.pars) ||
      any(!return.pars %in% c("posterior", "estimate", "X")))
    stop("ssrv.bootstrap return.pars supports only 'posterior', 'estimate', and 'X'")
  N <- ncol(Y); D <- nrow(Y)
  formula <- if (inherits(X, "formula")) X else NULL
  if (!is.null(formula)) {
    if (is.null(data)) stop("data should not be null if X is a formula")
    X <- t(model.matrix(formula, data))
  }
  if (!is.matrix(X) || !is.numeric(X) || ncol(X)!=N || !all(is.finite(X)) || nrow(X)<1)
    stop("X must be a finite numeric P x N model matrix")
  P <- nrow(X)
  if (is.null(rownames(X))) rownames(X) <- paste0("parameter_", seq_len(P))
  .ssrv_names(rownames(X))
  if (qr(t(X))$rank < P) stop("rank-deficient designs are not supported by method 'lm'")
  contrasts <- .aldex_validate_contrasts(contrasts, t(X))
  stream <- N*D*nsample*8/100000 > streamsize
  if (stream && floor(streamsize*100000/(N*D*8)) < 1) stop("streamsize too small")
  policy <- "paired rows"
  if (is.null(strata) && P==2L && any(apply(X, 1, function(x) all(x==1)))) {
    k <- which(!apply(X, 1, function(x) all(x==1)))
    if (length(k)==1L && length(unique(X[k, ]))==2L) {
      strata <- X[k, ]; policy <- "within two groups"
    }
  } else if (!is.null(strata)) policy <- "explicit strata"
  if (!is.null(strata) && (length(strata)!=N || !is.atomic(strata) || !is.null(dim(strata)) || anyNA(strata)))
    stop("bootstrap.strata must be a nonmissing length-N atomic vector")
  groups <- if (is.null(strata)) list(seq_len(N)) else split(seq_len(N), strata, drop=TRUE)
  if (!length(groups)) stop("bootstrap.strata has no groups")
  if (inherits(scale, "aldex_sparse_scale") || inherits(scale, "aldex_gaussian_scale"))
    scale <- compose.sm(scale)
  composable <- inherits(scale, "aldex_coefficient_scale")
  custom <- is.function(scale) && "A" %in% names(formals(scale))
  if (composable) {
    assigned <- unlist(lapply(scale, `[[`, "coefficients"), use.names=FALSE)
    if (anyDuplicated(assigned)) stop("duplicate or overlapping coefficient assignments")
    if (!setequal(assigned, rownames(X)))
      stop("scale assignments must match exact design coefficients; unknown: ",
           paste(setdiff(assigned, rownames(X)), collapse=", "), "; missing: ",
           paste(setdiff(rownames(X), assigned), collapse=", "))
  } else if (!(custom || is.function(scale) || is.matrix(scale))) {
    stop("supply a supported scale function, scale matrix, or compose.sm specification")
  }
  if (is.matrix(scale) && (!is.numeric(scale) || !identical(dim(scale), c(N, as.integer(nsample))) || any(!is.finite(scale))))
    stop("scale matrix must be finite N x nsample")
  entities <- rownames(Y)
  if (is.null(entities)) entities <- paste0("entity_", seq_len(D))
  A <- array(NA_real_, c(P,D,nsample), list(rownames(X), entities, NULL))
  shifts <- matrix(0, P, nsample, dimnames=list(rownames(X),NULL))
  laplace <- boot <- matrix(NA_real_, P, nsample, dimnames=list(rownames(X),NULL))
  invalid <- matrix(FALSE, P, nsample)
  laplace.max.variance <- setNames(rep(NA_real_,P),rownames(X))
  unbounded <- setNames(integer(P),rownames(X))
  rejected <- integer(nsample)
  sparse <- if (composable) unlist(lapply(scale, function(b)
    if (inherits(b,"aldex_sparse_scale")) match(b$coefficients, rownames(X)) else integer()), use.names=FALSE) else integer()
  for (s in seq_len(nsample)) {
    lc <- array(NA_real_, c(D,N,1L))
    for (i in seq_len(N)) lc[,i,1L] <- rLogDirichlet(1L, Y[,i]+0.5)
    for (attempt in seq_len(max.attempts)) {
      idx <- .ssrv_subject_indices(groups)
      xb <- t(X[,idx,drop=FALSE])
      if (qr(xb)$rank == P) break
      rejected[s] <- rejected[s]+1L
    }
    if (qr(xb)$rank < P)
      stop("rank-deficient bootstrap design at draw ", s, " after ", max.attempts,
           " attempts; rejected resamples so far: ", sum(rejected),
           ". Revise design or specify appropriate bootstrap.strata.")
    a <- fflm(aperm(lc[,idx,,drop=FALSE], c(2,1,3)), xb,
              .coefficients.only=TRUE)$estimate[,,1L,drop=FALSE]
    a <- matrix(a, P,D, dimnames=list(rownames(X),entities))
    A[,,s] <- a
    if (composable) {
      # One common feature bootstrap for all sparse coefficients in this draw.
      fi <- if (length(sparse)) sample.int(D,D,replace=TRUE) else integer()
      for (block in scale) {
        ks <- match(block$coefficients, rownames(X))
        if (inherits(block, "aldex_gaussian_scale")) {
          shifts[ks,s] <- if (!is.null(block$sd)) rnorm(length(ks), block$mean, block$sd) else
            as.vector(.draw_scale_mvn(1L,block$mean,block$covariance))
        } else for (k in ks) {
          if (block$uncertainty != "laplace") {
            kde <- .ssrv_kde(a[k,fi],block$bandwidth)
            boot[k,s] <- -kde$mode
          }
          if (block$uncertainty != "bootstrap") {
            kde <- .ssrv_kde(a[k,],block$bandwidth)
            v <- kde$variance
            if (!is.na(v)) laplace.max.variance[k] <-
              if(is.na(laplace.max.variance[k])) v else max(laplace.max.variance[k],v)
            unbounded[k] <- unbounded[k]+as.integer(is.infinite(v) && v>0)
            invalid[k,s] <- !is.finite(v) || v<=0
            laplace[k,s] <- if(invalid[k,s]) NA_real_ else rnorm(1L,-kde$mode,sqrt(v))
          }
        }
      }
    } else if (custom) {
      context <- c(list(A=a,X=X,logComp=lc,indices=idx,draw=s),args)
      context <- context[!duplicated(names(context))]
      if (!"..." %in% names(formals(scale))) context <- context[intersect(names(context),names(formals(scale)))]
      b <- do.call(scale,context)
      if (!is.numeric(b) || !is.null(dim(b)) || length(b)!=P || any(!is.finite(b)) ||
          is.null(names(b)) || anyDuplicated(names(b)) || !setequal(names(b),rownames(X)))
        stop("custom coefficient scale must return a finite named P-vector matching X")
      shifts[,s] <- b[rownames(X)]
    } else {
      if (is.matrix(scale)) ls <- scale[,s,drop=FALSE] else {
        context <- c(list(X=X,logComp=lc,Y=Y),args)
        context <- context[!duplicated(names(context))]
        # Match existing callback argument filtering, including deprecated arguments.
        context <- context[intersect(names(context), names(formals(scale)))]
        ls <- do.call(scale,context)
      }
      if (!is.matrix(ls) || !is.numeric(ls) || !identical(dim(ls),c(N,1L)) || any(!is.finite(ls)))
        stop("sample-scale callback must return a finite N x 1 matrix for each draw")
      shifts[,s] <- fflm(array(ls[idx,1L],c(N,1L,1L)),xb,
                          .coefficients.only=TRUE)$estimate[,1L,1L]
    }
  }
  branches <- setNames(character(length(sparse)),rownames(X)[sparse])
  if (composable) for (block in scale) if (inherits(block,"aldex_sparse_scale")) {
    for (k in match(block$coefficients,rownames(X))) {
      selected <- .ssrv_select(laplace[k,],boot[k,],block$uncertainty)
      shifts[k,] <- if(is.na(selected)) NA_real_ else
        if(selected=="bootstrap") boot[k,] else laplace[k,]
      branches[rownames(X)[k]] <- selected
    }
  }
  unavailable <- rownames(X)[apply(shifts,1,function(x) any(!is.finite(x)))]
  if(length(unavailable)) {
    reason <- if(any(unbounded>0)) "Unbounded Laplace variance" else "Undefined Laplace approximation"
    warning(reason," for: ",paste(unavailable,collapse=", "),
      "; affected results are unavailable, with no automatic bootstrap substitution. ",
      "Consider explicitly requesting uncertainty='bootstrap', or reconsider whether ",
      "the sparse scale model fits these data.",call.=FALSE)
    shifts[unavailable,] <- NA_real_
  }
  for (s in seq_len(nsample)) A[,,s] <- sweep(matrix(A[,,s],P,D),1L,shifts[,s],`+`)
  res <- list(X=X, method="lm", test="ssrv.bootstrap", streaming=stream,
              nsample=nsample, scale.basis="log2 composition",
              diagnostics=list(resampling=policy, strata=strata, max.attempts=max.attempts,
                               rejected=sum(rejected), rejections.by.draw=rejected,
                               scale.branch=branches, unavailable.coefficients=unavailable,
                               available.draws=setNames(ifelse(rownames(X)%in%unavailable,0L,nsample),rownames(X)),
                               laplace.unbounded.draws=unbounded,
                               laplace.max.conditional.variance=laplace.max.variance,
                               laplace.unavailable.draws=setNames(rowSums(invalid),rownames(X)),
                               laplace.variance=if (length(sparse)) apply(laplace[sparse,,drop=FALSE],1,var) else numeric(),
                               bootstrap.variance=if (length(sparse)) apply(boot[sparse,,drop=FALSE],1,var) else numeric()))
  if ("posterior" %in% return.pars) res$posterior <- .ssrv_posterior(A)
  if ("estimate" %in% return.pars) res$estimate <- A
  if (!is.null(contrasts)) {
    ca <- array(NA_real_,c(nrow(contrasts),D,nsample),list(rownames(contrasts),entities,NULL))
    for (r in seq_len(nrow(contrasts))) {
      active <- which(contrasts[r,]!=0)
      for (s in seq_len(nsample))
        ca[r,,s] <- contrasts[r,active,drop=FALSE] %*%
          matrix(A[active,,s,drop=FALSE],length(active),D)
    }
    res$diagnostics$unavailable.contrasts <- rownames(contrasts)[
      apply(contrasts[,unavailable,drop=FALSE]!=0,1,any)]
    res$contrasts <- list(matrix=contrasts)
    if ("posterior" %in% return.pars) res$contrasts$posterior <- .ssrv_posterior(ca)
    if ("estimate" %in% return.pars) res$contrasts$estimate <- ca
  }
  if (!is.null(formula)) res$formula <- formula
  if (!is.null(data)) res$data <- data
  class(res) <- "aldex"
  res
}
