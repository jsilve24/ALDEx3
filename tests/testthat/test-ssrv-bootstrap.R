ssrv_fixture <- function() {
  set.seed(101)
  data <- data.frame(treatment=rep(0:1,each=6), age=seq_len(12)/12)
  Y <- matrix(rpois(15*12,100),15,12,
              dimnames=list(paste0("taxon",1:15),paste0("sample",1:12)))
  list(Y=Y,data=data)
}
ssrv_scale <- function(mode="hybrid") compose.sm(
  gaussian.sm("(Intercept)",sd=0), sparse.sm("treatment",mode,bandwidth=0.2))

test_that("coefficients-only fflm preserves estimator and skips inference", {
  set.seed(2)
  x <- cbind(`(Intercept)`=1, x=rnorm(12))
  y <- array(rnorm(12*4*3),c(12,4,3))
  C <- matrix(c(0,2),1,dimnames=list("twice",colnames(x)))
  a <- fflm(y,x,contrasts=C,.coefficients.only=TRUE)
  b <- fflm(y,x,contrasts=C)
  expect_named(a,c("estimate","contrast.estimate"))
  expect_equal(a$estimate,b$estimate)
  expect_equal(unname(a$contrast.estimate),unname(b$contrast.estimate))
  expect_equal(a$estimate[,,1],unname(coef(lm(y[,,1]~x[,2]))))
})

test_that("stratified draws agree with independent resampled group means", {
  f <- ssrv_fixture(); seen <- list()
  capture <- function(A,X,logComp,indices,draw) {
    seen[[draw]] <<- list(A=A,lc=logComp,idx=indices)
    setNames(rep(0,nrow(X)),rownames(X))
  }
  set.seed(18)
  fit <- aldex(f$Y,~treatment,f$data,nsample=12,test="ssrv.bootstrap",scale=capture,
               bootstrap.strata=f$data$treatment)
  x <- f$data$treatment
  expected <- vapply(seen,function(z) {
    lc <- z$lc[,,1]; i <- z$idx
    rowMeans(lc[,i[x[i]==1],drop=FALSE])-rowMeans(lc[,i[x[i]==0],drop=FALSE])
  },numeric(15))
  expect_lt(max(abs(fit$estimate[2,,]-expected)),1e-12)
  expect_true(all(vapply(seen,function(z) identical(sort(x[z$idx]),sort(x)),logical(1))))
  expect_equal(fit$diagnostics$resampling,"explicit strata")
  seen <- list()
  set.seed(18)
  fit <- aldex(f$Y,~treatment,f$data,nsample=12,test="ssrv.bootstrap",scale=capture)
  expect_equal(fit$diagnostics$resampling,"paired rows")
  expect_true(any(vapply(seen,function(z) sum(x[z$idx]==1)!=6,logical(1))))
  expected <- simplify2array(lapply(seen,function(z)
    coef(lm.fit(t(fit$X[,z$idx]),t(z$lc[,z$idx,1])))))
  expect_lt(max(abs(simplify2array(lapply(seen,`[[`,"A"))-expected)),1e-10)
})

test_that("Laplace uses original-subject OLS for its mode and final effects", {
  f <- ssrv_fixture(); comps <- list(); profiles <- list()
  draw <- rLogDirichlet
  local_mocked_bindings(rLogDirichlet=function(...) {
    z <- draw(...); comps[[length(comps)+1L]] <<- z; z
  }, .ssrv_kde=function(a,bandwidth=NULL) {
    profiles[[length(profiles)+1L]] <<- a
    list(mode=mean(a),variance=0.04)
  }, .package="ALDEx3")
  set.seed(86)
  fit <- aldex(f$Y,~treatment,f$data,nsample=5,test="ssrv.bootstrap",
               scale=ssrv_scale("laplace"))
  expected <- vapply(1:5,function(s) {
    lc <- do.call(cbind,comps[((s-1)*12+1):(s*12)])
    coef(lm.fit(t(fit$X),t(lc)))[2,]
  },numeric(15))
  expect_lt(max(abs(do.call(cbind,profiles)-expected)),1e-10)
  shifts <- fit$estimate[2,,]-expected
  expect_lt(max(abs(sweep(shifts,2,shifts[1,]))),1e-10)
})

test_that("Laplace-only sparse models require no resampling", {
  f <- ssrv_fixture(); comps <- list(); draw <- rLogDirichlet
  local_mocked_bindings(rLogDirichlet=function(...) {
    z <- draw(...); comps[[length(comps)+1L]] <<- z; z
  }, .ssrv_subject_indices=function(...) stop("unexpected subject bootstrap"),
  .package="ALDEx3")
  fit <- aldex(f$Y,~treatment,f$data,nsample=3,test="ssrv.bootstrap",
               scale=ssrv_scale("laplace"),bootstrap.strata=f$data$treatment)
  expect_identical(fit$diagnostics$resampling,"original subjects")
  expect_identical(fit$diagnostics$laplace.resampling,"original subjects")
  expect_identical(fit$diagnostics$rejected,0L)
  expected <- vapply(1:3,function(s) {
    lc <- do.call(cbind,comps[((s-1)*12+1):(s*12)])
    coef(lm.fit(t(fit$X),t(lc)))[1,]
  },numeric(15))
  expect_lt(max(abs(fit$estimate[1,,]-expected)),1e-10)
  expect_error(aldex(f$Y,~treatment,f$data,nsample=3,test="ssrv.bootstrap",
                     scale=ssrv_scale()),"unexpected subject bootstrap")
})

test_that("hybrid selects and returns one complete candidate effect sequence", {
  f <- ssrv_fixture(); select <- .ssrv_select; candidates <- NULL
  draw <- rLogDirichlet; resample <- .ssrv_subject_indices
  comps <- list(); indices <- list()
  local_mocked_bindings(.ssrv_select=function(laplace,bootstrap,mode) {
    candidates <<- list(laplace=laplace,bootstrap=bootstrap)
    select(laplace,bootstrap,mode)
  },rLogDirichlet=function(...) {
    z <- draw(...); comps[[length(comps)+1L]] <<- z; z
  },.ssrv_subject_indices=function(groups) {
    z <- resample(groups); indices[[length(comps) %/% ncol(f$Y)]] <<- z; z
  },.package="ALDEx3")
  set.seed(52)
  fit <- aldex(f$Y,~treatment,f$data,nsample=12,test="ssrv.bootstrap",scale=ssrv_scale())
  widths <- vapply(candidates,function(z)
    median(apply(z,1,function(x) diff(quantile(x,c(.025,.975))))),numeric(1))
  branch <- if (widths["laplace"] >= widths["bootstrap"]) "laplace" else "bootstrap"
  expect_identical(unname(fit$diagnostics$scale.branch["treatment"]),branch)
  expect_equal(unname(fit$estimate[2,,]),unname(candidates[[branch]]),tolerance=1e-10)
  expect_equal(unname(fit$diagnostics$laplace.interval.width["treatment"]),unname(widths["laplace"]))
  expect_equal(unname(fit$diagnostics$bootstrap.interval.width["treatment"]),unname(widths["bootstrap"]))
  for (candidate in c("bootstrap","laplace")) {
    profiles <- vapply(1:12,function(s) {
      lc <- do.call(cbind,comps[((s-1)*12+1):(s*12)])
      i <- if (candidate=="bootstrap") indices[[s]] else seq_len(12)
      coef(lm.fit(t(fit$X[,i]),t(lc[,i])))[2,]
    },numeric(15))
    shifts <- candidates[[candidate]]-profiles
    expect_lt(max(abs(sweep(shifts,2,shifts[1,]))),1e-10)
  }
})

test_that("scales are shared scalar shifts and all coefficients retain subject uncertainty", {
  f <- ssrv_fixture(); f$Y <- round(f$Y*1e7)
  captured <- list()
  scale <- function(A,draw) {
    captured[[draw]] <<- A
    c(`(Intercept)`=2,treatment=-1,age=0.5)
  }
  fit <- aldex(f$Y,~treatment+age,f$data,nsample=25,scale=scale,test="ssrv.bootstrap")
  shifts <- fit$estimate-simplify2array(captured)
  expect_lt(max(abs(shifts-array(c(2,-1,0.5),dim(shifts)))),1e-10)
  expect_true(all(apply(fit$estimate,c(1,2),sd)>1e-4))
})

test_that("KDE mode, inverse log curvature, global peaks and translation behave", {
  z <- .ssrv_kde(rep(7,10),0.4)
  expect_equal(z$mode,7); expect_equal(z$variance,0.16)
  for (a in list(rep(0,10),c(-2,-1.9,-2.1,2,2.05,2.1,2.15,2.2))) {
    u <- .ssrv_kde(a,0.4); v <- .ssrv_kde(a+13,0.4)
    expect_equal(v$mode,u$mode+13,tolerance=1e-7)
    expect_equal(v$variance,u$variance,tolerance=1e-7)
  }
  a <- c(rep(-5,10),rep(3,30))
  z <- .ssrv_kde(a,0.2)
  expect_equal(z$mode,3,tolerance=1e-6)
  expect_equal(z$variance,0.04,tolerance=1e-6)
  expect_equal(.ssrv_kde(c(-2,2),0.1)$mode,-2,tolerance=1e-6)
  expect_error(sparse.sm("x",bandwidth=0),"positive")
  # Equal Gaussian kernels at +/-h have a flat maximum with zero log curvature.
  expect_identical(.ssrv_kde(c(-1,1),1)$variance,Inf)
  expect_true(is.finite(.ssrv_kde(c(-1,1),1.000000001)$variance))
  small <- c(-1,0,1)*1e-10
  expect_equal(.ssrv_kde(small)$bandwidth,bw.nrd0(small))
  expect_error(.ssrv_kde(rep(0,10)),"explicit positive bandwidth")
})

test_that("hybrid uses final-effect widths and never substitutes unavailable candidates", {
  l <- rbind(c(-1,1,-1,1),c(-1,1,-1,1),c(-1,1,-1,1))
  b <- rbind(c(-2,-2,2,2),c(-2,-2,2,2),c(-2,-2,2,2))
  expect_equal(.ssrv_select(l,b,"hybrid"),"bootstrap")
  expect_equal(.ssrv_select(l-b,b-b,"hybrid"),"laplace")
  expect_equal(.ssrv_select(l,l,"hybrid"),"laplace")
  # A rare extreme makes Laplace variance larger but its central interval narrower.
  rare <- matrix(c(rep(0,99),10),nrow=1)
  regular <- matrix(rep(c(-.5,.5),50),nrow=1)
  expect_gt(var(c(rare)),var(c(regular)))
  expect_equal(.ssrv_select(rare,regular,"hybrid"),"bootstrap")
  l[1,1] <- NA
  expect_true(is.na(.ssrv_select(l,b,"hybrid")))
  f <- ssrv_fixture()
  local_mocked_bindings(.ssrv_kde=function(a,bandwidth=NULL)
    list(mode=mean(a),variance=Inf),.package="ALDEx3")
  expect_warning(lap <- aldex(f$Y,~treatment,f$data,nsample=3,
    scale=ssrv_scale("laplace"),test="ssrv.bootstrap"),"Unbounded Laplace variance")
  C <- diag(2); dimnames(C) <- list(c("baseline","treatment"),c("(Intercept)","treatment"))
  expect_warning(b <- aldex(f$Y,~treatment,f$data,nsample=3,scale=ssrv_scale(),
    contrasts=C,test="ssrv.bootstrap"),"Unbounded Laplace variance")
  expect_true(is.na(b$diagnostics$scale.branch[["treatment"]]))
  expect_identical(b$diagnostics$laplace.max.conditional.variance[["treatment"]],Inf)
  expect_true(all(is.na(b$posterior$mean["treatment",])))
  expect_true(all(is.finite(b$contrasts$posterior$mean["baseline",])))
  expect_true(all(is.na(b$contrasts$posterior$mean["treatment",])))
  expect_true(all(summary(b,component="coefficients")$status=="unavailable"))
  expect_null(b$diagnostics$laplace.fallback)
  expect_true(all(is.na(lap$posterior$mean["treatment",])))
})

test_that("sparse corrections are translation equivariant under paired RNG", {
  f <- ssrv_fixture()
  original <- fflm
  run <- function() aldex(f$Y,~treatment,f$data,nsample=12,
                         scale=ssrv_scale(),test="ssrv.bootstrap")
  set.seed(33); a <- run()
  local_mocked_bindings(fflm=function(...) {
    z <- original(...); z$estimate[2,,] <- z$estimate[2,,]+4; z
  },.package="ALDEx3")
  set.seed(33); b <- run()
  expect_equal(a$estimate,b$estimate,tolerance=1e-7)
})

test_that("posterior summaries, contrasts, dimensions, retention and streaming agree", {
  f <- ssrv_fixture()
  x <- model.matrix(~treatment,f$data)
  C <- matrix(c(-1,2),1,dimnames=list("comparison",colnames(x)))
  run <- function(stream=Inf,pars=c("posterior","estimate")) {
    set.seed(16)
    aldex(f$Y,~treatment,f$data,nsample=20,test="ssrv.bootstrap",scale=ssrv_scale(),
          contrasts=C,streamsize=stream,return.pars=pars)
  }
  a <- run(); b <- run(0.03); compact <- run(0.03,"posterior")
  expect_equal(a$estimate,b$estimate)
  expect_equal(a$posterior,compact$posterior)
  expect_equal(a$diagnostics,b$diagnostics)
  expect_null(compact$estimate)
  expect_equal(a$posterior$mean,apply(a$estimate,1:2,mean))
  expect_equal(a$posterior$lower,apply(a$estimate,1:2,quantile,0.025))
  expect_equal(a$posterior$upper,apply(a$estimate,1:2,quantile,0.975))
  expected <- simplify2array(lapply(1:20,function(s) C%*%a$estimate[,,s]))
  expect_lt(max(abs(a$contrasts$estimate-expected)),1e-10)
  expect_equal(summary(a),summary(compact))
  expect_equal(summary(a),summary(run(Inf,"estimate")))
  expect_equal(a$p.val, compact$p.val)
  expect_equal(a$p.val.adj, compact$p.val.adj)
  score <- function(z) min(1, 2*min(mean(z <= 0),mean(z >= 0)))
  expect_equal(a$p.val, apply(a$estimate,1:2,score))
  expect_equal(a$contrasts$p.val, apply(a$contrasts$estimate,1:2,score))
  expect_equal(a$contrasts$p.val.adj[1,], p.adjust(a$contrasts$p.val[1,],"BH"))
  old.compact <- compact; old.compact$p.val <- old.compact$p.val.adj <- NULL
  expect_warning(old.summary <- summary(old.compact, component="coefficients"), "probabilities are unavailable")
  expect_true(all(is.na(old.summary$p.val)))
  expect_equal(attr(summary(a),"interval"),"equal-tailed")
  expect_equal(names(summary(a,component="coefficients")),
               c("parameter","entity","estimate","p2.5","p97.5","p.val","p.val.adj"))
  expect_true(all(c("coefficient","contrast") %in% summary(a)$type))
  expect_equal(summary(a,component="contrasts")$p2.5,
               as.vector(a$contrasts$posterior$lower))
  single <- aldex(f$Y[1,,drop=FALSE],t(x),nsample=5,test="ssrv.bootstrap",
                  scale=tss.sm,gamma=0,contrasts=C)
  expect_equal(dim(single$estimate),c(2,1,5))
  expect_equal(dim(single$contrasts$posterior$mean),c(1,1))
})

test_that("Gaussian blocks, legacy scales and joint custom scales keep their meaning", {
  f <- ssrv_fixture()
  sc <- compose.sm(gaussian.sm(c("(Intercept)","age"),mean=c(2,0.5),
                               covariance=matrix(c(1,1,1,1),2)),sparse.sm("treatment","bootstrap"))
  base <- fflm; captured <- list()
  local_mocked_bindings(fflm=function(...) {
    z <- base(...); captured[[length(captured)+1L]] <<- z$estimate; z
  },.package="ALDEx3")
  fit <- aldex(f$Y,~treatment+age,f$data,nsample=15,test="ssrv.bootstrap",scale=sc)
  shifts <- fit$estimate-simplify2array(lapply(captured,function(z) z[,,1]))
  expect_lt(max(abs((shifts[1,,]-2)-(shifts[3,,]-0.5))),1e-7)
  expect_true(sd(fit$estimate[3,1,])>0)
  x <- t(model.matrix(~treatment,f$data)); coeff <- c(1.5,-0.75)
  set.seed(12); prior <- aldex(f$Y,x,nsample=10,test="ssrv.bootstrap",scale=coefficient.sm,
                               c.mu=coeff,c.sd=c(0,0))
  set.seed(12); adapter <- aldex(f$Y,x,nsample=10,test="ssrv.bootstrap",
                                scale=matrix(t(x)%*%coeff,12,10))
  # Zero-SD Gaussian draws do not consume RNG, so both paths align exactly.
  expect_equal(prior$estimate,adapter$estimate)
  set.seed(12); samplefit <- aldex(f$Y,x,nsample=10,test="ssrv.bootstrap",scale=sample.sm,
                                  s.mu=as.vector(t(x)%*%coeff),s.sd=rep(0,12))
  expect_equal(prior$estimate,samplefit$estimate)
})

test_that("design identities and finite rank redraw policy are enforced", {
  f <- ssrv_fixture()
  f$data$treatment <- factor(f$data$treatment,levels=c(1,0))
  fit <- aldex(f$Y,~treatment,f$data,nsample=10,test="ssrv.bootstrap",scale=tss.sm,gamma=0)
  expect_equal(rownames(fit$X),c("(Intercept)","treatment0"))
  expect_equal(dimnames(fit$estimate)[[2]],rownames(f$Y))
  expect_error(aldex(f$Y,rbind(a=rep(1,12),b=rep(1,12)),nsample=3,
                     test="ssrv.bootstrap",scale=tss.sm),"rank-deficient")
  local_mocked_bindings(.ssrv_subject_indices=function(groups) rep(1L,12),.package="ALDEx3")
  # Explicit single stratum allows the controlled resample to be singular.
  expect_error(aldex(f$Y,~treatment,f$data,nsample=3,test="ssrv.bootstrap",scale=tss.sm,
                     bootstrap.strata=rep(1,12),bootstrap.max.attempts=2),"after 2 attempts")
})

test_that("scope checks and downstream guards fail early without breaking legacy objects", {
  f <- ssrv_fixture()
  for (m in c("blmm","lme4","nlme")) expect_error(
    aldex(f$Y,~treatment,f$data,method=m,test="ssrv.bootstrap"),"only method")
  expect_error(aldex(f$Y,~treatment+(1|age),f$data,test="ssrv.bootstrap"),"random-effects")
  expect_error(aldex(f$Y,~treatment,f$data,test="ssrv.bootstrap",correlation=1),"unsupported")
  expect_error(compose.sm(gaussian.sm("x",sd=1),sparse.sm("x")),"overlapping")
  expect_error(aldex(f$Y,~treatment,f$data,nsample=3,test="ssrv.bootstrap",scale=sparse.sm("treatment")),"missing")
  expect_error(aldex(f$Y,~treatment,f$data,nsample=1,test="ssrv.bootstrap",scale=tss.sm),"nsample")
  fit <- aldex(f$Y,~treatment,f$data,nsample=5,test="ssrv.bootstrap",scale=tss.sm)
  expect_false(is.null(fit$p.val)); expect_null(fit$std.error); expect_null(fit$df)
  for (fun in list(aldex.plot,aldex.effect,aldex.p.adjust,cohensd))
    expect_error(fun(fit),"not yet available")
  expect_error(aldex(f$Y,~treatment,f$data,nsample=3,test="ssrv.bootstrap",scale=tss.sm,
                     onesided="upper"),"one-sided")
  old <- aldex(f$Y,~treatment,f$data,nsample=5,scale=tss.sm)
  old$test <- NULL
  expect_s3_class(summary(old),"data.frame")
})

test_that("arbitrary sample scales use matched rows and exact projection", {
  f <- ssrv_fixture(); indices <- list(); comps <- list()
  draw <- .ssrv_subject_indices
  local_mocked_bindings(.ssrv_subject_indices=function(groups) {
    idx <- draw(groups); indices[[length(indices)+1L]] <<- idx; idx
  },.package="ALDEx3")
  values <- seq_len(12)^2/100
  callback <- function(X,logComp) {
    comps[[length(comps)+1L]] <<- logComp[,,1]
    matrix(values,12,1)
  }
  fit <- aldex(f$Y,~treatment,f$data,nsample=8,test="ssrv.bootstrap",scale=callback)
  expected <- simplify2array(lapply(1:8,function(s) {
    idx <- indices[[s]]
    coef(lm.fit(t(fit$X[,idx]),t(comps[[s]][,idx])+values[idx]))
  }))
  expect_lt(max(abs(fit$estimate-expected)),1e-10)
})

test_that("shared feature indices preserve relationships between sparse coefficients", {
  f <- ssrv_fixture(); original <- fflm
  local_mocked_bindings(fflm=function(...) {
    z <- original(...)
    z$estimate[3,,] <- z$estimate[2,,]+2
    z
  },.package="ALDEx3")
  sc <- compose.sm(gaussian.sm("(Intercept)",sd=0),
                   sparse.sm(c("treatment","age"),"bootstrap",bandwidth=0.2))
  fit <- aldex(f$Y,~treatment+age,f$data,nsample=8,test="ssrv.bootstrap",scale=sc)
  expect_equal(fit$estimate[2,,],fit$estimate[3,,],tolerance=1e-7)
})

test_that("matrix scale draws retain their global draw alignment across streaming", {
  f <- ssrv_fixture(); x <- t(model.matrix(~treatment,f$data))
  b <- rbind(seq_len(8)/10,seq_len(8)/20)
  run <- function(sc,stream=Inf) {
    set.seed(28)
    aldex(f$Y,x,nsample=8,test="ssrv.bootstrap",scale=sc,streamsize=stream)
  }
  a <- run(matrix(0,12,8)); shifted <- run(t(x)%*%b,0.03)
  expected <- simplify2array(lapply(1:8,function(s) matrix(b[,s],2,15)))
  expect_lt(max(abs(shifted$estimate-a$estimate-expected)),1e-10)
})

test_that("rank redraw diagnostics count rejected attempts and realize every draw", {
  f <- ssrv_fixture(); draw <- .ssrv_subject_indices; calls <- 0L
  local_mocked_bindings(.ssrv_subject_indices=function(groups) {
    calls <<- calls+1L
    if (calls==1L) rep(1L,12) else draw(groups)
  },.package="ALDEx3")
  fit <- aldex(f$Y,~treatment,f$data,nsample=4,test="ssrv.bootstrap",scale=tss.sm)
  expect_equal(fit$diagnostics$rejections.by.draw,c(1L,0L,0L,0L))
  expect_equal(fit$nsample,4)
  expect_equal(dim(fit$estimate)[3],4)
})

test_that("Gaussian blocks and custom coefficient callbacks use original log proportions", {
  f <- ssrv_fixture(); x <- t(model.matrix(~treatment+age,f$data))
  mu <- c(3,-0.5,0.2)
  run <- function(scale, ...) {
    set.seed(82)
    aldex(f$Y,x,nsample=12,test="ssrv.bootstrap",scale=scale,...)
  }
  for (covariance in c(FALSE,TRUE)) {
    if (covariance) {
      sigma <- matrix(c(1,0.1,0.02,0.1,0.2,0.01,0.02,0.01,0.05),3)
      sc <- compose.sm(gaussian.sm(rownames(x),mean=mu,covariance=sigma))
      legacy <- run(coefficient.sm,c.mu=mu,c.cov=sigma)
    } else {
      sc <- compose.sm(gaussian.sm(rownames(x),mean=mu,sd=c(0.1,0.2,0.05)))
      legacy <- run(coefficient.sm,c.mu=mu,c.sd=c(0.1,0.2,0.05))
    }
    composed <- run(sc)
    expect_equal(composed$estimate,legacy$estimate,tolerance=1e-10)
    expect_equal(composed$posterior,legacy$posterior,tolerance=1e-10)
    expect_identical(composed$scale.basis,"log2 composition")
  }
  custom <- run(function(A) setNames(rep(0,nrow(A)),rownames(A)))
  zero <- run(matrix(0,ncol(x),12))
  expect_equal(custom$estimate,zero$estimate,tolerance=1e-10)
  expect_lt(max(custom$estimate[1,,]),-1) # Detect accidental CLR centering.
})


test_that("zero-inclusive sign probabilities use all features in each family", {
  draws <- array(c(0,0,0,0, -1,0,1,1, 1,1,1,1, -1,-1,NA,-1),
                 c(1,4,4), dimnames=list("x",paste0("f",1:4),NULL))
  # Arrange four complete feature sequences explicitly.
  draws[1,,] <- rbind(c(0,0,0,0),c(-1,0,1,1),c(1,1,1,1),c(-1,-1,NA,-1))
  p <- .ssrv_probabilities(draws,"BH")
  expect_equal(unname(p$p.val[1,]),c(1,1,0,NA))
  expect_equal(unname(p$p.val.adj[1,]),p.adjust(c(1,1,0,NA),"BH",n=4))
  draws[1,2,] <- c(-1,1,1,1)
  expect_equal(.ssrv_probabilities(draws,"BH")$p.val[1,2],0.5)
  draws[1,2,1] <- Inf
  expect_true(is.na(.ssrv_probabilities(draws,"BH")$p.val[1,2]))
  family <- array(1, c(1,4,8))
  family[1,1,1] <- -1; family[1,2,1:2] <- -1
  family[1,3,] <- 0; family[1,4,1] <- NA
  scores <- .ssrv_probabilities(family,"BH")
  expect_equal(unname(scores$p.val[1,]), c(.25,.5,1,NA))
  expect_equal(unname(scores$p.val.adj[1,]), c(1,1,1,NA))
  f <- ssrv_fixture()
  set.seed(16)
  fit <- aldex(f$Y,~treatment,f$data,nsample=8,test="ssrv.bootstrap",scale=tss.sm,
               p.adjust.method="bonferroni",return.pars=c("estimate","p.val"))
  expect_identical(fit$p.adjust.method,"bonferroni")
  expect_null(fit$p.val.adj)
  expect_equal(summary(fit,ignore.intercept=FALSE)$p.val.adj,
               as.vector(t(apply(fit$p.val,1,p.adjust,method="bonferroni"))))
})
