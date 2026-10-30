#!/usr/bin/env Rscript
# From the package root: Rscript developer/benchmark-mixed-engines.R NEW_OUTPUT [REPS=3]
# Optional ALDEX3_DREAM_LIBRARY; ALDEX3_BENCH_INPUTS_ONLY=1 generates inputs without fits.
args <- commandArgs(TRUE)
if (length(args) < 1L || length(args) > 2L) stop("Usage: benchmark-mixed-engines.R NEW_OUTPUT [REPS=3]")
reps <- if (length(args)==2L) as.integer(args[2]) else 3L
stopifnot(length(reps)==1L, is.finite(reps), reps>0L, file.exists("DESCRIPTION"))
out <- args[1]
if (dir.exists(out)) stop("Output directory already exists; use a new path")
dir.create(out, recursive=TRUE)
out <- normalizePath(out)
input_only <- identical(Sys.getenv("ALDEX3_BENCH_INPUTS_ONLY"),"1")
dream_lib <- Sys.getenv("ALDEX3_DREAM_LIBRARY")
if (nzchar(dream_lib)) .libPaths(c(dream_lib,.libPaths()))
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
if (requireNamespace("RhpcBLASctl",quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
}
started <- proc.time()[[3]]
# pkgload's development rebuild defaults to -O0. Build explicitly with release
# flags so native-engine comparisons do not measure a debug/release difference.
if (!input_only) pkgbuild::compile_dll(".", force=TRUE, debug=FALSE, quiet=TRUE)
pkgload::load_all(".", quiet=TRUE, compile=FALSE)
ns <- asNamespace("ALDEx3")
loading <- proc.time()[[3]]-started
use_dream <- nzchar(dream_lib) && requireNamespace("variancePartition",quietly=TRUE)
if (nzchar(dream_lib) && !use_dream) stop("variancePartition unavailable in requested library")
commit <- "275c054b1b2d314a0a5ff6b592d335bd8e76b827"
old_source <- system2("git",c("show",paste0(commit,":R/sr-mem.R")),stdout=TRUE)
if (!is.null(attr(old_source,"status"))) stop("Unable to retrieve pinned fresh-fit source")
dir.create(file.path(out,"source")); dir.create(file.path(out,"fits"))
writeLines(old_source,file.path(out,"source","sr-mem-baseline.R"))
old_env <- new.env(parent=ns); eval(parse(text=old_source),envir=old_env)
sources <- c(list.files("R",pattern="\\.R$",full.names=TRUE),"src/blmm.cpp","DESCRIPTION","data/oral_mouthwash_data.rda",
             "developer/benchmark-mixed-engines.R")
for (f in sources) {
  target <- file.path(out,"source",f); dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  file.copy(f,target,overwrite=FALSE)
}
source_md5 <- tools::md5sum(sources)
write.csv(data.frame(file=names(source_md5),md5=unname(source_md5)),file.path(out,"source-md5.csv"),row.names=FALSE)
writeLines(c(paste("baseline",commit),system2("git",c("rev-parse","HEAD"),stdout=TRUE),
             system2("git",c("status","--short"),stdout=TRUE)),file.path(out,"git.txt"))
if (file.exists("/proc/cpuinfo")) invisible(file.copy("/proc/cpuinfo",file.path(out,"cpuinfo.txt")))
capture.output(sessionInfo(),file=file.path(out,"session.txt"))
packages <- c("ALDEx3","lme4","lmerTest","Matrix","TMB","pkgload",
              if(use_dream)c("variancePartition","limma","BiocParallel"))
write.csv(data.frame(package=packages,version=vapply(packages,function(p)as.character(packageVersion(p)),character(1))),
          file.path(out,"versions.csv"),row.names=FALSE)
if (use_dream) writeLines(capture.output(packageDescription("variancePartition")),file.path(out,"dream-description.txt"))
if (use_dream) for (f in c("dream","run_lmm","run_lmm_on_gene","as_lmerModLmerTest2"))
  writeLines(deparse(get(f,envir=asNamespace("variancePartition"))),file.path(out,"source",paste0("dream-",f,".R")))
snapshots <- list.files(file.path(out,"source"),recursive=TRUE,full.names=TRUE)
write.csv(data.frame(file=snapshots,md5=unname(tools::md5sum(snapshots))),file.path(out,"snapshot-md5.csv"),row.names=FALSE)
phases <- list(data.frame(fixture="all",rep=NA_integer_,method="all",phase="loading",seconds=loading))
warnings <- errors <- timings <- comparisons <- list()
record_phase <- function(fixture,rep,method,phase,seconds) {
  phases[[length(phases)+1L]] <<- data.frame(fixture=fixture,rep=rep,method=method,phase=phase,seconds=seconds)
}
make_synthetic <- function(kind,N=240L,G=24L,D=3L,S=8L,seed=20261030L) {
  set.seed(seed); size <- N/G
  dat <- data.frame(g=factor(rep(seq_len(G),each=size)),x=rep(seq(-1,1,length.out=size),G),z=rnorm(N))
  rownames(dat) <- paste0("sample_",seq_len(N)); X <- model.matrix(~x+z,dat)
  Y <- array(NA_real_,c(N,D,S),dimnames=list(rownames(dat),paste0("feature_",seq_len(D)),paste0("draw_",seq_len(S))))
  for(d in seq_len(D)) {
    bi <- rnorm(G,sd=c(.8,1.3,2)[d]); bs <- .35*bi+rnorm(G,sd=c(.5,.8,1.2)[d])
    mu <- drop(X%*%c(2+d,.25*d,-.3))+bi[dat$g]
    if(kind=="correlated_slope") mu <- mu+bs[dat$g]*dat$x
    base <- mu+rnorm(N,sd=.7+.15*d)
    for(s in seq_len(S)) Y[,d,s] <- base+rnorm(N,sd=.12+.03*s)
  }
  C <- rbind(intercept_plus_half_x=c(1,.5,0),x_minus_z=c(0,1,-1)); colnames(C) <- colnames(X)
  list(logW=Y,data=dat,formula=if(kind=="correlated_slope")~x+z+(1+x|g) else ~x+z+(1|g),
       contrasts=C,provenance=list(type="synthetic Gaussian common-design response draws",seed=seed))
}
make_oral <- function() {
  e <- new.env(); load("data/oral_mouthwash_data.rda",envir=e); dat <- e$oral_mouthwash_data
  counts <- dat$counts; keep <- rowMeans(counts==0)<=.75
  pooled <- rbind(counts[keep,,drop=FALSE],other=colSums(counts[!keep,,drop=FALSE]))
  meta <- dat$metadata; form <- ~treat*timec+(1|participant_id)
  X <- model.matrix(reformulas::nobars(form),meta)
  valid_scale <- all(c("log_scale","log_scale_sd")%in%names(meta)) &&
    all(is.finite(meta$log_scale)) && all(is.finite(meta$log_scale_sd)) && all(meta$log_scale_sd>=0)
  set.seed(20261031L)
  sampled <- ns$aldex.sampler(pooled,t(X),16L,
    scale=if(valid_scale)ns$sample.sm else ns$clr.sm,
    scale.args=if(valid_scale)list(s.mu=meta$log_scale,s.sd=meta$log_scale_sd) else list(gamma=.5),return.pars="logW")
  selected <- c(head(setdiff(rownames(pooled)[order(rowSums(pooled),decreasing=TRUE)],"other"),8L),"other")
  Y <- aperm(sampled$logW[match(selected,rownames(pooled)),,,drop=FALSE],c(2,1,3))
  dimnames(Y) <- list(colnames(pooled),selected,paste0("draw_",1:16))
  C <- matrix(0,2,ncol(X),dimnames=list(c("antiseptic_minus_water_2h","antiseptic_change_2h"),colnames(X)))
  C[1,c("treatantiseptic-mouthwash","treatantiseptic-mouthwash:timec2hrs-later")] <- 1
  C[2,c("timec2hrs-later","treatantiseptic-mouthwash:timec2hrs-later")] <- 1
  list(logW=Y,data=meta,formula=form,contrasts=C,
    provenance=list(type="public packaged oral microbiome",seed=20261031L,zero_pool_threshold=.75,
      publication="Marotz et al. (2021), Quantifying live microbial load in human saliva samples over time reveals stable composition and dynamic load",
      doi="10.1128/mSystems.01182-20",publication_url="https://journals.asm.org/doi/10.1128/msystems.01182-20",
      qiita_study="11899 (acute perturbation)",ena_accession="ERP117149",
      upstream_analysis="https://github.com/knightlab-analyses/Saliva_quantification_study",
      input_note="Uses ALDEx3's packaged derived table; does not reprocess raw reads from Qiita/ENA",
      sampled_features=rownames(pooled),selected_features=selected,
      scale=if(valid_scale)"sample.sm: packaged log_scale means and log_scale_sd SDs" else "clr.sm gamma=.5"))
}
started <- proc.time()[[3]]
fixtures <- list(random_intercept=make_synthetic("RI"),correlated_slope=make_synthetic("correlated_slope"),
                 oral=make_oral(),large_random_intercept=make_synthetic("RI",N=2000L,G=200L,S=16L,seed=20261032L))
for(z in fixtures) stopifnot(all(is.finite(z$logW)),dim(z$logW)[1]==nrow(z$data),
  identical(colnames(z$contrasts),colnames(model.matrix(reformulas::nobars(z$formula),z$data))))
record_phase("all",NA_integer_,"all","input_generation",proc.time()[[3]]-started)
saveRDS(fixtures,file.path(out,"inputs.rds"))
write.csv(do.call(rbind,lapply(names(fixtures),function(n) {
  z <- fixtures[[n]]; data.frame(fixture=n,N=dim(z$logW)[1],D=dim(z$logW)[2],S=dim(z$logW)[3],formula=deparse(z$formula))
})),file.path(out,"inputs.csv"),row.names=FALSE)
writeLines(unname(tools::md5sum(file.path(out,"inputs.rds"))),file.path(out,"inputs-md5.txt"))
if(input_only) { write.csv(do.call(rbind,phases),file.path(out,"phases.csv"),row.names=FALSE); message("Inputs generated; no fits run: ",out); quit(status=0L) }
methods <- c("fresh_baseline","exact_refit","adaptive_blmm",if(use_dream)"dream")
fields <- c("estimate","std.error","df","p.lower","p.upper")
run_engine <- function(z,method) {
  if(method=="dream") {
    N <- dim(z$logW)[1]; D <- dim(z$logW)[2]; S <- dim(z$logW)[3]
    expr <- t(matrix(z$logW,N,D*S)); colnames(expr) <- rownames(z$data)
    rownames(expr) <- paste0("response_",seq_len(D*S))
    raw <- variancePartition::dream(expr,z$formula,z$data,L=t(z$contrasts),useWeights=FALSE,
      REML=TRUE,ddf="Satterthwaite",BPPARAM=BiocParallel::SerialParam())
    if(length(attr(raw,"errors")) || length(attr(raw,"error.initial"))) stop("dream response fitting errors: ",paste(attr(raw,"errors"),collapse="; "))
    names_all <- c(colnames(z$contrasts),rownames(z$contrasts))
    reshape <- function(v) array(t(v[rownames(expr),names_all,drop=FALSE]),c(length(names_all),D,S),dimnames=list(names_all,NULL,NULL))
    est <- reshape(raw$coefficients); se <- reshape(raw$stdev.unscaled*raw$sigma); df <- reshape(raw$df.residual)
    arrays <- list(estimate=est,std.error=se,df=df,p.lower=pt(est/se,df),p.upper=pt(est/se,df,lower.tail=FALSE))
    P <- ncol(z$contrasts); fit <- list()
    for(f in fields) { fit[[f]] <- arrays[[f]][seq_len(P),,,drop=FALSE]; fit[[paste0("contrast.",f)]] <- arrays[[f]][P+seq_len(nrow(z$contrasts)),,,drop=FALSE] }
  } else {
    fit <- if(method=="adaptive_blmm") ns$blmm(z$logW,z$formula,z$data,n.cores=1L,inference="adaptive",contrasts=z$contrasts) else
      (if(method=="fresh_baseline")old_env$sr.mem else ns$sr.mem)(z$logW,z$formula,z$data,n.cores=1L,method="lme4",mem.args=list(),contrasts=z$contrasts)
    raw <- fit
  }
  list(raw=raw,fit=fit)
}
run <- function(z,fixture,rep,method,phase="fit") {
  gc(); started <- proc.time()[[3]]
  answer <- tryCatch(withCallingHandlers({
    value <- run_engine(z,method); tfit <- proc.time()[[3]]
    record_phase(fixture,rep,method,paste0(phase,"_engine_and_extraction"),tfit-started)
    value$final <- lapply(c(coefficients="",contrasts="contrast."),function(pre)
      ns$aldex.pvals(value$fit[[paste0(pre,"p.lower")]],value$fit[[paste0(pre,"p.upper")]],"BH",FALSE))
    record_phase(fixture,rep,method,paste0(phase,"_aggregation"),proc.time()[[3]]-tfit)
    value
  },warning=function(w) {
    warnings[[length(warnings)+1L]] <<- data.frame(fixture=fixture,rep=rep,method=method,phase=phase,message=conditionMessage(w))
    invokeRestart("muffleWarning")
  }),error=function(e)e)
  elapsed <- proc.time()[[3]]-started
  if(inherits(answer,"error")) errors[[length(errors)+1L]] <<- data.frame(fixture=fixture,rep=rep,method=method,phase=phase,message=conditionMessage(answer))
  if(phase=="fit") timings[[length(timings)+1L]] <<- data.frame(fixture=fixture,rep=rep,method=method,seconds=elapsed,failed=inherits(answer,"error"))
  else record_phase(fixture,rep,method,phase,elapsed)
  answer
}
compare <- function(a,b,fixture,rep,method) {
  add <- function(actual,reference,component,quantity,scale=NULL,decisions=FALSE) {
    stopifnot(identical(dim(actual),dim(reference)),identical(dimnames(actual)[[1]],dimnames(reference)[[1]]))
    ok <- is.finite(actual)&is.finite(reference); dif <- abs(actual-reference)
    stat <- function(v,fun) if(length(v)) fun(v) else NA_real_
    scaled <- if(is.null(scale)) rep(NA_real_,sum(ok)) else dif[ok]/pmax(abs(scale[ok]),.Machine$double.eps)
    comparisons[[length(comparisons)+1L]] <<- data.frame(fixture=fixture,rep=rep,method=method,component=component,quantity=quantity,
      pairs=length(ok),finite_pairs=sum(ok),nonfinite_pairs=sum(!ok),reference_nonfinite=sum(!is.finite(reference)),actual_nonfinite=sum(!is.finite(actual)),
      median_abs=stat(dif[ok],median),p95_abs=stat(dif[ok],function(x)unname(quantile(x,.95))),max_abs=stat(dif[ok],max),
      max_scaled=if(is.null(scale))NA_real_ else stat(scaled,max),decisions_disagree=if(decisions)sum((actual[ok]<=.05)!=(reference[ok]<=.05)) else NA_integer_)
  }
  for(component in c("coefficients","contrasts")) {
    pre <- if(component=="coefficients")"" else "contrast."
    for(f in fields) {
      key <- paste0(pre,f); scale <- if(f=="estimate")b$fit[[paste0(pre,"std.error")]] else if(f%in%c("std.error","df"))b$fit[[key]] else NULL
      add(a$fit[[key]],b$fit[[key]],component,f,scale)
    }
    for(f in c("p.res","p.adj.res")) add(a$final[[component]][[f]],b$final[[component]][[f]],component,f,decisions=TRUE)
  }
  if(!is.null(a$fit$random.eff)) add(a$fit$random.eff,b$fit$random.eff,"random_covariance","random.eff")
}
warm <- fixtures$random_intercept; warm$logW <- warm$logW[,1,1:2,drop=FALSE]
for(method in methods) { message("Warmup ",method); run(warm,"warmup",0L,method,"warmup") }
for(fixture in names(fixtures)) for(rep in seq_len(reps)) {
  order <- if(rep%%2L==0L)rev(methods) else methods
  order <- order[((seq_along(order)+rep-2L)%%length(order))+1L]
  answers <- list()
  for(position in seq_along(order)) {
    method <- order[position]; message(fixture," repetition ",rep,"/",reps,": ",method)
    answer <- run(fixtures[[fixture]],fixture,rep,method)
    timings[[length(timings)]]$position <- position
    saveRDS(answer,file.path(out,"fits",paste0(fixture,"-",rep,"-",method,".rds")))
    if(!inherits(answer,"error")) answers[[method]] <- answer
  }
  reference <- answers$fresh_baseline
  if(!is.null(reference)) for(method in setdiff(names(answers),"fresh_baseline")) compare(answers[[method]],reference,fixture,rep,method)
}
write_table <- function(x,name) write.csv(if(length(x))do.call(rbind,x) else data.frame(),file.path(out,name),row.names=FALSE)
write_table(timings,"timings.csv"); write_table(comparisons,"comparison.csv")
write_table(warnings,"warnings.csv"); write_table(errors,"errors.csv"); write_table(phases,"phases.csv")
tt <- do.call(rbind,timings); good <- tt[!tt$failed,,drop=FALSE]
if(nrow(good)) write.csv(aggregate(seconds~fixture+method,good,function(x)c(median=median(x),min=min(x),max=max(x))),file.path(out,"timing-summary.csv"),row.names=FALSE)
unchanged <- identical(source_md5,tools::md5sum(sources))
writeLines(c(paste("Source files unchanged throughout:",unchanged),
  "Times include engine fitting, inference, common result extraction and directional BH aggregation.",
  "Excluded phases: package/data loading, input generation, warmup, gc, serialization and comparisons.",
  "Within-process repetitions on one shared machine; independent jobs must be paused externally.",
  "Boundary and other nonfinite/large discrepancies are retained in comparisons and saved arrays.",
  "Reference warnings are rows with method=fresh_baseline in warnings.csv."),file.path(out,"README.txt"))
if(!unchanged) stop("Source changed during benchmarking; timings are not a frozen-source comparison")
message("Finished: ",out)
