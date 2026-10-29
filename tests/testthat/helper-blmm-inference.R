# Independent datasets; treatment is assigned at group level when requested.
inference_fixture <- function(seed = 731, groups = 8, size = 5,
                              design = "between", noise = 0.02, S = 4,
                              effect = 0.6, weak = FALSE) {
  set.seed(seed)
  sizes <- if (design == "unbalanced") rep(c(3, 5, 8), length.out = groups) else rep(size, groups)
  g <- rep(seq_len(groups), sizes)
  x <- if (design %in% c("between", "unbalanced", "nested"))
    rep(rep(c(-0.5, 0.5), length.out = groups), sizes) else
      unlist(lapply(sizes, function(n) seq(-1, 1, length.out = n)))
  dat <- data.frame(g = factor(g), x = x, block = factor(ceiling(g / 2)))
  a <- rnorm(groups, sd = if (weak) 0.05 else 1)
  b <- 0.4 * a + rnorm(groups, sd = 0.6)
  y <- effect * x + a[g] + rnorm(length(g), sd = 0.6)
  form <- ~ x + (1 | g)
  if (design == "slope") {
    y <- y + b[g] * x
    form <- ~ x + (1 + x | g)
  }
  if (design == "nested") {
    y <- y + rnorm(nlevels(dat$block), sd = 0.7)[dat$block]
    form <- ~ x + (1 | block) + (1 | g)
  }
  Y <- matrix(y, length(y), S) + matrix(rnorm(length(y) * S, sd = noise), length(y), S)
  list(Y = Y, logW = array(Y, c(length(y), 1, S)), data = dat, formula = form)
}

feature_fit <- function(z, ...) {
  p <- ALDEx3:::blmm_parse_formula(z$formula, z$data)
  do.call(ALDEx3:::blmm_fit_feature, c(list(d = 1, logW = z$logW,
    X = p$X, basis = p$basis, is_log = p$is_log, phi_init = p$phi_init,
    lower = p$lower, reTrms = p$reTrms, formula = z$formula, data = z$data), list(...)))
}

# Test oracle: the average per-draw profiled REML objective, which the batched
# TMB kernel minimises for the shared anchor.
blmm_anchor_objective_from_report <- function(report, n, p) {
  mean(ALDEx3:::blmm_profiled_objectives(report, n, p))
}

# aldex.mem.sim() fixtures: true log scales, sampled logW (N x D x S), and
# a BLMM/exact-lme4 pair fitted to the same draws.
blmm_true_scale <- function(sim, nsample) {
  matrix(rep(log2(colSums(sim$W)), nsample), nrow = ncol(sim$Y))
}

blmm_sample_logW <- function(Y, formula, data, nsample, scale, scale.args = list()) {
  X <- t(model.matrix(reformulas::nobars(formula), data))
  draws <- ALDEx3:::aldex.sampler(
    Y, X, nsample, scale, scale.args, return.pars = "logW"
  )
  aperm(draws$logW, c(2, 1, 3))
}

blmm_fit_pair <- function(logW, formula, data, inference = "adaptive") {
  list(
    blmm = ALDEx3:::blmm(logW, formula, data, n.cores = 1L,
                         inference = inference),
    lme4 = ALDEx3:::sr.mem(logW, formula, data, n.cores = 1L,
                           method = "lme4", mem.args = list())
  )
}
