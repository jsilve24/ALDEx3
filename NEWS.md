# ALDEx3 1.5.0

## Inference and linear hypotheses

- `method = "blmm"` now defaults to `test = "t.satterthwaite"`, the
  shared-anchor Satterthwaite approximation. Coefficient-specific degrees of
  freedom can change p-values even when estimates and standard errors barely
  change. `test = "t.residual"` retains the previous BLMM convention, including
  Satterthwaite inference for draws rescued with `lmerTest`.
- `test = NULL` resolves by engine: `"t.HC3"` for ordinary linear models,
  `"t.satterthwaite"` for BLMM, and native `"t"` inference for `lme4` and
  `nlme`. Explicit incompatible choices now fail before sampling. The two
  new selectors are available only with BLMM. Results record `method` and `test`.
- `contrasts` accepts a named numeric matrix of scalar linear hypotheses at
  fitting time. Contrasts use the full within-fit coefficient covariance and
  the original composition and scale draws. Results are stored separately in
  `fit$contrasts`; use `summary(fit, component = "contrasts")` to summarize them.
  The `return.pars` fields control coefficient and contrast results together;
  degrees of freedom are now included by default.
- Each coefficient and contrast keeps its own adjustment across features,
  within each draw and direction before Monte Carlo averaging. There is no
  joint FDR claim for combined discovery lists. Contrast weights need not sum
  to zero, and contrast rows need not be independent.
- Rewrote the mixed-effects vignette around a runnable age-by-genotype example,
  interpretation, scale assumptions, a compact evidence summary, and a separate
  mathematical explanation. The maintainer's adoption decision supersedes the
  experimental-status recommendations in the preserved investigation reports.
  The full publication benchmark remains incomplete; no completed IBD aggregate
  validation or universal accuracy, calibration, or acceleration claim is made.
- Preserve singleton array axes when combining streaming results, and preserve
  the coefficient dimension in intercept-only mixed-model fits.

# ALDEx3 1.4.0

## Maintenance

- Calculate both t-distribution tails directly in BLMM and the exact mixed-model
  wrapper, preserving small upper-tail probabilities previously lost by subtraction.
- Add internal experimental shared-anchor Satterthwaite inference for all BLMM
  coefficients, matched-input regression tests, and a reproducible development
  study. Production defaults remain unchanged: accepted BLMM draws still use
  residual df. Shared-anchor inference is not promoted pending resolution of
  accepted-draw errors and stronger calibration evidence. Clarify this distinction
  from the `lmerTest` Satterthwaite inference used by the `lme4` engine and rescue.
- Made BLMM reject invalid or ill-conditioned raw anchor Hessians instead of
  ridge-regularising and inverting them. Draw-specific Newton updates are now
  checked against their profiled REML objective, and rejected draws fall back
  selectively to exact `lme4` rather than returning pathological uncertainty.

# ALDEx3 1.3.1

## Maintenance

- Registered TMB native routines and disabled dynamic symbol lookup for
  portable compiled-code checks.
- Made the ALDEx2 compatibility regression test robust to the R 4.7.0
  binomial RNG change by checking overall p-value agreement and matching
  significance decisions instead of a strict elementwise tolerance.

# ALDEx3 1.3.0

## User Facing Changes

- Harmonized independent scale-uncertainty inputs on standard deviations:
  `sample.sm()` now accepts `s.sd`, and `coefficient.sm()` now accepts `c.sd`.
- Added explicitly named `s.cov` and `c.cov` arguments for covariance matrices.
  The legacy `s.var`, `s.cor`, and `c.cor` arguments retain their original
  semantics with deprecation warnings until ALDEx3 2.0.0.
- Clarified that `gamma` is a standard deviation in `clr.sm()` and `tss.sm()`,
  and documented that sample-level and coefficient-level priors have different
  covariance structures even when their SD inputs are numerically equal.
- Added a two-group coefficient-scale example explaining how factor reference
  levels determine the sign of a group offset.

## New Features

- Added `method = "blmm"`: an ALDEx3-specific approximate mixed-effects engine using a batched profiled REML anchor fit per feature, draw-specific local covariance updates, and exact conditional GLS fixed-effect solves. Falls back to `lme4` for features where the approximation cannot be evaluated cleanly.
- Initial small-dataset benchmarks showed substantial speedups. Later measurements
  found that dense large-design computations and fallback can eliminate them;
  runtime gains should not be extrapolated from Monte Carlo draw count alone.
- Added feature-level parallelism (`n.cores`) for `method = "blmm"` and `method = "lme4"`.
- Added dedicated mixed-effects vignette covering model setup, BLMM formulation, validation guidance, and runtime comparison with exact `lme4`.

# ALDEx3 1.0.3

## User Facing Changes

- Added `aldex.plot()` for pairwise ALDEx3 model contrasts, with volcano,
  effect, MA, and waterfall plot modes.
- Added `aldex.effect()` for ALDEx2-inspired effect diagnostics on a single
  binary ALDEx3 contrast. The function reports the mean contrast estimate,
  mean pooled within-group standard deviation, mean Cohen's d, and an
  ALDEx2-style directional `overlap` diagnostic.
- Updated `cohensd()` documentation and preserved its legacy return value: a
  feature by Monte Carlo sample matrix of Cohen's d draws.
- Effect diagnostics and MA plots now use reconstructed log abundance,
  `logComp + logScale`, and clearly report that they require stored
  Monte Carlo arrays.
- Fixed the Quickstart vignette setup so package data such as
  `gut_crohns_data` is available during vignette rebuilds.

# ALDEx3 0.1.1

## User Facing Changes

- 1000x (approx) speed up in HC3 and HC0 standard error calculations. Note this was rate limiting before so this is a huge performance boost. 
- Fixed error with how the streamsize variable was calibrated. Leads to massive improvement in memory management. 
