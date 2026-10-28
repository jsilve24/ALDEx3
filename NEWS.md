# Experimental branch: adaptive boundary handling

- Keep the default/shared-anchor and residual BLMM paths unchanged. Adaptive
  refinement alone uses bounded natural relative Cholesky coordinates.
- For a stationary all-zero covariance fit with validated raw curvature and
  feasible objective checks, use the cached QR calculation and the existing
  first-order Satterthwaite limit, N-p, instead of a full lmerTest refit.
- At a validated partial boundary, remove complete zero Cholesky columns from
  the Satterthwaite information calculation while retaining uncertainty from
  residual variance and nonzero covariance parameters. Check decoupling,
  curvature, stationarity and predicted inference changes before acceptance.
- Keep exact rescue for unsupported, unresolved and invalid-curvature fits.
  This option remains experimental; numerical equivalence does not establish
  statistical calibration at a data-selected variance boundary.

# ALDEx3 1.6.2

## Review fixes

- Add opt-in BLMM `test="t.satterthwaite.adaptive"`: per-draw centering,
  at most eight safeguarded covariance iterations, local Satterthwaite
  inference, and exact rescue when numerical checks fail. The shared-anchor
  default remains unchanged.
- Avoid duplicate TMB report evaluations and cache invariant design work;
  batch adaptive exact rescues while retaining individual exact fits.
- Avoid repeated wide QR decompositions when validating contrasts against a
  full-rank design, substantially reducing overhead for cell-level models.
- Keep coefficient plots separate from stored contrasts with the same name.
- Prevent failed exact refits of rejected BLMM draws from retrying the whole
  feature. Exact fitting errors identify the affected feature and draw indices.
- Retain valid fitted quantities when lmerTest returns invalid numerical
  contrast inference, warning and returning `NA` p-values. Missing inference
  propagates through Monte Carlo aggregation; fitting errors still stop.
- Include the shared-anchor Satterthwaite derivation at the end of the
  mixed-effects vignette.
- Update release metadata and documentation, add repository and issue-tracker
  URLs, and require testthat >= 3.2.0 for the mocking helpers used by the tests.

# ALDEx3 1.6.1

## Grouped multiple-testing adjustment

- Add `aldex.p.adjust(fit, groups = c("age", "genotype"))` to adjust jointly
  across all features of selected coefficients or stored contrasts, without
  refitting. A list of character vectors specifies several disjoint families;
  unlisted rows retain separate families. Coefficients and contrasts can share
  a family; qualified names distinguish rows with the same name.
- Reconstruct directional probabilities from stored per-draw estimates,
  standard errors, and degrees of freedom, then adjust within each draw and
  direction before Monte Carlo aggregation. Streaming fits are supported;
  fits that omitted these arrays cannot be readjusted. Existing fitting and
  default per-vector adjustments are unchanged.
- Record `p.adjust.method` and `onesided` in new fits. For older fits lacking
  these settings, the helper defaults to BH and two-sided inference; pass
  nondefault settings explicitly. The returned object records its adjustment
  families, and `summary()` uses the updated results.

# ALDEx3 1.6.0

## Mixed-model computation

- BLMM now caches the design and response cross-products used repeatedly during
  fitting, avoiding a large dense tensor containing covariance-parameter,
  random-effect, and observation basis terms. Sparse and dense factorizations
  use the same profiled REML objective and conditional GLS calculations,
  including correlated slopes and crossed effects. Inference defaults and the
  criteria for refitting failed approximations with `lmerTest` are unchanged.
- Sparse solves reuse the symbolic structure of the factorization and account
  for its permutations. Production anchor reports omit large random-effect
  matrices. `Matrix` is now an explicit runtime dependency (>= 1.6-0).
- Before allocating dense factors, BLMM estimates factorization and solve work
  from model size, random-effect factor fill-in (extra nonzeros created during
  factorization), and draw count. Internal `ALDEx3.blmm.backend` (`"auto"`,
  `"dense"`, or `"sparse"`) and `ALDEx3.blmm.memory` (estimated working bytes)
  options support reproducible validation. An optional
  `ALDEx3.blmm.diagnostics` callback receives the selected backend and its
  rationale once per BLMM invocation or streaming chunk. Memory estimates
  include worker and automatic-differentiation (AD) costs; they are not a hard
  process-memory limit.
- On the saved IBDMDB pilot (three taxa, eight draws each, 1,297 observations),
  automatic BLMM took a median 0.100 s versus 1.203 s for fresh `lmerTest`
  fitting, about 12 times faster. Three repetitions used varied engine order;
  each BLMM run refitted the same single draw for one taxon with `lmerTest`
  (one of the 24 fits).
  The IBD reference emitted convergence/Hessian warnings for Bacteroides; this
  timing comparison does not establish that those reference fits are reliable.
- On the saved large simulation pilot (one taxon, 300 draws, 6,000 observations,
  300 random-intercept groups), median automatic BLMM time was 0.680 s versus
  17.790 s for one fresh `lmerTest` run, about 26 times faster. Forced dense and
  sparse BLMM medians were 1.949 s and 0.601 s. Both benefit from cached matrix
  products; sparse factorization provides additional savings on this design.
- These Linux measurements use identical saved draws, one worker, and one BLAS
  thread; fitting includes time spent refitting rejected draws with `lmerTest` and excludes input
  loading, design parsing, and compilation. BLMM medians use three repetitions.
  They are pilot results, not universal speed guarantees or a completed IBDMDB
  reanalysis. High factor fill-in and frequent full `lmerTest` refits can reduce
  the benefit. Tests of numerical results and derivatives support equivalence to
  the existing method. They do not establish how well its approximate p-values
  are calibrated across different study designs.
- Correct the `streamsize` documentation to describe the existing heuristic of
  ten times the raw array size (the default corresponds to about 800 MB of response
  draws). Streaming behavior is unchanged; its threshold is not a hard memory
  limit, and changing BLMM chunk boundaries can change approximate results.

# ALDEx3 1.5.0

## Summary presentation

- `summary(fit)` now includes coefficients and any requested contrasts, with a
  `type` column when combining them. Explicit `component = "coefficients"` and
  `component = "contrasts"` filters remain available. Fits without requested
  contrasts retain the previous column schema. Intercept filtering applies only
  to coefficients; stored results and multiple-testing adjustments are unchanged.

## Inference and linear hypotheses

- `method = "blmm"` now defaults to `test = "t.satterthwaite"`, the
  shared-anchor Satterthwaite approximation. Coefficient-specific degrees of
  freedom can change p-values even when estimates and standard errors barely
  change. `test = "t.residual"` retains the previous BLMM convention, including
  Satterthwaite inference for draws refitted with `lmerTest`.
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
  from the `lmerTest` Satterthwaite inference used by the `lme4` engine and its
  fallback refits.
- Made BLMM reject invalid or ill-conditioned raw anchor Hessians instead of
  ridge-regularising and inverting them. Draw-specific Newton updates are now
  checked against their profiled REML objective, and rejected draws are
  selectively refitted with full `lmerTest` fits rather than returning
  pathological uncertainty.

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
