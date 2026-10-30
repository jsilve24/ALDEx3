# ALDEx3 1.8.1

* Adaptive BLMM reuses one response-updatable TMB objective across accepted
  draws of a feature instead of recording a new derivative tape for every
  draw. Covariance refinement, observed REML curvature, inference checks and
  exact-fit rescues are unchanged. Failed draws discard the reusable object.

* BLMM now computes residual and contrast-variance sensitivities analytically
  from the native conditional solves, replacing repeated finite differences.
  The observed REML Hessian, degrees-of-freedom formula and exact-fit rescue
  checks are unchanged.

* The `lme4` engine now reuses a fitted model with `lme4::refit()` for subsequent
  response draws. Covariance parameters are still optimized separately, and
  Satterthwaite inference is recomputed for every draw and requested contrast.
  Singular starts or results, and refits that raise warnings or errors, use
  fresh fits to retain the full engine's boundary handling.

# ALDEx3 1.8.0

* Add `aldex(test="ssrv.bootstrap", method="lm")`: subject bootstrap Scale Simulation Random Variables (SSRVs), final absolute coefficient draws, posterior means and
  pointwise equal-tailed 95% intervals, including joint linear contrasts.
* Bootstrap SSRV summaries include raw two-sided sign scores and adjusted scores
  across all supplied features per coefficient or contrast (BH by default).
  Interval columns are `p2.5` and `p97.5`; posterior-only fits retain the scores.
* Composable and custom coefficient scales use original log2-proportion
  coordinates, consistent with `coefficient.sm`, without automatic CLR centering.
* Add `sparse.sm`, `gaussian.sm`, and `compose.sm` for exact named coefficient
  assignments. Sparse Gaussian-KDE modes support Laplace, ordinary feature
  bootstrap, and hybrid uncertainty selected by median final-effect interval width
  across features, with Laplace on ties. Each scale shift is one scalar across taxa.
  Bootstrap uses global paired subject resampling by default; Laplace uses the
  original subjects.
* Adapt existing scale callbacks/matrices by the same bootstrap fitting operator;
  preserve correlated Gaussian blocks and all legacy inference/deprecated args.
* Add explicit strata, bounded rank redraws and numerical/selection diagnostics.
  Zero curvature is reported as infinite conditional Laplace variance; no
  curvature cutoff, variance cap, or automatic bootstrap substitution is used.
  Unavailable posterior results are warned about and marked explicitly.
  Default bandwidth has no floor; constant profiles require explicit bandwidth.
  Existing plotting, p-adjust and effect diagnostics explicitly reject this path.
  Mixed-effects bootstrap inference is not supported.

# ALDEx3 1.7.0

## Changes to mixed-model results

- `method = "blmm"` now defaults to `test = "t.satterthwaite.adaptive"`. Each
  Monte Carlo draw refines its own covariance fit from a shared starting
  point (at most eight steps) and gets its own Satterthwaite degrees of
  freedom. Draws that fail the numerical checks are refit with `lmerTest`.
  **BLMM calls that omit `test` will give different estimates, standard
  errors, df and p-values than in earlier versions.**
- To reproduce the BLMM default of 1.5.0–1.6.2, pass `test = "t.satterthwaite"`
  (the shared-anchor approximation; cheaper than the adaptive default).
- `test = "t.residual"`, which reproduces the BLMM default of 1.4.0 and
  earlier, is deprecated and warns once per `aldex()` call.
- Refits are less frequent when a variance component is near zero: the
  adaptive start no longer forces whole-feature refits, and covariance
  parameters that are numerically zero are set to the boundary.
- Boundary fits get local inference without a refit when every covariance
  parameter is zero, or when whole columns of the covariance Cholesky factor
  are zero. Other boundary cases are refit with `lmerTest`.
- The mixed-effects vignette (`vignette("ALDEx3-mixed-effects")`) now
  describes the adaptive method and its Satterthwaite df, with the
  `oral_mouthwash_data` example, contrasts and measured runtimes.

## Upgrading from 1.4.0 (last CRAN release)

Versions 1.5.0–1.6.2 were not released on CRAN. Compared with 1.4.0:

- The BLMM default changed as described above.
- New `contrasts` argument to `aldex()` for testing linear combinations of
  fixed effects; results are stored in `fit$contrasts`.
- `summary()` includes requested contrasts with the coefficients and adds a
  `type` column to tell them apart. Use `component = "coefficients"` or
  `"contrasts"` to select one.
- New `aldex.p.adjust()` adjusts several coefficients or contrasts as one
  multiple-testing family without refitting.
- `df` (per-draw degrees of freedom) is now included in `return.pars` by
  default.
- `test = NULL` resolves per engine, and incompatible `method`/`test`
  combinations now fail before sampling.
- `Matrix (>= 1.6-0)` is now a runtime dependency.

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
  mathematical explanation.
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
