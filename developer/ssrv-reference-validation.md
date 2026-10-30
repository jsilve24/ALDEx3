# Bootstrap SSRV reference validation — 2026-10-29

## Summary

Validated ALDEx3 1.8.0 on `sparse-pim`, starting at
`a15e54969e6cda2c4014ff52918d667e9ac0c98e`, against the supplied Won snapshot.
No demonstrated SSRV coding defect was found in the comparisons below. No R
implementation, statistical rule, legacy argument, or mixed-effects path was
changed. The vignette now describes the actual reference findings and corrects
its interpretation of an interval above zero. These checks establish numerical
agreement of specified stages, not consistency, coverage, or FDR control.

### Sources and reproducibility

The reference root was
`/home/jds6696/Downloads/sparse-DA-manuscript_code copy/manuscript_code copy`.
It was read only; no analysis driver was sourced or run there. Its README calls
this the September 2026 manuscript rebuild and records R 4.4.3. This validation
uses R 4.6.1 locally, rather than assuming the snapshot's reproduction claims.

SHA-256, with paths relative to that root:

| File | SHA-256 |
| --- | --- |
| `sparse_SSRV.R` | `23d31d4c87c5b171e0a4cb856271db51541ab197adfb3d8fe59d1f1f32732582` |
| `README.md` | `f90d734fd662187802ae76e35280e740dfbe9ccb0b4aa37982521e55e170fa35` |
| `count_input_policy.R` | `466279f4ab4bdf30c879a6f089b3479bb70a31ebfc7ef27ee10d2a5b618f3f84` |
| `2020_zemb_microOpen_spike/sparse_DA_09012026.R` | `8205796edd6d9891b9d5fa2358871f0b1d6248e9934ba151a04fc08905a7232d` |
| `2020_zemb_microOpen_spike/zemb_counts.csv` | `a316d75c18fafcbf6a9a50da2c88c19575be21c4e43e651f84ebbdd838cb1fd4` |
| `2020_zemb_microOpen_spike/zemb_metadata.csv` | `8cb53ce9d871a12997d154b1a5ac1c8a1ef70284a6e6dfc68c300bf1e614f381` |
| `simulation_rerun/src/sparse_SSRV_rerun.R` | `f013cf40cc314ce9c053f1a7af3aa07580a60fba90d1827ac34f256684c55ab8` |
| `simulation_rerun/src/run_methods_rerun.R` | `f387320d2c86e011263833bfb94d419ae9539cbf5b6b5614323f93cf10974d5b` |

Accompanying material was located separately under
`/home/jds6696/Dropbox/Research/sparse-da/docs/overleaf-aoas-resubmission`.
This is additional local manuscript evidence, not a claim that its revision is
identical to the supplied code snapshot. The inspected files were:

| File | SHA-256 |
| --- | --- |
| `main.tex` | `c83d1b7f57bd6073b2cb2a1442a7dea445859e6c1fb9b3a199e31b4962e50bf0` |
| `supplementary_body.tex` | `534284839309ba128c2b5392964037b5059b1b17411267359f65d0ce9adafa61` |
| `consistency_proof.tex` | `4514d87e79ee240b9022560982978d75e6cfd54a2f98a4238d2c505dc4c13497` |
| `algorithm_bootstrap.tex` | `a5e824beb276215afe66724db41eedac612b43ee4d6df822a95299a2750e7179` |

The small, temporary harness is in `/tmp/ssrv-reference-validation/`:

```sh
rtk proxy Rscript /tmp/ssrv-reference-validation/compare.R
rtk proxy Rscript /tmp/ssrv-reference-validation/end-to-end.R
rtk proxy Rscript /tmp/ssrv-reference-validation/boundary.R
```

It evaluates only reference function definitions in an isolated environment and
clones the active `sparse_auto` loop, returning its intermediates before selection
and reporting. The native Gaussian KDE/grid calculations remain intact. Wrappers
inject compositions, subject indices, feature indices, and standard normals;
`var_max=Inf` disables the reference's rejection threshold for these finite-variance
fixtures. This adaptation is explicit: it aligns resamples but does not validate
equality of the original resampling distributions or full hybrid procedures.

## Major issues: executable reference versus the agreed method

### Entry points and measurement model

* Real data: `run/01_fit_real_data.sh` dispatches six study scripts; these source
  `sparse_SSRV.R` and call `sparse_auto` and bandwidth sensitivity helpers. In the
  Zemb script see lines 120 and 793–814. Study A instead uses
  `simulation_rerun/src/run_methods_rerun.R:172–193` and the rerun's
  `sim_sparse_auto`. The legacy `sim_sparse_auto` in the top-level file is a
  different function, not the active Study A implementation.
* Core active computations require base/stats functions. The driver supplies
  `rdirichlet` and broad analysis dependencies; the README records `dirmult`,
  `LaplacesDemon`, `numDeriv`, ALDEx2, and numerous benchmark/plot packages.
  Active `kde_mode_var` uses analytic derivatives; old functions still require
  numerical Hessians. Full manuscript reproduction additionally needs the
  pinned Julia/MiMIX environment and other benchmark packages; it was not run.
* Raw example inputs include Zemb counts/metadata/qPCR CSVs, Stämmler study
  supplements, Prasse RDS/metadata, Jin spreadsheets, Vandeputte OTU/flow-cytometry
  inputs, and ALDEx2's SELEX data plus the supplied truth file. Study A generates
  SparseDOSSA2 data. The harness uses an explicit computational subset of raw
  Zemb data, not the manuscript's complete analysis or its feature filters.
* Drivers draw `Dirichlet(raw counts + 0.5)` for each original subject, store an
  S × D × N array of **log2 proportions**, and only then bootstrap subjects.
  Treatment is 1, control 0, and coefficients are treatment minus control.
  There is no CLR transform. The separate `Y=raw+1` matrix in study scripts is
  not the SSRV measurement input (`count_input_policy.R:1–47`).

### Bootstrap and KDE

`sparse_SSRV.R:608–638` samples one subject-index vector across all features,
resamples design labels with those subjects, and reuses a subject's stored
composition when its index repeats. It samples all N subjects globally, retrying
until both groups occur. ALDEx3 defaults to fixed-size within-group resampling
for a simple binary design; its explicit single-stratum option instead requests
paired global sampling with bounded rank redraws. These distributions differ.

The reference next resamples D coefficients to estimate one negative KDE mode,
then adds this scalar to the **original** feature coefficients. ALDEx3 has the
same feature-resampling role. The reference bootstrap uses the maximum of the
512-point `stats::density()` grid; ALDEx3 searches/refines the smooth Gaussian
KDE. Both default to `bw.nrd0`, re-estimated on the resampled profile. A numeric
`density(bw=...)` and ALDEx3 `bandwidth` are Gaussian kernel SDs in coefficient
units here; `adjust` or another library's dimensionless bandwidth factor would
need conversion. They cannot be interchanged by matching a numeric value alone.

### Laplace and hybrid

The active reference uses original-subject group means in its Laplace branch
(`sparse_SSRV.R:642–687`). `kde_mode_var:220–235` starts from the density grid
maximum and refines within two grid cells, using `kde_derivs:165–172` to compute
`-1 / (f''/f - (f'/f)^2)`. The Gaussian SD is its square root. Its active default
is **one scalar shared across taxa**. The old `sparse_auto_old`,
`sparse_auto_sens_multi`, and legacy simulation functions retain independent
taxon noise and/or Hessians of interpolated log density. Neither those functions
nor stale comments describe the active path reliably.

ALDEx3 deliberately includes subject resampling in **both** candidates. This is
an extension, not a missing attempt to reproduce the plug-in Laplace branch.
Absolute effects equal log-proportion coefficients plus total-scale coefficients
once, with one shift per coefficient/draw. Gaussian blocks retain the absolute
total meaning of `coefficient.sm`; no automatic CLR centering was introduced.

The snapshot's executable default is `scale_rule="width"`: select by the larger
median delivered 95% effect-interval width, with Laplace winning ties
(`sparse_SSRV.R:744–785`). Comments saying the default is `majority` are stale.
Other selectable rules compare variance contributions including covariance
terms; none is the agreed all-draw sample variance of the scalar shifts. A small
counterexample gives opposite choices: 99 bootstrap zeros plus one 10 have
sample variance 1 and 95% width 0, whereas 100 alternating Laplace values ±0.5
have sample variance 0.252525 and width 1. ALDEx3 correctly chooses bootstrap
under the agreed variance rule; the width rule chooses Laplace.

The reference rejects conditional variances at or above `var_max=100` and drops
failed Laplace draws. The real-data function stops with fewer than two surviving
draws; the rerun simulation can substitute bootstrap. Those behaviors are not
copied. ALDEx3 retains Inf in diagnostics, warns, and makes the affected entire
coefficient sequence unavailable; it neither samples an infinite-SD Gaussian nor
invents a variance or selects bootstrap automatically.

Reference effect tables report posterior means and 95% quantiles of their
retained draws, plus strict-sign tail counts, BH-adjusted values and decisions
(`ssrv_arm_table:258–303`). Separate scale summaries can use analytic Gaussian
mixture quantities. ALDEx3 reports actual final joint-draw means/quantiles and
contrasts, without adding those tail/BH quantities or replacing effect summaries
with analytic scale-mixture summaries.

## Suggested tests: focused checks actually performed

No permanent tests were added because the reference comparisons demonstrated no
new SSRV coding defect. Existing SSRV assertions cover shared shifts, scale
meaning, hybrid ties, joint contrasts, bounded rank redraws, streaming/retention,
and explicit unavailable results. Temporary comparisons were run independently
by the main agent after inspecting the harness.

| Comparison                                                                                                                         | Numerical evidence (maximum absolute error unless stated)          |
| ---                                                                                                                                | ---                                                                |
| Actual `rLogDirichlet` versus normalized Gamma(counts+0.5), shared RNG                                                             | `3.55e-15` in log2 proportions                                     |
| Synthetic fixture: shared latent compositions, subject and feature indices; Won group means versus `fflm`                          | `4.44e-15`                                                         |
| Synthetic native reference grid bootstrap final effects                                                                            | `4.89e-15`                                                         |
| Synthetic native grid versus ALDEx3 smooth bootstrap mode                                                                          | `0.00645821` log2 units                                            |
| Synthetic aligned smooth bootstrap modes, using Won `kde_mode_var`                                                                 | `1.90e-5`                                                          |
| Original plug-in Laplace: matched input, distinct numerical mode searches                                                          | mode `2.02e-5`; variance `5.05e-5`; shared-Z shifts `3.99e-5`      |
| Subject-bootstrap Laplace inputs aligned to the extension                                                                          | mode `2.03e-5`; variance `1.67e-4`; final shared-Z draws `4.50e-5` |
| Won analytic curvature evaluated at ALDEx3's same mode and bandwidth                                                               | variance `3.33e-16`, relative error `7.52e-16`                     |
| Zemb subset, native grid reference coefficients/shift/final draws                                                                  | at most `5.67e-15`                                                 |
| Zemb native grid versus smooth modes; aligned smooth modes                                                                         | `0.01188634`; `1.93e-5`                                            |
| Actual bootstrap entry point versus independent QR coefficients and assembled final draws/contrasts/summaries, both factor codings | at most `2.00e-14`                                                 |

The synthetic fixture has D=30, N=8, S=40. The real-count fixture takes the first
50 raw feature rows and tubes 1,21,22,23 versus 4,7,10,13, S=30; metadata confirms
the control/exposed order. Each original subject's alpha is checked before row
resampling. Repeated subject indices reuse that original latent draw. Feature
indices contain duplicates while output names/order retain the original 50
features. The same fixtures check a shared scalar shift, reversed factor coding,
and joint coefficient/contrast means and 0.025/0.975 quantiles.

On the Zemb fixture, candidate shift sample variances are 0.207733030 (bootstrap)
and 0.308033233 (Laplace); hybrid selects the entire Laplace sequence. Reversed
coding is also checked against its corresponding QR calculation, without
asserting pathwise equality under an unchanged Z sequence. Exact equal-variance
sequences select bootstrap. Expected final arrays use the package KDE helper;
the separate Won derivative/mode comparisons above independently check that
helper's formula and distinguish its numerical mode search.

One Zemb profile has `bw.nrd0=0.209275847166` in log2 units and
`0.145058963422` in natural-log units. Scaling bandwidth and mode by log(2), and
variance by log(2)^2, agrees to floating-point precision. The observed smooth-mode
differences reflect the two optimizers' tolerances/searches; at the same mode the
curvature formulas agree. This is not a general global-mode guarantee or a bound
on discrepancies for arbitrary multimodal data.

The boundary harness injects valid compositions (2/3,1/3) in controls and
(1/3,2/3) in treatment, obtaining the exact profile (-1,+1) through the actual
`fflm` and KDE. At h=1, pure Laplace and hybrid each warn once, preserve Inf and
four unbounded draws in diagnostics, and make the entire treatment sequence
unavailable. The intercept and a contrast with zero treatment weight remain
finite; hybrid's branch is NA. Explicit bootstrap stays finite without warning.
At h=1.000000001, conditional variance is 499999959.63 and is accepted without
warning. Controlled normal perturbations are zero in this boundary fixture;
random shared-Z sampling is checked separately above. No KDE/profile mock or
numerical cutoff is used.

Local package checks:

* Installed current 1.8.0 into `/tmp/ssrv-local-validation/library`; the user's
  installed 1.3.1 was not overwritten. Platform: Linux x86_64, R 4.6.1,
  lme4 2.0.6, Matrix 1.7.5, TMB 1.9.22, testthat 3.3.2.
* Focused SSRV tests: 213 assertions, zero failures, warnings or skips.
* Full `NOT_CRAN=true` suite: 1,632 assertions passed, no failures or skips.
  `test_local` emitted four existing BLMM PSOCK worker warnings that
  `'package:ALDEx3' may not be available when loading`. An earlier sandboxed
  run had three socket-creation errors; the rerun with local sockets enabled
  resolved those errors. These were outside the SSRV path.
* All vignettes rebuilt with `R CMD build`. The final bootstrap vignette was
  rendered again after the last prose correction (Pandoc warned that `--mathml`
  is deprecated). A staged final archive excludes only the unrelated existing
  `R/.agent-shell` transcripts; those files remain untouched in the checkout.
  All 65 R/test/manual files in that archive match the repository byte for byte.
  Archive: `/tmp/ssrv-final-validation/ALDEx3_1.8.0.tar.gz`, SHA-256
  `ac54361a98b830e76ae8a7f88991f507cbe26a012c6a912a97a9c1cd182c9f3f`.
* Final `NOT_CRAN=true R CMD check ALDEx3_1.8.0.tar.gz`: **Status: OK**, including
  installation, examples, all vignette rebuilds, and the PDF manual (no
  `--no-manual` option). Installed-package tests report
  `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1632 ]`. Logs are in
  `/tmp/ssrv-final-validation/check.log` and
  `/tmp/ssrv-final-validation/ALDEx3.Rcheck/`. The first, unstaged package check
  also passed all tests without warnings and completed the manual/vignettes;
  its sole NOTE concerned the preserved unrelated transcript directory.
* `git diff --check` passed. Branch/HEAD remain `sparse-pim` at the requested
  starting commit; no commit, push, PR, merge, tag, or release was made. Won's
  source hashes remained unchanged.

## Documentation updates

The bootstrap vignette now distinguishes initial development from subsequent
reference validation, identifies the grid/smooth mode difference, and states the
reference's sampling, selection and rejection differences. Its interval-above-zero
explanation now refers to posterior draws, rather than the simulated truth.
The earlier audit's historical checks remain historical; this report supersedes
its statement that original code and proof material had not yet been inspected.

## Nice-to-have improvements: methodological questions for author review

These remain questions, not proposed implementation changes:

1. **Reconcile the source versions.** `main.tex:407` and the active supplement's
   hybrid theorem (`supplementary_body.tex:1320–1397`) compare scale variances and
   give bootstrap ties. The executable snapshot instead defaults to delivered
   interval widths and Laplace ties. Which code/results should accompany that
   manuscript statement? The agreed ALDEx3 rule is retained.
2. **Scope of the extension.** The supplement proves a within-group bootstrap
   result for binary mean differences, and its Laplace theorem uses original
   subjects. Applying subject resampling to Laplace changes its limiting regime;
   general regression and rank-conditioned bounded redraws require separate
   assumptions. In particular, conditional unbiasedness of binary bootstrap
   means does not by itself establish unbiasedness after general rank rejection.
3. **Bandwidth and mode assumptions.** The supplement's KDE perturbation argument
   holds a common h_D fixed for fixed D, assumes h_D→0 and D h_D²→∞, and maximizes
   on a fixed compact interval (`supplementary_body.tex:645–744, 870–994`). Actual
   nrd0 is data dependent and re-estimated per feature bootstrap; an explicit
   fixed bandwidth need not shrink with D. ALDEx3 maximizes numerically over the
   kernel-center range, without a user-specified compact interval. Transfer of
   the proof, including boundedness/uniform-integrability steps, needs review.
4. **Means versus concentration and finite Monte Carlo.** The active theorem
   assumes integrability for Laplace and finite posterior variances for hybrid.
   Inverse log-density curvature measures density width and need not vanish as
   feature count grows; it is not automatically the sampling variance of the
   estimated mode. The theorem concerns posterior means with exact component
   moments, not a finite-S sample-variance decision, posterior concentration, or
   coverage/FDR. Infinite/undefined approximations are outside those assumptions.
5. **Pseudocode revisions.** Active S5 still prints independent vector Gaussian
   noise (`supplementary_body.tex:1712`), despite the author correction and active
   code using a shared scalar. Its bootstrap algorithm also needs explicit
   pairing of resampled labels/design rows; executable code does pair them.
   The separate `algorithm_bootstrap.tex` is older and even presents subject
   resampling inside a feature loop. The supplement's active proof uses one
   joint within-group subject resample instead.

No automatic centering, filtering, curvature threshold, bandwidth floor, variance
cap, substitute candidate, new significance rule, or consistency claim was added.
