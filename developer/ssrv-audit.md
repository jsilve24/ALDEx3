# Statistical-procedure audit: bootstrap SSRVs

Historical implementation audit. On 2026-10-29, the procedure was revised:
global paired subject bootstrap is the default for all designs, optional strata
apply only to bootstrap, and Laplace uses original subjects. Hybrid selects per
coefficient the entire final-effect candidate with the larger median across
features of all-draw 95% interval widths, with Laplace on exact ties. These
decisions supersede conflicting descriptions in the historical inventory below
and the earlier [reference validation report](ssrv-reference-validation.md).
Native end-to-end validation on three full datasets is recorded in
[ssrv-native-validation.md](ssrv-native-validation.md); the earlier reference
report records comparisons with Won's supplied code and locally available
manuscript and proof sources. Statements below about unavailable reference
material describe the earlier audit, not the current validation status.

## Future work

- [ ] Audit the full test suite for useful behavioral coverage, redundant
  assertions, brittle implementation coupling, independent reference calculations,
  and unnecessary runtime. Consolidate without weakening numerical or statistical
  checks. The cleanup so far covers SSRV tests and selected redundancies in three
  pre-existing test files; it is not a comprehensive suite audit.

This report supersedes the earlier audit's recommendation of a curvature cutoff
and forced bootstrap fallback. Those were changes to the procedure, not merely
numerical repairs. No consistency claim is made from passing tests. The original
proofs and reference implementation have not been checked against this extension.

## Restored hybrid rule

For available Laplace and feature-bootstrap candidate sequences, compare their
sample variances over all global draws. Use the entire larger-variance sequence
for that coefficient; bootstrap wins an exact tie. No chunkwise or drawwise
selection, curvature cutoff, variance cap, or forced branch substitution remains.

The curvature-based conditional Laplace variance and the empirical variance of a
completed candidate sequence are different quantities. Zero negative log-density
curvature produces an infinite conditional variance. That value is preserved in
`laplace.max.conditional.variance`, with affected iteration counts in
`laplace.unbounded.draws`. There is no proper infinite-SD Gaussian from which to
sample a finite candidate sequence. The program therefore warns and reports
unavailable results; it does not invent an empirical variance, manufacture draws,
or silently choose bootstrap. For hybrid the unresolved branch is NA. Explicit
pure Laplace also warns and leaves affected results unavailable.

Unavailable coefficient arrays/summaries are NA for the whole requested sequence;
no failed draws are deleted to summarize a shorter sequence. `summary()` adds a
status column, diagnostics name unavailable coefficients/contrasts, and
`available.draws` records zero for affected coefficients. Unaffected coefficients
and contrasts that do not involve unavailable coefficients remain available.
The warning suggests explicitly requesting bootstrap or reconsidering the sparse
model. This is failure reporting, not a replacement inferential procedure.

## Inventory of choices

| Choice | Origin | Action and statistical implications |
| --- | --- | --- |
| Automatic CLR centering | Assistant implementation choice | Removed earlier. All fits start from log2 proportions; Gaussian scale coefficients concern absolute totals. CLR changes Gaussian prior interpretation and can change hybrid scale-variance ordering. |
| Forced bootstrap on unavailable Laplace | Assistant choice permitted, but not required, by the initial request; explicitly rejected in the later discussion | Removed. No candidate substitution. The initial implementation permission was not evidence that the modified rule inherited any proof. |
| Fixed sqrt(machine epsilon) curvature cutoff | Assistant addition during the preceding audit | Removed. Every representable positive curvature-derived variance is used without this threshold. Zero gives Inf; undefined/negative precision remains undefined. |
| Fixed 1e-6 bandwidth floor on nondegenerate profiles | Assistant choice; the original request allowed an explicit bandwidth or documented degenerate-data floor, not a mandatory clamp on every estimate | Removed. Default nondegenerate bandwidth is bw.nrd0, re-estimated for feature bootstraps, with no fixed lower bound preventing bandwidth from shrinking. |
| Automatic bandwidth for constant profiles | Assistant fallback | Removed. Constant profiles, including a constant feature resample, require an explicit positive bandwidth. No default smoothing uncertainty is silently created. This may require explicit bandwidth with very few features. |
| Treating peaks within 1e-12 log density as ties | Assistant numerical policy | Removed. Only exact equality of computed densities is treated as a tie; the leftmost computed maximum wins. A near-equal but higher peak is not deliberately replaced by a lower one. |
| Fixed cap of 16385 regular grid points | Assistant numerical/resource choice | Removed. The grid targets spacing h/8 and includes every kernel center. An unrepresentable grid errors rather than being silently coarsened. Large range/bandwidth ratios can be expensive. |
| Bandwidth-unit root/optimizer tolerance | Numerical approximation, not a different prior | Retained and documented. Root refinement uses the analytic score when bracketed, otherwise bounded density optimization, in standardized KDE coordinates. Its location tolerance scales back with h rather than imposing a fixed coefficient-unit precision floor. Finite-precision mode-search agreement with reference code remains unverified. |
| Centering and scaling coordinates inside KDE | Algebraic numerical transformation | Retained. The mode is transformed back and curvature variance restored to original units. This is not CLR centering of the regression response. |
| bw.nrd0 bandwidth and re-estimation in feature resampling | Documented implementation choice authorized by the initial request, not verified against reference code | Retained, without a floor. Bandwidth changes the estimator and must be aligned in the later reference comparison; it is not a theorem-free implementation detail. |
| Bounded redraw of singular bootstrap designs | Explicit initial requirement | Retained. This conditions the successful resampling distribution on full rank; it is not identical to an unconditional row bootstrap on designs where singular resamples have appreciable probability. Errors after the bounded attempt count, with diagnostics; no pseudoinverse or coefficient deletion. Proof assumptions must account for this. |
| Within-group resampling for simple binary designs; paired rows otherwise | Explicit initial requirement | Retained with explicit optional strata. General-formula strata are not inferred. No cluster inference is claimed. |
| Subject resampling in both hybrid candidates | Explicit deliberate regression extension in the initial request | Retained and documented as different from the manuscript's original Laplace branch. Existing consistency proofs cannot automatically be asserted to cover this extension. |
| Gaussian block independence and shared sparse feature indices | Explicit composable-interface specification and documented choices | Retained. Joint Gaussian covariance within a block is preserved; no coordinate is overwritten to discard covariance. Custom joint coefficient draws remain supported. |
| Half-count measurement convention | Existing package behavior and explicit requirement | Retained. No feature filtering, reference-taxon substitution, or change to the multinomial/Dirichlet measurement model. |
| Legacy covariance symmetry/PSD tolerances and MVN sampling | Existing scale helpers, not introduced by this feature | Reused unchanged. They are not new statistical repairs; reference comparisons should still use the same covariance and sampling conventions. |
| Exact posterior quantiles and default 95% level | Explicit reporting specification | Retained over all completed final draws. No averaged chunk quantiles, extra residual t errors, p-values, BH on posterior probabilities, or marginal-endpoint contrast intervals. |

These choices are now explicit. Documenting them does not establish that the
paper's assumptions or consistency theorems cover them. There is no claimed
formal guarantee for the finite-precision global-mode search.

## Checks and scope

The exact synthetic boundary fixture uses two kernel centers -1 and +1 with
bandwidth 1. Analytic score refinement finds mode zero and zero precision, so
conditional variance is Inf. A nearby positive-precision fixture with bandwidth
1.000000001 has finite variance about 5e8 and is not thresholded away. This is a
constructed analytic example, not evidence that exact flatness commonly occurs
in real biological data. No real-data frequency estimate has been established.

Tests also cover larger-variance selection, bootstrap ties, absence of forced
fallback, raw-proportion Gaussian/legacy equivalence, shared shifts, exact
summaries, and unavailable versus unaffected contrasts. They remain in existing
test blocks; no new test files or general testing framework were added.

Earlier independent checks matched two-group bootstrap draws to separately
written gamma sampling, group means and direct-density mode optimization within
5.3e-9, and compared 100 randomized KDE fixtures with a dense-grid reference.
Those are computational checks, not consistency proofs or comparison against the
student's implementation. Updated checks on this revision:

* Independent two-group bootstrap draws and summaries matched within 4.443111e-9.
* 100 randomized KDE fixtures matched an independent dense-grid/refinement
  reference to a maximum log-density deficit of 2.220446e-16; standardized
  finite-difference curvature checks passed for resolvable positive precision.
* The real flat-KDE boundary was exercised through the inference path for pure
  Laplace and hybrid: Inf retained, one warning, no bootstrap substitution,
  unavailable results explicit, unaffected contrasts finite. Near-flat positive
  precision remained finite without warnings. Explicit bootstrap produced finite
  results, and streaming/retention choices agreed.
* A D=80, N=32, P=3, S=128 hybrid smoke fit took 0.862 seconds; returned size was
  290808 bytes with draws versus 38496 for summaries only. Summaries and
  diagnostics matched across streaming/retention choices. This is a smoke check,
  not a performance or memory-bound guarantee.
* R 4.5.0 on Debian 13: `NOT_CRAN=true R CMD check --no-manual` returned
  **Status: OK**, including installation, examples, and all vignette rebuilds.
  All tests were enabled: **1632 assertions passed, zero failures, warnings,
  or skipped tests**. No new test files or test blocks were introduced.
* The final source archive matches current R, test, manual and vignette files;
  `git diff --check` passed.

## Remaining scientific review and reference comparison

The user's methodological review and comparison with the forthcoming original
code remain necessary. Align log base, count prior, design coding, subject and
feature resamples, bandwidth choice/re-estimation, mode ties, curvature formula,
shared scalar perturbation, and all-draw variance comparison. Begin with the
bootstrap-only two-group case. Do not assume equality of the original full hybrid
with the deliberately extended regression hybrid.

The previous small simulation's observed coverage was imperfect and is not a
coverage/FDR guarantee. No release, push, tag, merge or PR submission was made.
