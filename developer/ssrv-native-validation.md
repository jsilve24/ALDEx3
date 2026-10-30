# Native SSRV result validation — 2026-10-29

ALDEx3 1.8.0 on `sparse-pim`, starting at
`a15e54969e6cda2c4014ff52918d667e9ac0c98e`, gives substantively similar
answers to Won's active two-group implementation on three full native inputs.
Both use the same basic estimator after the requested default changes. This
comparison concerns effects, uncertainty and selection, rather than identical
libraries or identical random numbers. No new coding defect was demonstrated.

This report supersedes the pending native-validation status and descriptions of
the old defaults in [the earlier stage report](ssrv-reference-validation.md).
That report retains the source inventory, manuscript/proof references and
controlled intermediate comparisons; its old subject-resampling and hybrid rules
are historical. The current method is described in [ssrv-bootstrap.md](ssrv-bootstrap.md).

## Inputs and comparison

The reference directory
`/home/jds6696/Downloads/sparse-DA-manuscript_code copy/manuscript_code copy`
was read only. Only function definitions were evaluated; analysis drivers with
file-writing side effects were not sourced. The active entry point was
`sparse_auto` in `sparse_SSRV.R`, SHA-256
`23d31d4c87c5b171e0a4cb856271db51541ab197adfb3d8fe59d1f1f32732582`.
The snapshot README records September 2026 and R 4.4.3; local validation used
R 4.6.1. The unavailable `dirmult` generator was replaced with
`LaplacesDemon::rdirichlet`, which implements the same normalized-Gamma
Dirichlet law. ALDEx3 was loaded from the current source, not the older installed
package.

| Native analysis | Features | Subjects (control / treatment) | Preprocessing reproduced |
| --- | ---: | ---: | --- |
| Zemb | 199 | 9 (5 / 4) | Remove no-feces sample; retain the 100 µl sample; mean raw count >= 5; remove internal standard Seq135 |
| Stämmler Pool1 versus Pool6 | 27 | 12 (6 / 6) | Native supplement parser and sample alignment; genus agglomeration; mean raw count >= 5 |
| Vandeputte control versus CD | 29 | 95 (66 / 29) | Remove first 40 input rows; genus agglomeration; count > 1 in strictly more than 20% of subjects |

These are complete retained inputs for those contrasts, not arbitrary feature
subsets. Counts, feature/sample order, labels and totals were independently
checked against the native parsers; Vandeputte aggregation was also reproduced
by a separate lineage-based summation. The Stämmler source parser emits 369
ragged-taxonomy-row warnings; these are input-parsing warnings, not fit failures.

Each dataset had three fits per implementation with 10,000 draws and separate
random streams. No latent compositions, bootstrap indices or normal draws were
forced to agree in these native comparisons. Each hybrid fit exposed both
completed candidates using observation-only wrappers, without changing RNG
consumption. Every fit used raw counts + 0.5, log2 proportions, treatment minus
control, global paired subject bootstrap, default `bw.nrd0`, and the native
median final-effect 95% width selector. Bootstrap and Laplace results were
compared separately as well as the selected hybrid. The reference's variance
cap and rejection behavior were left intact.

The earlier controlled harness separately verified composition reuse after
subject resampling, paired design rows, scalar feature-bootstrap shifts, scale
added once, log-KDE curvature and square-root Gaussian SD, and final-draw joint
contrasts and summaries. These native comparisons supplement those checks.

## Substantive results

Widths below are the range of median feature 95% interval widths across both
implementations and all three runs, in log2 units. The effect gap is the median
across features of the absolute difference between each implementation's
three-run average posterior means; it is a repeated-estimator diagnostic, not
a posterior interval formed by averaging quantiles.

| Analysis | Hybrid choice, both implementations | Bootstrap widths | Laplace widths | Hybrid effect gap |
| --- | --- | ---: | ---: | ---: |
| Zemb | Bootstrap, 3/3 runs | 2.009–2.040 | 1.924–1.948 | 0.00322 |
| Stämmler | Bootstrap, 3/3 runs | 1.734–1.777 | 1.424–1.457 | 0.00235 |
| Vandeputte | Laplace, 3/3 runs | 7.243–7.317 | 8.338–8.513 | 0.01639 |

Across all nine dataset/arm combinations, the largest feature mean gap after
averaging the three runs was 0.05865 log2 units. The largest feature-wise gap
divided by its average interval width was 0.752%. Bootstrap-only median mean gaps
were 0.00322, 0.00235 and 0.01163 for Zemb, Stämmler and Vandeputte; Laplace-only
gaps were 0.00361, 0.00571 and 0.01639. These small differences coexist with
ordinary between-run Monte Carlo variation and the documented numerical and
rejection differences; they do not establish exact distributional equality.

The designated Zemb spike (Seq1) and both Stämmler spike labels
(`AB076660XXXX`, `AB247615XXXX`) exclude zero in every arm, implementation and
run. Stämmler identifies exactly those two taxa throughout. For illustration,
the Zemb bootstrap spike intervals have lower endpoints 3.172–3.299 and upper
endpoints 8.899–9.035 across the six independent fits. Its Laplace endpoints
are 5.290–5.304 and 6.958–6.996.

Pointwise interval decisions are nearly, but not completely, identical:

- Zemb hybrid/bootstrap: 19 taxa in each reference run; 19, 18 and 19 in
  ALDEx3. The sole disagreement is Seq143 in run 2. The unused Laplace candidate
  has one or two disagreements per run (Seq151, Seq67, Seq180). All have an
  endpoint near zero: for example, Seq143's run-2 lower endpoints are +0.0437
  (reference) and -0.0068 (ALDEx3).
- Stämmler: no disagreements in any arm or run.
- Vandeputte bootstrap: the same eight taxa in every fit. Laplace/hybrid:
  reference 5, 5 and 6 versus ALDEx3 6, 6 and 6; the only disagreement is sp3.
  Its reference upper endpoints in runs 1 and 2 are +0.0377 and +0.0503;
  ALDEx3's are -0.0651 and -0.0865, against intervals about eight units wide.
  Won's own repeats therefore also change its exclusion decision.

These are pointwise 95% interval exclusions, not BH-adjusted discoveries. Won's
active code computes twice the smaller fraction of negative/positive effect
draws, applies BH across taxa, and uses q < 0.05 for discovery lists. The Zemb,
Stämmler and Vandeputte drivers all use those lists in their reported results
(`sparse_DA_09012026.R:942–955`, `sparse_DA_09012026.R:723–736`, and
`vandeputte_09012026.R:632–634,775–785`, respectively). These summaries are
used to construct the paper's discovery tables. ALDEx3 now reports these scores
as `p.val` and their adjusted values as `p.val.adj`, calculated after selection
from final absolute-effect draws. Each coefficient and each joint contrast has
a separate family across all supplied features; `aldex.p.adjust()` remains
unsupported for this method.

The approved zero-tie convention uses
`min(1, 2 * min(mean(theta <= 0), mean(theta >= 0)))`, so all-zero sequences
return one. Won uses strict `< 0` and `> 0`. No exact zeros occurred in the
tested treatment-effect draws, so this convention does not explain any native
comparison differences. Any nonfinite draw makes the score unavailable;
adjustment retains the declared family size rather than dropping missing taxa.

An additional Vandeputte comparison at the manuscript's 50,000 draws selected
Laplace in both implementations. Every one of the 29 feature interval decisions
agreed in each arm: eight excluded zero under bootstrap and six under
Laplace/hybrid. The borderline sp3 intervals were [-8.3899, -0.0155] (reference)
and [-8.3353, -0.0813] (ALDEx3), with means -4.18385 and -4.18360. The median
absolute feature-mean difference was 0.00591 for bootstrap and 0.00231 for
Laplace/hybrid; corresponding maximum differences were 0.02403 and 0.00894.
Median Laplace widths were 8.4890 and 8.4224. Actual scale means were 4.18261
and 4.18268. Neither implementation rejected a Laplace draw or reported a fit
warning/unavailable result. This run reused each implementation's replicate-1
seed; it increases draw count, rather than providing a fourth independent repeat.
It supports a Monte Carlo explanation for the borderline discrepancy without
establishing that all numerical/rejection differences are immaterial generally.

## Tail-score and BH validation

Feeding the same finite, nonzero completed draws to Won's `ssrv_arm_table`,
ALDEx3's reporting helper and an independent sign-count/sorted-BH calculation
matched within 1e-12 in every tested arm. Deterministic checks separately cover
the approved zero-inclusive rule, all-zero sequences, missing values and a fixed
full-family denominator. On small independent fixtures for all three uncertainty
modes, the reporting change preserved coefficient/contrast draws, interval
summaries, diagnostics and the final RNG state exactly.

Native comparisons used the same full inputs and each implementation's separate
streams. Seeds reuse the earlier comparisons; these are not additional independent
replicates. Zemb used its native 10,000 draws with three seeds; Stämmler's primary
Pool1/Pool6 contrast and Vandeputte used their native 50,000 draws.

| Analysis | Candidate | Reference BH discoveries | ALDEx3 BH discoveries |
| --- | --- | --- | --- |
| Zemb, three repeats | Bootstrap / selected hybrid | 2 / 4 / 3 | 4 / 4 / 4 |
| Zemb, three repeats | Laplace, not selected | 10 / 9 / 8 | 12 / 6 / 9 |
| Stämmler, 50,000 draws | Each of bootstrap, Laplace and hybrid | 2 | 2 |
| Vandeputte, 50,000 draws | Bootstrap, not selected | 5 | 5 |
| Vandeputte, 50,000 draws | Laplace / selected hybrid | 0 | 0 |

Stämmler and Vandeputte match feature identities as well as counts. The Stämmler
discoveries are the two designated spikes. All valid primary runs completed
without warnings, zero or nonfinite treatment-effect draws, or native variance-cap
exclusions. The Zemb seed-1 replay and seeds 2/3, plus both 50,000-draw analyses,
have clean process exits and passed the shared-draw checks.

Zemb's selected bootstrap lists agree on Seq1 and Seq47 throughout. In repeat 1,
ALDEx3 additionally calls Seq17 and Seq115; in repeat 3 it additionally calls
Seq115. Repeat 2 matches all four features. For Seq115, reference BH values are
0.05572 / 0.00995 / 0.07164, versus 0.02985 / 0 / 0.02985 in ALDEx3. In repeat 1,
this corresponds to seven versus three contrary-sign draws out of 10,000.
Monte Carlo variation demonstrably changes the reference decisions; the three
seeds do not establish that it is the sole cause of all differences.

A validation-only replay substituted the native 512-point grid mode for the
refined bootstrap mode on the same package stream. Selection and all four
bootstrap discoveries remained unchanged. Median feature width changed from
2.01508 to 2.02210, versus Laplace 1.93928. Thus the grid approximation does not
explain the selected-list discrepancy on that stream. Equal local raw scores can
still have different BH values because the rest of the family differs: Zemb
Laplace Seq47 has raw p=0.003 in both repeat-1 fits, but adjusted values 0.05427
and 0.04975.

The mode-location difference exceeded 0.1 log2 units in 2/10,000 draws and
0.5 in one draw. In the largest case (0.76267), direct Gaussian-mixture
evaluation and an independent derivative-root search confirmed that ALDEx3
found the higher smooth-density peak: modes -1.37481 and -0.60868 have log
densities -1.047436 and -1.047950. Won's grid selected near the second peak.
The inverse negative log-density curvature at the higher peak is 0.320115,
matching the reported conditional variance. This is evidence of closely
competing peaks and grid approximation, not a demonstrated refinement defect.

The temporary scripts, full score tables and process-exit records are under
`/tmp/ssrv-tail-validation/`. `run-v1.R` runs the native comparisons and
`run-grid-v3.R` records the grid-mode diagnostic and largest-difference profile;
`analyze-grid-profile.R` evaluates that profile. Primary comparisons use the
native draw counts and have recorded clean exits. Lower-draw Stämmler and
Vandeputte pilots are excluded from these primary results.

## Differences and scope

The requested corrections made global paired subject resampling the bootstrap
default, restored original-subject Laplace coefficients and KDE inputs, and
matched the native hybrid width rule including Laplace ties. Optional strata
apply only to subject bootstrap. Pure Laplace performs no subject or feature
bootstrap. These supersede the earlier deliberately different regression
Laplace branch and scale-variance selector.

Remaining differences are explicit:

- Won's feature-bootstrap mode is the maximum of the 512-point density grid;
  ALDEx3 refines the smooth Gaussian KDE. Both use Gaussian kernel SDs in log2
  coefficient units and re-estimate bandwidth after feature resampling.
- Won rejects conditional Laplace variances >= 100 and drops those draws.
  One Vandeputte run rejected 1/10,000; the other eight primary fits rejected
  none. ALDEx3 retains all finite conditional variances. No unbounded or
  unavailable package result occurred here. The prescribed Inf/warn/unavailable
  behavior remains covered by a constructed regression case, not inferred from
  these datasets.
- ALDEx3 bounds full-rank subject-bootstrap redraws. Primary fits rejected
  3/7/4 Stämmler and 57/74/60 Zemb resamples, and zero Vandeputte resamples.
  Each affected draw needed only one retry; no draw reached the attempt limit.
  For these binary designs, full rank means both groups occur, matching Won's
  successful-resample condition. General-design rank conditioning still needs
  consideration in any proof.
- The regression/composable-scale interface, Gaussian covariance blocks and
  joint contrasts extend the native two-group interface. Native results here
  validate the two-group case, not every possible regression design. Legacy
  inference and deprecated arguments remain supported; mixed effects are
  outside the new path.

There were no fit warnings or failures. Captured scale additions differed across
taxa by at most 3.56e-15, confirming a common scalar rather than independent
taxon perturbations on these fits. No new permanent comparison test suite or
numerical cap, floor, filtering rule or fallback was added.

## Scientific interpretation and remaining review

Both implementations estimate a positive Vandeputte Laplace scale shift
(means 4.148–4.194 log2 across the primary fits), while the supplied cytometry
contrast is -1.494. Applying the native read-share projection to the retained
taxa gives -1.628. This projection estimates retained totals from read fractions;
it is not a direct independent measurement of their absolute totals. The Zemb
projected measured contrast is +0.901 (whole-community +0.896).

Thus implementation agreement does not rescue a sparse assumption that appears
inappropriate for Vandeputte's dense changes. Release claims should distinguish
reproduction from calibration. These runs do not establish coverage, FDR control
or consistency. Author review remains needed for the following points; source
locations and hashes are in the earlier stage report:

- The manuscript/supplement's variance selector with bootstrap ties differs from
  the executable width selector with Laplace ties now followed here. The theorem
  and algorithm descriptions should be reconciled with the intended procedure.
- The proof uses within-group binary resampling, whereas the native code and
  current default use global paired resampling conditioned on both groups.
  General regression, bounded rank-conditioned redraws and composed blocks
  need appropriate assumptions; passing these comparisons does not supply them.
- The proof's bandwidth/mode assumptions need review for data-dependent nrd0
  re-estimated after feature resampling. S5 should also reflect the confirmed
  shared scalar Gaussian correction.

No assumptions were changed to resolve these questions. A full test-suite
usefulness audit remains a separate future TODO, not a completed audit. This
work did not regenerate every manuscript study, simulation or benchmark method.

## Reproducibility and package checks

The small temporary workspace is `/tmp/ssrv-native-validation/`: `data-prep.R`,
`verify-inputs.R`, `run-native.R`, `summarize-native.R`, native result RDS/CSV
files, independent summary tables and PNG/PDF scale-interval plots. Source/data
hashes and session information are recorded there. The reference source/input
SHA-256 manifest is `source-hashes.csv` (manifest hash
`18384430fcf40403e4b34afd4bae039e067b181f97bbcd3bccf9ca7d65ebcb85`).
The input RDS SHA-256 values are:

| Input | SHA-256 |
| --- | --- |
| Zemb | `d2aad8c9a7a6d6bb43c0cb26f278100048d704260aeabee2df460e3623783806` |
| Stämmler | `b87643f1d4eac61af388f725fd619a194348689062462ef965b11154a4c99e12` |
| Vandeputte | `f38b67c6749c4d2e2ab2e92f9e12120b2ec8fbd2a9ef9258df916e3011c6d7f5` |

The final current-source build and `NOT_CRAN=true R CMD check` returned
**Status: OK**: full tests (zero failures, warnings or skips), examples, all
vignette rebuilds and PDF manual. Logs are under
`/tmp/ssrv-tail-package-validation/`; source archive SHA-256 is
`834ce5dbf41098f22b0e822d9c34a61af5f63a4dff98dbc93c8ec77de1a6254f`.
All 75 packaged R/test/manual/vignette/native/metadata source files checked
matched the working tree byte for byte; DESCRIPTION differs only by standard
build formatting/generated fields. All 16 reference code/input hashes remained
unchanged. At completion of validation, branch/HEAD remained at the requested
starting commit, unrelated local changes were preserved, and no commit, push,
PR, merge, tag or release had been made.

The additional Vandeputte fit is saved as
`results/vandeputte_S50000_rep1.rds` in the temporary workspace, with matching
feature/scale CSVs and log. No release action was performed.
