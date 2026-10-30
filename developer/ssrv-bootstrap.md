# Bootstrap SSRV implementation (1.8.0)

`aldex(test="ssrv.bootstrap", method="lm")` dispatches to a private draw-major
procedure; legacy inference is unchanged. `fflm(..., .coefficients.only=TRUE)`
uses the existing solver but returns before residual uncertainty calculations.
The design is built once. Original-subject Dirichlet draws precede any joint row
resampling. Scale is applied exactly once; contrasts are linear transformations
of final absolute coefficient draws.

Public scale API: `sparse.sm(coefficients, uncertainty="hybrid", bandwidth=NULL)`;
`gaussian.sm(coefficients, mean=0, sd=NULL, covariance=NULL)`;
`compose.sm(...)`. Every exact design coefficient name must occur once. The
intercept must be assigned explicitly. Separate Gaussian blocks are independent;
within-block covariance is retained. Sparse blocks share feature indices per draw.
Gaussian block draws and Laplace perturbations are independent conditional on
their compositional coefficients; sparse feature resampling is shared. The
bootstrap and Laplace candidates use their respective resampled and original-subject
compositional baselines. No scale model reweights the compositional draws.
A joint custom `compose.sm` block is intentionally not needed: the full `scale`
may instead be a function with formals `A, X, logComp, indices, draw`, returning
one named coefficient vector; this retains joint dependence of arbitrary kind.
Existing sample-scale functions (without an `A` formal) retain their signature;
their original-sample draws are projected by coefficients-only fflm on the same
bootstrap rows. All models use raw log2-proportion coefficients as the explicit baseline;
legacy adapters add projected sample scale once. There is no automatic CLR
centering. Gaussian blocks describe total-scale coefficients as in coefficient.sm;
zero intercept scale shift fixes reference total abundance at one. Custom A
callbacks also receive log-proportion coefficients.

KDE uses Gaussian kernels and bw.nrd0 re-estimated for each feature bootstrap,
without a bandwidth floor. Constant profiles require explicit positive bandwidth.
Coordinate centering and bandwidth scaling are undone in the returned mode and
curvature variance; they do not change the log-proportion regression basis.
The mode grid targets h/8 spacing and includes all kernel centers without a fixed
coarsening cap. Analytic-score root refinement is used where bracketed; otherwise
bounded density optimization refines the peak. Only exact computed-density ties
choose the leftmost peak. Tolerances are in bandwidth units.

Conditional Laplace variance is the inverse negative log-density curvature:
positive precision is not thresholded, zero gives Inf, negative/undefined precision
is undefined. No variance cap or automatic bootstrap fallback is permitted.
An unbounded/undefined Laplace candidate has no completed finite draw sequence:
warn once and mark affected coefficient results unavailable, retaining Inf in
conditional-variance diagnostics. Never fabricate an empirical sequence variance
or compare it as if proper infinite-variance Gaussian draws had been generated.
For complete candidates, hybrid compares, per coefficient, the median across
features of each candidate’s all-draw q.975 minus q.025 final-effect widths.
It selects the entire candidate with the larger median; exact ties choose Laplace.
Unresolved hybrid branches are NA, not bootstrap. The
posterior summary helper never removes failed draws to summarize a smaller set.
For hybrid coefficients, `diagnostics$bootstrap.interval.width` and
`laplace.interval.width` record the two median feature widths, and
`selection.rule` identifies the width rule. The existing `bootstrap.variance`
and `laplace.variance` fields retain their meaning as scale-shift sample
variances; they no longer determine selection. `laplace.resampling` records
the original-subject policy when a Laplace candidate is requested.
Contrasts with zero weight on unavailable coefficients retain valid marginal
draws; affected contrasts are unavailable. Summary adds status for such fits.
Laplace uses original-subject compositions and regression for both scale and
final effects. The feature-bootstrap candidate uses subject-resampled regression
and feature resampling. If at least one sparse block is present and all use
Laplace, no subject or feature bootstrap is drawn, and accompanying Gaussian
coefficients use the original-subject baseline. Its resampling diagnostic is
`original subjects`. Gaussian-only, callback, and legacy scale paths use subject
bootstrap; mixed sparse modes and hybrid use both bootstrap and original-subject
Laplace candidates. See ssrv-audit.md for the historical choice inventory and
its supersession notice.

When a subject bootstrap is used, all designs default to global paired row
resampling; `bootstrap.strata` is an explicit length-N vector controlling only
actual bootstrap subject resampling. Rank-deficient draws are redrawn at most
`bootstrap.max.attempts` times per requested draw (default 100 total attempts). No draw is silently discarded.

Large sample arrays are processed one draw at a time for deterministic RNG order
regardless of streamsize. Coefficients require O(P D S) working memory even for
posterior-only results, to compute exact quantiles after hybrid selection.
Candidate shifts use P x S storage, with additional original-subject coefficient
baselines; overall working memory remains O(P D S). No residual t uncertainty is computed. Final draws yield two-sided sign scores
`min(1, 2 * min(mean(theta <= 0), mean(theta >= 0)))`; nonfinite draws make
the score unavailable. `p.adjust` uses all D features for each coefficient or
contrast, with `n=D` and BH by default. The reference implementation uses strict sign
tails (`< 0` and `> 0`), so its scores differ at exact zero ties; ALDEx3
counts zeros in both tails and reports 1 when every draw is zero. Posterior intervals are pointwise
equal-tailed 95% quantiles. Existing
inference-dependent helpers reject this path; summary and contrasts are supported.

Reference comparisons and scientific limitations are documented in
[ssrv-reference-validation.md](ssrv-reference-validation.md) and
[ssrv-native-validation.md](ssrv-native-validation.md).

Hybrid selection uses final-effect widths, incorporating covariance between each
compositional baseline and its scale shift. Gaussian priors refer to total-scale
coefficients on the original log-proportion basis; do not insert automatic CLR
centering.
