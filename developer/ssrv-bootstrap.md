# Bootstrap SSRV implementation (1.8.0)

`aldex(test="ssrv.bootstrap", method="lm")` dispatches to a private draw-major
procedure; legacy inference is unchanged. `fflm(..., .coefficients.only=TRUE)`
uses the existing solver but returns before residual uncertainty calculations.
The design is built once. Original-subject Dirichlet draws precede joint row
resampling. Scale is applied exactly once; contrasts are linear transformations
of final absolute coefficient draws.

Public scale API: `sparse.sm(coefficients, uncertainty="hybrid", bandwidth=NULL)`;
`gaussian.sm(coefficients, mean=0, sd=NULL, covariance=NULL)`;
`compose.sm(...)`. Every exact design coefficient name must occur once. The
intercept must be assigned explicitly. Separate Gaussian blocks are independent;
within-block covariance is retained. Sparse blocks share feature indices per draw.
Gaussian block draws and Laplace perturbations are independent conditional on
the shared compositional coefficients and subject resampling; sparse feature
resampling is shared. No scale model reweights the compositional draws.
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
For complete candidates, hybrid compares their all-draw sample variances; exact
ties choose bootstrap. Unresolved hybrid branches are NA, not bootstrap. The
posterior summary helper never removes failed draws to summarize a smaller set.
Contrasts with zero weight on unavailable coefficients retain valid marginal
draws; affected contrasts are unavailable. Summary adds status for such fits.
Both candidates include subject resampling, deliberately extending the original
Laplace branch. This extension and the numerical mode search are not certified
by an examined consistency proof. See ssrv-audit.md for the full choice inventory.

Simple intercept-plus-binary designs default to within-group resampling. Other
designs use paired row resampling; `bootstrap.strata` is an explicit length-N
vector. Rank-deficient draws are redrawn at most `bootstrap.max.attempts` times
per requested draw (default 100 total attempts). No draw is silently discarded.

Large sample arrays are processed one draw at a time for deterministic RNG order
regardless of streamsize. Coefficients require O(P D S) working memory even for
posterior-only results, to compute exact quantiles after hybrid selection. Only
P x S candidate shifts are stored. No residual t uncertainty, p-values, or BH are
computed. Posterior intervals are pointwise equal-tailed 95% quantiles. Existing
inference-dependent helpers reject this path; summary and contrasts are supported.

No manuscript PDFs or original reference code were supplied. Validation uses
independent algebraic fixtures, not claimed reproduction of unavailable code.

Hybrid reference dependence: a fixed offset per coefficient preserves both
candidate absolute effects and variance ordering. Draw-dependent offsets preserve
each candidate's absolute effects but can reverse scale-variance ordering through
different covariance terms. Selection is defined in original log-proportion
coordinates; do not insert an automatic CLR transform.
