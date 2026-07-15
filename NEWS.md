# ALDEx3 1.2.0

## New Features

- Added `method = "blmm"`: an ALDEx3-specific approximate mixed-effects engine using a batched profiled REML anchor fit per feature, draw-specific local covariance updates, and exact conditional GLS fixed-effect solves. Falls back to `lme4` for features where the approximation cannot be evaluated cleanly.
- On small datasets (with S\approx 20) blmm is approximately 40x faster than lme4, that factor should increase substantially as S increases to more realistic levels.
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

