# ALDEx3 1.1.0

## New Features

- Added `method = "blmm"`: an ALDEx3-specific approximate mixed-effects engine using a batched profiled REML anchor fit per feature, draw-specific local covariance updates, and exact conditional GLS fixed-effect solves. Falls back to `lme4` for features where the approximation cannot be evaluated cleanly.
- On small datasets (with S\approx 20) blmm is approximately 40x faster than lme4, that factor should increase substantially as S increases to more realistic levels. 
- Added feature-level parallelism (`n.cores`) for `method = "blmm"` and `method = "lme4"`.
- Added dedicated mixed-effects vignette covering model setup, BLMM formulation, validation guidance, and runtime comparison with exact `lme4`.


# ALDEx3 0.1.1

## User Facing Changes

- 1000x (approx) speed up in HC3 and HC0 standard error calculations. Note this was rate limiting before so this is a huge performance boost. 
- Fixed error with how the streamsize variable was calibrated. Leads to massive improvement in memory management. 


