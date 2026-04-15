# Unreleased

## Changes

- Added feature-level parallelism to the BLMM mixed-effects engine and aligned the user-facing docs with the implemented behavior.
- Removed stale rendered vignette outputs and source-side build artifacts from the repository.
- Tightened the BLMM anchor objective and fallback handling so the mixed-effects vignette and tests stay internally consistent.

# ALDEx3 0.1.1

## User Facing Changes

- 1000x (approx) speed up in HC3 and HC0 standard error calculations. Note this was rate limiting before so this is a huge performance boost. 
- Fixed error with how the streamsize variable was calibrated. Leads to massive improvement in memory management. 


