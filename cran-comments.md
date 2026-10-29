# ALDEx3 1.7.0 release candidate

This candidate updates CRAN version 1.4.0. NEWS records contrasts, selectable
mixed-model inference, sparse computation, grouped adjustment and refit/plot
fixes, together with bounded adaptive inference. BLMM now defaults to
the adaptive selector. Explicit
`t.satterthwaite` retains shared-anchor behavior. `t.residual` is deprecated:
it warns once per call but retains its fitting and inference behavior, with
no scheduled removal version. No df formula, clipping or rescue-rule change
was introduced by promotion.

## Local environment and exact archive

- Arch Linux x86_64, R 4.6.1 (2026-06-24), GCC/G++ 16.2.1, C++17.
- Fresh source build, including both vignettes, followed by fresh
  `R CMD check --as-cran` on that archive in an isolated directory/library.
- Standard `-g -O2` temporary Makevars avoid distribution-specific compiler
  flags. The package Makevars specifies only C++17.
- HTML Tidy 5.8.0 and V8 8.2.1 were provided from temporary check locations.
  These are validation tools, not package dependencies.
- Archive: `ALDEx3_1.7.0.tar.gz` (302,970 bytes).
- SHA256: `8fc59bf1c1c1f81d2f4efc0942a17dd23369cf74bdccf4787f6b8a802ef130c2`.
- Pre-promotion archive/logs are retained in
  `release/1.7.0-pre-adaptive-default-20261028/`.
- Archive inspection found no private investigation inputs, diagnostics,
  working libraries or compiled object files.

## Documentation review

The mixed-model vignette now presents the adaptive analysis first, then its
algorithm, local inference, warnings and resource guidance. Alternative selectors,
engines and historical comparisons are at the end. README, Quickstart and help
were reviewed together; mathematical definitions and stopping/boundary details
were checked against source. All runtime R expressions are unchanged from
`b2c93e7`. The previous checked archive is preserved in
`release/1.7.0-pre-docs-20261029/`.

## Check results

Final promoted-default `R CMD check --as-cran`: **OK**, zero errors, warnings
and notes. Tests reported 1,438 passing expectations with zero failures,
warnings or skips.
Examples including donttest, both vignette rebuilds, and PDF/HTML manuals passed.
The source archive and corresponding logs/checksums are preserved locally in
`release/`; earlier 1.6.2 artifacts were preserved before replacing current notes.
The ordinary installed user package was not overwritten.

## Numerical validation and default decision

All numerical solver function bodies and native source were verified against
the frozen, previously tested implementation. The public API change only
selects adaptive when `test` is omitted and warns for explicit residual mode. Matched real-input engine
comparisons use identical response arrays and full intended anchor blocks.
The independent null assessment uses explicitly synthetic counts on a retained
covariate design: 600 independent datasets, three variance regimes and 32 draws.
It is a limited calibration screen rather than universal validation.

Adaptive and full lmerTest agreed on the primary and family-level rejection
decisions in every paired dataset. The weak-variance primary rejection count
was 14/200 for both (7%; two-sided 95% interval 3.88–11.47%). Its uncertainty
did not clear the original absolute default-promotion screen. The initial
candidate retained the default. On 2026-10-29 the maintainer explicitly revised
the release judgment and authorized adaptive as default based on its paired
agreement with full fits and bounded computational cost. The absolute screen
remains failed; no study, outcome, cutoff or numerical method was changed to
make it pass. This result does not establish inflation or that the older
approximation is statistically safer. Calibration remains limited and unresolved;
boundary matching is not described as proof of calibration.
Convergence warnings from two full-reference calls and approximation-rescue
notifications were retained in the diagnostic evidence; they are not failures
or warnings from the package check above.

## Reverse dependencies

The CRAN package database queried during this preparation lists version 1.4.0
and no direct reverse dependencies through Depends, Imports, LinkingTo or
Suggests. No reverse-dependency check was therefore required for listed CRAN
packages. This does not cover unlisted downstream projects.

## Outstanding external checks — not submission-ready yet

Windows, macOS and atlas/Linux R-devel remain pending for this exact candidate.
The maintainer prefers established local checks plus the existing R-hub GitHub
workflow (`.github/workflows/rhub.yaml`). The candidate branch must be published
before that workflow can check it; no push or dispatch has been authorized here.
Verify the requested R versions and actual job logs, not just platform names.

No Win-builder upload occurred. That optional route was dropped in favor of
the existing workflow after the maintainer's clarification on 2026-10-29.
No GitHub push, release tag, CRAN submission or public check upload was made.
Update these notes with actual workflow results before tagging or submission.
