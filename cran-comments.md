# ALDEx3 1.7.0 release candidate

This candidate updates CRAN version 1.4.0. NEWS records contrasts, selectable
mixed-model inference, sparse computation, grouped adjustment and refit/plot
fixes, together with the new bounded adaptive option. Shared-anchor
Satterthwaite remains the maintained BLMM default; adaptive remains experimental
and opt-in. No new df method or df clipping was introduced.

## Local environment and exact archive

- Arch Linux x86_64, R 4.6.1 (2026-06-24), GCC/G++ 16.2.1, C++17.
- Fresh source build, including both vignettes, followed by fresh
  `R CMD check --as-cran` on that archive in an isolated directory/library.
- Standard `-g -O2` temporary Makevars avoid distribution-specific compiler
  flags. The package Makevars specifies only C++17.
- HTML Tidy 5.8.0 and V8 8.2.1 were provided from temporary check locations.
  These are validation tools, not package dependencies.
- Archive: `ALDEx3_1.7.0.tar.gz` (290,957 bytes).
- SHA256: `57bc30be119832b184888859648ceb2660ca69b7762904ff151857bf5ec77143`.
- Archive inspection found no private investigation inputs, diagnostics,
  working libraries or compiled object files.

## Check results

`R CMD check --as-cran`: **OK**, zero errors, warnings and notes.

Tests reported 1,431 passing expectations, zero failures, warnings and skips.
Examples including donttest, both vignette rebuilds, and PDF/HTML manuals passed.
The source archive and corresponding logs/checksums are preserved locally in
`release/`; earlier 1.6.2 artifacts were preserved before replacing current notes.
The ordinary installed user package was not overwritten.

## Numerical validation and default decision

All candidate R function bodies and native source were verified against the
frozen, previously tested numerical implementation. Matched real-input engine
comparisons use identical response arrays and full intended anchor blocks.
The independent null assessment uses explicitly synthetic counts on a retained
covariate design: 600 independent datasets, three variance regimes and 32 draws.
It is a limited calibration screen rather than universal validation.

Adaptive and full lmerTest agreed on the primary and family-level rejection
decisions in every paired dataset. The weak-variance primary rejection count
was 14/200 for both (7%; two-sided 95% interval 3.88–11.47%). Its uncertainty
did not clear the prespecified default-promotion criterion. Consequently the
default is unchanged. This does not establish that the older approximation is
statistically safer. Boundary matching is not described as proof of calibration.
Convergence warnings from two full-reference calls and approximation-rescue
notifications were retained in the diagnostic evidence; they are not failures
or warnings from the package check above.

## Reverse dependencies

The CRAN package database queried during this preparation lists version 1.4.0
and no direct reverse dependencies through Depends, Imports, LinkingTo or
Suggests. No reverse-dependency check was therefore required for listed CRAN
packages. This does not cover unlisted downstream projects.

## Outstanding external checks — not submission-ready yet

- Windows/R-devel: not run. Automatic approval review blocked uploading this
  source archive to Win-builder and its maintainer-email notification. Explicit
  authorization was requested; no upload occurred.
- macOS and atlas/Linux R-devel: not run. The current candidate is not on
  GitHub, and no verified R-hub upload identity is available locally. The live
  manifest confirms the requested `atlas` platform is R-devel; a check against
  older published main would not validate this archive.
- No GitHub push, release tag, CRAN submission or public R-hub upload was made.

Update these notes with actual platform results before tagging or submission.
