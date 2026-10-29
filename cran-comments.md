## Submission

This is an update from ALDEx3 1.4.0, the current CRAN version. Main changes:
a new default inference method for the `"blmm"` mixed-model engine, linear
contrasts (`contrasts` argument), grouped multiple-testing adjustment
(`aldex.p.adjust()`), and faster sparse mixed-model computation. Details are
in NEWS.md.

## Test environments

- Local: Arch Linux x86_64, R 4.6.1, GCC 16.2.1.
- TODO (pending): Windows, macOS and Linux R-devel via the R-hub GitHub
  Actions workflow. Results to be added here before submission.

## R CMD check results

Local `R CMD check --as-cran`:

0 errors | 0 warnings | 2 notes

Both NOTEs are specific to the local environment:

- Non-portable compilation flags, which come from the Arch Linux system R
  configuration, not from the package's `Makevars`.
- HTML validation of the manual was skipped because HTML Tidy and V8 are not
  available on this machine.

## Reverse dependencies

There are no reverse dependencies on CRAN.
