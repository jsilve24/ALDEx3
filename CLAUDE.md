# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

ALDEx3 is an R package for scalable generalized linear and mixed effects models for sequence count data (16S rRNA, RNA-seq). It models counts via a Dirichlet-multinomial sampling framework and introduces "scale models" to explicitly quantify uncertainty about scale rather than implicitly fixing scale via normalization. This generalizes the earlier ALDEx2 package.

Key papers:
- Nixon et al. (2025) Genome Biology — scale uncertainty framework
- McGovern et al. (2025) bioRxiv — SR-MEM method

## Common Commands

```r
# Run all tests
devtools::test()

# Run a single test file
devtools::test(filter = "aldex_lm")   # matches test-aldex_lm.R
devtools::test(filter = "fflm")       # matches test-fflm.R

# Regenerate documentation from roxygen2 comments
devtools::document()

# Full R CMD check (builds, checks examples, runs tests)
devtools::check()

# Install locally
devtools::install()
```

From the shell:
```bash
R CMD build ALDEx3/
R CMD check ALDEx3_*.tar.gz
```

## Architecture

### Core Data Flow

1. **Input**: D×N count matrix `Y` (D taxa, N samples) + formula/model matrix `X`
2. **Sampling** (`aldex-internals.R`): Draw `nsample` Monte Carlo replicates from a Dirichlet-multinomial distribution → 3D log-composition array (N × D × nsample)
3. **Scale models** (`scale.R`): Apply a scale function to each replicate → N × nsample log-scale matrix
4. **Model fitting**: Fit taxon-wise linear models across all Monte Carlo replicates
   - Fixed effects only → `fflm()` (`fflm.R`)
   - Mixed effects → `sr.mem()` (`sr-mem.R`)
5. **Aggregation** (`aldex-internals.R`): Average estimates across replicates, compute p-values
6. **Output**: List with estimates, std errors, p-values; `summary.aldex()` formats as a data frame

### Key Files

| File | Purpose |
|------|---------|
| `R/aldex.R` | Main `aldex()` entry point; routes to lm vs. mixed effects |
| `R/aldex-internals.R` | Dirichlet sampling, streaming, p-value aggregation |
| `R/scale.R` | Scale model functions: `clr.sm()`, `coefficient.sm()`, `sample.sm()`, `tss.sm()` |
| `R/fflm.R` | Fast vectorized fixed-effects linear model over 3D arrays; HC0/HC3 robust SEs |
| `R/sr-mem.R` | Scale-Reliant Mixed Effects Models via `lme4`/`nlme` backends |
| `R/aldex-summary.R` | `summary.aldex()` S3 method and `cohensd()` effect sizes |
| `R/simulation.R` | Simulation helpers used in tests |
| `R/dirichlet.R` | `rLogDirichlet()` sampling |

### Scale Models

Scale functions have a flexible API: they receive any subset of `{X, Y, logComp}` via dynamic parameter injection (`formalArgs()`), and return an N × nsample log-scale matrix.

- `clr.sm(gamma=0.5)` — CLR with scale uncertainty controlled by `gamma`
- `coefficient.sm()` — weak prior knowledge encoded via model coefficients
- `sample.sm()` — external scale measurements per sample
- `tss.sm()` — total-sum scaling

### Memory / Performance

- **Streaming**: When the estimated memory footprint exceeds `streamsize` (default 8000 MB), sampling is chunked; only aggregated estimates are stored, not intermediate arrays.
- **Parallelism**: `sr.mem()` uses `parallel::makeCluster()` (default: `detectCores()-1`) for taxon-wise mixed effects fitting.
- **Speed**: `fflm()` achieves ~1000× speedup over looped `lm()` for HC3/HC0 SEs via vectorized matrix ops on the 3D array.

### Testing

Tests use `testthat` edition 3. Simulation-based tests call `aldex.lm.sim.clr()` to generate synthetic count data with known ground truth.

- `test-aldex_lm.R` — main `aldex()` with fixed effects, streaming
- `test-aldex-mem.R` — mixed effects (lme4, nlme)
- `test-fflm.R` — fast linear model correctness
- `test-scale.R` — scale model outputs
- `test-dirichlet.R` — Dirichlet sampler
- `test-coda.R` — CODA utilities
- `test-aldex-summary.R` — summary output format

### CI/CD

R-hub via `.github/workflows/rhub.yaml`; triggered manually with configurable platforms (linux, windows, macos). Requires `RHUB_TOKEN` secret.
