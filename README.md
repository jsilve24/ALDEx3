# ALDEx3

## Overview

**ALDEx3** is the successor to **ALDEx2**, providing improved computational efficiency and support for both fixed- and mixed-effects models. It enables statistically principled inference on both **relative** and **absolute** abundances from sequencing count data (e.g., RNA-seq or 16S rRNA gene sequencing).

Like ALDEx2, ALDEx3 models sequencing counts using a **Dirichlet–multinomial** framework, which accounts for sampling variability and naturally accommodates zero counts through the multinomial likelihood (e.g., the model does not require treating zeros as true structural absences). 

High-throughput sequencing data tell us how reads are divided among features, but not how much biological material was present in each sample. We call that missing total the sample's **scale**. Depending on the experiment, scale might mean total microbial load, total cellular RNA, or total biomass. A normalization method fills in this missing information by assumption; ALDEx3 instead lets the user state the assumption and its uncertainty explicitly.

ALDEx3 addresses this limitation by explicitly modeling **scale uncertainty**. Rather than fixing scale through normalization, ALDEx3 uses **scale models**—probability distributions that represent uncertainty in the unobserved scale and propagate that uncertainty into downstream inference. This approach generalizes normalization while making its assumptions explicit and testable.

Choose the scale model that matches the information available:

1. Use `sample.sm` when each sample has an external measurement, such as flow cytometry, qPCR, or a spike-in estimate.
2. Use `coefficient.sm` when prior knowledge concerns a group difference or another model term, such as an expected treatment-related reduction in total microbial load.
3. Use `clr.sm` when CLR normalization is a reasonable center but its implied scale differences are uncertain.
4. Use `tss.sm` when total-sum scaling (TSS)—no systematic scale difference between groups—is a reasonable starting assumption.
5. Write a custom scale-model function when the built-in models do not match the study. See `?aldex`, under “Writing a custom scale model,” for the required function interface and output.

ALDEx3 carries this uncertainty through the analysis instead of treating one normalization as exact. This supports inference about relative or absolute changes and works with both fixed- and mixed-effects models.

Before using ALDEx3, users are **strongly encouraged** to read:

> **Nixon, Gloor, and Silverman (2025).**  
> *Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization*.  
> **Genome Biology**.  
> https://link.springer.com/article/10.1186/s13059-025-03609-3

This article provides the theoretical motivation for scale models, explains why normalization implicitly fixes unmeasured scale, and demonstrates how explicitly modeling scale uncertainty improves robustness and reproducibility in differential abundance and expression analyses.

## How it works

For each of `nsample` Monte Carlo draws, `aldex()`:

1. draws each sample's composition from a Dirichlet–multinomial posterior given its counts;
2. draws each sample's scale from the chosen scale model;
3. combines the two into log2 abundances (log composition plus log scale);
4. fits the linear model (`method = "lm"`) or mixed model (`"blmm"`, `"lme4"`, `"nlme"`) to every feature, giving estimates, standard errors and p-values for that draw;
5. adjusts the p-values across features within the draw. Estimates, standard errors and adjusted p-values are then averaged over draws.

## Installation

Install the released version from CRAN:

``` r
install.packages("ALDEx3")
```

Or install the development version from GitHub:

``` r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("jsilve24/ALDEx3")
```

## Preparing the examples

The examples below share one preprocessing step. Very sparse features are pooled into a single `other` row. This keeps each sample's total count, which the Dirichlet–multinomial model uses. The 75% threshold is a pragmatic default for speed and stability; adjust it to your sequencing depth and scientific question.

``` r
library(ALDEx3)

# Pool features that are zero in more than 75% of samples into an 'other' row.
filter_sparse <- function(Y, max.zero = 0.75) {
  sparse <- rowMeans(Y == 0) > max.zero
  rbind(Y[!sparse, , drop = FALSE], other = colSums(Y[sparse, , drop = FALSE]))
}

# Crohn's disease vs. control data
data(gut_crohns_data)
Y <- filter_sparse(gut_crohns_data$counts)
X <- gut_crohns_data$metadata   # samples x metadata values
X$Health.status <- factor(X$Health.status, levels = c("Control", "CD"))
```

## Quickstart (Inference with Relative Abundances)

The first example uses `clr.sm` with `gamma=0` to reproduce the fixed CLR assumption used by ALDEx2. This is a compatibility example, not a general recommendation. When no external scale measurements or biological prior information are available, `tss.sm` with `gamma=0.5` is a reasonable starting point. Results should still be checked across scientifically plausible scale models.

``` r
set.seed(42)
# Relative abundances, CD vs. Control
aldex.gut.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=clr.sm,   # CLR assumption
                       gamma=0)        # gamma=0: no scale uncertainty
head(summary(aldex.gut.raw))
```

## Mixed-Effects Models: BLMM, `lme4`, and `nlme`

ALDEx3 supports the published [SR-MEM method](https://link.springer.com/article/10.1186/s40168-026-02377-x) (Scale Reliant Mixed Effects Models). For a mixed-effects analysis, use `method = "blmm"` with its default adaptive Satterthwaite inference:

``` r
set.seed(42)
data(oral_mouthwash_data)
Y.oral <- filter_sparse(oral_mouthwash_data$counts)
X.oral <- oral_mouthwash_data$metadata

# Fixed effects: interaction between treatment and discrete time point
# Random effects: random intercept for each study participant
aldex.mouthwash.raw <- aldex(Y.oral,
                             ~treat*timec+(1|participant_id),
                             X.oral,
                             method="blmm",   # approximate mixed-effects engine
                             n.cores=4,       # parallelize across features
                             nsample=250,
                             scale=clr.sm,    # CLR assumption
                             gamma=0)         # gamma=0: no scale uncertainty
head(summary(aldex.mouthwash.raw))
```

For each feature, BLMM fits one shared covariance starting point (the **anchor**) across all draws, then refines each draw from there and computes Satterthwaite degrees of freedom per coefficient. Draws that fail its numerical checks are refit with `lmerTest`, and a warning reports how many.

To test linear combinations of fixed effects, pass a named contrast matrix via `contrasts`; to adjust several coefficients or contrasts as one multiple-testing family after fitting, use `aldex.p.adjust()`. See the [mixed-effects vignette](https://CRAN.R-project.org/package=ALDEx3/vignettes/ALDEx3-mixed-effects.html) (`vignette("ALDEx3-mixed-effects")`), `?aldex` and `?aldex.p.adjust` for a worked example, interpretation, practical guidance and the technical details.

### Other inference options and engines

| Selection | Use |
|:--|:--|
| `test = "t.satterthwaite"` | Shared-anchor approximation without per-draw refinement. Cheaper; reproduces the BLMM default of ALDEx3 1.5–1.6. |
| `method = "lme4"` | Fit every draw with `lmerTest` (Satterthwaite inference). Slower, useful as a comparison. |
| `method = "nlme"` | Native `nlme` fitting, including residual correlation structures (`correlation = ...`). |

`test = "t.residual"` is deprecated and warns when used. BLMM `test` values cannot be passed to the other engines.

## Inference with Absolute Abundances: Defining Scale Models 

The first analysis compared relative abundances in Crohn's disease (CD) and control samples. The following examples show how different kinds of scale information change that analysis. They reuse `Y` and `X` from [Preparing the examples](#preparing-the-examples).

### Prior Knowledge About a Group Difference

Suppose previous studies suggest that total microbial load is about fourfold lower in CD than in controls. A fourfold reduction is -2 on the log2 scale. We are not certain that the reduction is exactly fourfold, so we place an SD of 0.5 around it. Because `Health.status` uses Control as its reference level, the second model coefficient is the CD-minus-Control difference. We therefore set `c.mu = c(0, -2)` and `c.sd = c(0, 0.5)`: the first entries leave the intercept fixed, and the second entries describe the expected CD shift and its uncertainty.

``` r
set.seed(42)
# Prior for c(intercept, CD-minus-Control)
aldex.gut.coef <- aldex(Y,
                        ~Health.status,
                        X,
                        nsample=2000,
                        scale=coefficient.sm,
                        c.mu=c(0, -2),
                        c.sd=c(0, 0.5))
head(summary(aldex.gut.coef))
```

### External Measurements for Individual Samples

The Crohn's dataset also contains a flow-cytometry estimate of microbial load for every sample. This is more direct information than a group-level prior, so we use `sample.sm`. The mean scale for each sample is the log2 flow-cytometry measurement. We use `s.sd = 0.5` to allow measurement error around each value. These errors are independent across samples; if the measurement errors were correlated, we would supply their covariance matrix through `s.cov` instead.

``` r
set.seed(42)
# Log2 scale measurements, with SD 0.5 for technical error
s.mu <- log2(X$Average.cell.count..per.gram.of.frozen.feces.)
s.sd <- rep(0.5, length(s.mu))
aldex.gut.flow <- aldex(Y,
                        ~Health.status,
                        X,
                        nsample=2000,
                        scale=sample.sm,   # sample scale model
                        s.mu=s.mu,
                        s.sd=s.sd)
head(summary(aldex.gut.flow))
```

## CLR with Scale Uncertainty

When no external measurement or group-level prior is available, CLR can supply a starting point rather than a fixed answer. The `clr.sm` model centers the analysis on the CLR-implied sample scales. Here, `gamma=1` allows each scale-model coefficient to vary around that center with an SD of one log2 unit. This uncertainty can reduce false positives caused by treating CLR normalization as exact, as demonstrated for this dataset by [Nixon et al.](https://link.springer.com/article/10.1186/s13059-025-03609-3).

``` r
set.seed(42)
# Allow uncertainty around the CLR-implied scale differences
aldex.gut.clr1 <- aldex(Y,
                        ~Health.status,
                        X,
                        nsample=2000,
                        scale=clr.sm,  # CLR assumption
                        gamma=1)       # gamma=1: SD of 1 around CLR assumption
head(summary(aldex.gut.clr1))
```

## Model Outputs

`summary()` of an `aldex` fit is usually the most useful output. For the flow-cytometry (`sample.sm`) example above, `head(summary(aldex.gut.flow))` gives:

```
        parameter          entity   estimate std.error    p.val.adj
1 Health.statusCD    Agathobacter -1.8036338  1.704191 3.585030e-01
2 Health.statusCD     Akkermansia -6.1972896  1.171886 1.341654e-05
3 Health.statusCD       Alistipes -6.1010462  1.306783 4.760054e-05
4 Health.statusCD Anaerobutyricum -3.2366868  1.403536 4.306112e-02
5 Health.statusCD    Anaerostipes -6.2686881  1.171873 8.778170e-06
6 Health.statusCD     Bacteroides  0.1276413  1.499862 9.689928e-01
```

Exact values vary slightly with the random draws. The data.frame has one row per model coefficient and feature:

- `parameter`: the regression coefficient (or contrast) name, as in the output of `lm`, `lmer`, etc.
- `entity`: the taxon or gene name (the `row.names` of the input counts `Y`).
- `estimate`: the regression coefficient, averaged over Monte Carlo draws.
- `std.error`: the within-fit standard error, averaged over Monte Carlo draws.
- `p.val.adj`: the adjusted p-value (BH by default), adjusted within each draw and then averaged over draws.
- `type`: present only when contrasts were requested; marks each row as `"coefficient"` or `"contrast"`.

The per-draw arrays (`estimate`, `std.error`, `df`) and the raw `p.val` are stored in the fitted object; see `?aldex`.
