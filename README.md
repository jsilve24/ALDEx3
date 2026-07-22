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

## Installation

``` r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("jsilve24/ALDEx3")
```

## Quickstart (Inference with Relative Abundances)

The first example uses `clr.sm` with `gamma=0` to reproduce the fixed CLR assumption used by ALDEx2. This is a compatibility example, not a general recommendation. When no external scale measurements or biological prior information are available, `tss.sm` with `gamma=0.5` is a reasonable starting point. Results should still be checked across scientifically plausible scale models.

``` r
library(ALDEx3)
set.seed(42)
# Load Crohn's Disease vs. Control Data
data(gut_crohns_data)

# NOTE: We recommend a ~75% sparsity filter of counts!
# This is a pragmatic default for speed and stability in very sparse tables;
# adjust based on your depth and scientific question. Always keep an aggregated
# ‘other’ category so each sample’s total count is preserved under the count
Y <- gut_crohns_data$counts
keep_names <- row.names(Y[((rowSums(Y==0))/ncol(Y))<=0.75,])
other <- colSums(Y[((rowSums(Y==0))/ncol(Y))>0.75,])
Y <- Y[keep_names,]
Y <- rbind(Y, other)

# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))

# Relative Abundances CD vs. Control
aldex.gut.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=clr.sm,   # CLR assumption
                       gamma=0)        # Gamma=0 no scale uncertainty
aldex.gut.summary <- summary(aldex.gut.raw)
head(aldex.gut.summary)
```

## Mixed-Effects Models: BLMM, `lme4`, and `nlme`

ALDEx3 supports the published [SR-MEM method](https://link.springer.com/article/10.1186/s40168-026-02377-x) (Scale Reliant Mixed Effects Models) through three engines. The new `method = "blmm"` engine is an ALDEx3-specific approximation that avoids repeated exact variance-component optimization across Monte Carlo draws. It parallelizes across features with `n.cores` and falls back to exact `lme4` fits when the approximation cannot be evaluated cleanly. The exact `lme4` and `nlme` engines remain available for reference analyses.

See the **[ALDEx3 Mixed-Effects Engines vignette](vignettes/ALDEx3-mixed-effects.Rmd)** for model setup, the BLMM formulation, validation guidance, and a runtime and accuracy comparison with exact `lme4` fits.

``` r
library(ALDEx3)
set.seed(42)
# Load Data
data(oral_mouthwash_data)

# 75% Sparsity Filter
Y <- oral_mouthwash_data$counts
keep_names <- row.names(Y[((rowSums(Y==0))/ncol(Y))<=0.75,])
other <- colSums(Y[((rowSums(Y==0))/ncol(Y))>0.75,])
Y <- Y[keep_names,]
Y <- rbind(Y, other)
X <- oral_mouthwash_data$metadata

# Fixed Effects: Interaction between Treatment and Discrete Time Point
# Random Effects: Random intercept for each study participant
aldex.mouthwash.raw <- aldex(Y,
                             ~treat*timec+(1|participant_id),
                             X,
                             method="blmm",   # Fast approximate engine
                             n.cores=4,        # Parallelize across features
                             nsample=250,
                             scale=clr.sm,    # CLR assumption
                             gamma=0)         # Gamma=0 no scale uncertainty
aldex.mouthwash.summary <- summary(aldex.mouthwash.raw)
head(aldex.mouthwash.summary)
```

To run the same model with exact variance-component optimization, replace `method = "blmm"` with `method = "lme4"`.

## Inference with Absolute Abundances: Defining Scale Models 

The first analysis compared relative abundances in Crohn's disease (CD) and control samples. The following examples show how different kinds of scale information change that analysis.

### Prior Knowledge About a Group Difference

Suppose previous studies suggest that total microbial load is about fourfold lower in CD than in controls. A fourfold reduction is -2 on the log2 scale. We are not certain that the reduction is exactly fourfold, so we place an SD of 0.5 around it. Because `Health.status` uses Control as its reference level, the second model coefficient is the CD-minus-Control difference. We therefore set `c.mu = c(0, -2)` and `c.sd = c(0, 0.5)`: the first entries leave the intercept fixed, and the second entries describe the expected CD shift and its uncertainty.

``` r
library(ALDEx3)
set.seed(42)
# Load Data
data(gut_crohns_data)

# NOTE: We recommend a ~75% sparsity filter of counts!
# Genera w/ sparsity > 75% amalgamated into 'other' row
# Keep other counts for multinomial-Dirichlet accuracy
Y <- gut_crohns_data$counts
keep_names <- row.names(Y[((rowSums(Y==0))/ncol(Y))<=0.75,])
other <- colSums(Y[((rowSums(Y==0))/ncol(Y))>0.75,])
Y <- Y[keep_names,]
Y <- rbind(Y, other)

# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))
# Prior for c(intercept, CD-minus-Control)
aldex.gut.abs.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=coefficient.sm,
                       c.mu=c(0, -2),
                       c.sd=c(0, 0.5))

aldex.gut.abs.summary <- summary(aldex.gut.abs.raw)
head(aldex.gut.abs.summary)
```

### External Measurements for Individual Samples

The Crohn's dataset also contains a flow-cytometry estimate of microbial load for every sample. This is more direct information than a group-level prior, so we use `sample.sm`. The mean scale for each sample is the log2 flow-cytometry measurement. We use `s.sd = 0.5` to allow measurement error around each value. These errors are independent across samples; if the measurement errors were correlated, we would supply their covariance matrix through `s.cov` instead.

``` r
library(ALDEx3)
set.seed(42)
# Load Data
data(gut_crohns_data)

# NOTE: We recommend a ~75% sparsity filter of counts!
# Genera w/ sparsity > 75% amalgamated into 'other' row
# Keep other counts for multinomial-Dirichlet accuracy
Y <- gut_crohns_data$counts
keep_names <- row.names(Y[((rowSums(Y==0))/ncol(Y))<=0.75,])
other <- colSums(Y[((rowSums(Y==0))/ncol(Y))>0.75,])
Y <- Y[keep_names,]
Y <- rbind(Y, other)

# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))
# Log2 Scale Measurements
s.mu <- log2(X$Average.cell.count..per.gram.of.frozen.feces.)
# Account For Potential Technical Error in Measurements
s.sd <- rep(0.5, length(s.mu))
# Estimate Absolute Abundance Changes
aldex.gut.abs.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=sample.sm,         # sample scale model
                       s.mu=s.mu,
                       s.sd=s.sd)
aldex.gut.abs.summary <- summary(aldex.gut.abs.raw)
head(aldex.gut.abs.summary)
```

## CLR with Scale Uncertainty

When no external measurement or group-level prior is available, CLR can supply a starting point rather than a fixed answer. The `clr.sm` model centers the analysis on the CLR-implied sample scales. Here, `gamma=1` allows each scale-model coefficient to vary around that center with an SD of one log2 unit. This uncertainty can reduce false positives caused by treating CLR normalization as exact, as demonstrated for this dataset by [Nixon et al.](https://link.springer.com/article/10.1186/s13059-025-03609-3).

``` r
library(ALDEx3)
set.seed(42)
# Load Data
data(gut_crohns_data)

# NOTE: We recommend a ~75% sparsity filter of counts!
# Genera w/ sparsity > 75% amalgamated into 'other' row
# Keep other counts for multinomial-Dirichlet accuracy
Y <- gut_crohns_data$counts
keep_names <- row.names(Y[((rowSums(Y==0))/ncol(Y))<=0.75,])
other <- colSums(Y[((rowSums(Y==0))/ncol(Y))>0.75,])
Y <- Y[keep_names,]
Y <- rbind(Y, other)

# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))
# Allow uncertainty around the CLR-implied scale differences
aldex.gut.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=clr.sm,  # CLR assumption
                       gamma=1)       # Gamma=1: SD around CLR assumption
aldex.gut.summary <- summary(aldex.gut.raw)
head(aldex.gut.summary)
```

## Model Outputs

The output of `summary.aldex` is typically most useful. It looks like:
```
        parameter          entity   estimate std.error    p.val.adj
1 Health.statusCD    Agathobacter -1.7951101  1.702948 3.665291e-01
2 Health.statusCD     Akkermansia -6.1827824  1.168809 9.867942e-06
3 Health.statusCD       Alistipes -6.1123999  1.307723 4.871389e-05
4 Health.statusCD Anaerobutyricum -3.2153272  1.397338 4.276630e-02
5 Health.statusCD    Anaerostipes -6.2694178  1.177097 9.477381e-06
6 Health.statusCD     Bacteroides  0.1164949  1.503918 9.687594e-01
```

The output is a data.frame which includes:
- parameter: The name of the regression parameter for the result, typical of the summary output of `lm`, `lmer`, etc.
- entity: The taxon or gene name for that result (corresponding to `row.names` of input counts `Y`)
- estimate: The average regression coefficient across Monte Carlo draws
- std.err: The average regression standard error across Monte Carlo draws
- p.val.adj: The BH adjusted p-values combined from Monte Carlo draws using special procedure
