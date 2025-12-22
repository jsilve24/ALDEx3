# ALDEx3

## Overview

ALDEx3 is the successor to ALDEx2. ALDEx3 enables both fixed and mixed effects modeling of relative and absolute abundances from sequence count data (e.g., RNA-seq or 16S rRNA-seq data). ALDEx3 additionally offers improved computational efficiency over ALDEx2. Like ALDEx2, ALDEx3 accounts for zeros by modeling sequence count data with a multinomial-Dirichlet model. This model treats zeros as low abundance observations, and accounts for counting uncertainty (i.e., the fact there is more uncertainty in the true proportion of each taxon in each sample when a count of $1$ is observed vs. a count of $1000$). ALDEx3 can model either relative or absolute abundances.

For modeling relative abundances, ALDEx3 uses the Centered Log Ration (CLR) transformation because it preserves distances when mapping compositions (i.e., proportions) to the Euclidean space. Absolute abundances are unmeasured by sequence count data because the _scale_ is unmeasured. The term _scale_ refers, for instance, to the total microbial load (or total cellular transcription, etc.) from the biological system from which a sequencing sample was collected. The scale depends on the experiment and scientific question of interest. For example, in an oral microbiome study, the scale might be the total microbial load per $\mu$L of saliva, or in RNA-seq the total mRNA produced by a cell.

ALDEx3 enables modeling of absolute abundances using _scale models_. A _scale model_ is a probability distribution over the unmeasured scale. The advantage of using a probability distribution rather than typical normalization or bias-correction approaches is that it can explicitly account for uncertainty in the unmeasured scale. Scale models can be defined in a variety of ways such as distributions representing uncertainty in normalizations (e.g., normalizing to housekeeping genes, TSS normalization, etc.), distributions based on secondary measurements of scale (e.g., flow cytometry or qPCR measurements), or just weak biological knowledge from previous studies (e.g., assuming an antibiotic kills between 20-80\% of all microbes). See Defining Scale Models for full details.

## Installation

``` r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("jsilve24/ALDEx3")
```

## Quickstart (Inference with Relative Abundances)

``` r
library(ALDEx3)
# Load Data
data(gut_crohns_data)
# Genus Counts x Samples
Y <- gut_crohns_data$counts
# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))
# Fit Crohns Disease (CD) vs. Control
# relative abundance changes in gut
aldex.gut.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=clr.sm,   # CLR transform
                       gamma=0)        # Gamma=0 no uncertainty in CLR
aldex.gut.summary <- summary.aldex(aldex.gut.raw)
head(aldex.gut.summary)
```

## Mixed Effects Models (SR-MEM)

Here is an example data analysis using the published [SR-MEM method](https://www.biorxiv.org/content/10.1101/2025.08.05.668734v1) (Scale Reliant Mixed Effects Models). Mixed effects modeling can be peformed in ALDEx3 with either the lme4 or nlme interfaces.

``` r
library(ALDEx3)
# Load Data
data(oral_mouthwash_data)
Y <- oral_mouthwash_data$counts
X <- oral_mouthwash_data$metadata
# Fixed Effects: Interaction between Treatment and Discrete Time Point
# Random Effects: Random intercept for each study participant
aldex.mouthwash.raw <- aldex(Y,
                             ~treat*timec+(1|participant_id),
                             X,
                             method="lme4",   # Required
                             nsample=250, 
                             scale=clr.sm,    # Relative Abundances
                             gamma=0)
aldex.mouthwash.summary <- summary.aldex(aldex.mouthwash.raw)
head(aldex.mouthwash.summary)
```

## Inference with Absolute Abundances: Defining Scale Models 

In the Quickstart we considered changes in relative abundances of taxa in the guts of patients with Crohn's Disease (CD) vs. Control. Here we will infer changes in **absolute abundances** using three different approaches to **scale models**. Everything here can also be applied to mixed effects modeling.

### Weak Prior Biological Knowledge

Consider based on prior studies we believe individuals with Crohns (CD) have, on average, 2 log2-fold lower total microbial load (i.e., scale) on average. However, we are not certain about this belief, and want to account for potential error. Let θ<sup>tot</sup><sub>CD</sub> be the average log2-fold change in gut microbial load in samples CD vs. Control. We believe that θ<sup>tot</sup><sub>CD</sub>=-2. However to account for uncertainty in the assumption we will define our scale model using a standard deviation of $0.5$ around this assumption: θ<sup>tot</sup><sub>CD</sub> ~ *N* (-2, 0.5).
We can do this using the `coef.sm` scale model and the parameters `c.mu` and `c.cor`:

``` r
# Load Data
data(gut_crohns_data)
# Genus Counts x Samples
Y <- gut_crohns_data$counts
# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))
# Fit Before vs. After teeth brushing
# relative abundance changes
aldex.gut.abs.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=coef.sm,         # coef scale model
                       c.mu=c(0, -2),
                       c.cor=rbind(c(0, 0), c(0, 0.5**2)))

aldex.gut.abs.summary <- summary.aldex(aldex.gut.abs.raw)
head(aldex.gut.abs.summary)
```

### External Scale Measurements (Ideal)

In the Crohn's dataset, we actually have flow cytometry measurements of microbial load from the collected fecal samples. If this data exists this is usually the ideal approach. To include these measurements in ALDEx3, here we use the `sample.sm` scale model which takes `s.mu` and `s.var` (it can also techincally take a correlation matrix `s.cor`) as arguments. As an example of `s.mu`, consider we have three samples with qPCR $\log_2$ microbial load measurements of `s.mu=c(10.6, 9.5, 11.2)`. We might assume a variance of $0.25$ per sample to account for technical noise. For this dataset we can do:

``` r
# Load Data
data(gut_crohns_data)
# Genus Counts x Samples
Y <- gut_crohns_data$counts
# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))
# Log2 Scale Measurements
s.mu <- log2(X$Average.cell.count..per.gram.of.frozen.feces.)
# Account For Potential Technical Error in Measurements
s.var <- rep(0.25, length(s.mu))
# Estimate Absolute Abundance Changes
aldex.gut.abs.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=sample.sm,         # coef scale model
                       s.mu=s.mu,
                       s.var=s.var)
aldex.gut.abs.summary <- summary.aldex(aldex.gut.abs.raw)
head(aldex.gut.abs.summary)
```

## Normalization Uncertainty

In the Quickstart the CLR transform was used to infer changes in absolute abundances. While not the most powerful approach, we can actually just define a scale model to account for uncertainty in how closely CLR transformed abudances match absolute abundances. [Nixon et al. showed for this dataset that accounting for scale uncertainty can substantially reduce false positives](https://link.springer.com/article/10.1186/s13059-025-03609-3). 

``` r
# Load Data
data(gut_crohns_data)
# Genus Counts x Samples
Y <- gut_crohns_data$counts
# Metadata Samples x Metadata Values
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                           levels=c("Control", "CD"))
# Fit Before vs. After teeth brushing; absolute abundances
aldex.gut.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=clr.sm,  # CLR transform
                       gamma=1)       # Gamma=1: variance in CLR assumption 
aldex.gut.summary <- summary.aldex(aldex.gut.raw)
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
