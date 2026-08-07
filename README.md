# Contrasting performance of institutionally distinct empirical harvest control rules under identical operating models

This repository contains R code and data for the publication:

> Yagi and Ichinokawa (2026)

*Contrasting performance of institutionally distinct empirical harvest control rules under identical operating models*

## Overview

The R codes implement simulations for:

- age-structured management strategy evaluation (MSE) simulations;
- implementation of five empirical harvest control rules (HCRs);
- genetic-algorithm-based parameter optimization;
- sensitivity analyses;
- generation of figures and tables presented in the manuscript.

------------------------------------------------------------------------

## Requirements

The analyses were conducted in R version 4.4.2 (2024-10-31).

Required packages include:

``` r
GA
tidyverse
ggplot2
doParallel
FLCore
FLBRP
frasyr23
remotes
cat3advice
lemon
patchwork
ggh4x
```

------------------------------------------------------------------------

## Main scripts

### `default_simulations.R`

Runs the main simulations for four stocks:

- Pollack (*Pollachius pollachius*; k = 0.19)
- Thornback ray (*Raja clavata*; k = 0.09)
- European plaice (*Pleuronectes platessa*; k = 0.23)
- Anchovy (*Engraulis encrasicolus*; k = 0.44)

Five HCRs are evaluated:

1.  rfb rule
2.  type-2 rule
3.  chr rule
4.  rfb rule with $\bar{C}$
5.  type-2 rule with $f$

The simulations consist of:

- a 100-year historical period;
- a 30-year management period.

Outputs are used to produce figures 1–8 and S1-S8.

------------------------------------------------------------------------

### `functions.R`

Contains functions used throughout the analyses, including:

- operating model calculations;
- HCR implementations;
- performance indicator calculations;
- utility functions.

------------------------------------------------------------------------

### `all_ga_do.R`

Performs parameter optimization using a genetic algorithm (GA) for 4 stocks and 5 HCRs.

Optimization is conducted for:

- 11 fishing-histories for 4 stocks;

Outputs include optimal and semi-optimal parameter sets used in the manuscript.

------------------------------------------------------------------------

## Figure scripts

### `Fig1.R` – `Fig8.R`

Scripts used to generate Figures 1–8 in the manuscript.

Each script corresponds to a single figure (e.g., `Fig1.R` generates Figure 1).

------------------------------------------------------------------------

## Data files

### `brps_new.rds`

Reference-point data (`Fcrash`) used in this study.

This file was adapted from the GitHub repository accompanying:

> Fischer et al., *An exploration of the ICES approach for categories 4–6*

<https://github.com/shfischer/GA_MSE_cat456>

------------------------------------------------------------------------

### `stocks.csv`

Biological parameters for the study species.

This file was adapted from the GitHub repository accompanying:

> Fischer et al. (2020)

<https://github.com/shfischer/GA_MSE>

------------------------------------------------------------------------

## Sensitivity analyses

The folder:

``` text
sensitivity_analysis/
```

contains scripts used for all sensitivity analyses reported in the manuscript.
implementing sensitivity analyses by varying steepness to 0.6 and 0.9, process error to 0.3 and 0.9, and observation error to 0.4 and 0.6

------------------------------------------------------------------------

## Reproducing the analyses

A recommended execution order is:

1.  `default_simulations.R`
2.  `all_ga_do.R`
3.  scripts in `sensitivity_analysis/`

------------------------------------------------------------------------