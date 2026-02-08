<p align="center">
  <img src="docs/assets/logo.svg" alt="dlm logo" width="120">
</p>

<h1 align="center">dlm</h1>

<p align="center">
  <strong>Distributed Lag Models for R</strong><br>
  Generalize event studies to continuous treatments
</p>

<p align="center">
  <a href="https://tlcaputi.github.io/dlm/">Documentation</a> &middot;
  <a href="https://github.com/tlcaputi/dlm-stata">Stata version</a>
</p>

---

R implementation of the distributed lag model (DLM) framework from [Schmidheiny and Siegloch (2023, *Journal of Applied Econometrics*)](https://doi.org/10.1002/jae.2971).

DLMs generalize the canonical event study to settings with **continuous treatments** that can change in magnitude, sign, and timing throughout the study period. When applied to a binary absorbing treatment, the DLM produces estimates that are numerically identical to a binned-endpoint event study.

## Installation

```r
# install.packages("devtools")
devtools::install_github("tlcaputi/dlm")
```

## Quick Start

```r
library(dlm)
library(dplyr)

# Generate test data (true treatment effect = -3)
df <- generate_data(seed = 42, n_groups = 500, n_times = 20, treat_prob = 0.4)

outcome_data <- df %>% select(group, time, outcome)
exposure_data <- df %>% select(group, time, post) %>% distinct()

# Estimate DLM
mod <- distributed_lags_model(
  data = outcome_data,
  exposure_data = exposure_data,
  from_rt = -3, to_rt = 3,
  outcome = "outcome", exposure = "post",
  unit = "group", time = "time"
)

# View results
mod$betas
mod$plot
```

## Why DLMs?

The canonical event study uses event-time dummies, which only work for binary treatments that turn on once and stay on. Many empirical settings involve **continuous treatments** — tax rates, minimum wages, policy dosages — where event-time dummies don't apply. The DLM replaces these dummies with leads and lags of the treatment variable itself, naturally handling treatments that are continuous, change sign, vary in magnitude, and occur multiple times per unit.

See the [documentation](https://tlcaputi.github.io/dlm/) for a full explanation with concrete examples.

## Citation

> Schmidheiny, K. and S. Siegloch (2023). "On event studies and distributed-lag models: Equivalence, generalization and practical implications." *Journal of Applied Econometrics*, 38(5): 695-713.
