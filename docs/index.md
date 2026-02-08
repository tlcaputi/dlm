# dlm — Distributed Lag Models

Implementations of the distributed lag model (DLM) framework from [Schmidheiny and Siegloch (2023, *Journal of Applied Econometrics*)](https://doi.org/10.1002/jae.2971) for **R** and **Stata**.

DLMs generalize the canonical event study to settings with **continuous treatments** that can change in magnitude, sign, and timing throughout the study period. When applied to a binary treatment, the DLM produces estimates that are numerically identical to a binned-endpoint event study.

## Why DLMs?

The canonical event study uses event-time indicator variables (dummies), which only make sense when treatment is a one-time binary switch. But many empirical settings involve **continuous treatments** — tax rates that change by different amounts across states, minimum wages that increase and decrease at varying magnitudes over time, policy shocks of different sizes hitting different units at different times. In these settings, event-time dummies don't apply. Researchers are often forced to dichotomize continuous treatments into "big changes" versus no change, which fundamentally alters the analysis by discarding variation and potentially biasing estimates.

The DLM solves this by replacing event-time dummies with **leads and lags of the treatment variable itself**. Because it operates on the treatment variable directly, it naturally handles treatments that are continuous, change sign, vary in magnitude, and occur multiple times per unit.

**Key properties:**

- **Generalizes event studies to continuous treatments.** The DLM can produce event-study-style dynamic treatment effect plots for treatments that increase and decrease at different magnitudes over many time periods — something canonical event studies simply cannot do.
- **Equivalent to binned event studies for binary treatments.** When applied to the special case of a binary absorbing treatment, the DLM produces betas that are numerically identical to a properly specified event study with binned endpoints. The DLM is therefore a strict generalization of the canonical event study (see [Theory](theory/background.md)).
- **More data hungry than canonical event studies.** The DLM requires observing the treatment variable over a wider window than the outcome. Leads and lags at panel edges create missing values that reduce the estimation sample, so wider event windows require substantially more data.

## How It Works

1. **Regresses** the outcome on leads and lags of the treatment variable → produces **gamma** (γ) coefficients measuring incremental treatment effects
2. **Transforms** gammas into cumulative **beta** (β) coefficients via cumulative summation
3. **Propagates** standard errors correctly through the variance-covariance matrix

The resulting betas have the same interpretation as event-study coefficients: dynamic treatment effects at each time horizon relative to the reference period.

## Quick Example

=== "R"

    ```r
    # Install
    devtools::install_github("tlcaputi/dlm")

    # Estimate
    library(dlm)
    mod <- distributed_lags_model(
      data = outcome_data,
      exposure_data = treatment_data,
      from_rt = -3, to_rt = 3,
      outcome = "outcome", exposure = "treated",
      unit = "id", time = "year"
    )

    # View results
    mod$betas
    mod$plot
    ```

=== "Stata"

    ```stata
    * Install
    net install dlm, from("https://raw.githubusercontent.com/tlcaputi/dlm-stata/main/")

    * Estimate
    dlm outcome, exposure(treated) unit(id) time(year) from(-3) to(3)

    * View results
    matrix list e(betas)
    ```

## Features

| Feature | R | Stata |
|---|---|---|
| High-dimensional FE (unit + time) | fixest | reghdfe |
| Clustered SEs | Unit-level | Unit-level |
| Additional covariates | ✓ | ✓ |
| Additional fixed effects | ✓ | ✓ |
| Custom reference period | ✓ | ✓ |
| Event-study plot | Built-in | Via `e(betas)` matrix |
| Multiple outcomes per call | ✓ | Loop externally |
| Weighted regression | ✓ | — |

## GitHub

| | |
|---|---|
| **R package** | [github.com/tlcaputi/dlm](https://github.com/tlcaputi/dlm) |
| **Stata package** | [github.com/tlcaputi/dlm-stata](https://github.com/tlcaputi/dlm-stata) |
| **Developer** | [Theodore Caputi](https://www.theodorecaputi.com) |

## Citation

If you use this package, please cite:

> Schmidheiny, K. and S. Siegloch (2023). "On event studies and distributed-lag models: Equivalence, generalization and practical implications." *Journal of Applied Econometrics*, 38(5): 695–713.
