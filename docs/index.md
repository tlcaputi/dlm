# dlm — Distributed Lag Models

Implementations of the distributed lag model (DLM) framework from [Schmidheiny and Siegloch (2023, *Journal of Applied Econometrics*)](https://doi.org/10.1002/jae.2971) for **R** and **Stata**.

DLMs are mathematically equivalent to event-study regressions with binned endpoints — but often faster and more numerically stable for large panels.

## How It Works

Instead of creating event-time dummies and running a standard TWFE regression, the DLM approach:

1. **Regresses** the outcome on leads and lags of the treatment variable → produces **gamma** (γ) coefficients
2. **Transforms** gammas into cumulative **beta** (β) coefficients via cumulative summation
3. **Propagates** standard errors correctly through the variance-covariance matrix

The resulting betas are identical to what you'd get from a binned event-study regression — proven both mathematically and empirically (see [Theory](theory/background.md)).

## Quick Example

=== "Stata"

    ```stata
    * Install
    net install dlm, from("https://raw.githubusercontent.com/tlcaputi/dlm-stata/main/")

    * Estimate
    dlm outcome, exposure(treated) unit(id) time(year) from(-3) to(3)

    * View results
    matrix list e(betas)
    ```

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

## Citation

If you use this package, please cite:

> Schmidheiny, K. and S. Siegloch (2023). "On event studies and distributed-lag models: Equivalence, generalization and practical implications." *Journal of Applied Econometrics*, 38(5): 695–713.
