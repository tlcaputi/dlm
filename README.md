# dlm

Distributed lag models (DLMs) equivalent to event studies with binned endpoints, following [Schmidheiny and Siegloch (2023, *Journal of Applied Econometrics*)](https://doi.org/10.1002/jae.2971).

A DLM estimates the same parameters as a standard event-study regression with binned leads/lags but uses a different parameterization that can be more convenient for continuous or multi-valued treatments. This package implements the method of Schmidheiny and Siegloch and is based in part on their replication code.

## Installation

```r
# install.packages("devtools")
devtools::install_github("tlcaputi/dlm")
```

## Quick start

```r
library(dlm)

# Generate example panel data (true effect = -3)
df <- generate_data(seed = 1234)

# Run distributed lag model
results <- distributed_lags_models(
  data          = df,
  exposure_data = df,
  from_rt       = -5,
  to_rt         = 5,
  outcomes      = "outcome",
  exposure      = "post",
  unit          = "group",
  time          = "time",
  ref_period    = -1
)

# View coefficients
results[[1]]$betas

# Plot
results[[1]]$plot
```

## Equivalence to event studies

The key result from Schmidheiny and Siegloch (2023) is that DLM coefficients, obtained by cumulatively summing the distributed lag parameters, are numerically identical to event-study coefficients estimated with `fixest::i(..., bin = .)`. This package implements the cumulative summation and corresponding standard error calculations using the full variance-covariance matrix.

## References

Schmidheiny, K., and S. Siegloch. 2023. "On Event Studies and Distributed-Lags in Two-Way Fixed Effects Models: Identification, Equivalence, and Generalization." *Journal of Applied Econometrics* 38(5): 695-713.
