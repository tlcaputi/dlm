# R Guide

## Basic Usage

```r
library(dlm)

mod <- distributed_lags_model(
  data = outcome_data,
  exposure_data = treatment_data,
  from_rt = -3, to_rt = 3,
  outcome = "y", exposure = "treated",
  unit = "id", time = "year"
)
```

The R version takes **separate data frames** for outcomes and exposure. This is useful when the exposure is defined at a different level (e.g., pharmacy-level treatment merged onto beneficiary-level outcomes).

## Data Setup

The R DLM expects two data frames:

1. **`data`** — contains unit, time, outcome, and any covariates
2. **`exposure_data`** — contains unit, time, and the exposure variable (one row per unit-time)

```r
library(dlm)
library(dplyr)

# Example: separate outcome and exposure data
outcome_data <- full_data %>% select(id, year, y, x1, x2)
exposure_data <- full_data %>% select(id, year, treated) %>% distinct()
```

If your data is a single data frame, just select the relevant columns:

```r
df <- generate_data(seed = 42)

data <- df %>% select(group, time, outcome)
exposure_data <- df %>% select(group, time, post) %>% distinct()
```

## Examples

### Basic estimation

```r
library(dlm)
library(dplyr)

df <- generate_data(seed = 42, n_groups = 500, n_times = 20, treat_prob = 0.4)

data <- df %>% select(group, time, outcome)
exposure_data <- df %>% select(group, time, post) %>% distinct()

mod <- distributed_lags_model(
  data = data,
  exposure_data = exposure_data,
  from_rt = -3, to_rt = 3,
  outcome = "outcome", exposure = "post",
  unit = "group", time = "time"
)

mod$betas
```

### Custom reference period

```r
mod <- distributed_lags_model(
  data = data,
  exposure_data = exposure_data,
  from_rt = -3, to_rt = 3,
  outcome = "outcome", exposure = "post",
  unit = "group", time = "time",
  ref_period = -2
)
```

### With covariates

```r
data_with_covs <- df %>% select(group, time, outcome, x1, x2)

mod <- distributed_lags_model(
  data = data_with_covs,
  exposure_data = exposure_data,
  from_rt = -3, to_rt = 3,
  outcome = "outcome", exposure = "post",
  unit = "group", time = "time",
  covariates = c("x1", "x2")
)
```

### Additional fixed effects

```r
mod <- distributed_lags_model(
  data = data,
  exposure_data = exposure_data,
  from_rt = -3, to_rt = 3,
  outcome = "outcome", exposure = "post",
  unit = "group", time = "time",
  addl_fes = c("region")
)
```

### Multiple outcomes

Use `distributed_lags_models` (plural) to estimate multiple outcomes:

```r
mods <- distributed_lags_models(
  data = data,
  exposure_data = exposure_data,
  from_rt = -3, to_rt = 3,
  outcomes = c("outcome1", "outcome2"),
  exposure = "post",
  unit = "group", time = "time"
)

# Access each model
mods[[1]]$betas
mods[[2]]$betas
```

### Weighted regression

```r
mod <- distributed_lags_model(
  data = data,
  exposure_data = exposure_data,
  from_rt = -3, to_rt = 3,
  outcome = "outcome", exposure = "post",
  unit = "group", time = "time",
  weights = "pop_weight"
)
```

## Working with Results

The returned object is a list:

```r
names(mod)
# [1] "betas" "plot" "model" "vcov" "data_periods_included"
# [6] "fmla_str" "from_rt" "to_rt" "exposure" "outcome"
```

### Beta coefficients

```r
mod$betas
#   time_to_event       coef        se
# 1            -3  0.1234567 0.4567890
# 2            -2 -0.0987654 0.4321098
# 3             0 -2.8765432 0.4012345
# ...
```

### Built-in event-study plot

```r
mod$plot
```

The plot is a `ggplot2` object. Customize it:

```r
mod$plot +
  ggplot2::labs(title = "My Event Study", x = "Years to Treatment") +
  ggplot2::theme_minimal()
```

### Access the underlying fixest model

```r
summary(mod$model)
fixest::coeftable(mod$model)
```

### Aggregate treatment effects

The DLM's beta coefficients can be averaged directly:

```r
# Mean post-treatment effect
mean(mod$betas$coef[mod$betas$time_to_event >= 0])

# Mean pre-treatment effect (parallel-trends check)
mean(mod$betas$coef[mod$betas$time_to_event < 0])
```

For standard event studies estimated with `fixest::i()`, you can also use `aggr_es()`:

```r
# Only works with fixest models estimated via i(), NOT with DLM models
library(fixest)
est <- feols(y ~ i(event_time, treat, ref = -1) | unit + time, data = df)
aggr_es(est, period = "post")
aggr_es(est, period = "pre")
aggr_es(est, aggregation = "cumulative")
```

## Test Data Generator

```r
df <- generate_data(seed = 42, n_groups = 500, n_times = 20, treat_prob = 0.4)
```

| Argument | Default | Description |
|---|---|---|
| `seed` | 1234 | Random seed |
| `n_groups` | 676 | Number of panel units |
| `n_times` | 20 | Number of time periods |
| `treat_prob` | 0.4 | Treatment probability |

Key variables: `group`, `time`, `treat`, `treatment_time`, `years_to_treatment`, `post`, `outcome`.
