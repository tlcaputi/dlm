# R API Reference

## `distributed_lags_model()`

Estimate a single-outcome distributed lag model.

### Usage

```r
distributed_lags_model(
  data,
  exposure_data,
  from_rt,
  to_rt,
  outcome,
  exposure,
  unit,
  time,
  covariates = NULL,
  addl_fes = NULL,
  ref_period = -1,
  weights = NULL,
  dd = FALSE,
  n = 2,
  dict = NULL
)
```

### Arguments

| Argument | Description |
|---|---|
| `data` | Data frame with unit, time, outcome, and covariates |
| `exposure_data` | Data frame with unit, time, and exposure variable |
| `from_rt` | Starting relative period (negative integer) |
| `to_rt` | Ending relative period (positive integer) |
| `outcome` | String: name of the outcome variable |
| `exposure` | String: name of the exposure variable |
| `unit` | String: name of the unit identifier |
| `time` | String: name of the time variable |
| `covariates` | Character vector of covariate names (default: `NULL`) |
| `addl_fes` | Character vector of additional FE variables (default: `NULL`) |
| `ref_period` | Reference period (default: −1) |
| `weights` | String: name of weight variable (default: `NULL`) |
| `dd` | Logical: include DD companion estimate in plot (default: `FALSE`) |
| `n` | Integer: number of digits for DD companion (default: 2) |
| `dict` | Named vector for renaming axes in the plot (default: `NULL`) |

### Return Value

A list with:

| Element | Type | Description |
|---|---|---|
| `betas` | data.frame | Columns: `time_to_event`, `coef`, `se` |
| `plot` | ggplot | Event-study plot with confidence intervals |
| `model` | fixest | The estimated `feols` model object |
| `vcov` | matrix | Variance-covariance matrix of gamma coefficients |
| `data_periods_included` | vector | Calendar periods included in estimation |
| `fmla_str` | string | Formula string used in `feols` |
| `from_rt` | numeric | From period (echo) |
| `to_rt` | numeric | To period (echo) |
| `exposure` | string | Exposure variable name (echo) |
| `outcome` | string | Outcome variable name (echo) |

---

## `distributed_lags_models()`

Estimate DLMs for multiple outcomes in a single call.

### Usage

```r
distributed_lags_models(
  data, exposure_data, from_rt, to_rt,
  outcomes, exposure, unit, time, ...
)
```

Takes the same arguments as `distributed_lags_model()`, except `outcomes` (plural) is a character vector. Returns a list of model results, one per outcome.

---

## `generate_data()`

Generate a balanced panel with staggered treatment for testing.

### Usage

```r
generate_data(seed = 1234, n_groups = 676, n_times = 20, treat_prob = 0.4)
```

### Return Value

A data frame with columns: `group`, `time`, `treat`, `treatment_time`, `years_to_treatment`, `post`, `outcome`.

---

## `standard_twfe_for_comparison()`

Estimate a standard binned-endpoint TWFE event study for comparison with DLM results.

### Usage

```r
standard_twfe_for_comparison(
  data, from_rt, to_rt, outcome, time, unit,
  time_to_treatment, treat, covariates = NULL,
  ref_period = -1, weights = NULL
)
```

---

## `aggr_es()`

Aggregate event-study treatment effects.

### Usage

```r
aggr_es(mod, period = "post", agg = "mean")
```

| Argument | Description |
|---|---|
| `mod` | A DLM model object |
| `period` | `"pre"` or `"post"` |
| `agg` | `"mean"` or `"cumulative"` |

---

## Helper Functions

| Function | Description |
|---|---|
| `revcumsum(x)` | Reverse cumulative sum of a vector |
| `secumsum(cov)` | SE of forward cumulative sum from a VCV matrix |
| `serevcumsum(cov)` | SE of reverse cumulative sum from a VCV matrix |
| `iplot_data(mod)` | Extract plot data from a fixest model |
| `ggiplot(mod, ...)` | ggplot2-based event-study plot from a fixest model |
