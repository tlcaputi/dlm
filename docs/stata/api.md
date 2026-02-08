# Stata API Reference

## `dlm`

Estimate a distributed lag model.

### Syntax

```
dlm depvar [if] [in], exposure(varname) unit(varname) time(varname)
    from(#) to(#) [ref(#) covariates(varlist) addl_fes(varlist) verbose]
```

### Required Options

| Option | Description |
|---|---|
| `exposure(varname)` | Treatment/exposure variable (numeric, typically 0/1) |
| `unit(varname)` | Panel unit identifier |
| `time(varname)` | Time period variable (numeric, evenly spaced) |
| `from(#)` | Earliest relative period to estimate (negative integer, e.g. −3) |
| `to(#)` | Latest relative period to estimate (positive integer, e.g. 3) |

### Optional

| Option | Default | Description |
|---|---|---|
| `ref(#)` | −1 | Reference (omitted) period |
| `covariates(varlist)` | — | Control variables included in the regression |
| `addl_fes(varlist)` | — | Additional fixed effects passed to `reghdfe absorb()` |
| `verbose` | off | Display progress information |

### Stored Results

#### Scalars

| Name | Description |
|---|---|
| `e(N)` | Number of observations |
| `e(N_clust)` | Number of clusters |
| `e(from)` | Starting relative period |
| `e(to)` | Ending relative period |
| `e(ref_period)` | Reference period |

#### Macros

| Name | Description |
|---|---|
| `e(cmd)` | `"dlm"` |
| `e(outcome)` | Outcome variable name |
| `e(exposure)` | Exposure variable name |
| `e(unit)` | Unit variable name |
| `e(time)` | Time variable name |

#### Matrices

| Name | Dimensions | Description |
|---|---|---|
| `e(betas)` | (to − from + 1) × 5 | Columns: time_to_event, coef, se, ci_lo, ci_hi. Includes reference period row (coef = 0, se = 0). |
| `e(gamma)` | 1 × (|from| + to) | Raw DLM gamma coefficients from the lead/lag regression |
| `e(gamma_V)` | (|from| + to) × (|from| + to) | Variance-covariance matrix of gamma coefficients |

---

## `dlm_gen_data`

Generate balanced panel test data with staggered treatment.

### Syntax

```
dlm_gen_data [, n_groups(#) n_times(#) treat_prob(#) seed(#) clear]
```

### Options

| Option | Default | Description |
|---|---|---|
| `n_groups(#)` | 676 | Number of panel units |
| `n_times(#)` | 20 | Number of time periods |
| `treat_prob(#)` | 0.4 | Fraction of units randomly assigned to treatment (0 < p < 1) |
| `seed(#)` | 12345 | Random number seed |
| `clear` | — | Clear data in memory before generating |

### Variables Created

| Variable | Type | Description |
|---|---|---|
| `unit` | long | Unit identifier (1 to `n_groups`) |
| `time` | long | Time period (1 to `n_times`) |
| `treat` | byte | 1 if unit ever treated, 0 otherwise |
| `treatment_time` | long | Period when treatment begins (7, 8, or 9; missing if never treated) |
| `years_to_treatment` | long | Event time (−1000 for never-treated units) |
| `post` | byte | 1 if post-treatment, 0 otherwise |
| `outcome` | double | N(0, 5) noise + (−3) × post |

### Data Generating Process

Treatment timing is drawn uniformly from {7, 8, 9}. The outcome equals:

$$y_{it} = \varepsilon_{it} - 3 \cdot \mathbf{1}[t \geq t^*_i]$$

where $\varepsilon_{it} \sim N(0, 25)$ and $t^*_i$ is unit $i$'s treatment time.
