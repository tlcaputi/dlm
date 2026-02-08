# Stata Guide

## Basic Usage

```stata
dlm depvar, exposure(varname) unit(varname) time(varname) from(#) to(#)
```

The `dlm` command estimates a distributed lag model on `depvar`, using `exposure` as the treatment indicator. Unit and time fixed effects are absorbed via `reghdfe`, and standard errors are clustered at the unit level.

## Options

| Option | Required | Description |
|---|---|---|
| `exposure(varname)` | Yes | Treatment/exposure variable (typically 0/1) |
| `unit(varname)` | Yes | Panel unit identifier |
| `time(varname)` | Yes | Time period variable (numeric, evenly spaced) |
| `from(#)` | Yes | Earliest relative period (negative integer) |
| `to(#)` | Yes | Latest relative period (positive integer) |
| `ref(#)` | No | Reference period (default: −1) |
| `covariates(varlist)` | No | Additional control variables |
| `addl_fes(varlist)` | No | Additional fixed effects beyond unit and time |
| `verbose` | No | Display progress information |

## Examples

### Basic estimation

```stata
dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear

dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
```

### Wider event window

```stata
dlm outcome, exposure(post) unit(unit) time(time) from(-5) to(5)
```

### Custom reference period

```stata
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) ref(-2)
```

### With covariates

If your dataset includes control variables, pass them with `covariates()`:

```stata
* Generate data and add a covariate
dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
gen double x1 = rnormal()

dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) covariates(x1)
```

### Additional fixed effects

To absorb fixed effects beyond the default unit and time FEs (e.g., region-by-year):

```stata
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) addl_fes(region)
```

### Subsetting with if/in

```stata
dlm outcome if unit <= 250, exposure(post) unit(unit) time(time) from(-3) to(3)
```

### Multiple outcomes

Stata convention is one outcome per call. Loop for multiple:

```stata
foreach var of varlist outcome1 outcome2 outcome3 {
    dlm `var', exposure(post) unit(unit) time(time) from(-3) to(3)
    matrix betas_`var' = e(betas)
}
```

## Working with Results

### Extract the betas matrix

```stata
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
matrix b = e(betas)

* Access individual cells: b[row, col]
* Columns: 1=time_to_event, 2=coef, 3=se, 4=ci_lo, 5=ci_hi
display "Effect at t=0: " b[4, 2]  // 4th row = period 0 (for from=-3)
```

### Export to CSV

```stata
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
matrix b = e(betas)

* Convert to dataset and export
preserve
clear
local nrows = rowsof(b)
set obs `nrows'
gen time_to_event = .
gen coef = .
gen se = .
gen ci_lo = .
gen ci_hi = .
forvalues i = 1/`nrows' {
    replace time_to_event = b[`i', 1] in `i'
    replace coef = b[`i', 2] in `i'
    replace se = b[`i', 3] in `i'
    replace ci_lo = b[`i', 4] in `i'
    replace ci_hi = b[`i', 5] in `i'
}
export delimited using "dlm_results.csv", replace
restore
```

### Plot results

```stata
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
matrix b = e(betas)

* Convert betas to variables for plotting
preserve
clear
local nrows = rowsof(b)
set obs `nrows'
gen time_to_event = .
gen coef = .
gen ci_lo = .
gen ci_hi = .
forvalues i = 1/`nrows' {
    replace time_to_event = b[`i', 1] in `i'
    replace coef = b[`i', 2] in `i'
    replace ci_lo = b[`i', 4] in `i'
    replace ci_hi = b[`i', 5] in `i'
}

twoway (rcap ci_lo ci_hi time_to_event, lcolor(navy)) ///
       (scatter coef time_to_event, mcolor(navy) msymbol(circle)), ///
       yline(0, lpattern(dash) lcolor(gray)) ///
       xline(-0.5, lpattern(dash) lcolor(gray)) ///
       xtitle("Time to Treatment") ytitle("Coefficient") ///
       title("Event Study (DLM)") legend(off)
restore
```

## Data Preservation

`dlm` uses `preserve`/`restore` internally. Your original dataset is **never modified** — no temporary variables leak into your data, and the observation count remains unchanged.

## Test Data Generator

`dlm_gen_data` creates balanced panel data for testing:

```stata
dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
```

| Option | Default | Description |
|---|---|---|
| `n_groups(#)` | 676 | Number of panel units |
| `n_times(#)` | 20 | Number of time periods |
| `treat_prob(#)` | 0.4 | Fraction of units assigned to treatment |
| `seed(#)` | 12345 | Random seed for reproducibility |
| `clear` | — | Clear existing data in memory |
