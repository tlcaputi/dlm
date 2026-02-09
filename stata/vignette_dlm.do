/*
    vignette_dlm.do
    ================
    Distributed Lag Models in Stata — Tutorial

    This vignette demonstrates the dlm package for Stata.
    It covers installation, basic usage, interpretation,
    and comparison with standard event-study regressions.

    Prerequisites:
        net install dlm, from("https://raw.githubusercontent.com/tlcaputi/dlm-stata/main/")
        ssc install reghdfe
*/

clear all
set more off

// ============================================================================
// 1. Generate Example Data
// ============================================================================
// dlm_gen_data creates a balanced panel with staggered treatment.
// Units are randomly assigned to treatment with probability treat_prob.
// Treated units start treatment at time 7, 8, or 9 (uniformly).
// The true treatment effect is -3 (constant, immediate, persistent).

dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear

summarize outcome post years_to_treatment

// ============================================================================
// 2. Basic DLM Estimation
// ============================================================================
// Estimate the DLM with an event window from -3 to 3, reference period -1.
// This absorbs unit and time fixed effects and clusters SEs at the unit level.

dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)

// The output shows beta coefficients at each event time.
// Pre-treatment betas (t=-3, t=-2) should be near zero.
// Post-treatment betas (t=0, t=1, t=2, t=3) should be near -3.

// ============================================================================
// 3. Accessing Results
// ============================================================================
// All results are stored in e():

display "N observations:  " e(N)
display "N clusters:      " e(N_clust)
display "Event window:    [" e(from) ", " e(to) "]"
display "Reference period: " e(ref_period)

// The key output is the betas matrix:
matrix list e(betas), format(%10.4f)

// Raw gamma coefficients (from the lead/lag regression):
matrix list e(gamma), format(%10.4f)

// ============================================================================
// 4. Different Event Windows
// ============================================================================
// Wider window: from=-5 to 5
dlm outcome, exposure(post) unit(unit) time(time) from(-5) to(5)

// Asymmetric window: more post-treatment periods
dlm outcome, exposure(post) unit(unit) time(time) from(-2) to(6)

// ============================================================================
// 5. Custom Reference Period
// ============================================================================
// Use ref(-2) instead of the default ref(-1):
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) ref(-2)

// ============================================================================
// 6. Adding Covariates
// ============================================================================
// Generate a covariate
gen double x1 = rnormal()

dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) covariates(x1)

drop x1

// ============================================================================
// 7. Subsetting with if
// ============================================================================
// Run only on the first 250 units:
dlm outcome if unit <= 250, exposure(post) unit(unit) time(time) from(-3) to(3)

// ============================================================================
// 8. Verify: DLM = Event Study
// ============================================================================
// The DLM produces betas that are mathematically identical to a binned
// event-study regression. Let's verify.

// 8a. Run the DLM
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) ref(-1)
matrix dlm_b = e(betas)

// 8b. Run the equivalent event study
// Time window: restrict to periods where all leads/lags are observed
summarize time, meanonly
local tlo = r(min) + 3   // 4
local thi = r(max) - 2   // 18

// Create binned event-time dummies
gen byte es_m3 = (years_to_treatment <= -3) & (years_to_treatment != -1000)
gen byte es_m2 = (years_to_treatment == -2)
// ref = -1 is omitted
gen byte es_0  = (years_to_treatment == 0)
gen byte es_1  = (years_to_treatment == 1)
gen byte es_2  = (years_to_treatment == 2)
gen byte es_3  = (years_to_treatment >= 3) & (years_to_treatment != -1000)

reghdfe outcome es_m3 es_m2 es_0 es_1 es_2 es_3 ///
    if (time >= `tlo') & (time <= `thi'), ///
    absorb(unit time) vce(cluster unit)

// 8c. Compare
display ""
display "DLM vs Event Study comparison:"
display "  Period  DLM beta     ES beta      Difference"
display "  ------  ----------   ----------   ----------"

local periods "-3 -2 0 1 2 3"
local es_vars "es_m3 es_m2 es_0 es_1 es_2 es_3"
local dlm_rows "1 2 4 5 6 7"

forvalues j = 1/6 {
    local p : word `j' of `periods'
    local v : word `j' of `es_vars'
    local r : word `j' of `dlm_rows'
    local dlm_coef = dlm_b[`r', 2]
    local es_coef = _b[`v']
    local diff = abs(`dlm_coef' - `es_coef')
    display "  " %5.0f `p' "  " %10.6f `dlm_coef' "   " %10.6f `es_coef' "   " %10.2e `diff'
}

drop es_*

// ============================================================================
// 9. Plotting Results
// ============================================================================
// Use dlm_plot after estimation for publication-quality event study plot
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
dlm_plot, title("Event Study (DLM)") xtitle("Time to Treatment") ytitle("Coefficient")

// ============================================================================
// 10. Looping Over Multiple Outcomes
// ============================================================================
// Create a second outcome variable
gen double outcome2 = outcome + rnormal(0, 2)

foreach var in outcome outcome2 {
    dlm `var', exposure(post) unit(unit) time(time) from(-3) to(3)
    matrix betas_`var' = e(betas)
    display "Betas for `var':"
    matrix list betas_`var', format(%10.4f)
}

drop outcome2

display ""
display "Vignette complete."
