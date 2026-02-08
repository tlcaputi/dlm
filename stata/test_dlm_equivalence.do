/*
    test_dlm_equivalence.do
    =======================
    Verifies that DLM beta coefficients exactly match event-study (binned TWFE)
    coefficients to machine precision. This is the mathematical equivalence
    proven in Schmidheiny & Siegloch (2023, JAE).

    The test:
    1. Generates staggered treatment panel data
    2. Runs the DLM estimator
    3. Runs the equivalent binned event-study regression
    4. Compares coefficients — must match to < 1e-6

    Required: reghdfe, dlm.ado, dlm_gen_data.ado in adopath
*/

clear all
set more off

// Add current directory to adopath so Stata finds dlm.ado and dlm_gen_data.ado
adopath ++ "."

// ============================================================================
// Configuration
// ============================================================================
local from_rt = -3
local to_rt   = 3
local ref_period = -1
local tol = 1e-6
local seed = 12345

// ============================================================================
// Step 1: Generate test data
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 1: Generating test data"
display as text "{hline 72}"

dlm_gen_data, n_groups(676) n_times(20) treat_prob(0.4) seed(`seed') clear

summarize
display as text "Data generated: " _N " observations"

// ============================================================================
// Step 2: Run DLM
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 2: Running DLM estimator"
display as text "{hline 72}"

dlm outcome, exposure(post) unit(unit) time(time) ///
    from(`from_rt') to(`to_rt') ref(`ref_period') verbose

// Store DLM betas
matrix dlm_betas = e(betas)
scalar dlm_N = e(N)

display as text ""
display as text "DLM betas:"
matrix list dlm_betas, format(%12.8f)

// ============================================================================
// Step 3: Run equivalent binned event study
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 3: Running equivalent event-study regression"
display as text "{hline 72}"

// Time-window restriction to match DLM's implicit sample
// DLM drops obs where leads/lags are missing at panel edges
// Effective window: [min_time + to_rt, max_time - (abs(from_rt) - 1)]
summarize time, meanonly
local min_time = r(min)
local max_time = r(max)
local es_time_lo = `min_time' + `to_rt'
local es_time_hi = `max_time' - (abs(`from_rt') - 1)
display as text "  Time restriction: [`es_time_lo', `es_time_hi']"

// Create binned event-time dummies
// For from=-3, to=3, ref=-1:
//   es_m3: years_to_treatment <= -3 (binned endpoint, excluding never-treated)
//   es_m2: years_to_treatment == -2
//   es_0:  years_to_treatment == 0
//   es_1:  years_to_treatment == 1
//   es_2:  years_to_treatment == 2
//   es_3:  years_to_treatment >= 3 (binned endpoint, excluding never-treated)
// ref = -1 is omitted

// Dynamically create ES dummies for any from/to/ref
local es_vars ""

forvalues t = `from_rt'/`to_rt' {
    if (`t' == `ref_period') continue

    // Variable name: es_m2 for t=-2, es_0 for t=0, es_3 for t=3, etc.
    if (`t' < 0) {
        local absval = abs(`t')
        local vname "es_m`absval'"
    }
    else {
        local vname "es_`t'"
    }

    // Binned endpoints
    if (`t' == `from_rt') {
        gen byte `vname' = (years_to_treatment <= `from_rt') & ///
                           (years_to_treatment != -1000)
    }
    else if (`t' == `to_rt') {
        gen byte `vname' = (years_to_treatment >= `to_rt') & ///
                           (years_to_treatment != -1000)
    }
    else {
        gen byte `vname' = (years_to_treatment == `t')
    }

    local es_vars "`es_vars' `vname'"
}

display as text "  ES variables: `es_vars'"

// Run reghdfe on the time-windowed sample
reghdfe outcome `es_vars' if (time >= `es_time_lo') & (time <= `es_time_hi'), ///
    absorb(unit time) vce(cluster unit)

scalar es_N = e(N)

// Extract ES coefficients in order
local n_es_vars : word count `es_vars'
display as text ""
display as text "Event-study coefficients:"
matrix es_betas = J(`n_es_vars', 2, .)
local row = 0
forvalues t = `from_rt'/`to_rt' {
    if (`t' == `ref_period') continue
    local ++row

    if (`t' < 0) {
        local absval = abs(`t')
        local vname "es_m`absval'"
    }
    else {
        local vname "es_`t'"
    }

    matrix es_betas[`row', 1] = `t'
    matrix es_betas[`row', 2] = _b[`vname']
    display as text "  t=" %3.0f `t' "  ES beta = " %12.8f _b[`vname']
}

// ============================================================================
// Step 4: Compare DLM vs ES coefficients
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 4: Comparing DLM vs Event-Study coefficients"
display as text "{hline 72}"

// Verify sample sizes match
display as text "  DLM N = " dlm_N "  |  ES N = " es_N
if (dlm_N != es_N) {
    display as error "FAIL: Sample sizes differ (DLM=" dlm_N ", ES=" es_N ")"
    display as error "This suggests the time-windowing is not aligned."
    error 9
}

// Compare coefficients
local n_tests = 0
local n_pass = 0
local max_diff = 0

local es_row = 0
local nrows_dlm = rowsof(dlm_betas)

forvalues i = 1/`nrows_dlm' {
    local t = dlm_betas[`i', 1]
    local dlm_coef = dlm_betas[`i', 2]

    // Skip reference period
    if (`t' == `ref_period') continue

    local ++es_row
    local es_coef = es_betas[`es_row', 2]
    local diff = abs(`dlm_coef' - `es_coef')

    local ++n_tests
    if (`diff' < `tol') {
        local ++n_pass
        local status "PASS"
    }
    else {
        local status "FAIL"
    }

    if (`diff' > `max_diff') {
        local max_diff = `diff'
    }

    display as text "  t=" %3.0f `t' ///
        "  DLM=" %14.10f `dlm_coef' ///
        "  ES=" %14.10f `es_coef' ///
        "  diff=" %12.2e `diff' ///
        "  [`status']"
}

// ============================================================================
// Final verdict
// ============================================================================
display as text ""
display as text "{hline 72}"
if (`n_pass' == `n_tests') {
    display as result "ALL `n_tests' TESTS PASSED (max diff = " %10.2e `max_diff' ")"
    display as result "DLM coefficients match event-study to < `tol'"
}
else {
    local n_fail = `n_tests' - `n_pass'
    display as error "`n_fail' of `n_tests' TESTS FAILED (max diff = " %10.2e `max_diff' ")"
    display as error "DLM coefficients DO NOT match event-study"
    error 9
}
display as text "{hline 72}"

// ============================================================================
// Additional tests: different from/to/ref combinations
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "ADDITIONAL TEST: from=-2, to=4, ref=-1"
display as text "{hline 72}"

// Regenerate data (same seed)
dlm_gen_data, n_groups(676) n_times(20) treat_prob(0.4) seed(`seed') clear

// DLM
dlm outcome, exposure(post) unit(unit) time(time) ///
    from(-2) to(4) ref(-1)
matrix dlm_betas2 = e(betas)
scalar dlm_N2 = e(N)

// ES
summarize time, meanonly
local min_time = r(min)
local max_time = r(max)
local es_time_lo = `min_time' + 4
local es_time_hi = `max_time' - 1

// Create ES dummies for from=-2, to=4, ref=-1
gen byte es2_m2 = (years_to_treatment <= -2) & (years_to_treatment != -1000)
gen byte es2_0  = (years_to_treatment == 0)
gen byte es2_1  = (years_to_treatment == 1)
gen byte es2_2  = (years_to_treatment == 2)
gen byte es2_3  = (years_to_treatment == 3)
gen byte es2_4  = (years_to_treatment >= 4) & (years_to_treatment != -1000)

local es2_vars "es2_m2 es2_0 es2_1 es2_2 es2_3 es2_4"

reghdfe outcome `es2_vars' if (time >= `es_time_lo') & (time <= `es_time_hi'), ///
    absorb(unit time) vce(cluster unit)

scalar es_N2 = e(N)

// Compare
display as text ""
display as text "  DLM N = " dlm_N2 "  |  ES N = " es_N2

local n_tests2 = 0
local n_pass2 = 0
local max_diff2 = 0

local es_row = 0
local nrows_dlm2 = rowsof(dlm_betas2)

forvalues i = 1/`nrows_dlm2' {
    local t = dlm_betas2[`i', 1]
    local dlm_coef = dlm_betas2[`i', 2]

    if (`t' == -1) continue

    local ++es_row

    // Map es_row to the correct variable
    if (`es_row' == 1) local vname "es2_m2"
    if (`es_row' == 2) local vname "es2_0"
    if (`es_row' == 3) local vname "es2_1"
    if (`es_row' == 4) local vname "es2_2"
    if (`es_row' == 5) local vname "es2_3"
    if (`es_row' == 6) local vname "es2_4"

    local es_coef = _b[`vname']
    local diff = abs(`dlm_coef' - `es_coef')

    local ++n_tests2
    if (`diff' < `tol') {
        local ++n_pass2
        local status "PASS"
    }
    else {
        local status "FAIL"
    }

    if (`diff' > `max_diff2') {
        local max_diff2 = `diff'
    }

    display as text "  t=" %3.0f `t' ///
        "  DLM=" %14.10f `dlm_coef' ///
        "  ES=" %14.10f `es_coef' ///
        "  diff=" %12.2e `diff' ///
        "  [`status']"
}

display as text ""
if (`n_pass2' == `n_tests2') {
    display as result "ALL `n_tests2' TESTS PASSED (max diff = " %10.2e `max_diff2' ")"
}
else {
    local n_fail2 = `n_tests2' - `n_pass2'
    display as error "`n_fail2' of `n_tests2' TESTS FAILED (max diff = " %10.2e `max_diff2' ")"
    error 9
}
display as text "{hline 72}"

display as text ""
display as text "All equivalence tests completed successfully."
