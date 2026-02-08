/*
    test_cross_language.do
    ======================
    Cross-language equivalence test: Stata DLM vs R DLM on shared data.

    Strategy:
    1. Generate data in Stata, export as CSV
    2. Run R dlm package on same data via Rscript, export R results as CSV
    3. Run Stata DLM on same data
    4. Compare Stata betas vs R betas — must match to < 1e-6

    This proves the Stata implementation produces identical results to the
    proven R implementation on identical input data.
*/

clear all
set more off

// Use package files from this directory
adopath ++ "."

local tol = 1e-6
local stata_dir = "`c(pwd)'"

// ============================================================================
// Step 1: Generate data in Stata and export as CSV
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 1: Generate shared test data in Stata"
display as text "{hline 72}"

dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear

// Export to CSV for R
export delimited using "`stata_dir'/test_shared_data.csv", replace
display as text "Exported test data to test_shared_data.csv"

// ============================================================================
// Step 2: Run Stata DLM
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 2: Run Stata DLM (from=-3, to=3, ref=-1)"
display as text "{hline 72}"

dlm outcome, exposure(post) unit(unit) time(time) ///
    from(-3) to(3) ref(-1) verbose

matrix stata_betas = e(betas)
display as text ""
display as text "Stata DLM betas:"
matrix list stata_betas, format(%14.10f)

// Save Stata betas to CSV for comparison
local nrows = rowsof(stata_betas)
quietly {
    preserve
    clear
    set obs `nrows'
    gen double time_to_event = .
    gen double coef = .
    gen double se = .
    forvalues i = 1/`nrows' {
        replace time_to_event = stata_betas[`i', 1] in `i'
        replace coef = stata_betas[`i', 2] in `i'
        replace se = stata_betas[`i', 3] in `i'
    }
    export delimited using "`stata_dir'/test_stata_results.csv", replace
    restore
}
display as text "Exported Stata results to test_stata_results.csv"

// ============================================================================
// Step 3: Run R DLM on same data via Rscript
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 3: Run R DLM on same data"
display as text "{hline 72}"

// Write the R script
tempname fh
file open `fh' using "`stata_dir'/test_r_runner.R", write replace
file write `fh' `"library(dlm)"' _n
file write `fh' `"library(dplyr)"' _n
file write `fh' `""' _n
file write `fh' `"# Read Stata-generated data"' _n
file write `fh' `"args <- commandArgs(trailingOnly = TRUE)"' _n
file write `fh' `"data_path <- args[1]"' _n
file write `fh' `"out_path <- args[2]"' _n
file write `fh' `""' _n
file write `fh' `"df <- read.csv(data_path)"' _n
file write `fh' `"cat('R: Read', nrow(df), 'rows\n')"' _n
file write `fh' `""' _n
file write `fh' `"# Prepare data and exposure_data (R DLM expects separate)"' _n
file write `fh' `"data <- df %>% select(unit, time, outcome)"' _n
file write `fh' `"exposure_data <- df %>% select(unit, time, post) %>% distinct()"' _n
file write `fh' `""' _n
file write `fh' `"cat('R: data rows =', nrow(data), '\n')"' _n
file write `fh' `"cat('R: exposure_data rows =', nrow(exposure_data), '\n')"' _n
file write `fh' `""' _n
file write `fh' `"# Run R DLM"' _n
file write `fh' `"mod <- distributed_lags_model("' _n
file write `fh' `"    data = data,"' _n
file write `fh' `"    exposure_data = exposure_data,"' _n
file write `fh' `"    from_rt = -3,"' _n
file write `fh' `"    to_rt = 3,"' _n
file write `fh' `"    outcome = 'outcome',"' _n
file write `fh' `"    exposure = 'post',"' _n
file write `fh' `"    unit = 'unit',"' _n
file write `fh' `"    time = 'time',"' _n
file write `fh' `"    ref_period = -1"' _n
file write `fh' `")"' _n
file write `fh' `""' _n
file write `fh' `"cat('R DLM betas:\n')"' _n
file write `fh' `"print(mod[['betas']])"' _n
file write `fh' `""' _n
file write `fh' `"# Export results"' _n
file write `fh' `"write.csv(mod[['betas']], out_path, row.names = FALSE)"' _n
file write `fh' `"cat('R: Results written to', out_path, '\n')"' _n
file close `fh'

// Run the R script
display as text "Running R script..."
shell Rscript "`stata_dir'/test_r_runner.R" "`stata_dir'/test_shared_data.csv" "`stata_dir'/test_r_results.csv"

// ============================================================================
// Step 4: Compare Stata vs R results
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 4: Compare Stata vs R DLM betas"
display as text "{hline 72}"

// Load R results
preserve
import delimited using "`stata_dir'/test_r_results.csv", clear varnames(1)
display as text "R results loaded: " _N " rows"
list, sep(0)

// Store R results in matrix
local r_nrows = _N
matrix r_betas = J(`r_nrows', 3, .)
forvalues i = 1/`r_nrows' {
    matrix r_betas[`i', 1] = time_to_event[`i']
    matrix r_betas[`i', 2] = coef[`i']
    matrix r_betas[`i', 3] = se[`i']
}
restore

// Compare
display as text ""
display as text %14s "Time" %16s "Stata coef" %16s "R coef" %14s "Coef diff" %16s "Stata SE" %16s "R SE" %14s "SE diff" %8s "Status"
display as text "{hline 120}"

local n_tests = 0
local n_pass = 0
local max_coef_diff = 0
local max_se_diff = 0

// R results don't include the ref period row, so we need to map carefully
local r_row = 0
forvalues i = 1/`nrows' {
    local t = stata_betas[`i', 1]
    local stata_coef = stata_betas[`i', 2]
    local stata_se = stata_betas[`i', 3]

    // Skip ref period (not in R output)
    if (`t' == -1) {
        display as text %14.0f `t' %16s "(ref)" %16s "(ref)" %14s "." %16s "." %16s "." %14s "." %8s "---"
        continue
    }

    local ++r_row
    local r_coef = r_betas[`r_row', 2]
    local r_se = r_betas[`r_row', 3]

    local coef_diff = abs(`stata_coef' - `r_coef')
    local se_diff = abs(`stata_se' - `r_se')

    local ++n_tests
    local pass = (`coef_diff' < `tol') & (`se_diff' < `tol')
    if (`pass') {
        local ++n_pass
        local status "PASS"
    }
    else {
        local status "FAIL"
    }

    if (`coef_diff' > `max_coef_diff') local max_coef_diff = `coef_diff'
    if (`se_diff' > `max_se_diff') local max_se_diff = `se_diff'

    display as text %14.0f `t' %16.10f `stata_coef' %16.10f `r_coef' %14.2e `coef_diff' %16.10f `stata_se' %16.10f `r_se' %14.2e `se_diff' %8s "  [`status']"
}

// Verdict
display as text ""
display as text "{hline 72}"
if (`n_pass' == `n_tests') {
    display as result "CROSS-LANGUAGE: ALL `n_tests' TESTS PASSED"
    display as result "  Max coef diff = " %10.2e `max_coef_diff'
    display as result "  Max SE diff   = " %10.2e `max_se_diff'
    display as result "Stata DLM matches R DLM to < `tol'"
}
else {
    local n_fail = `n_tests' - `n_pass'
    display as error "CROSS-LANGUAGE: `n_fail' of `n_tests' TESTS FAILED"
    display as error "  Max coef diff = " %10.2e `max_coef_diff'
    display as error "  Max SE diff   = " %10.2e `max_se_diff'
    error 9
}
display as text "{hline 72}"

// ============================================================================
// Step 5: Additional test with from=-2, to=4, ref=-1
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "STEP 5: Cross-language test with from=-2, to=4, ref=-1"
display as text "{hline 72}"

// Reload data
import delimited using "`stata_dir'/test_shared_data.csv", clear varnames(1)

// Stata DLM
dlm outcome, exposure(post) unit(unit) time(time) from(-2) to(4) ref(-1)
matrix stata_betas2 = e(betas)

// R script for this config
tempname fh2
file open `fh2' using "`stata_dir'/test_r_runner2.R", write replace
file write `fh2' `"library(dlm)"' _n
file write `fh2' `"library(dplyr)"' _n
file write `fh2' `"args <- commandArgs(trailingOnly = TRUE)"' _n
file write `fh2' `"df <- read.csv(args[1])"' _n
file write `fh2' `"data <- df %>% select(unit, time, outcome)"' _n
file write `fh2' `"exposure_data <- df %>% select(unit, time, post) %>% distinct()"' _n
file write `fh2' `"mod <- distributed_lags_model("' _n
file write `fh2' `"    data = data, exposure_data = exposure_data,"' _n
file write `fh2' `"    from_rt = -2, to_rt = 4,"' _n
file write `fh2' `"    outcome = 'outcome', exposure = 'post',"' _n
file write `fh2' `"    unit = 'unit', time = 'time', ref_period = -1)"' _n
file write `fh2' `"write.csv(mod[['betas']], args[2], row.names = FALSE)"' _n
file close `fh2'

shell Rscript "`stata_dir'/test_r_runner2.R" "`stata_dir'/test_shared_data.csv" "`stata_dir'/test_r_results2.csv"

// Load R results
preserve
import delimited using "`stata_dir'/test_r_results2.csv", clear varnames(1)
local r_nrows2 = _N
matrix r_betas2 = J(`r_nrows2', 3, .)
forvalues i = 1/`r_nrows2' {
    matrix r_betas2[`i', 1] = time_to_event[`i']
    matrix r_betas2[`i', 2] = coef[`i']
    matrix r_betas2[`i', 3] = se[`i']
}
restore

// Compare
local n_tests2 = 0
local n_pass2 = 0
local max_coef_diff2 = 0
local max_se_diff2 = 0
local nrows2 = rowsof(stata_betas2)

local r_row = 0
forvalues i = 1/`nrows2' {
    local t = stata_betas2[`i', 1]
    local stata_coef = stata_betas2[`i', 2]
    local stata_se = stata_betas2[`i', 3]

    if (`t' == -1) continue

    local ++r_row
    local r_coef = r_betas2[`r_row', 2]
    local r_se = r_betas2[`r_row', 3]

    local coef_diff = abs(`stata_coef' - `r_coef')
    local se_diff = abs(`stata_se' - `r_se')

    local ++n_tests2
    local pass = (`coef_diff' < `tol') & (`se_diff' < `tol')
    if (`pass') {
        local ++n_pass2
        local status "PASS"
    }
    else {
        local status "FAIL"
    }

    if (`coef_diff' > `max_coef_diff2') local max_coef_diff2 = `coef_diff'
    if (`se_diff' > `max_se_diff2') local max_se_diff2 = `se_diff'

    display as text "  t=" %3.0f `t' ///
        "  coef_diff=" %12.2e `coef_diff' ///
        "  se_diff=" %12.2e `se_diff' ///
        "  [`status']"
}

display as text ""
if (`n_pass2' == `n_tests2') {
    display as result "CROSS-LANGUAGE (from=-2, to=4): ALL `n_tests2' TESTS PASSED"
    display as result "  Max coef diff = " %10.2e `max_coef_diff2'
    display as result "  Max SE diff   = " %10.2e `max_se_diff2'
}
else {
    local n_fail2 = `n_tests2' - `n_pass2'
    display as error "CROSS-LANGUAGE (from=-2, to=4): `n_fail2' of `n_tests2' TESTS FAILED"
    error 9
}
display as text "{hline 72}"

// Cleanup temp files
erase "`stata_dir'/test_shared_data.csv"
erase "`stata_dir'/test_stata_results.csv"
erase "`stata_dir'/test_r_results.csv"
erase "`stata_dir'/test_r_results2.csv"
erase "`stata_dir'/test_r_runner.R"
erase "`stata_dir'/test_r_runner2.R"

display as text ""
display as text "All cross-language equivalence tests completed successfully."
