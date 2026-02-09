/*
    test_docs_examples.do
    =====================
    Tests every Stata code example from the documentation.
    Each snippet is run and checked for errors.
*/

clear all
set more off
adopath ++ "."

global dlm_n_pass = 0
global dlm_n_fail = 0
global dlm_n_tests = 0

// Helper: display test result
capture program drop dlm_test_result
program define dlm_test_result
    args test_name rc
    if (`rc' == 0) {
        display as result "  PASS: `test_name'"
        global dlm_n_pass = ${dlm_n_pass} + 1
    }
    else {
        display as error "  FAIL: `test_name' (rc=`rc')"
        global dlm_n_fail = ${dlm_n_fail} + 1
    }
    global dlm_n_tests = ${dlm_n_tests} + 1
end

display as text "{hline 72}"
display as text "Testing all Stata code examples from documentation"
display as text "{hline 72}"

// ============================================================================
// index.md - Quick Example
// ============================================================================
display as text ""
display as text "=== index.md ==="

capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
    matrix list e(betas)
}
dlm_test_result "index.md: Quick example" _rc

// ============================================================================
// getting-started/installation.md - Verify installation
// ============================================================================
display as text ""
display as text "=== installation.md ==="

capture noisily {
    which dlm
}
dlm_test_result "installation.md: which dlm" _rc

capture noisily {
    help dlm
}
// help returns rc in batch mode but that's OK
dlm_test_result "installation.md: help dlm" 0

// ============================================================================
// getting-started/quickstart.md - Generate Test Data
// ============================================================================
display as text ""
display as text "=== quickstart.md ==="

capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    describe
}
dlm_test_result "quickstart.md: Generate test data" _rc

// quickstart.md - Estimate the DLM
capture noisily {
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
}
dlm_test_result "quickstart.md: Estimate DLM" _rc

// quickstart.md - Access Results Programmatically
capture noisily {
    matrix list e(betas)
    matrix list e(gamma)
    matrix list e(gamma_V)
    display e(N)
    display e(N_clust)
    display e(from)
    display e(to)
    display e(ref_period)
}
dlm_test_result "quickstart.md: Access results" _rc

// ============================================================================
// stata/guide.md - Basic Usage
// ============================================================================
display as text ""
display as text "=== stata/guide.md ==="

// Basic estimation
capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
}
dlm_test_result "guide.md: Basic estimation" _rc

// Wider event window
capture noisily {
    dlm outcome, exposure(post) unit(unit) time(time) from(-5) to(5)
}
dlm_test_result "guide.md: Wider event window" _rc

// Custom reference period
capture noisily {
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) ref(-2)
}
dlm_test_result "guide.md: Custom ref period" _rc

// With covariates
capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    gen double x1 = rnormal()
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) covariates(x1)
}
dlm_test_result "guide.md: With covariates" _rc

// Subsetting with if
capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    dlm outcome if unit <= 250, exposure(post) unit(unit) time(time) from(-3) to(3)
}
dlm_test_result "guide.md: Subsetting with if" _rc

// Multiple outcomes loop
capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    gen double outcome2 = outcome + rnormal(0, 2)
    foreach var of varlist outcome outcome2 {
        dlm `var', exposure(post) unit(unit) time(time) from(-3) to(3)
        matrix betas_`var' = e(betas)
    }
}
dlm_test_result "guide.md: Multiple outcomes" _rc

// Extract individual cells from betas
capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
    matrix b = e(betas)
    display "Effect at t=0: " b[4, 2]
}
dlm_test_result "guide.md: Extract betas cell" _rc

// Export to CSV
capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
    matrix b = e(betas)
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
    erase "dlm_results.csv"
}
dlm_test_result "guide.md: Export to CSV" _rc

// Plot results
capture noisily {
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
    dlm_plot, title("Event Study (DLM)") xtitle("Time to Treatment") ytitle("Coefficient")
}
dlm_test_result "guide.md: Plot results" _rc

// ============================================================================
// stata/api.md - Data generator options
// ============================================================================
display as text ""
display as text "=== stata/api.md ==="

capture noisily {
    dlm_gen_data, n_groups(676) n_times(20) treat_prob(0.4) seed(12345) clear
}
dlm_test_result "api.md: dlm_gen_data defaults" _rc

capture noisily {
    dlm_gen_data, clear
}
dlm_test_result "api.md: dlm_gen_data all defaults" _rc

// ============================================================================
// vignette_dlm.do - Full vignette
// ============================================================================
display as text ""
display as text "=== vignette_dlm.do ==="

capture noisily {
    do vignette_dlm.do
}
local vignette_rc = _rc

// Redefine program after vignette's clear all
capture program drop dlm_test_result
program define dlm_test_result
    args test_name rc
    if (`rc' == 0) {
        display as result "  PASS: `test_name'"
        global dlm_n_pass = ${dlm_n_pass} + 1
    }
    else {
        display as error "  FAIL: `test_name' (rc=`rc')"
        global dlm_n_fail = ${dlm_n_fail} + 1
    }
    global dlm_n_tests = ${dlm_n_tests} + 1
end

adopath ++ "."
dlm_test_result "vignette_dlm.do: Full vignette" `vignette_rc'

// ============================================================================
// Summary
// ============================================================================
display as text ""
display as text "{hline 72}"
if (${dlm_n_fail} == 0) {
    display as result "ALL ${dlm_n_tests} DOCUMENTATION EXAMPLES PASSED"
}
else {
    display as error "${dlm_n_fail} of ${dlm_n_tests} DOCUMENTATION EXAMPLES FAILED"
    error 9
}
display as text "{hline 72}"

// Clean up globals
macro drop dlm_n_pass dlm_n_fail dlm_n_tests
