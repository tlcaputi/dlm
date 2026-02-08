/*
    test_package_install.do
    =======================
    Tests the dlm Stata package installation workflow:
    1. net install from local directory
    2. help dlm works
    3. dlm command works after install (from any directory)
    4. e() returns are correct
    5. Error handling works
    6. net uninstall works cleanly
*/

clear all
set more off

// Store the package source directory
local pkg_dir = "`c(pwd)'"
display as text "Package source: `pkg_dir'"

// ============================================================================
// TEST 1: net install from local directory
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 1: net install dlm from local directory"
display as text "{hline 72}"

// First uninstall if already installed (clean slate)
capture ado uninstall dlm
display as text "  (any previous install removed)"

// Install from local directory
net install dlm, from("`pkg_dir'") replace
display as result "  PASS: net install dlm succeeded"

// ============================================================================
// TEST 2: Verify installed files
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 2: Verify installed files exist"
display as text "{hline 72}"

// Check dlm.ado is findable
capture which dlm
if _rc {
    display as error "  FAIL: dlm.ado not found after install"
    error 9
}
display as result "  PASS: which dlm -> `r(fn)'"

// Check dlm_gen_data.ado is findable
capture which dlm_gen_data
if _rc {
    display as error "  FAIL: dlm_gen_data.ado not found after install"
    error 9
}
display as result "  PASS: which dlm_gen_data -> `r(fn)'"

// ============================================================================
// TEST 3: help dlm works
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 3: help dlm renders without error"
display as text "{hline 72}"

capture noisily help dlm
if _rc {
    display as error "  FAIL: help dlm failed with rc=" _rc
    error 9
}
display as result "  PASS: help dlm rendered successfully"

// ============================================================================
// TEST 4: dlm command works post-install (no adopath hack needed)
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 4: dlm works after install (no adopath needed)"
display as text "{hline 72}"

// Reset adopath to defaults (remove any "." entries)
adopath - "."

// Generate data and run DLM — this should work purely from the installed package
dlm_gen_data, n_groups(200) n_times(15) treat_prob(0.3) seed(99) clear
display as text "  Data generated: " _N " obs"

dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) ref(-1)
display as result "  PASS: dlm command ran successfully post-install"

// ============================================================================
// TEST 5: e() return values are correct
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 5: e() return values"
display as text "{hline 72}"

local tests_passed = 0
local tests_total = 0

// e(cmd) should be "dlm"
local ++tests_total
if ("`e(cmd)'" == "dlm") {
    local ++tests_passed
    display as result `"  PASS: e(cmd) = "`e(cmd)'""'
}
else {
    display as error `"  FAIL: e(cmd) = "`e(cmd)'" (expected "dlm")"'
}

// e(N) should be positive
local ++tests_total
if (e(N) > 0) {
    local ++tests_passed
    display as result "  PASS: e(N) = " e(N)
}
else {
    display as error "  FAIL: e(N) = " e(N)
}

// e(N_clust) should be positive
local ++tests_total
if (e(N_clust) > 0) {
    local ++tests_passed
    display as result "  PASS: e(N_clust) = " e(N_clust)
}
else {
    display as error "  FAIL: e(N_clust) = " e(N_clust)
}

// e(from) should be -3
local ++tests_total
if (e(from) == -3) {
    local ++tests_passed
    display as result "  PASS: e(from) = " e(from)
}
else {
    display as error "  FAIL: e(from) = " e(from) " (expected -3)"
}

// e(to) should be 3
local ++tests_total
if (e(to) == 3) {
    local ++tests_passed
    display as result "  PASS: e(to) = " e(to)
}
else {
    display as error "  FAIL: e(to) = " e(to) " (expected 3)"
}

// e(ref_period) should be -1
local ++tests_total
if (e(ref_period) == -1) {
    local ++tests_passed
    display as result "  PASS: e(ref_period) = " e(ref_period)
}
else {
    display as error "  FAIL: e(ref_period) = " e(ref_period) " (expected -1)"
}

// e(outcome) should be "outcome"
local ++tests_total
if ("`e(outcome)'" == "outcome") {
    local ++tests_passed
    display as result `"  PASS: e(outcome) = "`e(outcome)'""'
}
else {
    display as error `"  FAIL: e(outcome) = "`e(outcome)'" (expected "outcome")"'
}

// e(exposure) should be "post"
local ++tests_total
if ("`e(exposure)'" == "post") {
    local ++tests_passed
    display as result `"  PASS: e(exposure) = "`e(exposure)'""'
}
else {
    display as error `"  FAIL: e(exposure) = "`e(exposure)'" (expected "post")"'
}

// e(betas) matrix should be 7x5 (from=-3 to=3 = 7 rows, 5 columns)
local ++tests_total
matrix betas_check = e(betas)
local brows = rowsof(betas_check)
local bcols = colsof(betas_check)
if (`brows' == 7 & `bcols' == 5) {
    local ++tests_passed
    display as result "  PASS: e(betas) is `brows'x`bcols'"
}
else {
    display as error "  FAIL: e(betas) is `brows'x`bcols' (expected 7x5)"
}

// Reference period row should have coef=0, se=0
local ++tests_total
// ref=-1 is row 3 (from=-3 is row 1, -2 is row 2, -1 is row 3)
local ref_coef = betas_check[3, 2]
local ref_se = betas_check[3, 3]
if (`ref_coef' == 0 & `ref_se' == 0) {
    local ++tests_passed
    display as result "  PASS: ref period has coef=0, se=0"
}
else {
    display as error "  FAIL: ref period coef=`ref_coef', se=`ref_se' (expected 0,0)"
}

// e(gamma) should be 1x6 (abs(-3)+3 = 6 gamma coefficients)
local ++tests_total
matrix gamma_check = e(gamma)
local gcols = colsof(gamma_check)
if (`gcols' == 6) {
    local ++tests_passed
    display as result "  PASS: e(gamma) has `gcols' columns"
}
else {
    display as error "  FAIL: e(gamma) has `gcols' columns (expected 6)"
}

// e(gamma_V) should be 6x6
local ++tests_total
matrix gv_check = e(gamma_V)
local gvrows = rowsof(gv_check)
local gvcols = colsof(gv_check)
if (`gvrows' == 6 & `gvcols' == 6) {
    local ++tests_passed
    display as result "  PASS: e(gamma_V) is `gvrows'x`gvcols'"
}
else {
    display as error "  FAIL: e(gamma_V) is `gvrows'x`gvcols' (expected 6x6)"
}

display as text ""
if (`tests_passed' == `tests_total') {
    display as result "  ALL `tests_total' e() TESTS PASSED"
}
else {
    local n_fail = `tests_total' - `tests_passed'
    display as error "  `n_fail' of `tests_total' e() TESTS FAILED"
    error 9
}

// ============================================================================
// TEST 6: Error handling
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 6: Error handling"
display as text "{hline 72}"

dlm_gen_data, n_groups(100) n_times(20) treat_prob(0.4) seed(1) clear

// from >= to should error
capture dlm outcome, exposure(post) unit(unit) time(time) from(3) to(-3)
if _rc {
    display as result "  PASS: from >= to correctly rejected (rc=" _rc ")"
}
else {
    display as error "  FAIL: from >= to was not rejected"
    error 9
}

// from >= 0 should error
capture dlm outcome, exposure(post) unit(unit) time(time) from(0) to(3)
if _rc {
    display as result "  PASS: from >= 0 correctly rejected (rc=" _rc ")"
}
else {
    display as error "  FAIL: from >= 0 was not rejected"
    error 9
}

// to <= 0 should error
capture dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(0)
if _rc {
    display as result "  PASS: to <= 0 correctly rejected (rc=" _rc ")"
}
else {
    display as error "  FAIL: to <= 0 was not rejected"
    error 9
}

// ref out of range should error
capture dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) ref(-5)
if _rc {
    display as result "  PASS: ref out of range correctly rejected (rc=" _rc ")"
}
else {
    display as error "  FAIL: ref out of range was not rejected"
    error 9
}

// ============================================================================
// TEST 7: Data preservation (dlm should not modify original data)
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 7: Data preservation"
display as text "{hline 72}"

dlm_gen_data, n_groups(100) n_times(20) treat_prob(0.4) seed(1) clear
local orig_N = _N
local orig_vars : char _dta[_orig_varcount]
describe, short
local nvars_before = r(k)

dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)

local post_N = _N
describe, short
local nvars_after = r(k)

if (`orig_N' == `post_N') {
    display as result "  PASS: N unchanged (`orig_N' -> `post_N')"
}
else {
    display as error "  FAIL: N changed (`orig_N' -> `post_N')"
    error 9
}

if (`nvars_before' == `nvars_after') {
    display as result "  PASS: variable count unchanged (`nvars_before' -> `nvars_after')"
}
else {
    display as error "  FAIL: variable count changed (`nvars_before' -> `nvars_after')"
    error 9
}

// Check no _dlm_ variables leaked into original data
capture confirm variable _dlm_lead1
if _rc {
    display as result "  PASS: no _dlm_ temporary variables in data"
}
else {
    display as error "  FAIL: _dlm_ temporary variables leaked into original data"
    error 9
}

// ============================================================================
// TEST 8: if/in support
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 8: if/in subsetting"
display as text "{hline 72}"

dlm_gen_data, n_groups(200) n_times(20) treat_prob(0.4) seed(1) clear

// Run with if condition
capture noisily dlm outcome if unit <= 100, exposure(post) unit(unit) time(time) from(-3) to(3)
if !_rc {
    display as result "  PASS: dlm with if condition works (N=" e(N) ")"
}
else {
    display as error "  FAIL: dlm with if condition failed (rc=" _rc ")"
    error 9
}

// ============================================================================
// TEST 9: net uninstall
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 9: net uninstall dlm"
display as text "{hline 72}"

ado uninstall dlm
display as text "  Uninstalled dlm"

// Verify it's gone
capture which dlm
if _rc {
    display as result "  PASS: dlm correctly not found after uninstall"
}
else {
    display as error "  FAIL: dlm still found after uninstall"
    error 9
}

// ============================================================================
// TEST 10: Reinstall for user
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "TEST 10: Reinstall for continued use"
display as text "{hline 72}"

net install dlm, from("`pkg_dir'") replace
display as result "  PASS: Reinstalled dlm"

// Final verify
dlm_gen_data, n_groups(100) n_times(15) seed(42) clear
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
display as result "  PASS: Post-reinstall run successful"

// ============================================================================
// Summary
// ============================================================================
display as text ""
display as text "{hline 72}"
display as result "ALL PACKAGE INSTALL TESTS PASSED"
display as text "{hline 72}"
display as text ""
display as text "  1. net install from local directory     - PASS"
display as text "  2. Installed files verified             - PASS"
display as text "  3. help dlm works                      - PASS"
display as text "  4. dlm works post-install               - PASS"
display as text "  5. e() return values correct            - PASS"
display as text "  6. Error handling                       - PASS"
display as text "  7. Data preservation                   - PASS"
display as text "  8. if/in subsetting                    - PASS"
display as text "  9. net uninstall works                 - PASS"
display as text " 10. Reinstall + verify                  - PASS"
display as text "{hline 72}"
