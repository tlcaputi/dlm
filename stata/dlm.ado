*! version 1.0.0
*! Distributed Lag Models (Schmidheiny & Siegloch 2023, JAE)
*! Estimates DLM and converts gamma coefficients to cumulative beta coefficients
*!
*! Syntax:
*!   dlm outcome, exposure(varname) unit(varname) time(varname) ///
*!       from(int) to(int) [ref(int) covariates(varlist) addl_fes(varlist) verbose]
*!
*! Required: reghdfe (ssc install reghdfe)

program define dlm, eclass
    version 15.0

    // =========================================================================
    // Parse syntax
    // =========================================================================
    syntax varname(numeric) [if] [in], ///
        EXPosure(varname numeric) ///
        UNIT(varname) ///
        TIME(varname numeric) ///
        FROM(integer) ///
        TO(integer) ///
        [REF(integer -1) ///
         COVariates(varlist numeric) ///
         ADDL_fes(varlist) ///
         VERBose]

    local outcome `varlist'
    local exposure `exposure'
    local unit `unit'
    local time `time'
    local from_rt `from'
    local to_rt `to'
    local ref_period `ref'
    local covariates `covariates'
    local addl_fes `addl_fes'
    local verbose = ("`verbose'" != "")

    // =========================================================================
    // Input validation
    // =========================================================================
    if (`from_rt' >= `to_rt') {
        display as error "from() must be less than to()"
        exit 198
    }
    if (`from_rt' >= 0) {
        display as error "from() must be negative"
        exit 198
    }
    if (`to_rt' <= 0) {
        display as error "to() must be positive"
        exit 198
    }
    if (`ref_period' < `from_rt' | `ref_period' > `to_rt') {
        display as error "ref() must be between from() and to()"
        exit 198
    }

    // Check reghdfe is installed
    capture which reghdfe
    if _rc {
        display as error "reghdfe is required. Install with: ssc install reghdfe"
        exit 199
    }

    // Mark the estimation sample
    marksample touse
    markout `touse' `outcome' `exposure' `unit' `time' `covariates' `addl_fes'

    // =========================================================================
    // Step 1: Create lead and lag variables
    // =========================================================================
    if (`verbose') display as text "Creating lead and lag variables..."

    // Preserve the original data — we'll restore after estimation
    preserve

    // Keep estimation sample only
    quietly keep if `touse'

    // Sort by unit and time for by-group operations
    sort `unit' `time'

    // Track the variables we create
    local lead_vars ""
    local lag_vars ""
    local all_dlm_vars ""

    // Create leads: for from=-3, we need lead2 and lead1
    // Lead k means: exposure at time t+k (future value)
    local num_leads = abs(`from_rt') - 1
    if (`num_leads' > 0) {
        forvalues k = `num_leads'(-1)1 {
            quietly by `unit': gen double _dlm_lead`k' = `exposure'[_n + `k']
            local lead_vars "`lead_vars' _dlm_lead`k'"
        }
    }

    // Create lags: for to=3, we need lag0, lag1, lag2, lag3
    // Lag k means: exposure at time t-k (past value)
    forvalues k = 0/`to_rt' {
        quietly by `unit': gen double _dlm_lag`k' = `exposure'[_n - `k']
        local lag_vars "`lag_vars' _dlm_lag`k'"
    }

    local all_dlm_vars "`lead_vars' `lag_vars'"

    // Count: total gamma coefficients = abs(from) + to
    //   leads: abs(from)-1 variables
    //   lags:  to+1 variables
    //   total: abs(from)-1 + to+1 = abs(from) + to
    local num_vars = abs(`from_rt') + `to_rt'

    if (`verbose') {
        display as text "  Lead variables: `lead_vars'"
        display as text "  Lag variables:  `lag_vars'"
        display as text "  Total gamma coefficients: `num_vars'"
    }

    // =========================================================================
    // Step 2: Run regression with reghdfe
    // =========================================================================
    if (`verbose') display as text "Running reghdfe regression..."

    // Build absorb() list: always unit + time, plus any additional FEs
    local absorb_list "`unit' `time'"
    if ("`addl_fes'" != "") {
        local absorb_list "`absorb_list' `addl_fes'"
    }

    // Run reghdfe
    // Variable order: leads (highest to lowest), then lags (0 to highest)
    // This matches the R package ordering
    quietly reghdfe `outcome' `all_dlm_vars' `covariates', ///
        absorb(`absorb_list') vce(cluster `unit')

    local N_obs = e(N)
    local N_clust = e(N_clust)

    if (`verbose') {
        display as text "  N observations: `N_obs'"
        display as text "  N clusters: `N_clust'"
    }

    // Verify all DLM coefficients were estimated
    // e(b) includes our DLM vars + covariates + _cons
    // Check that the first num_vars columns of e(b) match our variable names
    tempname check_b
    matrix `check_b' = e(b)
    local n_cols = colsof(`check_b')
    if (`n_cols' < `num_vars') {
        display as error "Only `n_cols' total coefficients but expected at least `num_vars' DLM coefficients."
        display as error "Check your data — panel may be too short for the requested window."
        restore
        exit 498
    }

    // =========================================================================
    // Step 3: Extract gamma coefficients and VCV, transform to betas (Mata)
    // =========================================================================
    if (`verbose') display as text "Transforming gamma -> beta coefficients..."

    // Store gamma coefficients in a Mata-accessible matrix
    tempname gamma_vec gamma_vcv betas_mat

    // Extract gamma vector (first num_vars coefficients from e(b))
    matrix `gamma_vec' = e(b)
    matrix `gamma_vec' = `gamma_vec'[1, 1..`num_vars']

    // Extract gamma VCV (top-left num_vars x num_vars block of e(V))
    matrix `gamma_vcv' = e(V)
    matrix `gamma_vcv' = `gamma_vcv'[1..`num_vars', 1..`num_vars']

    // Mata block: gamma -> beta transformation
    mata: _dlm_transform("`gamma_vec'", "`gamma_vcv'", `from_rt', `to_rt', `ref_period', "`betas_mat'")

    // =========================================================================
    // Step 4: Display results
    // =========================================================================
    display as text ""
    display as text "{hline 72}"
    display as text "Distributed Lag Model (Schmidheiny & Siegloch 2023)"
    display as text "{hline 72}"
    display as text "  Outcome:    `outcome'"
    display as text "  Exposure:   `exposure'"
    display as text "  Window:     [`from_rt', `to_rt'], ref = `ref_period'"
    display as text "  N obs:      `N_obs'"
    display as text "  N clusters: `N_clust'"
    display as text "{hline 72}"
    display as text ""
    display as text %12s "Time" %12s "Coef" %12s "SE" %14s "95% CI lo" %14s "95% CI hi"
    display as text "{hline 72}"

    local nrows = rowsof(`betas_mat')
    forvalues i = 1/`nrows' {
        local t = `betas_mat'[`i', 1]
        local b = `betas_mat'[`i', 2]
        local s = `betas_mat'[`i', 3]
        local lo = `betas_mat'[`i', 4]
        local hi = `betas_mat'[`i', 5]
        if (`t' == `ref_period') {
            display as text %12.0f `t' %12s "(ref)" %12s "." %14s "." %14s "."
        }
        else {
            display as text %12.0f `t' %12.6f `b' %12.6f `s' %14.6f `lo' %14.6f `hi'
        }
    }
    display as text "{hline 72}"

    // =========================================================================
    // Step 5: Post results to e()
    // =========================================================================

    // Label the betas matrix
    matrix colnames `betas_mat' = time_to_event coef se ci_lo ci_hi

    // Label the gamma vector and VCV with variable names
    matrix colnames `gamma_vec' = `all_dlm_vars'
    matrix rownames `gamma_vcv' = `all_dlm_vars'
    matrix colnames `gamma_vcv' = `all_dlm_vars'

    // Post to e()
    ereturn clear
    ereturn post, esample(`touse')

    ereturn matrix betas = `betas_mat'
    ereturn matrix gamma = `gamma_vec'
    ereturn matrix gamma_V = `gamma_vcv'

    ereturn scalar N = `N_obs'
    ereturn scalar N_clust = `N_clust'
    ereturn scalar from = `from_rt'
    ereturn scalar to = `to_rt'
    ereturn scalar ref_period = `ref_period'

    ereturn local cmd "dlm"
    ereturn local outcome "`outcome'"
    ereturn local exposure "`exposure'"
    ereturn local unit "`unit'"
    ereturn local time "`time'"

    // Restore original data
    restore

end


// =============================================================================
// Mata: Gamma-to-Beta transformation
// =============================================================================
// This implements the cumulative sum transformation from DLM gammas to
// event-study betas, with proper variance propagation through the VCV.

mata:
void _dlm_transform(string scalar gamma_name, string scalar vcv_name, ///
                     real scalar from_rt, real scalar to_rt, ///
                     real scalar ref_period, string scalar result_name)
{
    real matrix gamma, V
    real scalar num_vars, num_before, num_after
    real vector time_to_event, coefs, ses
    real matrix betas
    real scalar i, j, idx
    real rowvector a
    real scalar w

    // Load Stata matrices into Mata
    gamma = st_matrix(gamma_name)  // 1 x num_vars
    V = st_matrix(vcv_name)        // num_vars x num_vars
    num_vars = cols(gamma)

    // Build the time_to_event vector (all periods except ref)
    // Total output rows = to - from + 1 (including ref period row with 0)
    real scalar total_periods
    total_periods = to_rt - from_rt + 1

    // Initialize output: [time_to_event, coef, se, ci_lo, ci_hi]
    betas = J(total_periods, 5, 0)

    // Fill time_to_event column (including ref period)
    for (i = 1; i <= total_periods; i++) {
        betas[i, 1] = from_rt + i - 1
    }

    // Number of periods before ref (excluding ref itself)
    num_before = ref_period - from_rt
    // Number of periods after ref (excluding ref itself)
    num_after = to_rt - ref_period

    // === Before-reference periods ===
    // Gamma ordering: lead(abs(from)-1), ..., lead1, lag0, ..., lag(to)
    // Before-ref gammas are indices 1..num_before (the leads + possibly lag0 up to ref-1)
    // Actually: num_before = ref_period - from_rt
    //   For from=-3, ref=-1: num_before = -1 - (-3) = 2
    //   Gamma indices 1..2 are the "before" gammas
    //   Beta[i] = -sum(gamma[i..num_before])  (negative reverse cumulative sum)

    // Before periods: beta[i] = -sum(gamma[i..num_before])
    for (i = 1; i <= num_before; i++) {
        // Coefficient: negative sum from i to num_before
        coefs = 0
        for (j = i; j <= num_before; j++) {
            coefs = coefs + gamma[1, j]
        }
        coefs = -coefs

        // SE: sqrt(1' * V[i..num_before, i..num_before] * 1)
        real scalar len
        len = num_before - i + 1
        a = J(1, len, 1)
        w = a * V[i..num_before, i..num_before] * a'
        ses = sqrt(w)

        // Row index in betas: period from_rt corresponds to row 1,
        // period from_rt+1 to row 2, etc.
        idx = i  // row i = period from_rt + i - 1
        betas[idx, 2] = coefs
        betas[idx, 3] = ses
        betas[idx, 4] = coefs - 1.96 * ses
        betas[idx, 5] = coefs + 1.96 * ses
    }

    // === Reference period ===
    // Row for ref period: ref_period - from_rt + 1
    idx = ref_period - from_rt + 1
    betas[idx, 2] = 0
    betas[idx, 3] = 0
    betas[idx, 4] = 0
    betas[idx, 5] = 0

    // === After-reference periods ===
    // After-ref gammas are indices (num_before+1)..num_vars
    // Beta[j] = sum(gamma[num_before+1..num_before+j]) for j=1,2,...,num_after
    for (i = 1; i <= num_after; i++) {
        // Coefficient: sum from (num_before+1) to (num_before+i)
        coefs = 0
        for (j = num_before + 1; j <= num_before + i; j++) {
            coefs = coefs + gamma[1, j]
        }

        // SE: sqrt(1' * V[start..start+i-1, start..start+i-1] * 1)
        real scalar start
        start = num_before + 1
        a = J(1, i, 1)
        w = a * V[start..(start + i - 1), start..(start + i - 1)] * a'
        ses = sqrt(w)

        // Row index: ref_period + i corresponds to row (ref_period - from_rt + 1 + i)
        idx = ref_period - from_rt + 1 + i
        betas[idx, 2] = coefs
        betas[idx, 3] = ses
        betas[idx, 4] = coefs - 1.96 * ses
        betas[idx, 5] = coefs + 1.96 * ses
    }

    // Store result back to Stata
    st_matrix(result_name, betas)
}
end
