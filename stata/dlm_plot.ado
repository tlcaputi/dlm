*! version 1.0.0
*! Plot event study results from dlm estimation
*! Creates publication-quality event study plot from e(betas)
*!
*! Syntax:
*!   dlm_plot [, title(string) xtitle(string) ytitle(string) ///
*!       from_label(string) to_label(string) saving(string)]
*!
*! Must be run after dlm estimation (uses e(betas), e(N), e(ref_period))

program define dlm_plot
    version 15.0

    syntax [, TItle(string) XTItle(string) YTItle(string) ///
             FROM_label(string) TO_label(string) ///
             SAVing(string)]

    // =========================================================================
    // Validate that dlm was run
    // =========================================================================
    if ("`e(cmd)'" != "dlm") {
        display as error "dlm_plot requires prior dlm estimation"
        exit 301
    }

    // =========================================================================
    // Set defaults
    // =========================================================================
    if (`"`title'"' == "") local title "Event Study (DLM)"
    if (`"`xtitle'"' == "") {
        local expo_var "`e(exposure)'"
        local pretty_expo = proper(subinstr("`expo_var'", "_", " ", .))
        local xtitle "Time to Unit Change in `pretty_expo'"
    }
    if (`"`ytitle'"' == "") local ytitle "Coefficient"

    // Build caption with N and optional from/to labels
    local N_obs = e(N)
    local caption `"N=`N_obs'"'
    if (`"`from_label'"' != "" & `"`to_label'"' != "") {
        local caption `"`caption' | From `from_label' To `to_label'"'
    }

    // =========================================================================
    // Extract betas into variables for plotting
    // =========================================================================
    tempname b
    matrix `b' = e(betas)
    local ref_period = e(ref_period)

    preserve
    clear
    local nr = rowsof(`b')
    quietly set obs `nr'
    quietly gen double time_to_event = .
    quietly gen double coef = .
    quietly gen double ci_lo = .
    quietly gen double ci_hi = .
    forvalues i = 1/`nr' {
        quietly replace time_to_event = `b'[`i', 1] in `i'
        quietly replace coef = `b'[`i', 2] in `i'
        quietly replace ci_lo = `b'[`i', 4] in `i'
        quietly replace ci_hi = `b'[`i', 5] in `i'
    }

    // =========================================================================
    // Create plot with nc700-style formatting
    // =========================================================================
    // Colors: black line/points, grey20 ribbon at 20% opacity
    // Red dashed zero line, gray dashed reference line
    local vline_x = `ref_period' + 0.5

    twoway (rarea ci_lo ci_hi time_to_event, ///
                fcolor(gs3%20) lcolor(none) sort) ///
           (line coef time_to_event, ///
                lcolor(black) lwidth(medthick) sort) ///
           (scatter coef time_to_event, ///
                mcolor(black) msymbol(circle) msize(medsmall)), ///
           yline(0, lpattern(dash) lcolor(red)) ///
           xline(`vline_x', lpattern(dash) lcolor(gs8)) ///
           xtitle(`"`xtitle'"') ytitle(`"`ytitle'"') ///
           title(`"`title'"', size(medium) position(11)) ///
           note(`"`caption'"', size(vsmall) position(5)) ///
           legend(off) ///
           graphregion(color(white)) plotregion(color(white)) ///
           scheme(s2color)

    if (`"`saving'"' != "") {
        graph export `"`saving'"', replace
    }

    restore
end
