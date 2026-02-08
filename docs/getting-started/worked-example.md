# Worked Example

This page presents a complete, self-contained example that generates test data, estimates a DLM, runs the equivalent canonical event study, and verifies that the two produce identical results. You can copy and paste the code below directly into R or Stata.

## The Setup

We generate a balanced panel of 500 units observed over 20 time periods. About 40% of units are randomly assigned to treatment. Treated units begin treatment at time 7, 8, or 9 (uniformly drawn). The true treatment effect is −3: after treatment onset, the outcome drops by 3 units on average.

Because treatment is binary and absorbing (turns on and never turns off), both the DLM and a canonical event study with binned endpoints should produce identical estimates. This example verifies that.

## Step 1: Generate Data and Estimate the DLM

=== "R"

    ```r
    library(dlm)
    library(dplyr)

    # Generate test data
    df <- generate_data(seed = 42, n_groups = 500, n_times = 20, treat_prob = 0.4)

    # Separate outcome and exposure data (required by the R interface)
    outcome_data <- df %>% select(group, time, outcome)
    exposure_data <- df %>% select(group, time, post) %>% distinct()

    # Estimate DLM with event window [-3, 3] and reference period -1
    mod <- distributed_lags_model(
      data = outcome_data,
      exposure_data = exposure_data,
      from_rt = -3, to_rt = 3,
      outcome = "outcome", exposure = "post",
      unit = "group", time = "time"
    )

    # View beta coefficients
    mod$betas
    #   time_to_event       coef        se
    # 1            -3 -0.0639590 0.4835454
    # 2            -2  0.0946686 0.4804242
    # 3             0 -2.7639731 0.4619074
    # 4             1 -3.0942824 0.5204137
    # 5             2 -2.7076914 0.5549396
    # 6             3 -3.2569207 0.4262002
    ```

=== "Stata"

    ```stata
    * Generate test data
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear

    * Estimate DLM with event window [-3, 3] and reference period -1
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
    ```

    Output:

    ```
    ------------------------------------------------------------------------
    Distributed Lag Model (Schmidheiny & Siegloch 2023)
    ------------------------------------------------------------------------
      Outcome:    outcome
      Exposure:   post
      Window:     [-3, 3], ref = -1
      N obs:      7500
      N clusters: 500
    ------------------------------------------------------------------------

            Time        Coef          SE     95% CI lo     95% CI hi
    ------------------------------------------------------------------------
              -3   -0.063959    0.483545     -1.011708      0.883790
              -2    0.094669    0.480424     -0.846962      1.036299
              -1       (ref)           .             .             .
               0   -2.763973    0.461907     -3.669310     -1.858636
               1   -3.094282    0.520414     -4.114294     -2.074271
               2   -2.707691    0.554940     -3.795373     -1.620010
               3   -3.256921    0.426200     -4.092272     -2.421569
    ------------------------------------------------------------------------
    ```

**Interpretation.** The pre-treatment betas at $t = -3$ and $t = -2$ are close to zero (−0.06 and 0.09), consistent with the parallel trends assumption. The post-treatment betas at $t = 0$ through $t = 3$ cluster around −3, matching the true treatment effect of −3. The reference period $t = -1$ is normalized to zero by construction.

## Step 2: Run the Equivalent Event Study

Now we run a standard binned-endpoint event study on the same data and verify that the coefficients are identical.

=== "R"

    ```r
    # Run the equivalent binned-endpoint event study
    es <- standard_twfe_for_comparison(
      data = df,
      from_rt = -3, to_rt = 3,
      outcome = "outcome",
      time = "time", unit = "group",
      time_to_treatment = "years_to_treatment",
      treat = "treat",
      ref_period = -1
    )

    # Compare DLM betas to event-study betas
    comparison <- data.frame(
      time = mod$betas$time_to_event,
      dlm_beta = mod$betas$coef,
      es_beta = es$betas$coef,
      difference = abs(mod$betas$coef - es$betas$coef)
    )
    print(comparison)
    #   time    dlm_beta     es_beta   difference
    # 1   -3 -0.06395904 -0.06395904 1.110223e-15
    # 2   -2  0.09466863  0.09466863 4.440892e-16
    # 3    0 -2.76397307 -2.76397307 4.440892e-16
    # 4    1 -3.09428243 -3.09428243 8.881784e-16
    # 5    2 -2.70769138 -2.70769138 4.440892e-16
    # 6    3 -3.25692067 -3.25692067 0.000000e+00

    # Maximum difference (should be ~1e-14 or smaller)
    cat("Max difference:", max(comparison$difference), "\n")
    ```

=== "Stata"

    ```stata
    * Save DLM results
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3) ref(-1)
    matrix dlm_b = e(betas)

    * Determine time window (same observations the DLM uses)
    summarize time, meanonly
    local tlo = r(min) + 3   // drop first 3 periods (lags need history)
    local thi = r(max) - 2   // drop last 2 periods (leads need future)

    * Create binned event-time dummies
    gen byte es_m3 = (years_to_treatment <= -3) & (years_to_treatment != -1000)
    gen byte es_m2 = (years_to_treatment == -2)
    * ref = -1 is omitted
    gen byte es_0  = (years_to_treatment == 0)
    gen byte es_1  = (years_to_treatment == 1)
    gen byte es_2  = (years_to_treatment == 2)
    gen byte es_3  = (years_to_treatment >= 3) & (years_to_treatment != -1000)

    * Run the event study on the same time window
    reghdfe outcome es_m3 es_m2 es_0 es_1 es_2 es_3 ///
        if (time >= `tlo') & (time <= `thi'), ///
        absorb(unit time) vce(cluster unit)

    * Compare coefficients
    display ""
    display "DLM vs Event Study Comparison:"
    display "  Period  DLM beta     ES beta      Difference"
    display "  ------  ----------   ----------   ----------"

    local periods  "-3 -2 0 1 2 3"
    local es_vars  "es_m3 es_m2 es_0 es_1 es_2 es_3"
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
    ```

    Output:

    ```
    DLM vs Event Study Comparison:
      Period  DLM beta     ES beta      Difference
      ------  ----------   ----------   ----------
         -3   -0.063959    -0.063959     0.00e+00
         -2    0.094669     0.094669     0.00e+00
          0   -2.763973    -2.763973     0.00e+00
          1   -3.094282    -3.094282     0.00e+00
          2   -2.707691    -2.707691     0.00e+00
          3   -3.256921    -3.256921     0.00e+00
    ```

**The DLM and event-study estimates match to machine precision.** This confirms the theoretical equivalence result from Schmidheiny & Siegloch (2023): for a binary absorbing treatment, the DLM is a numerically identical reparametrization of the canonical binned-endpoint event study.

## Step 3: View the Event-Study Plot

=== "R"

    ```r
    # Built-in plot (ggplot2 object)
    mod$plot

    # Customize
    mod$plot +
      ggplot2::labs(title = "Treatment Effect Dynamics", x = "Periods to Treatment") +
      ggplot2::theme_minimal()
    ```

=== "Stata"

    ```stata
    dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)
    matrix b = e(betas)

    preserve
    clear
    local nr = rowsof(b)
    set obs `nr'
    gen time_to_event = .
    gen coef = .
    gen ci_lo = .
    gen ci_hi = .
    forvalues i = 1/`nr' {
        replace time_to_event = b[`i', 1] in `i'
        replace coef = b[`i', 2] in `i'
        replace ci_lo = b[`i', 4] in `i'
        replace ci_hi = b[`i', 5] in `i'
    }

    twoway (rcap ci_lo ci_hi time_to_event, lcolor(navy)) ///
           (scatter coef time_to_event, mcolor(navy) msymbol(circle)), ///
           yline(0, lpattern(dash) lcolor(gray)) ///
           xline(-0.5, lpattern(dash) lcolor(gray)) ///
           xtitle("Periods to Treatment") ytitle("Coefficient") ///
           title("Treatment Effect Dynamics") legend(off)
    restore
    ```

## Why This Matters

This example uses a binary absorbing treatment, which is the special case where both approaches work. The DLM's advantage appears when the treatment is **continuous** — e.g., a tax rate, minimum wage, or policy dosage that changes in magnitude over time. In those settings, event-time dummies don't apply, but the DLM's lead/lag approach works identically. See the [Theory](../theory/background.md) page for a concrete illustration with continuous treatment data.
