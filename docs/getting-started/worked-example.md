# Worked Example

This page walks through a complete analysis: we generate panel data, run both a **canonical event study** and a **DLM**, compare the results, and plot them. You can copy and paste the code directly into R or Stata.

## The Setup

We generate a balanced panel of 500 units observed over 20 time periods. About 40% of units receive treatment starting at time 7, 8, or 9 (randomly drawn). The true treatment effect is **−3**: after treatment onset, the outcome drops by 3 units on average.

Because treatment is binary and absorbing (turns on once and stays on), both the canonical event study and the DLM should produce **identical estimates**. This example demonstrates that.

## Step 1: Generate the Data

=== "R"

    ```r
    library(dlm)
    library(dplyr)
    library(ggplot2)

    # Generate test data (true effect = -3)
    df <- generate_data(seed = 42, n_groups = 500, n_times = 20, treat_prob = 0.4)

    # Separate outcome and exposure data (required by the R interface)
    outcome_data <- df %>% select(group, time, outcome)
    exposure_data <- df %>% select(group, time, post) %>% distinct()
    ```

=== "Stata"

    ```stata
    * Generate test data (true effect = -3)
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    ```

The generated data has columns `unit` (or `group`), `time`, `outcome`, `post` (binary treatment indicator), and `years_to_treatment` (event time, with −1000 for never-treated units).

## Step 2: Run the Canonical Event Study

First, we estimate a standard event study using binned-endpoint event-time dummies. This is the approach most researchers are familiar with: create dummy variables for each event-time period (e.g., $D_{-3}, D_{-2}, D_0, D_1, D_2, D_3$), omit the reference period ($t = -1$), and regress the outcome on these dummies with unit and time fixed effects.

The endpoints are "binned" — the $D_{-3}$ dummy equals 1 for all units at event time $-3$ **or earlier**, and $D_3$ equals 1 for all units at event time $3$ **or later**.

=== "R"

    ```r
    # Run a binned-endpoint event study
    es <- standard_twfe_for_comparison(
      data = df,
      from_rt = -3, to_rt = 3,
      outcome = "outcome",
      time = "time", unit = "group",
      time_to_treatment = "years_to_treatment",
      treat = "treat",
      ref_period = -1
    )

    # Event-study coefficients
    es$betas
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
    * Determine the estimation window
    * The DLM uses observations in [min_time + to, max_time - (abs(from) - 1)]
    * With from=-3, to=3, times 1-20: estimation window is [4, 18]
    summarize time, meanonly
    local tlo = r(min) + 3
    local thi = r(max) - 2

    * Create binned event-time dummies
    gen byte es_m3 = (years_to_treatment <= -3) & (years_to_treatment != -1000)
    gen byte es_m2 = (years_to_treatment == -2)
    * ref = -1 is omitted
    gen byte es_0  = (years_to_treatment == 0)
    gen byte es_1  = (years_to_treatment == 1)
    gen byte es_2  = (years_to_treatment == 2)
    gen byte es_3  = (years_to_treatment >= 3) & (years_to_treatment != -1000)

    * Run the event study
    reghdfe outcome es_m3 es_m2 es_0 es_1 es_2 es_3 ///
        if (time >= `tlo') & (time <= `thi'), ///
        absorb(unit time) vce(cluster unit)
    ```

    Output (abridged):

    ```
    HDFE Linear regression
    Absorbing 2 HDFE groups

          outcome |      Coef.   Std. Err.      t    P>|t|
    --------------+--------------------------------------------
            es_m3 |  -0.063959    0.483545    -0.13   0.895
            es_m2 |   0.094669    0.480424     0.20   0.844
             es_0 |  -2.763973    0.461907    -5.98   0.000
             es_1 |  -3.094282    0.520414    -5.95   0.000
             es_2 |  -2.707691    0.554940    -4.88   0.000
             es_3 |  -3.256921    0.426200    -7.64   0.000
    ```

**Interpretation.** The pre-treatment coefficients at $t = -3$ and $t = -2$ are near zero (−0.06 and 0.09), consistent with parallel trends. The post-treatment coefficients cluster around −3, matching the true effect.

## Step 3: Run the DLM

Now we estimate the same relationship using the DLM. Instead of event-time dummies, the DLM regresses the outcome on **leads and lags of the treatment variable itself**, then transforms the resulting "gamma" coefficients into "beta" coefficients via cumulative summation. The betas have the same interpretation as event-study coefficients.

=== "R"

    ```r
    # Estimate DLM with event window [-3, 3]
    mod <- distributed_lags_model(
      data = outcome_data,
      exposure_data = exposure_data,
      from_rt = -3, to_rt = 3,
      outcome = "outcome", exposure = "post",
      unit = "group", time = "time"
    )

    # DLM beta coefficients
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
    * Estimate DLM with event window [-3, 3]
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

The DLM betas are identical to the event-study coefficients from Step 2.

## Step 4: Verify Equivalence

Let's confirm the equivalence numerically.

=== "R"

    ```r
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

    cat("Max difference:", max(comparison$difference), "\n")
    # Max difference: 1.110223e-15
    ```

=== "Stata"

    ```stata
    * Save DLM results
    matrix dlm_b = e(betas)

    * Compare
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

**The estimates match to machine precision** (~$10^{-15}$). This confirms the theoretical result from Schmidheiny & Siegloch (2023): for a binary absorbing treatment, the DLM is a numerically identical reparametrization of the canonical binned-endpoint event study.

## Step 5: Plot the Results

Both approaches produce the same event-study plot — coefficients by event time, with the reference period normalized to zero.

=== "R"

    ```r
    # The DLM object includes a ready-made event-study plot
    mod$plot

    # Or customize it
    mod$plot +
      labs(title = "Event-Study Plot (from DLM)", x = "Periods to Treatment") +
      theme_minimal()

    # You can also build the plot manually from the event-study results
    ggplot(es$betas, aes(x = time_to_event, y = coef)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +
      geom_pointrange(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se)) +
      labs(
        title = "Event-Study Plot (from canonical ES)",
        x = "Periods to Treatment", y = "Coefficient"
      ) +
      theme_minimal()
    ```

=== "Stata"

    ```stata
    * Plot the event-study coefficients from the DLM
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
           title("Event-Study Plot") legend(off)
    restore

    drop es_*
    ```

The plot shows the classic event-study shape: flat pre-trends (coefficients near zero before treatment) and a sharp drop to around −3 at treatment onset, consistent with the true effect of −3.

## Why This Matters

This example uses a binary absorbing treatment — the special case where the canonical event study works. Both approaches give the same answer, so why bother with DLMs?

The DLM's advantage appears when the treatment is **continuous** — e.g., a tax rate that changes from 0% to 5% to 8.5%, a minimum wage that increases in steps, or a policy dosage that varies across units and over time. In those settings, you can't create event-time dummies (there is no single "event"), but the DLM's lead/lag framework handles it naturally. See the [Theory](../theory/background.md) page for a concrete illustration with continuous treatment data.
