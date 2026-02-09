# Quick Start

This page walks through a complete example in both R and Stata using generated test data. Both packages include a built-in test dataset generator so you can try the package immediately and verify that the DLM produces estimates identical to a canonical event study with binned endpoints.

## Generate Test Data

Both packages include a data generator that creates a balanced panel with staggered treatment adoption. The generated data has:

- A configurable number of units and time periods
- Treatment assigned randomly to a fraction of units
- Treatment onset at time 7, 8, or 9 (uniform)
- A treatment effect of −3 on the outcome (post-treatment)
- Gaussian noise (σ = 5)

=== "R"

    ```r
    library(dlm)
    df <- generate_data(seed = 42, n_groups = 500, n_times = 20, treat_prob = 0.4)
    str(df)
    ```

    Output:

    ```
    'data.frame':	10000 obs. of  7 variables:
     $ group             : Factor w/ 500 levels "group1","group2",..: 1 1 1 1 1 ...
     $ time              : int  1 2 3 4 5 6 7 8 9 10 ...
     $ treatment_time    : int  7 7 7 7 7 7 7 7 7 7 ...
     $ treat             : num  1 1 1 1 1 1 1 1 1 1 ...
     $ years_to_treatment: num  -6 -5 -4 -3 -2 -1 0 1 2 3 ...
     $ post              : num  0 0 0 0 0 0 1 1 1 1 ...
     $ outcome           : num  -4.68 3.52 2.48 8.83 -2.54 ...
    ```

=== "Stata"

    ```stata
    dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear
    describe
    ```

    Output:

    ```
    Generated balanced panel:
      Units:          500
      Time periods:   20
      Total obs:      10000
      Treated units:  195 (39%)
      Treatment time: uniformly in {7, 8, 9}
      Treatment effect: -3 (post-treatment)

    Contains data
     Observations:        10,000
        Variables:             7
    ---------------------------------------------------------------
    Variable      Storage   Display    Value
        name         type    format    label      Variable label
    ---------------------------------------------------------------
    unit            long    %12.0g                Unit identifier
    time            long    %12.0g                Time period
    treat           byte    %8.0g                 Treatment group indicator
    treatment_time  long    %12.0g                Period when treatment begins
    years_to_trea~t long    %12.0g                Periods relative to treatment
                                                    (-1000 = never treated)
    post            byte    %8.0g                 Post-treatment indicator
    outcome         double  %10.0g                Outcome (noise + treatment
                                                    effect)
    ---------------------------------------------------------------
    Sorted by: unit  time
    ```

!!! note
    R and Stata use different random number generators, so `seed = 42` produces different data in each language. The coefficients differ numerically but the DLM works identically in both.

## Estimate the DLM

=== "R"

    ```r
    library(dlm)
    library(dplyr)

    df <- generate_data(seed = 42, n_groups = 500, n_times = 20, treat_prob = 0.4)

    outcome_data <- df %>% select(group, time, outcome)
    exposure_data <- df %>% select(group, time, post) %>% distinct()

    mod <- distributed_lags_model(
      data = outcome_data,
      exposure_data = exposure_data,
      from_rt = -3, to_rt = 3,
      outcome = "outcome", exposure = "post",
      unit = "group", time = "time"
    )

    mod$betas
    ```

    Output:

    ```
               time_to_event        coef        se
    post_lead2            -3 -0.11846256 0.4729260
    post_lead1            -2 -0.06327287 0.5139965
    post_lag0              0 -2.64104056 0.5438700
    post_lag1              1 -2.26029265 0.5234023
    post_lag2              2 -3.04210098 0.5674556
    post_lag3              3 -2.61751913 0.4214689
    ```

=== "Stata"

    ```stata
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

## Interpret Results

The output table shows **beta coefficients** — cumulative treatment effects at each event time relative to the reference period (default: −1).

- **Pre-treatment betas** (t < ref): Should be near zero if parallel trends hold. These test the "no pre-trend" assumption.
- **Post-treatment betas** (t ≥ 0): Estimate the cumulative causal effect at each horizon.
- **Reference period** (t = −1): Normalized to zero by construction.

In both examples, the true treatment effect is −3. The estimated post-treatment betas cluster around −3, while pre-treatment betas are close to zero — exactly as expected.

Because the test data uses a binary absorbing treatment, these betas are numerically identical to what a canonical binned-endpoint event study would produce on the same data. Both the R and Stata packages include equivalence tests that verify this match to machine precision (~10⁻¹³).

## Plot the Results

=== "R"

    The returned model object includes a ready-made event-study plot:

    ```r
    mod$plot
    ```

    ![R Event-Study Plot](../assets/plot_quickstart_r_v3.png){ width="600" }

=== "Stata"

    After estimation, use `dlm_plot` for a publication-quality event study plot:

    ```stata
    dlm_plot, title("Event-Study Plot (DLM)") xtitle("Periods to Treatment") ytitle("Coefficient")
    ```

    ![Stata Event-Study Plot](../assets/plot_quickstart_stata_v3.png){ width="600" }

## Access Results Programmatically

=== "R"

    ```r
    # Beta coefficients (data.frame: time_to_event, coef, se)
    mod$betas
    #            time_to_event        coef        se
    # post_lead2            -3 -0.11846256 0.4729260
    # post_lead1            -2 -0.06327287 0.5139965
    # post_lag0              0 -2.64104056 0.5438700
    # post_lag1              1 -2.26029265 0.5234023
    # post_lag2              2 -3.04210098 0.5674556
    # post_lag3              3 -2.61751913 0.4214689

    # Event-study plot (ggplot2 object)
    mod$plot

    # Underlying fixest model object
    summary(mod$model)

    # Gamma variance-covariance matrix
    mod$vcov

    # Number of observations
    nobs(mod$model)
    # [1] 7500
    ```

=== "Stata"

    ```stata
    * Beta coefficients (time, coef, se, ci_lo, ci_hi)
    matrix list e(betas)
    ```

    ```
    e(betas)[7,5]
        time_to_ev~t          coef            se         ci_lo         ci_hi
    r1            -3    -.06395882     .48354546    -1.0117079     .88379028
    r2            -2     .09466869     .48042378    -.84696192     1.0362993
    r3            -1             0             0             0             0
    r4             0    -2.7639728     .46190658    -3.6693097    -1.8586358
    r5             1    -3.0942824     .52041411    -4.1142941    -2.0742708
    r6             2    -2.7076914     .55493963    -3.7953731    -1.6200097
    r7             3    -3.2569205     .42619968    -4.0922719    -2.4215691
    ```

    ```stata
    * Raw gamma coefficients from the lead/lag regression
    matrix list e(gamma)
    ```

    ```
    e(gamma)[1,6]
        _dlm_lead2  _dlm_lead1   _dlm_lag0   _dlm_lag1   _dlm_lag2   _dlm_lag3
    y1    .1586275  -.09466869  -2.7639728  -.33030967   .38659102  -.54922912
    ```

    ```stata
    * Scalars
    display e(N)          // number of observations
    display e(N_clust)    // number of clusters
    display e(from)       // from period
    display e(to)         // to period
    display e(ref_period) // reference period
    ```

    ```
    7500
    500
    -3
    3
    -1
    ```
