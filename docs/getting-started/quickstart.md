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
    ```

    Variables created:

    | Variable | Description |
    |---|---|
    | `unit` | Unit identifier (1 to N) |
    | `time` | Time period (1 to T) |
    | `treat` | 1 if unit ever treated, 0 otherwise |
    | `treatment_time` | Period when treatment begins (7, 8, or 9) |
    | `years_to_treatment` | Periods relative to treatment (−1000 if never treated) |
    | `post` | 1 if post-treatment, 0 otherwise |
    | `outcome` | Simulated outcome |

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

In the example above, the true treatment effect is −3. The estimated post-treatment betas (−2.76, −3.09, −2.71, −3.26) cluster around −3, while pre-treatment betas (−0.06, 0.09) are close to zero — exactly as expected.

Because the test data uses a binary absorbing treatment, these betas are numerically identical to what a canonical binned-endpoint event study would produce on the same data. Both the R and Stata packages include equivalence tests that verify this match to machine precision (~10⁻¹³).

## Access Results Programmatically

=== "R"

    ```r
    # Beta coefficients (data.frame: time_to_event, coef, se)
    mod$betas

    # Event-study plot (ggplot2 object)
    mod$plot

    # Underlying fixest model object
    summary(mod$model)

    # Gamma variance-covariance matrix
    mod$vcov

    # Number of observations
    nobs(mod$model)
    ```

=== "Stata"

    ```stata
    * Beta coefficients (time, coef, se, ci_lo, ci_hi)
    matrix list e(betas)

    * Raw gamma coefficients from the lead/lag regression
    matrix list e(gamma)

    * Gamma variance-covariance matrix
    matrix list e(gamma_V)

    * Scalars
    display e(N)          // number of observations
    display e(N_clust)    // number of clusters
    display e(from)       // from period
    display e(to)         // to period
    display e(ref_period) // reference period
    ```
