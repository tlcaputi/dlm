*! version 1.0.0
*! Generate balanced panel test data for the dlm package
*!
*! Syntax:
*!   dlm_gen_data, n_groups(int) n_times(int) treat_prob(real) seed(int) clear

program define dlm_gen_data
    version 15.0

    syntax , ///
        [N_groups(integer 676) ///
         N_times(integer 20) ///
         TREAT_prob(real 0.4) ///
         SEED(integer 12345) ///
         CLEAR]

    // Validate inputs
    if (`n_groups' < 2) {
        display as error "n_groups() must be at least 2"
        exit 198
    }
    if (`n_times' < 5) {
        display as error "n_times() must be at least 5"
        exit 198
    }
    if (`treat_prob' <= 0 | `treat_prob' >= 1) {
        display as error "treat_prob() must be between 0 and 1"
        exit 198
    }

    if ("`clear'" != "") {
        clear
    }
    else if (_N > 0) {
        display as error "data in memory; use clear option"
        exit 4
    }

    // Set seed for reproducibility
    set seed `seed'

    // =========================================================================
    // Step 1: Create balanced panel (1 row per unit-time)
    // =========================================================================
    local total_obs = `n_groups' * `n_times'
    quietly set obs `total_obs'

    // Generate unit IDs: 1..n_groups, each appearing n_times rows
    quietly gen long unit = 1 + floor((_n - 1) / `n_times')

    // Generate time within each unit
    quietly gen long time = mod(_n - 1, `n_times') + 1
    sort unit time

    // =========================================================================
    // Step 2: Assign treatment groups and treatment timing
    // =========================================================================
    // Each unit gets treat=1 with probability treat_prob
    // Treated units get treatment_time uniformly in {7, 8, 9}
    quietly {
        // Generate a uniform random for treatment assignment (one per unit)
        gen double _u_treat = .
        by unit: replace _u_treat = runiform() if _n == 1
        by unit: replace _u_treat = _u_treat[1]

        gen byte treat = (_u_treat < `treat_prob')

        // Treatment time: uniform in {7, 8, 9} for treated units
        gen long treatment_time = .
        gen double _u_time = .
        by unit: replace _u_time = runiform() if _n == 1
        by unit: replace _u_time = _u_time[1]
        replace treatment_time = 7 + floor(_u_time * 3) if treat == 1
        // Ensure treatment_time is in {7, 8, 9}
        replace treatment_time = min(treatment_time, 9) if treat == 1

        // Years to treatment
        gen long years_to_treatment = -1000 if treat == 0
        replace years_to_treatment = time - treatment_time if treat == 1

        // Post indicator
        gen byte post = 0
        replace post = (time >= treatment_time) if treat == 1

        // Outcome: random noise + treatment effect of -3 for post-treatment
        gen double outcome = rnormal(0, 5) + (-3) * (years_to_treatment >= 0 & treat == 1)

        // Clean up temporary variables
        drop _u_treat _u_time
    }

    // =========================================================================
    // Step 3: Sort and label
    // =========================================================================
    sort unit time
    order unit time treat treatment_time years_to_treatment post outcome

    label variable unit "Unit identifier"
    label variable time "Time period"
    label variable treat "Treatment group indicator"
    label variable treatment_time "Period when treatment begins"
    label variable years_to_treatment "Periods relative to treatment (-1000 = never treated)"
    label variable post "Post-treatment indicator"
    label variable outcome "Outcome (noise + treatment effect)"

    display as text "Generated balanced panel:"
    display as text "  Units:          `n_groups'"
    display as text "  Time periods:   `n_times'"
    display as text "  Total obs:      `total_obs'"
    quietly count if treat == 1 & time == 1
    local n_treated = r(N)
    display as text "  Treated units:  `n_treated' (" ///
        string(round(`n_treated' / `n_groups' * 100, 0.1)) "%)"
    display as text "  Treatment time: uniformly in {7, 8, 9}"
    display as text "  Treatment effect: -3 (post-treatment)"

end
