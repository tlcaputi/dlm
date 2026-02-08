# Mathematical Background

This page summarizes the key theoretical results from Schmidheiny and Siegloch (2023), "On event studies and distributed-lag models: Equivalence, generalization and practical implications," *Journal of Applied Econometrics*, 38(5): 695–713.

## Event-Study Regression

The standard event-study specification is:

$$y_{it} = \alpha_i + \lambda_t + \sum_{\substack{k = \underline{k} \\ k \neq k^*}}^{\bar{k}} \beta_k \cdot D_{it}^k + \varepsilon_{it}$$

where:

- $\alpha_i$ and $\lambda_t$ are unit and time fixed effects
- $D_{it}^k$ are event-time indicators (1 if unit $i$ is $k$ periods from treatment at time $t$)
- $k^*$ is the omitted reference period (typically −1)
- $\underline{k}$ and $\bar{k}$ are the endpoints of the event window
- The endpoint dummies are **binned**: $D_{it}^{\underline{k}} = \mathbf{1}[\text{event time} \leq \underline{k}]$ and $D_{it}^{\bar{k}} = \mathbf{1}[\text{event time} \geq \bar{k}]$

## Distributed Lag Model

The DLM regresses the outcome on leads and lags of the treatment indicator:

$$y_{it} = \alpha_i + \lambda_t + \sum_{j=0}^{|\underline{k}|-1} \gamma_j^{\text{lead}} \cdot x_{i,t+j} + \sum_{j=0}^{\bar{k}} \gamma_j^{\text{lag}} \cdot x_{i,t-j} + \varepsilon_{it}$$

where $x_{it}$ is the binary treatment variable (typically: 1 if post-treatment, 0 otherwise).

This produces **gamma** ($\gamma$) coefficients on the individual leads and lags.

## The Equivalence

**Theorem (Schmidheiny & Siegloch 2023):** The gamma coefficients from the DLM and the beta coefficients from the binned event study are linked by cumulative summation:

### Before the reference period

For $k < k^*$:

$$\beta_k = -\sum_{j=k}^{k^*-1} \gamma_j$$

(negative reverse cumulative sum of gammas from $k$ to the period just before the reference)

### After the reference period

For $k > k^*$:

$$\beta_k = \sum_{j=k^*+1}^{k} \gamma_j$$

(forward cumulative sum of gammas from just after the reference to $k$)

## Standard Error Propagation

Because betas are cumulative sums of gammas, the standard errors must account for covariances. For a cumulative sum $\beta = \sum_{j=a}^{b} \gamma_j$:

$$\text{SE}(\beta) = \sqrt{\mathbf{1}' \cdot V_{a:b, a:b} \cdot \mathbf{1}}$$

where $V$ is the variance-covariance matrix of the gammas and $\mathbf{1}$ is a conformable vector of ones. This correctly accounts for the correlation between adjacent gamma estimates.

!!! warning "Not just the sum of SEs"
    A common mistake is computing $\text{SE}(\beta) = \sqrt{\sum \text{SE}(\gamma_j)^2}$, which ignores off-diagonal covariances and produces incorrect confidence intervals.

## Why Use the DLM?

The DLM and event-study regressions produce **identical point estimates and standard errors** — they are algebraically equivalent. So why use the DLM formulation?

1. **Fewer indicator variables.** The DLM uses $|\underline{k}| + \bar{k}$ continuous leads/lags instead of $|\underline{k}| + \bar{k}$ indicator variables. With large event windows, this can be faster.

2. **No multicollinearity from binning.** The event-study approach requires careful construction of binned endpoints; errors in binning cause subtle bugs. The DLM avoids this entirely.

3. **Natural time-window restriction.** Missing leads/lags at panel edges automatically restrict the sample to the correct time window — no manual `if` conditions needed.

4. **Cleaner extension to continuous treatments.** The DLM generalizes naturally to continuous (non-binary) exposure variables, where event-time dummies don't apply.

## Verification

Both the R and Stata implementations include equivalence tests that verify DLM betas match event-study betas to machine precision (~10⁻¹³). The cross-language test confirms the R and Stata implementations agree to ~10⁻⁸ (the small difference is due to different HDFE solvers: `fixest` in R vs. `reghdfe` in Stata).

## Reference

Schmidheiny, K. and S. Siegloch (2023). "On event studies and distributed-lag models: Equivalence, generalization and practical implications." *Journal of Applied Econometrics*, 38(5): 695–713. [DOI: 10.1002/jae.2971](https://doi.org/10.1002/jae.2971)
