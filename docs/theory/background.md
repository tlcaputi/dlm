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

where $x_{it}$ is the treatment variable. In the binary case, $x_{it}$ is an indicator (1 if post-treatment, 0 otherwise). In the generalized case, $x_{it}$ can be a continuous measure of treatment intensity (e.g., a tax rate).

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

### Generalization to continuous treatments

The canonical event study relies on event-time indicator variables — dummies for "3 periods before treatment," "2 periods before," etc. These indicators only make sense when treatment is a one-time binary switch (e.g., a policy turns on and stays on). But many empirical settings involve treatments that are continuous, vary in intensity, change sign, or occur multiple times per unit. For instance:

- **Multiple tax reforms** of different magnitudes hitting different states at different times (Fuest et al., 2018; Suárez Serrato & Zidar, 2016)
- **Minimum wage changes** that increase by different amounts across jurisdictions (Cengiz et al., 2019)
- **Policy shocks** where the treatment variable is inherently continuous (e.g., benefit duration in weeks, tax rates in percentage points)

In these settings, event-time dummies do not apply. Researchers who want to use a canonical event study must typically dichotomize the continuous treatment into "big changes" vs. no change — discarding treatment variation, potentially biasing estimates, and fundamentally changing the analysis (Schmidheiny & Siegloch 2023, Section 4.3).

The DLM avoids this problem entirely. Because it regresses on leads and lags of the treatment variable directly, the treatment variable $T_{i,t}$ can be continuous. Treatment effects are interpreted as the effect of a one-unit increase, just as in a generalized difference-in-differences model. This allows researchers to produce event-study-style dynamic treatment effect plots for treatments that increase and decrease at different magnitudes over many time periods.

### Equivalence to binned event studies

When applied to the special case of a one-time binary absorbing treatment, the DLM produces point estimates and standard errors that are **numerically identical** to a properly specified event study with binned endpoints. The DLM is a reparametrization of the event study: the gamma coefficients measure *incremental* changes in treatment effects, while the beta coefficients (recovered by cumulative summation) measure *cumulative* treatment effects (Schmidheiny & Siegloch 2023, Remark 6). The DLM is therefore a strict generalization of the canonical event study.

### Data requirements (trade-off)

The DLM is more data hungry than a canonical event study. Because the model includes leads and lags of the treatment variable, researchers must observe the treatment variable over a **wider time window** than the outcome variable. Specifically, for a balanced panel observed from period $t$ to $\bar{t}$ with an event window $[\underline{k}, \bar{k}]$, the treatment variable must be observed from $t - \bar{k}$ to $\bar{t} + |\underline{k}| - 1$ (Schmidheiny & Siegloch 2023, Remark 4). Observations at the edges of the panel where leads or lags cannot be constructed are automatically dropped from the estimation sample.

This means that wider event windows require substantially more data. At some point, the estimation sample becomes too small for precise estimation. Researchers should experiment with different event window lengths and verify that estimates leading up to the endpoints converge, as recommended by Schmidheiny & Siegloch (2023, Remark 3).

## Verification

Both the R and Stata implementations include equivalence tests that verify DLM betas match event-study betas to machine precision (~10⁻¹³). The cross-language test confirms the R and Stata implementations agree to ~10⁻⁸ (the small difference is due to different HDFE solvers: `fixest` in R vs. `reghdfe` in Stata).

## Reference

Schmidheiny, K. and S. Siegloch (2023). "On event studies and distributed-lag models: Equivalence, generalization and practical implications." *Journal of Applied Econometrics*, 38(5): 695–713. [DOI: 10.1002/jae.2971](https://doi.org/10.1002/jae.2971)
