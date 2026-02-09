# Mathematical Background

This page summarizes the key theoretical results from Schmidheiny and Siegloch (2023), "On event studies and distributed-lag models: Equivalence, generalization and practical implications," *Journal of Applied Econometrics*, 38(5): 695–713.

## Intuition: From Event Studies to DLMs

Most researchers are familiar with the canonical event study. This section uses that familiarity as a starting point and shows, step by step, how the distributed lag model is a reparametrization of the same idea — and why that reparametrization unlocks continuous treatments.

### Step 1: The canonical event study (what you already know)

In a canonical event study, you observe a panel of units over time. Some units receive a binary treatment that turns on at a known date and stays on. You want to trace out the treatment effect *dynamically* — how does the outcome change at 2 periods before treatment, 1 period before, the period of treatment, 1 period after, etc.?

To do this, you create **event-time dummy variables**. Each dummy equals 1 if a unit is exactly $k$ periods away from its treatment onset, and 0 otherwise. You pick a reference period (typically $k = -1$) and omit it. Then you regress the outcome on these dummies plus unit and time fixed effects.

**Concrete example.** Suppose Unit A is treated at time 5 and you want a window from $-2$ to $+2$ with reference period $-1$. Your data looks like:

| unit | time | outcome | event\_time | $D_{-2}$ | $D_{0}$ | $D_{1}$ | $D_{2}$ |
|------|------|---------|------------|-----------|----------|----------|----------|
| A    | 3    | 12.3    | −2         | 1         | 0        | 0        | 0        |
| A    | 4    | 11.9    | −1 (ref)   | 0         | 0        | 0        | 0        |
| A    | 5    | 8.7     | 0          | 0         | 1        | 0        | 0        |
| A    | 6    | 9.1     | 1          | 0         | 0        | 1        | 0        |
| A    | 7    | 8.4     | 2          | 0         | 0        | 0        | 1        |

You then estimate:

$$y_{it} = \alpha_i + \lambda_t + \beta_{-2} D_{it}^{-2} + \beta_0 D_{it}^{0} + \beta_1 D_{it}^{1} + \beta_2 D_{it}^{2} + \varepsilon_{it}$$

The $\beta$ coefficients are your event-study estimates — each one measures the treatment effect at a particular time horizon relative to the reference period.

### Step 2: The DLM (same data, different columns)

The DLM takes the **same data** but constructs different regressors. Instead of asking "what event-time period is this unit in?", it asks: "what is this unit's treatment status at various time offsets?"

Specifically, it creates **leads** (future values) and **lags** (current and past values) of the treatment variable. For the same window ($-2$ to $+2$), it creates:

- $x_{t+1}$: treatment value 1 period in the **future** (1 lead, because $|-2| - 1 = 1$)
- $x_{t}$: treatment value **now** (lag 0)
- $x_{t-1}$: treatment value 1 period in the **past** (lag 1)
- $x_{t-2}$: treatment value 2 periods in the **past** (lag 2)

Here is the same data for Unit A (treated at time 5, so `post` switches from 0 to 1 at time 5):

| unit | time | outcome | post | $x_{t+1}$ | $x_{t}$ | $x_{t-1}$ | $x_{t-2}$ |
|------|------|---------|------|-----------|----------|-----------|-----------|
| A    | 3    | 12.3    | 0    | 0         | 0        | 0         | 0         |
| A    | 4    | 11.9    | 0    | 1         | 0        | 0         | 0         |
| A    | 5    | 8.7     | 1    | 1         | 1        | 0         | 0         |
| A    | 6    | 9.1     | 1    | 1         | 1        | 1         | 0         |
| A    | 7    | 8.4     | 1    | 1         | 1        | 1         | 1         |

You then estimate:

$$y_{it} = \alpha_i + \lambda_t + \gamma_1^{\text{lead}} \cdot x_{i,t+1} + \gamma_0^{\text{lag}} \cdot x_{it} + \gamma_1^{\text{lag}} \cdot x_{i,t-1} + \gamma_2^{\text{lag}} \cdot x_{i,t-2} + \varepsilon_{it}$$

The $\gamma$ coefficients measure **incremental** changes in the treatment effect at each time offset.

### Step 3: The equivalence (betas = cumulative sums of gammas)

The key result from Schmidheiny & Siegloch (2023) is that the event-study $\beta$ coefficients are just **cumulative sums** of the DLM $\gamma$ coefficients:

| Event-study $\beta$ | Equals | Intuition |
|---------------------|--------|-----------|
| $\beta_{-2}$        | $-\gamma_1^{\text{lead}}$ | Negative of the single lead gamma |
| $\beta_{-1}$        | 0 (reference) | Normalized to zero |
| $\beta_{0}$         | $\gamma_0^{\text{lag}}$ | Just the contemporaneous gamma |
| $\beta_{1}$         | $\gamma_0^{\text{lag}} + \gamma_1^{\text{lag}}$ | Cumulative sum of two gammas |
| $\beta_{2}$         | $\gamma_0^{\text{lag}} + \gamma_1^{\text{lag}} + \gamma_2^{\text{lag}}$ | Cumulative sum of three gammas |

In other words, **gammas are incremental** (the marginal change at each step) and **betas are cumulative** (the total effect up to that point). The two representations contain exactly the same information — you can always convert between them.

When applied to binary treatment data, the DLM betas are numerically identical to the event-study betas. Both the R and Stata packages verify this to machine precision (~10⁻¹³).

### Step 4: Why this matters — continuous treatments

With a binary treatment, both approaches give identical answers, and event-time dummies are perfectly natural. But now consider a **continuous treatment** — say, a tax rate that changes by different amounts at different times:

| unit | time | outcome | tax\_rate |
|------|------|---------|----------|
| A    | 3    | 12.3    | 0        |
| A    | 4    | 11.9    | 0        |
| A    | 5    | 8.7     | 5.0      |
| A    | 6    | 7.1     | 8.5      |
| A    | 7    | 9.4     | 3.0      |

**The canonical event study breaks down.** There is no single "treatment onset" to anchor event-time dummies. When did treatment start — when the rate went to 5? And what is "2 periods after treatment" when the rate keeps changing? Researchers who want to use event-time dummies here must typically dichotomize the treatment into "big changes" vs. no change, discarding variation and fundamentally altering the analysis (Schmidheiny & Siegloch 2023, Section 4.3).

**The DLM works perfectly.** Just create leads and lags of `tax_rate`:

| unit | time | outcome | tax\_rate | $x_{t+1}$ | $x_t$ | $x_{t-1}$ | $x_{t-2}$ |
|------|------|---------|----------|-----------|--------|-----------|-----------|
| A    | 3    | 12.3    | 0        | 0         | 0      | 0         | 0         |
| A    | 4    | 11.9    | 0        | 5.0       | 0      | 0         | 0         |
| A    | 5    | 8.7     | 5.0      | 8.5       | 5.0    | 0         | 0         |
| A    | 6    | 7.1     | 8.5      | 3.0       | 8.5    | 5.0       | 0         |
| A    | 7    | 9.4     | 3.0      | ...       | 3.0    | 8.5       | 5.0       |

The same regression specification works. The gammas now measure the effect of a **one-unit increase** in the treatment variable at each time offset, and the betas (cumulative sums) trace out a dynamic treatment effect plot — exactly the kind of plot researchers are used to seeing from event studies, but now for a continuous treatment.

This is the core contribution of the DLM framework: it **generalizes the canonical event study to continuous treatments** while being **exactly equivalent** in the binary case.

---

## Formal Specification

The sections below present the formal econometric specification from Schmidheiny & Siegloch (2023).

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
