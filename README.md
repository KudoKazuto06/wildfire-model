# Stochastic Modeling of Wildfire Intensity Under Changing Environmental Conditions

## STAT 471 — Probability I

University of Alberta  
Fall 2026

### Project Members

- Cong Thinh Nguyen
- Akhiruzzaman Habib

---

## 1. Project Overview

This project investigates wildfire intensity under changing environmental
conditions using stochastic differential equations (SDEs), stochastic
processes, mathematical analysis, and Monte Carlo simulation.

Wildfire intensity is influenced by several environmental factors,
including wind speed, relative humidity, fuel characteristics, and
unobserved sources of randomness.

Rather than treating environmental conditions as deterministic quantities,
this project represents wind speed and relative humidity as stochastic
processes that evolve continuously over time.

These environmental processes are coupled with a stochastic wildfire
intensity model to investigate:

- The evolution of wildfire intensity.
- The probability distribution of future intensity.
- The probability of reaching dangerous intensity thresholds.
- The expected time required to reach these thresholds.
- Sensitivity to environmental and stochastic parameters.
- Differences between seasonal environmental conditions.
- The consequences of environmental dependence.
- The importance of nonlinear wind-humidity interactions.

The project combines theoretical probability with numerical simulation
to investigate the mathematical properties of the proposed model.

This is a mathematical modeling and simulation study, rather than a
calibrated operational wildfire forecasting system.

---

## 2. Research Questions

| ID | Research Question | Status |
|---|---|---|
| Q1 | How can wildfire intensity and environmental conditions be represented as a coupled stochastic system? | Completed |
| Q2 | What is the probability distribution of wildfire intensity at a future time? | Completed |
| Q3 | What is the probability of reaching a dangerous intensity threshold, and when does this occur? | Completed |
| Q4 | How sensitive is wildfire intensity to changes in model parameters? | Completed |
| Q5 | How does wildfire intensity differ under seasonal environmental conditions? | Completed |
| Q6 | How do environmental dependence and alternative stochastic assumptions affect wildfire intensity? | Completed |

Detailed mathematical derivations, numerical experiments, and
interpretations are provided in the `research_questions/` directory.

---

## 3. Project Structure

    stat471-wildfire-model/
    |
    |-- data/
    |   |-- raw/
    |   `-- processed/
    |
    |-- src/
    |   |-- __init__.py
    |   |-- config.py
    |   |-- environment.py
    |   |-- intensity.py
    |   |-- simulation.py
    |   |-- theoretical.py
    |   |-- monte_carlo.py
    |   |-- stopping_times.py
    |   |-- sensitivity.py
    |   |-- seasonal.py
    |   `-- dependence.py
    |
    |-- experiments/
    |   |-- __init__.py
    |   |-- baseline.py
    |   |-- validate_ou.py
    |   |-- wildfire_baseline.py
    |   |-- intensity_distribution.py
    |   |-- threshold_crossing.py
    |   |-- sensitivity_analysis.py
    |   |-- plot_sensitivity.py
    |   |-- seasonal_comparison.py
    |   |-- plot_seasonal.py
    |   |-- validate_dependence_engine.py
    |   |-- validate_correlated_ou.py
    |   |-- dependence_analysis.py
    |   `-- plot_dependence.py
    |
    |-- tests/
    |   |-- test_environment.py
    |   |-- test_stopping_times.py
    |   |-- test_seasonal.py
    |   `-- test_dependence.py
    |
    |-- research_questions/
    |   |-- Q1_coupled_system.md
    |   |-- Q2_intensity_distribution.md
    |   |-- Q3_threshold_crossing.md
    |   |-- Q4_sensitivity_analysis.md
    |   |-- Q5_seasonal_comparison.md
    |   `-- Q6_dependence_assumptions.md
    |
    |-- results/
    |   |-- figures/
    |   |-- CSV experiment results
    |   `-- NPZ simulation samples
    |
    |-- requirements.txt
    `-- README.md

The `src/` directory contains reusable mathematical and numerical
implementations.

The `experiments/` directory contains reproducible numerical experiments
and visualization scripts.

The `tests/` directory contains automated mathematical and numerical
consistency tests.

The `research_questions/` directory contains detailed explanations of
each research question.

The `results/` directory contains generated figures, numerical summaries,
and saved simulation samples.

---

## 4. Historical Wildfire Data

### 4.1 Data Source

The empirical component of this project uses Alberta historical wildfire
records covering the period 2006–2025.

The dataset provides information about wildfire events, including:

- Wind speed.
- Relative humidity.
- Temperature.
- Recorded fire spread rate.
- Fuel classification.
- Seasonal occurrence.

### 4.2 Data Cleaning

| Description | Number of Records |
|---|---:|
| Original wildfire observations | 27,828 |
| Complete observations retained | 24,859 |
| Observations excluded | 2,969 |
| Retention rate | 89.33% |

Complete observations are used to construct environmental normalization
constants and seasonal environmental scenarios.

### 4.3 Environmental Summary Statistics

| Variable | Mean | Standard Deviation |
|---|---:|---:|
| Wind speed (km/h) | 8.920230 | 8.415991 |
| Relative humidity (%) | 45.203347 | 18.674584 |
| Temperature (°C) | 17.88 | 7.78 |
| Fire spread rate (m/min) | 0.91 | 2.61 |

Wind speed and relative humidity are the primary environmental variables
included in the stochastic model.

Temperature is not explicitly included in the current model.

### 4.4 Fuel Classification

The model considers a fixed fuel category:

    C2 — Boreal Spruce

C2 was selected because it is the most frequently recorded fuel
classification in the historical dataset.

Fixing the fuel category allows the project to focus on environmental
and stochastic effects without introducing additional fuel-dependent
parameters.

### 4.5 Data Limitations

The historical dataset consists primarily of wildfire event-level
observations rather than continuous longitudinal trajectories.

Consequently:

- Environmental normalization constants are data-derived.
- Seasonal initial conditions are data-derived.
- OU mean-reversion and diffusion parameters are assumed.
- Wildfire intensity growth coefficients are assumed.
- Brownian correlations are investigated through scenarios rather than
  directly estimated from the dataset.

The dataset does not provide all quantities required to reconstruct
continuous head fire intensity trajectories.

The assumed stochastic parameters should therefore not be interpreted
as statistically calibrated estimates of real wildfire dynamics.

---

## 5. Mathematical Model

### 5.1 Environmental Standardization

Let:

- W(t) denote wind speed.
- H(t) denote relative humidity.
- I(t) denote wildfire intensity.

The standardized environmental variables are:

    Z_W(t) = (W(t) - mu_W) / s_W

    Z_H(t) = (H(t) - mu_H) / s_H

where:

    mu_W = 8.920230
    s_W  = 8.415991

    mu_H = 45.203347
    s_H  = 18.674584

Positive standardized values represent conditions above the historical
sample mean, while negative values represent conditions below it.

---

### 5.2 Ornstein-Uhlenbeck Environmental Processes

Wind speed and relative humidity are modeled using Ornstein-Uhlenbeck
(OU) processes:

    dZ_W(t) = -kappa_W Z_W(t) dt + sigma_W dB_W(t)

    dZ_H(t) = -kappa_H Z_H(t) dt + sigma_H dB_H(t)

The OU process is selected because it provides:

- Continuous stochastic evolution.
- Mean-reverting behavior.
- Gaussian marginal distributions.
- Explicit analytical solutions.
- Explicit mean and variance.
- Exact numerical transitions.

For a generic OU process:

    dZ(t) = -kappa Z(t) dt + sigma dB(t)

the solution is:

    Z(t) = Z(0) exp(-kappa t)
           + sigma integral_0^t exp[-kappa(t-s)] dB(s)

Its expectation is:

    E[Z(t)] = Z(0) exp(-kappa t)

Its variance is:

    Var[Z(t)] =
        sigma^2 / (2 kappa) * [1 - exp(-2 kappa t)]

The stationary distribution is:

    Z(infinity) ~ Normal(0, sigma^2 / (2 kappa))

---

### 5.3 Wildfire Intensity Process

The wildfire intensity process is modeled as:

    dI(t) = I(t) r(Z_W(t), Z_H(t)) dt
            + sigma_I I(t) dB_I(t)

where the environmental growth rate is:

    r(w,h) = alpha_C
             + beta_W w
             - beta_H h
             + gamma w h

The parameters represent:

| Parameter | Interpretation |
|---|---|
| alpha_C | Baseline fuel-dependent growth rate |
| beta_W | Wind contribution to intensity growth |
| beta_H | Humidity contribution to intensity growth |
| gamma | Wind-humidity interaction coefficient |
| sigma_I | Intrinsic intensity volatility |

The interaction term allows the effect of one environmental variable
to depend on the other.

In particular:

    partial r / partial w = beta_W + gamma h

    partial r / partial h = -beta_H + gamma w

The additive model is obtained by setting:

    gamma = 0

The interaction model uses:

    gamma = -0.04

---

### 5.4 Log-Intensity Transformation

Applying Ito's formula to:

    Y(t) = log(I(t))

gives:

    dY(t) =
        [r(Z_W(t), Z_H(t)) - sigma_I^2 / 2] dt
        + sigma_I dB_I(t)

Therefore:

    Y(t) = log(I(0))
           + integral_0^t r(Z_W(s), Z_H(s)) ds
           - sigma_I^2 t / 2
           + sigma_I B_I(t)

The intensity process can be recovered through:

    I(t) = exp(Y(t))

This representation guarantees strictly positive intensity for the
continuous-time model.

It also motivates the positivity-preserving numerical scheme used
throughout the project.

---

### 5.5 Markov Property

Define the joint stochastic state:

    X(t) = (I(t), Z_W(t), Z_H(t))

The joint process is Markov because its future evolution depends on
the current state and future Brownian increments.

The intensity process alone is generally not Markov because its growth
rate depends on the environmental state.

The joint representation therefore provides the appropriate state space
for mathematical analysis.

---

### 5.6 Infinitesimal Generator

Under independent Brownian drivers, the generator of the joint process
is:

    Lf =
        i r(w,h) f_i
        - kappa_W w f_w
        - kappa_H h f_h
        + (sigma_I^2 i^2 / 2) f_ii
        + (sigma_W^2 / 2) f_ww
        + (sigma_H^2 / 2) f_hh

where the subscripts denote partial derivatives.

Under correlated environmental Brownian motions, an additional
mixed-derivative term appears:

    rho sigma_W sigma_H f_wh

The generator connects the stochastic model with martingale methods
and Kolmogorov equations.

---

### 5.7 Martingale Representation

For sufficiently regular test functions f, the process:

    M_f(t) =
        f(X(t)) - f(X(0))
        - integral_0^t Lf(X(s)) ds

is a local martingale under the appropriate regularity conditions.

This provides a theoretical connection between the stochastic
differential equations and the probability concepts studied in STAT 471.

---

## 6. Model Parameters

### 6.1 Environmental Parameters

| Parameter | Baseline Value |
|---|---:|
| kappa_W | 0.8 |
| sigma_W | 1.0 |
| kappa_H | 0.5 |
| sigma_H | 0.7 |

### 6.2 Wildfire Intensity Parameters

| Parameter | Baseline Value |
|---|---:|
| alpha_C | 0.05 |
| beta_W | 0.12 |
| beta_H | 0.10 |
| gamma | -0.04 |
| sigma_I | 0.15 |
| Initial intensity | 500 kW/m |

### 6.3 Initial Environmental Conditions

The baseline environmental scenario is:

    Z_W(0) = 1.5

    Z_H(0) = -1.0

This represents relatively windy and dry initial conditions.

Seasonal experiments use alternative initial environmental conditions
derived from historical seasonal averages.

### 6.4 Simulation Parameters

| Parameter | Value |
|---|---:|
| Simulation horizon | 10 |
| Time step | 0.01 |
| Number of time steps | 1,000 |
| Number of observation times | 1,001 |
| Monte Carlo paths per main configuration | 10,000 |
| Random seed | 42 |
| Default Brownian correlation | 0 |

The simulation time scale is illustrative and has not been calibrated
to a specific physical time unit.

---

## 7. Numerical Methods

### 7.1 Exact OU Transitions

The environmental processes are simulated using exact OU transitions.

For a time step dt:

    Z(n+1) =
        exp(-kappa dt) Z(n)
        + sigma sqrt(
            [1 - exp(-2 kappa dt)] / (2 kappa)
          ) epsilon(n)

where:

    epsilon(n) ~ Normal(0,1)

Unlike the Euler-Maruyama approximation, this transition reproduces
the exact conditional Gaussian distribution of the OU process at
the numerical observation times.

---

### 7.2 Log-Euler Intensity Simulation

The wildfire intensity process is simulated in logarithmic space:

    Y(n+1) =
        Y(n)
        + [r(Z_W(n), Z_H(n)) - sigma_I^2 / 2] dt
        + sigma_I sqrt(dt) epsilon_I(n)

Physical intensity is recovered using:

    I(n+1) = exp(Y(n+1))

This numerical scheme preserves positivity and avoids negative
intensity values that can arise from direct Euler discretization.

The environmental transitions are exact, but the intensity drift
integral is numerically approximated.

---

### 7.3 Vectorized Monte Carlo Simulation

The Monte Carlo engine simulates multiple trajectories simultaneously
using NumPy arrays.

For the baseline configuration, the intensity trajectory matrix has
dimensions:

    (10000, 1001)

The vectorized implementation supports:

- Large-scale simulation.
- Parameter overrides.
- Alternative initial conditions.
- Additive and interaction models.
- Correlated environmental processes.
- Reproducible random number generation.

---

### 7.4 Common Random Numbers

Common random numbers are used when comparing model configurations.

The same underlying random streams are reused across paired
simulations.

For example, additive and interaction models are simulated using
the same environmental trajectories and intensity innovations.

This allows differences between paired trajectories to be attributed
to changes in model assumptions and reduces the sampling variability
of estimated differences.

Trajectories are independent across Monte Carlo sample indices,
while corresponding trajectories across configurations are deliberately
paired.

---

### 7.5 Monte Carlo Confidence Intervals

For simulated observations X_1, ..., X_N, the sample mean is:

    X_bar = (1/N) sum_i X_i

The estimated standard error is:

    SE = s / sqrt(N)

An approximate pointwise 95% confidence interval is:

    X_bar +/- 1.96 SE

For paired model comparisons, confidence intervals are calculated
from the sample of paired differences.

For threshold probabilities, Bernoulli-based standard errors are used.

Normal-approximation intervals should be interpreted cautiously for
very rare events and small crossing counts.

These intervals quantify Monte Carlo sampling uncertainty, not
uncertainty in the assumed physical model parameters.

---

## 8. Q1 — Coupled Stochastic System

### Objective

Construct and validate a coupled stochastic model for environmental
conditions and wildfire intensity.

### Implementation

- `src/environment.py`
- `src/intensity.py`
- `src/simulation.py`
- `experiments/baseline.py`
- `experiments/validate_ou.py`
- `experiments/wildfire_baseline.py`

### Mathematical Validation

The numerical OU simulations are compared against their exact
theoretical means and variances.

At T = 10:

| Quantity | Theoretical | Monte Carlo |
|---|---:|---:|
| Wind mean | 0.00050 | 0.02007 |
| Wind variance | 0.62500 | 0.61631 |
| Humidity mean | -0.00674 | -0.02274 |
| Humidity variance | 0.48998 | 0.49803 |

The empirical results agree closely with the analytical OU moments.

### Main Findings

- Wind and humidity exhibit mean-reverting stochastic behavior.
- Environmental conditions influence wildfire intensity growth.
- The log-Euler method preserves positive intensity.
- The joint process provides a Markov representation.
- The simulated environmental moments agree with OU theory.

### Generated Figures

- Baseline environmental trajectories.
- Wind OU validation.
- Humidity OU validation.
- Coupled wildfire intensity trajectory.

### Detailed Documentation

See:

    research_questions/Q1_coupled_system.md

---

## 9. Q2 — Terminal Wildfire Intensity Distribution

### Objective

Investigate the probability distribution of wildfire intensity
at a future time T.

### Analytical Additive Model

When gamma = 0, define:

    J_W(T) = integral_0^T Z_W(t) dt

    J_H(T) = integral_0^T Z_H(t) dt

The terminal log-intensity is:

    log(I(T)) =
        log(I(0))
        + (alpha_C - sigma_I^2 / 2) T
        + beta_W J_W(T)
        - beta_H J_H(T)
        + sigma_I B_I(T)

Under deterministic initial conditions and independent Brownian
drivers, this is a Gaussian random variable.

Therefore:

    log(I(T)) ~ Normal(m, v)

and:

    I(T) ~ Lognormal(m, v)

The analytical moments are:

    E[I(T)] = exp(m + v/2)

    Var[I(T)] = [exp(v) - 1] exp(2m + v)

    Median[I(T)] = exp(m)

### Numerical Results

| Statistic | Additive Theory | Additive MC | Interaction MC |
|---|---:|---:|---:|
| Mean | 1478.03 | 1485.39 | 1596.97 |
| Standard deviation | 1258.99 | 1222.71 | 1438.32 |
| Median | 1125.17 | 1132.66 | 1166.65 |
| 2.5th percentile | 264.55 | 264.03 | 277.32 |
| 97.5th percentile | 4785.53 | 4727.54 | 5434.03 |
| Log mean | 7.02569 | 7.03215 | 7.07871 |
| Log variance | 0.54556 | 0.55024 | 0.58505 |

### Paired Interaction Effect

    Mean difference = 111.57796 kW/m

    Standard error = 2.93821 kW/m

    95% confidence interval:
        [105.81906, 117.33686]

    Percentage increase = 7.51%

    Fraction of trajectories with increased intensity = 72.95%

### Main Findings

- The additive model has an exact continuous-time lognormal
  terminal distribution.
- Monte Carlo estimates agree closely with the analytical distribution.
- The interaction model generally does not retain exact lognormality.
- Under the baseline conditions, the interaction increases mean
  intensity and upper quantiles.
- The interaction increases empirical intensity variance by
  approximately 38.38%.

### Implementation

- `src/theoretical.py`
- `src/monte_carlo.py`
- `experiments/intensity_distribution.py`

### Generated Figures

- Terminal intensity distributions.
- Log-intensity distributions.
- Normal Q-Q plots.
- Empirical cumulative distribution functions.

### Detailed Documentation

See:

    research_questions/Q2_intensity_distribution.md

---

## 10. Q3 — Threshold Crossing and Stopping Times

### Objective

Investigate the probability that wildfire intensity reaches a
specified threshold before the simulation horizon.

### First Passage Time

For an intensity threshold K:

    tau_K = inf {t >= 0 : I(t) >= K}

Define the running maximum:

    M(T) = sup_{0 <= t <= T} I(t)

Then:

    P(tau_K <= T) = P(M(T) >= K)

This probability differs from:

    P(I(T) >= K)

because intensity may cross a threshold and subsequently fall below it.

### Numerical Thresholds

    K = 1000, 2000, 3000, 4000, 5000 kW/m

These thresholds are illustrative numerical risk levels.

### Crossing Probability Results

| Threshold | Additive | Interaction | Difference (pp) |
|---|---:|---:|---:|
| 1000 | 73.77% | 75.59% | +1.82 |
| 2000 | 29.24% | 31.97% | +2.73 |
| 3000 | 12.79% | 15.41% | +2.62 |
| 4000 | 5.60% | 7.80% | +2.20 |
| 5000 | 2.72% | 4.08% | +1.36 |

### Conditional Mean First Passage Times

| Threshold | Additive | Interaction |
|---|---:|---:|
| 1000 | 3.89479 | 3.62587 |
| 2000 | 6.74927 | 6.35238 |
| 3000 | 7.71501 | 7.32248 |
| 4000 | 8.13120 | 7.83714 |
| 5000 | 8.50261 | 8.10963 |

First passage time statistics are calculated conditionally on
trajectories that cross the threshold within the simulation horizon.

### Main Findings

- Higher intensity thresholds have lower crossing probabilities.
- Terminal exceedance probabilities underestimate pathwise crossing
  probabilities.
- The interaction increases crossing probabilities under the
  baseline conditions.
- At K = 5000, the interaction increases crossing probability
  by 50% relative to the additive model.
- Conditional passage times are lower in the interaction model
  under the baseline configuration.

Conditional passage-time comparisons do not imply that every
individual trajectory crosses earlier.

### Implementation

- `src/stopping_times.py`
- `experiments/threshold_crossing.py`
- `tests/test_stopping_times.py`

### Generated Figures

- Threshold crossing probabilities.
- Terminal versus pathwise exceedance probabilities.
- Cumulative crossing probabilities.
- First passage time distributions.
- Paired interaction effects on crossing risk.

### Detailed Documentation

See:

    research_questions/Q3_threshold_crossing.md

---

## 11. Q4 — Parameter Sensitivity Analysis

### Objective

Determine how changes in environmental and intensity parameters
affect wildfire behavior.

### Method

One-at-a-time sensitivity analysis is performed across 11 parameters:

    alpha_C
    beta_W
    beta_H
    gamma
    sigma_I
    kappa_W
    kappa_H
    sigma_W
    sigma_H
    Z_W(0)
    Z_H(0)

Five values are considered for each parameter.

The experiment contains 55 scenario entries, with repeated baseline
configurations reused where applicable.

The primary response variables include:

- Mean terminal intensity.
- Standard deviation.
- Median.
- Upper intensity quantiles.
- Threshold crossing probability.

### Analytical Sensitivity

Selected derivatives of expected log-intensity are derived
analytically and compared with paired numerical estimates.

| Parameter | Analytical | Numerical |
|---|---:|---:|
| alpha_C | 10.00000 | 10.00000 |
| beta_W | 1.87437 | 1.83324 |
| beta_H | 1.98652 | 2.07059 |
| gamma | -1.15384 | -1.16399 |
| sigma_I | -1.50000 | -1.48006 |

The estimates agree closely overall, although the continuous-time
beta_H reference lies slightly outside its approximate numerical
confidence interval.

### Sensitivity of Mean Intensity

| Parameter | Mean Range (kW/m) |
|---|---:|
| kappa_W | 1384.02–2557.57 |
| kappa_H | 1416.41–2223.61 |
| alpha_C | 1243.72–2050.55 |
| beta_W | 1326.86–2008.66 |
| beta_H | 1359.48–1940.46 |
| Initial wind | 915.15–1596.97 |
| Initial humidity | 959.77–1596.97 |

### Main Findings

- Mean-reversion rates substantially affect terminal intensity.
- Initial environmental conditions strongly influence crossing risk.
- Wind and humidity coefficients affect wildfire growth.
- Environmental volatility increases uncertainty and can increase
  upper intensity quantiles.
- Intrinsic intensity volatility has little effect on the estimated
  arithmetic mean but substantially affects dispersion and risk.
- Parameter importance depends on the response variable and the
  selected parameter ranges.

The sensitivity rankings are local to the chosen scenarios and should
not be interpreted as global causal importance rankings.

### Implementation

- `src/sensitivity.py`
- `src/theoretical.py`
- `src/monte_carlo.py`
- `experiments/sensitivity_analysis.py`
- `experiments/plot_sensitivity.py`

### Generated Figures

- Mean intensity sensitivity panels.
- Crossing probability sensitivity panels.
- Mean intensity tornado diagram.
- Crossing probability tornado diagram.
- Multi-metric sensitivity heatmap.

### Detailed Documentation

See:

    research_questions/Q4_sensitivity_analysis.md

---

## 12. Q5 — Seasonal Environmental Comparison

### Objective

Investigate how wildfire intensity changes under representative
seasonal environmental conditions.

### Historical Seasonal Conditions

| Season | Observations | Mean Wind | Mean Humidity |
|---|---:|---:|---:|
| Winter | 273 | 12.176 | 55.121 |
| Spring | 9457 | 9.966 | 38.515 |
| Summer | 12492 | 7.894 | 49.064 |
| Fall | 2637 | 9.692 | 49.872 |

The seasonal environmental averages are transformed into standardized
initial conditions.

All dynamic model parameters are held fixed across seasons.

### Experimental Design

- Four seasonal scenarios.
- Two wildfire intensity models.
- 10,000 trajectories per configuration.
- 80,000 total simulated trajectories.

### Terminal Intensity Results

| Season | Additive Mean | Interaction Mean |
|---|---:|---:|
| Winter | 926.06 | 926.73 |
| Spring | 1062.62 | 1078.97 |
| Summer | 915.00 | 922.89 |
| Fall | 936.81 | 943.93 |

### Interaction Model Crossing Probabilities

| Season | K=1000 | K=2000 | K=5000 |
|---|---:|---:|---:|
| Winter | 41.36% | 11.04% | 0.55% |
| Spring | 50.84% | 15.70% | 1.01% |
| Summer | 40.96% | 10.86% | 0.53% |
| Fall | 42.32% | 11.52% | 0.59% |

### Main Findings

- Spring produces the highest simulated mean terminal intensity.
- Summer produces the lowest simulated mean terminal intensity.
- Spring produces the highest crossing probability at all
  investigated thresholds.
- Seasonal initial conditions alter the effect of the interaction term.
- The interaction decreases crossing probability at K = 1000
  across all four seasons.
- The interaction increases crossing probability at thresholds
  K = 2000 through K = 5000.
- All eight theoretical log-mean comparisons agree with the
  corresponding approximate Monte Carlo confidence intervals.

The seasonal results describe the proposed model under different
initial conditions, rather than a calibrated prediction of actual
seasonal wildfire severity.

### Implementation

- `src/seasonal.py`
- `experiments/seasonal_comparison.py`
- `experiments/plot_seasonal.py`
- `tests/test_seasonal.py`

### Generated Figures

- Seasonal environmental conditions.
- Seasonal terminal intensity statistics.
- Seasonal intensity distributions.
- Seasonal crossing probabilities.
- Seasonal first passage times.
- Seasonal interaction effects.
- Seasonal theoretical validation.

### Detailed Documentation

See:

    research_questions/Q5_seasonal_comparison.md

---

## 13. Q6 — Environmental Dependence

### Objective

Investigate how relaxing the environmental independence assumption
changes wildfire intensity and threshold crossing risk.

### Correlated Brownian Motions

The environmental Brownian motions satisfy:

    d<B_W, B_H>(t) = rho dt

where:

    -1 <= rho <= 1

The intrinsic intensity Brownian motion remains independent of
the environmental Brownian motions.

The original model is recovered when:

    rho = 0

### Environmental Covariance

For deterministic initial conditions:

    Cov(Z_W(t), Z_H(t)) =
        rho sigma_W sigma_H
        / (kappa_W + kappa_H)
        * [1 - exp(-(kappa_W + kappa_H)t)]

The stationary covariance is:

    Cov_infinity =
        rho sigma_W sigma_H
        / (kappa_W + kappa_H)

The stationary environmental correlation is:

    Corr_infinity =
        2 rho sqrt(kappa_W kappa_H)
        / (kappa_W + kappa_H)

Therefore, the correlation between environmental states is not
necessarily equal to the correlation between their Brownian drivers.

### Correlated Additive Model

The additive log-intensity remains Gaussian under correlated
environmental Brownian motions.

Its variance becomes:

    Var(log I(T)) =
        sigma_I^2 T
        + beta_W^2 Var(J_W)
        + beta_H^2 Var(J_H)
        - 2 beta_W beta_H Cov(J_W,J_H)

Negative environmental dependence increases log-intensity variance,
while positive dependence decreases it under the selected parameters.

The expected additive log-intensity remains independent of rho.

### Experimental Design

    rho = -1.00, -0.75, -0.50, -0.25, 0.00,
           0.25, 0.50, 0.75, 1.00

- Nine correlation scenarios.
- Two intensity models.
- 10,000 trajectories per configuration.
- 180,000 total simulated trajectories.

### Terminal Intensity Results

| Rho | Additive Mean | Interaction Mean |
|---|---:|---:|
| -1.00 | 1711.14 | 2355.28 |
| -0.75 | 1655.39 | 2132.14 |
| -0.50 | 1597.50 | 1932.73 |
| -0.25 | 1540.83 | 1755.56 |
| 0.00 | 1485.39 | 1596.97 |
| 0.25 | 1431.01 | 1454.09 |
| 0.50 | 1377.48 | 1324.70 |
| 0.75 | 1324.41 | 1206.88 |
| 1.00 | 1268.33 | 1096.28 |

### Interaction Effects

| Rho | Mean Interaction Effect |
|---|---:|
| -1.00 | +37.64% |
| -0.50 | +20.98% |
| 0.00 | +7.51% |
| 0.50 | -3.83% |
| 1.00 | -13.56% |

The interaction increases mean intensity under negative dependence
but decreases mean intensity under sufficiently positive dependence.

### Critical Correlations

Two distinct sign-change points are identified.

| Quantity | Critical Correlation | Method |
|---|---:|---|
| Expected log-interaction effect | 0.2321 | Analytical |
| Arithmetic mean interaction effect | 0.3261 | Numerical interpolation |

These values differ because the arithmetic mean of intensity
depends on the entire distribution of log-intensity.

### Threshold Crossing Risk

For K = 5000:

| Rho | Additive | Interaction |
|---|---:|---:|
| -1.00 | 6.42% | 12.25% |
| -0.50 | 4.73% | 8.31% |
| 0.00 | 2.72% | 4.08% |
| 0.50 | 1.18% | 1.18% |
| 1.00 | 0.14% | 0.03% |

Strong negative dependence substantially increases simulated
high-intensity crossing risk.

The smallest probabilities are based on very few observed crossings
and therefore have considerable relative sampling uncertainty.

### Main Findings

- Environmental dependence changes wildfire intensity distributions.
- Negative correlation amplifies high-intensity outcomes.
- Positive correlation generally reduces high-intensity risk.
- The additive expected log-intensity is invariant to correlation.
- The interaction expected log-intensity depends linearly on
  correlation.
- The interaction effect can reverse direction.
- Independence can underestimate or overestimate risk depending
  on the actual dependence structure.
- All 18 theoretical log-mean comparisons agree with their
  approximate pointwise Monte Carlo confidence intervals.
- All nine theoretical paired log-interaction comparisons agree
  with their approximate Monte Carlo confidence intervals.

### Implementation

- `src/dependence.py`
- `src/monte_carlo.py`
- `experiments/validate_dependence_engine.py`
- `experiments/validate_correlated_ou.py`
- `experiments/dependence_analysis.py`
- `experiments/plot_dependence.py`
- `tests/test_dependence.py`

### Generated Figures

Environmental validation:

- Correlated OU covariance.
- Correlated OU correlation.
- Integrated OU covariance.
- Marginal OU moments under dependence.

Wildfire dependence analysis:

- Terminal intensity statistics.
- Theoretical and numerical validation.
- Terminal intensity distributions.
- Threshold crossing probabilities.
- Conditional first passage times.
- Paired interaction effects.
- Interaction crossing-risk heatmap.
- Consequences of environmental independence.
- Correlation-dependent interaction sign changes.

### Detailed Documentation

See:

    research_questions/Q6_dependence_assumptions.md

---

## 14. Overall Mathematical and Numerical Validation

The project uses several complementary validation approaches.

| Validation | Purpose |
|---|---|
| Exact OU moments | Validate environmental simulation |
| Additive lognormal distribution | Validate intensity simulation |
| Integrated OU moments | Validate analytical terminal distributions |
| Stopping-time identities | Validate threshold calculations |
| Analytical derivatives | Validate sensitivity calculations |
| Seasonal theoretical expectations | Validate seasonal scenarios |
| Correlated OU covariance | Validate environmental dependence |
| Correlated additive distribution | Validate dependence extension |
| Paired log-interaction expectation | Validate nonlinear interaction effects |
| Backward compatibility | Preserve previous independent simulations |
| Automated tests | Verify mathematical and numerical consistency |

The validation results support the correctness of the numerical
implementation under the proposed stochastic assumptions.

They do not establish predictive accuracy against real wildfire
intensity trajectories.

---

## 15. Automated Testing

The automated test suite includes:

- OU environmental process tests.
- Stopping-time tests.
- Seasonal analysis tests.
- Correlation parameter validation.
- Boundary correlation tests.
- Reproducibility tests.
- Backward compatibility tests.
- Analytical reference identities.
- Numerical agreement with theoretical variance.
- Numerical agreement with expected interaction effects.

### Final Test Results

    50 passed in 6.15s

### Test Summary

| Test Group | Passed |
|---|---:|
| Environmental process | 2 |
| Stopping times | 10 |
| Seasonal analysis | 10 |
| Environmental dependence | 28 |
| Total | 50 |

All implemented automated tests passed.

---

## 16. Generated Figures

| Research Question | Number of Figures |
|---|---:|
| Q1 — Coupled stochastic system | 4 |
| Q2 — Terminal distribution | 4 |
| Q3 — Threshold crossing | 5 |
| Q4 — Sensitivity analysis | 5 |
| Q5 — Seasonal comparison | 7 |
| Q6 — Environmental dependence | 13 |
| Total | 38 |

The figures are stored in:

    results/figures/

Numerical summaries are stored as CSV files.

Selected simulation samples are stored as compressed NumPy NPZ files,
allowing visualization scripts to reproduce figures without repeating
the full Monte Carlo experiments.

---

## 17. Installation

### Requirements

- Python 3.13.5
- NumPy
- Pandas
- SciPy
- Matplotlib
- Pytest

### Create a Virtual Environment

From the project root directory:

    python -m venv .venv

Activate the virtual environment on Windows PowerShell:

    .venv\Scripts\Activate.ps1

On Linux or macOS:

    source .venv/bin/activate

Install the required dependencies:

    python -m pip install -r requirements.txt

Alternatively:

    python -m pip install numpy pandas scipy matplotlib pytest

---

## 18. Reproducing the Experiments

All commands should be executed from the project root directory.

### Q1 — Coupled Stochastic System

    python -m experiments.baseline

    python -m experiments.validate_ou

    python -m experiments.wildfire_baseline

### Q2 — Terminal Intensity Distribution

    python -m experiments.intensity_distribution

### Q3 — Threshold Crossing

    python -m experiments.threshold_crossing

### Q4 — Parameter Sensitivity

    python -m experiments.sensitivity_analysis

    python -m experiments.plot_sensitivity

### Q5 — Seasonal Comparison

    python -m src.seasonal

    python -m experiments.seasonal_comparison

    python -m experiments.plot_seasonal

### Q6 — Environmental Dependence

    python -m src.dependence

    python -m experiments.validate_dependence_engine

    python -m experiments.validate_correlated_ou

    python -m experiments.dependence_analysis

    python -m experiments.plot_dependence

### Run All Automated Tests

    python -m pytest -q

For detailed test output:

    python -m pytest -v

---

## 19. Example Usage

The Monte Carlo engine can be used directly.

    from src.monte_carlo import simulate_many

    times, paths = simulate_many(
        n_paths=10000,
        wind_z0=1.5,
        humidity_z0=-1.0,
        gamma=-0.04,
        rho=-0.5,
        horizon=10.0,
        dt=0.01,
        seed=42,
    )

The returned objects are:

    times.shape = (1001,)

    paths.shape = (10000, 1001)

Terminal intensity samples can be extracted using:

    terminal_intensity = paths[:, -1]

Maximum intensity samples can be extracted using:

    maximum_intensity = paths.max(axis=1)

The estimated probability of crossing K = 2000 is:

    crossing_probability = (
        maximum_intensity >= 2000
    ).mean()

Setting:

    rho = 0

recovers the independent environmental model.

Setting:

    gamma = 0

recovers the additive wildfire intensity model.

---

## 20. Model Assumptions and Limitations

### 20.1 Assumed Dynamic Parameters

The OU and wildfire intensity parameters are illustrative rather
than fully estimated from longitudinal observations.

### 20.2 Gaussian Environmental Processes

OU processes have unbounded Gaussian distributions.

Consequently, transformed environmental values can occasionally
represent physically unrealistic wind speeds or humidity levels.

Bounded or transformed stochastic processes could address this
limitation in future work.

### 20.3 Fixed Fuel Classification

Only one fuel category is explicitly represented.

Changes in fuel type, fuel availability, and fuel consumption are
not modeled dynamically.

### 20.4 Simplified Environmental Drivers

Temperature, precipitation, topography, fuel moisture, and spatial
fire spread are not explicitly included.

### 20.5 Constant Environmental Dependence

The Brownian correlation parameter is assumed constant throughout
each simulation.

Real environmental dependence may vary over time and across
meteorological conditions.

### 20.6 Discrete Threshold Monitoring

Threshold crossings are observed on a numerical time grid.

Crossings occurring between observation times may be missed.

Therefore, simulated grid-based crossing probabilities can
underestimate continuous-time crossing probabilities.

### 20.7 Monte Carlo Sampling Error

All numerical estimates contain finite-sample uncertainty.

Very rare threshold events require substantially larger sample
sizes or specialized rare-event simulation methods.

### 20.8 Seasonal Interpretation

Seasonal scenarios differ through initial environmental conditions.

The model does not estimate separate seasonal stochastic dynamics.

### 20.9 Nonlinear Interaction Distribution

The interaction model generally does not possess a simple
lognormal terminal distribution.

Its distribution is investigated numerically, and analytical
expectations of log-intensity are used for validation.

### 20.10 Predictive Interpretation

The model is intended to investigate stochastic mechanisms and
mathematical relationships.

It should not be used directly for operational wildfire prediction,
emergency response, or real-world risk classification without
additional calibration and validation.

---

## 21. Potential Extensions

Future work could investigate:

- Parameter estimation from longitudinal weather observations.
- Calibration using measured head fire intensity trajectories.
- Time-dependent environmental parameters.
- Time-varying environmental dependence.
- Bounded humidity processes.
- Nonnegative wind speed processes.
- Additional environmental drivers.
- Fuel-dependent stochastic dynamics.
- Spatial wildfire propagation.
- Stochastic control and wildfire suppression.
- Rare-event simulation.
- Numerical convergence analysis.
- Global sensitivity analysis.

---

## 22. Overall Conclusion

This project develops a coupled stochastic framework for investigating
wildfire intensity under changing environmental conditions.

The Ornstein-Uhlenbeck processes provide mathematically tractable
representations of environmental fluctuations, while the multiplicative
intensity SDE represents stochastic wildfire growth.

The additive model permits an analytical lognormal terminal
distribution, providing an important reference for numerical
validation.

Introducing nonlinear wind-humidity interactions changes terminal
intensity distributions and threshold crossing probabilities.

Sensitivity analysis demonstrates that wildfire behavior depends
strongly on environmental initial conditions, mean-reversion rates,
and growth parameters.

Seasonal comparisons demonstrate that different initial environmental
conditions can produce substantially different simulated outcomes.

Environmental dependence analysis further demonstrates that the
independence assumption can materially affect the conclusions of
the model.

In particular, the nonlinear interaction effect can reverse direction
as environmental Brownian correlation changes.

Together, the six research questions demonstrate the importance of
stochastic environmental dynamics, nonlinear interactions, stopping
times, and dependence structures in mathematical wildfire modeling.

The theoretical derivations, numerical validation experiments,
38 generated figures, and 50 passing automated tests provide a
reproducible framework for investigating these relationships.

---

## 23. Project Status

- [x] Historical data exploration and preprocessing
- [x] Environmental normalization
- [x] Mathematical model development
- [x] Ornstein-Uhlenbeck simulation
- [x] Analytical OU validation
- [x] Wildfire intensity simulation
- [x] Vectorized Monte Carlo engine
- [x] Q1 — Coupled stochastic system
- [x] Q2 — Terminal intensity distribution
- [x] Q3 — Threshold crossing and stopping times
- [x] Q4 — Parameter sensitivity analysis
- [x] Q5 — Seasonal environmental comparison
- [x] Q6 — Environmental dependence
- [x] Correlated OU validation
- [x] Analytical and numerical comparison
- [x] Numerical experiment visualization
- [x] Automated testing
- [x] Research question documentation
- [ ] Final repository review
- [ ] Code demonstration video
- [ ] Final project report
- [ ] Individual presentation video

---

## 24. References

- Government of Alberta. Alberta Historical Wildfire Data, 2006–2025.
- STAT 471 — Probability I, University of Alberta, Fall 2026.
- Project mathematical modeling document.
- Project historical wildfire data analysis document.
- NumPy documentation: https://numpy.org/doc/
- SciPy documentation: https://docs.scipy.org/doc/scipy/
- Matplotlib documentation: https://matplotlib.org/stable/
- Pytest documentation: https://docs.pytest.org/

---

Developed as part of STAT 471 — Probability I,
University of Alberta, Fall 2026.