# Q6 — Environmental Dependence and Alternative Stochastic Assumptions

## STAT 471 — Stochastic Wildfire Intensity Modeling

---

## 1. Research Question

How do environmental dependence and alternative stochastic assumptions
affect the distribution of wildfire intensity, threshold crossing
probabilities, and the influence of nonlinear environmental interactions?

The original stochastic model assumes that the Brownian motions
driving wind speed, relative humidity, and wildfire intensity are
mutually independent.

However, environmental conditions may exhibit statistical dependence.

For example, changes in wind speed may occur together with changes
in atmospheric humidity.

This research question investigates the consequences of relaxing
the independence assumption between the environmental Brownian
motions while preserving the remaining model structure.

---

## 2. Objectives

The objectives of this research question are:

1. Extend the independent environmental model to correlated
   Brownian motions.

2. Derive the covariance structure of correlated OU processes.

3. Construct exact correlated Gaussian OU transitions.

4. Derive the distribution of additive log-intensity under
   environmental dependence.

5. Derive the expected nonlinear interaction effect.

6. Investigate terminal intensity distributions across different
   correlation scenarios.

7. Investigate threshold crossing probabilities and conditional
   first passage times.

8. Determine whether environmental dependence can reverse the
   effect of the wind-humidity interaction.

9. Quantify the consequences of assuming environmental independence.

10. Validate the theoretical and numerical implementations.

---

## 3. Original Independence Assumption

The original environmental processes are:

    dZ_W(t) = -kappa_W Z_W(t) dt + sigma_W dB_W(t)

    dZ_H(t) = -kappa_H Z_H(t) dt + sigma_H dB_H(t)

The wildfire intensity process is:

    dI(t) =
        I(t) [
            alpha_C
            + beta_W Z_W(t)
            - beta_H Z_H(t)
            + gamma Z_W(t) Z_H(t)
        ] dt
        + sigma_I I(t) dB_I(t)

The original model assumes:

    B_W, B_H, B_I are mutually independent.

Under deterministic initial environmental conditions, this implies:

    Cov(Z_W(t), Z_H(t)) = 0

The environmental processes remain coupled to wildfire intensity
through the intensity growth rate.

However, the environmental processes themselves are independent.

---

## 4. Correlated Brownian Motion Extension

We introduce a Brownian correlation parameter:

    rho in [-1,1]

such that:

    d<B_W, B_H>(t) = rho dt

Equivalently:

    E[dB_W(t) dB_H(t)] = rho dt

The intrinsic intensity Brownian motion remains independent of
both environmental Brownian motions:

    d<B_W, B_I>(t) = 0

    d<B_H, B_I>(t) = 0

The Brownian correlation matrix is:

        [ 1    rho   0 ]
    R = [ rho  1     0 ]
        [ 0    0     1 ]

The eigenvalues are:

    1 + rho
    1 - rho
    1

Therefore, the matrix is positive semidefinite whenever:

    -1 <= rho <= 1

This ensures a mathematically valid Brownian covariance structure.

---

## 5. Interpretation of Environmental Correlation

### Negative Correlation

When rho < 0, positive wind innovations tend to occur with
negative humidity innovations.

Since wind increases the additive growth rate while humidity
decreases it, these innovations tend to reinforce each other.

This can increase variability in wildfire growth.

### Zero Correlation

When rho = 0, the original independent environmental model
is recovered.

### Positive Correlation

When rho > 0, wind and humidity innovations tend to move in
the same direction.

Because their additive effects on intensity have opposite signs,
their contributions can partially offset each other.

This can reduce variability in wildfire growth.

These interpretations concern stochastic innovations rather than
deterministic environmental trajectories.

---

## 6. Historical Data and Correlation

The historical wildfire dataset exhibits a negative Pearson
correlation between recorded wind speed and relative humidity:

    Historical Pearson correlation = -0.17945

The corresponding Spearman correlation is:

    Historical Spearman correlation = -0.22996

These values describe cross-sectional associations between
recorded environmental levels.

They are not direct estimates of the correlation between
continuous-time Brownian innovations.

Estimating Brownian correlation would require longitudinal
environmental observations and an appropriate stochastic
parameter estimation procedure.

Therefore, rho is investigated through hypothetical scenarios
rather than calibrated directly from the historical correlation.

---

## 7. Analytical Solutions of Correlated OU Processes

The environmental OU solutions are:

    Z_W(t) =
        z_W exp(-kappa_W t)
        + sigma_W integral_0^t
          exp[-kappa_W(t-s)] dB_W(s)

    Z_H(t) =
        z_H exp(-kappa_H t)
        + sigma_H integral_0^t
          exp[-kappa_H(t-s)] dB_H(s)

The marginal expectations remain:

    E[Z_W(t)] = z_W exp(-kappa_W t)

    E[Z_H(t)] = z_H exp(-kappa_H t)

The marginal variances remain:

    Var[Z_W(t)] =
        sigma_W^2 / (2 kappa_W)
        * [1 - exp(-2 kappa_W t)]

    Var[Z_H(t)] =
        sigma_H^2 / (2 kappa_H)
        * [1 - exp(-2 kappa_H t)]

Introducing environmental correlation does not change the
individual OU marginal distributions.

It changes their joint distribution.

---

## 8. Environmental Covariance Derivation

Let:

    C(t) = Cov(Z_W(t), Z_H(t))

Using the explicit OU solutions and the covariance identity
for correlated stochastic integrals:

    E[
        integral f(s) dB_W(s)
        integral g(s) dB_H(s)
    ]
    =
        rho integral f(s) g(s) ds

we obtain:

    C(t) =
        rho sigma_W sigma_H
        integral_0^t
        exp[-(kappa_W + kappa_H)(t-s)] ds

Define:

    a = kappa_W + kappa_H

Then:

    C(t) =
        rho sigma_W sigma_H
        integral_0^t exp[-a(t-s)] ds

Evaluating the integral:

    C(t) =
        rho sigma_W sigma_H / a
        * [1 - exp(-at)]

Therefore:

    Cov(Z_W(t), Z_H(t)) =
        rho sigma_W sigma_H
        / (kappa_W + kappa_H)
        * [1 - exp(-(kappa_W + kappa_H)t)]

This is the exact continuous-time covariance for deterministic
initial environmental conditions.

---

## 9. Alternative Derivation Using Ito's Product Rule

Consider the product:

    Z_W(t) Z_H(t)

Applying Ito's product rule:

    d(Z_W Z_H) =
        Z_H dZ_W
        + Z_W dZ_H
        + d<Z_W,Z_H>

Substituting the environmental SDEs:

    d(Z_W Z_H) =
        -(kappa_W + kappa_H) Z_W Z_H dt
        + sigma_W Z_H dB_W
        + sigma_H Z_W dB_H
        + rho sigma_W sigma_H dt

Taking expectations:

    d/dt E[Z_W(t) Z_H(t)] =
        -(kappa_W + kappa_H)
        E[Z_W(t) Z_H(t)]
        + rho sigma_W sigma_H

Define:

    q(t) = E[Z_W(t) Z_H(t)]

Then:

    q'(t) = -a q(t) + rho sigma_W sigma_H

with:

    q(0) = z_W z_H

Solving the linear ordinary differential equation:

    q(t) =
        z_W z_H exp(-at)
        + rho sigma_W sigma_H / a
          * [1 - exp(-at)]

Since:

    E[Z_W(t)] E[Z_H(t)] =
        z_W z_H exp(-at)

we recover:

    Cov(Z_W(t), Z_H(t)) =
        rho sigma_W sigma_H / a
        * [1 - exp(-at)]

This provides an independent derivation of the covariance formula
using Ito calculus.

---

## 10. Stationary Environmental Covariance

Taking the limit as t approaches infinity:

    Cov_infinity =
        rho sigma_W sigma_H
        / (kappa_W + kappa_H)

The stationary marginal variances are:

    Var_W,infinity = sigma_W^2 / (2 kappa_W)

    Var_H,infinity = sigma_H^2 / (2 kappa_H)

Therefore, the stationary environmental correlation is:

    Corr_infinity =
        Cov_infinity
        / sqrt(Var_W,infinity Var_H,infinity)

After simplification:

    Corr_infinity =
        2 rho sqrt(kappa_W kappa_H)
        / (kappa_W + kappa_H)

assuming both environmental volatilities are positive.

Using the baseline parameters:

    kappa_W = 0.8
    kappa_H = 0.5
    sigma_W = 1.0
    sigma_H = 0.7

we obtain:

    Cov_infinity = 0.53846154 rho

    Corr_infinity = 0.97300851 rho

Consequently, Brownian correlation and stationary environmental
correlation are not identical.

Even when:

    rho = 1

the stationary environmental correlation is approximately:

    Corr_infinity = 0.97301

This occurs because the two OU processes apply different
mean-reversion kernels to the Brownian innovations.

---

## 11. Analytical Environmental Dependence Results

| Rho | Covariance at T=10 | Stationary Correlation |
|---|---:|---:|
| -1.00 | -0.53846 | -0.97301 |
| -0.75 | -0.40385 | -0.72976 |
| -0.50 | -0.26923 | -0.48650 |
| -0.25 | -0.13462 | -0.24325 |
| 0.00 | 0.00000 | 0.00000 |
| 0.25 | 0.13462 | 0.24325 |
| 0.50 | 0.26923 | 0.48650 |
| 0.75 | 0.40385 | 0.72976 |
| 1.00 | 0.53846 | 0.97301 |

The environmental covariance is linear in rho.

At T = 10, the environmental covariance is already close to
its stationary value.

---

## 12. Exact Correlated OU Transitions

For a numerical time step dt, the exact OU transitions are:

    Z_W(n+1) =
        exp(-kappa_W dt) Z_W(n) + eta_W(n)

    Z_H(n+1) =
        exp(-kappa_H dt) Z_H(n) + eta_H(n)

The Gaussian transition innovations satisfy:

    E[eta_W] = E[eta_H] = 0

Their variances are:

    V_W =
        sigma_W^2 / (2 kappa_W)
        * [1 - exp(-2 kappa_W dt)]

    V_H =
        sigma_H^2 / (2 kappa_H)
        * [1 - exp(-2 kappa_H dt)]

Their covariance is:

    C_dt =
        rho sigma_W sigma_H
        / (kappa_W + kappa_H)
        * [1 - exp(-(kappa_W + kappa_H)dt)]

The exact transition covariance matrix is:

           [ V_W    C_dt ]
    Q_dt = [ C_dt   V_H  ]

The correlation between the standardized transition innovations is:

    r_dt = C_dt / sqrt(V_W V_H)

when both transition variances are positive.

This correlation is generally not exactly equal to rho.

Using rho directly as the correlation between exact OU transition
innovations would generally produce an incorrect finite-step
covariance structure.

---

## 13. Generating Correlated Gaussian Innovations

Let:

    epsilon_W ~ Normal(0,1)

    epsilon_H ~ Normal(0,1)

be independent standard Gaussian random variables.

Define:

    epsilon_H,corr =
        r_dt epsilon_W
        + sqrt(1 - r_dt^2) epsilon_H

Then:

    Var(epsilon_H,corr) = 1

and:

    Cov(epsilon_W, epsilon_H,corr) = r_dt

The exact OU transitions can therefore be implemented as:

    Z_W(n+1) =
        exp(-kappa_W dt) Z_W(n)
        + sqrt(V_W) epsilon_W

    Z_H(n+1) =
        exp(-kappa_H dt) Z_H(n)
        + sqrt(V_H) epsilon_H,corr

The intensity innovation remains independent.

If either environmental diffusion coefficient is zero, its
transition noise vanishes and the corresponding correlation
does not affect the joint innovation covariance.

When rho = 0:

    r_dt = 0

and:

    epsilon_H,corr = epsilon_H

Therefore, the original independent numerical implementation
is recovered.

---

## 14. Exact Transition Covariance Example

Using:

    rho = -0.5

    dt = 0.01

the calculated transition covariance matrix is:

        [  0.00992042   -0.00347735 ]
    Q = [ -0.00347735    0.00487558 ]

This matrix is used to generate correlated environmental
innovations with the correct exact OU transition distribution.

---

## 15. Joint Markov Process and Generator

The correlated joint state remains:

    X(t) = (I(t), Z_W(t), Z_H(t))

Introducing correlation does not destroy the Markov property.

The diffusion covariance matrix becomes:

    A(i,w,h) =

        [ sigma_I^2 i^2       0                    0              ]
        [ 0                   sigma_W^2            rho sigma_W sigma_H ]
        [ 0                   rho sigma_W sigma_H  sigma_H^2      ]

The generator is:

    L_rho f =
        i r(w,h) f_i
        - kappa_W w f_w
        - kappa_H h f_h
        + (sigma_I^2 i^2 / 2) f_ii
        + (sigma_W^2 / 2) f_ww
        + (sigma_H^2 / 2) f_hh
        + rho sigma_W sigma_H f_wh

The additional mixed partial derivative represents environmental
dependence.

When rho = 0, the original independent generator is recovered.

---

## 16. Integrated Environmental Processes

Define:

    J_W(T) = integral_0^T Z_W(t) dt

    J_H(T) = integral_0^T Z_H(t) dt

For a generic OU process:

    J(T) =
        z_0 [1 - exp(-kappa T)] / kappa
        + sigma / kappa
          integral_0^T
          [1 - exp(-kappa(T-s))] dB(s)

This follows by substituting the OU solution and applying
stochastic Fubini's theorem.

The expectation is:

    E[J(T)] =
        z_0 [1 - exp(-kappa T)] / kappa

The variance is:

    Var[J(T)] =
        sigma^2 / kappa^2
        * [
            T
            - 2[1 - exp(-kappa T)] / kappa
            + [1 - exp(-2 kappa T)] / (2 kappa)
          ]

These marginal integrated moments are unchanged by rho.

---

## 17. Integrated OU Covariance Derivation

Using the integrated OU representations:

    Cov(J_W(T), J_H(T)) =

        rho sigma_W sigma_H
        / (kappa_W kappa_H)

        * integral_0^T
          [1 - exp(-kappa_W(T-s))]
          [1 - exp(-kappa_H(T-s))] ds

Expanding the product:

    [1 - exp(-kappa_W u)]
    [1 - exp(-kappa_H u)]

    = 1
      - exp(-kappa_W u)
      - exp(-kappa_H u)
      + exp(-(kappa_W + kappa_H)u)

Integrating each term gives:

    Cov(J_W(T), J_H(T)) =

        rho sigma_W sigma_H
        / (kappa_W kappa_H)

        * [
            T
            - [1 - exp(-kappa_W T)] / kappa_W
            - [1 - exp(-kappa_H T)] / kappa_H
            + [1 - exp(-(kappa_W + kappa_H)T)]
              / (kappa_W + kappa_H)
          ]

For the baseline parameters and T = 10:

    Cov(J_W(T), J_H(T)) = approximately 13.18297 rho

### Analytical Results

| Rho | Integrated Covariance |
|---|---:|
| -1.00 | -13.18297 |
| -0.75 | -9.88723 |
| -0.50 | -6.59148 |
| -0.25 | -3.29574 |
| 0.00 | 0.00000 |
| 0.25 | 3.29574 |
| 0.50 | 6.59148 |
| 0.75 | 9.88723 |
| 1.00 | 13.18297 |

This covariance is particularly important because terminal
log-intensity depends on the integrated environmental processes,
rather than only their terminal values.

---

## 18. Correlated Additive Intensity Model

When:

    gamma = 0

the terminal log-intensity is:

    Y(T) =
        log(I_0)
        + (alpha_C - sigma_I^2 / 2) T
        + beta_W J_W(T)
        - beta_H J_H(T)
        + sigma_I B_I(T)

The integrated environmental processes are jointly Gaussian.

Since the intrinsic intensity Brownian motion remains independent,
Y(T) is also Gaussian.

Therefore:

    Y(T) ~ Normal(m, v(rho))

and:

    I(T) ~ Lognormal(m, v(rho))

The lognormal distribution remains exact for the continuous-time
additive model even when the environmental Brownian motions
are correlated.

---

## 19. Expected Additive Log-Intensity

Taking expectations:

    m =
        log(I_0)
        + (alpha_C - sigma_I^2 / 2) T
        + beta_W E[J_W(T)]
        - beta_H E[J_H(T)]

The marginal integrated environmental expectations do not
depend on rho.

Therefore:

    E[log I(T)] is independent of rho.

Using the baseline parameters:

    m = 7.02569

for every investigated correlation scenario.

---

## 20. Additive Log-Intensity Variance

Using the variance formula for a linear combination:

    Var(aX + bY) =
        a^2 Var(X)
        + b^2 Var(Y)
        + 2ab Cov(X,Y)

we obtain:

    v(rho) =
        sigma_I^2 T
        + beta_W^2 Var(J_W)
        + beta_H^2 Var(J_H)
        - 2 beta_W beta_H Cov(J_W,J_H)

The negative sign appears because humidity enters the intensity
growth rate with coefficient -beta_H.

For the baseline parameters:

    v(rho) = approximately 0.54556 - 0.31639 rho

Consequently:

    v(-1) > v(0) > v(1)

Negative correlation increases additive log-intensity variance.

Positive correlation decreases additive log-intensity variance.

---

## 21. Additive Terminal Intensity Moments

Since:

    I(T) ~ Lognormal(m,v(rho))

the arithmetic mean is:

    E[I(T)] = exp(m + v(rho)/2)

The variance is:

    Var[I(T)] =
        [exp(v(rho)) - 1]
        exp(2m + v(rho))

The median is:

    Median[I(T)] = exp(m)

The p-th quantile is:

    Q_p =
        exp(m + sqrt(v(rho)) Phi_inverse(p))

Although the expected log-intensity remains constant, the
arithmetic mean changes with rho because of the exponential
transformation.

The median remains constant because it depends only on m.

### Theoretical Results

| Rho | Log Mean | Log Variance | Mean | Median | Q97.5 |
|---|---:|---:|---:|---:|---:|
| -1.00 | 7.02569 | 0.86195 | 1731.36 | 1125.17 | 6941.97 |
| -0.75 | 7.02569 | 0.78285 | 1664.22 | 1125.17 | 6373.09 |
| -0.50 | 7.02569 | 0.70375 | 1599.69 | 1125.17 | 5824.92 |
| -0.25 | 7.02569 | 0.62466 | 1537.65 | 1125.17 | 5296.21 |
| 0.00 | 7.02569 | 0.54556 | 1478.03 | 1125.17 | 4785.53 |
| 0.25 | 7.02569 | 0.46646 | 1420.71 | 1125.17 | 4291.11 |
| 0.50 | 7.02569 | 0.38736 | 1365.62 | 1125.17 | 3810.58 |
| 0.75 | 7.02569 | 0.30827 | 1312.67 | 1125.17 | 3340.55 |
| 1.00 | 7.02569 | 0.22917 | 1261.77 | 1125.17 | 2875.40 |

This provides an exact analytical benchmark for the correlated
Monte Carlo implementation.

---

## 22. Nonlinear Interaction Model

When gamma is nonzero, the terminal log-intensity becomes:

    Y_gamma(T) =
        Y_0(T)
        + gamma integral_0^T Z_W(t) Z_H(t) dt

Define the integrated interaction:

    J_WH(T) =
        integral_0^T Z_W(t) Z_H(t) dt

Then:

    Y_gamma(T) - Y_0(T) = gamma J_WH(T)

Unlike the additive model, J_WH(T) contains a product of
Gaussian processes.

Therefore, the interaction log-intensity is generally not
Gaussian, and the terminal intensity is generally not exactly
lognormal.

However, its expected log-intensity can still be derived
analytically.

---

## 23. Expected Log-Interaction Effect

Taking expectations:

    Delta_log(rho) =
        E[Y_gamma(T) - Y_0(T)]

    = gamma integral_0^T
      E[Z_W(t) Z_H(t)] dt

Using:

    E[Z_W(t) Z_H(t)] =
        E[Z_W(t)] E[Z_H(t)]
        + Cov(Z_W(t), Z_H(t))

we obtain:

    E[Z_W(t) Z_H(t)] =

        z_W z_H exp(-at)

        + rho sigma_W sigma_H / a
          * [1 - exp(-at)]

where:

    a = kappa_W + kappa_H

Substituting and integrating:

    Delta_log(rho) =

        gamma z_W z_H
        * [1 - exp(-aT)] / a

        + gamma rho sigma_W sigma_H / a
          * [
              T - [1 - exp(-aT)] / a
            ]

The expected interaction effect can therefore be separated into
two components.

### Initial Environmental Effect

    Delta_initial =
        gamma z_W z_H
        * [1 - exp(-aT)] / a

### Environmental Covariance Effect

    Delta_covariance(rho) =
        gamma rho sigma_W sigma_H / a
        * [
            T - [1 - exp(-aT)] / a
          ]

Hence:

    Delta_log(rho) =
        Delta_initial + Delta_covariance(rho)

The expected log-interaction effect is linear in rho.

---

## 24. Baseline Expected Interaction Effect

Using:

    z_W = 1.5
    z_H = -1.0
    gamma = -0.04
    kappa_W = 0.8
    kappa_H = 0.5
    sigma_W = 1.0
    sigma_H = 0.7
    T = 10

we obtain:

    Delta_initial = approximately 0.04615

    Delta_covariance(rho) =
        approximately -0.19881 rho

Therefore:

    Delta_log(rho) =
        approximately 0.04615 - 0.19881 rho

### Analytical Results

| Rho | Initial Effect | Covariance Effect | Total Effect |
|---|---:|---:|---:|
| -1.00 | 0.04615 | 0.19882 | 0.24497 |
| -0.75 | 0.04615 | 0.14911 | 0.19527 |
| -0.50 | 0.04615 | 0.09941 | 0.14556 |
| -0.25 | 0.04615 | 0.04970 | 0.09586 |
| 0.00 | 0.04615 | 0.00000 | 0.04615 |
| 0.25 | 0.04615 | -0.04970 | -0.00355 |
| 0.50 | 0.04615 | -0.09941 | -0.05325 |
| 0.75 | 0.04615 | -0.14911 | -0.10296 |
| 1.00 | 0.04615 | -0.19882 | -0.15266 |

Negative dependence reinforces the positive initial interaction
effect.

Sufficiently positive dependence offsets and eventually reverses
the initial interaction effect.

---

## 25. Analytical Critical Correlation

The critical correlation is obtained by solving:

    Delta_log(rho_star) = 0

Using:

    A = [1 - exp(-aT)] / a

we have:

    gamma [
        z_W z_H A
        + rho_star sigma_W sigma_H / a * (T-A)
    ] = 0

For nonzero gamma:

    rho_star =
        -a z_W z_H A
        / [sigma_W sigma_H (T-A)]

Equivalently:

    rho_star =
        -z_W z_H [1 - exp(-aT)]
        / [
            sigma_W sigma_H
            * (T - [1 - exp(-aT)] / a)
          ]

For the baseline configuration:

    rho_star = approximately 0.2321

Therefore:

    rho < 0.2321:
        Expected log-interaction effect is positive.

    rho = 0.2321:
        Expected log-interaction effect is zero.

    rho > 0.2321:
        Expected log-interaction effect is negative.

This is an exact analytical sign-change condition for the
continuous-time expected log-interaction effect.

The numerical value depends on the initial conditions,
environmental parameters, and simulation horizon.

It is not a universal physical wildfire correlation threshold.

---

## 26. Expected Log-Intensity and Geometric Mean

The expected interaction log-intensity is:

    E[Y_gamma(T)] = m + Delta_log(rho)

The geometric mean is:

    G_gamma(rho) = exp(E[Y_gamma(T)])

The ratio between interaction and additive geometric means is:

    G_gamma(rho) / G_0(rho) =
        exp(Delta_log(rho))

### Analytical Results

| Rho | Interaction Log Mean | Geometric Mean |
|---|---:|---:|
| -1.00 | 7.27066 | 1437.49 |
| -0.75 | 7.22095 | 1367.79 |
| -0.50 | 7.17125 | 1301.47 |
| -0.25 | 7.12154 | 1238.36 |
| 0.00 | 7.07184 | 1178.31 |
| 0.25 | 7.02213 | 1121.18 |
| 0.50 | 6.97243 | 1066.81 |
| 0.75 | 6.92273 | 1015.08 |
| 1.00 | 6.87302 | 965.86 |

For the interaction model, the geometric mean should not
automatically be interpreted as the median.

---

## 27. Experimental Design

The following correlation scenarios are considered:

    rho = -1.00
    rho = -0.75
    rho = -0.50
    rho = -0.25
    rho =  0.00
    rho =  0.25
    rho =  0.50
    rho =  0.75
    rho =  1.00

For each correlation, two models are simulated:

    Additive model:
        gamma = 0

    Interaction model:
        gamma = -0.04

### Simulation Settings

| Setting | Value |
|---|---:|
| Correlation scenarios | 9 |
| Models per scenario | 2 |
| Total configurations | 18 |
| Trajectories per configuration | 10,000 |
| Total trajectories | 180,000 |
| Simulation horizon | 10 |
| Time step | 0.01 |
| Random seed | 42 |
| Initial intensity | 500 |
| Initial wind anomaly | 1.5 |
| Initial humidity anomaly | -1.0 |

All other model parameters remain fixed.

Common random numbers are used across configurations to
support paired comparisons.

---

## 28. Numerical Implementation

The main mathematical implementation is contained in:

    src/dependence.py

This module provides analytical calculations for:

- Correlation validation.
- Environmental covariance.
- Stationary environmental covariance.
- Stationary environmental correlation.
- Exact OU transition covariance.
- Exact OU innovation correlation.
- Integrated environmental covariance.
- Correlated additive log-intensity moments.
- Expected nonlinear interaction effects.

The Monte Carlo implementation is contained in:

    src/monte_carlo.py

The original simulation engine is extended with:

    rho = 0.0

as an optional parameter.

The default value preserves the original independent model.

The simulation engine:

1. Validates the correlation parameter.
2. Constructs the numerical time grid.
3. Initializes independent underlying random streams.
4. Computes exact OU transition variances.
5. Computes the exact innovation correlation.
6. Generates correlated environmental innovations.
7. Updates log-intensity using the current environmental state.
8. Updates the environmental processes.
9. Recovers physical intensity through exponentiation.
10. Validates the resulting intensity trajectories.

---

## 29. Backward Compatibility Validation

Before running the main dependence experiments, the extended
Monte Carlo engine is tested against the original independent
simulation.

The following configurations are compared:

    simulate_many()

    simulate_many(rho=0.0)

The resulting trajectories are exactly equal.

The previous independent interaction baseline is also reproduced:

    Mean = 1596.96718

    Standard deviation = 1438.32228

This confirms that the environmental dependence extension
preserves the numerical results from Q2–Q5 when rho = 0.

Invalid correlations outside [-1,1] are rejected.

---

## 30. Correlated OU Validation

The correlated environmental processes are independently
validated using 10,000 simulated trajectories per scenario.

The following correlations are investigated:

    rho = -1.0, -0.5, 0.0, 0.5, 1.0

### 30.1 Environmental Covariance Validation

| Rho | Theoretical | Empirical | Error | Z-score |
|---|---:|---:|---:|---:|
| -1.0 | -0.53846 | -0.52302 | 0.01544 | 1.99921 |
| -0.5 | -0.26923 | -0.25995 | 0.00928 | 1.50862 |
| 0.0 | 0.00000 | 0.00181 | 0.00181 | 0.32691 |
| 0.5 | 0.26923 | 0.26308 | -0.00615 | -0.99946 |
| 1.0 | 0.53846 | 0.52303 | -0.01543 | -1.99880 |

The empirical covariance follows the theoretical relationship
closely.

The two extreme correlations have standardized errors close
to two standard errors and slightly outside the conventional
1.96 reference boundary.

These differences are compatible with finite-sample variation
and do not indicate a substantial implementation discrepancy.

### 30.2 Environmental Correlation Validation

| Rho | Theoretical | Empirical |
|---|---:|---:|
| -1.0 | -0.97303 | -0.97232 |
| -0.5 | -0.48651 | -0.48283 |
| 0.0 | 0.00000 | 0.00335 |
| 0.5 | 0.48651 | 0.48755 |
| 1.0 | 0.97303 | 0.97233 |

The empirical environmental correlations closely reproduce
their analytical values.

### 30.3 Marginal Variance Validation

| Rho | Wind Theory | Wind MC | Humidity Theory | Humidity MC |
|---|---:|---:|---:|---:|
| -1.0 | 0.62500 | 0.60427 | 0.48998 | 0.47884 |
| -0.5 | 0.62500 | 0.60427 | 0.48998 | 0.47967 |
| 0.0 | 0.62500 | 0.60427 | 0.48998 | 0.48139 |
| 0.5 | 0.62500 | 0.60427 | 0.48998 | 0.48184 |
| 1.0 | 0.62500 | 0.60427 | 0.48998 | 0.47884 |

The wind simulation is identical across correlation scenarios
because the same underlying wind random stream is preserved.

Humidity trajectories change with rho, but their theoretical
marginal distribution remains unchanged.

The numerical results support this property.

### 30.4 Integrated Covariance Validation

| Rho | Theoretical | Empirical | Error |
|---|---:|---:|---:|
| -1.0 | -13.18297 | -13.10308 | 0.07989 |
| -0.5 | -6.59148 | -6.68395 | -0.09247 |
| 0.0 | 0.00000 | -0.15297 | -0.15297 |
| 0.5 | 6.59148 | 6.41899 | -0.17249 |
| 1.0 | 13.18297 | 13.10281 | -0.08016 |

The empirical integrated covariance follows the analytical
linear relationship with rho.

The largest absolute discrepancy is approximately 0.17249.

The differences arise from finite-sample variation and numerical
approximation of the environmental integrals.

---

## 31. Environmental Validation Figures

### Figure 1 — Integrated OU Covariance

The theoretical integrated covariance is a straight line in rho.

The empirical values closely follow this line.

This validates the covariance structure of the integrated
environmental processes used in the terminal intensity model.

### Figure 2 — Environmental Correlation Over Time

The empirical environmental correlations closely follow the
theoretical time-dependent correlations.

For rho = +/-1, the environmental state correlation gradually
approaches approximately +/-0.97301.

This demonstrates the distinction between Brownian correlation
and environmental state correlation.

At exactly t = 0, the environmental state correlation is
undefined because both initial states are deterministic.

The initial plotted value should be interpreted as the
right-hand limiting correlation.

### Figure 3 — Correlated OU Covariance

The environmental covariance starts at zero and approaches
its stationary value.

Negative correlations produce negative covariance.

Positive correlations produce positive covariance.

The independent scenario remains close to zero.

### Figure 4 — Marginal OU Moments

The marginal means and variances remain consistent with
their theoretical OU expressions across all correlations.

This demonstrates that environmental dependence changes the
joint distribution without changing the marginal distributions.

---

## 32. Terminal Intensity Statistics

### Additive Model

| Rho | Mean | Std | Median | Q97.5 |
|---|---:|---:|---:|---:|
| -1.00 | 1711.14 | 1958.51 | 1117.78 | 6646.28 |
| -0.75 | 1655.39 | 1721.81 | 1133.05 | 6214.46 |
| -0.50 | 1597.50 | 1536.46 | 1131.97 | 5682.87 |
| -0.25 | 1540.83 | 1373.36 | 1129.77 | 5162.10 |
| 0.00 | 1485.39 | 1222.71 | 1132.66 | 4727.54 |
| 0.25 | 1431.01 | 1078.48 | 1134.20 | 4242.22 |
| 0.50 | 1377.48 | 936.78 | 1132.60 | 3816.80 |
| 0.75 | 1324.41 | 794.63 | 1132.74 | 3376.68 |
| 1.00 | 1268.33 | 646.63 | 1132.61 | 2922.35 |

### Interaction Model

| Rho | Mean | Std | Median | Q97.5 |
|---|---:|---:|---:|---:|
| -1.00 | 2355.28 | 3624.93 | 1371.47 | 10268.79 |
| -0.75 | 2132.14 | 2769.96 | 1321.01 | 8915.79 |
| -0.50 | 1932.73 | 2180.59 | 1271.81 | 7555.91 |
| -0.25 | 1755.56 | 1759.88 | 1209.00 | 6359.37 |
| 0.00 | 1596.97 | 1438.32 | 1166.65 | 5434.03 |
| 0.25 | 1454.09 | 1175.57 | 1116.11 | 4536.52 |
| 0.50 | 1324.70 | 950.73 | 1068.76 | 3809.52 |
| 0.75 | 1206.88 | 753.09 | 1021.44 | 3163.96 |
| 1.00 | 1096.28 | 574.30 | 975.83 | 2542.04 |

### Interpretation

Increasing rho reduces mean terminal intensity in both models.

The interaction model is more sensitive to environmental
dependence than the additive model.

The additive median remains approximately constant, consistent
with the theoretical invariance of expected log-intensity.

The interaction median decreases as rho increases because
the expected log-intensity also decreases.

Negative environmental dependence substantially increases
dispersion and upper intensity quantiles.

At rho = -1, the interaction model has a 97.5th percentile
of approximately 10268.79 kW/m.

At rho = 1, this decreases to approximately 2542.04 kW/m.

---

## 33. Theoretical Log-Intensity Validation

| Rho | Additive Theory | Additive MC | Interaction Theory | Interaction MC |
|---|---:|---:|---:|---:|
| -1.00 | 7.02569 | 7.01842 | 7.27066 | 7.26229 |
| -0.75 | 7.02569 | 7.02510 | 7.22095 | 7.21990 |
| -0.50 | 7.02569 | 7.02817 | 7.17125 | 7.17362 |
| -0.25 | 7.02569 | 7.03044 | 7.12154 | 7.12646 |
| 0.00 | 7.02569 | 7.03215 | 7.07184 | 7.07871 |
| 0.25 | 7.02569 | 7.03335 | 7.02213 | 7.03041 |
| 0.50 | 7.02569 | 7.03400 | 6.97243 | 6.98151 |
| 0.75 | 7.02569 | 7.03384 | 6.92273 | 6.93173 |
| 1.00 | 7.02569 | 7.03007 | 6.87302 | 6.87807 |

All 18 theoretical expected log-intensity values fall within
their corresponding approximate pointwise 95% Monte Carlo
confidence intervals.

The largest absolute standardized log-mean error is:

    Maximum absolute Z-score = 1.57347

This provides strong numerical support for the analytical
expected log-intensity expressions.

The comparisons are correlated across scenarios because
common random numbers are used.

They should not be interpreted as 18 independent statistical tests.

---

## 34. Additive Log-Variance Validation

| Rho | Theoretical | Empirical |
|---|---:|---:|
| -1.00 | 0.86195 | 0.86153 |
| -0.75 | 0.78285 | 0.78650 |
| -0.50 | 0.70375 | 0.70851 |
| -0.25 | 0.62466 | 0.62965 |
| 0.00 | 0.54556 | 0.55024 |
| 0.25 | 0.46646 | 0.47043 |
| 0.50 | 0.38736 | 0.39038 |
| 0.75 | 0.30827 | 0.31023 |
| 1.00 | 0.22917 | 0.23049 |

The empirical log-variance follows the theoretical linear
relationship with rho.

The agreement supports the exact correlated OU innovation
implementation and the integrated covariance derivation.

---

## 35. Terminal Intensity Distribution Interpretation

The empirical intensity distributions remain positively skewed.

Negative correlation produces more dispersed distributions
with substantially larger upper intensity quantiles.

Positive correlation produces more concentrated distributions.

The interaction model further increases upper intensity outcomes
under negative dependence.

Under sufficiently positive dependence, the interaction
distribution shifts toward lower intensities.

The empirical cumulative distribution functions illustrate
these changes in dispersion and exceedance probabilities.

However, the results do not establish a universal stochastic
ordering across every possible intensity threshold.

---

## 36. Threshold Crossing Analysis

Define the first passage time:

    tau_K = inf {t >= 0 : I(t) >= K}

The crossing probability is:

    p_K(rho) = P(tau_K <= T)

The numerical estimate is:

    p_hat_K(rho) =
        (1/N) sum_i 1{max_n I_i(t_n) >= K}

The following thresholds are investigated:

    K = 1000, 2000, 3000, 4000, 5000

The same numerical time grid and simulation horizon are used
for every correlation scenario.

---

## 37. Additive Threshold Crossing Probabilities

All values are expressed as percentages.

| Rho | K=1000 | K=2000 | K=3000 | K=4000 | K=5000 |
|---|---:|---:|---:|---:|---:|
| -1.00 | 72.29 | 33.97 | 18.22 | 10.79 | 6.42 |
| -0.75 | 72.61 | 33.06 | 17.39 | 9.78 | 5.72 |
| -0.50 | 72.87 | 32.11 | 16.17 | 8.47 | 4.73 |
| -0.25 | 73.26 | 30.93 | 14.62 | 7.12 | 3.70 |
| 0.00 | 73.77 | 29.24 | 12.79 | 5.60 | 2.72 |
| 0.25 | 74.48 | 27.29 | 10.44 | 4.18 | 1.83 |
| 0.50 | 75.21 | 25.00 | 8.23 | 2.81 | 1.18 |
| 0.75 | 76.53 | 21.61 | 6.02 | 1.67 | 0.53 |
| 1.00 | 77.81 | 17.02 | 3.25 | 0.68 | 0.14 |

For thresholds K >= 2000, crossing probability decreases
as rho increases.

However, the lowest threshold K = 1000 exhibits the
opposite relationship.

Its crossing probability increases from 72.29% to 77.81%.

This demonstrates that the effect of environmental dependence
on crossing probability depends on the threshold level.

Reducing intensity variability does not necessarily reduce
the probability of crossing every possible threshold.

The additive terminal log-mean remains fixed while the
distribution becomes more concentrated.

This helps explain why lower and higher threshold risks
can respond differently, although pathwise crossing also
depends on the full temporal covariance structure.

---

## 38. Interaction Threshold Crossing Probabilities

All values are expressed as percentages.

| Rho | K=1000 | K=2000 | K=3000 | K=4000 | K=5000 |
|---|---:|---:|---:|---:|---:|
| -1.00 | 78.25 | 42.87 | 26.20 | 17.37 | 12.25 |
| -0.75 | 77.57 | 40.56 | 23.87 | 15.23 | 10.47 |
| -0.50 | 76.99 | 38.07 | 21.36 | 13.30 | 8.31 |
| -0.25 | 76.22 | 35.18 | 18.52 | 10.68 | 6.19 |
| 0.00 | 75.59 | 31.97 | 15.41 | 7.80 | 4.08 |
| 0.25 | 74.73 | 28.32 | 11.99 | 5.10 | 2.43 |
| 0.50 | 73.69 | 23.35 | 7.92 | 2.90 | 1.18 |
| 0.75 | 72.60 | 17.79 | 4.62 | 1.28 | 0.40 |
| 1.00 | 71.17 | 11.56 | 1.78 | 0.34 | 0.03 |

Unlike the additive model, the interaction model exhibits
decreasing empirical crossing probabilities across all
investigated thresholds as rho increases.

This is consistent with the combined effects of decreasing
expected log-intensity and decreasing log-intensity variance.

For K = 5000, the estimated crossing probability decreases
from 12.25% to 0.03% across the full correlation range.

The final estimate is based on only three observed crossings
and should be interpreted cautiously.

---

## 39. Conditional Mean First Passage Times

The conditional mean first passage time is:

    E[tau_K | tau_K <= T]

It is estimated using only trajectories that reach the
specified threshold within the simulation horizon.

### Selected Additive Results

| Rho | K=1000 | K=2000 | K=3000 | K=4000 | K=5000 |
|---|---:|---:|---:|---:|---:|
| -1.0 | 3.73387 | 6.26835 | 7.14486 | 7.73677 | 7.99977 |
| -0.5 | 3.80230 | 6.44039 | 7.39896 | 7.93834 | 8.25512 |
| 0.0 | 3.89479 | 6.74927 | 7.71501 | 8.13120 | 8.50261 |
| 0.5 | 4.02871 | 7.08891 | 8.02662 | 8.44473 | 8.88966 |
| 1.0 | 4.21782 | 7.49562 | 8.49483 | 8.97397 | 9.04929 |

### Selected Interaction Results

| Rho | K=1000 | K=2000 | K=3000 | K=4000 | K=5000 |
|---|---:|---:|---:|---:|---:|
| -1.0 | 3.42042 | 5.82382 | 6.68633 | 7.11295 | 7.49457 |
| -0.5 | 3.50568 | 6.05842 | 6.92331 | 7.49419 | 7.77839 |
| 0.0 | 3.62587 | 6.35238 | 7.32248 | 7.83714 | 8.10963 |
| 0.5 | 3.78399 | 6.79652 | 7.79292 | 8.27234 | 8.57390 |
| 1.0 | 4.01519 | 7.36026 | 8.42466 | 8.97588 | 8.64667 |

### Interpretation

Conditional mean passage times generally increase with rho.

This indicates that trajectories reaching the investigated
thresholds tend to do so later under positive dependence.

However, different correlation scenarios produce different
sets of crossing trajectories.

Therefore, conditional passage-time differences do not
establish a trajectory-by-trajectory ordering.

At high thresholds and strong positive correlation, very
few trajectories cross the threshold.

This makes conditional passage-time estimates unstable.

---

## 40. Rare-Event Interpretation

For K = 5000:

| Rho | Additive Crossings | Interaction Crossings |
|---|---:|---:|
| -1.00 | 642 | 1225 |
| -0.75 | 572 | 1047 |
| -0.50 | 473 | 831 |
| -0.25 | 370 | 619 |
| 0.00 | 272 | 408 |
| 0.25 | 183 | 243 |
| 0.50 | 118 | 118 |
| 0.75 | 53 | 40 |
| 1.00 | 14 | 3 |

At rho = 1, the interaction conditional mean first passage
time is calculated from only three trajectories.

The apparent decrease in conditional passage time at the
largest correlation should therefore not be interpreted as
a reliable reversal of the overall passage-time relationship.

The visualization uses hollow markers to identify estimates
based on fewer than 30 observed crossings.

This is a practical warning criterion rather than a formal
statistical guarantee.

Larger Monte Carlo samples or specialized rare-event methods
would be required for more precise estimates.

---

## 41. Paired Interaction Analysis

For each correlation scenario, the additive and interaction
models use common random numbers.

Define the paired terminal difference:

    D_i(rho) =
        I_gamma,i(T;rho) - I_0,i(T;rho)

The estimated mean interaction effect is:

    D_bar(rho) =
        (1/N) sum_i D_i(rho)

The estimated standard error is:

    SE_D =
        Std(D_i) / sqrt(N)

An approximate pointwise 95% confidence interval is:

    D_bar +/- 1.96 SE_D

This paired approach reduces uncertainty relative to comparing
independently generated samples.

---

## 42. Arithmetic Mean Interaction Effects

| Rho | Mean Difference | Percentage Change | CI Lower | CI Upper |
|---|---:|---:|---:|---:|
| -1.00 | 644.14 | 37.64% | 607.85 | 680.43 |
| -0.75 | 476.75 | 28.80% | 453.21 | 500.29 |
| -0.50 | 335.23 | 20.98% | 320.46 | 350.00 |
| -0.25 | 214.73 | 13.94% | 205.51 | 223.94 |
| 0.00 | 111.58 | 7.51% | 105.82 | 117.34 |
| 0.25 | 23.08 | 1.61% | 19.34 | 26.82 |
| 0.50 | -52.78 | -3.83% | -55.66 | -49.91 |
| 0.75 | -117.53 | -8.87% | -120.32 | -114.74 |
| 1.00 | -172.05 | -13.56% | -174.98 | -169.12 |

All nine approximate pointwise paired confidence intervals
exclude zero.

The estimated interaction effect is positive for correlations
from -1 through 0.25.

The estimated interaction effect is negative for correlations
from 0.50 through 1.

Therefore, the arithmetic mean interaction effect reverses
direction between rho = 0.25 and rho = 0.50.

---

## 43. Fraction of Trajectories With Increased Intensity

The fraction of paired trajectories for which the interaction
increases terminal intensity is:

    P_hat(I_gamma(T) > I_0(T))

| Rho | Fraction Increased |
|---|---:|
| -1.00 | 100.00% |
| -0.75 | 99.96% |
| -0.50 | 97.81% |
| -0.25 | 89.38% |
| 0.00 | 72.95% |
| 0.25 | 51.43% |
| 0.50 | 27.87% |
| 0.75 | 9.51% |
| 1.00 | 0.25% |

Under strong negative dependence, nearly every simulated
trajectory experiences increased terminal intensity due to
the interaction.

Under strong positive dependence, the interaction decreases
terminal intensity for almost every simulated trajectory.

The 100% result at rho = -1 describes the finite simulated
sample and does not establish an almost-sure mathematical result.

---

## 44. Paired Log-Interaction Validation

Define:

    D_log,i =
        log(I_gamma,i(T)) - log(I_0,i(T))

Under common environmental trajectories and intensity innovations:

    D_log,i =
        gamma integral_0^T Z_W,i(t) Z_H,i(t) dt

up to numerical approximation of the drift integral.

The empirical paired log difference can therefore be compared
directly against the analytical expected interaction effect.

| Rho | Theory | Empirical | Error |
|---|---:|---:|---:|
| -1.00 | 0.24497 | 0.24388 | -0.00110 |
| -0.75 | 0.19527 | 0.19480 | -0.00046 |
| -0.50 | 0.14556 | 0.14544 | -0.00012 |
| -0.25 | 0.09586 | 0.09602 | 0.00017 |
| 0.00 | 0.04615 | 0.04656 | 0.00041 |
| 0.25 | -0.00355 | -0.00294 | 0.00061 |
| 0.50 | -0.05325 | -0.05249 | 0.00076 |
| 0.75 | -0.10296 | -0.10210 | 0.00086 |
| 1.00 | -0.15266 | -0.15200 | 0.00067 |

All nine analytical expectations fall within their corresponding
approximate pointwise paired Monte Carlo confidence intervals.

The largest absolute discrepancy is approximately 0.00110.

This provides direct numerical validation of the analytical
interaction effect.

---

## 45. Numerical Arithmetic Mean Critical Correlation

The analytical expected log-interaction effect changes sign at:

    rho_log_star = 0.2321

However, the arithmetic mean interaction effect changes sign
at a different correlation.

The empirical arithmetic mean differences are:

    Delta_mean(0.25) = 23.08018

    Delta_mean(0.50) = -52.78144

Using linear interpolation:

    rho_mean_star =
        rho_1
        - Delta_1 (rho_2 - rho_1)
          / (Delta_2 - Delta_1)

Substituting:

    rho_mean_star =
        0.25
        - 23.08018(0.50 - 0.25)
          / (-52.78144 - 23.08018)

Therefore:

    rho_mean_star = approximately 0.3261

This is a numerical interpolation between two simulated
correlation scenarios.

It is not an exact analytical root or a confidence interval
for the true arithmetic mean sign-change point.

A finer correlation grid would be required to estimate the
arithmetic mean critical correlation more precisely.

---

## 46. Why the Two Critical Correlations Differ

The expected log-interaction effect is:

    E[log(I_gamma(T)) - log(I_0(T))]

The arithmetic mean interaction effect is:

    E[I_gamma(T) - I_0(T)]

These are different mathematical quantities.

Using the paired log representation:

    I_gamma(T) =
        I_0(T) exp(gamma J_WH(T))

Therefore:

    E[I_gamma(T) - I_0(T)] =

        E[
            I_0(T)
            (exp(gamma J_WH(T)) - 1)
        ]

The arithmetic mean interaction effect depends on:

- The distribution of the integrated interaction.
- The exponential transformation.
- The dependence between baseline intensity and the interaction.
- The frequency and magnitude of extreme intensity outcomes.

Consequently:

    E[gamma J_WH(T)] = 0

does not imply:

    E[I_gamma(T) - I_0(T)] = 0

This explains why:

    rho_log_star  = approximately 0.2321

    rho_mean_star = approximately 0.3261

The distinction illustrates the importance of nonlinear
transformations and distributional properties in stochastic
modeling.

---

## 47. Paired Interaction Effects on Crossing Risk

Define:

    Delta_p(K,rho) =
        P_gamma(tau_K <= T)
        - P_0(tau_K <= T)

The paired estimator is constructed from:

    D_i =
        1{M_gamma,i(T) >= K}
        - 1{M_0,i(T) >= K}

Each paired difference belongs to:

    {-1,0,1}

The estimated mean paired difference is the difference
between the two crossing probabilities.

### Numerical Results

All differences are expressed in percentage points.

| Rho | K=1000 | K=2000 | K=3000 | K=4000 | K=5000 |
|---|---:|---:|---:|---:|---:|
| -1.00 | +5.96 | +8.90 | +7.98 | +6.58 | +5.83 |
| -0.75 | +4.96 | +7.50 | +6.48 | +5.45 | +4.75 |
| -0.50 | +4.12 | +5.96 | +5.19 | +4.83 | +3.58 |
| -0.25 | +2.96 | +4.25 | +3.90 | +3.56 | +2.49 |
| 0.00 | +1.82 | +2.73 | +2.62 | +2.20 | +1.36 |
| 0.25 | +0.25 | +1.03 | +1.55 | +0.92 | +0.60 |
| 0.50 | -1.52 | -1.65 | -0.31 | +0.09 | 0.00 |
| 0.75 | -3.93 | -3.82 | -1.40 | -0.39 | -0.13 |
| 1.00 | -6.64 | -5.46 | -1.47 | -0.34 | -0.11 |

The interaction increases crossing risk under negative
correlation.

Under sufficiently positive correlation, the interaction
decreases crossing risk.

The correlation at which the interaction changes direction
depends on the threshold.

There is therefore no single universal critical correlation
for all wildfire risk metrics.

---

## 48. Statistical Interpretation of Crossing Effects

The interaction heatmap identifies approximate pointwise
95% paired confidence intervals that exclude zero.

Of the 45 correlation-threshold comparisons, 42 have
intervals excluding zero.

The three comparisons whose intervals include zero are:

    rho = 0.25, K = 1000

    rho = 0.50, K = 4000

    rho = 0.50, K = 5000

These scenarios occur near the transition between positive
and negative interaction effects.

The intervals are pointwise rather than simultaneous.

No multiple-comparison correction is applied.

Normal-approximation intervals may also be unreliable for
very rare crossings or small numbers of discordant pairs.

The heatmap therefore provides exploratory evidence of
correlation-dependent interaction effects rather than a
formal simultaneous significance guarantee.

---

## 49. Consequences of Environmental Independence

The original model assumes:

    rho = 0

To investigate the consequences of this assumption, every
correlation scenario is compared with the corresponding
independent model.

Define:

    Delta_independence(rho) =
        E[I(T;rho)] - E[I(T;0)]

A positive value indicates that the correlated model produces
a larger mean intensity than the independent model.

A negative value indicates that the correlated model produces
a smaller mean intensity.

---

## 50. Dependence Effects Relative to Independence

### Additive Model

| Rho | Mean Difference | Percentage Change |
|---|---:|---:|
| -1.00 | 225.75 | +15.20% |
| -0.75 | 170.00 | +11.44% |
| -0.50 | 112.11 | +7.55% |
| -0.25 | 55.44 | +3.73% |
| 0.00 | 0.00 | 0.00% |
| 0.25 | -54.38 | -3.66% |
| 0.50 | -107.91 | -7.26% |
| 0.75 | -160.98 | -10.84% |
| 1.00 | -217.06 | -14.61% |

### Interaction Model

| Rho | Mean Difference | Percentage Change |
|---|---:|---:|
| -1.00 | 758.31 | +47.48% |
| -0.75 | 535.17 | +33.51% |
| -0.50 | 335.76 | +21.02% |
| -0.25 | 158.59 | +9.93% |
| 0.00 | 0.00 | 0.00% |
| 0.25 | -142.87 | -8.95% |
| 0.50 | -272.27 | -17.05% |
| 0.75 | -390.09 | -24.43% |
| 1.00 | -500.69 | -31.35% |

The interaction model is substantially more sensitive to
environmental dependence.

For example, relative to independence:

    rho = -0.5:
        Interaction mean increases by 21.02%.

    rho = 0.5:
        Interaction mean decreases by 17.05%.

Therefore, even moderate environmental dependence can
materially alter the model's conclusions.

These differences are scenario comparisons, not estimates
of the actual error of a real-world forecasting system.

---

## 51. Dependence Effects on High-Intensity Risk

Consider the interaction model at K = 5000.

Under independence:

    P(tau_5000 <= T) = 4.08%

Under rho = -0.5:

    P(tau_5000 <= T) = 8.31%

The relative increase is:

    (8.31 - 4.08) / 4.08
        = approximately 103.68%

Thus, moderate negative dependence more than doubles the
estimated high-intensity crossing probability relative to
independence.

Under rho = 0.5:

    P(tau_5000 <= T) = 1.18%

The relative decrease is:

    (4.08 - 1.18) / 4.08
        = approximately 71.08%

The independence assumption can therefore produce either
underestimation or overestimation of simulated risk.

The direction depends on the environmental dependence structure
and the chosen risk metric.

---

## 52. Wildfire Dependence Figures

### Figure 5 — Terminal Wildfire Intensity Under Dependence

This figure contains:

- Mean terminal intensity.
- Standard deviation.
- Median.
- 97.5th percentile.
- Expected log-intensity.
- Log-intensity variance.

The additive expected log-intensity and median remain
approximately constant.

The interaction expected log-intensity decreases with rho.

The mean, standard deviation, and upper quantiles decrease
in both models.

The figure demonstrates the increased sensitivity introduced
by the nonlinear interaction.

### Figure 6 — Theoretical and Numerical Validation

This figure compares analytical and numerical results for:

- Expected log-intensity.
- Additive log-intensity variance.
- Additive arithmetic mean.
- Additive standard deviation.
- Additive 97.5th percentile.
- Standardized log-mean errors.

The empirical results closely follow the theoretical curves.

All standardized log-mean errors remain within the approximate
95% reference range.

### Figure 7 — Terminal Intensity Distributions

This figure compares empirical intensity distributions across
selected correlation scenarios.

Negative dependence produces greater dispersion and larger
upper intensity outcomes.

Positive dependence produces more concentrated distributions.

The empirical cumulative distribution functions further
illustrate the changing exceedance probabilities.

The displayed intensity range is truncated at the pooled
99th percentile for visualization purposes.

Histogram densities are normalized using the full samples.

### Figure 8 — Threshold Crossing Risk

This figure displays crossing probabilities for all five
intensity thresholds.

The additive model exhibits increasing crossing probability
at K = 1000 but decreasing crossing probability at the
higher thresholds.

The interaction model exhibits decreasing crossing probability
at every investigated threshold.

The logarithmic panel for K = 5000 highlights the substantial
changes in rare high-intensity events.

### Figure 9 — Conditional First Passage Times

This figure displays conditional mean first passage times.

Passage times generally increase as rho increases.

The final panel displays the number of observed crossings
at K = 5000.

Hollow markers identify conditional estimates based on
fewer than 30 crossings.

This prevents unstable rare-event estimates from being
interpreted without considering their sample sizes.

### Figure 10 — Paired Wind-Humidity Interaction Effects

This figure displays:

- Arithmetic mean interaction differences.
- Percentage changes in mean intensity.
- Fraction of trajectories with increased intensity.
- Interaction effects on crossing risk.

The interaction changes from intensity-amplifying under
negative dependence to intensity-reducing under sufficiently
positive dependence.

### Figure 11 — Interaction Crossing-Risk Heatmap

The heatmap summarizes the interaction effect across nine
correlations and five intensity thresholds.

Red cells indicate increased crossing risk.

Blue cells indicate decreased crossing risk.

Asterisks identify approximate pointwise paired confidence
intervals excluding zero.

The figure demonstrates that interaction effects depend
jointly on environmental correlation and threshold severity.

### Figure 12 — Consequences of Environmental Independence

This figure compares each correlation scenario with rho = 0.

The interaction model exhibits substantially larger changes
in mean intensity and crossing risk than the additive model.

The figure demonstrates the sensitivity of model conclusions
to the independence assumption.

### Figure 13 — Correlation-Dependent Sign Changes

This figure compares two critical correlations:

    Analytical expected log-interaction root:
        rho = 0.2321

    Interpolated arithmetic mean interaction root:
        rho = 0.3261

The left panel demonstrates the exact linear relationship
between expected log-interaction and rho.

The right panel demonstrates the numerical arithmetic mean
interaction effect.

The difference between the roots illustrates the importance
of the nonlinear transformation from log-intensity to intensity.

---

## 53. Automated Testing

The dependence test suite is contained in:

    tests/test_dependence.py

The tests cover:

| Category | Number of Tests |
|---|---:|
| Correlation parameter validation | 9 |
| Backward compatibility | 1 |
| Boundary correlation simulation | 3 |
| Invalid simulation correlation | 1 |
| Reproducibility and dependence | 2 |
| Analytical reference identities | 6 |
| Numerical agreement | 6 |
| Total | 28 |

The tests verify:

- Valid correlations are accepted.
- Invalid correlations are rejected.
- Default simulation equals explicit independence.
- Simulations remain finite and positive at rho = +/-1.
- Identical random seeds reproduce identical trajectories.
- Changing rho affects the simulated trajectories.
- Analytical integrated covariance satisfies the expected identities.
- Additive log-variance decreases with rho.
- The analytical critical correlation is approximately 0.2321.
- Monte Carlo additive log-variance agrees with theory.
- Paired log-interaction differences agree with theory.

### Dependence Test Results

    28 passed in 6.13s

### Full Project Test Results

    50 passed in 6.15s

All implemented automated tests passed.

The analytical reference identity tests verify the mathematical
reference calculations, while the numerical agreement tests
connect these calculations to the Monte Carlo implementation.

---

## 54. Files Used in Q6

### Mathematical Implementation

    src/dependence.py

### Monte Carlo Engine

    src/monte_carlo.py

### Validation Experiments

    experiments/validate_dependence_engine.py

    experiments/validate_correlated_ou.py

### Main Dependence Experiment

    experiments/dependence_analysis.py

### Visualization

    experiments/plot_dependence.py

### Automated Testing

    tests/test_dependence.py

### Generated Results

    results/

### Generated Figures

    results/figures/

---

## 55. Reproducing Q6

Execute the following commands from the project root directory.

### Analytical Calculations

    python -m src.dependence

### Validate the Extended Monte Carlo Engine

    python -m experiments.validate_dependence_engine

### Validate Correlated OU Processes

    python -m experiments.validate_correlated_ou

### Run the Main Dependence Experiment

    python -m experiments.dependence_analysis

### Generate Dependence Figures

    python -m experiments.plot_dependence

### Run Dependence Tests

    python -m pytest tests/test_dependence.py -v

### Run the Entire Project Test Suite

    python -m pytest -q

---

## 56. Limitations

### 56.1 Brownian Correlation Is Not Calibrated

The correlation scenarios are hypothetical.

The historical cross-sectional wind-humidity correlation is
not treated as a direct estimate of Brownian innovation
correlation.

### 56.2 Constant Correlation

The model assumes a constant rho throughout each trajectory.

Real environmental dependence may vary with time, season,
location, and meteorological conditions.

### 56.3 Gaussian Dependence

The dependence structure is represented through correlated
Gaussian Brownian innovations.

Non-Gaussian dependence, asymmetric dependence, and more
complex environmental interactions are not investigated.

### 56.4 Independent Intrinsic Intensity Noise

The intrinsic intensity Brownian motion remains independent
of the environmental Brownian motions.

Alternative dependence between intensity noise and environmental
noise could be considered in future work.

### 56.5 Illustrative Parameters

The numerical results depend on the assumed model parameters,
initial conditions, and simulation horizon.

The identified critical correlations are specific to the
investigated configuration.

### 56.6 Rare-Event Sampling

Some high-intensity crossing probabilities are estimated from
very few trajectories.

Their relative sampling uncertainty is substantial.

### 56.7 Numerical Time Discretization

Although the environmental transitions are exact, the
intensity drift integral and first passage times are
numerically approximated.

Between-grid threshold crossings may be missed.

### 56.8 Arithmetic Mean Critical Correlation

The arithmetic mean sign-change point is obtained using
linear interpolation between two simulated scenarios.

It is not an exact analytical result.

### 56.9 Alternative Stochastic Assumptions

This research question focuses on relaxing environmental
independence and comparing additive versus nonlinear
interaction structures.

Alternative marginal stochastic processes, such as bounded
humidity diffusions or non-Gaussian environmental processes,
are left for future investigation.

### 56.10 Model Validation Versus Predictive Validation

The analytical comparisons and automated tests validate
mathematical and numerical consistency.

They do not establish the predictive accuracy of the
proposed model against actual wildfire intensity trajectories.

---

## 57. Main Findings

### Finding 1 — Dependence Changes the Joint Environmental Distribution

Introducing Brownian correlation changes environmental
covariance without changing the individual OU marginal
distributions.

### Finding 2 — Negative Dependence Amplifies Additive Variability

Negative correlation increases additive log-intensity
variance because wind and humidity enter the growth rate
with opposite signs.

### Finding 3 — Additive Log-Mean Remains Constant

The expected additive log-intensity is independent of rho.

The arithmetic mean changes through the variance-dependent
exponential transformation.

### Finding 4 — Interaction Effects Depend on Environmental Covariance

The expected nonlinear interaction effect depends on both
the initial environmental conditions and the evolving
environmental covariance.

### Finding 5 — Interaction Effects Can Reverse Direction

The interaction increases mean terminal intensity by 37.64%
at rho = -1 but decreases it by 13.56% at rho = 1.

### Finding 6 — Different Metrics Have Different Critical Correlations

The expected log-interaction effect changes sign at
approximately rho = 0.2321.

The arithmetic mean interaction effect changes sign near
rho = 0.3261 based on numerical interpolation.

Threshold crossing effects exhibit additional
threshold-dependent sign changes.

### Finding 7 — Independence Can Materially Affect Risk Estimates

For the interaction model, moderate negative dependence
at rho = -0.5 more than doubles the estimated probability
of crossing K = 5000 relative to independence.

Moderate positive dependence substantially reduces this
probability.

### Finding 8 — The Interaction Model Is More Dependence-Sensitive

The interaction model exhibits larger changes in terminal
intensity and high-intensity risk across correlation scenarios.

### Finding 9 — Numerical Results Agree With Analytical Theory

The environmental covariance, integrated covariance,
additive lognormal moments, and expected interaction effects
are closely reproduced by Monte Carlo simulation.

All 28 dependence tests and all 50 project tests passed.

---

## 58. Conclusion

This research question extends the original stochastic wildfire
model by introducing correlated environmental Brownian motions.

The extension preserves the marginal Ornstein-Uhlenbeck
distributions while introducing a mathematically tractable
environmental covariance structure.

The correlated additive model retains an exact lognormal
terminal intensity distribution.

Its expected log-intensity remains independent of environmental
correlation, while its variance depends explicitly on the
integrated environmental covariance.

The nonlinear interaction model introduces an additional
dependence mechanism through the product of wind and humidity.

The expected log-interaction effect can be derived analytically
and is linear in the Brownian correlation parameter.

The numerical experiments demonstrate that environmental
dependence can substantially alter wildfire intensity
distributions, threshold crossing probabilities, and
conditional first passage times.

Most importantly, the nonlinear interaction effect is not
universally positive or negative.

Its direction depends on the environmental dependence structure,
initial conditions, simulation horizon, and chosen response metric.

The analytical expected log-interaction effect changes sign at
approximately rho = 0.2321, while the arithmetic mean interaction
effect changes sign near rho = 0.3261 in the numerical experiment.

These findings demonstrate that assumptions about environmental
dependence are an important component of stochastic wildfire
modeling.

Ignoring dependence can materially change the conclusions
obtained from the proposed model.

The theoretical derivations, environmental validation experiments,
180,000 simulated wildfire trajectories, 13 generated figures,
and 28 passing dependence tests provide a comprehensive
investigation of the consequences of relaxing environmental
independence.

---

## 59. Completion Status

- [x] Correlated Brownian motion formulation
- [x] Environmental covariance derivation
- [x] Ito product-rule derivation
- [x] Stationary covariance and correlation
- [x] Exact correlated OU transitions
- [x] Correlated Monte Carlo implementation
- [x] Backward compatibility validation
- [x] Integrated OU covariance derivation
- [x] Correlated additive lognormal distribution
- [x] Expected nonlinear interaction derivation
- [x] Analytical critical correlation
- [x] Correlated OU numerical validation
- [x] Terminal intensity distribution experiments
- [x] Threshold crossing analysis
- [x] Conditional first passage time analysis
- [x] Paired interaction analysis
- [x] Independence assumption comparison
- [x] Numerical critical correlation estimation
- [x] Environmental validation figures
- [x] Wildfire dependence figures
- [x] Automated dependence testing
- [x] Mathematical interpretation
- [x] Research question documentation

Q6 Status: COMPLETED

All six research questions have now been implemented,
numerically investigated, and documented.