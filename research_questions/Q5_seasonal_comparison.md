# Q5 — Seasonal Environmental Comparison

## 1. Research Question

How do seasonal differences in environmental conditions influence
the stochastic evolution of wildfire intensity?

Specifically, we investigate:

- Expected wildfire intensity across seasons.
- Terminal intensity distributions.
- Variability and extreme intensity outcomes.
- Threshold crossing probabilities.
- Conditional first-passage times.
- Seasonal differences in environmental interaction effects.
- Agreement between analytical and numerical results.

The objective is to determine how historical seasonal wind and
humidity conditions influence the behavior of the coupled
stochastic wildfire model.

---

## 2. Historical Seasonal Data

The historical Alberta wildfire dataset contains 24,859 complete
observations used in the environmental analysis.

The observations are divided into four seasons.

| Season | Observations | Wind (km/h) | Humidity (%) | Spread (m/min) |
|---|---:|---:|---:|---:|
| Winter | 273 | 12.176 | 55.121 | 0.353 |
| Spring | 9457 | 9.966 | 38.515 | 0.893 |
| Summer | 12492 | 7.894 | 49.064 | 1.044 |
| Fall | 2637 | 9.692 | 49.872 | 0.431 |

The seasonal definitions are:

- Winter: December–February.
- Spring: March–May.
- Summer: June–August.
- Fall: September–November.

### 2.1 Historical observations

Winter has the highest recorded mean wind speed.

Spring has the lowest recorded mean relative humidity.

Summer has the highest historical mean spread rate.

However, these observations alone cannot determine which season
produces the highest simulated wildfire intensity.

Wildfire intensity depends on the combined environmental effects
and their stochastic evolution over time.

### 2.2 Interpretation of seasonal averages

The seasonal averages describe environmental conditions recorded
during historical wildfire assessments.

They are not general seasonal climate normals.

The unequal sample sizes also indicate that some seasonal averages,
particularly winter, are based on fewer observations.

---

## 3. Standardization of Seasonal Conditions

The mathematical model uses standardized environmental anomalies
rather than raw environmental measurements.

The historical normalization constants are:

| Variable | Mean | Standard deviation |
|---|---:|---:|
| Wind speed | 8.920230 | 8.415991 |
| Relative humidity | 45.203347 | 18.674584 |

For season $s$, define:

$$
z_{0,s}^W=
\frac{\overline W_s-\mu_W}{s_W}
$$

and:

$$
z_{0,s}^H=
\frac{\overline H_s-\mu_H}{s_H}
$$

### 3.1 Example: Spring

The spring mean wind speed is:

$$
\overline W_{\text{spring}}=9.966
$$

Therefore:

$$
z_{0,\text{spring}}^W
=
\frac{9.966-8.920230}{8.415991}
\approx0.12426
$$

Similarly:

$$
z_{0,\text{spring}}^H
=
\frac{38.515-45.203347}{18.674584}
\approx-0.35815
$$

Spring therefore represents above-average wind speed combined
with below-average humidity.

### 3.2 Standardized seasonal results

| Season | Wind anomaly | Humidity anomaly |
|---|---:|---:|
| Winter | 0.38686 | 0.53108 |
| Spring | 0.12426 | -0.35815 |
| Summer | -0.12194 | 0.20673 |
| Fall | 0.09170 | 0.25000 |

Positive anomalies indicate conditions above the historical mean.

Negative anomalies indicate conditions below the historical mean.

Spring is the only investigated season combining positive wind
anomalies with negative humidity anomalies.

---

## 4. Seasonal Stochastic Model

For each season, the environmental processes satisfy:

$$
dZ_t^W=
-\kappa_WZ_t^Wdt+
\sigma_WdB_t^{(W)}
$$

$$
dZ_t^H=
-\kappa_HZ_t^Hdt+
\sigma_HdB_t^{(H)}
$$

with seasonal initial conditions:

$$
Z_0^W=z_{0,s}^W
$$

$$
Z_0^H=z_{0,s}^H
$$

Wildfire intensity satisfies:

$$
dI_t=
I_t
\left(
\alpha_c+
\beta_WZ_t^W-
\beta_HZ_t^H+
\gamma Z_t^WZ_t^H
\right)dt
+
\sigma_I I_t dB_t^{(I)}
$$

The Brownian drivers are assumed independent throughout Q5.

### 4.1 Controlled seasonal comparison

All dynamic model parameters remain identical across seasons.

Only the initial environmental conditions change.

This isolates the influence of seasonal initial wind and humidity
conditions from changes in other model parameters.

### 4.2 Important distinction

The seasonal scenarios are not four independently calibrated
seasonal stochastic processes.

All environmental processes revert toward the same equilibrium:

$$
\lim_{t\to\infty}E[Z_t^W]=0
$$

$$
\lim_{t\to\infty}E[Z_t^H]=0
$$

Consequently, seasonal differences in expected environmental
conditions gradually diminish over time.

---

## 5. Model Configurations

The following parameters are held fixed.

| Parameter | Symbol | Value |
|---|---|---:|
| Wind mean reversion | $\kappa_W$ | 0.8 |
| Wind volatility | $\sigma_W$ | 1.0 |
| Humidity mean reversion | $\kappa_H$ | 0.5 |
| Humidity volatility | $\sigma_H$ | 0.7 |
| Baseline growth | $\alpha_c$ | 0.05 |
| Wind coefficient | $\beta_W$ | 0.12 |
| Humidity coefficient | $\beta_H$ | 0.10 |
| Intensity volatility | $\sigma_I$ | 0.15 |
| Initial intensity | $I_0$ | 500 |
| Simulation horizon | $T$ | 10 |
| Time step | $\Delta t$ | 0.01 |
| Trajectories per configuration | $N$ | 10000 |
| Random seed | | 42 |
| Fuel category | | C2 |

Two interaction configurations are investigated.

### Additive model

$$
\gamma=0
$$

### Interaction model

$$
\gamma=-0.04
$$

The dynamic coefficients are illustrative and are not estimated
directly from longitudinal wildfire intensity measurements.

---

## 6. Initial Seasonal Growth Rates

The proportional intensity growth function is:

$$
r(w,h)=
\alpha_c+\beta_Ww-\beta_Hh+\gamma wh
$$

For season $s$:

$$
r_s=
\alpha_c+
\beta_Wz_{0,s}^W-
\beta_Hz_{0,s}^H+
\gamma z_{0,s}^Wz_{0,s}^H
$$

The initial interaction contribution is:

$$
\Delta r_s=
\gamma z_{0,s}^Wz_{0,s}^H
$$

### 6.1 Example: Spring

For the additive model:

$$
r_{\text{spring}}^{(0)}
=
0.05+
0.12(0.12426)-
0.10(-0.35815)
$$

Therefore:

$$
r_{\text{spring}}^{(0)}
\approx0.10073
$$

The interaction contribution is:

$$
\Delta r_{\text{spring}}
=
-0.04(0.12426)(-0.35815)
$$

$$
\Delta r_{\text{spring}}
\approx0.00178
$$

Hence:

$$
r_{\text{spring}}^{(\gamma)}
\approx0.10251
$$

### 6.2 Seasonal growth results

| Season | Additive | Interaction | Interaction contribution |
|---|---:|---:|---:|
| Winter | 0.04331 | 0.03510 | -0.00822 |
| Spring | 0.10073 | 0.10251 | +0.00178 |
| Summer | 0.01469 | 0.01570 | +0.00101 |
| Fall | 0.03600 | 0.03509 | -0.00092 |

### 6.3 Interpretation

Spring has the highest initial proportional intensity growth rate.

Winter has the highest wind anomaly, but its positive humidity
anomaly counteracts the wind contribution.

Summer combines below-average wind with above-average humidity,
resulting in the lowest initial growth rate.

The interaction contribution is negative when both anomalies
have the same sign and positive when their signs differ,
because the interaction coefficient is negative.

These are initial growth rates, not constant growth rates over
the complete stochastic trajectory.

---

## 7. Analytical Seasonal Environmental Expectations

The explicit solution of the wind OU process is:

$$
Z_t^W=
z_{0,s}^We^{-\kappa_Wt}
+
\sigma_W
\int_0^t
e^{-\kappa_W(t-u)}dB_u^{(W)}
$$

Taking expectations:

$$
E[Z_t^W]
=
z_{0,s}^We^{-\kappa_Wt}
$$

Similarly:

$$
E[Z_t^H]
=
z_{0,s}^He^{-\kappa_Ht}
$$

Define the integrated environmental processes:

$$
J_T^W=\int_0^T Z_t^Wdt
$$

$$
J_T^H=\int_0^T Z_t^Hdt
$$

Using the OU expectations:

$$
E[J_T^W]
=
\int_0^T z_{0,s}^We^{-\kappa_Wt}dt
$$

Therefore:

$$
E[J_T^W]
=
z_{0,s}^W
\frac{1-e^{-\kappa_WT}}{\kappa_W}
$$

Similarly:

$$
E[J_T^H]
=
z_{0,s}^H
\frac{1-e^{-\kappa_HT}}{\kappa_H}
$$

For convenience, define:

$$
A_W(T)=
\frac{1-e^{-\kappa_WT}}{\kappa_W}
$$

$$
A_H(T)=
\frac{1-e^{-\kappa_HT}}{\kappa_H}
$$

Then:

$$
E[J_T^W]=z_{0,s}^WA_W(T)
$$

$$
E[J_T^H]=z_{0,s}^HA_H(T)
$$

These expressions describe the expected accumulated influence
of the initial seasonal environmental anomalies.

---

## 8. Analytical Seasonal Log-Intensity

Applying Itô's formula to:

$$
Y_t=\log I_t
$$

gives:

$$
dY_t=
\left(
\alpha_c+
\beta_WZ_t^W-
\beta_HZ_t^H+
\gamma Z_t^WZ_t^H-
\frac{\sigma_I^2}{2}
\right)dt
+
\sigma_I dB_t^{(I)}
$$

Integrating from zero to $T$:

$$
\begin{aligned}
Y_T={}&
\log I_0+
\left(\alpha_c-\frac{\sigma_I^2}{2}\right)T\\
&+\beta_WJ_T^W-\beta_HJ_T^H\\
&+\gamma\int_0^T Z_t^WZ_t^Hdt\\
&+\sigma_I B_T^{(I)}.
\end{aligned}
$$

Taking expectations and using:

$$
E[B_T^{(I)}]=0
$$

gives:

$$
\begin{aligned}
E[Y_T]={}&
\log I_0+
\left(\alpha_c-\frac{\sigma_I^2}{2}\right)T\\
&+\beta_WE[J_T^W]
-\beta_HE[J_T^H]\\
&+\gamma\int_0^T E[Z_t^WZ_t^H]dt.
\end{aligned}
$$

This expression applies to both the additive and interaction models.

---

## 9. Analytical Additive Seasonal Distribution

When:

$$
\gamma=0
$$

the terminal log-intensity is a linear combination of Gaussian
random variables.

Consequently:

$$
Y_{T,s}^{(0)}
\sim N(m_s,v_T)
$$

where:

$$
\begin{aligned}
m_s={}&
\log I_0+
\left(\alpha_c-\frac{\sigma_I^2}{2}\right)T\\
&+\beta_Wz_{0,s}^WA_W(T)\\
&-\beta_Hz_{0,s}^HA_H(T).
\end{aligned}
$$

The log-intensity variance is:

$$
v_T=
\beta_W^2V_W(T)
+
\beta_H^2V_H(T)
+
\sigma_I^2T
$$

where the integrated OU variance is:

$$
V(T)=
\frac{\sigma^2}{\kappa^2}
\left[
T-
\frac{2(1-e^{-\kappa T})}{\kappa}
+
\frac{1-e^{-2\kappa T}}{2\kappa}
\right]
$$

### 9.1 Constant additive log variance

The integrated OU variance does not depend on the deterministic
initial condition.

Therefore, all four seasonal additive models have the same
theoretical log-intensity variance.

For the selected configuration:

$$
v_T\approx0.54556
$$

The seasonal initial conditions change the log-intensity mean
but not its variance.

### 9.2 Lognormal intensity distribution

Since:

$$
I_{T,s}^{(0)}=e^{Y_{T,s}^{(0)}}
$$

we obtain:

$$
I_{T,s}^{(0)}
\sim\operatorname{Lognormal}(m_s,v_T)
$$

The theoretical arithmetic mean is:

$$
E[I_{T,s}^{(0)}]
=
e^{m_s+v_T/2}
$$

The variance is:

$$
\operatorname{Var}(I_{T,s}^{(0)})
=
(e^{v_T}-1)e^{2m_s+v_T}
$$

The median is:

$$
\operatorname{Median}(I_{T,s}^{(0)})
=
e^{m_s}
$$

The theoretical quantile is:

$$
q_p=
e^{m_s+\sqrt{v_T}\Phi^{-1}(p)}
$$

### 9.3 Theoretical results

| Season | Log mean | Mean | Standard deviation | Median |
|---|---:|---:|---:|---:|
| Winter | 6.55462 | 922.78515 | 786.03327 | 702.48010 |
| Spring | 6.69189 | 1058.56345 | 901.68994 | 805.84278 |
| Summer | 6.54276 | 911.90414 | 776.76477 | 694.19681 |
| Fall | 6.56620 | 933.53186 | 795.18738 | 710.66115 |

The theoretical additive mean ranking is:

$$
\text{Spring}>
\text{Fall}>
\text{Winter}>
\text{Summer}
$$

Because all seasons have the same log variance, their theoretical
additive means, medians and quantiles share the same ranking.

---

## 10. Analytical Seasonal Interaction Effect

Under independent environmental Brownian motions and deterministic
initial conditions:

$$
E[Z_t^WZ_t^H]
=
E[Z_t^W]E[Z_t^H]
$$

Substituting the OU expectations:

$$
E[Z_t^WZ_t^H]
=
z_{0,s}^Wz_{0,s}^H
e^{-(\kappa_W+\kappa_H)t}
$$

Integrating:

$$
\int_0^T E[Z_t^WZ_t^H]dt
=
z_{0,s}^Wz_{0,s}^H
\frac{
1-e^{-(\kappa_W+\kappa_H)T}
}{
\kappa_W+\kappa_H
}
$$

Define:

$$
A_{WH}(T)=
\frac{
1-e^{-(\kappa_W+\kappa_H)T}
}{
\kappa_W+\kappa_H
}
$$

The expected log-interaction effect is therefore:

$$
\boxed{
\Delta m_s=
\gamma z_{0,s}^Wz_{0,s}^HA_{WH}(T)
}
$$

The complete expected interaction log-intensity becomes:

$$
m_s^{(\gamma)}
=
m_s^{(0)}+\Delta m_s
$$

### 10.1 Theoretical interaction results

| Season | Additive log mean | Interaction log mean | Difference |
|---|---:|---:|---:|
| Winter | 6.55462 | 6.54830 | -0.00632 |
| Spring | 6.69189 | 6.69326 | +0.00137 |
| Summer | 6.54276 | 6.54353 | +0.00078 |
| Fall | 6.56620 | 6.56549 | -0.00071 |

The expected interaction effect is negative in winter and fall.

It is positive in spring and summer.

### 10.2 Geometric mean

The geometric mean of terminal intensity is:

$$
G_s=\exp(E[\log I_{T,s}])
$$

For the interaction model:

| Season | Geometric mean |
|---|---:|
| Winter | 698.05336 |
| Spring | 806.94702 |
| Summer | 694.73547 |
| Fall | 710.16002 |

The geometric mean is not generally equal to the arithmetic mean.

It is also not necessarily equal to the median when the
interaction term is present.

### 10.3 Why the interaction distribution is different

The interaction model contains:

$$
\int_0^T Z_t^WZ_t^Hdt
$$

Although both environmental processes are Gaussian, their product
is not generally Gaussian.

Therefore, the terminal interaction log-intensity is not
generally normally distributed.

The additive lognormal formulas cannot automatically be applied
to the interaction model.

Monte Carlo simulation is used to investigate its distribution.

---

## 11. Monte Carlo Experimental Methodology

### 11.1 Experimental design

Each season is simulated under both model configurations.

| Property | Value |
|---|---:|
| Seasons | 4 |
| Models per season | 2 |
| Configurations | 8 |
| Trajectories per configuration | 10000 |
| Total simulated trajectories | 80000 |
| Time steps per trajectory | 1000 |
| Recorded time points | 1001 |
| Thresholds per configuration | 5 |
| Threshold evaluations | 40 |

### 11.2 Simulation procedure

For each seasonal configuration:

1. Obtain the historical seasonal mean wind and humidity.
2. Standardize the seasonal environmental conditions.
3. Initialize the environmental OU processes.
4. Set the initial wildfire intensity to 500 kW/m.
5. Simulate the additive model.
6. Simulate the interaction model using corresponding random numbers.
7. Extract terminal intensity samples.
8. Calculate distributional statistics.
9. Calculate threshold crossing probabilities.
10. Calculate conditional first-passage times.
11. Compare analytical and numerical log-intensity expectations.
12. Calculate paired interaction effects.
13. Save numerical results for visualization.

### 11.3 Common random numbers

The same random seed is used across corresponding configurations.

This allows differences between model outcomes to be investigated
using paired trajectories rather than unrelated samples.

Common random numbers reduce comparison noise when corresponding
outcomes are positively correlated.

The configurations are therefore not statistically independent
of one another, although trajectories within each configuration
are independently generated.

---

## 12. Terminal Intensity Results

| Season | Model | Mean | SD | Median | Q97.5 |
|---|---|---:|---:|---:|---:|
| Winter | Additive | 926.06 | 762.29 | 706.15 | 2947.35 |
| Winter | Interaction | 926.73 | 798.73 | 689.32 | 3061.59 |
| Spring | Additive | 1062.62 | 874.70 | 810.28 | 3381.99 |
| Spring | Interaction | 1078.97 | 945.19 | 797.05 | 3618.75 |
| Summer | Additive | 915.00 | 753.19 | 697.72 | 2912.18 |
| Summer | Interaction | 922.89 | 796.90 | 686.60 | 3057.41 |
| Fall | Additive | 936.81 | 771.14 | 714.35 | 2981.57 |
| Fall | Interaction | 943.93 | 815.97 | 701.88 | 3127.90 |

All intensity statistics are measured in kW/m.

### 12.1 Seasonal mean comparison

Both models reproduce the theoretical seasonal mean ranking:

$$
\text{Spring}>
\text{Fall}>
\text{Winter}>
\text{Summer}
$$

Under the interaction model, spring has approximately 16.91%
higher mean terminal intensity than summer.

### 12.2 Interaction-induced distributional changes

| Season | Mean change | SD change | Median change | Q97.5 change |
|---|---:|---:|---:|---:|
| Winter | +0.07% | +4.78% | -2.38% | +3.88% |
| Spring | +1.54% | +8.06% | -1.63% | +7.00% |
| Summer | +0.86% | +5.80% | -1.59% | +4.99% |
| Fall | +0.76% | +5.81% | -1.75% | +4.91% |

The interaction increases the empirical arithmetic mean in
every seasonal scenario.

However, the increase is relatively small.

The changes in variability and upper quantiles are substantially
larger than the changes in the arithmetic mean.

The empirical median decreases in every seasonal scenario.

### 12.3 Mean versus median

The arithmetic mean is sensitive to large intensity outcomes.

The median represents the central ranked observation.

An increase in large intensity outcomes can raise the arithmetic
mean even when the median decreases.

This explains why an interaction model may exhibit:

- Higher arithmetic mean.
- Lower median.
- Higher standard deviation.
- Higher upper quantiles.

These findings indicate changes in the distribution rather than
a uniform increase in all wildfire trajectories.

---

## 13. Analytical and Numerical Validation

The analytical expected log-intensity is compared with the
Monte Carlo estimate:

$$
\widehat m_s=
\frac{1}{N}
\sum_{j=1}^N
\log I_{T,s}^{(j)}
$$

The estimated standard error is:

$$
SE(\widehat m_s)
=
\frac{s_{\log I}}{\sqrt N}
$$

Define the validation statistic:

$$
z_s=
\frac{
\widehat m_s-m_s^{\mathrm{theory}}
}{
SE(\widehat m_s)
}
$$

### 13.1 Expected log-intensity validation

| Season | Model | Theory | Numerical | Error | z-score |
|---|---|---:|---:|---:|---:|
| Winter | Additive | 6.55462 | 6.55965 | 0.00503 | 0.67838 |
| Winter | Interaction | 6.54830 | 6.55404 | 0.00575 | 0.77689 |
| Spring | Additive | 6.69189 | 6.69721 | 0.00532 | 0.71670 |
| Spring | Interaction | 6.69326 | 6.69885 | 0.00559 | 0.74597 |
| Summer | Additive | 6.54276 | 6.54764 | 0.00489 | 0.65892 |
| Summer | Interaction | 6.54353 | 6.54895 | 0.00542 | 0.73111 |
| Fall | Additive | 6.56620 | 6.57119 | 0.00499 | 0.67332 |
| Fall | Interaction | 6.56549 | 6.57105 | 0.00556 | 0.75029 |

All eight absolute z-scores are below one.

All theoretical expectations fall inside their corresponding
approximate pointwise 95% Monte Carlo confidence intervals.

### 13.2 Analytical interaction validation

For paired trajectories:

$$
L_{s,j}
=
\log I_{T,s,j}^{(\gamma)}
-
\log I_{T,s,j}^{(0)}
$$

Using common environmental and intensity innovations:

$$
L_{s,j}
=
\gamma
\int_0^T
Z_{t,s,j}^WZ_{t,s,j}^Hdt
$$

The numerical implementation approximates this integral using
the simulation grid.

The Monte Carlo expected log-interaction effect is:

$$
\widehat{\Delta m_s}
=
\frac{1}{N}
\sum_{j=1}^N L_{s,j}
$$

The analytical and numerical results are:

| Season | Theory | Numerical | Error | z-score |
|---|---:|---:|---:|---:|
| Winter | -0.00632 | -0.00561 | 0.00071 | 0.91138 |
| Spring | 0.00137 | 0.00164 | 0.00027 | 0.35224 |
| Summer | 0.00078 | 0.00131 | 0.00053 | 0.68863 |
| Fall | -0.00071 | -0.00014 | 0.00057 | 0.73966 |

All four absolute z-scores are below one.

### 13.3 Interpretation of validation

The numerical estimates are consistent with the derived
continuous-time analytical expectations.

The small discrepancies can arise from:

- Finite Monte Carlo sampling.
- Numerical approximation of the intensity drift integral.
- Rounding of historical seasonal inputs.

The validation supports the mathematical implementation.

It does not demonstrate that the model accurately predicts
real-world wildfire intensity.

Because common random numbers are used, the validation results
across seasons should not be interpreted as independent tests.

---

## 14. Seasonal Threshold Crossing Probabilities

Define the first-passage time:

$$
\tau_K=
\inf\{t\geq0:I_t\geq K\}
$$

The seasonal threshold crossing probability is:

$$
p_s(K)=P_s(\tau_K\leq T)
$$

Equivalently:

$$
p_s(K)=
P_s\left(
\sup_{0\leq t\leq T}I_t\geq K
\right)
$$

The Monte Carlo estimator is:

$$
\widehat p_s(K)
=
\frac{1}{N}
\sum_{j=1}^N
\mathbf1
\left\{
\max_n I_{s,n}^{(j)}\geq K
\right\}
$$

Five thresholds are investigated:

$$
K\in\{1000,2000,3000,4000,5000\}
$$

### 14.1 Additive model results

| Season | 1000 | 2000 | 3000 | 4000 | 5000 |
|---|---:|---:|---:|---:|---:|
| Winter | 43.01% | 10.67% | 3.07% | 1.10% | 0.40% |
| Spring | 51.71% | 15.31% | 4.78% | 1.83% | 0.76% |
| Summer | 42.04% | 10.33% | 2.91% | 1.06% | 0.39% |
| Fall | 43.57% | 10.99% | 3.21% | 1.12% | 0.42% |

### 14.2 Interaction model results

| Season | 1000 | 2000 | 3000 | 4000 | 5000 |
|---|---:|---:|---:|---:|---:|
| Winter | 41.36% | 11.04% | 3.42% | 1.30% | 0.55% |
| Spring | 50.84% | 15.70% | 5.49% | 2.33% | 1.01% |
| Summer | 40.96% | 10.86% | 3.40% | 1.30% | 0.53% |
| Fall | 42.32% | 11.52% | 3.60% | 1.45% | 0.59% |

### 14.3 Seasonal comparison

Spring produces the highest estimated threshold crossing
probability at every investigated threshold under both models.

For the interaction model at $K=2000$:

$$
p_{\text{spring}}=0.1570
$$

$$
p_{\text{summer}}=0.1086
$$

The absolute difference is:

$$
0.1570-0.1086=0.0484
$$

This corresponds to 4.84 percentage points.

The relative difference is:

$$
\frac{0.1570-0.1086}{0.1086}
\times100
\approx44.57\%
$$

At $K=5000$, the spring crossing probability is approximately
1.91 times the summer crossing probability.

These are comparisons between simulated seasonal scenarios,
not estimates of actual seasonal wildfire occurrence rates.

### 14.4 Threshold monotonicity

For:

$$
K_1<K_2
$$

we have:

$$
\{\tau_{K_2}\leq T\}
\subseteq
\{\tau_{K_1}\leq T\}
$$

Therefore:

$$
P(\tau_{K_2}\leq T)
\leq
P(\tau_{K_1}\leq T)
$$

All simulated seasonal configurations satisfy this mathematical
consistency property.

---

## 15. Conditional First-Passage Times

The conditional mean first-passage time is:

$$
E[\tau_K\mid\tau_K\leq T]
$$

Its Monte Carlo estimator is:

$$
\widehat{\mu}_{\tau_K}
=
\frac{
\sum_{j:\tau_K^{(j)}\leq T}
\tau_K^{(j)}
}{
N_K
}
$$

where $N_K$ is the number of trajectories that cross the
threshold during the simulation.

### 15.1 Additive model

| Season | 1000 | 2000 | 3000 | 4000 | 5000 |
|---|---:|---:|---:|---:|---:|
| Winter | 5.78786 | 7.77555 | 8.41756 | 8.85636 | 8.98450 |
| Spring | 5.31005 | 7.54157 | 8.17098 | 8.64710 | 8.97263 |
| Summer | 5.86284 | 7.80680 | 8.41821 | 8.89396 | 9.03462 |
| Fall | 5.76357 | 7.75637 | 8.41467 | 8.83580 | 8.99643 |

### 15.2 Interaction model

| Season | 1000 | 2000 | 3000 | 4000 | 5000 |
|---|---:|---:|---:|---:|---:|
| Winter | 5.79946 | 7.71680 | 8.26050 | 8.69385 | 8.89109 |
| Spring | 5.23993 | 7.34912 | 8.00064 | 8.42923 | 8.74089 |
| Summer | 5.82932 | 7.71873 | 8.28468 | 8.70400 | 8.86434 |
| Fall | 5.74090 | 7.66552 | 8.22417 | 8.68366 | 8.83729 |

### 15.3 Interpretation

Spring has the lowest estimated conditional mean first-passage
time at every investigated threshold under both models.

For the interaction model at $K=1000$:

$$
5.82932-5.23993=0.58939
$$

Among trajectories crossing the threshold within the horizon,
the estimated spring conditional mean passage time is
approximately 0.59 time units earlier than summer's.

Higher thresholds produce later conditional mean crossing times
in all investigated configurations.

However, conditional passage-time statistics are calculated
from different sets of successful trajectories.

They should not be interpreted as proof that every corresponding
trajectory crosses earlier under one configuration.

The unrestricted expected first-passage time is not estimated.

---

## 16. Paired Seasonal Interaction Effects

For season $s$ and trajectory $j$, define:

$$
D_{s,j}
=
I_{T,s,j}^{(\gamma)}
-
I_{T,s,j}^{(0)}
$$

The estimated mean interaction effect is:

$$
\overline D_s
=
\frac{1}{N}
\sum_{j=1}^N D_{s,j}
$$

The estimated standard error is:

$$
SE(\overline D_s)
=
\frac{s_{D_s}}{\sqrt N}
$$

An approximate pointwise 95% confidence interval is:

$$
\overline D_s
\pm1.96SE(\overline D_s)
$$

### 16.1 Terminal intensity interaction results

| Season | Mean difference | Relative change | 95% CI |
|---|---:|---:|---|
| Winter | 0.67189 | 0.07255% | [-1.42971, 2.77350] |
| Spring | 16.34830 | 1.53849% | [13.63865, 19.05795] |
| Summer | 7.88709 | 0.86197% | [5.78879, 9.98539] |
| Fall | 7.12458 | 0.76052% | [4.96127, 9.28788] |

Spring exhibits the largest estimated positive arithmetic
mean interaction effect.

Winter is the only seasonal scenario whose approximate
paired mean-difference confidence interval includes zero.

Consequently, the winter mean interaction effect cannot be
clearly distinguished from zero using this interval.

This does not imply that the interaction has no effect on
winter intensity variability or upper quantiles.

### 16.2 Fraction of trajectories with increased intensity

| Season | Fraction increased |
|---|---:|
| Winter | 46.67% |
| Spring | 50.92% |
| Summer | 50.81% |
| Fall | 49.71% |

Approximately half of the paired seasonal trajectories
experience increased terminal intensity.

In winter and fall, fewer than half experience an increase,
despite the positive estimated arithmetic mean difference.

This illustrates how the magnitude of changes matters in
addition to their frequency.

### 16.3 Comparison with Q2

Under the original windy and dry baseline:

- Mean interaction effect: +111.58 kW/m.
- Relative mean increase: 7.51%.
- Fraction increased: 72.95%.

The seasonal interaction effects are considerably smaller.

The original baseline initial conditions were:

$$
z_0^W=1.5
$$

$$
z_0^H=-1.0
$$

Their product is:

$$
z_0^Wz_0^H=-1.5
$$

The seasonal anomaly products have substantially smaller
absolute magnitudes.

This helps explain why the expected log-interaction effects
are smaller in the seasonal experiment.

---

## 17. Threshold-Dependent Interaction Effects

For each season and threshold, define:

$$
\Delta p_s(K)
=
p_s^{(\gamma)}(K)-p_s^{(0)}(K)
$$

For paired trajectories, define:

$$
D_{s,j}(K)
=
\mathbf1\{\tau_{K,s,j}^{(\gamma)}\leq T\}
-
\mathbf1\{\tau_{K,s,j}^{(0)}\leq T\}
$$

The paired estimator is:

$$
\widehat{\Delta p_s}(K)
=
\frac{1}{N}
\sum_{j=1}^N D_{s,j}(K)
$$

### 17.1 Results in percentage points

| Season | 1000 | 2000 | 3000 | 4000 | 5000 |
|---|---:|---:|---:|---:|---:|
| Winter | -1.65 | +0.37 | +0.35 | +0.20 | +0.15 |
| Spring | -0.87 | +0.39 | +0.71 | +0.50 | +0.25 |
| Summer | -1.08 | +0.53 | +0.49 | +0.24 | +0.14 |
| Fall | -1.25 | +0.53 | +0.39 | +0.33 | +0.17 |

### 17.2 Lower-threshold interaction effect

At $K=1000$, introducing the interaction reduces the estimated
crossing probability in all four seasonal scenarios.

Winter exhibits the largest absolute reduction:

$$
41.36\%-43.01\%=-1.65\text{ pp}
$$

### 17.3 Higher-threshold interaction effect

At every investigated threshold of 2000 kW/m or higher,
the interaction increases the estimated crossing probability.

For spring at $K=3000$:

$$
5.49\%-4.78\%=0.71\text{ pp}
$$

The relative increase is:

$$
\frac{0.0549-0.0478}{0.0478}
\times100
\approx14.85\%
$$

At $K=5000$, spring's relative increase is approximately 32.89%.

### 17.4 Mathematical interpretation

The interaction changes the distribution of the accumulated
environmental drift.

The observed increases in variability and upper quantiles are
consistent with increased crossing risk at higher thresholds.

However, terminal variability alone does not mathematically
determine first-passage probabilities.

The pathwise interaction effects depend on the complete
evolution of intensity.

The results establish an empirical change in the direction
of the interaction effect between the investigated thresholds
of 1000 and 2000 kW/m.

They do not identify an exact transition threshold or prove
that the same pattern holds for every possible configuration.

### 17.5 Comparison with Q3

In Q3, the interaction increased crossing probabilities at
every investigated threshold.

In Q5, the interaction reduces crossing probabilities at the
lowest investigated threshold while increasing them at the
higher investigated thresholds.

Therefore, the interaction does not have a universally
positive effect on wildfire threshold crossing risk.

Its effect depends on the initial environmental conditions
and the selected threshold.

---

## 18. Figure Analysis

Seven figures are generated for the seasonal investigation.

### 18.1 Figure 1 — Historical Seasonal Environmental Conditions

The figure displays standardized initial wind and humidity
anomalies for the four seasonal scenarios.

Main observations:

- Winter has the highest positive wind anomaly.
- Winter also has the highest positive humidity anomaly.
- Spring combines positive wind and negative humidity.
- Summer combines negative wind and positive humidity.
- Fall has moderately positive wind and humidity anomalies.

The horizontal zero line represents the historical dataset mean.

Spring has the most favorable initial combination for intensity
growth among the four scenarios under the selected coefficients.

The figure explains the initial conditions used throughout Q5.

### 18.2 Figure 2 — Seasonal Terminal Intensity Statistics

The figure contains four panels:

1. Mean terminal intensity.
2. Standard deviation.
3. Median terminal intensity.
4. 97.5th percentile.

Spring has the highest value in all four panels under both models.

The interaction increases empirical arithmetic means,
standard deviations and upper quantiles in every season.

However, it decreases the empirical median in every season.

The differences in mean intensity are relatively small compared
with the differences in variability and upper quantiles.

This demonstrates why expected intensity alone is insufficient
for characterizing the interaction effect.

### 18.3 Figure 3 — Seasonal Terminal Intensity Distributions

The figure displays the empirical terminal intensity
distributions for each season.

The additive and interaction histograms overlap substantially.

All distributions are visually right-skewed.

Most simulated outcomes are concentrated at lower and
moderate intensity levels, while a smaller proportion reach
much larger intensity values.

Spring's distribution is shifted toward higher intensity
relative to the other seasonal scenarios.

The interaction distributions exhibit higher empirical
upper quantiles, although the histogram differences are subtle.

The histogram does not establish the exact distributional
family of the interaction model.

### 18.4 Figure 4 — Seasonal Threshold Crossing Probabilities

The figure compares threshold crossing probabilities for
the additive and interaction models.

Spring consistently exhibits the highest estimated crossing
probability at every investigated threshold.

All curves decrease as the threshold increases.

This agrees with the mathematical nesting property of
threshold crossing events.

The seasonal differences are most visible at lower thresholds
in absolute percentage-point terms.

At higher thresholds, crossing events become less frequent.

The interaction decreases crossing probability at the
lowest investigated threshold but increases it at the
higher investigated thresholds.

### 18.5 Figure 5 — Seasonal Conditional First-Passage Times

The figure compares conditional mean first-passage times
across seasons and thresholds.

Spring generally separates from the other seasonal curves
by exhibiting earlier conditional crossing times.

All curves increase over the investigated thresholds.

The conditional mean crossing times approach the simulation
horizon for the highest thresholds.

This reflects the observed timing of successful crossings
within the finite simulation window.

The figure does not include trajectories that fail to cross
the corresponding threshold.

### 18.6 Figure 6 — Paired Seasonal Interaction Effects

The left panel displays paired arithmetic mean terminal
intensity differences and approximate 95% confidence intervals.

Spring has the largest positive estimated mean difference.

Summer and fall also exhibit positive mean differences whose
pointwise confidence intervals exclude zero.

Winter's confidence interval includes zero.

The right panel displays paired threshold crossing probability
differences in percentage points.

All four seasonal curves are negative at 1000 kW/m.

All four curves become positive at 2000 kW/m and remain
positive at the other investigated thresholds.

This is one of the main findings of Q5.

The interaction can reduce moderate-threshold crossing risk
while increasing higher-threshold crossing risk under the
same seasonal configuration.

### 18.7 Figure 7 — Analytical and Numerical Consistency

The left panel displays:

$$
\widehat m_s-m_s^{\mathrm{theory}}
$$

for each season and model.

The right panel displays:

$$
\widehat{\Delta m_s}-\Delta m_s^{\mathrm{theory}}
$$

for each season.

The horizontal zero line represents exact agreement between
the numerical estimate and analytical expectation.

All displayed approximate confidence intervals contain zero.

All standardized validation errors have absolute values
below one.

The figure provides evidence of consistency between the
analytical derivations and numerical implementation.

It does not establish empirical calibration or predictive
accuracy against observed wildfire intensity.

---

## 19. Historical Spread Rate Versus Simulated Intensity

The historical mean spread rate ranking is:

$$
\text{Summer}>
\text{Spring}>
\text{Fall}>
\text{Winter}
$$

However, the simulated mean intensity ranking is:

$$
\text{Spring}>
\text{Fall}>
\text{Winter}>
\text{Summer}
$$

| Season | Historical spread | Simulated interaction mean |
|---|---:|---:|
| Winter | 0.353 | 926.73 |
| Spring | 0.893 | 1078.97 |
| Summer | 1.044 | 922.89 |
| Fall | 0.431 | 943.93 |

This discrepancy should be discussed rather than ignored.

### 19.1 Different physical quantities

Spread rate measures the rate at which a wildfire propagates.

Fire intensity measures energy release per unit time per
unit length of the fire front.

A physical relationship is:

$$
I=300CR
$$

where:

- $I$ is intensity in kW/m.
- $C$ is fuel consumption in kg/m².
- $R$ is spread rate in m/min.

Without fuel consumption measurements, historical intensity
cannot be reconstructed from spread rate alone.

### 19.2 Omitted environmental factors

The model does not explicitly incorporate:

- Temperature.
- Seasonal fuel moisture.
- Fuel consumption.
- Precipitation.
- Vegetation changes.
- Topography.
- Fire suppression.

These factors may contribute to historical seasonal differences.

### 19.3 Controlled experimental interpretation

Q5 investigates the isolated effect of seasonal initial wind
and humidity conditions while holding other model parameters fixed.

It does not attempt to reproduce the complete historical
seasonal wildfire intensity distribution.

Therefore, the simulated spring ranking should not be
interpreted as evidence that spring wildfires are universally
more intense than summer wildfires in Alberta.

---

## 20. Mathematical Consistency Checks

The seasonal experiment verifies selected mathematical
properties of the model.

These include:

- Positive simulated wildfire intensity.
- Consistent seasonal simulation grids.
- Valid threshold crossing probabilities.
- Non-increasing crossing probabilities as thresholds increase.
- Agreement between paired and marginal probability differences.
- Consistency of analytical and numerical log-intensity expectations.

The experiment reports:

    All seasonal mathematical consistency checks passed.

These checks support the internal consistency of the experiment.

---

## 21. Automated Testing

The seasonal test suite is implemented in:

    tests/test_seasonal.py

Ten automated tests are included.

| Test | Purpose |
|---|---|
| test_seasonal_data | Validate seasonal historical data |
| test_standardization_at_mean | Confirm zero anomaly at the mean |
| test_standardization_one_std | Confirm standardization by one SD |
| test_seasonal_anomaly_signs | Verify seasonal anomaly directions |
| test_initial_interaction_contribution | Validate the interaction term |
| test_expected_log_intensity_at_zero | Verify the initial log-intensity |
| test_additive_theoretical_consistency | Compare additive analytical formulas |
| test_analytical_interaction_effect | Validate expected log-interaction effects |
| test_additive_log_variance_constant | Verify identical seasonal log variance |
| test_theoretical_seasonal_ranking | Verify theoretical mean ranking |

### 21.1 Seasonal test results

Command:

    python -m pytest tests/test_seasonal.py -v

Result:

    10 passed in 1.49s

### 21.2 Complete project test results

Command:

    python -m pytest -q

Result:

    ...................... [100%]
    22 passed in 1.74s

All current project tests pass.

---

## 22. Implementation Documentation

### 22.1 `src/seasonal.py`

Responsible for preparing seasonal scenarios and calculating
their analytical properties.

Main responsibilities:

- Store historical seasonal summary statistics.
- Standardize seasonal environmental conditions.
- Construct seasonal initial conditions.
- Calculate initial growth rates.
- Calculate analytical additive terminal statistics.
- Calculate expected seasonal log-intensity.
- Calculate expected log-interaction effects.

This module implements the mathematical foundation of Q5.

### 22.2 `src/theoretical.py`

Provides reusable analytical calculations developed in Q2
and subsequent research questions.

Relevant calculations include:

- Integrated OU means.
- Integrated OU variances.
- Additive log-intensity moments.
- Lognormal terminal intensity statistics.
- Expected interaction contributions.

### 22.3 `src/monte_carlo.py`

Provides the vectorized stochastic simulation engine.

For Q5, the initial wind and humidity anomalies and the
interaction coefficient are varied between configurations.

The remaining dynamic parameters are held fixed.

### 22.4 `src/stopping_times.py`

Provides reusable stopping-time calculations.

Responsibilities include:

- First-passage time estimation.
- Threshold crossing probability estimation.
- Conditional passage-time statistics.
- Monte Carlo uncertainty calculations.

### 22.5 `experiments/seasonal_comparison.py`

Conducts the complete seasonal Monte Carlo experiment.

Responsibilities include:

- Simulate all eight configurations.
- Calculate terminal intensity statistics.
- Calculate threshold crossing statistics.
- Calculate conditional first-passage times.
- Perform paired interaction comparisons.
- Validate analytical expectations.
- Save numerical results.

### 22.6 `experiments/plot_seasonal.py`

Generates the seven seasonal figures using the saved
numerical experiment results.

Separating simulation and visualization allows figures
to be regenerated without repeating the full experiment.

### 22.7 `tests/test_seasonal.py`

Contains automated tests for the seasonal mathematical
calculations and scenario construction.

---

## 23. Reproducing Q5

Run all commands from the project root directory.

### Step 1 — Analytical seasonal calculations

    python -m src.seasonal

This displays:

- Historical seasonal data.
- Standardized seasonal conditions.
- Initial growth rates.
- Theoretical additive statistics.
- Theoretical interaction effects.

### Step 2 — Seasonal Monte Carlo experiment

    python -m experiments.seasonal_comparison

This simulates:

- Four seasonal scenarios.
- Two models per season.
- 10,000 trajectories per configuration.
- 80,000 total trajectories.
- 40 threshold evaluations.

### Step 3 — Generate figures

    python -m experiments.plot_seasonal

This generates the seven seasonal visualization figures.

### Step 4 — Seasonal automated tests

    python -m pytest tests/test_seasonal.py -v

### Step 5 — Complete project tests

    python -m pytest -q

---

## 24. Limitations

### 24.1 Seasonal initial conditions only

Seasonal differences are represented through deterministic
initial environmental anomalies.

All seasons share the same long-term OU equilibrium and
dynamic parameters.

The model therefore investigates transient seasonal effects
rather than persistent seasonal environmental regimes.

### 24.2 Historical sample imbalance

Winter contains only 273 observations, compared with
12,492 summer observations.

The reliability of historical seasonal averages may differ
between seasons.

Uncertainty in the seasonal averages is not propagated
through the current simulation.

### 24.3 Within-season variability

Each seasonal scenario begins at its historical mean
wind and humidity conditions.

The model does not sample initial conditions from the full
historical within-season environmental distribution.

### 24.4 Uncalibrated dynamic parameters

The OU and intensity parameters are illustrative.

They are not estimated from longitudinal seasonal
environmental and wildfire intensity measurements.

### 24.5 Independent Brownian drivers

Wind, humidity and intensity Brownian drivers are assumed
independent.

Actual environmental fluctuations may exhibit dependence.

This assumption is investigated in Q6.

### 24.6 Numerical discretization

The environmental OU transitions are exact at grid points.

The intensity drift integral is approximated numerically.

Threshold crossings between grid points may be missed.

### 24.7 Rare threshold events

The highest thresholds are crossed by relatively few
trajectories.

For example, under the interaction model, only 53 summer
trajectories cross 5000 kW/m.

Consequently, rare-event estimates have greater relative
Monte Carlo uncertainty.

### 24.8 Confidence intervals

The reported confidence intervals quantify finite
Monte Carlo sampling uncertainty.

They do not account for parameter uncertainty, historical
measurement errors, model misspecification or numerical
discretization bias.

They are pointwise rather than simultaneous intervals.

### 24.9 Physical environmental bounds

Gaussian OU processes are unbounded.

Consequently, the simulated environmental anomalies can
correspond to physically unrealistic wind or humidity values.

### 24.10 Simulation horizon

The simulation horizon and dynamic parameters are illustrative.

Conditional first-passage times are expressed in model
time units rather than calibrated physical hours or days.

---

## 25. Final Answer to Research Question 5

Seasonal differences in initial wind and humidity conditions
produce measurable differences in the simulated stochastic
behavior of wildfire intensity.

Under the selected model parameters, spring exhibits the
highest initial growth rate, expected terminal intensity,
upper intensity quantiles and threshold crossing probabilities.

Spring also exhibits the lowest estimated conditional mean
first-passage times across the investigated thresholds.

The analytical additive distribution predicts the seasonal
mean intensity ranking:

$$
\text{Spring}>
\text{Fall}>
\text{Winter}>
\text{Summer}
$$

The Monte Carlo experiments reproduce this ranking under
both the additive and interaction models.

Introducing the environmental interaction increases empirical
standard deviations and upper quantiles in every seasonal
scenario, while reducing empirical medians.

The interaction also produces threshold-dependent effects.

At 1000 kW/m, the interaction reduces estimated crossing
probabilities in every seasonal scenario.

At the investigated thresholds of 2000 kW/m and above,
the interaction increases estimated crossing probabilities.

Therefore, the interaction cannot be characterized simply
as increasing or decreasing wildfire risk.

Its effect depends on the environmental initial conditions,
the intensity threshold and the distributional properties
of the stochastic system.

The analytical expected log-intensity and expected
log-interaction effects agree reasonably well with their
Monte Carlo counterparts.

These findings demonstrate the importance of investigating
complete probability distributions and stopping-time behavior
rather than relying exclusively on arithmetic mean intensity.

The results should be interpreted as controlled mathematical
comparisons of seasonal initial environmental conditions,
not calibrated predictions of actual seasonal wildfire behavior.

---

## 26. Connection to Research Question 6

The analytical derivations in Q5 rely on:

$$
E[Z_t^WZ_t^H]
=
E[Z_t^W]E[Z_t^H]
$$

This equality follows from the independence of the
environmental processes under the current assumptions.

If wind and humidity become dependent, then:

$$
E[Z_t^WZ_t^H]
=
E[Z_t^W]E[Z_t^H]
+
\operatorname{Cov}(Z_t^W,Z_t^H)
$$

Consequently, the expected interaction contribution becomes:

$$
\begin{aligned}
E\left[
\gamma\int_0^T Z_t^WZ_t^Hdt
\right]
={}&
\gamma\int_0^T
E[Z_t^W]E[Z_t^H]dt\\
&+
\gamma\int_0^T
\operatorname{Cov}(Z_t^W,Z_t^H)dt.
\end{aligned}
$$

This introduces an additional covariance contribution
to the expected log-intensity.

Dependence also changes the variance of the additive
log-intensity distribution.

Under correlated environmental processes:

$$
\begin{aligned}
\operatorname{Var}(\log I_T^{(0)})
={}&
\beta_W^2\operatorname{Var}(J_T^W)\\
&+\beta_H^2\operatorname{Var}(J_T^H)\\
&-2\beta_W\beta_H
\operatorname{Cov}(J_T^W,J_T^H)\\
&+\sigma_I^2T.
\end{aligned}
$$

Q6 will investigate these additional covariance terms
and determine how environmental dependence changes
terminal intensity distributions and threshold crossing risk.

---

## 27. Mathematical Model References

| Component | Reference |
|---|---|
| Environmental standardization | Equation 2.1 |
| Joint state vector | Equation 2.2 |
| Environmental OU processes | Equations 2.3–2.4 |
| Additive growth function | Equation 2.5 |
| Interaction growth function | Equation 2.6 |
| Intensity SDE | Equation 2.7 |
| Wind OU solution | Equations 3.1–3.2 |
| Wind mean and variance | Equations 3.3–3.4 |
| Humidity OU solution | Equations 3.5–3.6 |
| Humidity mean and variance | Equations 3.7–3.8 |
| Intensity growth function | Equation 3.9 |
| Marginal wind effect | Equation 3.10 |
| Intensity SDE | Equation 3.11 |
| Itô transformation | Equations 3.12–3.13 |
| Integrated log-intensity | Equation 3.14 |
| Exponential intensity solution | Equations 3.15–3.16 |
| Stopping times | Section 4.4 |

Additional mathematical results developed during Q5 include:

- Seasonal integrated OU expectations.
- Seasonal additive lognormal distributions.
- Expected seasonal log-interaction effects.
- Paired seasonal interaction comparisons.
- Analytical and numerical consistency validation.

---

## 28. Completion Checklist

### Historical data preparation

- [x] Obtain historical seasonal statistics.
- [x] Validate seasonal observation counts.
- [x] Standardize seasonal wind conditions.
- [x] Standardize seasonal humidity conditions.
- [x] Construct four seasonal scenarios.

### Mathematical analysis

- [x] Derive seasonal environmental expectations.
- [x] Calculate initial seasonal growth rates.
- [x] Derive expected seasonal log-intensity.
- [x] Derive additive seasonal terminal distributions.
- [x] Derive expected seasonal interaction effects.
- [x] Establish constant additive log variance.

### Numerical experiments

- [x] Simulate four seasonal additive configurations.
- [x] Simulate four seasonal interaction configurations.
- [x] Generate 80,000 total trajectories.
- [x] Calculate terminal intensity statistics.
- [x] Calculate threshold crossing probabilities.
- [x] Calculate conditional first-passage times.
- [x] Calculate paired terminal interaction effects.
- [x] Calculate paired threshold interaction effects.

### Validation

- [x] Validate analytical expected log-intensity.
- [x] Validate analytical interaction expectations.
- [x] Verify mathematical consistency properties.
- [x] Implement ten seasonal automated tests.
- [x] Pass all 22 project automated tests.

### Visualization and documentation

- [x] Generate seven seasonal figures.
- [x] Interpret all seven figures.
- [x] Discuss historical spread-rate differences.
- [x] Discuss mathematical and physical limitations.
- [x] Update general README.
- [x] Complete Q5 research question documentation.

---

## 29. Status

Q5 is completed.

The seasonal mathematical derivations, numerical experiments,
analytical validation, visualization and automated testing
have been completed successfully.

The next research question investigates environmental
dependence and alternative stochastic assumptions.