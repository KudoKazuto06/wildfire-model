import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.stats import lognorm, norm, skew, probplot

from src.config import GAMMA, I0, T, SEED, N_PATHS

from src.monte_carlo import simulate_many

from src.theoretical import (theoretical_intensity_moments, 
                             theoretical_interaction_log_effect)


# ==========================================================
# 1. EMPIRICAL DISTRIBUTION STATISTICS
# ==========================================================

def calculate_statistics(samples):
    """
    Calculate empirical statistics of future wildfire intensity.

    Includes:
    - Mean and variance
    - Standard deviation
    - Median
    - Quantiles
    - Skewness
    - Log-intensity statistics
    """

    log_samples = np.log(samples)

    return {
        "mean": np.mean(samples),
        "variance": np.var(samples, ddof=1),
        "std": np.std(samples, ddof=1),
        "median": np.median(samples),

        "q025": np.quantile(samples, 0.025),
        "q975": np.quantile(samples, 0.975),

        "skewness": skew(samples, bias=False),

        "log_mean": np.mean(log_samples),
        "log_variance": np.var(log_samples, ddof=1),
        "log_skewness": skew(log_samples, bias=False)
    }


# ==========================================================
# 2. THEORETICAL ADDITIVE MODEL STATISTICS
# ==========================================================

def calculate_theoretical_statistics():
    """
    Calculate theoretical statistics for gamma = 0.

    Under the additive model:
        ln(I_T) ~ Normal(m_T, v_T)
        I_T ~ Lognormal(m_T, v_T)
    """

    stats = theoretical_intensity_moments()

    m = stats["log_mean"]
    v = stats["log_variance"]

    # Theoretical lognormal quantiles
    stats["q025"] = np.exp(
        m + np.sqrt(v) * norm.ppf(0.025)
    )

    stats["q975"] = np.exp(
        m + np.sqrt(v) * norm.ppf(0.975)
    )

    # Theoretical lognormal skewness
    stats["skewness"] = (
        np.exp(v) + 2
    ) * np.sqrt(np.exp(v) - 1)

    # Gaussian log-intensity has zero skewness
    stats["log_skewness"] = 0.0

    return stats


# ==========================================================
# 3. PRINT COMPARISON TABLE
# ==========================================================

def print_statistics(theory, additive, interaction):
    """
    Compare theoretical additive statistics with
    both Monte Carlo models.
    """

    print("\nFUTURE INTENSITY DISTRIBUTION")
    print("=" * 80)

    print(
        f"{'Statistic':<18}"
        f"{'Theory (g=0)':>20}"
        f"{'MC (g=0)':>20}"
        f"{'MC (g!=0)':>20}"
    )

    print("-" * 80)

    keys = [
        "mean",
        "variance",
        "std",
        "median",
        "q025",
        "q975",
        "skewness",
        "log_mean",
        "log_variance",
        "log_skewness"
    ]

    for key in keys:

        print(
            f"{key:<18}"
            f"{theory[key]:>20.5f}"
            f"{additive[key]:>20.5f}"
            f"{interaction[key]:>20.5f}"
        )

    print("=" * 80)


# ==========================================================
# 4. PAIRED INTERACTION ANALYSIS
# ==========================================================

def analyze_interaction(additive, interaction):
    """
    Compare corresponding trajectories generated
    using common random numbers.
    """

    differences = interaction - additive

    mean_difference = np.mean(differences)

    std_difference = np.std(
        differences,
        ddof=1
    )

    standard_error = (
        std_difference / np.sqrt(len(differences))
    )

    # Approximate 95% Monte Carlo confidence interval
    ci_lower = mean_difference - 1.96 * standard_error

    ci_upper = mean_difference + 1.96 * standard_error

    # Percentage change in the expected intensity
    percentage_change = (
        mean_difference / np.mean(additive)
    ) * 100

    # Fraction of paired trajectories where interaction
    # produces greater final intensity.
    probability_increase = np.mean(
        differences > 0
    )

    print("\nPAIRED INTERACTION ANALYSIS")
    print("=" * 55)

    print(
        f"Mean difference: "
        f"{mean_difference:.5f} kW/m"
    )

    print(
        f"Standard error: "
        f"{standard_error:.5f} kW/m"
    )

    print(
        f"95% confidence interval: "
        f"[{ci_lower:.5f}, {ci_upper:.5f}]"
    )

    print(
        f"Change in mean intensity: "
        f"{percentage_change:.2f}%"
    )

    print(
        f"Fraction with increased intensity: "
        f"{probability_increase:.4f}"
    )

    return differences


# ==========================================================
# 5. DISTRIBUTION VISUALIZATION
# ==========================================================

def plot_distributions(
    additive,
    interaction,
    theory,
    output_dir
):
    """
    Compare empirical intensity distributions against
    the theoretical additive lognormal distribution.
    """

    m = theory["log_mean"]
    v = theory["log_variance"]

    # Restrict visualization to the central 99.5% range
    # to prevent extreme outcomes from dominating the plot.
    # All observations remain included in the statistics.

    x_max = max(
        np.quantile(additive, 0.995),
        np.quantile(interaction, 0.995)
    )

    x = np.linspace(1e-6, x_max, 1000)

    # scipy lognormal parameterization:
    # shape = sqrt(v)
    # scale = exp(m)

    theoretical_pdf = lognorm.pdf(
        x,
        s=np.sqrt(v),
        scale=np.exp(m)
    )

    bins = np.linspace(0, x_max, 81)

    # Calculate histogram densities using the complete
    # sample size so omitted tail observations do not
    # artificially inflate the displayed density.

    additive_counts, _ = np.histogram(
        additive,
        bins=bins
    )

    interaction_counts, _ = np.histogram(
        interaction,
        bins=bins
    )

    widths = np.diff(bins)

    additive_density = (
        additive_counts / (len(additive) * widths)
    )

    interaction_density = (
        interaction_counts / (len(interaction) * widths)
    )

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.stairs(
        additive_density,
        bins,
        fill=True,
        alpha=0.35,
        label="Monte Carlo: gamma = 0"
    )

    ax.stairs(
        interaction_density,
        bins,
        fill=True,
        alpha=0.35,
        label=f"Monte Carlo: gamma = {GAMMA}"
    )

    ax.plot(
        x,
        theoretical_pdf,
        color="black",
        linewidth=2,
        linestyle="--",
        label="Theoretical lognormal"
    )

    ax.axvline(
        I0,
        color="gray",
        linestyle=":",
        label="Initial intensity"
    )

    ax.set_xlabel("Final wildfire intensity (kW/m)")
    ax.set_ylabel("Probability density")

    ax.set_title(
        f"Future Wildfire Intensity Distribution (T = {T})"
    )

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        output_dir / "intensity_distribution.png",
        dpi=200,
        bbox_inches="tight"
    )


# ==========================================================
# 6. LOG-INTENSITY DISTRIBUTION
# ==========================================================

def plot_log_distributions(
    additive,
    interaction,
    theory,
    output_dir
):
    """
    Investigate whether log-intensity follows
    a Gaussian distribution.
    """

    log_additive = np.log(additive)

    log_interaction = np.log(interaction)

    m = theory["log_mean"]
    v = theory["log_variance"]

    x_min = min(
        log_additive.min(),
        log_interaction.min()
    )

    x_max = max(
        log_additive.max(),
        log_interaction.max()
    )

    x = np.linspace(x_min, x_max, 1000)

    bins = np.linspace(x_min, x_max, 81)

    theoretical_pdf = norm.pdf(
        x,
        loc=m,
        scale=np.sqrt(v)
    )

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.hist(
        log_additive,
        bins=bins,
        density=True,
        alpha=0.4,
        label="Monte Carlo: gamma = 0"
    )

    ax.hist(
        log_interaction,
        bins=bins,
        density=True,
        alpha=0.4,
        label=f"Monte Carlo: gamma = {GAMMA}"
    )

    ax.plot(
        x,
        theoretical_pdf,
        color="black",
        linestyle="--",
        linewidth=2,
        label="Theoretical Gaussian"
    )

    ax.set_xlabel("Log wildfire intensity")
    ax.set_ylabel("Probability density")
    ax.set_title("Distribution of Log-Intensity")

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        output_dir / "log_intensity_distribution.png",
        dpi=200,
        bbox_inches="tight"
    )


# ==========================================================
# 7. NORMAL Q-Q PLOTS
# ==========================================================

def plot_qq(additive, interaction, output_dir):
    """
    Compare empirical log-intensity quantiles with
    theoretical Gaussian quantiles.
    """

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12, 5)
    )

    probplot(
        np.log(additive),
        dist="norm",
        plot=axes[0]
    )

    axes[0].set_title(
        "Additive Model: gamma = 0"
    )

    probplot(
        np.log(interaction),
        dist="norm",
        plot=axes[1]
    )

    axes[1].set_title(
        f"Interaction Model: gamma = {GAMMA}"
    )

    for ax in axes:
        ax.grid(alpha=0.3)

    fig.suptitle(
        "Normal Q-Q Plots of Log-Intensity"
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / "log_intensity_qq.png",
        dpi=200,
        bbox_inches="tight"
    )


# ==========================================================
# 8. EMPIRICAL CUMULATIVE DISTRIBUTION
# ==========================================================

def plot_cdf(
    additive,
    interaction,
    theory,
    output_dir
):
    """
    Compare the probability that final wildfire intensity
    is below a specified intensity value.
    """

    additive_sorted = np.sort(additive)

    interaction_sorted = np.sort(interaction)

    additive_cdf = (
        np.arange(1, len(additive) + 1)
        / len(additive)
    )

    interaction_cdf = (
        np.arange(1, len(interaction) + 1)
        / len(interaction)
    )

    x_max = max(
        np.quantile(additive, 0.995),
        np.quantile(interaction, 0.995)
    )

    x = np.linspace(0, x_max, 1000)

    theoretical_cdf = lognorm.cdf(
        x,
        s=np.sqrt(theory["log_variance"]),
        scale=np.exp(theory["log_mean"])
    )

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.step(
        additive_sorted,
        additive_cdf,
        where="post",
        label="Monte Carlo: gamma = 0"
    )

    ax.step(
        interaction_sorted,
        interaction_cdf,
        where="post",
        label=f"Monte Carlo: gamma = {GAMMA}"
    )

    ax.plot(
        x,
        theoretical_cdf,
        color="black",
        linestyle="--",
        linewidth=2,
        label="Theoretical lognormal CDF"
    )

    ax.set_xlim(0, x_max)
    ax.set_ylim(0, 1.02)

    ax.set_xlabel("Final wildfire intensity (kW/m)")
    ax.set_ylabel("Cumulative probability")

    ax.set_title(
        "Cumulative Distribution of Future Intensity"
    )

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        output_dir / "intensity_cdf.png",
        dpi=200,
        bbox_inches="tight"
    )


# ==========================================================
# 9. MAIN EXPERIMENT
# ==========================================================

def main():

    output_dir = Path("results/figures")

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # ------------------------------------------------------
    # Simulate additive model
    # ------------------------------------------------------

    print("Simulating additive model...")

    _, additive_paths = simulate_many(
        n_paths=N_PATHS,
        gamma=0.0,
        seed=SEED
    )

    additive_final = additive_paths[:, -1].copy()

    # Release complete trajectories because this experiment
    # only requires final intensity values.
    del additive_paths

    # ------------------------------------------------------
    # Simulate interaction model
    # ------------------------------------------------------

    print("Simulating interaction model...")

    _, interaction_paths = simulate_many(
        n_paths=N_PATHS,
        gamma=GAMMA,
        seed=SEED
    )

    interaction_final = interaction_paths[:, -1].copy()

    del interaction_paths

    # ------------------------------------------------------
    # Calculate statistical quantities
    # ------------------------------------------------------

    theory = calculate_theoretical_statistics()

    additive_stats = calculate_statistics(
        additive_final
    )

    interaction_stats = calculate_statistics(
        interaction_final
    )

    print_statistics(
        theory,
        additive_stats,
        interaction_stats
    )

    # ------------------------------------------------------
    # Paired comparison
    # ------------------------------------------------------

    analyze_interaction(
        additive_final,
        interaction_final
    )
    
    # ------------------------------------------------------
    # Validate the expected log-intensity interaction effect
    # ------------------------------------------------------

    theoretical_log_effect = (
        theoretical_interaction_log_effect(gamma=GAMMA)
    )

    empirical_log_effect = np.mean(
        np.log(interaction_final)
        - np.log(additive_final)
    )

    print("\nLOG-INTENSITY INTERACTION VALIDATION")
    print("=" * 55)

    print(
        f"Theoretical log effect: "
        f"{theoretical_log_effect:.5f}"
    )

    print(
        f"Empirical log effect:    "
        f"{empirical_log_effect:.5f}"
    )

    print(
        f"Absolute difference:   "
        f"{abs(empirical_log_effect - theoretical_log_effect):.5f}"
    )

    # ------------------------------------------------------
    # Generate figures
    # ------------------------------------------------------

    plot_distributions(
        additive_final,
        interaction_final,
        theory,
        output_dir
    )

    plot_log_distributions(
        additive_final,
        interaction_final,
        theory,
        output_dir
    )

    plot_qq(
        additive_final,
        interaction_final,
        output_dir
    )

    plot_cdf(
        additive_final,
        interaction_final,
        theory,
        output_dir
    )

    print("\nAll distribution figures saved successfully.")

    plt.show()


if __name__ == "__main__":
    main()
    
    
# ==========================================================
# Q2: FUTURE WILDFIRE INTENSITY DISTRIBUTION
# ==========================================================
#
# EXPERIMENT CONFIGURATION
# ----------------------------------------------------------
# Number of trajectories per model: 10000
# Simulation horizon: T = 10
# Time-step size: dt = 0.01
# Initial intensity: 500 kW/m
# Initial wind anomaly: 1.5
# Initial humidity anomaly: -1.0
# Additive model: gamma = 0
# Interaction model: gamma = -0.04
# Random seed: 42
#
# Common random numbers are used to isolate the effect
# of the wind-humidity interaction.
#
# NUMERICAL RESULTS
# ----------------------------------------------------------
# Statistic       Theory         Additive       Interaction
#
# Mean            1478.02852     1485.38922     1596.96718
# Variance        1585062.13     1495024.21     2068770.98
# Std             1258.99251     1222.71183     1438.32228
# Median          1125.16507     1132.65977     1166.65167
# Q025            264.54661      264.02823      277.31807
# Q975            4785.53261     4727.54330     5434.03406
# Skewness        3.17346        2.53473        3.02021
# Log mean        7.02569        7.03215        7.07871
# Log variance    0.54556        0.55024        0.58505
# Log skewness    0.00000       -0.03885        0.08614
#
# PAIRED INTERACTION ANALYSIS
# ----------------------------------------------------------
# Mean difference: 111.57796 kW/m
# Standard error: 2.93821 kW/m
# 95% MC confidence interval: [105.81906, 117.33686]
# Increase in mean intensity: 7.51%
# Increase in variance: 38.38%
# Increase in standard deviation: 17.63%
# Increase in upper 97.5th percentile: 14.94%
# Increase in prediction interval width: 15.53%
# Fraction with increased intensity: 72.95%
#
# ANALYTICAL INTERACTION EFFECT
# ----------------------------------------------------------
# The expected difference in log-intensity is
#
# E[ln(I_gamma) - ln(I_additive)]
#
# = gamma * wind_z0 * humidity_z0
#   / (KAPPA_W + KAPPA_H)
#   * (1 - exp(-(KAPPA_W + KAPPA_H) * T))
#
# Theoretical log effect: approximately 0.04615
# Empirical log effect: approximately 0.04656
#
# This explains the positive expected log-intensity effect
# under the selected initial environmental conditions.
#
# FIGURE ANALYSIS
# ----------------------------------------------------------
# Figure 1: intensity_distribution.png
# - Additive Monte Carlo distribution closely follows the
#   theoretical lognormal density.
# - Both intensity distributions are positively skewed.
# - Interaction increases variability and upper quantiles.
#
# Figure 2: log_intensity_distribution.png
# - Additive log-intensity closely follows a Gaussian.
# - Interaction log-intensity remains approximately Gaussian
#   under the selected parameter configuration.
#
# Figure 3: log_intensity_qq.png
# - Both models show approximately linear normal Q-Q plots.
# - Small departures occur at extreme quantiles.
# - Approximate normality does not establish exact normality.
#
# Figure 4: intensity_cdf.png
# - Additive empirical CDF closely matches theory.
# - Interaction CDF is lower over much of the plotted range,
#   indicating increased terminal exceedance probabilities.
#
# CONCLUSION
# ----------------------------------------------------------
# The additive model produces an exactly lognormal
# continuous-time intensity distribution under the stated
# independence and deterministic initial-state assumptions.
#
# The interaction model is generally not exactly lognormal,
# but remains approximately lognormal in this experiment.
#
# Including the interaction increases expected intensity,
# variability, upper quantiles and prediction interval width.
#
# The interaction increases final intensity in 72.95%
# of paired simulated trajectories.
#
# These conclusions apply to the selected illustrative
# parameters and initial conditions. They are not calibrated
# predictions of real wildfire behavior.
#
# Q2 STATUS
# ----------------------------------------------------------
# [x] Derive integrated OU moments.
# [x] Derive theoretical additive intensity distribution.
# [x] Implement theoretical calculations.
# [x] Implement vectorized Monte Carlo simulation.
# [x] Validate additive simulation against theory.
# [x] Simulate the interaction model.
# [x] Compare empirical distributions.
# [x] Analyze prediction intervals and skewness.
# [x] Perform paired interaction analysis.
# [x] Derive analytical expected log-interaction effect.
# [x] Generate and analyze distribution figures.
# [ ] Confirm analytical interaction validation output.
# [ ] Update README.md.
# [ ] Complete Q2_intensity_distribution.md.
# ==========================================================