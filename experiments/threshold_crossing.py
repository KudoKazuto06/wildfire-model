import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

from src.config import GAMMA, SEED
from src.monte_carlo import simulate_many
from src.stopping_times import (
    first_passage_times,
    threshold_statistics
)


# ============================================================
# EXPERIMENT CONFIGURATION
# ============================================================

N_PATHS = 10000

THRESHOLDS = [1000, 2000, 3000, 4000, 5000]

OUTPUT_DIR = Path("results/figures")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# THRESHOLD ANALYSIS
# ============================================================

def analyze_thresholds(times, paths):
    """
    Calculate threshold crossing statistics for
    every selected intensity threshold.
    """

    results = []

    for threshold in THRESHOLDS:

        stats = threshold_statistics(
            times,
            paths,
            threshold
        )

        results.append(stats)

    return pd.DataFrame(results)


# ============================================================
# PAIRED INTERACTION ANALYSIS
# ============================================================

def paired_analysis(additive_paths, interaction_paths):
    """
    Estimate the effect of the interaction on
    threshold crossing probabilities.

    Common random numbers allow us to compare
    corresponding trajectories directly.
    """

    additive_max = np.max(additive_paths, axis=1)
    interaction_max = np.max(interaction_paths, axis=1)

    results = []

    for threshold in THRESHOLDS:

        additive_crossed = additive_max >= threshold
        interaction_crossed = interaction_max >= threshold

        # Possible differences:
        # +1: Only interaction model crosses.
        #  0: Both models have the same crossing outcome.
        # -1: Only additive model crosses.

        differences = (
            interaction_crossed.astype(float)
            - additive_crossed.astype(float)
        )

        mean_difference = np.mean(differences)

        standard_error = (
            np.std(differences, ddof=1)
            / np.sqrt(len(differences))
        )

        ci_lower = mean_difference - 1.96 * standard_error
        ci_upper = mean_difference + 1.96 * standard_error

        results.append({
            "threshold": threshold,
            "probability_difference": mean_difference,
            "standard_error": standard_error,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "interaction_only": np.mean(differences == 1),
            "additive_only": np.mean(differences == -1)
        })

    return pd.DataFrame(results)


# ============================================================
# FIGURE 1: CROSSING PROBABILITY VS THRESHOLD
# ============================================================

def plot_crossing_probability(additive, interaction):

    fig, ax = plt.subplots(figsize=(10, 6))

    for data, label in [
        (additive, "Additive model"),
        (interaction, "Interaction model")
    ]:

        ax.errorbar(
            data["threshold"],
            data["crossing_probability"],
            yerr=1.96 * data["standard_error"],
            marker="o",
            capsize=4,
            linewidth=2,
            label=label
        )

    ax.set_xlabel("Intensity threshold (kW/m)")
    ax.set_ylabel("Threshold crossing probability")
    ax.set_title("Probability of Reaching an Intensity Threshold")

    ax.set_ylim(0, 1.05)
    ax.set_xticks(THRESHOLDS)

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "threshold_crossing_probability.png",
        dpi=200,
        bbox_inches="tight"
    )


# ============================================================
# FIGURE 2: TERMINAL VS PATHWISE EXCEEDANCE
# ============================================================

def plot_terminal_vs_pathwise(additive, interaction):

    fig, axes = plt.subplots(
        1, 2,
        figsize=(14, 5),
        sharey=True
    )

    for ax, data, title in [
        (axes[0], additive, "Additive model"),
        (axes[1], interaction, "Interaction model")
    ]:

        ax.plot(
            data["threshold"],
            data["crossing_probability"],
            marker="o",
            linewidth=2,
            label="Pathwise crossing"
        )

        ax.plot(
            data["threshold"],
            data["terminal_probability"],
            marker="s",
            linestyle="--",
            linewidth=2,
            label="Terminal exceedance"
        )

        # Highlight the probability of crossing before T
        # but finishing below the threshold.

        ax.fill_between(
            data["threshold"],
            data["terminal_probability"],
            data["crossing_probability"],
            alpha=0.15
        )

        ax.set_title(title)
        ax.set_xlabel("Intensity threshold (kW/m)")
        ax.set_xticks(THRESHOLDS)
        ax.set_ylim(0, 1.05)

        ax.legend()
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Probability")

    fig.suptitle(
        "Terminal Exceedance vs Pathwise Crossing Probability"
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "terminal_vs_pathwise.png",
        dpi=200,
        bbox_inches="tight"
    )


# ============================================================
# FIGURE 3: CUMULATIVE FIRST-PASSAGE PROBABILITY
# ============================================================

def cumulative_crossing_probability(times, passage_times):
    """
    Estimate P(tau_K <= t) at every simulation time.

    Non-crossing trajectories remain in the denominator.
    """

    sorted_passage_times = np.sort(passage_times)

    cumulative_probability = np.searchsorted(
        sorted_passage_times,
        times,
        side="right"
    ) / len(passage_times)

    return cumulative_probability


def plot_cumulative_crossing(
    times,
    additive_paths,
    interaction_paths
):

    fig, axes = plt.subplots(
        1, 2,
        figsize=(14, 5),
        sharey=True
    )

    for ax, paths, title in [
        (axes[0], additive_paths, "Additive model"),
        (axes[1], interaction_paths, "Interaction model")
    ]:

        for threshold in THRESHOLDS:

            passage_times = first_passage_times(
                times,
                paths,
                threshold
            )

            cumulative_probability = (
                cumulative_crossing_probability(
                    times,
                    passage_times
                )
            )

            ax.plot(
                times,
                cumulative_probability,
                linewidth=2,
                label=f"K = {threshold}"
            )

        ax.set_title(title)
        ax.set_xlabel("Time (illustrative units)")
        ax.set_ylim(0, 1.05)

        ax.legend()
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Cumulative crossing probability")

    fig.suptitle("Cumulative First-Passage Probability")

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "cumulative_crossing_probability.png",
        dpi=200,
        bbox_inches="tight"
    )


# ============================================================
# FIGURE 4: CONDITIONAL FIRST-PASSAGE DISTRIBUTIONS
# ============================================================

def plot_first_passage_distributions(
    times,
    additive_paths,
    interaction_paths
):

    selected_thresholds = [1000, 3000, 5000]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(16, 5),
        sharex=True
    )

    bins = np.linspace(times[0], times[-1], 31)

    for ax, threshold in zip(axes, selected_thresholds):

        for paths, label in [
            (additive_paths, "Additive model"),
            (interaction_paths, "Interaction model")
        ]:

            passage_times = first_passage_times(
                times,
                paths,
                threshold
            )

            # Only include trajectories that crossed
            # within the observation horizon.

            successful_times = passage_times[
                np.isfinite(passage_times)
            ]

            if len(successful_times) > 0:

                ax.hist(
                    successful_times,
                    bins=bins,
                    density=True,
                    alpha=0.5,
                    label=label
                )

        ax.set_title(f"K = {threshold} kW/m")
        ax.set_xlabel("First-passage time")
        ax.legend()
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Conditional probability density")

    fig.suptitle(
        "First-Passage Time Distributions Given Crossing"
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "first_passage_distributions.png",
        dpi=200,
        bbox_inches="tight"
    )


# ============================================================
# FIGURE 5: PAIRED INTERACTION EFFECT
# ============================================================

def plot_interaction_effect(paired_results):

    fig, ax = plt.subplots(figsize=(10, 6))

    # Convert probability differences into
    # percentage-point differences.

    differences = (
        paired_results["probability_difference"] * 100
    )

    errors = (
        1.96 * paired_results["standard_error"] * 100
    )

    ax.errorbar(
        paired_results["threshold"],
        differences,
        yerr=errors,
        marker="o",
        linewidth=2,
        capsize=5
    )

    ax.axhline(
        0,
        color="black",
        linestyle="--",
        label="No interaction effect"
    )

    ax.set_xlabel("Intensity threshold (kW/m)")

    ax.set_ylabel(
        "Change in crossing probability (percentage points)"
    )

    ax.set_title("Paired Effect of Wind-Humidity Interaction")

    ax.set_xticks(THRESHOLDS)

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "threshold_interaction_effect.png",
        dpi=200,
        bbox_inches="tight"
    )


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def main():

    print("Simulating additive model...")

    times, additive_paths = simulate_many(
        n_paths=N_PATHS,
        gamma=0.0,
        seed=SEED
    )

    print("Simulating interaction model...")

    interaction_times, interaction_paths = simulate_many(
        n_paths=N_PATHS,
        gamma=GAMMA,
        seed=SEED
    )

    # Both models must use the same time grid.

    np.testing.assert_array_equal(
        times,
        interaction_times
    )

    # Calculate threshold statistics.

    additive_results = analyze_thresholds(
        times,
        additive_paths
    )

    interaction_results = analyze_thresholds(
        times,
        interaction_paths
    )

    paired_results = paired_analysis(
        additive_paths,
        interaction_paths
    )

    # ========================================================
    # MATHEMATICAL CONSISTENCY CHECKS
    # ========================================================

    for results in [additive_results, interaction_results]:

        # Terminal probability cannot exceed
        # pathwise crossing probability.

        assert np.all(
            results["terminal_probability"]
            <= results["crossing_probability"]
        )

        # Crossing probability must be non-increasing
        # as the intensity threshold increases.

        assert np.all(
            np.diff(results["crossing_probability"]) <= 0
        )

    np.testing.assert_allclose(
        paired_results["probability_difference"],
        interaction_results["crossing_probability"]
        - additive_results["crossing_probability"]
    )

    print("\nAll mathematical consistency checks passed.")

    # ========================================================
    # PRINT RESULTS
    # ========================================================

    columns = [
        "threshold",
        "n_crossed",
        "crossing_probability",
        "terminal_probability",
        "standard_error",
        "ci_lower",
        "ci_upper",
        "mean_passage_time",
        "median_passage_time"
    ]

    print("\nADDITIVE MODEL")
    print("=" * 100)

    print(
        additive_results[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}"
        )
    )

    print("\nINTERACTION MODEL")
    print("=" * 100)

    print(
        interaction_results[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}"
        )
    )

    print("\nPAIRED INTERACTION ANALYSIS")
    print("=" * 100)

    print(
        paired_results.to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}"
        )
    )

    # ========================================================
    # SAVE NUMERICAL RESULTS
    # ========================================================

    results_dir = Path("results")

    results_dir.mkdir(parents=True, exist_ok=True)

    additive_results.to_csv(
        results_dir / "threshold_additive.csv",
        index=False
    )

    interaction_results.to_csv(
        results_dir / "threshold_interaction.csv",
        index=False
    )

    paired_results.to_csv(
        results_dir / "threshold_paired_comparison.csv",
        index=False
    )

    # ========================================================
    # GENERATE FIGURES
    # ========================================================

    plot_crossing_probability(
        additive_results,
        interaction_results
    )

    plot_terminal_vs_pathwise(
        additive_results,
        interaction_results
    )

    plot_cumulative_crossing(
        times,
        additive_paths,
        interaction_paths
    )

    plot_first_passage_distributions(
        times,
        additive_paths,
        interaction_paths
    )

    plot_interaction_effect(paired_results)

    print("\nAll threshold crossing results saved successfully.")

    plt.show()


if __name__ == "__main__":
    main()
    

# ============================================================
# Q3 EXPERIMENTAL RESULTS AND INTERPRETATION
# ============================================================
#
# Experiment:
#   Number of trajectories per model = 10000
#   Initial intensity = 500 kW/m
#   Time horizon = 10
#   Time step = 0.01
#   Additive model gamma = 0
#   Interaction model gamma = -0.04
#
# Both models use common random numbers.
#
# ------------------------------------------------------------
# 1. THRESHOLD CROSSING PROBABILITIES
# ------------------------------------------------------------
#
# Threshold    Additive    Interaction    Difference
# 1000         0.73770     0.75590        +0.01820
# 2000         0.29240     0.31970        +0.02730
# 3000         0.12790     0.15410        +0.02620
# 4000         0.05600     0.07800        +0.02200
# 5000         0.02720     0.04080        +0.01360
#
# Crossing probabilities decrease as the threshold increases.
# This agrees with the theoretical monotonicity property.
#
# The interaction increases crossing probabilities at all
# five investigated thresholds.
#
# Relative probability increases:
#   K = 1000:  2.47%
#   K = 2000:  9.34%
#   K = 3000: 20.48%
#   K = 4000: 39.29%
#   K = 5000: 50.00%
#
# The largest absolute increase occurs at K = 2000.
# The largest relative increase occurs at K = 5000.
#
# ------------------------------------------------------------
# 2. TERMINAL VS PATHWISE EXCEEDANCE
# ------------------------------------------------------------
#
# Terminal exceedance probabilities are smaller than
# pathwise crossing probabilities for every threshold.
#
# Mathematical relationship:
#
#   P(I_T >= K) <= P(tau_K <= T)
#
# The difference represents trajectories that reached
# the threshold but finished below it.
#
# For the interaction model at K = 1000:
#
#   Pathwise probability = 0.75590
#   Terminal probability = 0.58030
#   Difference           = 0.17560
#
# Therefore, 1756 out of 10000 trajectories crossed
# the threshold but finished below it.
#
# Terminal intensity analysis alone would miss these
# observed high-intensity events.
#
# ------------------------------------------------------------
# 3. CUMULATIVE FIRST-PASSAGE PROBABILITY
# ------------------------------------------------------------
#
# F_K(t) = P(tau_K <= t)
#
# The cumulative crossing curves increase over time.
# Higher thresholds have lower cumulative probabilities.
#
# Lower thresholds tend to be reached earlier.
# Higher thresholds accumulate crossings later.
#
# The interaction model exhibits earlier accumulation
# of crossing probability under the selected configuration.
#
# ------------------------------------------------------------
# 4. CONDITIONAL FIRST-PASSAGE TIMES
# ------------------------------------------------------------
#
# Threshold    Additive    Interaction
# 1000         3.89479     3.62587
# 2000         6.74927     6.35238
# 3000         7.71501     7.32248
# 4000         8.13120     7.83714
# 5000         8.50261     8.10963
#
# Conditional mean crossing times are lower in the
# interaction model for every investigated threshold.
#
# These statistics estimate:
#
#   E[tau_K | tau_K <= T]
#
# They do not estimate the unrestricted expected
# first-passage time E[tau_K].
#
# The two models condition on potentially different
# successful trajectories. Consequently, the difference
# in conditional means is not a paired timing effect.
#
# ------------------------------------------------------------
# 5. PAIRED INTERACTION ANALYSIS
# ------------------------------------------------------------
#
# Threshold    Difference    95% Monte Carlo CI
# 1000         0.01820       [0.01480, 0.02160]
# 2000         0.02730       [0.02362, 0.03098]
# 3000         0.02620       [0.02289, 0.02951]
# 4000         0.02200       [0.01907, 0.02493]
# 5000         0.01360       [0.01133, 0.01587]
#
# All five pointwise confidence intervals exclude zero.
#
# This provides Monte Carlo statistical evidence of
# a positive interaction effect on threshold crossing
# probabilities under the selected model parameters.
#
# ------------------------------------------------------------
# 6. VALIDATION
# ------------------------------------------------------------
#
# All 10 stopping-time unit tests passed.
# All 12 existing project tests passed.
#
# Mathematical consistency checks passed:
#
#   1. Terminal probability <= crossing probability.
#   2. Crossing probability decreases or remains constant
#      as the threshold increases.
#   3. Paired probability differences agree with the
#      differences between empirical crossing probabilities.
#
# ------------------------------------------------------------
# 7. LIMITATIONS
# ------------------------------------------------------------
#
# Model parameters are illustrative rather than calibrated.
#
# Thresholds are illustrative experimental values, not
# validated operational wildfire danger classifications.
#
# Crossings between numerical time points may be missed.
#
# Trajectories that do not cross within T are right-censored
# and represented computationally using infinity.
#
# Confidence intervals quantify Monte Carlo sampling error,
# not parameter uncertainty or discretization error.
#
# The results are conditional on the selected environmental
# initial conditions and model configuration.
#
# ------------------------------------------------------------
# 8. FINAL CONCLUSION
# ------------------------------------------------------------
#
# The wind-humidity interaction increases the probability
# of reaching high-intensity thresholds and reduces the
# observed conditional mean first-passage times.
#
# The relative effect on crossing probability becomes
# larger at higher intensity thresholds in this experiment.
#
# Pathwise risk analysis provides additional information
# that cannot be obtained from terminal intensity alone.
#
# Q3: Numerical implementation and analysis completed.
#
# ============================================================