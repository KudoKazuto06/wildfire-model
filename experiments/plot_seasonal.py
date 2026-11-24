"""
Q5: Visualization of seasonal wildfire simulations.

Uses previously saved numerical results.
No additional Monte Carlo simulation is required.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


RESULTS_DIR = Path("results")

OUTPUT_DIR = RESULTS_DIR / "figures"

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SEASONS = [
    "Winter",
    "Spring",
    "Summer",
    "Fall",
]

MODELS = [
    "Additive",
    "Interaction",
]

MODEL_COLORS = {
    "Additive": "tab:blue",
    "Interaction": "tab:orange",
}

SEASON_COLORS = {
    "Winter": "tab:blue",
    "Spring": "tab:green",
    "Summer": "tab:red",
    "Fall": "tab:orange",
}


def save_figure(fig, filename):

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / filename,
        dpi=200,
        bbox_inches="tight",
    )
    
def plot_environment(scenarios):

    scenarios = (
        scenarios
        .set_index("season")
        .loc[SEASONS]
    )

    x = np.arange(len(SEASONS))

    width = 0.35

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )

    ax.bar(
        x - width / 2,
        scenarios["wind_z0"],
        width,
        label="Wind anomaly",
        color="tab:blue",
    )

    ax.bar(
        x + width / 2,
        scenarios["humidity_z0"],
        width,
        label="Humidity anomaly",
        color="tab:green",
    )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    ax.set_xticks(x, SEASONS)

    ax.set_ylabel(
        "Standardized initial anomaly"
    )

    ax.set_title(
        "Historical Seasonal Environmental Conditions"
    )

    ax.legend()

    ax.grid(
        axis="y",
        alpha=0.25,
    )

    save_figure(
        fig,
        "seasonal_environment.png",
    )

def plot_terminal_statistics(terminal):

    metrics = {
        "mean": "Mean intensity",
        "std": "Standard deviation",
        "median": "Median intensity",
        "q975": "97.5th percentile",
    }

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13, 9),
    )

    x = np.arange(len(SEASONS))

    width = 0.35

    for ax, (metric, label) in zip(
        axes.flat,
        metrics.items(),
    ):

        for i, model in enumerate(MODELS):

            data = (
                terminal[
                    terminal["model"] == model
                ]
                .set_index("season")
                .loc[SEASONS]
            )

            offset = (
                i - 0.5
            ) * width

            ax.bar(
                x + offset,
                data[metric],
                width,
                label=model,
                color=MODEL_COLORS[model],
            )

        ax.set_xticks(x, SEASONS)

        ax.set_title(label)

        ax.set_ylabel("Intensity (kW/m)")

        ax.grid(
            axis="y",
            alpha=0.25,
        )

        ax.legend()

    fig.suptitle(
        "Seasonal Terminal Intensity Statistics",
        fontsize=16,
    )

    save_figure(
        fig,
        "seasonal_terminal_statistics.png",
    )

def plot_distributions(samples):

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13, 9),
        sharex=True,
        sharey=True,
    )

    # Common plotting range across all seasons.
    upper_limit = max(
        np.quantile(samples[key], 0.995)
        for key in samples.files
    )

    bins = np.linspace(
        0,
        upper_limit,
        75,
    )

    for ax, season in zip(
        axes.flat,
        SEASONS,
    ):

        for model in MODELS:

            values = samples[
                f"{season}_{model}"
            ]

            counts, edges = np.histogram(
                values,
                bins=bins,
            )

            # Normalize by the full sample size,
            # including observations outside
            # the displayed plotting range.
            density = counts / (
                len(values) * np.diff(edges)
            )

            ax.stairs(
                density,
                edges,
                fill=True,
                alpha=0.35,
                label=model,
                color=MODEL_COLORS[model],
            )

        ax.set_title(season)

        ax.set_xlabel(
            "Terminal intensity (kW/m)"
        )

        ax.set_ylabel("Density")

        ax.legend()

        ax.grid(alpha=0.2)

    fig.suptitle(
        "Seasonal Terminal Intensity Distributions",
        fontsize=16,
    )

    save_figure(
        fig,
        "seasonal_distributions.png",
    )
    
def plot_crossing_probabilities(thresholds):

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14, 6),
        sharey=True,
    )

    for ax, model in zip(
        axes,
        MODELS,
    ):

        for season in SEASONS:

            data = thresholds[
                (thresholds["season"] == season)
                &
                (thresholds["model"] == model)
            ].sort_values("threshold")

            ax.errorbar(
                data["threshold"],
                100 * data["crossing_probability"],
                yerr=(
                    100
                    * 1.96
                    * data["standard_error"]
                ),
                marker="o",
                capsize=3,
                label=season,
                color=SEASON_COLORS[season],
            )

        ax.set_title(model)

        ax.set_xlabel(
            "Threshold (kW/m)"
        )

        ax.set_ylabel(
            "Crossing probability (%)"
        )

        ax.legend()

        ax.grid(alpha=0.25)

    fig.suptitle(
        "Seasonal Threshold Crossing Probabilities",
        fontsize=16,
    )

    save_figure(
        fig,
        "seasonal_crossing_probability.png",
    )
    
def plot_first_passage(thresholds):

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14, 6),
        sharey=True,
    )

    for ax, model in zip(
        axes,
        MODELS,
    ):

        for season in SEASONS:

            data = thresholds[
                (thresholds["season"] == season)
                &
                (thresholds["model"] == model)
            ].sort_values("threshold")

            ax.plot(
                data["threshold"],
                data["mean_passage_time"],
                marker="o",
                label=season,
                color=SEASON_COLORS[season],
            )

        ax.set_title(model)

        ax.set_xlabel(
            "Threshold (kW/m)"
        )

        ax.set_ylabel(
            "Conditional mean first-passage time"
        )

        ax.legend()

        ax.grid(alpha=0.25)

    fig.suptitle(
        "Seasonal Conditional First-Passage Times",
        fontsize=16,
    )

    save_figure(
        fig,
        "seasonal_first_passage.png",
    )
    
def plot_interaction_effects(
    paired_terminal,
    paired_thresholds,
):

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15, 6),
    )

    # Panel A: Paired terminal mean differences.
    data = (
        paired_terminal
        .set_index("season")
        .loc[SEASONS]
    )

    x = np.arange(len(SEASONS))

    axes[0].errorbar(
        x,
        data["mean_difference"],
        yerr=1.96 * data["standard_error"],
        fmt="o",
        capsize=5,
        color="tab:blue",
    )

    axes[0].axhline(
        0,
        color="black",
        linestyle="--",
    )

    axes[0].set_xticks(x, SEASONS)

    axes[0].set_title(
        "Effect on Mean Terminal Intensity"
    )

    axes[0].set_ylabel(
        "Mean difference (kW/m)"
    )

    axes[0].grid(alpha=0.25)

    # Panel B: Paired crossing probability differences.
    for season in SEASONS:

        data = paired_thresholds[
            paired_thresholds["season"] == season
        ].sort_values("threshold")

        axes[1].errorbar(
            data["threshold"],
            100 * data["probability_difference"],
            yerr=(
                100
                * 1.96
                * data["standard_error"]
            ),
            marker="o",
            capsize=3,
            label=season,
            color=SEASON_COLORS[season],
        )

    axes[1].axhline(
        0,
        color="black",
        linestyle="--",
    )

    axes[1].set_title(
        "Effect on Threshold Crossing Risk"
    )

    axes[1].set_xlabel(
        "Threshold (kW/m)"
    )

    axes[1].set_ylabel(
        "Probability difference (pp)"
    )

    axes[1].legend()

    axes[1].grid(alpha=0.25)

    fig.suptitle(
        "Paired Seasonal Interaction Effects",
        fontsize=16,
    )

    save_figure(
        fig,
        "seasonal_interaction_effects.png",
    )
    
def plot_validation(
    validation,
    interaction_validation,
):

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14, 6),
    )

    # Panel A: Expected log-intensity validation.
    x = np.arange(len(SEASONS))

    for i, model in enumerate(MODELS):

        data = (
            validation[
                validation["model"] == model
            ]
            .set_index("season")
            .loc[SEASONS]
        )

        offset = (
            i - 0.5
        ) * 0.18

        axes[0].errorbar(
            x + offset,
            data["log_mean_error"],
            yerr=(
                1.96
                * data["log_standard_error"]
            ),
            fmt="o",
            capsize=4,
            label=model,
            color=MODEL_COLORS[model],
        )

    axes[0].axhline(
        0,
        color="black",
        linestyle="--",
    )

    axes[0].set_xticks(x, SEASONS)

    axes[0].set_title(
        "Expected Log-Intensity Validation"
    )

    axes[0].set_ylabel(
        "Monte Carlo minus theoretical"
    )

    axes[0].legend()

    axes[0].grid(alpha=0.25)

    # Panel B: Analytical interaction validation.
    data = (
        interaction_validation
        .set_index("season")
        .loc[SEASONS]
    )

    axes[1].errorbar(
        x,
        data["log_effect_error"],
        yerr=(
            1.96
            * data["log_standard_error"]
        ),
        fmt="o",
        capsize=4,
        color="tab:green",
    )

    axes[1].axhline(
        0,
        color="black",
        linestyle="--",
    )

    axes[1].set_xticks(x, SEASONS)

    axes[1].set_title(
        "Expected Log-Interaction Validation"
    )

    axes[1].set_ylabel(
        "Monte Carlo minus theoretical"
    )

    axes[1].grid(alpha=0.25)

    fig.suptitle(
        "Analytical and Numerical Consistency",
        fontsize=16,
    )

    save_figure(
        fig,
        "seasonal_validation.png",
    )
    
def main():

    scenarios = pd.read_csv(
        RESULTS_DIR / "seasonal_scenarios.csv"
    )

    terminal = pd.read_csv(
        RESULTS_DIR / "seasonal_terminal.csv"
    )

    thresholds = pd.read_csv(
        RESULTS_DIR / "seasonal_thresholds.csv"
    )

    paired_terminal = pd.read_csv(
        RESULTS_DIR / "seasonal_paired_terminal.csv"
    )

    paired_thresholds = pd.read_csv(
        RESULTS_DIR / "seasonal_paired_thresholds.csv"
    )

    validation = pd.read_csv(
        RESULTS_DIR / "seasonal_validation.csv"
    )

    interaction_validation = pd.read_csv(
        RESULTS_DIR
        / "seasonal_interaction_validation.csv"
    )

    plot_environment(scenarios)

    plot_terminal_statistics(terminal)

    with np.load(
        RESULTS_DIR / "seasonal_terminal_samples.npz"
    ) as samples:

        plot_distributions(samples)

    plot_crossing_probabilities(thresholds)

    plot_first_passage(thresholds)

    plot_interaction_effects(
        paired_terminal,
        paired_thresholds,
    )

    plot_validation(
        validation,
        interaction_validation,
    )

    print(
        "All seven seasonal figures saved successfully."
    )

    plt.show()


if __name__ == "__main__":
    main()
    
# ============================================================
# Q5: SEASONAL COMPARISON - RESULTS AND INTERPRETATION
# ============================================================
#
# EXPERIMENT DESIGN
# ------------------------------------------------------------
# Four seasonal scenarios were constructed using historical
# mean wind speed and relative humidity.
#
# Each season was simulated under both the additive and
# interaction models, with 10,000 trajectories per model.
#
# Total configurations: 8
# Total simulated trajectories: 80,000
#
# All dynamic parameters were held constant to isolate
# the effects of seasonal initial environmental conditions.
#
#
# FIGURE 1: HISTORICAL ENVIRONMENTAL CONDITIONS
# ------------------------------------------------------------
# Winter has the highest initial wind and humidity anomalies.
#
# Spring combines above-average wind with below-average
# humidity, producing the highest initial intensity drift.
#
# Summer combines below-average wind with above-average
# humidity, producing the lowest initial intensity drift.
#
# These results demonstrate the importance of considering
# wind and humidity jointly rather than independently.
#
#
# FIGURE 2: TERMINAL INTENSITY STATISTICS
# ------------------------------------------------------------
# Both models produce the following mean intensity ranking:
#
# Spring > Fall > Winter > Summer
#
# Introducing the interaction term:
#   - Slightly increases empirical mean intensity.
#   - Increases standard deviation in every season.
#   - Decreases median intensity in every season.
#   - Increases the 97.5th percentile in every season.
#
# Therefore, the interaction affects distributional
# variability and upper quantiles more strongly than
# arithmetic mean intensity.
#
#
# FIGURE 3: TERMINAL INTENSITY DISTRIBUTIONS
# ------------------------------------------------------------
# All eight empirical distributions are right-skewed.
#
# The additive model has an analytical lognormal distribution.
# The interaction model is not generally exactly lognormal.
#
# Despite substantial overlap between the histograms,
# the interaction model exhibits higher empirical upper
# quantiles and standard deviations.
#
# The histograms do not establish asymptotic tail behavior.
#
#
# FIGURE 4: THRESHOLD CROSSING PROBABILITIES
# ------------------------------------------------------------
# Crossing probabilities decrease as thresholds increase,
# consistent with the mathematical nesting of crossing events.
#
# Spring has the highest estimated crossing probability
# at every investigated threshold under both models.
#
# Under the interaction model:
#
# K = 1000:
#     Spring: 50.84%
#     Summer: 40.96%
#
# K = 5000:
#     Spring: 1.01%
#     Summer: 0.53%
#
# Relative seasonal differences become larger at higher
# thresholds, although absolute differences become smaller.
#
#
# FIGURE 5: CONDITIONAL FIRST-PASSAGE TIMES
# ------------------------------------------------------------
# Higher thresholds have larger empirical conditional
# mean first-passage times in all configurations.
#
# Spring has the lowest estimated conditional mean
# first-passage time at every investigated threshold.
#
# These statistics are conditional on successful crossing
# before the simulation horizon.
#
# They do not represent unrestricted expected stopping times
# or establish that individual trajectories cross earlier.
#
#
# FIGURE 6: PAIRED INTERACTION EFFECTS
# ------------------------------------------------------------
# Mean terminal intensity changes:
#
# Winter: +0.67189 kW/m
# Spring: +16.34830 kW/m
# Summer: +7.88709 kW/m
# Fall:   +7.12458 kW/m
#
# Winter's paired confidence interval includes zero.
# The remaining three intervals exclude zero.
#
# Most importantly, the interaction effect on crossing
# probability changes direction across thresholds:
#
# At K = 1000:
#     All four seasons exhibit negative estimated effects.
#
# At K = 2000, 3000, 4000, 5000:
#     All four seasons exhibit positive estimated effects.
#
# This contrasts with Q3, where the interaction increased
# crossing probability at all five investigated thresholds.
#
# Therefore, interaction effects depend on both initial
# environmental conditions and the selected threshold.
#
#
# FIGURE 7: ANALYTICAL VALIDATION
# ------------------------------------------------------------
# All eight theoretical expected log-intensities fall
# within their corresponding approximate Monte Carlo
# confidence intervals.
#
# All four analytical expected log-interaction effects
# are also consistent with the paired numerical estimates.
#
# All absolute validation z-scores are below one.
#
# Common random numbers reduce uncertainty in paired
# comparisons but introduce dependence across scenarios.
#
# These results support numerical and analytical consistency.
# They do not establish real-world predictive accuracy.
#
#
# LIMITATIONS
# ------------------------------------------------------------
# Seasonal historical averages are used only as initial
# environmental conditions.
#
# All seasons share the same long-term OU equilibrium
# and dynamic model parameters.
#
# Temperature, seasonal fuel moisture, and fuel consumption
# are not explicitly modeled.
#
# Historical spread rate and simulated fire intensity are
# different physical quantities.
#
# Confidence intervals represent Monte Carlo uncertainty,
# not parameter uncertainty or real-world prediction error.
#
# Threshold crossings are evaluated on a discrete time grid.
#
# ============================================================