from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from experiments.sensitivity_analysis import (
    PARAMETER_LABELS
)


# ======================================================
# 1. Output configuration
# ======================================================

RESULTS_DIR = Path("results")

OUTPUT_DIR = RESULTS_DIR / "figures"

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ======================================================
# 2. Parameter response curves
# ======================================================

def plot_response_curves(
    results,
    metric,
    ylabel,
    title,
    filename,
    scale=1.0
):
    """
    Generate individual parameter response curves.

    Each subplot changes one parameter while
    holding all other parameters fixed.
    """

    parameters = list(PARAMETER_LABELS.keys())

    fig, axes = plt.subplots(
        3,
        4,
        figsize=(18, 12)
    )

    axes = axes.flatten()

    for ax, parameter in zip(axes, parameters):

        group = results[
            results["parameter"] == parameter
        ].sort_values("value")

        baseline = group[
            group["is_baseline"]
        ].iloc[0]

        # Plot parameter response.
        ax.plot(
            group["value"],
            group[metric] * scale,
            marker="o",
            linewidth=2,
            color="steelblue"
        )

        # Highlight baseline configuration.
        ax.scatter(
            baseline["value"],
            baseline[metric] * scale,
            color="red",
            s=65,
            zorder=5,
            label="Baseline"
        )

        # Baseline parameter value.
        ax.axvline(
            baseline["value"],
            color="gray",
            linestyle="--",
            alpha=0.6
        )

        # Baseline outcome value.
        ax.axhline(
            baseline[metric] * scale,
            color="gray",
            linestyle=":",
            alpha=0.6
        )

        ax.set_title(
            PARAMETER_LABELS[parameter]
        )

        ax.set_xlabel(parameter)

        ax.set_ylabel(ylabel)

        ax.grid(alpha=0.25)

    # Hide unused subplot.
    for ax in axes[len(parameters):]:

        ax.set_visible(False)

    fig.suptitle(
        title,
        fontsize=18
    )

    fig.tight_layout(
        rect=[0, 0, 1, 0.95]
    )

    fig.savefig(
        OUTPUT_DIR / filename,
        dpi=300,
        bbox_inches="tight"
    )


# ======================================================
# 3. Tornado sensitivity diagrams
# ======================================================

def plot_tornado(
    summary,
    lower_column,
    upper_column,
    xlabel,
    title,
    filename
):
    """
    Display minimum and maximum changes relative
    to the common baseline.

    Parameters are ordered by their observed
    outcome range.
    """

    data = summary.copy()

    data["span"] = (
        data[upper_column]
        - data[lower_column]
    )

    data = data.sort_values(
        "span",
        ascending=True
    )

    labels = data["parameter_label"]

    lower = data[lower_column].to_numpy()

    upper = data[upper_column].to_numpy()

    y = np.arange(len(data))

    fig, ax = plt.subplots(
        figsize=(11, 7)
    )

    # Negative changes relative to baseline.
    ax.barh(
        y,
        lower,
        color="steelblue",
        label="Decrease from baseline"
    )

    # Positive changes relative to baseline.
    ax.barh(
        y,
        upper,
        color="darkorange",
        label="Increase from baseline"
    )

    ax.axvline(
        0,
        color="black",
        linewidth=1
    )

    ax.set_yticks(y)

    ax.set_yticklabels(labels)

    ax.set_xlabel(xlabel)

    ax.set_title(title)

    ax.grid(
        axis="x",
        alpha=0.25
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / filename,
        dpi=300,
        bbox_inches="tight"
    )


# ======================================================
# 4. Multi-metric sensitivity heatmap
# ======================================================

def plot_sensitivity_heatmap(results):
    """
    Compare parameter sensitivity across multiple
    wildfire intensity and risk metrics.

    Each cell represents the maximum minus minimum
    outcome, divided by the baseline outcome.
    """

    metrics = {
        "mean": "Mean intensity",
        "std": "Standard deviation",
        "median": "Median intensity",
        "q975": "97.5th percentile",
        "crossing_probability": "Crossing probability"
    }

    parameters = list(PARAMETER_LABELS.keys())

    matrix = np.zeros(
        (len(parameters), len(metrics))
    )

    for i, parameter in enumerate(parameters):

        group = results[
            results["parameter"] == parameter
        ]

        baseline = group[
            group["is_baseline"]
        ].iloc[0]

        for j, metric in enumerate(metrics):

            minimum = group[metric].min()

            maximum = group[metric].max()

            baseline_value = baseline[metric]

            if baseline_value != 0:

                matrix[i, j] = (
                    100
                    * (maximum - minimum)
                    / baseline_value
                )

            else:

                matrix[i, j] = np.nan

    fig, ax = plt.subplots(
        figsize=(12, 8)
    )

    image = ax.imshow(
        matrix,
        cmap="YlOrRd",
        aspect="auto"
    )

    ax.set_xticks(
        np.arange(len(metrics))
    )

    ax.set_xticklabels(
        list(metrics.values()),
        rotation=25,
        ha="right"
    )

    ax.set_yticks(
        np.arange(len(parameters))
    )

    ax.set_yticklabels([
        PARAMETER_LABELS[p]
        for p in parameters
    ])

    # Annotate every heatmap cell.
    maximum_value = np.nanmax(matrix)

    for i in range(len(parameters)):

        for j in range(len(metrics)):

            value = matrix[i, j]

            text_color = (
                "white"
                if value > maximum_value * 0.6
                else "black"
            )

            ax.text(
                j,
                i,
                f"{value:.1f}%",
                ha="center",
                va="center",
                color=text_color,
                fontsize=10
            )

    colorbar = fig.colorbar(
        image,
        ax=ax
    )

    colorbar.set_label(
        "Outcome range relative to baseline (%)"
    )

    ax.set_title(
        "Multi-Metric Parameter Sensitivity"
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "sensitivity_heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )


# ======================================================
# 5. Main visualization workflow
# ======================================================

def main():

    results = pd.read_csv(
        RESULTS_DIR / "sensitivity_results.csv"
    )

    summary = pd.read_csv(
        RESULTS_DIR / "sensitivity_summary.csv"
    )

    # --------------------------------------------------
    # Figure 1: Mean intensity response
    # --------------------------------------------------

    plot_response_curves(
        results=results,
        metric="mean",
        ylabel="Mean intensity (kW/m)",
        title="Sensitivity of Mean Terminal Intensity",
        filename="sensitivity_mean_intensity.png"
    )

    # --------------------------------------------------
    # Figure 2: Threshold crossing response
    # --------------------------------------------------

    plot_response_curves(
        results=results,
        metric="crossing_probability",
        ylabel="Crossing probability (%)",
        title="Sensitivity of Threshold Crossing Probability",
        filename="sensitivity_crossing_probability.png",
        scale=100
    )

    # --------------------------------------------------
    # Figure 3: Mean intensity tornado diagram
    # --------------------------------------------------

    plot_tornado(
        summary=summary,
        lower_column="mean_change_min_pct",
        upper_column="mean_change_max_pct",
        xlabel="Change in mean intensity (%)",
        title="Tornado Diagram: Mean Intensity Sensitivity",
        filename="sensitivity_tornado_mean.png"
    )

    # --------------------------------------------------
    # Figure 4: Crossing probability tornado diagram
    # --------------------------------------------------

    plot_tornado(
        summary=summary,
        lower_column="crossing_change_min_pp",
        upper_column="crossing_change_max_pp",
        xlabel="Change in crossing probability (pp)",
        title="Tornado Diagram: Threshold Crossing Sensitivity",
        filename="sensitivity_tornado_crossing.png"
    )

    # --------------------------------------------------
    # Figure 5: Multi-metric sensitivity heatmap
    # --------------------------------------------------

    plot_sensitivity_heatmap(results)

    print("All sensitivity figures saved successfully.")

    print(f"Output directory: {OUTPUT_DIR}")

    plt.show()


if __name__ == "__main__":
    main()
    

# ======================================================
# Q4: SENSITIVITY ANALYSIS - RESULTS AND INTERPRETATION
# ======================================================
#
# Experimental configuration:
#
# - 11 model inputs investigated.
# - 5 settings per input.
# - 55 recorded scenarios.
# - 45 unique simulated configurations.
# - 10,000 trajectories per configuration.
# - Common random numbers across configurations.
# - Threshold K = 2000 kW/m.
# - Simulation horizon T = 10.
#
# Baseline results:
#
# Mean terminal intensity = 1596.96718 kW/m
# Standard deviation      = 1438.32228 kW/m
# Upper 97.5th percentile = 5434.03406 kW/m
# Crossing probability    = 0.31970
#
# ------------------------------------------------------
# FIGURE 1: MEAN TERMINAL INTENSITY
# ------------------------------------------------------
#
# Mean intensity increases with:
# - Baseline growth coefficient.
# - Wind coefficient.
# - Humidity coefficient.
# - Wind volatility.
# - Humidity volatility.
# - Initial wind anomaly.
#
# Mean intensity decreases with:
# - Interaction coefficient increasing toward zero.
# - Wind mean-reversion rate.
# - Humidity mean-reversion rate.
# - Initial humidity anomaly.
#
# These directions apply to the investigated scenarios,
# not necessarily every possible model configuration.
#
# Intensity volatility produces almost no variation
# in mean terminal intensity.
#
# Minimum mean = 1596.79883
# Maximum mean = 1598.31008
# Relative range = 0.09463%
#
# This supports the theoretical result that independent
# multiplicative intensity noise does not change the
# expected terminal intensity, when the mean is finite.
#
# ------------------------------------------------------
# FIGURE 2: THRESHOLD CROSSING PROBABILITY
# ------------------------------------------------------
#
# Baseline crossing probability = 31.97%.
#
# Increasing intensity volatility from 0.075 to 0.225
# increases crossing probability from 27.82% to 35.23%.
#
# Absolute change = 7.41 percentage points.
#
# This demonstrates that mean intensity and extreme
# event risk are distinct stochastic quantities.
#
# Initial wind and humidity conditions also produce
# substantial changes in threshold crossing probability.
#
# Initial wind:
# - Minimum probability = 10.76%.
# - Maximum probability = 31.97%.
#
# Initial humidity:
# - Minimum probability = 11.99%.
# - Maximum probability = 31.97%.
#
# ------------------------------------------------------
# FIGURE 3: MEAN INTENSITY TORNADO DIAGRAM
# ------------------------------------------------------
#
# Largest observed mean intensity sensitivity ranges:
#
# 1. Wind mean reversion:     73.49%
# 2. Humidity mean reversion: 50.55%
# 3. Baseline growth:         50.52%
#
# Wind mean reversion produces asymmetric sensitivity.
#
# Reducing kappa_w from 0.8 to 0.4 increases mean
# intensity by approximately 60.15%.
#
# Increasing kappa_w from 0.8 to 1.2 decreases mean
# intensity by approximately 13.33%.
#
# Slower mean reversion increases persistence of
# initial anomalies and stationary OU variance when
# environmental volatility is held constant.
#
# Initial condition bars are one-sided because the
# baseline uses the highest tested initial wind and
# lowest tested initial humidity.
#
# ------------------------------------------------------
# FIGURE 4: CROSSING PROBABILITY TORNADO DIAGRAM
# ------------------------------------------------------
#
# Largest observed crossing probability ranges:
#
# 1. Initial wind anomaly:  21.21 percentage points.
# 2. Wind mean reversion:   21.17 percentage points.
# 3. Baseline growth:       21.08 percentage points.
#
# These three sensitivity ranges are very similar.
# Their ordering should not be interpreted as a
# statistically established ranking.
#
# Wind mean reversion produces the largest increase
# above baseline crossing probability:
#
# 46.40% - 31.97% = 14.43 percentage points.
#
# ------------------------------------------------------
# FIGURE 5: MULTI-METRIC SENSITIVITY HEATMAP
# ------------------------------------------------------
#
# Baseline growth produces identical relative ranges
# for mean, standard deviation, median and quantiles.
#
# Theoretical explanation:
#
# I_T(alpha + delta) = exp(delta * T) * I_T(alpha)
#
# Therefore, all four terminal distribution metrics
# scale by the same positive multiplicative factor.
#
# The theoretical relative range is:
#
# 100 * (exp(0.25) - exp(-0.25)) = 50.52246%.
#
# Wind mean reversion produces the largest observed
# sensitivity for:
#
# - Standard deviation: 175.6%.
# - Upper 97.5th percentile: 138.2%.
#
# Environmental volatility produces relatively small
# changes in median intensity but substantial changes
# in standard deviation and upper quantiles.
#
# Intensity volatility produces:
#
# - Mean sensitivity: 0.1%.
# - Standard deviation sensitivity: 46.5%.
# - Median sensitivity: 21.4%.
# - Upper quantile sensitivity: 31.7%.
# - Relative crossing probability sensitivity: 23.2%.
#
# Heatmap values are relative outcome ranges, whereas
# the crossing tornado diagram uses percentage points.
#
# ------------------------------------------------------
# MATHEMATICAL VALIDATION
# ------------------------------------------------------
#
# Analytical versus numerical log-intensity derivatives:
#
# Parameter     Theory      Numerical
# alpha_c       10.00000    10.00000
# beta_w         1.87437     1.83324
# beta_h         1.98652     2.07059
# gamma         -1.15384    -1.16399
# sigma_i       -1.50000    -1.48006
#
# The humidity coefficient derivative differs from
# the continuous-time theoretical value by about 4.23%.
#
# Its theoretical value falls slightly outside the
# approximate numerical 95% confidence interval.
#
# Possible sources of discrepancy include Monte Carlo
# sampling variability and time discretization.
#
# All mathematical consistency checks passed.
# All 12 existing automated tests passed.
#
# ------------------------------------------------------
# LIMITATIONS
# ------------------------------------------------------
#
# - Model parameters are illustrative, not calibrated.
# - Sensitivity ranges depend on the chosen inputs.
# - One-at-a-time analysis does not measure joint
#   parameter uncertainty or parameter interactions.
# - Initial condition ranges are not symmetric
#   around the baseline configuration.
# - Rankings are descriptive, not statistically
#   established parameter importance rankings.
# - Threshold sensitivity depends on K and T.
# - Discrete observation can miss crossings between
#   consecutive numerical time points.
# - Sensitivity ranges do not include uncertainty
#   intervals for the estimated parameter rankings.
#
# ------------------------------------------------------
# CONCLUSION
# ------------------------------------------------------
#
# Environmental persistence, initial environmental
# conditions and baseline growth have substantial
# effects on simulated wildfire intensity and risk.
#
# Environmental mean reversion strongly influences
# dispersion and upper-tail intensity.
#
# Independent intensity volatility has little effect
# on expected intensity but substantially affects
# uncertainty and threshold crossing probability.
#
# Parameter importance depends on the outcome metric
# and the investigated parameter range.
#
# Q4 numerical sensitivity analysis completed.
# ======================================================