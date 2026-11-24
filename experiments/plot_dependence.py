import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import PercentFormatter


# ==========================================================
# CONFIGURATION
# ==========================================================

RESULTS_DIR = Path("results")

FIGURES_DIR = RESULTS_DIR / "figures"

FIGURES_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

CORRELATIONS = [
    -1.0,
    -0.75,
    -0.50,
    -0.25,
    0.0,
    0.25,
    0.50,
    0.75,
    1.0,
]

THRESHOLDS = [
    1000,
    2000,
    3000,
    4000,
    5000,
]

MODELS = [
    "Additive",
    "Interaction",
]

MODEL_COLORS = {
    "Additive": "tab:blue",
    "Interaction": "tab:orange",
}

MODEL_MARKERS = {
    "Additive": "o",
    "Interaction": "s",
}

SELECTED_CORRELATIONS = [
    -1.0,
    -0.5,
    0.0,
    0.5,
    1.0,
]

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 10,
    "figure.titlesize": 16,
    "axes.grid": True,
    "grid.alpha": 0.25,
})


# ==========================================================
# 1. GENERAL PLOTTING UTILITIES
# ==========================================================

def save_figure(fig, filename):
    """
    Save a figure and release its memory.
    """

    path = FIGURES_DIR / filename

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {path}")


def configure_correlation_axis(ax):
    """
    Apply consistent formatting to correlation axes.
    """

    ax.set_xlabel("Brownian correlation")

    ax.set_xticks(CORRELATIONS)

    ax.axvline(
        0,
        color="gray",
        linestyle=":",
        alpha=0.7,
    )

    ax.grid(
        True,
        alpha=0.25,
    )


def confidence_errors(mean, lower, upper):
    """
    Convert confidence interval endpoints into
    asymmetric error-bar lengths.
    """

    return np.vstack([
        np.maximum(0, mean - lower),
        np.maximum(0, upper - mean),
    ])


# ==========================================================
# 2. TERMINAL INTENSITY STATISTICS
# ==========================================================

def plot_terminal_statistics(terminal_df):
    """
    Compare terminal intensity statistics across
    correlation scenarios and intensity models.
    """

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(17, 10),
    )

    metrics = [
        ("mean", "Mean terminal intensity", "Intensity (kW/m)"),
        ("std", "Standard deviation", "Intensity (kW/m)"),
        ("median", "Median terminal intensity", "Intensity (kW/m)"),
        ("q975", "97.5th percentile", "Intensity (kW/m)"),
        ("log_mean", "Expected log-intensity", "Log-intensity"),
        ("log_variance", "Log-intensity variance", "Variance"),
    ]

    for ax, (metric, title, ylabel) in zip(
        axes.flat,
        metrics,
    ):

        for model in MODELS:

            subset = terminal_df[
                terminal_df["model"] == model
            ].sort_values("rho")

            ax.plot(
                subset["rho"],
                subset[metric],
                color=MODEL_COLORS[model],
                marker=MODEL_MARKERS[model],
                linewidth=2,
                label=model,
            )

        ax.set_title(title)

        ax.set_ylabel(ylabel)

        configure_correlation_axis(ax)

        ax.legend()

    fig.suptitle(
        "Terminal Wildfire Intensity Under Environmental Dependence"
    )

    save_figure(
        fig,
        "dependence_terminal_statistics.png",
    )


# ==========================================================
# 3. THEORETICAL INTENSITY VALIDATION
# ==========================================================

def plot_theoretical_validation(
    terminal_df,
    theoretical_df,
    validation_df,
):
    """
    Compare analytical and empirical intensity results.

    Both models have analytical expected log-intensity.

    Only the additive model has the exact lognormal
    distribution used for the remaining comparisons.
    """

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(17, 10),
    )

    # --------------------------------------------------
    # Panel A: Expected log-intensity
    # --------------------------------------------------

    ax = axes[0, 0]

    for model in MODELS:

        empirical = terminal_df[
            terminal_df["model"] == model
        ].sort_values("rho")

        theory = theoretical_df[
            theoretical_df["model"] == model
        ].sort_values("rho")

        color = MODEL_COLORS[model]

        ax.plot(
            theory["rho"],
            theory["log_mean"],
            color=color,
            linewidth=2,
            label=f"{model}: theory",
        )

        ax.errorbar(
            empirical["rho"],
            empirical["log_mean"],
            yerr=1.96 * empirical["log_standard_error"],
            color=color,
            marker=MODEL_MARKERS[model],
            linestyle="none",
            capsize=3,
            label=f"{model}: Monte Carlo",
        )

    ax.set_title("Expected log-intensity")

    ax.set_ylabel("Log-intensity")

    configure_correlation_axis(ax)

    ax.legend(fontsize=8)

    # --------------------------------------------------
    # Panels B-E: Additive model statistics
    # --------------------------------------------------

    additive_empirical = terminal_df[
        terminal_df["model"] == "Additive"
    ].sort_values("rho")

    additive_theory = theoretical_df[
        theoretical_df["model"] == "Additive"
    ].sort_values("rho")

    metrics = [
        ("log_variance", "Additive log variance", "Variance"),
        ("mean", "Additive arithmetic mean", "Intensity (kW/m)"),
        ("std", "Additive standard deviation", "Intensity (kW/m)"),
        ("q975", "Additive 97.5th percentile", "Intensity (kW/m)"),
    ]

    selected_axes = [
        axes[0, 1],
        axes[0, 2],
        axes[1, 0],
        axes[1, 1],
    ]

    for ax, (metric, title, ylabel) in zip(
        selected_axes,
        metrics,
    ):

        ax.plot(
            additive_theory["rho"],
            additive_theory[metric],
            color="tab:blue",
            linewidth=2,
            label="Theoretical",
        )

        ax.plot(
            additive_empirical["rho"],
            additive_empirical[metric],
            color="tab:orange",
            marker="s",
            linestyle="--",
            linewidth=1.5,
            label="Monte Carlo",
        )

        ax.set_title(title)

        ax.set_ylabel(ylabel)

        configure_correlation_axis(ax)

        ax.legend()

    # --------------------------------------------------
    # Panel F: Standardized log-mean errors
    # --------------------------------------------------

    ax = axes[1, 2]

    ax.axhspan(
        -1.96,
        1.96,
        color="gray",
        alpha=0.12,
        label="Approximate 95% reference range",
    )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    ax.axhline(
        1.96,
        color="gray",
        linestyle="--",
    )

    ax.axhline(
        -1.96,
        color="gray",
        linestyle="--",
    )

    for model in MODELS:

        subset = validation_df[
            validation_df["model"] == model
        ].sort_values("rho")

        ax.plot(
            subset["rho"],
            subset["log_mean_z"],
            color=MODEL_COLORS[model],
            marker=MODEL_MARKERS[model],
            label=model,
        )

    ax.set_title("Standardized log-mean errors")

    ax.set_ylabel("Z-score")

    configure_correlation_axis(ax)

    ax.legend(fontsize=8)

    fig.suptitle(
        "Theoretical and Numerical Validation Under Dependence"
    )

    save_figure(
        fig,
        "dependence_theoretical_validation.png",
    )


# ==========================================================
# 4. TERMINAL INTENSITY DISTRIBUTIONS
# ==========================================================

def plot_distributions(terminal_samples):
    """
    Compare terminal intensity distributions for
    selected correlation scenarios.

    Histograms use the full sample size when calculating
    density, even when the displayed range is truncated.
    """

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(17, 10),
    )

    # Determine a common plotting range.
    all_samples = np.concatenate([
        terminal_samples[f"{model}_{rho:.2f}"]
        for model in MODELS
        for rho in SELECTED_CORRELATIONS
    ])

    upper_limit = np.quantile(
        all_samples,
        0.99,
    )

    bins = np.linspace(
        0,
        upper_limit,
        65,
    )

    bin_widths = np.diff(bins)

    for ax, rho in zip(
        axes.flat,
        SELECTED_CORRELATIONS,
    ):

        for model in MODELS:

            samples = terminal_samples[
                f"{model}_{rho:.2f}"
            ]

            counts, edges = np.histogram(
                samples,
                bins=bins,
            )

            density = (
                counts
                / (len(samples) * bin_widths)
            )

            ax.stairs(
                density,
                edges,
                color=MODEL_COLORS[model],
                linewidth=1.8,
                label=model,
            )

            ax.stairs(
                density,
                edges,
                color=MODEL_COLORS[model],
                fill=True,
                alpha=0.12,
            )

        ax.set_title(
            f"Brownian correlation = {rho:.2f}"
        )

        ax.set_xlabel("Terminal intensity (kW/m)")

        ax.set_ylabel("Probability density")

        ax.set_xlim(0, upper_limit)

        ax.legend()

    # --------------------------------------------------
    # Final panel: Interaction model empirical CDFs
    # --------------------------------------------------

    ax = axes.flat[-1]

    cmap = plt.get_cmap("coolwarm")

    for rho in SELECTED_CORRELATIONS:

        samples = np.sort(
            terminal_samples[
                f"Interaction_{rho:.2f}"
            ]
        )

        probabilities = (
            np.arange(1, len(samples) + 1)
            / len(samples)
        )

        ax.plot(
            samples,
            probabilities,
            color=cmap((rho + 1) / 2),
            linewidth=1.8,
            label=f"rho = {rho:.2f}",
        )

    ax.set_title("Interaction model empirical CDFs")

    ax.set_xlabel("Terminal intensity (kW/m)")

    ax.set_ylabel("Cumulative probability")

    ax.set_xlim(0, upper_limit)

    ax.set_ylim(0, 1.02)

    ax.legend()

    fig.suptitle(
        "Terminal Intensity Distributions Under Dependence"
    )

    fig.text(
        0.5,
        0.01,
        "Displayed intensity range is truncated at the pooled "
        "99th percentile; densities use all simulated samples.",
        ha="center",
        fontsize=10,
    )

    fig.tight_layout(
        rect=[0, 0.035, 1, 0.96]
    )

    path = (
        FIGURES_DIR
        / "dependence_distributions.png"
    )

    fig.savefig(
        path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {path}")


# ==========================================================
# 5. THRESHOLD CROSSING PROBABILITIES
# ==========================================================

def plot_crossing_probabilities(threshold_df):
    """
    Compare pathwise threshold crossing probabilities
    across correlation scenarios.

    Separate panels allow the effect of dependence
    to be examined at each intensity threshold.
    """

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(17, 10),
    )

    for ax, threshold in zip(
        axes.flat,
        THRESHOLDS,
    ):

        for model in MODELS:

            subset = threshold_df[
                (threshold_df["model"] == model)
                & (threshold_df["threshold"] == threshold)
            ].sort_values("rho")

            ax.plot(
                subset["rho"],
                subset["crossing_probability"],
                color=MODEL_COLORS[model],
                marker=MODEL_MARKERS[model],
                linewidth=2,
                label=model,
            )

        ax.set_title(
            f"Threshold K = {threshold} kW/m"
        )

        ax.set_ylabel("Crossing probability")

        ax.yaxis.set_major_formatter(
            PercentFormatter(1.0)
        )

        ax.set_ylim(bottom=0)

        configure_correlation_axis(ax)

        ax.legend()

    # --------------------------------------------------
    # Final panel: High-threshold crossing probabilities
    # on a logarithmic scale.
    # --------------------------------------------------

    ax = axes.flat[-1]

    subset = threshold_df[
        threshold_df["threshold"] == 5000
    ]

    for model in MODELS:

        model_data = subset[
            subset["model"] == model
        ].sort_values("rho")

        probabilities = model_data[
            "crossing_probability"
        ].to_numpy()

        # Zero probabilities cannot be plotted
        # on a logarithmic axis.
        probabilities = np.where(
            probabilities > 0,
            probabilities,
            np.nan,
        )

        ax.plot(
            model_data["rho"],
            probabilities,
            color=MODEL_COLORS[model],
            marker=MODEL_MARKERS[model],
            linewidth=2,
            label=model,
        )

    ax.set_yscale("log")

    ax.set_title("K = 5000: logarithmic scale")

    ax.set_ylabel("Crossing probability")

    configure_correlation_axis(ax)

    ax.legend()

    fig.suptitle(
        "Threshold Crossing Risk Under Environmental Dependence"
    )

    save_figure(
        fig,
        "dependence_crossing_probability.png",
    )


# ==========================================================
# 6. CONDITIONAL FIRST PASSAGE TIMES
# ==========================================================

def plot_first_passage_times(threshold_df):
    """
    Plot conditional mean first passage times.

    Estimates based on fewer than 30 observed crossings
    are marked with hollow symbols.
    """

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(17, 10),
    )

    for ax, threshold in zip(
        axes.flat,
        THRESHOLDS,
    ):

        for model in MODELS:

            subset = threshold_df[
                (threshold_df["model"] == model)
                & (threshold_df["threshold"] == threshold)
            ].sort_values("rho")

            rho = subset["rho"].to_numpy()

            passage_times = subset[
                "mean_passage_time"
            ].to_numpy()

            n_crossed = subset[
                "n_crossed"
            ].to_numpy()

            sufficient = n_crossed >= 30

            insufficient = (
                (n_crossed > 0)
                & (n_crossed < 30)
            )

            color = MODEL_COLORS[model]

            marker = MODEL_MARKERS[model]

            ax.plot(
                rho,
                passage_times,
                color=color,
                linewidth=1.8,
                label=model,
            )

            ax.scatter(
                rho[sufficient],
                passage_times[sufficient],
                color=color,
                marker=marker,
                s=35,
                zorder=3,
            )

            ax.scatter(
                rho[insufficient],
                passage_times[insufficient],
                facecolors="none",
                edgecolors=color,
                marker=marker,
                s=65,
                linewidths=1.8,
                zorder=4,
            )

        ax.set_title(
            f"Threshold K = {threshold} kW/m"
        )

        ax.set_ylabel("Conditional mean passage time")

        configure_correlation_axis(ax)

        ax.legend()

    # --------------------------------------------------
    # Final panel: Number of observed crossings
    # for the highest threshold.
    # --------------------------------------------------

    ax = axes.flat[-1]

    for model in MODELS:

        subset = threshold_df[
            (threshold_df["model"] == model)
            & (threshold_df["threshold"] == 5000)
        ].sort_values("rho")

        ax.plot(
            subset["rho"],
            subset["n_crossed"],
            color=MODEL_COLORS[model],
            marker=MODEL_MARKERS[model],
            linewidth=2,
            label=model,
        )

    ax.axhline(
        30,
        color="gray",
        linestyle="--",
        label="30 crossings",
    )

    ax.set_title("Observed crossings at K = 5000")

    ax.set_ylabel("Number of crossings")

    configure_correlation_axis(ax)

    ax.legend()

    fig.suptitle(
        "Conditional First Passage Times Under Dependence"
    )

    fig.text(
        0.5,
        0.01,
        "Hollow markers indicate estimates based on fewer "
        "than 30 observed crossings.",
        ha="center",
        fontsize=10,
    )

    fig.tight_layout(
        rect=[0, 0.035, 1, 0.96]
    )

    path = (
        FIGURES_DIR
        / "dependence_first_passage.png"
    )

    fig.savefig(
        path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {path}")


# ==========================================================
# 7. PAIRED INTERACTION EFFECTS
# ==========================================================

def plot_interaction_effects(
    paired_terminal_df,
    paired_threshold_df,
):
    """
    Compare paired additive-interaction differences
    under different environmental correlations.
    """

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(14, 10),
    )

    data = paired_terminal_df.sort_values(
        "rho"
    )

    # --------------------------------------------------
    # Panel A: Arithmetic mean difference
    # --------------------------------------------------

    ax = axes[0, 0]

    errors = confidence_errors(
        data["mean_difference"],
        data["ci_lower"],
        data["ci_upper"],
    )

    ax.errorbar(
        data["rho"],
        data["mean_difference"],
        yerr=errors,
        color="tab:purple",
        marker="o",
        linewidth=2,
        capsize=4,
        label="Paired difference with 95% CI",
    )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    ax.set_title("Arithmetic mean interaction effect")

    ax.set_ylabel("Mean difference (kW/m)")

    configure_correlation_axis(ax)

    ax.legend()

    # --------------------------------------------------
    # Panel B: Percentage change
    # --------------------------------------------------

    ax = axes[0, 1]

    changes = data[
        "percentage_change"
    ].to_numpy()

    colors = [
        "tab:red" if change > 0 else "tab:blue"
        for change in changes
    ]

    ax.bar(
        data["rho"],
        changes,
        width=0.18,
        color=colors,
        alpha=0.8,
    )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    ax.set_title("Percentage change in mean intensity")

    ax.set_ylabel("Percentage change (%)")

    configure_correlation_axis(ax)

    # --------------------------------------------------
    # Panel C: Fraction of trajectories with increased
    # terminal intensity.
    # --------------------------------------------------

    ax = axes[1, 0]

    ax.plot(
        data["rho"],
        data["fraction_increased"],
        color="tab:green",
        marker="o",
        linewidth=2,
    )

    ax.axhline(
        0.5,
        color="gray",
        linestyle="--",
        label="50% of trajectories",
    )

    ax.set_ylim(-0.02, 1.02)

    ax.yaxis.set_major_formatter(
        PercentFormatter(1.0)
    )

    ax.set_title("Trajectories with increased intensity")

    ax.set_ylabel("Fraction of trajectories")

    configure_correlation_axis(ax)

    ax.legend()

    # --------------------------------------------------
    # Panel D: Paired crossing probability differences
    # --------------------------------------------------

    ax = axes[1, 1]

    for threshold in THRESHOLDS:

        subset = paired_threshold_df[
            paired_threshold_df["threshold"] == threshold
        ].sort_values("rho")

        ax.plot(
            subset["rho"],
            100 * subset["mean_difference"],
            marker="o",
            linewidth=1.5,
            label=f"K = {threshold}",
        )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    ax.set_title("Interaction effect on crossing risk")

    ax.set_ylabel("Probability difference (percentage points)")

    configure_correlation_axis(ax)

    ax.legend()

    fig.suptitle(
        "Paired Wind-Humidity Interaction Effects"
    )

    save_figure(
        fig,
        "dependence_interaction_effects.png",
    )


# ==========================================================
# 8. INTERACTION CROSSING HEATMAP
# ==========================================================

def plot_interaction_heatmap(paired_threshold_df):
    """
    Visualize paired crossing probability differences.

    Positive values indicate higher crossing probability
    under the interaction model.

    An asterisk indicates that the approximate pointwise
    95% paired confidence interval excludes zero.
    """

    matrix = paired_threshold_df.pivot(
        index="rho",
        columns="threshold",
        values="mean_difference",
    ).reindex(
        index=CORRELATIONS,
        columns=THRESHOLDS,
    )

    lower = paired_threshold_df.pivot(
        index="rho",
        columns="threshold",
        values="ci_lower",
    ).reindex(
        index=CORRELATIONS,
        columns=THRESHOLDS,
    )

    upper = paired_threshold_df.pivot(
        index="rho",
        columns="threshold",
        values="ci_upper",
    ).reindex(
        index=CORRELATIONS,
        columns=THRESHOLDS,
    )

    values = 100 * matrix.to_numpy()

    significant = (
        (lower.to_numpy() > 0)
        | (upper.to_numpy() < 0)
    )

    limit = max(
        np.max(np.abs(values)),
        1e-10,
    )

    norm = TwoSlopeNorm(
        vmin=-limit,
        vcenter=0,
        vmax=limit,
    )

    fig, ax = plt.subplots(
        figsize=(10, 8)
    )

    image = ax.imshow(
        values,
        cmap="RdBu_r",
        norm=norm,
        aspect="auto",
    )

    ax.set_xticks(
        np.arange(len(THRESHOLDS))
    )

    ax.set_xticklabels(THRESHOLDS)

    ax.set_yticks(
        np.arange(len(CORRELATIONS))
    )

    ax.set_yticklabels([
        f"{rho:.2f}"
        for rho in CORRELATIONS
    ])

    for i in range(len(CORRELATIONS)):

        for j in range(len(THRESHOLDS)):

            value = values[i, j]

            suffix = (
                "*"
                if significant[i, j]
                else ""
            )

            text_color = (
                "white"
                if abs(value) > 0.6 * limit
                else "black"
            )

            ax.text(
                j,
                i,
                f"{value:+.2f}{suffix}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=10,
            )

    colorbar = fig.colorbar(
        image,
        ax=ax,
    )

    colorbar.set_label(
        "Crossing probability difference (percentage points)"
    )

    ax.set_xlabel(
        "Intensity threshold (kW/m)"
    )

    ax.set_ylabel(
        "Brownian correlation"
    )

    ax.set_title(
        "Interaction Effect on Threshold Crossing Risk"
    )

    ax.grid(False)

    fig.text(
        0.5,
        0.01,
        "* Approximate pointwise 95% paired confidence "
        "interval excludes zero.",
        ha="center",
        fontsize=10,
    )

    fig.tight_layout(
        rect=[0, 0.04, 1, 1]
    )

    path = (
        FIGURES_DIR
        / "dependence_interaction_heatmap.png"
    )

    fig.savefig(
        path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {path}")


# ==========================================================
# 9. DEPENDENCE RELATIVE TO INDEPENDENCE
# ==========================================================

def plot_independence_comparison(
    dependence_terminal_df,
    dependence_threshold_df,
):
    """
    Compare correlated models against their corresponding
    independent environmental baseline.
    """

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(15, 10),
    )

    # --------------------------------------------------
    # Panel A: Mean intensity differences
    # --------------------------------------------------

    ax = axes[0, 0]

    for model in MODELS:

        subset = dependence_terminal_df[
            dependence_terminal_df["model"] == model
        ].sort_values("rho")

        errors = confidence_errors(
            subset["mean_difference"],
            subset["ci_lower"],
            subset["ci_upper"],
        )

        ax.errorbar(
            subset["rho"],
            subset["mean_difference"],
            yerr=errors,
            color=MODEL_COLORS[model],
            marker=MODEL_MARKERS[model],
            linewidth=2,
            capsize=3,
            label=model,
        )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    ax.set_title("Mean difference relative to independence")

    ax.set_ylabel("Mean difference (kW/m)")

    configure_correlation_axis(ax)

    ax.legend()

    # --------------------------------------------------
    # Panel B: Percentage change
    # --------------------------------------------------

    ax = axes[0, 1]

    for model in MODELS:

        subset = dependence_terminal_df[
            dependence_terminal_df["model"] == model
        ].sort_values("rho")

        ax.plot(
            subset["rho"],
            subset["percentage_change"],
            color=MODEL_COLORS[model],
            marker=MODEL_MARKERS[model],
            linewidth=2,
            label=model,
        )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    ax.set_title("Percentage change relative to independence")

    ax.set_ylabel("Percentage change (%)")

    configure_correlation_axis(ax)

    ax.legend()

    # --------------------------------------------------
    # Panels C-D: Threshold risk differences
    # --------------------------------------------------

    for ax, model in zip(
        axes[1],
        MODELS,
    ):

        for threshold in THRESHOLDS:

            subset = dependence_threshold_df[
                (dependence_threshold_df["model"] == model)
                & (
                    dependence_threshold_df["threshold"]
                    == threshold
                )
            ].sort_values("rho")

            ax.plot(
                subset["rho"],
                100 * subset["mean_difference"],
                marker="o",
                linewidth=1.5,
                label=f"K = {threshold}",
            )

        ax.axhline(
            0,
            color="black",
            linewidth=1,
        )

        ax.set_title(
            f"{model}: crossing risk relative to independence"
        )

        ax.set_ylabel(
            "Probability difference (percentage points)"
        )

        configure_correlation_axis(ax)

        ax.legend()

    fig.suptitle(
        "Consequences of the Environmental Independence Assumption"
    )

    save_figure(
        fig,
        "dependence_vs_independence.png",
    )


# ==========================================================
# 10. INTERACTION SIGN CHANGE
# ==========================================================

def plot_interaction_sign_change(
    paired_terminal_df,
):
    """
    Compare the analytical sign change in expected
    log-intensity with the empirical sign change in
    arithmetic mean intensity.
    """

    data = paired_terminal_df.sort_values(
        "rho"
    )

    rho = data["rho"].to_numpy()

    theoretical_effect = data[
        "theoretical_log_effect"
    ].to_numpy()

    empirical_effect = data[
        "log_mean_difference"
    ].to_numpy()

    # The theoretical expected log-interaction effect
    # is exactly affine in rho.
    slope, intercept = np.polyfit(
        rho,
        theoretical_effect,
        1,
    )

    critical_rho_log = -intercept / slope

    # Estimate the arithmetic mean sign change by
    # interpolating between adjacent scenarios.
    arithmetic_effect = data[
        "mean_difference"
    ].to_numpy()

    critical_rho_mean = np.nan

    for i in range(len(rho) - 1):

        left = arithmetic_effect[i]

        right = arithmetic_effect[i + 1]

        if left == 0:

            critical_rho_mean = rho[i]

            break

        if left * right < 0:

            critical_rho_mean = (
                rho[i]
                - left
                * (rho[i + 1] - rho[i])
                / (right - left)
            )

            break

    if arithmetic_effect[-1] == 0:

        critical_rho_mean = rho[-1]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15, 6),
    )

    # --------------------------------------------------
    # Panel A: Expected log-interaction effect
    # --------------------------------------------------

    ax = axes[0]

    fine_rho = np.linspace(
        -1,
        1,
        200,
    )

    ax.plot(
        fine_rho,
        intercept + slope * fine_rho,
        color="tab:blue",
        linewidth=2,
        label="Analytical expectation",
    )

    errors = confidence_errors(
        empirical_effect,
        data["log_ci_lower"],
        data["log_ci_upper"],
    )

    ax.errorbar(
        rho,
        empirical_effect,
        yerr=errors,
        color="tab:orange",
        marker="s",
        linestyle="none",
        capsize=4,
        label="Monte Carlo with 95% CI",
    )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    if -1 <= critical_rho_log <= 1:

        ax.axvline(
            critical_rho_log,
            color="tab:red",
            linestyle="--",
            label=(
                f"Analytical root = "
                f"{critical_rho_log:.4f}"
            ),
        )

    ax.set_title("Expected log-interaction effect")

    ax.set_ylabel("Expected log-intensity difference")

    configure_correlation_axis(ax)

    ax.legend()

    # --------------------------------------------------
    # Panel B: Arithmetic mean interaction effect
    # --------------------------------------------------

    ax = axes[1]

    errors = confidence_errors(
        arithmetic_effect,
        data["ci_lower"],
        data["ci_upper"],
    )

    ax.errorbar(
        rho,
        arithmetic_effect,
        yerr=errors,
        color="tab:purple",
        marker="o",
        linewidth=2,
        capsize=4,
        label="Paired mean difference with 95% CI",
    )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
    )

    if np.isfinite(critical_rho_mean):

        ax.axvline(
            critical_rho_mean,
            color="tab:red",
            linestyle="--",
            label=(
                f"Interpolated root = "
                f"{critical_rho_mean:.4f}"
            ),
        )

    ax.set_title("Arithmetic mean interaction effect")

    ax.set_ylabel("Mean intensity difference (kW/m)")

    configure_correlation_axis(ax)

    ax.legend()

    fig.suptitle(
        "Correlation-Dependent Sign Changes in Interaction Effects"
    )

    save_figure(
        fig,
        "dependence_interaction_sign_change.png",
    )

    print(
        f"Analytical log-effect root: "
        f"{critical_rho_log:.5f}"
    )

    print(
        f"Interpolated arithmetic mean root: "
        f"{critical_rho_mean:.5f}"
    )


# ==========================================================
# 11. MAIN PLOTTING WORKFLOW
# ==========================================================

def main():

    print("Q6 — DEPENDENCE ANALYSIS FIGURES")
    print("=" * 80)

    # --------------------------------------------------
    # Load previously generated numerical results
    # --------------------------------------------------

    terminal_df = pd.read_csv(
        RESULTS_DIR / "dependence_terminal.csv"
    )

    threshold_df = pd.read_csv(
        RESULTS_DIR / "dependence_thresholds.csv"
    )

    theoretical_df = pd.read_csv(
        RESULTS_DIR / "dependence_theoretical.csv"
    )

    validation_df = pd.read_csv(
        RESULTS_DIR / "dependence_validation.csv"
    )

    paired_terminal_df = pd.read_csv(
        RESULTS_DIR / "dependence_paired_terminal.csv"
    )

    paired_threshold_df = pd.read_csv(
        RESULTS_DIR / "dependence_paired_thresholds.csv"
    )

    dependence_terminal_df = pd.read_csv(
        RESULTS_DIR
        / "dependence_vs_independence_terminal.csv"
    )

    dependence_threshold_df = pd.read_csv(
        RESULTS_DIR
        / "dependence_vs_independence_thresholds.csv"
    )

    print("All numerical result tables loaded.")

    # --------------------------------------------------
    # Generate all dependence figures
    # --------------------------------------------------

    plot_terminal_statistics(
        terminal_df
    )

    plot_theoretical_validation(
        terminal_df,
        theoretical_df,
        validation_df,
    )

    with np.load(
        RESULTS_DIR
        / "dependence_terminal_samples.npz"
    ) as terminal_samples:

        plot_distributions(
            terminal_samples
        )

    plot_crossing_probabilities(
        threshold_df
    )

    plot_first_passage_times(
        threshold_df
    )

    plot_interaction_effects(
        paired_terminal_df,
        paired_threshold_df,
    )

    plot_interaction_heatmap(
        paired_threshold_df
    )

    plot_independence_comparison(
        dependence_terminal_df,
        dependence_threshold_df,
    )

    plot_interaction_sign_change(
        paired_terminal_df
    )

    print("\nAll 9 dependence analysis figures generated.")

    print("No additional simulations were required.")

    print("Q6 plotting completed.")


if __name__ == "__main__":
    main()