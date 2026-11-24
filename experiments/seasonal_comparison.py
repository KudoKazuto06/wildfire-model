"""
Q5: Seasonal environmental comparison.

Investigate how historical seasonal initial
environmental conditions influence:

1. Terminal intensity distributions.
2. Expected log-intensity.
3. Threshold crossing probabilities.
4. Conditional first-passage times.
5. Wind-humidity interaction effects.

All dynamic model parameters remain fixed.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import skew

from src.config import (
    GAMMA,
    N_PATHS,
    SEED,
    T,
    DT,
)

from src.seasonal import (
    build_seasonal_scenarios,
    theoretical_seasonal_comparison,
)

from src.monte_carlo import simulate_many

from src.stopping_times import threshold_statistics


# Experiment configuration.
THRESHOLDS = [1000, 2000, 3000, 4000, 5000]

RESULTS_DIR = Path("results")

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

def terminal_statistics(samples):
    """
    Calculate empirical terminal intensity statistics.
    """

    samples = np.asarray(samples)

    log_samples = np.log(samples)

    n_paths = len(samples)

    mean = np.mean(samples)

    std = np.std(samples, ddof=1)

    standard_error = std / np.sqrt(n_paths)

    log_mean = np.mean(log_samples)

    log_std = np.std(log_samples, ddof=1)

    log_standard_error = (
        log_std / np.sqrt(n_paths)
    )

    return {
        "n_paths": n_paths,

        "mean": mean,

        "variance": np.var(
            samples,
            ddof=1,
        ),

        "std": std,

        "median": np.median(samples),

        "q025": np.quantile(samples, 0.025),

        "q975": np.quantile(samples, 0.975),

        "skewness": skew(
            samples,
            bias=False,
        ),

        "log_mean": log_mean,

        "log_variance": np.var(
            log_samples,
            ddof=1,
        ),

        "geometric_mean": np.exp(log_mean),

        "mean_standard_error": standard_error,

        "mean_ci_lower": (
            mean - 1.96 * standard_error
        ),

        "mean_ci_upper": (
            mean + 1.96 * standard_error
        ),

        "log_standard_error": log_standard_error,

        "log_ci_lower": (
            log_mean - 1.96 * log_standard_error
        ),

        "log_ci_upper": (
            log_mean + 1.96 * log_standard_error
        ),
    }
    
def simulate_season(
    season,
    wind_z0,
    humidity_z0,
    model,
    gamma,
):
    """
    Simulate one seasonal model configuration.

    Returns:
        terminal intensity statistics,
        threshold crossing statistics,
        terminal intensity samples,
        maximum intensity samples.
    """

    times, paths = simulate_many(
        n_paths=N_PATHS,
        wind_z0=wind_z0,
        humidity_z0=humidity_z0,
        gamma=gamma,
        horizon=T,
        dt=DT,
        seed=SEED,
    )

    # Terminal intensity of every trajectory.
    terminal = paths[:, -1].copy()

    # Maximum intensity reached by each trajectory.
    maximum = np.max(paths, axis=1)

    # Calculate terminal distribution statistics.
    terminal_result = terminal_statistics(
        terminal
    )

    terminal_result.update(
        {
            "season": season,
            "model": model,
            "gamma": gamma,
            "wind_z0": wind_z0,
            "humidity_z0": humidity_z0,
        }
    )

    threshold_results = []

    # Reuse the stopping-time implementation
    # developed in Q3.
    for threshold in THRESHOLDS:

        statistics = threshold_statistics(
            times,
            paths,
            threshold,
        )

        statistics.update(
            {
                "season": season,
                "model": model,
                "gamma": gamma,
            }
        )

        threshold_results.append(statistics)

    # Full trajectories are no longer needed.
    del paths

    return (
        terminal_result,
        threshold_results,
        terminal,
        maximum,
    )
    
def paired_interaction_analysis(
    season,
    additive_terminal,
    interaction_terminal,
    additive_maximum,
    interaction_maximum,
):
    """
    Compare additive and interaction models
    using common random numbers.
    """

    n_paths = len(additive_terminal)

    # Paired arithmetic intensity differences.
    differences = (
        interaction_terminal
        - additive_terminal
    )

    mean_difference = np.mean(differences)

    standard_error = (
        np.std(differences, ddof=1)
        / np.sqrt(n_paths)
    )

    # Paired log-intensity differences.
    log_differences = (
        np.log(interaction_terminal)
        - np.log(additive_terminal)
    )

    log_difference = np.mean(log_differences)

    log_standard_error = (
        np.std(log_differences, ddof=1)
        / np.sqrt(n_paths)
    )

    terminal_result = {
        "season": season,

        "mean_difference": mean_difference,

        "standard_error": standard_error,

        "ci_lower": (
            mean_difference
            - 1.96 * standard_error
        ),

        "ci_upper": (
            mean_difference
            + 1.96 * standard_error
        ),

        "mean_change_pct": (
            100
            * mean_difference
            / np.mean(additive_terminal)
        ),

        "fraction_increased": np.mean(
            differences > 0
        ),

        "log_difference": log_difference,

        "log_standard_error": log_standard_error,

        "log_ci_lower": (
            log_difference
            - 1.96 * log_standard_error
        ),

        "log_ci_upper": (
            log_difference
            + 1.96 * log_standard_error
        ),
    }

    threshold_results = []

    for threshold in THRESHOLDS:

        additive_crossed = (
            additive_maximum >= threshold
        )

        interaction_crossed = (
            interaction_maximum >= threshold
        )

        # Values:
        #  1: interaction model only crossed
        # -1: additive model only crossed
        #  0: both or neither crossed
        paired_difference = (
            interaction_crossed.astype(float)
            - additive_crossed.astype(float)
        )

        probability_difference = np.mean(
            paired_difference
        )

        probability_se = (
            np.std(paired_difference, ddof=1)
            / np.sqrt(n_paths)
        )

        threshold_results.append(
            {
                "season": season,

                "threshold": threshold,

                "additive_probability": np.mean(
                    additive_crossed
                ),

                "interaction_probability": np.mean(
                    interaction_crossed
                ),

                "probability_difference": (
                    probability_difference
                ),

                "standard_error": probability_se,

                "ci_lower": (
                    probability_difference
                    - 1.96 * probability_se
                ),

                "ci_upper": (
                    probability_difference
                    + 1.96 * probability_se
                ),

                "interaction_only": np.mean(
                    paired_difference == 1
                ),

                "additive_only": np.mean(
                    paired_difference == -1
                ),
            }
        )

    return terminal_result, threshold_results

def validate_results(
    terminal_df,
    threshold_df,
    paired_terminal_df,
    paired_threshold_df,
):
    """
    Validate mathematical consistency
    of seasonal simulation results.
    """

    # Four seasons and two models.
    assert len(terminal_df) == 8

    # Five thresholds per configuration.
    assert len(threshold_df) == 40

    assert len(paired_terminal_df) == 4

    assert len(paired_threshold_df) == 20

    # All terminal intensity means must be positive.
    assert (terminal_df["mean"] > 0).all()

    # Geometric mean cannot exceed arithmetic mean.
    assert (
        terminal_df["geometric_mean"]
        <= terminal_df["mean"] + 1e-10
    ).all()

    # Quantile ordering.
    assert (
        terminal_df["q025"]
        <= terminal_df["median"]
    ).all()

    assert (
        terminal_df["median"]
        <= terminal_df["q975"]
    ).all()

    # Crossing probabilities must be valid.
    assert threshold_df[
        "crossing_probability"
    ].between(0, 1).all()

    # Terminal exceedance cannot exceed
    # pathwise crossing probability.
    assert (
        threshold_df["terminal_probability"]
        <= threshold_df["crossing_probability"]
        + 1e-12
    ).all()

    # Crossing probability must decrease
    # or remain unchanged as threshold increases.
    for _, group in threshold_df.groupby(
        ["season", "model"]
    ):

        probabilities = group.sort_values(
            "threshold"
        )["crossing_probability"].to_numpy()

        assert np.all(
            np.diff(probabilities) <= 1e-12
        )

    # Paired probability differences must equal
    # differences between marginal probabilities.
    np.testing.assert_allclose(
        paired_threshold_df[
            "probability_difference"
        ],
        (
            paired_threshold_df[
                "interaction_probability"
            ]
            - paired_threshold_df[
                "additive_probability"
            ]
        ),
        atol=1e-12,
    )

    print(
        "\nAll seasonal mathematical "
        "consistency checks passed."
    )
    
def main():

    scenarios = build_seasonal_scenarios()

    theoretical = theoretical_seasonal_comparison()

    terminal_results = []

    threshold_results = []

    paired_terminal_results = []

    paired_threshold_results = []

    terminal_samples = {}

    maximum_samples = {}

    print("SEASONAL MONTE CARLO EXPERIMENT")
    print("=" * 90)

    for _, row in scenarios.iterrows():

        season = row["season"]

        wind_z0 = row["wind_z0"]

        humidity_z0 = row["humidity_z0"]

        print(f"\nInvestigating {season}...")

        # Additive and interaction configurations.
        for model, gamma in [
            ("Additive", 0.0),
            ("Interaction", GAMMA),
        ]:

            (
                terminal_result,
                threshold_result,
                terminal,
                maximum,
            ) = simulate_season(
                season=season,
                wind_z0=wind_z0,
                humidity_z0=humidity_z0,
                model=model,
                gamma=gamma,
            )

            terminal_results.append(
                terminal_result
            )

            threshold_results.extend(
                threshold_result
            )

            terminal_samples[
                f"{season}_{model}"
            ] = terminal

            maximum_samples[
                f"{season}_{model}"
            ] = maximum

            print(
                f"  {model:12s}"
                f" | Mean = "
                f"{terminal_result['mean']:10.2f}"
                f" | Median = "
                f"{terminal_result['median']:10.2f}"
                f" | Q97.5 = "
                f"{terminal_result['q975']:10.2f}"
            )

        # Compare both models within this season.
        (
            paired_terminal,
            paired_threshold,
        ) = paired_interaction_analysis(
            season=season,

            additive_terminal=terminal_samples[
                f"{season}_Additive"
            ],

            interaction_terminal=terminal_samples[
                f"{season}_Interaction"
            ],

            additive_maximum=maximum_samples[
                f"{season}_Additive"
            ],

            interaction_maximum=maximum_samples[
                f"{season}_Interaction"
            ],
        )

        paired_terminal_results.append(
            paired_terminal
        )

        paired_threshold_results.extend(
            paired_threshold
        )

    # Convert numerical results into DataFrames.
    terminal_df = pd.DataFrame(
        terminal_results
    )

    threshold_df = pd.DataFrame(
        threshold_results
    )

    paired_terminal_df = pd.DataFrame(
        paired_terminal_results
    )

    paired_threshold_df = pd.DataFrame(
        paired_threshold_results
    )

    validate_results(
        terminal_df,
        threshold_df,
        paired_terminal_df,
        paired_threshold_df,
    )

    # Compare analytical and numerical results.
    validation = terminal_df.merge(
        theoretical,
        on="season",
        suffixes=("", "_theory"),
    )

    validation["theoretical_log_mean"] = np.where(
        validation["model"] == "Additive",
        validation["additive_log_mean"],
        validation["interaction_log_mean"],
    )

    validation["log_mean_error"] = (
        validation["log_mean"]
        - validation["theoretical_log_mean"]
    )

    validation["log_mean_z_score"] = (
        validation["log_mean_error"]
        / validation["log_standard_error"]
    )

    validation["theory_inside_ci"] = (
        (
            validation["theoretical_log_mean"]
            >= validation["log_ci_lower"]
        )
        &
        (
            validation["theoretical_log_mean"]
            <= validation["log_ci_upper"]
        )
    )

    # Validate paired interaction effects.
    paired_validation = paired_terminal_df.merge(
        theoretical[
            [
                "season",
                "log_interaction_effect",
            ]
        ],
        on="season",
    )

    paired_validation["log_effect_error"] = (
        paired_validation["log_difference"]
        - paired_validation["log_interaction_effect"]
    )

    paired_validation["log_effect_z_score"] = (
        paired_validation["log_effect_error"]
        / paired_validation["log_standard_error"]
    )

    # Print terminal intensity results.
    print("\nTERMINAL INTENSITY COMPARISON")
    print("=" * 110)

    print(
        terminal_df[
            [
                "season",
                "model",
                "mean",
                "std",
                "median",
                "q975",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # Print analytical validation.
    print("\nANALYTICAL LOG-INTENSITY VALIDATION")
    print("=" * 110)

    print(
        validation[
            [
                "season",
                "model",
                "theoretical_log_mean",
                "log_mean",
                "log_mean_error",
                "log_mean_z_score",
                "theory_inside_ci",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # Print threshold crossing comparison.
    print("\nTHRESHOLD CROSSING COMPARISON")
    print("=" * 110)

    crossing_table = threshold_df.pivot(
        index=["season", "model"],
        columns="threshold",
        values="crossing_probability",
    )

    print(
        crossing_table.to_string(
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # Print conditional first-passage times.
    print("\nCONDITIONAL MEAN FIRST-PASSAGE TIMES")
    print("=" * 110)

    passage_table = threshold_df.pivot(
        index=["season", "model"],
        columns="threshold",
        values="mean_passage_time",
    )

    print(
        passage_table.to_string(
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # Print paired interaction results.
    print("\nPAIRED SEASONAL INTERACTION EFFECT")
    print("=" * 110)

    print(
        paired_terminal_df[
            [
                "season",
                "mean_difference",
                "mean_change_pct",
                "ci_lower",
                "ci_upper",
                "fraction_increased",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    print("\nANALYTICAL INTERACTION VALIDATION")
    print("=" * 110)

    print(
        paired_validation[
            [
                "season",
                "log_interaction_effect",
                "log_difference",
                "log_effect_error",
                "log_effect_z_score",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    print("\nPAIRED CROSSING PROBABILITY DIFFERENCES")
    print("=" * 110)

    interaction_table = paired_threshold_df.pivot(
        index="season",
        columns="threshold",
        values="probability_difference",
    )

    print(
        interaction_table.to_string(
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # Save numerical results.
    outputs = {
        "seasonal_scenarios.csv": scenarios,

        "seasonal_theoretical.csv": theoretical,

        "seasonal_terminal.csv": terminal_df,

        "seasonal_thresholds.csv": threshold_df,

        "seasonal_paired_terminal.csv":
            paired_terminal_df,

        "seasonal_paired_thresholds.csv":
            paired_threshold_df,

        "seasonal_validation.csv": validation,

        "seasonal_interaction_validation.csv":
            paired_validation,
    }

    for filename, dataframe in outputs.items():

        dataframe.to_csv(
            RESULTS_DIR / filename,
            index=False,
        )

    # Save raw samples for later visualization.
    # This avoids repeating the Monte Carlo experiment.
    np.savez_compressed(
        RESULTS_DIR / "seasonal_terminal_samples.npz",
        **terminal_samples,
    )

    np.savez_compressed(
        RESULTS_DIR / "seasonal_maximum_samples.npz",
        **maximum_samples,
    )

    print("\nEXPERIMENT SUMMARY")
    print("=" * 90)

    print("Seasons:", len(scenarios))

    print("Configurations:", len(terminal_df))

    print("Trajectories per configuration:", N_PATHS)

    print(
        "Total simulated trajectories:",
        len(terminal_df) * N_PATHS,
    )

    print(
        "Threshold evaluations:",
        len(threshold_df),
    )

    print(
        "\nAll seasonal results saved successfully."
    )


if __name__ == "__main__":
    main()