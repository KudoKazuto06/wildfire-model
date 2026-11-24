from pathlib import Path

import numpy as np
import pandas as pd

from src.config import N_PATHS, SEED

from src.sensitivity import (
    BASELINE_PARAMETERS,
    evaluate_parameter
)


# ======================================================
# 1. Experiment configuration
# ======================================================

THRESHOLD = 2000

RELATIVE_FACTORS = [
    0.50,
    0.75,
    1.00,
    1.25,
    1.50
]

PARAMETER_LABELS = {
    "alpha_c": "Baseline growth",
    "beta_w": "Wind coefficient",
    "beta_h": "Humidity coefficient",
    "gamma": "Interaction coefficient",
    "sigma_i": "Intensity volatility",
    "kappa_w": "Wind mean reversion",
    "kappa_h": "Humidity mean reversion",
    "sigma_w": "Wind volatility",
    "sigma_h": "Humidity volatility",
    "wind_z0": "Initial wind anomaly",
    "humidity_z0": "Initial humidity anomaly"
}


# ======================================================
# 2. Construct parameter experiment ranges
# ======================================================

def build_parameter_ranges():
    """
    Construct one-at-a-time sensitivity scenarios.

    Positive dynamic parameters use relative
    perturbations around their baseline values.

    Interaction and initial environmental conditions
    use explicitly defined scenario ranges.
    """

    ranges = {}

    for parameter, baseline in BASELINE_PARAMETERS.items():

        if parameter == "gamma":

            ranges[parameter] = [
                -0.08,
                -0.06,
                -0.04,
                -0.02,
                0.0
            ]

        elif parameter == "wind_z0":

            ranges[parameter] = [
                -1.5,
                -0.75,
                0.0,
                0.75,
                1.5
            ]

        elif parameter == "humidity_z0":

            ranges[parameter] = [
                -1.0,
                -0.5,
                0.0,
                0.5,
                1.0
            ]

        else:

            ranges[parameter] = [
                baseline * factor
                for factor in RELATIVE_FACTORS
            ]

    return ranges


# ======================================================
# 3. Calculate relative outcome changes
# ======================================================

def relative_change(value, baseline):
    """
    Calculate percentage change relative to baseline.
    """

    if baseline == 0:
        return np.nan

    return 100 * (value - baseline) / baseline


def add_relative_changes(results, baseline):
    """
    Compare simulation outcomes against the common
    baseline configuration.
    """

    results = results.copy()

    metrics = [
        "mean",
        "std",
        "median",
        "q975",
        "crossing_probability"
    ]

    for metric in metrics:

        results[f"{metric}_change_pct"] = (
            results[metric].apply(
                lambda value: relative_change(
                    value,
                    baseline[metric]
                )
            )
        )

    # Absolute change in crossing probability,
    # expressed in percentage points.
    results["crossing_change_pp"] = (
        100 * (
            results["crossing_probability"]
            - baseline["crossing_probability"]
        )
    )

    # Absolute change in expected log-intensity.
    results["log_mean_change"] = (
        results["log_mean"]
        - baseline["log_mean"]
    )

    return results


# ======================================================
# 4. Run one-at-a-time sensitivity experiments
# ======================================================

def run_sensitivity_experiments():

    parameter_ranges = build_parameter_ranges()

    results = []

    # --------------------------------------------------
    # Simulate the baseline configuration once
    # --------------------------------------------------

    print("Simulating common baseline configuration...")

    baseline = evaluate_parameter(
        parameter="alpha_c",
        value=BASELINE_PARAMETERS["alpha_c"],
        n_paths=N_PATHS,
        threshold=THRESHOLD,
        seed=SEED
    )

    print(
        f"Baseline mean intensity: "
        f"{baseline['mean']:.5f}"
    )

    print(
        f"Baseline crossing probability: "
        f"{baseline['crossing_probability']:.5f}"
    )

    # --------------------------------------------------
    # Iterate through all parameters and scenarios
    # --------------------------------------------------

    for parameter, values in parameter_ranges.items():

        print(
            f"\nInvestigating "
            f"{PARAMETER_LABELS[parameter]}..."
        )

        baseline_value = BASELINE_PARAMETERS[parameter]

        for value in values:

            is_baseline = np.isclose(
                value,
                baseline_value,
                rtol=0,
                atol=1e-12
            )

            # Reuse baseline results rather than
            # simulating the same configuration again.
            if is_baseline:

                result = baseline.copy()

                result["parameter"] = parameter

                result["value"] = value

            else:

                result = evaluate_parameter(
                    parameter=parameter,
                    value=value,
                    n_paths=N_PATHS,
                    threshold=THRESHOLD,
                    seed=SEED
                )

            result["baseline_value"] = baseline_value

            result["is_baseline"] = is_baseline

            result["parameter_label"] = (
                PARAMETER_LABELS[parameter]
            )

            results.append(result)

            print(
                f"  Value = {value:9.4f} | "
                f"Mean = {result['mean']:10.2f} | "
                f"P(Crossing) = "
                f"{result['crossing_probability']:.4f}"
            )

    results_df = pd.DataFrame(results)

    results_df = add_relative_changes(
        results_df,
        baseline
    )

    return results_df, baseline


# ======================================================
# 5. Summarize sensitivity across parameter ranges
# ======================================================

def summarize_sensitivity(results):

    summaries = []

    for parameter, group in results.groupby(
        "parameter",
        sort=False
    ):

        min_mean = group["mean"].min()

        max_mean = group["mean"].max()

        min_probability = (
            group["crossing_probability"].min()
        )

        max_probability = (
            group["crossing_probability"].max()
        )

        min_mean_change = (
            group["mean_change_pct"].min()
        )

        max_mean_change = (
            group["mean_change_pct"].max()
        )

        min_crossing_change = (
            group["crossing_change_pp"].min()
        )

        max_crossing_change = (
            group["crossing_change_pp"].max()
        )

        summaries.append({
            "parameter": parameter,

            "parameter_label":
                PARAMETER_LABELS[parameter],

            "min_mean": min_mean,

            "max_mean": max_mean,

            "mean_range": max_mean - min_mean,

            "mean_change_min_pct": min_mean_change,

            "mean_change_max_pct": max_mean_change,

            "mean_change_span_pct":
                max_mean_change - min_mean_change,

            "min_crossing_probability":
                min_probability,

            "max_crossing_probability":
                max_probability,

            "crossing_range_pp":
                100 * (
                    max_probability
                    - min_probability
                ),

            "crossing_change_min_pp":
                min_crossing_change,

            "crossing_change_max_pp":
                max_crossing_change
        })

    summary_df = pd.DataFrame(summaries)

    summary_df = summary_df.sort_values(
        by="crossing_range_pp",
        ascending=False
    )

    return summary_df


# ======================================================
# 6. Validate mathematical consistency
# ======================================================

def validate_results(results):

    expected_rows = (
        len(BASELINE_PARAMETERS) * 5
    )

    assert len(results) == expected_rows

    assert results["is_baseline"].sum() == (
        len(BASELINE_PARAMETERS)
    )

    # Intensity must remain positive.
    assert (results["mean"] > 0).all()

    assert (results["median"] > 0).all()

    assert (results["std"] >= 0).all()

    # Probabilities must lie between zero and one.
    assert results["crossing_probability"].between(
        0, 1
    ).all()

    assert results["terminal_probability"].between(
        0, 1
    ).all()

    # Terminal exceedance cannot exceed
    # first-passage probability.
    assert (
        results["terminal_probability"]
        <= results["crossing_probability"]
    ).all()

    # Quantile ordering.
    assert (
        results["q025"] <= results["median"]
    ).all()

    assert (
        results["median"] <= results["q975"]
    ).all()

    # Every baseline configuration should have
    # zero change relative to the common baseline.
    baseline_rows = results[
        results["is_baseline"]
    ]

    assert np.allclose(
        baseline_rows["mean_change_pct"],
        0
    )

    assert np.allclose(
        baseline_rows["crossing_change_pp"],
        0
    )

    # Increasing alpha_c must increase intensity
    # pathwise when common random numbers are used.
    alpha_results = results[
        results["parameter"] == "alpha_c"
    ].sort_values("value")

    assert (
        np.diff(alpha_results["mean"]) >= -1e-10
    ).all()

    assert (
        np.diff(
            alpha_results["crossing_probability"]
        ) >= -1e-10
    ).all()

    print("\nAll sensitivity consistency checks passed.")


# ======================================================
# 7. Main experiment
# ======================================================

def main():

    output_dir = Path("results")

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    results, baseline = (
        run_sensitivity_experiments()
    )

    validate_results(results)

    summary = summarize_sensitivity(results)

    # Save complete numerical results.
    results.to_csv(
        output_dir / "sensitivity_results.csv",
        index=False
    )

    # Save parameter-level summary.
    summary.to_csv(
        output_dir / "sensitivity_summary.csv",
        index=False
    )

    # --------------------------------------------------
    # Print sensitivity summary
    # --------------------------------------------------

    print("\nSENSITIVITY SUMMARY")
    print("=" * 120)

    display_columns = [
        "parameter",
        "min_mean",
        "max_mean",
        "mean_change_span_pct",
        "min_crossing_probability",
        "max_crossing_probability",
        "crossing_range_pp"
    ]

    print(
        summary[display_columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}"
        )
    )

    print("\nTotal scenarios:", len(results))

    print("Trajectories per scenario:", N_PATHS)

    print(
        "Total unique simulated configurations:",
        1 + int((~results["is_baseline"]).sum())
    )

    print("\nAll sensitivity results saved successfully.")


if __name__ == "__main__":
    main()