import numpy as np
import pandas as pd

from pathlib import Path

from scipy.stats import skew, norm

from src.config import (
    KAPPA_W,
    SIGMA_W,
    KAPPA_H,
    SIGMA_H,
    ALPHA_C,
    BETA_W,
    BETA_H,
    GAMMA,
    SIGMA_I,
    I0,
    T,
    DT,
    SEED,
    N_PATHS,
)

from src.monte_carlo import simulate_many

from src.stopping_times import threshold_statistics

from src.theoretical import integrated_ou_moments

from src.dependence import integrated_ou_covariance


# ==========================================================
# EXPERIMENT CONFIGURATION
# ==========================================================

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

MODELS = {
    "Additive": 0.0,
    "Interaction": GAMMA,
}

WIND_Z0 = 1.5

HUMIDITY_Z0 = -1.0

RESULTS_DIR = Path("results")


# ==========================================================
# 1. THEORETICAL INTENSITY STATISTICS
# ==========================================================

def theoretical_statistics(rho, gamma):
    """
    Calculate theoretical terminal intensity statistics
    under correlated environmental Brownian motions.

    The additive model has an exact lognormal terminal
    distribution.

    For the interaction model, calculate the exact
    expected log-intensity and geometric mean.
    """

    # Integrated marginal OU moments.
    wind_mean, wind_variance = integrated_ou_moments(
        WIND_Z0,
        KAPPA_W,
        SIGMA_W,
        T,
    )

    humidity_mean, humidity_variance = (
        integrated_ou_moments(
            HUMIDITY_Z0,
            KAPPA_H,
            SIGMA_H,
            T,
        )
    )

    # Covariance between integrated OU processes.
    integrated_covariance = integrated_ou_covariance(
        rho=rho,
        horizon=T,
    )

    # Expected additive log-intensity.
    additive_log_mean = (
        np.log(I0)
        + (ALPHA_C - 0.5 * SIGMA_I**2) * T
        + BETA_W * wind_mean
        - BETA_H * humidity_mean
    )

    # Exact additive log-intensity variance.
    additive_log_variance = (
        BETA_W**2 * wind_variance
        + BETA_H**2 * humidity_variance
        - 2 * BETA_W * BETA_H
        * integrated_covariance
        + SIGMA_I**2 * T
    )

    # --------------------------------------------------
    # Expected interaction contribution
    # --------------------------------------------------

    combined_kappa = KAPPA_W + KAPPA_H

    decay_integral = (
        -np.expm1(-combined_kappa * T)
        / combined_kappa
    )

    # Contribution from deterministic initial states.
    initial_effect = (
        gamma
        * WIND_Z0
        * HUMIDITY_Z0
        * decay_integral
    )

    # Contribution from environmental covariance.
    covariance_effect = (
        gamma
        * rho
        * SIGMA_W
        * SIGMA_H
        / combined_kappa
        * (T - decay_integral)
    )

    interaction_effect = (
        initial_effect + covariance_effect
    )

    log_mean = (
        additive_log_mean
        + interaction_effect
    )

    result = {
        "rho": rho,
        "gamma": gamma,
        "integrated_covariance": integrated_covariance,
        "log_mean": log_mean,
        "geometric_mean": np.exp(log_mean),
        "initial_effect": initial_effect,
        "covariance_effect": covariance_effect,
        "interaction_effect": interaction_effect,
    }

    # --------------------------------------------------
    # Exact lognormal moments for the additive model
    # --------------------------------------------------

    if gamma == 0.0:

        log_variance = additive_log_variance

        mean = np.exp(
            log_mean + 0.5 * log_variance
        )

        variance = (
            np.expm1(log_variance)
            * np.exp(
                2 * log_mean + log_variance
            )
        )

        result.update({
            "log_variance": log_variance,

            "mean": mean,

            "variance": variance,

            "std": np.sqrt(variance),

            "median": np.exp(log_mean),

            "q025": np.exp(
                log_mean
                + np.sqrt(log_variance)
                * norm.ppf(0.025)
            ),

            "q975": np.exp(
                log_mean
                + np.sqrt(log_variance)
                * norm.ppf(0.975)
            ),

            "skewness": (
                np.exp(log_variance) + 2
            ) * np.sqrt(
                np.expm1(log_variance)
            ),
        })

    else:

        # The interaction model is generally not
        # lognormal. Do not apply additive lognormal
        # formulas to its terminal intensity.

        result.update({
            "log_variance": np.nan,
            "mean": np.nan,
            "variance": np.nan,
            "std": np.nan,
            "median": np.nan,
            "q025": np.nan,
            "q975": np.nan,
            "skewness": np.nan,
        })

    return result


# ==========================================================
# 2. EMPIRICAL TERMINAL INTENSITY STATISTICS
# ==========================================================

def terminal_statistics(samples):
    """
    Calculate empirical terminal intensity statistics
    from Monte Carlo trajectories.
    """

    n = len(samples)

    log_samples = np.log(samples)

    mean = np.mean(samples)

    variance = np.var(
        samples,
        ddof=1,
    )

    std = np.sqrt(variance)

    mean_se = std / np.sqrt(n)

    log_mean = np.mean(log_samples)

    log_variance = np.var(
        log_samples,
        ddof=1,
    )

    log_se = np.sqrt(
        log_variance / n
    )

    return {
        "n_paths": n,

        "mean": mean,

        "variance": variance,

        "std": std,

        "median": np.median(samples),

        "q025": np.quantile(samples, 0.025),

        "q975": np.quantile(samples, 0.975),

        "skewness": skew(
            samples,
            bias=False,
        ),

        "log_mean": log_mean,

        "log_variance": log_variance,

        "log_skewness": skew(
            log_samples,
            bias=False,
        ),

        "geometric_mean": np.exp(log_mean),

        "mean_standard_error": mean_se,

        "mean_ci_lower": mean - 1.96 * mean_se,

        "mean_ci_upper": mean + 1.96 * mean_se,

        "log_standard_error": log_se,

        "log_ci_lower": log_mean - 1.96 * log_se,

        "log_ci_upper": log_mean + 1.96 * log_se,
    }


# ==========================================================
# 3. SIMULATE ONE MODEL CONFIGURATION
# ==========================================================

def simulate_configuration(rho, model, gamma):
    """
    Simulate one correlation scenario and one
    wildfire intensity model.
    """

    times, paths = simulate_many(
        n_paths=N_PATHS,
        wind_z0=WIND_Z0,
        humidity_z0=HUMIDITY_Z0,
        gamma=gamma,
        horizon=T,
        dt=DT,
        seed=SEED,
        rho=rho,
    )

    # Extract terminal intensity and maximum intensity.
    terminal = paths[:, -1].copy()

    maximum = np.max(
        paths,
        axis=1,
    )

    # Calculate terminal statistics.
    terminal_result = {
        "rho": rho,
        "model": model,
        "gamma": gamma,
    }

    terminal_result.update(
        terminal_statistics(terminal)
    )

    # Calculate threshold crossing statistics.
    threshold_results = []

    for threshold in THRESHOLDS:

        result = threshold_statistics(
            times,
            paths,
            threshold,
        )

        result.update({
            "rho": rho,
            "model": model,
            "gamma": gamma,
        })

        threshold_results.append(result)

    # Release the complete trajectory matrix.
    # Only terminal and maximum samples are retained.
    del paths

    return (
        terminal_result,
        threshold_results,
        terminal,
        maximum,
    )


# ==========================================================
# 4. PAIRED DIFFERENCE STATISTICS
# ==========================================================

def paired_statistics(differences):
    """
    Calculate the mean, standard error, and approximate
    95% confidence interval of paired differences.
    """

    mean = np.mean(differences)

    standard_error = (
        np.std(differences, ddof=1)
        / np.sqrt(len(differences))
    )

    return {
        "mean_difference": mean,

        "standard_error": standard_error,

        "ci_lower": (
            mean - 1.96 * standard_error
        ),

        "ci_upper": (
            mean + 1.96 * standard_error
        ),

        "fraction_increased": np.mean(
            differences > 0
        ),
    }


# ==========================================================
# 5. PAIRED INTERACTION ANALYSIS
# ==========================================================

def paired_interaction_analysis(
    rho,
    additive_terminal,
    interaction_terminal,
    additive_maximum,
    interaction_maximum,
):
    """
    Compare the interaction model against the additive
    model under the same environmental correlation.
    """

    # --------------------------------------------------
    # Arithmetic terminal intensity differences
    # --------------------------------------------------

    differences = (
        interaction_terminal
        - additive_terminal
    )

    terminal_result = {
        "rho": rho,
    }

    terminal_result.update(
        paired_statistics(differences)
    )

    terminal_result["percentage_change"] = (
        100
        * np.mean(differences)
        / np.mean(additive_terminal)
    )

    # --------------------------------------------------
    # Log-intensity differences
    # --------------------------------------------------

    log_differences = (
        np.log(interaction_terminal)
        - np.log(additive_terminal)
    )

    log_result = paired_statistics(
        log_differences
    )

    for key, value in log_result.items():

        terminal_result[f"log_{key}"] = value

    # Exact expected log-interaction effect.
    theory = theoretical_statistics(
        rho,
        GAMMA,
    )

    theoretical_effect = theory[
        "interaction_effect"
    ]

    terminal_result["theoretical_log_effect"] = (
        theoretical_effect
    )

    terminal_result["log_effect_error"] = (
        np.mean(log_differences)
        - theoretical_effect
    )

    terminal_result["log_theory_inside_ci"] = (
        log_result["ci_lower"]
        <= theoretical_effect
        <= log_result["ci_upper"]
    )

    # --------------------------------------------------
    # Paired threshold crossing differences
    # --------------------------------------------------

    threshold_results = []

    for threshold in THRESHOLDS:

        additive_crossed = (
            additive_maximum >= threshold
        )

        interaction_crossed = (
            interaction_maximum >= threshold
        )

        # Convert Boolean indicators to integers
        # before subtracting.
        differences = (
            interaction_crossed.astype(float)
            - additive_crossed.astype(float)
        )

        result = {
            "rho": rho,
            "threshold": threshold,

            "additive_probability": np.mean(
                additive_crossed
            ),

            "interaction_probability": np.mean(
                interaction_crossed
            ),

            "interaction_only": np.mean(
                interaction_crossed
                & ~additive_crossed
            ),

            "additive_only": np.mean(
                additive_crossed
                & ~interaction_crossed
            ),
        }

        result.update(
            paired_statistics(differences)
        )

        threshold_results.append(result)

    return terminal_result, threshold_results


# ==========================================================
# 6. EFFECT OF CORRELATION RELATIVE TO INDEPENDENCE
# ==========================================================

def compare_against_independence(
    terminal_samples,
    maximum_samples,
):
    """
    Compare every correlation scenario against rho = 0
    using common random numbers.

    Perform the comparison separately for each model.
    """

    terminal_results = []

    threshold_results = []

    for model in MODELS:

        baseline_key = (model, 0.0)

        baseline_terminal = terminal_samples[
            baseline_key
        ]

        baseline_maximum = maximum_samples[
            baseline_key
        ]

        for rho in CORRELATIONS:

            key = (model, rho)

            terminal = terminal_samples[key]

            maximum = maximum_samples[key]

            # Paired terminal intensity differences.
            differences = (
                terminal - baseline_terminal
            )

            result = {
                "model": model,
                "rho": rho,
                "baseline_mean": np.mean(
                    baseline_terminal
                ),
                "scenario_mean": np.mean(
                    terminal
                ),
            }

            result.update(
                paired_statistics(differences)
            )

            result["percentage_change"] = (
                100
                * np.mean(differences)
                / np.mean(baseline_terminal)
            )

            terminal_results.append(result)

            # Paired threshold crossing differences.
            for threshold in THRESHOLDS:

                baseline_crossed = (
                    baseline_maximum >= threshold
                )

                scenario_crossed = (
                    maximum >= threshold
                )

                differences = (
                    scenario_crossed.astype(float)
                    - baseline_crossed.astype(float)
                )

                result = {
                    "model": model,
                    "rho": rho,
                    "threshold": threshold,

                    "baseline_probability": np.mean(
                        baseline_crossed
                    ),

                    "scenario_probability": np.mean(
                        scenario_crossed
                    ),
                }

                result.update(
                    paired_statistics(differences)
                )

                threshold_results.append(result)

    return terminal_results, threshold_results


# ==========================================================
# 7. THEORETICAL VALIDATION
# ==========================================================

def build_theoretical_validation(
    terminal_df,
    theoretical_df,
):
    """
    Compare theoretical and empirical log-intensity
    means for both models.

    Also compare log-intensity variances and arithmetic
    means for the additive model.
    """

    rows = []

    for _, empirical in terminal_df.iterrows():

        rho = empirical["rho"]

        model = empirical["model"]

        theory = theoretical_df.loc[
            (theoretical_df["rho"] == rho)
            & (theoretical_df["model"] == model)
        ].iloc[0]

        log_error = (
            empirical["log_mean"]
            - theory["log_mean"]
        )

        log_z = (
            log_error
            / empirical["log_standard_error"]
        )

        result = {
            "rho": rho,

            "model": model,

            "theoretical_log_mean":
                theory["log_mean"],

            "empirical_log_mean":
                empirical["log_mean"],

            "log_mean_error": log_error,

            "log_mean_z": log_z,

            "log_theory_inside_ci": (
                empirical["log_ci_lower"]
                <= theory["log_mean"]
                <= empirical["log_ci_upper"]
            ),

            "theoretical_log_variance":
                theory["log_variance"],

            "empirical_log_variance":
                empirical["log_variance"],

            "theoretical_mean":
                theory["mean"],

            "empirical_mean":
                empirical["mean"],
        }

        if model == "Additive":

            result["log_variance_error"] = (
                empirical["log_variance"]
                - theory["log_variance"]
            )

            result["mean_error"] = (
                empirical["mean"]
                - theory["mean"]
            )

            result["mean_z"] = (
                result["mean_error"]
                / (
                    theory["std"]
                    / np.sqrt(N_PATHS)
                )
            )

        else:

            result["log_variance_error"] = np.nan

            result["mean_error"] = np.nan

            result["mean_z"] = np.nan

        rows.append(result)

    return pd.DataFrame(rows)


# ==========================================================
# 8. MATHEMATICAL CONSISTENCY CHECKS
# ==========================================================

def validate_results(
    terminal_df,
    threshold_df,
    theoretical_df,
    paired_terminal_df,
    paired_threshold_df,
    dependence_terminal_df,
    dependence_threshold_df,
):
    """
    Verify mathematical and structural consistency
    of the complete dependence experiment.
    """

    # Expected number of model configurations.
    assert len(terminal_df) == 18

    assert len(theoretical_df) == 18

    assert len(threshold_df) == 90

    assert len(paired_terminal_df) == 9

    assert len(paired_threshold_df) == 45

    assert len(dependence_terminal_df) == 18

    assert len(dependence_threshold_df) == 90

    # Positive terminal intensity statistics.
    assert np.all(
        terminal_df["mean"] > 0
    )

    assert np.all(
        terminal_df["median"] > 0
    )

    assert np.all(
        terminal_df["variance"] >= 0
    )

    # Quantile ordering.
    assert np.all(
        terminal_df["q025"]
        <= terminal_df["median"]
    )

    assert np.all(
        terminal_df["median"]
        <= terminal_df["q975"]
    )

    # Arithmetic mean >= geometric mean.
    assert np.all(
        terminal_df["mean"]
        >= terminal_df["geometric_mean"] - 1e-10
    )

    # Probability bounds.
    assert np.all(
        threshold_df["crossing_probability"]
        .between(0, 1)
    )

    assert np.all(
        threshold_df["terminal_probability"]
        .between(0, 1)
    )

    # Terminal exceedance implies pathwise crossing.
    assert np.all(
        threshold_df["terminal_probability"]
        <= threshold_df["crossing_probability"]
    )

    # Crossing probability must decrease or remain
    # constant as the threshold increases.
    for _, group in threshold_df.groupby(
        ["rho", "model"]
    ):

        group = group.sort_values(
            "threshold"
        )

        probabilities = group[
            "crossing_probability"
        ].to_numpy()

        assert np.all(
            np.diff(probabilities) <= 1e-12
        )

    # Paired interaction probability identity.
    expected_difference = (
        paired_threshold_df[
            "interaction_probability"
        ]
        - paired_threshold_df[
            "additive_probability"
        ]
    )

    assert np.allclose(
        expected_difference,
        paired_threshold_df["mean_difference"],
    )

    # Paired correlation probability identity.
    expected_difference = (
        dependence_threshold_df[
            "scenario_probability"
        ]
        - dependence_threshold_df[
            "baseline_probability"
        ]
    )

    assert np.allclose(
        expected_difference,
        dependence_threshold_df[
            "mean_difference"
        ],
    )

    # Independence compared against itself must
    # produce exactly zero differences.
    independent_rows = dependence_terminal_df[
        dependence_terminal_df["rho"] == 0.0
    ]

    assert np.all(
        independent_rows["mean_difference"] == 0
    )

    assert np.all(
        independent_rows["standard_error"] == 0
    )

    # Recover original Q2 baseline.
    baseline = terminal_df.loc[
        (terminal_df["rho"] == 0.0)
        & (terminal_df["model"] == "Interaction")
    ].iloc[0]

    assert np.isclose(
        baseline["mean"],
        1596.96718,
        rtol=0,
        atol=1e-5,
    )

    assert np.isclose(
        baseline["std"],
        1438.32228,
        rtol=0,
        atol=1e-5,
    )

    print("\nAll mathematical consistency checks passed.")


# ==========================================================
# 9. MAIN EXPERIMENT
# ==========================================================

def main():

    RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("Q6 — ENVIRONMENTAL DEPENDENCE ANALYSIS")
    print("=" * 100)

    print(f"Correlation scenarios: {len(CORRELATIONS)}")

    print(f"Models: {len(MODELS)}")

    print(f"Paths per configuration: {N_PATHS}")

    print(
        f"Total trajectories: "
        f"{len(CORRELATIONS) * len(MODELS) * N_PATHS}"
    )

    terminal_rows = []

    threshold_rows = []

    theoretical_rows = []

    paired_terminal_rows = []

    paired_threshold_rows = []

    terminal_samples = {}

    maximum_samples = {}

    # --------------------------------------------------
    # Run all correlation scenarios
    # --------------------------------------------------

    for rho in CORRELATIONS:

        print(f"\nSimulating rho = {rho:.2f}")

        for model, gamma in MODELS.items():

            print(f"  {model} model ...")

            (
                terminal_result,
                threshold_results,
                terminal,
                maximum,
            ) = simulate_configuration(
                rho,
                model,
                gamma,
            )

            terminal_rows.append(
                terminal_result
            )

            threshold_rows.extend(
                threshold_results
            )

            theory = theoretical_statistics(
                rho,
                gamma,
            )

            theory["model"] = model

            theoretical_rows.append(theory)

            terminal_samples[
                (model, rho)
            ] = terminal

            maximum_samples[
                (model, rho)
            ] = maximum

        # Paired additive-interaction comparison.
        paired_terminal, paired_thresholds = (
            paired_interaction_analysis(
                rho,

                terminal_samples[
                    ("Additive", rho)
                ],

                terminal_samples[
                    ("Interaction", rho)
                ],

                maximum_samples[
                    ("Additive", rho)
                ],

                maximum_samples[
                    ("Interaction", rho)
                ],
            )
        )

        paired_terminal_rows.append(
            paired_terminal
        )

        paired_threshold_rows.extend(
            paired_thresholds
        )

    # --------------------------------------------------
    # Compare correlation scenarios against rho = 0
    # --------------------------------------------------

    (
        dependence_terminal_rows,
        dependence_threshold_rows,
    ) = compare_against_independence(
        terminal_samples,
        maximum_samples,
    )

    # --------------------------------------------------
    # Convert results to DataFrames
    # --------------------------------------------------

    terminal_df = pd.DataFrame(
        terminal_rows
    )

    threshold_df = pd.DataFrame(
        threshold_rows
    )

    theoretical_df = pd.DataFrame(
        theoretical_rows
    )

    paired_terminal_df = pd.DataFrame(
        paired_terminal_rows
    )

    paired_threshold_df = pd.DataFrame(
        paired_threshold_rows
    )

    dependence_terminal_df = pd.DataFrame(
        dependence_terminal_rows
    )

    dependence_threshold_df = pd.DataFrame(
        dependence_threshold_rows
    )

    validation_df = build_theoretical_validation(
        terminal_df,
        theoretical_df,
    )

    # --------------------------------------------------
    # Validate numerical results
    # --------------------------------------------------

    validate_results(
        terminal_df,
        threshold_df,
        theoretical_df,
        paired_terminal_df,
        paired_threshold_df,
        dependence_terminal_df,
        dependence_threshold_df,
    )

    # --------------------------------------------------
    # Print terminal intensity statistics
    # --------------------------------------------------

    print("\nTERMINAL INTENSITY STATISTICS")
    print("=" * 110)

    columns = [
        "rho",
        "model",
        "mean",
        "std",
        "median",
        "q975",
        "log_mean",
        "log_variance",
    ]

    print(
        terminal_df[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Print theoretical validation
    # --------------------------------------------------

    print("\nTHEORETICAL LOG-INTENSITY VALIDATION")
    print("=" * 110)

    columns = [
        "rho",
        "model",
        "theoretical_log_mean",
        "empirical_log_mean",
        "log_mean_error",
        "log_mean_z",
        "log_theory_inside_ci",
    ]

    print(
        validation_df[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Print threshold crossing probabilities
    # --------------------------------------------------

    print("\nTHRESHOLD CROSSING PROBABILITIES")
    print("=" * 110)

    crossing_table = threshold_df.pivot(
        index=["rho", "model"],
        columns="threshold",
        values="crossing_probability",
    )

    print(
        crossing_table.to_string(
            float_format=lambda x: f"{x:.5f}"
        )
    )

    # --------------------------------------------------
    # Print conditional first passage times
    # --------------------------------------------------

    print("\nCONDITIONAL MEAN FIRST PASSAGE TIMES")
    print("=" * 110)

    passage_table = threshold_df.pivot(
        index=["rho", "model"],
        columns="threshold",
        values="mean_passage_time",
    )

    print(
        passage_table.to_string(
            float_format=lambda x: f"{x:.5f}"
        )
    )

    # --------------------------------------------------
    # Print paired interaction effects
    # --------------------------------------------------

    print("\nPAIRED INTERACTION EFFECTS")
    print("=" * 110)

    columns = [
        "rho",
        "mean_difference",
        "percentage_change",
        "ci_lower",
        "ci_upper",
        "fraction_increased",
    ]

    print(
        paired_terminal_df[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Print analytical interaction validation
    # --------------------------------------------------

    print("\nEXPECTED LOG-INTERACTION EFFECT")
    print("=" * 110)

    columns = [
        "rho",
        "theoretical_log_effect",
        "log_mean_difference",
        "log_effect_error",
        "log_ci_lower",
        "log_ci_upper",
        "log_theory_inside_ci",
    ]

    print(
        paired_terminal_df[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Print paired interaction crossing differences
    # --------------------------------------------------

    print("\nPAIRED INTERACTION CROSSING EFFECTS")
    print("=" * 110)

    interaction_table = (
        paired_threshold_df.pivot(
            index="rho",
            columns="threshold",
            values="mean_difference",
        )
    )

    print(
        interaction_table.to_string(
            float_format=lambda x: f"{x:.5f}"
        )
    )

    # --------------------------------------------------
    # Print correlation effects relative to independence
    # --------------------------------------------------

    print("\nDEPENDENCE EFFECT RELATIVE TO RHO = 0")
    print("=" * 110)

    columns = [
        "rho",
        "model",
        "mean_difference",
        "percentage_change",
        "ci_lower",
        "ci_upper",
    ]

    print(
        dependence_terminal_df[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Save numerical results
    # --------------------------------------------------

    tables = {
        "dependence_terminal": terminal_df,

        "dependence_thresholds": threshold_df,

        "dependence_theoretical": theoretical_df,

        "dependence_validation": validation_df,

        "dependence_paired_terminal":
            paired_terminal_df,

        "dependence_paired_thresholds":
            paired_threshold_df,

        "dependence_vs_independence_terminal":
            dependence_terminal_df,

        "dependence_vs_independence_thresholds":
            dependence_threshold_df,
    }

    for filename, dataframe in tables.items():

        dataframe.to_csv(
            RESULTS_DIR / f"{filename}.csv",
            index=False,
        )

    # --------------------------------------------------
    # Save terminal and maximum intensity samples
    # --------------------------------------------------

    terminal_archive = {}

    maximum_archive = {}

    for model, rho in terminal_samples:

        key = f"{model}_{rho:.2f}"

        terminal_archive[key] = terminal_samples[
            (model, rho)
        ]

        maximum_archive[key] = maximum_samples[
            (model, rho)
        ]

    np.savez_compressed(
        RESULTS_DIR
        / "dependence_terminal_samples.npz",
        **terminal_archive,
    )

    np.savez_compressed(
        RESULTS_DIR
        / "dependence_maximum_samples.npz",
        **maximum_archive,
    )

    print("\nAll dependence analysis results saved.")

    print("Q6 dependence experiment completed.")


if __name__ == "__main__":
    main()
    
    
# ==========================================================
# Q6 — ENVIRONMENTAL DEPENDENCE EXPERIMENT RESULTS
# ==========================================================
#
# Experimental configuration:
#
# Correlation scenarios: 9
# Wildfire intensity models: 2
# Trajectories per configuration: 10000
# Total simulated trajectories: 180000
# Simulation horizon: 10.0
# Time step: 0.01
# Random seed: 42
#
# Correlation values:
# -1.00, -0.75, -0.50, -0.25, 0.00,
#  0.25,  0.50,  0.75,  1.00
#
# Thresholds:
# 1000, 2000, 3000, 4000, 5000 kW/m
#
# ----------------------------------------------------------
# 1. Terminal intensity results
# ----------------------------------------------------------
#
# Rho      Additive Mean     Interaction Mean
# -1.00    1711.13756        2355.27896
# -0.75    1655.39002        2132.13870
# -0.50    1597.49894        1932.72588
# -0.25    1540.83117        1755.55721
#  0.00    1485.38922        1596.96718
#  0.25    1431.01229        1454.09247
#  0.50    1377.48175        1324.70031
#  0.75    1324.40881        1206.87803
#  1.00    1268.32763        1096.28050
#
# Increasing environmental correlation decreases
# empirical mean intensity in both models.
#
# The interaction model exhibits substantially greater
# sensitivity to environmental dependence.
#
# From rho = -1 to rho = 1:
#
# Additive mean reduction:       25.88%
# Interaction mean reduction:    53.45%
# Interaction SD reduction:      84.16%
# Interaction Q97.5 reduction:   75.25%
#
# ----------------------------------------------------------
# 2. Theoretical validation
# ----------------------------------------------------------
#
# Additive expected log-intensity:
#
# E[log(I_T)] = 7.02569
#
# Additive log-intensity variance:
#
# Var(log(I_T)) =
#     0.54556 - 0.31639 * rho
#
# Interaction expected log-intensity:
#
# E[log(I_T)] =
#     7.07184 - 0.19882 * rho
#
# All 18 theoretical log-intensity expectations fall
# within their empirical pointwise 95% confidence
# intervals.
#
# Maximum absolute standardized error:
#     1.57347
#
# All 9 theoretical paired log-interaction effects
# also fall within their corresponding confidence
# intervals.
#
# These results support numerical consistency with
# the analytical continuous-time model.
#
# ----------------------------------------------------------
# 3. Threshold crossing probabilities
# ----------------------------------------------------------
#
# Threshold K = 1000:
#
# Rho      Additive       Interaction
# -1.00    0.72290        0.78250
#  0.00    0.73770        0.75590
#  1.00    0.77810        0.71170
#
# Threshold K = 2000:
#
# Rho      Additive       Interaction
# -1.00    0.33970        0.42870
#  0.00    0.29240        0.31970
#  1.00    0.17020        0.11560
#
# Threshold K = 5000:
#
# Rho      Additive       Interaction
# -1.00    0.06420        0.12250
#  0.00    0.02720        0.04080
#  1.00    0.00140        0.00030
#
# Increasing correlation decreases empirical crossing
# probabilities for the interaction model at every
# investigated threshold.
#
# For the additive model, the lowest threshold
# exhibits the opposite pattern:
#
# P(tau_1000 <= T) increases with rho.
#
# Higher additive thresholds exhibit decreasing
# crossing probabilities as rho increases.
#
# Therefore, the effect of environmental dependence
# on crossing risk is threshold-dependent.
#
# Very small crossing probabilities have substantial
# relative Monte Carlo uncertainty.
#
# ----------------------------------------------------------
# 4. Paired interaction effects
# ----------------------------------------------------------
#
# Rho      Mean Difference     Percentage Change
# -1.00     644.14140           37.64405%
# -0.75     476.74868           28.79978%
# -0.50     335.22694           20.98449%
# -0.25     214.72604           13.93573%
#  0.00     111.57796            7.51170%
#  0.25      23.08018            1.61286%
#  0.50     -52.78144           -3.83173%
#  0.75    -117.53078           -8.87421%
#  1.00    -172.04713          -13.56488%
#
# All nine pointwise paired confidence intervals
# for arithmetic mean differences exclude zero.
#
# Negative environmental correlation amplifies
# the positive interaction effect.
#
# Sufficiently positive environmental correlation
# reverses the arithmetic mean interaction effect.
#
# Theoretical expected log-interaction effect:
#
# Delta_log(rho) =
#     0.04615 - 0.19882 * rho
#
# Analytical log-effect sign change:
#     rho approximately 0.2321
#
# Empirical arithmetic mean effect changes sign
# between rho = 0.25 and rho = 0.50.
#
# Linear interpolation estimates the arithmetic
# mean sign change near rho = 0.3261.
#
# This interpolation is not an exact theoretical
# critical correlation.
#
# The two sign changes differ because arithmetic
# mean intensity depends on the distribution of
# log-intensity, not only its expectation.
#
# ----------------------------------------------------------
# 5. Conditional first passage times
# ----------------------------------------------------------
#
# Threshold K = 2000:
#
# Rho      Additive       Interaction
# -1.00    6.26835        5.82382
#  0.00    6.74927        6.35238
#  1.00    7.49562        7.36026
#
# Increasing correlation increases the conditional
# mean first passage time at this threshold.
#
# These expectations are conditional on crossing
# before the simulation horizon.
#
# They do not imply that every individual trajectory
# crosses earlier or later.
#
# ----------------------------------------------------------
# 6. Effect relative to environmental independence
# ----------------------------------------------------------
#
# Rho      Additive Change    Interaction Change
# -1.00     15.19792%          47.48449%
# -0.50      7.54750%          21.02477%
#  0.00      0.00000%           0.00000%
#  0.50     -7.26459%         -17.04900%
#  1.00    -14.61311%         -31.35235%
#
# Environmental independence is a consequential
# modeling assumption.
#
# Moderate environmental dependence can substantially
# change simulated wildfire intensity and risk.
#
# ----------------------------------------------------------
# 7. Limitations
# ----------------------------------------------------------
#
# Brownian correlation is a scenario parameter,
# not a calibrated historical estimate.
#
# Historical cross-sectional environmental correlation
# is not equivalent to Brownian innovation correlation.
#
# Common random numbers introduce dependence between
# estimates across correlation scenarios.
#
# Confidence intervals are pointwise rather than
# simultaneous confidence intervals.
#
# Numerical threshold crossing probabilities may
# underestimate continuous-time crossing probabilities
# because crossings between observation times can
# be missed.
#
# Conditional first passage time estimates may be
# unreliable when very few trajectories cross.
#
# Interaction model distribution statistics are
# empirical Monte Carlo estimates rather than
# exact analytical moments.
#
# ----------------------------------------------------------
# 8. Conclusion
# ----------------------------------------------------------
#
# Environmental dependence substantially affects
# terminal wildfire intensity, variability, extreme
# quantiles, and threshold crossing probabilities.
#
# The nonlinear interaction term increases sensitivity
# to environmental dependence and can change from
# amplifying intensity to suppressing intensity.
#
# The numerical results reproduce the analytical
# expected log-intensity relationships and the
# original independent-model baseline.
#
# All mathematical consistency checks passed.
#
# The experiment supports the internal mathematical
# and numerical consistency of the assumed model,
# not its predictive accuracy for real wildfires.
# ==========================================================