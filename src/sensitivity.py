import numpy as np

from src.config import (
    ALPHA_C,
    BETA_W,
    BETA_H,
    GAMMA,
    SIGMA_I,
    KAPPA_W,
    KAPPA_H,
    SIGMA_W,
    SIGMA_H,
    N_PATHS,
    SEED
)

from src.monte_carlo import simulate_many
from src.stopping_times import threshold_statistics
from src.theoretical import theoretical_log_sensitivity


# ======================================================
# 1. Baseline parameter configuration
# ======================================================

BASELINE_PARAMETERS = {
    "alpha_c": ALPHA_C,
    "beta_w": BETA_W,
    "beta_h": BETA_H,
    "gamma": GAMMA,
    "sigma_i": SIGMA_I,
    "kappa_w": KAPPA_W,
    "kappa_h": KAPPA_H,
    "sigma_w": SIGMA_W,
    "sigma_h": SIGMA_H,
    "wind_z0": 1.5,
    "humidity_z0": -1.0
}


# ======================================================
# 2. Calculate simulation outcome metrics
# ======================================================

def calculate_metrics(
    times,
    intensity_paths,
    threshold=2000
):
    """
    Calculate terminal intensity and first-passage
    statistics for a given simulation configuration.
    """

    final_intensity = intensity_paths[:, -1]

    log_intensity = np.log(final_intensity)

    # Reuse stopping-time calculations from Q3
    crossing_stats = threshold_statistics(
        times,
        intensity_paths,
        threshold
    )

    return {
        "mean": np.mean(final_intensity),

        "std": np.std(
            final_intensity,
            ddof=1
        ),

        "median": np.median(final_intensity),

        "q025": np.quantile(
            final_intensity,
            0.025
        ),

        "q975": np.quantile(
            final_intensity,
            0.975
        ),

        "log_mean": np.mean(log_intensity),

        "log_variance": np.var(
            log_intensity,
            ddof=1
        ),

        "crossing_probability":
            crossing_stats["crossing_probability"],

        "terminal_probability":
            crossing_stats["terminal_probability"],

        "mean_passage_time":
            crossing_stats["mean_passage_time"]
    }


# ======================================================
# 3. Evaluate one parameter configuration
# ======================================================

def evaluate_parameter(
    parameter,
    value,
    n_paths=N_PATHS,
    threshold=2000,
    seed=SEED
):
    """
    Change one parameter while holding all other
    model parameters at their baseline values.
    """

    if parameter not in BASELINE_PARAMETERS:
        raise ValueError(
            f"Unknown parameter: {parameter}"
        )

    # Copy baseline parameters to avoid modifying
    # the original configuration.
    parameters = BASELINE_PARAMETERS.copy()

    # Override only the selected parameter.
    parameters[parameter] = value

    # Simulate wildfire trajectories.
    times, paths = simulate_many(
        n_paths=n_paths,
        seed=seed,
        **parameters
    )

    # Calculate outcome statistics.
    metrics = calculate_metrics(
        times,
        paths,
        threshold
    )

    return {
        "parameter": parameter,
        "value": value,
        **metrics
    }


# ======================================================
# 4. Estimate numerical log-intensity sensitivity
# ======================================================

def numerical_log_sensitivity(
    parameter,
    h=0.001,
    n_paths=N_PATHS,
    seed=SEED
):
    """
    Estimate the derivative of expected terminal
    log-intensity using paired central finite
    differences and common random numbers.
    """

    if parameter not in BASELINE_PARAMETERS:
        raise ValueError(
            f"Unknown parameter: {parameter}"
        )

    if not np.isfinite(h) or h <= 0:
        raise ValueError(
            "Step size h must be finite and positive"
        )

    if n_paths < 2:
        raise ValueError(
            "At least two trajectories are required"
        )

    baseline = BASELINE_PARAMETERS[parameter]

    # Construct upper and lower configurations.
    parameters_plus = BASELINE_PARAMETERS.copy()

    parameters_minus = BASELINE_PARAMETERS.copy()

    parameters_plus[parameter] = baseline + h

    parameters_minus[parameter] = baseline - h

    # --------------------------------------------------
    # Simulate upper parameter configuration
    # --------------------------------------------------

    _, paths_plus = simulate_many(
        n_paths=n_paths,
        seed=seed,
        **parameters_plus
    )

    log_plus = np.log(paths_plus[:, -1])

    # Release the full trajectory matrix.
    del paths_plus

    # --------------------------------------------------
    # Simulate lower parameter configuration
    # --------------------------------------------------

    _, paths_minus = simulate_many(
        n_paths=n_paths,
        seed=seed,
        **parameters_minus
    )

    log_minus = np.log(paths_minus[:, -1])

    del paths_minus

    # --------------------------------------------------
    # Calculate paired finite differences
    # --------------------------------------------------

    differences = (
        log_plus - log_minus
    ) / (2 * h)

    derivative = np.mean(differences)

    standard_error = (
        np.std(differences, ddof=1)
        / np.sqrt(n_paths)
    )

    ci_lower = derivative - 1.96 * standard_error

    ci_upper = derivative + 1.96 * standard_error

    return {
        "parameter": parameter,
        "baseline": baseline,
        "step_size": h,
        "derivative": derivative,
        "standard_error": standard_error,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    }


# ======================================================
# 5. Compare analytical and numerical derivatives
# ======================================================

if __name__ == "__main__":

    theoretical = theoretical_log_sensitivity()

    print("ANALYTICAL VS NUMERICAL SENSITIVITY")
    print("=" * 100)

    print(
        f"{'Parameter':<12}"
        f"{'Theory':>12}"
        f"{'Numerical':>12}"
        f"{'Abs. Error':>14}"
        f"{'Std. Error':>14}"
        f"{'CI Lower':>14}"
        f"{'CI Upper':>14}"
    )

    for parameter, theoretical_value in theoretical.items():

        result = numerical_log_sensitivity(
            parameter=parameter,
            h=0.001
        )

        numerical_value = result["derivative"]

        absolute_error = abs(
            numerical_value - theoretical_value
        )

        print(
            f"{parameter:<12}"
            f"{theoretical_value:12.5f}"
            f"{numerical_value:12.5f}"
            f"{absolute_error:14.5f}"
            f"{result['standard_error']:14.5f}"
            f"{result['ci_lower']:14.5f}"
            f"{result['ci_upper']:14.5f}"
        )

    # --------------------------------------------------
    # Preliminary one-at-a-time sensitivity experiment
    # --------------------------------------------------

    print("\nWIND COEFFICIENT SENSITIVITY")
    print("=" * 80)

    print(
        f"{'Beta W':>10}"
        f"{'Mean':>15}"
        f"{'Std':>15}"
        f"{'Q97.5':>15}"
        f"{'P(Crossing)':>15}"
    )

    for value in [0.06, 0.09, 0.12, 0.15, 0.18]:

        result = evaluate_parameter(
            parameter="beta_w",
            value=value
        )

        print(
            f"{value:10.2f}"
            f"{result['mean']:15.5f}"
            f"{result['std']:15.5f}"
            f"{result['q975']:15.5f}"
            f"{result['crossing_probability']:15.5f}"
        )