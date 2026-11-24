import numpy as np


def first_passage_times(times, intensity_paths, threshold):
    """
    Calculate the discrete first-passage time for each
    simulated wildfire intensity trajectory.

    Mathematical definition:
        tau_K = inf{t >= 0 : I_t >= K}

    Numerical approximation:
        tau_K = min{t_n : I(t_n) >= K}

    Parameters
    ----------
    times : np.ndarray
        Array of simulation time points.

    intensity_paths : np.ndarray
        Matrix of simulated intensity trajectories.
        Shape: (number of paths, number of time points).

    threshold : float
        Wildfire intensity threshold K.

    Returns
    -------
    passage_times : np.ndarray
        First observed crossing time for each trajectory.
        np.inf indicates no observed crossing.
    """

    times = np.asarray(times)
    intensity_paths = np.asarray(intensity_paths)

    # Validate input dimensions.
    if times.ndim != 1 or len(times) == 0:
        raise ValueError("Times must be a nonempty 1D array")

    if intensity_paths.ndim != 2:
        raise ValueError("Intensity paths must be a 2D array")

    if intensity_paths.shape[1] != len(times):
        raise ValueError("Number of time points must match")

    if intensity_paths.shape[0] == 0:
        raise ValueError("At least one trajectory is required")

    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("Threshold must be finite and positive")

    if not np.all(np.isfinite(times)):
        raise ValueError("Times must be finite")

    if not np.all(np.diff(times) > 0):
        raise ValueError("Times must be strictly increasing")

    if not np.all(np.isfinite(intensity_paths)):
        raise ValueError("Intensity paths must be finite")

    if np.any(intensity_paths <= 0):
        raise ValueError("Intensity must remain positive")

    # Identify every time point at which intensity
    # reaches or exceeds the threshold.
    crossed = intensity_paths >= threshold

    # Determine whether each trajectory crosses
    # the threshold at least once.
    has_crossed = np.any(crossed, axis=1)

    # Find the index of the first crossing.
    first_indices = np.argmax(crossed, axis=1)

    # Initially assume that no trajectory crosses.
    passage_times = np.full(
        intensity_paths.shape[0],
        np.inf
    )

    # Assign crossing times only to trajectories
    # that actually reached the threshold.
    passage_times[has_crossed] = times[
        first_indices[has_crossed]
    ]

    return passage_times


def threshold_statistics(times, intensity_paths, threshold):
    """
    Calculate threshold crossing and first-passage
    statistics for a collection of wildfire trajectories.
    """

    passage_times = first_passage_times(
        times,
        intensity_paths,
        threshold
    )

    n_paths = len(passage_times)

    # A finite first-passage time means that the
    # threshold was reached during the simulation.
    crossed = np.isfinite(passage_times)

    n_crossed = np.sum(crossed)

    # Estimate P(tau_K <= T).
    crossing_probability = np.mean(crossed)

    # Estimate P(I_T >= K).
    terminal_probability = np.mean(
        intensity_paths[:, -1] >= threshold
    )

    # Estimate Monte Carlo sampling uncertainty.
    standard_error = np.sqrt(
        crossing_probability
        * (1 - crossing_probability)
        / n_paths
    )

    # Approximate 95% Monte Carlo confidence interval.
    ci_lower = max(
        0.0,
        crossing_probability - 1.96 * standard_error
    )

    ci_upper = min(
        1.0,
        crossing_probability + 1.96 * standard_error
    )

    # First-passage statistics are conditional on
    # reaching the threshold within the time horizon.
    if n_crossed > 0:

        successful_times = passage_times[crossed]

        mean_passage_time = np.mean(successful_times)
        median_passage_time = np.median(successful_times)

        q025 = np.quantile(successful_times, 0.025)
        q975 = np.quantile(successful_times, 0.975)

    else:

        mean_passage_time = np.nan
        median_passage_time = np.nan
        q025 = np.nan
        q975 = np.nan

    return {
        "threshold": threshold,
        "n_paths": n_paths,
        "n_crossed": n_crossed,
        "crossing_probability": crossing_probability,
        "terminal_probability": terminal_probability,
        "standard_error": standard_error,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "mean_passage_time": mean_passage_time,
        "median_passage_time": median_passage_time,
        "q025": q025,
        "q975": q975
    }


if __name__ == "__main__":

    # Simple deterministic example for validation.
    times = np.array([0, 1, 2, 3, 4])

    intensity_paths = np.array([
        [500, 800, 1200, 900, 700],
        [500, 600, 800, 900, 950],
        [500, 1100, 1300, 1500, 1600],
        [500, 700, 900, 1000, 1200]
    ])

    threshold = 1000

    passage_times = first_passage_times(
        times,
        intensity_paths,
        threshold
    )

    print("First-passage times:")
    print(passage_times)

    print("\nThreshold statistics")
    print("-" * 45)

    statistics = threshold_statistics(
        times,
        intensity_paths,
        threshold
    )

    for name, value in statistics.items():
        print(f"{name:25s}: {value}")