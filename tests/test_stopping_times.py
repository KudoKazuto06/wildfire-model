import numpy as np
import pytest

from src.stopping_times import (
    first_passage_times,
    threshold_statistics
)


def test_first_passage_times():
    """
    Verify first-passage times for trajectories
    with different crossing behaviors.
    """

    times = np.array([0, 1, 2, 3, 4])

    paths = np.array([
        [500, 800, 1200, 900, 700],
        [500, 600, 800, 900, 950],
        [500, 1100, 1300, 1500, 1600],
        [500, 700, 900, 1000, 1200]
    ])

    result = first_passage_times(
        times, paths, threshold=1000
    )

    expected = np.array([2.0, np.inf, 1.0, 3.0])

    np.testing.assert_array_equal(result, expected)


def test_initial_threshold_crossing():
    """
    If initial intensity already reaches the threshold,
    the first-passage time must be zero.
    """

    times = np.array([0, 1, 2])

    paths = np.array([
        [1000, 800, 600],
        [1200, 900, 700]
    ])

    result = first_passage_times(
        times, paths, threshold=1000
    )

    np.testing.assert_array_equal(
        result, np.array([0.0, 0.0])
    )


def test_no_threshold_crossing():
    """
    If no trajectory crosses the threshold,
    all observed first-passage times must be infinite.
    """

    times = np.array([0, 1, 2])

    paths = np.array([
        [500, 600, 700],
        [500, 800, 900]
    ])

    result = first_passage_times(
        times, paths, threshold=1000
    )

    assert np.all(np.isinf(result))


def test_exact_threshold_crossing():
    """
    Reaching the threshold exactly counts as a crossing
    because the mathematical condition is I_t >= K.
    """

    times = np.array([0, 1, 2])

    paths = np.array([
        [500, 1000, 800]
    ])

    result = first_passage_times(
        times, paths, threshold=1000
    )

    assert result[0] == 1.0


def test_threshold_statistics():
    """
    Validate crossing probabilities and conditional
    first-passage time statistics.
    """

    times = np.array([0, 1, 2, 3, 4])

    paths = np.array([
        [500, 800, 1200, 900, 700],
        [500, 600, 800, 900, 950],
        [500, 1100, 1300, 1500, 1600],
        [500, 700, 900, 1000, 1200]
    ])

    stats = threshold_statistics(
        times, paths, threshold=1000
    )

    assert stats["n_paths"] == 4
    assert stats["n_crossed"] == 3

    assert np.isclose(
        stats["crossing_probability"], 0.75
    )

    assert np.isclose(
        stats["terminal_probability"], 0.50
    )

    assert np.isclose(
        stats["mean_passage_time"], 2.0
    )

    assert np.isclose(
        stats["median_passage_time"], 2.0
    )

    assert np.isclose(stats["q025"], 1.05)
    assert np.isclose(stats["q975"], 2.95)


def test_terminal_probability_bound():
    """
    Terminal exceedance probability cannot exceed
    pathwise threshold crossing probability.
    """

    times = np.array([0, 1, 2, 3])

    paths = np.array([
        [500, 1500, 800, 600],
        [500, 800, 1200, 1500],
        [500, 600, 700, 800]
    ])

    stats = threshold_statistics(
        times, paths, threshold=1000
    )

    assert (
        stats["terminal_probability"]
        <= stats["crossing_probability"]
    )


def test_threshold_monotonicity():
    """
    Higher thresholds cannot have higher crossing
    probabilities when evaluated on the same paths.
    """

    times = np.array([0, 1, 2, 3])

    paths = np.array([
        [500, 1200, 2500, 800],
        [500, 800, 1500, 1800],
        [500, 600, 700, 900]
    ])

    thresholds = [1000, 2000, 3000]

    probabilities = [
        threshold_statistics(
            times, paths, threshold
        )["crossing_probability"]
        for threshold in thresholds
    ]

    assert np.all(np.diff(probabilities) <= 0)


def test_no_crossing_statistics():
    """
    Conditional first-passage statistics are undefined
    when no trajectory reaches the threshold.
    """

    times = np.array([0, 1, 2])

    paths = np.array([
        [500, 600, 700],
        [500, 800, 900]
    ])

    stats = threshold_statistics(
        times, paths, threshold=1000
    )

    assert stats["n_crossed"] == 0
    assert stats["crossing_probability"] == 0.0
    assert stats["terminal_probability"] == 0.0

    assert np.isnan(stats["mean_passage_time"])
    assert np.isnan(stats["median_passage_time"])
    assert np.isnan(stats["q025"])
    assert np.isnan(stats["q975"])


def test_invalid_threshold():
    """
    Wildfire intensity thresholds must be positive.
    """

    times = np.array([0, 1, 2])

    paths = np.array([
        [500, 600, 700]
    ])

    with pytest.raises(ValueError):
        first_passage_times(
            times, paths, threshold=-100
        )


def test_invalid_dimensions():
    """
    The number of observations in each trajectory
    must match the number of time points.
    """

    times = np.array([0, 1, 2])

    paths = np.array([
        [500, 600]
    ])

    with pytest.raises(ValueError):
        first_passage_times(
            times, paths, threshold=1000
        )