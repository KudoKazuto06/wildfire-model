"""
Unit tests for seasonal environmental scenarios
and analytical intensity expectations.
"""

import numpy as np

from src.config import (
    WIND_MEAN,
    WIND_STD,
    HUMIDITY_MEAN,
    HUMIDITY_STD,
    KAPPA_W,
    KAPPA_H,
    GAMMA,
    I0,
    T,
)

from src.seasonal import (
    load_seasonal_data,
    standardize_environment,
    build_seasonal_scenarios,
    expected_log_intensity,
    theoretical_seasonal_comparison,
)

from src.theoretical import (
    theoretical_intensity_moments,
)


def test_seasonal_data():

    data = load_seasonal_data()

    assert len(data) == 4

    assert set(data["season"]) == {
        "Winter",
        "Spring",
        "Summer",
        "Fall",
    }

    assert data["n_observations"].sum() == 24859


def test_standardization_at_mean():

    wind_z, humidity_z = standardize_environment(
        WIND_MEAN,
        HUMIDITY_MEAN,
    )

    np.testing.assert_allclose(
        [wind_z, humidity_z],
        [0.0, 0.0],
    )


def test_standardization_one_std():

    wind_z, humidity_z = standardize_environment(
        WIND_MEAN + WIND_STD,
        HUMIDITY_MEAN - HUMIDITY_STD,
    )

    np.testing.assert_allclose(
        [wind_z, humidity_z],
        [1.0, -1.0],
    )


def test_seasonal_anomaly_signs():

    scenarios = (
        build_seasonal_scenarios()
        .set_index("season")
    )

    assert scenarios.loc["Spring", "wind_z0"] > 0

    assert scenarios.loc["Spring", "humidity_z0"] < 0

    assert scenarios.loc["Summer", "wind_z0"] < 0

    assert scenarios.loc["Summer", "humidity_z0"] > 0


def test_initial_interaction_contribution():

    scenarios = build_seasonal_scenarios()

    expected = (
        GAMMA
        * scenarios["wind_z0"]
        * scenarios["humidity_z0"]
    )

    actual = (
        scenarios["interaction_growth_rate"]
        - scenarios["additive_growth_rate"]
    )

    np.testing.assert_allclose(
        actual,
        expected,
        atol=1e-12,
    )


def test_expected_log_intensity_at_zero():

    result = expected_log_intensity(
        wind_z0=1.5,
        humidity_z0=-1.0,
        gamma=GAMMA,
        horizon=0.0,
    )

    np.testing.assert_allclose(
        result,
        np.log(I0),
    )


def test_additive_theoretical_consistency():

    scenarios = build_seasonal_scenarios()

    for _, row in scenarios.iterrows():

        wind_z0 = row["wind_z0"]

        humidity_z0 = row["humidity_z0"]

        actual = expected_log_intensity(
            wind_z0=wind_z0,
            humidity_z0=humidity_z0,
            gamma=0.0,
        )

        expected = theoretical_intensity_moments(
            wind_z0=wind_z0,
            humidity_z0=humidity_z0,
        )["log_mean"]

        np.testing.assert_allclose(
            actual,
            expected,
            atol=1e-12,
        )


def test_analytical_interaction_effect():

    wind_z0 = 1.5

    humidity_z0 = -1.0

    combined_kappa = KAPPA_W + KAPPA_H

    expected = (
        GAMMA
        * wind_z0
        * humidity_z0
        * (
            1 - np.exp(-combined_kappa * T)
        )
        / combined_kappa
    )

    actual = (
        expected_log_intensity(
            wind_z0,
            humidity_z0,
            gamma=GAMMA,
        )
        -
        expected_log_intensity(
            wind_z0,
            humidity_z0,
            gamma=0.0,
        )
    )

    np.testing.assert_allclose(
        actual,
        expected,
        atol=1e-12,
    )


def test_additive_log_variance_constant():

    results = theoretical_seasonal_comparison()

    variances = results[
        "additive_log_variance"
    ].to_numpy()

    np.testing.assert_allclose(
        variances,
        variances[0],
        atol=1e-12,
    )


def test_theoretical_seasonal_ranking():

    results = theoretical_seasonal_comparison()

    ranking = results.sort_values(
        "additive_mean",
        ascending=False,
    )["season"].tolist()

    assert ranking == [
        "Spring",
        "Fall",
        "Winter",
        "Summer",
    ]