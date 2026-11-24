"""
Seasonal environmental scenarios for Q5.

Historical seasonal averages are converted into
standardized initial environmental conditions.

All dynamic model parameters remain unchanged.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from src.config import (
    WIND_MEAN,
    WIND_STD,
    HUMIDITY_MEAN,
    HUMIDITY_STD,
    ALPHA_C,
    BETA_W,
    BETA_H,
    GAMMA,
)

from src.intensity import growth_rate

from src.config import (
    KAPPA_W,
    KAPPA_H,
    SIGMA_I,
    I0,
    T,
)

from src.theoretical import theoretical_intensity_moments


# Resolve the project directory independently
# of the current working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[1]

SEASONAL_DATA_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "seasonal_summary.csv"
)


def load_seasonal_data(path=SEASONAL_DATA_PATH):
    """
    Load and validate historical seasonal summaries.
    """

    data = pd.read_csv(path)

    required_columns = {
        "season",
        "n_observations",
        "wind_mean",
        "humidity_mean",
        "spread_mean",
    }

    if not required_columns.issubset(data.columns):
        raise ValueError(
            "Missing required seasonal data columns"
        )

    expected_seasons = {
        "Winter",
        "Spring",
        "Summer",
        "Fall",
    }

    if (
        len(data) != 4
        or set(data["season"]) != expected_seasons
    ):
        raise ValueError(
            "Exactly four unique seasons are required"
        )

    numeric_columns = [
        "n_observations",
        "wind_mean",
        "humidity_mean",
        "spread_mean",
    ]

    if not np.all(
        np.isfinite(data[numeric_columns].to_numpy())
    ):
        raise ValueError(
            "Seasonal statistics must be finite"
        )

    if (data["n_observations"] <= 0).any():
        raise ValueError(
            "Observation counts must be positive"
        )

    if (data["wind_mean"] < 0).any():
        raise ValueError(
            "Wind speed cannot be negative"
        )

    if not data["humidity_mean"].between(0, 100).all():
        raise ValueError(
            "Relative humidity must be between 0 and 100"
        )

    if (data["spread_mean"] < 0).any():
        raise ValueError(
            "Spread rate cannot be negative"
        )

    return data


def standardize_environment(wind, humidity):
    """
    Convert physical environmental measurements
    into standardized anomalies.

    Mathematical reference: Equation (2.1).
    """

    wind_z = (wind - WIND_MEAN) / WIND_STD

    humidity_z = (
        humidity - HUMIDITY_MEAN
    ) / HUMIDITY_STD

    return wind_z, humidity_z


def build_seasonal_scenarios(path=SEASONAL_DATA_PATH):
    """
    Construct standardized seasonal initial conditions
    and calculate their initial intensity growth rates.
    """

    data = load_seasonal_data(path).copy()

    # Convert seasonal physical averages into anomalies.
    data["wind_z0"], data["humidity_z0"] = (
        standardize_environment(
            data["wind_mean"],
            data["humidity_mean"],
        )
    )

    # Calculate the initial additive growth rate.
    data["additive_growth_rate"] = growth_rate(
        data["wind_z0"],
        data["humidity_z0"],
        ALPHA_C,
        BETA_W,
        BETA_H,
        0.0,
    )

    # Calculate the initial interaction growth rate.
    data["interaction_growth_rate"] = growth_rate(
        data["wind_z0"],
        data["humidity_z0"],
        ALPHA_C,
        BETA_W,
        BETA_H,
        GAMMA,
    )

    # Isolate the initial interaction contribution.
    data["interaction_contribution"] = (
        GAMMA
        * data["wind_z0"]
        * data["humidity_z0"]
    )

    return data


def expected_log_intensity(
    wind_z0,
    humidity_z0,
    gamma=GAMMA,
    horizon=T,
):
    """
    Calculate the analytical expected terminal
    log-intensity under independent OU processes.

    Valid for both additive and interaction models.

    Mathematical references:
        Section 3.1: OU process expectations
        Section 3.2: Log-intensity representation
    """

    if horizon < 0:
        raise ValueError("Horizon cannot be negative")

    # Initial log-intensity.
    initial_term = np.log(I0)

    # Baseline growth with Ito correction.
    baseline_term = (
        ALPHA_C - 0.5 * SIGMA_I**2
    ) * horizon

    # Expected integrated wind contribution.
    wind_term = (
        BETA_W
        * wind_z0
        * (1 - np.exp(-KAPPA_W * horizon))
        / KAPPA_W
    )

    # Expected integrated humidity contribution.
    humidity_term = (
        -BETA_H
        * humidity_z0
        * (1 - np.exp(-KAPPA_H * horizon))
        / KAPPA_H
    )

    # Combined environmental mean reversion.
    combined_kappa = KAPPA_W + KAPPA_H

    # Expected integrated wind-humidity interaction.
    interaction_term = (
        gamma
        * wind_z0
        * humidity_z0
        * (1 - np.exp(-combined_kappa * horizon))
        / combined_kappa
    )

    # Brownian motion has zero expectation.
    log_mean = (
        initial_term
        + baseline_term
        + wind_term
        + humidity_term
        + interaction_term
    )

    return log_mean


def theoretical_seasonal_comparison(
    horizon=T,
):
    """
    Calculate analytical seasonal statistics.

    Additive model:
        Exact lognormal distribution.

    Interaction model:
        Exact expected log-intensity, but not
        necessarily a lognormal distribution.
    """

    scenarios = build_seasonal_scenarios()

    results = []

    for _, row in scenarios.iterrows():

        season = row["season"]

        wind_z0 = row["wind_z0"]
        humidity_z0 = row["humidity_z0"]

        # Exact additive lognormal statistics
        # derived previously in Q2.
        additive = theoretical_intensity_moments(
            wind_z0=wind_z0,
            humidity_z0=humidity_z0,
            horizon=horizon,
        )

        # Expected log-intensity with interaction.
        interaction_log_mean = expected_log_intensity(
            wind_z0=wind_z0,
            humidity_z0=humidity_z0,
            gamma=GAMMA,
            horizon=horizon,
        )

        # Difference in expected log-intensity.
        log_difference = (
            interaction_log_mean
            - additive["log_mean"]
        )

        results.append(
            {
                "season": season,
                "wind_z0": wind_z0,
                "humidity_z0": humidity_z0,

                # Exact additive statistics.
                "additive_log_mean": additive["log_mean"],
                "additive_log_variance": additive["log_variance"],
                "additive_mean": additive["mean"],
                "additive_std": additive["std"],
                "additive_median": additive["median"],

                # Interaction model statistics.
                "interaction_log_mean": interaction_log_mean,

                # Geometric mean, not arithmetic mean.
                "interaction_geometric_mean": np.exp(
                    interaction_log_mean
                ),

                # Analytical interaction effect.
                "log_interaction_effect": log_difference,
            }
        )

    return pd.DataFrame(results)


if __name__ == "__main__":

    scenarios = build_seasonal_scenarios()

    print("HISTORICAL SEASONAL DATA")
    print("=" * 85)

    print(
        scenarios[
            [
                "season",
                "n_observations",
                "wind_mean",
                "humidity_mean",
                "spread_mean",
            ]
        ].to_string(index=False)
    )

    print("\nSTANDARDIZED SEASONAL CONDITIONS")
    print("=" * 85)

    print(
        scenarios[
            [
                "season",
                "wind_z0",
                "humidity_z0",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    print("\nINITIAL INTENSITY GROWTH RATES")
    print("=" * 100)

    print(
        scenarios[
            [
                "season",
                "additive_growth_rate",
                "interaction_growth_rate",
                "interaction_contribution",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    print(
        "\nTotal historical observations:",
        scenarios["n_observations"].sum(),
    )
    
    print("\nTHEORETICAL SEASONAL COMPARISON")
    print("=" * 110)

    theoretical = theoretical_seasonal_comparison()

    print(
        theoretical[
            [
                "season",
                "additive_log_mean",
                "additive_mean",
                "additive_std",
                "additive_median",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    print("\nTHEORETICAL INTERACTION EFFECT")
    print("=" * 110)

    print(
        theoretical[
            [
                "season",
                "additive_log_mean",
                "interaction_log_mean",
                "interaction_geometric_mean",
                "log_interaction_effect",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # Validate the new analytical expression against
    # the existing additive benchmark from Q2.
    for _, row in theoretical.iterrows():

        calculated = expected_log_intensity(
            wind_z0=row["wind_z0"],
            humidity_z0=row["humidity_z0"],
            gamma=0.0,
        )

        np.testing.assert_allclose(
            calculated,
            row["additive_log_mean"],
            rtol=1e-10,
            atol=1e-10,
        )

    # At time zero, log-intensity must equal log(I0).
    np.testing.assert_allclose(
        expected_log_intensity(
            wind_z0=1.0,
            humidity_z0=-1.0,
            horizon=0.0,
        ),
        np.log(I0),
    )

    print("\nAll analytical consistency checks passed.")