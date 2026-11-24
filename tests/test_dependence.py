import numpy as np
import pytest

from src.config import (
    KAPPA_W,
    KAPPA_H,
    SIGMA_W,
    SIGMA_H,
    SIGMA_I,
    BETA_W,
    BETA_H,
    GAMMA,
    T,
    SEED,
)

from src.dependence import validate_correlation
from src.monte_carlo import simulate_many


# ============================================================
# 1. Correlation parameter validation
# ============================================================

@pytest.mark.parametrize("rho", [-1.0, -0.5, 0.0, 0.5, 1.0])
def test_valid_correlation(rho):
    """Valid Brownian correlations should be accepted."""
    validate_correlation(rho)


@pytest.mark.parametrize("rho", [-1.1, 1.1, -2.0, 2.0])
def test_invalid_correlation(rho):
    """Correlations outside [-1, 1] should be rejected."""
    with pytest.raises(ValueError):
        validate_correlation(rho)


# ============================================================
# 2. Backward compatibility
# ============================================================

def test_default_equals_independence():
    """Default simulation must equal explicit rho=0."""

    times_default, paths_default = simulate_many(
        n_paths=100,
        horizon=1.0,
        seed=SEED,
    )

    times_independent, paths_independent = simulate_many(
        n_paths=100,
        horizon=1.0,
        seed=SEED,
        rho=0.0,
    )

    np.testing.assert_array_equal(
        times_default,
        times_independent,
    )

    np.testing.assert_array_equal(
        paths_default,
        paths_independent,
    )


# ============================================================
# 3. Boundary correlations and simulation validity
# ============================================================

@pytest.mark.parametrize("rho", [-1.0, 0.0, 1.0])
def test_simulation_boundary_correlations(rho):
    """Simulations must remain valid at extreme correlations."""

    times, paths = simulate_many(
        n_paths=100,
        horizon=1.0,
        dt=0.01,
        seed=SEED,
        rho=rho,
    )

    assert paths.shape == (100, len(times))
    assert np.all(np.isfinite(paths))
    assert np.all(paths > 0)
    assert np.all(np.diff(times) > 0)


def test_simulation_rejects_invalid_correlation():
    """The Monte Carlo engine must reject invalid rho."""

    with pytest.raises(ValueError):
        simulate_many(
            n_paths=10,
            horizon=1.0,
            rho=1.5,
        )


# ============================================================
# 4. Reproducibility
# ============================================================

def test_correlated_simulation_reproducibility():
    """Identical parameters and seeds produce identical paths."""

    times_1, paths_1 = simulate_many(
        n_paths=100,
        horizon=1.0,
        rho=-0.5,
        seed=SEED,
    )

    times_2, paths_2 = simulate_many(
        n_paths=100,
        horizon=1.0,
        rho=-0.5,
        seed=SEED,
    )

    np.testing.assert_array_equal(times_1, times_2)
    np.testing.assert_array_equal(paths_1, paths_2)


def test_correlation_changes_trajectories():
    """Changing environmental dependence must affect intensity."""

    _, independent_paths = simulate_many(
        n_paths=100,
        horizon=1.0,
        rho=0.0,
        seed=SEED,
    )

    _, correlated_paths = simulate_many(
        n_paths=100,
        horizon=1.0,
        rho=-0.5,
        seed=SEED,
    )

    assert not np.array_equal(
        independent_paths,
        correlated_paths,
    )


# ============================================================
# 5. Independent analytical reference calculations
# ============================================================

def integrated_covariance_reference(rho, horizon):
    """Continuous-time covariance of integrated OU processes."""

    kw = KAPPA_W
    kh = KAPPA_H

    aw = -np.expm1(-kw * horizon) / kw
    ah = -np.expm1(-kh * horizon) / kh

    awh = (
        -np.expm1(-(kw + kh) * horizon)
        / (kw + kh)
    )

    return (
        rho * SIGMA_W * SIGMA_H / (kw * kh)
        * (horizon - aw - ah + awh)
    )


def integrated_variance_reference(kappa, sigma, horizon):
    """Continuous-time variance of an integrated OU process."""

    return (
        sigma**2 / kappa**2
        * (
            horizon
            - 2 * (-np.expm1(-kappa * horizon)) / kappa
            + (-np.expm1(-2 * kappa * horizon))
            / (2 * kappa)
        )
    )


def additive_log_variance_reference(rho, horizon):
    """Theoretical log-intensity variance for gamma=0."""

    vw = integrated_variance_reference(
        KAPPA_W, SIGMA_W, horizon
    )

    vh = integrated_variance_reference(
        KAPPA_H, SIGMA_H, horizon
    )

    covariance = integrated_covariance_reference(
        rho, horizon
    )

    return (
        SIGMA_I**2 * horizon
        + BETA_W**2 * vw
        + BETA_H**2 * vh
        - 2 * BETA_W * BETA_H * covariance
    )


def log_interaction_reference(rho, horizon):
    """Expected paired log-intensity interaction effect."""

    wind_z0 = 1.5
    humidity_z0 = -1.0

    k = KAPPA_W + KAPPA_H

    a = -np.expm1(-k * horizon) / k

    initial_effect = wind_z0 * humidity_z0 * a

    dependence_effect = (
        rho * SIGMA_W * SIGMA_H / k
        * (horizon - a)
    )

    return GAMMA * (
        initial_effect + dependence_effect
    )


# ============================================================
# 6. Analytical reference identities
# ============================================================

def test_integrated_covariance_at_independence():
    assert integrated_covariance_reference(0.0, T) == 0.0


def test_integrated_covariance_at_initial_time():
    assert integrated_covariance_reference(0.5, 0.0) == 0.0


def test_integrated_covariance_symmetry():
    negative = integrated_covariance_reference(-0.5, T)
    positive = integrated_covariance_reference(0.5, T)

    assert negative < 0
    assert positive > 0

    assert negative == pytest.approx(-positive)


def test_additive_log_variance_ordering():
    """Negative correlation increases additive log variance."""

    negative = additive_log_variance_reference(-1.0, T)
    independent = additive_log_variance_reference(0.0, T)
    positive = additive_log_variance_reference(1.0, T)

    assert negative > independent > positive > 0


def test_analytical_critical_correlation():
    """Expected log-interaction effect vanishes at its root."""

    k = KAPPA_W + KAPPA_H

    a = -np.expm1(-k * T) / k

    rho_star = (
        1.5 * a * k
        / (SIGMA_W * SIGMA_H * (T - a))
    )

    assert rho_star == pytest.approx(0.2321, abs=1e-4)

    assert log_interaction_reference(
        rho_star, T
    ) == pytest.approx(0.0, abs=1e-10)


def test_log_interaction_sign_change():
    """The interaction effect reverses with correlation."""

    assert log_interaction_reference(-1.0, T) > 0
    assert log_interaction_reference(0.0, T) > 0
    assert log_interaction_reference(1.0, T) < 0


# ============================================================
# 7. Numerical agreement with analytical references
# ============================================================

@pytest.mark.parametrize("rho", [-1.0, 0.0, 1.0])
def test_additive_log_variance_against_theory(rho):
    """Monte Carlo variance should agree with OU theory."""

    _, paths = simulate_many(
        n_paths=4000,
        horizon=T,
        gamma=0.0,
        rho=rho,
        seed=SEED,
    )

    terminal_log = np.log(paths[:, -1])

    empirical = np.var(terminal_log, ddof=1)

    theoretical = additive_log_variance_reference(
        rho, T
    )

    assert empirical == pytest.approx(
        theoretical,
        rel=0.08,
    )


@pytest.mark.parametrize("rho", [-1.0, 0.0, 1.0])
def test_paired_log_interaction_against_theory(rho):
    """Paired log differences should agree with theory."""

    _, additive_paths = simulate_many(
        n_paths=4000,
        horizon=T,
        gamma=0.0,
        rho=rho,
        seed=SEED,
    )

    _, interaction_paths = simulate_many(
        n_paths=4000,
        horizon=T,
        gamma=GAMMA,
        rho=rho,
        seed=SEED,
    )

    differences = (
        np.log(interaction_paths[:, -1])
        - np.log(additive_paths[:, -1])
    )

    empirical = np.mean(differences)

    standard_error = (
        np.std(differences, ddof=1)
        / np.sqrt(len(differences))
    )

    theoretical = log_interaction_reference(rho, T)

    # Allow for Monte Carlo sampling error and
    # time discretization error.
    tolerance = 5 * standard_error + 0.002

    assert abs(empirical - theoretical) <= tolerance