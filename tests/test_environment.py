import numpy as np

from src.environment import simulate_ou


def test_initial_condition():
    rng = np.random.default_rng(42)

    times, z = simulate_ou(
        z0=1.5,
        kappa=0.8,
        sigma=1.0,
        T=10.0,
        dt=0.01,
        rng=rng
    )

    assert z[0] == 1.5
    assert len(times) == len(z)
    assert np.all(np.isfinite(z))


def test_deterministic_solution():
    rng = np.random.default_rng(42)

    times, z = simulate_ou(
        z0=2.0,
        kappa=0.8,
        sigma=0.0,
        T=10.0,
        dt=0.01,
        rng=rng
    )

    expected = 2.0 * np.exp(-0.8 * times)

    assert np.allclose(z, expected)
