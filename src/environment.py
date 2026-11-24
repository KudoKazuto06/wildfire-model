
import numpy as np


def simulate_ou(z0, kappa, sigma, T, dt, rng):
    """
    Simulate an Ornstein-Uhlenbeck process.

    dZ_t = -kappa * Z_t dt + sigma * dB_t

    Uses the exact OU transition formula.
    """

    if kappa <= 0:
        raise ValueError("kappa must be positive")

    if sigma < 0:
        raise ValueError("sigma cannot be negative")

    if T <= 0 or dt <= 0:
        raise ValueError("T and dt must be positive")

    n_steps = int(np.ceil(T / dt))
    dt = T / n_steps

    times = np.linspace(0, T, n_steps + 1)

    z = np.zeros(n_steps + 1)
    z[0] = z0

    decay = np.exp(-kappa * dt)

    noise_std = sigma * np.sqrt(
        -np.expm1(-2 * kappa * dt) / (2 * kappa)
    )

    for n in range(n_steps):
        epsilon = rng.normal()

        z[n + 1] = (
            decay * z[n]
            + noise_std * epsilon
        )

    return times, z
