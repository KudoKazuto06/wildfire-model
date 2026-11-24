
import numpy as np

from src.environment import simulate_ou
from src.intensity import simulate_intensity

from src.config import (
    KAPPA_W, SIGMA_W,
    KAPPA_H, SIGMA_H,
    ALPHA_C, BETA_W, BETA_H,
    GAMMA, SIGMA_I, I0,
    T, DT, SEED
)


def simulate_wildfire(
    wind_z0=1.5,
    humidity_z0=-1.0,
    seed=SEED,
    gamma=GAMMA
):
    """
    Simulate the coupled wildfire system.

    Returns time, wind anomaly, humidity anomaly,
    and wildfire intensity trajectories.
    """

    # Independent random generators for the three noises
    streams = np.random.SeedSequence(seed).spawn(3)

    wind_rng = np.random.default_rng(streams[0])
    humidity_rng = np.random.default_rng(streams[1])
    intensity_rng = np.random.default_rng(streams[2])

    # 1. Simulate wind
    times, wind = simulate_ou(
        z0=wind_z0,
        kappa=KAPPA_W,
        sigma=SIGMA_W,
        T=T,
        dt=DT,
        rng=wind_rng
    )

    # 2. Simulate humidity
    _, humidity = simulate_ou(
        z0=humidity_z0,
        kappa=KAPPA_H,
        sigma=SIGMA_H,
        T=T,
        dt=DT,
        rng=humidity_rng
    )

    # 3. Simulate wildfire intensity
    intensity = simulate_intensity(
        times=times,
        wind=wind,
        humidity=humidity,
        I0=I0,
        alpha_c=ALPHA_C,
        beta_w=BETA_W,
        beta_h=BETA_H,
        gamma=gamma,
        sigma_i=SIGMA_I,
        rng=intensity_rng
    )

    return times, wind, humidity, intensity
