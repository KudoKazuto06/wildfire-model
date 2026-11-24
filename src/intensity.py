
import numpy as np


def growth_rate(wind, humidity, alpha_c,
                beta_w, beta_h, gamma):
    """
    Calculate wildfire intensity growth rate.

    r(w, h) = alpha_c + beta_w*w
              - beta_h*h + gamma*w*h
    """

    return (
        alpha_c
        + beta_w * wind
        - beta_h * humidity
        + gamma * wind * humidity
    )


def simulate_intensity(
    times,
    wind,
    humidity,
    I0,
    alpha_c,
    beta_w,
    beta_h,
    gamma,
    sigma_i,
    rng
):
    """
    Simulate wildfire intensity using log-Euler.

    Wind and humidity are standardized OU trajectories.
    """

    times = np.asarray(times)
    wind = np.asarray(wind)
    humidity = np.asarray(humidity)

    if I0 <= 0:
        raise ValueError("Initial intensity must be positive")

    if sigma_i < 0:
        raise ValueError("Intensity volatility cannot be negative")

    if len(times) != len(wind) or len(times) != len(humidity):
        raise ValueError("Input arrays must have equal lengths")

    if len(times) < 2 or not np.all(np.diff(times) > 0):
        raise ValueError("Times must be strictly increasing")

    if not (
        np.all(np.isfinite(times))
        and np.all(np.isfinite(wind))
        and np.all(np.isfinite(humidity))
    ):
        raise ValueError("Inputs must be finite")

    log_intensity = np.empty(len(times))
    log_intensity[0] = np.log(I0)

    for n in range(len(times) - 1):
        dt = times[n + 1] - times[n]

        r = growth_rate(
            wind[n],
            humidity[n],
            alpha_c,
            beta_w,
            beta_h,
            gamma
        )

        epsilon = rng.normal()

        log_intensity[n + 1] = (
            log_intensity[n]
            + (r - 0.5 * sigma_i**2) * dt
            + sigma_i * np.sqrt(dt) * epsilon
        )

    intensity = np.exp(log_intensity)

    if not np.all(np.isfinite(intensity)) or np.any(intensity <= 0):
        raise FloatingPointError(
            "Intensity overflowed or underflowed. "
            "Check parameters and simulation horizon."
        )

    return intensity
