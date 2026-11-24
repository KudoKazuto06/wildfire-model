import numpy as np

from src.intensity import growth_rate

from src.config import (
    KAPPA_W,
    SIGMA_W,
    KAPPA_H,
    SIGMA_H,
    ALPHA_C,
    BETA_W,
    BETA_H,
    GAMMA,
    SIGMA_I,
    I0,
    T,
    DT,
    SEED,
    N_PATHS,
)

from src.dependence import (
    validate_correlation,
    ou_innovation_correlation,
)


def simulate_many(
    n_paths=N_PATHS,
    wind_z0=1.5,
    humidity_z0=-1.0,
    gamma=GAMMA,
    horizon=T,
    dt=DT,
    seed=SEED,
    *,
    kappa_w=KAPPA_W,
    sigma_w=SIGMA_W,
    kappa_h=KAPPA_H,
    sigma_h=SIGMA_H,
    alpha_c=ALPHA_C,
    beta_w=BETA_W,
    beta_h=BETA_H,
    sigma_i=SIGMA_I,
    initial_intensity=I0,
    rho=0.0,
):
    """
    Simulate multiple independent wildfire trajectories
    using exact OU transitions and log-Euler intensity
    updates.

    All dynamic model parameters can be overridden
    for sensitivity analysis.

    Environmental Brownian motions may be correlated
    through the parameter rho.

    rho = 0:
        Independent environmental Brownian motions.

    rho < 0:
        Negatively correlated environmental Brownian motions.

    rho > 0:
        Positively correlated environmental Brownian motions.

    The intensity Brownian motion remains independent
    of both environmental Brownian motions.

    Default arguments reproduce the original baseline
    simulation used in Q2 and Q3.

    Returns
    -------
    times : ndarray
        Numerical observation times.

    intensity_paths : ndarray
        Intensity trajectory matrix with shape
        (n_paths, n_steps + 1).
    """

    # --------------------------------------------------
    # 1. Validate simulation inputs
    # --------------------------------------------------

    if not isinstance(n_paths, (int, np.integer)):
        raise ValueError("n_paths must be an integer")

    if n_paths <= 0:
        raise ValueError("n_paths must be positive")

    parameters = [
        wind_z0,
        humidity_z0,
        gamma,
        horizon,
        dt,
        kappa_w,
        sigma_w,
        kappa_h,
        sigma_h,
        alpha_c,
        beta_w,
        beta_h,
        sigma_i,
        initial_intensity,
    ]

    if not np.all(np.isfinite(parameters)):
        raise ValueError("Model parameters must be finite")

    if horizon <= 0 or dt <= 0:
        raise ValueError(
            "horizon and dt must be positive"
        )

    if initial_intensity <= 0:
        raise ValueError(
            "Initial intensity must be positive"
        )

    if kappa_w <= 0 or kappa_h <= 0:
        raise ValueError(
            "Mean-reversion rates must be positive"
        )

    if sigma_w < 0 or sigma_h < 0 or sigma_i < 0:
        raise ValueError(
            "Volatility parameters cannot be negative"
        )

    if beta_w < 0 or beta_h < 0:
        raise ValueError(
            "Environmental coefficients cannot be negative"
        )

    # Q6: Validate Brownian correlation.
    validate_correlation(rho)

    # --------------------------------------------------
    # 2. Construct numerical time grid
    # --------------------------------------------------

    n_steps = int(np.ceil(horizon / dt))

    # Adjust dt so the final observation is exactly
    # at the specified simulation horizon.
    dt = horizon / n_steps

    times = np.linspace(
        0,
        horizon,
        n_steps + 1,
    )

    # --------------------------------------------------
    # 3. Initialize independent random streams
    # --------------------------------------------------

    streams = np.random.SeedSequence(seed).spawn(3)

    wind_rng = np.random.default_rng(streams[0])

    humidity_rng = np.random.default_rng(streams[1])

    intensity_rng = np.random.default_rng(streams[2])

    # The underlying random streams remain independent.
    # Environmental dependence is introduced later
    # through a linear transformation of the innovations.

    # --------------------------------------------------
    # 4. Initialize stochastic state variables
    # --------------------------------------------------

    wind = np.full(
        n_paths,
        wind_z0,
        dtype=float,
    )

    humidity = np.full(
        n_paths,
        humidity_z0,
        dtype=float,
    )

    log_intensity = np.full(
        n_paths,
        np.log(initial_intensity),
        dtype=float,
    )

    # Store complete intensity trajectories.
    intensity_paths = np.empty(
        (n_paths, n_steps + 1)
    )

    intensity_paths[:, 0] = initial_intensity

    # --------------------------------------------------
    # 5. Precompute exact OU transition coefficients
    # --------------------------------------------------

    # Exact wind mean-reversion coefficient.
    wind_decay = np.exp(-kappa_w * dt)

    # Exact wind innovation standard deviation.
    wind_noise_std = sigma_w * np.sqrt(
        -np.expm1(-2 * kappa_w * dt)
        / (2 * kappa_w)
    )

    # Exact humidity mean-reversion coefficient.
    humidity_decay = np.exp(-kappa_h * dt)

    # Exact humidity innovation standard deviation.
    humidity_noise_std = sigma_h * np.sqrt(
        -np.expm1(-2 * kappa_h * dt)
        / (2 * kappa_h)
    )

    # --------------------------------------------------
    # 5.1 Q6: Exact environmental innovation correlation
    # --------------------------------------------------

    # Brownian correlation rho is not necessarily equal
    # to the correlation between exact OU innovations.
    #
    # Calculate the correct innovation correlation
    # using the exact OU transition covariance matrix.

    innovation_corr = ou_innovation_correlation(
        dt=dt,
        rho=rho,
        kappa_w=kappa_w,
        sigma_w=sigma_w,
        kappa_h=kappa_h,
        sigma_h=sigma_h,
    )

    # Weight of the independent humidity innovation.
    #
    # epsilon_h_correlated =
    #     innovation_corr * epsilon_w
    #     + independent_weight * epsilon_h
    #
    # This transformation preserves unit variance.

    independent_weight = np.sqrt(
        max(0.0, 1.0 - innovation_corr**2)
    )

    # --------------------------------------------------
    # 5.2 Intensity transition coefficients
    # --------------------------------------------------

    # Intensity diffusion coefficient.
    intensity_noise_std = sigma_i * np.sqrt(dt)

    # Ito correction from transforming I_t to log(I_t).
    ito_correction = 0.5 * sigma_i**2

    # --------------------------------------------------
    # 6. Simulate all trajectories simultaneously
    # --------------------------------------------------

    for n in range(n_steps):

        # Generate independent standard Gaussian
        # innovations for all trajectories.

        epsilon_w = wind_rng.normal(size=n_paths)

        epsilon_h = humidity_rng.normal(size=n_paths)

        epsilon_i = intensity_rng.normal(size=n_paths)

        # ----------------------------------------------
        # 6.1 Q6: Introduce environmental dependence
        # ----------------------------------------------

        # Transform the humidity innovation so that
        # its correlation with the wind innovation
        # matches the exact OU transition correlation.
        #
        # When rho = 0, no transformation is applied.
        # This preserves the original independent
        # simulation and its numerical reproducibility.

        if innovation_corr != 0.0:

            epsilon_h = (
                innovation_corr * epsilon_w
                + independent_weight * epsilon_h
            )

        # ----------------------------------------------
        # 6.2 Evaluate the current intensity growth rate
        # ----------------------------------------------

        # Evaluate growth using the environmental
        # conditions at the beginning of the time step.

        r = growth_rate(
            wind,
            humidity,
            alpha_c,
            beta_w,
            beta_h,
            gamma,
        )

        # ----------------------------------------------
        # 6.3 Log-Euler intensity update
        # ----------------------------------------------

        # Update log-intensity rather than intensity
        # directly to preserve positivity.
        #
        # The intensity innovation remains independent
        # of both environmental innovations.

        log_intensity += (
            (r - ito_correction) * dt
            + intensity_noise_std * epsilon_i
        )

        # ----------------------------------------------
        # 6.4 Exact wind OU transition
        # ----------------------------------------------

        wind = (
            wind_decay * wind
            + wind_noise_std * epsilon_w
        )

        # ----------------------------------------------
        # 6.5 Exact humidity OU transition
        # ----------------------------------------------

        # The humidity innovation now incorporates
        # environmental dependence when rho != 0.

        humidity = (
            humidity_decay * humidity
            + humidity_noise_std * epsilon_h
        )

        # ----------------------------------------------
        # 6.6 Recover physical wildfire intensity
        # ----------------------------------------------

        intensity_paths[:, n + 1] = np.exp(
            log_intensity
        )

    # --------------------------------------------------
    # 7. Validate numerical outputs
    # --------------------------------------------------

    if not np.all(np.isfinite(intensity_paths)):
        raise FloatingPointError(
            "Non-finite intensity values detected"
        )

    if np.any(intensity_paths <= 0):
        raise FloatingPointError(
            "Non-positive intensity values detected"
        )

    return times, intensity_paths


if __name__ == "__main__":

    # --------------------------------------------------
    # Baseline regression check
    # --------------------------------------------------

    print("BASELINE ADDITIVE MODEL")
    print("=" * 55)

    times, paths = simulate_many(
        gamma=0.0,
        rho=0.0,
    )

    final_intensity = paths[:, -1]

    print(f"Number of paths: {paths.shape[0]}")
    print(f"Number of time points: {paths.shape[1]}")
    print(f"Mean: {np.mean(final_intensity):.5f}")
    print(f"Variance: {np.var(final_intensity, ddof=1):.5f}")
    print(f"Median: {np.median(final_intensity):.5f}")

    del final_intensity, paths

    # --------------------------------------------------
    # Preliminary parameter sensitivity check
    # --------------------------------------------------

    print("\nWIND COEFFICIENT SENSITIVITY CHECK")
    print("=" * 55)

    print(
        f"{'Beta W':>10}"
        f"{'Mean':>15}"
        f"{'Std':>15}"
        f"{'P(Max >= 2000)':>18}"
    )

    for beta_w_value in [0.06, 0.12, 0.18]:

        times, paths = simulate_many(
            beta_w=beta_w_value,
            seed=SEED,
            rho=0.0,
        )

        final_intensity = paths[:, -1]

        mean = np.mean(final_intensity)

        std = np.std(
            final_intensity,
            ddof=1,
        )

        crossing_probability = np.mean(
            np.max(paths, axis=1) >= 2000
        )

        print(
            f"{beta_w_value:10.2f}"
            f"{mean:15.5f}"
            f"{std:15.5f}"
            f"{crossing_probability:18.5f}"
        )

        del final_intensity, paths

    # --------------------------------------------------
    # Q6: Preliminary dependence sensitivity check
    # --------------------------------------------------

    print("\nENVIRONMENTAL DEPENDENCE CHECK")
    print("=" * 85)

    print(
        f"{'Rho':>10}"
        f"{'Mean':>15}"
        f"{'Std':>15}"
        f"{'Log Mean':>15}"
        f"{'P(Max >= 2000)':>18}"
    )

    for rho_value in [-1.0, -0.5, 0.0, 0.5, 1.0]:

        times, paths = simulate_many(
            rho=rho_value,
            gamma=GAMMA,
            seed=SEED,
        )

        final_intensity = paths[:, -1]

        mean = np.mean(final_intensity)

        std = np.std(
            final_intensity,
            ddof=1,
        )

        log_mean = np.mean(
            np.log(final_intensity)
        )

        crossing_probability = np.mean(
            np.max(paths, axis=1) >= 2000
        )

        print(
            f"{rho_value:10.2f}"
            f"{mean:15.5f}"
            f"{std:15.5f}"
            f"{log_mean:15.5f}"
            f"{crossing_probability:18.5f}"
        )

        del final_intensity, paths


# ==========================================================
# MONTE CARLO VALIDATION — ADDITIVE MODEL (gamma = 0)
# ==========================================================
#
# Simulation settings:
# Number of trajectories: 10000
# Number of time steps: 1000
# Number of time points: 1001
# Intensity matrix shape: (10000, 1001)
# Execution time: 1.074 seconds
# Environmental Brownian correlation: rho = 0
#
# Statistic       Theoretical       Monte Carlo
# ----------------------------------------------------------
# Log mean        7.02569           7.03215
# Log variance    0.54556           0.55024
# Mean            1478.02852        1485.38922
# Variance        1585062.13318     1495024.21055
# Std             1258.99251        1222.71183
# Median          1125.16507        1132.65977
#
# Observations:
# - The empirical mean differs from theory by approximately
#   0.50%.
# - The empirical median differs from theory by approximately
#   0.67%.
# - The empirical log variance differs by approximately 0.86%.
# - The larger intensity variance error is consistent with
#   the sensitivity of variance estimation to extreme outcomes.
#
# Monte Carlo standard error of the mean:
# SE = 1258.99251 / sqrt(10000)
#    = 12.58993
#
# Standardized difference:
# z = (1485.38922 - 1478.02852) / 12.58993
#   = 0.585
#
# Conclusion:
# The empirical statistics agree closely with the theoretical
# additive model. This supports the correctness of the
# vectorized Monte Carlo implementation.
#
# Remaining differences may arise from Monte Carlo sampling
# variation and discretization of the intensity drift integral.
#
# This validates the numerical implementation, not the
# predictive accuracy of the model against real wildfire data.
# ==========================================================


# ==========================================================
# Q6 — EXTENSION TO CORRELATED ENVIRONMENTAL PROCESSES
# ==========================================================
#
# The original model assumes independent environmental
# Brownian motions:
#
#     d<B_W, B_H>_t = 0
#
# The extended model allows:
#
#     d<B_W, B_H>_t = rho * dt
#
# where -1 <= rho <= 1.
#
# Exact OU transition covariance:
#
# q_WH =
#     rho * sigma_w * sigma_h
#     * (1 - exp(-(kappa_w + kappa_h) * dt))
#     / (kappa_w + kappa_h)
#
# Exact innovation correlation:
#
#     r_dt = q_WH / sqrt(q_W * q_H)
#
# Correlated Gaussian construction:
#
#     epsilon_W = U
#
#     epsilon_H =
#         r_dt * U + sqrt(1 - r_dt**2) * V
#
# where U and V are independent standard normal variables.
#
# Important properties:
#
# - Individual OU marginal distributions are preserved.
# - Environmental innovations have the required covariance.
# - Intensity Brownian motion remains independent.
# - Different Monte Carlo trajectories remain independent.
# - Intensity remains positive through the log-Euler method.
# - rho = 0 preserves the original simulation algorithm.
# - Existing Q1-Q5 experiments require no modifications.
# - Common random numbers support paired comparisons
#   across different correlation scenarios.
#
# The historical correlation between observed wind and
# humidity levels is not treated as a calibrated estimate
# of the underlying Brownian correlation.
#
# Instead, rho is investigated through scenario analysis.
#
# Exact environmental transitions do not make the entire
# wildfire intensity simulation exact. The intensity drift
# integral is still approximated using the log-Euler method.
# ==========================================================