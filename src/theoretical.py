import numpy as np

from src.config import (
    KAPPA_W, SIGMA_W,
    KAPPA_H, SIGMA_H,
    ALPHA_C, BETA_W, BETA_H,
    SIGMA_I, I0, T
)


def integrated_ou_moments(z0, kappa, sigma, T):
    """
    Calculate the theoretical mean and variance of

        J_T = integral from 0 to T of Z_t dt

    where Z_t follows an Ornstein-Uhlenbeck process.

    Derived from the analytical OU solution in
    Section 3.1 of the mathematical model.
    """

    if kappa <= 0:
        raise ValueError("kappa must be positive")

    if sigma < 0:
        raise ValueError("sigma cannot be negative")

    if T < 0:
        raise ValueError("T cannot be negative")

    # Expected value of the integrated OU process
    mean = (
        z0 / kappa
        * (1 - np.exp(-kappa * T))
    )

    # Variance of the integrated OU process
    variance = (
        sigma**2 / kappa**2
        * (
            T
            - 2 * (1 - np.exp(-kappa * T)) / kappa
            + (1 - np.exp(-2 * kappa * T)) / (2 * kappa)
        )
    )

    # Protect against small negative rounding errors
    variance = max(0.0, variance)

    return mean, variance


def theoretical_log_intensity(
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T
):
    """
    Calculate the theoretical mean and variance
    of log wildfire intensity when gamma = 0.

    Assumptions:
    - Fixed initial conditions.
    - Independent Brownian motions.
    - Additive environmental growth rate.
    """

    # Integrated wind OU moments
    wind_mean, wind_var = integrated_ou_moments(
        z0=wind_z0,
        kappa=KAPPA_W,
        sigma=SIGMA_W,
        T=horizon
    )

    # Integrated humidity OU moments
    humidity_mean, humidity_var = integrated_ou_moments(
        z0=humidity_z0,
        kappa=KAPPA_H,
        sigma=SIGMA_H,
        T=horizon
    )

    # Mean of log intensity
    log_mean = (
        np.log(I0)
        + (ALPHA_C - 0.5 * SIGMA_I**2) * horizon
        + BETA_W * wind_mean
        - BETA_H * humidity_mean
    )

    # Variance of log intensity
    log_variance = (
        BETA_W**2 * wind_var
        + BETA_H**2 * humidity_var
        + SIGMA_I**2 * horizon
    )

    return log_mean, log_variance


def theoretical_intensity_moments(
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T
):
    """
    Calculate theoretical wildfire intensity statistics
    under the additive model (gamma = 0).
    """

    log_mean, log_variance = theoretical_log_intensity(
        wind_z0=wind_z0,
        humidity_z0=humidity_z0,
        horizon=horizon
    )

    # Lognormal expectation
    mean = np.exp(log_mean + 0.5 * log_variance)

    # Lognormal variance
    variance = (
        np.exp(log_variance) - 1
    ) * np.exp(2 * log_mean + log_variance)

    # Standard deviation
    std = np.sqrt(variance)

    # Lognormal median
    median = np.exp(log_mean)

    return {
        "log_mean": log_mean,
        "log_variance": log_variance,
        "mean": mean,
        "variance": variance,
        "std": std,
        "median": median
    }

def theoretical_interaction_log_effect(
    gamma,
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T
):
    """
    Calculate the theoretical expected difference
    in log-intensity between the interaction
    and additive models.

    E[ln(I_gamma) - ln(I_additive)]

    Assumptions:
    - Independent environmental Brownian motions.
    - Deterministic initial environmental states.
    - Identical parameters except gamma.
    """

    combined_kappa = KAPPA_W + KAPPA_H

    effect = (
        gamma
        * wind_z0
        * humidity_z0
        / combined_kappa
        * (1 - np.exp(-combined_kappa * horizon))
    )

    return effect

def theoretical_log_sensitivity(
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T,
    kappa_w=KAPPA_W,
    kappa_h=KAPPA_H,
    sigma_i=SIGMA_I
):
    """
    Calculate analytical derivatives of expected
    terminal log-intensity with respect to the
    intensity model parameters.

    Assumptions:
    - Independent Brownian drivers.
    - Deterministic initial environmental conditions.
    - Constant model coefficients.
    """

    if horizon < 0:
        raise ValueError("Horizon cannot be negative")

    if kappa_w <= 0 or kappa_h <= 0:
        raise ValueError(
            "Mean-reversion rates must be positive"
        )

    if sigma_i < 0:
        raise ValueError(
            "Volatility cannot be negative"
        )

    # --------------------------------------------------
    # Sensitivity to baseline growth coefficient
    # --------------------------------------------------

    sensitivity_alpha = horizon

    # --------------------------------------------------
    # Sensitivity to wind coefficient
    # --------------------------------------------------

    sensitivity_beta_w = (
        wind_z0 / kappa_w
        * (1 - np.exp(-kappa_w * horizon))
    )

    # --------------------------------------------------
    # Sensitivity to humidity coefficient
    # --------------------------------------------------

    sensitivity_beta_h = (
        -humidity_z0 / kappa_h
        * (1 - np.exp(-kappa_h * horizon))
    )

    # --------------------------------------------------
    # Sensitivity to wind-humidity interaction
    # --------------------------------------------------

    combined_kappa = kappa_w + kappa_h

    sensitivity_gamma = (
        wind_z0 * humidity_z0 / combined_kappa
        * (1 - np.exp(-combined_kappa * horizon))
    )

    # --------------------------------------------------
    # Sensitivity to intensity volatility
    # --------------------------------------------------

    sensitivity_sigma_i = -sigma_i * horizon

    return {
        "alpha_c": sensitivity_alpha,
        "beta_w": sensitivity_beta_w,
        "beta_h": sensitivity_beta_h,
        "gamma": sensitivity_gamma,
        "sigma_i": sensitivity_sigma_i
    }

if __name__ == "__main__":

    results = theoretical_intensity_moments()

    print("Theoretical Additive Model (gamma = 0)")
    print("-" * 45)

    for name, value in results.items():
        print(f"{name:15s}: {value:.5f}")
        
    print("\nANALYTICAL LOG-INTENSITY SENSITIVITY")
    print("=" * 55)

    sensitivities = theoretical_log_sensitivity()

    for parameter, value in sensitivities.items():

        print(f"{parameter:15s}: {value:12.5f}")
        
        
# ==========================================================
# THEORETICAL RESULTS — ADDITIVE MODEL (gamma = 0)
# ==========================================================
#
# Initial intensity: 500 kW/m
# Initial wind anomaly: 1.5
# Initial humidity anomaly: -1.0
# Simulation horizon: T = 10
#
# Log-intensity distribution:
# ln(I_T) ~ Normal(mean=7.02569, variance=0.54556)
#
# Intensity distribution:
# I_T ~ Lognormal(meanlog=7.02569, varlog=0.54556)
#
# Theoretical statistics:
# E[I_T]   = 1478.02852 kW/m
# Var(I_T) = 1585062.13318 (kW/m)^2
# SD(I_T)  = 1258.99251 kW/m
# Median   = 1125.16507 kW/m
#
# Observations:
# - The mean exceeds the median because the lognormal
#   distribution is positively skewed.
# - The large standard deviation indicates substantial
#   uncertainty in future wildfire intensity.
# - These results provide an analytical benchmark for
#   validating the Monte Carlo simulation.
# - These calculations assume gamma = 0, fixed initial
#   conditions, and independent Brownian motions.
# - Model parameters are illustrative, not calibrated.
#
# Status: Theoretical calculations implemented.
# Next: Monte Carlo simulation and numerical validation.
# ==========================================================