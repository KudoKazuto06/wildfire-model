import numpy as np

from scipy.stats import norm

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
)

from src.theoretical import integrated_ou_moments


def validate_correlation(rho):
    """
    Validate the Brownian correlation parameter.
    """
    if not np.isfinite(rho) or not -1 <= rho <= 1:
        raise ValueError(
            "Correlation must be between -1 and 1"
        )


def environmental_covariance(
    t,
    rho,
    kappa_w=KAPPA_W,
    sigma_w=SIGMA_W,
    kappa_h=KAPPA_H,
    sigma_h=SIGMA_H,
):
    """
    Calculate Cov(Z_t^W, Z_t^H) analytically.

    Assumes deterministic initial environmental conditions
    and correlated environmental Brownian motions.
    """

    validate_correlation(rho)

    t = np.asarray(t, dtype=float)

    if np.any(t < 0) or not np.all(np.isfinite(t)):
        raise ValueError("Time must be finite and non-negative")

    if kappa_w <= 0 or kappa_h <= 0:
        raise ValueError("Mean reversion must be positive")

    if sigma_w < 0 or sigma_h < 0:
        raise ValueError("Volatility cannot be negative")

    kappa_sum = kappa_w + kappa_h

    covariance = (
        rho
        * sigma_w
        * sigma_h
        * (-np.expm1(-kappa_sum * t))
        / kappa_sum
    )

    return covariance


def stationary_covariance(
    rho,
    kappa_w=KAPPA_W,
    sigma_w=SIGMA_W,
    kappa_h=KAPPA_H,
    sigma_h=SIGMA_H,
):
    """
    Calculate the limiting environmental covariance
    as time approaches infinity.
    """

    validate_correlation(rho)

    if kappa_w <= 0 or kappa_h <= 0:
        raise ValueError("Mean reversion must be positive")

    if sigma_w < 0 or sigma_h < 0:
        raise ValueError("Volatility cannot be negative")

    return (
        rho * sigma_w * sigma_h
        / (kappa_w + kappa_h)
    )


def stationary_correlation(
    rho,
    kappa_w=KAPPA_W,
    kappa_h=KAPPA_H,
):
    """
    Calculate the stationary correlation between the
    environmental OU processes.

    Assumes both environmental volatilities are positive.
    """

    validate_correlation(rho)

    if kappa_w <= 0 or kappa_h <= 0:
        raise ValueError("Mean reversion must be positive")

    return (
        2 * rho * np.sqrt(kappa_w * kappa_h)
        / (kappa_w + kappa_h)
    )


def ou_transition_covariance(
    dt,
    rho,
    kappa_w=KAPPA_W,
    sigma_w=SIGMA_W,
    kappa_h=KAPPA_H,
    sigma_h=SIGMA_H,
):
    """
    Calculate the covariance matrix of the exact
    OU transition innovations over one time step.
    """

    validate_correlation(rho)

    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Time step must be positive")

    if kappa_w <= 0 or kappa_h <= 0:
        raise ValueError("Mean reversion must be positive")

    if sigma_w < 0 or sigma_h < 0:
        raise ValueError("Volatility cannot be negative")

    wind_variance = (
        sigma_w**2
        * (-np.expm1(-2 * kappa_w * dt))
        / (2 * kappa_w)
    )

    humidity_variance = (
        sigma_h**2
        * (-np.expm1(-2 * kappa_h * dt))
        / (2 * kappa_h)
    )

    covariance = environmental_covariance(
        t=dt,
        rho=rho,
        kappa_w=kappa_w,
        sigma_w=sigma_w,
        kappa_h=kappa_h,
        sigma_h=sigma_h,
    )

    return np.array([
        [wind_variance, covariance],
        [covariance, humidity_variance],
    ])


def integrated_decay(kappa, horizon):
    """
    Calculate the integral of exponential decay:

        A(kappa, T) = (1 - exp(-kappa * T)) / kappa
    """

    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("Mean reversion must be positive")

    if not np.isfinite(horizon) or horizon < 0:
        raise ValueError("Horizon must be non-negative")

    return -np.expm1(-kappa * horizon) / kappa


def integrated_ou_covariance(
    rho,
    horizon=T,
    kappa_w=KAPPA_W,
    sigma_w=SIGMA_W,
    kappa_h=KAPPA_H,
    sigma_h=SIGMA_H,
):
    """
    Calculate the covariance between the integrated
    wind and humidity OU processes:

        Cov(J_W, J_H)

    where:

        J_W = integral of Z_W over [0, T]
        J_H = integral of Z_H over [0, T]
    """

    validate_correlation(rho)

    if (
        not np.isfinite(sigma_w)
        or not np.isfinite(sigma_h)
        or sigma_w < 0
        or sigma_h < 0
    ):
        raise ValueError("Volatility must be non-negative")

    # Calculate integrated exponential decay terms.
    a_w = integrated_decay(kappa_w, horizon)

    a_h = integrated_decay(kappa_h, horizon)

    a_wh = integrated_decay(
        kappa_w + kappa_h,
        horizon,
    )

    # Evaluate the analytical covariance formula.
    covariance = (
        rho
        * sigma_w
        * sigma_h
        / (kappa_w * kappa_h)
        * (
            horizon
            - a_w
            - a_h
            + a_wh
        )
    )

    return covariance


def correlated_additive_log_moments(
    rho,
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T,
):
    """
    Calculate the analytical mean and variance of
    terminal log-intensity under the correlated
    additive model (gamma = 0).
    """

    validate_correlation(rho)

    # Individual integrated OU moments.
    wind_mean, wind_variance = integrated_ou_moments(
        wind_z0,
        KAPPA_W,
        SIGMA_W,
        horizon,
    )

    humidity_mean, humidity_variance = (
        integrated_ou_moments(
            humidity_z0,
            KAPPA_H,
            SIGMA_H,
            horizon,
        )
    )

    # New covariance term introduced by dependence.
    covariance = integrated_ou_covariance(
        rho=rho,
        horizon=horizon,
    )

    # Correlation does not affect the additive log mean.
    log_mean = (
        np.log(I0)
        + (ALPHA_C - 0.5 * SIGMA_I**2) * horizon
        + BETA_W * wind_mean
        - BETA_H * humidity_mean
    )

    # Correlation changes the additive log variance.
    log_variance = (
        BETA_W**2 * wind_variance
        + BETA_H**2 * humidity_variance
        - 2 * BETA_W * BETA_H * covariance
        + SIGMA_I**2 * horizon
    )

    # Protect against small negative rounding errors.
    if log_variance < -1e-12:
        raise FloatingPointError(
            "Negative theoretical log variance"
        )

    log_variance = max(0.0, log_variance)

    return log_mean, log_variance


def correlated_additive_statistics(
    rho,
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T,
):
    """
    Calculate analytical terminal intensity statistics
    for the correlated additive model.

    Terminal intensity remains lognormally distributed.
    """

    log_mean, log_variance = (
        correlated_additive_log_moments(
            rho=rho,
            wind_z0=wind_z0,
            humidity_z0=humidity_z0,
            horizon=horizon,
        )
    )

    log_std = np.sqrt(log_variance)

    # Lognormal arithmetic mean.
    mean = np.exp(
        log_mean + 0.5 * log_variance
    )

    # Lognormal variance.
    variance = (
        np.expm1(log_variance)
        * np.exp(2 * log_mean + log_variance)
    )

    # Lognormal median.
    median = np.exp(log_mean)

    # Central 95% distribution quantiles.
    q025 = np.exp(
        log_mean + log_std * norm.ppf(0.025)
    )

    q975 = np.exp(
        log_mean + log_std * norm.ppf(0.975)
    )

    return {
        "rho": rho,
        "log_mean": log_mean,
        "log_variance": log_variance,
        "mean": mean,
        "variance": variance,
        "std": np.sqrt(variance),
        "median": median,
        "q025": q025,
        "q975": q975,
    }


def correlated_interaction_log_effect(
    rho,
    gamma=GAMMA,
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T,
):
    """
    Calculate the expected log-intensity contribution
    of the environmental interaction term.

    Separates the contribution into:

    1. Initial environmental conditions.
    2. Environmental Brownian dependence.
    """

    validate_correlation(rho)

    kappa_sum = KAPPA_W + KAPPA_H

    a_wh = integrated_decay(
        kappa_sum,
        horizon,
    )

    # Interaction contribution from deterministic
    # initial environmental conditions.
    initial_effect = (
        gamma
        * wind_z0
        * humidity_z0
        * a_wh
    )

    # Additional interaction contribution caused
    # by environmental dependence.
    covariance_effect = (
        gamma
        * rho
        * SIGMA_W
        * SIGMA_H
        / kappa_sum
        * (horizon - a_wh)
    )

    total_effect = (
        initial_effect + covariance_effect
    )

    return {
        "initial_effect": initial_effect,
        "covariance_effect": covariance_effect,
        "total_effect": total_effect,
    }


def correlated_interaction_log_mean(
    rho,
    gamma=GAMMA,
    wind_z0=1.5,
    humidity_z0=-1.0,
    horizon=T,
):
    """
    Calculate the exact expected terminal log-intensity
    for the correlated interaction model.

    This does not assume that the interaction model
    has a lognormal terminal distribution.
    """

    additive_log_mean, _ = (
        correlated_additive_log_moments(
            rho=rho,
            wind_z0=wind_z0,
            humidity_z0=humidity_z0,
            horizon=horizon,
        )
    )

    interaction_effect = (
        correlated_interaction_log_effect(
            rho=rho,
            gamma=gamma,
            wind_z0=wind_z0,
            humidity_z0=humidity_z0,
            horizon=horizon,
        )
    )

    return (
        additive_log_mean
        + interaction_effect["total_effect"]
    )
    
def analyze_correlated_intensity():
    """
    Investigate analytical terminal intensity statistics
    across different environmental correlations.
    """

    correlations = [
        -1.0,
        -0.75,
        -0.50,
        -0.25,
        0.0,
        0.25,
        0.50,
        0.75,
        1.0,
    ]

    print("\nINTEGRATED ENVIRONMENTAL COVARIANCE")
    print("=" * 85)

    print(
        f"{'Rho':>8}"
        f"{'Cov(J_W, J_H)':>22}"
    )

    for rho in correlations:

        covariance = integrated_ou_covariance(rho)

        print(
            f"{rho:8.2f}"
            f"{covariance:22.5f}"
        )

    print("\nCORRELATED ADDITIVE MODEL")
    print("=" * 110)

    print(
        f"{'Rho':>8}"
        f"{'Log Mean':>15}"
        f"{'Log Variance':>15}"
        f"{'Mean':>15}"
        f"{'Std':>15}"
        f"{'Median':>15}"
        f"{'Q97.5':>15}"
    )

    results = []

    for rho in correlations:

        stats = correlated_additive_statistics(rho)

        results.append(stats)

        print(
            f"{rho:8.2f}"
            f"{stats['log_mean']:15.5f}"
            f"{stats['log_variance']:15.5f}"
            f"{stats['mean']:15.5f}"
            f"{stats['std']:15.5f}"
            f"{stats['median']:15.5f}"
            f"{stats['q975']:15.5f}"
        )

    print("\nCORRELATED INTERACTION MODEL")
    print("=" * 110)

    print(
        f"{'Rho':>8}"
        f"{'Initial Effect':>18}"
        f"{'Covariance Effect':>20}"
        f"{'Total Effect':>18}"
        f"{'Log Mean':>15}"
        f"{'Geometric Mean':>18}"
    )

    interaction_log_means = []

    for rho in correlations:

        effect = correlated_interaction_log_effect(rho)

        log_mean = correlated_interaction_log_mean(rho)

        interaction_log_means.append(log_mean)

        geometric_mean = np.exp(log_mean)

        print(
            f"{rho:8.2f}"
            f"{effect['initial_effect']:18.5f}"
            f"{effect['covariance_effect']:20.5f}"
            f"{effect['total_effect']:18.5f}"
            f"{log_mean:15.5f}"
            f"{geometric_mean:18.5f}"
        )

    print("\nINTENSITY CONSISTENCY CHECKS")
    print("=" * 85)

    # Check 1:
    # Independence must produce zero integrated covariance.
    assert np.isclose(
        integrated_ou_covariance(0),
        0,
    )

    # Check 2:
    # Zero horizon must produce zero integrated covariance.
    assert np.isclose(
        integrated_ou_covariance(0.5, horizon=0),
        0,
    )

    # Check 3:
    # Reversing correlation reverses covariance.
    assert np.isclose(
        integrated_ou_covariance(-0.5),
        -integrated_ou_covariance(0.5),
    )

    # Check 4:
    # Recover the original independent additive model.
    from src.theoretical import (
        theoretical_intensity_moments,
    )

    independent = theoretical_intensity_moments()

    correlated_zero = correlated_additive_statistics(0)

    for key in [
        "log_mean",
        "log_variance",
        "mean",
        "variance",
        "std",
        "median",
    ]:

        assert np.isclose(
            independent[key],
            correlated_zero[key],
        )

    # Check 5:
    # Additive log means remain constant.
    log_means = [
        result["log_mean"]
        for result in results
    ]

    assert np.allclose(log_means, log_means[0])

    # Check 6:
    # Additive medians remain constant.
    medians = [
        result["median"]
        for result in results
    ]

    assert np.allclose(medians, medians[0])

    # Check 7:
    # With positive wind and humidity coefficients,
    # log variance decreases as correlation increases.
    log_variances = [
        result["log_variance"]
        for result in results
    ]

    assert np.all(
        np.diff(log_variances) <= 1e-12
    )

    # Check 8:
    # Arithmetic mean decreases with correlation.
    means = [
        result["mean"]
        for result in results
    ]

    assert np.all(
        np.diff(means) <= 1e-12
    )

    # Check 9:
    # Independence eliminates the covariance contribution
    # to the expected interaction effect.
    zero_effect = correlated_interaction_log_effect(0)

    assert np.isclose(
        zero_effect["covariance_effect"],
        0,
    )

    # Check 10:
    # Setting gamma to zero eliminates all interaction.
    for rho in correlations:

        effect = correlated_interaction_log_effect(
            rho=rho,
            gamma=0,
        )

        assert np.isclose(
            effect["total_effect"],
            0,
        )

    # Check 11:
    # At T = 0, intensity must equal its initial value.
    for rho in correlations:

        stats = correlated_additive_statistics(
            rho=rho,
            horizon=0,
        )

        assert np.isclose(stats["mean"], I0)

        assert np.isclose(stats["variance"], 0)

        assert np.isclose(stats["median"], I0)

    # Check 12:
    # With the selected negative gamma, expected
    # interaction log-intensity decreases with rho.
    assert np.all(
        np.diff(interaction_log_means) <= 1e-12
    )

    print("All intensity consistency checks passed.")


def ou_innovation_correlation(
    dt,
    rho,
    kappa_w=KAPPA_W,
    sigma_w=SIGMA_W,
    kappa_h=KAPPA_H,
    sigma_h=SIGMA_H,
):
    """
    Calculate the correlation between the exact
    wind and humidity OU transition innovations.

    This is generally different from the underlying
    Brownian correlation rho.
    """

    covariance_matrix = ou_transition_covariance(
        dt=dt,
        rho=rho,
        kappa_w=kappa_w,
        sigma_w=sigma_w,
        kappa_h=kappa_h,
        sigma_h=sigma_h,
    )

    wind_variance = covariance_matrix[0, 0]

    humidity_variance = covariance_matrix[1, 1]

    covariance = covariance_matrix[0, 1]

    # If either process is deterministic, the
    # innovation covariance is zero.
    # No correlation transformation is necessary.
    if wind_variance == 0 or humidity_variance == 0:
        return 0.0

    correlation = (
        covariance
        / np.sqrt(
            wind_variance * humidity_variance
        )
    )

    # Protect against floating-point rounding errors.
    return float(
        np.clip(correlation, -1.0, 1.0)
    )


def main():

    correlations = [
        -1.0,
        -0.75,
        -0.50,
        -0.25,
        0.0,
        0.25,
        0.50,
        0.75,
        1.0,
    ]

    print("ANALYTICAL ENVIRONMENTAL DEPENDENCE")
    print("=" * 85)

    print(
        f"{'Rho':>8}"
        f"{'Covariance at T':>22}"
        f"{'Stationary Cov':>22}"
        f"{'Stationary Corr':>22}"
    )

    for rho in correlations:

        covariance = environmental_covariance(T, rho)

        limiting_covariance = stationary_covariance(rho)

        limiting_correlation = stationary_correlation(rho)

        print(
            f"{rho:8.2f}"
            f"{covariance:22.5f}"
            f"{limiting_covariance:22.5f}"
            f"{limiting_correlation:22.5f}"
        )

    print("\nEXACT OU TRANSITION COVARIANCE")
    print("=" * 85)

    rho = -0.5

    covariance_matrix = ou_transition_covariance(
        dt=DT,
        rho=rho,
    )

    print(f"Brownian correlation: {rho}")
    print(f"Time step: {DT}")
    print("\nCovariance matrix:")
    print(covariance_matrix)

    print("\nMATHEMATICAL CONSISTENCY CHECKS")
    print("=" * 85)

    # Deterministic initial conditions imply zero covariance.
    assert np.isclose(
        environmental_covariance(0, 0.5),
        0,
    )

    # Independent Brownian drivers imply zero covariance.
    assert np.isclose(
        environmental_covariance(T, 0),
        0,
    )

    # Changing the sign of rho changes the covariance sign.
    assert np.isclose(
        environmental_covariance(T, -0.5),
        -environmental_covariance(T, 0.5),
    )

    # Long-run covariance approaches its stationary value.
    assert np.isclose(
        environmental_covariance(100, 0.5),
        stationary_covariance(0.5),
    )

    # All transition covariance matrices must be
    # positive semidefinite.
    for rho in correlations:

        matrix = ou_transition_covariance(DT, rho)

        eigenvalues = np.linalg.eigvalsh(matrix)

        assert np.all(eigenvalues >= -1e-12)

    print("All analytical consistency checks passed.")
    
    analyze_correlated_intensity()
    

if __name__ == "__main__":
    main()