import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

from src.config import (
    KAPPA_W,
    SIGMA_W,
    KAPPA_H,
    SIGMA_H,
    T,
    DT,
    SEED,
    N_PATHS,
)

from src.dependence import (
    validate_correlation,
    environmental_covariance,
    stationary_correlation,
    integrated_ou_covariance,
    ou_innovation_correlation,
)


# ----------------------------------------------------------
# Experiment configuration
# ----------------------------------------------------------

CORRELATIONS = [-1.0, -0.5, 0.0, 0.5, 1.0]

WIND_Z0 = 1.5

HUMIDITY_Z0 = -1.0

RESULTS_DIR = Path("results")

FIGURES_DIR = RESULTS_DIR / "figures"


def simulate_correlated_environment(
    rho,
    n_paths=N_PATHS,
    horizon=T,
    dt=DT,
    seed=SEED,
):
    """
    Simulate correlated wind and humidity OU processes.

    Store complete environmental trajectories to validate
    marginal moments, covariance, correlation, and
    integrated environmental covariance.
    """

    validate_correlation(rho)

    # --------------------------------------------------
    # 1. Construct numerical time grid
    # --------------------------------------------------

    n_steps = int(np.ceil(horizon / dt))

    dt = horizon / n_steps

    times = np.linspace(
        0,
        horizon,
        n_steps + 1,
    )

    # --------------------------------------------------
    # 2. Initialize independent random streams
    # --------------------------------------------------

    # Use the same environmental random streams as
    # the main Monte Carlo simulation engine.

    streams = np.random.SeedSequence(seed).spawn(3)

    wind_rng = np.random.default_rng(streams[0])

    humidity_rng = np.random.default_rng(streams[1])

    # --------------------------------------------------
    # 3. Initialize environmental states
    # --------------------------------------------------

    wind = np.full(
        n_paths,
        WIND_Z0,
        dtype=float,
    )

    humidity = np.full(
        n_paths,
        HUMIDITY_Z0,
        dtype=float,
    )

    wind_paths = np.empty(
        (n_paths, n_steps + 1)
    )

    humidity_paths = np.empty(
        (n_paths, n_steps + 1)
    )

    wind_paths[:, 0] = wind

    humidity_paths[:, 0] = humidity

    # --------------------------------------------------
    # 4. Precompute exact OU transition coefficients
    # --------------------------------------------------

    wind_decay = np.exp(-KAPPA_W * dt)

    humidity_decay = np.exp(-KAPPA_H * dt)

    wind_noise_std = SIGMA_W * np.sqrt(
        -np.expm1(-2 * KAPPA_W * dt)
        / (2 * KAPPA_W)
    )

    humidity_noise_std = SIGMA_H * np.sqrt(
        -np.expm1(-2 * KAPPA_H * dt)
        / (2 * KAPPA_H)
    )

    # Exact correlation between OU innovations.
    innovation_corr = ou_innovation_correlation(
        dt=dt,
        rho=rho,
    )

    independent_weight = np.sqrt(
        max(0.0, 1.0 - innovation_corr**2)
    )

    # --------------------------------------------------
    # 5. Simulate correlated environmental trajectories
    # --------------------------------------------------

    for n in range(n_steps):

        epsilon_w = wind_rng.normal(size=n_paths)

        epsilon_h = humidity_rng.normal(size=n_paths)

        # Introduce environmental dependence.
        if innovation_corr != 0.0:

            epsilon_h = (
                innovation_corr * epsilon_w
                + independent_weight * epsilon_h
            )

        # Exact OU transitions.
        wind = (
            wind_decay * wind
            + wind_noise_std * epsilon_w
        )

        humidity = (
            humidity_decay * humidity
            + humidity_noise_std * epsilon_h
        )

        wind_paths[:, n + 1] = wind

        humidity_paths[:, n + 1] = humidity

    return times, wind_paths, humidity_paths


def theoretical_environment_statistics(times, rho):
    """
    Calculate theoretical environmental means,
    variances, covariance, and correlation.
    """

    # OU marginal means.
    wind_mean = WIND_Z0 * np.exp(
        -KAPPA_W * times
    )

    humidity_mean = HUMIDITY_Z0 * np.exp(
        -KAPPA_H * times
    )

    # OU marginal variances.
    wind_variance = (
        SIGMA_W**2
        * (-np.expm1(-2 * KAPPA_W * times))
        / (2 * KAPPA_W)
    )

    humidity_variance = (
        SIGMA_H**2
        * (-np.expm1(-2 * KAPPA_H * times))
        / (2 * KAPPA_H)
    )

    # Theoretical environmental covariance.
    covariance = environmental_covariance(
        t=times,
        rho=rho,
    )

    denominator = np.sqrt(
        wind_variance * humidity_variance
    )

    # Correlation is undefined when either variance
    # equals zero, particularly at t = 0.
    correlation = np.full(
        len(times),
        np.nan,
    )

    np.divide(
        covariance,
        denominator,
        out=correlation,
        where=denominator > 0,
    )

    return {
        "wind_mean": wind_mean,
        "humidity_mean": humidity_mean,
        "wind_variance": wind_variance,
        "humidity_variance": humidity_variance,
        "covariance": covariance,
        "correlation": correlation,
    }


def empirical_environment_statistics(
    wind_paths,
    humidity_paths,
):
    """
    Estimate environmental statistics across Monte
    Carlo trajectories at every observation time.
    """

    n_paths = wind_paths.shape[0]

    # Empirical marginal means.
    wind_mean = np.mean(
        wind_paths,
        axis=0,
    )

    humidity_mean = np.mean(
        humidity_paths,
        axis=0,
    )

    # Empirical marginal variances.
    wind_variance = np.var(
        wind_paths,
        axis=0,
        ddof=1,
    )

    humidity_variance = np.var(
        humidity_paths,
        axis=0,
        ddof=1,
    )

    # Center each process around its empirical mean.
    centered_wind = wind_paths - wind_mean

    centered_humidity = humidity_paths - humidity_mean

    # Sample covariance at each observation time.
    covariance = np.einsum(
        "ij,ij->j",
        centered_wind,
        centered_humidity,
    ) / (n_paths - 1)

    denominator = np.sqrt(
        wind_variance * humidity_variance
    )

    correlation = np.full(
        wind_paths.shape[1],
        np.nan,
    )

    np.divide(
        covariance,
        denominator,
        out=correlation,
        where=denominator > 0,
    )

    # Initial environmental conditions are deterministic.
    correlation[0] = np.nan

    return {
        "wind_mean": wind_mean,
        "humidity_mean": humidity_mean,
        "wind_variance": wind_variance,
        "humidity_variance": humidity_variance,
        "covariance": covariance,
        "correlation": correlation,
    }


def validate_integrated_covariance(
    times,
    wind_paths,
    humidity_paths,
    rho,
):
    """
    Compare empirical covariance of integrated OU
    processes with the continuous-time analytical result.

    Use left-endpoint integration to match the intensity
    drift discretization in the Monte Carlo engine.
    """

    dt = times[1] - times[0]

    # J_W = integral of wind over [0, T].
    integrated_wind = np.sum(
        wind_paths[:, :-1],
        axis=1,
    ) * dt

    # J_H = integral of humidity over [0, T].
    integrated_humidity = np.sum(
        humidity_paths[:, :-1],
        axis=1,
    ) * dt

    empirical_covariance = np.cov(
        integrated_wind,
        integrated_humidity,
        ddof=1,
    )[0, 1]

    theoretical_covariance = integrated_ou_covariance(
        rho=rho,
        horizon=times[-1],
    )

    return {
        "integrated_cov_theory": theoretical_covariance,
        "integrated_cov_empirical": empirical_covariance,
        "integrated_cov_error": (
            empirical_covariance
            - theoretical_covariance
        ),
    }


def covariance_standard_error(
    wind_variance,
    humidity_variance,
    covariance,
    n_paths,
):
    """
    Calculate the standard error of the sample covariance
    for jointly Gaussian observations.
    """

    return np.sqrt(
        (
            wind_variance * humidity_variance
            + covariance**2
        )
        / (n_paths - 1)
    )


def build_validation_summary(
    rho,
    theory,
    empirical,
    integrated_results,
):
    """
    Construct a numerical validation summary at the
    final observation time.
    """

    theoretical_cov = theory["covariance"][-1]

    empirical_cov = empirical["covariance"][-1]

    covariance_error = (
        empirical_cov - theoretical_cov
    )

    standard_error = covariance_standard_error(
        wind_variance=theory["wind_variance"][-1],
        humidity_variance=theory["humidity_variance"][-1],
        covariance=theoretical_cov,
        n_paths=N_PATHS,
    )

    covariance_z = (
        covariance_error / standard_error
    )

    result = {
        "rho": rho,

        "wind_mean_theory":
            theory["wind_mean"][-1],

        "wind_mean_empirical":
            empirical["wind_mean"][-1],

        "wind_variance_theory":
            theory["wind_variance"][-1],

        "wind_variance_empirical":
            empirical["wind_variance"][-1],

        "humidity_mean_theory":
            theory["humidity_mean"][-1],

        "humidity_mean_empirical":
            empirical["humidity_mean"][-1],

        "humidity_variance_theory":
            theory["humidity_variance"][-1],

        "humidity_variance_empirical":
            empirical["humidity_variance"][-1],

        "covariance_theory":
            theoretical_cov,

        "covariance_empirical":
            empirical_cov,

        "covariance_error":
            covariance_error,

        "covariance_standard_error":
            standard_error,

        "covariance_z":
            covariance_z,

        "correlation_theory":
            theory["correlation"][-1],

        "correlation_empirical":
            empirical["correlation"][-1],

        "stationary_correlation":
            stationary_correlation(rho),
    }

    result.update(integrated_results)

    return result


def plot_covariance_validation(results):
    """
    Compare empirical and theoretical environmental
    covariance over time for each correlation scenario.
    """

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(16, 9),
        sharex=True,
        sharey=True,
    )

    axes = axes.flatten()

    for ax, rho in zip(axes, CORRELATIONS):

        result = results[rho]

        times = result["times"]

        theory = result["theory"]

        empirical = result["empirical"]

        ax.plot(
            times,
            theory["covariance"],
            label="Theoretical",
            linewidth=2.5,
        )

        ax.plot(
            times,
            empirical["covariance"],
            linestyle="--",
            label="Monte Carlo",
            linewidth=1.5,
        )

        ax.axhline(
            0,
            color="gray",
            linewidth=0.8,
        )

        ax.set_title(f"Brownian correlation = {rho}")

        ax.set_xlabel("Time")

        ax.set_ylabel("Environmental covariance")

        ax.grid(alpha=0.25)

        ax.legend()

    # Remove unused subplot.
    axes[-1].axis("off")

    fig.suptitle(
        "Validation of Correlated OU Covariance",
        fontsize=16,
    )

    fig.tight_layout()

    fig.savefig(
        FIGURES_DIR / "dependence_ou_covariance.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_correlation_validation(results):
    """
    Compare empirical environmental correlations
    against their theoretical values.
    """

    fig, ax = plt.subplots(
        figsize=(11, 7)
    )

    colors = plt.cm.coolwarm(
        np.linspace(0, 1, len(CORRELATIONS))
    )

    for rho, color in zip(CORRELATIONS, colors):

        result = results[rho]

        times = result["times"]

        ax.plot(
            times,
            result["theory"]["correlation"],
            color=color,
            linewidth=2.5,
            label=f"Theory: rho = {rho}",
        )

        ax.plot(
            times,
            result["empirical"]["correlation"],
            color=color,
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
        )

    ax.axhline(
        0,
        color="gray",
        linewidth=0.8,
    )

    ax.set_xlabel("Time")

    ax.set_ylabel("Environmental correlation")

    ax.set_title(
        "Theoretical and Empirical OU Correlation"
    )

    ax.set_ylim(-1.05, 1.05)

    ax.grid(alpha=0.25)

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        FIGURES_DIR / "dependence_ou_correlation.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_marginal_validation(results):
    """
    Verify that environmental dependence does not
    change the marginal OU means and variances.
    """

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(14, 10),
    )

    statistics = [
        ("wind_mean", "Wind mean"),
        ("humidity_mean", "Humidity mean"),
        ("wind_variance", "Wind variance"),
        ("humidity_variance", "Humidity variance"),
    ]

    for ax, (key, title) in zip(
        axes.flatten(),
        statistics,
    ):

        # All correlation scenarios have the same
        # theoretical marginal moments.
        baseline = results[0.0]

        ax.plot(
            baseline["times"],
            baseline["theory"][key],
            color="black",
            linewidth=2.5,
            label="Theoretical",
        )

        for rho in CORRELATIONS:

            result = results[rho]

            ax.plot(
                result["times"],
                result["empirical"][key],
                linestyle="--",
                linewidth=1.2,
                label=f"rho = {rho}",
            )

        ax.set_title(title)

        ax.set_xlabel("Time")

        ax.grid(alpha=0.25)

        ax.legend()

    fig.suptitle(
        "Marginal OU Moments Under Dependence",
        fontsize=16,
    )

    fig.tight_layout()

    fig.savefig(
        FIGURES_DIR / "dependence_ou_marginals.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_integrated_covariance(summary):
    """
    Compare theoretical and empirical covariance
    between integrated environmental processes.
    """

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        summary["rho"],
        summary["integrated_cov_theory"],
        marker="o",
        linewidth=2.5,
        label="Theoretical",
    )

    ax.plot(
        summary["rho"],
        summary["integrated_cov_empirical"],
        marker="s",
        linestyle="--",
        linewidth=1.8,
        label="Monte Carlo",
    )

    ax.axhline(
        0,
        color="gray",
        linewidth=0.8,
    )

    ax.set_xlabel("Brownian correlation")

    ax.set_ylabel("Integrated environmental covariance")

    ax.set_title(
        "Validation of Integrated OU Covariance"
    )

    ax.grid(alpha=0.25)

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        FIGURES_DIR
        / "dependence_integrated_covariance.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def main():

    RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("CORRELATED OU VALIDATION")
    print("=" * 85)

    print(f"Number of paths per scenario: {N_PATHS}")

    print(f"Simulation horizon: {T}")

    print(f"Time step: {DT}")

    results = {}

    summary_rows = []

    # --------------------------------------------------
    # Run environmental dependence scenarios
    # --------------------------------------------------

    for rho in CORRELATIONS:

        print(f"\nSimulating rho = {rho:.2f} ...")

        times, wind_paths, humidity_paths = (
            simulate_correlated_environment(rho)
        )

        theory = theoretical_environment_statistics(
            times,
            rho,
        )

        empirical = empirical_environment_statistics(
            wind_paths,
            humidity_paths,
        )

        integrated_results = validate_integrated_covariance(
            times,
            wind_paths,
            humidity_paths,
            rho,
        )

        summary = build_validation_summary(
            rho,
            theory,
            empirical,
            integrated_results,
        )

        summary_rows.append(summary)

        results[rho] = {
            "times": times,
            "theory": theory,
            "empirical": empirical,
        }

        # Release complete trajectory matrices before
        # running the next correlation scenario.
        del wind_paths, humidity_paths

    summary = pd.DataFrame(summary_rows)

    # --------------------------------------------------
    # Print covariance validation
    # --------------------------------------------------

    print("\nENVIRONMENTAL COVARIANCE VALIDATION")
    print("=" * 100)

    columns = [
        "rho",
        "covariance_theory",
        "covariance_empirical",
        "covariance_error",
        "covariance_z",
    ]

    print(
        summary[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Print correlation validation
    # --------------------------------------------------

    print("\nENVIRONMENTAL CORRELATION VALIDATION")
    print("=" * 100)

    columns = [
        "rho",
        "correlation_theory",
        "correlation_empirical",
        "stationary_correlation",
    ]

    print(
        summary[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Print marginal variance validation
    # --------------------------------------------------

    print("\nMARGINAL VARIANCE VALIDATION")
    print("=" * 100)

    columns = [
        "rho",
        "wind_variance_theory",
        "wind_variance_empirical",
        "humidity_variance_theory",
        "humidity_variance_empirical",
    ]

    print(
        summary[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Print integrated covariance validation
    # --------------------------------------------------

    print("\nINTEGRATED OU COVARIANCE VALIDATION")
    print("=" * 100)

    columns = [
        "rho",
        "integrated_cov_theory",
        "integrated_cov_empirical",
        "integrated_cov_error",
    ]

    print(
        summary[columns].to_string(
            index=False,
            float_format=lambda x: f"{x:.5f}",
        )
    )

    # --------------------------------------------------
    # Save numerical results
    # --------------------------------------------------

    summary.to_csv(
        RESULTS_DIR / "dependence_ou_validation.csv",
        index=False,
    )

    # --------------------------------------------------
    # Generate validation figures
    # --------------------------------------------------

    plot_covariance_validation(results)

    plot_correlation_validation(results)

    plot_marginal_validation(results)

    plot_integrated_covariance(summary)

    print("\nValidation results and figures saved.")

    print("Correlated OU validation completed.")


if __name__ == "__main__":
    main()
    
# ==========================================================
# Q6 — CORRELATED OU NUMERICAL VALIDATION RESULTS
# ==========================================================
#
# Simulation configuration:
#
# Number of trajectories per scenario: 10000
# Simulation horizon: 10.0
# Time step: 0.01
# Number of correlation scenarios: 5
# Total environmental trajectory pairs: 50000
#
# Brownian correlation scenarios:
#     rho = -1.0, -0.5, 0.0, 0.5, 1.0
#
# ----------------------------------------------------------
# 1. Environmental covariance validation
# ----------------------------------------------------------
#
# Rho     Theoretical    Empirical      Z-score
# -1.0    -0.53846       -0.52302       1.99921
# -0.5    -0.26923       -0.25995       1.50862
#  0.0     0.00000        0.00181       0.32691
#  0.5     0.26923        0.26308      -0.99946
#  1.0     0.53846        0.52303      -1.99880
#
# The empirical covariance closely follows the theoretical
# linear relationship with Brownian correlation.
#
# The extreme correlation scenarios have standardized
# errors of approximately two standard errors.
#
# These discrepancies are compatible with ordinary
# Monte Carlo variation and do not independently establish
# an implementation error.
#
# Common random numbers introduce dependence between
# estimates obtained under different correlation scenarios.
#
# ----------------------------------------------------------
# 2. Environmental correlation validation
# ----------------------------------------------------------
#
# Rho     Theoretical    Empirical
# -1.0    -0.97303       -0.97232
# -0.5    -0.48651       -0.48283
#  0.0     0.00000        0.00335
#  0.5     0.48651        0.48755
#  1.0     0.97303        0.97233
#
# All absolute correlation errors are below 0.004.
#
# The stationary environmental correlation is:
#
# Corr(Z_W, Z_H) =
#     2 * rho * sqrt(kappa_w * kappa_h)
#     / (kappa_w + kappa_h)
#
# Under our parameters:
#
#     Stationary correlation = 0.97301 * rho
#
# Perfect Brownian correlation does not imply perfect
# environmental state correlation because the two OU
# processes have different mean-reversion rates.
#
# ----------------------------------------------------------
# 3. Marginal variance validation
# ----------------------------------------------------------
#
# Theoretical wind variance:     0.62500
# Empirical wind variance:       0.60427
#
# Theoretical humidity variance: 0.48998
#
# Rho     Empirical humidity variance
# -1.0    0.47884
# -0.5    0.47967
#  0.0    0.48139
#  0.5    0.48184
#  1.0    0.47884
#
# Wind trajectories are identical across correlation
# scenarios because the same wind random stream is used.
#
# Dependence is introduced by transforming the humidity
# innovations while preserving their unit variance.
#
# The theoretical marginal OU distributions remain
# unchanged under different Brownian correlations.
#
# The empirical wind variance is approximately 3.32%
# below theory, which is about 2.35 standard errors.
#
# ----------------------------------------------------------
# 4. Integrated OU covariance validation
# ----------------------------------------------------------
#
# Rho     Theoretical    Empirical      Error
# -1.0    -13.18297      -13.10308       0.07989
# -0.5     -6.59148       -6.68395      -0.09247
#  0.0      0.00000       -0.15297      -0.15297
#  0.5      6.59148        6.41899      -0.17249
#  1.0     13.18297       13.10281      -0.08016
#
# Relative absolute errors for nonzero correlations
# are below 3%.
#
# The independent scenario has a theoretical covariance
# of zero, so a relative error is not defined.
#
# Integrated environmental covariance is particularly
# important because it contributes directly to the
# variance of terminal log-intensity in the additive
# wildfire model.
#
# Differences between analytical and empirical integrated
# covariance may arise from both Monte Carlo sampling
# variation and left-endpoint integration error.
#
# ----------------------------------------------------------
# 5. Conclusion
# ----------------------------------------------------------
#
# The numerical results support the correctness of the
# correlated OU transition construction.
#
# Environmental dependence changes the joint distribution
# of wind and humidity while preserving their theoretical
# marginal distributions.
#
# The results also support the analytical covariance
# formulas for the integrated environmental processes.
#
# The validation concerns numerical consistency with
# the assumed stochastic model, not predictive accuracy
# against observed wildfire behavior.
#
# The Brownian correlation parameter remains a scenario
# parameter rather than a calibrated historical estimate.
# ==========================================================