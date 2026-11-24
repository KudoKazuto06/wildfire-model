 
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from src.environment import simulate_ou
from src.config import (
    KAPPA_W, SIGMA_W,
    KAPPA_H, SIGMA_H,
    T, DT, SEED
)


# Number of independent simulated trajectories
N_SIMULATIONS = 5000


def theoretical_ou(times, z0, kappa, sigma):
    """
    Calculate theoretical OU mean and variance.
    """

    mean = z0 * np.exp(-kappa * times)

    variance = (
        sigma**2 / (2 * kappa)
        * (1 - np.exp(-2 * kappa * times))
    )

    return mean, variance


def validate_process(name, z0, kappa, sigma, rng):
    """
    Compare Monte Carlo estimates with OU theory.
    """

    # Generate one trajectory to determine time grid
    times, first_path = simulate_ou(
        z0=z0,
        kappa=kappa,
        sigma=sigma,
        T=T,
        dt=DT,
        rng=rng
    )

    # Store all simulated trajectories
    paths = np.empty((N_SIMULATIONS, len(times)))
    paths[0] = first_path

    for i in range(1, N_SIMULATIONS):
        _, paths[i] = simulate_ou(
            z0=z0,
            kappa=kappa,
            sigma=sigma,
            T=T,
            dt=DT,
            rng=rng
        )

    # Empirical Monte Carlo statistics
    empirical_mean = np.mean(paths, axis=0)
    empirical_variance = np.var(paths, axis=0, ddof=1)

    # Theoretical statistics
    theoretical_mean, theoretical_variance = theoretical_ou(
        times, z0, kappa, sigma
    )

    # Plot the comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].plot(
        times, empirical_mean,
        label="Monte Carlo mean"
    )
    axes[0].plot(
        times, theoretical_mean,
        "--", label="Theoretical mean"
    )
    axes[0].set_title(f"{name}: Mean")
    axes[0].set_xlabel("Time")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].plot(
        times, empirical_variance,
        label="Monte Carlo variance"
    )
    axes[1].plot(
        times, theoretical_variance,
        "--", label="Theoretical variance"
    )
    axes[1].set_title(f"{name}: Variance")
    axes[1].set_xlabel("Time")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    fig.tight_layout()

    output_dir = Path("results/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_dir / f"ou_validation_{name.lower()}.png",
        dpi=200,
        bbox_inches="tight"
    )

    plt.show()

    # Numerical comparison at final time
    print(f"\n{name} OU Validation at T = {T}")
    print("-" * 45)

    print(
        f"Mean:     empirical={empirical_mean[-1]:.5f}, "
        f"theoretical={theoretical_mean[-1]:.5f}"
    )

    print(
        f"Variance: empirical={empirical_variance[-1]:.5f}, "
        f"theoretical={theoretical_variance[-1]:.5f}"
    )


def main():
    rng = np.random.default_rng(SEED)

    validate_process(
        name="Wind",
        z0=1.5,
        kappa=KAPPA_W,
        sigma=SIGMA_W,
        rng=rng
    )

    validate_process(
        name="Humidity",
        z0=-1.0,
        kappa=KAPPA_H,
        sigma=SIGMA_H,
        rng=rng
    )


if __name__ == "__main__":
    main()


# Monte Carlo OU Validation Results
#
# Number of trajectories: 5000
# Final time: T = 10.0
#
# Wind:
#   Empirical mean:         0.02007
#   Theoretical mean:       0.00050
#   Empirical variance:     0.61631
#   Theoretical variance:   0.62500
#
# Humidity:
#   Empirical mean:         -0.02274
#   Theoretical mean:       -0.00674
#   Empirical variance:     0.49803
#   Theoretical variance:   0.48998
#
# Intepretation:
#   The empirical means and variances are consistent with
#   theoretical OU predictions within Monte Carlo error
#
#   This supports the correctness of the OU simulation
#   It does not establish that the model fits real weather