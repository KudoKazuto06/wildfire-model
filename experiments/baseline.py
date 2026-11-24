import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from src.environment import simulate_ou
from src.config import (
    KAPPA_W, SIGMA_W,
    KAPPA_H, SIGMA_H,
    T, DT, SEED
)


def main():
    rng = np.random.default_rng(SEED)

    # Simulate wind
    times, wind = simulate_ou(
        z0=1.5,
        kappa=KAPPA_W,
        sigma=SIGMA_W,
        T=T,
        dt=DT,
        rng=rng
    )

    # Simulate humidity
    _, humidity = simulate_ou(
        z0=-1.0,
        kappa=KAPPA_H,
        sigma=SIGMA_H,
        T=T,
        dt=DT,
        rng=rng
    )

    # Create figure
    plt.figure(figsize=(10, 5))

    plt.plot(times, wind, label="Wind anomaly")
    plt.plot(times, humidity, label="Humidity anomaly")

    plt.axhline(0, color="black", linestyle="--")

    plt.xlabel("Time (illustrative units)")
    plt.ylabel("Standardized anomaly")
    plt.title("OU Simulation: Wind and Humidity")

    plt.legend()
    plt.grid(alpha=0.3)

    # Save the figure
    output_dir = Path("results/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(
        output_dir / "environment_baseline.png",
        dpi=200,
        bbox_inches="tight"
    )

    plt.show()


if __name__ == "__main__":
    main()


# Insights from the baseline OU simulation
#
# 1. Both wind and humidity exhibit random fluctuations
#    around their reference means (zero anomaly)
#
# 2. The wind anomaly begins at +1.5 and generally moves
#    toward zero initially, illustrating mean reversion
#
# 3. The humidity anomaly begins at -1.0 and also
#    fluctuates toward and away from zero
#
# 4. Mean reversion does not guarantee that every
#    trajectory stays close to zero. Random fluctuations
#    can temporarily dominate the drift
#
# 5. The two trajectories are simulated using independent
#    Brownian noise, as assumed in our mathematical model