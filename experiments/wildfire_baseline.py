
import matplotlib.pyplot as plt
from pathlib import Path

from src.simulation import simulate_wildfire
from src.config import GAMMA, I0


def main():
    times, wind, humidity, intensity = simulate_wildfire()

    fig, axes = plt.subplots(
        3, 1,
        figsize=(10, 9),
        sharex=True
    )

    axes[0].plot(times, wind)
    axes[0].axhline(0, color="black", linestyle="--")
    axes[0].set_ylabel("Wind anomaly")
    axes[0].grid(alpha=0.3)

    axes[1].plot(times, humidity)
    axes[1].axhline(0, color="black", linestyle="--")
    axes[1].set_ylabel("Humidity anomaly")
    axes[1].grid(alpha=0.3)

    axes[2].plot(times, intensity)
    axes[2].axhline(
        I0,
        color="black",
        linestyle="--",
        label="Initial intensity"
    )
    axes[2].set_ylabel("Intensity (kW/m)")
    axes[2].set_xlabel("Time (illustrative units)")
    axes[2].legend()
    axes[2].grid(alpha=0.3)

    fig.suptitle(
        f"Coupled Wildfire SDE Simulation (gamma={GAMMA})"
    )

    fig.tight_layout()

    output_dir = Path("results/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_dir / "wildfire_baseline.png",
        dpi=200,
        bbox_inches="tight"
    )

    print(f"Initial intensity: {intensity[0]:.2f} kW/m")
    print(f"Final intensity:   {intensity[-1]:.2f} kW/m")
    print(f"Maximum intensity: {intensity.max():.2f} kW/m")

    plt.show()


if __name__ == "__main__":
    main()


# =====================================================
# Insights from the coupled wildfire simulation
# =====================================================
#
# 1. Wildfire intensity remains positive throughout
#    the simulation, consistent with the exponential
#    formulation of the intensity process.
#
# 2. During the early simulation, above-average wind
#    and below-average humidity favor intensity growth.
#
# 3. Intensity increases from 500 kW/m to a maximum
#    of approximately 890 kW/m.
#
# 4. As wind weakens and humidity increases, intensity
#    generally decreases.
#
# 5. The final intensity is approximately 260 kW/m.
#
# 6. Intensity does not change monotonically because
#    both environmental conditions and intensity
#    contain stochastic components.
#
# 7. These observations are consistent with the
#    mathematical structure of our coupled SDE model.
#
# Limitations:
# - This is only one simulated realization.
# - Intensity parameters are illustrative, not fitted.
# - We cannot draw general statistical conclusions
#   from a single trajectory.
# - The interaction effect cannot be isolated without
#   comparing against a simulation with gamma = 0.
