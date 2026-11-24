import numpy as np

from src.config import GAMMA, N_PATHS, SEED

from src.monte_carlo import simulate_many


def main():

    print("DEPENDENCE ENGINE VALIDATION")
    print("=" * 85)

    # -------------------------------------------------
    # Test 1: Default and explicit independence
    # -------------------------------------------------

    times_default, paths_default = simulate_many(
        n_paths=N_PATHS,
        seed=SEED,
    )

    times_zero, paths_zero = simulate_many(
        n_paths=N_PATHS,
        seed=SEED,
        rho=0.0,
    )

    assert np.array_equal(
        times_default,
        times_zero,
    )

    assert np.array_equal(
        paths_default,
        paths_zero,
    )

    print(
        "PASS: Default simulation equals explicit rho=0."
    )

    # -------------------------------------------------
    # Test 2: Recover previous Q2 baseline
    # -------------------------------------------------

    terminal = paths_zero[:, -1]

    mean = np.mean(terminal)

    std = np.std(terminal, ddof=1)

    print("\nINDEPENDENT INTERACTION BASELINE")
    print("=" * 85)

    print(f"Mean: {mean:.5f}")

    print(f"Standard deviation: {std:.5f}")

    assert np.isclose(
        mean,
        1596.96718,
        rtol=0,
        atol=1e-5,
    )

    assert np.isclose(
        std,
        1438.32228,
        rtol=0,
        atol=1e-5,
    )

    print("PASS: Previous Q2 baseline reproduced.")

    del paths_default, paths_zero

    # -------------------------------------------------
    # Test 3: Investigate different correlations
    # -------------------------------------------------

    correlations = [
        -1.0,
        -0.5,
        0.0,
        0.5,
        1.0,
    ]

    print("\nCORRELATED SIMULATION RESULTS")
    print("=" * 85)

    print(
        f"{'Rho':>8}"
        f"{'Mean':>15}"
        f"{'Std':>15}"
        f"{'Median':>15}"
        f"{'Log Mean':>15}"
        f"{'Log Variance':>15}"
    )

    for rho in correlations:

        times, paths = simulate_many(
            n_paths=N_PATHS,
            gamma=GAMMA,
            seed=SEED,
            rho=rho,
        )

        # Check simulation dimensions.
        assert paths.shape == (
            N_PATHS,
            len(times),
        )

        # Check mathematical positivity.
        assert np.all(paths > 0)

        # Check numerical stability.
        assert np.all(np.isfinite(paths))

        terminal = paths[:, -1]

        log_terminal = np.log(terminal)

        print(
            f"{rho:8.2f}"
            f"{np.mean(terminal):15.5f}"
            f"{np.std(terminal, ddof=1):15.5f}"
            f"{np.median(terminal):15.5f}"
            f"{np.mean(log_terminal):15.5f}"
            f"{np.var(log_terminal, ddof=1):15.5f}"
        )

        del terminal, log_terminal, paths

    # -------------------------------------------------
    # Test 4: Invalid correlation
    # -------------------------------------------------

    try:

        simulate_many(
            n_paths=10,
            rho=1.5,
        )

    except ValueError:

        print(
            "\nPASS: Invalid correlation rejected."
        )

    else:

        raise AssertionError(
            "Invalid correlation was not rejected."
        )

    print("\nAll dependence engine checks passed.")


if __name__ == "__main__":
    main()