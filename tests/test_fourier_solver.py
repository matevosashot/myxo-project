from pathlib import Path

import h5py
import numpy as np
import pytest

from myxo import FourierSolver
from myxo.fourier_solver import apply_fft, apply_ifft


REFERENCE_H5 = Path(__file__).parent / "biharmonic_solution.h5"

# Parameters fixed to match draft_p_fluct_visualisation.ipynb so that the
# reference output in biharmonic_solution.h5 (produced by the old
# draft.BiharmonicSolverFourier code) is the ground truth here.
N = 10
L = 5
MODEL_PARAMS = dict(S0=1.0, l=1.0, lM=7.0 / 3, lm=0.7 / 3, R=10)


def test_rfft_and_ifft_match():
    n = 10
    L = 10.0
    dx = 2.0 * L / n
    pts = np.linspace(-L, L - dx, n)
    x = np.random.randn(n, n, n, n)
    x_hat = apply_fft(x)
    x_hat_r = apply_ifft(x_hat)
    np.testing.assert_allclose(x, x_hat_r)


@pytest.fixture(scope="module")
def reference():
    with h5py.File(REFERENCE_H5, "r") as f:
        return {"f": f["f"][:], "fdiag": f["fdiag"][:], "grid": f["grid"][:]}


@pytest.fixture(scope="module")
def solver():
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32)
    s.precompute()
    return s


def test_grid_matches_reference(solver, reference):
    np.testing.assert_allclose(solver.grid_1d, reference["grid"], rtol=1e-6)


def test_P_contribution_matches_old_solver(solver, reference):
    # The old draft.BiharmonicSolverFourier solver and the new FourierSolver
    # use opposite sign conventions, so P_contrib ≈ -f_old.
    P = solver.solve_P_contribution()
    assert P.shape == reference["f"].shape
    np.testing.assert_allclose(P, -reference["f"], atol=1e-5)


def test_P_contribution_diag_matches_old_solver(solver, reference):
    P = solver.solve_P_contribution()
    idx = np.arange(N)
    P_diag = P[idx[:, None], idx[None, :], idx[:, None], idx[None, :]]
    assert P_diag.shape == reference["fdiag"].shape
    np.testing.assert_allclose(P_diag, -reference["fdiag"], atol=1e-5)
