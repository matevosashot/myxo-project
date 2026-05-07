import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from myxo import FourierSolver
from myxo.fourier_solver import apply_fft, apply_ifft


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
def solver():
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32, mmap=False)
    s.precompute()
    return s



def test_Q_contribution_symmetry_optimized_matches_naive():
    # Use float64 so the naive and optimized paths agree to tight tolerance;
    # any disagreement here is from index-symmetry bookkeeping, not roundoff.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s.precompute()
    Q_naive = s.solve_Q_contribution()
    Q_opt = s.solve_Q_contribution__symmetry_optimized()
    assert Q_opt.shape == Q_naive.shape
    np.testing.assert_allclose(Q_opt, Q_naive, rtol=1e-10, atol=1e-12)


def test_P_contribution_symmetry_optimized_matches_naive():
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s.precompute()
    P_naive = s.solve_P_contribution()
    P_opt = s.solve_P_contribution__symmetry_optimized()
    assert P_opt.shape == P_naive.shape
    np.testing.assert_allclose(P_opt, P_naive, rtol=1e-10, atol=1e-12)


MMAP_TEST_PATH = "/scratch/phi_mmap_for_test.mmap"

requires_scratch = pytest.mark.skipif(
    not os.path.isdir(os.path.dirname(MMAP_TEST_PATH)),
    reason=f"{os.path.dirname(MMAP_TEST_PATH)} not present on this host",
)


@pytest.fixture
def mmap_test_path():
    # Distinct from the production default (/scratch/phi_temp.mmap) so tests
    # never clobber a real run sharing the same node.
    if os.path.exists(MMAP_TEST_PATH):
        os.unlink(MMAP_TEST_PATH)
    yield MMAP_TEST_PATH
    # Defensive: precompute() unlinks immediately on success, so this only
    # fires if the test failed before the unlink ran.
    if os.path.exists(MMAP_TEST_PATH):
        os.unlink(MMAP_TEST_PATH)


@requires_scratch
def test_phi_is_memmap_when_enabled(mmap_test_path):
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32,
                      mmap=True, mmap_path=mmap_test_path)
    s.precompute()
    assert isinstance(s.gcalc.phi, np.memmap)
    assert s.gcalc.phi.shape == (N, N, N, N)
    assert s.gcalc.phi.dtype == np.float32


def test_phi_is_not_memmap_when_disabled():
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32, mmap=False)
    s.precompute()
    assert not isinstance(s.gcalc.phi, np.memmap)


@requires_scratch
def test_mmap_phi_matches_in_memory_phi(mmap_test_path):
    s_mmap = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32,
                           mmap=True, mmap_path=mmap_test_path)
    s_mmap.precompute()
    s_mem = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32, mmap=False)
    s_mem.precompute()
    np.testing.assert_array_equal(np.asarray(s_mmap.gcalc.phi), s_mem.gcalc.phi)


@requires_scratch
def test_mmap_solve_P_matches_in_memory(mmap_test_path):
    s_mmap = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                           mmap=True, mmap_path=mmap_test_path)
    s_mmap.precompute()
    P_mmap = s_mmap.solve_P_contribution()

    s_mem = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s_mem.precompute()
    P_mem = s_mem.solve_P_contribution()

    np.testing.assert_array_equal(P_mmap, P_mem)


@requires_scratch
def test_mmap_solve_Q_matches_in_memory(mmap_test_path):
    s_mmap = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                           mmap=True, mmap_path=mmap_test_path)
    s_mmap.precompute()
    Q_mmap = s_mmap.solve_Q_contribution()

    s_mem = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s_mem.precompute()
    Q_mem = s_mem.solve_Q_contribution()

    np.testing.assert_array_equal(Q_mmap, Q_mem)


@requires_scratch
def test_mmap_file_is_unlinked_after_precompute(mmap_test_path):
    # The memmap is created at mmap_path then immediately unlinked; the
    # kernel keeps the inode alive while the np.memmap exists, but the path
    # itself should no longer resolve.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32,
                      mmap=True, mmap_path=mmap_test_path)
    s.precompute()
    assert not os.path.exists(mmap_test_path)
    # phi is still usable despite the path being gone.
    assert s.gcalc.phi.sum() != 0
