import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from myxo import FourierSolver
from myxo.fourier_solver import apply_fft, apply_ifft, apply_fft_spatial, apply_ifft_spatial


N = 11
L = 5
MODEL_PARAMS = dict(S0=1.0, l=1.0, lM=7.0 / 3, lm=0.7 / 3, R=10)


def test_rfft_and_ifft_match():
    n = 10
    L = 10.0
    dx = 2.0 * L / n
    pts = np.linspace(-L, L - dx, n)
    x = np.random.randn(n, n, n, n)
    x_hat = apply_fft(x, fft_type="rfft")
    x_hat_r = apply_ifft(x_hat, fft_type="rfft", s=(n, n, n, n))
    np.testing.assert_allclose(x, x_hat_r)


@pytest.mark.parametrize("fft_type", ["rfft", "fft"])
def test_apply_fft_ifft_round_trip(fft_type):
    # Both modes must be a faithful inverse pair on real input. Use a fixed
    # seed so the test is deterministic; copy x before passing to apply_fft
    # because overwrite_x=True may clobber it.
    x = np.random.default_rng(0).standard_normal((N, N, N, N))
    x_hat = apply_fft(x.copy(), fft_type=fft_type)
    x_back = apply_ifft(x_hat, fft_type=fft_type, s=(N, N, N, N))
    np.testing.assert_allclose(x_back, x, rtol=1e-12, atol=1e-13)


@pytest.fixture(scope="module")
def solver():
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32, mmap=False)
    s.precompute()
    return s



def _reference_pressure_field(solver):
    """Independent implementation of solve_pressure_field for cross-checks."""
    Q = solver.gcalc_Q.Q
    Q_fourier = apply_fft_spatial(Q, fft_type=solver.fft_type)

    if solver.fft_type == "rfft":
        k_comp = (solver.qfull[:, None], solver.qhalf[None, :])
    else:
        k_comp = (solver.qfull[:, None], solver.qfull[None, :])

    k_inv_squared = 1.0 / (k_comp[0] ** 2 + k_comp[1] ** 2 + solver.eps)
    k_outer = np.empty(k_inv_squared.shape + (2, 2), dtype=solver.dtype)
    for a in range(2):
        for b in range(2):
            k_outer[..., a, b] = k_comp[a] * k_comp[b]

    kernel = k_outer * k_inv_squared[..., None, None]
    p_fourier = -np.sum(kernel * Q_fourier, axis=(2, 3))
    p_fourier[0, 0] = 0.0
    return apply_ifft_spatial(
        p_fourier, fft_type=solver.fft_type, s=(solver.n, solver.n),
    )


def test_solve_pressure_field_shape_and_finiteness():
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s.precompute()
    p = s.solve_pressure_field()
    assert p.shape == (N, N)
    assert np.all(np.isfinite(p))


def test_solve_pressure_field_matches_reference():
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s.precompute()
    p = s.solve_pressure_field()
    p_ref = _reference_pressure_field(s)
    np.testing.assert_allclose(p, p_ref, rtol=0, atol=0)


def test_solve_pressure_field_zero_dc_mode():
    # The solver pins the k=0 Fourier mode to zero, so the spatial mean
    # vanishes (up to float64 roundoff from the inverse transform).
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s.precompute()
    p = s.solve_pressure_field()
    np.testing.assert_allclose(p.mean(), 0.0, atol=1e-12)


def test_solve_pressure_field_rfft_matches_fft():
    s_r = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS,
                        dtype=np.float64, mmap=False, fft_type="rfft")
    s_f = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS,
                        dtype=np.float64, mmap=False, fft_type="fft")
    s_r.precompute()
    s_f.precompute()
    p_r = s_r.solve_pressure_field()
    p_f = s_f.solve_pressure_field()
    np.testing.assert_allclose(p_f, p_r, rtol=1e-10, atol=1e-12)


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


def test_gcalc_Q_default_aliases_gcalc():
    # Without gcalc_Q_cls, the Q-side reuses the P-side Gcalc — no extra
    # phi buffer, no extra precompute.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s.precompute()
    assert s.gcalc_Q is s.gcalc


def test_drop_phi_releases_inactive_gcalc_phi_during_solve():
    # If gcalc_Q is a distinct instance whose calc_phi materialises a full
    # N^4 phi (i.e. inherits v1 without overriding), the inactive side's
    # phi must be released for the duration of the active solve. Otherwise
    # the inner-loop peak goes back to 4*N^4.
    from myxo.correlation_def import Gcalc

    class _FullPhiQGcalc(Gcalc):
        """A trivial subclass with a v1-style full N^4 phi."""
        pass

    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                      mmap=False, drop_phi=True, gcalc_Q_cls=_FullPhiQGcalc)
    s.precompute()
    # Both Gcalcs have full phi resident after precompute.
    assert s.gcalc is not s.gcalc_Q
    assert s.gcalc.phi is not None
    assert s.gcalc_Q.phi is not None

    s.solve_P_contribution()
    # solve_P should have released gcalc_Q.phi (inactive side) at entry
    # and drop_phi loop should leave gcalc.phi None at exit.
    assert s.gcalc.phi is None
    assert s.gcalc_Q.phi is None

    s.solve_Q_contribution()
    # solve_Q rebuilds gcalc_Q.phi inside the loop, then drops it on the
    # last iteration. gcalc.phi should remain None (it was already None
    # and solve_Q doesn't touch it except via _drop_inactive_phi).
    assert s.gcalc.phi is None
    assert s.gcalc_Q.phi is None


def test_drop_inactive_phi_is_no_op_when_aliased():
    # No gcalc_Q_cls => gcalc_Q is gcalc; _drop_inactive_phi must not
    # touch the active side.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                      mmap=False, drop_phi=True)
    s.precompute()
    assert s.gcalc_Q is s.gcalc
    s._drop_inactive_phi(s.gcalc)
    # phi must still be present after the no-op.
    assert s.gcalc.phi is not None


def test_gcalc_Q_cls_constructs_separate_instance():
    # With gcalc_Q_cls, FourierSolver builds a second Gcalc instance for
    # the Q-side and routes solve_Q_* through it. v2 uses a different C_Q
    # formula, so we just verify the routing wiring (separate instance, of
    # the right type, P-side untouched) and that the Q-side solver runs to
    # completion with finite output.
    from myxo.correlation_def_v2 import GcalcForQ

    s_v2 = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                         mmap=False, gcalc_Q_cls=GcalcForQ)
    s_v2.precompute()
    assert s_v2.gcalc_Q is not s_v2.gcalc
    assert isinstance(s_v2.gcalc_Q, GcalcForQ)
    # P-side untouched: still the v1 Gcalc.
    assert type(s_v2.gcalc).__name__ == "Gcalc"
    assert s_v2.gcalc.__class__.__module__.endswith("correlation_def")

    Q_v2 = s_v2.solve_Q_contribution()
    assert Q_v2.shape == (N, N, N, N)
    assert np.all(np.isfinite(Q_v2))


@pytest.mark.parametrize("indices", [
    # Off-diagonal C_Q classes (ab != mn) — these are the call sites in
    # solve_Q_contribution__symmetry_optimized that route through the fused
    # path. Diagonal classes use plain calc_C_Q_index and aren't relevant.
    (0, 0, 0, 1),
    (0, 0, 1, 1),
    (0, 1, 1, 1),
])
def test_calc_C_Q_index_symmetrized_matches_two_step(indices):
    a, b, m, n = indices
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64, mmap=False)
    s.precompute()

    c = s.gcalc.calc_C_Q_index(a, b, m, n)
    expected = c + c.transpose(2, 3, 0, 1)

    actual = s.gcalc.calc_C_Q_index_symmetrized(a, b, m, n)

    # ne.evaluate processes in chunks; the fused expression and the two-step
    # path differ in operation order, so allow a few ULPs of float64 slack.
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)


def test_gcalc_drop_phi_clears_phi_and_transpose():
    # Both attributes must be set to None: phi_transpose is a numpy view of
    # phi, so leaving it bound would keep the underlying buffer alive.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32, mmap=False)
    s.precompute()
    assert s.gcalc.phi is not None
    assert s.gcalc.phi_transpose is not None

    s.gcalc.drop_phi()

    assert s.gcalc.phi is None
    assert s.gcalc.phi_transpose is None


def test_drop_phi_false_keeps_phi_resident_after_solve():
    # The default flag must leave the cached phi in place so callers can run
    # multiple solves without paying the rebuild cost.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float32,
                      mmap=False, drop_phi=False)
    s.precompute()
    s.solve_P_contribution()
    assert s.gcalc.phi is not None
    assert s.gcalc.phi_transpose is not None


def test_drop_phi_true_two_consecutive_solves():
    # After solve_P, drop_phi=True leaves phi=None. The next solve_*
    # must rebuild via _ensure_phi on its first iteration.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                      mmap=False, drop_phi=True)
    s.precompute()

    p_first = s.solve_P_contribution()
    assert s.gcalc.phi is None

    q = s.solve_Q_contribution()
    assert s.gcalc.phi is None
    assert q.shape == (N, N, N, N)

    # Re-running solve_P after a Q in between must give bit-identical output.
    p_second = s.solve_P_contribution()
    np.testing.assert_array_equal(p_first, p_second)


@requires_scratch
def test_drop_phi_with_mmap(mmap_test_path):
    # drop_phi + mmap: phi is rebuilt as a fresh np.memmap each iteration
    # (the previous one is GC'd when drop_phi() nulls the attribute, which
    # closes the fd and releases the unlinked inode). Output must still match
    # the cached-RAM path bit-for-bit.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                      mmap=True, mmap_path=mmap_test_path, drop_phi=True)
    s.precompute()
    out = s.solve_P_contribution()

    s_ref = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                          mmap=False, drop_phi=False)
    s_ref.precompute()
    out_ref = s_ref.solve_P_contribution()

    np.testing.assert_array_equal(out, out_ref)
    assert s.gcalc.phi is None  # dropped after the last iteration


@pytest.mark.parametrize("solve_method", [
    "solve_P_contribution",
    "solve_P_contribution__symmetry_optimized",
    "solve_Q_contribution",
    "solve_Q_contribution__symmetry_optimized",
])
def test_drop_phi_matches_no_drop(solve_method):
    # drop_phi recomputes phi each iteration; result must be bit-for-bit
    # identical to the cached-phi path. float64 keeps the comparison exact.
    s_keep = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                           mmap=False, drop_phi=False)
    s_keep.precompute()
    out_keep = getattr(s_keep, solve_method)()

    s_drop = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS, dtype=np.float64,
                           mmap=False, drop_phi=True)
    s_drop.precompute()
    out_drop = getattr(s_drop, solve_method)()

    np.testing.assert_array_equal(out_drop, out_keep)
    # phi was rebuilt at least once, then dropped at the end of the last
    # iteration — final state must be None for both attributes so the
    # backing buffer can be reclaimed.
    assert s_drop.gcalc.phi is None
    assert s_drop.gcalc.phi_transpose is None


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


# ---------------------------------------------------------------------------
# fft_type: rfft (default, half-spectrum on the last axis) vs fft (full
# spectrum on every axis). Both modes are mathematically the same transform
# on real input, so every solve_* entry point must agree to tight tolerance.
# ---------------------------------------------------------------------------


def test_invalid_fft_type_raises():
    with pytest.raises(ValueError):
        FourierSolver(n=N, L=L, model_params=MODEL_PARAMS,
                      dtype=np.float64, mmap=False, fft_type="bogus")


@pytest.mark.parametrize("fft_type, expected_last", [
    ("rfft", N // 2 + 1),
    ("fft",  N),
])
def test_precompute_qyp_shape_matches_fft_type(fft_type, expected_last):
    # Last axis of the q' grid must shrink to N//2+1 in rfft mode and stay
    # at N in fft mode; q2_inv_sq broadcasts off qyp so it must follow.
    s = FourierSolver(n=N, L=L, model_params=MODEL_PARAMS,
                      dtype=np.float64, mmap=False, fft_type=fft_type)
    s.precompute()
    assert s.qyp.shape[-1] == expected_last
    assert s.q2_inv_sq.shape[-1] == expected_last


@pytest.mark.parametrize("solve_method", [
    "solve_P_contribution",
    "solve_P_contribution__symmetry_optimized",
    "solve_Q_contribution",
    "solve_Q_contribution__symmetry_optimized",
])
def test_fft_type_rfft_matches_fft(solve_method):
    # Cross-mode equivalence: rfft and fft are the same transform on real
    # input, just with different storage layouts in Fourier space. Every
    # solve_* method must yield the same real-space output regardless.
    #
    # Use an odd grid size so there is no Nyquist degeneracy on any axis
    # (np.fft.fftfreq has a Nyquist mode at index n/2 only when n is even,
    # and rfftfreq's positive-Nyquist convention there disagrees with
    # fftfreq's negative-Nyquist for kernels that are odd in qyp). For odd
    # n the two grids agree pointwise on the half-spectrum, so the kernel
    # values match and the inverse transforms are bit-for-bit equivalent
    # up to float64 roundoff.
    n_odd = 11
    s_r = FourierSolver(n=n_odd, L=L, model_params=MODEL_PARAMS,
                        dtype=np.float64, mmap=False, fft_type="rfft")
    s_f = FourierSolver(n=n_odd, L=L, model_params=MODEL_PARAMS,
                        dtype=np.float64, mmap=False, fft_type="fft")
    s_r.precompute()
    s_f.precompute()
    out_r = getattr(s_r, solve_method)()
    out_f = getattr(s_f, solve_method)()
    assert out_r.shape == out_f.shape == (n_odd, n_odd, n_odd, n_odd)
    np.testing.assert_allclose(out_f, out_r, rtol=1e-10, atol=1e-12)
