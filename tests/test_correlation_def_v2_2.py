"""
Tests pinning the three Q-side methods of the v2.2 GcalcForQ (v2.1's
Gaussian-kernel formula additionally multiplied by the scalar prefactor
S(r1) * S(r2)) against each other so any future change to one must be
matched in the others.

Cross-consistency invariants exercised here:
  (1) calc_C_Q_index runs and produces finite output of the expected
      shape/dtype.
  (2) calc_C_Q_index_symmetrized == calc_C_Q_index + calc_C_Q_index.T.
      Like v2.1, c is NOT symmetric under r1<->r2 here -- the `-P_ab*PT_mn`
      cross-terms break that symmetry; the S(r1)*S(r2) prefactor is
      symmetric so it does not restore it.
  (3) calc_C_Q[..., a, b, m, n] == calc_C_Q_index(a, b, m, n) for every
      (a, b, m, n).
"""
import numpy as np
import pytest

from myxo.correlation_def_v2_2 import GcalcForQ


N = 8
L = 5.0
MODEL_PARAMS = dict(S0=1.0, l=1.0, lM=7.0 / 3, lm=0.7 / 3, R=10)


@pytest.fixture(scope="module")
def gcalc_Q():
    grid_1d = np.linspace(-L, L - 2 * L / N, N).astype(np.float64)
    g = GcalcForQ(grid_1d, **MODEL_PARAMS, dtype=np.dtype(np.float64))
    g.precompute()
    return g


def test_S_field_is_populated(gcalc_Q):
    # v2.2's calcQ override must materialise self.S as the scalar field
    # on the (N, N) grid; calc_C_Q_index and friends rely on it.
    assert hasattr(gcalc_Q, "S")
    assert gcalc_Q.S.shape == (N, N)
    assert gcalc_Q.S.dtype == np.float64
    assert np.all(np.isfinite(gcalc_Q.S))


@pytest.mark.parametrize("indices", [
    (0, 0, 0, 0),
    (0, 0, 0, 1),
    (0, 1, 0, 1),
    (0, 1, 1, 0),
    (1, 1, 1, 1),
])
def test_calc_C_Q_index_runs(gcalc_Q, indices):
    c = gcalc_Q.calc_C_Q_index(*indices)
    assert c.shape == (N, N, N, N)
    assert c.dtype == np.float64
    assert np.all(np.isfinite(c))


@pytest.mark.parametrize("indices", [
    (0, 0, 0, 0),
    (0, 0, 0, 1),
    (0, 1, 0, 1),
    (0, 1, 1, 0),
    (1, 1, 1, 1),
])
def test_calc_C_Q_index_symmetrized_matches_two_step(gcalc_Q, indices):
    # The v2.2 calc_C_Q_index is NOT symmetric under r1<->r2 (the
    # `-P_ab*PT_mn` cross-term in both halves breaks the swap symmetry;
    # S(r1)*S(r2) is symmetric and does not restore it), so the symmetrised
    # value must be checked against explicit c + c.transpose(2,3,0,1).
    c = gcalc_Q.calc_C_Q_index(*indices)
    expected = c + c.transpose(2, 3, 0, 1)
    actual = gcalc_Q.calc_C_Q_index_symmetrized(*indices)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)


def test_calc_C_Q_full_matches_index_loop(gcalc_Q):
    full = gcalc_Q.calc_C_Q()
    assert full.shape == (N, N, N, N, 2, 2, 2, 2)
    assert full.dtype == np.float64
    for a in range(2):
        for b in range(2):
            for m in range(2):
                for n in range(2):
                    single = gcalc_Q.calc_C_Q_index(a, b, m, n)
                    np.testing.assert_allclose(
                        full[..., a, b, m, n], single,
                        rtol=1e-12, atol=1e-14,
                        err_msg=f"mismatch at (a,b,m,n)=({a},{b},{m},{n})",
                    )


def test_calc_C_Q_index_out_parameter_writes_in_place(gcalc_Q):
    buf = np.empty((N, N, N, N), dtype=np.float64)
    result = gcalc_Q.calc_C_Q_index(0, 1, 0, 1, out=buf)
    # ne.evaluate with out= returns the buffer; the result must be the
    # same object so callers can reuse a pre-allocated scratch array.
    assert result is buf
    expected = gcalc_Q.calc_C_Q_index(0, 1, 0, 1)
    np.testing.assert_array_equal(buf, expected)


def test_calc_C_Q_index_symmetrized_out_parameter_writes_in_place(gcalc_Q):
    buf = np.empty((N, N, N, N), dtype=np.float64)
    result = gcalc_Q.calc_C_Q_index_symmetrized(0, 1, 0, 1, out=buf)
    assert result is buf
    expected = gcalc_Q.calc_C_Q_index_symmetrized(0, 1, 0, 1)
    np.testing.assert_array_equal(buf, expected)


def test_phi_kernel_matches_expected_gaussian_with_S_prefactor(gcalc_Q):
    # Spot-check at r1 == r2 (dx = dy = 0, so G = 1 and PT collapses to P).
    # For a = b = m = n = 0 the v2.2 index formula reduces to
    #     G * S(r1) * S(r2) * [(3*P_00^2 - P_00^2) + (3*P_00^2 - P_00^2)]
    #   = 1 * S(r)^2          *  4 * P_00^2.
    c = gcalc_Q.calc_C_Q_index(0, 0, 0, 0)
    P00 = gcalc_Q.P[..., 0, 0]
    S = gcalc_Q.S
    expected_diag = (S * S) * 4.0 * (P00 * P00)
    diag = np.einsum("ijij->ij", c)
    np.testing.assert_allclose(diag, expected_diag, rtol=1e-13, atol=1e-15)
