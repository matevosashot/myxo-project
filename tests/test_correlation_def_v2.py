"""
Tests pinning the three Q-side methods of GcalcForQ against each other
so any future change to one must be matched in the others.

Cross-consistency invariants exercised here:
  (1) calc_C_Q_index runs and produces finite output of the expected
      shape/dtype.
  (2) c[r1,r2] == c[r2,r1] -- v2's phi and P-product structure are both
      symmetric under r1 <-> r2.
  (3) calc_C_Q_index_symmetrized == calc_C_Q_index + calc_C_Q_index.T,
      which (by 2) reduces to 2 * calc_C_Q_index.
  (4) calc_C_Q[..., a, b, m, n] == calc_C_Q_index(a, b, m, n) for every
      (a, b, m, n).
"""
import numpy as np
import pytest

from myxo.correlation_def_v2 import GcalcForQ


N = 8
L = 5.0
MODEL_PARAMS = dict(S0=1.0, l=1.0, lM=7.0 / 3, lm=0.7 / 3, R=10)


@pytest.fixture(scope="module")
def gcalc_Q():
    grid_1d = np.linspace(-L, L - 2 * L / N, N).astype(np.float64)
    g = GcalcForQ(grid_1d, **MODEL_PARAMS, dtype=np.dtype(np.float64))
    g.precompute()
    return g


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
def test_calc_C_Q_index_is_symmetric_under_r1_r2_swap(gcalc_Q, indices):
    # v2's phi(r1,r2) = 2*exp(-|r2-r1|/lM) is symmetric, and the
    # P-product sum is too (swapping P<->PT only permutes the four terms).
    # The symmetry is analytically exact; the tolerance accounts for
    # ne.evaluate's chunked-order rounding.
    c = gcalc_Q.calc_C_Q_index(*indices)
    np.testing.assert_allclose(c, c.transpose(2, 3, 0, 1), rtol=1e-13, atol=1e-15)

@pytest.mark.parametrize("indices", [
    (0, 0, 0, 0),
    (0, 0, 0, 1),
    (0, 1, 0, 1),
    (0, 1, 1, 0),
    (1, 1, 1, 1),
])
def test_calc_C_Q_index_symmetrized_matches_two_step(gcalc_Q, indices):
    # By (r1,r2) symmetry the symmetrised value is 2*calc_C_Q_index, but
    # the test still pins it against the explicit two-step computation so
    # any future asymmetric override is caught immediately.
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
