"""
Tests pinning the three Q-side methods of the v3 GcalcForQ against each
other so any future change to one must be matched in the others.

v3 keeps v2.2's S(r1)*S(r2) scalar prefactor but replaces the Gaussian
kernel inside calc_C_Q_index by an exponential exp(-|r2-r1|/lm) and the
simple P-product P_ab*P_mn structure by a wrapped-Gaussian form built
from
    sigma(r)   = sqrt(-2 log S(r))     -- wrapped-Gaussian std-dev
    varphi(r)  = arctan2(y, x)         -- nematic-axis identifier
    G(r1,r2)   = exp(-|r2-r1| / lm)
combined as
    c[r1,r2; a,b,m,n] = factor * factor_T * S(r1) * S(r2)
        * ( T_ab(varphi)   * T_mn(varphi_T)   * cosh(G * sigma * sigma_T)
          + U_ab(varphi)   * U_mn(varphi_T)   * sinh(G * sigma * sigma_T) )
where T_xy = cos if x == y else sin (with a `cossin -> -sin*cos`
sign convention on the sinh term), and
    factor   = -1 if a == 1 and b == 1 else 1
    factor_T = -1 if m == 1 and n == 1 else 1.

Cross-consistency invariants exercised here:
  (1) calc_C_Q_index runs and produces finite output of the expected
      shape/dtype.
  (2) calc_C_Q_index_symmetrized == calc_C_Q_index + calc_C_Q_index.T.
      Like v2.1 / v2.2, c is NOT symmetric under r1<->r2 here -- the
      factor/factor_T (a==1,b==1) sign rule and the asymmetric cossin /
      sincos pairings break the swap symmetry; the S(r1)*S(r2) and
      G*sigma*sigma_T pieces are themselves symmetric so they do not
      restore it.
  (3) calc_C_Q[..., a, b, m, n] == calc_C_Q_index(a, b, m, n) for every
      (a, b, m, n).
"""
import numpy as np
import pytest

from myxo.correlation_def_v3 import GcalcForQ


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
    # v3's calcQ override (inherited from v2.2) must materialise self.S
    # as the scalar nematic-strength field on the (N, N) grid;
    # calc_C_Q_index and friends rely on it.
    assert hasattr(gcalc_Q, "S")
    assert gcalc_Q.S.shape == (N, N)
    assert gcalc_Q.S.dtype == np.float64
    assert np.all(np.isfinite(gcalc_Q.S))


def test_varphi_sigma_field_is_populated(gcalc_Q):
    # varphi_sigma(r) = sqrt(-2 log max(S(r), eps)) is the wrapped-Gaussian
    # standard deviation; v3's calcQ override stashes it on the (N, N) grid
    # under the name `varphi_sigma` (NOT `sigma`, which the parent
    # Gcalc.__init__ already binds to a scalar covariance knob).
    assert hasattr(gcalc_Q, "varphi_sigma")
    sigma = gcalc_Q.varphi_sigma
    assert isinstance(sigma, np.ndarray), (
        "self.varphi_sigma must be the (N, N) wrapped-Gaussian std-dev "
        "grid; the parent Gcalc.__init__ stores a scalar `sigma` "
        "covariance knob under self.sigma, hence the rename.")
    assert sigma.shape == (N, N)
    assert sigma.dtype == np.float64
    assert np.all(np.isfinite(sigma))
    np.testing.assert_allclose(
        sigma, np.sqrt(-2.0 * np.log(np.maximum(gcalc_Q.S, 1e-15))),
        rtol=1e-13, atol=1e-15,
    )


def test_varphi_field_is_populated(gcalc_Q):
    # varphi(r) = arctan2(y, x) is the nematic-axis identifier
    # (varphi = 2*theta); v3's calcQ stashes it on the (N, N) grid so
    # calc_C_Q_index can pick it up without recomputing per call.
    assert hasattr(gcalc_Q, "varphi")
    assert gcalc_Q.varphi.shape == (N, N)
    assert gcalc_Q.varphi.dtype == np.float64
    assert np.all(np.isfinite(gcalc_Q.varphi))
    grid = gcalc_Q.grid
    x = grid[:, None]
    y = grid[None, :]
    np.testing.assert_allclose(
        gcalc_Q.varphi, np.arctan2(y, x),
        rtol=1e-13, atol=1e-15,
    )


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
    # v3's calc_C_Q_index is NOT symmetric under r1<->r2: the
    # factor/factor_T (a==1,b==1) sign rule and the cossin / sincos
    # off-diagonal pairings break the swap symmetry, while the S(r1)*S(r2)
    # and G*sigma*sigma_T pieces are themselves symmetric and do not
    # restore it. So the symmetrised value must be checked against an
    # explicit c + c.transpose(2,3,0,1) -- not against 2*c.
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


def test_phi_kernel_matches_expected_form_at_diagonal(gcalc_Q):
    # Spot-check at r1 == r2: dx = dy = 0, so G = exp(0) = 1 and
    # varphi_T == varphi, sigma_T == sigma. For (a, b, m, n) = (0, 0, 0, 0)
    # we have factor = factor_T = +1 and the expr is "coscos", so the
    # v3 index formula reduces to
    #     S * ST * ( cos(varphi)*cos(varphi_T) * cosh(G * sigma * sigma_T)
    #              + sin(varphi)*sin(varphi_T) * sinh(G * sigma * sigma_T) )
    #   = S(r)^2 * ( cos^2(varphi) * cosh(sigma^2)
    #              + sin^2(varphi) * sinh(sigma^2) ).
    c = gcalc_Q.calc_C_Q_index(0, 0, 0, 0)
    S = gcalc_Q.S
    sigma = gcalc_Q.varphi_sigma
    varphi = gcalc_Q.varphi
    expected_diag = (S * S) * (
        np.cos(varphi) ** 2 * np.cosh(sigma * sigma)
        + np.sin(varphi) ** 2 * np.sinh(sigma * sigma)
    )
    diag = np.einsum("ijij->ij", c)
    np.testing.assert_allclose(diag, expected_diag, rtol=1e-13, atol=1e-15)


def symmetry_test_1(gcalc_Q):
    CQ = gcalc_Q.calc_C_Q()

    pass