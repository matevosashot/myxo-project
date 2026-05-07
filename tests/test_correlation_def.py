import numpy as np
import pytest

from myxo.correlation_def import Gcalc, S, f, _mvn_pdf_2d, _mvn_pdf_2d_fast


DTYPE = np.dtype(np.float64)


@pytest.fixture(scope="module")
def gcalc():
    grid = np.linspace(-3.0, 3.0, 6, dtype=DTYPE)
    g = Gcalc(grid, S0=1.0, l=1.0, lM=7.0 / 3, lm=0.7 / 3, R=10.0, dtype=DTYPE)
    g.precompute()
    return g


@pytest.fixture(scope="module")
def small_gcalc():
    grid = np.linspace(-2.0, 2.0, 4, dtype=DTYPE)
    g = Gcalc(grid, S0=0.8, l=1.2, lM=2.0, lm=0.5, R=8.0, dtype=DTYPE)
    g.precompute()
    return g


# ---------- helpers ----------

def test_f_at_zero():
    assert f(np.array(0.0)) == 0.0


def test_S_at_zero():
    assert S(np.array(0.0), S0=1.0, l=1.0, R=10.0) == 0.0


def test_mvn_pdf_fast_matches_reference():
    rng = np.random.default_rng(0)
    dx = rng.normal(size=(3, 4)).astype(np.float64)
    dy = rng.normal(size=(3, 4)).astype(np.float64)
    cov = np.array([[1.5, 0.3], [0.3, 0.9]], dtype=np.float64)
    cov_b = np.broadcast_to(cov, (3, 4, 2, 2))

    ref = _mvn_pdf_2d(dx, dy, cov_b)
    fast = _mvn_pdf_2d_fast(dx, dy, cov_b)
    np.testing.assert_allclose(fast, ref, rtol=1e-12, atol=1e-14)


# ---------- Q tensor ----------

def test_Q_symmetric(gcalc):
    Q = gcalc.Q
    np.testing.assert_allclose(Q[..., 0, 1], Q[..., 1, 0])


def test_Q_traceless(gcalc):
    Q = gcalc.Q
    np.testing.assert_allclose(Q[..., 0, 0] + Q[..., 1, 1], 0.0, atol=1e-14)


def test_Q_shape(gcalc):
    N = len(gcalc.grid)
    assert gcalc.Q.shape == (N, N, 2, 2)


# ---------- P tensor ----------

def test_P_equals_half_I_plus_Q(gcalc):
    eye = np.eye(2)
    expected = 0.5 * (eye + gcalc.Q)
    np.testing.assert_allclose(gcalc.P, expected)


def test_P_trace_is_one(gcalc):
    tr = gcalc.P[..., 0, 0] + gcalc.P[..., 1, 1]
    np.testing.assert_allclose(tr, 1.0)


def test_P_symmetric(gcalc):
    P = gcalc.P
    np.testing.assert_allclose(P[..., 0, 1], P[..., 1, 0])


# ---------- Sigma ----------

def test_Sigma_symmetric(gcalc):
    Sigma = gcalc.Sigma
    np.testing.assert_allclose(Sigma[..., 0, 1], Sigma[..., 1, 0])


def test_Sigma_positive_definite(gcalc):
    Sigma = gcalc.Sigma
    det = Sigma[..., 0, 0] * Sigma[..., 1, 1] - Sigma[..., 0, 1] * Sigma[..., 1, 0]
    assert (Sigma[..., 0, 0] > 0).all()
    assert (det > 0).all()


def test_Sigma_far_from_origin_isotropic():
    # At grid points outside the cutoff R the order parameter S is exponentially
    # damped, so Q ≈ 0 and Sigma collapses to 0.5*(lM^2 + lm^2)*I.
    grid = np.linspace(-50.0, 50.0, 4, dtype=DTYPE)
    g = Gcalc(grid, S0=1.0, l=1.0, lM=2.0, lm=0.5, R=5.0, dtype=DTYPE, sigma=1.0)
    g.precompute()
    expected_diag = 0.5 * (g.lM**2 + g.lm**2)
    # S has a soft Gaussian cutoff exp(-r/R); at r=50, R=5 the residual is
    # ~exp(-10) ≈ 5e-5, so Sigma deviates from the isotropic value at that order.
    np.testing.assert_allclose(g.Sigma[0, 0, 0, 0], expected_diag, atol=1e-3)
    np.testing.assert_allclose(g.Sigma[0, 0, 0, 1], 0.0, atol=1e-3)


# ---------- phi ----------

def test_phi_nonnegative(gcalc):
    assert (gcalc.phi >= 0).all()


def test_phi_transpose(gcalc):
    np.testing.assert_array_equal(gcalc.phi_transpose, gcalc.phi.transpose(2, 3, 0, 1))


def test_phi_shape(gcalc):
    N = len(gcalc.grid)
    assert gcalc.phi.shape == (N, N, N, N)


# ---------- C_P consistency (mirrors selection lines 282-286) ----------

def test_calc_C_P_matches_indexed(gcalc):
    value = gcalc.calc_C_P()
    value_2 = np.empty_like(value)
    for a in range(2):
        for b in range(2):
            value_2[..., a, b] = gcalc.calc_C_P_index(a, b)
    np.testing.assert_allclose(value, value_2, rtol=1e-12, atol=1e-14)


def test_calc_C_P_matches_indexed_alt_params(small_gcalc):
    value = small_gcalc.calc_C_P()
    value_2 = np.empty_like(value)
    for a in range(2):
        for b in range(2):
            value_2[..., a, b] = small_gcalc.calc_C_P_index(a, b)
    np.testing.assert_allclose(value, value_2, rtol=1e-12, atol=1e-14)


def test_C_P_symmetric_in_ab(gcalc):
    # P is symmetric so C_P(a,b) should equal C_P(b,a).
    C = gcalc.calc_C_P()
    np.testing.assert_allclose(C[..., 0, 1], C[..., 1, 0], rtol=1e-12, atol=1e-14)


def test_C_P_shape(gcalc):
    N = len(gcalc.grid)
    assert gcalc.calc_C_P().shape == (N, N, N, N, 2, 2)


# ---------- C_Q consistency (mirrors selection lines 289-298) ----------

def test_calc_C_Q_matches_indexed(gcalc):
    value = gcalc.calc_C_Q()
    value_2 = np.empty_like(value)
    for a in range(2):
        for b in range(2):
            for m in range(2):
                for n in range(2):
                    value_2[..., a, b, m, n] = gcalc.calc_C_Q_index(a, b, m, n)
    np.testing.assert_allclose(value, value_2, rtol=1e-12, atol=1e-14)


def test_calc_C_Q_matches_indexed_alt_params(small_gcalc):
    value = small_gcalc.calc_C_Q()
    value_2 = np.empty_like(value)
    for a in range(2):
        for b in range(2):
            for m in range(2):
                for n in range(2):
                    value_2[..., a, b, m, n] = small_gcalc.calc_C_Q_index(a, b, m, n)
    np.testing.assert_allclose(value, value_2, rtol=1e-12, atol=1e-14)


def test_C_Q_shape(gcalc):
    N = len(gcalc.grid)
    assert gcalc.calc_C_Q().shape == (N, N, N, N, 2, 2, 2, 2)


def test_C_Q_swap_r1_r2_symmetry(gcalc):
    # The construction in calc_C_Q is part_left + part_right where part_right
    # is the (r1<->r2, ab<->mn) swap of part_left. So swapping the spatial
    # indices and the (ab,mn) tensor pairs must leave the result invariant.
    C = gcalc.calc_C_Q()
    swapped = C.transpose(2, 3, 0, 1, 6, 7, 4, 5)
    np.testing.assert_allclose(C, swapped, rtol=1e-12, atol=1e-14)
