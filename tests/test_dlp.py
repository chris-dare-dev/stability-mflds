"""Algorithm 2: the Drezet-Le Potier curve delta(mu) on P^2 (corrected)."""

from fractions import Fraction as F

from bridgeland_stability.dlp import delta, dlp_curve, moduli_nonempty
from bridgeland_stability.exceptional import Bundle, enumerate_exceptional


def _bundles():
    return enumerate_exceptional(-3, 3, R_max=60)


def test_delta_at_integers_is_one():
    bs = _bundles()
    for n in (-1, 0, 1, 2):
        assert delta(n, bs) == 1


def test_delta_canonical_values():
    bs = _bundles()
    # CORRECTED values (CH convention). Brief claimed nu(1/2)=3/4, nu(1/3)=8/9 -- both wrong.
    assert delta(F(1, 2), bs) == F(5, 8)  # cusp of T(-1): 1 - 3/8
    assert delta(F(1, 3), bs) == F(5, 9)  # controlled by O at distance 1/3: P(-1/3)
    assert delta(F(2, 3), bs) == F(5, 9)  # symmetric
    assert delta(F(1, 4), bs) == F(21, 32)  # P(-1/4)
    assert delta(F(2, 5), bs) == F(13, 25)  # cusp of rank-5 bundle: 1 - 12/25


def test_delta_is_at_least_one_half():
    bs = _bundles()
    for k in range(1, 40):
        mu = F(k, 40)
        assert delta(mu, bs) >= F(1, 2)


def test_exceptional_bundles_lie_below_the_curve():
    # An exceptional bundle E_alpha sits at Delta_alpha, strictly below delta(alpha).
    bs = _bundles()
    T = Bundle.T_minus1()
    assert T.discriminant == F(3, 8)
    assert delta(T.slope, bs) == F(5, 8)
    assert T.discriminant < delta(T.slope, bs)


def test_dlp_curve_structure():
    curve = dlp_curve(0, 1, R_max=30)
    # cusp tips lie ON the curve at height 1 - Delta_alpha
    cusp_dict = dict(curve.cusps)
    assert cusp_dict[F(1, 2)] == F(5, 8)
    # the exceptional point itself is below, at Delta_alpha
    exc_dict = dict(curve.exceptional_points)
    assert exc_dict[F(1, 2)] == F(3, 8)
    # sampled curve never dips below 1/2
    assert all(dval >= F(1, 2) for dval in curve.deltas)


def test_moduli_nonempty_exceptional_is_isolated_point():
    res = moduli_nonempty(Bundle.T_minus1())
    assert res["exceptional"] is True
    assert res["nonempty"] is True
    assert res["positive_dimensional"] is False  # 3/8 < delta(1/2)=5/8
    assert "isolated point below" in res["reason"]


def test_moduli_nonempty_above_curve():
    # rank 2, c1=1, ch2=-5/2 -> Delta = 1/8 + 5/4 = 11/8 > delta(1/2)=5/8
    E = Bundle(2, 1, F(-5, 2))
    res = moduli_nonempty(E)
    assert res["discriminant"] == F(11, 8)
    assert res["positive_dimensional"] is True
    assert res["nonempty"] is True


def test_moduli_empty_between_exceptional_and_curve():
    # (4,2,-1) = ch of T(-1)^{+2}: Delta=3/8 < delta(1/2)=5/8, not exceptional -> empty
    E = Bundle(4, 2, F(-1))
    res = moduli_nonempty(E)
    assert res["discriminant"] == F(3, 8)
    assert res["exceptional"] is False
    assert res["positive_dimensional"] is False
    assert res["nonempty"] is False
    assert "EMPTY" in res["reason"]
