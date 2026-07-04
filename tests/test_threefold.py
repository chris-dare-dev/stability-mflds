"""Algorithms 4 & 5: BG inequality (surface + threefold tilt-BG boundary)."""

import math
from fractions import Fraction as F

from bridgeland_stability.bg_check import check_bg_surface
from bridgeland_stability.chern import ChernChar
from bridgeland_stability.exceptional import Bundle
from bridgeland_stability.threefold import (
    ThreefoldChern,
    alpha_crit,
    bg_boundary_curve,
    bmt_Q_at,
    check_bg_threefold,
)
from bridgeland_stability.varieties import P2, P3, QUINTIC, BLOWUP_P3_POINT


# ---- surface BG (Algorithm 4) --------------------------------------------
def test_bg_surface_exceptional_ok():
    res = check_bg_surface(Bundle.T_minus1(), P2)
    assert res.satisfies
    assert res.discriminant == F(3, 8)
    assert res.discriminant_brief == F(3, 4)
    assert res.bogomolov_integer == 3  # c1^2 - 2 r ch2 = 1 + 2 = 3


def test_bg_surface_equality_case():
    res = check_bg_surface(ChernChar(2, F(0), F(0)), P2)  # O + O
    assert res.satisfies
    assert res.discriminant == 0
    assert res.equality
    assert res.bogomolov_integer == 0


def test_bg_surface_violation():
    res = check_bg_surface(ChernChar(1, F(0), F(1, 4)), P2)
    assert not res.satisfies
    assert res.discriminant == F(-1, 4)  # CH convention
    assert res.discriminant_brief == F(-1, 2)  # brief convention


# ---- threefold tilt-BG (Algorithm 5) -------------------------------------
def test_p3_null_correlation_alpha_crit():
    # CORRECTED values. Brief: beta=1/2 -> sqrt(29)/4 (wrong); beta=1 -> sqrt(2)/2 (wrong).
    v = ThreefoldChern(2, F(0), F(1), F(0))  # (r, ch1.H^2, ch2.H, ch3)
    assert bmt_Q_at(v, F(1, 2), P3.d3) == 3
    assert abs(alpha_crit(v, F(1, 2), P3.d3) - math.sqrt(3)) < 1e-12
    assert bmt_Q_at(v, F(1), P3.d3) == 0
    assert alpha_crit(v, F(1), P3.d3) == 0.0
    # the brief's bad value would require ch3b = -7/6 (dropped rank factor)
    tw = v.twist(F(1), P3.d3)
    assert tw.a3 == F(-4, 3)  # correct; NOT -7/6


def test_p3_null_correlation_degenerate_at_zero():
    v = ThreefoldChern(2, F(0), F(1), F(0))
    # beta=0 gives ch1b = 0 -> degenerate vertical wall
    assert v.twist(F(0), P3.d3).a1 == 0
    assert alpha_crit(v, F(0), P3.d3) is None


def test_quintic_structure_sheaf_on_bg_boundary():
    O = ThreefoldChern(1, F(0), F(0), F(0))
    for beta in (F(1), F(2, 3), F(-3, 7), F(5)):
        assert bmt_Q_at(O, beta, QUINTIC.d3) == 0  # Q vanishes identically
        if beta != 0:
            assert alpha_crit(O, beta, QUINTIC.d3) == 0.0


def test_threefold_bg_check_threshold():
    v = ThreefoldChern(2, F(0), F(1), F(0))
    ok = check_bg_threefold(v, alpha=1, beta=F(1, 2), d3=P3.d3)
    assert ok["Q"] == 3 and ok["threshold"] == 1 and ok["satisfies"]
    bad = check_bg_threefold(v, alpha=2, beta=F(1, 2), d3=P3.d3)
    assert bad["threshold"] == 4 and not bad["satisfies"]  # 3 < 4


def test_bg_boundary_curve_flags_proven_status():
    v = ThreefoldChern(2, F(0), F(1), F(0))
    proven = bg_boundary_curve(v, P3, beta_range=(-2, 2), N=200)
    assert proven.bg_proven is True
    assert len(proven.alphas) == len(proven.betas) > 0
    unproven = bg_boundary_curve(v, BLOWUP_P3_POINT, beta_range=(-2, 2), N=200)
    assert unproven.bg_proven is False
