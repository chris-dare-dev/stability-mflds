"""Algorithm 4: dedicated Bogomolov-Gieseker (surface) tests, plus the
threefold tilt-BG re-export.  E1-M4 (test-only; no core formula changes).

Surface BG uses the CH-normalized discriminant  Delta = (1/2) mu^2 - ch2/(r d);
BG (Bogomolov 1978, Gieseker 1977) requires  Delta >= 0  for a mu-semistable
sheaf (a necessary condition).  The threefold tilt-BG form
Q = 4 ch2b^2 - 6 ch1b ch3b  lives in bridgeland_stability.threefold and is
RE-EXPORTED from bg_check as check_bg_threefold; this module pins that
re-export and its sign agreement with threefold.bmt_Q on the P^3
null-correlation class (2,0,1,0).  All asserted values are exact Fractions,
re-derived by hand (see docs/CORRECTIONS.md sections 4 and 5).
"""

from fractions import Fraction as F

import pytest

from bridgeland_stability import bg_check
from bridgeland_stability import threefold as threefold_mod
from bridgeland_stability.bg_check import (
    check_bg_surface,
    check_bg_threefold,
    check_existence_k3,
    check_existence_abelian,
    delta_min,
    ExistenceResult,
)
from bridgeland_stability.chern import ChernChar
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.threefold import ThreefoldChern, bmt_Q, bmt_Q_at
from bridgeland_stability.varieties import P1xP1, K3, abelian_surface


def _sign(x):
    return (x > 0) - (x < 0)


def test_bg_surface_boundary():
    # Only surface.d (= H^2 = 2 here) is used by check_bg_surface -- the
    # rank-independent Bogomolov baseline.  P^1 x P^1 (a rational surface with
    # H^2 = 2) supplies d = 2; no K3 / existence semantics are involved (that is
    # the distinct later milestone E3-M1, check_existence_k3).
    d2 = P1xP1

    # ACCEPT: ChernChar(2,0,-3), mu = 0, Delta = 0 - (-3)/(2*2) = 3/4 >= 0.
    acc = check_bg_surface(ChernChar(2, 0, F(-3)), d2)
    assert acc.satisfies is True
    assert acc.discriminant == F(3, 4)            # CH convention
    assert acc.discriminant_brief == F(3, 2)      # doubled = 2 * Delta
    assert acc.bogomolov_integer == 12            # c^2/d - 2 r e = 0 - 2*2*(-3)
    assert acc.equality is False
    assert acc.discriminant_brief == 2 * acc.discriminant

    # REJECT: sign-flipped ChernChar(2,0,3), Delta = -3/4 < 0.
    rej = check_bg_surface(ChernChar(2, 0, F(3)), d2)
    assert rej.satisfies is False
    assert rej.discriminant == F(-3, 4)
    assert rej.discriminant_brief == F(-3, 2)
    assert rej.bogomolov_integer == -12
    assert rej.discriminant_brief == 2 * rej.discriminant


def test_bg_threefold_reexport():
    # check_bg_threefold is the SAME object re-exported from threefold.py.
    assert bg_check.check_bg_threefold is threefold_mod.check_bg_threefold

    # P^3 null-correlation bundle (r, ch1.H^2, ch2.H, ch3) = (2,0,1,0), d3 = 1.
    # Q^beta = 4(1 - beta^2) exactly, so Q > 0 for |beta| < 1, = 0 at beta = 1,
    # < 0 for |beta| > 1.  The dict Q from check_bg_threefold must equal bmt_Q
    # of the beta-twist (and bmt_Q_at) and agree with it in sign.
    v = ThreefoldChern(2, F(0), F(1), F(0))
    d3 = 1
    for beta, expected_Q in ((F(1, 2), 3), (F(1), 0), (F(2), -12)):
        res = check_bg_threefold(v, alpha=1, beta=beta, d3=d3)
        q_direct = bmt_Q(v.twist(beta, d3))
        assert res["Q"] == expected_Q
        assert q_direct == expected_Q
        assert res["Q"] == bmt_Q_at(v, beta, d3)
        assert _sign(res["Q"]) == _sign(q_direct)


# --------------------------------------------------------------------------
# E3-M1 / G6: K3 and abelian moduli non-emptiness with integrality guard.
# Every asserted numeric value is re-derived exactly (Fraction) and confirmed;
# citations: Yoshioka / Bayer-Macri (arXiv:1301.6968), Yanagida-Yoshioka
# (arXiv:1203.0884).  These are additions only; the two tests above are unchanged.
# --------------------------------------------------------------------------


def test_k3_structure_sheaf():
    # v=(1,0,0) on K3(2): Mukai (1,0,1), Delta=0, delta_min=0, v^2=-2 -> ACCEPT (a point).
    res = check_existence_k3(ChernChar(1, 0, 0), K3(2))
    assert res.mukai == (1, 0, 1)
    assert res.discriminant == F(0)
    assert res.delta_min == F(0)          # (1 - 1/1^2)/2 = 0
    assert res.v_squared == -2
    assert res.moduli_dim == 0
    assert res.nonempty is True
    assert res.satisfies_bound is True
    assert res.primitive is True


def test_k3_reject_rank2():
    # v=(2,2,0) on K3(2): Delta=1/8 < 3/8=delta_min, v^2=-6 < -2 -> REJECT.
    res = check_existence_k3(ChernChar(2, 2, 0), K3(2))
    assert res.discriminant == F(1, 8)
    assert res.delta_min == F(3, 8)       # (1 - 1/4)/2
    assert res.discriminant < res.delta_min
    assert res.mukai == (2, 1, 2)
    assert res.v_squared == -6
    assert res.nonempty is False


def test_k3_accept_rank2():
    # v=(2,0,-3) on K3(2): Delta=3/4 >= 3/8, v^2=4 -> ACCEPT, moduli_dim=6.
    res = check_existence_k3(ChernChar(2, 0, -3), K3(2))
    assert res.discriminant == F(3, 4)
    assert res.delta_min == F(3, 8)
    assert res.mukai == (2, 0, -1)
    assert res.v_squared == 4
    assert res.moduli_dim == 6            # v^2 + 2
    assert res.nonempty is True


def test_abelian_boundary():
    # v=(1,0,0) on abelian_surface(2): BARE Mukai (1,0,0), v^2=0 -> ACCEPT (no +r shift).
    res = check_existence_abelian(ChernChar(1, 0, 0), abelian_surface(2))
    assert res.mukai == (1, 0, 0)         # bare, NOT (1,0,1)
    assert res.v_squared == 0
    assert res.delta_min == F(0)
    assert res.nonempty is True
    assert res.primitive is True


def test_integrality_guard_k3():
    # c=1, d=2: 1 % 2 != 0 -> must RAISE (not silently floor to l=0 and return v^2=-8).
    with pytest.raises(ValueError):
        check_existence_k3(ChernChar(2, 1, 0), K3(2))


def test_integrality_guard_abelian():
    # c=1, d=2: l=1/2 not in Z -> RAISE (guard lives in both checks).
    with pytest.raises(ValueError):
        check_existence_abelian(ChernChar(1, 1, 0), abelian_surface(2))


def test_delta_min_helper():
    assert delta_min(1, 2) == F(0)        # (1-1)/2
    assert delta_min(2, 2) == F(3, 8)     # (4-1)/(4*2)
    assert delta_min(2, 1) == F(3, 4)
    assert delta_min(3, 1) == F(8, 9)
    assert delta_min(5, 1) == F(24, 25)


def test_k3_bound_equivalence():
    # The two formulations agree exactly: (Delta >= delta_min) iff (v^2 >= -2).
    for ch in (ChernChar(1, 0, 0), ChernChar(2, 2, 0), ChernChar(2, 0, -3)):
        res = check_existence_k3(ch, K3(2))
        assert (res.discriminant >= res.delta_min) == (res.v_squared >= -2)
        assert res.nonempty == (res.v_squared >= -2)


def test_imprimitive_not_asserted():
    # v=(2,0,-2) -> Mukai (2,0,0) = 2*(1,0,0), imprimitive. v^2=0 satisfies the bound
    # but Yoshioka's biconditional needs a primitive vector -> nonempty NOT asserted.
    res = check_existence_k3(ChernChar(2, 0, -2), K3(2))
    assert res.mukai == (2, 0, 0)
    assert res.v_squared == 0
    assert res.satisfies_bound is True
    assert res.primitive is False
    assert res.nonempty is None
    assert "IMPRIMITIVE" in res.note


def test_check_bg_surface_untouched():
    # Regression: the rank-independent Bogomolov baseline still ACCEPTS (2,2,0) (Delta=1/8>=0)
    # even though the sharper Yoshioka existence check REJECTS it (Delta=1/8 < 3/8). The two
    # are distinct and check_bg_surface is unchanged.
    bg = check_bg_surface(ChernChar(2, 2, 0), K3(2))
    assert bg.discriminant == F(1, 8)
    assert bg.satisfies is True
    assert check_existence_k3(ChernChar(2, 2, 0), K3(2)).nonempty is False


def test_existence_certificate_and_kind_guards():
    # Certificates: PROVEN inequality, primitivity/genericity hypotheses, correct citations.
    k3 = check_existence_k3(ChernChar(2, 0, -3), K3(2))
    assert k3.certificate.rigor == Rigor.PROVEN
    assert any("primitive" in h for h in k3.certificate.hypotheses)
    assert any("generic" in h for h in k3.certificate.hypotheses)
    assert "arXiv:1301.6968" in k3.certificate.citations
    ab = check_existence_abelian(ChernChar(1, 0, 0), abelian_surface(2))
    assert ab.certificate.rigor == Rigor.PROVEN
    assert "arXiv:1203.0884" in ab.certificate.citations
    # kind guards: K3 shift must not touch abelian input and vice versa.
    with pytest.raises(ValueError):
        check_existence_k3(ChernChar(1, 0, 0), abelian_surface(2))
    with pytest.raises(ValueError):
        check_existence_abelian(ChernChar(1, 0, 0), K3(2))
