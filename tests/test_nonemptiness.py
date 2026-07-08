"""E11-M3 / G18a: the in-architecture ``Delta >= delta_H`` non-emptiness evaluator.

Every asserted numeric value is hand-derived AND cross-checked two ways:

* the ``delta_H`` values regress to the already-pinned Drezet-Le Potier
  ``dlp.delta`` (``delta(1/2)=5/8``, ``delta(1/3)=5/9``, ``delta(2/5)=13/25``,
  ``delta(1/4)=21/32`` -- see ``tests/test_dlp.py`` / ``docs/CORRECTIONS.md``);
* the ``discriminant_H`` values equal the pinned CH-normalized
  ``ChernChar.discriminant`` on ``P^2`` (``d=1``).

Hand derivations (mu = <c1,H>/(r d), d=1 on P^2):
  SurfaceBundle(2,(1,),0)   -> mu=1/2, delta_H=5/8,  Delta=1/2*1/4 - 0      = 1/8
  SurfaceBundle(3,(1,),0)   -> mu=1/3, delta_H=5/9
  SurfaceBundle(5,(2,),0)   -> mu=2/5, delta_H=13/25
  SurfaceBundle(4,(1,),0)   -> mu=1/4, delta_H=21/32
  SurfaceBundle(2,(1,),-5/2)-> mu=1/2, Delta=1/8 - (-5/2)/2 = 1/8 + 5/4     = 11/8  (>= 5/8 -> True)
  SurfaceBundle(3,(2,),1)   -> mu=2/3, Delta=1/2*4/9 - 1/3  = 2/9 - 3/9     = -1/9  (= ChernChar(3,2,1).discriminant(1))
  P1xP1 SurfaceBundle(2,(1,0),0) -> mu=<(1,0),(1,1)>/(2*2)=1/4, Delta=1/2*1/16=1/32; delta_H floor=0 -> True but HEURISTIC.
"""

from fractions import Fraction as F
import subprocess
import sys

import pytest

from bridgeland_stability.nonemptiness_rational import (
    delta_H,
    discriminant_H,
    moduli_nonempty,
    NonemptinessVerdict,
    HNMode,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi
from bridgeland_stability.varieties import P2, P1xP1, enriques
from bridgeland_stability.dlp import delta as dlp_delta
from bridgeland_stability.dlp import moduli_nonempty as dlp_moduli_nonempty
from bridgeland_stability.exceptional import Bundle, enumerate_exceptional
from bridgeland_stability.chern import ChernChar
from bridgeland_stability.rigor import Rigor


# --------------------------------------------------------------------------
# Acceptance (2): P^2-limit / single-ray case regresses to the pinned DLP delta.
# --------------------------------------------------------------------------
def test_p2_limit_regresses_to_dlp():
    # delta_H on P^2 IS the closed-form DLP curve; two-way: pinned literature
    # value + an independent dlp.delta recompute over the same bundle window.
    for r, c, expected in [(2, 1, F(5, 8)), (3, 1, F(5, 9)),
                           (5, 2, F(13, 25)), (4, 1, F(21, 32))]:
        xi = SurfaceBundle(r, (c,), 0)
        mu = F(c, r)  # d = 1 on P^2
        assert delta_H(xi, P2) == expected
        bundles = enumerate_exceptional(mu - 3, mu + 3, 60)
        assert dlp_delta(mu, bundles) == expected  # independent regression


# --------------------------------------------------------------------------
# discriminant_H is the CH-normalized Delta (never the doubled brief); two-way
# against the pinned ChernChar.discriminant on P^2 (d=1).
# --------------------------------------------------------------------------
def test_discriminant_H_matches_chern_convention():
    assert discriminant_H(SurfaceBundle(2, (1,), 0), P2) == F(1, 8)
    assert discriminant_H(SurfaceBundle(2, (1,), F(-5, 2)), P2) == F(11, 8)
    assert discriminant_H(SurfaceBundle(3, (2,), 1), P2) == F(-1, 9)
    # independent recompute via ChernChar (rank-1, d=1)
    assert ChernChar(2, 1, 0).discriminant(1) == F(1, 8)
    assert ChernChar(2, 1, F(-5, 2)).discriminant(1) == F(11, 8)
    assert ChernChar(3, 2, 1).discriminant(1) == F(-1, 9)


# --------------------------------------------------------------------------
# Acceptance (1): moduli_nonempty computes Delta >= delta_H exactly, returns a
# NonemptinessVerdict with a G5 Certificate.
# --------------------------------------------------------------------------
def test_moduli_nonempty_p2_positive():
    v = moduli_nonempty(2, (1,), F(-5, 2), P2)
    assert isinstance(v, NonemptinessVerdict)
    assert v.discriminant == F(11, 8)
    assert v.delta_H == F(5, 8)
    assert v.nonempty is True                 # 11/8 >= 5/8
    assert v.mode is HNMode.DLP
    assert v.certified is True
    assert v.certificate.rigor == Rigor.PROVEN
    # every core quantity stays exact Fraction
    assert isinstance(v.discriminant, F) and isinstance(v.delta_H, F)


def test_moduli_nonempty_p2_below_curve():
    # SurfaceBundle(2,(1,),0): Delta=1/8 < delta_H=5/8 -> the inequality is False.
    # NOT an exceptional bundle (c2=1/2 non-integral, chi(E,E)=3), so the DLP
    # exceptional disjunct does NOT rescue it: it stays EMPTY.
    v = moduli_nonempty(2, (1,), 0, P2)
    assert v.discriminant == F(1, 8)
    assert v.delta_H == F(5, 8)
    assert v.nonempty is False
    assert v.exceptional is False
    assert v.mode is HNMode.DLP
    assert v.certificate.rigor == Rigor.PROVEN


# --------------------------------------------------------------------------
# MUST-FIX (E11-M3 STAGE-4): the DLP exceptional-bundle disjunct.  On P^2 the
# DLP theorem gives non-emptiness iff  Delta >= delta(mu)  OR  xi is exceptional
# (a non-empty ISOLATED POINT that may sit STRICTLY BELOW the curve).  Dropping
# the second disjunct made moduli_nonempty return a PROVEN nonempty=False for a
# class that EXISTS -- contradicting dlp.moduli_nonempty and the theorem.  These
# tests pin the corrected behaviour and its agreement with the sibling.
# --------------------------------------------------------------------------
def test_p2_exceptional_T_minus1_is_nonempty():
    # T_{P^2}(-1) = (2, 1, -1/2): Delta = 1/8 + 1/4 = 3/8 < 5/8 = delta(1/2),
    # yet it is a genuine exceptional bundle -> non-empty isolated point.
    v = moduli_nonempty(2, (1,), F(-1, 2), P2)
    assert v.discriminant == F(3, 8)           # two-way: ChernChar recompute below
    assert ChernChar(2, 1, F(-1, 2)).discriminant(1) == F(3, 8)
    assert v.delta_H == F(5, 8)                 # pinned DLP delta(1/2)
    assert v.discriminant < v.delta_H          # STRICTLY below the curve
    assert v.nonempty is True                   # rescued by the exceptional disjunct
    assert v.exceptional is True
    assert v.mode is HNMode.DLP
    assert v.certified is True
    assert v.certificate.rigor == Rigor.PROVEN  # PROVEN NON-empty, never provably-empty
    assert "exceptional" in v.reason
    # agrees bit-for-bit with the sibling dlp.moduli_nonempty on the same class
    assert dlp_moduli_nonempty(Bundle(2, 1, F(-1, 2)))["nonempty"] is True


def test_p2_exceptional_below_curve_matches_dlp():
    # Every genuine exceptional bundle below the curve must read non-empty AND
    # agree with dlp.moduli_nonempty.  Includes the rank-5 slope-2/5 bundle whose
    # Delta=12/25 sits just below the pinned delta(2/5)=13/25.
    cases = [
        (1, 0, F(0), F(0), F(1)),          # O:        Delta 0   < delta(0)=1
        (2, 1, F(-1, 2), F(3, 8), F(5, 8)),  # T(-1):  Delta 3/8 < delta(1/2)=5/8
        (5, 2, F(-2), F(12, 25), F(13, 25)),  # rk5:    Delta 12/25 < delta(2/5)=13/25
    ]
    for r, c1, ch2, exp_disc, exp_dH in cases:
        v = moduli_nonempty(r, (c1,), ch2, P2)
        assert v.discriminant == exp_disc
        assert v.delta_H == exp_dH
        assert v.discriminant < v.delta_H       # strictly below the DLP curve
        assert v.exceptional is True
        assert v.nonempty is True
        assert v.certificate.rigor == Rigor.PROVEN
        assert dlp_moduli_nonempty(Bundle(r, c1, ch2))["nonempty"] is True


def test_p2_below_curve_nonexceptional_stays_empty():
    # An INTEGRAL, non-exceptional class strictly below the curve is genuinely
    # empty: (2,1,1/2) has Delta=-1/8 < 5/8, integral c2=0, chi(E,E)=5 (not 1).
    # The exceptional disjunct must NOT over-correct -- verdict stays False and
    # agrees with dlp.moduli_nonempty.
    v = moduli_nonempty(2, (1,), F(1, 2), P2)
    assert v.discriminant == F(-1, 8)
    assert v.exceptional is False
    assert v.nonempty is False
    assert v.certificate.rigor == Rigor.PROVEN
    assert dlp_moduli_nonempty(Bundle(2, 1, F(1, 2)))["nonempty"] is False


def test_exceptional_disjunct_is_p2_only():
    # The exceptional disjunct is a P^2 (Pic = Z.H) fact; off P^2 the flag never
    # fires (the sharp delta_H is delegated -- E11-M4/M5), so a HEURISTIC verdict
    # is never silently upgraded via a spurious "exceptional" claim.
    v = moduli_nonempty(2, (1, 0), 0, P1xP1)
    assert v.exceptional is False
    assert v.mode is HNMode.HEURISTIC


# --------------------------------------------------------------------------
# Acceptance (4): THIN slice off P^2 -- Bogomolov floor, HEURISTIC verdict.
# --------------------------------------------------------------------------
def test_p1xp1_heuristic_floor():
    v = moduli_nonempty(2, (1, 0), 0, P1xP1)
    assert v.discriminant == F(1, 32)         # mu=1/4 -> 1/2*(1/4)^2
    assert v.delta_H == F(0)                  # Bogomolov floor (sharp delta_H delegated)
    assert v.nonempty is True                 # 1/32 >= 0
    assert v.mode is HNMode.HEURISTIC
    assert v.certified is False
    assert v.certificate.rigor == Rigor.HEURISTIC


# --------------------------------------------------------------------------
# Acceptance (3): delta_H >= 0 sanity floor + the verdict reports its mode.
# --------------------------------------------------------------------------
def test_delta_h_nonnegative():
    cases = [
        (SurfaceBundle(2, (1,), 0), P2),
        (SurfaceBundle(3, (1,), 0), P2),
        (SurfaceBundle(5, (2,), 0), P2),
        (SurfaceBundle(4, (1,), 0), P2),
        (SurfaceBundle(2, (1, 0), 0), P1xP1),
        (SurfaceBundle(3, (1, 1), 0), P1xP1),
    ]
    for xi, surf in cases:
        assert delta_H(xi, surf) >= 0


def test_verdict_reports_mode():
    vp = moduli_nonempty(2, (1,), 0, P2)
    assert vp.mode is HNMode.DLP and vp.mode.value == "dlp"
    assert vp.certificate.rigor == Rigor.PROVEN
    assert "dlp" in vp.reason

    vh = moduli_nonempty(2, (1, 0), 0, P1xP1)
    assert vh.mode is HNMode.HEURISTIC
    assert vh.certificate.rigor == Rigor.HEURISTIC
    assert "HEURISTIC" in vh.reason


# --------------------------------------------------------------------------
# Certified external HN datum (E11-M4/M5 hook): PROVEN only with a certified
# hn_source; an uncertified source (or none) is rejected.
# --------------------------------------------------------------------------
def test_certified_external_target_is_proven():
    v = moduli_nonempty(2, (1, 0), 0, P1xP1,
                        delta_H_target=F(1, 64), hn_source=HNMode.PAPER)
    assert v.delta_H == F(1, 64)
    assert v.discriminant == F(1, 32)
    assert v.mode is HNMode.PAPER
    assert v.certified is True
    assert v.certificate.rigor == Rigor.PROVEN
    assert v.nonempty is True                 # 1/32 >= 1/64


def test_uncertified_target_rejected():
    with pytest.raises(ValueError):
        moduli_nonempty(2, (1, 0), 0, P1xP1,
                        delta_H_target=F(1, 64), hn_source=HNMode.HEURISTIC)
    with pytest.raises(ValueError):
        moduli_nonempty(2, (1, 0), 0, P1xP1, delta_H_target=F(1, 64))


# --------------------------------------------------------------------------
# G12 guard: torsion-canonical surfaces (Enriques) are refused, not mis-modelled.
# --------------------------------------------------------------------------
def test_torsion_canonical_refused():
    with pytest.raises(NotImplementedError):
        delta_H(SurfaceBundle(2, (1,), 0), enriques())
    with pytest.raises(NotImplementedError):
        moduli_nonempty(2, (1,), 0, enriques())


# --------------------------------------------------------------------------
# E11 M1<->M3 INTERACTION.  M3 consumes M1's SurfaceBundle and decides P^2
# exceptionality via the P^2 Bundle's exceptional.chi, while M1's chi(E,E,P2) is
# an INDEPENDENT Euler form (exceptional_surface.chi built from the RR formula in
# NS coordinates).  The DLP exceptional disjunct of M3 is correct ONLY if these
# two Euler forms agree on the exact classes M3 relies on -- pin that agreement
# so a sign / normalization drift in either implementation fails loudly here.
# --------------------------------------------------------------------------
def test_m1_euler_agrees_with_m3_exceptional_disjunct():
    # For every P^2 class, M1's self-Euler chi(E,E)==1 (exceptional_surface.chi)
    # must agree bit-for-bit with M3's `.exceptional` flag (exceptional.chi path
    # inside _is_p2_exceptional).  Genuine exceptional bundles read chi==1 AND
    # exceptional; a non-exceptional class below the curve reads neither.
    cases = [
        (2, 1, F(-1, 2), True),   # T_{P^2}(-1): exceptional, chi(E,E)=1
        (5, 2, F(-2), True),      # rank-5 slope-2/5 exceptional bundle, chi(E,E)=1
        (1, 0, F(0), True),       # O_{P^2}: line bundle, chi(E,E)=1
        (2, 1, F(1, 2), False),   # integral non-exceptional below curve, chi(E,E)=5
    ]
    for r, c1, ch2, is_exc in cases:
        xi = SurfaceBundle(r, (c1,), ch2)
        m1_says_exceptional = chi(xi, xi, P2) == 1        # M1 Euler self-pairing
        m3_says_exceptional = moduli_nonempty(r, (c1,), ch2, P2).exceptional
        assert m1_says_exceptional is is_exc              # M1 matches ground truth
        assert m3_says_exceptional is is_exc              # M3 matches ground truth
        assert m1_says_exceptional == m3_says_exceptional  # the two modules agree


def test_discriminant_H_offp2_matches_chern_convention():
    # M3's discriminant_H is only two-way cross-checked on P^2 (d=1) elsewhere;
    # pin the CH normalization on a GENUINE d=2 surface (P^1 x P^1, H=(1,1)) too.
    # discriminant_H uses c = <c1,H> and d = H^2 = 2, exactly ChernChar(r,c,e)'s
    # H-numerical convention -- so the two must agree bit-for-bit.
    cases = [
        (SurfaceBundle(2, (1, 0), 0), F(1, 32)),      # c=<(1,0),(1,1)>=1: 1/2*(1/4)^2
        (SurfaceBundle(3, (1, 1), F(1, 2)), F(-1, 36)),  # c=<(1,1),(1,1)>=2
    ]
    for xi, expected in cases:
        c = P1xP1.lattice.pairing(xi.c1, P1xP1.H)     # H-degree, an integer here
        assert c.denominator == 1
        got = discriminant_H(xi, P1xP1)
        assert got == expected                         # pinned hand value
        assert got == ChernChar(xi.r, int(c), xi.ch2).discriminant(P1xP1.d)  # two-way


# --------------------------------------------------------------------------
# Zero-runtime-dependency invariant: nonemptiness_rational is stdlib-only at
# import time (no matplotlib / plotly / numpy / sympy pulled in).
# --------------------------------------------------------------------------
def test_stdlib_only():
    code = (
        "import sys\n"
        "import bridgeland_stability.nonemptiness_rational as m\n"
        "heavy = {'matplotlib', 'plotly', 'numpy', 'sympy'}\n"
        "bad = sorted(n for n in sys.modules if n.split('.')[0] in heavy)\n"
        "assert not bad, bad\n"
        "print('OK')\n"
    )
    out = subprocess.run([sys.executable, "-c", code],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "OK" in out.stdout
