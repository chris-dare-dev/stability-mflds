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
    hirzebruch_with_polarization,
    moduli_nonempty,
    NonemptinessVerdict,
    HNMode,
    PaperDeltaHTarget,
    paper_delta_H_targets,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi
from bridgeland_stability.oracle.m2 import find_m2  # skip-guard only (matches test_oracle.py)
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
# E11-M4 / G18b [RESEARCH]: paper-tabulated delta_H targets.  Each verified entry
# of paper_delta_H_targets() is fed in as a certified PAPER target and must
# reproduce the primary source's yes/no verdict with a PROVEN HN-length-one
# hypothesis.  Every numeric value is traceable to arXiv:1907.06739
# (Coskun-Huizenga) or arXiv:1910.14060 (Levine-Zhang) -- see each entry's
# .citation / .note (R1/R3).
# --------------------------------------------------------------------------
def test_paper_verdicts():
    # Acceptance (1)+(2)+(3): moduli_nonempty reproduces the paper yes/no for
    # every tabulated class, returns the exact delta_H target, and carries
    # rigor=PROVEN ONLY because the HN-length-one datum came from the PAPER table.
    targets = paper_delta_H_targets()
    assert len(targets) >= 1
    for e in targets:
        v = moduli_nonempty(e.r, e.c1, e.ch2, e.surface,
                            delta_H_target=e.delta_H, hn_source=HNMode.PAPER)
        assert v.nonempty is e.paper_nonempty            # paper's yes/no reproduced
        assert v.mode is HNMode.PAPER
        assert v.certified is True
        assert v.certificate.rigor == Rigor.PROVEN       # PROVEN via the certified table
        assert v.delta_H == e.delta_H                     # exact tabulated target returned
        assert isinstance(v.delta_H, F)


def test_paper_table_is_exact_fractions():
    # Exact-arithmetic invariant: every tabulated numeric is a Fraction and every
    # entry cites an arXiv:19xx primary source.
    targets = paper_delta_H_targets()
    assert len(targets) > 0
    for e in targets:
        assert isinstance(e, PaperDeltaHTarget)
        assert isinstance(e.delta_H, F)
        assert isinstance(e.delta_H_paper, F)
        assert isinstance(e.ch2, F)
        assert "arXiv:19" in e.citation


def test_paper_p2_targets_match_native_dlp():
    # Independent two-way anchor: on P^2 the PAPER target must equal the package's
    # OWN native DLP delta-curve, and discriminant_H must equal the pinned
    # ChernChar.discriminant (d=1).  Ties each paper datum to in-package machinery.
    seen = 0
    for e in paper_delta_H_targets():
        if e.surface is not P2:
            continue
        seen += 1
        xi = SurfaceBundle(e.r, e.c1, e.ch2)
        assert delta_H(xi, P2) == e.delta_H                       # native DLP == PAPER target
        assert discriminant_H(xi, P2) == \
            ChernChar(e.r, int(e.c1[0]), e.ch2).discriminant(1)   # two-way discriminant
    assert seen == 4


def test_paper_p1xp1_normalization_conversion():
    # The load-bearing normalization conversion (R3): for the c1||H entries on F_0
    # the paper's full-NS Delta = 1/2 nu^2 - ch2/r equals d*discriminant_H, and the
    # CH-package target is delta_H_paper/d.  The exceptional disjunct is P^2-only.
    seen = 0
    for e in paper_delta_H_targets():
        if e.surface is not P1xP1:
            continue
        seen += 1
        nu = tuple(F(x, e.r) for x in e.c1)                        # nu = c1/r, exact
        Delta_paper = F(1, 2) * P1xP1.lattice.self_pairing(nu) - e.ch2 / e.r
        xi = SurfaceBundle(e.r, e.c1, e.ch2)
        assert Delta_paper == P1xP1.d * discriminant_H(xi, P1xP1)  # c1||H identity
        assert e.delta_H == e.delta_H_paper / P1xP1.d              # /d conversion
        assert e.delta_H_paper == 1                                # controlling l.b.: chi(O_{F_0})
        v = moduli_nonempty(e.r, e.c1, e.ch2, P1xP1,
                            delta_H_target=e.delta_H, hn_source=HNMode.PAPER)
        assert v.exceptional is False                             # disjunct is P^2-only
    assert seen == 2


def test_paper_exceptional_coexists_with_target():
    # E4: the rank-5 exceptional bundle sits STRICTLY BELOW its paper delta_H yet is
    # nonempty via the exceptional disjunct, chi(E,E)=1, and agrees with dlp.
    e = next(t for t in paper_delta_H_targets()
             if t.surface is P2 and t.r == 5 and t.c1 == (2,))
    xi = SurfaceBundle(e.r, e.c1, e.ch2)
    v = moduli_nonempty(e.r, e.c1, e.ch2, P2,
                        delta_H_target=e.delta_H, hn_source=HNMode.PAPER)
    assert discriminant_H(xi, P2) < e.delta_H                     # strictly below the curve
    assert v.nonempty is True
    assert v.exceptional is True
    assert chi(xi, xi, P2) == 1                                   # genuine exceptional bundle
    assert dlp_moduli_nonempty(Bundle(5, 2, e.ch2))["nonempty"] is True


# --------------------------------------------------------------------------
# E11-M5 / G18c [PROVEN]: F_n polarization-dependence witness.
#
# HONEST SCOPING (rules H1/H2).  delta_H(xi, surface) is polarization-INDEPENDENT
# off P^2: it returns the Bogomolov FLOOR 0 for every non-P^2 surface regardless
# of H (verified below: delta_H == 0 for BOTH polarizations).  A literal
# "delta_H differs between two ample H" witness would therefore be FALSE (0 == 0);
# the sharp polarization-dependent delta_H envelope on F_n (the Hirzebruch
# exceptional-bundle DLP limiting procedure) is the DEFERRED G18 remainder and is
# NOT computed here.  So the [PROVEN] witness is at the level of the quantities
# the package GENUINELY computes on F_n for a FIXED class under two DIFFERENT
# ample H: mu = <c1,H>/(r d), d = H^2, discriminant_H, and hence the
# moduli_nonempty VERDICT -- all of which straddle and flip.  This is a real,
# exact, machine-checkable polarization dependence that the P^2 model
# (Pic(P^2) = Z.H, a single ample ray) STRUCTURALLY cannot express.
# --------------------------------------------------------------------------
def test_fn_polarization_dependence():
    # Fixed F_1 class with NON-diagonal c1=(1,1) (not proportional to either ample
    # H); a genuine effective integral class: c1^2 = <(1,1),(1,1)> = 2*1*1 - 1 = 1,
    # so c2 = c1^2/2 - ch2 = 1/2 - 1/2 = 0 (realized e.g. by O + O(f+s)).  Two
    # STRICTLY-ample H in different regions of the F_1 ample cone; the factory
    # H=(1,1) is only nef-and-big (a-nb=0), so we build our own strictly-ample H.
    HA = hirzebruch_with_polarization(1, (2, 1))   # d = <(2,1),(2,1)> = 2*2 - 1 = 3
    HB = hirzebruch_with_polarization(1, (4, 1))   # d = <(4,1),(4,1)> = 2*4 - 1 = 7
    assert HA.d == 3 and HB.d == 7                  # d = H^2 genuinely differs
    xi = SurfaceBundle(2, (1, 1), F(1, 2))
    assert xi.c1 == (F(1), F(1))                    # c2 = <c1,c1>/2 - ch2 = 1/2 - 1/2 = 0

    # mu = <c1,H>/(r d) differs between the two ample H (exact Fractions):
    muA = HA.lattice.pairing(xi.c1, HA.H) / (xi.r * HA.d)   # <(1,1),(2,1)>/(2*3)=2/6
    muB = HB.lattice.pairing(xi.c1, HB.H) / (xi.r * HB.d)   # <(1,1),(4,1)>/(2*7)=4/14
    assert muA == F(1, 3) and muB == F(2, 7) and muA != muB

    # discriminant_H = 1/2 mu^2 - ch2/(r d) differs AND straddles 0:
    dA = discriminant_H(xi, HA)   # 1/2*(1/3)^2 - (1/2)/6  = 1/18 - 1/12 = -1/36
    dB = discriminant_H(xi, HB)   # 1/2*(2/7)^2 - (1/2)/14 = 2/49  - 1/28 =  1/196
    assert dA == F(-1, 36) and dB == F(1, 196) and dA < 0 < dB

    # HONEST SCOPING (rules H1/H2): delta_H is the Bogomolov FLOOR 0 for BOTH H --
    # polarization-INDEPENDENT off P^2.  The sharp polarization-dependent delta_H
    # envelope on F_n (Hirzebruch exceptional-bundle DLP) is the DEFERRED G18
    # remainder, NOT computed here.  We assert EQUALITY (never fake a difference):
    assert delta_H(xi, HA) == 0 and delta_H(xi, HB) == 0

    # ...therefore the moduli_nonempty VERDICT genuinely differs (empty vs nonempty)
    # -- a polarization dependence P^2 (Pic = Z.H, one ray) structurally cannot show.
    vA = moduli_nonempty(xi.r, xi.c1, xi.ch2, HA)
    vB = moduli_nonempty(xi.r, xi.c1, xi.ch2, HB)
    assert vA.nonempty is False and vB.nonempty is True     # verdict flips with H
    assert vA.mode is HNMode.HEURISTIC and vB.mode is HNMode.HEURISTIC
    assert vA.discriminant == F(-1, 36) and vB.discriminant == F(1, 196)
    assert vA.delta_H == 0 and vB.delta_H == 0              # floor unchanged both ways


# --------------------------------------------------------------------------
# E11-M5 / G16 optional cross-check.  oracle.moduli_nonempty_by_construction is
# P^2-only and M2-gated; on this host M2 is ABSENT, so this test SKIPS cleanly
# (exactly like the E10 @requires_m2 tests in tests/test_oracle.py).  When a user
# provisions M2 it asserts the sufficient-only constructive witness AGREES with
# the formula-layer moduli_nonempty on a shared P^2 class list.  The oracle import
# is kept LOCAL to the test body, preserving the core-never-imports-oracle
# invariant (nonemptiness_rational never imports oracle at module top level).
# --------------------------------------------------------------------------
@pytest.mark.skipif(find_m2() is None, reason="Macaulay2 (M2) not installed")
def test_m2_cross_check():
    from bridgeland_stability.oracle.m2 import moduli_nonempty_by_construction
    # Shared P^2 class list inside the oracle's rank-1 ideal-sheaf witness family.
    # moduli_nonempty_by_construction is SUFFICIENT-ONLY: returns True or None,
    # NEVER False.  Agreement = "when the oracle CONSTRUCTS a witness (True), the
    # formula layer also says nonempty".
    shared = [(1, 0, F(-2)),   # I_Z, l=2    -> formula disc 2   >= delta(0)=1  -> True
              (1, 0, F(0)),    # O,   l=0    -> exceptional bundle             -> True
              (1, 2, F(1))]    # I_Z(2), l=1 -> disc 1 >= delta(2)=1 (boundary) -> True
    for r, c1, ch2 in shared:
        oracle = moduli_nonempty_by_construction(r, c1, ch2, P2)
        formula = moduli_nonempty(r, (c1,), ch2, P2)
        assert oracle is not False               # witness never claims emptiness
        if oracle is True:                       # only assert on constructed classes
            assert formula.nonempty is True      # oracle AGREES with formula layer


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
