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
    discriminant,
    discriminant_H,
    hirzebruch_with_polarization,
    moduli_nonempty,
    NonemptinessVerdict,
    HNMode,
    PaperDeltaHTarget,
    paper_delta_H_targets,
)
from bridgeland_stability.dlp_hirzebruch import exceptional_discriminant, is_stable_exceptional
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi
from bridgeland_stability.oracle.m2 import find_m2  # skip-guard only (matches test_oracle.py)
from bridgeland_stability.varieties import P2, P1xP1, enriques, hirzebruch
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


def test_exceptional_disjunct_on_F0_needs_a_genuine_exceptional_bundle():
    # (2,(1,0)) is NOT potentially exceptional on F_0: Lemma "excFacts" (2) forbids even
    # rank when e is even.  So the disjunct must not fire, and the verdict is decided by
    # the envelope alone -- never silently upgraded via a spurious "exceptional" claim.
    v = moduli_nonempty(2, (1, 0), 0, P1xP1)
    assert v.exceptional is False
    assert v.mode is HNMode.DLP_ANTICANONICAL     # F_0 ships H = f+s, the -K ray


# --------------------------------------------------------------------------
# E11-M6 / G18: off P^2 the Bogomolov floor is REPLACED by the natively computed
# Coskun-Huizenga envelope.  This test previously pinned the deferred behaviour
# (delta_H == 0, verdict "nonempty" with HEURISTIC rigor); the envelope now proves
# the opposite, and the old verdict was wrong.
# --------------------------------------------------------------------------
def test_p1xp1_envelope_proves_emptiness():
    # xi = (2, f, 0) has full-NS Delta = 1/2 <nu,nu> - ch2/r = 0 (nu = f/2, <f,f> = 0),
    # NOT the H-projected 1/32.  The line bundle O is a mu_H-stable exceptional bundle
    # with w = nu - 0 = (1/2, 0), w.H = 1/2 (strictly inside the strip |w.H| <= 2), so any
    # mu_H-semistable sheaf of this slope obeys Delta >= P(-w) - 0 = 1 + 1/2(0 - 1) = 1/2.
    # Since 0 < 1/2 there is no such sheaf: M_H(xi) is PROVABLY EMPTY.
    v = moduli_nonempty(2, (1, 0), 0, P1xP1)
    assert v.discriminant == F(0)                 # full-NS Delta (polarization-independent)
    assert discriminant_H(SurfaceBundle(2, (1, 0), 0), P1xP1) == F(1, 32)   # the old surrogate
    assert v.delta_H == F(1, 2)                   # DLP_{-K}(1/2, 0), floored by Cor. "K1/2"
    assert v.nonempty is False
    assert v.mode is HNMode.DLP_ANTICANONICAL
    assert v.certificate.rigor == Rigor.PROVEN    # certified emptiness, not a floor guess
    assert "PROVEN empty" in v.reason


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

    # F_0 with H = f+s is the anticanonical ray -> the sharp envelope mode.
    vk = moduli_nonempty(2, (1, 0), 0, P1xP1)
    assert vk.mode is HNMode.DLP_ANTICANONICAL and vk.mode.value == "dlp_anticanonical"

    # An F_n surface whose only polarization is the nef-and-big factory H has no
    # envelope at all: the Bogomolov floor + HEURISTIC path is still reachable.
    vh = moduli_nonempty(2, (1, 0), 0, hirzebruch(3))
    assert vh.mode is HNMode.HEURISTIC
    assert vh.certificate.rigor == Rigor.HEURISTIC
    assert vh.delta_H == 0
    assert "HEURISTIC" in vh.reason


# --------------------------------------------------------------------------
# Certified external HN datum (E11-M4/M5 hook): PROVEN only with a certified
# hn_source; an uncertified source (or none) is rejected.
# --------------------------------------------------------------------------
def test_certified_external_target_is_proven():
    # An externally supplied target overrides the native envelope entirely (that is the
    # E11-M4 paper hook / E11-M5 oracle hook).  Compared against the full-NS Delta = 2.
    v = moduli_nonempty(2, (0, 0), F(-4), P1xP1,
                        delta_H_target=F(1), hn_source=HNMode.PAPER)
    assert v.delta_H == F(1)
    assert v.discriminant == F(2)
    assert v.mode is HNMode.PAPER
    assert v.certified is True
    assert v.certificate.rigor == Rigor.PROVEN
    assert v.nonempty is True                 # 2 >= 1

    # ...and the target really does override: an absurd target flips the verdict.
    w = moduli_nonempty(2, (0, 0), F(-4), P1xP1,
                        delta_H_target=F(5), hn_source=HNMode.PAPER)
    assert w.nonempty is False and w.delta_H == F(5)


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


def test_paper_p1xp1_targets_match_native_envelope():
    # E11-M6 supersedes the M4 "/d" conversion: moduli_nonempty now compares against the
    # paper's own full-NS Delta, so delta_H == delta_H_paper.  The c1||H identity
    # Delta = d * discriminant_H is retained as a documented relation (it is exactly why
    # the M4 rescaling worked for these two entries and nowhere else).  Crucially the
    # target is no longer a transcription: dlp_envelope computes it natively.
    from bridgeland_stability.dlp_hirzebruch import dlp_envelope, total_slope

    seen = 0
    for e in paper_delta_H_targets():
        if e.surface is not P1xP1:
            continue
        seen += 1
        xi = SurfaceBundle(e.r, e.c1, e.ch2)
        nu = total_slope(xi)
        Delta = F(1, 2) * P1xP1.lattice.self_pairing(nu) - e.ch2 / e.r
        assert discriminant(xi, P1xP1) == Delta                    # full-NS, exact
        assert Delta == P1xP1.d * discriminant_H(xi, P1xP1)        # c1||H identity (relic)
        assert e.delta_H == e.delta_H_paper == 1                   # controlling l.b.: chi(O_{F_0})

        env = dlp_envelope(nu, P1xP1, 8)                           # NATIVE computation
        assert env.value == e.delta_H_paper                        # regression vs the paper
        assert env.exact and env.sharp                             # certified sharp truncation

        v = moduli_nonempty(e.r, e.c1, e.ch2, P1xP1,
                            delta_H_target=e.delta_H, hn_source=HNMode.PAPER)
        assert v.nonempty is e.paper_nonempty
        assert v.exceptional is False                              # not an exceptional bundle
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
# E11-M5 / G18c: F_n polarization-dependence witness -- REDONE for E11-M6 / G18.
#
# The original M5 witness was built on ``discriminant_H``, which depends on H by
# construction (it is 1/2 mu_H^2 - ch2/(r d)).  Under the primary sources' actual
# discriminant Delta = 1/2 <nu,nu> - ch2/r that class, xi = (2,(1,1),1/2), has
# Delta = -1/8 < 0: it violates Bogomolov and is empty for EVERY polarization.  Its
# apparent polarization dependence was an artefact of the H-projected surrogate.
#
# Delta is polarization-independent; the polarization dependence lives entirely in
# delta_H.  The honest witness is therefore a class whose delta_H -- and hence whose
# non-emptiness -- genuinely moves with H, both verdicts PROVEN.  See
# docs/CORRECTIONS.md (C7) and tests/test_dlp_hirzebruch.py.
# --------------------------------------------------------------------------
def test_old_m5_class_violates_bogomolov_under_the_true_discriminant():
    HA = hirzebruch_with_polarization(1, (2, 1))   # d = <(2,1),(2,1)> = 2*2 - 1 = 3
    HB = hirzebruch_with_polarization(1, (4, 1))   # d = <(4,1),(4,1)> = 2*4 - 1 = 7
    assert HA.d == 3 and HB.d == 7
    xi = SurfaceBundle(2, (1, 1), F(1, 2))

    # The H-projected surrogate straddles 0 and moves with H (this is the artefact):
    assert discriminant_H(xi, HA) == F(-1, 36) and discriminant_H(xi, HB) == F(1, 196)

    # The true discriminant is the SAME for both, and negative: empty for every H.
    assert discriminant(xi, HA) == discriminant(xi, HB) == F(-1, 8)
    for S in (HA, HB):
        v = moduli_nonempty(xi.r, xi.c1, xi.ch2, S)
        assert v.discriminant == F(-1, 8)
        assert v.nonempty is False
        assert v.certificate.rigor == Rigor.PROVEN          # Bogomolov, not a floor guess
        assert "Bogomolov" in v.certificate.hypotheses[0]


def test_fn_polarization_dependence():
    # THE witness.  xi = (2, f+s, -1/2) on F_1 is the rank-2 exceptional character:
    # Delta = 3/8 = 1/2 - 1/(2*2^2), c1^2 = <(1,1),(1,1)> = 2 - 1 = 1, c2 = 1/2 + 1/2 = 1.
    # Paper Table 2 gives its stability interval I_V = (0, 1) in the H_m parametrization
    # H_m = E + (e+m)F, which in the package (f, s) basis is the ray of (e+m, 1).
    #
    #   * H = (3,2)  <-> m = 3/2 - 1 = 1/2, the ANTICANONICAL polarization (-K = 3f+2s),
    #     and 1/2 in I_V: V is a mu_H-stable exceptional bundle, so M_H(xi) is a single
    #     reduced point -- NON-EMPTY -- sitting strictly BELOW its own envelope
    #     delta_H = 5/8 (exactly as T_{P^2}(-1) sits below the P^2 DLP curve).
    #   * H = (3,1)  <-> m = 3 - 1 = 2, OUTSIDE I_V: no mu_H-stable exceptional bundle of
    #     this character exists.  The line bundle O(f) then certifies emptiness:
    #     w = nu - nu_{O(f)} = (-1/2, 1/2), w.H = 1/2 > 0 (strict branch, inside the strip
    #     since -1/2 K.H = 7/2), so any mu_H-semistable sheaf obeys
    #     Delta >= P(-w) = 1 + 1/2((-w)^2 - (-w).K) = 1 + 1/2(-3/4 + 1/2) = 7/8 > 3/8.
    #
    # Delta is identical for both; delta_H is not; both verdicts are PROVEN.
    xi = SurfaceBundle(2, (1, 1), F(-1, 2))
    HK = hirzebruch_with_polarization(1, (3, 2))    # anticanonical, d = 2*6 - 4 = 8
    HA = hirzebruch_with_polarization(1, (3, 1))    # m = 2,          d = 2*3 - 1 = 5
    assert HK.d == 8 and HA.d == 5

    # Delta does NOT depend on the polarization:
    assert discriminant(xi, HK) == discriminant(xi, HA) == F(3, 8) == exceptional_discriminant(2)

    # ...but the exceptional bundle exists only for the anticanonical H:
    assert is_stable_exceptional(2, (1, 1), HK) is True
    assert is_stable_exceptional(2, (1, 1), HA) is False

    vK = moduli_nonempty(xi.r, xi.c1, xi.ch2, HK)
    vA = moduli_nonempty(xi.r, xi.c1, xi.ch2, HA)

    # delta_H genuinely differs -- the real polarization dependence:
    assert vK.delta_H == F(5, 8) and vA.delta_H == F(7, 8)
    assert vK.delta_H != vA.delta_H

    # ...and so does the verdict, with BOTH sides PROVEN:
    assert vK.nonempty is True and vK.exceptional is True
    assert vK.mode is HNMode.DLP_ANTICANONICAL
    assert vK.certificate.rigor == Rigor.PROVEN
    assert vK.discriminant < vK.delta_H              # strictly below its own envelope

    assert vA.nonempty is False and vA.exceptional is False
    assert vA.mode is HNMode.DLP_LOWER
    assert vA.certificate.rigor == Rigor.PROVEN
    assert "PROVEN empty" in vA.reason

    # delta_H is sharp only on the anticanonical del Pezzo ray; off it, a lower bound.
    assert delta_H(xi, HK) == F(5, 8) and delta_H(xi, HA) == F(7, 8)


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
