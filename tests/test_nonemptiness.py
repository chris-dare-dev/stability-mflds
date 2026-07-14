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
    SharpBoundEvidence,
    HNMode,
    PaperDeltaHTarget,
    paper_delta_H_targets,
)
from bridgeland_stability.dlp_hirzebruch import (
    dlp_envelope,
    emptiness_bound,
    exceptional_ch2,
    exceptional_discriminant,
    is_semiexceptional,
    is_stable_exceptional,
    total_slope,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi
from bridgeland_stability.oracle.m2 import find_m2  # skip-guard only (matches test_oracle.py)
from bridgeland_stability.varieties import K3, P2, P1xP1, abelian_surface, enriques, hirzebruch
from bridgeland_stability.dlp import delta as dlp_delta
from bridgeland_stability.dlp import moduli_nonempty as dlp_moduli_nonempty
from bridgeland_stability.exceptional import Bundle, enumerate_exceptional
from bridgeland_stability.chern import ChernChar
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.nonemptiness_rational import (
    validate_character,
    VerdictStatus,
    is_semiexceptional_p2,
)


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


# --------------------------------------------------------------------------
# E12-M1 (A4): the certified rank cutoff R_max = denominator(mu).  The old
# hard-coded R_max=60 missed the rank-89 DLP cusp at mu=34/89, spuriously
# lifting a class over its (too-low) truncated delta and returning a PROVEN
# non-empty verdict for a class that is EMPTY.  See docs/CORRECTIONS.md sec.8.
# --------------------------------------------------------------------------
def test_delta_H_certified_rank_cutoff_at_rank_89():
    # mu = 34/89 (denom 89 > old R_max=60); the default R_max now bumps to 89.
    xi = SurfaceBundle(89, (34,), F(0))
    assert delta_H(xi, P2) == F(3961, 7921)           # 1 - (1 - 1/89^2)/2 = 3961/7921
    v = moduli_nonempty(8010, (3060,), F(-3421), P2)
    assert v.discriminant == F(356489, 712890)         # 1/2(34/89)^2 + 3421/8010
    assert v.delta_H == F(3961, 7921)                  # = 356490/712890
    assert v.nonempty is False                         # 356489 < 356490 -> EMPTY
    assert v.exceptional is False                      # not (semi)exceptional either
    assert v.certificate.rigor == Rigor.PROVEN


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


# --------------------------------------------------------------------------
# E12-M2 (A1, A3): the P^2 semiexceptional disjunct + Thm 2.2 character validation,
# and the branch-derived VerdictStatus.  See docs/CORRECTIONS.md sec.8 (E12-M2).
# --------------------------------------------------------------------------
def test_p2_semiexceptional_direct_sum_nonempty():   # A1 (nonemptiness_rational side)
    v = moduli_nonempty(4, (2,), F(-1), P2)
    assert v.discriminant == F(3, 8)                 # 1/8 + 1/4
    assert v.delta_H == F(5, 8)
    assert v.discriminant < v.delta_H
    assert v.nonempty is True
    assert v.exceptional is False                    # rank 4 != denom 2: not a single exc bundle
    assert v.certificate.rigor == Rigor.PROVEN
    assert v.status is VerdictStatus.PROVEN_NONEMPTY
    assert "semiexceptional" in v.reason
    assert dlp_moduli_nonempty(Bundle(4, 2, F(-1)))["nonempty"] is True     # dlp twin agrees
    # O^{+2}: Delta=0 < delta(0)=1, semiexceptional via O=(1,0,0)
    w = moduli_nonempty(2, (0,), F(0), P2)
    assert w.discriminant == F(0) and w.delta_H == F(1)
    assert w.nonempty is True and w.exceptional is False
    assert w.status is VerdictStatus.PROVEN_NONEMPTY
    assert dlp_moduli_nonempty(Bundle(2, 0, F(0)))["nonempty"] is True


def test_p2_invalid_character_is_empty():            # A3
    v = moduli_nonempty(1, (0,), F(-3, 2), P2)       # chi = -1/2 not in Z
    assert v.nonempty is False
    assert v.status is VerdictStatus.PROVEN_EMPTY
    assert "invalid" in v.reason
    assert validate_character(1, (0,), F(-3, 2), P2) is False    # c2 = 3/2 not integral
    assert validate_character(2, (1,), F(0), P2) is False        # chi = 7/2 not integral
    assert validate_character(2, (1,), F(-5, 2), P2) is True     # chi = 1
    assert validate_character(4, (2,), F(-1), P2) is True        # chi = 6


def test_non_integral_rank_is_invalid():             # E13 re-audit R3
    """validate_character documents r in Z but only tested r < 1: Fraction(3,2)
    passed, and downstream theorem APIs int()-truncated it -- theorem-level answers
    for a DIFFERENT character than the caller supplied."""
    assert validate_character(F(3, 2), (0, 0), 0, P1xP1) is False
    assert validate_character(F(3, 2), (0,), F(0), P2) is False
    assert validate_character(F(4, 2), (0,), F(0), P2) is True    # 4/2 == 2 in Z: fine
    assert validate_character(0, (0,), F(0), P2) is False         # r < 1 still refused
    # end-to-end: the verdict is the invalid-character PROVEN_EMPTY, never a
    # silently truncated verdict for the r=1 character.
    v = moduli_nonempty(F(3, 2), (0, 0), 0, P1xP1)
    assert v.status is VerdictStatus.PROVEN_EMPTY and "invalid" in v.reason


def test_fe_invalid_character_is_empty_not_forged_nonempty():   # A3 off P^2 (E12-M2 IMPROVE)
    # The native F_e envelope path must reject a non-integral character BEFORE any
    # exceptional-bundle branch can forge a PROVEN non-empty verdict.  Each class below has a
    # non-integral c2 = 1/2 <c1,c1> - ch2, so no coherent sheaf has it and M_H(xi) is empty.
    # Before the fix these tripped _hirzebruch_verdict's "Delta > delta_H sharp" and
    # "exceptional bundle" branches -> a forged PROVEN_NONEMPTY (invariant-7's worst outcome).
    lat = P1xP1.ns_lattice                          # F_0 NS Gram [[0,1],[1,0]]
    forge = [
        (2, (0, 0), F(-7, 2)),   # <c1,c1>=0 -> c2 = 0 + 7/2 = 7/2   (was: Delta=7/4 > delta_H=1)
        (3, (1, 1), F(-9, 2)),   # <c1,c1>=2 -> c2 = 1 + 9/2 = 11/2  (was: exceptional branch)
        (2, (2, 2), F(-3, 2)),   # <c1,c1>=8 -> c2 = 4 + 3/2 = 11/2  (was: Delta=7/4 > delta_H=1)
    ]
    for r, c1, ch2 in forge:
        c2 = F(1, 2) * lat.self_pairing(tuple(F(x) for x in c1)) - ch2
        assert c2.denominator != 1                             # genuinely non-integral c2
        # rejected OFF P^2 -- the c2 clause is K_X-independent, no K_H repair needed:
        assert validate_character(r, c1, ch2, P1xP1) is False
        v = moduli_nonempty(r, c1, ch2, P1xP1)
        assert v.nonempty is False                             # the forged non-empty is gone
        assert v.status is VerdictStatus.PROVEN_EMPTY          # honest: emptiness is a theorem
        assert v.certificate.rigor == Rigor.PROVEN
        assert v.exceptional is False                          # no exceptional-branch forgery
        assert "invalid Chern character" in v.reason

    # Positive control: a VALID F_1 character (integral c2) is untouched -- the rank-2
    # exceptional bundle at the anticanonical polarization stays PROVEN non-empty.
    HK = hirzebruch_with_polarization(1, (3, 2))               # anticanonical ray on F_1
    xi = SurfaceBundle(2, (1, 1), F(-1, 2))                    # c2 = 1/2 + 1/2 = 1 (integral)
    assert validate_character(xi.r, xi.c1, xi.ch2, HK) is True
    vK = moduli_nonempty(xi.r, xi.c1, xi.ch2, HK)
    assert vK.nonempty is True and vK.exceptional is True
    assert vK.status is VerdictStatus.PROVEN_NONEMPTY
    assert vK.certificate.rigor == Rigor.PROVEN


def test_k3_with_hyperbolic_ns_never_mints_a_proof():   # E13 re-audit R1
    """A projective K3 with NS = U, H = (2,1), ch = (5,(-3,1),-3): the NS Gram is
    F_0-shaped, and Gram-only dispatch used to classify it as F_0 and return
    PROVEN_NONEMPTY via the exceptional-bundle branch.  The verdict is FALSE on a K3:
    the Mukai vector is v = (5,(-3,1),2), v^2 = <c1,c1> - 2 r s = -6 - 20 = -26 < -2,
    and c1.H = -1 is coprime to r = 5 (semistable => stable), so a stable sheaf would
    need v^2 >= -2 -- the moduli space is EMPTY.  Post-R1 the surface never enters the
    CH F_e theory: the verdict falls back to the honest HEURISTIC Bogomolov floor
    (UNKNOWN), and hirzebruch_index refuses the surface outright."""
    from bridgeland_stability.nslattice import NSLattice
    from bridgeland_stability.varieties import Surface
    from bridgeland_stability.dlp_hirzebruch import hirzebruch_index

    U = NSLattice(2, ((0, 1), (1, 0)))
    fake = Surface(name="K3 (NS = U)", d=4, K=(0, 0), chi_O=2, picard_rank=2,
                   kind="K3", H=(2, 1), ns_lattice=U)
    # Two-way: the Mukai arithmetic asserted above, recomputed exactly.
    c1 = (F(-3), F(1))
    assert U.self_pairing(c1) == F(-6)
    assert U.self_pairing(c1) - 2 * 5 * (5 + (-3)) == -26        # v^2, s = r + ch2
    assert U.pairing(c1, (F(2), F(1))) == F(-1)                  # c1.H, gcd(5,1)=1
    with pytest.raises(NotImplementedError):
        hirzebruch_index(fake)
    v = moduli_nonempty(5, (-3, 1), -3, fake)
    assert v.status is VerdictStatus.UNKNOWN                    # never PROVEN_NONEMPTY
    assert v.exceptional is False
    assert v.certificate.rigor is not Rigor.PROVEN


def test_verdict_status_is_branch_derived():
    assert moduli_nonempty(2, (1,), F(-5, 2), P2).status is VerdictStatus.PROVEN_NONEMPTY  # 11/8>=5/8
    assert moduli_nonempty(2, (1,), F(1, 2), P2).status is VerdictStatus.PROVEN_EMPTY       # -1/8<5/8, valid
    assert moduli_nonempty(1, (0,), F(-3, 2), P2).status is VerdictStatus.PROVEN_EMPTY      # invalid
    # off P^2: strict > gives PROVEN_NONEMPTY (target 1 == native sharp bound 1).
    above = moduli_nonempty(2, (0, 0), F(-4), P1xP1, delta_H_target=F(1), hn_source=HNMode.PAPER)
    assert above.discriminant == F(2) and above.status is VerdictStatus.PROVEN_NONEMPTY     # 2 > 1
    # E12-M4: target 2 != the package's native sharp bound 1 for this class -> a forged target,
    # refused, never a boundary verdict-flip that stays certified.  (The native Delta==delta_H
    # UNKNOWN band is pinned by test_certified_target_off_p2_band_is_unknown_not_proven_empty.)
    with pytest.raises(ValueError):
        moduli_nonempty(2, (0, 0), F(-4), P1xP1, delta_H_target=F(2), hn_source=HNMode.PAPER)


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

    # E12-M4 (A5): a target != the package's OWN certified sharp bound (native delta_H = 1
    # here) is a forge -> ValueError, never a verdict-flip that stays PROVEN.  Retires the
    # old "an absurd target overrides and flips the verdict" behaviour, which pinned the bug.
    with pytest.raises(ValueError):
        moduli_nonempty(2, (0, 0), F(-4), P1xP1, delta_H_target=F(5), hn_source=HNMode.PAPER)


def test_uncertified_target_rejected():
    with pytest.raises(ValueError):
        moduli_nonempty(2, (1, 0), 0, P1xP1,
                        delta_H_target=F(1, 64), hn_source=HNMode.HEURISTIC)
    with pytest.raises(ValueError):
        moduli_nonempty(2, (1, 0), 0, P1xP1, delta_H_target=F(1, 64))


# --------------------------------------------------------------------------
# E12-M4 / A5: class-bound sharp-bound evidence retires the forgeable
# (delta_H_target, hn_source) pair.  A target is honoured (-> PROVEN-eligible) only
# when it equals the package's OWN certified sharp bound for that exact class; a
# forged / mis-derived / mis-attributed target is a ValueError, never PROVEN.  See
# docs/CORRECTIONS.md "E12-M4".
# --------------------------------------------------------------------------
def test_sharp_bound_evidence_is_class_bound():
    # crit. 1: evidence built for one class cannot certify another.
    ev = SharpBoundEvidence(surface=P1xP1, r=2, c1=(0, 0), ch2=F(-4),
                            sharp_bound=F(1), sharp_bound_source=HNMode.PAPER,
                            hn_length_one_source=HNMode.PAPER, citation="test")
    assert ev.matches(2, (0, 0), F(-4), P1xP1) is True
    v = moduli_nonempty(2, (0, 0), F(-4), P1xP1, evidence=ev)     # right class, bound == native 1
    assert v.status is VerdictStatus.PROVEN_NONEMPTY and v.discriminant == F(2)
    assert v.delta_H == F(1) and v.mode is HNMode.PAPER
    with pytest.raises(ValueError):                               # wrong class (r matches, c1 differs)
        moduli_nonempty(2, (2, 2), F(-2), P1xP1, evidence=ev)


def test_sharp_bound_evidence_wrong_value_refused():
    # crit. 2 (the value gate): even for the RIGHT class, a sharp_bound != the package's own
    # certified sharp bound (native = 1) cannot mint a verdict -- it is a mis-derived claim.
    ev = SharpBoundEvidence(surface=P1xP1, r=2, c1=(0, 0), ch2=F(-4),
                            sharp_bound=F(5), sharp_bound_source=HNMode.PAPER,
                            hn_length_one_source=HNMode.PAPER, citation="wrong value")
    assert ev.matches(2, (0, 0), F(-4), P1xP1) is True           # class matches...
    with pytest.raises(ValueError):                              # ...but the value is forged
        moduli_nonempty(2, (0, 0), F(-4), P1xP1, evidence=ev)


def test_oracle_evidence_is_mint_guarded():
    # crit. 3: ORACLE-sourced evidence cannot be forged outside oracle/ (no capability token).
    with pytest.raises(TypeError):
        SharpBoundEvidence(surface=P1xP1, r=2, c1=(0, 0), ch2=F(-4), sharp_bound=F(1),
                           sharp_bound_source=HNMode.ORACLE, hn_length_one_source=HNMode.ORACLE)
    # a PAPER-sourced object is freely constructible (needed by the legacy path)
    ok = SharpBoundEvidence(surface=P1xP1, r=2, c1=(0, 0), ch2=F(-4), sharp_bound=F(1),
                            sharp_bound_source=HNMode.PAPER, hn_length_one_source=HNMode.PAPER)
    assert ok.hn_length_one_source is HNMode.PAPER


def test_raw_oracle_target_is_refused():
    # crit. 3: even a CORRECT raw ORACLE target is refused -- ORACLE must be a capability object
    # minted after a real construction, never a raw (delta_H_target, hn_source=ORACLE) pair.
    with pytest.raises(ValueError):
        moduli_nonempty(2, (0, 0), F(-4), P1xP1, delta_H_target=F(1), hn_source=HNMode.ORACLE)


def test_p2_forged_target_refused_both_paths():
    # roadmap crit. 5: the P^2 PAPER and P^2 ORACLE forges.  native delta(0|P2)=1 != 0, so a
    # target of 0 is a forge on BOTH source paths -- the same "fixed one path, missed its twin"
    # shape as A4b.  Without the fix, (3,(0),-2, target=0, ORACLE) minted PROVEN nonempty=True
    # for a class that is natively PROVEN empty.
    for src in (HNMode.PAPER, HNMode.ORACLE):
        with pytest.raises(ValueError):
            moduli_nonempty(3, (0,), F(-2), P2, delta_H_target=F(0), hn_source=src)
    # native truth is unchanged: PROVEN empty (Delta = 2/3 < delta(0) = 1, rank 3 not Markov).
    assert moduli_nonempty(3, (0,), F(-2), P2).status is VerdictStatus.PROVEN_EMPTY


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


def test_certified_target_off_p2_keeps_exceptional_disjunct():
    # E12-M2 regression (defect B: the certificate forger).  The OFF-P^2 analogue of
    # test_paper_exceptional_coexists_with_target: a GENUINE mu_H-stable exceptional bundle
    # fed its OWN correct sharp delta_H via the certified PAPER-target path must stay
    # non-empty (its moduli space is a single reduced point).  Pre-fix the certified-target
    # branch dropped the exceptional/semiexceptional disjunct off P^2, so it returned
    # nonempty=False with Rigor.PROVEN -> a FALSE PROVEN_EMPTY contradicting the SAME
    # function's native verdict for the identical class.
    #
    # xi = (3, (1,1), -1) on F_0 = P^1 x P^1.  nu = (1/3, 1/3), so the full-NS
    # Delta = 1/2<nu,nu> - ch2/r = 1/2*(2/9) + 1/3 = 4/9 == Delta_V(3) = (1 - 1/9)/2: a
    # rank-3 mu_H-stable exceptional bundle.  Its sharp anticanonical envelope is
    # delta_H = 5/9 (exact, sharp), and 4/9 < 5/9 -- it sits STRICTLY BELOW the envelope.
    r, c1, ch2 = 3, (1, 1), F(-1)
    xi = SurfaceBundle(r, c1, ch2)
    assert validate_character(r, c1, ch2, P1xP1) is True
    assert discriminant(xi, P1xP1) == F(4, 9) == exceptional_discriminant(3)   # Delta == Delta_V
    assert is_stable_exceptional(r, c1, P1xP1) is True                         # genuine exc bundle
    sharp = delta_H(xi, P1xP1)
    assert sharp == F(5, 9)                                                    # sharp envelope
    assert discriminant(xi, P1xP1) < sharp                                    # strictly below

    native = moduli_nonempty(r, c1, ch2, P1xP1)
    assert native.status is VerdictStatus.PROVEN_NONEMPTY                      # native: correct

    # The certified-target path, fed the class's OWN sharp delta_H, must AGREE with native.
    tgt = moduli_nonempty(r, c1, ch2, P1xP1, delta_H_target=sharp, hn_source=HNMode.PAPER)
    assert tgt.nonempty is True                                               # NOT a false empty
    assert tgt.exceptional is True                                            # via the disjunct
    assert tgt.certificate.rigor == Rigor.PROVEN
    assert tgt.status is VerdictStatus.PROVEN_NONEMPTY                        # == native
    assert tgt.mode is HNMode.PAPER

    # A6 guard: the same (r, c1) with a WRONG ch2 (Delta != Delta_V) is NOT an exceptional
    # bundle, so the disjunct must NOT rescue it -- the target verdict stays PROVEN_EMPTY.
    impostor = moduli_nonempty(3, (1, 1), F(0), P1xP1,
                               delta_H_target=sharp, hn_source=HNMode.PAPER)
    assert discriminant(SurfaceBundle(3, (1, 1), F(0)), P1xP1) != exceptional_discriminant(3)
    assert impostor.nonempty is False and impostor.exceptional is False
    assert impostor.status is VerdictStatus.PROVEN_EMPTY


def test_A6_native_ch2_guard_is_proven_empty():
    # E12-M3.  (3,(1,1),0) on F_0 = P^1xP^1 (Gram [[0,1],[1,0]]): nu=(1/3,1/3).
    # <c1,c1> = 2*1*1 = 2;  Delta = 1/2*(2/9) - 0/3 = 1/9.
    # exceptional_ch2(3,(1,1)) = <c1,c1>/(2r) - r*Delta_V = 2/6 - 3*(4/9) = 1/3 - 4/3 = -1 != 0,
    #   so (r,c1) carries a rank-3 mu_H-stable exceptional bundle but THIS ch2 (=0) is not it.
    # emptiness_bound = 5/9 (theorem-branch max) > Delta = 1/9  ->  PROVEN_EMPTY.
    # Pre-fix _hirzebruch_verdict stamped exceptional=True off (r,c1) alone (the A6 bug);
    # the ch2 guard now requires ch2 == exceptional_ch2, so the class falls through to the
    # certified emptiness branch.  This pins the NATIVE path (no delta_H_target); the
    # certified-target twin is pinned by test_certified_target_off_p2_keeps_exceptional_disjunct.
    r, c1, ch2 = 3, (1, 1), F(0)
    xi = SurfaceBundle(r, c1, ch2)
    assert discriminant(xi, P1xP1) == F(1, 9)
    assert exceptional_discriminant(3) == F(4, 9)
    assert exceptional_ch2(3, c1, P1xP1) == F(-1)          # ch2 guard: -1 != 0
    assert is_stable_exceptional(r, c1, P1xP1) is True     # the (r,c1)-only test still fires
    assert is_semiexceptional(xi, P1xP1) is False          # Delta != Delta_V
    assert emptiness_bound(xi, P1xP1) == F(5, 9)
    assert discriminant(xi, P1xP1) < emptiness_bound(xi, P1xP1)
    v = moduli_nonempty(r, c1, ch2, P1xP1)
    assert v.nonempty is False
    assert v.exceptional is False
    assert v.status is VerdictStatus.PROVEN_EMPTY
    assert v.certificate.rigor is Rigor.PROVEN             # theorem-backed emptiness, not UNKNOWN


def test_hirzebruch_branch_disjointness_no_double_fire():
    # E12-M3.  The native _hirzebruch_verdict has two PROVEN-empty branch predicates
    # (Delta<0 Bogomolov; Delta<emptiness_bound) and two PROVEN-nonempty ones
    # ((semi)exceptional; env.certified_sharp and Delta>delta_H).  A6 was exactly a class
    # -- (3,(1,1),0) on F_0 -- that satisfied BOTH a PROVEN-empty predicate (Delta=1/9<5/9)
    # and the (buggy) PROVEN-nonempty exceptional predicate; only the source ORDER of the
    # branches decided the verdict.  With the ch2 guard the two families are provably disjoint,
    # which this crisp invariant pins: no character is reported PROVEN_NONEMPTY while sitting
    # strictly below its own certified emptiness_bound (nor below the Bogomolov floor 0).
    #
    # Runtime: moduli_nonempty off P^2 enumerates the DLP envelope (~65 ms/call), so the
    # ROADMAP index range (r<=8, |c1|<=6, a ch2 box) is *hours* of wall time.  This in-suite
    # box is bounded for runtime yet still (a) contains the A6 witness (3,(1,1),0) on F_0 and
    # (b) fires BOTH proven branches many times (asserted non-vacuous below).  A wider box
    # (r<=6, |c1|<=4, c2 in [-4,4], ~8.7k characters) was verified OFFLINE with zero
    # double-fires; see docs/CORRECTIONS.md section 8, subsection "E12-M3".
    F0 = hirzebruch_with_polarization(0, (1, 1))     # -K, del Pezzo, e=0
    F1 = hirzebruch_with_polarization(1, (3, 2))     # -K = 3f+2s, del Pezzo, e=1
    checked = n_nonempty = n_empty = 0
    saw_a6 = False
    for surf in (F0, F1):
        for r in range(1, 7):
            for a in range(-1, 2):
                for b in range(-1, 2):
                    c1 = (a, b)
                    c1c1 = surf.lattice.self_pairing((F(a), F(b)))
                    for c2 in range(0, 3):
                        ch2 = F(c1c1, 2) - c2               # integral-c2 => valid character
                        xi = SurfaceBundle(r, c1, ch2)
                        v = moduli_nonempty(r, c1, ch2, surf)
                        checked += 1
                        if surf is F0 and r == 3 and c1 == (1, 1) and ch2 == 0:
                            saw_a6 = True                  # the A6 witness is in-box
                            assert v.nonempty is False
                            assert v.status is VerdictStatus.PROVEN_EMPTY
                        if v.status is VerdictStatus.PROVEN_NONEMPTY:
                            n_nonempty += 1
                            # a PROVEN-nonempty class must NOT also satisfy a PROVEN-empty predicate
                            disc = discriminant(xi, surf)
                            eb = emptiness_bound(xi, surf)
                            assert disc >= 0, (surf.name, r, c1, ch2, "Bogomolov")
                            assert disc >= eb, (surf.name, r, c1, ch2, disc, eb)
                        elif v.status is VerdictStatus.PROVEN_EMPTY:
                            n_empty += 1
    assert checked > 0
    assert saw_a6                                          # A6 cell (3,(1,1),0)/F_0 was swept
    assert n_nonempty > 0 and n_empty > 0                  # non-vacuous: both families fire


def test_certified_target_off_p2_band_is_unknown_not_proven_empty():
    # E12-M2 regression (IMPROVE round 3, the certificate forger).  A NON-(semi)exceptional class
    # whose full-NS Delta lands in the gap emptiness_bound <= Delta < delta_H, fed its OWN correct
    # sharp delta_H through the certified PAPER-target path, must return UNKNOWN -- the SAME value
    # the function's native envelope verdict returns for the identical class.  Emptiness is NOT a
    # theorem in that band: emptiness_bound is strictly weaker than the envelope (it drops the
    # non-theorem (nu - nu_V).H = 0 branch), so a PROVEN_EMPTY there is a non-theorem-backed false
    # verdict (invariant-7's worst outcome).  Pre-fix the certified-target tail read emptiness off
    # the flat Delta < delta_H and forged PROVEN_EMPTY across the whole gap.
    #
    # xi = (3, (1,2), -1) on F_0 = P^1 x P^1.  nu = (1/3, 2/3); full-NS Delta = 5/9,
    # emptiness_bound = 2/9, sharp delta_H = 8/9, so 2/9 <= 5/9 < 8/9 -- strictly in the gap.
    r, c1, ch2 = 3, (1, 2), F(-1)
    xi = SurfaceBundle(r, c1, ch2)
    assert validate_character(r, c1, ch2, P1xP1) is True
    assert is_stable_exceptional(r, c1, P1xP1) is False        # not a single exceptional bundle
    assert is_semiexceptional(xi, P1xP1) is False              # nor a sum of copies of one
    eb, disc, sharp = emptiness_bound(xi, P1xP1), discriminant(xi, P1xP1), delta_H(xi, P1xP1)
    assert (eb, disc, sharp) == (F(2, 9), F(5, 9), F(8, 9))
    assert eb <= disc < sharp                                  # strictly inside the gap

    native = moduli_nonempty(r, c1, ch2, P1xP1)                # native envelope verdict
    assert native.status is VerdictStatus.UNKNOWN             # honest: existence undecided here
    assert native.certificate.rigor is not Rigor.PROVEN

    # The certified-target path, fed the class's OWN sharp delta_H, must AGREE with native --
    # not forge a PROVEN_EMPTY off the flat Delta < delta_H comparison.
    tgt = moduli_nonempty(r, c1, ch2, P1xP1, delta_H_target=sharp, hn_source=HNMode.PAPER)
    assert tgt.delta_H == sharp == native.delta_H             # identical bound; only the path differs
    assert tgt.status is VerdictStatus.UNKNOWN               # FIXED: was a forged PROVEN_EMPTY
    assert tgt.certificate.rigor is not Rigor.PROVEN
    assert tgt.exceptional is False

    # The band downgrade is NOT a blanket "always UNKNOWN": the two theorem-backed boundaries of
    # the certified-target path are untouched.  Below emptiness_bound the verdict stays PROVEN
    # empty (a mu_H-stable exceptional bundle forces chi <= 0), and strictly above delta_H it
    # stays PROVEN nonempty (Thm "deltaSurface" (1), which needs a STRICT >).
    below_xi = SurfaceBundle(3, (1, 1), F(0))                  # nu = (1/3,1/3): Delta=1/9 < eb=5/9
    below = moduli_nonempty(3, (1, 1), F(0), P1xP1, delta_H_target=F(5, 9), hn_source=HNMode.PAPER)
    assert (discriminant(below_xi, P1xP1), emptiness_bound(below_xi, P1xP1)) == (F(1, 9), F(5, 9))
    assert below.status is VerdictStatus.PROVEN_EMPTY and below.exceptional is False
    above_xi = SurfaceBundle(3, (1, 2), F(-3))                 # nu = (1/3,2/3): Delta=11/9 > dH=8/9
    above = moduli_nonempty(3, (1, 2), F(-3), P1xP1, delta_H_target=F(8, 9), hn_source=HNMode.PAPER)
    assert discriminant(above_xi, P1xP1) == F(11, 9)
    assert above.status is VerdictStatus.PROVEN_NONEMPTY


def test_target_less_certified_source_off_p2_is_not_forged_proven():
    # E12-M2 regression (IMPROVE round 4, the certificate forger, target-LESS variant).  Every
    # earlier Defect-B round concerned the delta_H_target-SUPPLIED path; this one is the
    # documented E11-M5 hook used WITHOUT a target.  Off P^2 with no delta_H_target the package's
    # only bound is its OWN native envelope.  A bare certified hn_source (ORACLE/PAPER/DLP)
    # certifies the HN-length-one hypothesis, NOT a sharp delta_H -- so it must not raise the
    # verdict to PROVEN when that envelope is a certified LOWER bound (env.sharp False) or absent
    # (K3 / abelian / nef-big F_n, delta_H = the Bogomolov floor 0).  Pre-fix, the tail stamped
    # _MODE_CERT[hn_source] = PROVEN over env.value, forging a PROVEN_NONEMPTY for every class in
    # the gap [env.value, sharp delta_H): a class that is EMPTY reported as PROVEN non-empty
    # (invariant-7's worst outcome).  The verdict must instead EQUAL the native (no-source) one.
    certified = (HNMode.ORACLE, HNMode.PAPER, HNMode.DLP)

    # (1) Lower-bound envelope: F_0 with an ample, NON-anticanonical H = 2f + s (the finding's
    #     repro).  env.sharp is False, so DLP_H(nu) is only a certified lower bound (= the 1/2
    #     floor); disc = 2 >= 1/2 does NOT imply disc >= the (unknown, larger) sharp delta_H.
    S = hirzebruch_with_polarization(0, (2, 1))
    r, c1, ch2 = 2, (-3, -2), F(-1)
    env = dlp_envelope(total_slope(SurfaceBundle(r, c1, ch2)), S, 20)
    assert env.sharp is False and env.certified_sharp is False and env.value == F(1, 2)
    native = moduli_nonempty(r, c1, ch2, S)
    assert native.status is VerdictStatus.UNKNOWN and native.certificate.rigor is Rigor.HEURISTIC
    for src in certified:
        v = moduli_nonempty(r, c1, ch2, S, hn_source=src)      # certified source, NO target
        assert v.status is VerdictStatus.UNKNOWN               # FIXED: was a forged PROVEN_NONEMPTY
        assert v.certificate.rigor is not Rigor.PROVEN
        # a bare label adds no sharp bound: the verdict is IDENTICAL to the native one
        assert (v.status, v.certificate.rigor, v.nonempty, v.delta_H) == (
            native.status, native.certificate.rigor, native.nonempty, native.delta_H)
        assert v.delta_H == env.value                          # compared against the LOWER bound

    # (2) No envelope at all (env is None): K3, abelian, and a nef-and-big (non-ample) F_2 factory
    #     polarization.  delta_H is the Bogomolov floor 0 -- never sharp -- so a bare certified
    #     source cannot certify existence against it either.
    for surface, cc in ((K3(2), (1,)), (abelian_surface(2), (1,)), (hirzebruch(2), (1, 0))):
        base = moduli_nonempty(2, cc, F(-1), surface)
        assert base.status is VerdictStatus.UNKNOWN and base.delta_H == 0
        for src in certified:
            v = moduli_nonempty(2, cc, F(-1), surface, hn_source=src)
            assert v.status is VerdictStatus.UNKNOWN and v.certificate.rigor is not Rigor.PROVEN
            assert (v.status, v.nonempty, v.delta_H) == (base.status, base.nonempty, base.delta_H)

    # (3) The fix is NOT over-broad: where the native envelope IS certified sharp (F_0 = P^1 x P^1,
    #     the anticanonical del Pezzo) a class strictly ABOVE its sharp delta_H still reaches
    #     PROVEN_NONEMPTY (Thm "deltaSurface" (1)) -- with or without a bare certified source,
    #     which is honoured only because the envelope itself is sharp, per-branch, never by label.
    rs, cs, es = 2, (0, 0), F(-4)                              # disc = 2 > sharp delta_H = 1
    envK = dlp_envelope(total_slope(SurfaceBundle(rs, cs, es)), P1xP1, 20)
    assert envK.certified_sharp is True and envK.value == F(1)
    for kw in ({}, {"hn_source": HNMode.ORACLE}, {"hn_source": HNMode.PAPER}):
        v = moduli_nonempty(rs, cs, es, P1xP1, **kw)
        assert v.status is VerdictStatus.PROVEN_NONEMPTY and v.discriminant == F(2)


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
