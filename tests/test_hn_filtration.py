"""E13-M3a / G18: the HN-length-one existence criterion + Thm 1.13 structure.

Primary source: arXiv:1907.06739 (Coskun-Huizenga, "Existence of semistable sheaves on
Hirzebruch surfaces") Sec. 1.6 (the HN-length-one theorem), Sec. 5 (the generic HN
filtration; Thm 1.6), Sec. 7 (Thm 1.13 = Cor 7.7), Example 1.9, Example 1.14.  See
``docs/CORRECTIONS.md`` Sec. 11 for the two-way evidence.

The load-bearing theorem (Sec. 1, verbatim): "there exists an H_m-semistable sheaf with
Chern character v if and only if the generic H_m-Harder-Narasimhan filtration has
length 1."

Every asserted value is anchored two ways: (1) the exact ``Fraction`` discriminant
``Delta = 1/2 <nu,nu> - ch2/r`` (the full-NS CH discriminant, CLAUDE.md invariant 2 --
NEVER ``discriminant_H``), hand-recomputed in ``docs/CORRECTIONS.md`` Sec. 11 and pinned
below, and (2) the tri-state verdict, which must agree bit-for-bit with the shipped
:func:`bridgeland_stability.nonemptiness_rational.moduli_nonempty` status (the
no-fabrication delegation guarantee) and, on the anticanonical del Pezzo ray, with the
certified-sharp :func:`bridgeland_stability.dlp_hirzebruch.dlp_envelope`.
"""

import os
import sys
from collections import namedtuple
from fractions import Fraction as F

import pytest

from bridgeland_stability import (
    HNRegion,
    HNVerdict,
    hn_verdict,
    hn_region,
    semistable_exists,
    generic_hn_length,
    THM_1_13_MIN_DELTA,
    THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS,
)
from bridgeland_stability.dlp_hirzebruch import (
    total_slope,
    discriminant,
    dlp_envelope,
    exceptional_discriminant,
    is_del_pezzo_anticanonical,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle
from bridgeland_stability.nonemptiness_rational import (
    moduli_nonempty,
    VerdictStatus,
    hirzebruch_with_polarization,
)
from bridgeland_stability.varieties import P2, P1xP1, K3, hirzebruch

# --- oracle import (tests/oracle/; put tests/ on sys.path) ------------------- #
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _TESTS_DIR not in sys.path:
    sys.path.insert(0, _TESTS_DIR)

from oracle.dlp_reference import reference_semistable_exists
from oracle.corpus import CORPUS


# --------------------------------------------------------------------------- #
# Surfaces.  F_0 = P1xP1 with H=(1,1) and F_1 with H=(3,2) are both ample       #
# ANTICANONICAL del Pezzo, where dlp_envelope is certified sharp.               #
# --------------------------------------------------------------------------- #
F0 = P1xP1
F1 = hirzebruch_with_polarization(1, (3, 2))
F2 = hirzebruch_with_polarization(2, (5, 2))          # e=2: out of M3a scope

_STATUS_TO_EXISTS = {
    VerdictStatus.PROVEN_NONEMPTY: True,
    VerdictStatus.PROVEN_EMPTY: False,
    VerdictStatus.UNKNOWN: None,
}


def se(r, c1, ch2, S):
    """Mirror a user: derive ``(nu, Delta)`` from the character, then ask the criterion."""
    xi = SurfaceBundle(r, c1, ch2)
    return semistable_exists(r, total_slope(xi), discriminant(xi, S), S)


def region(r, c1, ch2, S):
    xi = SurfaceBundle(r, c1, ch2)
    return hn_region(r, total_slope(xi), discriminant(xi, S), S)


def verdict(r, c1, ch2, S):
    xi = SurfaceBundle(r, c1, ch2)
    return hn_verdict(r, total_slope(xi), discriminant(xi, S), S)


# --------------------------------------------------------------------------- #
# The pinned anchor table (all values probe-confirmed against shipped code and  #
# hand-recomputed in docs/CORRECTIONS.md Sec. 11).  Each row carries the exact  #
# full-NS Delta, the tri-state existence verdict, and the HN region.            #
# --------------------------------------------------------------------------- #
Anchor = namedtuple("Anchor", "label r c1 ch2 surface Delta exists region")

ANCHORS = (
    # ---- P^2 (d=1): region S / exceptional disjunct / certified empty. --------
    Anchor("P2_S_ideal_sheaf", 1, (0,), F(-5), P2, F(5), True, HNRegion.S),
    Anchor("P2_S_above_curve", 2, (1,), F(-5, 2), P2, F(11, 8), True, HNRegion.S),
    Anchor("P2_exc_T(-1)", 2, (1,), F(-1, 2), P2, F(3, 8), True, HNRegion.S),
    Anchor("P2_exc_rank5", 5, (2,), F(-2), P2, F(12, 25), True, HNRegion.S),
    Anchor("P2_empty_rank3", 3, (0,), F(-2), P2, F(2, 3), False, HNRegion.EMPTY),
    # ---- F_0 = P1xP1, H=(1,1), d=2, Gram [[0,1],[1,0]]. -----------------------
    Anchor("F0_S_slope00", 2, (0, 0), F(-4), F0, F(2), True, HNRegion.S),
    Anchor("F0_S_diag", 2, (2, 2), F(-2), F0, F(2), True, HNRegion.S),
    Anchor("F0_exc_rank3", 3, (1, 1), F(-1), F0, F(4, 9), True, HNRegion.S),
    Anchor("F0_exc_rank5", 5, (1, 2), F(-2), F0, F(12, 25), True, HNRegion.S),
    Anchor("F0_band_flagship", 2, (1, 1), F(0), F0, F(1, 4), True, HNRegion.S),   # decided by M3b
    Anchor("F0_empty_band", 2, (1, 1), F(1, 2), F0, F(0), False, HNRegion.EMPTY),
    Anchor("F0_empty_below_eb", 3, (0, 0), F(-1), F0, F(1, 3), False, HNRegion.EMPTY),
    Anchor("F0_bogomolov", 2, (0, 0), F(1, 2), F0, F(-1, 4), False, HNRegion.EMPTY),
    # ---- F_1, H=(3,2), d=8, Gram [[0,1],[1,-1]]. ------------------------------
    Anchor("F1_S", 2, (0, 0), F(-4), F1, F(2), True, HNRegion.S),
    Anchor("F1_bogomolov", 2, (0, 0), F(1), F1, F(-1, 2), False, HNRegion.EMPTY),
    Anchor("F1_boundary", 2, (0, 0), F(-2), F1, F(1), True, HNRegion.S),          # decided by M3b
    Anchor("F1_empty", 2, (1, 1), F(0), F1, F(1, 8), False, HNRegion.EMPTY),
)

_P2_ANCHORS = tuple(a for a in ANCHORS if a.surface is P2)
_F0_ANCHORS = tuple(a for a in ANCHORS if a.surface is F0)
_F1_ANCHORS = tuple(a for a in ANCHORS if a.surface is F1)


# --------------------------------------------------------------------------- #
# 0. Two-way exact-Fraction anchor: the recomputed full-NS Delta.              #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("a", ANCHORS, ids=[a.label for a in ANCHORS])
def test_anchor_discriminant_is_exact(a):
    """Delta(xi) matches the hand-recomputed full-NS value (docs/CORRECTIONS.md Sec. 11)."""
    xi = SurfaceBundle(a.r, a.c1, a.ch2)
    assert discriminant(xi, a.surface) == a.Delta


# --------------------------------------------------------------------------- #
# 1-3. Per-anchor tri-state existence + region (P^2, F_0, F_1).                #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("a", ANCHORS, ids=[a.label for a in ANCHORS])
def test_anchor_semistable_exists(a):
    """semistable_exists reproduces the pinned tri-state (True / False / None)."""
    assert se(a.r, a.c1, a.ch2, a.surface) is a.exists


@pytest.mark.parametrize("a", ANCHORS, ids=[a.label for a in ANCHORS])
def test_anchor_region(a):
    """The HN region label (Example 1.14: S / K / EMPTY) matches the pinned value."""
    assert region(a.r, a.c1, a.ch2, a.surface) is a.region


def test_p2_is_a_total_verdict_over_a_grid():
    """P^2: the Drezet-Le Potier closed form is sharp everywhere, so the verdict is
    total -- every character is decided True/False, NEVER None (no K region on P^2)."""
    for r in range(1, 7):
        for c1 in range(-8, 9):
            for c2 in range(0, 7):
                ch2 = F(c1 * c1, 2) - c2          # integral-c2 character
                v = se(r, (c1,), ch2, P2)
                assert v is not None, (r, c1, ch2)
                assert isinstance(v, bool)


def test_flagship_decides_to_S_with_computed_length_one():
    """The flagship: (2,(1,1),0) on F_0 IS O(1,0) (+) O(0,1) -- polystable (the summands
    share a reduced Hilbert polynomial), so a semistable sheaf EXISTS and the Sec. 1.6
    criterion gives generic HN length one.  The E13 re-audit R2 class: M3a honestly
    reported UNKNOWN here; E13-M3b now COMPUTES the generic HN filtration and decides
    True with length 1 -- matching the polystable truth (the paper exhibits exactly
    this bundle as -K-semistable with Delta = 1/4 < DLP_{-K}, arXiv:1907.06739 Sec. 7,
    the example after cor-delPezzoKss)."""
    # Both rank-1 summands O(1,0) and O(0,1) exist (region S).
    assert se(1, (1, 0), F(0), F0) is True
    assert se(1, (0, 1), F(0), F0) is True
    assert region(1, (1, 0), F(0), F0) is HNRegion.S
    # The direct sum ch = (2,(1,1),0) is DECIDED: exists, computed length 1.
    assert se(2, (1, 1), F(0), F0) is True
    assert region(2, (1, 1), F(0), F0) is HNRegion.S
    v = verdict(2, (1, 1), F(0), F0)
    assert v.generic_hn_length == 1 and v.factors == ((2, (1, 1), F(0)),)
    # The band structure that made M3a defer: emptiness_bound <= Delta <= DLP_{-K}.
    xi = SurfaceBundle(2, (1, 1), F(0))
    assert discriminant(xi, F0) == F(1, 4)
    assert dlp_envelope(total_slope(xi), F0).value == F(3, 4)


def test_F1_boundary_is_decided_by_the_computed_filtration():
    """F_1 (2,(0,0),-2): Delta=1 lands exactly ON the sharp envelope DLP_{-K}(0,0)=1 and
    on emptiness_bound=1, where CH Thm 'deltaSurface' (1) is silent (strict inequality).
    E13-M3b decides the boundary by the Sec. 5 algorithm: no maximal destabilizing
    character passes the cor-algorithm iff (hand-checks in tests/test_generic_hn.py),
    so semistable sheaves EXIST at the boundary -- closing this instance of the E11-M6
    open question O2."""
    xi = SurfaceBundle(2, (0, 0), F(-2))
    assert discriminant(xi, F1) == F(1)
    assert dlp_envelope(total_slope(xi), F1).value == F(1)
    assert se(2, (0, 0), F(-2), F1) is True
    assert region(2, (0, 0), F(-2), F1) is HNRegion.S


# --------------------------------------------------------------------------- #
# 4. Anticanonical-ray regression: bit-for-bit tie to dlp_envelope.            #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("a", _F0_ANCHORS + _F1_ANCHORS,
                         ids=[a.label for a in _F0_ANCHORS + _F1_ANCHORS])
def test_sharp_bound_ties_to_certified_envelope(a):
    """On the anticanonical del Pezzo ray the verdict's sharp_bound IS the certified-sharp
    dlp_envelope value -- the regression tying M3a to the shipped sharp theory."""
    assert is_del_pezzo_anticanonical(a.surface)      # F_0/F_1 anchors are all on the -K ray
    xi = SurfaceBundle(a.r, a.c1, a.ch2)
    env = dlp_envelope(total_slope(xi), a.surface)
    v = verdict(a.r, a.c1, a.ch2, a.surface)
    assert v.sharp_bound == env.value
    # And the exists mapping never contradicts the moduli_nonempty status: a PROVEN
    # status binds; an UNKNOWN status is decided by the M3b filtration (bool).
    mv = moduli_nonempty(a.r, a.c1, a.ch2, a.surface)
    expected = _STATUS_TO_EXISTS[mv.status]
    if expected is None:
        assert isinstance(v.exists, bool)
    else:
        assert v.exists is expected


# --------------------------------------------------------------------------- #
# 5. Delegation invariance: semistable_exists == status-map(moduli_nonempty).  #
#    The no-fabrication guarantee -- M3a re-derives NO verdict logic.           #
# --------------------------------------------------------------------------- #
def _delegation_grid(surface, ranks, coords, ch2s):
    dim = 1 if surface.is_p2 else 2
    if dim == 1:
        c1_iter = [(c,) for c in coords]
    else:
        c1_iter = [(a, b) for a in coords for b in coords]
    for r in ranks:
        for c1 in c1_iter:
            for ch2 in ch2s:
                yield r, c1, ch2


def test_delegation_invariance_P2():
    for r, c1, ch2 in _delegation_grid(
        P2, range(1, 6), range(-5, 6), [F(x, 2) for x in range(-10, 5)]
    ):
        v = moduli_nonempty(r, c1, ch2, P2)
        assert se(r, c1, ch2, P2) is _STATUS_TO_EXISTS[v.status], (r, c1, ch2)


def test_delegation_invariance_F0():
    # F_e moduli_nonempty enumerates the envelope (expensive) and this grid hits it
    # TWICE per point (directly + via se), so it is kept small; it still spans both c1
    # parities and the S/K/EMPTY bands.  Since E13-M3b the delegation contract is:
    # a PROVEN moduli status BINDS the verdict; an UNKNOWN status is DECIDED by the
    # computed generic HN filtration (a bool, never None) -- and never contradicted.
    for r, c1, ch2 in _delegation_grid(
        F0, range(1, 3), range(-1, 2), [F(x, 2) for x in range(-2, 3)]
    ):
        v = moduli_nonempty(r, c1, ch2, F0)
        got = se(r, c1, ch2, F0)
        expected = _STATUS_TO_EXISTS[v.status]
        if expected is None:
            assert isinstance(got, bool), (r, c1, ch2)     # M3b decided the band
        else:
            assert got is expected, (r, c1, ch2)


def test_delegation_invariance_F1():
    for r, c1, ch2 in _delegation_grid(
        F1, range(1, 3), range(-1, 2), [F(x, 2) for x in range(-2, 3)]
    ):
        v = moduli_nonempty(r, c1, ch2, F1)
        got = se(r, c1, ch2, F1)
        expected = _STATUS_TO_EXISTS[v.status]
        if expected is None:
            assert isinstance(got, bool), (r, c1, ch2)     # M3b decided the band
        else:
            assert got is expected, (r, c1, ch2)


# --------------------------------------------------------------------------- #
# 6. generic_hn_length: 1 iff semistable_exists is True, else None.            #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("a", ANCHORS, ids=[a.label for a in ANCHORS])
def test_generic_hn_length_matches_existence(a):
    xi = SurfaceBundle(a.r, a.c1, a.ch2)
    nu, D = total_slope(xi), discriminant(xi, a.surface)
    length = generic_hn_length(a.r, nu, D, a.surface)
    if a.exists is True:
        assert length == 1
    else:
        assert length is None


def test_generic_hn_length_over_grid():
    """length == 1 exactly where sheaves exist; M3b-decided empties carry their
    exact computed length; envelope-decided empties carry None."""
    # P^2 (fast, closed form) keeps a wide grid; the F_e sweeps (expensive envelope
    # enumeration) are kept tight but still span both parities and the S/K/EMPTY bands.
    for surface, coords, r_hi in ((P2, range(-4, 5), 4), (F0, range(-1, 2), 3), (F1, range(-1, 2), 3)):
        dim = 1 if surface.is_p2 else 2
        c1s = [(c,) for c in coords] if dim == 1 else [(a, b) for a in coords for b in coords]
        for r in range(1, r_hi):
            for c1 in c1s:
                for ch2 in [F(x, 2) for x in range(-6, 3)] if surface.is_p2 else [F(x, 2) for x in range(-4, 2)]:
                    xi = SurfaceBundle(r, c1, ch2)
                    nu, D = total_slope(xi), discriminant(xi, surface)
                    length = generic_hn_length(r, nu, D, surface)
                    exists = semistable_exists(r, nu, D, surface)
                    assert (length == 1) is (exists is True), (surface.name, r, c1, ch2)
                    # exact computed lengths (2..4) appear where M3b decided emptiness;
                    # None only on the envelope-decided empty paths.
                    assert length in (1, 2, 3, 4, None)


# --------------------------------------------------------------------------- #
# 7. Thm 1.13 = Cor 7.7 / Example 1.14 structural pins.                        #
# --------------------------------------------------------------------------- #
def test_thm_1_13_min_delta_is_three_eighths_two_ways():
    """The Thm 1.13 threshold Delta >= 3/8, pinned two ways (CLAUDE.md invariant 3):
    3/8 == exceptional_discriminant(2) == (1 - 1/2^2)/2."""
    assert THM_1_13_MIN_DELTA == F(3, 8)
    assert THM_1_13_MIN_DELTA == exceptional_discriminant(2)
    assert exceptional_discriminant(2) == F(1, 2) - F(1, 2 * 2 * 2)   # (1 - 1/4)/2


def test_thm_1_13_at_most_one_non_semiexceptional_factor():
    """Thm 1.13 shape: at most one HN factor of the general prioritary sheaf is not a
    semiexceptional bundle -- so region K has exactly one Kronecker-module factor."""
    assert THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS == 1


def test_verdict_region_labels_on_P2_and_F0():
    """The package's verdict-region labels where the DLP curve is known: S (PROVEN
    length 1) / UNCLASSIFIED (honest UNKNOWN band) / EMPTY (certified).  These label
    the VERDICT, not the sheaf's generic-HN shape (E13 re-audit R2): Example 1.14's
    S/K/R/empty shapes occur on P^2 too, but the P^2 existence boolean is total, so
    UNCLASSIFIED never fires there."""
    assert region(1, (0,), F(-5), P2) is HNRegion.S
    assert region(3, (0,), F(-2), P2) is HNRegion.EMPTY
    assert region(2, (0, 0), F(-4), F0) is HNRegion.S
    assert region(2, (1, 1), F(0), F0) is HNRegion.S       # decided by M3b
    assert region(2, (0, 0), F(1, 2), F0) is HNRegion.EMPTY


def test_P2_existence_verdict_is_total_over_grid():
    """P^2's existence boolean is total (DLP sharp everywhere), so the UNCLASSIFIED
    band never occurs there -- and region K is reserved for M3b, so it never occurs
    anywhere yet.  This asserts verdict totality only, NOT 'P^2 has no Kronecker
    generic-HN shapes' (it does -- E13 re-audit R2)."""
    for r in range(1, 6):
        for c1 in range(-6, 7):
            for c2 in range(0, 6):
                ch2 = F(c1 * c1, 2) - c2
                reg = region(r, (c1,), ch2, P2)
                assert reg in (HNRegion.S, HNRegion.EMPTY), (r, c1, ch2)


def test_hn_verdict_fields_are_consistent():
    """The HNVerdict dataclass exposes a consistent tri-state (exists / length / region)."""
    v = verdict(2, (1, 1), F(0), F0)                 # the flagship: M3b-decided
    assert isinstance(v, HNVerdict)
    assert v.exists is True and v.generic_hn_length == 1
    assert v.region is HNRegion.S
    assert v.factors == ((2, (1, 1), F(0)),)
    assert v.discriminant == F(1, 4) and v.sharp_bound == F(3, 4)
    v2 = verdict(1, (0,), F(-5), P2)                  # region S
    assert v2.exists is True and v2.generic_hn_length == 1 and v2.region is HNRegion.S


# --------------------------------------------------------------------------- #
# 8. Scope guards: M3a is del Pezzo e in {0,1} plus P^2 only.                  #
# --------------------------------------------------------------------------- #
def test_scope_e_ge_2_is_unlocked():
    """E13-M3c: any strictly ample F_e is in verdict scope -- the envelope's
    PROVEN branches bind and the M3b algorithm decides the rest (see
    tests/test_generic_hn.py for the pi-equivariance gate)."""
    assert semistable_exists(1, (F(0), F(0)), F(0), F2) is True     # O_{F_2}
    v = hn_verdict(2, (F(0), F(0)), F(1), F2)
    assert v.exists in (True, False)
    assert v.certificate.rigor.name == "PROVEN"


def test_scope_non_hirzebruch_raises():
    """A K3 (no F_e NS lattice) is refused -- NotImplementedError from hirzebruch_index."""
    with pytest.raises(NotImplementedError):
        semistable_exists(1, (F(0),), F(0), K3(2))


def test_non_integral_rank_is_never_truncated():
    """E13 re-audit R3: hn_verdict(r=3/2, ...) used to int()-truncate to r=1 and answer
    for a DIFFERENT character.  A non-integral rank is not the Chern character of any
    sheaf, so the verdict must be the invalid-character PROVEN_EMPTY -- and must differ
    from the truncated r=1 verdict wherever that one is region S."""
    # r=1, nu=(1,1), Delta=2 -> the truncated character (1,(1,1),...) is region S...
    assert hn_verdict(1, (F(1), F(1)), F(2), F0).region is HNRegion.S
    # ...but r=3/2 with the same (nu, Delta) is no character at all: certified EMPTY.
    v = hn_verdict(F(3, 2), (F(1), F(1)), F(2), F0)
    assert v.exists is False and v.region is HNRegion.EMPTY


def test_scope_nef_and_big_factory_raises():
    """The hirzebruch(1) FACTORY polarization H=(1,1) is only nef-and-big, not ample;
    it is honestly refused (a ValueError), directing the caller to
    hirzebruch_with_polarization with a strictly ample H."""
    with pytest.raises(ValueError):
        semistable_exists(1, (F(0), F(0)), F(0), hirzebruch(1))


# --------------------------------------------------------------------------- #
# 9. Oracle differential: package se on P^2 == independent reference.          #
# --------------------------------------------------------------------------- #
def test_oracle_differential_on_E12_corpus():
    """On the E12 falsification corpus, the package's P^2 semistable_exists agrees with
    the independent, package-free oracle reference_semistable_exists (which imports
    nothing from the package and uses no float)."""
    for row in CORPUS:
        assert row.surface == "P^2", row
        pkg = se(row.r, (row.c1,), row.ch2, P2)
        ref = reference_semistable_exists(row.r, row.c1, row.ch2)
        assert pkg is ref, (row.r, row.c1, row.ch2, "pkg", pkg, "ref", ref)


def test_oracle_differential_over_grid():
    """A differential sweep on P^2: package se == oracle over an integral-c2 box.

    Both compare the SAME full-NS discriminant against the SAME DLP curve (plus the
    (semi)exceptional disjunct), so they must agree on every class; a divergence names
    it.  Note se returns bool on P^2 (never None), matching the oracle's True/False."""
    for r in range(1, 7):
        for c1 in range(-8, 9):
            for c2 in range(0, 8):
                ch2 = F(c1 * c1, 2) - c2          # integral-c2 => oracle never INVALID
                pkg = se(r, (c1,), ch2, P2)
                ref = reference_semistable_exists(r, c1, ch2)
                assert pkg is ref, (r, c1, ch2, "pkg", pkg, "ref", ref)
