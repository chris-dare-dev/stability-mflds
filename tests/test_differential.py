"""Differential sweep: the package under test vs the independent oracle (E12-M0).

This is the one E12-M0 module that DOES import the package -- that is its whole
job.  It compares :mod:`bridgeland_stability` against ``tests/oracle/dlp_reference``
(which imports nothing from the package) on ``(status, Delta, delta)`` and pins the
six confirmed audit defects A1..A6 as ``xfail(strict=True)`` tripwires that flip when
their milestone lands.

The suite is GREEN at E12-M0 *because* A1..A6 are xfails: each asserts the
theorem-correct behaviour, which the package does not yet provide, so it is reported
``xfailed``.  When the fixing milestone lands, the assertion passes -> ``XPASS`` ->
strict-xfail failure -> the implementer must remove the marker (the "flip").

Faithful fast path
------------------
The package's P^2 verdict is, verbatim from ``nonemptiness_rational.moduli_nonempty``:

    nonempty = validate_character(xi)  AND  ( discriminant(xi, P2) >= delta_H(xi, P2, R_max)
               OR is_exceptional(xi)  OR  is_semiexceptional_p2(xi) )

(see ``_is_p2_exceptional`` / ``_is_p2_semiexceptional`` and ``validate_character``).  ``delta_H``
depends only on the slope, so we cache it per ``mu``.  Every in-box character has integral c2, so
``validate_character`` is always True there and the fast path drops it; the exhaustive box sweep
evaluates the remaining expression directly (fast), and a strided sub-sweep re-checks it against
the REAL ``moduli_nonempty`` so a future change to the wiring is caught here.
"""

import os
import random
import sys
from fractions import Fraction as F

import pytest

# --- oracle import (tests/oracle/; put tests/ on sys.path) ------------------- #
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _TESTS_DIR not in sys.path:
    sys.path.insert(0, _TESTS_DIR)

from oracle.dlp_reference import (
    Status,
    reference_nonempty,
    reference_delta,
    reference_discriminant,
    reference_is_exceptional,
    reference_is_semiexceptional,
)

# --- package under test (the differential module MAY import it) ------------- #
from bridgeland_stability.nonemptiness_rational import (
    moduli_nonempty,
    delta_H as pkg_delta_H,
    discriminant as pkg_discriminant,
    validate_character,
    SurfaceBundle,
    HNMode,
)
from bridgeland_stability.exceptional import (
    Bundle,
    is_exceptional as pkg_is_exceptional,
    is_semiexceptional_p2,
)
from bridgeland_stability.varieties import P2, P1xP1
from bridgeland_stability import dlp as pkg_dlp

# moduli_nonempty's default truncation.
_R_MAX = 60

# Box definition (E12-M0 acceptance criterion).
_BOX_R = range(1, 21)       # 1 <= r <= 20
_BOX_C1 = range(-40, 41)    # |c1| <= 40
_BOX_C2 = range(0, 61)      # 0 <= c2 <= 60

# Frozen divergence counts for the box (see the module tests below).  These pin the
# CURRENT (buggy) package; a fix changes the count and this assertion forces the
# fixing milestone to update it in lockstep with removing the matching xfail.
_FROZEN_A1 = 0    # E12-M2 flipped this from 99: the semiexceptional disjunct now lands, so the
                  # package matches the oracle on every box class (docs/CORRECTIONS.md sec.8 E12-M2)
_FROZEN_A2 = 0    # no epsilon-impostor survives the fixed is_exceptional (E12-M1 flipped this from 6)


# --------------------------------------------------------------------------- #
# Package helpers (faithful fast path + real verdict).                         #
# --------------------------------------------------------------------------- #
_pkg_delta_cache = {}


def _pkg_delta(mu):
    """Package ``delta_H`` at slope ``mu``, cached (delta_H is slope-only).

    Called exactly as ``moduli_nonempty`` calls it: ``delta_H(xi, P2, 60)``.
    """
    if mu not in _pkg_delta_cache:
        xi = SurfaceBundle(mu.denominator, (mu.numerator,), F(0))
        _pkg_delta_cache[mu] = pkg_delta_H(xi, P2, _R_MAX)
    return _pkg_delta_cache[mu]


def _pkg_fast(r, c1, ch2):
    """The package's P^2 verdict via its documented wiring; returns (disc, dH, status, exc).

    Box characters are all integral-c2, so ``validate_character`` is always True there and is a
    no-op in this fast path (the box never produces an INVALID class); it is exercised through the
    REAL ``moduli_nonempty`` by ``test_box_fast_path_matches_real_moduli_nonempty``.
    """
    disc = pkg_discriminant(SurfaceBundle(r, (c1,), ch2), P2)
    dH = _pkg_delta(F(c1, r))
    exc = pkg_is_exceptional(Bundle(r, c1, ch2))
    semi = is_semiexceptional_p2(r, c1, ch2)            # E12-M2 (A1): m>=1 disjunct
    status = Status.NONEMPTY if (disc >= dH or exc or semi) else Status.EMPTY
    return disc, dH, status, exc


def _pkg_real_status(r, c1, ch2):
    return Status.NONEMPTY if moduli_nonempty(r, (c1,), ch2, P2).nonempty else Status.EMPTY


# --------------------------------------------------------------------------- #
# 4a. Exhaustive box sweep.                                                    #
# --------------------------------------------------------------------------- #
def test_box_delta_and_discriminant_agree_everywhere():
    """Delta and delta agree on EVERY integral character in the box.

    Both compute the same CH discriminant, and every in-box slope has denominator
    <= 20 <= R_max = 60, so the package's truncated delta_H equals the oracle's
    exact DLP delta.  A single disagreement names the class.
    """
    for r in _BOX_R:
        for c1 in _BOX_C1:
            mu = F(c1, r)
            odelta = reference_delta(mu)
            pdelta = _pkg_delta(mu)
            assert pdelta == odelta, ("delta", r, c1, mu, "pkg", pdelta, "ref", odelta)
            for c2 in _BOX_C2:
                ch2 = F(c1 * c1, 2) - c2  # integral-c2 character
                disc = pkg_discriminant(SurfaceBundle(r, (c1,), ch2), P2)
                assert disc == reference_discriminant(r, c1, ch2), (
                    "disc", r, c1, ch2, "pkg", disc, "ref", reference_discriminant(r, c1, ch2)
                )


def test_box_status_divergences_are_exactly_A1_and_A2():
    """Package vs oracle non-emptiness over the whole box.

    Every divergence is attributed to A1 (missing semiexceptional disjunct) or A2
    (epsilon-impostor accepted as exceptional); an unattributable divergence names
    the class and fails.  The A1/A2 counts are frozen so a fix in E12-M1/M2 changes
    the count and forces this test (and the matching xfail) to be updated.
    """
    n_a1 = n_a2 = 0
    unexpected = []
    for r in _BOX_R:
        for c1 in _BOX_C1:
            mu = F(c1, r)
            odelta = reference_delta(mu)   # oracle delta once per slope (naive; reused)
            for c2 in _BOX_C2:
                ch2 = F(c1 * c1, 2) - c2
                disc, dH, pkg, exc = _pkg_fast(r, c1, ch2)
                # Oracle verdict, reusing odelta (integral c2 in-box => never INVALID).
                # This is reference_nonempty's exact composition; the composed public
                # function is pinned separately by tests/test_oracle_integrity.py.
                if reference_discriminant(r, c1, ch2) >= odelta:
                    ref = Status.NONEMPTY
                elif reference_is_semiexceptional(r, c1, ch2):
                    ref = Status.NONEMPTY
                else:
                    ref = Status.EMPTY
                if pkg is ref:
                    continue
                if ref is Status.NONEMPTY and pkg is Status.EMPTY:
                    # Package missed a non-empty class.  Below the curve (disc < delta),
                    # so the only theorem branch it could be is semiexceptionality (A1).
                    assert disc < odelta, (r, c1, ch2)
                    assert reference_is_semiexceptional(r, c1, ch2), ("not semiexc", r, c1, ch2)
                    assert not reference_is_exceptional(r, c1, ch2), (r, c1, ch2)
                    n_a1 += 1
                elif ref is Status.EMPTY and pkg is Status.NONEMPTY:
                    # Package called it non-empty; oracle says empty.  Must be an
                    # impostor accepted as exceptional (A2).
                    assert exc and not reference_is_exceptional(r, c1, ch2), ("not impostor", r, c1, ch2)
                    n_a2 += 1
                else:
                    unexpected.append((r, c1, ch2, "ref", ref, "pkg", pkg))
    assert unexpected == [], f"UNEXPECTED divergence (possible new bug): {unexpected[:20]}"
    assert (n_a1, n_a2) == (_FROZEN_A1, _FROZEN_A2), (
        f"divergence counts changed: A1={n_a1} (frozen {_FROZEN_A1}), "
        f"A2={n_a2} (frozen {_FROZEN_A2}); update the frozen counts and the matching xfail"
    )


def test_box_fast_path_matches_real_moduli_nonempty():
    """The documented fast-path wiring equals the REAL moduli_nonempty on a strided
    sub-sweep; a future change to the verdict wiring is caught here."""
    stride = 29
    i = 0
    for r in _BOX_R:
        for c1 in _BOX_C1:
            for c2 in _BOX_C2:
                i += 1
                if i % stride:
                    continue
                ch2 = F(c1 * c1, 2) - c2
                _, _, fast, _ = _pkg_fast(r, c1, ch2)
                real = _pkg_real_status(r, c1, ch2)
                assert fast is real, (r, c1, ch2, "fast", fast, "real", real)


# --------------------------------------------------------------------------- #
# 4b. Seeded high-rank sweep (reaches rank >= 610) + truncation (A4).          #
# --------------------------------------------------------------------------- #
def _classify_divergence(r, c1, ch2, disc, dH, exc, ref, pkg, odelta):
    """Return 'A1'/'A2'/'A4' for an expected divergence, else None."""
    if ref is Status.NONEMPTY and pkg is Status.EMPTY:
        if reference_is_semiexceptional(r, c1, ch2):
            return "A1"
    if ref is Status.EMPTY and pkg is Status.NONEMPTY:
        if exc and not reference_is_exceptional(r, c1, ch2):
            return "A2"
        # Above the package's TRUNCATED curve but below the TRUE curve, not
        # exceptional: the R_max=60 cutoff (A4).
        if (not exc) and disc >= dH and disc < odelta:
            return "A4"
    return None


def test_high_rank_random_sweep_reaches_610_and_only_A1_A2_A4():
    """Deterministic (seed 0) high-rank sweep.

    * Delta always agrees.
    * delta agrees exactly when denominator(mu) <= 60; for larger denominators the
      package's R_max=60 enumeration sees a SUBSET of exceptional slopes, so its
      delta_H can only UNDER-estimate the true delta (A4 truncation).
    * Every status divergence is A1/A2/A4-shaped; an unclassified one names the class.
    """
    rng = random.Random(0)
    saw_high_denom = False
    saw_a4 = False
    for _ in range(400):
        r = rng.randint(1, 700)
        c1 = rng.randint(-2 * r, 2 * r)
        c2 = rng.randint(0, 300)
        ch2 = F(c1 * c1, 2) - c2
        mu = F(c1, r)
        disc, dH, pkg, exc = _pkg_fast(r, c1, ch2)
        assert reference_discriminant(r, c1, ch2) == disc, (r, c1, ch2)
        odelta = reference_delta(mu)
        if mu.denominator <= 60:
            assert dH == odelta, ("delta", r, c1, mu, dH, odelta)
        else:
            saw_high_denom = True
            assert dH <= odelta, ("truncation must under-estimate", r, c1, mu, dH, odelta)
        ref = reference_nonempty(r, c1, ch2)   # integral c2 => never INVALID
        if pkg is not ref:
            label = _classify_divergence(r, c1, ch2, disc, dH, exc, ref, pkg, odelta)
            assert label is not None, ("UNCLASSIFIED divergence", r, c1, ch2, ref, pkg)
            saw_a4 = saw_a4 or label == "A4"
    assert saw_high_denom, "sweep never reached a slope of denominator > 60"


# Rank-610/985/1325 epsilon-impostors: their rank IS a Markov number and their ch2
# matches the exceptional formula, yet the slope is not in the epsilon-image, so the
# oracle rejects them.  (The package accepts them -- that is A2, pinned in test_A2.)
_IMPOSTOR_SLOPES = [
    (133, 610), (477, 610), (183, 985), (802, 985), (182, 1325), (1143, 1325),
]


def test_reference_rejects_markov_rank_impostors():
    for p, q in _IMPOSTOR_SLOPES:
        ch2 = F(p * p - q * q + 1, 2 * q)          # the exceptional ch2 at slope p/q
        assert reference_is_exceptional(q, p, ch2) is False, (p, q)


def test_A4_cutoff_now_captures_rank_89_cusp():
    """The concrete A4 witness, post-fix: mu = 34/89, denominator 89 > the old R_max = 60.

    E12-M1 bumps the P^2 enumeration cutoff to R_max = max(R_max, denominator(mu)), so the
    package now sees the rank-89 cusp and its delta_H equals the oracle's exact delta (the
    STRONGER corrected value, pinned here; see docs/CORRECTIONS.md sec.8).  The oracle's
    verdict for (8010,3060,-3421) is still EMPTY (Delta = delta - 1/712890), now reproduced
    by the package (test_A4_truncation_flips_to_empty)."""
    r, c1, ch2 = 8010, 3060, F(-3421)
    mu = F(c1, r)
    assert mu == F(34, 89)
    odelta = reference_delta(mu)
    assert odelta == F(3961, 7921)
    assert _pkg_delta(mu) == odelta                 # E12-M1: R_max bumps to denom(mu)=89, cusp now seen
    # And the oracle's verdict is EMPTY (Delta = delta - 1/712890 < delta).
    assert reference_nonempty(r, c1, ch2) is Status.EMPTY
    assert reference_discriminant(r, c1, ch2) == F(356489, 712890)


# --------------------------------------------------------------------------- #
# 4c. dlp.moduli_nonempty cross-check (secondary).                            #
# --------------------------------------------------------------------------- #
_DLP_CHECK = [
    (4, 2, F(-1)),        # A1: 2*T(-1) semiexceptional
    (2, 0, F(0)),         # A1: O^2 semiexceptional
    (10, 3, F(-9, 2)),    # A2: 3/10 impostor
    (2, 1, F(-5, 2)),     # above curve
    (2, 1, F(-1, 2)),     # T(-1) exceptional
    (5, 2, F(-2)),        # 2/5 exceptional
    (1, 0, F(0)),         # O exceptional
    (2, 1, F(1, 2)),      # empty
    (1, 0, F(-1)),        # P^2[1] boundary
    (1, 0, F(-5)),        # P^2[5]
]


def test_dlp_moduli_nonempty_cross_check():
    """dlp.moduli_nonempty now matches the oracle on every class in _DLP_CHECK.

    E12-M1 closed A2 (impostor rejection) and E12-M2 closed A1 (the semiexceptional disjunct)
    and A4b (the R_max cutoff) in the P^2-only twin, so both former divergence classes are gone:
    the once-permissive A1/A2 tripwire is now a hard equality check (docs/CORRECTIONS.md sec.8)."""
    for r, c1, ch2 in _DLP_CHECK:
        res = pkg_dlp.moduli_nonempty(Bundle(r, c1, ch2))
        pkg = Status.NONEMPTY if res["nonempty"] else Status.EMPTY
        ref = reference_nonempty(r, c1, ch2)
        assert ref is not Status.INVALID   # this list is all valid characters
        assert pkg is ref, (r, c1, ch2, "ref", ref, "pkg", pkg)


def test_box_dlp_moduli_nonempty_matches_reference():
    """The dlp twin, strided over the box, equals the oracle on every sampled class.

    A4b existed because the differential sweep only ever exercised nonemptiness_rational; this
    strided sweep pulls dlp.moduli_nonempty into the same lens.  The box is all integral-c2, so
    the oracle never returns INVALID, and (post-M2) both disjuncts + the cutoff agree with it."""
    stride, i = 29, 0
    for r in _BOX_R:
        for c1 in _BOX_C1:
            for c2 in _BOX_C2:
                i += 1
                if i % stride:
                    continue
                ch2 = F(c1 * c1, 2) - c2
                pkg = (Status.NONEMPTY if pkg_dlp.moduli_nonempty(Bundle(r, c1, ch2))["nonempty"]
                       else Status.EMPTY)
                ref = reference_nonempty(r, c1, ch2)      # box is integral-c2 => never INVALID
                assert pkg is ref, (r, c1, ch2, "pkg", pkg, "ref", ref)


# --------------------------------------------------------------------------- #
# 4d. The A1..A6 strict-xfail tripwires (each names the milestone that flips it). #
# --------------------------------------------------------------------------- #
def test_A1_semiexceptional_nonempty():
    # xi = 2*ch(T(-1)); Delta=3/8 < delta(1/2)=5/8 but semiexceptional => NONEMPTY.  (E12-M2:
    # xfail removed -- the P^2 semiexceptional disjunct now lands, so this passes normally.)
    assert moduli_nonempty(4, (2,), F(-1), P2).nonempty is True


def test_A2_impostor_not_exceptional():
    # chi(E,E)=1 and c2 integral, but the slope is not an epsilon-slope.  (E12-M1: xfail
    # removed -- the fixed is_exceptional now rejects these, so this passes normally.)
    assert pkg_is_exceptional(Bundle(10, 3, F(-9, 2))) is False
    assert pkg_is_exceptional(Bundle(610, 133, F(-581, 2))) is False   # rank IS Markov


def test_A3_invalid_character_is_empty():
    # chi = -1/2 not in Z (c2 = 3/2): not a Chern character.  (E12-M2: xfail removed --
    # validate_character now rejects non-integral chi, so this passes normally.)
    assert moduli_nonempty(1, (0,), F(-3, 2), P2).nonempty is False


def test_A4_truncation_flips_to_empty():
    # mu = 34/89; Delta = delta(34/89) - 1/712890 < delta, not (semi)exceptional.  (E12-M1:
    # xfail removed -- R_max now bumps to denom(mu)=89, so the package reports EMPTY.)
    assert moduli_nonempty(8010, (3060,), F(-3421), P2).nonempty is False


def test_A4b_dlp_truncation_flips_to_empty():
    # A4b (E12-M2): dlp.moduli_nonempty must carry E12-M1's cutoff too.  R_max bumps to
    # denom(34/89)=89, so the rank-89 cusp is seen and delta = 3961/7921 = 356490/712890 >
    # Delta = 356489/712890 -> EMPTY.  Delivered as a direct passing test (like A2/A4), not an xfail.
    assert pkg_dlp.moduli_nonempty(Bundle(8010, 3060, F(-3421)))["nonempty"] is False


def test_A5_forged_target_rejected():
    # E12-M4: a certified target != the package's own certified sharp bound is a forge ->
    # ValueError on EVERY surface, PAPER and ORACLE alike.  native delta_H(2,(0,0),-4|P1xP1)=1,
    # native delta(0|P2)=1 != 0.  A raw ORACLE target is refused outright (an ORACLE datum is a
    # capability object minted after a real construction, never a bare pair).
    with pytest.raises(ValueError):                                 # original off-P^2 ORACLE forge
        moduli_nonempty(2, (0, 0), F(-4), P1xP1,
                        delta_H_target=F(10 ** 6), hn_source=HNMode.ORACLE)
    with pytest.raises(ValueError):                                 # P^2 PAPER forge (delta(0)=1 != 0)
        moduli_nonempty(3, (0,), F(-2), P2, delta_H_target=F(0), hn_source=HNMode.PAPER)
    with pytest.raises(ValueError):                                 # P^2 ORACLE forge (raw ORACLE refused)
        moduli_nonempty(3, (0,), F(-2), P2, delta_H_target=F(0), hn_source=HNMode.ORACLE)


def test_A6_F0_ch2_guard_empty():
    # (3,(1,1),0) on F_0: nu=(1/3,1/3), Delta=1/9; emptiness_bound 5/9 > 1/9 => EMPTY.  (E12-M3:
    # xfail removed -- _hirzebruch_verdict now guards `exceptional` on ch2 == exceptional_ch2
    # (3,(1,1)) = -1 != 0, so the (r,c1)-only shortcut no longer forges non-empty -> PROVEN_EMPTY.)
    assert moduli_nonempty(3, (1, 1), F(0), P1xP1).nonempty is False
