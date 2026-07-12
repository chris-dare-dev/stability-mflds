"""E13-M1 / G18: the Coskun-Huizenga F_e -> F_{e-2} reduction map pi + Lemma 11.3.

Primary source: arXiv:1907.06739 (Coskun-Huizenga, "Existence of semistable sheaves on
Hirzebruch surfaces") Sec. 11.1 and Lemma 11.3.

The map, in the package's ``(f, s)`` basis (``f = F`` the fiber, ``s = E`` the section), is the
unimodular ``SL_2(Z)`` NS isometry ``M = [[1,-1],[0,1]]`` acting on ``c1``, fixing ``r`` and ``ch2``:

    pi(r, (x, y), ch2) = (r, (x - y, y), ch2),        M^T G_{e-2} M = G_e.

Every value below is hand-derived AND independently confirmed against the package's own exact
``Fraction`` invariants (:mod:`bridgeland_stability.dlp_hirzebruch`,
:func:`bridgeland_stability.exceptional_surface.chi`).  Two worked characters anchor the pins:

    (2,(3,1),-1) / F_2, H=A_1=(3,1), d=4   ->  (2,(2,1),-1) / F_0, H'=A'_1=(2,1), d=4
        Delta=1, chi(O,.)=4, mu_H=1/2, hilbert_P(nu)=3, c2=3  (all equal on both sides)
    (2,(2,1),-1/2) / F_3, H=(4,1), d=5     ->  (2,(1,1),-1/2) / F_1, H'=(3,1), d=5
        Delta=3/8 (= Delta_V(2)), c2=1

and the telescope (3,(3,1),0) / F_4 H=(5,1) -> F_2 -> F_0 has Delta=1/9 at both ends.
"""

import math
import os
import sys
from fractions import Fraction as F

import pytest

from bridgeland_stability.reduction import (
    reduce,
    reduce_to_del_pezzo,
    pi_c1,
    REDUCTION_MATRIX,
)
from bridgeland_stability.nonemptiness_rational import hirzebruch_with_polarization
from bridgeland_stability.dlp_hirzebruch import (
    hirzebruch_index,
    discriminant,
    total_slope,
    hilbert_P,
    dlp_envelope,
    emptiness_bound,
    exceptional_discriminant,
    is_del_pezzo_anticanonical,
)
from bridgeland_stability.exceptional_surface import (
    SurfaceBundle,
    canonical_class,
    chi as surface_chi,
)
from bridgeland_stability.varieties import P2, P1xP1, hirzebruch, NSLattice

# --- oracle import (tests/oracle/; put tests/ on sys.path) ------------------- #
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _TESTS_DIR not in sys.path:
    sys.path.insert(0, _TESTS_DIR)

from oracle.dlp_reference import reference_reduce_pi


# --------------------------------------------------------------------------- #
# Small exact helpers (all Fraction / int; no float ever).                     #
# --------------------------------------------------------------------------- #
def _O(surface):
    """The structure sheaf ch(O_X) = (1, (0,...,0), 0) on ``surface``."""
    return SurfaceBundle(1, (0,) * surface.lattice.rank, F(0))


def _mu_H(xi, surface):
    """The A_m-slope mu_H = <c1, H>/(r d) (exact Fraction)."""
    return surface.lattice.pairing(xi.c1, surface.H) / (xi.r * surface.d)


def _c2(xi, surface):
    """c2 = 1/2 <c1, c1> - ch2 (exact Fraction)."""
    return F(1, 2) * surface.lattice.self_pairing(xi.c1) - xi.ch2


def _Fe(e, H):
    """A strictly-ample F_e = (F_e, H)."""
    return hirzebruch_with_polarization(e, H)


# A sweep of ample F_e over both parities: H=(e+1, 1) (Nakai: a-e b = 1 > 0).
_SWEEP_ES = (2, 3, 4, 5, 6)


def _ample_Fe(e):
    return _Fe(e, (e + 1, 1))


# --------------------------------------------------------------------------- #
# The map + the matrix + the isometry.                                         #
# --------------------------------------------------------------------------- #
def test_reduction_matrix_is_unimodular():
    """M = [[1,-1],[0,1]] is in SL_2(Z) (det = 1)."""
    (a, b), (c, d) = REDUCTION_MATRIX
    assert (a, b, c, d) == (1, -1, 0, 1)
    assert a * d - b * c == 1


def test_pi_c1_matches_matrix_action():
    """pi(x, y) = (x - y, y) is the column action M . (x, y), exact and type-preserving."""
    assert pi_c1((5, 2)) == (3, 2)
    assert pi_c1((F(3, 2), F(1, 2))) == (F(1), F(1, 2))
    # matrix action agrees
    (m00, m01), (m10, m11) = REDUCTION_MATRIX
    for x, y in [(3, 1), (-2, 4), (7, 7)]:
        assert pi_c1((x, y)) == (m00 * x + m01 * y, m10 * x + m11 * y)
    with pytest.raises(ValueError):
        pi_c1((1, 2, 3))


def test_pi_is_an_isometry_M_T_G_M_equals_G():
    """M^T G_{e-2} M = G_e for both parities: pi is an NS isometry F_e -> F_{e-2}
    (arXiv:1907.06739 Lemma 11.3(1))."""
    (m00, m01), (m10, m11) = REDUCTION_MATRIX
    M = [[m00, m01], [m10, m11]]
    MT = [[M[j][i] for j in range(2)] for i in range(2)]

    def mul(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

    for e in _SWEEP_ES:
        G_e = [[0, 1], [1, -e]]
        G_em2 = [[0, 1], [1, -(e - 2)]]
        assert mul(mul(MT, G_em2), M) == G_e


# --------------------------------------------------------------------------- #
# reduce(): the two worked characters (table in the module docstring).         #
# --------------------------------------------------------------------------- #
def test_reduce_worked_character_1_F2_to_F0():
    """(2,(3,1),-1)/F_2 H=(3,1) -> (2,(2,1),-1)/F_0 H=(2,1); all Lemma 11.3 invariants pinned."""
    S2 = _Fe(2, (3, 1))
    v = SurfaceBundle(2, (3, 1), F(-1))
    vr, S0 = reduce(v, S2)

    # the map itself: r, ch2 fixed; c1 -> (x-y, y); index drops by 2; H -> A'_1 = (2,1)
    assert (vr.r, vr.c1, vr.ch2) == (2, (F(2), F(1)), F(-1))
    assert hirzebruch_index(S2) == 2 and hirzebruch_index(S0) == 0
    assert S0.H == (2, 1) and S2.d == 4 and S0.d == 4
    assert S0.lattice.gram == ((0, 1), (1, 0))            # genuine rank-2 F_0 = P^1xP^1 lattice

    # (2) discriminant, (4) chi(O,.), (5) mu_H + hilbert_P, and c2 -- all equal both sides
    assert discriminant(v, S2) == discriminant(vr, S0) == F(1)
    assert surface_chi(_O(S2), v, S2) == surface_chi(_O(S0), vr, S0) == 4
    assert _mu_H(v, S2) == _mu_H(vr, S0) == F(1, 2)
    assert hilbert_P(total_slope(v), S2) == hilbert_P(total_slope(vr), S0) == F(3)
    assert _c2(v, S2) == _c2(vr, S0) == F(3)


def test_reduce_worked_character_2_F3_to_F1():
    """(2,(2,1),-1/2)/F_3 H=(4,1) -> (2,(1,1),-1/2)/F_1 H=(3,1); Delta=3/8=Delta_V(2), c2=1."""
    S3 = _Fe(3, (4, 1))
    w = SurfaceBundle(2, (2, 1), F(-1, 2))
    wr, S1 = reduce(w, S3)

    assert (wr.r, wr.c1, wr.ch2) == (2, (F(1), F(1)), F(-1, 2))
    assert hirzebruch_index(S1) == 1 and S1.H == (3, 1)
    assert S3.d == 5 and S1.d == 5
    assert discriminant(w, S3) == discriminant(wr, S1) == F(3, 8) == exceptional_discriminant(2)
    assert _c2(w, S3) == _c2(wr, S1) == F(1)


def test_reduce_fixes_r_and_ch2_and_stays_exact():
    """r and ch2 are untouched; every coordinate is an exact Fraction (no float)."""
    S4 = _Fe(4, (5, 1))
    xi = SurfaceBundle(3, (3, 1), F(7, 3))
    xir, S2 = reduce(xi, S4)
    assert xir.r == xi.r == 3
    assert xir.ch2 == xi.ch2 == F(7, 3)
    assert all(isinstance(c, F) for c in xir.c1)
    assert isinstance(discriminant(xir, S2), F)


# --------------------------------------------------------------------------- #
# Guards: e < 2, non-ample H, non-F_e surface.                                 #
# --------------------------------------------------------------------------- #
def test_reduce_raises_for_e_below_2():
    """F_0 / F_1 are already del Pezzo -- nothing to reduce (ValueError)."""
    for e in (0, 1):
        S = _Fe(e, (e + 2, 1))
        with pytest.raises(ValueError):
            reduce(SurfaceBundle(2, (1, 1), F(0)), S)


def test_reduce_raises_for_nonample_H():
    """The nef-and-big factory polarization H=(n,1) is refused; reduce needs strict ampleness."""
    with pytest.raises(ValueError):
        reduce(SurfaceBundle(2, (1, 1), F(0)), hirzebruch(3))


def test_reduce_raises_off_hirzebruch():
    """A non-F_e surface has no fiber/section Gram -> NotImplementedError (via hirzebruch_index)."""
    with pytest.raises(NotImplementedError):
        reduce(SurfaceBundle(2, (1,), F(0)), P2)


# --------------------------------------------------------------------------- #
# Lemma 11.3 (1)-(6), each as an exact identity over a sweep + both parities.  #
# --------------------------------------------------------------------------- #
def test_lemma_113_1_pairing_preserved():
    """(1) <pi u, pi v>_{e-2} = <u, v>_e.  Pinned example u=(3,1),v=(1,2),e=2 -> 3."""
    S2, S0 = _Fe(2, (3, 1)), _Fe(0, (2, 1))
    assert S2.lattice.pairing((3, 1), (1, 2)) == S0.lattice.pairing(pi_c1((3, 1)), pi_c1((1, 2))) == 3
    for e in _SWEEP_ES:
        Se, Sem2 = _ample_Fe(e), _ample_Fe(e - 2)   # pairing is H-independent; any lattice rep
        for u in [(1, 0), (0, 1), (2, -3), (-1, 4), (5, 5)]:
            for v in [(1, 1), (3, -2), (0, 2)]:
                assert Se.lattice.pairing(u, v) == Sem2.lattice.pairing(pi_c1(u), pi_c1(v))


def test_lemma_113_2_discriminant_preserved():
    """(2) Delta(pi xi) = Delta(xi) over a full integral sweep, both parities."""
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        for r in (1, 2, 3):
            for x in range(-3, 4):
                for y in range(-3, 4):
                    xi = SurfaceBundle(r, (x, y), F(1, 2))
                    xir, Sem2 = reduce(xi, Se)
                    assert discriminant(xir, Sem2) == discriminant(xi, Se)


def test_lemma_113_3_pi_of_canonical_and_structure_sheaf():
    """(3) pi(K_{F_e}) = K_{F_{e-2}} = (-e,-2), pi(-K ray) = (e,2), pi(ch O) = ch O; K^2=8."""
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        Sem2 = reduce(_O(Se), Se)[1]
        # canonical class transports to the reduced canonical class K_{F_{e-2}} = (-e, -2)
        assert tuple(pi_c1(canonical_class(Se))) == tuple(canonical_class(Sem2)) == (-e, -2)
        assert canonical_class(Se) == (-(e + 2), -2)
        # the -K ray (e+2, 2) -> (e, 2)
        assert pi_c1((e + 2, 2)) == (e, 2)
        # K^2 = 8 on every F_e (so ch(O(K)) has ch2 = 4, transported unchanged)
        assert Se.lattice.self_pairing(canonical_class(Se)) == 8
        # ch(O) is fixed
        Or, _ = reduce(_O(Se), Se)
        assert (Or.r, Or.c1, Or.ch2) == (1, (F(0), F(0)), F(0))


def test_lemma_113_4_euler_pairings_preserved():
    """(4) chi(v) := chi(O_X, v) and chi(v, w), chi(w, v) all preserved (RR is isometry-built).

    v, w are genuine integral K-theory classes on each F_e (c2 = 0, i.e. ch2 = 1/2<c1,c1>, which is
    isometry-preserved so the reductions are genuine too); chi() would raise on a non-class.
    """
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        lat = Se.lattice
        cv, cw = (3, 1), (1, 2)                       # genuine section-carrying c1 (transport is nontrivial)
        v = SurfaceBundle(2, cv, F(1, 2) * lat.self_pairing(cv))    # c2 = 0
        w = SurfaceBundle(3, cw, F(1, 2) * lat.self_pairing(cw))    # c2 = 0
        vr, Sem2 = reduce(v, Se)
        wr, _ = reduce(w, Se)
        # chi(O, v) and the two-argument Euler forms are all preserved under pi
        assert surface_chi(_O(Se), v, Se) == surface_chi(_O(Sem2), vr, Sem2)
        assert surface_chi(v, w, Se) == surface_chi(vr, wr, Sem2)
        assert surface_chi(w, v, Se) == surface_chi(wr, vr, Sem2)


def test_lemma_113_4_integral_and_primitive_preserved():
    """(4) M unimodular: integral -> integral, and gcd primitivity is preserved (e = 2..6)."""
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        for x in range(-4, 5):
            for y in range(-4, 5):
                xi = SurfaceBundle(2, (x, y), F(0))
                xir, _ = reduce(xi, Se)
                # image is integral
                assert all(c.denominator == 1 for c in xir.c1)
                px, py = int(xir.c1[0]), int(xir.c1[1])
                # gcd(x - y, y) == gcd(x, y): primitive -> primitive
                assert math.gcd(abs(px), abs(py)) == math.gcd(abs(x), abs(y))


def test_lemma_113_5_slope_and_hilbert_transport():
    """(5) A_m -> A'_m: mu_H(v) = mu_{H'}(pi v) and hilbert_P(nu) = hilbert_P(pi nu)."""
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        # A_m identity: pi(A_m) = A'_m.  A_m = -1/2 K + m F = ((e+2)/2 + m, 1); H=(e+1,1) is A_{m},
        # m = (e+1) - (e+2)/2 = e/2; pi(H) = (e, 1) = ((e-2+2)/2 + e/2, 1) = A'_{e/2} on F_{e-2}.
        assert pi_c1(Se.H) == (e, 1)
        for r in (1, 2):
            for c1 in [(2, 1), (3, 2), (-1, 1)]:
                xi = SurfaceBundle(r, c1, F(-1))
                xir, Sem2 = reduce(xi, Se)
                assert _mu_H(xi, Se) == _mu_H(xir, Sem2)
                assert hilbert_P(total_slope(xi), Se) == hilbert_P(total_slope(xir), Sem2)


def test_lemma_113_6_additive_on_direct_sums():
    """(6) pi is additive: pi(A (+) B (+) C) = pi(A) (+) pi(B) (+) pi(C).

    The theorem-backed proxy for the paper's named Sec.11.1(6) three-term direct-sum character
    (the exact named example is open question O1; additivity is what the paper's identity rests on).
    """
    S2 = _Fe(2, (3, 1))
    A = SurfaceBundle(1, (1, 0), F(0))
    B = SurfaceBundle(2, (3, 1), F(-1))
    C = SurfaceBundle(1, (0, 1), F(0))

    def dsum(*xs):
        r = sum(x.r for x in xs)
        c1 = tuple(sum(x.c1[i] for x in xs) for i in range(2))
        ch2 = sum(x.ch2 for x in xs)
        return SurfaceBundle(r, c1, ch2)

    lhs, _ = reduce(dsum(A, B, C), S2)
    parts = [reduce(x, S2)[0] for x in (A, B, C)]
    rhs = dsum(*parts)
    assert (lhs.r, lhs.c1, lhs.ch2) == (rhs.r, rhs.c1, rhs.ch2)
    # and the named example is a genuine direct sum: r and ch2 add, c1 adds componentwise.
    # A(+)B(+)C = (4, (1+3+0, 0+1+1), -1) = (4, (4,2), -1); pi(4,2) = (4-2, 2) = (2,2).
    assert lhs.r == 4 and lhs.c1 == (F(2), F(2)) and lhs.ch2 == F(-1)


def test_lemma_113_6_named_direct_sum_character():
    """(6) FAITHFUL to the paper's *named* character (arXiv:1907.06739 Sec.11.1(6)), closing O1.

    pi carries  O(-E+(n-1)F)^A (+) O^B (+) O(-F)^C  to  O(-E'+nF')^A (+) O'^B (+) O(-F')^C.
    In the (f,s) basis (f=F, s=E):  O(-E+(n-1)F) has c1=(n-1,-1) -> pi=(n,-1)=O(-E'+nF');
    O(-F)=(-1,0) is fixed = O(-F');  O=(0,0) is fixed.  This pins WHICH line bundle each summand
    maps to (stronger than generic additivity), including the exact ch2 that the isometry forces.
    """
    def named_sum(surf, twist_f, A, B, C):
        # twist_f is the fiber(f)-coefficient of the twisted summand: n-1 on F_e, n on F_{e-2}.
        L1 = SurfaceBundle.line_bundle((twist_f, -1), surf)   # O(-E + twist_f * F)
        L2 = SurfaceBundle.line_bundle((0, 0), surf)          # O
        L3 = SurfaceBundle.line_bundle((-1, 0), surf)         # O(-F)
        r = A * L1.r + B * L2.r + C * L3.r
        c1 = tuple(A * L1.c1[i] + B * L2.c1[i] + C * L3.c1[i] for i in range(2))
        ch2 = A * L1.ch2 + B * L2.ch2 + C * L3.ch2
        return SurfaceBundle(r, c1, ch2)

    for e, H, n, A, B, C in [(2, (3, 1), 2, 1, 1, 1),
                             (4, (5, 1), 3, 2, 1, 1),
                             (2, (5, 2), 0, 1, 2, 3)]:
        Se = _Fe(e, H)
        src = named_sum(Se, n - 1, A, B, C)                   # O(-E+(n-1)F)^A (+) O^B (+) O(-F)^C
        Sem2 = hirzebruch_with_polarization(e - 2, pi_c1(H))
        tgt = named_sum(Sem2, n, A, B, C)                     # O(-E'+nF')^A (+) O'^B (+) O(-F')^C
        red, _ = reduce(src, Se)
        assert (red.r, red.c1, red.ch2) == (tgt.r, tgt.c1, tgt.ch2), (e, H, n, A, B, C)


# --------------------------------------------------------------------------- #
# reduce_to_del_pezzo: telescoping to F_0 / F_1, discriminant invariant.       #
# --------------------------------------------------------------------------- #
def test_reduce_to_del_pezzo_telescope_even():
    """(3,(3,1),0)/F_4 H=(5,1) -> F_2 -> F_0; Delta = 1/9 at both ends (property (2) telescopes)."""
    S4 = _Fe(4, (5, 1))
    t = SurfaceBundle(3, (3, 1), F(0))
    tdp, Sdp = reduce_to_del_pezzo(t, S4)
    assert hirzebruch_index(Sdp) == 0
    assert tdp.c1 == (F(1), F(1))
    assert discriminant(t, S4) == discriminant(tdp, Sdp) == F(1, 9)


def test_reduce_to_del_pezzo_lands_on_correct_parity():
    """Even e -> F_0, odd e -> F_1; the composed discriminant always equals Delta(xi)."""
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        xi = SurfaceBundle(2, (2, 1), F(-1, 3))
        xdp, Sdp = reduce_to_del_pezzo(xi, Se)
        assert hirzebruch_index(Sdp) == (0 if e % 2 == 0 else 1)
        assert discriminant(xdp, Sdp) == discriminant(xi, Se)


def test_reduce_to_del_pezzo_is_noop_on_del_pezzo():
    """e in {0, 1}: nothing to reduce -- the character and surface are returned unchanged."""
    for e in (0, 1):
        Se = _Fe(e, (e + 2, 1))
        xi = SurfaceBundle(2, (1, 1), F(0))
        xdp, Sdp = reduce_to_del_pezzo(xi, Se)
        assert xdp is xi and Sdp is Se


# --------------------------------------------------------------------------- #
# Cross-check against the E11-M6 envelope (SCOPED at M1; see the docstring).    #
# --------------------------------------------------------------------------- #
def test_cross_check_envelope_lower_bound_equality():
    """The reduced character's DLP envelope on F_{e-2} equals the direct one on F_e, and BOTH are
    non-sharp lower bounds -- the theorem-backed EQUALITY that realizes ``lower_bound <= reduced``.

    Honest scope (E13-M1).  Reducing a strictly ample F_e (e >= 2) can never reach an anticanonical
    (sharp) F_0/F_1 ray (pi(H) ~ -K_{F_{e-2}} iff H ~ -K_{F_e}, and -K_{F_e} is not ample for
    e >= 2), so both envelopes have ``sharp=False`` and the reduced value is a certified lower bound
    EQUAL to the direct one, not a strictly-sharper delta_H.  The 'strictly sharper delta_H /
    certified_sharp on F_0/F_1' sub-clause needs the sharp non-anticanonical delta_H theory and is
    DEFERRED to E13-M2/M3 (docs/CORRECTIONS.md Sec. 9, open question O2)."""
    # flagship: worked character 1, slope nu = (3/2, 1/2) on F_2 -> (1, 1/2) on F_0
    S2 = _Fe(2, (3, 1))
    v = SurfaceBundle(2, (3, 1), F(-1))
    vr, S0 = reduce(v, S2)
    e_src = dlp_envelope(total_slope(v), S2)
    e_red = dlp_envelope(total_slope(vr), S0)
    assert e_src.value == e_red.value == F(1)
    assert e_src.sharp is False and e_red.sharp is False
    assert not is_del_pezzo_anticanonical(S2) and not is_del_pezzo_anticanonical(S0)
    # emptiness_bound (the certified-empty theorem) is transported too
    assert emptiness_bound(v, S2) == emptiness_bound(vr, S0) == F(1, 2)

    # a genuine cusp slope nu = (2/3, 1/3) on F_2 -> (1/3, 1/3) on F_0, value 10/9 both sides
    nu_src = (F(2, 3), F(1, 3))
    env_src = dlp_envelope(nu_src, S2)
    env_red = dlp_envelope(pi_c1(nu_src), S0)
    assert env_src.value == env_red.value == F(10, 9)
    assert env_src.sharp is False and env_red.sharp is False


def test_cross_check_no_contradiction_over_sweep():
    """For every ample F_e (e = 2..6), reduce_to_del_pezzo's envelope equals the direct F_e
    envelope value (so ``lower_bound <= reduced`` holds with equality; never a contradiction)."""
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        for c1 in [(2, 1), (3, 1), (1, 2)]:
            xi = SurfaceBundle(2, c1, F(-1))
            xdp, Sdp = reduce_to_del_pezzo(xi, Se)
            direct = dlp_envelope(total_slope(xi), Se).value
            reduced = dlp_envelope(total_slope(xdp), Sdp).value
            assert direct == reduced
            assert direct <= reduced        # the acceptance inequality, here an exact equality


# --------------------------------------------------------------------------- #
# Differential regression against the independent oracle reference.            #
# --------------------------------------------------------------------------- #
def test_differential_against_oracle_reference():
    """The package pi agrees with the independent tests/oracle reference_reduce_pi (append-only,
    no frozen-row mutation): same r, c1, ch2, and reduced index over a sweep and both parities."""
    for e in _SWEEP_ES:
        Se = _ample_Fe(e)
        for r in (1, 2, 3):
            for x in range(-3, 4):
                for y in range(-3, 4):
                    xi = SurfaceBundle(r, (x, y), F(-1, 2))
                    xir, Sem2 = reduce(xi, Se)
                    ref_r, ref_c1, ref_ch2, ref_e = reference_reduce_pi(r, (x, y), F(-1, 2), e)
                    assert xir.r == ref_r
                    assert tuple(int(c) for c in xir.c1) == ref_c1
                    assert xir.ch2 == ref_ch2
                    assert hirzebruch_index(Sem2) == ref_e


def test_oracle_reference_rejects_del_pezzo():
    """reference_reduce_pi refuses e < 2, matching reduce()'s guard."""
    with pytest.raises(ValueError):
        reference_reduce_pi(2, (1, 1), F(0), 1)
