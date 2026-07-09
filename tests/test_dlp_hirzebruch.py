"""E11-M6 / G18: the polarization-dependent Drezet-Le Potier envelope on F_e.

Ground truth is Coskun-Huizenga, arXiv:1907.06739 ("Existence of semistable sheaves on
Hirzebruch surfaces").  Every pinned value below is either

  (a) transcribed from a numbered statement / table of that paper, or
  (b) re-derived here by exact ``Fraction`` arithmetic from a formula of that paper,

per the two-way rule of ``docs/CORRECTIONS.md``.  The strongest anchors are the paper's
**Tables 1 and 2** (stability intervals of every exceptional bundle of rank <= 19 on F_0
and <= 20 on F_1): reproducing them bit-for-bit exercises the potentially-exceptional
congruence, the DLP_{H,V} branches, the strip, the bounded search region, and the
rank-induction of Cor. "DLPExceptional" all at once.
"""
from fractions import Fraction as F

import pytest

from bridgeland_stability.dlp_hirzebruch import (
    DLPEnvelope,
    discriminant,
    dlp_bundle,
    dlp_envelope,
    dlp_restricted,
    emptiness_bound,
    exceptional_ch2,
    exceptional_discriminant,
    hilbert_P,
    hirzebruch_index,
    is_ample,
    is_anticanonical,
    is_del_pezzo_anticanonical,
    is_potentially_exceptional,
    is_semiexceptional,
    is_stable_exceptional,
    total_slope,
)
from bridgeland_stability.exceptional import P as P2_hilbert
from bridgeland_stability.exceptional_surface import SurfaceBundle, canonical_class
from bridgeland_stability.exceptional_surface import chi as surface_chi
from bridgeland_stability.nonemptiness_rational import (
    discriminant_H,
    hirzebruch_with_polarization,
)
from bridgeland_stability.varieties import P2, P1xP1, hirzebruch

# The anticanonical polarizations.  -K_{F_e} = (e+2) f + 2 s in the package's (f, s)
# basis; F_0 already ships H = f + s = (1,1), the same ray as -K_{F_0} = (2,2).
F1_K = hirzebruch_with_polarization(1, (3, 2))     # -K_{F_1}, strictly ample, d = 8


# --------------------------------------------------------------------------
# Paper Tables 1 and 2 -- the load-bearing regression.
#
# Coordinates.  The paper writes c1 = a E + b F with E the (-e)-section and F the
# fiber; this package stores NS-vectors in the basis (f, s) = (F, E).  So a paper
# entry (r, (a, b)) is the package character (r, c1 = (b, a)).  On F_0 the paper uses
# the two fiber classes (F_1, F_2) with c1 = a F_1 + b F_2, and the Gram [[0,1],[1,0]]
# is symmetric in them, so (a, b) may be read directly.
# --------------------------------------------------------------------------
TABLE_F0 = [  # (r, c1(paper)) for exceptional bundles on F_0, 0 <= a < r/2, a <= b < r
    (1, (0, 0)), (3, (1, 1)), (5, (1, 2)), (7, (1, 3)), (9, (1, 4)),
    (11, (1, 5)), (11, (4, 4)), (13, (1, 6)), (15, (1, 7)),
    (17, (1, 8)), (17, (5, 5)), (19, (1, 9)), (19, (4, 7)),
]
TABLE_F1 = [  # (r, c1 = a E + b F) on F_1, 0 <= a <= r/2, 0 <= b < r
    (1, (0, 0)), (2, (1, 1)), (4, (1, 2)), (5, (2, 2)), (6, (1, 3)),
    (8, (1, 4)), (10, (1, 5)), (11, (3, 5)), (12, (1, 6)), (13, (5, 5)),
    (14, (1, 7)), (16, (1, 8)), (18, (1, 9)), (19, (5, 10)), (20, (1, 10)),
]


def test_table_1_exceptional_bundles_on_F0():
    got = [
        (r, (a, b))
        for r in range(1, 20)
        for a in range(0, (r + 1) // 2) if 2 * a < r
        for b in range(a, r)
        if is_stable_exceptional(r, (b, a), P1xP1)          # package order (f, s)
    ]
    assert got == TABLE_F0


def test_table_2_exceptional_bundles_on_F1():
    got = [
        (r, (a, b))
        for r in range(1, 21)
        for a in range(0, r // 2 + 1)
        for b in range(0, r)
        if is_stable_exceptional(r, (b, a), F1_K)
    ]
    assert got == TABLE_F1


def test_ex_Kgraphs_contributing_ranks():
    # arXiv:1907.06739 Example "Kgraphs": DLP^{<8}_{-K_{F_0}} receives contributions from
    # exceptional bundles of ranks 1,3,5,7; DLP^{<7}_{-K_{F_1}} from ranks 1,2,4,5,6.
    assert sorted({r for r, _ in TABLE_F0 if r < 8}) == [1, 3, 5, 7]
    assert sorted({r for r, _ in TABLE_F1 if r < 7}) == [1, 2, 4, 5, 6]


def test_potentially_exceptional_is_not_exceptional():
    # (9,(2,2)) on F_0 satisfies the congruence of Lemma "excFacts" (2) -- it IS
    # potentially exceptional -- yet it is absent from Table 1: the rank induction of
    # Cor. "DLPExceptional" rejects it (Delta = 40/81 < DLP^{<9}_{-K}(2/9,2/9)).
    # Forced ch2 = <c1,c1>/(2r) - r Delta_V = 8/18 - 9*(40/81) = -4; then c2 = 4+4 = 8.
    assert is_potentially_exceptional(9, (2, 2), P1xP1) is True
    assert is_stable_exceptional(9, (2, 2), P1xP1) is False
    assert discriminant(SurfaceBundle(9, (2, 2), F(-4)), P1xP1) == exceptional_discriminant(9)
    assert exceptional_discriminant(9) == F(40, 81)


def test_no_even_rank_exceptionals_when_e_even():
    # Lemma "excFacts" (2): r is odd when e is even.  Exhaustive for small (r, c1) on F_0.
    for r in range(2, 13, 2):
        for a in range(-4, 5):
            for b in range(-4, 5):
                assert is_potentially_exceptional(r, (a, b), P1xP1) is False


# --------------------------------------------------------------------------
# The full-NS discriminant (G18a) and its relation to the H-projected surrogate.
# --------------------------------------------------------------------------
def test_discriminant_reduces_to_discriminant_H_on_p2():
    # On P^2, d = H^2 = 1 and Pic = Z.H, so the full-NS Delta and the H-projected
    # discriminant_H agree bit-for-bit.  All the pinned P^2 values are therefore stable.
    for xi in [SurfaceBundle(2, (1,), 0), SurfaceBundle(2, (1,), F(-5, 2)),
               SurfaceBundle(3, (2,), 1), SurfaceBundle(5, (2,), F(-2))]:
        assert discriminant(xi, P2) == discriminant_H(xi, P2)


def test_discriminant_equals_d_times_discriminant_H_iff_c1_parallel_H():
    # Delta = d * Delta_H exactly when c1 || H (in particular on every Picard-rank-1
    # surface).  The M4 F_0 entries are of this type; a non-diagonal c1 is not.
    diag = SurfaceBundle(2, (2, 2), F(-2))               # c1 = 2H, H = (1,1)
    assert discriminant(diag, P1xP1) == P1xP1.d * discriminant_H(diag, P1xP1) == 2

    off = SurfaceBundle(2, (1, 0), 0)                    # c1 = f, NOT parallel to H
    assert discriminant(off, P1xP1) == 0                 # nu^2 = 0, ch2 = 0
    assert discriminant_H(off, P1xP1) == F(1, 32)        # the lossy H-projected surrogate
    assert discriminant(off, P1xP1) != P1xP1.d * discriminant_H(off, P1xP1)


def test_discriminant_is_polarization_independent():
    # Delta = 1/2 <nu,nu> - ch2/r never mentions H.  Two different ample polarizations
    # on F_1 give the same Delta -- unlike discriminant_H, which is built from mu_H.
    A = hirzebruch_with_polarization(1, (2, 1))
    B = hirzebruch_with_polarization(1, (4, 1))
    xi = SurfaceBundle(2, (1, 1), F(1, 2))
    assert discriminant(xi, A) == discriminant(xi, B) == F(-1, 8)
    assert discriminant_H(xi, A) == F(-1, 36) != F(1, 196) == discriminant_H(xi, B)


def test_hilbert_P_reduces_to_p2_polynomial():
    # P(nu) = chi(O) + 1/2 (nu^2 - nu.K) must reduce on P^2 to (m^2 + 3m + 2)/2.
    for m in [F(0), F(1), F(-3), F(1, 2), F(-5, 3)]:
        assert hilbert_P((m,), P2) == P2_hilbert(m)


def test_hilbert_P_on_Fe_matches_paper_closed_form():
    # arXiv:1907.06739 Sec. 2.1: P(a E + b F) = (a+1)(b+1 - a e / 2) on F_e.
    for surf, e in [(P1xP1, 0), (F1_K, 1)]:
        for a in [F(0), F(1), F(-2), F(3, 2)]:
            for b in [F(0), F(2), F(-1), F(1, 3)]:
                nu = (b, a)                              # package (f, s) = (F, E)
                assert hilbert_P(nu, surf) == (a + 1) * (b + 1 - a * e / 2)


# --------------------------------------------------------------------------
# DLP_{H,V}, the strip, and the envelope.
# --------------------------------------------------------------------------
def test_dlp_bundle_outside_strip_is_none():
    # The bound is only proved on |(nu - nu_V).H| <= -1/2 K.H; outside it, None.
    assert dlp_bundle((F(50), F(50)), 1, (0, 0), P1xP1) is None
    assert dlp_bundle((F(0), F(0)), 1, (0, 0), P1xP1) == 1        # w = 0: P(0) - 0 = chi(O)


def test_cor_38_line_bundle_envelope_floor():
    # arXiv:1907.06739 Cor. "38": for e = 0 or 1 and ANY polarization, DLP^1_H(nu) >= 3/8.
    # DLP^{<2} is exactly the line-bundle envelope DLP^1.
    surfaces = [P1xP1, hirzebruch_with_polarization(0, (3, 1)), F1_K,
                hirzebruch_with_polarization(1, (3, 1))]
    for surf in surfaces:
        for i in range(0, 5):
            for j in range(0, 5):
                nu = (F(i, 4), F(j, 4))
                assert dlp_restricted(nu, surf, 2) >= F(3, 8)


def test_delta_K_at_line_bundle_slopes_is_one():
    # At nu = nu(O(D)) the controlling bundle is the line bundle itself: Delta_V = 0 and
    # DLP_{-K,O(D)}(nu) = P(0) - 0 = chi(O_{F_e}) = 1.  These are exactly the two F_0
    # entries the E11-M4 paper table could only TRANSCRIBE (delta_H_paper = 1); they are
    # now computed natively, and the truncation certifies itself as exact.
    for nu in [(F(0), F(0)), (F(1), F(1)), (F(2), F(-1))]:
        env = dlp_envelope(nu, P1xP1, 8)
        assert env.value == 1 and env.exact and env.sharp
        assert env.witness == (1, (int(nu[0]), int(nu[1])))
    env = dlp_envelope((F(0), F(0)), F1_K, 8)
    assert env.value == 1 and env.exact and env.sharp


def test_delta_K_at_exceptional_slopes_is_the_cusp():
    # At the slope of a mu_{-K}-stable exceptional bundle V, DLP_{-K}(nu_V) = 1 - Delta_V
    # = 1/2 + 1/(2 r_V^2) -- the cusp.  (Derived in the module docstring; the same shape
    # as the P^2 curve's cusps.)
    env = dlp_envelope((F(1, 3), F(1, 3)), P1xP1, 8)     # rank-3 exceptional (3,(1,1))
    assert env.value == F(5, 9) == 1 - exceptional_discriminant(3)
    assert env.exact and env.sharp and env.witness == (3, (1, 1))

    env = dlp_envelope((F(1, 2), F(1, 2)), F1_K, 8)      # rank-2 exceptional (2,(1,1))
    assert env.value == F(5, 8) == 1 - exceptional_discriminant(2)
    assert env.exact and env.sharp and env.witness == (2, (1, 1))


def test_envelope_floor_is_one_half_on_anticanonical_del_pezzo():
    # arXiv:1907.06739 Cor. "K1/2": DLP_{-K_{F_e}}(nu) >= 1/2 for e = 0, 1.  (The F_e
    # analogue of the 1/2 clamp of the P^2 Drezet-Le Potier curve.)
    for surf in (P1xP1, F1_K):
        for i in range(0, 6):
            for j in range(0, 6):
                assert dlp_envelope((F(i, 5), F(j, 5)), surf, 9).value >= F(1, 2)


def test_envelope_is_a_lower_bound_off_the_anticanonical_ray():
    # For a non-anticanonical H the envelope is only a CERTIFIED LOWER BOUND on
    # delta_H^{mu-s} (Cor. "deltaDLPe"): never claim exact/sharp.
    S = hirzebruch_with_polarization(1, (3, 1))          # m = 2, not anticanonical
    assert not is_anticanonical(S) and not is_del_pezzo_anticanonical(S)
    env = dlp_envelope((F(1, 2), F(1, 2)), S, 8)
    assert env.sharp is False and env.exact is False
    assert env.value > 0


# --------------------------------------------------------------------------
# The counterexample the paper states just before Cor. "K1/2" -- why the third
# branch of DLP_{H,V} must NOT be used as an emptiness certificate.
# --------------------------------------------------------------------------
def test_third_branch_is_not_an_emptiness_certificate():
    # On F_0 with the (non-generic) anticanonical H, the bundle O(F_1) + O(F_2) is
    # Gieseker-semistable with Delta = 1/4, yet DLP_{-K}(1/2,1/2) = 3/4 -- computed by
    # the excluded branch (nu - nu_V).H = 0 at V = O(F_1).  emptiness_bound must drop
    # that branch and conclude nothing, or the package would call an existing sheaf empty.
    xi = SurfaceBundle(2, (1, 1), 0)                     # ch(O(f) + O(s)) = (2, f+s, 0)
    assert discriminant(xi, P1xP1) == F(1, 4)
    env = dlp_envelope((F(1, 2), F(1, 2)), P1xP1, 8)
    assert env.value == F(3, 4) and env.witness == (1, (1, 0))   # via the excluded branch
    assert emptiness_bound(xi, P1xP1, 8) == F(1, 4)              # strict branches only
    assert not discriminant(xi, P1xP1) < emptiness_bound(xi, P1xP1, 8)   # no emptiness claim


def test_emptiness_bound_at_a_line_bundle_slope():
    # nu = nu(O) with Delta != Delta_O: Gieseker semistability forces Delta >= 1 - 0 = 1.
    # This is the branch that reproduces the M4 F_0 targets (delta_H_paper = 1).
    xi = SurfaceBundle(2, (0, 0), F(-4))                 # Delta = 2
    assert discriminant(xi, P1xP1) == 2
    assert emptiness_bound(xi, P1xP1, 8) == 1


# --------------------------------------------------------------------------
# Polarization dependence: the stability interval I_V of Table 2 moves a bundle
# in and out of existence.  This is the E11-M5 witness done for real.
# --------------------------------------------------------------------------
def test_stability_interval_of_the_rank_2_bundle_on_F1():
    # Table 2: V = (2, E + F) on F_1 has I_V = (0, 1).  H_m = E + (e+m) F is, in the
    # package (f, s) basis, (e + m, 1) up to scale.  m = 1/2 is anticanonical.
    inside = [hirzebruch_with_polarization(1, (3, 2)),   # m = 3/2 - 1 = 1/2
              hirzebruch_with_polarization(1, (19, 10))]  # m = 19/10 - 1 = 9/10 < 1
    outside = [hirzebruch_with_polarization(1, (2, 1)),  # m = 2 - 1 = 1   (endpoint, open)
               hirzebruch_with_polarization(1, (3, 1))]  # m = 3 - 1 = 2
    for S in inside:
        assert is_stable_exceptional(2, (1, 1), S) is True
    for S in outside:
        assert is_stable_exceptional(2, (1, 1), S) is False


def test_stability_intervals_from_table_1():
    # F_0 Table 1: V = (3,(1,1)) has I_V = (1/2, 2); V = (17,(5,5)) has I_V = (8/9, 9/8).
    # H_m on F_0 is (m, 1) up to scale in the (f, s) basis.
    def Fm(num, den):                                    # H_{num/den} = (num, den)
        return hirzebruch_with_polarization(0, (num, den))

    for num, den, want in [(2, 5, False), (1, 2, False), (3, 5, True), (1, 1, True),
                           (19, 10, True), (2, 1, False), (21, 10, False)]:
        assert is_stable_exceptional(3, (1, 1), Fm(num, den)) is want

    for num, den, want in [(4, 5, False), (8, 9, False), (1, 1, True),
                           (9, 8, False), (6, 5, False)]:
        assert is_stable_exceptional(17, (5, 5), Fm(num, den)) is want


def test_polarization_dependence_of_nonemptiness():
    # THE witness.  xi = (2, f + s, -1/2) on F_1 is the rank-2 exceptional character
    # (Delta = 3/8 = 1/2 - 1/(2*4), c2 = 1).  Delta does NOT depend on H.  Yet:
    #   * H = (3,2) = -K (m = 1/2, inside I_V = (0,1)):  V is mu_H-stable exceptional,
    #     so M_H(xi) is a single reduced point -- NON-EMPTY, sitting strictly below its
    #     own envelope delta_H = 5/8 (exactly like T_{P^2}(-1) below the DLP curve).
    #   * H = (3,1) (m = 2, outside I_V):  no mu_H-stable exceptional bundle of this
    #     character; the line bundle O(f) certifies Delta = 3/8 < 7/8 = emptiness_bound,
    #     so M_H(xi) is PROVABLY EMPTY.
    # Hand check of the emptiness bound: w = nu - nu_{O(f)} = (-1/2, 1/2); with H = (3,1)
    # and Gram [[0,1],[1,-1]], w.H = 1/2 > 0 (strict branch, inside the strip since
    # -1/2 K.H = 7/2), so the bound is P(-w) = 1 + 1/2((-w)^2 - (-w).K)
    #                                        = 1 + 1/2(-3/4 + 1/2) = 7/8.
    xi = SurfaceBundle(2, (1, 1), F(-1, 2))
    K = hirzebruch_with_polarization(1, (3, 2))
    A = hirzebruch_with_polarization(1, (3, 1))

    assert discriminant(xi, K) == discriminant(xi, A) == F(3, 8) == exceptional_discriminant(2)
    assert is_potentially_exceptional(2, (1, 1), K)

    assert is_stable_exceptional(2, (1, 1), K) is True
    assert dlp_envelope(total_slope(xi), K, 8).value == F(5, 8)
    assert emptiness_bound(xi, K, 8) == F(3, 8)                  # == Delta: no emptiness claim

    assert is_stable_exceptional(2, (1, 1), A) is False
    assert emptiness_bound(xi, A, 8) == F(7, 8)                  # hand-derived above
    assert discriminant(xi, A) < emptiness_bound(xi, A, 8)       # PROVEN empty


def test_semiexceptional_detection():
    # V^{+m} is Gieseker-semistable, so it is non-empty even though Delta = Delta_V lies
    # strictly below the envelope.  On F_1, V = (2, f+s, -1/2); V^{+2} = (4, 2f+2s, -1).
    V = SurfaceBundle(2, (1, 1), F(-1, 2))
    V2 = SurfaceBundle(4, (2, 2), F(-1))
    assert discriminant(V2, F1_K) == discriminant(V, F1_K) == F(3, 8)
    assert is_semiexceptional(V, F1_K) is True
    assert is_semiexceptional(V2, F1_K) is True
    assert is_semiexceptional(SurfaceBundle(2, (1, 1), 0), F1_K) is False   # Delta = 1/8


# --------------------------------------------------------------------------
# Guards: the machinery refuses what it cannot honestly model.
# --------------------------------------------------------------------------
def test_nef_but_not_ample_factory_polarization_is_refused():
    # varieties.hirzebruch(n) carries H = n f + s, which is nef-and-big but NOT ample
    # (a - n b = 0).  The CH theory needs an ample H, so the envelope refuses it.
    for n in (1, 2, 3):
        S = hirzebruch(n)
        assert is_ample(S) is False
        with pytest.raises(ValueError):
            dlp_envelope((F(0), F(0)), S, 4)


def test_non_hirzebruch_surface_is_refused():
    with pytest.raises(NotImplementedError):
        hirzebruch_index(P2)
    with pytest.raises(NotImplementedError):
        dlp_envelope((F(0),), P2, 4)


def test_index_and_polarization_predicates():
    assert hirzebruch_index(P1xP1) == 0
    assert hirzebruch_index(F1_K) == 1
    assert is_anticanonical(P1xP1) and is_del_pezzo_anticanonical(P1xP1)
    assert is_anticanonical(F1_K) and is_del_pezzo_anticanonical(F1_K)
    S = hirzebruch_with_polarization(1, (3, 1))
    assert is_ample(S) and not is_anticanonical(S)
    # F_e is a del Pezzo surface only for e = 0, 1.  For e >= 2 the anticanonical class
    # -K = (e+2) f + 2 s fails Nakai (a - e b = 2 - e <= 0), so it is not even ample and
    # cannot be handed to the envelope at all -- there is no anticanonical F_2 to test.
    with pytest.raises(ValueError):
        hirzebruch_with_polarization(2, (4, 2))
    S2 = hirzebruch_with_polarization(2, (5, 1))          # ample, but not anticanonical
    assert is_ample(S2) and not is_anticanonical(S2) and not is_del_pezzo_anticanonical(S2)


# --------------------------------------------------------------------------
# Two-way verification against an INDEPENDENT Euler form.
#
# exceptional_surface.chi implements Riemann-Roch directly from the Gram, K_X and
# chi(O_X) -- it shares no code with hilbert_P / discriminant.  Pinning agreement
# makes a sign or normalization drift in EITHER implementation fail loudly, exactly as
# test_nonemptiness.py::test_m1_euler_agrees_with_m3_exceptional_disjunct does on P^2.
# --------------------------------------------------------------------------
def _bundle(r, c1_paper, surface):
    c1 = (c1_paper[1], c1_paper[0])                  # paper (a=E, b=F) -> package (f, s)
    return SurfaceBundle(r, c1, exceptional_ch2(r, c1, surface))


def test_every_tabulated_bundle_has_self_euler_one():
    # chi(V,V) = 1 is the DEFINITION of potentially exceptional; verify it for all 28
    # bundles of Tables 1-2 through the independent RR of exceptional_surface.chi.
    for surf, table in ((P1xP1, TABLE_F0), (F1_K, TABLE_F1)):
        for r, c1p in table:
            V = _bundle(r, c1p, surf)
            assert surface_chi(V, V, surf) == 1
            assert discriminant(V, surf) == exceptional_discriminant(r)
            c2 = F(1, 2) * surf.lattice.self_pairing(V.c1) - V.ch2
            assert c2.denominator == 1               # integral character


def test_riemann_roch_pairing_matches_independent_euler_form():
    # arXiv:1907.06739 Sec. 2.1: chi(V,W) = r_V r_W (P(nu_W - nu_V) - Delta_V - Delta_W).
    # LHS from exceptional_surface.chi; RHS from this module's hilbert_P + discriminant.
    for surf, table in ((P1xP1, TABLE_F0[:6]), (F1_K, TABLE_F1[:6])):
        bundles = [_bundle(r, c1p, surf) for r, c1p in table]
        for V in bundles:
            for W in bundles:
                w = tuple(x - y for x, y in zip(total_slope(W), total_slope(V)))
                rhs = V.r * W.r * (hilbert_P(w, surf)
                                   - discriminant(V, surf) - discriminant(W, surf))
                assert surface_chi(V, W, surf) == rhs


def test_canonical_class_of_Fe():
    # K_{F_e} = -(e+2) f - 2 s; K^2 = 8 for every e (F_e is a rational surface).
    for surf, e in ((P1xP1, 0), (F1_K, 1)):
        K = canonical_class(surf)
        assert K == (F(-(e + 2)), F(-2))
        assert surf.lattice.self_pairing(K) == 8


def test_everything_is_exact_fractions():
    env = dlp_envelope((F(1, 3), F(1, 3)), P1xP1, 6)
    assert isinstance(env.value, F) and isinstance(env, DLPEnvelope)
    assert isinstance(discriminant(SurfaceBundle(2, (1, 1), F(1, 3)), P1xP1), F)
    assert isinstance(hilbert_P((F(1, 3), F(1, 5)), P1xP1), F)
    assert isinstance(emptiness_bound(SurfaceBundle(2, (0, 0), F(-4)), P1xP1, 6), F)
