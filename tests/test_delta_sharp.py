"""E14-M1: the sharp Bogomolov function ``delta_m^{mu-s}(nu)`` as a certified sandwich.

Pins the two Sec. 8 Kronecker values of arXiv:1907.06739 as SANDWICH facts (the inf
is not attained at finite rank -- see ``docs/CORRECTIONS.md`` Sec. 17), the paper's
two ``DLP^{<r}`` computer values, the anticanonical pinches, and the certified
decision procedure ``mu_stable_exists`` on both sides of the mu-stable/Gieseker gap.

Every asserted Fraction was hand-derived from the paper (Ex. "KroneckerF0" /
"KroneckerF1", Table "stabilityInterval1") and re-computed exactly before pinning.
"""

from fractions import Fraction as F

import pytest

from bridgeland_stability.delta_sharp import (
    DeltaSharp,
    delta_mu_stable,
    mu_stable_exists,
    polarization_index,
    surface_with_index,
)
from bridgeland_stability.dlp_hirzebruch import (
    discriminant,
    dlp_envelope,
    dlp_restricted,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi
from bridgeland_stability.hn_filtration import hn_verdict
from bridgeland_stability.reduction import reduce_to_del_pezzo
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)

# Paper slopes in package (f, s) coordinates: paper eps*E + phi*F -> (phi, eps).
NU_F0 = (F(1, 3), F(1, 5))       # Ex. KroneckerF0: nu = 1/5 E + 1/3 F, m = 25/9
NU_F1 = (F(6, 13), F(3, 13))     # Ex. KroneckerF1: nu = 3/13 E + 6/13 F, m = 12/7
M_F0 = F(25, 9)
M_F1 = F(12, 7)
DELTA_F0 = F(3, 5)               # delta_{25/9}^{mu-s}(NU_F0), thm-deltaKronecker
DELTA_F1 = F(98, 169)            # delta_{12/7}^{mu-s}(NU_F1)


# ---------------------------------------------------------------------------
# The Sec. 8 Kronecker pins: wall classes refuse, the sandwich brackets the value
# ---------------------------------------------------------------------------

def test_f0_kronecker_wall_class_has_no_mu_stable_sheaf():
    # At m = 25/9 the general sheaf of character (15, NU_F0, 3/5) is STRICTLY
    # mu-semistable (thm-intervalKronecker); mu-stable sheaves of that exact
    # character do not exist even though Delta equals delta^{mu-s} -- the inf
    # is approached, not attained.  One lattice step up, they exist
    # (thm-deltaSurface (1)).
    S = surface_with_index(0, M_F0)
    assert mu_stable_exists(15, NU_F0, DELTA_F0, S) is False
    assert mu_stable_exists(15, NU_F0, F(2, 3), S) is True


def test_f0_kronecker_sandwich_brackets_the_paper_value():
    r = delta_mu_stable(NU_F0, M_F0, F0, max_rank=15)
    assert isinstance(r, DeltaSharp)
    # lower: the DLP envelope truncation; equals the paper's computer value
    # DLP^{<15}_{H_{25/9}}(NU_F0) = 19/35 (Ex. KroneckerF0) -- see the direct
    # dlp_restricted pin below.
    assert r.lower == F(19, 35)
    # upper: the least rank-15 lattice discriminant > 3/5 is 2/3, and it
    # certifies; 3/5 itself refuses (previous test), so upper = 2/3 exactly.
    assert r.upper == F(2, 3)
    assert r.lower <= DELTA_F0 <= r.upper       # the sandwich brackets the truth
    assert r.upper > DELTA_F0                    # ... and cannot attain it (Sec. 17)
    assert r.exact is False
    assert r.witness[0] == 15
    assert r.certificate.rigor == Rigor.PROVEN


def test_f1_kronecker_wall_class_and_sandwich():
    S = surface_with_index(1, M_F1)
    # The M3b paper-pin class (13, (6,3), -13/2) sits exactly at Delta = 98/169.
    assert mu_stable_exists(13, NU_F1, DELTA_F1, S) is False
    r = delta_mu_stable(NU_F1, M_F1, F1, max_rank=13)
    # lower = the paper's computer value DLP^{<13}_{H_{12/7}}(NU_F1) = 523/1014.
    assert r.lower == F(523, 1014)
    # first rank-13 lattice step above 98/169: ch2 = t/2 (t odd), Delta =
    # (27 - 13 t)/338; t = -13 gives 98/169 (refused), t = -15 gives 111/169.
    assert r.upper == F(111, 169)
    assert r.lower <= DELTA_F1 <= r.upper
    assert r.exact is False


def test_paper_dlp_restricted_computer_values_pin_bit_for_bit():
    # arXiv:1907.06739 Ex. KroneckerF0 / KroneckerF1: "computer calculations show"
    # DLP^{<15}_{H_{25/9}}(1/5 E + 1/3 F) = 19/35 and
    # DLP^{<13}_{H_{12/7}}(3/13 E + 6/13 F) = 523/1014.  First literature
    # cross-check of dlp_restricted off the anticanonical ray.
    S0 = surface_with_index(0, M_F0)
    S1 = surface_with_index(1, M_F1)
    assert dlp_restricted(NU_F0, S0, 15) == F(19, 35)
    assert dlp_restricted(NU_F1, S1, 13) == F(523, 1014)


# ---------------------------------------------------------------------------
# Anticanonical pinches: the sandwich closes where the sharp theory is shipped
# ---------------------------------------------------------------------------

def test_anticanonical_pinch_nu_zero():
    for surf, e in ((F0, 0), (F1, 1)):
        m0 = F(2 - e, 2)
        r = delta_mu_stable((0, 0), m0, surf, max_rank=4)
        assert r.lower == 1 and r.upper == 1 and r.exact is True
        # the witness is the rank-one ideal sheaf of a point: (1, (0,0), c2 = 1)
        assert r.witness[0] == 1
        # and the lower bound came from the sharp envelope
        env = dlp_envelope((0, 0), surface_with_index(e, m0))
        assert env.value == 1 and env.sharp is True


def test_anticanonical_pinch_half_half_on_f0():
    r = delta_mu_stable((F(1, 2), F(1, 2)), F(1), F0, max_rank=4)
    assert r.lower == F(3, 4) and r.upper == F(3, 4) and r.exact is True
    assert r.witness == (2, (1, 1), F(-1))


def test_flagship_cousin_separates_mu_stable_from_semistable():
    # (2, (1,1), 0) on anticanonical F_0 -- the E13 flagship: Gieseker-semistable
    # sheaves EXIST (O(F1) + O(F2) is polystable; hn_verdict says so), but a
    # mu-stable sheaf would be simple with chi(v,v) = 4(1 - 2*(1/4)) = 2 >= 2 --
    # impossible.  The two notions genuinely separate on this class.
    S = surface_with_index(0, F(1))
    nu = (F(1, 2), F(1, 2))
    assert hn_verdict(2, nu, F(1, 4), S).exists is True
    assert mu_stable_exists(2, nu, F(1, 4), S) is False


# ---------------------------------------------------------------------------
# The decision procedure: rank one, invalid characters, the exceptional branch,
# the honest Delta = 1/2 edge, e >= 2 transport
# ---------------------------------------------------------------------------

def test_rank_one_is_the_ideal_sheaf_test():
    S = surface_with_index(0, F(3, 2))
    assert mu_stable_exists(1, (0, 0), F(2), S) is True     # I_Z, len 2
    assert mu_stable_exists(1, (0, 0), F(0), S) is True     # O
    assert mu_stable_exists(1, (0, 0), F(-1), S) is False   # c2 < 0
    assert mu_stable_exists(1, (0, 0), F(1, 2), S) is False  # c2 non-integral
    assert mu_stable_exists(1, (F(1, 2), 0), F(1), S) is False  # c1 non-integral


def test_invalid_characters_and_bogomolov_refuse():
    S = surface_with_index(0, F(1))
    assert mu_stable_exists(3, (F(1, 2), 0), F(7, 10), S) is False   # c1 not integral
    assert mu_stable_exists(2, (0, 0), F(-3), S) is False            # Delta < 0
    with pytest.raises(ValueError):
        mu_stable_exists(F(3, 2), (0, 0), F(1), S)                   # never truncate rank


def test_delta_half_edge_is_honestly_undecided():
    S = surface_with_index(0, F(1))
    assert mu_stable_exists(2, (0, 0), F(1, 2), S) is None


def test_exceptional_branch_follows_the_stability_interval():
    # (2, (1,1)) on F_1 is exceptional with Delta = 3/8 and stability interval
    # I_V = (0, 1) (Table "stabilityInterval1").  At the anticanonical index
    # m = 1/2 in I_V it is mu-stable; at m = 3/2 outside the closure it is not.
    nu = (F(1, 2), F(1, 2))
    assert mu_stable_exists(2, nu, F(3, 8), surface_with_index(1, F(1, 2))) is True
    assert mu_stable_exists(2, nu, F(3, 8), surface_with_index(1, F(3, 2))) is False


def test_chi_identity_underpinning_the_small_delta_branch():
    # chi(v,v) = r^2 (1 - 2 Delta) on F_e (Lemma "excFacts" (1)) -- the identity
    # the Delta < 1/2 branch trips on loudly if violated.
    for surf in (F0, F1):
        for (r, c1, ch2) in ((2, (1, 1), F(0)), (3, (1, 2), F(-1)),
                             (5, (2, 3), F(1)), (4, (0, 1), F(-2))):
            x = SurfaceBundle(r, c1, ch2)
            D = discriminant(x, surf)
            assert chi(x, x, surf) == r * r * (1 - 2 * D)


def test_e_ge_2_transport_matches_manual_reduction():
    # cor-highermus: the decision on F_2 equals the decision for the pi-reduced
    # character on F_0 with the transported polarization -- exercised through the
    # public API on both sides.
    S2 = surface_with_index(2, F(3, 4))          # H = (11, 4) on F_2, ample
    for (r, nu, D) in ((2, (F(1, 2), F(1, 2)), F(3, 4)),
                       (2, (F(1, 2), F(1, 2)), F(5, 4)),
                       (3, (F(1, 3), F(2, 3)), F(8, 9)),
                       (1, (0, 0), F(1))):
        c1 = tuple(int(r * x) for x in nu)
        ch2 = r * (F(1, 2) * S2.lattice.self_pairing(nu) - D)
        xi_red, S_red = reduce_to_del_pezzo(SurfaceBundle(r, c1, ch2), S2)
        nu_red = tuple(F(c) / r for c in xi_red.c1)
        assert mu_stable_exists(r, nu, D, S2) == mu_stable_exists(
            r, nu_red, D, S_red)


# ---------------------------------------------------------------------------
# Structural gates: sandwich ordering, scan monotonicity, mu-s => Gieseker-ss
# ---------------------------------------------------------------------------

def test_lower_never_exceeds_upper_on_a_sweep():
    # A certified lower bound above a certified upper bound would falsify one of
    # the two theorem routes -- hard failure by construction.
    for surf, e in ((F0, 0), (F1, 1)):
        for m in (F(1), F(3, 2), F(5, 7)):
            for nu in ((0, 0), (F(1, 2), F(1, 2)), (F(1, 3), F(1, 3))):
                r = delta_mu_stable(nu, m, surf, max_rank=6)
                if r.upper is not None:
                    assert r.lower <= r.upper
                    assert r.exact == (r.lower == r.upper)


def test_upper_is_nonincreasing_in_max_rank():
    uppers = [delta_mu_stable((F(1, 2), F(1, 2)), F(1), F0, max_rank=k).upper
              for k in (2, 4, 6)]
    assert all(u is not None for u in uppers)
    assert uppers[0] >= uppers[1] >= uppers[2]


def test_mu_stable_implies_gieseker_semistable_at_m():
    # mu-stable => Gieseker-stable => Gieseker-semistable AT the same m: whenever
    # the certifier says True, hn_verdict at m itself must agree sheaves exist.
    S = surface_with_index(0, F(3, 2))
    for (r, nu, D) in ((2, (F(1, 2), F(1, 2)), F(3, 4)),
                       (2, (0, 0), F(1)),
                       (3, (F(1, 3), F(2, 3)), F(7, 9))):
        if mu_stable_exists(r, nu, D, S) is True:
            assert hn_verdict(r, nu, D, S).exists is True


def test_polarization_index_round_trip():
    for e in (0, 1, 2, 3):
        for m in (F(1), F(1, 2), F(25, 9), F(12, 7)):
            assert polarization_index(surface_with_index(e, m)) == m
    with pytest.raises(ValueError):
        surface_with_index(0, F(0))
    with pytest.raises(ValueError):
        surface_with_index(1, F(-1, 2))


# ===========================================================================
# E14-M2: thm-deltaKronecker -- the closed formula on the Kronecker triangle
# ===========================================================================

from bridgeland_stability.delta_sharp import (   # noqa: E402
    delta_kronecker,
    kronecker_data,
    _psi_gt,
)

F2 = hirzebruch(2)
F3 = hirzebruch(3)


def test_kronecker_formula_reproduces_both_paper_values():
    # Ex. KroneckerF0: delta_{25/9}^{mu-s}(1/5 E + 1/3 F) = 3/5 on F_0;
    # Ex. KroneckerF1: delta_{12/7}^{mu-s}(3/13 E + 6/13 F) = 98/169 on F_1.
    # Hand re-derivation of the full RR simplification: CORRECTIONS Sec. 18.
    assert delta_kronecker(NU_F0, M_F0, F0, l=3) == DELTA_F0
    assert delta_kronecker(NU_F1, M_F1, F1, l=3) == DELTA_F1
    # the l=None window scan finds the same (single) window
    assert delta_kronecker(NU_F0, M_F0, F0) == DELTA_F0
    assert delta_kronecker(NU_F1, M_F1, F1) == DELTA_F1


def test_kronecker_internals_pin_ex_triangle_data():
    # ex-triangle / ex-KroneckerF1 (e = 1, l = 3, m = 12/7): the line of slope
    # -12/7 meets L_K at nu1 = (1/2, 0) and L_L at nu2 = (2/11, 6/11), with
    # exponent ratios b/a = 1 and d/c = 13/2.  PAPER ERRATUM (recorded in
    # CORRECTIONS Sec. 18): ex-triangle prints nu = (2/13, 6/13), but the point
    # of the stated chord is (3/13, 6/13) -- the ex-KroneckerF1 slope; (2/13,
    # 6/13) fails the collinearity check y = -12/7 (x - 1/2) exactly.
    d = kronecker_data(NU_F1, M_F1, F1, 3)
    assert d is not None
    assert (d.x0, d.y0) == (F(3, 13), F(6, 13))
    assert d.x1 == F(1, 2) and d.x2 == F(2, 11)
    assert d.b_over_a == 1 and d.d_over_c == F(13, 2)
    assert d.lam == F(2, 13)                      # nu = lam*nu1 + (1-lam)*nu2
    assert d.lam * d.x1 + (1 - d.lam) * d.x2 == d.x0
    assert d.value == DELTA_F1
    assert (d.k, d.N, d.M) == (2, 3, 7)


def test_psi_comparisons_are_exact():
    # psi_3 = (3 + sqrt 5)/2 ~ 2.618, psi_4 = 2 + sqrt 3 ~ 3.732: integer-square
    # comparisons, no float.
    assert _psi_gt(3, F(5, 2)) is True
    assert _psi_gt(3, F(8, 3)) is False
    assert _psi_gt(4, F(7, 2)) is True
    assert _psi_gt(4, F(15, 4)) is False


def test_kronecker_predicate_rejects_outside_the_window():
    # m out of the geometric range (1 - e/2, k)
    assert kronecker_data(NU_F0, F(3), F0, 3) is None          # m = k
    assert kronecker_data(NU_F0, F(7, 2), F0, 3) is None       # m > k
    assert kronecker_data(NU_F0, F(-1), F0, 3) is None         # m <= 0
    # the vertex P4 = (1/(2l-e), l/(2l-e)) itself: x2 hits the OPEN right
    # endpoint, refused (this is where the formula denominator would vanish)
    assert kronecker_data((F(1, 2), F(1, 6)), F(2), F0, 3) is None
    # a slope far outside the triangle
    assert kronecker_data((F(5), F(5)), F(2), F0, 3) is None
    with pytest.raises(ValueError):
        kronecker_data(NU_F0, M_F0, F0, 2)                     # l >= 3
    with pytest.raises(ValueError):
        kronecker_data(NU_F0, M_F0, F2, 3)                     # del Pezzo statement


def test_kronecker_pins_admit_exactly_one_window():
    # Empirically the (nu, m) windows of distinct l tile without overlap (a
    # ~50k-probe sweep found none -- CORRECTIONS Sec. 18); where they WOULD
    # overlap, delta_kronecker(l=None) asserts equal values (both equal
    # delta_m^{mu-s} by the theorem).  The paper points admit exactly l = 3.
    for surf, nu, m in ((F0, NU_F0, M_F0), (F1, NU_F1, M_F1)):
        admitting = [l for l in range(3, 16)
                     if kronecker_data(nu, m, surf, l) is not None]
        assert admitting == [3]


def test_kronecker_formula_vs_sandwich_differential():
    # Wherever the closed formula applies, the M1 sandwich must bracket it --
    # with STRICT upper (the inf is not attained at these wall polarizations:
    # the general sheaf at Delta = delta is strictly mu-semistable).  Two
    # independent theorem routes to the same number; a violation falsifies one.
    for surf, nu, m, mr in ((F0, NU_F0, M_F0, 15),
                            (F0, NU_F0, F(5, 2), 15),
                            (F1, NU_F1, M_F1, 13),
                            (F1, NU_F1, F(3, 2), 13)):
        dk = delta_kronecker(nu, m, surf)
        assert dk is not None
        r = delta_mu_stable(nu, m, surf, max_rank=mr)
        assert r.lower <= dk < r.upper


def test_kronecker_values_increase_with_m():
    # cor-deltaMonotone preview on the closed formula: for m >= m0 = 1 - e/2 the
    # value is nondecreasing in m (here strictly increasing on the sampled
    # chords through each triangle).
    f0_vals = [delta_kronecker(NU_F0, m, F0)
               for m in (F(5, 2), F(8, 3), F(11, 4), M_F0, F(14, 5))]
    f1_vals = [delta_kronecker(NU_F1, m, F1)
               for m in (F(3, 2), F(5, 3), M_F1, F(7, 4), F(9, 5))]
    for vals in (f0_vals, f1_vals):
        assert all(v is not None for v in vals)
        assert all(a < b for a, b in zip(vals, vals[1:]))
    assert f0_vals == [F(26, 45), F(62, 105), F(242, 405), F(3, 5), F(298, 495)]
    assert f1_vals == [F(379, 676), F(1653, 2873), DELTA_F1, F(2169, 3718),
                       F(895, 1521)]


def test_kronecker_transports_to_e_ge_2():
    # delta_{m, F_e} = delta_{m+1, F_{e-2}} o pi (cor-highermus + Lemma 11.3):
    # the F_0 pin lifts to F_2 at slope pi^{-1}(nu) = (8/15, 1/5), m = 25/9 - 1
    # = 16/9, and the F_1 pin to F_3 at (9/13, 3/13), m = 12/7 - 1 = 5/7 (the
    # index SHRINKS by one per lift: the reduction adds one).
    assert delta_kronecker((F(8, 15), F(1, 5)), F(16, 9), F2) == DELTA_F0
    assert delta_kronecker((F(9, 13), F(3, 13)), F(5, 7), F3) == DELTA_F1
    # and the M1 decision procedure agrees on F_2 directly: the transported wall
    # class refuses, one lattice step up certifies.
    S2 = surface_with_index(2, F(16, 9))
    assert mu_stable_exists(15, (F(8, 15), F(1, 5)), DELTA_F0, S2) is False
    assert mu_stable_exists(15, (F(8, 15), F(1, 5)), F(2, 3), S2) is True
