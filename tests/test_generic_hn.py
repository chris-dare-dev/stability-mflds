"""E13-M3b / G18: the generic Harder-Narasimhan filtration (arXiv:1907.06739 Sec. 5).

Ground truth is the paper's Sec. 5 algorithm (thm-HNcriterion / cor-algorithm) and
its two EXPLICIT worked multi-factor examples (Sec. 1.5, "Harder-Narasimhan
filtrations from orthogonal Kronecker modules") -- the strongest possible pins:
full rank-13 / rank-15 characters whose generic HN factors the paper states
exactly, with NEITHER factor semiexceptional.  Coordinates: the paper writes
``nu = eps E + phi F``; the package stores ``(f, s)``-vectors, so a paper slope
``eps E + phi F`` of rank ``r`` is the package ``c1 = (r*phi, r*eps)``.

Every pinned value is verified two ways (CLAUDE.md invariant 3): the exact
package-coordinate conversion is recomputed in comments/assertions here (sums,
discriminants, chi-orthogonality), and the algorithm output must reproduce the
paper tuple bit-for-bit -- at TWO distinct small perturbations ``eps`` of the
polarization, so a wall-crossing accident at one ``eps`` cannot fake agreement.
"""

from fractions import Fraction as F

import pytest

import bridgeland_stability.generic_hn as ghn
from bridgeland_stability.dlp_hirzebruch import (
    discriminant,
    hilbert_P,
    is_semiexceptional,
    total_slope,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi as surface_chi
from bridgeland_stability.generic_hn import generic_hn_factors, semistable_exists_hn
from bridgeland_stability.hn_filtration import HNRegion, hn_verdict
from bridgeland_stability.nonemptiness_rational import (
    VerdictStatus,
    hirzebruch_with_polarization,
    moduli_nonempty,
)
from bridgeland_stability.prioritary import delta_prioritary
from bridgeland_stability.varieties import P2, P1xP1, hirzebruch

F0 = P1xP1                                       # H = (1,1): m = 1 (the -K ray)
F1 = hirzebruch_with_polarization(1, (3, 2))     # H = (3,2): m = 1/2 (the -K ray)


def _surface_at_m(e, m):
    """The ample F_e whose polarization ray is H_m (a/b = m + e)."""
    q = m + e
    return hirzebruch_with_polarization(e, (q.numerator, q.denominator))


# --------------------------------------------------------------------------- #
# 1. The paper's two orthogonal-Kronecker pins (Sec. 1.5) -- neither factor     #
#    semiexceptional; reproduced at two distinct small eps.                     #
# --------------------------------------------------------------------------- #
#: F_1, m = 12/7 + eps:  v = (13, 3/13 E + 6/13 F, 98/169)
#:   -> (2, 1/2 E, 5/8)  and  (11, 2/11 E + 6/11 F, 65/121).
#: Package conversion: v = (13, (6, 3), -13/2); factors (2, (0, 1), -3/2) and
#: (11, (6, 2), -5).  Sum check: 2+11=13, (0,1)+(6,2)=(6,3), -3/2-5=-13/2.
_F1_PIN_V = (13, (6, 3), F(-13, 2))
_F1_PIN_FACTORS = ((2, (0, 1), F(-3, 2)), (11, (6, 2), F(-5)))
#: F_0, m = 25/9 + eps:  v = (15, 1/5 E + 1/3 F, 3/5)
#:   -> (2, 1/2 E - 1/2 F, 3/4)  and  (13, 2/13 E + 6/13 F, 90/169).
#: Package: v = (15, (5, 3), -8); factors (2, (-1, 1), -2) and (13, (6, 2), -6).
_F0_PIN_V = (15, (5, 3), F(-8))
_F0_PIN_FACTORS = ((2, (-1, 1), F(-2)), (13, (6, 2), F(-6)))


@pytest.mark.parametrize("eps", [F(1, 70), F(1, 700)])
def test_paper_kronecker_pin_F1(eps):
    S = _surface_at_m(1, F(12, 7) + eps)
    r, c1, ch2 = _F1_PIN_V
    got = generic_hn_factors(r, c1, ch2, S)
    assert got == _F1_PIN_FACTORS
    _verify_filtration_identities(got, (r, c1, ch2), S)
    # the paper's point: NEITHER factor is semiexceptional (Kronecker pair)
    for (rf, c1f, ch2f) in got:
        assert not is_semiexceptional(SurfaceBundle(rf, c1f, ch2f), S)
    # paper discriminants: Delta = 5/8 and 65/121
    assert discriminant(SurfaceBundle(*got[0]), S) == F(5, 8)
    assert discriminant(SurfaceBundle(*got[1]), S) == F(65, 121)


@pytest.mark.parametrize("eps", [F(1, 90), F(1, 900)])
def test_paper_kronecker_pin_F0(eps):
    S = _surface_at_m(0, F(25, 9) + eps)
    r, c1, ch2 = _F0_PIN_V
    got = generic_hn_factors(r, c1, ch2, S)
    assert got == _F0_PIN_FACTORS
    _verify_filtration_identities(got, (r, c1, ch2), S)
    for (rf, c1f, ch2f) in got:
        assert not is_semiexceptional(SurfaceBundle(rf, c1f, ch2f), S)
    assert discriminant(SurfaceBundle(*got[0]), S) == F(3, 4)
    assert discriminant(SurfaceBundle(*got[1]), S) == F(90, 169)


def _verify_filtration_identities(factors, v, S):
    """The thm-HNcriterion identities, re-verified from the OUTSIDE with the
    package's own general-purpose chi/discriminant (not the module's inlined
    integer forms): sum = v; chi(w_i, w_j) = 0 for i < j; Delta_i >= 0."""
    r, c1, ch2 = v
    assert sum(w[0] for w in factors) == r
    assert tuple(sum(w[1][k] for w in factors) for k in (0, 1)) == tuple(c1)
    assert sum((w[2] for w in factors), F(0)) == ch2
    bundles = [SurfaceBundle(w[0], w[1], w[2]) for w in factors]
    for i in range(len(bundles)):
        assert discriminant(bundles[i], S) >= 0            # Bogomolov per factor
        for j in range(i + 1, len(bundles)):
            assert surface_chi(bundles[i], bundles[j], S) == 0   # lem-HNorthogonal


# --------------------------------------------------------------------------- #
# 2. The flagship and the boundary: the band decides.                          #
# --------------------------------------------------------------------------- #
def test_flagship_decides_to_length_one():
    """(2,(1,1),0) on F_0 = O(1,0) (+) O(0,1): the E13 re-audit R2 class.  The
    paper exhibits it as -K-semistable (Sec. 7, the example after
    cor-delPezzoKss) even though Delta = 1/4 < DLP_{-K} = 3/4.  The computed
    generic HN filtration must have length 1 -- existence, not UNKNOWN."""
    assert generic_hn_factors(2, (1, 1), F(0), F0) == ((2, (1, 1), F(0)),)
    assert semistable_exists_hn(2, (1, 1), F(0), F0) is True


def test_F1_boundary_decides():
    """(2,(0,0),-2) on F_1 (-K ray, m=1/2): Delta = 1 sits exactly ON the sharp
    envelope, the strict-inequality open question of E11-M6.  The algorithm
    decides EXISTENCE: no valid maximal destabilizing character passes the
    cor-algorithm iff.  Hand-verified rejections of the natural candidates:
    w1 = O has chi(O, I_{Z_2}) = 1 - 2 = -1 != 0 (condition (4) fails);
    w1 = (1,(0,0),-1) ties the reduced Hilbert polynomial with v itself
    (q_1 > q_v fails, both slope 0 with P(0)-Delta = 0)."""
    assert generic_hn_factors(2, (0, 0), F(-2), F1) == ((2, (0, 0), F(-2)),)


def test_rank_one_base_case():
    # c2 >= 0: the ideal-sheaf twist exists and is stable (length 1)
    assert generic_hn_factors(1, (0, 0), F(-2), F0) == ((1, (0, 0), F(-2)),)
    assert generic_hn_factors(1, (2, -1), F(-2), F0) == ((1, (2, -1), F(-2)),)
    # c2 < 0: not the character of any torsion-free sheaf: empty prioritary stack
    assert generic_hn_factors(1, (0, 0), F(1), F0) is None
    assert semistable_exists_hn(1, (0, 0), F(1), F0) is False


# --------------------------------------------------------------------------- #
# 3. The inlined exact forms, cross-pinned against the package's general ones. #
# --------------------------------------------------------------------------- #
def test_inlined_chi_and_qkey_match_package_forms():
    ctx = ghn._Ctx(F1)
    # all integral characters ON F_1 (c2 in Z with the [[0,1],[1,-1]] Gram)
    chars = [(2, (0, 1), F(-3, 2)), (11, (6, 2), F(-5)), (1, (0, 0), F(0)),
             (3, (1, 1), F(-1, 2)), (2, (1, 1), F(-1, 2)), (5, (-2, 3), F(-7, 2))]
    for w1 in chars:
        b1 = SurfaceBundle(w1[0], w1[1], w1[2])
        # q_key second component == P(nu) - Delta via the package general forms
        mu, tail = ctx.q_key(w1)
        nu = total_slope(b1)
        assert tail == hilbert_P(nu, F1) - discriminant(b1, F1)
        assert mu == F1.lattice.pairing(nu, F1.H) / F1.H[1]
        for w2 in chars:
            b2 = SurfaceBundle(w2[0], w2[1], w2[2])
            assert ctx.chi(w1, w2) == surface_chi(b1, b2, F1)


# --------------------------------------------------------------------------- #
# 4. Differential: the algorithm vs the envelope verdicts (independent          #
#    theorem routes -- CH Sec. 5 vs the DLP envelope machinery).               #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("S", [F0, F1], ids=["F0", "F1"])
def test_algorithm_never_contradicts_the_envelope(S):
    """Over an integral-c2 grid: wherever moduli_nonempty is PROVEN, the
    computed filtration must agree (length 1 <-> PROVEN_NONEMPTY; length >= 2
    or empty stack <-> PROVEN_EMPTY).  Where it is UNKNOWN the algorithm
    supplies the decision -- there is nothing to contradict."""
    lat = S.lattice
    checked = unknown_decided = 0
    for r in range(1, 4):
        for a in range(-2, 3):
            for b in range(-2, 3):
                c1 = (a, b)
                for c2 in range(-2, 5):
                    ch2 = F(1, 2) * lat.self_pairing(c1) - c2
                    v = moduli_nonempty(r, c1, ch2, S)
                    factors = generic_hn_factors(r, c1, ch2, S)
                    exists = factors is not None and len(factors) == 1
                    if v.status is VerdictStatus.PROVEN_NONEMPTY:
                        assert exists, (r, c1, ch2)
                    elif v.status is VerdictStatus.PROVEN_EMPTY:
                        assert not exists, (r, c1, ch2)
                    else:
                        unknown_decided += 1
                    checked += 1
    assert checked > 300
    assert unknown_decided > 0        # the band is genuinely exercised


def test_exists_implies_prioritary_iff(surface=F1):
    """Existence (computed length 1) implies Delta >= delta^p_n (Cor 4.17's
    necessary direction for the semistable subcase) over a spot grid."""
    n = 1                                          # ceil(1/2)
    lat = surface.lattice
    for r in range(2, 4):
        for a in range(-1, 2):
            for b in range(-1, 2):
                for c2 in range(0, 4):
                    ch2 = F(1, 2) * lat.self_pairing((a, b)) - c2
                    if not semistable_exists_hn(r, (a, b), ch2, surface):
                        continue
                    xi = SurfaceBundle(r, (a, b), ch2)
                    assert discriminant(xi, surface) >= delta_prioritary(
                        total_slope(xi), n, surface)


def test_paranoid_uniqueness_sweep():
    """cor-algorithm's iff makes the maximal destabilizing character UNIQUE:
    with the full-sweep flag on, any second distinct winner raises the tripwire.
    Run over a small grid rich in length->=2 classes."""
    old = ghn.PARANOID_UNIQUENESS
    ghn.PARANOID_UNIQUENESS = True
    try:
        lat = F0.lattice
        for r in range(2, 4):
            for a in range(-1, 2):
                for b in range(-1, 2):
                    for c2 in range(0, 3):
                        ch2 = F(1, 2) * lat.self_pairing((a, b)) - c2
                        generic_hn_factors(r, (a, b), ch2, F0)   # must not raise
    finally:
        ghn.PARANOID_UNIQUENESS = old


# --------------------------------------------------------------------------- #
# 5. hn_verdict integration: the band is decided, labels are earned.           #
# --------------------------------------------------------------------------- #
def test_hn_verdict_flagship_flips_to_S():
    v = hn_verdict(2, (F(1, 2), F(1, 2)), F(1, 4), F0)
    assert v.exists is True and v.generic_hn_length == 1
    assert v.region is HNRegion.S
    assert v.factors == ((2, (1, 1), F(0)),)
    assert v.certificate.rigor.name == "PROVEN"


def test_hn_verdict_earns_region_K_on_the_kronecker_pin():
    """The F_1 orthogonal-Kronecker class through the verdict layer: PROVEN
    empty, computed length 2, region K EARNED (a non-semiexceptional factor),
    the factors exhibited."""
    S = _surface_at_m(1, F(12, 7) + F(1, 70))
    xi = SurfaceBundle(*_F1_PIN_V)
    v = hn_verdict(13, total_slope(xi), discriminant(xi, S), S)
    assert v.exists is False
    assert v.generic_hn_length == 2
    assert v.region is HNRegion.K
    assert v.factors == _F1_PIN_FACTORS
    assert v.certificate.rigor.name == "PROVEN"
    assert "factors" in v.reason


@pytest.mark.parametrize("S", [F0, F1], ids=["F0", "F1"])
def test_unclassified_never_returned_on_scope(S):
    """E13-M3b totality: the verdict is True/False with a PROVEN certificate
    everywhere on the del Pezzo scope; UNCLASSIFIED never fires."""
    lat = S.lattice
    for r in range(1, 4):
        for a in range(-1, 2):
            for b in range(-1, 2):
                for c2 in range(-1, 4):
                    ch2 = F(1, 2) * lat.self_pairing((a, b)) - c2
                    xi = SurfaceBundle(r, (a, b), ch2)
                    v = hn_verdict(r, total_slope(xi), discriminant(xi, S), S)
                    assert v.exists in (True, False), (r, (a, b), ch2)
                    assert v.region is not HNRegion.UNCLASSIFIED
                    assert v.certificate.rigor.name == "PROVEN"


# --------------------------------------------------------------------------- #
# 6. Scope and domain guards.                                                  #
# --------------------------------------------------------------------------- #
def test_scope_guards():
    with pytest.raises(NotImplementedError):
        generic_hn_factors(2, (1, 1), F(0), P2)          # not an F_e
    with pytest.raises(ValueError):
        generic_hn_factors(2, (1, 1), F(0), hirzebruch(1))   # nef-and-big factory
    with pytest.raises(ValueError):
        generic_hn_factors(2, (1, 1), F(1, 2), F0)       # c2 = 1/2 not integral
    with pytest.raises(ValueError):
        generic_hn_factors(F(3, 2), (1, 1), F(0), F0)    # non-integral rank
    assert semistable_exists_hn(2, (1, 1), F(1, 2), F0) is False   # invalid: empty


def test_e_ge_2_direct_smoke():
    """The Sec. 5 algorithm is uniform in e; a light e = 2 smoke (full M3c will
    differential this against the E13-M1 reduction pi)."""
    S2 = hirzebruch_with_polarization(2, (5, 2))
    assert generic_hn_factors(1, (0, 0), F(-1), S2) == ((1, (0, 0), F(-1)),)
    f = generic_hn_factors(2, (0, 0), F(-2), S2)
    assert f is not None and 1 <= len(f) <= 4
