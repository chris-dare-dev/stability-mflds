"""E13-M2 / G18: the prioritary sharp bound ``delta^p_n(nu)`` on ``F_e``.

Primary source: arXiv:1907.06739 (Coskun-Huizenga, "Existence of semistable sheaves on
Hirzebruch surfaces") Sec. 2.4 (prioritary sheaves), Sec. 4.1-4.3 (the L-Gaeta
resolution + the direct-sum construction), Prop 4.15, Cor 4.17, Cor 4.18, Remark 1.4,
Remark 4.13.  See ``docs/CORRECTIONS.md`` Sec. 10 for the two-way evidence.

Every asserted ``delta^p`` value is anchored two ways: (1) hand-computed from the
Prop 4.15 master formula ``delta^p_n = max{1/2 lam1 (lam2 (e+2n-2) + lam3 (e+2n)), 0}``
on the triangle ``T``, and (2) cross-checked against
:func:`bridgeland_stability.dlp_hirzebruch.discriminant` of the assembled direct-sum
witness ``V`` (which the theorem says realizes the bound).  The package basis is
``(f, s) = (F, E)``, so ``eps = nu[1]`` (the ``s = E`` coeff) and ``phi = nu[0]`` (the
``f = F`` coeff) -- a swapped mapping fails the two-way ``discriminant`` pins loudly.
"""

import os
import subprocess
import sys
from fractions import Fraction as F

import pytest

from bridgeland_stability.prioritary import (
    delta_prioritary,
    prioritary_nonempty,
    generic_prioritary_index,
    delta_prioritary_bundle,
    _L0_and_psi,
)
from bridgeland_stability.dlp_hirzebruch import (
    discriminant,
    total_slope,
    dlp_envelope,
    is_del_pezzo_anticanonical,
)
from bridgeland_stability.nonemptiness_rational import hirzebruch_with_polarization
from bridgeland_stability.exceptional_surface import SurfaceBundle
from bridgeland_stability.varieties import P1xP1, P2

# --- oracle import (tests/oracle/; put tests/ on sys.path) ------------------- #
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _TESTS_DIR not in sys.path:
    sys.path.insert(0, _TESTS_DIR)

from oracle.dlp_reference import reference_delta_prioritary


# --------------------------------------------------------------------------- #
# Surfaces.  F_1 = the ample anticanonical del Pezzo (H = (3,2) = -K ray);      #
# P1xP1 = F_0 with H = (1,1) (also anticanonical del Pezzo).                    #
# --------------------------------------------------------------------------- #
F1 = hirzebruch_with_polarization(1, (3, 2))


def _ample_Fe(e):
    """A strictly-ample F_e; delta^p is H-independent, so any ample H works."""
    return P1xP1 if e == 0 else hirzebruch_with_polarization(e, (e + 1, 1))


# --------------------------------------------------------------------------- #
# 1-2. Prop 4.15 on the triangle T, pinned two ways (formula + discriminant).   #
# --------------------------------------------------------------------------- #
def test_delta_p_on_triangle_F0():
    """F_0, nu = -1/2 E - 1/4 F (package (f,s)=(-1/4,-1/2)), n=1: lam=(1/2,1/4,1/4),
    delta^p_1 = 1/2 * 1/2 * (1/4*0 + 1/4*2) = 1/8.  Witness V = O(-E)^2 (+) O (+) O(-F)
    = (4, (-1,-2), 0) has discriminant 1/8 (the two-way check)."""
    nu = (F(-1, 4), F(-1, 2))
    assert delta_prioritary(nu, 1, P1xP1) == F(1, 8)
    assert discriminant(SurfaceBundle(4, (-1, -2), F(0)), P1xP1) == F(1, 8)
    assert delta_prioritary(nu, 1, P1xP1) == discriminant(
        SurfaceBundle(4, (-1, -2), F(0)), P1xP1)


def test_delta_p_on_triangle_F1():
    """Same slope on F_1 (e=1), n=1: 1/2 * 1/2 * (1/4*1 + 1/4*3) = 1/4.  Witness
    V = (4, (-1,-2), -1) (ch2 = 1/2 * 2 * (-e) = -1) has discriminant 1/4."""
    nu = (F(-1, 4), F(-1, 2))
    assert delta_prioritary(nu, 1, F1) == F(1, 4)
    assert discriminant(SurfaceBundle(4, (-1, -2), F(-1)), F1) == F(1, 4)
    assert delta_prioritary(nu, 1, F1) == discriminant(
        SurfaceBundle(4, (-1, -2), F(-1)), F1)


# --------------------------------------------------------------------------- #
# 3. Remark 4.13: twist + dual invariance.                                     #
# --------------------------------------------------------------------------- #
def test_delta_p_twist_dual_invariant():
    """delta^p_1 = 1/8 on F_0 for the base slope (-1/4,-1/2), two integer twists of it,
    and its dual (1/4,1/2) (Remark 4.13: integer twists and the dual preserve delta^p)."""
    for nu in [(F(-1, 4), F(-1, 2)),      # base (eps,phi) = (-1/2, -1/4)
               (F(-13, 4), F(3, 2)),      # twist by (a,b)=(2,-3) in (eps,phi)
               (F(19, 4), F(-3, 2)),      # twist by (-1,5)
               (F(1, 4), F(1, 2))]:       # dual of base
        assert delta_prioritary(nu, 1, P1xP1) == F(1, 8), nu


# --------------------------------------------------------------------------- #
# 4. Def 4.11: integer section slope eps in Z => delta^p = 0 (any n).           #
# --------------------------------------------------------------------------- #
def test_delta_p_eps_integer_is_zero():
    """nu = 1/3 F + 2 E (package (f,s)=(1/3, 2)) has eps = 2 in Z, so delta^p_n = 0
    for every n (Def 4.11 / Remark 4.7)."""
    nu = (F(1, 3), F(2))
    for n in (-1, 0, 1, 2, 3, 5):
        assert delta_prioritary(nu, n, P1xP1) == 0, n
        assert delta_prioritary(nu, n, F1) == 0, n


# --------------------------------------------------------------------------- #
# 5. Example 4.12: n <= -e => delta^p_n(nu) = 0 for ALL nu.                     #
# --------------------------------------------------------------------------- #
def test_delta_p_zero_for_n_le_minus_e():
    """For n <= -e the coefficient e+2n <= 0, so the master formula is <= 0 everywhere
    and delta^p_n = 0 for every slope (Example 4.12)."""
    for r in (2, 3):
        for a in range(-3, 4):
            for b in range(-3, 4):
                nu = (F(b, r), F(a, r))
                # F_0 (e=0): n in {0, -1}
                assert delta_prioritary(nu, 0, P1xP1) == 0, nu
                assert delta_prioritary(nu, -1, P1xP1) == 0, nu
                # F_1 (e=1): n in {-1, -2}
                assert delta_prioritary(nu, -1, F1) == 0, nu
                assert delta_prioritary(nu, -2, F1) == 0, nu


# --------------------------------------------------------------------------- #
# 6. Remark 1.4 floor: delta^p >= 0 (strong Bogomolov).                         #
# --------------------------------------------------------------------------- #
def test_delta_p_nonnegative_sweep():
    """delta^p_n(nu) >= 0 over r in {2,3,4}, a,b in [-4,4], n in {1,2,3}, e in {0,1}."""
    for r in (2, 3, 4):
        for a in range(-4, 5):
            for b in range(-4, 5):
                nu = (F(a, r), F(b, r))
                for n in (1, 2, 3):
                    for S in (P1xP1, F1):
                        assert delta_prioritary(nu, n, S) >= 0, (nu, n, S.name)


# --------------------------------------------------------------------------- #
# 7. Monotonicity in n (self-consistency) + the exact F_1 sequence pin.         #
# --------------------------------------------------------------------------- #
def test_delta_p_monotonic_in_n():
    """delta^p_n(nu) <= delta^p_{n+1}(nu) over a sweep, and the exact sequence on F_1 at
    nu = -1/2 E - 1/4 F: [delta^p_n]_{n=-2..4} = [0, 0, 0, 1/4, 1/2, 3/4, 1]."""
    for r in (2, 3, 4):
        for a in range(-3, 4):
            for b in range(-3, 4):
                nu = (F(a, r), F(b, r))
                for S, ns in [(P1xP1, range(-2, 5)), (F1, range(-3, 5))]:
                    vals = [delta_prioritary(nu, n, S) for n in ns]
                    for lo, hi in zip(vals, vals[1:]):
                        assert lo <= hi, (nu, S.name, lo, hi)
    nu = (F(-1, 4), F(-1, 2))
    seq = [delta_prioritary(nu, n, F1) for n in range(-2, 5)]
    assert seq == [F(0), F(0), F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]


# --------------------------------------------------------------------------- #
# 8-9. Remark 1.4 (the load-bearing acceptance): on the certified-sharp del      #
# Pezzo anticanonical ray, delta^{mu-s} (= dlp_envelope) >= delta^p_2.           #
# For F_0 (H=H_1, m=1) and F_1 (H=H_{1/2}, m=1/2) both use ceil(m)+1 = 2.        #
# --------------------------------------------------------------------------- #
_SWEEP_F0 = [
    (F(0), F(0)), (F(1), F(0)), (F(0), F(1)), (F(1), F(1)), (F(1, 2), F(1, 2)),
    (F(-1, 2), F(1, 2)), (F(1, 2), F(-1, 2)), (F(3, 2), F(1, 2)), (F(-1), F(1, 2)),
    (F(2, 3), F(1, 3)), (F(1, 3), F(2, 3)), (F(-2, 3), F(1, 3)),
]
# Package (f,s)=(phi,eps).  The eps = +-1/2 rows are the TIGHT cases where
# delta^p_2 = 5/8 > 1/2 and dlp_envelope = 5/8 exactly (the bound is sharp there).
_SWEEP_F1 = [
    (F(-1), F(-1, 2)), (F(-1), F(1, 2)), (F(0), F(-1, 2)), (F(0), F(1, 2)),
    (F(1), F(-1, 2)), (F(1), F(1, 2)), (F(0), F(0)), (F(1), F(1)),
    (F(1, 2), F(1, 2)), (F(1, 3), F(1, 3)), (F(2, 3), F(1, 3)), (F(-1, 2), F(3, 2)),
]


def test_remark_1_4_vs_certified_sharp_envelope_F0():
    """On F_0 with H = (1,1) (= H_1, the anticanonical del Pezzo, m=1), the package's
    mu-stable delta^{mu-s}_1 = dlp_envelope.value.  Remark 1.4 gives
    delta^{mu-s}_1 >= delta^p_{ceil(1)+1} = delta^p_2 -- asserted over the sweep."""
    assert is_del_pezzo_anticanonical(P1xP1)
    for nu in _SWEEP_F0:
        env = dlp_envelope(nu, P1xP1)
        assert env.value >= delta_prioritary(nu, 2, P1xP1), (nu, env.value)


def test_remark_1_4_vs_certified_sharp_envelope_F1():
    """On F_1 with H = (3,2) (= H_{1/2}, anticanonical del Pezzo, m=1/2),
    delta^{mu-s}_{1/2} = dlp_envelope.value >= delta^p_{ceil(1/2)+1} = delta^p_2
    (Remark 1.4).  The eps = +-1/2 slopes make this tight (delta^p_2 = 5/8 > 1/2), so
    the assertion is not vacuously satisfied by the 1/2 floor.

    Scope note: off the anticanonical ray dlp_envelope is only a LOWER bound for
    delta^{mu-s}, so ``delta^p <= dlp_envelope`` is NOT a theorem there and is not
    asserted -- both F_0 and F_1 here are on the certified-sharp -K ray."""
    assert is_del_pezzo_anticanonical(F1)
    saw_tight = False
    for nu in _SWEEP_F1:
        env = dlp_envelope(nu, F1)
        dp = delta_prioritary(nu, 2, F1)
        assert env.value >= dp, (nu, env.value, dp)
        if dp > F(1, 2):
            saw_tight = True
    assert saw_tight, "the F_1 sweep must include a slope with delta^p_2 > 1/2 (teeth)"


# --------------------------------------------------------------------------- #
# 10. Cor 4.17: P_{F,H_n}(v) != empty  <=>  Delta >= delta^p_n(nu).             #
# --------------------------------------------------------------------------- #
def test_prioritary_nonempty_characterization():
    """Cor 4.17 on F_1, nu = -1/2 E - 1/4 F, n=1 (delta^p = 1/4).  The boundary is
    inclusive (Delta = delta^p -> nonempty); above -> nonempty; below -> empty.  The
    Delta values are chosen so the character (r, nu, Delta) is a genuine integral class
    (c2 in Z); a non-integral class is out of Cor 4.17's domain and raises."""
    nu = (F(-1, 4), F(-1, 2))
    assert delta_prioritary(nu, 1, F1) == F(1, 4)
    # Delta in (1/4)Z gives an integral character here (c2 = 4*Delta).
    assert prioritary_nonempty(4, nu, F(1, 4), 1, F1) is True    # Delta = delta^p (inclusive)
    assert prioritary_nonempty(4, nu, F(1, 2), 1, F1) is True    # Delta > delta^p
    assert prioritary_nonempty(4, nu, F(0), 1, F1) is False      # Delta < delta^p
    # non-integral character (c2 = 4*3/8 = 3/2 not in Z) is refused
    with pytest.raises(ValueError):
        prioritary_nonempty(4, nu, F(3, 8), 1, F1)
    # non-integral c1 (3*nu = (-3/4,-3/2) not in NS) is refused
    with pytest.raises(ValueError):
        prioritary_nonempty(3, nu, F(1, 4), 1, F1)
    # negative Delta violates the Bogomolov hypothesis of Cor 4.17
    with pytest.raises(ValueError):
        prioritary_nonempty(4, nu, F(-1), 1, F1)
    # consistency with delta_prioritary over a sweep of genuine integral characters
    for r in (2, 3, 4):
        for a in range(-3, 4):
            for b in range(-3, 4):
                for c2 in range(-2, 5):
                    c1 = (a, b)
                    ch2 = F(1, 2) * F1.lattice.self_pairing(c1) - c2
                    xi = SurfaceBundle(r, c1, ch2)
                    nu2, D = total_slope(xi), discriminant(xi, F1)
                    if D < 0:
                        continue
                    for n in (1, 2):
                        assert prioritary_nonempty(r, nu2, D, n, F1) is (
                            D >= delta_prioritary(nu2, n, F1))


# --------------------------------------------------------------------------- #
# 11. Cor 4.18 / Example 4.9 / Figure 2: the generic prioritary index (anchored). #
# --------------------------------------------------------------------------- #
def test_generic_prioritary_index_figure2():
    """arXiv:1907.06739 Example 4.9 / Figure 2: nu = 1/2 E + 1/3 F (package
    (f,s)=(1/3,1/2)), Delta = 11/10, e = 1.  psi = -97/60, L_0 = ceil(eps) E +
    ceil(psi) F = (1,-1) (matching the Figure 2 caption (a0,b0)=(1,-1)), rho_gen = 4.
    (psi and rho_gen=4 are stated in Example 4.9's text; the caption carries only
    nu, Delta, e, (a0,b0).)"""
    nu = (F(1, 3), F(1, 2))
    assert generic_prioritary_index(nu, F(11, 10), F1) == 4
    L0, psi = _L0_and_psi(F(1, 2), F(1, 3), F(11, 10), 1)
    assert psi == F(-97, 60)
    assert L0 == (1, -1)                       # (a0, b0) in (E, F) coords
    # eps in Z is out of Cor 4.18's domain
    with pytest.raises(ValueError):
        generic_prioritary_index((F(1, 3), F(2)), F(1), F1)
    # Delta < 0 is out of Cor 4.18's domain too: the prioritary stack P_F(v) is
    # nonempty iff Delta >= 0 (Walter), so rho_gen is undefined below the Bogomolov
    # floor -- the formula used to return a meaningless -4 here (E13 re-audit R3).
    with pytest.raises(ValueError):
        generic_prioritary_index((F(1, 3), F(1, 2)), F(-1), P1xP1)


def test_non_integral_rank_is_never_truncated():
    """E13 re-audit R3: r = Fraction(3,2) used to be int()-truncated to r = 1, making
    prioritary_nonempty answer Cor 4.17 -- and delta_prioritary_bundle mint a Prop 4.15
    witness -- for a DIFFERENT character than the caller supplied."""
    with pytest.raises(ValueError, match="positive integer"):
        prioritary_nonempty(F(3, 2), (F(0), F(0)), F(0), 1, P1xP1)
    with pytest.raises(ValueError, match="positive integer"):
        delta_prioritary_bundle((F(0), F(-1, 2)), 1, F(3, 2), P1xP1)
    # An integral-valued Fraction is fine (4/2 == 2 in Z).
    assert prioritary_nonempty(F(4, 2), (F(0), F(0)), F(0), 1, P1xP1) is True


# --------------------------------------------------------------------------- #
# 12. Generalized two-way: discriminant of the Prop 4.15 witness == delta^p.     #
# --------------------------------------------------------------------------- #
def test_delta_prioritary_bundle_matches_discriminant():
    """Over an in-T sweep with integral multiplicities A,B,C = r*lam_i, the assembled
    direct-sum witness V = delta_prioritary_bundle(nu,n,r,S) has
    discriminant(V, S) == delta_prioritary(nu, n, S) whenever the latter is > 0."""
    checked = 0
    for S in (P1xP1, F1):
        for r in (2, 3, 4, 5, 6):
            for n in (1, 2, 3):
                for A in range(1, r):               # A = r*lam1 = -r*eps
                    for C in range(0, r):           # C = r*lam3
                        if A + C > r:
                            continue                # B = r-A-C would be negative
                        eps = F(-A, r)
                        phi = -(n - 1) * eps - F(C, r)   # lam3 = C/r
                        nu = (phi, eps)             # package (f,s)
                        dp = delta_prioritary(nu, n, S)
                        if dp > 0:
                            V = delta_prioritary_bundle(nu, n, r, S)
                            assert discriminant(V, S) == dp, (S.name, r, n, nu)
                            assert V.r == r
                            checked += 1
    assert checked > 100                            # the sweep actually exercised the path


# --------------------------------------------------------------------------- #
# 13. Differential regression against the independent oracle reference.          #
# --------------------------------------------------------------------------- #
def test_differential_against_oracle_reference():
    """The package delta_prioritary agrees with the append-only, no-package-import
    oracle reference_delta_prioritary over e in {0,1,2,3}, r in {1,2,3}, a,b in [-3,3],
    n in {1,2,3}.  The oracle takes (eps,phi) = (a/r, b/r); the package takes
    (f,s) = (phi,eps) = (b/r, a/r) -- the basis flip is part of the assertion."""
    for e in (0, 1, 2, 3):
        S = _ample_Fe(e)
        for r in (1, 2, 3):
            for a in range(-3, 4):
                for b in range(-3, 4):
                    for n in (1, 2, 3):
                        pkg = delta_prioritary((F(b, r), F(a, r)), n, S)
                        ref = reference_delta_prioritary(F(a, r), F(b, r), n, e)
                        assert pkg == ref, (e, r, a, b, n, pkg, ref)


# --------------------------------------------------------------------------- #
# 14-15. Scope + exactness guards.                                             #
# --------------------------------------------------------------------------- #
def test_delta_prioritary_off_hirzebruch_raises():
    """P^2 has no F_e fiber/section Gram, so delta^p is undefined there
    (NotImplementedError via hirzebruch_index)."""
    with pytest.raises(NotImplementedError):
        delta_prioritary((F(0),), 1, P2)
    with pytest.raises(NotImplementedError):
        prioritary_nonempty(2, (F(0),), F(0), 1, P2)
    with pytest.raises(NotImplementedError):
        generic_prioritary_index((F(1, 2),), F(1), P2)


def test_prioritary_module_is_stdlib_only():
    """prioritary.py must not pull matplotlib/plotly at import time (subprocess, so the
    check is not polluted by other tests that loaded viz)."""
    child = r'''
import sys
class _Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in {"matplotlib", "plotly"}:
            raise ImportError("blocked: " + name)
        return None
sys.meta_path.insert(0, _Blocker())
import bridgeland_stability.prioritary
assert "matplotlib" not in sys.modules
assert "plotly" not in sys.modules
print("STDLIB_ONLY_OK")
'''
    r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "STDLIB_ONLY_OK" in r.stdout


def test_delta_prioritary_returns_exact_fraction():
    """Every delta_prioritary return is an exact ``Fraction`` -- no float leaks into the
    core value (delta^p is purely rational, unlike the sqrt-bearing Gieseker bound)."""
    for r in (2, 3):
        for a in range(-3, 4):
            for b in range(-3, 4):
                for n in (1, 2):
                    for S in (P1xP1, F1):
                        v = delta_prioritary((F(a, r), F(b, r)), n, S)
                        assert isinstance(v, F) and not isinstance(v, float)
