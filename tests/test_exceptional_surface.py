"""E11-M1 / G14: SurfaceBundle + the Riemann-Roch Euler pairing chi.

Two-way verification for every asserted integer:

  (a) the closed form  chi(O(a,b), O(c,d)) = (c-a+1)(d-b+1)  for line bundles on
      P^1 x P^1 [Beilinson, PROVEN -- the exceptional collection
      <O, O(1,0), O(0,1), O(1,1)>], and
  (b) an independent re-derivation by the Riemann-Roch Euler formula
      chi(E,F) = r_E r_F chi(O_X) - (1/2)<r_E c1_F - r_F c1_E, K_X>
                 + (r_E ch2_F + r_F ch2_E - <c1_E, c1_F>),
      cross-checked on P^2 against the pinned exceptional.chi fixtures
      chi(O,O(3)) = 10, chi(O,O(-3)) = 1, chi(T(-1),T(-1)) = 1.

All arithmetic is exact Fraction; chi returns an exact int (never float).
"""

import subprocess
import sys

from fractions import Fraction as F

import pytest

from bridgeland_stability import exceptional
from bridgeland_stability.exceptional_surface import (
    SurfaceBundle,
    canonical_class,
    chi,
    euler_gram,
    p1xp1_collection,
    hirzebruch_collection,
    del_pezzo_collection,
    is_exceptional_collection,
)
from bridgeland_stability.varieties import P1xP1, P2, hirzebruch


# --------------------------------------------------------------------------
# The pinned P^1 x P^1 Euler Gram in the Beilinson basis.
# --------------------------------------------------------------------------
def test_p1xp1_euler_gram():
    basis = [SurfaceBundle.line_bundle(c, P1xP1) for c in [(0, 0), (1, 0), (0, 1), (1, 1)]]
    assert euler_gram(basis, P1xP1) == (
        (1, 2, 2, 4),
        (0, 1, 0, 2),
        (0, 0, 1, 2),
        (0, 0, 0, 1),
    )


def test_diagonal_unit():
    # Each generator is a line bundle: B = 0 and <c1,c1> is cancelled by ch2, so
    # chi(E,E) = chi(O_X) = 1.
    for c in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        E = SurfaceBundle.line_bundle(c, P1xP1)
        assert chi(E, E, P1xP1) == 1


def test_offdiagonal_witness():
    O10 = SurfaceBundle.line_bundle((1, 0), P1xP1)
    O01 = SurfaceBundle.line_bundle((0, 1), P1xP1)
    assert chi(O10, O01, P1xP1) == 0
    assert chi(O01, O10, P1xP1) == 0


def test_line_bundle_ch2():
    # ch2 = (1/2)<D,D>;  <(1,1),(1,1)> = 2 and <(1,0),(1,0)> = 0 (Gram [[0,1],[1,0]]).
    assert SurfaceBundle.line_bundle((1, 1), P1xP1).ch2 == F(1)
    assert SurfaceBundle.line_bundle((1, 0), P1xP1).ch2 == F(0)


def test_canonical_class_p1xp1():
    K = canonical_class(P1xP1)
    assert K == (F(-2), F(-2))
    # K.H reproduces the pinned scalar datum P1xP1.K_H = -4.
    assert P1xP1.lattice.pairing(K, P1xP1.H) == P1xP1.K_H == -4


def test_chi_matches_p2_exceptional():
    # P^2 line bundles as SurfaceBundle(1, (n,), n^2/2).
    O = SurfaceBundle(1, (F(0),), F(0))
    O3 = SurfaceBundle(1, (F(3),), F(9, 2))
    Om3 = SurfaceBundle(1, (F(-3),), F(9, 2))
    Tm1 = SurfaceBundle(2, (F(1),), F(-1, 2))  # T_{P^2}(-1): rank 2, c1=1, ch2=-1/2

    # Beilinson/Euler value chi(O,O(3)) = 10, tied to the pinned exceptional.chi.
    assert chi(O, O3, P2) == 10
    assert int(exceptional.chi(exceptional.Bundle.O(0), exceptional.Bundle.O(3))) == 10
    # chi(O, O(-3)) = P(-3) = 1.
    assert chi(O, Om3, P2) == 1
    assert int(exceptional.chi(exceptional.Bundle.O(0), exceptional.Bundle.O(-3))) == 1
    # Exceptional bundle: chi(T(-1), T(-1)) = 1.
    assert chi(Tm1, Tm1, P2) == 1
    assert int(exceptional.chi(exceptional.Bundle.T_minus1(), exceptional.Bundle.T_minus1())) == 1


def test_chi_returns_exact_int():
    O10 = SurfaceBundle.line_bundle((1, 0), P1xP1)
    O11 = SurfaceBundle.line_bundle((1, 1), P1xP1)
    v = chi(O10, O11, P1xP1)
    assert type(v) is int          # exact int -- not bool, not float, not Fraction
    assert not isinstance(v, bool)
    assert v == 2


def test_canonical_class_p2():
    # Picard-rank-1 shim path: K = (K.H / d) . H = (-3/1).(1,) = (-3,).
    assert canonical_class(P2) == (F(-3),)


def test_surface_bundle_frozen():
    E = SurfaceBundle(1, (0, 0), F(1, 3))
    with pytest.raises(Exception):
        E.r = 5                     # frozen dataclass
    # __post_init__ coerces c1 -> tuple(Fraction) and ch2 -> Fraction.
    assert isinstance(E.c1, tuple)
    assert all(isinstance(x, F) for x in E.c1)
    assert isinstance(E.ch2, F)


def test_chi_non_integral_raises():
    # A non-lattice class: c1=(0,0) but ch2=1/3 -> chi(E,E) = 1 + 2/3 = 5/3, not an
    # integer.  A genuine K-theory class always has integral chi, so this guards
    # the -> int contract honestly.
    bad = SurfaceBundle(1, (F(0), F(0)), F(1, 3))
    with pytest.raises(ValueError):
        chi(bad, bad, P1xP1)


def test_citation_provenance():
    # Guards the miscited-arXiv defect class: arXiv:1611.02674 is Coskun-Huizenga,
    # "Weak Brill-Noether for rational surfaces" -- NOT the (arXiv-less) Gokova
    # survey "The birational geometry of the moduli spaces of sheaves on P^2".
    # Verified on arXiv (title verbatim: "Weak Brill-Noether for rational surfaces").
    import bridgeland_stability.exceptional_surface as es

    doc = es.__doc__ or ""
    assert "arXiv:1611.02674" in doc
    assert "Weak Brill-Noether" in doc
    # The survey title must never be the one paired with this ID (the bug that was).
    assert "birational geometry" not in doc.lower()


def test_stdlib_only():
    # exceptional_surface.py must not pull matplotlib/plotly at import time.
    child = r'''
import sys
class _Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in {"matplotlib", "plotly"}:
            raise ImportError("blocked: " + name)
        return None
sys.meta_path.insert(0, _Blocker())
import bridgeland_stability.exceptional_surface
assert "matplotlib" not in sys.modules
assert "plotly" not in sys.modules
print("STDLIB_ONLY_OK")
'''
    r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "STDLIB_ONLY_OK" in r.stdout


# --------------------------------------------------------------------------
# E11-M2 / G14: per-family exceptional generators + necessary-only check.
#
# Every asserted value below is a LENGTH (= chi_top = rank K_0) or the
# NECESSARY numerical Euler-Gram condition (upper-triangular unipotent).  The
# F_n Grams are exactly re-derived: for line bundles chi(O(D),O(D')) =
# 1 + a(b+1) + b - (n/2) b(b+1) with (a,b) = D' - D on Gram [[0,1],[1,-n]],
# K = (-(n+2),-2); the collection D in {(0,0),(1,0),(0,1),(1,1)} makes every
# strictly-lower entry vanish for ALL n (diagonal = 1).
# --------------------------------------------------------------------------
def test_p1xp1_length():
    # chi_top(P^1 x P^1) = 4 = rank K_0 (Beilinson collection, transcription-free).
    assert len(p1xp1_collection()) == 4


def test_del_pezzo_length():
    # chi_top(dP_deg) = 12 - deg; for deg = 3 (cubic surface, Bl_6 P^2): 12 - 3 = 9.
    assert len(del_pezzo_collection(3)) == 9


def test_is_exceptional_collection_necessary():
    # P^1 x P^1 Beilinson Gram ((1,2,2,4),(0,1,0,2),(0,0,1,2),(0,0,0,1)) is
    # upper-triangular unipotent => the NECESSARY numerical condition holds.
    assert is_exceptional_collection(p1xp1_collection(), P1xP1) is True


def test_hirzebruch_collection_exceptional():
    # Exactly re-derived Grams (upper-triangular unipotent for all n):
    #   F_2 = ((1,2,0,2),(0,1,-2,0),(0,0,1,2),(0,0,0,1))
    #   F_3 = ((1,2,-1,1),(0,1,-3,-1),(0,0,1,2),(0,0,0,1))
    expected = {
        2: ((1, 2, 0, 2), (0, 1, -2, 0), (0, 0, 1, 2), (0, 0, 0, 1)),
        3: ((1, 2, -1, 1), (0, 1, -3, -1), (0, 0, 1, 2), (0, 0, 0, 1)),
    }
    for n in (1, 2, 3):
        coll = hirzebruch_collection(n)
        assert len(coll) == 4
        assert is_exceptional_collection(coll, hirzebruch(n)) is True
        if n in expected:
            assert euler_gram(coll, hirzebruch(n)) == expected[n]


def test_is_exceptional_collection_rejects_reversed():
    # Reversed Beilinson Gram ((1,0,0,0),(2,1,0,0),(2,0,1,0),(4,2,2,1)) is
    # LOWER-triangular: G[1][0] = 2 != 0 => necessary condition fails.
    reversed_coll = list(reversed(p1xp1_collection()))
    assert is_exceptional_collection(reversed_coll, P1xP1) is False


def test_del_pezzo_length_family():
    # (deg, length) with length = 12 - deg across the del Pezzo family.
    for deg, length in [(9, 3), (8, 4), (6, 6), (3, 9), (1, 11)]:
        assert len(del_pezzo_collection(deg)) == length


def test_hirzebruch_collection_requires_positive_n():
    with pytest.raises(ValueError):
        hirzebruch_collection(0)


def test_is_exceptional_collection_necessary_only_doc():
    # Delegation-honesty guard: the docstring must flag NECESSARY-only and G16.
    doc = is_exceptional_collection.__doc__ or ""
    assert "NECESSARY" in doc
    assert "G16" in doc
