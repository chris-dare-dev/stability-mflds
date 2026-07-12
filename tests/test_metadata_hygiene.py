"""E12-M6 -- metadata + input-hygiene defects A7..A11 (permanent regressions).

Each defect is a metadata/provenance bug (not a P^2 class->verdict fact), so it is
recorded in ``tests/oracle/corpus.py`` as a ``METADATA_DEFECTS`` entry rather than
as a differential-sweep row, and is pinned here by a direct passing test:

* A7  ``dlp.moduli_nonempty`` reported ``moduli_dim`` in the doubled (brief)
      convention -- ``r^2(2*Delta_brief - 1)+1 = 4n`` on P^2[n] instead of the CHW
      ``r^2(2*Delta_CH - 1)+1 = 2n`` (arXiv:1401.1613 sec.2).
* A8  ``Surface.K_H`` was a hard-coded ``-2`` placeholder on every F_n; it is now
      the derived intersection number ``K.H = <K, H>``, correct on any polarization.
* A9  ``hirzebruch_with_polarization`` silently ``int()``-truncated a rational H
      (``(3/2,1) -> (1,1)``); it now raises.
* A10 ``SurfaceBundle.line_bundle`` accepted a fractional divisor class D (``O(D)``
      is undefined for non-integral D); it now raises.
* A11 ``canonical_class`` keyed on the Gram matrix, so a K3 with ``NS = U`` (same
      Gram as F_0) got the false ``K=(-2,-2)``; it now keys on ``surface.kind``.
"""

from fractions import Fraction as F

import pytest

from bridgeland_stability.dlp import moduli_nonempty
from bridgeland_stability.exceptional import Bundle
from bridgeland_stability.exceptional_surface import SurfaceBundle, canonical_class, chi
from bridgeland_stability.nonemptiness_rational import hirzebruch_with_polarization
from bridgeland_stability.nslattice import NSLattice
from bridgeland_stability.varieties import (
    P2,
    P1xP1,
    Surface,
    K3,
    abelian_surface,
    hirzebruch,
)


# --------------------------------------------------------------------------- #
# A7 -- moduli_dim in the CH convention: dim M(P^2[n]) = 2n, not 4n.           #
# --------------------------------------------------------------------------- #
def test_A7_moduli_dim_p2n_is_2n():
    # P^2[n] = (1,0,-n): mu=0, Delta_CH = n, delta(0)=1 so positive-dimensional
    # for n>=1.  dim = r^2(2 Delta_CH - 1)+1 = 1*(2n-1)+1 = 2n (was 4n, doubled).
    for n in range(1, 6):
        res = moduli_nonempty(Bundle(1, 0, F(-n)))
        assert res["discriminant"] == n
        assert res["positive_dimensional"] is True
        assert res["moduli_dim"] == 2 * n
        assert res["moduli_dim"] != 4 * n          # the old doubled value is gone


# --------------------------------------------------------------------------- #
# A8 -- Surface.K_H computed as <K, H> from the lattice, not a -2 placeholder. #
# --------------------------------------------------------------------------- #
def test_A8_K_H_computed_from_lattice():
    # Rank-1 shim rows: K.H = <K, H>.
    assert P2.K_H == -3                              # <(-3,),(1,)> on [[1]]
    assert P1xP1.K_H == -4                           # <(-2,-2),(1,1)> on [[0,1],[1,0]]
    # F_n via the factory (its degree-n polarization H = n f + s): K.H = -(n+2).
    assert hirzebruch(3).K_H == -5                   # was -2 placeholder
    assert hirzebruch(5).K_H == -7
    # F_n with an explicit strictly-ample H = a f + b s: K.H = (n-2)b - 2a.
    assert hirzebruch_with_polarization(3, (4, 1)).K_H == -7     # (3-2)*1 - 2*4 = -7
    assert hirzebruch_with_polarization(0, (1, 1)).K_H == -4     # (0-2)*1 - 2*1 = -4
    assert hirzebruch_with_polarization(0, (3, 2)).K_H == -10    # (0-2)*2 - 2*3 = -10

    # K numerically trivial families: K.H = 0.
    assert K3().K_H == 0
    assert abelian_surface().K_H == 0

    # K_H is LIVE: trivial_canonical reads it.
    assert hirzebruch(3).trivial_canonical is False
    assert K3().trivial_canonical is True
    assert abelian_surface().trivial_canonical is True


def test_A8_K_H_matches_pairing_directly():
    # Cross-check the property against a hand-evaluated <K, H>.
    S = hirzebruch_with_polarization(3, (4, 1))
    assert S.K_H == S.lattice.pairing(S.K, S.H)
    assert S.K == (-5, -2)                           # -(3+2) f - 2 s


def test_A8_K_validated_on_construction():
    # A mismatched-length or non-integral K is rejected by __post_init__.
    with pytest.raises(ValueError):
        Surface(name="bad-K-len", d=2, K=(0,), chi_O=1, picard_rank=2,
                kind="rational", H=(1, 1),
                ns_lattice=NSLattice(2, ((0, 1), (1, 0))))
    with pytest.raises(ValueError):
        Surface(name="bad-K-frac", d=1, K=(F(1, 2),), chi_O=1, kind="P2")


# --------------------------------------------------------------------------- #
# A9 -- hirzebruch_with_polarization rejects a non-integral polarization.     #
# --------------------------------------------------------------------------- #
def test_A9_non_integral_polarization_raises():
    with pytest.raises(ValueError):
        hirzebruch_with_polarization(0, (F(3, 2), 1))
    # The integral ray representative is accepted; d = 2ab - n b^2 = 2*3*2 = 12.
    assert hirzebruch_with_polarization(0, (3, 2)).d == 12


# --------------------------------------------------------------------------- #
# A10 -- line_bundle rejects a fractional divisor class.                       #
# --------------------------------------------------------------------------- #
def test_A10_line_bundle_rejects_fractional_divisor():
    with pytest.raises(ValueError):
        SurfaceBundle.line_bundle((F(1, 2),), P2)        # O_{P^2}(H/2): no such bundle
    with pytest.raises(ValueError):
        SurfaceBundle.line_bundle((F(1, 2), 0), P1xP1)
    # Integral divisors still build the pinned classes.
    assert SurfaceBundle.line_bundle((1, 1), P1xP1).ch2 == F(1)
    assert SurfaceBundle.line_bundle((1, 0), P1xP1).ch2 == F(0)


# --------------------------------------------------------------------------- #
# A11 -- canonical_class keys on surface.kind, not the Gram matrix.            #
# --------------------------------------------------------------------------- #
def test_A11_canonical_class_returns_stored_K_not_gram_inferred():
    # A K3 carrying NS = U (Gram [[0,1],[1,0]], identical to F_0).  The OLD
    # Gram-pattern match returned a FALSE K=(-2,-2) inferred from the shared Gram.
    # Since K is now a stored field (A8), canonical_class returns the TRUE stored
    # K=(0,0) and never Gram-infers -- which is the real content of A11.  (An
    # interim E12-M6 revision keyed on kind and raised here; that broke chi() /
    # euler_gram() on K3/abelian, restored in the E12 code review -- see
    # docs/CORRECTIONS.md sec.8.)
    k3U = Surface(name="K3 with NS=U", d=2, K=(0, 0), chi_O=2,
                  picard_rank=2, kind="K3", H=(1, 1),
                  ns_lattice=NSLattice(2, ((0, 1), (1, 0))))
    assert canonical_class(k3U) == (F(0), F(0))     # stored K, NOT the Gram-inferred (-2,-2)


def test_A11_canonical_class_works_on_k3_abelian():
    # The general-purpose Riemann-Roch chi()/euler_gram() must stay defined on
    # numerically-trivial (K-trivial) surfaces: canonical_class returns (0,...).
    from bridgeland_stability.varieties import K3, abelian_surface
    assert canonical_class(K3(2)) == (F(0),)
    assert canonical_class(abelian_surface(2)) == (F(0),)
    # chi(O,O,K3) = chi(O_K3) = 2 (K-trivial: canonical_term vanishes).
    O = SurfaceBundle(1, (0,), F(0))
    assert chi(O, O, K3(2)) == 2
    assert chi(O, O, abelian_surface(2)) == 0        # chi(O_abelian) = 0


def test_A11_canonical_class_unchanged_on_rational_surfaces():
    # Every surface that currently reaches canonical_class gets its OLD value.
    assert canonical_class(P2) == (F(-3),)
    assert canonical_class(P1xP1) == (F(-2), F(-2))
    assert canonical_class(hirzebruch(3)) == (F(-5), F(-2))     # -(3+2) f - 2 s
