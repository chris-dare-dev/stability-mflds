"""E11-M1 / G14: ``SurfaceBundle`` + the Riemann-Roch Euler pairing ``chi``.

The *numerical* layer of the rational-surface exceptional-collection program:
a K-theory class on a polarized surface ``(X, H)`` recorded as ``(r, c1, ch2)``
with ``c1`` an NS-vector in the coordinates of ``Surface.lattice``, plus the
exact Euler form

    chi(E, F) = integral over X of  ch(E)^v . ch(F) . td(X)

computed from the intersection form ``<.,.>`` (``Surface.lattice.pairing``), the
canonical class ``K_X in NS``, and ``chi(O_X)``.  With
``ch(E)^v = (r_E, -c1_E, ch2_E)`` and ``td(X) = (1, -K_X/2, chi(O_X))`` the
degree-two part integrates to the closed form

    chi(E,F) = r_E r_F chi(O_X)
             - (1/2) < r_E c1_F - r_F c1_E , K_X >
             + ( r_E ch2_F + r_F ch2_E - <c1_E, c1_F> ).                 (RR)

This *reduces exactly* to the pinned P^2 pairing :func:`bridgeland_stability.
exceptional.chi`: on P^2 (``chi(O_X)=1``, ``K_X=-3H``, ``d=H^2=1``) the middle
term becomes ``+(3/2)(r_E c1_F - r_F c1_E)``, matching ``exceptional.py`` line
for line (cross-checked in ``tests/test_exceptional_surface.py``).

Scope (numerical only).  This module ships the Euler *table* (the "uniform
deliverable" of G14).  Genuine sheaf-level exceptionality
(``Ext^{>0}(E_i, E_i) = 0``) is NOT asserted here -- it is delegated to G16
(E10).  ``euler_gram`` is the exact, citable Gram; the pinned P^1 x P^1 basis
``(O, O(1,0), O(0,1), O(1,1))`` is the Beilinson exceptional collection.

References
----------
* Coskun-Huizenga, "Weak Brill-Noether for Rational Surfaces", arXiv:1611.02674
  -- weak Brill-Noether on rational (Hirzebruch / del Pezzo) surfaces; the RR
  Euler form in NS coordinates.
* Levine-Zhang (Coskun-Huizenga circle), arXiv:1910.14060 -- non-emptiness on
  rational surfaces; the RR Euler pairing in NS coordinates.
* Beilinson: ``<O, O(1,0), O(0,1), O(1,1)>`` is a full exceptional collection on
  P^1 x P^1 [PROVEN]; the closed form ``chi(O(a,b),O(c,d)) = (c-a+1)(d-b+1)``.
  This transcription-free Beilinson collection is the ONLY per-family generator
  set pinned here without recourse to a secondary transcription.
* Rudakov et al., "Helices and Vector Bundles" (LMS Lecture Note Series 148,
  1990) -- exceptional collections and helices on Hirzebruch surfaces F_n; the
  length-4 collection ``<O, O(f), O(s), O(f+s)>`` used by
  :func:`hirzebruch_collection` [RESEARCH-light: read off Rudakov, not proven
  here -- :func:`is_exceptional_collection` checks only the NECESSARY numerical
  Euler-Gram condition].
* Kuleshov-Orlov, "Exceptional sheaves on Del Pezzo surfaces" (Izv. Ross. Akad.
  Nauk Ser. Mat. 58 (1994)) -- exceptional sheaves / collections on del Pezzo
  surfaces; the full collection has length ``12 - deg`` on
  ``dP_deg = Bl_{9-deg}(P^2)`` [RESEARCH-light].
* Orlov, "Projective bundles, monoidal transformations, and derived categories
  of coherent sheaves" (Izv. Math. 41 (1993)) -- the blow-up semiorthogonal
  decomposition ``D^b(Bl_p X) = <D^b(X), O_E>`` underpinning the del Pezzo
  length count (``pi^* <O,O(1),O(2)>`` from P^2 plus one ``O_{E_i}`` per
  blown-up point); the generator identities of :func:`del_pezzo_collection` are
  [RESEARCH-light] and only the LENGTH ``12 - deg`` is pinned as a test value.

Exact arithmetic.  Every quantity stays ``fractions.Fraction``; the ONLY
narrowing is the final ``int(val)`` in :func:`chi`, taken *after* an integrality
assertion (an Euler characteristic of a genuine K-theory class is an integer).
``int`` is exact -- no float ever appears in this module.  Stdlib-only at import
time (``fractions``, ``dataclasses``, ``typing`` + the stdlib-only intra-package
``nslattice`` / ``varieties``), preserving the zero-runtime-dependency invariant.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence, Tuple, Union

from .nslattice import NSLattice  # noqa: F401  (re-exported type; documents the pairing source)
from .varieties import Surface, P1xP1, hirzebruch

Number = Union[int, Fraction]
NSVector = Tuple[Fraction, ...]  # coordinates in surface.lattice's basis


def _Q(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


@dataclass(frozen=True)
class SurfaceBundle:
    """A K-theory class ``(r, c1, ch2)`` on a polarized surface.

    ``c1`` is an NS-vector in the coordinates of ``Surface.lattice`` (coerced to
    a tuple of ``Fraction`` in ``__post_init__``); ``ch2`` is a ``Fraction``.
    Frozen -- build new instances, never mutate.
    """

    r: int
    c1: NSVector
    ch2: Fraction

    def __post_init__(self) -> None:
        object.__setattr__(self, "c1", tuple(_Q(x) for x in self.c1))
        object.__setattr__(self, "ch2", _Q(self.ch2))

    @classmethod
    def line_bundle(cls, c1: Sequence[Number], surface: Surface) -> "SurfaceBundle":
        """The line bundle ``O(D)`` for ``D = c1``:  ``ch = (1, D, D^2/2)``.

        ``ch2 = (1/2) <D, D>`` from the surface's intersection form, so on
        P^1 x P^1 ``O(1,1)`` has ``ch2 = (1/2)*2 = 1`` and ``O(1,0)`` has
        ``ch2 = (1/2)*0 = 0``.
        """
        d = tuple(_Q(x) for x in c1)
        ch2 = Fraction(1, 2) * surface.lattice.self_pairing(d)
        return cls(1, d, ch2)

    def dual(self) -> "SurfaceBundle":
        """The dual class ``E^v = (r, -c1, ch2)`` (ch2 is even in the dual)."""
        return SurfaceBundle(self.r, tuple(-x for x in self.c1), self.ch2)


def canonical_class(surface: Surface) -> NSVector:
    """The canonical class ``K_X`` as an NS-vector in ``surface.lattice``'s basis.

    Two exact paths cover everything E11-M1/M2 needs:

    * Picard-rank-1 rows (``ns_lattice is None``): ``Pic(X) tensor Q = Q.H``, so
      ``K_X = (K.H / H^2) . H`` exactly -- ``P^2 -> (-3,)``, ``K3/abelian ->
      (0,)`` -- read off the genuine scalar datum ``surface.K_H`` and ``d``.
    * Explicit rank-2 fiber/section lattices with Gram ``[[0,1],[1,-n]]``
      (``P^1 x P^1``: ``n=0``; ``F_n``: ``n = -gram[1][1]``):
      ``K_{F_n} = -2 s - (n+2) f``, i.e. coordinate ``(-(n+2), -2)`` in the basis
      ``(f, s)``.  This path does NOT use ``surface.K_H`` -- ``hirzebruch``'s
      ``K_H = -2`` is a documented PLACEHOLDER (not ``K.H``; the true
      ``K_{F_n}.H = -(n+2)``).

    Torsion-canonical surfaces (Enriques / bielliptic, ``canonical_order != 0``)
    and any unrecognized lattice raise ``NotImplementedError``: their ``K_X`` is
    not determined by the numerical ``(K.H, d)`` data alone (needs the
    torsion-aware NS lattice, G12).
    """
    if surface.canonical_order != 0:
        raise NotImplementedError(
            f"{surface.name}: K_X is torsion (canonical_order="
            f"{surface.canonical_order}); its NS-vector is not determined by the "
            "numerical (K.H, d) data (needs the torsion-aware NS lattice, G12)."
        )
    if surface.ns_lattice is None:
        # Picard-rank-1: K_X = (K.H / d) . H  (exact on P^2, K3, abelian).
        scalar = Fraction(surface.K_H, surface.d)
        return tuple(scalar * _Q(h) for h in surface.H)
    gram = surface.ns_lattice.gram
    if surface.ns_lattice.rank == 2 and gram[0][0] == 0 and gram[0][1] == 1:
        # fiber/section normalization f^2=0, f.s=1, s^2=-n:  K = -(n+2) f - 2 s.
        n = -gram[1][1]
        return (Fraction(-(n + 2)), Fraction(-2))
    raise NotImplementedError(
        f"{surface.name}: canonical_class is implemented for Picard-rank-1 "
        "surfaces (K_X = (K.H/d).H) and rank-2 fiber/section lattices "
        f"[[0,1],[1,-n]]; got Gram {gram!r}."
    )


def chi(E: SurfaceBundle, F: SurfaceBundle, surface: Surface) -> int:
    """Exact Riemann-Roch Euler pairing ``chi(E,F)`` on ``surface`` -- see (RR).

    All arithmetic is ``Fraction``; the ``1/2`` is ``Fraction(1, 2)`` and every
    intersection number comes from ``surface.lattice.pairing`` (a ``Fraction``).
    The Euler characteristic of a genuine K-theory class is an integer, so the
    result is asserted integral (``ValueError`` otherwise -- this guards the
    ``-> int`` contract honestly) and returned as an exact ``int``.
    """
    lattice = surface.lattice
    K = canonical_class(surface)
    # B = r_E c1_F - r_F c1_E  (an NS-vector in the lattice basis).
    B = tuple(E.r * cf - F.r * ce for ce, cf in zip(E.c1, F.c1))
    rank_term = E.r * F.r * surface.chi_O
    canonical_term = -Fraction(1, 2) * lattice.pairing(B, K)
    deg2 = E.r * F.ch2 + F.r * E.ch2 - lattice.pairing(E.c1, F.c1)
    val = rank_term + canonical_term + deg2
    if val.denominator != 1:
        raise ValueError(
            f"chi(E,F) = {val} is not an integer; (r, c1, ch2) is not a genuine "
            "K-theory class on this surface (its c1/ch2 are not in the lattice)."
        )
    return int(val)


def euler_gram(
    bundles: Sequence[SurfaceBundle], surface: Surface
) -> Tuple[Tuple[int, ...], ...]:
    """The Euler Gram ``G[i][j] = chi(bundles[i], bundles[j], surface)``.

    The exact, citable "uniform deliverable" of G14.  For the Beilinson basis
    ``(O, O(1,0), O(0,1), O(1,1))`` on P^1 x P^1 this is the pinned
    upper-triangular unipotent matrix ``((1,2,2,4),(0,1,0,2),(0,0,1,2),
    (0,0,0,1))``.
    """
    return tuple(tuple(chi(E, F, surface) for F in bundles) for E in bundles)


# --------------------------------------------------------------------------
# E11-M2 / G14: per-family exceptional generators + necessary-only check.
#
# Each family's generators are HARD-CODED from its source (NO uniform
# recursion): the Beilinson collection on P^1 x P^1 is transcription-free
# [PROVEN]; the F_n and del Pezzo collections are read off Rudakov and
# Kuleshov-Orlov / Orlov respectively [RESEARCH-light].  The only pinned test
# values are the collection LENGTHS (= chi_top = rank K_0) and the NECESSARY
# numerical Euler-Gram condition of :func:`is_exceptional_collection`.
# --------------------------------------------------------------------------


def p1xp1_collection() -> "list[SurfaceBundle]":
    """Beilinson full exceptional collection ``<O, O(1,0), O(0,1), O(1,1)>`` on
    ``P^1 x P^1``.

    Transcription-free [PROVEN]; length ``4 = chi_top(P^1 x P^1) = rank K_0``.
    Generators are hard-coded (no uniform recursion).
    """
    return [SurfaceBundle.line_bundle(c, P1xP1) for c in [(0, 0), (1, 0), (0, 1), (1, 1)]]


def hirzebruch_collection(n: int) -> "list[SurfaceBundle]":
    """Rudakov-style full exceptional collection ``<O, O(f), O(s), O(f+s)>`` on
    ``F_n`` (``n >= 1``).

    Length ``4 = chi_top(F_n)``; specializes to the Beilinson collection at
    ``n = 0`` (``= P^1 x P^1``).  Generators are hard-coded from Rudakov (no
    uniform recursion); genuine sheaf-level exceptionality is NOT asserted --
    :func:`is_exceptional_collection` verifies only the NECESSARY Euler-Gram
    condition.
    """
    if n < 1:
        raise ValueError("hirzebruch_collection needs n>=1; use p1xp1_collection() for F_0")
    S = hirzebruch(n)
    return [SurfaceBundle.line_bundle(c, S) for c in [(0, 0), (1, 0), (0, 1), (1, 1)]]


def del_pezzo_collection(deg: int) -> "list[SurfaceBundle]":
    """Full exceptional collection of length ``12 - deg`` on
    ``dP_deg = Bl_{9-deg}(P^2)`` (Orlov blow-up decomposition; Kuleshov-Orlov
    exceptional sheaves).

    Only the LENGTH ``12 - deg`` is a pinned test value here (PROVEN =
    ``chi_top``); the exact generator identities are [RESEARCH-light]
    (Kuleshov-Orlov).  Returned as bundle DATA in the ``I_{1,k}`` basis
    ``(L, E_1..E_k)``; computing its Euler Gram on a dP surface is future work
    (needs a dP NS lattice + ``canonical_class(I_{1,k})`` extension, out of
    E11-M2 scope).
    """
    if not (1 <= deg <= 9):
        raise ValueError("del Pezzo degree must be in 1..9")
    k = 9 - deg
    rank = k + 1                              # NS = I_{1,k}: L^2=1, E_i^2=-1, L.E_i=0
    Z = (Fraction(0),) * rank
    # Orlov pullbacks pi*<O, O(1), O(2)> from P^2 :  ch2(O(mL)) = (1/2)<mL,mL> = m^2/2
    coll = [
        SurfaceBundle(1, Z, Fraction(0)),                                        # O_X
        SurfaceBundle(1, (Fraction(1),) + (Fraction(0),) * k, Fraction(1, 2)),   # O_X(L)
        SurfaceBundle(1, (Fraction(2),) + (Fraction(0),) * k, Fraction(2)),      # O_X(2L)
    ]
    for i in range(k):                        # exceptional-curve structure sheaves O_{E_i}
        c1 = tuple(Fraction(1) if j == i + 1 else Fraction(0) for j in range(rank))
        coll.append(SurfaceBundle(0, c1, Fraction(1, 2)))  # ch(O_{E_i})=(0,E_i,-E_i^2/2)=(0,E_i,1/2)
    return coll


def is_exceptional_collection(
    bundles: Sequence[SurfaceBundle], surface: Surface
) -> bool:
    """NECESSARY numerical condition ONLY: the Euler Gram
    ``G[i][j] = chi(E_i, E_j, surface)`` is upper-triangular unipotent (unit
    diagonal, strictly-lower triangle zero).

    Genuine sheaf-level exceptionality ``Ext^{>0}(E_i, E_i) = 0`` is delegated
    to G16 (E10) and is NEVER asserted by the core (E11-M2 acceptance).  A
    ``True`` return is therefore a numerical *necessary* certificate, not a
    proof of exceptionality.
    """
    G = euler_gram(bundles, surface)
    for i in range(len(G)):
        if G[i][i] != 1:
            return False
        for j in range(i):                    # strictly-lower triangle must vanish
            if G[i][j] != 0:
                return False
    return True
