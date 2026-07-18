"""E15-M3 / G18+: the Conjecture A falsification harness.

arXiv:1907.06739's Sec. 1.5 conjecture ("orthogonal Kronecker only"): if the
generic ``H_m``-Harder-Narasimhan filtration of a character on ``F_e`` has MORE
than one non-semiexceptional factor, then the length is 2 and the two factors
are (integer) linear combinations of the blocks ``(E_3, E_4)`` resp.
``(E_1, E_2)`` of a full exceptional collection.

Two falsification targets, both observable since E13-M3b exhibits factors:

1. a generic filtration of length >= 3 with two or more non-semiexceptional
   factors -- an immediate counterexample (:func:`classify_generic_filtration`
   reports the pattern; the sweep script scans for it);
2. a length-2 both-non-semiexceptional pair admitting NO block decomposition
   within the searched collection family (:func:`block_decomposition`) -- a
   RANKED CANDIDATE only, since the search family is bounded: the paper's
   Sec. 8 collections

       O(-E - l F),  O,  O(F),  O(E - (l - 1 - e) F)        (l >= 3)

   twisted by an arbitrary line bundle (twisting preserves full exceptional
   collections).  The paper proves this family is a full exceptional
   collection; ``is_exceptional_collection`` (the necessary chi-orthogonality
   test, E11-M2) is asserted on every family member as a tripwire.

The membership test is exact: solve ``alpha ch E_i + beta ch E_j = v`` on the
``(r, c1)`` components over ``Q`` and verify the ``ch_2`` component, then
require ``alpha, beta`` INTEGERS (the paper's constructions realize the factors
as extensions/cokernels of integer multiples; a Q-but-not-Z solution is
reported as a near-miss, never as a match).

Exact ``Fraction`` arithmetic; stdlib-only at import.  Record: the E15-M3
section of ``docs/CORRECTIONS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Optional, Sequence, Tuple

from .chern import Number
from .varieties import Surface
from .exceptional_surface import SurfaceBundle, is_exceptional_collection
from .dlp_hirzebruch import hirzebruch_index, is_semiexceptional
from .generic_hn import generic_hn_factors
from .rigor import Certificate, Rigor

__all__ = [
    "BlockWitness",
    "FiltrationPattern",
    "block_decomposition",
    "classify_generic_filtration",
]

_HALF = Fraction(1, 2)

_CITATIONS = (
    "arXiv:1907.06739 Sec. 1.5 (the conjecture), Sec. 8 (the collection family and the "
    "Kronecker constructions), Ex. KroneckerF0/KroneckerF1",
    "E13-M3b generic_hn (factors exhibited); E11-M2 is_exceptional_collection",
)


def _Q(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


def _line(surface: Surface, c1: Tuple[int, int]) -> SurfaceBundle:
    """``ch O(D) = (1, D, <D,D>/2)``."""
    return SurfaceBundle(1, c1, _HALF * surface.lattice.self_pairing(c1))


def _collection(surface: Surface, e: int, l: int, twist: Tuple[int, int],
                swapped: bool = False) -> Tuple[SurfaceBundle, ...]:
    """The Sec. 8 full exceptional collection, twisted by ``O(twist)`` --
    package ``(f, s)`` coordinates (paper ``aE + bF`` -> ``(b, a)``):

        E1 = O(-l, -1),  E2 = O,  E3 = O(1, 0),  E4 = O(-(l - 1 - e), 1).

    For ``e = 0`` the ruling-swap automorphism of ``P^1 x P^1`` (exchange the
    two fibre classes, i.e. the package coordinates) carries full exceptional
    collections to full exceptional collections; ``swapped=True`` returns that
    image family, which the search must also cover.
    """
    tf, ts = twist
    coords = (
        (-l + tf, -1 + ts),
        (tf, ts),
        (1 + tf, ts),
        (-(l - 1 - e) + tf, 1 + ts),
    )
    if swapped:
        coords = tuple((y, x) for (x, y) in coords)
    return tuple(_line(surface, c) for c in coords)


def _in_z_span(v: Tuple[int, Tuple[int, int], Fraction], A: SurfaceBundle,
               B: SurfaceBundle) -> Optional[Tuple[Fraction, Fraction, bool]]:
    """Solve ``alpha ch A + beta ch B = v`` exactly, using the rank and the
    ``s``-component of ``c1`` (checking the ``f``-component and ``ch_2``), or
    the ``f``-component where the first system is singular.  Returns
    ``(alpha, beta, integral)`` for a Q-solution, ``None`` if none exists.
    """
    r, c1, ch2 = v
    for idx in (1, 0):                            # s-component first, then f
        det = A.r * B.c1[idx] - B.r * A.c1[idx]
        if det == 0:
            continue
        alpha = Fraction(r * B.c1[idx] - B.r * c1[idx], det)
        beta = Fraction(A.r * c1[idx] - r * A.c1[idx], det)
        other = 1 - idx
        if alpha * A.c1[other] + beta * B.c1[other] != c1[other]:
            return None
        if alpha * A.ch2 + beta * B.ch2 != ch2:
            return None
        return alpha, beta, alpha.denominator == 1 and beta.denominator == 1
    return None                                   # both systems singular


@dataclass(frozen=True)
class BlockWitness:
    """A found block decomposition: the collection parameters and the integer
    coefficients ``v1 = a ch E3 + b ch E4``, ``v2 = c ch E1 + d ch E2``."""

    l: int
    twist: Tuple[int, int]
    coeffs_v1: Tuple[Fraction, Fraction]
    coeffs_v2: Tuple[Fraction, Fraction]
    certificate: Certificate


def block_decomposition(v1: Tuple[int, Sequence[int], Number],
                        v2: Tuple[int, Sequence[int], Number],
                        surface: Surface,
                        l_max: int = 10,
                        twist_max: int = 8) -> Optional[BlockWitness]:
    """Search the (bounded) Sec. 8 collection family for a block decomposition
    of the ordered factor pair ``(v1, v2)``: ``v1`` in the Z-span of
    ``(ch E3, ch E4)`` and ``v2`` in the Z-span of ``(ch E1, ch E2)``.

    Returns the first witness found, or ``None`` -- which, per the module
    docstring, is a SEARCH-BOUNDED statement (a ranked candidate for the
    conjecture's falsification, never a counterexample claim).
    """
    e = hirzebruch_index(surface)
    v1n = (int(v1[0]), tuple(int(x) for x in v1[1]), _Q(v1[2]))
    v2n = (int(v2[0]), tuple(int(x) for x in v2[1]), _Q(v2[2]))
    # l ranges BELOW the paper's l >= 3 as well: their bound served the Sec. 8
    # stability analysis, not collection-fullness.  Quadruples failing the
    # chi-orthogonality test are silently skipped (they are not exceptional
    # collections, hence outside the conjecture's quantifier); an orthogonal
    # length-4 collection on a del Pezzo is full (Kuleshov-Orlov).  The paper
    # range l >= 3 is orthogonal by the paper itself -- pinned in the tests.
    variants = (False, True) if e == 0 else (False,)
    for l in range(-2, l_max + 1):
        for tf in range(-twist_max, twist_max + 1):
            for ts in range(-twist_max, twist_max + 1):
              for swapped in variants:
                E1, E2, E3, E4 = _collection(surface, e, l, (tf, ts), swapped)
                s1 = _in_z_span(v1n, E3, E4)
                if s1 is None or not s1[2]:
                    continue
                s2 = _in_z_span(v2n, E1, E2)
                if s2 is None or not s2[2]:
                    continue
                if not is_exceptional_collection([E1, E2, E3, E4], surface):
                    continue                       # not an exceptional collection
                return BlockWitness(
                    l=l, twist=(tf, ts),
                    coeffs_v1=(s1[0], s1[1]), coeffs_v2=(s2[0], s2[1]),
                    certificate=Certificate(
                        rigor=Rigor.PROVEN, hypotheses=(
                            "collection is the Sec. 8 family (full, per the paper) "
                            "twisted by a line bundle",
                        ), citations=_CITATIONS,
                        note=f"block decomposition at l={l}, twist={(tf, ts)}"
                             + (", ruling-swapped" if swapped else "")))
    return None


@dataclass(frozen=True)
class FiltrationPattern:
    """The semiexceptionality pattern of a computed generic HN filtration."""

    length: Optional[int]                          # None: prioritary stack empty
    semiexceptional: Tuple[bool, ...]
    non_semiexceptional_count: int
    factors: Optional[tuple]
    conjecture_a_violation: bool                   # length >= 3 with >= 2 non-semiexc


def classify_generic_filtration(r: int, c1: Sequence[int], ch2: Number,
                                surface: Surface) -> FiltrationPattern:
    """Compute the generic HN factors (Sec. 5) and classify each factor's
    semiexceptionality -- the observable Conjecture A tests: a length >= 3
    filtration with two or more non-semiexceptional factors is an immediate
    counterexample (``conjecture_a_violation``); a length-2 both-non-semiexc
    pair feeds :func:`block_decomposition`."""
    facts = generic_hn_factors(r, c1, _Q(ch2), surface)
    if facts is None:
        return FiltrationPattern(None, (), 0, None, False)
    flags = tuple(
        bool(is_semiexceptional(SurfaceBundle(ri, c1i, ch2i), surface))
        for (ri, c1i, ch2i) in facts)
    non = sum(1 for f in flags if not f)
    return FiltrationPattern(
        length=len(facts), semiexceptional=flags,
        non_semiexceptional_count=non, factors=facts,
        conjecture_a_violation=(len(facts) >= 3 and non >= 2))
