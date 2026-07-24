"""Independent del Pezzo (F_0 / F_1) reference: exceptional CLASSES by K-theoretic mutation.

CHARTER -- why this oracle exists
=================================
The E15-M2 sweep (CORRECTIONS Sec. 24) enumerates the del Pezzo exceptional
characters by SCANNING every ``c1`` in ``[0,r)^2`` and testing each candidate with
the package's Drezet-Le Potier rank induction (``dlp_hirzebruch.is_stable_exceptional``).
The E15/E16 adversarial audit (Sec. 26.4) verified that scan only through rank 30
and flagged the rank-130 totals as never independently rerun -- a missed family
member is exactly the failure mode that would falsify the sweep's uniqueness
claim, so the enumeration needed a second, non-shared derivation.

This module is that derivation.  It generates the same objects from a DIFFERENT
THEOREM, and shares no code with the package:

  * it imports NOTHING from ``bridgeland_stability`` -- the NS Gram, the canonical
    class and the Riemann-Roch Euler pairing are all re-derived here (see below),
    so a transcription error in the package's ``chi`` or lattice cannot propagate;
  * generation is by MUTATION of full exceptional collections -- pure integer /
    ``Fraction`` linear algebra, ``L_A B = chi(A,B).A - B`` and
    ``R_B A = chi(A,B).B - A`` ([GR]) -- seeded on the line-bundle collection
    ``<O, O(1,0), O(0,1), O(1,1)>``;
  * completeness is [KO]/[R]'s transitivity theorem (on a del Pezzo surface the
    braid group acts transitively on full exceptional collections, and every
    exceptional bundle belongs to one), NOT the DLP envelope description;
  * NO envelope, NO ``delta``, NO rank induction, and no notion of stability
    appears anywhere in the generator.

The two enumerators can therefore agree only if both theorems are transcribed
correctly.  ``tests/test_delpezzo_mutation_oracle.py`` runs that differential.

THE EULER FORM (re-derived, not imported)
=========================================
On a surface ``X``, ``chi(E,F) = int ch(E)^v . ch(F) . td(X)`` with
``ch(E^v) = (r, -c1, ch2)`` and ``td(X) = (1, -K/2, chi(O_X))``.  Expanding,

    chi(E,F) = r r' chi(O_X) - (1/2) K.(r c1' - r' c1) + r ch2' + r' ch2 - c1.c1'

On ``F_e`` take the NS basis ``(F, E)`` -- fibre and ``(-e)``-section -- so
``F^2 = 0``, ``F.E = 1``, ``E^2 = -e``, giving the Gram ``[[0,1],[1,-e]]``.  Then
``K = -2E - (e+2)F = (-(e+2), -2)`` in those coordinates and ``K^2 = 8``, the
degree-8 del Pezzo value for ``e in {0,1}``.  ``chi(O_{F_e}) = 1``.  These are
asserted, not assumed: :func:`_self_pairing` of ``K`` is checked to be 8.

WHY THE HELIX IS OMITTED
========================
A helix shift twists one member by ``+-K``, moving its ``c1`` by ``+-r.K`` and so
leaving the twist class ``c1 mod r`` FIXED.  It therefore contributes no new
class while inflating the collection orbit enormously (measured: the orbit does
not terminate either way, but dropping the helix roughly halves the frontier).
Only the six pairwise mutations are used.

TERMINATION AND COMPLETENESS IN PRACTICE
========================================
The braid orbit of collections is infinite even modulo twist and under a rank
bound, so the walk cannot be run to exhaustion.  It is instead run to CLASS
SATURATION: it stops once ``sat_patience`` consecutive newly-discovered
collections have yielded no new class.  Two empirical guards back that up, both
recorded in CORRECTIONS Sec. 26.4:

  * the result is invariant under the intermediate-rank cap (caps 160 and 220
    give bit-identical rank-<=130 sets, with the front running out to rank 219);
  * it matches the production scan MEMBER FOR MEMBER on every range where that
    scan is affordable -- rank <= 15, <= 30 and <= 40.

Like the other oracles this uses exact ``int``/``Fraction`` arithmetic only: no
float, no ``math`` import, and it is written to be obviously-correct, not fast.

PRIMARY SOURCES (quoted, not recalled)
======================================
[GR]  Gorodentsev-Rudakov, "Exceptional vector bundles on projective spaces",
      Duke Math. J. 54 (1987), 115-130.  Mutation of exceptional pairs.
[KO]  Kuleshov-Orlov, "Exceptional sheaves on del Pezzo surfaces",
      Russian Acad. Sci. Izv. Math. 44 (1995), 479-513; also
      arXiv:alg-geom/9511016.  Exceptional sheaves on del Pezzo surfaces and
      transitivity of mutation on exceptional collections.
[R]   Rudakov, "Exceptional vector bundles on a quadric",
      Math. USSR Izv. 33 (1989), 115-138.  The ``P^1 x P^1 = F_0`` case.
[CH]  Coskun-Huizenga, arXiv:1907.06739, Sec. 11 -- the conjecture whose
      verified range the E15-M2 sweep extends.
"""

from collections import deque
from fractions import Fraction

__all__ = [
    "euler_pairing",
    "line_bundle",
    "standard_collection",
    "exceptional_classes",
]

_HALF = Fraction(1, 2)


# --------------------------------------------------------------------------
# NS lattice of F_e in the basis (F, E):  F^2 = 0, F.E = 1, E^2 = -e
# --------------------------------------------------------------------------

def _pairing(a, b, e):
    """Intersection pairing of two NS vectors on ``F_e``."""
    return a[0] * b[1] + a[1] * b[0] - e * a[1] * b[1]


def _self_pairing(a, e):
    return 2 * a[0] * a[1] - e * a[1] * a[1]


def _canonical(e):
    """``K = -2E - (e+2)F`` in the (F, E) basis."""
    return (-(e + 2), -2)


def euler_pairing(A, B, e):
    """``chi(A, B)`` on ``F_e`` by Riemann-Roch; classes are ``(r, c1, ch2)``."""
    r, c1, ch2 = A
    r2, c12, ch22 = B
    lin = (r * c12[0] - r2 * c1[0], r * c12[1] - r2 * c1[1])
    return (r * r2
            - Fraction(_pairing(_canonical(e), lin, e), 2)
            + r * ch22 + r2 * ch2
            - _pairing(c1, c12, e))


def line_bundle(c1, e):
    """The class of ``O(c1)``: rank 1, ``ch2 = c1^2 / 2``."""
    return (1, tuple(c1), Fraction(_self_pairing(c1, e), 2))


def standard_collection(e):
    """``<O, O(1,0), O(0,1), O(1,1)>`` -- a full exceptional collection on F_0/F_1."""
    return tuple(line_bundle(c, e) for c in ((0, 0), (1, 0), (0, 1), (1, 1)))


# --------------------------------------------------------------------------
# mutation + twist normalization
# --------------------------------------------------------------------------

def _twist(w, D, e):
    """Tensor ``w`` by the line bundle ``O(D)``."""
    r, c1, ch2 = w
    return (r,
            (c1[0] + r * D[0], c1[1] + r * D[1]),
            ch2 + _pairing(c1, D, e) + r * _HALF * _self_pairing(D, e))


def _normalize(coll, e):
    """Canonical representative of the common line-twist orbit of a collection."""
    for r1, c11, _ in coll:
        if r1:
            a = abs(r1)
            D = tuple(((c % a) - c) // r1 for c in c11)
            break
    else:
        raise ValueError("cannot twist-normalize an all-rank-zero collection")
    if D == (0, 0):
        return coll
    return tuple(_twist(w, D, e) for w in coll)


def _mutations(coll, e):
    """The six pairwise mutations.  (The helix is omitted -- see the charter.)"""
    for i in range(3):
        A, B = coll[i], coll[i + 1]
        x = euler_pairing(A, B, e)
        LAB = (x * A[0] - B[0],
               (x * A[1][0] - B[1][0], x * A[1][1] - B[1][1]),
               x * A[2] - B[2])
        RBA = (x * B[0] - A[0],
               (x * B[1][0] - A[1][0], x * B[1][1] - A[1][1]),
               x * B[2] - A[2])
        yield coll[:i] + (LAB, A) + coll[i + 2:]
        yield coll[:i] + (B, RBA) + coll[i + 2:]


# --------------------------------------------------------------------------
# the enumeration
# --------------------------------------------------------------------------

def exceptional_classes(e, rank_max, rank_cap=None, sat_patience=25_000,
                        node_cap=20_000_000, progress=None):
    """Twist classes ``(r, c1 mod r)`` of the exceptional bundles on ``F_e``,
    ``2 <= r <= rank_max``.

    ``rank_cap`` bounds the ranks allowed in an INTERMEDIATE collection (default
    ``3 * rank_max``); the result must be invariant under raising it.  The walk
    stops after ``sat_patience`` consecutive new collections yield no new class.

    Returns ``(classes, info)`` -- a ``frozenset`` and a telemetry dict.
    """
    if e not in (0, 1):
        raise ValueError("del Pezzo bases only: e in {0, 1}")
    assert _self_pairing(_canonical(e), e) == 8, "K^2 must be 8 on a degree-8 del Pezzo"
    if rank_cap is None:
        rank_cap = 3 * rank_max

    start = _normalize(standard_collection(e), e)
    seen = {start}
    queue = deque([start])
    classes = set()

    def collect(coll):
        for r, c1, _ in coll:
            if 2 <= r <= rank_max:
                classes.add((r, (c1[0] % r, c1[1] % r)))

    collect(start)
    nodes = 0
    last_growth = 0
    max_rank = 0
    while queue:
        coll = queue.popleft()
        nodes += 1
        if nodes > node_cap:
            status = "NODE_CAP"
            break
        before = len(classes)
        for nc in _mutations(coll, e):
            if any(abs(w[0]) > rank_cap for w in nc):
                continue
            nc = _normalize(nc, e)
            if nc in seen:
                continue
            seen.add(nc)
            collect(nc)
            m = max(abs(w[0]) for w in nc)
            if m > max_rank:
                max_rank = m
            queue.append(nc)
        if len(classes) > before:
            last_growth = nodes
        if nodes - last_growth > sat_patience:
            status = "SATURATED"
            break
        if progress is not None and nodes % 50_000 == 0:
            progress(nodes, len(seen), len(classes), max_rank)
    else:
        status = "EXHAUSTED"

    return frozenset(classes), {
        "status": status, "nodes": nodes, "collections": len(seen),
        "classes": len(classes), "max_rank": max_rank,
        "rank_cap": rank_cap, "sat_gap": nodes - last_growth,
    }
