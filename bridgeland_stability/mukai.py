"""K3 Mukai lattice, Mukai vectors, pairing, and wall-type classification.

For a sheaf E on a K3 surface, the Mukai vector is ``v(E) = ch(E) sqrt(td) =
(r, c1, ch2 + r)`` (since ``sqrt(td_K3) = (1, 0, 1)``).  We work with Picard
rank one, ``c1 = l H`` (``l`` in Z), so a Mukai vector is the integral triple
``(r, l, s)`` with ``s = ch2 + r``, and the Mukai pairing (with ``d = H^2``) is

    <(r,l,s), (r',l',s')> = d l l' - r s' - r' s ,    <v,w> = -chi(E,F).

Self-intersection ``v^2 = <v,v> = d l^2 - 2 r s``.  A class is *spherical* iff
``v^2 = -2`` (e.g. v(O_K3) = (1,0,1), v^2 = -2 = -chi(O,O)); *isotropic* iff
``v^2 = 0``.  For a generic stability condition the moduli space M(v) has
dimension ``v^2 + 2`` (Bayer-Macri, arXiv:1301.6968, Thm 2.15): a point if
``v^2 = -2``, dimension 2 if ``v^2 = 0``, positive-dimensional if ``v^2 >= 2``.

Wall classification (Bayer-Macri Thm 5.7) for the rank-2 hyperbolic lattice
``H`` spanned by ``v`` and a wall class ``w`` (spherical s^2=-2, isotropic
w^2=0):

  * DIVISORIAL if there is a spherical ``s in H`` with ``(s,v)=0`` (Brill-Noether),
    or an isotropic ``u in H`` with ``(u,v)=1`` (Hilbert-Chow), or ``(u,v)=2``
    (Li-Gieseker-Uhlenbeck);
  * else FLOPPING if ``v = a + b`` for two positive classes or there is a
    spherical ``s in H`` with ``0 < (s,v) <= v^2/2``;
  * otherwise a FAKE wall (totally semistable) or not a wall.

(The brief's "delta^2 = -2/0/2" trichotomy is garbled: there is no "+2" class
type; +2 is only the wrong-sign self-pairing of the mislabeled (1,0,-1).  See
docs/CORRECTIONS.md.)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass(frozen=True)
class MukaiVector:
    """An algebraic Mukai vector ``(r, l, s)`` on a Picard-rank-1 K3 (``c1 = l H``)."""

    r: int
    l: int
    s: int

    @classmethod
    def from_chern(cls, r: int, l: int, ch2: int) -> "MukaiVector":
        """Build ``v(E) = (r, l, ch2 + r)`` from Chern data ``(r, c1 = lH, ch2)``."""
        return cls(r, l, ch2 + r)

    def __iter__(self):
        return iter((self.r, self.l, self.s))


def pairing(v: MukaiVector, w: MukaiVector, d: int) -> int:
    """Mukai pairing ``<v,w> = d l l' - r s' - r' s`` (``= -chi(E,F)``)."""
    return d * v.l * w.l - v.r * w.s - w.r * v.s


def self_pairing(v: MukaiVector, d: int) -> int:
    """``v^2 = <v,v> = d l^2 - 2 r s``."""
    return d * v.l * v.l - 2 * v.r * v.s


def moduli_dim(v: MukaiVector, d: int) -> int:
    """``dim M_sigma(v) = v^2 + 2`` (Bayer-Macri Thm 2.15); negative => empty."""
    return self_pairing(v, d) + 2


def is_spherical(v: MukaiVector, d: int) -> bool:
    return self_pairing(v, d) == -2


def is_isotropic(v: MukaiVector, d: int) -> bool:
    return self_pairing(v, d) == 0


def _combo(a: int, b: int, v: MukaiVector, w: MukaiVector) -> MukaiVector:
    return MukaiVector(a * v.r + b * w.r, a * v.l + b * w.l, a * v.s + b * w.s)


@dataclass
class WallClassification:
    """Result of classifying a numerical wall of ``v`` defined by class ``w``."""

    wall_type: str  # 'divisorial', 'flopping', 'fake-or-none'
    subtype: Optional[str]  # 'brill-noether', 'hilbert-chow', 'li-gieseker-uhlenbeck', ...
    v_squared: int
    witness: Optional[Tuple[int, int, int]]  # the spherical/isotropic class found
    detail: str


def classify_wall(
    v: MukaiVector, w: MukaiVector, d: int, search: int = 8
) -> WallClassification:
    """Classify the wall of ``v`` through class ``w`` per Bayer-Macri Thm 5.7.

    Searches integer combinations ``a v + b w`` (``|a|,|b| <= search``) for
    spherical / isotropic witnesses and applies the theorem.  This is exact
    given the search bound; widen ``search`` for large invariants.
    """
    v2 = self_pairing(v, d)

    spherical: List[MukaiVector] = []
    isotropic: List[MukaiVector] = []
    for a in range(-search, search + 1):
        for b in range(-search, search + 1):
            if a == 0 and b == 0:
                continue
            x = _combo(a, b, v, w)
            x2 = self_pairing(x, d)
            if x2 == -2:
                spherical.append(x)
            elif x2 == 0 and (x.r, x.l, x.s) != (0, 0, 0):
                isotropic.append(x)

    # --- divisorial sub-cases ---
    for s in spherical:
        if pairing(s, v, d) == 0:
            return WallClassification(
                "divisorial", "brill-noether", v2, (s.r, s.l, s.s),
                "spherical s with (s,v)=0",
            )
    for u in isotropic:
        if pairing(u, v, d) == 1:
            return WallClassification(
                "divisorial", "hilbert-chow", v2, (u.r, u.l, u.s),
                "isotropic u with (u,v)=1",
            )
    for u in isotropic:
        if pairing(u, v, d) == 2:
            return WallClassification(
                "divisorial", "li-gieseker-uhlenbeck", v2, (u.r, u.l, u.s),
                "isotropic u with (u,v)=2",
            )

    # --- flopping ---
    if v2 > 0:
        for s in spherical:
            sv = pairing(s, v, d)
            # 0 < (s,v) <= v^2/2  with (s,v) integer  <=>  sv > 0 and 2*sv <= v^2
            if sv > 0 and 2 * sv <= v2:
                return WallClassification(
                    "flopping", "spherical", v2, (s.r, s.l, s.s),
                    "spherical s with 0 < (s,v) <= v^2/2",
                )

    return WallClassification(
        "fake-or-none", None, v2, None,
        "no Brill-Noether/Hilbert-Chow/LGU/flopping witness found in search box",
    )
