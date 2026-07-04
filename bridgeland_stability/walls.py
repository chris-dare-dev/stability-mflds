"""Algorithm 3 - Bridgeland numerical walls in the (s, t) upper half-plane.

For a fixed Chern character ``v`` on a surface with ``d = H^2`` and a potential
sub/quotient object ``w``, the wall ``W(v, w)`` is the locus where
``arg Z_{s,t}(v) = arg Z_{s,t}(w)``.  Writing ``v = (r,c,e)``, ``w = (r',c',e')``
and the three 2x2 "Mukai-like" minors

    W_rc = r c' - r' c,   W_re = r e' - r' e,   W_ce = c e' - c' e,

the wall is, when ``W_rc != 0``, the semicircle centered on the s-axis at

    s0  = W_re / W_rc,
    rho^2 = s0^2 - 2 * W_ce / (d * W_rc) = (s0 - mu_v)^2 - 2 * Delta_v(CH),

and it is a genuine (real) wall iff ``rho^2 >= 0``.  When ``W_rc = 0`` (equal
Mumford slopes) the wall degenerates to the vertical line ``s = mu_v``.

This matches Coskun-Huizenga (survey, section 5: ``s = (mu1+mu2)/2 -
(Delta1-Delta2)/(mu1-mu2)``, ``rho^2 = (s-mu1)^2 - 2 Delta1``) and ABCH
(arXiv:1203.0316).  Cross-check: for P^2[2], ``v = ch(I_Z) = (1,0,-2)`` is
destabilized by ``O(-1) = (1,-1,1/2)`` along the unique wall with center
``-5/2`` and radius ``3/2`` -- reproduced exactly by ``numerical_wall``.

All centers/radii are exact (``Fraction``); ``radius`` casts to float only at
the very end for geometry/plotting.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import List, Optional, Tuple

from .chern import ChernChar, Number, Q
from .varieties import Surface


@dataclass(frozen=True)
class Wall:
    """A semicircular numerical wall: center ``(s0, 0)`` and ``radius^2 = rho^2``."""

    center: Fraction
    radius_sq: Fraction
    subobject: ChernChar
    v: ChernChar

    @property
    def is_real(self) -> bool:
        return self.radius_sq >= 0

    @property
    def radius(self) -> float:
        return math.sqrt(float(self.radius_sq)) if self.radius_sq >= 0 else float("nan")

    @property
    def s_left(self) -> float:
        return float(self.center) - self.radius

    @property
    def s_right(self) -> float:
        return float(self.center) + self.radius

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"Wall(center={self.center}, radius_sq={self.radius_sq}, "
            f"R={self.radius:.4f}, w={self.subobject})"
        )


@dataclass(frozen=True)
class VerticalWall:
    """The degenerate wall ``s = s_value`` (equal Mumford slopes, ``W_rc = 0``)."""

    s_value: Optional[Fraction]
    subobject: ChernChar
    v: ChernChar

    is_real = True
    radius = float("inf")


def _minors(v: ChernChar, w: ChernChar) -> Tuple[Fraction, Fraction, Fraction]:
    W_rc = v.r * w.c - w.r * v.c
    W_re = v.r * w.e - w.r * v.e
    W_ce = v.c * w.e - w.c * v.e
    return W_rc, W_re, W_ce


def numerical_wall(v: ChernChar, w: ChernChar, d: int):
    """The wall ``W(v, w)`` as a :class:`Wall` (semicircle) or :class:`VerticalWall`.

    The returned :class:`Wall` may have ``radius_sq < 0`` (``is_real`` False),
    meaning ``w`` does not give a genuine wall for ``v`` (excluded by BG).
    """
    W_rc, W_re, W_ce = _minors(v, w)
    if W_rc == 0:
        s_val = v.slope(d) if v.r != 0 else None
        return VerticalWall(s_val, w, v)
    s0 = Fraction(W_re, W_rc)
    rho_sq = s0 * s0 - Fraction(2 * W_ce, d * W_rc)
    return Wall(center=s0, radius_sq=rho_sq, subobject=w, v=v)


def _half_integers(lo: Fraction, hi: Fraction):
    """Yield all half-integers (elements of (1/2)Z) in ``[lo, hi]``."""
    start = math.ceil(2 * lo)
    end = math.floor(2 * hi)
    for k in range(start, end + 1):
        yield Fraction(k, 2)


def _bg_ok(ch: ChernChar, d: int) -> bool:
    """Bogomolov necessary condition for a (sub/quotient) class to bound a real wall.

    Rank 0 (torsion) classes are admissible iff nonzero; positive/negative rank
    classes must satisfy ``Delta >= 0``.
    """
    if ch.r == 0:
        return ch.c != 0 or ch.e != 0
    return ch.discriminant(d) >= 0


def compute_walls(
    v: ChernChar,
    surface: Surface,
    s_range: Tuple[Number, Number] = (-6, 6),
    rank_bound: int = 6,
    degree_bound: int = 6,
    include_torsion: bool = True,
    require_bg: bool = True,
) -> List[Wall]:
    """Enumerate real numerical walls for ``v`` with center in ``s_range``.

    Strategy (exact, bounded).  For each candidate sub-rank ``r'`` in
    ``[0, min(rank_bound, rank v)]`` (``include_torsion`` adds the ``r'=0``
    curve-class destabilizers), the Mukai-degree ``W_rc = r_v c' - r' c_v`` is
    swept over ``[-degree_bound, degree_bound]`` (the natural finiteness
    parameter, since the wall center is ``W_re/W_rc``); each value fixes ``c'``,
    and the ``e'`` placing the center inside ``s_range`` form a finite interval
    of half-integers.  A wall is kept iff ``rho^2 >= 0``, its center is in
    ``s_range``, and (when ``require_bg``) BOTH the subobject ``w`` and the
    quotient ``v - w`` satisfy Bogomolov.  Walls are de-duplicated by
    ``(center, radius_sq)`` and returned sorted by radius (largest first).

    This returns *numerical/potential* walls (the BG-on-both necessary
    condition); certifying which are *actual* destabilizing walls requires the
    ABCH/Maciocia subobject analysis.  ``numerical_wall`` is the exact primitive
    for a single subobject; widen ``degree_bound``/``rank_bound`` for breadth.
    """
    if v.r <= 0:
        raise ValueError("compute_walls expects a positive-rank v")
    d = surface.d
    s_min, s_max = Q(s_range[0]), Q(s_range[1])
    if s_min > s_max:
        raise ValueError("s_range must be (s_min, s_max) with s_min <= s_max")
    seen = set()
    walls: List[Wall] = []

    r_lo = 0 if include_torsion else 1
    r_hi = min(rank_bound, v.r)
    for rp in range(r_lo, r_hi + 1):
        # Bound the Mukai-degree W_rc = v.r*c' - rp*v.c by degree_bound -- the
        # natural finiteness parameter (the wall center is W_re/W_rc).  For each
        # admissible W_rc, recover the unique integer c' (when one exists).
        for W_int in range(-degree_bound, degree_bound + 1):
            if W_int == 0:
                continue  # equal slopes -> vertical wall (use numerical_wall directly)
            cp = (Fraction(W_int) + rp * v.c) / v.r
            if cp.denominator != 1:
                continue  # no integer c' realizes this W_rc
            if rp == 0 and cp <= 0:
                continue  # a torsion subobject needs effective degree c' > 0
            W_rc = Fraction(W_int)
            # e' placing the center s0 = (v.r e' - rp v.e)/W_rc inside [s_min, s_max]
            e_a = Fraction(s_min * W_rc + rp * v.e, v.r)
            e_b = Fraction(s_max * W_rc + rp * v.e, v.r)
            e_lo, e_hi = (e_a, e_b) if e_a <= e_b else (e_b, e_a)
            for ep in _half_integers(e_lo, e_hi):
                w = ChernChar(rp, cp, ep)
                if require_bg:
                    q = ChernChar(v.r - rp, v.c - cp, v.e - ep)
                    if not (_bg_ok(w, d) and _bg_ok(q, d)):
                        continue
                wall = numerical_wall(v, w, d)
                if isinstance(wall, VerticalWall):
                    continue
                if not wall.is_real:
                    continue
                if not (s_min <= wall.center <= s_max):
                    continue
                key = (wall.center, wall.radius_sq)
                if key in seen:
                    continue
                seen.add(key)
                walls.append(wall)

    walls.sort(key=lambda W: (-W.radius_sq, W.center))
    return walls


def walls_from_subobjects(
    v: ChernChar, subobjects: List[ChernChar], surface: Surface
) -> List[Wall]:
    """Exact walls of ``v`` for an explicit list of destabilizing classes.

    The honest way to reproduce a published wall diagram: supply the known
    destabilizers (e.g. ``O(-1) = ChernChar(1,-1,1/2)`` for P^2[2]) and get the
    exact semicircles, with no enumeration heuristic.  Non-real walls
    (``radius_sq < 0``) and vertical walls are included as returned by
    :func:`numerical_wall`; filter with ``[w for w in result if getattr(w,
    'is_real', False)]`` as needed.
    """
    return [numerical_wall(v, w, surface.d) for w in subobjects]


# --------------------------------------------------------------------------
# Actual (potential) walls -- the certified-necessary refinement
# --------------------------------------------------------------------------
def _bogomolov_integer(ch: ChernChar, d: int) -> Fraction:
    """``c^2/d - 2 r e``: >= 0 is the Bogomolov inequality (rank 0 -> always >= 0)."""
    return Fraction(ch.c * ch.c, d) - 2 * ch.r * ch.e


def _is_integral_chern(ch: ChernChar) -> bool:
    """True iff ``(r, c1=c, ch2=e)`` lies in the Chern-character lattice of P^2,
    i.e. ``c2 = c1^2/2 - ch2`` is an integer (so it is the class of an actual object)."""
    return (Fraction(ch.c * ch.c, 2) - ch.e).denominator == 1


def _heart_ordering_ok(v: ChernChar, w: ChernChar, s0: Fraction, d: int) -> bool:
    """Heart/phase condition at the wall: ``Im Z(w) > 0`` and ``Im Z(v-w) > 0``.

    Along the wall ``Z(w)`` is parallel to ``Z(v)``; for ``w`` to be a genuine
    sub-object of an object of class ``v`` in the tilted heart, both ``w`` and
    the quotient must have strictly positive imaginary central charge there.
    ``Im Z(x) = t (c_x - s_0 r_x d)``, so the sign is that of ``c_x - s_0 r_x d``.
    """
    im_w = w.c - s0 * w.r * d
    im_u = (v.c - w.c) - s0 * (v.r - w.r) * d
    return im_w > 0 and im_u > 0


@dataclass(frozen=True)
class ActualWall:
    """A wall satisfying the rigorous necessary conditions for being an actual wall.

    Same geometry as :class:`Wall`, plus the certified destabilizing class and
    the quotient.
    """

    center: Fraction
    radius_sq: Fraction
    subobject: ChernChar
    quotient: ChernChar
    v: ChernChar

    is_real = True

    @property
    def radius(self) -> float:
        return math.sqrt(float(self.radius_sq))

    @property
    def s_left(self) -> float:
        return float(self.center) - self.radius

    @property
    def s_right(self) -> float:
        return float(self.center) + self.radius

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"ActualWall(center={self.center}, R={self.radius:.4f}, "
            f"sub={self.subobject}, quot={self.quotient})"
        )


def actual_walls(
    v: ChernChar,
    surface: Surface,
    *,
    rank_bound: Optional[int] = None,
    degree_bound: int = 16,
    center_window: Number = 12,
    include_torsion: bool = True,
) -> List[ActualWall]:
    """Enumerate the walls of ``v`` that satisfy the necessary conditions for
    being **actual** destabilizing walls, returning each with its destabilizer.

    A numerical wall ``W(v, w)`` is kept iff ALL of:

    * **rank reduction** -- ``0 <= rank(w) <= rank(v)`` (the first destabilizing
      object has rank at most ``rank(v)``; Coskun-Huizenga survey, sec. 6);
    * **Bogomolov on both pieces** -- ``Delta(w) >= 0`` and ``Delta(v-w) >= 0``;
    * **real semicircle** -- ``radius^2 > 0`` (strict; a zero-radius "wall" is a
      point, not a wall);
    * **heart/phase ordering** -- ``Im Z(w) > 0`` and ``Im Z(v-w) > 0`` on the
      wall, i.e. ``w`` is a genuine sub-object in the tilted heart there.

    These conditions are NECESSARY for an actual wall and -- crucially -- finite
    (the dense numerical-wall set returned by :func:`compute_walls` is cut down
    to the destabilizers that can really occur).  For the Hilbert scheme
    ``P^2[n]`` and the coprime / small-rank cases covered by the ABCH /
    Coskun-Huizenga theorems this is exactly the set of actual walls; e.g. for
    ``P^2[2]`` (``v = (1,0,-2)``) it returns the single wall ``(-5/2, 3/2)``
    realized by ``O(-1)``.  For general ``v`` it is the certified-necessary
    candidate set; use :func:`actual_walls_complete` to check stability of the
    set under widening the search bounds.

    Walls are de-duplicated by ``(center, radius_sq)`` (collecting all
    destabilizers found for each) and sorted by radius, largest first (the
    largest is the Gieseker wall).
    """
    if v.r <= 0:
        raise ValueError("actual_walls expects a positive-rank v")
    d = surface.d
    if v.discriminant(d) < 0:
        raise ValueError("v violates Bogomolov (Delta < 0): not a semistable class")
    if not _is_integral_chern(v):
        raise ValueError(
            "v is not in the Chern-character lattice (c2 not integral): it is not "
            "the class of an actual object, so 'actual walls' are undefined"
        )
    mu_v = v.slope(d)
    S = Q(center_window)
    s_lo, s_hi = mu_v - S, mu_v + S
    r_cap = v.r if rank_bound is None else min(rank_bound, v.r)
    r_lo = 0 if include_torsion else 1

    best: dict = {}  # (center, radius_sq) -> ActualWall
    for rp in range(r_lo, r_cap + 1):
        for W_int in range(-degree_bound, degree_bound + 1):
            if W_int == 0:
                continue  # vertical wall
            cp = (Fraction(W_int) + rp * v.c) / v.r
            if cp.denominator != 1:
                continue
            if rp == 0 and cp <= 0:
                continue  # torsion sub needs effective degree > 0
            W_rc = Fraction(W_int)
            # e' placing the center in [s_lo, s_hi]
            e_a = Fraction(s_lo * W_rc + rp * v.e, v.r)
            e_b = Fraction(s_hi * W_rc + rp * v.e, v.r)
            e_lo, e_hi = (e_a, e_b) if e_a <= e_b else (e_b, e_a)
            for ep in _half_integers(e_lo, e_hi):
                w = ChernChar(rp, cp, ep)
                u = ChernChar(v.r - rp, v.c - cp, v.e - ep)
                if not _is_integral_chern(w):
                    continue  # w (hence u) must be the class of an actual object
                if _bogomolov_integer(w, d) < 0:
                    continue
                if _bogomolov_integer(u, d) < 0:
                    continue
                wall = numerical_wall(v, w, d)
                if isinstance(wall, VerticalWall) or wall.radius_sq <= 0:
                    continue
                if not (s_lo <= wall.center <= s_hi):
                    continue
                if not _heart_ordering_ok(v, w, wall.center, d):
                    continue
                key = (wall.center, wall.radius_sq)
                if key not in best:
                    best[key] = ActualWall(wall.center, wall.radius_sq, w, u, v)
    out = list(best.values())
    out.sort(key=lambda W: (-W.radius_sq, W.center))
    return out


def actual_walls_complete(
    v: ChernChar, surface: Surface, **kwargs
) -> Tuple[List[ActualWall], bool]:
    """Run :func:`actual_walls` and check the set is stable under doubled bounds.

    Returns ``(walls, complete)`` where ``complete`` is True iff doubling both
    ``degree_bound`` and ``center_window`` yields the same set of semicircles
    (an empirical completeness certificate for the chosen bounds).
    """
    deg = kwargs.pop("degree_bound", 16)
    win = kwargs.pop("center_window", 12)
    base = actual_walls(v, surface, degree_bound=deg, center_window=win, **kwargs)
    wide = actual_walls(v, surface, degree_bound=2 * deg, center_window=2 * Q(win), **kwargs)
    base_keys = {(w.center, w.radius_sq) for w in base}
    wide_keys = {(w.center, w.radius_sq) for w in wide}
    return wide, base_keys == wide_keys
