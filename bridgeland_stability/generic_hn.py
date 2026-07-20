"""E13-M3b / G18: the generic Harder-Narasimhan filtration on ``F_e`` (arXiv:1907.06739 Sec. 5).

This module ships the Kronecker-module datum E13-M3a deferred: for a character
``v = (r, c1, ch2)`` on an ample-polarized Hirzebruch surface it COMPUTES the
characters of the ``H_m``-Harder-Narasimhan factors of the general prioritary
sheaf, by the paper's finite inductive procedure.  Since (Sec. 1.6, verbatim)

    "there exists an H_m-semistable sheaf with Chern character v if and only if
     the generic H_m-Harder-Narasimhan filtration has length 1,"

this makes non-emptiness DECIDABLE for every strictly ample polarization on
``F_e`` -- including the band ``emptiness_bound <= Delta <= delta_H`` that
E13-M3a honestly reported UNKNOWN, and the boundary ``Delta == delta_H``.

The mathematics (all citations are \\label names in arXiv:1907.06739)
---------------------------------------------------------------------
Fix ``m in Q_{>0}`` and the polarization ray ``H_m = E + (m+e)F``.  A package
surface with strictly ample ``H = a f + b s`` (Nakai: ``b > 0``, ``a > e b``)
carries ``m = a/b - e`` and ``H_m = H / b`` in ``NS otimes Q``.  The standing
hypothesis of Sec. 5 ("sec-genHN") is that the stack of ``F``- and
``H_{ceil(m)}``-prioritary sheaves of character ``v`` is nonempty; then the
general such sheaf has an ``H_m``-HN filtration with quotient characters
``gr_1, ..., gr_l`` depending only on ``(m, v)``:

* ``lem-HNclose``:   ``0 <= (nu_1 - nu_l).H_m <= 1`` and ``|(nu_i - nu).H_m| < 1``.
* ``lem-HNorthogonal``: ``chi(gr_i, gr_j) = 0`` for ``i < j``.
* ``thm-HNcriterion`` (the load-bearing iff): positive-rank characters
  ``w_1, ..., w_k`` are exactly the ``gr_i`` iff (1) they sum to ``v``,
  (2) their reduced ``H_m``-Hilbert polynomials strictly decrease,
  (3) ``mu_{H_m}(w_1) - mu_{H_m}(w_k) <= 1``, (4) ``chi(w_i, w_j) = 0`` for
  ``i < j``, and (5) every ``M_{H_m}(w_i)`` is nonempty.
* ``lem-slopeQuad``: each ``nu_i`` lies in the parallelogram
  ``|(nu_i - nu).H_m| < 1``, ``|(nu_i - nu).F| < max{1, 2/(e+2m)}`` -- the
  bound that makes the search finite.  (The filtration length is ``<= 4``.)
* ``cor-algorithm`` (the inductive form implemented here): ``w_1 = gr_1`` iff
  (1) there are ``H_{ceil(m)}``-prioritary sheaves of character ``u = v - w_1``
  -- whose generic HN factors ``w_2, ..., w_k`` are known by induction on the
  rank -- (2) ``q_1 > q_2``, (3) ``mu_{H_m}(w_1) - mu_{H_m}(w_k) <= 1``,
  (4) ``chi(w_1, w_i) = 0`` for ``i >= 2``, and (5) ``M_{H_m}(w_1)`` is
  nonempty.  If NO such ``w_1`` exists, the filtration has length one and
  ``M_{H_m}(v)`` is nonempty.

Supporting facts used by the implementation:

* ``prop-ssPrior``: a ``mu_{H_m}``-semistable sheaf is ``H_n``-prioritary for
  every integer ``n < m + 2`` -- in particular ``n = ceil(m)``.  So if the
  ``F``-and-``H_{ceil(m)}``-prioritary stack is EMPTY then ``M_{H_m}(v)`` is
  empty (this module returns ``None`` and the verdict layer reports a PROVEN
  emptiness).
* ``cor-prioritaryDelta`` (= the shipped E13-M2 ``prioritary_nonempty``): for
  ``Delta >= 0`` the ``F``-and-``H_n``-prioritary stack is nonempty iff
  ``Delta >= delta^p_n(nu)``.
* For ``Delta < 0`` there are NO ``F``-and-``H_n``-prioritary sheaves once
  ``n >= 1``: by ``thm-prioritaryNecessary`` prioritary existence forces
  ``chi(v(-L_0 - H_n)) <= 0``, and the Riemann-Roch expansion in
  ``cor-equivalentInequality`` evaluates that inequality as
  ``n <= Delta/((ceil(eps)-eps)(eps-floor(eps))) - e/2 + 1 - (ceil(psi)-psi)``
  for ``eps`` non-integral -- whose right side is ``<= 4*Delta + 1 < 1`` when
  ``Delta < 0`` -- and as ``-r*Delta <= 0``, i.e. ``Delta >= 0``, when ``eps``
  is integral (``rem-epInteger``: ``ceil(eps) - eps = 0`` collapses the
  expansion to ``chi/r = -Delta``).  Either way ``Delta < 0`` refuses.
* Rank-one base case: a rank-one torsion-free sheaf is ``I_Z(L)`` with
  ``c_2 = length(Z) >= 0``; it is stable for every polarization.  So
  ``M_{H_m}(1, c1, ch2)`` is nonempty iff the character is integral with
  ``c_2 >= 0``, and there is no rank-one character with ``c_2 < 0``.

Reduced ``H_m``-Hilbert polynomials, exactly
--------------------------------------------
``q_w(t) = P(nu_w + t H_m) - Delta_w`` is a degree-2 polynomial in ``t`` whose
leading coefficient ``H_m^2/2`` is common to every character, so the total
order ``q_i > q_j`` is the lexicographic comparison of
``( nu.H_m , P(nu) - Delta )`` -- both exact ``Fraction``s here (the
``-K.H_m/2`` contribution to the linear coefficient is also common and
cancels).  No polynomial is ever materialized.

Honest scope
------------
* Any ``F_e`` (``e >= 0``) carrying a STRICTLY AMPLE ``H`` -- authenticated by
  :func:`~bridgeland_stability.dlp_hirzebruch.hirzebruch_index` (E13 R1) and
  :func:`~bridgeland_stability.dlp_hirzebruch.is_ample`.  The nef-and-big
  factory polarization is refused.  (The E13-M3a verdict layer exposes only
  the del Pezzo ``e in {0, 1}``; unlocking ``e >= 2`` there is E13-M3c, where
  this direct computation will be cross-checked against the E13-M1 reduction.)
* NOT ``P^2`` -- the Drezet-Le Potier verdict there is already total.
* Exact ``Fraction`` arithmetic throughout; stdlib-only at import time.
* The search enumerates CLOSED versions of the strict paper windows (a
  superset of candidates -- extra candidates simply fail the iff conditions),
  and asserts the ``cor-algorithm`` uniqueness: every passing ``w_1`` must
  yield the same filtration, or the transcription is broken (a tripwire, like
  the Markov assertion of the mutation oracle).

References
----------
* Coskun-Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
  arXiv:1907.06739 -- Sec. 1.6, Sec. 5 (lem-HNclose, lem-HNorthogonal,
  thm-HNcriterion, lem-slopeQuad, lem-discBound, cor-algorithm), Sec. 4
  (thm-prioritaryNecessary, def-L0, cor-equivalentInequality, rem-epInteger,
  prop-ssPrior, cor-prioritaryDelta).
"""

from __future__ import annotations

from fractions import Fraction
from math import ceil, floor
from typing import Dict, Optional, Sequence, Tuple, Union

from .dlp_hirzebruch import (
    discriminant,
    emptiness_bound,
    hirzebruch_index,
    is_ample,
    is_semiexceptional,
)
from .exceptional_surface import SurfaceBundle
from .nonemptiness_rational import validate_character
from .prioritary import delta_prioritary
from .varieties import Surface

__all__ = [
    "generic_hn_factors",
    "semistable_exists_hn",
    "Character",
]

Number = Union[int, Fraction]
#: A factor character ``(r, c1, ch2)`` -- ``c1`` an integer NS-vector in the
#: package ``(f, s)`` basis, ``ch2`` an exact ``Fraction``.
Character = Tuple[int, Tuple[int, int], Fraction]

_CH = "arXiv:1907.06739"


# --------------------------------------------------------------------------
# Exact per-surface context.
# --------------------------------------------------------------------------
#: E15-M1b: cross-call persistent caches, keyed by (e, H_a, H_b) -- see _Ctx.
_PERSISTENT_CACHES: Dict[tuple, tuple] = {}

#: E15-M1b: optional progress hook for long computations (the 67-CPU-hour
#: post-mortem's telemetry requirement).  Set via :func:`set_progress`; called
#: as ``fn(event, depth, payload)`` with events ``"decide"`` (a character's
#: filtration was computed, payload = (character, length-or-None)) and
#: ``"gr1"`` (the TOP-LEVEL first factor was pinned, payload = the factor) --
#: the rigid-factor test (CORRECTIONS Sec. 21) needs only gr_1, so a caller
#: can act the moment it fires.  The hook must be cheap; rate-limit inside it.
_PROGRESS = None


def set_progress(fn) -> None:
    """Install (or clear, with ``None``) the module progress hook."""
    global _PROGRESS
    _PROGRESS = fn


class _Ctx:
    """Precomputed exact data for one (surface, H) pair: ``e``, ``m``, ``H_m = H/b``."""

    def __init__(self, surface: Surface):
        self.surface = surface
        self.e = hirzebruch_index(surface)          # authenticates the family (R1)
        if not is_ample(surface):
            raise ValueError(
                f"{surface.name}: the generic-HN algorithm needs a strictly ample H "
                "(Nakai b>0, a>e*b); the nef-and-big factory polarization is refused")
        a, b = surface.H
        self.b = b
        self.m = Fraction(a, b) - self.e             # H = b * H_m, m in Q_{>0}
        self.n = ceil(self.m)                        # the prioritary index ceil(m) >= 1
        self.lat = surface.lattice
        # E15-M1b: the caches PERSIST across calls, keyed by the exact (e, H)
        # pair -- the results are pure functions of (surface family, H,
        # character), and the Sec. 5 recursion revisits the same sub-characters
        # across calls (the Sec. 19 memoized-induction lesson).  Unbounded by
        # design; a long-running sweep that must bound memory can clear
        # _PERSISTENT_CACHES explicitly.
        # PARANOID_UNIQUENESS must see fresh caches: the tripwire's value is in
        # RECOMPUTING with the full-sweep assertions on, which a warm persistent
        # cache would silently skip.
        if PARANOID_UNIQUENESS:
            store = ({}, {}, {}, {})
        else:
            key = (self.e, surface.H[0], surface.H[1])
            store = _PERSISTENT_CACHES.setdefault(key, ({}, {}, {}, {}))
        self.cache: Dict[Character, Optional[Tuple[Character, ...]]] = store[0]
        self._qcache: Dict[Character, Tuple[Fraction, Fraction]] = store[1]
        self._mucache: Dict[Character, Fraction] = store[2]
        self._priocache: Dict[Character, bool] = store[3]

    # -- exact numerical invariants of a character -------------------------
    def nu(self, w: Character) -> Tuple[Fraction, Fraction]:
        r, c1, _ = w
        return (Fraction(c1[0], r), Fraction(c1[1], r))

    def mu_Hm(self, w: Character) -> Fraction:
        """``mu_{H_m}(w) = <nu_w, H_m> = <nu_w, H>/b`` (exact, memoized)."""
        mu = self._mucache.get(w)
        if mu is None:
            # <(x,y)/r, (a,b)> / b with Gram [[0,1],[1,-e]]:
            #   <(x,y),(a,b)> = x*b + y*a - e*y*b.
            r, (x, y), _ = w
            a, b = self.surface.H
            mu = Fraction(x * b + y * a - self.e * y * b, r * b)
            self._mucache[w] = mu
        return mu

    def delta(self, w: Character) -> Fraction:
        r, (x, y), ch2 = w
        # <c1,c1> = 2xy - e y^2 (Gram [[0,1],[1,-e]]); Delta = <c1,c1>/(2 r^2) - ch2/r.
        return Fraction(2 * x * y - self.e * y * y, 2 * r * r) - Fraction(ch2, 1) / r

    def q_key(self, w: Character) -> Tuple[Fraction, Fraction]:
        """The reduced ``H_m``-Hilbert polynomial as a lexicographic sort key
        ``(nu.H_m, P(nu) - Delta)`` -- see the module docstring.  Memoized.

        ``P(nu) - Delta`` is inlined in integer arithmetic:
        ``P(nu) = chi(O) + (1/2)(<nu,nu> - <nu,K>)`` with ``chi(O) = 1``,
        ``<nu,nu> = (2xy - e y^2)/r^2`` and ``<nu,K> = (-2x + (e-2)y)/r``
        (``K = (-(e+2), -2)``, Gram ``[[0,1],[1,-e]]``), and
        ``Delta = <nu,nu>/2 - ch2/r`` -- so ``P(nu) - Delta = 1 - <nu,K>/2 +
        ch2/r``.  Cross-pinned against :func:`~bridgeland_stability.
        dlp_hirzebruch.hilbert_P` in the tests.
        """
        q = self._qcache.get(w)
        if q is None:
            r, (x, y), ch2 = w
            nuK = Fraction(-2 * x + (self.e - 2) * y, r)
            q = (self.mu_Hm(w), 1 - nuK / 2 + Fraction(ch2, 1) / r)
            self._qcache[w] = q
        return q

    def chi(self, w1: Character, w2: Character) -> int:
        """Integer Riemann-Roch Euler pairing, inlined for the hot path.

        With ``K = (-(e+2), -2)`` in the (f, s) basis, ``chi(O) = 1``, and the
        Gram ``[[0,1],[1,-e]]``, the RR form
        ``chi(A,B) = rA rB chi(O) - <rA c1B - rB c1A, K>/2
        + (rA ch2B + rB ch2A - <c1A, c1B>)`` expands (in doubled units, all
        integers) to the closed form below.  Pinned against the package
        :func:`bridgeland_stability.exceptional_surface.chi` in the tests.

        E15-M1b: the whole computation runs over ``int`` (the ch2 components,
        half-integers with denominator 1 or 2, enter as twice-ch2 integers) --
        the 67-CPU-hour post-mortem (CORRECTIONS Sec. 21) sampled 100% of the
        runtime in Fraction arithmetic under this pairing.
        """
        rA, c1A, chA = w1
        rB, c1B, chB = w2
        two_chi = _two_chi(self.e, rA, c1A, chA, rB, c1B, chB)
        if two_chi & 1:
            raise ValueError(
                f"chi({w1!r}, {w2!r}) = {two_chi}/2 is not integral")
        return two_chi >> 1

    def c2(self, w: Character) -> Fraction:
        _, (x, y), ch2 = w
        return Fraction(2 * x * y - self.e * y * y, 2) - ch2

    # -- the prioritary standing hypothesis (condition (1)) ----------------
    def prioritary_nonempty(self, w: Character) -> bool:
        """``P_{F, H_n}(w)`` nonempty, ``n = ceil(m) >= 1`` -- the classification.

        * rank 1: always nonempty ... provided the character is that of a
          torsion-free sheaf, i.e. ``c2 >= 0`` (thm-prioritaryChi's proof:
          "If r = 1, then the inequality automatically holds and the stack is
          nonempty").
        * ``Delta < 0``: EMPTY -- PROVEN via thm-prioritaryNecessary +
          cor-equivalentInequality / rem-epInteger (module docstring).
        * ``Delta >= 0``: cor-prioritaryDelta, i.e. ``Delta >= delta^p_n(nu)``
          (the shipped E13-M2 sharp bound).  Memoized.
        """
        got = self._priocache.get(w)
        if got is not None:
            return got
        r = w[0]
        if r == 1:
            ok = self.c2(w) >= 0
        else:
            d = self.delta(w)
            if d < 0:
                ok = False
            else:
                ok = d >= delta_prioritary(self.nu(w), self.n, self.surface)
        self._priocache[w] = ok
        return ok


def _twice(ch2) -> int:
    """``2 * ch2`` as a pure ``int`` for the half-integer ``ch2`` of an integral
    character (denominator 1 or 2) -- no Fraction arithmetic on the hot path."""
    if isinstance(ch2, int):
        return 2 * ch2
    d = ch2.denominator
    if d == 1:
        return 2 * ch2.numerator
    if d == 2:
        return ch2.numerator
    raise ValueError(f"ch2 = {ch2!r} is not a half-integer")


def _two_chi(e, rA, c1A, chA, rB, c1B, chB):
    """``2 * chi(A, B)`` on ``F_e`` in doubled integer-friendly units.

    RR with ``chi(O) = 1``, ``K = (-(e+2), -2)``, Gram ``[[0,1],[1,-e]]``:
    ``chi(A,B) = rA rB - <rA c1B - rB c1A, K>/2 + rA ch2B + rB ch2A - <c1A, c1B>``
    where ``<D, K> = -2 D0 + (e-2) D1`` and
    ``<c1A, c1B> = xA yB + yA xB - e yA yB``.  All-int (E15-M1b): the ch2
    terms enter as ``rA * (2 ch2B) + rB * (2 ch2A)`` via :func:`_twice`.
    """
    xA, yA = c1A
    xB, yB = c1B
    d0 = rA * xB - rB * xA
    d1 = rA * yB - rB * yA
    return (2 * rA * rB
            - (-2 * d0 + (e - 2) * d1)
            + rA * _twice(chB) + rB * _twice(chA)
            - 2 * (xA * yB + yA * xB - e * yA * yB))


# --------------------------------------------------------------------------
# The inductive algorithm (cor-algorithm), memoized, rank-decreasing.
# --------------------------------------------------------------------------
def _decide(ctx: _Ctx, v: Character) -> Optional[Tuple[Character, ...]]:
    """The generic ``H_m``-HN factor characters of ``v``, or ``None`` when the
    ``F``-and-``H_{ceil(m)}``-prioritary stack of character ``v`` is empty
    (in which case ``M_{H_m}(v)`` is empty by prop-ssPrior)."""
    if v in ctx.cache:
        return ctx.cache[v]
    r = v[0]
    if r == 1:
        result: Optional[Tuple[Character, ...]] = (v,) if ctx.c2(v) >= 0 else None
    elif not ctx.prioritary_nonempty(v):
        result = None
    else:
        found = _search_gr1(ctx, v)
        result = found if found is not None else (v,)
    ctx.cache[v] = result
    if _PROGRESS is not None:
        _PROGRESS("decide", len(ctx.cache), (v, None if result is None else len(result)))
    return result


#: Opt-in paranoia: sweep ALL candidates and assert the cor-algorithm uniqueness
#: (every passing ``w_1`` is gr_1, so all filtrations must coincide).  Default
#: off -- the search returns on the first winner, which the iff makes correct;
#: the dedicated test flips this on over small-rank grids as a tripwire.
PARANOID_UNIQUENESS = False


def _delta1_cap(ctx: _Ctx, Bf: Fraction) -> Fraction:
    """A PROVEN upper bound for ``Delta_1`` when the filtration has length >= 2.

    lem-HNorthogonal gives ``chi(gr_1, gr_j) = 0`` for any ``j >= 2``, so by
    Riemann-Roch ``Delta_1 = P(nu_j - nu_1) - Delta_j <= P(nu_j - nu_1)``
    (Bogomolov ``Delta_j >= 0`` for the semistable factor).  The difference
    ``d = nu_j - nu_1`` lies in the doubled lem-HNclose/lem-slopeQuad window
    ``|d.H_m| <= 2``, ``|d.F| <= 2 Bf``.  Writing ``d = eps E + phi F``:
    ``d.F = eps`` and ``d.H_m = eps m + phi``, so ``|eps| <= 2 Bf`` and
    ``|phi| <= 2 + 2 Bf m``; then

        P(d) = 1 + (1/2)(<d,d> - d.K),   <d,d> = 2 eps phi - e eps^2 <= 2|eps||phi|,
        -d.K = eps(2 - e) + 2 phi <= |eps|(2 + e) + 2|phi|,

    which yields the exact coarse cap returned here.  Over-estimation is
    harmless (a wider window); UNDER-estimation is impossible by the above.
    """
    E = 2 * Bf
    Phi = 2 + 2 * Bf * ctx.m
    return 1 + Fraction(1, 2) * (2 * E * Phi + E * (2 + ctx.e) + 2 * Phi)


def _search_gr1(ctx: _Ctx, v: Character) -> Optional[Tuple[Character, ...]]:
    """Search for ``w_1 = gr_1`` per cor-algorithm; return the full factor list
    ``(w_1, w_2, ..., w_k)`` or ``None`` if no valid ``w_1`` exists (length 1).

    The enumeration windows are CLOSED supersets of the strict lem-HNclose /
    lem-slopeQuad bounds; every candidate is then vetted against the exact
    cor-algorithm conditions, so over-enumeration cannot create a factor.
    Three PROVEN prunes bound the discriminant sweep:

    * the TOTAL orthogonality ``chi(w_1, u) = sum_i chi(w_1, w_i) = 0`` --
      forced by condition (4) and chi-bilinearity -- is LINEAR in ``ch2_1``
      with coefficient ``2(r_u - r_1)``; when ``r_1 != r_u`` it pins ``ch2_1``
      to at most ONE value per ``(r_1, c_1)``, no scan at all.  In the balanced
      case ``r_1 = r_u`` it is ``ch2_1``-independent: either identically
      violated (skip the cell) or identically satisfied (scan the window);
    * ``Delta_1 <= _delta1_cap`` (orthogonality + Bogomolov, see there);
    * lem-discBound: ``Delta_1`` is the SMALLEST discriminant of an
      ``H_m``-semistable sheaf of slope ``nu_1`` (in both of its cases), so at
      a fixed ``(r_1, c_1)`` the ascending-``Delta_1`` scan may stop after the
      first candidate whose moduli space is nonempty -- any larger ``Delta_1``
      at this slope is dominated by that semistable sheaf and cannot be
      minimal.
    """
    r, c1, ch2 = v
    lat, H, b = ctx.lat, ctx.surface.H, ctx.b
    nu_v = ctx.nu(v)
    mu_v = ctx.mu_Hm(v)
    q_v = ctx.q_key(v)
    # lem-slopeQuad F-window half-width: max{1, 2/(e+2m)} (exact Fraction).
    Bf = max(Fraction(1), Fraction(2, 1) / (ctx.e + 2 * ctx.m))
    cap = _delta1_cap(ctx, Bf)
    a_H, b_H = H
    winners = []

    for r1 in range(1, r):
        ru = r - r1
        # (nu_1)_s = y/r1 with |y/r1 - (nu_v)_s| <= Bf   [nu.F = the s-coefficient]
        y_lo = _ceil_frac(r1 * (nu_v[1] - Bf))
        y_hi = _floor_frac(r1 * (nu_v[1] + Bf))
        for y in range(y_lo, y_hi + 1):
            # 0 <= <(x,y), H>/(r1 b) - mu_v <= 1, with <(x,y),(a,b)> = x*b + y*a - e*y*b
            # (Gram [[0,1],[1,-e]]).  b > 0, so solve for x:
            const = y * a_H - ctx.e * y * b_H
            lo = r1 * b * mu_v - const                # x*b >= lo
            hi = r1 * b * (mu_v + 1) - const          # x*b <= hi
            x_lo = _ceil_frac(Fraction(lo, b))
            x_hi = _floor_frac(Fraction(hi, b))
            for x in range(x_lo, x_hi + 1):
                c11 = (x, y)
                c1u = (c1[0] - x, c1[1] - y)
                # ch2 window: 0 <= Delta_1 <= cap  AND  Delta(u) >= 0 (else
                # condition (1) provably fails -- module docstring).  ch2_1 runs
                # over the integral-c2 lattice ch2_1 = <c11,c11>/2 - c2_1.
                s11 = lat.self_pairing(c11)
                s1u = lat.self_pairing(c1u)
                ch2_hi = Fraction(s11, 2 * r1)                    # Delta_1 >= 0
                ch2_lo = max(ch2 - Fraction(s1u, 2 * ru),         # Delta_u >= 0
                             Fraction(s11, 2 * r1) - r1 * cap)    # Delta_1 <= cap
                c2_lo = _ceil_frac(Fraction(s11, 2) - ch2_hi)     # c2_1 = <c11,c11>/2 - ch2_1
                c2_hi = _floor_frac(Fraction(s11, 2) - ch2_lo)
                if c2_hi < c2_lo:
                    continue
                # Total orthogonality chi(w_1, u) = 0 (NECESSARY: condition (4)
                # summed over u's factors; chi is bilinear).  In doubled units
                # 2*chi(w_1, u) = A + 2*(ru - r1)*ch2_1 with A independent of
                # ch2_1 -- computed here at ch2_1 = 0, ch2_u = ch2.
                A = _two_chi(ctx.e, r1, c11, Fraction(0), ru, c1u, ch2)
                if r1 != ru:
                    # unique root; must land on the integral-c2 lattice + window
                    ch2_1_star = Fraction(A, 2 * (r1 - ru))
                    c2_star = Fraction(s11, 2) - ch2_1_star
                    if c2_star.denominator != 1:
                        continue
                    c2_range = ((int(c2_star),)
                                if c2_lo <= int(c2_star) <= c2_hi else ())
                elif A != 0:
                    continue                     # balanced case, identically violated
                else:
                    c2_range = range(c2_lo, c2_hi + 1)   # balanced, identically zero
                for c2_1 in c2_range:                    # Delta_1 ASCENDING in c2_1
                    ch2_1 = Fraction(s11, 2) - c2_1
                    w1: Character = (r1, c11, ch2_1)
                    u: Character = (ru, c1u, ch2 - ch2_1)
                    # PROVEN prune: q_1 strictly exceeds q_v when l >= 2
                    # (P additivity: q_v is the rank-weighted mean of the q_i).
                    if not ctx.q_key(w1) > q_v:
                        continue
                    # cheap first: (1)'s prioritary formula for u
                    if not ctx.prioritary_nonempty(u):
                        continue
                    # (5) M_{H_m}(w_1) nonempty  [rank r1 < r: induction, with
                    # the E15-M1c PROVEN envelope shortcuts before recursing]
                    if not _ss_exists(ctx, w1):
                        continue
                    passed = _check_tail(ctx, w1, u)
                    if passed is not None:
                        if not PARANOID_UNIQUENESS:
                            return passed
                        winners.append(passed)
                    # lem-discBound: w_1 was semistable, so any LARGER Delta_1 at
                    # this (r_1, c_1) is not minimal-at-slope: stop the scan.
                    break
    if not winners:
        return None
    # cor-algorithm uniqueness tripwire (PARANOID_UNIQUENESS): every passing
    # w_1 IS gr_1, so all candidate filtrations must coincide.
    first = winners[0]
    for other in winners[1:]:
        if other != first:
            raise AssertionError(
                f"cor-algorithm uniqueness violated for v={v!r}: {first!r} vs {other!r} "
                "(two distinct maximal destabilizing characters passed the iff)")
    return first


_MISS = object()


def _ss_exists(ctx: _Ctx, w: Character) -> bool:
    """The condition-(5) boolean: is ``M_{H_m}(w)`` nonempty (Sec. 1.6: the
    generic filtration has length one)?

    E15-M1c: the ``w_1``-validation of ``cor-algorithm`` consumes ONLY this
    boolean, yet the plain recursion computes the full sub-filtration -- the
    measured super-exponential wall of CORRECTIONS Sec. 21.  Two PROVEN
    envelope-side deciders answer many candidates in O(1) after per-surface
    cache warmup:

    * ``is_semiexceptional`` -- ``V^{+k}`` of a mu_H-stable exceptional ``V``
      is Gieseker-semistable (the E12-M3 disjunct), so ``M`` is NONEMPTY and
      the generic filtration is length one; caching ``(w,)`` is exactly what
      ``_decide`` would return, so the cache stays truthful;
    * ``discriminant < emptiness_bound`` -- a certified lower bound on the
      discriminant of EVERY Gieseker-semistable sheaf of that slope (only the
      theorem-backed branches; any ample ``F_e``), so ``M`` is EMPTY.  The
      factor list is NOT computed on this path, so nothing is cached.

    Everything else falls through to the full recursion, bit-identically.
    """
    cached = ctx.cache.get(w, _MISS)
    if cached is not _MISS:
        return cached is not None and len(cached) == 1
    r, c1, ch2 = w
    if r == 1:
        return ctx.c2(w) >= 0
    xi = SurfaceBundle(r, c1, ch2)
    if is_semiexceptional(xi, ctx.surface):
        ctx.cache[w] = (w,)                       # truthful: length-one filtration
        if _PROGRESS is not None:
            _PROGRESS("shortcut", len(ctx.cache), (w, "semiexceptional"))
        return True
    if discriminant(xi, ctx.surface) < emptiness_bound(xi, ctx.surface):
        if _PROGRESS is not None:
            _PROGRESS("shortcut", len(ctx.cache), (w, "emptiness_bound"))
        return False                              # factors unknown: never cached
    fw = _decide(ctx, w)
    return fw is not None and len(fw) == 1


def _check_tail(ctx: _Ctx, w1: Character, u: Character) -> Optional[Tuple[Character, ...]]:
    """Conditions (1)-(4) of cor-algorithm against u's recursive factors."""
    fu = _decide(ctx, u)
    if fu is None:
        return None
    # (2) q_1 > q_2
    if not ctx.q_key(w1) > ctx.q_key(fu[0]):
        return None
    # (3) mu(w_1) - mu(w_k) <= 1
    if ctx.mu_Hm(w1) - ctx.mu_Hm(fu[-1]) > 1:
        return None
    # (4) chi(w_1, w_i) = 0 for every factor of u
    if any(ctx.chi(w1, wi) != 0 for wi in fu):
        return None
    return (w1,) + fu


def _ceil_frac(x: Fraction) -> int:
    x = Fraction(x)
    return -((-x.numerator) // x.denominator)


def _floor_frac(x: Fraction) -> int:
    x = Fraction(x)
    return x.numerator // x.denominator


# --------------------------------------------------------------------------
# Public API.
# --------------------------------------------------------------------------
def _as_character(r: int, c1: Sequence[Number], ch2: Number) -> Character:
    return (int(r), (int(c1[0]), int(c1[1])), Fraction(ch2))


def generic_hn_factors(
    r: int, c1: Sequence[Number], ch2: Number, surface: Surface
) -> Optional[Tuple[Character, ...]]:
    """The characters of the generic ``H_m``-Harder-Narasimhan factors of ``v``.

    Returns the ordered tuple ``(gr_1, ..., gr_l)`` (each a ``(r, c1, ch2)``
    triple, ``1 <= l <= 4``), or ``None`` when the ``F``-and-``H_{ceil(m)}``-
    prioritary stack of character ``v`` is empty -- in which case
    ``M_{H_m}(v)`` is empty too (prop-ssPrior).  ``M_{H_m}(v)`` is nonempty
    iff the result is a length-one tuple (Sec. 1.6).

    Raises ``ValueError`` for a non-integral character (the K-theory machinery
    is stated for ``v in K(F_e)``), a non-ample polarization, or (via
    ``hirzebruch_index``) a non-Hirzebruch surface.
    """
    ctx = _Ctx(surface)
    rf = Fraction(r)
    if rf.denominator != 1 or rf < 1:
        raise ValueError(f"rank must be a positive integer; got {r!r}")
    if any(Fraction(x).denominator != 1 for x in c1):
        raise ValueError(f"c1 must be an integral NS-vector; got {tuple(c1)!r}")
    if not validate_character(int(rf), c1, ch2, surface):
        raise ValueError(
            f"(r, c1, ch2) = ({r}, {tuple(c1)}, {ch2}) is not an integral Chern "
            "character (c2 not in Z): out of the domain of the Sec. 5 algorithm")
    return _decide(ctx, _as_character(int(rf), c1, ch2))


def semistable_exists_hn(
    r: int, c1: Sequence[Number], ch2: Number, surface: Surface
) -> bool:
    """``True`` iff ``M_{H_m}(r, c1, ch2)`` is nonempty, by the Sec. 5 algorithm.

    An invalid (non-integral) character returns ``False`` -- it is not the
    Chern character of any sheaf, matching the package's invalid-character
    PROVEN-empty convention.
    """
    rf = Fraction(r)
    if rf.denominator != 1 or rf < 1:
        return False
    if any(Fraction(x).denominator != 1 for x in c1):
        return False
    if not validate_character(int(rf), c1, ch2, surface):
        return False
    factors = generic_hn_factors(int(rf), c1, ch2, surface)
    return factors is not None and len(factors) == 1
