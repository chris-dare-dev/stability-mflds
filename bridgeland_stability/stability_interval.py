"""E14-M3 / G18: stability intervals ``I_V`` of exceptional bundles on ``F_e``.

For a bundle ``V`` on ``F_e``, ``I_V = {m > 0 : V is mu_{H_m}-stable}`` (paper
Sec. 6.3).  An exceptional bundle is determined by its character
(``prop-excPrior`` (2)), and by ``cor-DLPexcdelPezzo`` / ``cor-DLPExceptional``

    I_V  =  {m > 0 : Delta(V) >= DLP_{H_m}^{<r}(nu)}
         =  {m > 0 : a mu_{H_m}-stable exceptional bundle of character V exists},

an OPEN interval containing the anticanonical index ``m_0 = 1 - e/2`` on a del
Pezzo (every exceptional bundle on ``F_0`` / ``F_1`` is ``mu_{-K}``-stable --
Gorodentsev).  This module computes the interval EXACTLY by the paper's rank
induction (``thm-stabilityInterval``):

    S_V = { m_{V,W} : W exceptional, r(W) < r(V), chi(W,V) > 0, m_{V,W} in I_W }

where ``m_{V,W}`` is the positive solution of ``mu_{H_m}(V) = mu_{H_m}(W)``
(``m = -b/a`` for the paper slope difference ``nu(V) - nu(W) = aE + bF``); then
``I_V`` is the connected component of ``R_{>0} minus S_V`` containing ``m_0``.
The membership test ``m in I_W`` is exactly
:func:`~bridgeland_stability.dlp_hirzebruch.is_stable_exceptional` at ``H_m``
(the same set by ``cor-DLPexcdelPezzo`` + uniqueness), so the induction
terminates through the shipped ``DLP^{<r}`` machinery.

Effectivity is ``rem-stabilityIntervalCompute``: a slope difference computing an
element of ``S_V`` must lie in the region ``R_3 . R_1`` (two unit-width strips;
for ``e = 1`` the horizontal strip is cut to the triangle ``(0,-1), (0,0),
(2,0)`` by ``P > 0``), and far-out LINE BUNDLES already compute elements
``M_1 > m_0`` (and ``M_0 < m_0`` when ``e = 0``; ``M_0 = 0`` when ``e = 1``),
which bounds the slope window ``m in (M_0, M_1)`` and hence the candidate set to
a finite lattice search.

For ``e >= 2`` the interval transports down the E13-M1 reduction: by
``cor-highermus``, ``mu_{A_m}``-stable existence corresponds along ``pi`` with
the H-index shifted by one, so ``I_v = {t > 0 : t + 1 in I_{pi(v)}}`` -- the
paper's Sec. 11 corollary ``I_V = (0, m_1 - 1)`` (the left endpoint clamps to
``0``).  This is how the paper decides that characters like
``(3, 1/3 E + F, 4/9)`` on ``F_4`` and the conjecture's first candidate
``(107, 25/107 E + 76/107 F)`` on ``F_3`` admit NO slope-stable exceptional
bundle for any polarization: the transported interval is empty.

Exact ``Fraction`` arithmetic; stdlib-only at import.  Primary source:
arXiv:1907.06739 Sec. 6.3 (``prop-interval``, ``thm-stabilityInterval``,
``rem-stabilityIntervalCompute``, ``rem-stabilityIntervalQuotient``,
``ex-stabilityIntervals`` Tables 1-2), Sec. 11 (``cor-highermus`` and the
exceptional-bundle corollaries).  Record: ``docs/CORRECTIONS.md`` Sec. 19.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Optional, Sequence, Tuple

from .chern import Number
from .varieties import Surface, P1xP1, hirzebruch
from .exceptional_surface import SurfaceBundle, chi
from .dlp_hirzebruch import (
    hirzebruch_index,
    is_potentially_exceptional,
    is_stable_exceptional,
)
from .delta_sharp import surface_with_index
from .reduction import pi_c1
from .rigor import Certificate, Rigor

__all__ = ["StabilityInterval", "stability_interval"]

_HALF = Fraction(1, 2)
_PROBE_CAP = 200

#: Session cache of computed intervals, keyed by (e, r, c1) -- the interval of a
#: character is a pure function of the character and the surface FAMILY, and the
#: rank induction revisits low-rank characters constantly.
_INTERVAL_MEMO: dict = {}

_CITATIONS = (
    "arXiv:1907.06739 thm-stabilityInterval, rem-stabilityIntervalCompute, "
    "cor-DLPexcdelPezzo, cor-DLPExceptional",
    "arXiv:1907.06739 Sec. 11 (cor-highermus; the exceptional-bundle transport)",
    "Gorodentsev (mu_{-K}-stability of exceptional bundles on del Pezzo surfaces)",
)


def _Q(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


@dataclass(frozen=True)
class StabilityInterval:
    """The open interval ``I_v = {m > 0 : a mu_{H_m}-stable exceptional bundle of
    character v exists}`` -- for the (unique) exceptional bundle ``V`` of that
    character, its stability interval.

    ``lo`` / ``hi`` are the open endpoints (``lo = 0`` and/or ``hi = None`` for
    an interval reaching ``0`` / ``infinity``); ``empty`` marks ``I_v`` empty
    (possible only for ``e >= 2``, where the transported window may close --
    on a del Pezzo an existing exceptional bundle is ``mu_{-K}``-stable).
    ``witness_lo`` / ``witness_hi`` are ``(rank, c1)`` of exceptional bundles
    computing the endpoints in the sense of ``thm-stabilityInterval`` (``None``
    at ``0`` / ``infinity`` endpoints, and for transported ``e >= 2`` left
    endpoints).
    """

    lo: Fraction
    hi: Optional[Fraction]
    empty: bool
    witness_lo: Optional[Tuple[int, Tuple[int, int]]]
    witness_hi: Optional[Tuple[int, Tuple[int, int]]]
    e: int
    r: int
    c1: Tuple[int, int]
    certificate: Certificate

    def contains(self, m: Number) -> bool:
        """Open-interval membership ``m in I_v``."""
        m = _Q(m)
        if self.empty or m <= self.lo:
            return False
        return self.hi is None or m < self.hi


def _exceptional_ch2(r: int, c1: Sequence[int], surface: Surface) -> Fraction:
    """The ``ch2`` forced by ``chi(v,v) = 1``: ``Delta = (1 - 1/r^2)/2``."""
    nu = tuple(Fraction(x, r) for x in c1)
    delta = _HALF * (1 - Fraction(1, r * r))
    return r * (_HALF * surface.lattice.self_pairing(nu) - delta)


def _cheap_member(r_w: int, c1_w: Tuple[int, int], xi_V: SurfaceBundle,
                  surface: Surface) -> bool:
    """The polarization-independent ``S_V`` conditions: ``W`` potentially
    exceptional (integral) and ``chi(W, V) > 0``."""
    if not is_potentially_exceptional(r_w, c1_w, surface):
        return False
    W = SurfaceBundle(r_w, c1_w, _exceptional_ch2(r_w, c1_w, surface))
    return chi(W, xi_V, surface) > 0


def _member(m: Fraction, e: int, r_w: int, c1_w: Tuple[int, int],
            xi_V: SurfaceBundle, surface: Surface) -> bool:
    """The three ``S_V`` conditions beyond the slope equation: ``W`` potentially
    exceptional (integral), ``chi(W, V) > 0``, and ``m in I_W`` (=
    ``is_stable_exceptional`` at ``H_m`` by ``cor-DLPexcdelPezzo``; trivial for
    rank one -- line bundles are mu-stable for every polarization)."""
    if not _cheap_member(r_w, c1_w, xi_V, surface):
        return False
    if r_w == 1:
        return True
    return bool(is_stable_exceptional(r_w, c1_w, surface_with_index(e, m)))


def _line_probe(nu_V: Tuple[Fraction, Fraction], xi_V: SurfaceBundle, e: int,
                surface: Surface, vertical: bool) -> Tuple[Fraction, Tuple[int, int]]:
    """A far-out line bundle computing an element of ``S_V`` on the requested
    side (``rem-stabilityIntervalCompute``): ``vertical=True`` walks the strip
    ``a in (-1,0), b -> +inf`` (elements ``m -> infinity``, so past ``m_0``);
    ``vertical=False`` walks ``a -> +inf, b in (-1,0)`` (elements ``m -> 0``).

    Paper coordinates: the slope difference is ``xi = nu(V) - c1(L) = aE + bF``
    with ``a`` the package s-component and ``b`` the package f-component;
    ``m_{V,L} = -b/a``.  The strip pins one coordinate of ``c1(L)`` to the
    unique integer inside the unit window (exceptional slopes of rank >= 2 have
    non-integral coordinates -- ``rem-notStableForever``), and Riemann-Roch
    makes ``chi(L, V) > 0`` for the far tail.  A loud error after
    ``_PROBE_CAP`` steps rather than a silent wrong window.
    """
    nu_f, nu_s = nu_V
    m0 = Fraction(2 - e, 2)
    if vertical:
        if nu_s.denominator == 1:
            raise AssertionError("exceptional slope with integral s-coordinate (rank >= 2)")
        s_L = -((-nu_s.numerator) // nu_s.denominator)      # ceil: a = nu_s - s_L in (-1,0)
        f_L = nu_f.numerator // nu_f.denominator            # floor: b = nu_f - f_L > 0 start
        for _ in range(_PROBE_CAP):
            a, b = nu_s - s_L, nu_f - f_L
            if b > 0:
                m = -b / a                                   # = b/|a| > 0
                if m > m0 and _member(m, e, 1, (f_L, s_L), xi_V, surface):
                    return m, (f_L, s_L)
            f_L -= 1
    else:
        if nu_f.denominator == 1:
            raise AssertionError("exceptional slope with integral f-coordinate (rank >= 2)")
        f_L = -((-nu_f.numerator) // nu_f.denominator)      # ceil: b = nu_f - f_L in (-1,0)
        s_L = nu_s.numerator // nu_s.denominator            # floor: a = nu_s - s_L > 0 start
        for _ in range(_PROBE_CAP):
            a, b = nu_s - s_L, nu_f - f_L
            if a > 0:
                m = -b / a                                   # = |b|/a > 0
                if m < m0 and _member(m, e, 1, (f_L, s_L), xi_V, surface):
                    return m, (f_L, s_L)
            s_L -= 1
    raise AssertionError(
        "no line-bundle S_V member found within the probe cap -- "
        "rem-stabilityIntervalCompute guarantees one exists")


def _int_range(lo: Fraction, hi: Fraction):
    """Integers in the OPEN interval (lo, hi)."""
    n = lo.numerator // lo.denominator + 1                   # floor(lo) + 1 > lo
    while n < hi:
        yield n
        n += 1


def stability_interval(r: int, c1: Sequence[Number], surface: Surface) -> StabilityInterval:
    """The exact stability interval ``I_v`` of the exceptional character
    ``(r, c1)`` on the ``F_e`` family of ``surface`` (carried polarization
    ignored).  Raises ``ValueError`` for a non-potentially-exceptional or
    non-integral character; rank 1 returns ``(0, infinity)``.
    """
    e = hirzebruch_index(surface)
    r_frac = _Q(r)
    if r_frac.denominator != 1 or r_frac <= 0:
        raise ValueError(f"rank must be a positive integer; got {r!r}")
    r = int(r_frac)
    c1_t = tuple(_Q(x) for x in c1)
    if len(c1_t) != 2 or any(x.denominator != 1 for x in c1_t):
        raise ValueError(f"c1 must be an integral length-2 NS-vector; got {c1!r}")
    c1_i: Tuple[int, int] = (int(c1_t[0]), int(c1_t[1]))
    memo_key = (e, r, c1_i)
    if memo_key in _INTERVAL_MEMO:
        return _INTERVAL_MEMO[memo_key]

    def _cert(note: str) -> Certificate:
        return Certificate(rigor=Rigor.PROVEN, hypotheses=(
            "I_v = {m : Delta >= DLP^{<r}_{H_m}(nu)} (cor-DLPexcdelPezzo / cor-DLPExceptional)",
            "endpoints by thm-stabilityInterval + rem-stabilityIntervalCompute",
        ), citations=_CITATIONS, note=note)

    def _done(res: StabilityInterval) -> StabilityInterval:
        _INTERVAL_MEMO[memo_key] = res
        return res

    if r == 1:                                               # line bundles: always mu-stable
        return _done(StabilityInterval(Fraction(0), None, False, None, None,
                                       e, r, c1_i, _cert("rank 1: I = (0, infinity)")))
    if not is_potentially_exceptional(r, c1_i, surface):
        raise ValueError(
            f"(r={r}, c1={c1_i}) is not a potentially exceptional character on F_{e}")

    if e >= 2:
        # cor-highermus transport: I_v = {t > 0 : t + 1 in I_{pi(v)}} on F_{e-2}.
        base = P1xP1 if e - 2 == 0 else hirzebruch(e - 2)
        inner = stability_interval(r, pi_c1(c1_i), base)
        lo = max(Fraction(0), inner.lo - 1)
        hi = None if inner.hi is None else inner.hi - 1
        empty = inner.empty or (hi is not None and hi <= lo)
        w_hi = None
        if not empty and inner.witness_hi is not None:
            rw, cw = inner.witness_hi                        # transport the witness up: pi^{-1}
            w_hi = (rw, (cw[0] + cw[1], cw[1]))
        return _done(StabilityInterval(lo, hi, empty, None, w_hi, e, r, c1_i,
                                       _cert(f"transported from F_{e - 2} (cor-highermus)")))

    m0 = Fraction(2 - e, 2)
    if not is_stable_exceptional(r, c1_i, surface_with_index(e, m0)):
        # no exceptional bundle of this character exists on the del Pezzo (or it
        # would be mu_{-K}-stable, Gorodentsev): DLP monotonicity empties I_v.
        return _done(StabilityInterval(Fraction(0), Fraction(0), True, None, None,
                                       e, r, c1_i, _cert("no -K-stable exceptional: I empty")))

    nu_V = (Fraction(c1_i[0], r), Fraction(c1_i[1], r))      # (f, s) components
    xi_V = SurfaceBundle(r, c1_i, _exceptional_ch2(r, c1_i, surface))

    # Far line-bundle members bracketing m0 (rem-stabilityIntervalCompute):
    M1, wit1 = _line_probe(nu_V, xi_V, e, surface, vertical=True)
    if e == 0:
        M0, wit0 = _line_probe(nu_V, xi_V, e, surface, vertical=False)
    else:
        M0, wit0 = Fraction(0), None                          # e = 1: left bound 0

    # Bounded candidate search (cheap, polarization-independent conditions only):
    # slope differences xi = aE + bF in the two strips with m = -b/a in (M0, M1).
    # Vertical strip: a in (-1,0), b in (0, M1); horizontal strip: b in (-1,0),
    # a in (0, A_max) with A_max = 1/M0 (e = 0, from m = |b|/a > M0, |b| < 1) or
    # 2 (e = 1, the P > 0 triangle).
    a_max = (1 / M0) if e == 0 else Fraction(2)
    cands = []                                               # (m, r_w, c1_w), cheap-validated
    for r_w in range(1, r):
        nf, ns = r_w * nu_V[0], r_w * nu_V[1]
        # vertical strip: c1_s in (ns, ns + r_w), c1_f in (nf - r_w*M1, nf)
        for c1_s in _int_range(ns, ns + r_w):
            for c1_f in _int_range(nf - r_w * M1, nf):
                a, b = nu_V[1] - Fraction(c1_s, r_w), nu_V[0] - Fraction(c1_f, r_w)
                if not (-1 < a < 0 and b > 0):
                    continue
                m = -b / a
                if not (M0 < m < M1) or m == m0:
                    continue
                if _cheap_member(r_w, (c1_f, c1_s), xi_V, surface):
                    cands.append((m, r_w, (c1_f, c1_s)))
        # horizontal strip: c1_s in (ns - r_w*a_max, ns), c1_f in (nf, nf + r_w)
        for c1_s in _int_range(ns - r_w * a_max, ns):
            for c1_f in _int_range(nf, nf + r_w):
                a, b = nu_V[1] - Fraction(c1_s, r_w), nu_V[0] - Fraction(c1_f, r_w)
                if not (a > 0 and -1 < b < 0):
                    continue
                if e == 1 and (b + 1 - a / 2) <= 0:          # P(xi) > 0 triangle cut
                    continue
                m = -b / a
                if not (M0 < m < M1) or m == m0:
                    continue
                if _cheap_member(r_w, (c1_f, c1_s), xi_V, surface):
                    cands.append((m, r_w, (c1_f, c1_s)))

    # Endpoint search: the endpoint is the S_V member NEAREST m0 on each side,
    # so walk the cheap-validated candidates inward-first and stop at the first
    # whose full membership (m in I_W) holds.  Two prunes keep the expensive
    # is_stable_exceptional calls to a handful:
    #
    # * a rank-1 candidate certifies FREE (line bundles are mu-stable for every
    #   polarization), ending the walk at its value;
    # * a rank >= 2 candidate W is first filtered by ITS OWN line-bundle probe
    #   window: I_W is the component of m0 in the complement of S_W, and W's
    #   far probes are members of S_W, so I_W is contained in the open window
    #   between them -- a sound necessary condition costing microseconds.
    #   (A probe assertion -- e.g. an integral slope coordinate, where no
    #   exceptional bundle exists -- falls through to the honest full test.)
    cands.append((M1, 1, wit1))
    if wit0 is not None:
        cands.append((M0, 1, wit0))

    _probe_window_cache = {}

    def _passes_probe_window(r_w: int, c1_w: Tuple[int, int], m: Fraction) -> bool:
        key = (r_w, c1_w)
        if key not in _probe_window_cache:
            try:
                nu_W = (Fraction(c1_w[0], r_w), Fraction(c1_w[1], r_w))
                xi_W = SurfaceBundle(r_w, c1_w, _exceptional_ch2(r_w, c1_w, surface))
                hi_W, _ = _line_probe(nu_W, xi_W, e, surface, vertical=True)
                lo_W = (_line_probe(nu_W, xi_W, e, surface, vertical=False)[0]
                        if e == 0 else Fraction(0))
                _probe_window_cache[key] = (lo_W, hi_W)
            except AssertionError:
                _probe_window_cache[key] = None              # no filter: decide honestly
        win = _probe_window_cache[key]
        return win is None or (win[0] < m < win[1])

    def _first_member(seq):
        for m, r_w, c1_w in seq:
            if r_w == 1:
                return m, (r_w, c1_w)                        # line bundles: free
            if not _passes_probe_window(r_w, c1_w, m):
                continue
            # m in I_W by the MEMOIZED rank induction (rank strictly descends,
            # so the recursion is well-founded and each character's interval is
            # computed once per session) -- the same set as
            # is_stable_exceptional at H_m by cor-DLPexcdelPezzo, but a lookup
            # instead of a fresh DLP evaluation per (W, m).
            if stability_interval(r_w, c1_w, surface).contains(m):
                return m, (r_w, c1_w)
        return None

    hi_m, w_hi = _first_member(sorted(c for c in cands if c[0] > m0))
    hit_lo = _first_member(sorted((c for c in cands if c[0] < m0), reverse=True))
    if hit_lo is None:
        lo, w_lo = Fraction(0), None                         # stable down to 0 (e = 1 rows)
        if e == 0:
            raise AssertionError("e = 0 must have a finite left endpoint (the M0 probe)")
    else:
        lo, w_lo = hit_lo
    return _done(StabilityInterval(lo, hi_m, False, w_lo, w_hi,
                                   e, r, c1_i, _cert("thm-stabilityInterval rank induction")))
