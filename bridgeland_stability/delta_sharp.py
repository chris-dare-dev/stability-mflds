"""E14-M1 / G18: the sharp Bogomolov function ``delta_m^{mu-s}(nu)`` as a computable object.

The paper's headline object (arXiv:1907.06739 ``def-deltass``) is

    delta_m^{mu-s}(nu) = inf{ Delta >= 1/2 : there is a mu_{H_m}-stable sheaf
                              of total slope nu and discriminant Delta }

on ``F_e`` with ``H_m = E + (e+m)F``.  E13 made per-character *Gieseker-semistable*
existence decidable (``hn_verdict``); this module upgrades that to the
*slope-stable* question and evaluates the function itself as a certified sandwich.

Two exact deliverables
----------------------

1. :func:`mu_stable_exists` -- decide "is there a mu_{H_m}-stable sheaf of character
   ``(r, nu, Delta)``?" at a RATIONAL ``m``.  Rational polarizations are never
   "generic" in the paper's sense (an integral class orthogonal to ``H_m`` always
   exists), so the generic-polarization bridges ``prop-ssIMPs`` / ``prop-sIMPmus``
   never apply at ``m`` directly.  The decision instead goes through the **generic
   stability interval** ``I_v = {m > 0 : the general sheaf of P_F(v) is
   mu_{H_m}-stable}``:

   * existence of a mu_{H_m}-stable sheaf of character ``v``  <=>  ``m in I_v``
     (slope stability is an open property in flat families and ``P_F(v)`` is
     irreducible by Walter's theorem, so one stable sheaf makes the general sheaf
     stable);
   * ``I_v`` is convex (slope stability for two ample classes implies it for every
     positive combination, and two dense-open loci of an irreducible stack meet),
     open, and -- on a del Pezzo ``F_e``, ``e in {0,1}`` -- contains the
     anticanonical index ``m_0 = 1 - e/2`` whenever it is nonempty
     (``cor-KstabilityEasy``);
   * for ``Delta > 1/2`` and an IRRATIONAL ``m'`` (hence generic),
     Gieseker-semistable existence at ``m'`` implies mu_{H_{m'}}-stable existence
     (``prop-ssIMPs`` then ``prop-sIMPmus``).

   The certifier therefore samples the FIRST WALL-FREE CHAMBER beside ``m`` (on the
   side away from ``m_0``): Gieseker-semistable existence is constant on
   ``(m, m + g)`` for the explicit gap ``g`` of :func:`_chamber_gap` (every
   condition of the Sec. 5 criterion ``thm-HNcriterion`` / ``cor-algorithm`` can
   only flip where a bounded-window slope relation crosses ``0`` or ``1``, where
   the ``lem-slopeQuad`` F-window boundary crosses a candidate, or at an integer
   prioritary index -- all rational points bounded away from ``m`` by a lattice
   argument; see ``docs/CORRECTIONS.md`` Sec. 17).  One exact ``hn_verdict`` call at
   the rational midpoint of that chamber then decides BOTH directions:

   * semistable there  =>  mu-stable at an irrational point of the chamber
     =>  (convexity through ``m_0``) ``m in I_v``  =>  exists at ``m``;
   * not semistable there  =>  no mu-stable sheaf anywhere in the open chamber
     =>  ``I_v`` misses the chamber  =>  ``m not in I_v`` (``I_v`` is open)
     =>  no mu-stable sheaf at ``m``.

   For ``e >= 2`` the question transports to the del Pezzo base first:
   ``cor-highermus`` makes mu_{A_m}-stable existence equivalent under the E13-M1
   reduction ``pi``, and every strictly ample ``H`` on ``F_e`` is an ``A_m`` in the
   transported range.  The ``Delta < 1/2`` band is decided without sampling:
   a mu_H-stable sheaf is simple with ``Ext^2(V,V) = 0`` (``K.H < 0``), so
   ``chi(v,v) = r^2(1 - 2 Delta) <= 1``; ``chi >= 2`` refuses, and ``chi = 1`` is
   ``cor-DLPExceptional`` (:func:`~bridgeland_stability.dlp_hirzebruch.
   is_stable_exceptional`).  ``Delta = 1/2`` exactly is returned as an honest
   ``None`` (out of certification scope), never guessed.

2. :func:`delta_mu_stable` -- the function value as a certified sandwich
   ``DeltaSharp(lower, upper, exact)``:

   * ``lower  = max(1/2, DLP_{H_m}(nu) truncation, delta^p_{ceil(m)+1}(nu))`` --
     certified by ``cor-deltaDLPe`` and the chain ``delta^p_{ceil(m)+1} <=
     delta^{mu-ss}_m <= delta^{mu-s}_m`` (Sec. 3.1 of the paper);
   * ``upper``  = the least lattice discriminant ``Delta > 1/2`` of a scanned rank
     ``r <= max_rank`` where :func:`mu_stable_exists` decides True -- an upper
     bound because the inf ranges over a superset;
   * ``exact = (upper == lower)``.

   **The inf need not be attained** (``def-deltass`` is an infimum): at the paper's
   Sec. 8 Kronecker values (``delta^{mu-s}_{25/9} = 3/5`` on ``F_0``,
   ``delta^{mu-s}_{12/7} = 98/169`` on ``F_1``) the general sheaf AT the sharp
   discriminant is strictly mu-semistable and the paper reaches the value as a
   limit ``Delta_{m+eps} -> delta`` of nearby-slope characters of growing rank.
   The scan's ``upper`` therefore converges to but NEVER equals such values at any
   finite ``max_rank``; closing that gap exactly is E14-M2 (``thm-deltaKronecker``).

Exact ``Fraction`` arithmetic throughout; stdlib-only at import (zero-dependency
invariant).  Primary source: arXiv:1907.06739 -- ``def-deltass``, Sec. 3.1 chain,
``thm-deltaGeneric``, ``prop-ssIMPs``, ``prop-sIMPmus``, ``cor-KstabilityEasy``,
``cor-deltaDLPe``, ``thm-deltaSurface``, ``cor-highermus``; Walter,
"Irreducibility of moduli spaces of vector bundles on birationally ruled surfaces".
Record: ``docs/CORRECTIONS.md`` Sec. 17.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Optional, Sequence, Tuple

from .chern import Number
from .varieties import Surface
from .exceptional_surface import SurfaceBundle, chi
from .dlp_hirzebruch import (
    DEFAULT_RANK_MAX,
    dlp_envelope,
    hirzebruch_index,
    is_ample,
    is_stable_exceptional,
)
from .hn_filtration import hn_verdict
from .nonemptiness_rational import hirzebruch_with_polarization
from .prioritary import delta_prioritary
from .reduction import reduce_to_del_pezzo
from .rigor import Certificate, Rigor

__all__ = [
    "DeltaSharp",
    "delta_mu_stable",
    "mu_stable_exists",
    "polarization_index",
    "surface_with_index",
]

_HALF = Fraction(1, 2)

_CITATIONS = (
    "arXiv:1907.06739 def-deltass, thm-deltaGeneric, prop-ssIMPs, prop-sIMPmus",
    "arXiv:1907.06739 cor-KstabilityEasy, cor-deltaDLPe, thm-deltaSurface, cor-highermus",
    "Walter, Irreducibility of moduli spaces of vector bundles on birationally ruled surfaces",
)


def _Q(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


def _ceil(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def polarization_index(surface: Surface) -> Fraction:
    """The Coskun-Huizenga index ``m`` of the carried polarization: ``H = (a, b)``
    is the ray of ``H_m = E + (e+m)F`` with ``m = a/b - e``.  Requires strict
    ampleness (``m > 0``)."""
    e = hirzebruch_index(surface)
    if not is_ample(surface):
        raise ValueError(
            f"{surface.name}: needs a strictly ample H (the factory H is nef-and-big; "
            "use nonemptiness_rational.hirzebruch_with_polarization)")
    a, b = surface.H
    return Fraction(a, b) - e


def surface_with_index(e: int, m: Fraction) -> Surface:
    """The ample ``F_e`` carrying (the primitive integral point of) the ray ``H_m``,
    ``m = p/q > 0`` rational: ``H = (p + e q, q)`` (coprime since ``gcd(p, q) = 1``)."""
    m = _Q(m)
    if m <= 0:
        raise ValueError(f"H_m is ample iff m > 0; got m = {m}")
    p, q = m.numerator, m.denominator
    return hirzebruch_with_polarization(e, (p + e * q, q))


def _chamber_gap(r: int, e: int, m: Fraction, side_left: bool) -> Fraction:
    """A PROVEN lower bound ``g`` on the distance from ``m`` to the nearest value
    ``m' != m`` (on the sampled side) where Gieseker-semistable existence of a
    rank-``r`` character can change.

    Every condition of the Sec. 5 generic-HN criterion is a comparison whose
    outcome, as ``m'`` varies, can only flip where one of the following holds for
    a slope difference ``xi = nu_w - nu_u`` of two characters of the recursion
    (all ranks ``<= r``, so both coordinates of ``xi`` lie in ``(1/r^2) Z``, and
    ``|xi . F| <= 8 Ymax`` -- the ``lem-slopeQuad`` window ``max{1, 2/(e+2m')}``
    stacked over the recursion depth ``<= 4`` of the ``ell <= 4`` lemma, doubled
    for pairwise differences):

    * ``xi . H_{m'} = 0`` or ``= 1``  (slope crossings; q-key first-component ties;
      the ``lem-HNclose`` window boundary; the slope-spread-``<= 1`` bound).  In
      package ``(f, s)`` coordinates ``xi . H_{m'} = x + y m'``, so the crossing is
      at ``m' = -x/y`` resp. ``(1-x)/y`` and ``|m' - m| = |t| / |y|`` with
      ``t = x + m y`` (resp. ``- 1``) a nonzero element of ``(1/(r^2 q)) Z``
      (``q = den(m)``); hence ``|m' - m| >= 1/(8 Ymax r^2 q)``.
    * ``|xi . F| = 2/(e + 2m')`` (the F-window boundary moving past a candidate):
      ``m' - m`` is then a nonzero rational with numerator a nonzero integer over
      ``2 q |y| r^2 <= 16 Ymax r^2 q`` -- same bound, up to the factor 2.
    * an integer (the prioritary index ``ceil(m')`` and the ``delta^p_n`` inputs
      jump): at distance ``>= 1/q`` from a non-integral ``m`` and ``>= 1`` from an
      integral ``m``.

    ``g = 1/(32 Ymax r^2 q)`` under-runs all three.  Left samples are additionally
    capped at ``m/2`` (the sample must stay ample), and ``Ymax`` is evaluated at
    the far end of the sampled side so the window bound is conservative.
    """
    q = m.denominator
    m_far = m / 2 if side_left else m           # smallest m' the sample may reach
    y_max = max(Fraction(1), Fraction(2, 1) / (e + 2 * m_far))
    g = Fraction(1, 32) / (y_max * r * r * q)
    if side_left:
        g = min(g, m / 2)
    return g


def _integral_character(r: int, nu: Tuple[Fraction, ...], Delta: Fraction,
                        surface: Surface) -> Optional[Tuple[Tuple[int, ...], Fraction]]:
    """``(c1, ch2)`` of the exact round-trip character, or ``None`` if no coherent
    sheaf carries it (non-integral ``c1``, or ``c2 = <c1,c1>/2 - ch2`` non-integral)."""
    c1 = tuple(r * x for x in nu)
    if any(x.denominator != 1 for x in c1):
        return None
    c1 = tuple(int(x) for x in c1)
    ch2 = r * (_HALF * surface.lattice.self_pairing(nu) - Delta)
    c2 = _HALF * surface.lattice.self_pairing(c1) - ch2
    if c2.denominator != 1:
        return None
    return c1, ch2


def mu_stable_exists(r: int, nu: Sequence[Number], Delta: Number,
                     surface: Surface) -> Optional[bool]:
    """Is there a ``mu_{H_m}``-stable sheaf of character ``(r, nu, Delta)`` on
    ``surface`` (an ``F_e`` carrying a strictly ample ``H``, whose index ``m`` is
    read off the carried polarization)?

    Returns ``True`` / ``False`` as a PROVEN verdict, or ``None`` exactly on the
    one undecided band ``Delta == 1/2`` with ``r >= 2`` (see the module docstring;
    never a guess).  Decision routes, in order:

    * invalid character (non-integral ``c1`` or ``c2``), or ``Delta < 0``
      (Bogomolov for the mu-semistable closure)  ->  ``False``;
    * ``r == 1``: rank-one torsion-free sheaves are always mu-stable, and exist
      iff ``c2 = Delta`` is a non-negative integer  ->  that integrality test;
    * ``e >= 2``: transport along the E13-M1 reduction (``cor-highermus``) and
      recurse on the del Pezzo base;
    * ``Delta < 1/2``: ``chi(v,v) = r^2(1 - 2 Delta) >= 2`` refuses (a stable
      sheaf is simple with ``Ext^2 = 0``, so ``chi <= 1``); ``chi(v,v) = 1`` is
      decided by ``cor-DLPExceptional``;
    * ``Delta > 1/2``: the chamber-sample certificate of the module docstring --
      one ``hn_verdict`` call at the midpoint of the first wall-free chamber on
      the side of ``m`` away from the anticanonical index ``m_0 = 1 - e/2``.
    """
    e = hirzebruch_index(surface)
    m = polarization_index(surface)
    r_frac = _Q(r)
    if r_frac.denominator != 1 or r_frac <= 0:
        raise ValueError(f"rank must be a positive integer; got {r!r}")
    r = int(r_frac)
    Delta = _Q(Delta)
    nu_t = tuple(_Q(x) for x in nu)
    if len(nu_t) != 2:
        raise ValueError("nu must be a length-2 F_e NS-vector (f, s)")

    ic = _integral_character(r, nu_t, Delta, surface)
    if ic is None:
        return False                             # no coherent sheaf has this character
    c1, ch2 = ic

    if Delta < 0:
        return False                             # Bogomolov: mu-stable => mu-ss => Delta >= 0
    if r == 1:
        # rank-one torsion-free = line bundle tensor ideal sheaf: always mu-stable;
        # exists iff c2 = Delta is a non-negative integer (already integral here).
        return Delta.denominator == 1 and Delta >= 0

    if e >= 2:
        # cor-highermus: mu_{A_m}-stable existence is equivalent under pi, and every
        # strictly ample H on F_e is an A_m with m in the transported range.  Delta,
        # r, ch2 are pi-invariant (Lemma 11.3); the polarization transports with it.
        xi_red, surf_red = reduce_to_del_pezzo(SurfaceBundle(r, c1, ch2), surface)
        return mu_stable_exists(xi_red.r, tuple(_Q(c) / r for c in xi_red.c1),
                                Delta, surf_red)

    if Delta < _HALF:
        # A mu_H-stable V is simple and has Ext^2(V,V) = Hom(V, V(K))* = 0 since
        # K.H < 0 on an ample F_e, so chi(v,v) = 1 - ext^1 <= 1.  The Euler pairing
        # is evaluated by the package's general RR chi -- not transcribed.
        x = SurfaceBundle(r, c1, ch2)
        chi_vv = chi(x, x, surface)
        if chi_vv != r * r * (1 - 2 * Delta):    # loud tripwire, never silent
            raise AssertionError(
                f"chi(v,v) = {chi_vv} != r^2(1-2Delta) = {r * r * (1 - 2 * Delta)}")
        if chi_vv >= 2:
            return False
        return bool(is_stable_exceptional(r, c1, surface))
    if Delta == _HALF:
        return None                              # honest: outside the certified scope

    # Delta > 1/2: chamber-sample certificate (module docstring; CORRECTIONS Sec. 17).
    m0 = Fraction(2 - e, 2)                      # anticanonical index 1 - e/2
    side_left = m < m0
    g = _chamber_gap(r, e, m, side_left)
    sample = m - g / 2 if side_left else m + g / 2
    verdict = hn_verdict(r, nu_t, Delta, surface_with_index(e, sample))
    if verdict.exists is None:                   # never happens on the total scope
        return None
    return bool(verdict.exists)


@dataclass(frozen=True)
class DeltaSharp:
    """A certified two-sided estimate of ``delta_m^{mu-s}(nu)``.

    Attributes
    ----------
    lower : Fraction
        ``max(1/2, dlp_envelope value, delta^p_{ceil(m)+1})`` -- a certified lower
        bound (``cor-deltaDLPe`` + the Sec. 3.1 chain).
    upper : Fraction | None
        The least scanned lattice discriminant ``> 1/2`` certified to carry a
        mu_{H_m}-stable sheaf (an upper bound for the inf), or ``None`` if the
        scan certified nothing within ``max_rank`` / ``delta_cap``.
    exact : bool
        ``True`` iff the sandwich pinches (``upper == lower``).  ``False`` is NOT
        a statement that the value is unknown to mathematics -- only that this
        scan did not close the gap (see the module docstring on unattained infs).
    witness : tuple | None
        ``(r, c1, ch2)`` realizing ``upper``.
    nu, m, e, rank_scanned :
        The query, for auditability.
    lower_dlp, lower_prioritary : Fraction
        The two nontrivial ingredients of ``lower``, for transparency.
    certificate : Certificate
        PROVEN, with the load-bearing citations.
    """

    lower: Fraction
    upper: Optional[Fraction]
    exact: bool
    witness: Optional[Tuple[int, Tuple[int, ...], Fraction]]
    nu: Tuple[Fraction, ...]
    m: Fraction
    e: int
    rank_scanned: int
    lower_dlp: Fraction
    lower_prioritary: Fraction
    certificate: Certificate


def _lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def delta_mu_stable(nu: Sequence[Number], m: Number, surface: Surface,
                    max_rank: int = 16, delta_cap: Optional[Fraction] = None,
                    envelope_rank_max: int = DEFAULT_RANK_MAX) -> DeltaSharp:
    """The certified sandwich for ``delta_m^{mu-s}(nu)`` on the ``F_e`` family of
    ``surface`` (the carried polarization is ignored; ``H_m`` is built from ``m``).

    The scan visits ranks ``r`` that are multiples of the slope denominator
    ``r0`` up to ``max_rank`` and, per rank, ascends the integral-``c2`` lattice
    of discriminants starting strictly above ``1/2``, stopping at the first
    :func:`mu_stable_exists` certification (per-rank first hits are minima:
    slope-stability survives elementary modifications, so per rank the existence
    set is upward closed).  ``thm-deltaSurface`` (1) + the certifier's totality on
    ``Delta > 1/2`` guarantee the per-rank walk terminates at latest one lattice
    step above the true value; ``delta_cap`` (default ``lower + 8``) is a loud
    safety net -- ranks whose walk exceeds it contribute nothing, and if NO rank
    certifies, ``upper`` is an honest ``None``.

    Practical cost: each certification is one ``hn_verdict`` call, whose Sec. 5
    recursion tree grows steeply with the rank -- ranks up to ~15 decide in
    fractions of a second to seconds, while rank 30 takes tens of minutes.
    Choose ``max_rank`` (and slope denominators) accordingly.
    """
    e = hirzebruch_index(surface)
    m = _Q(m)
    nu_t = tuple(_Q(x) for x in nu)
    if len(nu_t) != 2:
        raise ValueError("nu must be a length-2 F_e NS-vector (f, s)")
    S = surface_with_index(e, m)

    env = dlp_envelope(nu_t, S, rank_max=envelope_rank_max)
    dp = delta_prioritary(nu_t, _ceil(m) + 1, S)
    lower = max(_HALF, env.value, dp)
    if delta_cap is None:
        delta_cap = lower + 8

    r0 = 1
    for x in nu_t:
        r0 = _lcm(r0, x.denominator)

    upper: Optional[Fraction] = None
    witness = None
    scanned = 0
    r = r0
    while r <= max_rank:
        scanned = r
        base = _HALF * S.lattice.self_pairing(nu_t) \
            - S.lattice.self_pairing(tuple(r * x for x in nu_t)) / (2 * r)
        # Delta(c2) = base + c2/r, c2 integral: start strictly above 1/2.
        c2 = _ceil((_HALF - base) * r)
        if base + Fraction(c2, r) <= _HALF:
            c2 += 1
        while True:
            D = base + Fraction(c2, r)
            if upper is not None and D >= upper:
                break
            if D > delta_cap:
                break                            # this rank certifies nothing below the cap
            if mu_stable_exists(r, nu_t, D, S):
                if upper is None or D < upper:
                    upper, witness = D, (r, tuple(int(r * x) for x in nu_t), r * (_HALF * S.lattice.self_pairing(nu_t) - D))
                break
            c2 += 1
        r += r0

    cert = Certificate(
        rigor=Rigor.PROVEN,
        hypotheses=(
            "lower: cor-deltaDLPe + delta^p_{ceil(m)+1} <= delta^{mu-ss} <= delta^{mu-s}",
            "upper: mu_stable_exists chamber-sample certificate (CORRECTIONS Sec. 17)",
            "def-deltass is an infimum; upper == value need not be attained at finite rank",
        ),
        citations=_CITATIONS,
        note="delta_m^{mu-s}(nu) sandwich (E14-M1)",
    )
    return DeltaSharp(
        lower=lower, upper=upper, exact=(upper is not None and upper == lower),
        witness=witness, nu=nu_t, m=m, e=e, rank_scanned=scanned,
        lower_dlp=env.value, lower_prioritary=dp, certificate=cert,
    )
