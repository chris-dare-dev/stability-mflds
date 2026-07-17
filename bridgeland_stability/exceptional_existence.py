"""E15-M1 / G18+: existence obstructions for exceptional bundles on ``F_e``.

The Sec. 11 conjecture of arXiv:1907.06739 (every exceptional bundle on ``F_e``,
``e >= 2``, is slope-stable near ``H_0``) reduces, for a *potentially exceptional*
character ``v`` whose stability interval is EMPTY, to the question: does an
exceptional bundle of character ``v`` exist at all?  The paper's own first open
case is ``v_107 = (107, 25/107 E + 76/107 F, 5724/11449)`` on ``F_3``.

Two necessary conditions are implemented, both valid on every ``F_e``:

1. **The prioritary index** (the paper's F_4-example route): an exceptional
   bundle is simple, hence ``H_2``-prioritary (``lem-simple``), and rigid, hence
   its point is dense-open in Walter's irreducible stack ``P_F(v)`` -- so the
   GENERAL sheaf is ``H_2``-prioritary and ``rho_gen(v) >= 2``
   (``cor-prioritaryRho`` / ``prop-excPrior``).  ``rho_gen(v) = 1`` refutes.
   [Recorded fact: ``rho_gen(v_107) = 2`` -- this route is INCONCLUSIVE for the
   paper's candidate, which is presumably why it was left open.]

2. **The rigid-factor obstruction** (:func:`exceptional_refutation`, new here):
   if the exceptional bundle ``V`` exists, then ``V`` *is* the general sheaf of
   ``P_F(v)`` (rigidity is open and the stack is irreducible), so the generic
   ``H_m``-Harder-Narasimhan factors computed by the Sec. 5 algorithm are the
   factors of ``V`` itself, at every ``m``.  Along any Gieseker-HN filtration on
   ``F_e`` the hypotheses of ``prop-mukai``(2) hold unconditionally
   (``Hom(W, U) = 0`` by the reduced-Hilbert ordering; ``Ext^2(U, W) =
   Hom(W, U(K))* = 0`` since twisting by ``K`` strictly lowers ``q`` for ample
   ``H``), so every factor of the rigid ``V`` is RIGID, giving the numerical
   consequence

       chi(gr_i, gr_i) = hom + ext^2 >= 1,   i.e.   Delta_i <= (1 - 1/r_i^2)/2

   for every factor -- testable exactly.  Moreover, at a CHAMBER-GENERIC sample
   (an irrational-equivalent polarization, via the E14-M1 chamber gap) a
   length-ONE filtration also refutes when the slope denominator equals the
   rank: ``V`` would be Gieseker-semistable there, hence mu-semistable, hence
   mu-STABLE (no proper subsheaf can match a slope of exact denominator ``r``
   at a generic polarization), contradicting the empty stability interval.

   NOTE the del Pezzo theorem ``thm-rigidSplit`` (Kuleshov-Orlov: rigid bundles
   split into exceptionals) is NOT available on ``F_e`` for ``e >= 2`` -- the
   obstruction above deliberately uses only ``prop-mukai``, which holds on any
   smooth surface.

A refutation is PROVEN (the conjecture holds for that character); a pass is
honestly INCONCLUSIVE -- existence is not decided.  Exact ``Fraction``
arithmetic; stdlib-only at import.  Record: ``docs/CORRECTIONS.md`` Sec. 21.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Optional, Sequence, Tuple

from .chern import Number
from .varieties import Surface
from .exceptional_surface import SurfaceBundle
from .dlp_hirzebruch import (
    discriminant,
    hirzebruch_index,
    is_potentially_exceptional,
)
from .delta_sharp import _chamber_gap, surface_with_index
from .generic_hn import generic_hn_factors
from .prioritary import generic_prioritary_index
from .rigor import Certificate, Rigor
from .stability_interval import stability_interval

__all__ = ["ExceptionalRefutation", "exceptional_refutation"]

_HALF = Fraction(1, 2)

_CITATIONS = (
    "arXiv:1907.06739 prop-mukai (Mukai/Gorodentsev; any smooth surface), lem-simple, "
    "prop-excPrior, cor-prioritaryRho, Sec. 5 (the generic HN algorithm), Sec. 11 conjecture",
    "Walter (irreducibility of the prioritary stack); E14-M1 chamber gap (CORRECTIONS Sec. 17)",
)


def _Q(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


@dataclass(frozen=True)
class ExceptionalRefutation:
    """Outcome of the E15-M1 necessary-condition battery for "an exceptional
    bundle of character ``(r, c1)`` exists".

    ``refuted = True`` is PROVEN non-existence (``reason`` says which condition
    fired); ``refuted = False`` means every implemented condition passed --
    honestly INCONCLUSIVE, never evidence of existence.  ``rho_gen`` and
    ``samples`` (per sampled polarization: ``(anchor m, computed length, factor
    tuples with their discriminant bounds)``) record the evidence either way.
    """

    refuted: bool
    reason: str
    rho_gen: int
    interval_empty: bool
    samples: Tuple[tuple, ...]
    r: int
    c1: Tuple[int, int]
    e: int
    certificate: Certificate


def exceptional_refutation(r: int, c1: Sequence[Number], surface: Surface,
                           anchors: Sequence[Number] = (Fraction(1),),
                           ) -> ExceptionalRefutation:
    """Run the necessary-condition battery for existence of an exceptional
    bundle of character ``(r, c1)`` on the ``F_e`` family of ``surface``.

    Preconditions checked: the character must be potentially exceptional
    (otherwise no exceptional bundle by definition -- refuted trivially).  The
    length-one branch additionally requires the character's stability interval
    to be EMPTY and its slope denominator to equal ``r``; where those fail the
    branch is skipped (recorded), never guessed.
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

    def _cert(note: str) -> Certificate:
        return Certificate(rigor=Rigor.PROVEN, hypotheses=(
            "rigid => generic in the irreducible prioritary stack",
            "prop-mukai(2) along the Gieseker-HN filtration (any F_e, any m)",
            "rigid factor => chi(gr,gr) >= 1 => Delta <= (1 - 1/rank^2)/2",
        ), citations=_CITATIONS, note=note)

    if not is_potentially_exceptional(r, c1_i, surface):
        return ExceptionalRefutation(
            True, "not potentially exceptional (chi(v,v) != 1 or non-integral c2)",
            0, False, (), r, c1_i, e, _cert("refuted: invalid character"))

    # forced exceptional discriminant and the exact round-trip slope
    nu = (Fraction(c1_i[0], r), Fraction(c1_i[1], r))
    delta = _HALF * (1 - Fraction(1, r * r))

    # Condition 1: the prioritary index (the paper's F_4 route).
    rho = generic_prioritary_index(nu, delta, surface)
    if rho < 2:
        return ExceptionalRefutation(
            True, f"rho_gen = {rho} < 2 (exceptional => H_2-prioritary => rho_gen >= 2)",
            rho, False, (), r, c1_i, e,
            _cert("refuted by cor-prioritaryRho / prop-excPrior"))

    # Condition 2: the rigid-factor obstruction.
    interval = stability_interval(r, c1_i, surface)
    denom = 1
    for x in nu:
        denom = denom * x.denominator // gcd(denom, x.denominator)
    length_one_refutes = interval.empty and denom == r

    # the exact round-trip ch2 of the forced exceptional discriminant
    ch2 = r * (_HALF * surface.lattice.self_pairing(nu) - delta)

    samples = []
    for anchor in anchors:
        anchor = _Q(anchor)
        g = _chamber_gap(r, e, anchor, side_left=False)
        m = anchor + g / 2                       # chamber-generic sample
        S = surface_with_index(e, m)
        # The Sec. 5 filtration of the H_{ceil(m)}-prioritary stack directly
        # (hn_verdict early-exits through prop-ssPrior when P_{H_{ceil(m)+1}}
        # is empty, which is exactly the rho_gen = 2 regime of interest; the
        # generic filtration of P_{H_{ceil(m)}} still exists -- the paper's
        # Sec. 1.4 remark).  V exceptional would be simple, hence D-prioritary
        # for every D with -(K+D) effective nontrivial (lem-simple), in
        # particular H_1- and H_2-prioritary -- so for anchors < 2 the sheaf V
        # would be a GENERAL member of the filtered stack.
        if anchor >= 2:
            raise ValueError(
                f"anchor {anchor} >= 2: the sample must keep ceil(m) <= 2, the "
                "prioritary regime lem-simple guarantees for an exceptional V")
        facts_v = generic_hn_factors(r, c1_i, ch2, S)
        if facts_v is None:
            samples.append((anchor, None, None))
            continue
        facts = []
        for (r_i, c1_f, ch2_i) in facts_v:
            d_i = discriminant(SurfaceBundle(r_i, c1_f, ch2_i), S)
            bound = _HALF * (1 - Fraction(1, r_i * r_i))
            facts.append((r_i, tuple(c1_f), d_i, bound))
            if d_i > bound:
                return ExceptionalRefutation(
                    True, f"generic HN factor (r={r_i}, c1={tuple(c1_f)}) at the sample "
                          f"near m={anchor} has Delta = {d_i} > {bound}: it cannot be "
                          "rigid, but every factor of a rigid sheaf is rigid "
                          "(prop-mukai)",
                    rho, interval.empty, tuple(samples), r, c1_i, e,
                    _cert("refuted: non-rigid generic HN factor"))
        if len(facts_v) == 1 and length_one_refutes:
            return ExceptionalRefutation(
                True, f"generic HN length 1 at the generic sample near m={anchor} with "
                      "prime slope denominator: V would be mu-stable, but the "
                      "stability interval is empty",
                rho, interval.empty, tuple(samples), r, c1_i, e,
                _cert("refuted: length-one at a generic polarization"))
        samples.append((anchor, len(facts_v), tuple(facts)))

    return ExceptionalRefutation(
        False, "every implemented necessary condition passed -- existence UNDECIDED",
        rho, interval.empty, tuple(samples), r, c1_i, e,
        _cert("inconclusive: no condition fired"))
