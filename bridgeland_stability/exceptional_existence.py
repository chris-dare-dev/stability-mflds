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
from .exceptional_surface import SurfaceBundle, chi
from .dlp_hirzebruch import (
    discriminant,
    hirzebruch_index,
    is_potentially_exceptional,
)
from .delta_sharp import _chamber_gap, surface_with_index
from .generic_hn import generic_hn_factors
from .prioritary import general_betti, generic_prioritary_index
from .rigor import Certificate, Rigor
from .stability_interval import stability_interval

__all__ = ["ExceptionalRefutation", "chi_box_conditions", "exceptional_refutation",
           "gaeta_star_conditions"]

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


def chi_box_conditions(r: int, c1: Sequence[int], surface: Surface
                       ) -> Tuple[Tuple[Tuple[int, int], int], ...]:
    """E15-M1d: the simple-prioritary χ-box — a family of PROVEN necessary
    conditions for the existence of an exceptional bundle of character
    ``(r, c1)``, one per divisor ``D`` in the finite box

        D effective nontrivial,  -(K + D) effective nontrivial.

    For such ``D`` and an exceptional (hence simple) ``V``:
    ``Hom(V, V(-D)) = 0`` (a nonzero map composed with the section injection
    ``V(-D) -> V`` would be a non-scalar endomorphism) and
    ``Ext^2(V, V(-D)) = Hom(V, V(K+D))* = 0`` (``lem-simple`` verbatim), so

        chi(v, v(-D)) = -ext^1(V, V(-D)) <= 0.

    Returns ``((D, chi), ...)`` for every box divisor; any ``chi > 0`` REFUTES
    existence.  Pure integer Riemann-Roch through the package ``chi`` — no
    recursion, no polarization.  NOTE the honest track record (CORRECTIONS
    Sec. 21): on every refutation case known so far (the paper's F_4 example)
    and on ``v_107`` itself, all box values are <= 0 — the family is strictly
    weaker than the ``rho_gen`` route on those, and ships as a cheap widener
    of the battery, not as a decided improvement.
    """
    e = hirzebruch_index(surface)
    lat = surface.lattice
    nu = tuple(Fraction(x, r) for x in c1)
    delta = _HALF * (1 - Fraction(1, r * r))
    ch2 = r * (_HALF * lat.self_pairing(nu) - delta)
    v = SurfaceBundle(r, tuple(c1), ch2)
    kf, ks = e + 2, 2                            # -K = (e+2, 2) in (f, s)
    rows = []
    for df in range(kf + 1):
        for ds in range(ks + 1):
            if (df, ds) in ((0, 0), (kf, ks)):
                continue
            D = (df, ds)
            c1t = (c1[0] - r * df, c1[1] - r * ds)
            ch2t = ch2 - lat.pairing(tuple(c1), D) + r * _HALF * lat.self_pairing(D)
            rows.append((D, chi(v, SurfaceBundle(r, c1t, ch2t), surface)))
    return tuple(rows)


def gaeta_star_conditions(r: int, c1: Sequence[int], surface: Surface
                          ) -> Tuple[Tuple[Tuple[int, int], int, int], ...]:
    """E15-M1e: the Gaeta dimension inequality — generic-``D``-prioritariness
    obstructions at EVERY box divisor, including the ``s``-coefficient-2 ones
    the ``H_n``-ray theory cannot see.

    Apply ``Hom(-, V(-D'))`` to the ``L_0``-Gaeta resolution of the general
    ``V in P_F(v)`` (``prop-Gaeta``; ``alpha >= 0`` is automatic for ``e >= 2``
    and holds whenever computed here — asserted): if ``V`` is
    ``D'``-prioritary, exactness + ``Ext^3 = 0`` force
    ``Ext^2(middle) -> Ext^2(kernel)`` injective, i.e.

        beta h2(V(-L0+E+eF-D')) + gamma h2(V(-L0+F-D')) + delta h2(V(-L0-D'))
            <= alpha h2(V(-L0+E+(e+1)F-D'))                            (star)

    with every ``h^2`` a Betti number of a twisted GENERAL sheaf (thm-BN,
    :func:`~bridgeland_stability.prioritary.general_betti`; twisting is an
    isomorphism of the irreducible stack, so the twist of a general sheaf is
    general).  An exceptional bundle of the character would be rigid, hence
    the general sheaf, hence ``D'``-prioritary for every box ``D'``
    (``prop-excPrior``) — so ``LHS > RHS`` at any box divisor REFUTES
    existence.  Returns ``((D', LHS, RHS), ...)``.

    Cross-validated (tests): on ``v_107`` the inequality holds at ``H_2`` and
    is VIOLATED at ``H_3``/``H_4`` — reproducing ``rho_gen = 2`` through
    independent machinery (Gaeta exponents + thm-BN vs ``cor-prioritaryRho``).
    """
    e = hirzebruch_index(surface)
    lat = surface.lattice
    r = int(r)
    c1 = tuple(int(x) for x in c1)
    eps = Fraction(c1[1], r)                      # paper epsilon (E-coeff of nu)
    phi = Fraction(c1[0], r)                      # paper phi (F-coeff of nu)
    delta_v = _HALF * (1 - Fraction(1, r * r))    # forced exceptional Delta
    ch2 = r * (_HALF * lat.self_pairing((phi, eps)) - delta_v)
    v = (r, c1, ch2)

    def _ceil(x: Fraction) -> int:
        return -((-x.numerator) // x.denominator)

    def tw(w, D):
        rr, cc, hh = w
        return (rr, (cc[0] + rr * D[0], cc[1] + rr * D[1]),
                hh + lat.pairing(cc, D) + rr * _HALF * lat.self_pairing(D))

    def chi_of(w) -> int:
        rr, cc, hh = w
        K = (-(e + 2), -2)
        val = rr - Fraction(lat.pairing(cc, K), 2) + hh
        if val.denominator != 1:
            raise ValueError(f"non-integral chi for {w!r}")
        return int(val)

    if eps.denominator == 1:
        return ()                                 # rem-epInteger: no restriction

    ce = _ceil(eps)
    psi = phi + _HALF * e * (ce - eps) - delta_v / (1 - (ce - eps))
    L0 = (_ceil(psi), ce)                         # pkg (F-coeff, E-coeff)
    alpha = -chi_of(tw(v, (-L0[0] - 1, -L0[1] - 1)))
    beta = -chi_of(tw(v, (-L0[0], -L0[1] - 1)))
    gamma = -chi_of(tw(v, (-L0[0] - 1, -L0[1])))
    dlt = chi_of(tw(v, (-L0[0], -L0[1])))
    if alpha < 0:
        # lem-alpha<0: the general sheaf is a rigid line-bundle sum; the star
        # machinery does not apply -- return no conditions (honest skip).
        return ()
    if -alpha + beta + gamma + dlt != r:          # resolution rank tripwire
        raise AssertionError("Gaeta exponents violate the rank identity")

    rows = []
    kf, ks = e + 2, 2
    for df in range(kf + 1):
        for ds in range(ks + 1):
            if (df, ds) in ((0, 0), (kf, ks)):
                continue
            Dp = (df, ds)
            h2_K = general_betti(*tw(v, (-L0[0] + e + 1 - Dp[0], -L0[1] + 1 - Dp[1])),
                                 surface)[2]
            h2_b = general_betti(*tw(v, (-L0[0] + e - Dp[0], -L0[1] + 1 - Dp[1])),
                                 surface)[2]
            h2_g = general_betti(*tw(v, (-L0[0] + 1 - Dp[0], -L0[1] - Dp[1])),
                                 surface)[2]
            h2_d = general_betti(*tw(v, (-L0[0] - Dp[0], -L0[1] - Dp[1])),
                                 surface)[2]
            rows.append((Dp, beta * h2_b + gamma * h2_g + dlt * h2_d,
                         alpha * h2_K))
    return tuple(rows)


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

    # Condition 0 (E15-M1d, cheapest first): the simple-prioritary chi-box.
    for (D, x) in chi_box_conditions(r, c1_i, surface):
        if x > 0:
            return ExceptionalRefutation(
                True, f"chi(v, v(-D)) = {x} > 0 at D = {D}: a simple bundle has "
                      "hom = ext^2 = 0 there (lem-simple both ways), forcing "
                      "chi <= 0",
                0, False, (), r, c1_i, e, _cert("refuted by the chi-box (M1d)"))

    # Condition 0.5 (E15-M1e): the Gaeta dimension inequality at every box
    # divisor (generic-D-prioritariness beyond the H_n ray).
    for (Dp, lhs, rhs) in gaeta_star_conditions(r, c1_i, surface):
        if lhs > rhs:
            return ExceptionalRefutation(
                True, f"Gaeta star inequality violated at D' = {Dp}: "
                      f"{lhs} > {rhs} -- the general sheaf is not D'-prioritary, "
                      "but a rigid exceptional bundle would be (prop-excPrior)",
                0, False, (), r, c1_i, e, _cert("refuted by the Gaeta star (M1e)"))

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
