"""Algorithm 2 - the Drezet-Le Potier curve delta(mu) on P^2 (CORRECTED).

The moduli space M(r, c1, ch2) of Gieseker-semistable sheaves on P^2 is
non-empty iff the integrality conditions hold (c1 = r*mu in Z and
chi = r(P(mu) - Delta) in Z) and either

    Delta >= delta(mu)            (positive-dimensional moduli, ABOVE the curve)

or (r, c1, ch2) is the Chern character of an exceptional bundle (an isolated
point STRICTLY BELOW the curve).

The DLP curve itself is the fractal upper envelope (Coskun-Huizenga survey,
Theorem 4.15 / Figure 1):

    delta(mu) = sup over exceptional slopes alpha with |mu - alpha| < 3 of
                ( P(-|mu - alpha|) - Delta_alpha ),     clamped below by 1/2,

with P(m) = (1/2)(m^2 + 3m + 2).  Each exceptional bundle contributes an
upward CUSP of height ``1 - Delta_alpha`` at mu = alpha; between cusps the
curve dips to 1/2.  The control interval where alpha is the binding bundle is
``I_alpha = (alpha - x_alpha, alpha + x_alpha)`` with

    x_alpha = (3 - sqrt(5 + 8 Delta_alpha)) / 2.

WHY THE BRIEF'S ALGORITHM 2 IS WRONG (three independent errors):
  (1) it places exceptional bundles ON the curve (nu(alpha) = Delta_alpha);
      they actually sit at Delta_alpha which is BELOW delta(alpha)=1-Delta_alpha;
  (2) it models each arc as one parabola through a Farey "mediant"; the true
      local shape is a two-branch cusp of the single bundle's own parabola;
  (3) it uses the non-existent rank-3 mediant (see exceptional.py).
Consequently the brief's nu(1/2)=3/4, nu(1/3)=8/9 are wrong; the true CH
values are delta(1/2)=5/8 and delta(1/3)=5/9.  See docs/CORRECTIONS.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import List, Optional, Tuple

from .chern import ChernChar, Number, Q
from .exceptional import Bundle, P, chi, enumerate_exceptional, is_exceptional

# Re-export so callers can do ``from bridgeland_stability.dlp import P``.
__all__ = [
    "P",
    "delta",
    "control_interval_halfwidth",
    "DLPCurve",
    "dlp_curve",
    "moduli_nonempty",
]


def delta(mu: Number, bundles: List[Bundle]) -> Fraction:
    """Exact DLP curve value at ``mu`` from a list of exceptional ``bundles``.

    ``delta(mu) = max(1/2, sup_{|mu-alpha|<3} (P(-|mu-alpha|) - Delta_alpha))``.
    The result is exact (``Fraction``).  Pass enough bundles (rank up to some
    ``R_max`` over a window containing ``[mu-3, mu+3]``) to resolve the curve;
    unresolved thin gaps default to the correct limiting value ``1/2``.
    """
    mu = Q(mu)
    best = Fraction(1, 2)
    for b in bundles:
        dist = abs(mu - b.slope)
        if dist < 3:
            val = P(-dist) - b.discriminant
            if val > best:
                best = val
    return best


def control_interval_halfwidth(b: Bundle) -> float:
    """``x_alpha = (3 - sqrt(5 + 8 Delta_alpha)) / 2`` (float), the half-width of I_alpha."""
    da = float(b.discriminant)
    return (3.0 - math.sqrt(5.0 + 8.0 * da)) / 2.0


@dataclass
class DLPCurve:
    """A sampled Drezet-Le Potier curve over ``[mu_min, mu_max]`` at resolution ``R_max``.

    Attributes
    ----------
    mus, deltas : list[Fraction]
        Sampled curve points (exact).  ``deltas[i] = delta(mus[i])``.
    bundles : list[Bundle]
        The exceptional bundles used (rank <= R_max), sorted by slope.
    cusps : list[tuple[Fraction, Fraction]]
        ``(alpha, 1 - Delta_alpha)`` -- the cusp tips that lie ON the curve.
    exceptional_points : list[tuple[Fraction, Fraction]]
        ``(alpha, Delta_alpha)`` -- the exceptional bundles, BELOW the curve.
    R_max : int
    """

    mus: List[Fraction]
    deltas: List[Fraction]
    bundles: List[Bundle]
    cusps: List[Tuple[Fraction, Fraction]]
    exceptional_points: List[Tuple[Fraction, Fraction]]
    R_max: int

    def as_floats(self) -> Tuple[List[float], List[float]]:
        return [float(m) for m in self.mus], [float(d) for d in self.deltas]


def dlp_curve(
    mu_min: Number = 0,
    mu_max: Number = 1,
    R_max: int = 50,
    samples_per_unit: int = 400,
) -> DLPCurve:
    """Build a sampled DLP curve, injecting cusp points + interval endpoints for crispness."""
    mu_min, mu_max = Q(mu_min), Q(mu_max)
    # enumerate exceptional bundles over a window padded by 3 (the |mu-alpha|<3 reach)
    bundles = enumerate_exceptional(mu_min - 3, mu_max + 3, R_max)

    # build a sorted set of sample abscissae
    n = max(2, int(samples_per_unit * float(mu_max - mu_min)) + 1)
    xs = {mu_min + Fraction(i, n) * (mu_max - mu_min) for i in range(n + 1)}
    for b in bundles:
        a = b.slope
        if mu_min <= a <= mu_max:
            xs.add(a)  # cusp tip
            xh = control_interval_halfwidth(b)
            # add interval endpoints (approx, as rationals) to capture the dip to 1/2
            for sgn in (-1, 1):
                e = a + Fraction(sgn) * Fraction(round(xh * 10000), 10000)
                if mu_min <= e <= mu_max:
                    xs.add(e)
    mus = sorted(xs)
    deltas = [delta(m, bundles) for m in mus]

    cusps = []
    excpts = []
    for b in bundles:
        a = b.slope
        if mu_min <= a <= mu_max:
            cusps.append((a, 1 - b.discriminant))
            excpts.append((a, b.discriminant))
    cusps.sort()
    excpts.sort()
    return DLPCurve(mus, deltas, bundles, cusps, excpts, R_max)


def moduli_nonempty(
    E: Bundle, R_max: int = 50
) -> dict:
    """Decide whether the P^2 moduli space M(E) is non-empty, with full reasoning.

    Returns a dict with keys: ``mu``, ``discriminant`` (CH), ``delta``,
    ``exceptional`` (bool), ``positive_dimensional`` (Delta >= delta),
    ``integral`` (c1, chi integral), ``nonempty`` (bool), ``moduli_dim``
    (= r^2(2*Delta_brief - 1) + 1 when positive-dimensional, else 0), ``reason``.
    """
    mu = E.slope
    Dm = E.discriminant  # CH
    exceptional = is_exceptional(E)
    # enumerate around mu
    bundles = enumerate_exceptional(mu - 3, mu + 3, R_max)
    d_val = delta(mu, bundles)
    positive_dim = Dm >= d_val
    # integrality: c1 in Z (given) and chi(E) = r(P(mu) - Delta) in Z
    chi_E = E.r * (P(mu) - Dm)
    integral = (E.c2.denominator == 1) and (chi_E.denominator == 1)
    nonempty = integral and (positive_dim or exceptional)
    # Drezet-Le Potier dimension r^2(2 Delta - 1) + 1; Delta here in the brief's
    # doubled convention (= 2*CH) so dim = r^2(2*Delta_brief - 1)+1.
    dim = E.r * E.r * (2 * E.discriminant_brief - 1) + 1 if positive_dim else 0
    if exceptional and not positive_dim:
        reason = "exceptional bundle: isolated point below the DLP curve"
    elif positive_dim and integral:
        reason = f"Delta={Dm} >= delta({mu})={d_val}: positive-dimensional moduli"
    elif not integral:
        reason = "non-integral Chern data (c2 or chi not an integer)"
    else:
        reason = f"Delta={Dm} < delta({mu})={d_val} and not exceptional: EMPTY"
    return {
        "mu": mu,
        "discriminant": Dm,
        "delta": d_val,
        "exceptional": exceptional,
        "positive_dimensional": positive_dim,
        "integral": integral,
        "nonempty": nonempty,
        "moduli_dim": dim,
        "reason": reason,
    }
