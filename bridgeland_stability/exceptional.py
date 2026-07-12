"""Algorithm 1 - exceptional bundles on P^2 (CORRECTED).

An exceptional bundle on P^2 is a stable bundle E with Ext^1(E,E) = 0
(equivalently chi(E,E) = 1).  By Drezet-Le Potier / Rudakov:

  * the RANKS of exceptional bundles are exactly the MARKOV NUMBERS
    {1, 2, 5, 13, 29, 34, 89, 169, ...} (solutions feeding x^2+y^2+z^2=3xyz),
  * their SLOPES are the "Markov fractions": a slope p/q (lowest terms) is
    exceptional iff q is a Markov number, and then the bundle has rank r = q,
    c1 = p, and  ch2 = (p^2 - q^2 + 1)/(2q),  Delta = (1/2)(1 - 1/r^2).

The slopes are produced by the Drezet-Le Potier  epsilon : Z[1/2] -> slopes
recursion (Coskun-Huizenga survey, section 4.2):

    eps(n) = n,   eps((2p+1)/2^{q+1}) = eps(p/2^q) . eps((p+1)/2^q)

where the binary operation on two adjacent exceptional slopes alpha, beta is

    alpha . beta = (alpha + beta)/2 + (Delta_beta - Delta_alpha)/(3 + alpha - beta).

WHY THE BRIEF'S ALGORITHM 1 IS WRONG.  The brief used the naive Farey rule
``r_G = r_E + r_F`` (ranks 1,2,3,4,...).  That manufactures a "rank-3 bundle"
(3, 1, -7/6) at slope 1/3 -- but its c2 = ch1^2/2 - ch2 = 1/2 + 7/6 = 5/3 is
not an integer, so no such bundle exists.  3 is not a Markov number; the true
exceptional bundle between O (slope 0) and T(-1) (slope 1/2) has rank 5 at
slope 2/5.  See docs/CORRECTIONS.md.

All arithmetic is exact (``fractions.Fraction``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Set

from .chern import ChernChar, Number, Q


def P(m: Number) -> Fraction:
    """Hilbert/Euler polynomial of O_{P^2}:  P(m) = (1/2)(m^2 + 3m + 2) = chi(O,O(m))."""
    m = Q(m)
    return Fraction(1, 2) * (m * m + 3 * m + 2)


@dataclass(frozen=True)
class Bundle:
    """An exceptional bundle on P^2, recorded as ``(r, c1, ch2)``.

    ``c1`` is an integer (P^2 has Pic = Z.H), ``ch2`` is a half-integer.
    """

    r: int
    c1: int
    ch2: Fraction

    def __post_init__(self) -> None:
        object.__setattr__(self, "ch2", Q(self.ch2))

    @property
    def slope(self) -> Fraction:
        return Fraction(self.c1, self.r)

    @property
    def discriminant(self) -> Fraction:
        """CH-normalized discriminant on P^2 (d = 1):  (1/2) mu^2 - ch2/r."""
        mu = self.slope
        return Fraction(1, 2) * mu * mu - Fraction(self.ch2, self.r)

    @property
    def discriminant_brief(self) -> Fraction:
        return 2 * self.discriminant

    @property
    def c2(self) -> Fraction:
        """Second Chern class  c2 = ch1^2/2 - ch2  (must be an integer for a bundle)."""
        return Fraction(self.c1 * self.c1, 2) - self.ch2

    def chern_char(self) -> ChernChar:
        return ChernChar(self.r, Fraction(self.c1), self.ch2)

    # -- constructors -------------------------------------------------------
    @classmethod
    def O(cls, n: int) -> "Bundle":
        """The line bundle O_{P^2}(n):  rank 1, c1 = n, ch2 = n^2/2,  Delta = 0."""
        return cls(1, n, Fraction(n * n, 2))

    @classmethod
    def T_minus1(cls) -> "Bundle":
        """The exceptional bundle T_{P^2}(-1):  rank 2, c1 = 1, ch2 = -1/2,  Delta = 3/8."""
        return cls(2, 1, Fraction(-1, 2))

    @classmethod
    def from_slope(cls, alpha: Number) -> "Bundle":
        """The exceptional bundle of slope ``alpha = p/q`` (lowest terms): rank q.

        Uses  r = q,  c1 = p,  ch2 = (p^2 - q^2 + 1)/(2q).  This is a genuine
        bundle only when ``q`` is a Markov number; otherwise ``is_exceptional``
        on the result is ``False`` (its c2 is not an integer).
        """
        a = Q(alpha)
        p, q = a.numerator, a.denominator
        ch2 = Fraction(p * p - q * q + 1, 2 * q)
        return cls(q, p, ch2)

    # -- presentation (LaTeX) ----------------------------------------------
    def to_latex(self) -> str:
        """LaTeX tuple ``(r,c_1,\\mathrm{ch}_2)=(5,2,-\\frac{1}{5})``.  No ``$``."""
        from ._latex import latex_chern

        return latex_chern(
            self.r, self.c1, self.ch2, name="", components=("r", "c_1", r"\mathrm{ch}_2")
        )

    def slope_latex(self) -> str:
        """LaTeX for the slope ``\\mu = c_1/r`` as an exact fraction.  No ``$``."""
        from ._latex import latex_frac

        return latex_frac(self.slope)

    def __format__(self, spec: str) -> str:  # pragma: no cover - cosmetic
        if spec.lower() in {"latex", "l"}:
            return self.to_latex()
        return format(str(self), spec)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"Bundle(r={self.r}, c1={self.c1}, ch2={self.ch2})"


# --------------------------------------------------------------------------
# Riemann-Roch on P^2 (independent of the exceptionality shortcut)
# --------------------------------------------------------------------------
def chi(E: Bundle, F: Bundle) -> Fraction:
    """Euler pairing ``chi(E,F) = sum (-1)^i ext^i(E,F)`` on P^2 via Riemann-Roch.

    Computed as the degree-2 part of  ch(E^v) . ch(F) . td(P^2),
    with td(P^2) = (1, (3/2)H, H^2).  This is the ground-truth pairing; for an
    exceptional bundle it returns ``chi(E,E) = 1``.
    """
    # ch(E^v) = (rE, -c1E, ch2E);  ch(F) = (rF, c1F, ch2F)
    deg0 = E.r * F.r
    deg1 = E.r * F.c1 - F.r * E.c1
    deg2 = E.r * F.ch2 + F.r * E.ch2 - E.c1 * F.c1
    # integrate against td(P^2) = (1, 3/2, 1):  deg2*1 + deg1*(3/2) + deg0*1
    return deg2 + Fraction(3, 2) * deg1 + deg0


def is_exceptional(E: Bundle) -> bool:
    """True iff ``E`` is a genuine exceptional bundle: its slope lies in the Drezet-Le
    Potier epsilon-recursion image (Theoreme A), with rank = slope denominator and the
    exceptional ``ch2 = Bundle.from_slope(slope).ch2``.  ``chi(E,E)=1`` and integral
    ``c2`` are a cheap NECESSARY pre-filter only -- they are met by infinitely many
    non-exceptional integral classes (e.g. the Markov-rank impostor ``(610,133,-581/2)``,
    whose rank 610 = 2*5*61 IS a Markov number).  See docs/CORRECTIONS.md section 8.
    """
    if chi(E, E) != 1 or E.c2.denominator != 1:      # cheap necessary pre-filter
        return False
    mu = E.slope
    if mu.denominator != E.r:                         # rank must equal the reduced denominator
        return False
    if E.ch2 != Bundle.from_slope(mu).ch2:            # the unique exceptional ch2 at this slope
        return False
    return is_exceptional_slope(mu, E.r)              # sufficiency: epsilon-image membership


def is_semiexceptional_p2(r: int, c1: int, ch2: Number) -> bool:
    """True iff (r, c1, ch2) = m * ch(E) for an exceptional bundle E on P^2 and m >= 1.

    Coskun-Huizenga arXiv:1907.06739 Ex.1.14: a semiexceptional bundle is a direct sum of
    copies of an exceptional bundle.  m=1 makes every exceptional bundle semiexceptional.
    The P^2 analogue of dlp_hirzebruch.is_semiexceptional; agrees with the oracle's
    reference_is_semiexceptional (tests/oracle/dlp_reference.py) bit-for-bit.
    """
    ch2 = Q(ch2)
    if r < 1:
        return False
    for m in range(1, r + 1):
        if r % m or c1 % m:
            continue
        if is_exceptional(Bundle(r // m, c1 // m, ch2 / m)):
            return True
    return False


# --------------------------------------------------------------------------
# The Drezet-Le Potier epsilon recursion
# --------------------------------------------------------------------------
def _delta_of_slope(alpha: Fraction) -> Fraction:
    """Discriminant of the exceptional bundle of slope ``alpha`` (rank = denominator)."""
    q = alpha.denominator
    return Fraction(1, 2) * (1 - Fraction(1, q * q))


def exceptional_mediant(alpha: Number, beta: Number) -> Fraction:
    """The DLP operation ``alpha . beta`` producing the exceptional slope between them.

        alpha . beta = (alpha + beta)/2 + (Delta_beta - Delta_alpha)/(3 + alpha - beta).

    For (0, 1) -> 1/2 (rank 2);  for (0, 1/2) -> 2/5 (rank 5);  for (0, 2/5) -> 5/13.
    """
    a, b = Q(alpha), Q(beta)
    da, db = _delta_of_slope(a), _delta_of_slope(b)
    return (a + b) / 2 + (db - da) / (3 + a - b)


def enumerate_exceptional(
    mu_min: Number, mu_max: Number, R_max: int
) -> List[Bundle]:
    """All exceptional bundles with slope in ``[mu_min, mu_max]`` and rank <= ``R_max``.

    Generated by the epsilon recursion: integer slopes are the seeds (line
    bundles O(n)); each interval between adjacent exceptional slopes is
    subdivided by ``exceptional_mediant`` until the child rank exceeds R_max.
    Ranks strictly increase down the tree, so the recursion terminates and the
    pruning is exact.  The result is sorted by slope.
    """
    mu_min, mu_max = Q(mu_min), Q(mu_max)
    if mu_min > mu_max:
        raise ValueError("mu_min must be <= mu_max")
    found: Dict[Fraction, Bundle] = {}

    def add(alpha: Fraction) -> None:
        if mu_min <= alpha <= mu_max and alpha.denominator <= R_max:
            found[alpha] = Bundle.from_slope(alpha)

    lo = math.floor(mu_min) - 1
    hi = math.ceil(mu_max) + 1
    for n in range(lo, hi + 1):
        add(Fraction(n))

    stack = [(Fraction(n), Fraction(n + 1)) for n in range(lo, hi)]
    while stack:
        a, b = stack.pop()
        g = exceptional_mediant(a, b)
        if g.denominator > R_max:
            continue
        add(g)
        stack.append((a, g))
        stack.append((g, b))

    return [found[k] for k in sorted(found)]


def certified_rank_cutoff(mu: Number, R_max: int = 0) -> int:
    """The certified-exact rank cutoff for evaluating ``delta(mu)`` on ``P^2``.

    An exceptional slope ``alpha`` of rank ``rho`` lifts ``delta`` above ``1/2`` only
    inside its control interval of half-width ``x_rho <= 0.382/rho^2``; and if
    ``alpha != mu`` with ``q = denominator(mu)`` then ``|mu - alpha| >= 1/(q rho)``,
    which forces ``rho < 0.382 q``.  The only remaining case is ``alpha = mu``
    (``rho = q``).  So enumerating exceptional bundles up to rank ``denominator(mu)``
    resolves ``delta(mu)`` EXACTLY -- it is a theorem, not a truncation
    (docs/CORRECTIONS.md sec.8).  Returns ``max(R_max, denominator(mu))`` so a caller's
    larger explicit budget is still honoured.

    This is the SINGLE source of the cutoff shared by both P^2 decision procedures
    (:func:`bridgeland_stability.dlp.moduli_nonempty` and
    :func:`bridgeland_stability.nonemptiness_rational.delta_H`); they must agree by
    construction, and defect A4b was exactly a cutoff that landed in one and not the
    other.  Do not re-inline ``max(R_max, mu.denominator)`` at a call site.
    """
    return max(int(R_max), Q(mu).denominator)


def exceptional_slopes(mu: Number, r_max: int, reach: Number = 3) -> List[Fraction]:
    """The exceptional (DLP epsilon-recursion image) slopes within ``[mu-reach, mu+reach]``
    of rank (= slope denominator) ``<= r_max``, sorted.  Wraps :func:`enumerate_exceptional`
    -- no new recursion -- so it generates exactly the same epsilon-image.  ``reach=3``
    matches the ``|mu-alpha| < 3`` DLP window used by :func:`bridgeland_stability.dlp.delta`.
    """
    m = Q(mu)
    return [b.slope for b in enumerate_exceptional(m - reach, m + reach, r_max)]


def is_exceptional_slope(alpha: Number, r_max: "int | None" = None) -> bool:
    """``True`` iff ``alpha`` is an exceptional slope: a member of the Drezet-Le Potier
    epsilon-recursion image (Theoreme A).  ``r_max`` defaults to ``denominator(alpha)`` --
    the certified-exact cutoff (docs/CORRECTIONS.md section 8) -- and a slope whose
    denominator exceeds ``r_max`` is rejected outright.  Generated epsilon-slope
    denominators strictly increase down the tree, so a slope of denominator ``q`` is
    produced iff it is genuinely in the image and ``r_max >= q``.
    """
    a = Q(alpha)
    if r_max is None:
        r_max = a.denominator
    if a.denominator > r_max:
        return False
    lo = math.floor(a)
    return any(b.slope == a for b in enumerate_exceptional(lo, lo + 1, r_max))


# --------------------------------------------------------------------------
# Markov numbers (the set of exceptional-bundle ranks) - for labeling / tests
# --------------------------------------------------------------------------
def markov_numbers(max_value: int) -> Set[int]:
    """All Markov numbers <= ``max_value`` (the possible ranks of exceptional bundles).

    Generated from the Markov tree seeded at the triple (1, 1, 2) via the Vieta
    involution ``z -> 3xy - z``.
    """
    result: Set[int] = set()
    if max_value < 1:
        return result
    start = (1, 1, 2)
    stack = [start]
    seen = {start}
    while stack:
        x, y, z = stack.pop()
        for v in (x, y, z):
            if v <= max_value:
                result.add(v)
        # Vieta jumps on each coordinate; prune on the MAX entry (not the min):
        # a triple whose largest entry exceeds max_value can only lead to larger
        # triples, so it contributes no new Markov number <= max_value.
        for nx, ny, nz in ((3 * y * z - x, y, z), (x, 3 * x * z - y, z), (x, y, 3 * x * y - z)):
            if min(nx, ny, nz) < 1:
                continue
            if max(nx, ny, nz) > max_value:
                continue
            key = tuple(sorted((nx, ny, nz)))
            if key not in seen:
                seen.add(key)
                stack.append((nx, ny, nz))
    return result


def is_markov_number(n: int) -> bool:
    return n in markov_numbers(n)
