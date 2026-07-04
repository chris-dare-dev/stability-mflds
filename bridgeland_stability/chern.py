"""Chern characters on a polarized surface, and the discriminant conventions.

CANONICAL DISCRIMINANT CONVENTION (Coskun-Huizenga / Drezet-Le Potier)
----------------------------------------------------------------------
For a torsion-free sheaf E on a polarized surface (X, H) with d = H^2, we
record only the H-numerical data of the Chern character:

    r = ch0(E)            rank (integer)
    c = ch1(E) . H        H-degree
    e = ch2(E)            the number  int_X ch2

and define

    mu(E)    = c / (r * d)                      (Mumford H-slope)
    Delta(E) = (1/2) * mu^2 - e / (r * d)       (normalized discriminant)

This is the convention used throughout the modern literature (Coskun-Huizenga
surveys, Drezet-Le Potier 1985).  In it:

  * line bundles O(n) have Delta = 0,
  * an exceptional bundle of rank r has  Delta = (1/2)(1 - 1/r^2) < 1/2,
  * the DLP curve takes values in [1/2, 1],
  * Riemann-Roch reads chi(E,F) = r_E r_F (P(mu_F - mu_E) - Delta_E - Delta_F)
    with P(m) = (1/2)(m^2 + 3m + 2)  (see ``exceptional.chi`` and ``dlp.P``).

NOTE.  The original project brief used  Delta_brief = mu^2 - 2 e/(r d) = 2*Delta
(twice this one).  We deliberately use the Coskun-Huizenga normalization because
every DLP / wall / Bogomolov formula in the literature is stated in it.  Use
``ChernChar.discriminant_brief`` if the doubled convention is required.

References
----------
I. Coskun, J. Huizenga, "The birational geometry of the moduli spaces of
    sheaves on P^2", section 4 (definitions 4.1, logarithmic invariants).
J.-M. Drezet, J. Le Potier, "Fibres stables et fibres exceptionnels sur P_2",
    Ann. Sci. ENS 18 (1985), 193-243.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Tuple, Union

Number = Union[int, Fraction]


def Q(x: Number) -> Fraction:
    """Coerce ``int``/``Fraction`` (or a string like ``'-1/2'``) to ``Fraction``."""
    return x if isinstance(x, Fraction) else Fraction(x)


@dataclass(frozen=True)
class ChernChar:
    """H-numerical Chern character ``(r, c, e)`` of a sheaf on a surface.

    ``r`` is the rank (``ch0``), ``c = ch1 . H`` and ``e = ch2``.  Slope and
    discriminant require the surface's ``d = H^2`` because all H-normalized
    formulas carry that factor (a classic P^2-to-general-surface porting bug).
    """

    r: int
    c: Fraction
    e: Fraction

    def __post_init__(self) -> None:
        object.__setattr__(self, "c", Q(self.c))
        object.__setattr__(self, "e", Q(self.e))

    # -- invariants ---------------------------------------------------------
    def slope(self, d: int) -> Fraction:
        if self.r == 0:
            raise ValueError("Mumford slope is undefined for a rank-0 (torsion) object")
        return Fraction(self.c, self.r * d)

    def discriminant(self, d: int) -> Fraction:
        """Normalized discriminant  Delta = (1/2) mu^2 - e/(r d)  (CH convention)."""
        mu = self.slope(d)
        return Fraction(1, 2) * mu * mu - Fraction(self.e, self.r * d)

    def discriminant_brief(self, d: int) -> Fraction:
        """Doubled discriminant  mu^2 - 2 e/(r d)  used in the original brief."""
        return 2 * self.discriminant(d)

    def bogomolov_discriminant(self, d: int) -> Fraction:
        """Integer-flavored Bogomolov discriminant  (c^2/d - 2 r e) = 2 r^2 d Delta.

        Equals the classical ``ch1^2 - 2 ch0 ch2`` when Pic has rank 1 and
        ch1 = (c/d) H.  Non-negativity of this is equivalent to BG.
        """
        return Fraction(self.c * self.c, d) - 2 * self.r * self.e

    # -- twisted characters -------------------------------------------------
    def twist(self, s: Number, d: int) -> "ChernChar":
        """The ``B = sH`` twisted Chern character ``ch^s = ch . e^{-sH}``.

        ch^s_0 = r,  ch^s_1.H = c - s r d,  ch^s_2 = e - s c + (s^2/2) r d.
        The B-twist preserves the discriminant: ``twist(s).discriminant == discriminant``.
        """
        s = Q(s)
        return ChernChar(
            self.r,
            self.c - s * self.r * d,
            self.e - s * self.c + (s * s / 2) * self.r * d,
        )

    # -- Bridgeland central charge -----------------------------------------
    def central_charge(self, s: Number, t: Number, d: int) -> Tuple[float, float]:
        """The Bridgeland central charge ``Z_{s,t}(E) = -int e^{-(s+it)H} ch(E)``.

        Returns ``(Re Z, Im Z)``.  Accepts rational or float ``s, t``; the
        result is returned as ``float`` for plotting convenience.
        ``Re Z = -(e - s c + (s^2 - t^2)/2 r d)``,  ``Im Z = t (c - s r d)``.
        """
        sf, tf, df = float(s), float(t), float(d)
        rf, cf, ef = float(self.r), float(self.c), float(self.e)
        re = -(ef - sf * cf + (sf * sf - tf * tf) / 2 * rf * df)
        im = tf * (cf - sf * rf * df)
        return re, im

    # -- presentation (LaTeX) ----------------------------------------------
    def to_latex(self, name: str = "v") -> str:
        """LaTeX for this class, e.g. ``v=(r,c,\\mathrm{ch}_2)=(2,0,-\\frac{1}{4})``.

        Uses ``\\frac`` (not ``\\tfrac``) so it renders in matplotlib mathtext as
        well as plotly's MathJax.  No surrounding ``$``; the caller wraps it.
        Purely cosmetic -- changes no computed value.
        """
        from ._latex import latex_chern

        return latex_chern(self.r, self.c, self.e, name=name)

    def __format__(self, spec: str) -> str:  # pragma: no cover - cosmetic
        """``format(v, 'latex')`` -> :meth:`to_latex`; any other spec is the default."""
        if spec.lower() in {"latex", "l"}:
            return self.to_latex()
        return format(str(self), spec)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"ChernChar(r={self.r}, c={self.c}, e={self.e})"
