"""Dependency-free LaTeX-string helpers (mathtext- and MathJax-compatible).

Used by :meth:`ChernChar.to_latex` / :meth:`Bundle.to_latex` and by the
visualization style layer.  This module has **no third-party dependencies** so
that ``import bridgeland_stability`` never pulls in matplotlib/plotly.

Conventions
-----------
* Fractions are emitted with ``\\frac`` (which renders in *both* matplotlib
  mathtext and plotly's MathJax) -- never ``\\tfrac``/``\\dfrac`` (mathtext does
  not implement them).
* Returned strings carry **no** surrounding ``$`` delimiters; the caller wraps
  them (matplotlib wants ``$...$``; plotly wants ``$...$`` in labels too).
"""

from __future__ import annotations

from fractions import Fraction
from typing import Union

Number = Union[int, Fraction]


def latex_frac(x: Number) -> str:
    """LaTeX for an exact rational ``x`` (``int`` or :class:`~fractions.Fraction`).

    Integers render bare (``-3``); proper fractions render as ``\\frac{a}{b}``
    with the sign pulled out front (``-\\frac{1}{4}``).  No ``$`` delimiters.
    """
    fr = x if isinstance(x, Fraction) else Fraction(x)
    if fr.denominator == 1:
        return str(fr.numerator)
    sign = "-" if fr < 0 else ""
    return f"{sign}\\frac{{{abs(fr.numerator)}}}{{{fr.denominator}}}"


def latex_frac_from_float(value: float, max_denominator: int = 64) -> str:
    """LaTeX for the simplest rational approximating ``value`` (for tick labels).

    Uses :meth:`Fraction.limit_denominator`; exact for values that *are* simple
    rationals (tick positions like ``0.4`` -> ``\\frac{2}{5}``), and a tidy
    nearby fraction otherwise.  No ``$`` delimiters.
    """
    return latex_frac(Fraction(value).limit_denominator(max_denominator))


def latex_sqrt(x: Number) -> str:
    """LaTeX ``\\sqrt{...}`` of an exact rational, simplifying perfect squares.

    ``9/4 -> \\frac{3}{2}``; ``3 -> \\sqrt{3}``; ``1/2 -> \\frac{\\sqrt{2}}{2}``.
    Falls back to ``\\sqrt{\\frac{a}{b}}`` when no tidy form is available.
    """
    fr = x if isinstance(x, Fraction) else Fraction(x)
    if fr < 0:
        return f"\\sqrt{{{latex_frac(fr)}}}"
    n, d = fr.numerator, fr.denominator

    def _isqrt_exact(m: int):
        r = int(m ** 0.5)
        for cand in (r - 1, r, r + 1):
            if cand >= 0 and cand * cand == m:
                return cand
        return None

    rn, rd = _isqrt_exact(n), _isqrt_exact(d)
    if rn is not None and rd is not None:
        return latex_frac(Fraction(rn, rd))
    if rd is not None and rd != 0:  # sqrt(n)/rd, rationalize not needed
        inner = "" if n == 1 else f"\\sqrt{{{n}}}"
        if n == 1:
            return latex_frac(Fraction(1, rd))
        return f"\\frac{{\\sqrt{{{n}}}}}{{{rd}}}"
    return f"\\sqrt{{{latex_frac(fr)}}}"


def latex_chern(r, c, e, *, name: str = "v", components=("r", "c", r"\mathrm{ch}_2")) -> str:
    """LaTeX tuple for a numerical Chern character, e.g.
    ``v=(r,c,\\mathrm{ch}_2)=(2,0,-\\frac{1}{4})``.  No ``$`` delimiters.
    """
    comps = ",".join(components)
    body = f"({latex_frac(r)},\\,{latex_frac(c)},\\,{latex_frac(e)})"
    head = f"{name}=" if name else ""
    return f"{head}({comps})={body}"
