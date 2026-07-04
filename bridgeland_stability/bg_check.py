"""Algorithm 4 - the Bogomolov-Gieseker inequality checker (surfaces).

THEOREM (Bogomolov 1978, Gieseker 1977).  For a torsion-free mu_H-semistable
sheaf E on a smooth projective surface with ample H, ``Delta(E) >= 0`` (with
the CH-normalized ``Delta = (1/2)mu^2 - ch2/(rd)``; equivalently the integer
``c1^2/d - 2 r ch2 >= 0``).  Equality holds iff E is projectively flat.

This is a *necessary* condition for mu-semistability; ``check_bg_surface`` only
asserts the inequality (callers must supply a mu-semistable E -- see the gotcha
about O(1)+O(-1), which is NOT mu-semistable and for which BG does not apply).

The threefold tilt-BG checker lives in :mod:`bridgeland_stability.threefold`
(it requires a tilt slope and the Bayer-Macri-Toda form); it is re-exported
here as :func:`check_bg_threefold` for convenience.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Union

from .chern import ChernChar, Number, Q
from .threefold import ThreefoldChern, check_bg_threefold  # noqa: F401  (re-export)
from .varieties import Surface


@dataclass
class BGResult:
    satisfies: bool
    discriminant: Fraction  # CH convention
    discriminant_brief: Fraction  # doubled (brief) convention
    bogomolov_integer: Fraction  # c1^2/d - 2 r ch2
    equality: bool
    note: str


def check_bg_surface(E: Union[ChernChar, "BundleLike"], surface: Surface) -> BGResult:
    """Check BG for a surface sheaf: returns ``Delta`` (both conventions) and verdict.

    Accepts a :class:`~bridgeland_stability.chern.ChernChar` or anything with a
    ``chern_char()`` method (e.g. an exceptional :class:`Bundle`).
    """
    ch = E.chern_char() if hasattr(E, "chern_char") else E
    if not isinstance(ch, ChernChar):
        raise TypeError("expected a ChernChar or an object with .chern_char()")
    d = surface.d
    Dm = ch.discriminant(d)
    return BGResult(
        satisfies=Dm >= 0,
        discriminant=Dm,
        discriminant_brief=2 * Dm,
        bogomolov_integer=ch.bogomolov_discriminant(d),
        equality=(Dm == 0),
        note=(
            "BG holds"
            if Dm > 0
            else "BG saturated (projectively flat)"
            if Dm == 0
            else "BG VIOLATED (not the ch of a mu-semistable sheaf)"
        ),
    )


# Typing shim for the docstring/annotation above.
BundleLike = ChernChar
