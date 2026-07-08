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

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional, Tuple, Union

from .chern import ChernChar, Number, Q
from .rigor import Rigor, Certificate
from .mukai import MukaiVector, self_pairing
from .threefold import ThreefoldChern, check_bg_threefold  # noqa: F401  (re-export)
from .varieties import Surface, require_faithful_computation


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


# --------------------------------------------------------------------------
# E3-M1 / G6: K3 and abelian moduli non-emptiness (Yoshioka; Yanagida-Yoshioka)
# with an integrality guard.  NO new mathematics: these are the PROVEN existence
# inequalities re-verified in exact arithmetic; no value in check_bg_surface or
# any core formula changes.  The K3 verdict uses the ch2 -> ch2 + r Mukai shift;
# the abelian verdict uses the BARE (unshifted) Mukai vector.
# --------------------------------------------------------------------------


def delta_min(r: int, d: int) -> Fraction:
    """Yoshioka rank-dependent K3 existence threshold on Delta: (1 - 1/r^2)/d.

    A K3 Mukai class has v^2 >= -2 (Yoshioka nonemptiness) iff Delta >= delta_min(r, d);
    the two are equivalent because v^2 = 2 r^2 d (Delta - delta_min) - 2 (exact, re-derived).
    delta_min(1, d) = 0 (the Bogomolov boundary); it rises toward 1/d as r grows -- the sharp
    Mukai/O'Grady moduli-dimension bound, stronger than Bogomolov Delta >= 0 for r >= 2.
    Reference: Bayer-Macri arXiv:1301.6968 (Thm 2.15); Yoshioka moduli on K3.
    """
    if r <= 0:
        raise ValueError("delta_min is defined for positive rank r")
    if d <= 0:
        raise ValueError("delta_min requires positive degree d = H^2")
    return Fraction(r * r - 1, r * r * d)


@dataclass
class ExistenceResult:
    nonempty: Optional[bool]        # True/False when v primitive; None if imprimitive (not asserted)
    satisfies_bound: bool           # the raw inequality v^2 >= bound
    v_squared: int                  # Mukai self-pairing (K3: shifted; abelian: bare)
    mukai: Tuple[int, int, int]     # the Mukai vector used
    discriminant: Fraction          # CH Delta
    delta_min: Fraction             # Delta-threshold for the verdict (K3: (1-1/r^2)/d; abelian: 0)
    moduli_dim: int                 # v^2 + 2 (expected dim of M_H(v))
    primitive: bool
    certificate: Certificate
    note: str = ""


def _lattice_l(v: ChernChar, d: int, fn: str) -> int:
    """Assert v lies in the Picard/Mukai lattice, then return l = c/d. Raises otherwise.
    The c % d guard MUST precede the floor division (danger zone: silent floor -> wrong v^2)."""
    if v.c % d != 0:
        raise ValueError(
            f"{fn}: class (r={v.r}, c={v.c}, e={v.e}) is not in the Picard lattice on d={d}: "
            f"l = c/d = {v.c}/{d} is not an integer, so it has no integral Mukai vector "
            "(refusing to floor c//d and silently return a wrong v^2)."
        )
    if v.e.denominator != 1:
        raise ValueError(
            f"{fn}: class (r={v.r}, c={v.c}, e={v.e}) is not in the Mukai lattice: "
            f"ch2 = {v.e} is not an integer."
        )
    return int(v.c // d)


def _is_primitive(mv: MukaiVector) -> bool:
    return math.gcd(math.gcd(abs(mv.r), abs(mv.l)), abs(mv.s)) == 1


def check_existence_k3(v: ChernChar, surface: Surface) -> ExistenceResult:
    """Yoshioka non-emptiness for a Picard-rank-1 K3: M_H(v) != empty iff v^2 >= -2
    (equivalently Delta >= (1 - 1/r^2)/d). Uses the K3 ch2 -> ch2 + r Mukai shift.
    Raises if surface.kind != 'K3' or v is not in the Mukai lattice (c % d == 0, ch2 in Z)."""
    if surface.kind != "K3":
        raise ValueError(
            "check_existence_k3 requires a K3 surface (surface.kind == 'K3'); got "
            f"kind={surface.kind!r}. Use check_existence_abelian for abelian surfaces "
            "(the ch2 -> ch2 + r shift is K3-only)."
        )
    require_faithful_computation(surface)  # G12 guard (defensive; K3 rows are faithful)
    d = surface.d
    l = _lattice_l(v, d, "check_existence_k3")
    mv = MukaiVector.from_chern(int(v.r), l, int(v.e))     # (r, l, ch2 + r) -- K3 shift
    v2 = self_pairing(mv, d)
    dm = delta_min(int(v.r), d)
    Dm = v.discriminant(d)
    prim = _is_primitive(mv)
    satisfies = v2 >= -2
    nonempty = satisfies if prim else None
    cert = Certificate(
        Rigor.PROVEN,
        ("v primitive (Yoshioka nonemptiness requires a primitive Mukai vector)",
         "generic polarization",
         "Picard rank 1 K3 Mukai lattice"),
        ("arXiv:1301.6968",),
        "M_H(v) nonempty iff v^2 >= -2 (Yoshioka); moduli dim v^2 + 2 (Bayer-Macri Thm 2.15).",
    )
    note = (
        f"v^2 = {v2} {'>=' if satisfies else '<'} -2 (Delta = {Dm} "
        f"{'>=' if Dm >= dm else '<'} delta_min = {dm})"
    )
    if not prim:
        note += "; v is IMPRIMITIVE -- Yoshioka's biconditional needs a primitive Mukai vector, "\
                "so non-emptiness is NOT asserted (nonempty=None)."
    return ExistenceResult(nonempty, satisfies, v2, tuple(mv), Dm, dm, v2 + 2, prim, cert, note)


def check_existence_abelian(v: ChernChar, surface: Surface) -> ExistenceResult:
    """Yanagida-Yoshioka non-emptiness for a Picard-rank-1 abelian surface (arXiv:1203.0884):
    M_H(v) != empty iff v^2 >= 0 on the BARE Mukai vector (r, l, ch2) -- NO K3 +r shift.
    Raises if surface.kind != 'abelian' or v is not in the lattice."""
    if surface.kind != "abelian":
        raise ValueError(
            "check_existence_abelian requires an abelian surface (surface.kind == 'abelian'); got "
            f"kind={surface.kind!r}. Use check_existence_k3 for K3 surfaces. The K3 ch2 -> ch2 + r "
            "shift must NOT be applied to abelian input (it injects a spurious +2/d)."
        )
    require_faithful_computation(surface)  # G12 guard (defensive; abelian rows are faithful)
    d = surface.d
    l = _lattice_l(v, d, "check_existence_abelian")
    mv = MukaiVector(int(v.r), l, int(v.e))               # BARE -- no +r shift
    v2 = self_pairing(mv, d)
    Dm = v.discriminant(d)
    prim = _is_primitive(mv)
    satisfies = v2 >= 0
    nonempty = satisfies if prim else None
    cert = Certificate(
        Rigor.PROVEN,
        ("v primitive (Yanagida-Yoshioka nonemptiness requires a primitive Mukai vector)",
         "generic polarization",
         "Picard rank 1 abelian surface"),
        ("arXiv:1203.0884",),
        "M_H(v) nonempty iff v^2 >= 0 (Yanagida-Yoshioka); bare Mukai vector (no +r shift).",
    )
    note = f"v^2 = {v2} {'>=' if satisfies else '<'} 0 (bare Mukai self-pairing)"
    if not prim:
        note += "; v is IMPRIMITIVE -- non-emptiness NOT asserted (nonempty=None)."
    return ExistenceResult(nonempty, satisfies, v2, tuple(mv), Dm, Fraction(0), v2 + 2, prim, cert, note)
