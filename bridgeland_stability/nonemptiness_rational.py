"""E11-M3 / G18a: the in-architecture ``Delta >= delta_H`` non-emptiness evaluator.

Coskun-Huizenga's non-emptiness theory for moduli of semistable sheaves on a
rational surface ``(X, H)`` splits into two layers:

  1. an **out-of-architecture, sheaf-theoretic core** -- pass to the prioritary
     sheaf stack, compute the Harder-Narasimhan filtration of the *generic*
     prioritary sheaf; semistable sheaves exist **iff that HN filtration has
     length one**.  This is genuinely sheaf-theoretic (not pure Chern
     arithmetic) and is DELEGATED (M2/OSCAR oracle G16 -- E11-M5, a paper table
     E11-M4, or deferred G19).  Nothing in this module computes it.
  2. an **in-architecture numerical verdict** -- *given* the HN-length-one datum,
     the existence question reduces to the sharp exact inequality
     ``Delta(xi) >= delta_H(xi)`` in the polarization ``H``.  THIS is the thin
     slice shipped here: an exact-``Fraction`` comparison that returns yes/no
     with a G5 :class:`~bridgeland_stability.rigor.Certificate`.

Scope (THIN, stated as such).  This module ships ONLY layer (2).  On ``P^2`` the
sharp bound is the closed-form Drezet-Le Potier curve
:func:`bridgeland_stability.dlp.delta` -- a theorem (HN length one is implicit),
so the P^2 verdict is ``PROVEN``.  Off ``P^2`` there is **no closed
delta-curve**: the sharp, *polarization-dependent* ``delta_H`` is exactly the
paper-tabulated / oracle-supplied datum of E11-M4/M5.  Absent a certified HN
datum this module falls back to the **Bogomolov floor** ``delta_H = 0`` (the
sanity floor ``Delta >= 0``) and marks the verdict ``HEURISTIC``.  A verdict is
``PROVEN`` only when the HN-length-one hypothesis came from a certified source
(the P^2 DLP closed form, or a ``PAPER`` / ``ORACLE`` supplier).

The P^2 exceptional disjunct (do NOT drop it).  The Drezet-Le Potier theorem
states that ``M(r, c1, ch2)`` on ``P^2`` is non-empty **iff**
``Delta(xi) >= delta(mu)`` **OR** ``xi`` is the Chern character of an
*exceptional bundle* -- a NON-empty isolated point sitting STRICTLY BELOW the
curve (e.g. ``T_{P^2}(-1) = (2, 1, -1/2)`` has ``Delta = 3/8 < 5/8 = delta(1/2)``
yet its moduli space is a single reduced point).  The inequality
``Delta >= delta_H`` alone therefore under-reports: it would return a PROVEN
``nonempty=False`` for a class that exists.  ``moduli_nonempty`` ORs in the
exceptional disjunct on ``P^2`` (via :func:`_is_p2_exceptional`) so its verdict
matches :func:`bridgeland_stability.dlp.moduli_nonempty` and the DLP theorem
exactly; the verdict's :attr:`NonemptinessVerdict.exceptional` flag records which
disjunct fired.  A ``PROVEN`` cert here certifies moduli NON-emptiness, never
emptiness.

Convention (fixed).  Coskun-Huizenga normalization
``Delta = 1/2 mu^2 - ch2/(r d)`` with ``mu = <c1, H>/(r d)`` and ``d = H^2``
(the module never uses the doubled ``discriminant_brief``).

Exact arithmetic.  Every quantity is ``fractions.Fraction``; no float ever
appears.  Stdlib-only at import time (``fractions``, ``dataclasses``, ``enum``,
``typing`` + the stdlib-only intra-package ``dlp`` / ``exceptional`` /
``exceptional_surface`` / ``varieties`` / ``rigor``), preserving the
zero-runtime-dependency invariant.

References
----------
* Coskun-Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
  arXiv:1907.06739 [PROVEN] -- sharp Bogomolov inequalities ``Delta >= delta_H(c1/r)``
  with semistable sheaves existing iff the prioritary-sheaf HN filtration has length
  one; the sharp polarization-dependent ``delta_H`` is computed by a limiting
  procedure (via exceptional bundles), NOT a simple closed form / table.
  (The "nef cone ... strong Bogomolov inequalities" title is a DIFFERENT paper,
  arXiv:1512.02661 -- do not conflate.)
* Levine-Zhang (Coskun-Huizenga circle), arXiv:1910.14060 [PROVEN] -- del Pezzo
  deg >= 3 existence / anticanonical.
* Coskun-Huizenga weak Brill-Noether for rational surfaces, arXiv:1611.02674
  [PROVEN].
* Drezet-Le Potier, Ann. Sci. ENS 18 (1985) -- the P^2 closed-form
  ``delta``-curve (the ``DLP`` mode).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Optional, Sequence, Union

from .dlp import delta
from .exceptional import Bundle, enumerate_exceptional, is_exceptional
from .exceptional_surface import SurfaceBundle
from .rigor import Certificate, Rigor
from .varieties import Surface, require_faithful_computation

__all__ = [
    "HNMode",
    "NonemptinessVerdict",
    "discriminant_H",
    "delta_H",
    "moduli_nonempty",
]

Number = Union[int, Fraction]


class HNMode(str, Enum):
    """Which supplier certified the HN-length-one hypothesis behind a verdict."""

    DLP = "dlp"            # PROVEN: P^2 Drezet-Le Potier closed form (HN length one implicit)
    PAPER = "paper"        # PROVEN: paper-tabulated HN-length-one datum (E11-M4 hook)
    ORACLE = "oracle"      # PROVEN: M2/OSCAR-constructed HN filtration (E11-M5 hook)
    HEURISTIC = "heuristic"  # HEURISTIC: no certified HN datum (Bogomolov floor)


#: The modes whose HN-length-one hypothesis is certified -> verdict rigor PROVEN.
_CERTIFIED = frozenset({HNMode.DLP, HNMode.PAPER, HNMode.ORACLE})

_MODE_CERT = {
    HNMode.DLP: Certificate(
        Rigor.PROVEN,
        ("HN filtration has length one (implicit): the P^2 Drezet-Le Potier "
         "closed-form delta-curve supplies the sharp bound delta_H = delta(mu)",),
        ("Drezet-Le Potier, Ann. Sci. ENS 18 (1985)", "arXiv:1611.02674"),
        "P^2 DLP closed form: delta_H = delta(mu) is a theorem (HN length one).",
    ),
    HNMode.PAPER: Certificate(
        Rigor.PROVEN,
        ("HN-length-one datum supplied by a certified paper table (sharp "
         "polarization-dependent delta_H target)",),
        ("arXiv:1907.06739", "arXiv:1910.14060"),
        "paper-tabulated sharp delta_H (Coskun-Huizenga / Levine-Zhang; E11-M4).",
    ),
    HNMode.ORACLE: Certificate(
        Rigor.PROVEN,
        ("HN-length-one datum supplied by an M2/OSCAR-constructed "
         "prioritary-sheaf HN filtration (G16 oracle)",),
        ("arXiv:1907.06739",),
        "M2/OSCAR-constructed prioritary-sheaf HN filtration (E11-M5 oracle hook).",
    ),
    HNMode.HEURISTIC: Certificate(
        Rigor.HEURISTIC,
        ("no certified HN-length-one datum: delta_H defaults to the Bogomolov "
         "floor 0; the sharp polarization-dependent delta_H is delegated",),
        ("arXiv:1907.06739", "arXiv:1910.14060"),
        "HEURISTIC Bogomolov floor; sharp delta_H delegated to the paper table "
        "(E11-M4) or the M2/OSCAR oracle (E11-M5).",
    ),
}


@dataclass(frozen=True)
class NonemptinessVerdict:
    """The result of the in-architecture ``Delta >= delta_H`` evaluation.

    Attributes
    ----------
    nonempty : bool
        The full non-emptiness verdict: ``Delta(xi) >= delta_H(xi)`` **OR**
        (on ``P^2``) ``xi`` is an exceptional bundle (DLP theorem, both
        disjuncts).  Matches :func:`bridgeland_stability.dlp.moduli_nonempty`.
    discriminant : Fraction
        ``Delta(xi)`` in the Coskun-Huizenga normalization.
    delta_H : Fraction
        The ``delta_H(xi)`` bound actually used (sharp on P^2 / from a certified
        supplier; the Bogomolov floor ``0`` in the HEURISTIC fallback).
    mode : HNMode
        Which supplier produced the HN-length-one datum.
    certificate : Certificate
        The G5 provenance stamp (``PROVEN`` iff ``mode`` is certified).  It
        certifies the (non-)emptiness *verdict*; a ``PROVEN`` cert is never
        readable as "provably empty" for an existing class -- see ``exceptional``.
    reason : str
        A human-readable one-line summary of the comparison and the mode.
    exceptional : bool
        ``True`` iff non-emptiness holds via the DLP *exceptional-bundle*
        disjunct on ``P^2`` -- i.e. ``xi`` is a genuine exceptional bundle, a
        non-empty isolated point that may sit STRICTLY BELOW the curve
        (``Delta < delta_H``).  Distinguishes "above the curve" from "isolated
        exceptional point below it"; defaults to ``False`` (backward-compatible).
    """

    nonempty: bool
    discriminant: Fraction
    delta_H: Fraction
    mode: HNMode
    certificate: Certificate
    reason: str
    exceptional: bool = False

    @property
    def certified(self) -> bool:
        """``True`` iff the verdict's rigor is ``PROVEN`` (certified HN datum)."""
        return self.certificate.rigor == Rigor.PROVEN


def _mu(xi: SurfaceBundle, surface: Surface) -> Fraction:
    """Exact slope ``mu = <c1, H>/(r d)`` in ``surface``'s NS lattice."""
    return surface.lattice.pairing(xi.c1, surface.H) / (xi.r * surface.d)


def _is_p2_exceptional(xi: SurfaceBundle, surface: Surface) -> bool:
    """``True`` iff ``xi`` is a genuine exceptional bundle on ``P^2``.

    The Drezet-Le Potier non-emptiness theorem on ``P^2`` has TWO disjuncts:
    ``M(xi)`` is non-empty iff ``Delta >= delta(mu)`` OR ``xi`` is (the Chern
    character of) an exceptional bundle -- a non-empty isolated point that may
    sit STRICTLY BELOW the curve.  This detects the second disjunct, delegating
    to the pinned :func:`bridgeland_stability.exceptional.is_exceptional`
    (``chi(E,E)=1`` and integral ``c2``), so ``moduli_nonempty`` agrees
    bit-for-bit with :func:`bridgeland_stability.dlp.moduli_nonempty` and never
    reports a PROVEN ``nonempty=False`` for a class that exists.

    Guarded to ``P^2`` only (``Pic(P^2) = Z.H``): a non-``P^2`` surface, a
    higher-Picard class (``len(c1) != 1``), or a non-integral ``c1`` is never an
    exceptional bundle in this sense and returns ``False``.
    """
    if not surface.is_p2 or len(xi.c1) != 1:
        return False
    c1 = xi.c1[0]
    if c1.denominator != 1:
        return False
    return is_exceptional(Bundle(xi.r, int(c1), xi.ch2))


def discriminant_H(xi: SurfaceBundle, surface: Surface) -> Fraction:
    """Exact discriminant ``Delta(xi)`` on ``surface`` (Coskun-Huizenga norm.).

    ``Delta = 1/2 mu^2 - ch2/(r d)`` with ``mu = <c1, H>/(r d)`` and ``d = H^2``.
    Bit-for-bit the CH convention used everywhere in this package (see
    :meth:`bridgeland_stability.chern.ChernChar.discriminant`); never the doubled
    ``discriminant_brief``.
    """
    d = surface.d
    mu = _mu(xi, surface)
    return Fraction(1, 2) * mu * mu - xi.ch2 / (xi.r * d)


def delta_H(xi: SurfaceBundle, surface: Surface, R_max: int = 60) -> Fraction:
    """The sharp non-emptiness bound ``delta_H(xi)`` for ``surface``'s polarization.

    * On ``P^2`` (``surface.is_p2``) this is the PROVEN Drezet-Le Potier
      closed-form curve :func:`bridgeland_stability.dlp.delta` evaluated at
      ``mu`` -- regressing exactly to the pinned ``delta(1/2)=5/8``,
      ``delta(1/3)=5/9``, ... .
    * Off ``P^2`` there is no closed ``delta``-curve; the sharp,
      polarization-dependent bound is the paper-tabulated / oracle-supplied datum
      of E11-M4/M5.  This THIN slice returns the **Bogomolov floor** ``0`` (the
      sanity floor ``Delta >= 0``); callers wanting a certified sharp target pass
      it to :func:`moduli_nonempty` as ``delta_H_target``.

    The G12 faithful-computation guard is applied first: a torsion-canonical
    surface (Enriques / bielliptic) is refused with the NS-lattice-refactor
    error rather than silently mis-modelled.  Always exact (``Fraction``).
    """
    require_faithful_computation(surface)  # G12 guard: torsion-canonical rows refused
    if surface.is_p2:
        mu = _mu(xi, surface)
        bundles = enumerate_exceptional(mu - 3, mu + 3, R_max)
        return delta(mu, bundles)
    return Fraction(0)  # THIN slice: Bogomolov floor; sharp delta_H is E11-M4/M5.


def moduli_nonempty(
    r: int,
    c1: Sequence[Number],
    ch2: Number,
    surface: Surface,
    *,
    delta_H_target: Optional[Number] = None,
    hn_source: Optional[HNMode] = None,
    R_max: int = 60,
) -> NonemptinessVerdict:
    """Decide ``Delta(xi) >= delta_H(xi)`` for ``xi = (r, c1, ch2)`` on ``surface``.

    GIVEN the HN-length-one datum (implicit on ``P^2`` via the DLP closed form,
    or supplied by a certified ``PAPER`` / ``ORACLE`` source), evaluate the sharp
    inequality exactly in ``Fraction`` and return a :class:`NonemptinessVerdict`
    carrying a G5 :class:`~bridgeland_stability.rigor.Certificate`.

    Mode selection
    --------------
    * ``delta_H_target`` given -- use that certified external target; it MUST
      come with a certified ``hn_source`` (``DLP`` / ``PAPER`` / ``ORACLE``),
      else ``ValueError``.  (The E11-M4 paper-table / E11-M5 oracle entry point.)
    * else on ``P^2`` -- ``DLP`` mode, the closed-form ``delta_H`` (PROVEN).
    * else (off ``P^2``, no target) -- ``HEURISTIC`` Bogomolov floor unless a
      certified ``hn_source`` is passed (E11-M5 hook).

    The verdict's certificate is ``PROVEN`` only for a certified mode; otherwise
    ``HEURISTIC``.  The G12 faithful-computation guard runs first.

    P^2 exceptional disjunct.  The full DLP verdict is
    ``Delta >= delta_H`` **OR** (on ``P^2``) ``xi`` is an exceptional bundle
    (a non-empty isolated point that may sit strictly below the curve).  Both
    disjuncts are honoured so the result matches
    :func:`bridgeland_stability.dlp.moduli_nonempty` and the DLP theorem; the
    returned ``exceptional`` flag records the second disjunct.  A ``PROVEN`` cert
    therefore certifies NON-emptiness -- never emptiness -- for an existing class.
    """
    require_faithful_computation(surface)  # G12 guard
    xi = SurfaceBundle(r, c1, ch2)
    disc = discriminant_H(xi, surface)

    if delta_H_target is not None:
        if hn_source not in _CERTIFIED:
            raise ValueError(
                "delta_H_target must come with a certified hn_source "
                "(HNMode.DLP / HNMode.PAPER / HNMode.ORACLE): an externally "
                "supplied sharp delta_H is only PROVEN when its HN-length-one "
                "hypothesis is certified (E11-M4 paper table / E11-M5 oracle)."
            )
        dH, mode = Fraction(delta_H_target), hn_source
    elif surface.is_p2:
        dH, mode = delta_H(xi, surface, R_max), HNMode.DLP
    else:
        dH = delta_H(xi, surface, R_max)
        mode = hn_source if hn_source in _CERTIFIED else HNMode.HEURISTIC

    cert = _MODE_CERT[mode]
    above_curve = disc >= dH
    exceptional = _is_p2_exceptional(xi, surface)  # DLP second disjunct (P^2 only)
    nonempty = above_curve or exceptional
    reason = f"mode={mode.value}: Delta={disc} {'>=' if above_curve else '<'} delta_H={dH}"
    if exceptional and not above_curve:
        reason += "; exceptional bundle: non-empty isolated point below the DLP curve"
    if mode not in _CERTIFIED:
        reason += " (HEURISTIC: Bogomolov floor; sharp delta_H is E11-M4/M5)"
    return NonemptinessVerdict(nonempty, disc, dH, mode, cert, reason, exceptional)
