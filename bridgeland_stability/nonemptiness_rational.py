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
from .varieties import Surface, P2, P1xP1, require_faithful_computation

__all__ = [
    "HNMode",
    "NonemptinessVerdict",
    "PaperDeltaHTarget",
    "discriminant_H",
    "delta_H",
    "moduli_nonempty",
    "paper_delta_H_targets",
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


# --------------------------------------------------------------------------
# E11-M4 / G18b [RESEARCH]: the fixed finite paper-tabulated delta_H table.
#
# Each entry pins an EXACT ``Fraction`` delta_H target + the primary source's
# yes/no verdict for one class ``(r, c1, ch2)`` on one polarized surface, fed
# into ``moduli_nonempty(..., delta_H_target=e.delta_H, hn_source=HNMode.PAPER)``
# so the in-architecture ``Delta >= delta_H`` inequality reproduces the paper's
# verdict with a certified (PROVEN) HN-length-one hypothesis.
#
# NORMALIZATION (the load-bearing research finding, R3).  The primary sources use
# the FULL-NS discriminant ``Delta_paper = 1/2 nu^2 - ch2/r`` (nu = c1/r an NS
# class), verified verbatim in arXiv:1907.06739 (Delta = 1/2 nu^2 - ch2/r) and
# its RR pairing ``chi(V,W) = r_V r_W (P(nu_W - nu_V) - Delta_V - Delta_W)`` with
# ``P(nu) = chi(O_X) + 1/2(nu^2 - nu.K_X)``.  The package's ``discriminant_H`` is
# the H-PROJECTED scalar ``1/2 mu^2 - ch2/(r d)``, ``mu = <c1,H>/(r d)``,
# ``d = H^2``.  On P^2 (d=1) the two coincide; for a class with ``c1`` PROPORTIONAL
# to H one has the exact identity ``Delta_paper = d * discriminant_H`` (machine-
# checked).  Every entry stores the CH-PACKAGE-NORMALIZED target ``delta_H`` (what
# the hook compares against ``discriminant_H``): ``delta_H == delta_H_paper`` on
# P^2 (d=1) and ``delta_H == delta_H_paper / d`` for the c1 || H entries on F_0.
# Only classes where ``discriminant_H`` is the right object to compare (P^2, or
# c1 || H on anticanonical F_0) are pinnable; a general off-P^2 non-diagonal class
# is NOT (the scalar model cannot reproduce the paper's full-NS verdict) and is an
# open question, per research-integrity rules R2/R4 (a small verified table beats a
# larger one with an invented value).
#
# The DLP peak values delta(mu) = chi(O_X) - Delta_E at an exceptional slope
# (Delta_E = (1 - 1/r_E^2)/2, classical Drezet-Le Potier) are the package's own
# already-pinned dlp.delta machinery (delta(1/2)=5/8, delta(1/3)=5/9,
# delta(2/5)=13/25) -- every P^2 target below is independently regressed against
# delta_H(xi, P2), and every F_0 target uses the exceptional line bundle at its
# integral slope (Delta_E=0 -> delta_H_paper = chi(O_{F_0}) = 1).
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class PaperDeltaHTarget:
    """One primary-source-tabulated ``delta_H`` target + verdict for a class.

    Attributes
    ----------
    surface : Surface
        The polarized surface the class lives on (``P2`` or ``P1xP1`` here).
    r : int
        Rank of the class.
    c1 : tuple
        First Chern class as an NS-vector in ``surface.lattice``'s basis
        (``(int, ...)``; coerced to ``Fraction`` by ``SurfaceBundle`` at use).
    ch2 : Fraction
        Second Chern character (exact ``Fraction``).
    delta_H : Fraction
        The CH-PACKAGE-NORMALIZED sharp target compared against ``discriminant_H``
        (``== delta_H_paper`` on P^2; ``== delta_H_paper / d`` for a ``c1 || H``
        class on ``F_0``).
    paper_nonempty : bool
        The primary source's yes/no non-emptiness verdict for this class.
    delta_H_paper : Fraction
        The paper's full-NS ``delta_H`` BEFORE the ``/d`` conversion
        (``== delta_H`` when ``d == 1``).
    citation : str
        Primary source: arXiv id + exact statement (Thm/Cor number) + the class
        ``(r, c1, ch2)`` + the polarization ``H``.
    note : str
        The exact ``Fraction`` arithmetic (R3): the paper datum -> the CH-normalized
        ``delta_H``, the independent ``Delta`` recompute, and the two-way verdict.
    """

    surface: Surface
    r: int
    c1: tuple
    ch2: Fraction
    delta_H: Fraction
    paper_nonempty: bool
    delta_H_paper: Fraction
    citation: str
    note: str = ""


_PAPER_DELTA_H_TABLE = (
    # -- P^2 = Levine-Zhang X_0 (m=0), the degree-9 anticanonical del Pezzo -----
    PaperDeltaHTarget(
        surface=P2, r=2, c1=(1,), ch2=Fraction(-5, 2),
        delta_H=Fraction(5, 8), paper_nonempty=True, delta_H_paper=Fraction(5, 8),
        citation=(
            "Levine-Zhang arXiv:1910.14060 Thm 1.4 (=Thm 7.15), m=0 (X_0=P^2); "
            "class (r,c1,ch2)=(2,(1),-5/2), H=-K=3H_0 (d=H_0^2=1). Existence via "
            "LZ Thm 1.4 (Delta=11/8>1/2, DL condition holds). delta_H=delta(1/2)="
            "5/8 classical Drezet-Le Potier (Ann. Sci. ENS 18, 1985)."
        ),
        note=(
            "mu=<1*H0,H0>/(2*1)=1/2; discriminant_H=1/2*(1/2)^2-(-5/2)/2=1/8+5/4=11/8; "
            "delta_H=delta(1/2)=chi(O_{P2})-Delta_{rk2 exc}=1-(1-1/2^2)/2=1-3/8=5/8; "
            "c2=1/2-(-5/2)=3. Two-way: 11/8>=5/8 -> nonempty; matches dlp.moduli_nonempty."
        ),
    ),
    PaperDeltaHTarget(
        surface=P2, r=3, c1=(0,), ch2=Fraction(-2),
        delta_H=Fraction(1), paper_nonempty=False, delta_H_paper=Fraction(1),
        citation=(
            "Levine-Zhang arXiv:1910.14060 Thm 1.4 (=Thm 7.15), m=0 (X_0=P^2); "
            "class (3,(0),-2), H=-K=3H_0 (d=1). LZ Thm 1.4 applies (Delta=2/3>1/2): "
            "the DL existence condition FAILS (2/3 < delta(0)=1) and the class is "
            "non-exceptional (rank 3 is not a Markov number), so NO -K-(semi)stable "
            "sheaf exists. delta(0)=1=chi(O_{P2}) classical Drezet-Le Potier."
        ),
        note=(
            "mu=0; discriminant_H=1/2*0-(-2)/3=2/3; delta_H=delta(0)=chi(O_{P2})-Delta_O="
            "1-0=1; c2=0-(-2)=2; is_exceptional(3,0,-2)=False (rank 3 not Markov). "
            "Two-way: 2/3<1 and non-exceptional -> EMPTY; matches dlp.moduli_nonempty."
        ),
    ),
    PaperDeltaHTarget(
        surface=P2, r=3, c1=(1,), ch2=Fraction(-5, 2),
        delta_H=Fraction(5, 9), paper_nonempty=True, delta_H_paper=Fraction(5, 9),
        citation=(
            "Levine-Zhang arXiv:1910.14060 Thm 1.4 (=Thm 7.15), m=0 (X_0=P^2); "
            "class (3,(1),-5/2), H=-K=3H_0 (d=1). LZ Thm 1.4 applies (Delta=8/9>1/2): "
            "DL condition holds (8/9 >= delta(1/3)=5/9) -> exists. delta(1/3)=5/9 "
            "classical Drezet-Le Potier (Ann. Sci. ENS 18, 1985)."
        ),
        note=(
            "mu=1/3; discriminant_H=1/2*(1/3)^2-(-5/2)/3=1/18+15/18=8/9; "
            "delta_H=delta(1/3)=chi(O)-Delta_{rk3 exc}=1-(1-1/3^2)/2=1-4/9=5/9; "
            "c2=1/2-(-5/2)=3. Two-way: 8/9>=5/9 -> nonempty; matches dlp.moduli_nonempty."
        ),
    ),
    PaperDeltaHTarget(
        surface=P2, r=5, c1=(2,), ch2=Fraction(-2),
        delta_H=Fraction(13, 25), paper_nonempty=True, delta_H_paper=Fraction(13, 25),
        citation=(
            "Coskun-Huizenga arXiv:1907.06739 (abstract + Cor 9.13: exceptional "
            "bundles are -K-stable on an anticanonically polarized del Pezzo); "
            "class (5,(2),-2) is the rank-5 slope-2/5 exceptional bundle on P^2 "
            "(Delta=Delta_E=12/25), H=-K=3H_0 (d=1). Nonempty via the exceptional "
            "bundle itself -- NOT via Levine-Zhang arXiv:1910.14060 Thm 1.4, whose "
            "Delta>1/2 hypothesis FAILS here (12/25<1/2). delta_H=delta(2/5)=13/25 "
            "classical Drezet-Le Potier."
        ),
        note=(
            "mu=2/5; discriminant_H=1/2*(2/5)^2-(-2)/5=2/25+10/25=12/25; "
            "Delta_E(r=5)=(1-1/5^2)/2=12/25 -> the class IS the exceptional bundle "
            "(is_exceptional=True, chi(E,E)=1); delta_H=delta(2/5)=13/25; c2=2-(-2)=4. "
            "12/25 < 13/25 (strictly below the curve) yet nonempty=True via the DLP "
            "exceptional-bundle disjunct (P^2 only); matches dlp.moduli_nonempty."
        ),
    ),
    # -- P^1 x P^1 = F_0, the degree-8 anticanonical del Pezzo (c1 || H) --------
    PaperDeltaHTarget(
        surface=P1xP1, r=2, c1=(0, 0), ch2=Fraction(-4),
        delta_H=Fraction(1, 2), paper_nonempty=True, delta_H_paper=Fraction(1),
        citation=(
            "Coskun-Huizenga arXiv:1907.06739 Thm 1.8 / Cor 9.13: on F_0 "
            "anticanonical delta_H(nu)=DLP_{-K_{F_0}}(nu). Class (2,(0,0),-4), "
            "H=(1,1) (same ample ray as -K=O(2,2)), d=H^2=2, NS Gram [[0,1],[1,0]]. "
            "At the line-bundle slope nu=(0,0) the controlling exceptional bundle is "
            "O (Delta_E=0), so delta_H_paper=P(0)-Delta_O=chi(O_{F_0})=1."
        ),
        note=(
            "nu=c1/r=(0,0), nu^2=0; Delta_paper=1/2*0-(-4)/2=2; discriminant_H: "
            "mu=<(0,0),(1,1)>/(2*2)=0, =1/2*0-(-4)/(2*2)=1; c1||H identity "
            "Delta_paper=d*discriminant_H=2*1=2; delta_H_paper=1, CH target=1/2; "
            "c2=0-(-4)=4. Two-way: Delta_paper=2>=delta_H_paper=1 <=> "
            "discriminant_H=1>=1/2 -> nonempty; exceptional=False (P^2-only disjunct)."
        ),
    ),
    PaperDeltaHTarget(
        surface=P1xP1, r=2, c1=(2, 2), ch2=Fraction(-2),
        delta_H=Fraction(1, 2), paper_nonempty=True, delta_H_paper=Fraction(1),
        citation=(
            "Coskun-Huizenga arXiv:1907.06739 Thm 1.8 / Cor 9.13: on F_0 "
            "anticanonical delta_H(nu)=DLP_{-K_{F_0}}(nu). Class (2,(2,2),-2), "
            "H=(1,1), d=2, NS Gram [[0,1],[1,0]]. At the diagonal line-bundle slope "
            "nu=(1,1) the controlling exceptional bundle is O(1,1) (Delta_E=0), so "
            "delta_H_paper=P(0)-Delta_{O(1,1)}=chi(O_{F_0})=1."
        ),
        note=(
            "nu=c1/r=(1,1), nu^2=<(1,1),(1,1)>=2; Delta_paper=1/2*2-(-2)/2=1+1=2; "
            "discriminant_H: mu=<(2,2),(1,1)>/(2*2)=1, =1/2*1^2-(-2)/(2*2)=1/2+1/2=1; "
            "c1||H identity Delta_paper=d*discriminant_H=2*1=2; delta_H_paper=1, "
            "CH target=1/2; c2=<(2,2),(2,2)>/2-(-2)=4+2=6. Two-way: 2>=1 <=> 1>=1/2 "
            "-> nonempty; exceptional=False (P^2-only disjunct)."
        ),
    ),
)


def paper_delta_H_targets() -> "tuple[PaperDeltaHTarget, ...]":
    """The fixed finite table of primary-source ``delta_H`` targets (E11-M4).

    Each entry carries an exact ``Fraction`` CH-package-normalized ``delta_H``, the
    paper's yes/no verdict, and a per-entry arXiv citation (arXiv:1907.06739
    Coskun-Huizenga / arXiv:1910.14060 Levine-Zhang) naming the exact statement and
    the class.  Feed an entry into
    ``moduli_nonempty(e.r, e.c1, e.ch2, e.surface, delta_H_target=e.delta_H,
    hn_source=HNMode.PAPER)`` to reproduce the paper verdict with a PROVEN
    HN-length-one hypothesis.
    """
    return _PAPER_DELTA_H_TABLE


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
