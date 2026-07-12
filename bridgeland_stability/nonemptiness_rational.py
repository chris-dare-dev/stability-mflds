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
* Coskun-Huizenga-Woolf, arXiv:1401.1613 Sec.2 Thm 2.2 -- the P^2 non-emptiness criterion this
  module implements (integrality + Delta >= delta(mu) OR exceptional).
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

from dataclasses import dataclass, field
from enum import Enum
from fractions import Fraction
from typing import Optional, Sequence, Union

from .dlp import delta
from .dlp_hirzebruch import (
    DEFAULT_RANK_MAX,
    discriminant,
    dlp_envelope,
    emptiness_bound,
    exceptional_ch2,
    hirzebruch_index,
    is_ample,
    is_semiexceptional,
    is_stable_exceptional,
    total_slope,
)
from .exceptional import (
    Bundle, P, certified_rank_cutoff, enumerate_exceptional, is_exceptional,
    is_semiexceptional_p2,
)
from .exceptional_surface import SurfaceBundle
from .rigor import Certificate, Rigor
from .varieties import Surface, P2, P1xP1, hirzebruch, require_faithful_computation

__all__ = [
    "HNMode",
    "VerdictStatus",
    "NonemptinessVerdict",
    "SharpBoundEvidence",
    "PaperDeltaHTarget",
    "discriminant",
    "discriminant_H",
    "delta_H",
    "hirzebruch_with_polarization",
    "moduli_nonempty",
    "validate_character",
    "is_semiexceptional_p2",
    "paper_delta_H_targets",
]

Number = Union[int, Fraction]


class VerdictStatus(str, Enum):
    """Branch-derived tri-state, computed from (nonempty, certificate.rigor) -- never the mode."""

    PROVEN_NONEMPTY = "proven_nonempty"
    PROVEN_EMPTY = "proven_empty"
    UNKNOWN = "unknown"


class HNMode(str, Enum):
    """Which supplier certified the HN-length-one hypothesis behind a verdict."""

    DLP = "dlp"            # PROVEN: P^2 Drezet-Le Potier closed form (HN length one implicit)
    PAPER = "paper"        # PROVEN: paper-tabulated HN-length-one datum (E11-M4 hook)
    ORACLE = "oracle"      # PROVEN: M2-constructed rank-1 ideal-sheaf witness (E11-M5 hook)
    HEURISTIC = "heuristic"  # HEURISTIC: no certified HN datum (Bogomolov floor)
    # -- E11-M6 / G18: the natively-computed Coskun-Huizenga envelope on F_e --------
    DLP_ANTICANONICAL = "dlp_anticanonical"  # e in {0,1}, H = -K: delta_H = DLP_{-K} (sharp)
    DLP_LOWER = "dlp_lower"                  # any other ample H on F_e: certified LOWER bound


#: The modes whose HN-length-one hypothesis is certified -> verdict rigor PROVEN.
#: (The two ``DLP_*`` envelope modes carry a per-verdict rigor instead -- see
#: :func:`_hirzebruch_verdict` -- because whether the verdict is PROVEN depends on which
#: side of the envelope ``Delta`` falls, not on the mode alone.)
_CERTIFIED = frozenset({HNMode.DLP, HNMode.PAPER, HNMode.ORACLE})

_CH = "arXiv:1907.06739"          # Coskun-Huizenga, Existence of semistable sheaves on F_e

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
        ("existence witnessed by a Macaulay2-constructed rank-1 ideal sheaf I_Z(c1) on P^2 "
         "(length l = c1^2/2 - ch2 = c2 >= 0), torsion-free of rank 1 hence mu-stable: a "
         "SUFFICIENT-ONLY construction returning True or None, never False; it does NOT compute "
         "the Harder-Narasimhan filtration of a prioritary sheaf and does NOT handle F_n",),
        ("arXiv:1907.06739",),
        "Macaulay2-constructed rank-1 ideal sheaf I_Z(c1) on P^2 (sufficient-only witness, "
        "returns True|None; see oracle/m2.py::moduli_nonempty_by_construction).",
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

#: Boundary certificate for an external (PAPER/ORACLE) target off P^2 with Delta == delta_H.
#: CH Thm "deltaSurface" (1) requires a STRICT inequality off P^2, so the boundary is UNKNOWN,
#: not PROVEN.  Rigor is non-PROVEN, so ``NonemptinessVerdict.status`` maps it to UNKNOWN;
#: HEURISTIC matches ``_hirzebruch_verdict``'s own boundary handling.
_BOUNDARY_CERT = Certificate(
    Rigor.HEURISTIC,
    ("Delta == delta_H off P^2: CH Thm 'deltaSurface' (1) requires a STRICT inequality; "
     "the boundary is undecided",),
    (_CH,),
    "boundary Delta == delta_H off P^2 is UNKNOWN (a strict inequality is required).",
)

#: Certificate for a character that is not the Chern character of any coherent sheaf, so
#: ``M(xi)`` is empty on EVERY surface for EVERY polarization.  A coherent sheaf has integral
#: Chern classes -- ``c1 in NS`` and ``c2 = 1/2 <c1,c1> - ch2 in Z`` (``c2 in H^4(X,Z) = Z``) --
#: and by Riemann-Roch (``c1.(c1 - K_X)`` is even on any surface, Wu) these two force ``chi in Z``,
#: so the ``(c1, c2)``-integrality clause is the WHOLE Thm 2.2 integrality condition and needs no
#: canonical class -- exact off ``P^2`` too.  ``PROVEN`` empty is honest and unfalsifiable here:
#: the emptiness is a theorem (no counterexample can exist), matching the ``P^2`` path's verdict
#: for the same defect (A3).
_INVALID_CHARACTER_CERT = Certificate(
    Rigor.PROVEN,
    ("a coherent sheaf has integral Chern classes: c1 in NS and "
     "c2 = 1/2 <c1,c1> - ch2 in Z; Riemann-Roch (c1.(c1 - K) even) then forces chi in Z",),
    (_CH,),
    "invalid Chern character (non-integral c1 or c2): no coherent sheaf has this character, "
    "so M(xi) is empty for every polarization.",
)


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

    @property
    def status(self) -> VerdictStatus:
        """PROVEN_NONEMPTY / PROVEN_EMPTY / UNKNOWN, derived from the branch (rigor), not the mode."""
        if self.certificate.rigor is not Rigor.PROVEN:
            return VerdictStatus.UNKNOWN
        return VerdictStatus.PROVEN_NONEMPTY if self.nonempty else VerdictStatus.PROVEN_EMPTY


# --------------------------------------------------------------------------
# E12-M4 / A5: class-bound sharp-bound evidence (retires the forgeable pair).
#
# The (delta_H_target, hn_source) pair let a caller assert ANY sharp bound for ANY
# class and still mint Rigor.PROVEN (defect A5: delta_H_target=10**6, ORACLE ->
# PROVEN).  A ``SharpBoundEvidence`` fixes both leaks:
#   * it is BOUND to a specific class ``(surface, r, c1, ch2)`` via ``matches`` --
#     evidence derived for one class can never certify another (crit. 1); and
#   * it separates the two independent claims the audit called out (crit. 2): the
#     VALUE claim "the sharp bound equals this" (``sharp_bound`` + its provenance
#     ``sharp_bound_source``) from the SHEAF-THEORETIC claim "the generic
#     prioritary HN filtration has length one" (``hn_length_one_source``).
# The verdict is still PROVEN only when the supplied ``sharp_bound`` EQUALS the
# package's OWN independently-certified sharp bound for that exact class
# (:func:`_certified_sharp_bound`): a mismatch is a forged / mis-derived target
# and raises ``ValueError``, never a PROVEN verdict (invariant 7).
#
# ORACLE-sourced evidence is a CAPABILITY object: it carries a module-private
# token that only :func:`bridgeland_stability.oracle.mint_oracle_evidence` holds,
# and that mint runs only after a real construction returned ``True`` (crit. 3).
# A raw ``(delta_H_target, hn_source=ORACLE)`` pair can no longer forge one.
# --------------------------------------------------------------------------
#: Module-private capability token.  Only ``oracle.mint_oracle_evidence`` imports and
#: passes it; a direct ORACLE-sourced construction lacks it and is refused (crit. 3).
_ORACLE_TOKEN = object()


@dataclass(frozen=True)
class SharpBoundEvidence:
    """A class-bound certificate that a sharp ``delta_H`` and an HN-length-one datum hold.

    Retires the forgeable ``(delta_H_target, hn_source)`` pair (A5).  The evidence is
    tied to one Chern character on one polarized surface (:meth:`matches`), and it
    separates the two distinct claims the audit conflated into a single enum:

    Attributes
    ----------
    surface, r, c1, ch2 :
        The exact class this evidence was derived for.  :meth:`matches` refuses to let
        it certify any other class (crit. 1).
    sharp_bound : Fraction
        CLAIM 1 (the VALUE claim): "the sharp non-emptiness bound ``delta_H`` for this
        class equals this value".  :func:`moduli_nonempty` honours it only if it equals
        the package's own :func:`_certified_sharp_bound` for the class.
    sharp_bound_source : HNMode
        Provenance of claim 1 (``DLP`` / ``PAPER`` / ``ORACLE``).
    hn_length_one_source : HNMode
        CLAIM 2 (the SHEAF-THEORETIC claim): the supplier of the "generic prioritary HN
        filtration has length one" hypothesis.  This is the ``mode`` the verdict reports;
        an ``ORACLE`` value requires the capability token (see ``__post_init__``).
    citation : str
        Free-text provenance note.
    """

    surface: Surface
    r: int
    c1: tuple
    ch2: Fraction
    sharp_bound: Fraction                 # CLAIM 1: "the sharp bound equals this value"
    sharp_bound_source: HNMode            # provenance of claim 1 (DLP / PAPER / ORACLE)
    hn_length_one_source: HNMode          # CLAIM 2: "the generic prioritary HN filtration has length one"
    citation: str = ""
    #: Kept out of equality/repr so value semantics are preserved; only the oracle mint sets it.
    _oracle_token: object = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.hn_length_one_source is HNMode.ORACLE and self._oracle_token is not _ORACLE_TOKEN:
            raise TypeError(
                "ORACLE-sourced evidence must be minted via "
                "bridgeland_stability.oracle.mint_oracle_evidence after a construction "
                "actually returned True (E12-M4 crit. 3); it cannot be built directly."
            )

    def matches(self, r, c1, ch2, surface) -> bool:
        """``True`` iff this evidence was derived for exactly the queried class."""
        return (self.r == r
                and tuple(Fraction(x) for x in self.c1) == tuple(Fraction(x) for x in c1)
                and Fraction(self.ch2) == Fraction(ch2)
                and self.surface == surface)


# --------------------------------------------------------------------------
# E11-M4 / G18b [RESEARCH]; relabelled E12-M5 (A12): REGRESSION FIXTURE of hand-derived
# delta_H targets -- NOT a verbatim paper table.  Each delta_H is DERIVED from general
# theorems (the Drezet-Le Potier closed form delta(mu) on P^2; DLP_{-K} on the del Pezzo
# F_e) and independently REGRESSED against the package's own delta_H / dlp_envelope (see
# tests/test_nonemptiness.py::test_paper_p2_targets_match_native_dlp and
# ::test_paper_p1xp1_targets_match_native_envelope).  The per-entry arXiv citation names
# the primary source for the EXISTENCE verdict, not for the numeric value.
#
# Each entry pins an EXACT ``Fraction`` delta_H target + the primary source's
# yes/no verdict for one class ``(r, c1, ch2)`` on one polarized surface, fed
# into ``moduli_nonempty(..., delta_H_target=e.delta_H, hn_source=HNMode.PAPER)``
# so the in-architecture ``Delta >= delta_H`` inequality reproduces the paper's
# verdict with a certified (PROVEN) HN-length-one hypothesis.
#
# NORMALIZATION (E11-M4's load-bearing research finding R3, and its E11-M6 RESOLUTION).
# The primary sources use the FULL-NS discriminant ``Delta = 1/2 nu^2 - ch2/r``
# (nu = c1/r an NS class), verified verbatim in arXiv:1907.06739 Sec. 2.1, together with
# ``chi(V,W) = r_V r_W (P(nu_W - nu_V) - Delta_V - Delta_W)`` and
# ``P(nu) = chi(O_X) + 1/2(nu^2 - nu.K_X)``.  M4 could not compare against that Delta --
# the package then only had the H-PROJECTED scalar ``discriminant_H = 1/2 mu^2 -
# ch2/(r d)`` -- so each entry stored a rescaled target ``delta_H = delta_H_paper / d``,
# and only classes with ``c1`` PROPORTIONAL to ``H`` (where ``Delta = d * discriminant_H``)
# were pinnable.
#
# E11-M6 / G18 ships the genuine full-NS ``discriminant`` (dlp_hirzebruch.discriminant),
# and ``moduli_nonempty`` now compares against IT.  So the rescaling is gone: every entry
# stores the paper's own target, ``delta_H == delta_H_paper``, on every surface.  The
# ``delta_H_paper`` field is kept (equal to ``delta_H``) so the historical distinction
# stays legible and the /d identity remains testable.  See docs/CORRECTIONS.md (C7).
#
# The two F_0 targets below are no longer merely transcribed: ``dlp_envelope`` computes
# ``DLP_{-K_{F_0}}(nu) = 1`` natively at both slopes and certifies the rank truncation as
# exact (tests/test_dlp_hirzebruch.py::test_delta_K_at_line_bundle_slopes_is_one), so the
# table is now a REGRESSION of the in-package envelope against the primary source rather
# than the package's only source of truth.
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
    """One hand-derived, regression-checked ``delta_H`` target + verdict for a class.

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
        The sharp target, in the paper's own (full-NS) normalization, compared against
        :func:`~bridgeland_stability.dlp_hirzebruch.discriminant`.  Since E11-M6 this
        equals ``delta_H_paper`` on every surface (the old ``/d`` rescaling is gone).
    paper_nonempty : bool
        The primary source's yes/no non-emptiness verdict for this class.
    delta_H_paper : Fraction
        The paper's full-NS ``delta_H``.  Retained (now always ``== delta_H``) so the
        historical M4 ``/d`` normalization finding stays legible and testable.
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
            "delta(1/3)=P(-1/3)-Delta_O=5/9-0=5/9, attained by O (rank 1, slope 0): "
            "P(-1/3)=((1/3)^2-3*(1/3)+2)/2=(1/9+1)/2=5/9 and Delta_O=(1-1/1^2)/2=0. "
            "There is NO rank-3 exceptional bundle (rank 3 is not Markov; Bundle.from_slope(1/3) "
            "has c2=5/3 not integral, see test_rank3_pseudobundle_does_not_exist), so the earlier "
            "chi(O)-minus-a-rank-3-exceptional-discriminant derivation was fictitious -- it agreed "
            "only by the numerical coincidence P(-1/3)=1-4/9. c2=1/2-(-5/2)=3. Two-way: 8/9>=5/9 "
            "-> nonempty; matches dlp.moduli_nonempty."
        ),
    ),
    PaperDeltaHTarget(
        surface=P2, r=5, c1=(2,), ch2=Fraction(-2),
        delta_H=Fraction(13, 25), paper_nonempty=True, delta_H_paper=Fraction(13, 25),
        citation=(
            "Drezet-Le Potier, Ann. Sci. ENS 18 (1985), Thm A: the rank-5 slope-2/5 "
            "exceptional bundle on P^2 EXISTS (2/5 is in the image of epsilon; rank 5 is "
            "Markov), Delta=Delta_E=12/25; its moduli space is a single reduced point, so "
            "M(5,(2),-2) is nonempty via the exceptional bundle itself -- NOT via Levine-Zhang "
            "arXiv:1910.14060 Thm 1.4, whose Delta>1/2 hypothesis FAILS here (12/25<1/2). "
            "delta(2/5)=P(0)-Delta_E=1-12/25=13/25 (classical DLP). Provenance (E12-M5/A12): CH "
            "arXiv:1907.06739 Cor 9.13 states delta^{mu-s}_{1-e/2}(nu)=DLP_{-K}(nu) on the del "
            "Pezzo F_e (e in {0,1}); the -K-stability of exceptional bundles is attributed to "
            "Gorodentsev, NOT a statement of Cor 9.13. H=-K=3H_0 (d=1)."
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
        delta_H=Fraction(1), paper_nonempty=True, delta_H_paper=Fraction(1),
        citation=(
            "Coskun-Huizenga arXiv:1907.06739 Cor 'deltaDLP': on F_0 with the "
            "anticanonical polarization delta_H^{mu-s}(nu)=DLP_{-K_{F_0}}(nu). Class "
            "(2,(0,0),-4), H=(1,1) (same ample ray as -K=O(2,2)), d=H^2=2, NS Gram "
            "[[0,1],[1,0]]. At the line-bundle slope nu=(0,0) the controlling exceptional "
            "bundle is O (Delta_E=0), so delta_H=P(0)-Delta_O=chi(O_{F_0})=1."
        ),
        note=(
            "nu=c1/r=(0,0), nu^2=0; full-NS Delta=1/2*0-(-4)/2=2; delta_H=1; c2=0-(-4)=4. "
            "Two-way: 2>=1 -> nonempty.  The /d relic: discriminant_H=1 and c1||H gives "
            "Delta=d*discriminant_H=2*1=2.  E11-M6 regression: dlp_envelope((0,0),P1xP1) "
            "computes delta_H=1 natively, exact=True, sharp=True, witness=O."
        ),
    ),
    PaperDeltaHTarget(
        surface=P1xP1, r=2, c1=(2, 2), ch2=Fraction(-2),
        delta_H=Fraction(1), paper_nonempty=True, delta_H_paper=Fraction(1),
        citation=(
            "Coskun-Huizenga arXiv:1907.06739 Cor 'deltaDLP': on F_0 with the "
            "anticanonical polarization delta_H^{mu-s}(nu)=DLP_{-K_{F_0}}(nu). Class "
            "(2,(2,2),-2), H=(1,1), d=2, NS Gram [[0,1],[1,0]]. At the diagonal "
            "line-bundle slope nu=(1,1) the controlling exceptional bundle is O(1,1) "
            "(Delta_E=0), so delta_H=P(0)-Delta_{O(1,1)}=chi(O_{F_0})=1."
        ),
        note=(
            "nu=c1/r=(1,1), nu^2=<(1,1),(1,1)>=2; full-NS Delta=1/2*2-(-2)/2=1+1=2; "
            "delta_H=1; c2=<(2,2),(2,2)>/2-(-2)=4+2=6. Two-way: 2>=1 -> nonempty.  The "
            "/d relic: discriminant_H=1 and c1||H gives Delta=d*discriminant_H=2. "
            "E11-M6 regression: dlp_envelope((1,1),P1xP1)=1, exact=True, witness=O(1,1)."
        ),
    ),
)


def paper_delta_H_targets() -> "tuple[PaperDeltaHTarget, ...]":
    """Regression fixture of hand-derived ``delta_H`` targets (E11-M4; relabelled E12-M5/A12).

    **Not a verbatim paper table**: each ``delta_H`` is derived from general theorems plus the
    package's own DLP machinery and regressed against it.  Each entry carries an exact
    ``Fraction`` CH-package-normalized ``delta_H``, the paper's yes/no verdict, and a per-entry
    arXiv citation (arXiv:1907.06739 Coskun-Huizenga / arXiv:1910.14060 Levine-Zhang) naming the
    primary source for the EXISTENCE verdict, not for the numeric value.  Feed an entry into
    ``moduli_nonempty(e.r, e.c1, e.ch2, e.surface, delta_H_target=e.delta_H,
    hn_source=HNMode.PAPER)`` to reproduce the paper verdict with a PROVEN HN-length-one
    hypothesis.
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
    to the pinned :func:`bridgeland_stability.exceptional.is_exceptional`, whose
    test is membership of the slope in the Drezet-Le Potier epsilon-recursion
    image (DLP Theoreme A; ``chi(E,E)=1`` and integral ``c2`` are only a NECESSARY
    pre-filter -- E12-M1 corrected the old ``chi==1``-only test that admitted
    Markov-rank impostors such as ``(610,133,-581/2)``).  So ``moduli_nonempty``
    agrees bit-for-bit with :func:`bridgeland_stability.dlp.moduli_nonempty` and
    never reports a PROVEN ``nonempty=False`` for a class that exists.

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


def _is_p2_semiexceptional(xi: SurfaceBundle, surface: Surface) -> bool:
    """True iff xi = m*ch(E) for an exceptional bundle E on P^2 (m>=1).  Guarded to P^2.

    The m>=1 disjunct of the Drezet-Le Potier / Coskun-Huizenga non-emptiness theorem on
    P^2 (arXiv:1907.06739 Ex.1.14: a semiexceptional bundle is a direct sum of copies of an
    exceptional bundle).  Delegates to :func:`bridgeland_stability.exceptional.
    is_semiexceptional_p2`, so it agrees with the oracle's ``reference_is_semiexceptional``.
    A non-``P^2`` surface, a higher-Picard class, or a non-integral ``c1`` returns ``False``.
    """
    if not surface.is_p2 or len(xi.c1) != 1:
        return False
    c1 = xi.c1[0]
    if c1.denominator != 1:
        return False
    return is_semiexceptional_p2(xi.r, int(c1), xi.ch2)


def _exceptional_disjunct(xi: SurfaceBundle, surface: Surface) -> "tuple[bool, bool]":
    """``(exceptional, semiexceptional)``: does ``xi`` satisfy the Drezet-Le Potier
    *exceptional-bundle* non-emptiness disjunct on ``surface``?

    A genuine (semi)exceptional class is a NON-empty moduli point (an exceptional bundle's
    moduli space is a single reduced point; ``V^{+m}`` is Gieseker-semistable) that sits
    STRICTLY BELOW the sharp envelope, so a ``Delta``-vs-``delta_H`` comparison alone
    under-reports it.  :func:`moduli_nonempty` must OR this disjunct in on **every** surface
    before it may conclude a PROVEN emptiness -- exactly as :func:`_hirzebruch_verdict` does
    natively -- so the certified-target path can never forge a PROVEN "empty" for a class
    fed its OWN correct sharp ``delta_H`` (E12-M2; the off-``P^2`` analogue of the pinned
    ``test_paper_exceptional_coexists_with_target``).

    * ``P^2`` -- the pinned P^2 detectors (:func:`_is_p2_exceptional` /
      :func:`_is_p2_semiexceptional`).
    * an **ample-polarized** ``F_e`` -- the surface-native
      :func:`bridgeland_stability.dlp_hirzebruch.is_semiexceptional`.  It requires
      ``Delta == Delta_V``, hence already subsumes the ``m = 1`` pure-exceptional case with
      the correct ``ch2 == exceptional_ch2`` -- so it does NOT carry the A6 raw-
      ``is_stable_exceptional`` bug (which is a function of ``(r, c1)`` only).  The ``m = 1``
      sub-flag is re-derived under the SAME ``ch2`` guard, purely for an honest reason string.
    * any surface with no wired exceptional-bundle theory (K3 / abelian / an ``F_n`` carrying
      only the nef-and-big factory polarization) -- ``(False, False)``, unchanged.
    """
    if surface.is_p2:
        return _is_p2_exceptional(xi, surface), _is_p2_semiexceptional(xi, surface)
    try:                                    # only ample-polarized F_e carry the CH theory
        hirzebruch_index(surface)
    except NotImplementedError:
        return False, False
    if not is_ample(surface) or xi.r < 1 or not _is_integral_c1(xi):
        return False, False
    semiexceptional = is_semiexceptional(xi, surface)      # ch2-guarded (Delta == Delta_V)
    c1i = tuple(int(c) for c in xi.c1)
    # m = 1 iff xi is itself a mu_H-stable exceptional bundle AND its ch2 is the exceptional
    # value; guarding on `semiexceptional` short-circuits and keeps A6's (r,c1)-only test out.
    exceptional = bool(
        semiexceptional
        and is_stable_exceptional(xi.r, c1i, surface)
        and xi.ch2 == exceptional_ch2(xi.r, c1i, surface)
    )
    return exceptional, semiexceptional


def hirzebruch_with_polarization(n: int, H: Sequence[Number]) -> Surface:
    """Build F_n (n>=0) carrying a GIVEN strictly-ample polarization H = a f + b s.

    Constructs a fresh frozen ``Surface`` on the FIXED F_n NS lattice (Gram
    ``[[0,1],[1,-n]]``, basis f, s), with ``d = <H,H>`` computed exactly so
    ``discriminant_H`` / ``_mu`` stay self-consistent.  ``n = 0`` gives the genuine
    rank-2 ``F_0 = P^1 x P^1`` lattice (the :func:`~bridgeland_stability.varieties.
    hirzebruch` factory falls back to the rank-1 shim there, so we take the lattice
    from :data:`~bridgeland_stability.varieties.P1xP1`).

    Strict ampleness on F_n (Nakai): ``H.f = b > 0`` and ``H.C_0 = a - n b > 0``.  The
    factory polarization ``H = n f + s = (n,1)`` is only nef-and-big (``a - n b = 0``)
    and is honestly REFUSED here.  Note that for ``n >= 2`` the anticanonical class
    ``-K = (n+2) f + 2 s`` is itself NOT ample (``a - n b = 2 - n <= 0``) -- F_n is a
    del Pezzo surface only for ``n = 0, 1``.

    Every ample class on ``F_n`` lies on the ray of some ``H`` with integer coordinates,
    and both ``DLP_{H,V}`` and ``mu_H``-stability depend on ``H`` only through its ray,
    so integer coordinates lose no generality: the polarization ``H_m = E + (n+m) F``
    of arXiv:1907.06739 is the ray of ``(n + m) f + s``, i.e. ``(p, q)`` with
    ``m = p/q - n``.
    """
    if n < 0:
        raise ValueError("hirzebruch_with_polarization needs n>=0")
    lat = P1xP1.ns_lattice if n == 0 else hirzebruch(n).ns_lattice
    if lat is None or lat.rank != 2:            # defensive
        raise ValueError("expected the rank-2 F_n NS lattice")
    if any(int(x) != x for x in H):
        raise ValueError(
            f"hirzebruch_with_polarization needs an integral polarization vector "
            f"(pass the primitive integer point on the ample ray, e.g. (3,2) not "
            f"(3/2,1)); got {tuple(H)!r}")                     # A9: no silent int() truncation
    Hv = tuple(int(x) for x in H)
    if len(Hv) != 2:
        raise ValueError("H must be a length-2 NS vector (a f + b s)")
    a, b = Hv
    if not (b > 0 and a - n * b > 0):
        raise ValueError(
            f"H={Hv} is not strictly ample on F_{n} (Nakai needs b>0 and a>n*b; the "
            f"factory H=({n},1) is only nef-and-big, and -K is not ample for n>=2)")
    d = lat.self_pairing(Hv)                     # = 2ab - n b^2, an exact integer
    return Surface(
        # K_{F_n} = -(n+2) f - 2 s (A8); derived K.H = (n-2)b - 2a on this polarization.
        name=f"F_{n} (H={a}f+{b}s, d={d})", d=int(d), K=(-(n + 2), -2), chi_O=1,
        picard_rank=2, kind="hirzebruch",
        note=f"F_{n} strictly-ample H={Hv}; carries the CH Drezet-Le Potier envelope.",
        H=Hv, ns_lattice=lat,
    )


def discriminant_H(xi: SurfaceBundle, surface: Surface) -> Fraction:
    """The **H-projected** scalar discriminant ``1/2 mu^2 - ch2/(r d)``, ``d = H^2``.

    .. warning::
       This is **not** the Coskun-Huizenga discriminant when ``rho(X) >= 2``.  The
       primary sources use the *full-NS* ``Delta = 1/2 <nu,nu> - ch2/r`` with
       ``nu = c1/r`` -- see :func:`bridgeland_stability.dlp_hirzebruch.discriminant`,
       which is what :func:`moduli_nonempty` compares against.  The two are related by

           ``Delta = d * discriminant_H``   **iff**  ``c1`` is proportional to ``H``,

       which holds automatically on every Picard-rank-1 surface, and on ``P^2``
       (``d = 1``) they are equal.  For a non-diagonal ``c1`` on ``P^1 x P^1`` or ``F_n``
       they genuinely differ, and ``discriminant_H`` -- being built from ``mu_H`` --
       spuriously *depends on the polarization*, whereas ``Delta`` does not.

    Retained as the H-numerical scalar of the ``(r, ch1.H, ch2)`` model (it agrees
    bit-for-bit with :meth:`bridgeland_stability.chern.ChernChar.discriminant`) and for
    comparison, exactly as ``discriminant_brief`` is retained.  See ``docs/CORRECTIONS.md``
    (C7).
    """
    d = surface.d
    mu = _mu(xi, surface)
    return Fraction(1, 2) * mu * mu - xi.ch2 / (xi.r * d)


def _hirzebruch_envelope(xi: SurfaceBundle, surface: Surface, rank_max: int):
    """The CH envelope on an ample-polarized ``F_e``, or ``None`` if that theory
    does not apply to ``surface`` (not a Hirzebruch surface, or ``H`` not ample)."""
    try:
        hirzebruch_index(surface)
    except NotImplementedError:
        return None                       # P^2, K3, abelian, ... : no F_e NS lattice
    if not is_ample(surface):
        return None                       # the hirzebruch(n) factory H is only nef-and-big
    return dlp_envelope(total_slope(xi), surface, rank_max)


def _fe_emptiness_bound(
    xi: SurfaceBundle, surface: Surface, rank_max: int
) -> Optional[Fraction]:
    """``emptiness_bound(xi, surface)`` on an ample-polarized ``F_e``, else ``None``.

    Mirrors :func:`_hirzebruch_envelope`'s guard: the certified emptiness theorem
    (arXiv:1907.06739 Sec. 5.4 -- a ``mu_H``-stable exceptional bundle forcing ``chi <= 0``)
    is defined only on a Hirzebruch surface carrying a strictly ample ``H``.  ``None`` signals
    "no ``F_e`` emptiness theory here" to the certified-target tail of
    :func:`moduli_nonempty`, which then falls back to the boundary-only downgrade.
    """
    try:
        hirzebruch_index(surface)
    except NotImplementedError:
        return None                       # P^2, K3, abelian, ... : no F_e emptiness bound
    if not is_ample(surface):
        return None                       # the hirzebruch(n) factory H is only nef-and-big
    return emptiness_bound(xi, surface, rank_max)


def delta_H(
    xi: SurfaceBundle,
    surface: Surface,
    R_max: int = 60,
    rank_max: int = DEFAULT_RANK_MAX,
) -> Fraction:
    """The sharp non-emptiness bound ``delta_H(xi)`` for ``surface``'s polarization.

    * On ``P^2`` this is the PROVEN Drezet-Le Potier closed-form curve
      :func:`bridgeland_stability.dlp.delta` evaluated at ``mu`` -- regressing exactly to
      the pinned ``delta(1/2)=5/8``, ``delta(1/3)=5/9``, ... .  ``R_max`` is the rank
      cutoff of the ``P^2`` exceptional-bundle enumeration.
    * On a Hirzebruch surface ``F_e`` with a **strictly ample** ``H`` this is the
      Coskun-Huizenga Drezet-Le Potier envelope ``DLP_H(nu)``, computed natively from the
      ``mu_H``-stable exceptional bundles of rank ``<= rank_max``
      (:func:`bridgeland_stability.dlp_hirzebruch.dlp_envelope`).  It is the **sharp**
      bound ``delta_H^{mu-s}(nu)`` when ``e in {0,1}`` and ``H`` is anticanonical
      (arXiv:1907.06739 Cor. "deltaDLP"), and a **certified lower bound** for every other
      ample ``H`` (Cor. "deltaDLPe").  This retires the G18 remainder.
    * Otherwise (K3, abelian, or an ``F_n`` carrying only the nef-and-big factory
      polarization) there is no envelope and this returns the **Bogomolov floor** ``0``.

    The G12 faithful-computation guard is applied first: a torsion-canonical surface
    (Enriques / bielliptic) is refused with the NS-lattice-refactor error rather than
    silently mis-modelled.  Always exact (``Fraction``).
    """
    require_faithful_computation(surface)  # G12 guard: torsion-canonical rows refused
    if surface.is_p2:
        mu = _mu(xi, surface)
        R_max = certified_rank_cutoff(mu, R_max)   # shared with dlp.moduli_nonempty (A4b)
        bundles = enumerate_exceptional(mu - 3, mu + 3, R_max)
        return delta(mu, bundles)
    env = _hirzebruch_envelope(xi, surface, rank_max)
    return Fraction(0) if env is None else env.value


def _certified_sharp_bound(
    xi: SurfaceBundle, surface: Surface, R_max: int, rank_max: int
) -> Optional[Fraction]:
    """The package's OWN theorem-certified sharp ``delta_H`` for ``xi``, or ``None``.

    A sharp bound is a THEOREM only where the package can prove one (E12-M4 gate):

    * **P^2** -- the Drezet-Le Potier closed form :func:`delta` is always the sharp curve
      ``delta(mu)`` (HN length one implicit), so this returns it unconditionally.
    * an **ample F_e with ``env.certified_sharp``** -- ``e in {0,1}`` and ``H``
      anticanonical, where CH Cor. "deltaDLP" gives ``delta_H^{mu-s}(nu) = DLP_{-K}(nu)``
      and the rank truncation is certified exact.  Returns ``env.value``.
    * everywhere else (a non-anticanonical ample ``F_e``, ``e >= 2``, K3, abelian, a
      nef-and-big ``F_n``) -- **``None``**: no theorem gives a sharp bound, so an external
      target cannot be verified against one and must be refused (invariant 7).

    Byte-identical to the ``delta_H`` the accepted certified-target path already uses (both
    read the same DLP curve / envelope value), so honouring a target that equals this value
    changes no verdict.  Always exact (``Fraction``); never a float.
    """
    if surface.is_p2:
        return delta_H(xi, surface, R_max)              # DLP closed form: always sharp
    env = _hirzebruch_envelope(xi, surface, rank_max)
    if env is not None and env.certified_sharp:
        return env.value                                # CH Cor. deltaDLP: sharp on the -K del Pezzo ray
    return None


def _is_integral_c1(xi: SurfaceBundle) -> bool:
    return all(c.denominator == 1 for c in xi.c1)


def validate_character(
    r: int, c1: Sequence[Number], ch2: Number, surface: Surface
) -> bool:
    """Thm 2.2 integrality: the character is integral on EVERY surface, or ``M(xi)`` is empty.

    A coherent sheaf has integral Chern classes -- ``r in Z``, ``c1 in NS`` and
    ``c2 = 1/2 <c1,c1> - ch2 in Z`` (``c2 in H^4(X,Z) = Z``) -- so a character failing any of
    these is not the Chern character of any sheaf and ``M(xi)`` is trivially empty.  All three
    are checked here, on every surface:

    * ``r >= 1`` and every ``c1`` coordinate integral (``c1 = r*mu in Z``, Thm 2.2);
    * ``c2 = 1/2 <c1,c1> - ch2 in Z`` via the NS self-pairing (``surface.lattice`` -- a rank-1
      shim on ``P^2``).  This is **K_X-independent** (it never touches the canonical class), so
      it is exact off ``P^2`` and does NOT wait for the ``Surface.K_H`` repair (A8 / E12-M6).
      Together with ``c1``-integrality it forces ``chi in Z`` by Riemann-Roch, because
      ``c1.(c1 - K_X)`` is even on any surface (Wu); so ``(c1, c2)``-integrality is the WHOLE
      Thm 2.2 integrality clause, not a fragment of it.

    On ``P^2`` the equivalent literal ``chi = r(P(mu) - Delta) = ch2 + 3 c1/2 + r`` is
    additionally checked, reproducing the oracle's ``_chi`` bit-for-bit
    (tests/oracle/dlp_reference.py); there ``chi in Z <=> c2 in Z``, so it is redundant with the
    ``c2`` test above but kept as the theorem's own form.
    """
    if r < 1:
        return False
    c1f = tuple(Fraction(x) for x in c1)
    if any(x.denominator != 1 for x in c1f):
        return False
    # c2 = 1/2 <c1,c1> - ch2 must be an integer (Chern classes are integral).  K_X-independent,
    # so this catches the off-P^2 A3 forge -- a non-integral c2 the P^2-only chi clause below
    # never sees -- without the deferred K_H repair.
    c2 = Fraction(1, 2) * surface.lattice.self_pairing(c1f) - Fraction(ch2)
    if c2.denominator != 1:
        return False
    if surface.is_p2:
        xi = SurfaceBundle(r, c1, ch2)
        mu = _mu(xi, surface)                      # = c1/r on P^2 (d = 1)
        chi = r * (P(mu) - discriminant(xi, surface))
        if chi.denominator != 1:
            return False
    return True


def _hirzebruch_verdict(
    xi: SurfaceBundle, surface: Surface, disc: Fraction, env, rank_max: int
) -> NonemptinessVerdict:
    """The E11-M6 / G18 verdict on an ample-polarized Hirzebruch surface.

    An invalid character (non-integral ``c1`` or ``c2``) short-circuits FIRST to a PROVEN empty
    verdict: it is not the Chern character of any sheaf, so no exceptional-bundle branch below
    may forge a PROVEN non-empty for it (A3, off ``P^2``).  The check is ``K_X``-independent
    (:func:`validate_character`'s ``c2`` clause), hence exact on ``F_e`` today.  Then four
    disjoint certified regimes, then an honest HEURISTIC remainder.

    1. ``Delta < 0`` -- Bogomolov: no ``mu_H``-semistable sheaf.  **PROVEN empty.**
    2. ``xi`` is a ``mu_H``-stable exceptional bundle, or ``V^{+m}`` for one.  Then
       ``M_H(xi)`` contains that (semi)stable sheaf.  **PROVEN non-empty**, even though
       such a class sits strictly BELOW the envelope -- the F_e analogue of the
       Drezet-Le Potier exceptional disjunct on ``P^2``.
    3. ``Delta < emptiness_bound(xi)`` -- some ``mu_H``-stable exceptional bundle ``V``
       forces ``chi(V, W) <= 0`` (or ``chi(W, V) <= 0``) on any semistable ``W`` of this
       slope, contradicting Riemann-Roch.  **PROVEN empty.**  Note ``emptiness_bound``
       is strictly weaker than the envelope: it drops the ``(nu - nu_V).H = 0``,
       ``nu != nu_V`` branch, which is not a theorem (see
       :func:`bridgeland_stability.dlp_hirzebruch.emptiness_bound`).
    4. ``Delta > delta_H`` with a **sharp and exact** envelope (``e in {0,1}``, ``H``
       anticanonical, truncation certified): ``mu_H``-stable sheaves exist by
       arXiv:1907.06739 Thm. "deltaSurface" (1) + Cor. "deltaDLP".  **PROVEN non-empty.**

    Everything else -- notably the boundary ``Delta == delta_H``, and every verdict under
    a merely lower-bounding envelope -- returns the inequality's truth value with
    ``HEURISTIC`` rigor.  The boundary case is a recorded open question: Thm.
    "deltaSurface" (1) needs a strict inequality.
    """
    mode = HNMode.DLP_ANTICANONICAL if env.sharp else HNMode.DLP_LOWER
    dH = env.value

    # A3 (off P^2): reject a non-integral character BEFORE any exceptional-bundle branch can
    # forge a PROVEN non-empty verdict for a class that is not the Chern character of any sheaf.
    # validate_character's c2 clause is K_X-independent, so this is exact on F_e today.
    if not validate_character(xi.r, xi.c1, xi.ch2, surface):
        return NonemptinessVerdict(
            False, disc, dH, mode, _INVALID_CHARACTER_CERT,
            f"mode={mode.value}: invalid Chern character "
            f"(c1 not in NS, or c2 = 1/2<c1,c1> - ch2 not integral): no coherent sheaf has "
            f"this character -> PROVEN empty", False)

    integral = _is_integral_c1(xi)
    c1i = tuple(int(c) for c in xi.c1) if integral else None

    # A6 (E12-M3): `is_stable_exceptional` is a function of (r, c1) ONLY; a class is a genuine
    # mu_H-stable exceptional bundle only if its ch2 is ALSO the forced exceptional value
    # (ch2 == exceptional_ch2, i.e. Delta == Delta_V).  Without this guard the `exceptional or`
    # short-circuit forged a PROVEN non-empty for an impostor with the right (r, c1) but wrong ch2.
    exceptional = bool(
        integral
        and is_stable_exceptional(xi.r, c1i, surface)
        and xi.ch2 == exceptional_ch2(xi.r, c1i, surface)
    )
    semiexceptional = bool(integral and (exceptional or is_semiexceptional(xi, surface)))

    def _cert(rigor, hyps, note):
        return Certificate(rigor, tuple(hyps), (_CH,), note)

    if disc < 0:
        return NonemptinessVerdict(
            False, disc, dH, mode,
            _cert(Rigor.PROVEN, ("Bogomolov inequality: Delta >= 0 for any mu_H-semistable sheaf",),
                  "Delta < 0: no semistable sheaf of this character exists, for any polarization."),
            f"mode={mode.value}: Delta={disc} < 0 violates Bogomolov -> PROVEN empty", False)

    if semiexceptional:
        return NonemptinessVerdict(
            True, disc, dH, mode,
            _cert(Rigor.PROVEN,
                  ("xi is (a direct sum of copies of) a mu_H-stable exceptional bundle",
                   "Cor. 'DLPExceptional': Delta >= DLP_H^{<r}(nu) certifies mu_H-stability"),
                  "exceptional/semiexceptional bundle: a non-empty point that may sit "
                  "strictly BELOW the envelope (F_e analogue of the DLP disjunct)."),
            f"mode={mode.value}: Delta={disc}, delta_H={dH}; exceptional bundle: "
            f"non-empty {'point' if exceptional else 'semiexceptional class'} below the envelope",
            True)

    eb = emptiness_bound(xi, surface, rank_max)
    if disc < eb:
        return NonemptinessVerdict(
            False, disc, dH, mode,
            _cert(Rigor.PROVEN,
                  ("a mu_H-stable exceptional bundle V with 0 < |(nu-nu_V).H| <= -1/2 K.H "
                   "(or nu = nu_V with Delta != Delta_V) forces Delta >= DLP_{H,V}(nu)",
                   "Sec. 5.4: Hom and Ext^2 vanish by stability + Serre duality"),
                  "Delta below the certified emptiness bound: no Gieseker-semistable sheaf."),
            f"mode={mode.value}: Delta={disc} < emptiness_bound={eb} -> PROVEN empty", False)

    if env.certified_sharp and disc > dH:
        return NonemptinessVerdict(
            True, disc, dH, mode,
            _cert(Rigor.PROVEN,
                  ("e in {0,1} and H anticanonical: delta_H^{mu-s}(nu) = DLP_{-K}(nu) "
                   "(Cor. 'deltaDLP')",
                   "the rank truncation is certified exact (DLP_{-K,V} <= 1/2 + 1/(2 r_V^2))",
                   "Thm. 'deltaSurface' (1): Delta > delta_H^{mu-s} gives mu_H-stable sheaves"),
                  "Delta strictly above the sharp anticanonical envelope."),
            f"mode={mode.value}: Delta={disc} > delta_H={dH} (sharp, exact) -> PROVEN nonempty",
            False)

    above = disc >= dH
    why = ("the envelope is a certified LOWER bound only (H is not an anticanonical del "
           "Pezzo polarization; delta_H may be computed by Kronecker modules)"
           if not env.sharp else
           "Delta sits ON the sharp envelope (Thm. 'deltaSurface' (1) needs a strict "
           "inequality) or the rank truncation is not certified exact")
    return NonemptinessVerdict(
        above, disc, dH, mode,
        Certificate(Rigor.HEURISTIC, (why,), (_CH,),
                    "envelope comparison without a certified disjunct: HEURISTIC."),
        f"mode={mode.value}: Delta={disc} {'>=' if above else '<'} delta_H={dH} "
        f"(HEURISTIC: {why})", False)


def moduli_nonempty(
    r: int,
    c1: Sequence[Number],
    ch2: Number,
    surface: Surface,
    *,
    delta_H_target: Optional[Number] = None,
    hn_source: Optional[HNMode] = None,
    evidence: Optional["SharpBoundEvidence"] = None,
    R_max: int = 60,
    rank_max: int = DEFAULT_RANK_MAX,
) -> NonemptinessVerdict:
    """Decide ``Delta(xi) >= delta_H(xi)`` for ``xi = (r, c1, ch2)`` on ``surface``.

    GIVEN the HN-length-one datum (implicit on ``P^2`` via the DLP closed form,
    or supplied by a certified ``PAPER`` / ``ORACLE`` source), evaluate the sharp
    inequality exactly in ``Fraction`` and return a :class:`NonemptinessVerdict`
    carrying a G5 :class:`~bridgeland_stability.rigor.Certificate`.

    Mode selection
    --------------
    * ``evidence`` or ``delta_H_target`` given -- an externally supplied sharp bound
      (the E11-M4 paper-table / E11-M5 oracle entry point).  E12-M4 gate (A5): the
      bound is honoured (-> PROVEN-eligible) **only** when it equals the package's OWN
      theorem-certified sharp bound for this exact class (:func:`_certified_sharp_bound`
      -- P^2, or an ample ``F_e`` with ``env.certified_sharp``).  A ``delta_H_target``
      still needs a certified ``hn_source`` (``DLP`` / ``PAPER``); a raw
      ``hn_source=ORACLE`` pair is refused outright (an ORACLE datum is the capability
      object :func:`bridgeland_stability.oracle.mint_oracle_evidence` returns, minted
      only after a real construction succeeded).  A bound bound to a different class, or
      unequal to the certified one, or on a surface with no certified sharp bound, is a
      ``ValueError`` -- never a forged PROVEN verdict (invariant 7).
    * else on ``P^2`` -- ``DLP`` mode, the closed-form ``delta_H`` (PROVEN).
    * else (off ``P^2``, no target) -- the native Coskun-Huizenga envelope verdict
      (:func:`_hirzebruch_verdict` on an ample ``F_e``; the ``HEURISTIC`` Bogomolov
      floor ``0`` on a K3 / abelian / nef-and-big ``F_n``).  A bare certified
      ``hn_source`` (no ``delta_H_target``) does **not** upgrade this: it certifies
      only the HN-length-one hypothesis, not a sharp ``delta_H``, so it is honoured
      (``PROVEN``) solely where the native envelope is itself certified sharp
      (``e in {0,1}``, ``H`` anticanonical) -- gated per-branch inside
      ``_hirzebruch_verdict``, never by the source label.  The E11-M5 oracle must
      pass its sharp ``delta_H`` as ``delta_H_target`` to reach ``PROVEN`` off ``P^2``.

    The verdict's certificate is ``PROVEN`` only for a certified mode; otherwise
    ``HEURISTIC``.  The G12 faithful-computation guard runs first, and every path enforces
    Thm 2.2 character integrality (:func:`validate_character`) -- an invalid ``c1``/``c2`` is
    ``PROVEN`` empty on **every** surface, never a forged non-empty (A3).

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
    disc = discriminant(xi, surface)       # full-NS CH discriminant (G18a)

    if evidence is not None or delta_H_target is not None:
        # E12-M4 (A5): an external sharp bound is honoured only when it is (1) bound to
        # THIS class and (2) equal to the package's OWN certified sharp bound for it.
        if evidence is None:
            # Legacy (delta_H_target, hn_source) entry point -> wrap it as evidence.
            if hn_source not in _CERTIFIED:
                raise ValueError(
                    "delta_H_target must come with a certified hn_source "
                    "(HNMode.DLP / HNMode.PAPER / HNMode.ORACLE): an externally "
                    "supplied sharp delta_H is only PROVEN when its HN-length-one "
                    "hypothesis is certified (E11-M4 paper table / E11-M5 oracle)."
                )
            if hn_source is HNMode.ORACLE:
                # A5: a raw ORACLE target is a forge.  An ORACLE datum is a capability
                # object minted inside bridgeland_stability.oracle after a construction
                # actually returned True -- never a bare (delta_H_target, ORACLE) pair.
                raise ValueError(
                    "an ORACLE HN datum is a capability object minted via "
                    "bridgeland_stability.oracle.mint_oracle_evidence after a construction "
                    "actually returned True (E12-M4); it cannot be passed as a raw "
                    "(delta_H_target, hn_source=ORACLE) pair"
                )
            evidence = SharpBoundEvidence(
                surface=surface, r=r, c1=tuple(xi.c1), ch2=xi.ch2,
                sharp_bound=Fraction(delta_H_target),
                sharp_bound_source=hn_source, hn_length_one_source=hn_source,
                citation="caller-supplied legacy delta_H_target",
            )
        # (1) class-bound: evidence derived for another class cannot certify this one.
        if not evidence.matches(r, tuple(xi.c1), xi.ch2, surface):
            raise ValueError(
                "evidence was derived for a different class than the one queried "
                f"(evidence: r={evidence.r}, c1={tuple(evidence.c1)}, ch2={evidence.ch2}, "
                f"surface={evidence.surface.name!r}; queried: r={r}, c1={tuple(xi.c1)}, "
                f"ch2={xi.ch2}, surface={surface.name!r})"
            )
        # (2) value gate: the supplied bound must equal the package's OWN certified sharp
        # bound for this class.  Where no theorem certifies one, an external target is
        # unverifiable and is refused rather than trusted (invariant 7).
        native = _certified_sharp_bound(xi, surface, R_max, rank_max)
        if native is None:
            raise ValueError(
                "no theorem-certified sharp bound on this surface, so an external target "
                "cannot be verified and is refused (invariant 7): the package certifies a "
                "sharp bound only on P^2 (Drezet-Le Potier closed form) and on an ample "
                "F_e with e in {0,1} and H anticanonical (CH Cor. deltaDLP)."
            )
        if evidence.sharp_bound != native:
            raise ValueError(
                f"delta_H_target {evidence.sharp_bound} != the package's own certified "
                f"sharp bound {native} for this class: a forged / mis-derived target cannot "
                "mint a PROVEN verdict (A5)"
            )
        dH, mode = evidence.sharp_bound, evidence.hn_length_one_source
    elif surface.is_p2:
        dH, mode = delta_H(xi, surface, R_max), HNMode.DLP
    else:
        env = _hirzebruch_envelope(xi, surface, rank_max)
        if env is not None:
            # Native Coskun-Huizenga envelope verdict.  _hirzebruch_verdict validates the
            # character first (A3, so an invalid c1/c2 cannot reach a PROVEN non-empty branch)
            # and gates PROVEN on env.certified_sharp PER BRANCH -- so a merely lower-bounding
            # envelope (any ample H off the anticanonical del Pezzo ray, or e >= 2) never
            # certifies existence.
            #
            # A bare certified hn_source (ORACLE/PAPER/DLP with NO delta_H_target) is
            # deliberately NOT honoured on this path.  It certifies only the HN-length-one
            # hypothesis, not a sharp delta_H; without an externally supplied sharp target the
            # package has ONLY its own envelope to compare against.  Previously this path stamped
            # _MODE_CERT[hn_source] = PROVEN over env.value while env.sharp was False, forging a
            # false PROVEN_NONEMPTY for every class in the gap [env.value, sharp delta_H) -- a
            # class that is EMPTY reported as PROVEN non-empty (invariant 7's worst outcome),
            # reachable through the documented E11-M5 ORACLE hook with one public call.  The
            # oracle MUST carry its sharp delta_H as delta_H_target (the first branch above) to
            # reach PROVEN off P^2.  See docs/CORRECTIONS.md E12-M2 (IMPROVE round 4).
            return _hirzebruch_verdict(xi, surface, disc, env, rank_max)
        # env is None: K3 / abelian / a nef-and-big (non-ample) F_n -- no CH envelope at all,
        # so delta_H falls back to the Bogomolov floor 0, which is never sharp.  A bare certified
        # hn_source cannot certify existence against the floor either, so the verdict stays
        # HEURISTIC (the E11-M5 oracle hook reaches PROVEN only via a real delta_H_target).
        dH = delta_H(xi, surface, R_max, rank_max)   # = 0 (Bogomolov floor)
        mode = HNMode.HEURISTIC

    cert = _MODE_CERT[mode]
    valid = validate_character(r, xi.c1, xi.ch2, surface)
    above_curve = disc >= dH
    # DLP exceptional-bundle disjunct on EVERY surface (the pinned P^2 detectors on P^2, the
    # ch2-guarded F_e detectors off P^2) -- so an exceptional/semiexceptional class fed its
    # OWN correct sharp delta_H target is still reported non-empty here, never a false PROVEN
    # "empty" for a class whose moduli space is a single reduced point (E12-M2).  Mirrors
    # _hirzebruch_verdict's native disjunct, which the certified-target path previously
    # dropped off P^2 -> a PROVEN_EMPTY contradicting the same function's native verdict.
    exceptional, semiexceptional = _exceptional_disjunct(xi, surface)
    nonempty = valid and (above_curve or exceptional or semiexceptional)
    # Branch-derived rigor off P^2.  The CH non-emptiness theorem (Thm "deltaSurface" (1))
    # needs a STRICT inequality Delta > delta_H, and its converse "Delta < delta_H => empty"
    # is a theorem only BELOW the certified emptiness_bound -- which is strictly weaker than the
    # envelope (CLAUDE.md invariant; emptiness_bound drops the non-theorem (nu-nu_V).H=0 branch).
    # So, mirroring _hirzebruch_verdict, downgrade the certificate to HEURISTIC (-> UNKNOWN) for
    # the whole band emptiness_bound <= Delta <= delta_H when no exceptional disjunct fires.  A
    # flat "Delta < delta_H => PROVEN empty" over-claimed emptiness in the gap emptiness_bound <=
    # Delta < delta_H that the SAME function's native envelope verdict reports UNKNOWN -- a
    # non-theorem-backed false PROVEN_EMPTY (E12-M2, IMPROVE round 3).  Below emptiness_bound the
    # verdict stays PROVEN empty (theorem); strictly above delta_H it stays PROVEN nonempty; the
    # exceptional/semiexceptional disjunct proves non-emptiness independently of the envelope and
    # so is never downgraded.  Off a Hirzebruch F_e (K3, abelian, nef-and-big F_n) there is no
    # emptiness_bound theory, so only the boundary Delta == delta_H is downgraded, as before.
    # An invalid character is not the Chern character of ANY coherent sheaf, so M(xi)
    # is empty on every surface for every polarization -- a K_X-independent theorem
    # (validate_character's c2 clause).  Swap in the branch-derived PROVEN-empty cert
    # so the verdict is PROVEN_EMPTY uniformly, matching _hirzebruch_verdict; without
    # this the K3/abelian path (mode=HEURISTIC) under-claimed it as UNKNOWN while
    # P^2/F_e reported PROVEN_EMPTY -- an inconsistency the E12 code review flagged.
    band_unknown = False
    if not valid:
        cert = _INVALID_CHARACTER_CERT
    elif (mode in _CERTIFIED and not surface.is_p2
            and not (exceptional or semiexceptional)):
        eb = _fe_emptiness_bound(xi, surface, rank_max)
        in_band = (disc == dH) or (eb is not None and eb <= disc <= dH)
        if in_band:
            cert = _BOUNDARY_CERT
            band_unknown = True
    reason = f"mode={mode.value}: Delta={disc} {'>=' if above_curve else '<'} delta_H={dH}"
    if not valid:
        reason += "; invalid Chern character (c1 or chi not integral): no sheaf exists"
    elif exceptional and not above_curve:
        reason += "; exceptional bundle: non-empty isolated point below the DLP curve"
    elif semiexceptional and not above_curve:
        reason += "; semiexceptional (m*ch(E)): non-empty polystable point below the DLP curve"
    elif band_unknown:
        reason += ("; emptiness_bound <= Delta <= delta_H: emptiness is not a theorem in this "
                   "band (the envelope is strictly stronger than emptiness_bound) -> UNKNOWN")
    if mode not in _CERTIFIED:
        reason += " (HEURISTIC: Bogomolov floor; no Drezet-Le Potier envelope applies)"
    return NonemptinessVerdict(nonempty, disc, dH, mode, cert, reason, exceptional)
