"""E13-M3a / G18: the HN-length-one existence criterion + Thm 1.13 structure.

The Coskun-Huizenga non-emptiness program is built toward a single theorem
(arXiv:1907.06739 Sec. 1.6, verbatim):

    "there exists an H_m-semistable sheaf with Chern character v if and only if
     the generic H_m-Harder-Narasimhan filtration has length 1."

That HN-length-one datum is exactly what the E11-M3 numerical evaluator
*delegated* ("genuinely sheaf-theoretic, not pure Chern arithmetic").  This
module ships it as the sharp-verdict interface on the **del Pezzo Hirzebruch
surfaces** ``e in {0,1}`` (plus ``P^2``), where the generic HN length is already
determined by the shipped :func:`bridgeland_stability.nonemptiness_rational.
moduli_nonempty` verdict:

* generic HN length ``1`` (region **S**) iff the class is above the sharp
  envelope (``Delta > delta_H^{mu-s}``) OR is (semi)exceptional -- a non-empty
  point below the envelope, the ``F_e`` analogue of the Drezet-Le Potier
  disjunct.  ``moduli_nonempty`` reports this as ``PROVEN_NONEMPTY``.
* generic HN length ``>= 2`` with a certified obstruction (region **EMPTY**):
  ``Delta < 0`` (Bogomolov), or ``Delta`` below the certified ``emptiness_bound``.
  ``moduli_nonempty`` reports this as ``PROVEN_EMPTY``.
* the remaining band ``emptiness_bound <= Delta <= delta_H`` (and the boundary
  ``Delta == delta_H``) -- where ``moduli_nonempty`` honestly reports UNKNOWN --
  is now **DECIDED (E13-M3b)** by computing the generic ``H_m``-Harder-Narasimhan
  filtration (:mod:`bridgeland_stability.generic_hn`, arXiv:1907.06739 Sec. 5):
  length 1 -> region S (semistable sheaves exist); length >= 2 -> PROVEN empty,
  with the factor characters EXHIBITED in :attr:`HNVerdict.factors` and region
  **K** earned exactly when a factor is not semiexceptional (the Kronecker-type
  datum); an empty prioritary stack -> PROVEN empty (prop-ssPrior).  The E13
  re-audit R2 discipline stands: every structural label is now backed by a
  computed filtration, never asserted from an epistemic UNKNOWN.

Faithful reframing, no re-derived logic
---------------------------------------
This module does NOT re-implement the envelope / verdict logic.  It builds the
character ``(r, c1, ch2)`` from ``(r, nu, Delta)`` exactly, calls the native
``moduli_nonempty`` (no ``delta_H_target`` / ``evidence``), and maps its
branch-derived :class:`~bridgeland_stability.nonemptiness_rational.VerdictStatus`:

    ================  =================  ============  ==============================
    moduli status     semistable_exists region        generic HN length
    ================  =================  ============  ==============================
    PROVEN_NONEMPTY   True               S             1
    PROVEN_EMPTY      False              EMPTY         None (>= 2; not computed here)
    UNKNOWN           decided by M3b     S / K / EMPTY exact computed length (1..4)
    ================  =================  ============  ==============================

Delegating guarantees zero drift: it inherits ``moduli_nonempty``'s A3/A5/A6 /
E12-M2 soundness fixes (the strict ``>`` vs ``>=`` boundary handling and the
``emptiness_bound`` band) for free.  Re-deriving the ``>`` boundary would
re-introduce exactly the over-claim the E12 audit closed.

Thm 1.13 = Cor 7.7 structure (arXiv:1907.06739 Sec. 7; Example 1.14)
--------------------------------------------------------------------
For ``e in {0,1}``, ``Delta >= 3/8``, and ``H`` sufficiently close to ``-K``: if
there are no ``H``-semistable sheaves then **at most one Harder-Narasimhan factor
of the general prioritary sheaf is not a semiexceptional bundle**.  That one
non-semiexceptional factor (when present) is a **Kronecker module**.  Note the
theorem bounds the number of NON-semiexceptional factors -- it does NOT assert
an exactly length-two filtration (E13 re-audit R2).  The threshold
``3/8`` is pinned as :data:`THM_1_13_MIN_DELTA` (two-way: ``3/8 =
exceptional_discriminant(2) = (1 - 1/2^2)/2``), and the "at most one factor"
shape as :data:`THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS`.

Honest scope
------------
* ``P^2`` -- the Drezet-Le Potier closed form is sharp everywhere, so the
  **existence boolean** is total: every character is decided (``True`` or
  ``False``, never ``None``).  Totality does NOT extend to the generic-HN
  *shape*: Example 1.14's S/K/R/empty shapes all occur on ``P^2`` too, and this
  module's region label describes the package's verdict, never the sheaf's HN
  structure (E13 re-audit R2).
* an ample-polarized del Pezzo ``F_0`` / ``F_1`` -- **TOTAL since E13-M3b**:
  the envelope decides region S and the certified-empty regions; the remaining
  band is decided by the computed generic HN filtration.  The flagship
  ``(2,(1,1),0)`` on ``F_0`` = ``O(1,0) (+) O(0,1)`` -- the E13 re-audit R2
  counterexample that forced the honest UNCLASSIFIED label -- now DECIDES to
  ``exists=True`` with computed length 1, matching the polystable truth (the
  paper exhibits exactly this bundle as ``-K``-semistable with
  ``Delta = 1/4 < DLP_{-K}``, the Sec. 7 example after cor-delPezzoKss).
* ``e >= 2`` -- **not** in M3a/M3b verdict scope: assembled via the E13-M1
  reduction ``pi`` in E13-M3c (the Sec. 5 algorithm itself is uniform in ``e``
  and :mod:`bridgeland_stability.generic_hn` accepts any ample ``F_e``; M3c
  will cross-check it against the reduction).  A :exc:`NotImplementedError`.
* K3 / abelian / a nef-and-big factory ``F_n`` -- no del Pezzo CH theory; refused.

E13-M3a partially closed the E11-M6 open question O2; E13-M3b closes it on the
del Pezzo scope (docs/CORRECTIONS.md Sec. 11/14): even the boundary
``Delta == delta_H`` is decided by the computed filtration.

Exact arithmetic; full-NS discriminant
--------------------------------------
Every quantity is ``fractions.Fraction``; no float ever appears.  ``Delta`` is
the **full-NS** Coskun-Huizenga discriminant ``Delta = 1/2 <nu,nu> - ch2/r``
(CLAUDE.md invariant 2), NEVER the H-projected surrogate ``discriminant_H``.
Stdlib-only at import time (``fractions``, ``enum``, ``dataclasses``, ``typing``
+ the stdlib-only intra-package ``nonemptiness_rational`` / ``dlp_hirzebruch`` /
``exceptional_surface`` / ``varieties`` / ``rigor``), preserving the
zero-runtime-dependency invariant.

References
----------
* Coskun-Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
  arXiv:1907.06739 -- Sec. 1.6 (the HN-length-one theorem), Sec. 5 (the generic
  HN filtration; Thm 1.6), Sec. 7 (Thm 1.13 = Cor 7.7), Example 1.9, Example 1.14.
* Drezet-Le Potier, Ann. Sci. ENS 18 (1985) -- the ``P^2`` closed-form
  ``delta``-curve underlying the total ``P^2`` verdict.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Optional, Sequence, Union

from .nonemptiness_rational import moduli_nonempty, VerdictStatus
from .dlp_hirzebruch import hirzebruch_index, is_ample, is_semiexceptional
from .exceptional_surface import SurfaceBundle
from .generic_hn import generic_hn_factors
from .rigor import Certificate, Rigor
from .varieties import Surface

__all__ = [
    "HNRegion",
    "HNVerdict",
    "hn_verdict",
    "semistable_exists",
    "generic_hn_length",
    "hn_region",
    "THM_1_13_MIN_DELTA",
    "THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS",
]

Number = Union[int, Fraction]

# --- Thm 1.13 = Cor 7.7 structural constants (arXiv:1907.06739 Sec. 7) --------
#: The Thm 1.13 discriminant threshold ``Delta >= 3/8``.  Two-way (CLAUDE.md
#: invariant 3): ``3/8 == exceptional_discriminant(2) == (1 - 1/2^2)/2`` -- the
#: rank-2 exceptional discriminant, the smallest ``Delta_V`` above the rank-1
#: floor ``0``.  Pinned as a literal from the paper; the test cross-checks it
#: against the package's own ``exceptional_discriminant(2)``.
THM_1_13_MIN_DELTA: Fraction = Fraction(3, 8)
#: Thm 1.13 shape: "at most one HN factor of the general prioritary sheaf is not
#: a semiexceptional bundle".  An upper bound on the non-semiexceptional factor
#: count ONLY -- it neither fixes the filtration length at two nor guarantees a
#: Kronecker factor exists (E13 re-audit R2).
THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS: int = 1


class HNRegion(str, Enum):
    """The verdict regions on a del Pezzo ``e in {0,1}`` (plus ``P^2``).

    Since E13-M3b the verdict is TOTAL on an ample ``F_0`` / ``F_1``: the band
    E13-M3a honestly reported UNKNOWN is now decided by COMPUTING the generic
    ``H_m``-Harder-Narasimhan filtration (:mod:`bridgeland_stability.generic_hn`,
    arXiv:1907.06739 Sec. 5), and the region label is EARNED, never asserted
    (E13 re-audit R2 discipline):

    * ``S`` -- the computed filtration has length 1: semistable sheaves exist.
    * ``K`` -- PROVEN empty with a computed length->=2 filtration at least one
      of whose factors is NOT semiexceptional (the Kronecker-type datum of
      Ex. 1.14 / Sec. 1.5's orthogonal-Kronecker examples); the factors are
      exhibited in :attr:`HNVerdict.factors`.
    * ``EMPTY`` -- certified empty with no Kronecker datum: invalid character,
      ``Delta < 0``, below ``emptiness_bound``, an empty prioritary stack, or a
      computed filtration all of whose factors are semiexceptional (the rigid
      "R"-type decomposition of Ex. 1.9).
    * ``UNCLASSIFIED`` -- retained as the honest fallback label; never returned
      on the M3a/M3b scope (P^2 / ample del Pezzo ``F_0``, ``F_1``).
    """

    S = "S"          # generic HN filtration length 1: semistable sheaves exist
    K = "K"          # PROVEN length->=2 filtration with a non-semiexceptional factor (M3b, exhibited)
    UNCLASSIFIED = "unclassified"  # honest-fallback label; never returned on the M3a/M3b scope
    EMPTY = "empty"  # certified empty: invalid, Delta<0, emptiness_bound, no prioritary stack, or all-semiexceptional factors


@dataclass(frozen=True)
class HNVerdict:
    """The generic-HN-length verdict for a class ``(r, nu, Delta)`` on ``surface``.

    Attributes
    ----------
    exists : bool | None
        ``True`` (region S, HN length 1, semistable sheaves exist) or ``False``
        (certified empty).  TOTAL on the M3a/M3b scope since E13-M3b; the
        ``None`` state is retained in the type only as the honest fallback and
        is never returned here.  Never a fabricated verdict.
    generic_hn_length : int | None
        ``1`` where the generic HN filtration is PROVEN length-1; the exact
        computed length (``2..4``) where the M3b algorithm decided emptiness;
        ``None`` on the envelope-decided PROVEN_EMPTY paths (empty, but the
        filtration was not computed).
    region : HNRegion
        ``S`` / ``K`` / ``EMPTY`` -- see :class:`HNRegion`; ``K`` is EARNED by a
        computed filtration with a non-semiexceptional factor.
    discriminant : Fraction
        The full-NS Coskun-Huizenga ``Delta(xi)`` the verdict compared.
    sharp_bound : Fraction
        The ``delta_H`` the underlying ``moduli_nonempty`` verdict used (on the
        anticanonical del Pezzo ray this is the sharp ``dlp_envelope`` value).
    certificate : Certificate
        The G5 provenance stamp: PROVEN everywhere on the M3a/M3b scope (the
        M3b branches carry the Sec. 5 algorithm citations).
    reason : str
        A human-readable one-line summary.
    factors : tuple | None
        E13-M3b: the computed generic-HN factor characters, where the verdict
        came from the algorithm (see the field comment below).
    """

    exists: Optional[bool]
    generic_hn_length: Optional[int]
    region: HNRegion
    discriminant: Fraction
    sharp_bound: Fraction
    certificate: Certificate
    reason: str
    #: E13-M3b: the COMPUTED characters of the generic H_m-Harder-Narasimhan
    #: factors ((r, c1, ch2) triples, ordered), when the verdict came from the
    #: Sec. 5 algorithm; ``None`` on the envelope-decided paths (where the
    #: verdict predates the filtration computation) and on P^2.
    factors: Optional[tuple] = None


def _require_del_pezzo_scope(surface: Surface) -> None:
    """M3a scope guard: ``P^2`` (fully sharp DLP) or an ample ``F_0`` / ``F_1``.

    * ``P^2`` -- allowed.
    * ``F_e`` with ``e in {0,1}`` and a strictly ample ``H`` -- allowed.
    * ``F_e`` with ``e in {0,1}`` but only the nef-and-big factory ``H`` -- a
      :exc:`ValueError` (use :func:`~bridgeland_stability.nonemptiness_rational.
      hirzebruch_with_polarization`).
    * ``e >= 2`` -- a :exc:`NotImplementedError` (E13-M3c assembles it via the
      reduction ``pi``).
    * K3 / abelian / any non-``F_e`` surface -- ``hirzebruch_index`` raises
      :exc:`NotImplementedError`, propagated here.
    """
    if surface.is_p2:
        return
    e = hirzebruch_index(surface)            # NotImplementedError off F_e (K3/abelian/...)
    if e in (0, 1):
        if not is_ample(surface):
            raise ValueError(
                f"{surface.name}: needs a strictly ample H on F_{e} "
                "(use nonemptiness_rational.hirzebruch_with_polarization); the factory "
                "H is only nef-and-big")
        return
    raise NotImplementedError(
        f"e={e} >= 2 is E13-M3c (assemble via the reduction pi: reduce_to_del_pezzo); "
        "E13-M3a is del Pezzo e in {0,1} plus P^2")


def hn_verdict(r: int, nu: Sequence[Number], Delta: Number, surface: Surface) -> HNVerdict:
    """The generic-HN-length verdict for ``(r, nu, Delta)`` on a del Pezzo ``surface``.

    Builds the character ``(r, c1 = r*nu, ch2 = r*(1/2 <nu,nu> - Delta))`` exactly
    -- an exact round-trip, so ``discriminant(xi) == Delta`` -- and delegates to the
    native :func:`~bridgeland_stability.nonemptiness_rational.moduli_nonempty` (no
    ``delta_H_target`` / ``evidence``).  The verdict's branch-derived
    :class:`~bridgeland_stability.nonemptiness_rational.VerdictStatus` maps to the
    tri-state ``(exists, generic_hn_length, region)`` per the module table.

    Raises per :func:`_require_del_pezzo_scope` off M3a scope (``e >= 2`` /
    K3 / a nef-and-big factory ``F_n``).
    """
    _require_del_pezzo_scope(surface)
    # E13 re-audit R3: never int()-truncate the rank.  A non-integral r is not the
    # Chern character of any sheaf; passing it through lets moduli_nonempty return the
    # invalid-character PROVEN_EMPTY verdict -- a truncated r would silently answer for
    # a DIFFERENT character than the caller supplied.
    r_frac = Fraction(r)
    r = int(r_frac) if r_frac.denominator == 1 else r_frac
    Delta = Fraction(Delta)
    nu_f = tuple(Fraction(x) for x in nu)
    # Exact round-trip from (r, nu, Delta) back to a Chern character:
    #   c1 = r * nu ; ch2 = r * (1/2 <nu,nu> - Delta)  =>  discriminant(xi) == Delta.
    c1 = tuple(r * x for x in nu_f)
    ch2 = r * (Fraction(1, 2) * surface.lattice.self_pairing(nu_f) - Delta)
    v = moduli_nonempty(r, c1, ch2, surface)   # NATIVE path (no target / evidence)
    st = v.status
    factors: Optional[tuple] = None
    cert = v.certificate
    tail = v.reason
    if st is VerdictStatus.PROVEN_NONEMPTY:
        exists: Optional[bool] = True
        length: Optional[int] = 1
        region = HNRegion.S
    elif st is VerdictStatus.PROVEN_EMPTY:
        exists, length, region = False, None, HNRegion.EMPTY
    else:
        # UNKNOWN band -- E13-M3b: DECIDE it by computing the generic H_m-HN
        # filtration (arXiv:1907.06739 Sec. 5 / cor-algorithm).  UNKNOWN only
        # occurs off P^2, on the ample F_e the scope guard admitted, so the
        # algorithm's own preconditions hold; the character is valid (an
        # invalid one is PROVEN_EMPTY above).
        factors = generic_hn_factors(int(r), c1, ch2, surface)
        if factors is None:
            exists, length, region = False, None, HNRegion.EMPTY
            cert = Certificate(
                Rigor.PROVEN,
                ("the stack of F- and H_ceil(m)-prioritary sheaves of this character "
                 "is empty (thm-prioritaryNecessary / cor-prioritaryDelta)",
                 "a mu_{H_m}-semistable sheaf is H_ceil(m)-prioritary (prop-ssPrior), "
                 "so M_{H_m}(v) is empty"),
                ("arXiv:1907.06739",),
                "empty prioritary stack: no semistable sheaf (E13-M3b).")
            tail = "prioritary stack empty"
        elif len(factors) == 1:
            exists, length, region = True, 1, HNRegion.S
            cert = Certificate(
                Rigor.PROVEN,
                ("the generic H_m-Harder-Narasimhan filtration, computed by the "
                 "Sec. 5 inductive algorithm (thm-HNcriterion / cor-algorithm), "
                 "has length one",
                 "an H_m-semistable sheaf exists iff the generic HN filtration "
                 "has length one (Sec. 1.6)"),
                ("arXiv:1907.06739",),
                "computed generic HN length one: semistable sheaves exist (E13-M3b).")
            tail = "generic HN filtration computed: length 1"
        else:
            exists, length = False, len(factors)
            non_semiexc = [
                w for w in factors
                if not is_semiexceptional(SurfaceBundle(w[0], w[1], w[2]), surface)
            ]
            region = HNRegion.K if non_semiexc else HNRegion.EMPTY
            cert = Certificate(
                Rigor.PROVEN,
                ("the generic H_m-Harder-Narasimhan filtration, computed by the "
                 f"Sec. 5 inductive algorithm, has length {len(factors)} >= 2",
                 "an H_m-semistable sheaf exists iff the generic HN filtration "
                 "has length one (Sec. 1.6), so M_{H_m}(v) is empty",
                 f"factor characters (r, c1, ch2): {tuple(factors)!r}"),
                ("arXiv:1907.06739",),
                "computed generic HN length >= 2: no semistable sheaf; the factors "
                "are exhibited (E13-M3b).")
            tail = (f"generic HN filtration computed: length {len(factors)}, "
                    f"factors {tuple(factors)!r}"
                    + ("" if non_semiexc else " (all semiexceptional)"))
    shape = (
        "1 (semistable sheaves exist)" if exists
        else (f"{length} (destabilized)" if length is not None and length >= 2
              else (">=2 (destabilized)" if exists is False
                    else "undecided"))
    )
    reason = (
        f"region {region.value}: generic HN filtration length {shape}; {tail}"
    )
    return HNVerdict(exists, length, region, v.discriminant, v.delta_H, cert, reason, factors)


def semistable_exists(
    r: int, nu: Sequence[Number], Delta: Number, surface: Surface
) -> Optional[bool]:
    """``True`` (region S, HN length 1) / ``False`` (certified empty) / ``None`` (UNCLASSIFIED).

    The HN-length-one existence criterion (arXiv:1907.06739 Sec. 1.6): an
    ``H``-semistable sheaf of character ``(r, nu, Delta)`` exists iff the generic
    HN filtration has length 1.  TOTAL on the M3a/M3b scope: the former
    UNCLASSIFIED band is decided by the computed generic HN filtration
    (E13-M3b); ``None`` is never returned here.
    """
    return hn_verdict(r, nu, Delta, surface).exists


def generic_hn_length(
    r: int, nu: Sequence[Number], Delta: Number, surface: Surface
) -> Optional[int]:
    """The PROVEN generic HN length: ``1`` on region S; the exact computed
    length (``2..4``) where the E13-M3b algorithm decided emptiness; ``None``
    only on the envelope-decided PROVEN_EMPTY paths (empty, filtration not
    computed -- call :func:`bridgeland_stability.generic_hn.generic_hn_factors`
    for the factors there too).
    """
    return hn_verdict(r, nu, Delta, surface).generic_hn_length


def hn_region(r: int, nu: Sequence[Number], Delta: Number, surface: Surface) -> HNRegion:
    """The verdict region ``S`` / ``K`` / ``EMPTY`` for ``(r, nu, Delta)`` on ``surface``."""
    return hn_verdict(r, nu, Delta, surface).region
