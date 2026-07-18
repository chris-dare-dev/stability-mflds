"""E15-M4 / G18+: the Conjecture A-gated evaluator of ``delta_m^{mu-s}(nu)``.

arXiv:1907.06739 Sec. 1.5: an affirmative answer to the orthogonal-Kronecker
conjecture "will allow an exact inductive computation of delta_m^{mu-s}(nu)".
This module implements that computation CONDITIONALLY, carrying the G5
provenance lattice honestly:

    delta_conjectural(nu, m, surface) = max(
        the DLP part      -- dlp_envelope, certified machinery (E11-M6),
        the Kronecker part -- thm-deltaKronecker (E14-M2) evaluated over the
                              twist/swap orbit of the Sec. 8 collection family)

* ``delta_m^{mu-s}`` is TWIST-INVARIANT (tensoring by a line bundle preserves
  slope stability and Delta), so the twisted-collection Kronecker values at
  ``nu`` are the untwisted formula at translated slopes ``nu - L``;
* on ``F_0`` the ruling swap is an automorphism carrying ``H_m`` to the
  ``H_{1/m}`` ray, so the swapped-family values are the formula at
  ``(swap(nu), 1/m)``;
* the result's rigor is ``PROVEN`` exactly where the DLP part is certified
  sharp AND dominates (the anticanonical del Pezzo ray); everywhere else it is
  ``CONJECTURAL`` -- exact IF Conjecture A holds and the searched orbit
  suffices, an explicit hypothesis on the certificate.

Differential gates (tests): the Sec. 18 grid points reproduce the E14-M2
values; the anticanonical ray reproduces the sharp envelope with the Kronecker
part never exceeding it (a genuine consistency check: a Kronecker value above
the PROVEN sharp bound would falsify the formula or the machinery); the E14-M1
sandwich brackets the conjectural value wherever both apply.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Optional, Sequence, Tuple

from .chern import Number
from .varieties import Surface
from .dlp_hirzebruch import (
    DEFAULT_RANK_MAX,
    dlp_envelope,
    hirzebruch_index,
    is_del_pezzo_anticanonical,
)
from .delta_sharp import delta_kronecker, surface_with_index
from .rigor import Certificate, Rigor

__all__ = ["DeltaConjectural", "delta_conjectural"]

_CITATIONS = (
    "arXiv:1907.06739 Sec. 1.5 (the conjecture and the promised inductive computation)",
    "thm-deltaKronecker (E14-M2), dlp_envelope (E11-M6), cor-deltaDLP (sharp on -K)",
)


def _Q(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


@dataclass(frozen=True)
class DeltaConjectural:
    """The conditional value with its provenance split.

    ``value = max(dlp_part, kronecker_part)``; ``certificate.rigor`` is PROVEN
    only where the DLP part is certified sharp and dominates, else CONJECTURAL
    (exact under Conjecture A + the searched-orbit hypothesis, recorded).
    ``kronecker_witness`` is ``(translated nu, index m used, value)`` for the
    best Kronecker contribution, ``None`` if none applied.
    """

    value: Fraction
    dlp_part: Fraction
    kronecker_part: Optional[Fraction]
    kronecker_witness: Optional[Tuple[Tuple[Fraction, Fraction], Fraction, Fraction]]
    certificate: Certificate


def delta_conjectural(nu: Sequence[Number], m: Number, surface: Surface,
                      twist_max: int = 4, l_max: int = 10,
                      envelope_rank_max: int = DEFAULT_RANK_MAX) -> DeltaConjectural:
    """Evaluate the Conjecture A-conditional ``delta_m^{mu-s}(nu)`` on the
    ``F_e`` family of ``surface`` (carried polarization ignored)."""
    e = hirzebruch_index(surface)
    m = _Q(m)
    nu_t = tuple(_Q(x) for x in nu)
    if len(nu_t) != 2:
        raise ValueError("nu must be a length-2 F_e NS-vector (f, s)")

    S = surface_with_index(e, m)
    env = dlp_envelope(nu_t, S, rank_max=envelope_rank_max)
    dlp_part = env.value

    best_k: Optional[Fraction] = None
    witness = None
    orbits = [(nu_t, m)]
    if e == 0:
        orbits.append(((nu_t[1], nu_t[0]), 1 / m))          # ruling swap
    for base_nu, base_m in orbits:
        for tf in range(-twist_max, twist_max + 1):
            for ts in range(-twist_max, twist_max + 1):
                shifted = (base_nu[0] - tf, base_nu[1] - ts)
                dk = delta_kronecker(shifted, base_m, surface, l_max=l_max)
                if dk is not None and (best_k is None or dk > best_k):
                    best_k = dk
                    witness = (shifted, base_m, dk)

    value = dlp_part if best_k is None else max(dlp_part, best_k)
    sharp_dlp = is_del_pezzo_anticanonical(S) and env.exact
    proven = sharp_dlp and (best_k is None or best_k <= dlp_part)
    cert = Certificate(
        rigor=Rigor.PROVEN if proven else Rigor.CONJECTURAL,
        hypotheses=(() if proven else (
            "Conjecture A (arXiv:1907.06739 Sec. 1.5): orthogonal Kronecker pairs are "
            "the only non-DLP source of generic HN filtrations",
            f"the searched orbit suffices (Sec. 8 family, l <= {l_max}, "
            f"|twist| <= {twist_max}, ruling swap on F_0)",
        )),
        citations=_CITATIONS,
        note="delta_m^{mu-s} via the conjecture-gated max (E15-M4)")
    return DeltaConjectural(value=value, dlp_part=dlp_part, kronecker_part=best_k,
                            kronecker_witness=witness, certificate=cert)
