"""Surface and threefold data classes, plus a small catalog.

Only the *numerical* invariants needed by the stability algorithms are stored:
``d = H^2`` (surface) or ``d3 = H^3`` (threefold), the canonical-class pairing
``K . H``, and ``chi(O_X)``.  Exceptional-bundle / DLP machinery is P^2-only
(``Surface.is_p2``); walls and BG work for any surface; the threefold BG
boundary requires a (proven or conjectural) BG inequality flagged by
``Threefold.bg_proven``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List


@dataclass(frozen=True)
class Surface:
    """A polarized smooth projective surface ``(X, H)``.

    Parameters
    ----------
    name : str
    d : int
        ``H^2`` (the degree of the polarization).
    K_H : int
        ``K_X . H``.
    chi_O : int
        ``chi(O_X)``.
    picard_rank : int
    kind : str
        Free-form family tag ("rational", "K3", "abelian", "enriques", ...).
    """

    name: str
    d: int
    K_H: int
    chi_O: int
    picard_rank: int = 1
    kind: str = ""
    note: str = ""

    @property
    def is_p2(self) -> bool:
        return self.kind == "P2"

    @property
    def trivial_canonical(self) -> bool:
        """K3 / abelian: ``K_X`` numerically trivial (``K . H = 0``)."""
        return self.K_H == 0


@dataclass(frozen=True)
class Threefold:
    """A polarized smooth projective threefold ``(X, H)``.

    ``d3 = H^3``.  ``bg_proven`` records whether the Bayer-Macri-Toda
    Bogomolov-Gieseker-type inequality for tilt-stable objects is a *theorem*
    for this variety (vs. conjectural), which determines whether Algorithm 5's
    output is rigorous.
    """

    name: str
    d3: int
    chi_O: int = 1
    bg_proven: bool = False
    kind: str = ""
    note: str = ""
    references: List[str] = field(default_factory=list)


# --------------------------------------------------------------------------
# Surface catalog
# --------------------------------------------------------------------------
P2 = Surface(
    name="P^2",
    d=1,
    K_H=-3,  # K = -3H, K.H = -3
    chi_O=1,
    picard_rank=1,
    kind="P2",
    note="The Drezet-Le Potier surface; exceptional bundles + DLP curve fully explicit.",
)

P1xP1 = Surface(
    name="P^1 x P^1",
    d=2,  # H = O(1,1), H^2 = 2
    K_H=-4,  # K = O(-2,-2), K.H = -4
    chi_O=1,
    picard_rank=2,
    kind="rational",
    note="Hirzebruch F_0 with the (1,1) polarization.",
)


def hirzebruch(n: int) -> Surface:
    """Hirzebruch surface F_n with the polarization ``H`` of self-intersection ``n``.

    Records only ``H^2 = n``; the precise ``K.H`` depends on the chosen H, so
    the slope/discriminant machinery (which needs only ``d = H^2``) is exact,
    while ``K_H`` here is a placeholder for the standard fiber-class normalization.
    """
    return Surface(
        name=f"F_{n}",
        d=n if n > 0 else 1,
        K_H=-2,
        chi_O=1,
        picard_rank=2,
        kind="hirzebruch",
        note="Hirzebruch surface (Coskun-Huizenga non-emptiness algorithm applies).",
    )


def K3(h2: int = 2) -> Surface:
    """A K3 surface with ``H^2 = h2`` (``h2`` even).  ``K = 0``, ``chi(O) = 2``."""
    if h2 % 2 != 0 or h2 <= 0:
        raise ValueError("a K3 polarization has even positive self-intersection H^2 = 2n")
    return Surface(
        name=f"K3 (H^2={h2})",
        d=h2,
        K_H=0,
        chi_O=2,
        picard_rank=1,
        kind="K3",
        note="Algebraic Mukai lattice; see bridgeland_stability.mukai.",
    )


def abelian_surface(h2: int = 2) -> Surface:
    """An abelian surface with ``H^2 = h2``.  ``K = 0``, ``chi(O) = 0``."""
    return Surface(
        name=f"abelian surface (H^2={h2})",
        d=h2,
        K_H=0,
        chi_O=0,
        picard_rank=1,
        kind="abelian",
        note="K trivial, chi(O)=0.",
    )


# --------------------------------------------------------------------------
# Threefold catalog (bg_proven reflects the literature, see docs/CORRECTIONS.md)
# --------------------------------------------------------------------------
P3 = Threefold(
    name="P^3",
    d3=1,
    chi_O=1,
    bg_proven=True,
    kind="fano",
    note="BG proven: Bayer-Macri-Toda (1103.5010), Macri (1207.4980).",
    references=["arXiv:1103.5010", "arXiv:1207.4980"],
)

QUADRIC3 = Threefold(
    name="quadric 3-fold Q in P^4",
    d3=2,
    chi_O=1,
    bg_proven=True,
    kind="fano",
    note="BG proven (Schmidt; Bernardara-Macri-Schmidt-Zhao).",
    references=["arXiv:1607.07182", "arXiv:1705.04011"],
)


def fano_picard_one(name: str, d3: int) -> Threefold:
    """A Fano threefold of Picard rank one (BG proven, Li arXiv:1510.04089)."""
    return Threefold(
        name=name,
        d3=d3,
        chi_O=1,
        bg_proven=True,
        kind="fano",
        note="BG proven for all Fano 3-folds of Picard rank 1 (C. Li, 1510.04089).",
        references=["arXiv:1510.04089"],
    )


def abelian_threefold(d3: int) -> Threefold:
    """An abelian threefold with ``H^3 = d3`` (BG proven, BMS arXiv:1410.1585)."""
    return Threefold(
        name=f"abelian 3-fold (H^3={d3})",
        d3=d3,
        chi_O=0,
        bg_proven=True,
        kind="abelian",
        note="BG proven: Bayer-Macri-Stellari (1410.1585).",
        references=["arXiv:1410.1585"],
    )


QUINTIC = Threefold(
    name="quintic CY3",
    d3=5,
    chi_O=0,
    bg_proven=True,
    kind="calabi-yau",
    note="BG proven: C. Li (1810.03434); first geometric stability on a strict CY3.",
    references=["arXiv:1810.03434"],
)

BLOWUP_P3_POINT = Threefold(
    name="Bl_p(P^3)",
    d3=1,
    chi_O=1,
    bg_proven=False,
    kind="fano-blowup",
    note=(
        "Schmidt's counterexample (1602.05055): the STRONGER/generalized BMT "
        "inequality FAILS here. Algorithm 5 must warn that its boundary is not "
        "rigorous for this variety."
    ),
    references=["arXiv:1602.05055", "arXiv:1705.04011"],
)
