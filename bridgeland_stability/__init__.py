"""bridgeland_stability - DLP curves, BG boundaries, and Bridgeland walls.

A standalone, exact-arithmetic toolkit for three tightly-coupled objects on
projective surfaces and threefolds:

  * the Drezet-Le Potier curve delta(mu) on P^2 and the exceptional bundles
    that shape it (``exceptional``, ``dlp``);
  * the Bogomolov-Gieseker positivity boundary (``bg_check``, ``threefold``);
  * Bridgeland stability walls in the (s, t) upper half-plane (``walls``);
  * the K3 Mukai lattice and Bayer-Macri wall classification (``mukai``).

All core computations use ``fractions.Fraction`` (exact); floats appear only
for geometry/plotting.  See ``docs/CORRECTIONS.md`` for the (substantial) list
of corrections relative to the original project brief, with citations.

Convention: the discriminant is the Coskun-Huizenga normalization
``Delta = (1/2) mu^2 - ch2/(r d)`` (the brief's was twice this).
"""

from __future__ import annotations

from .chern import ChernChar, Q
from .varieties import (
    Surface,
    Threefold,
    P2,
    P1xP1,
    P3,
    QUADRIC3,
    QUINTIC,
    BLOWUP_P3_POINT,
    hirzebruch,
    K3,
    abelian_surface,
    abelian_threefold,
    fano_picard_one,
)
from .exceptional import (
    Bundle,
    P,
    chi,
    is_exceptional,
    exceptional_mediant,
    enumerate_exceptional,
    markov_numbers,
    is_markov_number,
)
from .dlp import (
    delta,
    dlp_curve,
    DLPCurve,
    control_interval_halfwidth,
    moduli_nonempty,
)
from .walls import (
    Wall,
    VerticalWall,
    ActualWall,
    numerical_wall,
    compute_walls,
    walls_from_subobjects,
    actual_walls,
    actual_walls_complete,
)
from .bg_check import check_bg_surface, check_bg_threefold, BGResult
from .threefold import (
    ThreefoldChern,
    bmt_Q,
    bmt_Q_at,
    alpha_crit,
    bg_boundary_curve,
    BGBoundary,
)
from .mukai import (
    MukaiVector,
    pairing,
    self_pairing,
    moduli_dim,
    is_spherical,
    is_isotropic,
    classify_wall,
    WallClassification,
)

__version__ = "0.1.0"

__all__ = [
    "ChernChar",
    "Q",
    "Surface",
    "Threefold",
    "P2",
    "P1xP1",
    "P3",
    "QUADRIC3",
    "QUINTIC",
    "BLOWUP_P3_POINT",
    "hirzebruch",
    "K3",
    "abelian_surface",
    "abelian_threefold",
    "fano_picard_one",
    "Bundle",
    "P",
    "chi",
    "is_exceptional",
    "exceptional_mediant",
    "enumerate_exceptional",
    "markov_numbers",
    "is_markov_number",
    "delta",
    "dlp_curve",
    "DLPCurve",
    "control_interval_halfwidth",
    "moduli_nonempty",
    "Wall",
    "VerticalWall",
    "ActualWall",
    "numerical_wall",
    "compute_walls",
    "walls_from_subobjects",
    "actual_walls",
    "actual_walls_complete",
    "check_bg_surface",
    "check_bg_threefold",
    "BGResult",
    "ThreefoldChern",
    "bmt_Q",
    "bmt_Q_at",
    "alpha_crit",
    "bg_boundary_curve",
    "BGBoundary",
    "MukaiVector",
    "pairing",
    "self_pairing",
    "moduli_dim",
    "is_spherical",
    "is_isotropic",
    "classify_wall",
    "WallClassification",
]
