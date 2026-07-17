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
from .nslattice import NSLattice
from .rigor import Rigor, Certificate, meet, UNKNOWN_CERTIFICATE
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
    enriques,
    bielliptic,
    abelian_threefold,
    fano_picard_one,
    koseki_product,
    require_faithful_computation,
)
from .exceptional import (
    Bundle,
    P,
    chi,
    is_exceptional,
    is_exceptional_slope,
    exceptional_slopes,
    exceptional_mediant,
    enumerate_exceptional,
    markov_numbers,
    is_markov_number,
)
from .exceptional_surface import (
    SurfaceBundle,
    canonical_class,
    euler_gram,
    p1xp1_collection,
    hirzebruch_collection,
    del_pezzo_collection,
    is_exceptional_collection,
)
from .exceptional_surface import chi as surface_chi
from .dlp import (
    delta,
    dlp_curve,
    DLPCurve,
    control_interval_halfwidth,
    moduli_nonempty,
)
from .nonemptiness_rational import (
    delta_H,
    discriminant_H,
    hirzebruch_with_polarization,
    NonemptinessVerdict,
    SharpBoundEvidence,
    HNMode,
    moduli_nonempty as moduli_nonempty_rational,
)
from .dlp_hirzebruch import (
    DLPEnvelope,
    discriminant,
    dlp_bundle,
    dlp_envelope,
    dlp_restricted,
    emptiness_bound,
    exceptional_discriminant,
    hilbert_P,
    is_semiexceptional,
    is_stable_exceptional,
    total_slope,
)
from .reduction import (
    reduce as reduce_character,   # avoid shadowing functools.reduce at the package top level
    reduce_to_del_pezzo,
    pi_c1,
)
from .prioritary import (
    delta_prioritary,
    prioritary_nonempty,
    generic_prioritary_index,
    delta_prioritary_bundle,
)
from .hn_filtration import (
    HNRegion,
    HNVerdict,
    hn_verdict,
    semistable_exists,
    generic_hn_length,
    hn_region,
    THM_1_13_MIN_DELTA,
    THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS,
)
from .generic_hn import (
    generic_hn_factors,
    semistable_exists_hn,
)
from .delta_sharp import (
    DeltaSharp,
    KroneckerData,
    delta_kronecker,
    delta_mu_stable,
    kronecker_data,
    mu_stable_exists,
    polarization_index,
    surface_with_index,
)
from .stability_interval import (
    StabilityInterval,
    stability_interval,
)
from .walls import (
    Wall,
    VerticalWall,
    ActualWall,
    numerical_wall,
    numerical_wall_ns,
    abelian_wall,
    compute_walls,
    walls_from_subobjects,
    actual_walls,
    actual_walls_complete,
    maciocia_wall_bound,
    enumerate_ns_walls,
)
from .bg_check import (
    check_bg_surface,
    check_bg_threefold,
    BGResult,
    ExistenceResult,
    check_existence_k3,
    check_existence_abelian,
    delta_min,
)
from .threefold import (
    ThreefoldChern,
    bmt_Q,
    bmt_Q_at,
    nu,
    alpha_crit,
    bg_boundary_curve,
    BGBoundary,
    ThreefoldTiltWall,
    ThreefoldVerticalTiltWall,
    numerical_tilt_wall,
    tilt_wall_coeffs,
    BridgelandWall,
    bridgeland_wall,
    bridgeland_wall_coeffs,
    enumerate_tilt_walls,
    BridgelandWallVerdict,
    bridgeland_wall_verdict,
    is_bridgeland_certified,
    TILT_SOLVER_RIGOR,
)
from .mukai import (
    MukaiVector,
    pairing,
    self_pairing,
    moduli_dim,
    is_spherical,
    is_isotropic,
    classify_wall,
    classify_wall_certified,
    WallClassification,
    k3_wall,
    k3_wall_classified,
)

__version__ = "0.1.0"

__all__ = [
    "ChernChar",
    "Q",
    "NSLattice",
    "Rigor",
    "Certificate",
    "meet",
    "UNKNOWN_CERTIFICATE",
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
    "enriques",
    "bielliptic",
    "abelian_threefold",
    "fano_picard_one",
    "koseki_product",
    "require_faithful_computation",
    "Bundle",
    "P",
    "chi",
    "is_exceptional",
    "is_exceptional_slope",
    "exceptional_slopes",
    "exceptional_mediant",
    "enumerate_exceptional",
    "markov_numbers",
    "is_markov_number",
    "SurfaceBundle",
    "canonical_class",
    "euler_gram",
    "p1xp1_collection",
    "hirzebruch_collection",
    "del_pezzo_collection",
    "is_exceptional_collection",
    "surface_chi",
    "delta",
    "dlp_curve",
    "DLPCurve",
    "control_interval_halfwidth",
    "moduli_nonempty",
    "delta_H",
    "discriminant_H",
    "hirzebruch_with_polarization",
    "NonemptinessVerdict",
    "SharpBoundEvidence",
    "HNMode",
    "moduli_nonempty_rational",
    # E11-M6 / G18: the polarization-dependent Drezet-Le Potier envelope on F_e
    "DLPEnvelope",
    "discriminant",
    "dlp_bundle",
    "dlp_envelope",
    "dlp_restricted",
    "emptiness_bound",
    "exceptional_discriminant",
    "hilbert_P",
    "is_semiexceptional",
    "is_stable_exceptional",
    "total_slope",
    # E13-M1 / G18: the F_e -> F_{e-2} reduction map pi (arXiv:1907.06739 Sec.11.1)
    "reduce_character",
    "reduce_to_del_pezzo",
    "pi_c1",
    # E13-M2 / G18: the prioritary sharp bound delta^p_n (arXiv:1907.06739 Prop 4.15 / Cor 4.17)
    "delta_prioritary",
    "prioritary_nonempty",
    "generic_prioritary_index",
    "delta_prioritary_bundle",
    # E13-M3a / G18: HN-length-one existence criterion + Thm 1.13 structure
    "HNRegion",
    "HNVerdict",
    "hn_verdict",
    "generic_hn_factors",
    "semistable_exists_hn",
    "semistable_exists",
    "generic_hn_length",
    "hn_region",
    "THM_1_13_MIN_DELTA",
    "THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS",
    # E14-M1 / G18: the sharp Bogomolov function delta_m^{mu-s} as a computable object
    "DeltaSharp",
    "delta_mu_stable",
    "mu_stable_exists",
    "polarization_index",
    "surface_with_index",
    # E14-M2 / G18: thm-deltaKronecker -- the closed formula on the Kronecker triangle
    "KroneckerData",
    "delta_kronecker",
    "kronecker_data",
    # E14-M3 / G18: stability intervals I_V of exceptional bundles
    "StabilityInterval",
    "stability_interval",
    "Wall",
    "VerticalWall",
    "ActualWall",
    "numerical_wall",
    "numerical_wall_ns",
    "abelian_wall",
    "compute_walls",
    "walls_from_subobjects",
    "actual_walls",
    "actual_walls_complete",
    "maciocia_wall_bound",
    "enumerate_ns_walls",
    "check_bg_surface",
    "check_bg_threefold",
    "BGResult",
    "ExistenceResult",
    "check_existence_k3",
    "check_existence_abelian",
    "delta_min",
    "ThreefoldChern",
    "bmt_Q",
    "bmt_Q_at",
    "nu",
    "alpha_crit",
    "bg_boundary_curve",
    "BGBoundary",
    "ThreefoldTiltWall",
    "ThreefoldVerticalTiltWall",
    "numerical_tilt_wall",
    "tilt_wall_coeffs",
    "BridgelandWall",
    "bridgeland_wall",
    "bridgeland_wall_coeffs",
    "enumerate_tilt_walls",
    "BridgelandWallVerdict",
    "bridgeland_wall_verdict",
    "is_bridgeland_certified",
    "TILT_SOLVER_RIGOR",
    "MukaiVector",
    "pairing",
    "self_pairing",
    "moduli_dim",
    "is_spherical",
    "is_isotropic",
    "classify_wall",
    "classify_wall_certified",
    "WallClassification",
    "k3_wall",
    "k3_wall_classified",
]
