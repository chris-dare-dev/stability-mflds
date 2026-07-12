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
from typing import List, Optional, Tuple

from .nslattice import NSLattice, rank1_shim, RANK1_AMPLE
from .rigor import Rigor, Certificate, UNKNOWN_CERTIFICATE


@dataclass(frozen=True)
class Surface:
    """A polarized smooth projective surface ``(X, H)``.

    Parameters
    ----------
    name : str
    d : int
        ``H^2`` (the degree of the polarization).
    K : Tuple[int, ...]
        The canonical class ``K_X`` as an integer NS-vector in the coordinates of
        ``lattice`` (rank-1 shim: ``(-3,)`` on P^2, ``(0,)`` on K3/abelian).  The
        intersection number ``K_X . H`` is the derived property :attr:`K_H`
        (A8) -- it is NO LONGER a stored scalar, so it stays correct on every
        F_n polarization instead of a ``-2`` placeholder.
    chi_O : int
        ``chi(O_X)``.
    picard_rank : int
    canonical_order : int
        Torsion order of ``K_X`` in ``Pic(X)``.  ``0`` means none / genuinely
        trivial (P^2 has ``K = -3H``; K3 and abelian have ``K = 0``); a positive
        value records torsion (Enriques ``2``; bielliptic ``2``/``3``/``4``/``6``).
        Distinguishes numerically-K-trivial families (Enriques vs. K3/abelian)
        without the full NS-lattice refactor (G12).
    kind : str
        Free-form family tag ("rational", "K3", "abelian", "enriques",
        "bielliptic", ...).
    """

    name: str
    d: int
    K: Tuple[int, ...]                          # canonical class as an NS-vector (A8); K_H is derived
    chi_O: int
    picard_rank: int = 1
    canonical_order: int = 0   # torsion order of K_X in Pic (0 = none/genuinely trivial)
    kind: str = ""
    note: str = ""
    certificate: Certificate = UNKNOWN_CERTIFICATE
    H: Tuple[int, ...] = RANK1_AMPLE           # polarization coordinate vector: ample, or
                                               # nef-and-big on the nef boundary (F_n); shim: (1,)
    ns_lattice: Optional[NSLattice] = None     # explicit NS lattice; None => rank-1 shim from d

    def __post_init__(self) -> None:
        # H must match the NS-lattice rank (rank-1 shim when ns_lattice is unset) and be
        # integral -- guards a mismatched polarization vector now that explicit rank-2
        # lattices exist (E8-M4; keyed on the lattice rank, NOT picard_rank).
        rank = self.ns_lattice.rank if self.ns_lattice is not None else 1
        if len(self.H) != rank:
            raise ValueError(
                f"H has length {len(self.H)} but the NS lattice has rank {rank}")
        if any(int(x) != x for x in self.H):
            raise ValueError(f"H must be an integer vector; got {self.H!r}")
        # K (canonical class as an NS-vector, A8): same rank as H, integral.  Validate
        # integrality FIRST (so a fractional K raises rather than being truncated), then
        # normalize to an int-tuple -- matching the `Tuple[int, ...]` annotation and making
        # equality/hashing independent of whether the caller passed ints or Fractions
        # (canonical_class re-coerces to Fraction via _Q where it needs them).
        # A clean error when K is a scalar -- e.g. a caller using the OLD positional
        # form Surface(name, d, -3, ...) where the 3rd field used to be the int K_H.
        # Without this, `tuple(self.K)` raises a cryptic "'int' object is not iterable".
        if isinstance(self.K, (int, float)) or not hasattr(self.K, "__iter__"):
            raise ValueError(
                f"K must be an NS-vector (a tuple like (-3,) or (0, 0)), not the scalar "
                f"{self.K!r}.  Note: the old scalar `K_H` field was replaced by the NS-vector "
                f"`K` (A8); `K_H` is now a derived property.")
        Kt = tuple(self.K)
        if len(Kt) != rank:
            raise ValueError(
                f"K has length {len(Kt)} but the NS lattice has rank {rank}")
        if any(int(x) != x for x in Kt):
            raise ValueError(f"K must be an integer NS-vector; got {self.K!r}")
        object.__setattr__(self, "K", tuple(int(x) for x in Kt))

    @property
    def is_p2(self) -> bool:
        return self.kind == "P2"

    @property
    def K_H(self) -> int:
        """``K_X . H`` computed from the canonical class and the polarization (A8).

        ``<K, H>`` via ``lattice.pairing`` -- so it is correct on every F_n
        polarization (``(n-2)b - 2a`` for ``H = a f + b s``), where the old stored
        ``K_H = -2`` was a placeholder.  Asserted integral (a canonical
        intersection number) and returned as an exact ``int`` -- no truncation.
        """
        val = self.lattice.pairing(self.K, self.H)
        if val.denominator != 1:
            raise ValueError(f"K.H = {val} is not integral for {self.name}")
        return int(val)

    @property
    def trivial_canonical(self) -> bool:
        """K3 / abelian: ``K_X`` numerically trivial (``K . H = 0``)."""
        return self.K_H == 0

    @property
    def faithful_computation_supported(self) -> bool:
        """Whether the scalar rank-1 (r, ch1.H, ch2) model faithfully REPRESENTS
        this surface's Chern classes.

        Keyed on canonical torsion (``canonical_order == 0``), NOT on Picard rank:
        - True for P^2, K3, abelian (rho=1) AND the rho>=2 rational rows
          (P^1xP^1, F_n): their H-projected computation is valid -- only
          HEURISTIC at rho>=2 (see E2-M5), never refused.
        - False for torsion-canonical surfaces (Enriques order 2; bielliptic
          order 2/3/4/6): the model lacks a torsion-aware NS invariant so the
          classes are unrepresentable until the NS-lattice refactor (G12).
        The E7-M2 surface-consuming-API guard keys its 'requires NS-lattice
        refactor (G12)' error on this flag.
        """
        return self.canonical_order == 0

    @property
    def lattice(self) -> NSLattice:
        """The Neron-Severi lattice with its symmetric intersection form.

        For a Picard-rank-1 row where ``ns_lattice`` is unset, this is the
        rank-1 backward-compat shim ``NSLattice(1, ((d,),))`` whose ample
        generator ``H=(1,)`` satisfies ``<H,H> = d = H^2``.  E8-M4 supplies the
        genuine rank-2 lattices for P^1xP^1 / F_n.
        """
        return self.ns_lattice if self.ns_lattice is not None else rank1_shim(self.d)


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
    certificate: Certificate = UNKNOWN_CERTIFICATE


# --------------------------------------------------------------------------
# Surface catalog
# --------------------------------------------------------------------------
P2 = Surface(
    name="P^2",
    d=1,
    K=(-3,),  # K = -3H (rank-1 shim: K.H = <(-3,),(1,)> = -3)
    chi_O=1,
    picard_rank=1,
    kind="P2",
    note="The Drezet-Le Potier surface; exceptional bundles + DLP curve fully explicit.",
    certificate=Certificate(
        Rigor.PROVEN,
        ("Picard rank 1", "P^2 exceptional-bundle / DLP structure fully explicit"),
        ("Drezet-Le Potier, Ann. Sci. ENS 18 (1985)", "arXiv:1203.0316"),
        "DLP curve and ABCH (s,t)-wall structure are theorems on P^2.",
    ),
)

P1xP1 = Surface(
    name="P^1 x P^1",
    d=2,  # H = O(1,1) = f + s, H^2 = <f+s, f+s> = 2
    K=(-2, -2),  # K = O(-2,-2); K.H = <(-2,-2),(1,1)> = -4 on Gram [[0,1],[1,0]]
    chi_O=1,
    picard_rank=2,
    kind="rational",
    note="Hirzebruch F_0 with the (1,1) polarization; NS = hyperbolic plane U "
         "(basis fiber classes f, s; Gram [[0,1],[1,0]]), H = f + s.",
    H=(1, 1),                                   # H = f + s, <H,H> = 2 = d (E8-M4/G12.4)
    ns_lattice=NSLattice(2, ((0, 1), (1, 0))),  # rank-2 hyperbolic NS lattice
)


def hirzebruch(n: int) -> Surface:
    """Hirzebruch surface F_n with the polarization ``H`` of self-intersection ``n``.

    NS(F_n) is the rank-2 lattice with Gram ``[[0,1],[1,-n]]`` in the basis
    ``(f, s)`` -- ``f`` the fiber class (``f^2 = 0``), ``s`` the negative section
    (``s^2 = -n``, ``f.s = 1``); ``det = -1 < 0`` is the Hodge-index signature
    ``(1,1)`` (E8-M4 / G12.4).  The degree-``n`` polarization ``H = n f + s``
    (coordinate ``(n, 1)``) reproduces ``<H,H> = n = d``.  (``H = C_0 + n f`` is nef
    and big -- it sits on the nef boundary and is semiample, NOT strictly ample, for
    ``n >= 1``; only ``<H,H> = n`` enters the slope/discriminant machinery.)  The canonical
    class is stored as the NS-vector ``K = -(n+2) f - 2 s = (-(n+2), -2)`` (A8), so the
    derived :attr:`Surface.K_H` is the true ``K.H = -(n+2)`` on this polarization, not the
    old ``-2`` placeholder.  For ``n <= 0`` (``F_0 = P^1 x P^1``, or a degenerate input)
    this keeps the rank-1 backward-compat shim unchanged (``K = (-2,)``, preserving
    ``K.H = -2`` there); use the named catalog row :data:`P1xP1` for the genuine rank-2
    ``F_0`` lattice.
    """
    if n >= 1:
        d, H, ns = n, (n, 1), NSLattice(2, ((0, 1), (1, -n)))
        Kv = (-(n + 2), -2)                       # K_{F_n} = -(n+2) f - 2 s; K.H = -(n+2)
        note = ("Hirzebruch surface (Coskun-Huizenga non-emptiness algorithm "
                "applies); NS Gram [[0,1],[1,-n]], H = n f + s.")
    else:  # F_0 = P^1 x P^1 (see P1xP1) or degenerate input: rank-1 shim, unchanged
        d, H, ns = 1, RANK1_AMPLE, None
        Kv = (-2,)                                # rank-1 shim: K.H = -2 (old placeholder)
        note = "Hirzebruch surface (Coskun-Huizenga non-emptiness algorithm applies)."
    return Surface(
        name=f"F_{n}", d=d, K=Kv, chi_O=1, picard_rank=2, kind="hirzebruch",
        note=note, H=H, ns_lattice=ns,
    )


def K3(h2: int = 2) -> Surface:
    """A K3 surface with ``H^2 = h2`` (``h2`` even).  ``K = 0``, ``chi(O) = 2``."""
    if h2 % 2 != 0 or h2 <= 0:
        raise ValueError("a K3 polarization has even positive self-intersection H^2 = 2n")
    return Surface(
        name=f"K3 (H^2={h2})",
        d=h2,
        K=(0,),  # K = 0 (rank-1 shim: K.H = 0)
        chi_O=2,
        picard_rank=1,
        kind="K3",
        note="Algebraic Mukai lattice; see bridgeland_stability.mukai.",
        certificate=Certificate(
            Rigor.PROVEN,
            ("Picard rank 1",
             "K3 surface carries Bridgeland stability conditions"),
            ("arXiv:math/0307164", "arXiv:1301.6968"),
            "K3 (s,t)-wall geometry is a theorem at Picard rank 1 "
            "(Bridgeland; Bayer-Macri).",
        ),
    )


def abelian_surface(h2: int = 2) -> Surface:
    """An abelian surface with ``H^2 = h2``.  ``K = 0``, ``chi(O) = 0``.

    Wall convention (see ``docs/CORRECTIONS.md`` §6).  An abelian surface has
    trivial tangent bundle, so ``sqrt(td) = (1, 0, 0)`` and the Mukai vector IS the
    *bare* Chern triple ``(r, c1, ch2)`` with NO shift.  Route abelian (s, t) walls
    through ``walls.abelian_wall(v, w, surface)`` (which returns exactly
    ``numerical_wall`` on the bare ``(r, c, e)``), NOT through the K3
    ``ch2 -> ch2 + r`` Mukai shim (``mukai.from_chern`` / ``mukai.k3_wall``, G3):
    that shim is K3-only and would inject a spurious ``+2/d`` on an abelian class.
    """
    return Surface(
        name=f"abelian surface (H^2={h2})",
        d=h2,
        K=(0,),  # K = 0 (rank-1 shim: K.H = 0)
        chi_O=0,
        picard_rank=1,
        kind="abelian",
        note="K trivial, chi(O)=0.",
        certificate=Certificate(
            Rigor.PROVEN,
            ("Picard rank 1",
             "abelian surface carries Bridgeland stability conditions"),
            ("arXiv:0708.2247", "arXiv:1203.0884"),
            "abelian (s,t)-wall geometry is a theorem at Picard rank 1 "
            "(Arcara-Bertram-Lieblich; Yanagida-Yoshioka).",
        ),
    )


def enriques(h2: int = 2) -> Surface:
    """An Enriques surface, representative polarization H^2 = h2 (even).

    chi(O)=1 (q=0, p_g=0); K_X is 2-torsion -> numerically trivial (K.H=0) but
    NOT trivial, recorded by canonical_order=2. Picard rank 10 (Num = U + E8(-1),
    even). Record-only row: faithful_computation_supported is False (needs G12).
    """
    if h2 % 2 != 0 or h2 <= 0:
        raise ValueError("an Enriques polarization has even positive H^2 "
                         "(the Enriques lattice U + E8(-1) is even)")
    return Surface(
        name=f"Enriques surface (H^2={h2})",
        d=h2, K=(0,), chi_O=1,  # K numerically trivial (2-torsion): K.H = 0
        picard_rank=10, canonical_order=2, kind="enriques",
        note="K_X 2-torsion (numerically trivial); record-only row -- faithful "
             "walls need the NS-lattice refactor (G12).",
        certificate=Certificate(
            Rigor.PROVEN,
            ("K_X is 2-torsion: numerically trivial (K.H=0) but non-trivial",
             "Enriques surfaces carry Bridgeland stability conditions; "
             "projective moduli / MMP via wall-crossing established",
             "scalar rank-1 model cannot faithfully compute walls "
             "(needs a torsion-aware NS lattice, G12)"),
            ("arXiv:1901.04848", "arXiv:1607.04946"),
            "Enriques moduli established (Nuer-Yoshioka; Yoshioka); recorded as a "
            "catalog row, faithful computation deferred to G12.",
        ),
    )


def bielliptic(h2: int = 2) -> Surface:
    """A bielliptic (hyperelliptic) surface -- order-2 Bagnera-de Franchis type.

    chi(O)=0 (q=1, p_g=0); K_X torsion, numerically trivial (K.H=0). This factory
    returns the order-2 representative (canonical_order=2, Picard rank 2); across
    the seven Bagnera-de Franchis types K_X has order 2, 3, 4, or 6, but the Picard
    rank is always 2 (p_g=0 => NS(X)_Q = H^2, so rho = b_2 = 2 for every type).
    Record-only row: faithful_computation_supported is False (needs G12).

    H^2 must be even and positive: chi(O)=0 with numerically-trivial K gives, by
    Riemann-Roch, chi(O(D)) = D^2/2 in Z, so D^2 (hence H^2) is even. H^2 = 1
    cannot occur (matching the Enriques even-lattice constraint).
    """
    if h2 % 2 != 0 or h2 <= 0:
        raise ValueError("a bielliptic polarization has even positive H^2 "
                         "(chi(O)=0 with numerically-trivial K forces "
                         "chi(O(D)) = D^2/2 in Z by Riemann-Roch)")
    return Surface(
        name=f"bielliptic surface (H^2={h2})",
        d=h2, K=(0,), chi_O=0,  # K torsion, numerically trivial: K.H = 0
        picard_rank=2, canonical_order=2, kind="bielliptic",
        note="K_X torsion (order 2 here; orders 2,3,4,6 occur), numerically "
             "trivial; record-only row -- faithful walls need the NS-lattice "
             "refactor (G12).",
        certificate=Certificate(
            Rigor.PROVEN,
            ("K_X is torsion (order 2 for this Bagnera-de Franchis type; orders "
             "2,3,4,6 occur across the seven types), numerically trivial",
             "Chern-character classification of stable sheaves and Bridgeland "
             "wall-crossing on bielliptic surfaces is complete",
             "Picard rank 2 (p_g=0 => rho = b_2 = 2 for every type) with torsion "
             "K_X; scalar rank-1 model cannot faithfully compute walls (needs the "
             "NS lattice, G12)"),
            ("arXiv:2107.13370",),
            "bielliptic stable-sheaf classification / wall-crossing complete "
            "(Nuer, 'Stable sheaves on bielliptic surfaces', arXiv:2107.13370); "
            "recorded as a catalog row, faithful computation deferred to G12.",
        ),
    )


# --------------------------------------------------------------------------
# Surface-consuming-API guard (E7-M2 / G11).  The scalar rank-1 (r, ch1.H, ch2)
# model cannot faithfully represent a torsion-canonical surface's Chern classes
# (Enriques / bielliptic: canonical_order != 0), so every Surface-consuming
# entry point (compute_walls, actual_walls, moduli_nonempty, abelian_wall,
# k3_wall_classified, check_existence_*) routes through this guard BEFORE
# computing.  It lives here in varieties.py -- imported by every consumer, and
# imports only .rigor, so there is no import cycle.  numerical_wall(v, w, d) is
# deliberately NOT guarded: its signature takes only a ChernChar and an int, so
# it never sees a Surface (the caveat there is docs-only).
# --------------------------------------------------------------------------
def require_faithful_computation(surface: Surface) -> None:
    """Raise ``NotImplementedError`` unless ``surface`` can be faithfully computed.

    A no-op for every ``faithful_computation_supported`` surface (P^2, K3,
    abelian, and the rho>=2 rational rows).  For a torsion-canonical surface
    (Enriques order 2; bielliptic order 2/3/4/6, where
    ``faithful_computation_supported is False``) it raises with a message
    containing the literal substring ``NS-lattice refactor (G12)``: the scalar
    rank-1 ``(r, ch1.H, ch2)`` model lacks a torsion-aware NS invariant, so the
    classes are unrepresentable until the NS-lattice refactor (G12).  The
    exception type is ``NotImplementedError`` (semantically "needs G12").
    """
    if not surface.faithful_computation_supported:
        raise NotImplementedError(
            f"{surface.name}: faithful Bridgeland-wall / moduli computation on a "
            f"torsion-canonical surface (canonical_order={surface.canonical_order}) "
            "requires the NS-lattice refactor (G12). The scalar rank-1 "
            "(r, ch1.H, ch2) model cannot represent its Chern classes; this variety "
            "is recorded as a catalog row only."
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
    certificate=Certificate(
        Rigor.PROVEN,
        ("threefold tilt-BG (BMT) inequality is a theorem here",),
        ("arXiv:1103.5010", "arXiv:1207.4980"),
        "tilt-BG boundary rigorous (bg_proven).",
    ),
)

QUADRIC3 = Threefold(
    name="quadric 3-fold Q in P^4",
    d3=2,
    chi_O=1,
    bg_proven=True,
    kind="fano",
    note=(
        "BG proven as a Fano 3-fold of Picard rank 1 (C. Li, 1510.04089); "
        "existence of Bridgeland stability conditions on all Fano 3-folds "
        "(Bernardara-Macri-Schmidt-Zhao [BMSZ], 1607.08199). Schmidt (1509.04608) "
        "develops threefold Bridgeland wall-crossing techniques with worked "
        "examples on P^3 (twisted cubics, complete intersections) -- it is NOT "
        "quadric-specific; the quadric's own generalized BG inequality (Schmidt, "
        "1309.4265) is subsumed by the Fano Picard-rank-1 theorem 1510.04089."
    ),
    references=["arXiv:1510.04089", "arXiv:1607.08199", "arXiv:1509.04608"],
    certificate=Certificate(
        Rigor.PROVEN,
        ("threefold tilt-BG (BMT) inequality is a theorem here",),
        ("arXiv:1510.04089", "arXiv:1607.08199", "arXiv:1509.04608"),
        "tilt-BG boundary rigorous (bg_proven).",
    ),
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
        certificate=Certificate(
            Rigor.PROVEN,
            ("threefold tilt-BG (BMT) inequality is a theorem here",),
            ("arXiv:1510.04089",),
            "tilt-BG boundary rigorous (bg_proven).",
        ),
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
        certificate=Certificate(
            Rigor.PROVEN,
            ("threefold tilt-BG (BMT) inequality is a theorem here",),
            ("arXiv:1410.1585",),
            "tilt-BG boundary rigorous (bg_proven).",
        ),
    )


QUINTIC = Threefold(
    name="quintic CY3",
    d3=5,
    chi_O=0,
    bg_proven=True,
    kind="calabi-yau",
    note="BG proven: C. Li (1810.03434); first geometric stability on a strict CY3.",
    references=["arXiv:1810.03434"],
    certificate=Certificate(
        Rigor.PROVEN,
        ("threefold tilt-BG (BMT) inequality is a theorem here",),
        ("arXiv:1810.03434",),
        "tilt-BG boundary rigorous (bg_proven).",
    ),
)

BLOWUP_P3_POINT = Threefold(
    name="Bl_p(P^3)",
    d3=1,
    chi_O=1,
    bg_proven=False,
    kind="fano-blowup",
    note=(
        "Schmidt's counterexample (1602.05055): the STRONGER/generalized "
        "(strong BMT) tilt-BG inequality FAILS here, so bg_proven=False and "
        "Algorithm 5's boundary is not rigorous for this variety. Bl_p(P^3) is "
        "nonetheless Fano and DOES carry Bridgeland stability conditions "
        "(BMSZ, 1607.08199); Piyaratne (1705.04011) established an optimal "
        "modified BG-type inequality. Only the strong BMT boundary fails; "
        "stability conditions still exist."
    ),
    references=["arXiv:1602.05055", "arXiv:1705.04011"],
    certificate=Certificate(
        Rigor.CONJECTURAL,
        ("threefold strong-BMT / tilt-BG inequality is NOT a theorem here",),
        ("arXiv:1602.05055", "arXiv:1607.08199", "arXiv:1705.04011"),
        "strong BMT boundary not rigorous (Schmidt 1602.05055); stability "
        "conditions nonetheless exist (Fano; BMSZ 1607.08199; Piyaratne 1705.04011).",
    ),
)


# --------------------------------------------------------------------------
# Koseki product threefolds (E7-M2 / G11) [RESEARCH/UNVERIFIED].
#
# Koseki proves the tilt-BG inequality (hence Bridgeland stability conditions)
# for the product threefolds P^1 x S (S abelian), P^2 x C and P^1 x P^1 x C
# (C elliptic): the product-specific primary source is "Stability conditions on
# product threefolds of projective spaces and Abelian varieties"
# (arXiv:1703.07042), generalized to all threefolds with nef tangent bundles in
# arXiv:1811.03267.  Koseki's EXACT hypotheses on the factors are NOT yet checked
# against the papers, so until KOSEKI_BG_VERIFIED is flipped to True after that
# literature sub-task, every koseki_product row stays bg_proven=False and its
# Certificate carries rigor CONJECTURAL (never PROVEN).
#
# This is DISTINCT from Koseki's separate weighted-hypersurface paper on
# Calabi-Yau double/triple solids (arXiv:2007.00044) -- do NOT conflate the
# two; koseki_product cites the product-specific arXiv:1703.07042 and the
# nef-tangent generalization arXiv:1811.03267.  (Note: arXiv:1510.04474 is
# NOT Koseki -- it is Chuang-Lai, "Bogomolov-Gieseker Type Inequality on
# Calabi-Yau and Fano 3-folds", WITHDRAWN -- and must never be cited here.)
# --------------------------------------------------------------------------
KOSEKI_BG_VERIFIED = False

_KOSEKI_PRODUCTS = {
    "P1xS":    "P^1 x S (S an abelian surface)",
    "P2xC":    "P^2 x C (C an elliptic curve)",
    "P1xP1xC": "P^1 x P^1 x C (C an elliptic curve)",
}


def koseki_product(product: str, d3: int) -> Threefold:
    """A Koseki product threefold with nef tangent bundle (arXiv:1811.03267).

    Mirrors :func:`abelian_threefold` / :func:`fano_picard_one`: the polarization
    degree ``d3 = H^3`` is caller-supplied (it depends on the product
    polarization).  ``product`` must be one of ``"P1xS"`` (``P^1 x S``, ``S`` an
    abelian surface), ``"P2xC"`` (``P^2 x C``, ``C`` elliptic) or ``"P1xP1xC"``
    (``P^1 x P^1 x C``, ``C`` elliptic); raises ``ValueError`` otherwise, or if
    ``d3 <= 0``.

    ``chi(O) = 0`` for all three products (re-derived, Kunneth): ``chi(O_{XxY})
    = chi(O_X) * chi(O_Y)``; ``chi(O_{P^1}) = chi(O_{P^2}) = chi(O_{P^1xP^1}) =
    1``, while an abelian surface (``h^{0,1,2} = 1,2,1 -> chi = 0``) or an
    elliptic curve (``chi = 1 - g = 0``) contributes ``0``.

    [RESEARCH/UNVERIFIED].  ``bg_proven`` and the Certificate rigor are gated on
    the module flag ``KOSEKI_BG_VERIFIED`` (currently ``False``): until Koseki's
    exact hypotheses are checked against arXiv:1811.03267, ``bg_proven`` is
    ``False`` and the rigor is ``CONJECTURAL`` (never ``PROVEN``).  Cites only
    arXiv:1811.03267 -- NOT Koseki's separate Calabi-Yau double/triple-solid
    (weighted-hypersurface) paper arXiv:2007.00044.

    (Guidance-only, non-asserted product-polarization degrees, from standard
    intersection theory: ``P^2 x C`` with ``H = h + f`` gives ``H^3 = 3``;
    ``P^1 x S`` with ``H = f + H_S`` gives ``H^3 = 3 * H_S^2``; ``P^1 x P^1 x C``
    with ``H = f1 + f2 + g`` gives ``H^3 = 6``.)
    """
    if product not in _KOSEKI_PRODUCTS:
        raise ValueError(
            f"unknown Koseki product {product!r}; expected one of "
            f"{tuple(_KOSEKI_PRODUCTS)}."
        )
    if d3 <= 0:
        raise ValueError("a Koseki product polarization has positive degree H^3 > 0")
    desc = _KOSEKI_PRODUCTS[product]
    rigor = Rigor.PROVEN if KOSEKI_BG_VERIFIED else Rigor.CONJECTURAL
    note = (
        "[RESEARCH/UNVERIFIED] tilt-BG / Bridgeland stability for this product "
        "threefold (Koseki, product-specific 'Stability conditions on product "
        "threefolds of projective spaces and Abelian varieties', arXiv:1703.07042; "
        "generalized to nef tangent bundles in arXiv:1811.03267). bg_proven is "
        "gated on KOSEKI_BG_VERIFIED (currently False): Koseki's exact hypotheses "
        "on the factors are not yet verified against the papers, so this row stays "
        "CONJECTURAL, never PROVEN. Distinct from Koseki's separate Calabi-Yau "
        "double/triple-solid (weighted-hypersurface) paper (do not conflate; the "
        "specific ID is in the module note). chi(O)=0 by Kunneth."
    )
    return Threefold(
        name=f"{desc} (H^3={d3})",
        d3=d3,
        chi_O=0,  # re-derived via Kunneth (see docstring)
        bg_proven=KOSEKI_BG_VERIFIED,
        kind="koseki-product",
        note=note,
        references=["arXiv:1703.07042", "arXiv:1811.03267"],
        certificate=Certificate(
            rigor,
            (f"{desc}: product threefold with nef tangent bundle",
             "threefold tilt-BG (BMT) inequality for these products proved by "
             "Koseki -- product-specific arXiv:1703.07042, generalized to nef "
             "tangent bundles in arXiv:1811.03267",
             "[RESEARCH/UNVERIFIED] Koseki's exact hypotheses are NOT yet checked "
             "against the papers; bg_proven gated on KOSEKI_BG_VERIFIED "
             "(currently False)"),
            ("arXiv:1703.07042", "arXiv:1811.03267"),
            note,
        ),
    )


# --------------------------------------------------------------------------
# Provenance guard (G1): curated allowlist of vetted algebraic-geometry arXiv
# IDs, plus the tuple of named catalog threefolds.
#
# The allowlist is deliberately maintenance-brittle: it MUST be appended
# whenever a legitimately new AG reference is added to the catalog.  A format
# regex such as ^arXiv:\d{4}\.\d{4,5}$ is NOT a sufficient guard -- the
# probability paper arXiv:1607.07182 ("Interlacing Diffusions", math.PR) that
# QUADRIC3 formerly miscited matches that regex perfectly.  Membership in this
# curated set, not the regex, is the mechanism that catches a stray non-AG ID.
# --------------------------------------------------------------------------
AG_ALLOWLIST = frozenset(
    {
        "arXiv:1103.5010",   # Bayer-Macri-Toda: tilt-stability / BMT inequality
        "arXiv:1207.4980",   # Macri: generalized BG inequality on P^3
        "arXiv:1410.1585",   # Bayer-Macri-Stellari: abelian 3-folds
        "arXiv:1509.04608",  # Schmidt: threefold wall crossings
        "arXiv:1510.04089",  # C. Li: Fano 3-folds of Picard number one
        "arXiv:1602.05055",  # Schmidt: counterexample to the strong BMT bound
        "arXiv:1607.08199",  # Bernardara-Macri-Schmidt-Zhao (BMSZ): Fano 3-folds
        "arXiv:1705.04011",  # Piyaratne: modified BG-type inequality / Fano 3-folds
        "arXiv:1810.03434",  # C. Li: quintic CY3
        "arXiv:1703.07042",  # Koseki: product threefolds (projective spaces x abelian varieties)
        "arXiv:1811.03267",  # Koseki: threefolds with nef tangent bundles (generalization)
    }
)

# Named threefold catalog rows.  Factory-produced rows (fano_picard_one,
# abelian_threefold) are covered by including their arXiv IDs in AG_ALLOWLIST
# and by instantiating them in the allowlist-membership test.
ALL_THREEFOLDS = (P3, QUADRIC3, QUINTIC, BLOWUP_P3_POINT)
