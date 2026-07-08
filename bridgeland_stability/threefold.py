"""Algorithm 5 - threefold tilt-stability Bogomolov-Gieseker boundary.

A Chern character on a polarized threefold ``(X, H)`` with ``d3 = H^3`` is
recorded by its H-degrees ``(r, a1, a2, a3) = (ch0, ch1.H^2, ch2.H, ch3)``.
The ``beta``-twist ``ch^beta = ch . e^{-beta H}`` gives

    ch1b = a1 - beta r d3
    ch2b = a2 - beta a1 + (beta^2/2) r d3
    ch3b = a3 - beta a2 + (beta^2/2) a1 - (beta^3/6) r d3

The Bayer-Macri-Toda quadratic form and tilt-BG inequality are

    Q = 4 ch2b^2 - 6 ch1b ch3b ,        Q >= alpha^2 (ch1b)^2   (for nu-semistable E),

so the critical-radius boundary curve in the (beta, alpha) upper half-plane is

    alpha_crit(beta) = sqrt(max(0, Q)) / |ch1b|     (undefined where ch1b = 0).

The inequality ``Q >= 0`` is a THEOREM for: P^3 and all Fano 3-folds of Picard
rank 1, abelian 3-folds, and the quintic; it FAILS in its stronger form on
Bl_p(P^3) (Schmidt).  ``bg_boundary_curve`` warns when ``threefold.bg_proven``
is False.

NOTE on the brief's test values: the brief's "alpha_crit(1/2)=sqrt(29)/4" and
"beta=1 => Q=2, alpha_crit=sqrt(2)/2" are arithmetic errors.  Correct, for the
P^3 null-correlation bundle (2,0,1,0): beta=1/2 -> Q=3, alpha_crit=sqrt(3);
beta=1 -> Q=0, alpha_crit=0.  The brief's beta=1 slip drops the rank factor in
the cubic term (ch3b=-7/6 instead of -4/3).  See docs/CORRECTIONS.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import List, Optional, Tuple

from .chern import ChernChar, Number, Q
from .rigor import Rigor, Certificate, UNKNOWN_CERTIFICATE, meet
from .varieties import Threefold
from .walls import numerical_wall, VerticalWall, _half_integers


@dataclass(frozen=True)
class ThreefoldChern:
    """H-degrees ``(r, a1, a2, a3) = (ch0, ch1.H^2, ch2.H, ch3)`` on a threefold."""

    r: int
    a1: Fraction  # ch1 . H^2
    a2: Fraction  # ch2 . H
    a3: Fraction  # ch3

    def __post_init__(self) -> None:
        object.__setattr__(self, "a1", Q(self.a1))
        object.__setattr__(self, "a2", Q(self.a2))
        object.__setattr__(self, "a3", Q(self.a3))

    def twist(self, beta: Number, d3: int) -> "ThreefoldChern":
        b = Q(beta)
        a1 = self.a1 - b * self.r * d3
        a2 = self.a2 - b * self.a1 + (b * b / 2) * self.r * d3
        a3 = self.a3 - b * self.a2 + (b * b / 2) * self.a1 - (b * b * b / 6) * self.r * d3
        return ThreefoldChern(self.r, a1, a2, a3)


def bmt_Q(ch: ThreefoldChern) -> Fraction:
    """The Bayer-Macri-Toda form ``Q = 4 ch2^2 - 6 ch1 ch3`` of a (twisted) character."""
    return 4 * ch.a2 * ch.a2 - 6 * ch.a1 * ch.a3


def bmt_Q_at(v: ThreefoldChern, beta: Number, d3: int) -> Fraction:
    """``Q`` of the ``beta``-twist of ``v`` (exact when ``beta`` is rational)."""
    return bmt_Q(v.twist(beta, d3))


def nu(v: ThreefoldChern, alpha_sq: Number, beta: Number, d3: int) -> Optional[Fraction]:
    r"""BMT tilt-slope ``nu_{alpha,beta}`` of ``v`` (exact ``Fraction``, or ``None`` = ``+inf``).

    The Bayer-Macri-Toda tilt slope on a polarized ``n``-fold ``(X, H)`` is

        nu_{alpha,beta}(E) = (H^{n-2}.ch2^b - 1/2 alpha^2 H^n.ch0^b)
                             / (H^{n-1}.ch1^b),

    where ``ch^b = ch . e^{-beta H}`` is the ``beta``-twist and ``alpha > 0`` is
    the tilt parameter.  In the ``ThreefoldChern(r, a1, a2, a3)`` H-degree
    encoding with ``d3 = H^3`` the numerator is ``a2^b - 1/2 alpha^2 * r * d3``
    and the denominator is ``a1^b``, where ``a1^b, a2^b`` are the ``twist(beta)``
    components.  It is ``+inf`` when the denominator ``H^{n-1}.ch1^b`` vanishes.

    Parameters
    ----------
    alpha_sq:
        This is ``alpha^2`` (NOT ``alpha``), taken as an exact ``Fraction`` so
        that the tilt radius -- which is naturally ``sqrt``-valued -- can be fed
        in without ever introducing a float.  ``Q(alpha_sq) == Fraction(alpha_sq)``.

    Returns
    -------
    A pure ``fractions.Fraction`` (no float anywhere in this function), or
    ``None`` to encode ``nu = +inf`` in the degenerate vertical case
    ``tw.a1 == 0`` (i.e. ``H^{n-1}.ch1^b = 0``).

    Notes
    -----
    ``nu`` depends only on ``(ch0, ch1, ch2)`` and NEVER on ``ch3 = a3``: the
    numerator and denominator use only ``twist``'s ``a2^b`` and ``a1^b``, and
    those are functions of ``(r, a1, a2)`` alone.  (This is the ``ch3``-freedom
    that ``test_nu_ch3_independent`` pins.)

    The pinned P^3 null-correlation values -- for ``v = ThreefoldChern(2,0,1,0)``,
    ``d3 = 1``: ``nu(v,0,1,1) = -1``, ``nu(v,1,1,1) = -1/2``, ``nu(v,2,1,1) = 0``,
    ``nu(v,0,0,1) = None`` -- are EXACTLY DERIVED from the formula above
    (self-consistent exact ``Fraction`` re-computations), NOT read off a source.
    They are consistent with the pinned twist ``v.twist(1,1).a3 = -4/3`` (see
    docs/CORRECTIONS.md sec. 5).  A genuine Schmidt cross-check would first
    require confirming that Schmidt's H-power normalization (num
    ``= H^{n-2}.ch2^b - 1/2 alpha^2 H^n.ch0^b``, den ``= H^{n-1}.ch1^b``) matches
    the convention used here; until that normalization is confirmed these are
    self-consistent derivations, not yet a Schmidt cross-check.

    References
    ----------
    Bayer-Macri-Toda, arXiv:1103.5010 (tilt slope / BMT inequality);
    Schmidt, arXiv:1509.04608 (threefold tilt-stability wall computations).
    """
    tw = v.twist(beta, d3)
    if tw.a1 == 0:
        return None  # H^{n-1}.ch1^b = 0 -> vertical/degenerate (nu = +infinity)
    return (tw.a2 - Q(alpha_sq) / 2 * v.r * d3) / tw.a1


def check_bg_threefold(
    v: ThreefoldChern, alpha: Number, beta: Number, d3: int
) -> dict:
    """Tilt-BG check ``Q >= alpha^2 (ch1b)^2`` at a given ``(alpha, beta)``.

    Meaningful only for objects already known/assumed nu_{alpha,beta}-semistable
    (it reports what BG *requires*, not whether an arbitrary sheaf satisfies it).
    """
    tw = v.twist(beta, d3)
    Qv = bmt_Q(tw)
    threshold = Q(alpha) * Q(alpha) * tw.a1 * tw.a1
    return {
        "ch_twisted": tw,
        "Q": Qv,
        "threshold": threshold,
        "satisfies": Qv >= threshold,
        "slack": Qv - threshold,
    }


def alpha_crit(v: ThreefoldChern, beta: Number, d3: int) -> Optional[float]:
    """Critical radius ``sqrt(max(0,Q))/|ch1b|`` at ``beta``; ``None`` if ``ch1b = 0``."""
    tw = v.twist(beta, d3)
    if tw.a1 == 0:
        return None  # degenerate vertical wall ch1b = 0
    Qv = bmt_Q(tw)
    return math.sqrt(max(0.0, float(Qv))) / abs(float(tw.a1))


@dataclass
class BGBoundary:
    """Sampled BG boundary ``alpha_crit(beta)`` for a threefold."""

    betas: List[float]
    alphas: List[float]
    bg_proven: bool
    threefold_name: str
    degenerate_betas: List[float]
    certificate: Certificate = UNKNOWN_CERTIFICATE


def bg_boundary_curve(
    v: ThreefoldChern,
    threefold: Threefold,
    beta_range: Tuple[float, float] = (-3.0, 3.0),
    N: int = 600,
) -> BGBoundary:
    """Sample ``alpha_crit(beta)`` over ``beta_range`` (float grid)."""
    b0, b1 = beta_range
    betas: List[float] = []
    alphas: List[float] = []
    degenerate: List[float] = []
    for i in range(N + 1):
        beta = b0 + (b1 - b0) * i / N
        # exact ch1b sign check to detect degeneracy robustly
        a1b = float(v.a1) - beta * float(v.r) * float(threefold.d3)
        if abs(a1b) < 1e-12:
            degenerate.append(beta)
            continue
        # float twist for the grid
        rd = float(v.r) * float(threefold.d3)
        ch1b = float(v.a1) - beta * rd
        ch2b = float(v.a2) - beta * float(v.a1) + (beta * beta / 2) * rd
        ch3b = (
            float(v.a3)
            - beta * float(v.a2)
            + (beta * beta / 2) * float(v.a1)
            - (beta ** 3 / 6) * rd
        )
        Qv = 4 * ch2b * ch2b - 6 * ch1b * ch3b
        betas.append(beta)
        alphas.append(math.sqrt(max(0.0, Qv)) / abs(ch1b))
    alg_cert = (
        Certificate(
            Rigor.PROVEN,
            ("threefold tilt-BG (BMT) inequality is a theorem here",),
            tuple(threefold.references),
            "tilt-BG boundary rigorous (bg_proven)",
        )
        if threefold.bg_proven else
        Certificate(
            Rigor.CONJECTURAL,
            ("threefold strong-BMT boundary is NOT a theorem here",),
            tuple(threefold.references),
            "strong BMT boundary not rigorous; stability conditions nonetheless exist",
        )
    )
    return BGBoundary(
        betas=betas,
        alphas=alphas,
        bg_proven=threefold.bg_proven,
        threefold_name=threefold.name,
        degenerate_betas=degenerate,
        certificate=meet(alg_cert, threefold.certificate),
    )


# ==========================================================================
# E4-M2 / G9a -- PROVEN tilt-layer (nu) wall solver `numerical_tilt_wall`
# ==========================================================================
#
# The tilt slope ``nu_{alpha,beta}`` (shipped as ``threefold.nu``, E4-M1)
# depends ONLY on the truncated triple ``(r, a1, a2) = (ch0, ch1.H^2, ch2.H)``
# and NEVER on ``a3 = ch3``: the numerator ``a2^b - 1/2 alpha^2 r d3`` and
# denominator ``a1^b`` of ``nu`` use only ``twist``'s ``a1^b, a2^b``, which are
# functions of ``(r, a1, a2)`` alone.  The tilt wall ``nu(v) = nu(w)`` is
# therefore EXACTLY the surface ``(s, t)`` semicircle of ``walls.numerical_wall``
# under the substitution ``(s, t) -> (beta, alpha)``, ``(r, c, e) -> (r, a1, a2)``,
# ``d -> d3``.  This introduces NO new mathematics: it re-reads the existing
# exact ``numerical_wall`` primitive on a truncated character.
#
# Two-way verification (twisted cubic on P^3, d3 = 1).
#   v = ThreefoldChern(1, 0, -3, 5)  (truncated (1, 0, -3)),
#   w = ThreefoldChern(1, -1, 1/2, -1/6)  (truncated (1, -1, 1/2), the tilt
#   image of O(-1) on P^3).
#   (a) numerical_wall minors on the truncated triples:
#       W_rc = 1*(-1) - 1*0    = -1
#       W_re = 1*(1/2) - 1*(-3) = 7/2
#       W_ce = 0*(1/2) - (-1)*(-3) = -3
#       center  = W_re/W_rc = (7/2)/(-1) = -7/2
#       radius^2 = center^2 - 2 W_ce/(d3 W_rc) = 49/4 - 2(-3)/(1(-1)) = 25/4
#   (b) clearing nu(v) = nu(w) directly with
#       nu(E) = (a2^b - 1/2 alpha^2 r d3)/a1^b,
#       a1^b = a1 - beta r d3, a2^b = a2 - beta a1 + 1/2 beta^2 r d3,
#       yields alpha^2 + beta^2 + 7 beta + 6 = 0, i.e. (beta + 7/2)^2 + alpha^2
#       = 25/4 -> center -7/2, radius^2 25/4.  Both routes agree exactly.
#
# General identity (same expansion): the raw tilt-wall conic is
#   -(d3/2) W_rc (alpha^2 + beta^2) + d3 W_re beta - W_ce = 0,
# whose alpha^2 coefficient (from the -1/2 alpha^2 r d3 tilt term) and beta^2
# coefficient (from the twist's +1/2 beta^2 r d3 plus the (-beta a1)(-beta r d3)
# cross term) are two independently-derived expressions that are both
# ``-(d3/2) W_rc`` -- hence equal (the beta-axis-centered-circle fact).  Center
# ``W_re/W_rc`` is d3-independent; radius^2 carries d3 via ``-2 W_ce/(d3 W_rc)``.
#
# References: Bayer-Macri-Toda arXiv:1103.5010 (tilt slope);
# Schmidt arXiv:1509.04608 (threefold tilt-stability walls).

_TILT_WALL_CERTIFICATE = Certificate(
    Rigor.PROVEN,
    ("tilt slope nu_{alpha,beta} depends only on (ch0, ch1.H^2, ch2.H): "
     "ch3-independence is an exact algebraic identity",
     "Picard rank 1 (H-numerical (r, ch1.H^2, ch2.H) reduction is exact)"),
    ("arXiv:1103.5010", "arXiv:1509.04608"),  # BMT tilt slope; Schmidt tilt walls
    "tilt-layer (nu) wall only, = walls.numerical_wall on the truncated triple; "
    "NOT the Bridgeland wall (needs 2nd tilt + BG: E4-M3/E5). Schmidt H-power "
    "normalization not yet cross-checked (see nu docstring).",
)


@dataclass(frozen=True)
class ThreefoldTiltWall:
    """A tilt-layer (nu) wall: a beta-axis-centered semicircle in the (beta, alpha)
    upper half-plane where ``nu_{alpha,beta}(v) = nu_{alpha,beta}(w)``.

    Mirrors :class:`walls.Wall`.  ``center`` and ``radius_sq`` are exact
    ``Fraction`` (no float in the invariant); ``radius`` casts to float only for
    display/geometry.  Because ``nu`` uses only the truncated triple
    ``(r, a1, a2)``, this wall is ``ch3``-independent by an exact algebraic
    identity and is computed by delegating to :func:`walls.numerical_wall` on the
    truncated triple with ``d -> d3``.

    ``bridgeland_certified`` is always ``False`` here: this is the tilt (nu) wall,
    not the full Bridgeland wall (the Bridgeland upgrade -- second tilt + BG --
    is E4-M3/E5).  The attached ``certificate`` is PROVEN for the tilt-layer
    reduction (Picard rank 1); it does NOT certify a Bridgeland wall.
    """

    center: Fraction
    radius_sq: Fraction
    subobject: ThreefoldChern
    v: ThreefoldChern
    bridgeland_certified: bool = False
    certificate: Certificate = UNKNOWN_CERTIFICATE

    @property
    def is_real(self) -> bool:
        return self.radius_sq >= 0

    @property
    def radius(self) -> float:
        return math.sqrt(float(self.radius_sq)) if self.radius_sq >= 0 else float("nan")

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (f"ThreefoldTiltWall(center={self.center}, radius_sq={self.radius_sq}, "
                f"R={self.radius:.4f}, w={self.subobject}, "
                f"bridgeland_certified={self.bridgeland_certified}) "
                f"[{self.certificate.rigor.name}]")


@dataclass(frozen=True)
class ThreefoldVerticalTiltWall:
    """The degenerate vertical tilt wall ``beta = beta_value`` (equal tilt-Mumford
    slopes, ``W_rc = v.r a1_w - w.r a1_v = 0``).

    The tilt-layer analogue of :class:`walls.VerticalWall`.  ``beta_value`` is the
    exact ``Fraction`` ``a1_v/(r_v d3)`` (or ``None`` when ``r_v = 0``).  Like
    :class:`walls.VerticalWall` it carries no ``certificate`` field.
    ``bridgeland_certified`` is always ``False`` (tilt wall, not Bridgeland).
    """

    beta_value: Optional[Fraction]
    subobject: ThreefoldChern
    v: ThreefoldChern
    bridgeland_certified: bool = False

    is_real = True
    radius = float("inf")


def tilt_wall_coeffs(v: ThreefoldChern, w: ThreefoldChern, d3: int
                     ) -> Tuple[Fraction, Fraction, Fraction, Fraction]:
    """Exact conic coefficients ``(coeff_alpha_sq, coeff_beta_sq, coeff_beta,
    coeff_const)`` of the raw tilt-wall equation

        coeff_alpha_sq * alpha^2 + coeff_beta_sq * beta^2
            + coeff_beta * beta + coeff_const = 0,

    derived by clearing ``nu_{alpha,beta}(v) = nu_{alpha,beta}(w)`` (see the
    module-section derivation).  With the three ``numerical_wall`` minors on the
    truncated triples ``(r, a1, a2)``,

        W_rc = r_v a1_w - r_w a1_v,
        W_re = r_v a2_w - r_w a2_v,
        W_ce = a1_v a2_w - a1_w a2_v,

    the coefficients are ``-(d3/2) W_rc`` (both alpha^2 and beta^2, from two
    independently-derived expressions that coincide -- the beta-axis-centered
    fact), ``d3 W_re`` (beta), and ``-W_ce`` (constant).  Pure ``Fraction``.
    ``coeff_alpha_sq == coeff_beta_sq == 0`` exactly when ``W_rc = 0`` (the
    degenerate vertical-wall case).
    """
    W_rc = v.r * w.a1 - w.r * v.a1
    W_re = v.r * w.a2 - w.r * v.a2
    W_ce = v.a1 * w.a2 - w.a1 * v.a2
    coeff_alpha_sq = -Fraction(d3, 2) * W_rc   # from the -1/2 alpha^2 r d3 tilt term
    coeff_beta_sq = -Fraction(d3, 2) * W_rc    # from the twist's beta^2 structure
    coeff_beta = Fraction(d3) * W_re
    coeff_const = -W_ce
    return coeff_alpha_sq, coeff_beta_sq, coeff_beta, coeff_const


def numerical_tilt_wall(v: ThreefoldChern, w: ThreefoldChern, d3: int):
    """The tilt (nu) wall ``nu_{alpha,beta}(v) = nu_{alpha,beta}(w)`` as a
    :class:`ThreefoldTiltWall` (semicircle) or :class:`ThreefoldVerticalTiltWall`.

    ``ch3``-INDEPENDENCE (exact identity).  ``nu`` uses only the truncated triple
    ``(r, a1, a2)``, so only ``(v.r, v.a1, v.a2)`` and ``(w.r, w.a1, w.a2)`` enter
    -- ``v.a3`` and ``w.a3`` are never read.  The wall is computed by delegating
    to :func:`walls.numerical_wall` on the truncated ``ChernChar(r, a1, a2)`` with
    the surface degree ``d`` replaced by ``d3``: it is the SAME three 2x2 minors,
    so ``center`` and ``radius_sq`` are bit-for-bit equal to the surface wall on
    the truncated triple (the reduction identity).

    Returns a :class:`ThreefoldTiltWall` (frozen; exact ``Fraction`` center and
    radius_sq, float only in the ``radius`` display property) with
    ``bridgeland_certified=False`` and the PROVEN tilt-layer certificate, or a
    :class:`ThreefoldVerticalTiltWall` when ``W_rc = 0`` (equal tilt-Mumford
    slopes).  This is the tilt (nu) wall only, NOT the Bridgeland wall (which
    needs the second tilt + BG; E4-M3/E5) -- hence ``bridgeland_certified`` is
    always ``False``.
    """
    vt = ChernChar(v.r, v.a1, v.a2)
    wt = ChernChar(w.r, w.a1, w.a2)
    geom = numerical_wall(vt, wt, d3)
    if isinstance(geom, VerticalWall):
        return ThreefoldVerticalTiltWall(beta_value=geom.s_value, subobject=w, v=v)
    return ThreefoldTiltWall(center=geom.center, radius_sq=geom.radius_sq,
                             subobject=w, v=v, bridgeland_certified=False,
                             certificate=_TILT_WALL_CERTIFICATE)


# ==========================================================================
# E4-M3 / G9b -- second-tilt / Bridgeland wall `bridgeland_wall` [RESEARCH-light]
# ==========================================================================
#
# The full second-tilt central charge on a threefold (H-normalized H-degrees
# ``(r, a1, a2, a3) = (ch0, ch1.H^2, ch2.H, ch3)``, ``A = alpha^2``) is
#
#     Z_{alpha,beta}(E) = (-a3^b + (A/2) a1^b) + i (a2^b - (A/2) r d3),
#
# with ``a1^b, a2^b, a3^b`` the beta-twist components (see ``ThreefoldChern.twist``).
# This is the SECOND tilt: unlike the tilt slope ``nu`` (E4-M1/E4-M2), whose wall
# reads only the truncated triple ``(r, a1, a2)`` and is ``ch3``-independent, the
# Bridgeland wall's real part carries ``a3 = ch3``.
#
# HAND DERIVATION (why the naive locus is NOT a beta-axis circle).  Expanding the
# wall ``F = Re Z(v) Im Z(w) - Im Z(v) Re Z(w) = 0`` in full by hand, with the
# five 2x2 truncated-triple minors
#
#     W_rc = r_v a1_w - r_w a1_v      # (r, a1)  -- also the tilt wall's alpha^2 coeff
#     W_re = r_v a2_w - r_w a2_v      # (r, a2)
#     W_ce = a1_v a2_w - a1_w a2_v    # (a1, a2)
#     W_rs = r_v a3_w - r_w a3_v      # (r, a3)   -- enters only with ch3
#     W_as = a1_v a3_w - a1_w a3_v    # (a1, a3)  -- enters only with ch3
#
# gives a DEGREE-4 curve in (alpha, beta): the ``A^2 = alpha^4`` coefficient is
# ``(d3/4) W_rc`` (nonzero), the ``A = alpha^2`` coefficient collapses (all
# beta-dependence cancels) to ``(1/2)(W_ce - d3 W_rs)``, the ``beta^2`` coefficient
# is ``(1/2)(W_ce + d3 W_rs)``, the ``beta^3`` coefficient is
# ``(d3/3)(a2_v r_w - r_v a2_w)`` and the ``beta^4`` coefficient is ``(d3/12) W_rc``.
# So once ``ch3`` enters (``W_rs != 0``) the naive ``Im(Z_v conj(Z_w)) = 0`` locus
# is a quartic whose ``alpha^2`` and ``beta^2`` coefficients DIFFER by ``d3 W_rs``:
# it is NOT the beta-axis-centered circle ``x(alpha^2+beta^2) + y beta + z = 0`` the
# roadmap describes.  That circle requires Schmidt's specific second-tilt
# reduction / normalization (arXiv:1509.04608 sec.3/sec.7) -- which is exactly why
# the roadmap marks the ``(y, z)`` transcription and the twisted-cubic numeric pin
# ``[RESEARCH-light]`` / DEFERRED.
#
# CANDIDATE (y, z) SHIPPED BEHIND THE FLAG.  We ship the STRUCTURE + gate +
# enumerator and a clearly-flagged CANDIDATE ``(y, z)`` chosen so that: (a) ``x``
# comes ONLY from the ``(r, a1)`` minor (as in the tilt wall, so the locus is a
# beta-axis circle by construction); (b) ``ch3`` enters ``y, z`` via the two
# real-part minors ``W_rs = (r, a3)`` and ``W_as = (a1, a3)``; and (c) it reduces
# EXACTLY to the E4-M2 tilt wall when the two ch3-minors vanish:
#
#     x = -(d3/2) W_rc
#     y =  d3 W_re + d3 W_rs
#     z = -W_ce - W_as
#
# ``W_rs = W_as = 0`` (e.g. ``a3_v = a3_w = 0``) => ``(x, y, z) =
# (-(d3/2) W_rc, d3 W_re, -W_ce)`` = the tilt wall, so ``center = -y/(2x)`` and
# ``radius_sq = (y^2 - 4 x z)/(4 x^2)`` collapse to the E4-M2 / Schmidt-Thm-3.3
# values.  Hand recompute (twisted cubic ``V=(1,0,-3,5)``, ``W=(1,-1,1/2,-1/6)``,
# ``d3=1``): ``x=1/2, y=-5/3, z=-2 -> center=5/3, radius_sq=61/9`` -- these are
# CANDIDATE values, deliberately NOT asserted against Schmidt; with ``a3=0`` both,
# ``center=-7/2, radius_sq=25/4`` = the E4-M2 pin.
#
# IMPORTANT (provenance).  Do NOT cite Schmidt Thm 3.3 for the ch3-dependent
# semicircle: the ``+d3 W_rs`` / ``-W_as`` ch3 entries of ``(y, z)`` are a
# CANDIDATE, and the exact ``(y, z)`` is DEFERRED to Schmidt sec.3/sec.7
# (arXiv:1509.04608).  The ``x`` coefficient (from the ``(r, a1)`` minor) is exact,
# as in the tilt wall.  ``bridgeland_certified`` reflects ``threefold.bg_proven``
# (a genuine Bridgeland wall exists only where the threefold BG inequality is a
# THEOREM; Piyaratne-Toda arXiv:1504.01177) -- NOT a claim that ``(y, z)`` is
# verified.
#
# References: Schmidt arXiv:1509.04608 (second-tilt walls); Bayer-Macri-Toda
# arXiv:1103.5010 (tilt slope / BMT); Piyaratne-Toda arXiv:1504.01177 (moduli /
# properness of Bridgeland stability on threefolds, conditional on BG).

_BRIDGELAND_WALL_CERTIFICATE = Certificate(
    Rigor.CONJECTURAL,
    ("x coeff from the (r,a1) minor is exact (as in the tilt wall)",
     "second-tilt Z_{alpha,beta}: ch3-dependent (y,z) coefficients are a CANDIDATE, "
     "NOT yet transcribed from Schmidt arXiv:1509.04608 sec.3/sec.7",
     "the beta-axis CIRCLE SHAPE is ITSELF a candidate hypothesis: the cubic "
     "ch3^beta term means the true (alpha,beta) locus is NOT automatically a "
     "circle and may be a higher-degree curve -- locus_verified=False until "
     "confirmed against Schmidt sec.3/sec.7"),
    ("arXiv:1509.04608", "arXiv:1103.5010", "arXiv:1504.01177"),
    "Bridgeland (second-tilt) wall LOCUS is a candidate pending Schmidt (y,z) "
    "transcription AND shape confirmation (the circle form may not be the true "
    "locus). bridgeland_certified reflects threefold.bg_proven (BG theorem; "
    "Piyaratne-Toda arXiv:1504.01177), NOT locus/shape verification. Not attributed to Thm 3.3.",
)


@dataclass(frozen=True)
class BridgelandWall:
    """A second-tilt (Bridgeland) wall: the beta-axis-centered semicircle

        x (alpha^2 + beta^2) + y beta + z = 0

    of the full central charge ``Z_{alpha,beta}`` in the ``(beta, alpha)`` upper
    half-plane.  Mirrors :class:`ThreefoldTiltWall`, adding the exact fitted
    coefficients ``(x, y, z)``.  ``center`` and ``radius_sq`` are exact
    ``Fraction`` (no float in the invariant); ``radius`` casts to float only for
    display/geometry.

    Unlike the tilt (nu) wall, this reads ``a3 = ch3``: the ``x`` coefficient is
    the exact ``(r, a1)`` minor (identical to the tilt wall, which is why the
    locus stays a beta-axis circle), while ``y, z`` carry ``ch3`` through the
    ``(r, a3)`` and ``(a1, a3)`` minors.  The ``ch3``-dependent ``(y, z)`` is a
    CANDIDATE pending transcription from Schmidt arXiv:1509.04608 sec.3/sec.7 (see
    the module-section derivation); the attached ``certificate`` is CONJECTURAL and
    is NOT attributed to Thm 3.3.  ``bridgeland_certified`` echoes
    ``threefold.bg_proven`` (a genuine Bridgeland wall only where the threefold BG
    inequality is a theorem).
    """

    center: Fraction
    radius_sq: Fraction
    x: Fraction            # alpha^2/beta^2 coeff, from the (r, a1) minor
    y: Fraction            # beta coeff (carries ch3)
    z: Fraction            # constant (carries ch3)
    subobject: ThreefoldChern
    v: ThreefoldChern
    bridgeland_certified: bool = False
    # locus_verified is about THIS computed wall's shape+coeffs (Schmidt-deferred),
    # distinct from bridgeland_certified (about the threefold's BG theorem).  It is
    # False until the (y,z) circle candidate is confirmed against Schmidt sec.3/sec.7
    # (the cubic ch3^beta term means the true locus is not automatically a circle).
    locus_verified: bool = False
    certificate: Certificate = UNKNOWN_CERTIFICATE

    @property
    def is_real(self) -> bool:
        return self.radius_sq >= 0

    @property
    def radius(self) -> float:
        return math.sqrt(float(self.radius_sq)) if self.radius_sq >= 0 else float("nan")

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        marker = "" if self.locus_verified else " CANDIDATE-LOCUS"
        return (f"BridgelandWall(center={self.center}, radius_sq={self.radius_sq}, "
                f"R={self.radius:.4f}, x={self.x}, y={self.y}, z={self.z}, "
                f"w={self.subobject}, bridgeland_certified={self.bridgeland_certified}, "
                f"locus_verified={self.locus_verified}) "
                f"[{self.certificate.rigor.name}]{marker}")


def bridgeland_wall_coeffs(v: ThreefoldChern, w: ThreefoldChern, d3: int
                           ) -> Tuple[Fraction, Fraction, Fraction]:
    """CANDIDATE ``(x, y, z)`` coefficients of the second-tilt Bridgeland wall UNDER
    THE HYPOTHESIS that its ``(alpha, beta)`` locus is the beta-axis circle

        x (alpha^2 + beta^2) + y beta + z = 0.

    [RESEARCH-light] WARNING.  This circle SHAPE is itself unverified: the cubic
    ``ch3^beta`` term in the second-tilt central charge means the true locus is NOT
    automatically a circle (it may be a higher-degree curve).  Confirm the shape and
    transcribe the exact ``(y, z)`` against Schmidt arXiv:1509.04608 sec.3/sec.7
    before trusting either.  With the five truncated-triple minors

        W_rc = r_v a1_w - r_w a1_v,   W_re = r_v a2_w - r_w a2_v,
        W_ce = a1_v a2_w - a1_w a2_v, W_rs = r_v a3_w - r_w a3_v,
        W_as = a1_v a3_w - a1_w a3_v,

    the coefficients are ``x = -(d3/2) W_rc`` (from the ``(r, a1)`` minor only --
    identical to the tilt wall's alpha^2/beta^2 coefficient, so ``x`` stays exact),
    ``y = d3 (W_re + W_rs)`` and ``z = -(W_ce + W_as)`` (the ``W_rs``/``W_as``
    entries carry ``ch3`` and are the CANDIDATE part; see the section derivation).
    Pure ``Fraction``.  ``x == 0`` exactly when ``W_rc = 0`` (vertical case).
    """
    W_rc = v.r * w.a1 - w.r * v.a1        # (r, a1)
    W_re = v.r * w.a2 - w.r * v.a2        # (r, a2)
    W_ce = v.a1 * w.a2 - w.a1 * v.a2      # (a1, a2)
    W_rs = v.r * w.a3 - w.r * v.a3        # (r, a3)  -- ch3
    W_as = v.a1 * w.a3 - w.a1 * v.a3      # (a1, a3) -- ch3
    x = -Fraction(d3, 2) * W_rc
    y = Fraction(d3) * W_re + Fraction(d3) * W_rs
    z = -W_ce - W_as
    return x, y, z


def bridgeland_wall(v: ThreefoldChern, w: ThreefoldChern, threefold: Threefold):
    """The second-tilt (Bridgeland) wall of ``v`` destabilized by ``w`` as a
    :class:`BridgelandWall` (semicircle) or :class:`ThreefoldVerticalTiltWall`.

    Uses the full central charge ``Z_{alpha,beta}`` (so ``a3 = ch3`` enters via the
    ``(y, z)`` coefficients of the CANDIDATE circle ``x(alpha^2 + beta^2) + y beta +
    z = 0``); the ``x`` coefficient still comes ONLY from the ``(r, a1)`` minor.
    [RESEARCH-light]: the beta-axis circle SHAPE is a candidate hypothesis -- the
    cubic ``ch3^beta`` term means the true locus is not automatically a circle (see
    ``bridgeland_wall_coeffs``); the returned :class:`BridgelandWall` carries
    ``locus_verified=False`` until confirmed against Schmidt sec.3/sec.7.  Returns

        center = -y/(2x),   radius_sq = (y^2 - 4 x z)/(4 x^2)

    as exact ``Fraction`` (``radius`` is float only via the display property), or a
    :class:`ThreefoldVerticalTiltWall` at ``beta = a1_v/(r_v d3)`` when ``x = 0``
    (i.e. ``W_rc = 0``, equal tilt-Mumford slopes -- read off the ``(r, a1)`` minor).

    ``bridgeland_certified = threefold.bg_proven``: BG is a THEOREM on this THREEFOLD
    (True for P^3 / QUADRIC3 / all Fano rho=1 / abelian / quintic; False for
    Bl_p(P^3)), so genuine Bridgeland walls EXIST here.  Do NOT conflate that with
    ``locus_verified`` (always ``False`` here), which is about whether THIS computed
    locus -- its shape AND its ``ch3``-dependent ``(y, z)`` -- is verified: it is a
    CANDIDATE pending Schmidt arXiv:1509.04608 sec.3/sec.7.  The attached
    ``certificate`` is CONJECTURAL and is NOT attributed to Thm 3.3 (see
    ``bridgeland_wall_coeffs`` and the module-section derivation).
    """
    d3 = threefold.d3
    x, y, z = bridgeland_wall_coeffs(v, w, d3)
    if x == 0:  # W_rc == 0 -> vertical beta = const (from the (r, a1) minor)
        beta = Fraction(v.a1, v.r * d3) if v.r != 0 else None
        return ThreefoldVerticalTiltWall(beta_value=beta, subobject=w, v=v,
                                         bridgeland_certified=threefold.bg_proven)
    center = -y / (2 * x)
    radius_sq = (y * y - 4 * x * z) / (4 * x * x)
    return BridgelandWall(center, radius_sq, x, y, z, w, v,
                          bridgeland_certified=threefold.bg_proven,
                          certificate=meet(_BRIDGELAND_WALL_CERTIFICATE,
                                           threefold.certificate))


def enumerate_tilt_walls(
    v: ThreefoldChern,
    threefold: Threefold,
    *,
    rank_bound: Optional[int] = None,
    degree_bound: int = 6,
    center_window: Number = 12,
    include_torsion: bool = True,
) -> List[ThreefoldTiltWall]:
    """Enumerate the PROVEN, ``ch3``-independent tilt (nu) walls of ``v``.

    The threefold analogue of :func:`walls.actual_walls` / :func:`walls.compute_walls`,
    working in the truncated ``(r, a1, a2)`` space with the surface degree ``d``
    replaced by ``d3 = H^3``.  For each candidate sub-rank ``r'`` in
    ``[r_lo, min(rank_bound, r_v)]`` (``include_torsion`` adds ``r' = 0``), the
    ``(r, a1)`` minor ``W_rc = r_v a1' - r' a1_v`` is swept over
    ``[-degree_bound, degree_bound] \\ {0}`` (each value fixing the integer ``a1'``),
    and the half-integer ``a2'`` placing the wall center ``W_re/W_rc`` inside the
    ``beta`` window ``[mu - center_window, mu + center_window]`` (``mu = a1_v/(r_v d3)``)
    is swept; the truncated subobject ``ThreefoldChern(r', a1', a2', 0)`` is fed to
    :func:`numerical_tilt_wall`.  Real semicircles (``radius_sq > 0``) are kept,
    de-duplicated by ``(center, radius_sq)``, and sorted by ``radius_sq`` (largest
    first, the outermost wall).

    This enumerates the tilt (nu) layer only, which is PROVEN and ``ch3``-independent
    (``numerical_tilt_wall`` reads only ``(r, a1, a2)``) -- hence every returned
    :class:`ThreefoldTiltWall` has ``bridgeland_certified = False``.  It can be
    upgraded to :func:`bridgeland_wall` once the Schmidt ``(y, z)`` transcription
    (arXiv:1509.04608 sec.3/sec.7) lands (E4-M3 [RESEARCH-light] / E5).
    """
    if v.r <= 0:
        raise ValueError("enumerate_tilt_walls expects a positive-rank v")
    d3 = threefold.d3
    mu_v = Fraction(v.a1, v.r * d3)
    S = Q(center_window)
    s_lo, s_hi = mu_v - S, mu_v + S
    r_cap = v.r if rank_bound is None else min(rank_bound, v.r)
    r_lo = 0 if include_torsion else 1

    best: dict = {}  # (center, radius_sq) -> ThreefoldTiltWall
    for rp in range(r_lo, r_cap + 1):
        for W_int in range(-degree_bound, degree_bound + 1):
            if W_int == 0:
                continue  # equal tilt-Mumford slopes -> vertical wall (skip)
            a1p = (Fraction(W_int) + rp * v.a1) / v.r
            if a1p.denominator != 1:
                continue  # no integer a1' realizes this W_rc
            if rp == 0 and a1p <= 0:
                continue  # a torsion subobject needs effective degree a1' > 0
            W_rc = Fraction(W_int)
            # a2' placing the center s0 = (v.r a2' - rp v.a2)/W_rc inside the window
            e_a = Fraction(s_lo * W_rc + rp * v.a2, v.r)
            e_b = Fraction(s_hi * W_rc + rp * v.a2, v.r)
            e_lo, e_hi = (e_a, e_b) if e_a <= e_b else (e_b, e_a)
            for a2p in _half_integers(e_lo, e_hi):
                w = ThreefoldChern(rp, a1p, a2p, Fraction(0))
                wall = numerical_tilt_wall(v, w, d3)
                if isinstance(wall, ThreefoldVerticalTiltWall):
                    continue
                if wall.radius_sq <= 0:
                    continue
                if not (s_lo <= wall.center <= s_hi):
                    continue
                key = (wall.center, wall.radius_sq)
                if key not in best:
                    best[key] = wall
    out = list(best.values())
    out.sort(key=lambda W: (-W.radius_sq, W.center))
    return out


# ------------------------------------------------------------------
# E5-M1 / G15 -- conjecture-gated Bridgeland-wall verdict + rigor propagation
# ------------------------------------------------------------------
# The tilt (nu) wall locus is PROVEN (numerical_tilt_wall: the reduction to
# walls.numerical_wall is an exact identity, Rigor.PROVEN).  Whether that tilt
# wall LIFTS to a genuine Bridgeland wall requires the threefold BMT/BG
# inequality Q >= 0 to be a THEOREM -- i.e. threefold.bg_proven.  So the verdict
# rigor is min(TILT_SOLVER_RIGOR, PROVEN if bg_proven else CONJECTURAL).  Pure
# additive labeling: NO Fraction value changes.
#
# References: BMT arXiv:1103.5010; Schmidt arXiv:1509.04608 (tilt walls);
# Piyaratne-Toda arXiv:1504.01177 (moduli properness, CONDITIONAL on BG -- this
# hypothesis must be carried even on bg_proven rows).

TILT_SOLVER_RIGOR: Rigor = _TILT_WALL_CERTIFICATE.rigor   # == Rigor.PROVEN

_PIYARATNE_TODA = "arXiv:1504.01177"   # conditional moduli properness (assuming BG)

_VERDICT_UNPROVEN_NOTE = (
    "tilt-stability wall only; Bridgeland upgrade unproven (threefold BG open here)"
)
_VERDICT_BLOWUP_EXTRA = (
    "strong BMT boundary FAILS (Schmidt arXiv:1602.05055); "
    "stability conditions nonetheless EXIST (Fano, BMSZ arXiv:1607.08199)"
)
_VERDICT_CERTIFIED_NOTE = (
    "Bridgeland wall certified: the threefold tilt-BG (BMT) inequality is a "
    "theorem here, so the PROVEN tilt (nu) wall lifts to a genuine Bridgeland "
    "wall; moduli properness is conditional on BG "
    "(Piyaratne-Toda arXiv:1504.01177)."
)


def is_bridgeland_certified(threefold: Threefold) -> bool:
    """True iff the threefold BMT/BG inequality is a THEOREM here, so a tilt wall
    lifts to a genuine Bridgeland wall.  A named accessor for ``threefold.bg_proven``
    (Piyaratne-Toda arXiv:1504.01177 gives conditional properness).
    ``bg_boundary_curve`` already attaches a matching PROVEN/CONJECTURAL
    ``Certificate``."""
    return threefold.bg_proven


@dataclass(frozen=True)
class BridgelandWallVerdict:
    """Rigor-tagged verdict on whether a tilt (nu) wall is a genuine Bridgeland wall.

    ``tilt_wall`` is the RAW PROVEN locus from :func:`numerical_tilt_wall`
    (unchanged: no ``Fraction`` value moves).  ``certified == threefold.bg_proven``;
    ``rigor`` is the meet ``min(TILT_SOLVER_RIGOR, PROVEN if bg_proven else
    CONJECTURAL)``.  Never emits 'no stability condition'.
    """

    tilt_wall: object          # ThreefoldTiltWall | ThreefoldVerticalTiltWall
    certified: bool
    rigor: Rigor
    citations: Tuple[str, ...]
    note: str

    @property
    def certificate(self) -> Certificate:
        return Certificate(self.rigor, (), self.citations, self.note)


def bridgeland_wall_verdict(
    v: ThreefoldChern, w: ThreefoldChern, threefold: Threefold
) -> BridgelandWallVerdict:
    """Wrap the PROVEN tilt (nu) wall of ``v`` destabilized by ``w`` with a
    Bridgeland verdict on ``threefold``.  Signature mirrors
    ``bridgeland_wall(v, w, threefold)``: ``d3 = threefold.d3``.  ``certified =
    threefold.bg_proven``; ``rigor = min(TILT_SOLVER_RIGOR, PROVEN if bg_proven
    else CONJECTURAL)``.  Pure metadata -- no ``Fraction`` value changes; never
    emits 'no stability condition'."""
    tilt_wall = numerical_tilt_wall(v, w, threefold.d3)   # raw PROVEN locus (unchanged)
    certified = threefold.bg_proven
    bg_rigor = Rigor.PROVEN if certified else Rigor.CONJECTURAL
    rigor = min(TILT_SOLVER_RIGOR, bg_rigor)               # the meet law, exactly
    citations = tuple(dict.fromkeys(                       # order-preserving union
        _TILT_WALL_CERTIFICATE.citations                  # arXiv:1103.5010, 1509.04608
        + tuple(threefold.references)                     # threefold BG refs
        + (_PIYARATNE_TODA,)                              # conditional properness (always)
    ))
    if certified:
        note = _VERDICT_CERTIFIED_NOTE
    else:
        note = _VERDICT_UNPROVEN_NOTE
        if threefold.kind == "fano-blowup":               # Bl_p(P^3): Fano, carries stab. cond.
            note = note + "; " + _VERDICT_BLOWUP_EXTRA
    return BridgelandWallVerdict(tilt_wall, certified, rigor, citations, note)
