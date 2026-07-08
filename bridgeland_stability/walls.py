"""Algorithm 3 - Bridgeland numerical walls in the (s, t) upper half-plane.

For a fixed Chern character ``v`` on a surface with ``d = H^2`` and a potential
sub/quotient object ``w``, the wall ``W(v, w)`` is the locus where
``arg Z_{s,t}(v) = arg Z_{s,t}(w)``.  Writing ``v = (r,c,e)``, ``w = (r',c',e')``
and the three 2x2 "Mukai-like" minors

    W_rc = r c' - r' c,   W_re = r e' - r' e,   W_ce = c e' - c' e,

the wall is, when ``W_rc != 0``, the semicircle centered on the s-axis at

    s0  = W_re / W_rc,
    rho^2 = s0^2 - 2 * W_ce / (d * W_rc) = (s0 - mu_v)^2 - 2 * Delta_v(CH),

and it is a genuine (real) wall iff ``rho^2 >= 0``.  When ``W_rc = 0`` (equal
Mumford slopes) the wall degenerates to the vertical line ``s = mu_v``.

This matches Coskun-Huizenga (survey, section 5: ``s = (mu1+mu2)/2 -
(Delta1-Delta2)/(mu1-mu2)``, ``rho^2 = (s-mu1)^2 - 2 Delta1``) and ABCH
(arXiv:1203.0316).  Cross-check: for P^2[2], ``v = ch(I_Z) = (1,0,-2)`` is
destabilized by ``O(-1) = (1,-1,1/2)`` along the unique wall with center
``-5/2`` and radius ``3/2`` -- reproduced exactly by ``numerical_wall``.

All centers/radii are exact (``Fraction``); ``radius`` casts to float only at
the very end for geometry/plotting.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product
from typing import List, Optional, Sequence, Tuple

from .chern import ChernChar, Number, Q
from .rigor import Rigor, Certificate, UNKNOWN_CERTIFICATE, meet
from .varieties import Surface, require_faithful_computation


@dataclass(frozen=True)
class Wall:
    """A semicircular numerical wall: center ``(s0, 0)`` and ``radius^2 = rho^2``."""

    center: Fraction
    radius_sq: Fraction
    subobject: ChernChar
    v: ChernChar
    certificate: Certificate = UNKNOWN_CERTIFICATE

    @property
    def is_real(self) -> bool:
        return self.radius_sq >= 0

    @property
    def radius(self) -> float:
        return math.sqrt(float(self.radius_sq)) if self.radius_sq >= 0 else float("nan")

    @property
    def s_left(self) -> float:
        return float(self.center) - self.radius

    @property
    def s_right(self) -> float:
        return float(self.center) + self.radius

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"Wall(center={self.center}, radius_sq={self.radius_sq}, "
            f"R={self.radius:.4f}, w={self.subobject}) "
            f"[{self.certificate.rigor.name}]"
        )


@dataclass(frozen=True)
class VerticalWall:
    """The degenerate wall ``s = s_value`` (equal Mumford slopes, ``W_rc = 0``)."""

    s_value: Optional[Fraction]
    subobject: ChernChar
    v: ChernChar

    is_real = True
    radius = float("inf")


def _minors(v: ChernChar, w: ChernChar, d: int) -> Tuple[Fraction, Fraction, Fraction]:
    """The three 2x2 "Mukai-like" minors of ``(v, w)`` with the H-degree
    ``c = <ch1, H>`` read through the NS-lattice pairing (E8-M3 / G12.3):
    ``cv = v.ch1_dot_H(d)``, ``cw = w.ch1_dot_H(d)``.

    On the Picard-rank-1 shim ``<ch1, H> == v.c`` exactly, so every minor is
    bit-for-bit identical to the closed form ``W_rc = r c' - r' c`` etc.  Routing
    through ``ch1_dot_H`` makes the H-projection an NS-lattice quantity; the
    full rho>=2 wall (which also uses the polarization-orthogonal component of
    ch1) is E9/G13, NOT this milestone.
    """
    cv = v.ch1_dot_H(d)          # <ch1_v, H> via the NS-lattice pairing (== v.c at rho=1)
    cw = w.ch1_dot_H(d)          # <ch1_w, H>
    W_rc = v.r * cw - w.r * cv
    W_re = v.r * w.e - w.r * v.e
    W_ce = cv * w.e - cw * v.e
    return W_rc, W_re, W_ce


def numerical_wall(v: ChernChar, w: ChernChar, d: int):
    """The wall ``W(v, w)`` as a :class:`Wall` (semicircle) or :class:`VerticalWall`.

    The returned :class:`Wall` may have ``radius_sq < 0`` (``is_real`` False),
    meaning ``w`` does not give a genuine wall for ``v`` (excluded by BG).
    """
    W_rc, W_re, W_ce = _minors(v, w, d)
    if W_rc == 0:
        s_val = v.slope(d) if v.r != 0 else None
        return VerticalWall(s_val, w, v)
    s0 = Fraction(W_re, W_rc)
    rho_sq = s0 * s0 - Fraction(2 * W_ce, d * W_rc)
    return Wall(center=s0, radius_sq=rho_sq, subobject=w, v=v)


def _ns_ch1(
    chern: ChernChar, surface: Surface, ch1: Optional[Sequence[Number]]
) -> Tuple[Fraction, ...]:
    """The ``ch1`` of ``chern`` as an NS-lattice vector on ``surface`` (E9-M1/G13).

    * ``ch1 is None`` -> the Picard-rank-1 shim ``chern.ch1(surface.d) = (c/d,)``
      (the H-coordinate of ``ch1`` in the ample-generator basis).  On a rank-1
      surface this is the genuine ``ch1``, so the ``gamma`` twist below collapses
      to ``0`` and :func:`numerical_wall_ns` reduces to :func:`numerical_wall`
      bit-for-bit; the ``ch1`` argument is only needed to distinguish classes
      that share ``<ch1, H>`` but differ in ``ch1^2`` (Picard rank >= 2).
    * an explicit vector -> coerced to ``Fraction`` and validated: it must have
      ``len == surface.lattice.rank`` and its H-projection
      ``<ch1, H> = surface.lattice.pairing(ch1, surface.H)`` must equal the
      stored scalar ``chern.c`` (the H-degree the scalar core carries), else
      ``ValueError`` -- the supplied bidegree has to be consistent with the class.
    """
    if ch1 is None:
        if surface.lattice.rank > 1:      # fail LOUD: the (c/d,) shim is rho=1 only
            raise ValueError(
                f"ch1 vector required on Picard-rank >= 2 surface {surface.name} "
                f"(NS rank {surface.lattice.rank}): the (c/d,) shim is only valid at "
                "rho = 1; pass v_ch1 / w_ch1 to distinguish classes sharing <ch1,H>"
            )
        return chern.ch1(surface.d)
    L = surface.lattice
    vec = tuple(Q(x) for x in ch1)
    if len(vec) != L.rank:
        raise ValueError(
            f"ch1 vector has length {len(vec)} but the NS lattice of "
            f"{surface.name} has rank {L.rank}"
        )
    proj = L.pairing(vec, surface.H)
    if proj != chern.c:
        raise ValueError(
            f"ch1 vector {vec} has H-projection <ch1,H> = {proj}, which does not "
            f"match the stored H-degree c = {chern.c} of the class; the supplied "
            "bidegree is inconsistent with the Chern character"
        )
    return vec


def numerical_wall_ns(
    v: ChernChar,
    w: ChernChar,
    surface: Surface,
    *,
    v_ch1: Optional[Sequence[Number]] = None,
    w_ch1: Optional[Sequence[Number]] = None,
):
    """Full Neron-Severi Bridgeland wall ``W(v, w)`` from the whole intersection
    form -- via ``<ch1, ch1>`` and the polarization-orthogonal ``gamma^2`` -- not
    merely the H-projection ``<ch1, H>`` (E9-M1 / G13).

    Why the scalar model is blind to ``ch1^2``.  On the pure ``(s, t)`` slice
    (``beta = sH``, ``omega = tH``) the central charge ``Z_{s,t}`` depends only on
    ``(r, c=<ch1,H>, e)``, so the ``(s, t)`` semicircle can never see ``ch1^2``.
    Two classes with equal ``<ch1, H>`` but different ``ch1^2`` are conflated by
    :func:`numerical_wall`.  To use the full form one must leave the ``RH`` slice.

    Maciocia ``beta = b*omega + gamma`` decomposition (arXiv:1202.4587).  Write
    ``beta = b*H + gamma`` with ``gamma`` H-orthogonal (``<gamma, H> = 0``, hence
    ``gamma^2 <= 0`` by the Hodge index theorem).  Substituting into ``Z_{beta,
    omega}`` leaves the imaginary part unchanged (``gamma _|_ H``) and shifts the
    real part by exactly a ``gamma``-twist of ``ch2``:

        e~(x) := ch2(x) - <ch1(x), gamma> + (r_x / 2) * gamma^2
               (= ch2 of ch(x) . e^{-gamma})

    Therefore the full-NS wall in this slice is the *scalar* wall computed on the
    ``gamma``-twisted characters ``(r, c, e~)``.  This function builds that twist
    and delegates the minor geometry to :func:`numerical_wall`, then re-wraps the
    result with the ORIGINAL ``v, w`` so ``wall.subobject`` / ``wall.v`` are the
    true classes.

    ``gamma`` convention (documented, [RESEARCH-light]).
    ``gamma := H^perp(r_v * ch1_w - r_w * ch1_v)`` -- the H-orthogonal part of the
    *relative* first Chern class.  ``gamma = 0`` exactly when the two ch1-slopes
    are H-parallel; nonzero for a genuine bidegree destabilizer.

    Picard-rank-1 reduction is bit-for-bit exact (automatic, not special-cased).
    At rho = 1 the H-orthogonal complement ``H^perp = {0}``, so ``gamma = 0``,
    ``e~ = e``, and the twisted characters equal the untwisted ones: the return is
    identical to ``numerical_wall(v, w, surface.d)`` -- same ``center``,
    ``radius_sq``, ``subobject``, ``v`` and default ``UNKNOWN_CERTIFICATE`` (see
    ``test_ns_reduces_to_scalar``).  At rho >= 2 two classes with equal ``<ch1,
    H>`` but different ``ch1^2`` induce different ``gamma`` and get DIFFERENT walls
    -- the gain the scalar solver cannot see (``test_ns_finds_bidegree_wall``).

    Parameters
    ----------
    v, w : ChernChar
        The reference class and the potential sub/quotient class.
    surface : Surface
        Supplies ``d = H^2``, the NS lattice, and the polarization vector ``H``.
        Torsion-canonical surfaces (Enriques / bielliptic) are refused by the G12
        guard :func:`require_faithful_computation`.
    v_ch1, w_ch1 : Sequence[Number], optional
        Explicit NS-lattice ``ch1`` vectors (length ``surface.lattice.rank``).
        When omitted, the Picard-rank-1 shim ``(c/d,)`` is used -- correct at rho
        = 1, where the result equals :func:`numerical_wall` exactly.  On a Picard-
        rank >= 2 surface a genuine bidegree destabilizer must be supplied here
        (its H-projection must equal the class's stored ``c``, else ``ValueError``).

    [RESEARCH-light].  The rho = 1 reduction is exact and convention-independent
    (``gamma = 0`` regardless of normalization).  The DEMONSTRATED, test-backed
    rho > 1 gain is that the full-NS solver distinguishes classes the scalar model
    conflates (different ``ch1^2`` -> different ``gamma`` -> different walls) and
    finds a real wall for a bidegree destabilizer that the ``beta = sH`` slice
    misses; the EXACT rho > 1 radius, however, depends on the normalization of the
    ``gamma`` convention above.  Until one exact rho > 1 wall value is read off a
    worked example in Maciocia arXiv:1202.4587, the constructed
    P^1 x P^1 presence/absence pair (``test_ns_finds_bidegree_wall``) is the
    acceptance evidence, and the rho > 1 value is pinned under THIS documented
    ``gamma = H^perp(relative-ch1)`` convention rather than treated as literature-
    anchored.  Everything stays exact ``fractions.Fraction`` (no float, no viz).
    """
    require_faithful_computation(surface)          # G12 guard (Enriques/bielliptic refused)
    tv, tw = _ns_twisted_pair(v, w, surface, v_ch1, w_ch1)
    res = numerical_wall(tv, tw, surface.d)
    if isinstance(res, VerticalWall):
        return VerticalWall(res.s_value, w, v)
    return Wall(center=res.center, radius_sq=res.radius_sq, subobject=w, v=v)


def _ns_twisted_pair(
    v: ChernChar,
    w: ChernChar,
    surface: Surface,
    v_ch1: Optional[Sequence[Number]],
    w_ch1: Optional[Sequence[Number]],
) -> Tuple[ChernChar, ChernChar]:
    """The Maciocia gamma-twisted SCALAR reps ``(r, c, e~)`` of ``v`` and ``w`` on
    ``surface`` (E9-M1/G13 body, factored out for reuse by :func:`enumerate_ns_walls`).

    ``e~(x) = ch2(x) - <ch1(x), gamma> + (r_x / 2) * gamma^2`` with the twist vector
    ``gamma = H^perp(r_v ch1_w - r_w ch1_v)`` (the H-orthogonal part of the relative
    ch1).  :func:`numerical_wall(tv, tw, d)` on the returned pair reproduces
    :func:`numerical_wall_ns`'s geometry bit-for-bit.  The twist preserves the scalar
    H-degree ``c`` and is additive in ``ch``, so ``t(v - w) = tv - tw`` componentwise
    (used by the ρ>=2 Bogomolov-on-both necessary conditions in
    :func:`enumerate_ns_walls`).  At Picard rank 1, ``H^perp = {0}`` so ``gamma = 0``,
    ``e~ = e`` and the pair is the untwisted ``(v, w)`` -- the exact reduction.
    """
    d, L = surface.d, surface.lattice
    ch1_v = _ns_ch1(v, surface, v_ch1)
    ch1_w = _ns_ch1(w, surface, w_ch1)
    rel = tuple(v.r * ch1_w[i] - w.r * ch1_v[i] for i in range(L.rank))
    c_rel = L.pairing(rel, surface.H)
    gamma = tuple(rel[i] - Fraction(c_rel, d) * surface.H[i] for i in range(L.rank))
    g2 = L.self_pairing(gamma)                     # gamma^2 <= 0 (Hodge index)
    e_v = v.e - L.pairing(ch1_v, gamma) + Fraction(v.r, 2) * g2
    e_w = w.e - L.pairing(ch1_w, gamma) + Fraction(w.r, 2) * g2
    return ChernChar(v.r, v.c, e_v), ChernChar(w.r, w.c, e_w)


def _scalar_wall_alg_certificate(surface: Surface) -> Certificate:
    """Algebraic-rigor certificate for the scalar-minor (s, t) wall geometry.

    The scalar-minor wall (three 2x2 ``(r, c, e)`` minors) is PROVEN-correct ONLY
    at Picard rank 1: there ``c = ch1.H`` fully determines ``ch1`` (``NS = Z H``),
    so :func:`numerical_wall` computes the genuine (s, t) semicircle.  At
    ``picard_rank >= 2`` the scalar ``c = ch1.H`` is only the H-*projection* of
    ``ch1`` -- the polarization-orthogonal component of ``ch1`` is discarded, so
    the wall is only the H-projected one (Maciocia arXiv:1202.4587), not reliable.
    The PROVEN tag is therefore gated on ``surface.picard_rank == 1``; a rho >= 2
    surface falls back to a HEURISTIC cert flagging the H-projection and the
    NS-lattice refactor (E8/G12) that would fix it.
    """
    if surface.picard_rank == 1:
        return Certificate(
            Rigor.PROVEN,
            ("Picard rank 1",
             "scalar-minor (s,t) semicircle exact at Picard rank 1"),
            ("arXiv:1203.0316", "arXiv:1202.4587"),  # ABCH semicircle; Maciocia rho=1 nested walls
            "single (s,t)-wall geometry proven at Picard rank 1",
        )
    return Certificate(
        Rigor.HEURISTIC,
        ("Picard rank >= 2: scalar c=ch1.H is only the H-projection of ch1",),
        ("arXiv:1202.4587",),
        "H-projected wall only; needs the NS-lattice refactor (E8/G12)",
    )


def abelian_wall(v: ChernChar, w: ChernChar, surface: Surface):
    """The abelian-surface Bridgeland (s, t) wall ``W(v, w)`` -- the BARE wall.

    An abelian surface has trivial tangent bundle, so ``sqrt(td) = (1, 0, 0)`` and
    the Mukai vector is the *bare* Chern triple ``(r, c, ch2)`` with NO shift.
    Feeding the bare ``(r, c, e)`` to :func:`numerical_wall` therefore already gives
    the genuine abelian (s, t) semicircle: this returns exactly
    ``numerical_wall(v, w, surface.d)`` with NO Mukai shift applied -- in contrast to
    the K3 ``ch2 -> ch2 + r`` shim (``mukai.from_chern`` / ``mukai.k3_wall``, G3),
    which would inject a spurious ``+2/d`` here.

    For the shared fixture ``v=(1,0,-2)``, ``w=(1,-1,1/2)`` on ``d=2`` this is center
    ``-5/2``, radius^2 ``17/4`` -- exactly ``2/d = 1`` LESS than the K3 value ``21/4``.

    Only the single-semicircle primitive is a free win.  Do NOT route abelian classes
    through :func:`compute_walls` / :func:`actual_walls`: their ``_is_integral_chern``
    filter tests the ``d = 1`` P^2 Chern-lattice condition ``c1^2/2 - ch2 in Z``, which
    does NOT transfer to an abelian surface (``d != 1``), and ``actual_walls`` layers
    ABCH / P^2-specific necessary conditions -- so they are NOT abelian completeness
    oracles.  Picard-rank >= 2 abelian surfaces (``E_1 x E_2``) are gated on the
    NS-lattice refactor (E8 / G12) and are not covered here.

    Provenance (E2-M5 / G5).  The returned :class:`Wall` carries
    ``meet(_scalar_wall_alg_certificate(surface), surface.certificate)``: PROVEN at
    ``surface.picard_rank == 1`` (the scalar minors are the genuine wall), downgraded
    to HEURISTIC with an "H-projected; needs the NS-lattice refactor (E8/G12)" note
    at ``picard_rank >= 2`` (the scalar ``c = ch1.H`` is only the H-projection).  A
    degenerate equal-slope result is a :class:`VerticalWall`, which has NO
    ``certificate`` field, so it stays untagged (documented limitation).

    Raises ValueError unless ``surface.kind == "abelian"`` (guarding against the K3
    ``+r`` shift being applied to abelian input, and vice versa).
    """
    if surface.kind != "abelian":
        raise ValueError(
            "abelian_wall requires an abelian surface (surface.kind == 'abelian'); "
            f"got kind={surface.kind!r}. Use the K3 shim for K3 surfaces; "
            "numerical_wall(v, w, d) is the bare, surface-agnostic primitive."
        )
    # G12 guard (defensive; catalog abelian rows have canonical_order=0, so this
    # never fires there -- it refuses a hand-built torsion-canonical abelian input).
    require_faithful_computation(surface)
    wall = numerical_wall(v, w, surface.d)
    if isinstance(wall, VerticalWall):
        return wall
    cert = meet(_scalar_wall_alg_certificate(surface), surface.certificate)
    return replace(wall, certificate=cert)


def _half_integers(lo: Fraction, hi: Fraction):
    """Yield all half-integers (elements of (1/2)Z) in ``[lo, hi]``."""
    start = math.ceil(2 * lo)
    end = math.floor(2 * hi)
    for k in range(start, end + 1):
        yield Fraction(k, 2)


def _bg_ok(ch: ChernChar, d: int) -> bool:
    """Bogomolov necessary condition for a (sub/quotient) class to bound a real wall.

    Rank 0 (torsion) classes are admissible iff nonzero; positive/negative rank
    classes must satisfy ``Delta >= 0``.
    """
    if ch.r == 0:
        return ch.c != 0 or ch.e != 0
    return ch.discriminant(d) >= 0


def compute_walls(
    v: ChernChar,
    surface: Surface,
    s_range: Tuple[Number, Number] = (-6, 6),
    rank_bound: int = 6,
    degree_bound: int = 6,
    include_torsion: bool = True,
    require_bg: bool = True,
) -> List[Wall]:
    """Enumerate real numerical walls for ``v`` with center in ``s_range``.

    Strategy (exact, bounded).  For each candidate sub-rank ``r'`` in
    ``[0, min(rank_bound, rank v)]`` (``include_torsion`` adds the ``r'=0``
    curve-class destabilizers), the Mukai-degree ``W_rc = r_v c' - r' c_v`` is
    swept over ``[-degree_bound, degree_bound]`` (the natural finiteness
    parameter, since the wall center is ``W_re/W_rc``); each value fixes ``c'``,
    and the ``e'`` placing the center inside ``s_range`` form a finite interval
    of half-integers.  A wall is kept iff ``rho^2 >= 0``, its center is in
    ``s_range``, and (when ``require_bg``) BOTH the subobject ``w`` and the
    quotient ``v - w`` satisfy Bogomolov.  Walls are de-duplicated by
    ``(center, radius_sq)`` and returned sorted by radius (largest first).

    This returns *numerical/potential* walls (the BG-on-both necessary
    condition); certifying which are *actual* destabilizing walls requires the
    ABCH/Maciocia subobject analysis.  ``numerical_wall`` is the exact primitive
    for a single subobject; widen ``degree_bound``/``rank_bound`` for breadth.
    """
    require_faithful_computation(surface)  # G12 guard: torsion-canonical rows refused
    if v.r <= 0:
        raise ValueError("compute_walls expects a positive-rank v")
    d = surface.d
    s_min, s_max = Q(s_range[0]), Q(s_range[1])
    if s_min > s_max:
        raise ValueError("s_range must be (s_min, s_max) with s_min <= s_max")
    seen = set()
    walls: List[Wall] = []

    r_lo = 0 if include_torsion else 1
    r_hi = min(rank_bound, v.r)
    for rp in range(r_lo, r_hi + 1):
        # Bound the Mukai-degree W_rc = v.r*c' - rp*v.c by degree_bound -- the
        # natural finiteness parameter (the wall center is W_re/W_rc).  For each
        # admissible W_rc, recover the unique integer c' (when one exists).
        for W_int in range(-degree_bound, degree_bound + 1):
            if W_int == 0:
                continue  # equal slopes -> vertical wall (use numerical_wall directly)
            cp = (Fraction(W_int) + rp * v.c) / v.r
            if cp.denominator != 1:
                continue  # no integer c' realizes this W_rc
            if rp == 0 and cp <= 0:
                continue  # a torsion subobject needs effective degree c' > 0
            W_rc = Fraction(W_int)
            # e' placing the center s0 = (v.r e' - rp v.e)/W_rc inside [s_min, s_max]
            e_a = Fraction(s_min * W_rc + rp * v.e, v.r)
            e_b = Fraction(s_max * W_rc + rp * v.e, v.r)
            e_lo, e_hi = (e_a, e_b) if e_a <= e_b else (e_b, e_a)
            for ep in _half_integers(e_lo, e_hi):
                w = ChernChar(rp, cp, ep)
                if require_bg:
                    q = ChernChar(v.r - rp, v.c - cp, v.e - ep)
                    if not (_bg_ok(w, d) and _bg_ok(q, d)):
                        continue
                wall = numerical_wall(v, w, d)
                if isinstance(wall, VerticalWall):
                    continue
                if not wall.is_real:
                    continue
                if not (s_min <= wall.center <= s_max):
                    continue
                key = (wall.center, wall.radius_sq)
                if key in seen:
                    continue
                seen.add(key)
                walls.append(wall)

    walls.sort(key=lambda W: (-W.radius_sq, W.center))
    alg_cert = Certificate(
        Rigor.HEURISTIC,
        ("dense numerical enumeration; BG-on-both necessary only",),
        (),
        "compute_walls returns potential/numerical walls, not the certified actual-wall set",
    )
    cert = meet(alg_cert, surface.certificate)  # meet only lowers -> never PROVEN
    return [replace(w, certificate=cert) for w in walls]


def walls_from_subobjects(
    v: ChernChar, subobjects: List[ChernChar], surface: Surface
) -> List[Wall]:
    """Exact walls of ``v`` for an explicit list of destabilizing classes.

    The honest way to reproduce a published wall diagram: supply the known
    destabilizers (e.g. ``O(-1) = ChernChar(1,-1,1/2)`` for P^2[2]) and get the
    exact semicircles, with no enumeration heuristic.  Non-real walls
    (``radius_sq < 0``) and vertical walls are included as returned by
    :func:`numerical_wall`; filter with ``[w for w in result if getattr(w,
    'is_real', False)]`` as needed.
    """
    return [numerical_wall(v, w, surface.d) for w in subobjects]


# --------------------------------------------------------------------------
# Actual (potential) walls -- the certified-necessary refinement
# --------------------------------------------------------------------------
def _bogomolov_integer(ch: ChernChar, d: int) -> Fraction:
    """``c^2/d - 2 r e``: >= 0 is the Bogomolov inequality (rank 0 -> always >= 0)."""
    return Fraction(ch.c * ch.c, d) - 2 * ch.r * ch.e


def _is_integral_chern(ch: ChernChar) -> bool:
    """True iff ``(r, c1=c, ch2=e)`` lies in the Chern-character lattice of P^2,
    i.e. ``c2 = c1^2/2 - ch2`` is an integer (so it is the class of an actual object)."""
    return (Fraction(ch.c * ch.c, 2) - ch.e).denominator == 1


def _heart_ordering_ok(v: ChernChar, w: ChernChar, s0: Fraction, d: int) -> bool:
    """Heart/phase condition at the wall: ``Im Z(w) > 0`` and ``Im Z(v-w) > 0``.

    Along the wall ``Z(w)`` is parallel to ``Z(v)``; for ``w`` to be a genuine
    sub-object of an object of class ``v`` in the tilted heart, both ``w`` and
    the quotient must have strictly positive imaginary central charge there.
    ``Im Z(x) = t (c_x - s_0 r_x d)``, so the sign is that of ``c_x - s_0 r_x d``.
    """
    im_w = w.c - s0 * w.r * d
    im_u = (v.c - w.c) - s0 * (v.r - w.r) * d
    return im_w > 0 and im_u > 0


@dataclass(frozen=True)
class ActualWall:
    """A wall satisfying the rigorous necessary conditions for being an actual wall.

    Same geometry as :class:`Wall`, plus the certified destabilizing class and
    the quotient.
    """

    center: Fraction
    radius_sq: Fraction
    subobject: ChernChar
    quotient: ChernChar
    v: ChernChar
    certificate: Certificate = UNKNOWN_CERTIFICATE

    is_real = True

    @property
    def radius(self) -> float:
        return math.sqrt(float(self.radius_sq))

    @property
    def s_left(self) -> float:
        return float(self.center) - self.radius

    @property
    def s_right(self) -> float:
        return float(self.center) + self.radius

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"ActualWall(center={self.center}, R={self.radius:.4f}, "
            f"sub={self.subobject}, quot={self.quotient}) "
            f"[{self.certificate.rigor.name}]"
        )


def actual_walls(
    v: ChernChar,
    surface: Surface,
    *,
    rank_bound: Optional[int] = None,
    degree_bound: int = 16,
    center_window: Number = 12,
    include_torsion: bool = True,
) -> List[ActualWall]:
    """Enumerate the walls of ``v`` that satisfy the necessary conditions for
    being **actual** destabilizing walls, returning each with its destabilizer.

    A numerical wall ``W(v, w)`` is kept iff ALL of:

    * **rank reduction** -- ``0 <= rank(w) <= rank(v)`` (the first destabilizing
      object has rank at most ``rank(v)``; Coskun-Huizenga survey, sec. 6);
    * **Bogomolov on both pieces** -- ``Delta(w) >= 0`` and ``Delta(v-w) >= 0``;
    * **real semicircle** -- ``radius^2 > 0`` (strict; a zero-radius "wall" is a
      point, not a wall);
    * **heart/phase ordering** -- ``Im Z(w) > 0`` and ``Im Z(v-w) > 0`` on the
      wall, i.e. ``w`` is a genuine sub-object in the tilted heart there.

    These conditions are NECESSARY for an actual wall and -- crucially -- finite
    (the dense numerical-wall set returned by :func:`compute_walls` is cut down
    to the destabilizers that can really occur).  For the Hilbert scheme
    ``P^2[n]`` and the coprime / small-rank cases covered by the ABCH /
    Coskun-Huizenga theorems this is exactly the set of actual walls; e.g. for
    ``P^2[2]`` (``v = (1,0,-2)``) it returns the single wall ``(-5/2, 3/2)``
    realized by ``O(-1)``.  For general ``v`` it is the certified-necessary
    candidate set; use :func:`actual_walls_complete` to check stability of the
    set under widening the search bounds.

    Walls are de-duplicated by ``(center, radius_sq)`` (collecting all
    destabilizers found for each) and sorted by radius, largest first (the
    largest is the Gieseker wall).
    """
    require_faithful_computation(surface)  # G12 guard: torsion-canonical rows refused
    if v.r <= 0:
        raise ValueError("actual_walls expects a positive-rank v")
    d = surface.d
    if v.discriminant(d) < 0:
        raise ValueError("v violates Bogomolov (Delta < 0): not a semistable class")
    if not _is_integral_chern(v):
        raise ValueError(
            "v is not in the Chern-character lattice (c2 not integral): it is not "
            "the class of an actual object, so 'actual walls' are undefined"
        )
    mu_v = v.slope(d)
    S = Q(center_window)
    s_lo, s_hi = mu_v - S, mu_v + S
    r_cap = v.r if rank_bound is None else min(rank_bound, v.r)
    r_lo = 0 if include_torsion else 1

    # Completeness of this finite set is a THEOREM only for the ABCH / Coskun-
    # Huizenga-covered classes on P^2 -- the Hilbert scheme P^2[n] and the coprime
    # cases gcd(rank, degree) = 1.  For a general (non-coprime) P^2 class the set is
    # only the certified-NECESSARY candidate set, so we must NOT stamp it PROVEN:
    # the whole point of the provenance lattice (G5) is honesty about rigor.
    _coprime_on_p2 = (
        surface.is_p2
        and v.c.denominator == 1
        and math.gcd(v.r, v.c.numerator) == 1
    )
    alg_cert = (
        Certificate(
            Rigor.PROVEN,
            ("P^2 Hilbert-scheme / ABCH-Coskun-Huizenga wall structure",
             "Picard rank 1", "coprime class: gcd(rank, degree) = 1"),
            ("arXiv:1203.0316", "arXiv:1202.4587"),
            "actual walls are the ABCH/CH-certified finite set (coprime P^2 class)",
        )
        if _coprime_on_p2 else
        Certificate(
            Rigor.HEURISTIC,
            ("certified-necessary candidate set; completeness proven only for "
             "P^2[n] / coprime / small-rank (ABCH-Coskun-Huizenga)",),
            ("arXiv:1202.4587",),
            "certified-necessary candidates; completeness only empirical outside the "
            "ABCH/CH-covered classes (see actual_walls_complete, doubling certificate)",
        )
    )
    cert = meet(alg_cert, surface.certificate)

    best: dict = {}  # (center, radius_sq) -> ActualWall
    for rp in range(r_lo, r_cap + 1):
        for W_int in range(-degree_bound, degree_bound + 1):
            if W_int == 0:
                continue  # vertical wall
            cp = (Fraction(W_int) + rp * v.c) / v.r
            if cp.denominator != 1:
                continue
            if rp == 0 and cp <= 0:
                continue  # torsion sub needs effective degree > 0
            W_rc = Fraction(W_int)
            # e' placing the center in [s_lo, s_hi]
            e_a = Fraction(s_lo * W_rc + rp * v.e, v.r)
            e_b = Fraction(s_hi * W_rc + rp * v.e, v.r)
            e_lo, e_hi = (e_a, e_b) if e_a <= e_b else (e_b, e_a)
            for ep in _half_integers(e_lo, e_hi):
                w = ChernChar(rp, cp, ep)
                u = ChernChar(v.r - rp, v.c - cp, v.e - ep)
                if not _is_integral_chern(w):
                    continue  # w (hence u) must be the class of an actual object
                if _bogomolov_integer(w, d) < 0:
                    continue
                if _bogomolov_integer(u, d) < 0:
                    continue
                wall = numerical_wall(v, w, d)
                if isinstance(wall, VerticalWall) or wall.radius_sq <= 0:
                    continue
                if not (s_lo <= wall.center <= s_hi):
                    continue
                if not _heart_ordering_ok(v, w, wall.center, d):
                    continue
                key = (wall.center, wall.radius_sq)
                if key not in best:
                    best[key] = ActualWall(
                        wall.center, wall.radius_sq, w, u, v, certificate=cert
                    )
    out = list(best.values())
    out.sort(key=lambda W: (-W.radius_sq, W.center))
    return out


def _abch_covered(v: ChernChar, surface: Surface) -> bool:
    """True where :func:`maciocia_wall_bound` is the PROVEN ABCH outermost wall
    (ABCH arXiv:1203.0316 sec 9): Picard rank 1, ``P^2``, coprime class
    ``gcd(rank, degree) = 1`` (includes every ``P^2[n]`` ideal sheaf
    ``v = (1, 0, -n)``: ``r = 1``).

    This is the SAME coprime-``P^2`` gate already inlined in :func:`actual_walls`
    (``_coprime_on_p2``); factored out so :func:`enumerate_ns_walls` and
    :func:`actual_walls_complete` share exactly one notion of "the family where the
    ``(bog-1)^2/4`` bound is a theorem".  PROVEN completeness holds ONLY here; at
    Picard rank >= 2 the bound is a documented heuristic estimate (E9-M2
    [RESEARCH-light]; the general Thm-3.11 constant is not transcribed).
    """
    return (
        surface.picard_rank == 1
        and surface.is_p2
        and v.c.denominator == 1
        and math.gcd(v.r, v.c.numerator) == 1
    )


def _ceil_sqrt_fraction(x: Fraction) -> int:
    """Smallest non-negative integer ``m`` with ``m*m >= x`` (exact; no float ``sqrt``).

    Used to size a center-search window from an exact ``radius^2``-style bound so the
    window provably reaches a target center distance while everything stays a
    ``Fraction`` (the exact-arithmetic invariant: no float in the core).
    """
    if x <= 0:
        return 0
    m = math.isqrt(x.numerator // x.denominator)   # <= floor(sqrt(x))
    while m * m < x:
        m += 1
    return m


def _covered_center_reach_sq(v: ChernChar, bound: Fraction, d: int) -> Fraction:
    """Exact ``radius^2``-style reach: the squared distance from ``mu_v`` to the CENTER
    of the outermost in-bound wall of ``v`` (no ``sqrt``).

    From the CH wall equation ``rho^2 = (s0 - mu_v)^2 - 2*Delta_v`` (module docstring),
    every real wall with ``radius_sq = rho^2 <= bound`` has its center within
    ``sqrt(bound + 2*Delta_v)`` of ``mu_v``; the outermost (Gieseker) wall, whose
    ``radius_sq == bound``, attains it.  A center-window ``W`` therefore provably
    reaches every in-bound wall CENTER iff ``W^2 >= bound + 2*Delta_v`` -- the sound
    replacement for the unsound ``W^2 >= bound`` (which ignores the ``+2*Delta_v`` and
    truncates the outermost wall while the bound-filter passes vacuously).
    """
    return bound + 2 * v.discriminant(d)


def enumerate_ns_walls(
    v: ChernChar,
    surface: Surface,
    *,
    v_ch1: Optional[Sequence[Number]] = None,
    rank_bound: Optional[int] = None,
    degree_bound: int = 16,
    bidegree_bound: int = 6,
    center_window: Number = 12,
    include_torsion: bool = True,
) -> List[Wall]:
    """The finitely many real Neron-Severi walls of ``v`` inside the Maciocia bound
    ``maciocia_wall_bound(v, surface)`` (E9-M2 / G13).

    Returns a ``List[Wall]`` (uniformly, regardless of Picard rank), de-duplicated by
    ``(center, radius_sq)`` and sorted by ``radius_sq`` descending (the outermost /
    Gieseker wall first).  Every returned wall satisfies
    ``radius_sq <= maciocia_wall_bound(v, surface)``.

    Rigor by Picard rank (E9-M2 [RESEARCH-light], resolved as in E3-M2):

    * **Picard rank 1** (``surface.lattice.rank == 1``).  Delegates to the certified
      :func:`actual_walls` scalar enumerator and filters by the bound.  On the
      ABCH-covered family (``_abch_covered``, e.g. every ``P^2[n]``) the bound is the
      PROVEN outermost wall.  Because :func:`actual_walls` keeps a wall by its CENTER
      and the outermost (Gieseker) wall's center sits at distance
      ``sqrt(bound + 2*Delta_v) > sqrt(bound)`` from ``mu_v``, the delegated
      ``center_window`` is FIRST GROWN (via ``_covered_center_reach_sq`` /
      ``_ceil_sqrt_fraction``) so it provably reaches every in-bound wall center -- else
      the PROVEN outermost wall would be silently dropped while still tagged complete
      (``P^2[n]``, ``n >= 12``).  With that window nothing in-bound is truncated, so the
      returned walls carry :func:`actual_walls`' PROVEN certificate and this is a
      COMPLETE enumeration by the theorem.  For a non-covered rank-1 class the filter
      still holds but completeness is only the certified-necessary/HEURISTIC guarantee
      inherited from :func:`actual_walls`.

    * **Picard rank >= 2** (``surface.lattice.rank >= 2``).  ``v_ch1`` is REQUIRED
      (else :func:`numerical_wall_ns` raises via ``_ns_ch1``).  Sweeps integer NS
      bidegrees ``ch1_w`` over the box ``[-bidegree_bound, bidegree_bound]^rank`` and
      half-integer ``e'`` placing the (gamma-twisted) wall center in
      ``[mu_v - center_window, mu_v + center_window]``, keeping each real wall with
      ``0 < radius_sq <= bound`` whose twisted sub/quotient reps both satisfy Bogomolov
      and the heart/phase ordering.  This is a HEURISTIC bounded search INSIDE the
      ``(bog-1)^2/4`` estimate: it is NOT certified complete, because the general
      Thm-3.11 bounding-semicircle constant is ``[RESEARCH]`` (Maciocia arXiv:1202.4587
      sec 3, not transcribed) and a genuine rank->=2 wall MAY lie OUTSIDE this estimate
      (witness: on ``P^1 x P^1`` the class ``v = (1,0,-2)`` has a real wall at
      ``radius_sq = 26`` while the estimate is only ``9/4`` -- so it is dropped here).
      Every wall returned at Picard rank >= 2 carries a HEURISTIC certificate; none is
      ever tagged PROVEN.

    Everything stays exact ``fractions.Fraction`` (no float; the semicircle-containment
    check used by :func:`actual_walls_complete` -- ``win^2 >= bound + 2*Delta_v`` on the
    center reach, ``left*left >= radius_sq`` on the edges -- avoids ``sqrt``).
    """
    require_faithful_computation(surface)          # G12 guard: torsion-canonical refused
    if v.r <= 0:
        raise ValueError("enumerate_ns_walls expects a positive-rank v")
    bound = maciocia_wall_bound(v, surface)

    # -- Picard rank 1: the certified scalar enumerator, filtered by the bound -----
    if surface.lattice.rank == 1:
        win = Q(center_window)
        if _abch_covered(v, surface):
            # SOUNDNESS (E9-M2): the delegated actual_walls filters candidates by their
            # CENTER, and the outermost (Gieseker) wall of the covered family sits at
            # center-distance sqrt(bound + 2*Delta_v) > sqrt(bound) from mu_v.  A window
            # sized only to sqrt(bound) would silently DROP that PROVEN wall while the
            # bound-filter passed vacuously (P^2[n], n >= 12 lost (-(2n+1)/2, (2n-1)^2/4)
            # though tagged Rigor.PROVEN -- an unsound "complete enumeration").  Grow the
            # window so it provably reaches every in-bound wall center before filtering,
            # keeping the enumeration genuinely COMPLETE for the covered family.
            need = _covered_center_reach_sq(v, bound, surface.d)
            win = max(win, Fraction(_ceil_sqrt_fraction(need)))
        aw = actual_walls(
            v, surface,
            rank_bound=rank_bound,
            degree_bound=degree_bound,
            center_window=win,
            include_torsion=include_torsion,
        )
        return [
            Wall(w.center, w.radius_sq, w.subobject, w.v, w.certificate)
            for w in aw
            if w.radius_sq <= bound
        ]

    # -- Picard rank >= 2: bounded NS-bidegree search (HEURISTIC, inside the bound) -
    d, L = surface.d, surface.lattice
    mu_v = v.slope(d)
    S = Q(center_window)
    s_lo, s_hi = mu_v - S, mu_v + S
    r_cap = v.r if rank_bound is None else min(rank_bound, v.r)
    r_lo = 0 if include_torsion else 1

    heuristic_cert = Certificate(
        Rigor.HEURISTIC,
        ("Picard rank >= 2: bounded NS-bidegree search inside the (bog-1)^2/4 estimate",
         "NOT certified complete: Thm 3.11 constant is [RESEARCH]"),
        ("arXiv:1202.4587",),
        "walls inside the heuristic (bog-1)^2/4 estimate; NOT certified complete at "
        "Picard rank >= 2 -- the Thm 3.11 constant is [RESEARCH] (Maciocia "
        "arXiv:1202.4587 sec 3); genuine walls may lie outside this estimate",
    )
    # Fold in the surface's provenance, but only when it actually asserts something:
    # an UNTAGGED (UNKNOWN) surface certificate carries no claim, so it must not drag
    # the deliberate HEURISTIC search stamp down to UNKNOWN (the rho>=2 catalog rows
    # P^1xP^1 / F_n ship an untagged certificate).  A surface WITH a real (<=HEURISTIC)
    # certificate still lowers the tag via meet -- never over-claiming.
    if surface.certificate.rigor == Rigor.UNKNOWN:
        cert = heuristic_cert
    else:
        cert = meet(heuristic_cert, surface.certificate)

    seen = set()
    walls: List[Wall] = []
    for rp in range(r_lo, r_cap + 1):
        for box in product(range(-bidegree_bound, bidegree_bound + 1), repeat=L.rank):
            ch1_w = tuple(Fraction(x) for x in box)
            cp = L.pairing(ch1_w, surface.H)       # c' = <ch1_w, H> (integer-valued)
            W_rc = v.r * cp - rp * v.c
            if W_rc == 0:
                continue                            # equal twisted slopes -> vertical
            if rp == 0 and cp <= 0:
                continue                            # torsion sub needs effective degree > 0
            # The gamma-twist is affine in e' with unit slope: e~_w = e' + K_w.
            # Compute the twist once (e'=0) to recover K_w and the constant e~_v.
            tv, tw0 = _ns_twisted_pair(
                v, ChernChar(rp, cp, Fraction(0)), surface, v_ch1, ch1_w
            )
            K_w = tw0.e                             # twist offset:  e~_w = e' + K_w
            # Twisted center s0 = (v.r*(e'+K_w) - rp*e~_v)/W_rc; solve s0 in [s_lo, s_hi].
            e_a = Fraction(s_lo * W_rc + rp * tv.e, v.r) - K_w
            e_b = Fraction(s_hi * W_rc + rp * tv.e, v.r) - K_w
            e_lo, e_hi = (e_a, e_b) if e_a <= e_b else (e_b, e_a)
            for ep in _half_integers(e_lo, e_hi):
                w = ChernChar(rp, cp, ep)
                tw = ChernChar(rp, w.c, ep + K_w)   # == _ns_twisted_pair(v, w, ...)[1]
                res = numerical_wall(tv, tw, d)     # == numerical_wall_ns(v, w, ...) geometry
                if isinstance(res, VerticalWall):
                    continue
                if not (0 < res.radius_sq <= bound):
                    continue                        # inside the bound (the roadmap's spec)
                if not (s_lo <= res.center <= s_hi):
                    continue
                # Necessary conditions on the twisted reps (t(v-w) = tv - tw).
                quot = ChernChar(tv.r - tw.r, tv.c - tw.c, tv.e - tw.e)
                if not (_bg_ok(tw, d) and _bg_ok(quot, d)):
                    continue
                if not _heart_ordering_ok(tv, tw, res.center, d):
                    continue
                key = (res.center, res.radius_sq)
                if key in seen:
                    continue
                seen.add(key)
                walls.append(Wall(res.center, res.radius_sq, w, v, cert))

    walls.sort(key=lambda W: (-W.radius_sq, W.center))
    return walls


def actual_walls_complete(
    v: ChernChar, surface: Surface, **kwargs
) -> Tuple[List[ActualWall], bool]:
    """Run :func:`actual_walls` and certify the set is complete.

    Returns ``(walls, complete)``.  For the ABCH-covered ``P^2`` family
    (:func:`_abch_covered`) completeness is a THEOREM, but the certificate is granted
    ONLY when the search window provably contains the WHOLE bounding semicircle of the
    enumerated outermost wall -- NOT merely ``center_window^2 >= maciocia_wall_bound``.
    Since :func:`actual_walls` keeps a wall by its CENTER and the outermost (Gieseker)
    wall's center sits at distance ``sqrt(bound + 2*Delta_v) = radius + 1 > radius`` from
    ``mu_v``, the naive ``win^2 >= bound`` check is UNSOUND: it can truncate that wall
    away while ``all(radius_sq <= bound)`` passes vacuously (``P^2[n]``, ``n >= 12`` used
    to return ``complete=True`` with the outermost wall ``(-(2n+1)/2, (2n-1)^2/4)``
    MISSING).  The fast path now fires only when the enumerated outermost wall really is
    the bound wall (``radius_sq == bound``) and both window edges provably cover its
    semicircle (``left*left >= radius_sq`` and ``right*right >= radius_sq``, exact -- no
    ``sqrt``); then by ABCH arXiv:1203.0316 sec 9 + Bertram's nested-wall theorem the set
    is complete and the walls keep their PROVEN :func:`actual_walls` certificate
    (E9-M2 / G7, replacing the empirical doubling for that family).

    For every OTHER class -- and for a covered class whose window does NOT provably cover
    the semicircle -- the fallback is the unchanged empirical doubling certificate:
    ``complete`` is True iff doubling both ``degree_bound`` and ``center_window`` yields
    the same set of semicircles.  Doubling returns the WIDENED set (so the outermost wall
    is included) with ``complete`` reporting stability honestly; it never over-claims for
    an uncovered class (Picard rank >= 2 stays HEURISTIC; the general Thm-3.11 constant is
    [RESEARCH], E9-M2 [RESEARCH-light]).
    """
    deg = kwargs.pop("degree_bound", 16)
    win = kwargs.pop("center_window", 12)
    base = actual_walls(v, surface, degree_bound=deg, center_window=win, **kwargs)
    if _abch_covered(v, surface) and base:
        bound = maciocia_wall_bound(v, surface)
        # SOUND completeness (E9-M2 fix).  ``win^2 >= bound`` only means the window is
        # wider than the outermost RADIUS -- but actual_walls keeps a wall by its
        # CENTER, and the outermost (Gieseker) wall's center sits at distance
        # sqrt(bound + 2*Delta_v) = radius + 1 > radius from mu_v.  So ``win^2 >= bound``
        # can truncate that wall away while ``all(radius_sq <= bound)`` passes VACUOUSLY
        # (P^2[n], n >= 12: it dropped (-(2n+1)/2, (2n-1)^2/4) yet stamped complete=True).
        # Certify complete ONLY when the window provably contains the whole bounding
        # semicircle of the enumerated outermost wall AND that wall really is the bound
        # wall (radius_sq == bound, so it was not truncated on either the center or the
        # degree axis).  Exact -- no sqrt.  Otherwise fall through to empirical doubling,
        # which returns the widened (wall-included) set with complete=False.
        g = base[0]                              # outermost enumerated wall (sorted first)
        mu_v = v.slope(surface.d)
        W = Q(win)
        left = g.center - (mu_v - W)             # window-left reach vs. semicircle left edge
        right = (mu_v + W) - g.center            # window-right reach vs. semicircle right edge
        if (
            g.radius_sq == bound
            and all(w.radius_sq <= bound for w in base)
            and left >= 0 and left * left >= g.radius_sq
            and right >= 0 and right * right >= g.radius_sq
        ):
            return base, True
    wide = actual_walls(v, surface, degree_bound=2 * deg, center_window=2 * Q(win), **kwargs)
    base_keys = {(w.center, w.radius_sq) for w in base}
    wide_keys = {(w.center, w.radius_sq) for w in wide}
    return wide, base_keys == wide_keys


def maciocia_wall_bound(v: ChernChar, surface_or_d) -> Fraction:
    """Radius^2 of the outermost (Gieseker) wall for a ``P^2[n]`` ideal-sheaf class.

    Accepts EITHER a :class:`~bridgeland_stability.varieties.Surface` (E9-M2 overload)
    OR the raw degree ``d = H^2`` (the original signature).  ``Surface`` is unpacked to
    its ``.d`` so every existing ``maciocia_wall_bound(v, 1)`` / ``(v, P2.d)`` call is
    unchanged bit-for-bit (``P2.d == 1``).

    For the ABCH / Coskun-Huizenga-covered rank-1 P^2 ideal-sheaf classes
    ``v = (1, 0, -n)`` (the Hilbert scheme ``P^2[n]``) the outermost wall is the
    Gieseker wall realized by ``O(-1)``:

        center = -(2n+1)/2 ,   radius^2 = ((2n-1)/2)^2

    (ABCH, arXiv:1203.0316 sec 9; already pinned in :func:`actual_walls`).  Writing
    the integer Bogomolov discriminant ``bog = c^2/d - 2 r e`` (= ``2n`` for that
    family), this is exactly ``radius^2 = (bog - 1)^2 / 4`` -- returned here as an
    exact ``Fraction``.  **For that ``P^2[n]`` family** every real wall of ``v`` then
    satisfies ``radius_sq <= maciocia_wall_bound(v, d)`` (the Gieseker wall is the
    outermost one).

    Scope / rigor (this is the [RESEARCH]-gated milestone -- do not over-claim).
    The closed form ``(bog-1)^2/4`` is EXACTLY-DERIVED and PROVEN as the outermost
    wall ONLY for the ABCH-covered rank-1 P^2 classes above (two-way verified: exact
    Fraction recompute + ABCH arXiv:1203.0316 sec 9).  It is NOT a certified bound for
    a general class: for arbitrary ``v`` the value is only a heuristic estimate and a
    real wall MAY exceed it.  The general finiteness / boundedness of the wall set of a
    fixed class is due to Maciocia (arXiv:1202.4587), who proves the walls are bounded
    *in suitable planes of stability conditions* and gives conditions for global
    finiteness, together with Bertram's nested-wall theorem -- but the EXACT general
    bounding constant is [RESEARCH] (not transcribed from arXiv:1202.4587 sec 3), and
    Maciocia's boundedness is conditional, NOT a blanket single-semicircle bound for
    every class [UNVERIFIED at the specific-theorem level].  Consequently
    For the covered ``P^2[n]`` / coprime-``P^2`` family :func:`actual_walls_complete`
    now certifies completeness FROM this bound (window covers the semicircle; E9-M2 /
    G7); outside that family it keeps the empirical doubling certificate.

    Picard rank >= 2 (E9-M2 [RESEARCH-light]).  On a rank->=2 surface ``v`` never
    carries a ``ch1`` vector here, so ``bog = v.bogomolov_discriminant(d)`` is only the
    **H-projected** ``bog = c^2/d - 2 r e`` (the polarization-orthogonal ``ch1^2`` is
    invisible).  The value is therefore only a HEURISTIC estimate at Picard rank >= 2,
    and a genuine full-NS wall MAY exceed it: on ``P^1 x P^1`` the class ``v=(1,0,-2)``
    has a real NS wall at ``radius^2 = 26`` (the E9-M1-pinned ``ns2 = (-11/2, 26)``)
    while ``maciocia_wall_bound(v, P1xP1) = 9/4``.  The general Thm-3.11
    bounding-semicircle constant is ``[RESEARCH]`` (Maciocia arXiv:1202.4587 sec 3, NOT
    transcribed here), so no PROVEN completeness tag is set at Picard rank >= 2 --
    :func:`enumerate_ns_walls` searches only INSIDE this estimate and stays HEURISTIC.
    Pure ``Fraction`` (no float).
    """
    d = surface_or_d.d if isinstance(surface_or_d, Surface) else surface_or_d
    if v.r <= 0:
        raise ValueError("maciocia_wall_bound expects a positive-rank v")
    bog = v.bogomolov_discriminant(d)  # c^2/d - 2 r e (exact Fraction)
    return (bog - 1) ** 2 / 4
