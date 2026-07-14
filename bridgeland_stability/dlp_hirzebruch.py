"""E11-M6 / G18: the polarization-dependent Drezet-Le Potier envelope on ``F_e``.

This module retires the deferred **G18 remainder**: the sharp, polarization-dependent
non-emptiness bound ``delta_H`` off ``P^2``.  It ships two things that the E11-M3
"thin slice" could only delegate:

1. **The full-NS discriminant** ``Delta = 1/2 <nu,nu> - ch2/r`` with ``nu = c1/r`` an
   NS-vector (:func:`discriminant`).  This is the discriminant of the primary sources
   -- verbatim in arXiv:1907.06739 Sec. 2.1 -- and it is *polarization independent*.
   The package's older :func:`bridgeland_stability.nonemptiness_rational.discriminant_H`
   is the H-PROJECTED scalar ``1/2 mu^2 - ch2/(r d)``; the two agree (up to the exact
   factor ``d``) only when ``c1`` is proportional to ``H`` -- in particular on every
   Picard-rank-1 surface, and on ``P^2`` (``d = 1``) they are equal.  See
   ``docs/CORRECTIONS.md`` (C7).

2. **The Drezet-Le Potier envelope** ``DLP_H(nu)`` built from the ``mu_H``-stable
   exceptional bundles of ``F_e``, which *is* the sharp bound on an anticanonically
   polarized del Pezzo Hirzebruch surface (``e = 0, 1``).

Mathematical content (all citations to Coskun-Huizenga, "Existence of semistable
sheaves on Hirzebruch surfaces", arXiv:1907.06739)
--------------------------------------------------------------------------------
Write ``F_e = P(O_{P^1} + O_{P^1}(e))`` with ``E`` a section of self-intersection
``-e`` and ``F`` the fiber class (``E^2 = -e``, ``E.F = 1``, ``F^2 = 0``).  This
package stores the NS lattice in the basis ``(f, s)`` with Gram ``[[0,1],[1,-n]]``,
i.e. ``f = F`` and ``s = E`` and ``n = e``.

* **Riemann-Roch** (Sec. 2.1).  ``P(nu) = chi(O_X) + 1/2 (nu^2 - nu.K_X)`` and
  ``chi(V,W) = r_V r_W (P(nu_W - nu_V) - Delta_V - Delta_W)``.
* **Potentially exceptional** characters (Sec. 5, Lemma "excFacts"): ``chi(v,v) = 1``,
  equivalently ``Delta = 1/2 - 1/(2 r^2)``.  Integrality (``c_2 in Z``) forces, in the
  ``(E,F)`` basis with ``c1 = a E + b F``, that ``gcd(r,a) = 1``, that ``r`` is odd
  when ``e`` is even, and the congruence
  ``2ab = a^2 e + a e r - r^2 - 1  (mod 2r)`` -- so ``b`` is determined **mod r** by
  ``(r, a)``.  That congruence is what makes the enumeration below finite and cheap.
* **The per-bundle bound** ``DLP_{H,V}`` (Sec. 5.4, "The Drezet-Le Potier surface").
  For a ``mu_H``-stable exceptional bundle ``V`` and ``w = nu - nu(V)``,

      DLP_{H,V}(nu) = P( w)       - Delta_V     if  1/2 K.H <= w.H < 0
                    = P(-w)       - Delta_V     if  0 < w.H <= -1/2 K.H
                    = max{P(w),P(-w)} - Delta_V if  w.H = 0

  on the strip ``|w.H| <= -1/2 K.H``.  The first two branches are exactly
  ``chi(V,W) <= 0`` resp. ``chi(W,V) <= 0``, which hold because ``Hom`` and ``Ext^2``
  vanish by stability + Serre duality.
* **The envelope** ``DLP_H(nu) = sup_V DLP_{H,V}(nu)`` over ``mu_H``-stable exceptional
  ``V`` in the strip, and the rank-truncated ``DLP_H^{<r}(nu)``, whose supremum is a
  **maximum** (bounded sublevel sets) and hence exactly computable.
* **Stability of exceptional bundles** (Cor. "DLPExceptional").  A potentially
  exceptional ``v = (r, nu, Delta)`` carries a ``mu_H``-stable exceptional bundle
  **iff** ``Delta >= DLP_H^{<r}(nu)``.  This is the induction on rank that makes the
  whole thing computable, and it is what :func:`is_stable_exceptional` implements.
* **Sharpness** (Cor. "deltaDLP" / "deltaDLPe").  For every ample ``H`` one has
  ``delta_H^{mu-s}(nu) >= DLP_H(nu)``; and if ``e = 0`` or ``1`` and ``H`` is the
  **anticanonical** polarization then ``delta_{1-e/2}^{mu-s}(nu) = DLP_{-K}(nu)``
  exactly.  Also ``DLP_{-K_{F_e}}(nu) >= 1/2`` (Cor. "K1/2") -- the ``F_e`` analogue
  of the ``1/2`` floor of the ``P^2`` curve.

Certified truncation (re-derived here; the paper only says "computable as a limit")
----------------------------------------------------------------------------------
On an **anticanonical** polarization (``K = -lambda H``, ``lambda > 0``) write
``w = alpha H + beta T`` with ``T.H = 0``.  Then ``<T,T> = -d`` and ``w.K = -lambda
alpha d``, so on the strip ``|alpha| <= lambda/2`` the value collapses to

    DLP_{-K,V}(nu) = 1 + 1/2 d (alpha^2 - beta^2) - 1/2 lambda d |alpha| - Delta_V
                   = 1 - Delta_V + 1/2 d |alpha| (|alpha| - lambda) - 1/2 d beta^2
                  <= 1 - Delta_V   =  1/2 + 1/(2 r_V^2),

with equality iff ``w = 0``.  Hence every exceptional bundle of rank ``> R``
contributes at most ``1/2 + 1/(2 (R+1)^2)``, so if the rank-``<= R`` maximum already
attains that value the truncation is **exact**.  :func:`dlp_envelope` reports this as
:attr:`DLPEnvelope.exact`; it is what upgrades a lower bound to the sharp ``delta_H``.

Scope and honesty
-----------------
* Only ``F_e`` (``e >= 0``) with a **strictly ample** ``H`` is handled.  The
  ``hirzebruch(n)`` factory polarization ``H = n f + s`` is nef-and-big but *not*
  ample and is refused; use
  :func:`bridgeland_stability.nonemptiness_rational.hirzebruch_with_polarization`.
* ``DLP_H(nu)`` is **sharp** (= ``delta_H^{mu-s}``) only for ``e in {0,1}`` with the
  anticanonical ``H``.  For any other ample ``H`` the returned value is a *certified
  lower bound* -- the paper exhibits characters where the inequality is strict and
  ``delta_H`` is computed by Kronecker modules, not exceptional bundles.
* Exact arithmetic throughout: every value is a ``Fraction``.  The only non-``Fraction``
  arithmetic is ``math.isqrt`` / ``floor`` / ``ceil`` used to bound the **search
  region** (never a mathematical value), and integer ``floor``/``ceil`` of a
  ``Fraction`` is exact.  Stdlib-only at import time.

References
----------
* Coskun-Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
  arXiv:1907.06739 (Adv. Math. 381 (2021)) -- everything above.
* Drezet-Le Potier, Ann. Sci. ENS 18 (1985) -- the ``P^2`` closed-form curve, which
  this generalizes (see :mod:`bridgeland_stability.dlp`).
* Rudakov, Gorodentsev, Kuleshov-Orlov -- exceptional bundles on del Pezzo surfaces.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import ceil, floor, gcd, isqrt
from typing import Optional, Sequence, Tuple, Union

from .exceptional_surface import SurfaceBundle, canonical_class
from .varieties import Surface

__all__ = [
    "DLPEnvelope",
    "discriminant",
    "total_slope",
    "hilbert_P",
    "exceptional_discriminant",
    "exceptional_ch2",
    "hirzebruch_index",
    "is_ample",
    "is_anticanonical",
    "is_del_pezzo_anticanonical",
    "is_potentially_exceptional",
    "is_stable_exceptional",
    "is_semiexceptional",
    "dlp_bundle",
    "dlp_restricted",
    "dlp_envelope",
    "emptiness_bound",
    "DEFAULT_RANK_MAX",
]

Number = Union[int, Fraction]
NSVector = Tuple[Fraction, ...]

#: Default rank cutoff for the exceptional-bundle enumeration.  Ranks above this are
#: not enumerated; the envelope then reports ``exact=False`` unless the certified
#: stopping rule (see module docstring) already fired at a lower rank.
DEFAULT_RANK_MAX = 12


# --------------------------------------------------------------------------
# Basic invariants (valid on ANY surface carrying a genuine NS lattice)
# --------------------------------------------------------------------------
def total_slope(xi: SurfaceBundle) -> NSVector:
    """The total slope ``nu = c1 / r`` as an exact NS-vector (arXiv:1907.06739 Sec. 2.1)."""
    if xi.r <= 0:
        raise ValueError("total_slope needs positive rank")
    return tuple(c / xi.r for c in xi.c1)


def discriminant(xi: SurfaceBundle, surface: Surface) -> Fraction:
    """The **full-NS** Coskun-Huizenga discriminant ``Delta = 1/2 <nu,nu> - ch2/r``.

    This is the discriminant of the primary sources, and unlike the H-projected
    surrogate :func:`bridgeland_stability.nonemptiness_rational.discriminant_H` it does
    **not** depend on the polarization.  The two are related by ``Delta = d * Delta_H``
    exactly when ``c1`` is proportional to ``H`` (in particular on every Picard-rank-1
    surface); on ``P^2`` (``d = 1``) they are equal.
    """
    nu = total_slope(xi)
    return Fraction(1, 2) * surface.lattice.self_pairing(nu) - Fraction(xi.ch2, xi.r)


def hilbert_P(nu: Sequence[Number], surface: Surface) -> Fraction:
    """``P(nu) = chi(O_X) + 1/2 (<nu,nu> - <nu, K_X>)`` (arXiv:1907.06739 Sec. 2.1).

    Reduces on ``P^2`` to :func:`bridgeland_stability.exceptional.P`, i.e.
    ``P(m) = (m^2 + 3m + 2)/2``.
    """
    lat = surface.lattice
    K = canonical_class(surface)
    nu = tuple(Fraction(x) for x in nu)
    return surface.chi_O + Fraction(1, 2) * (lat.self_pairing(nu) - lat.pairing(nu, K))


def exceptional_discriminant(r: int) -> Fraction:
    """``Delta_V = 1/2 - 1/(2 r^2)`` -- forced by ``chi(V,V) = 1`` when ``chi(O_X) = 1``."""
    if r < 1:
        raise ValueError("rank must be >= 1")
    return Fraction(1, 2) - Fraction(1, 2 * r * r)


# --------------------------------------------------------------------------
# F_e recognition + polarization predicates
# --------------------------------------------------------------------------
#: Family tags that are never a Hirzebruch surface, refused outright even when the
#: caller hands them an F_e-shaped Gram matrix (E13 re-audit R1: a projective K3 with
#: NS = U carries the F_0 Gram, and rebased the F_e Gram for every even e).
_NON_HIRZEBRUCH_KINDS = frozenset({"P2", "K3", "abelian", "enriques", "bielliptic"})


def hirzebruch_index(surface: Surface) -> int:
    """The index ``e`` of ``F_e``, AUTHENTICATED -- never read off the Gram matrix alone.

    The NS Gram ``[[0,1],[1,-e]]`` supplies the candidate ``e``, but the Gram does not
    identify the surface family: a projective K3 with ``NS(X) = U`` carries the ``F_0``
    Gram (and, in a rebased ``(f, s - k f)`` basis, the ``F_{2k}`` Gram), yet has
    ``K_X = 0`` and ``chi(O_X) = 2``.  Dispatching such a surface into the CH theory
    minted a false PROVEN_NONEMPTY -- ``ch = (5,(-3,1),-3)`` on that K3 has Mukai
    self-pairing ``v^2 = -26 < -2`` with ``gcd(r, c1.H) = 1``, so its moduli space is
    empty (E13 re-audit, defect R1).  So beyond the Gram shape this now requires the
    invariants a genuine ``F_e`` must carry:

    * ``e >= 0``;
    * ``surface.K == (-(e+2), -2)`` -- ``K_{F_e} = -(e+2) f - 2 s`` in the ``(f, s)``
      basis (the Lemma 11.3(3) normalization; ``K = 0`` fails this for every ``e``);
    * ``surface.chi_O == 1`` (``chi(O_{F_e}) = 1``; a K3 has ``chi(O) = 2``);
    * ``surface.kind`` is not a declared non-Hirzebruch family
      (:data:`_NON_HIRZEBRUCH_KINDS`).

    ``(Gram, K, chi(O))`` jointly pin the family: a smooth projective surface with this
    rank-2 NS lattice, ``K^2 = 8`` and ``chi(O) = 1`` is a minimal rational surface,
    hence a Hirzebruch surface.  ``P^1 x P^1`` is ``F_0``.  Raises
    ``NotImplementedError`` for anything else -- the CH theory is stated for Hirzebruch
    surfaces and this module never silently mis-models another surface.
    """
    lat = surface.ns_lattice
    if lat is None or lat.rank != 2:
        raise NotImplementedError(
            f"{surface.name}: the CH Drezet-Le Potier envelope needs the rank-2 F_e NS "
            "lattice (Gram [[0,1],[1,-e]]); got no explicit ns_lattice.")
    g = lat.gram
    if not (g[0][0] == 0 and g[0][1] == 1 and g[1][0] == 1):
        raise NotImplementedError(
            f"{surface.name}: expected the F_e fiber/section Gram [[0,1],[1,-e]]; got {g!r}.")
    e = -g[1][1]
    if e < 0:
        raise NotImplementedError(
            f"{surface.name}: Gram [[0,1],[1,{-e}]] has e = {e} < 0; not an F_e lattice "
            "in the fiber/section normalization.")
    # E13 re-audit R1: authenticate the surface FAMILY, not just the lattice shape.
    if surface.kind in _NON_HIRZEBRUCH_KINDS:
        raise NotImplementedError(
            f"{surface.name}: kind={surface.kind!r} is not a Hirzebruch surface; the CH "
            "F_e theory does not apply even though the NS Gram is F_e-shaped.")
    if tuple(surface.K) != (-(e + 2), -2):
        raise NotImplementedError(
            f"{surface.name}: K = {tuple(surface.K)} != (-(e+2), -2) = {(-(e + 2), -2)}, "
            f"the canonical class of F_{e} in the (f, s) basis: not a Hirzebruch surface "
            "(a K3/abelian surface with NS = U has K = 0 and MUST be refused here; "
            "treating it as F_e forges non-emptiness proofs).")
    if surface.chi_O != 1:
        raise NotImplementedError(
            f"{surface.name}: chi(O) = {surface.chi_O} != 1 = chi(O_{{F_e}}): not a "
            "Hirzebruch surface.")
    return e


def is_ample(surface: Surface) -> bool:
    """Nakai-Moishezon on ``F_e``: ``H = a f + b s`` is ample iff ``b > 0`` and ``a > e b``.

    The :func:`bridgeland_stability.varieties.hirzebruch` factory polarization
    ``H = e f + s`` sits on the nef boundary (``a - e b = 0``) and is **not** ample.
    """
    e = hirzebruch_index(surface)
    a, b = surface.H
    return b > 0 and a - e * b > 0


def is_anticanonical(surface: Surface) -> bool:
    """``True`` iff ``H`` lies on the ray of ``-K_{F_e} = (e+2) f + 2 s``.

    ``DLP_{H,V}`` and ``mu_H``-stability depend on ``H`` only through its ray, so this
    proportionality test (not equality) is the right notion.
    """
    e = hirzebruch_index(surface)
    a, b = surface.H
    return 2 * a == (e + 2) * b


def is_del_pezzo_anticanonical(surface: Surface) -> bool:
    """``True`` iff ``e in {0,1}`` and ``H`` is anticanonical -- the case where
    ``DLP_{-K}(nu)`` is *exactly* the sharp bound ``delta_H^{mu-s}(nu)``
    (arXiv:1907.06739 Cor. "deltaDLP")."""
    return hirzebruch_index(surface) in (0, 1) and is_anticanonical(surface)


def _require_ample_hirzebruch(surface: Surface) -> int:
    e = hirzebruch_index(surface)
    if not is_ample(surface):
        a, b = surface.H
        raise ValueError(
            f"{surface.name}: H={surface.H} is not strictly ample on F_{e} "
            f"(Nakai needs b>0 and a>e*b; the hirzebruch({e}) factory H=({e},1) is only "
            "nef-and-big).  Use hirzebruch_with_polarization(e, H) with an ample H.")
    return e


# --------------------------------------------------------------------------
# Potentially exceptional characters
# --------------------------------------------------------------------------
def exceptional_ch2(r: int, c1: Sequence[int], surface: Surface) -> Fraction:
    """The ``ch2`` forced on a potentially exceptional character of rank ``r``.

    ``chi(v,v) = 1`` pins ``Delta = 1/2 - 1/(2r^2)``, and inverting
    ``Delta = 1/2 <nu,nu> - ch2/r`` gives ``ch2 = <c1,c1>/(2r) - r Delta_V``.  The
    character ``(r, c1, ch2)`` is integral iff ``c2 = <c1,c1>/2 - ch2`` is an integer --
    which is what :func:`is_potentially_exceptional` tests.
    """
    c1c1 = surface.lattice.self_pairing(tuple(Fraction(x) for x in c1))
    return c1c1 / (2 * r) - r * exceptional_discriminant(r)


_ch2_of_exceptional = exceptional_ch2          # internal alias (kept for readability below)


def is_potentially_exceptional(r: int, c1: Sequence[int], surface: Surface) -> bool:
    """``True`` iff ``(r, c1)`` carries an integral character with ``chi(v,v) = 1``.

    ``chi(v,v) = r^2 (chi(O_X) - 2 Delta) = 1`` forces ``Delta = 1/2 - 1/(2 r^2)``
    (arXiv:1907.06739 Lemma "excFacts" (1)); the character is then integral iff
    ``c_2 = <c1,c1>/2 - ch2`` is an integer.  This is the basis-free form of the
    congruence ``2ab = a^2 e + aer - r^2 - 1 (mod 2r)`` used by the enumerator
    :func:`_b_residue` (their agreement is a pinned test).
    """
    if r < 1:
        return False
    if any(Fraction(x).denominator != 1 for x in c1):
        return False
    c1c1 = surface.lattice.self_pairing(tuple(Fraction(x) for x in c1))
    c2 = Fraction(1, 2) * c1c1 - _ch2_of_exceptional(r, c1, surface)
    return c2.denominator == 1


def _b_residue(e: int, r: int, a: int) -> Optional[int]:
    """The unique ``b mod r`` with ``2ab = a^2 e + a e r - r^2 - 1 (mod 2r)``.

    arXiv:1907.06739 Lemma "excFacts" (2).  Here ``a`` is the ``E = s`` coefficient and
    ``b`` the ``F = f`` coefficient of ``c1``.  Returns ``None`` when no exceptional
    character of rank ``r`` has this ``a`` (``gcd(r,a) != 1``, or ``e`` even and ``r``
    even).  This turns the enumeration from a 2-D box into a 1-D arithmetic progression.
    """
    if r == 1:
        return 0                                   # every line bundle is exceptional
    if gcd(r, a) != 1:
        return None
    if e % 2 == 0 and r % 2 == 0:
        return None
    R = a * a * e + a * e * r - r * r - 1
    if R % 2:                                      # defensive: R is always even here
        return None
    return ((R // 2) * pow(a, -1, r)) % r


# --------------------------------------------------------------------------
# DLP_{H,V} and the envelope
# --------------------------------------------------------------------------
def dlp_bundle(
    nu: Sequence[Number], r_V: int, c1_V: Sequence[int], surface: Surface
) -> Optional[Fraction]:
    """``DLP_{H,V}(nu)`` for the exceptional bundle ``V = (r_V, c1_V)``.

    Returns ``None`` when ``nu`` lies **outside** ``V``'s strip
    ``|(nu - nu_V).H| <= -1/2 K.H``, where the paper proves nothing.  See the module
    docstring for the three branches.
    """
    lat = surface.lattice
    K = canonical_class(surface)
    H = tuple(Fraction(x) for x in surface.H)
    kappa = -Fraction(1, 2) * lat.pairing(K, H)
    nu = tuple(Fraction(x) for x in nu)
    nuV = tuple(Fraction(c, r_V) for c in c1_V)
    w = tuple(x - y for x, y in zip(nu, nuV))
    wH = lat.pairing(w, H)
    if wH > kappa or wH < -kappa:
        return None
    dV = exceptional_discriminant(r_V)
    mw = tuple(-x for x in w)
    if wH < 0:
        return hilbert_P(w, surface) - dV
    if wH > 0:
        return hilbert_P(mw, surface) - dV
    return max(hilbert_P(w, surface), hilbert_P(mw, surface)) - dV


def _isqrt_ceil(x: Fraction) -> Fraction:
    """An exact ``Fraction`` over-estimate of ``sqrt(x)`` (``x >= 0``), via integer isqrt."""
    p, q = x.numerator, x.denominator
    return Fraction(isqrt(p * q) + 1, q)


@lru_cache(maxsize=None)
def _halfwidths(e: int, H: Tuple[int, int], gram: Tuple[Tuple[int, ...], ...],
                K: Tuple[Fraction, ...], chi_O: int) -> Tuple[Fraction, Fraction]:
    """Coordinate half-widths of the region of ``w = nu - nu_V`` where ``DLP_{H,V}(nu) >= 0``.

    Derivation.  Put ``T = (e h1 - h0, h1)``, the integral vector with ``T.H = 0``; one
    computes ``<T,T> = -d``.  Writing ``w = alpha H + beta T`` gives ``w.H = alpha d``
    (so ``|alpha| <= kappa/d``) and ``w^2 = d(alpha^2 - beta^2)``.  Since
    ``DLP_{H,V}(nu) <= 1 + w^2/2 + |w.K|/2`` and ``Delta_V >= 0``, the requirement
    ``DLP >= 0`` forces ``(d/2) beta^2 - (|T.K|/2) |beta| - A <= 0`` with
    ``A = 1 + 3 kappa^2/(2d)``, hence a bound on ``|beta|``.  Only a **search bound** --
    never a mathematical value -- so the integer ``isqrt`` over-estimate is harmless.
    """
    from .nslattice import NSLattice

    lat = NSLattice(2, gram)
    Hf = tuple(Fraction(x) for x in H)
    d = lat.self_pairing(Hf)
    kappa = -Fraction(1, 2) * lat.pairing(K, Hf)
    T = (Fraction(e * H[1] - H[0]), Fraction(H[1]))
    TK = lat.pairing(T, K)
    A = 1 + Fraction(3, 2) * kappa * kappa / d
    disc = TK * TK + 8 * d * A
    B_beta = (abs(TK) + _isqrt_ceil(disc)) / (2 * d)
    B_alpha = kappa / d
    W0 = B_alpha * abs(Hf[0]) + B_beta * abs(T[0])
    W1 = B_alpha * abs(Hf[1]) + B_beta * abs(T[1])
    return W0, W1


def _surface_key(surface: Surface) -> Tuple:
    """Hashable identity of everything the envelope depends on: ``(e, H, Gram, K, chi_O)``."""
    e = hirzebruch_index(surface)
    return (e, tuple(int(x) for x in surface.H), surface.lattice.gram,
            canonical_class(surface), surface.chi_O)


def _strip_candidates(key: Tuple, nu: Tuple[Fraction, ...], rank_max: int):
    """Yield every potentially-exceptional ``(rho, c1_V)`` with ``rho <= rank_max`` whose
    ``DLP_{H,V}(nu)`` could be ``>= 0``.  Finite by :func:`_halfwidths`; ``b`` walks the
    arithmetic progression fixed by :func:`_b_residue` rather than a 2-D box."""
    e = key[0]
    W0, W1 = _halfwidths(*key)
    for rho in range(1, rank_max + 1):
        if e % 2 == 0 and rho % 2 == 0:               # no even-rank exceptionals when e even
            continue
        a_lo, a_hi = floor(rho * (nu[1] - W1)), ceil(rho * (nu[1] + W1))
        b_lo, b_hi = floor(rho * (nu[0] - W0)), ceil(rho * (nu[0] + W0))
        for a in range(a_lo, a_hi + 1):               # a = E = s coefficient
            b0 = _b_residue(e, rho, a)
            if b0 is None:
                continue
            b = b_lo + ((b0 - b_lo) % rho)            # b = F = f coefficient
            while b <= b_hi:
                yield rho, (b, a)                     # package (f, s) coordinate order
                b += rho


@lru_cache(maxsize=None)
def _stable_exc_cached(key: Tuple, r: int, c1: Tuple[int, int]) -> bool:
    surface = _surface_from_key(key)
    if not is_potentially_exceptional(r, c1, surface):
        return False
    if r == 1:
        return True                                    # line bundles: no stability condition
    nu = tuple(Fraction(c, r) for c in c1)
    bound = _dlp_restricted_cached(key, nu, r)
    return True if bound is None else exceptional_discriminant(r) >= bound


@lru_cache(maxsize=None)
def _dlp_restricted_cached(key: Tuple, nu: Tuple[Fraction, ...], r: int) -> Optional[Fraction]:
    surface = _surface_from_key(key)
    best: Optional[Fraction] = None
    for rho, c1V in _strip_candidates(key, nu, r - 1):
        val = dlp_bundle(nu, rho, c1V, surface)
        if val is None or val < 0:                    # Bogomolov floor: negatives never bind
            continue
        if best is not None and val <= best:          # cheap value test before the recursion
            continue
        if not _stable_exc_cached(key, rho, c1V):
            continue
        best = val
    return best


@lru_cache(maxsize=None)
def _surface_from_key(key: Tuple) -> Surface:
    """Rebuild the minimal ``Surface`` the envelope needs from its hashable key."""
    from .nslattice import NSLattice

    e, H, gram, K, chi_O = key
    lat = NSLattice(2, gram)
    d = lat.self_pairing(tuple(Fraction(x) for x in H))
    # K is canonical_class(original) = (-(e+2), -2) (A8); reuse it directly so the
    # reconstructed surface's derived K.H matches the original, not the -2 placeholder.
    return Surface(name=f"F_{e} (H={H})", d=int(d), K=K, chi_O=chi_O,
                   picard_rank=2, kind="hirzebruch", H=H, ns_lattice=lat)


def is_stable_exceptional(r: int, c1: Sequence[int], surface: Surface) -> bool:
    """``True`` iff ``(r, c1)`` is the character of a ``mu_H``-stable exceptional bundle.

    arXiv:1907.06739 Cor. "DLPExceptional": a potentially exceptional character carries a
    ``mu_H``-stable exceptional bundle **iff** ``Delta >= DLP_H^{<r}(nu)``.  Since
    ``DLP_H^{<r}`` only involves strictly smaller ranks, this is a terminating induction
    (rank 1 -- line bundles -- is unconditional).

    Reproduces, bit-for-bit, Tables 1 and 2 of the paper (stability intervals for every
    exceptional bundle of rank <= 19 on ``F_0`` and <= 20 on ``F_1``); see
    ``tests/test_dlp_hirzebruch.py``.
    """
    _require_ample_hirzebruch(surface)
    return _stable_exc_cached(_surface_key(surface), int(r), tuple(int(x) for x in c1))


def is_semiexceptional(xi: SurfaceBundle, surface: Surface) -> bool:
    """``True`` iff ``xi`` is ``V^{+m}`` for a ``mu_H``-stable exceptional bundle ``V``.

    A semiexceptional character has ``r = m r_V``, ``c1 = m c1_V`` and ``Delta = Delta_V``;
    ``V^{+m}`` is Gieseker-semistable, so its moduli space is non-empty even though it
    sits strictly below the envelope.  This is the ``F_e`` analogue of the ``P^2``
    exceptional disjunct of the Drezet-Le Potier theorem, and it must be honoured or
    ``moduli_nonempty`` would report a PROVEN "empty" for a character that exists.
    """
    _require_ample_hirzebruch(surface)
    r = xi.r
    if r < 1 or any(c.denominator != 1 for c in xi.c1):
        return False
    delta = discriminant(xi, surface)
    for m in range(1, r + 1):
        if r % m:
            continue
        rV = r // m
        if any(int(c) % m for c in xi.c1):
            continue
        if delta != exceptional_discriminant(rV):
            continue
        c1V = tuple(int(c) // m for c in xi.c1)
        if is_stable_exceptional(rV, c1V, surface):
            return True
    return False


@dataclass(frozen=True)
class DLPEnvelope:
    """The rank-truncated Drezet-Le Potier envelope ``DLP_H^{<= R}(nu)``.

    Attributes
    ----------
    value : Fraction
        ``max(0, DLP_H^{<R+1}(nu))`` -- always a **certified lower bound** for
        ``DLP_H(nu)`` and hence (Cor. "deltaDLPe") for ``delta_H^{mu-s}(nu)``.
    rank_max : int
        The rank cutoff ``R`` actually used.
    exact : bool
        ``True`` iff the truncation provably equals ``DLP_H(nu)``.  Certified only on an
        anticanonical del Pezzo (``e in {0,1}``), via the stopping rule
        ``value >= 1/2 + 1/(2(R+1)^2)`` derived in the module docstring.
    sharp : bool
        ``True`` iff ``DLP_H(nu)`` equals the sharp bound ``delta_H^{mu-s}(nu)``,
        i.e. iff the surface is an anticanonically polarized del Pezzo Hirzebruch
        surface (Cor. "deltaDLP").  Independent of ``exact``.
    witness : tuple | None
        ``(r_V, c1_V)`` of an exceptional bundle attaining ``value``, if any.
    """

    value: Fraction
    rank_max: int
    exact: bool
    sharp: bool
    witness: Optional[Tuple[int, Tuple[int, int]]] = None

    @property
    def certified_sharp(self) -> bool:
        """``True`` iff ``value`` *is* ``delta_H^{mu-s}(nu)`` (up to the ``1/2`` floor)."""
        return self.exact and self.sharp


def dlp_restricted(
    nu: Sequence[Number], surface: Surface, r: int
) -> Optional[Fraction]:
    """``DLP_H^{<r}(nu)``: the exact max of ``DLP_{H,V}(nu)`` over ``mu_H``-stable
    exceptional ``V`` of rank ``< r`` in the strip.  ``None`` iff no such ``V`` exists.

    Values ``< 0`` are dropped (they never bind: Bogomolov gives ``Delta >= 0``), so a
    non-``None`` result is always ``>= 0``.
    """
    _require_ample_hirzebruch(surface)
    key = _surface_key(surface)
    return _dlp_restricted_cached(key, tuple(Fraction(x) for x in nu), int(r))


def dlp_envelope(
    nu: Sequence[Number], surface: Surface, rank_max: int = DEFAULT_RANK_MAX
) -> DLPEnvelope:
    """The truncated envelope ``DLP_H^{<= rank_max}(nu)`` with its exactness certificate.

    On an anticanonically polarized del Pezzo (``e in {0,1}``) the returned ``value`` is
    additionally floored at ``1/2`` by Cor. "K1/2" (``DLP_{-K}(nu) >= 1/2``), matching the
    ``1/2`` clamp of the ``P^2`` curve in :func:`bridgeland_stability.dlp.delta`.
    """
    e = _require_ample_hirzebruch(surface)
    if rank_max < 1:
        raise ValueError("rank_max must be >= 1")
    key = _surface_key(surface)
    nu_t = tuple(Fraction(x) for x in nu)
    raw = _dlp_restricted_cached(key, nu_t, rank_max + 1)

    sharp = is_del_pezzo_anticanonical(surface)
    witness = _envelope_witness(key, nu_t, rank_max, raw, surface)

    if sharp:
        value = max(Fraction(1, 2), raw if raw is not None else Fraction(0))
        # every V of rank > rank_max contributes at most 1/2 + 1/(2 (R+1)^2)
        exact = raw is not None and raw >= Fraction(1, 2) + Fraction(1, 2 * (rank_max + 1) ** 2)
    else:
        value = max(Fraction(0), raw if raw is not None else Fraction(0))
        exact = False                    # no certified stopping rule off the anticanonical ray
    return DLPEnvelope(value, rank_max, exact, sharp, witness)


def _envelope_witness(key, nu, rank_max, raw, surface):
    if raw is None:
        return None
    for rho, c1V in _strip_candidates(key, nu, rank_max):
        if dlp_bundle(nu, rho, c1V, surface) == raw and _stable_exc_cached(key, rho, c1V):
            return (rho, c1V)
    return None


def emptiness_bound(
    xi: SurfaceBundle, surface: Surface, rank_max: int = DEFAULT_RANK_MAX
) -> Fraction:
    """A lower bound on ``Delta`` valid for **every Gieseker-semistable** sheaf of slope
    ``nu(xi)``.  If ``discriminant(xi, surface) < emptiness_bound(xi, surface)`` then
    ``M_{X,H}(xi)`` is **provably empty**.

    This is deliberately *weaker* than :func:`dlp_envelope`, and the difference is not
    cosmetic.  Only two of the three branches of ``DLP_{H,V}`` are theorems about
    semistable sheaves (arXiv:1907.06739 Sec. 5.4):

    * ``0 < |(nu - nu_V).H| <= -1/2 K.H`` -- ``Hom(V,W) = 0`` and ``Ext^2(V,W) = 0`` by
      stability and Serre duality, so ``chi(V,W) <= 0`` and Riemann-Roch gives
      ``Delta >= DLP_{H,V}(nu)``.  Included.
    * ``nu = nu_V`` and ``Delta != Delta_V`` -- one of ``Hom(V,W)``, ``Hom(W,V)`` vanishes
      by Gieseker semistability (the reduced Hilbert polynomials differ), and
      Riemann-Roch gives ``Delta >= 1 - Delta_V = 1/2 + 1/(2 r_V^2)``.  Included.
    * ``(nu - nu_V).H = 0`` with ``nu != nu_V`` -- the paper calls this branch
      "somewhat arbitrary" and it is **not** a valid bound.  **Excluded.**

    The exclusion is load-bearing.  On ``F_0`` with the (non-generic) anticanonical
    polarization, ``O(F_1) + O(F_2)`` is Gieseker-semistable with ``Delta = 1/4``, yet
    ``DLP_{-K}(1/2, 1/2) = 3/4`` -- computed precisely by the excluded branch at
    ``V = O(F_1)``.  Using the full envelope as an emptiness test would declare an
    existing sheaf empty; this function returns ``1/4`` there and correctly concludes
    nothing.  (The paper states this counterexample immediately before Cor. "K1/2".)
    """
    e = _require_ample_hirzebruch(surface)
    key = _surface_key(surface)
    lat = surface.lattice
    H = tuple(Fraction(x) for x in surface.H)
    nu = total_slope(xi)
    delta = discriminant(xi, surface)
    best = Fraction(0)                                # Bogomolov: Delta >= 0 always
    for rho, c1V in _strip_candidates(key, nu, rank_max):
        nuV = tuple(Fraction(c, rho) for c in c1V)
        w = tuple(x - y for x, y in zip(nu, nuV))
        wH = lat.pairing(w, H)
        if wH == 0:
            if any(x != 0 for x in w):                # excluded branch -- proves nothing
                continue
            if delta == exceptional_discriminant(rho):
                continue                              # xi is semiexceptional at V: no bound
            cand = 1 - exceptional_discriminant(rho)
        else:
            cand = dlp_bundle(nu, rho, c1V, surface)
            if cand is None:                          # outside V's strip
                continue
        if cand <= best:
            continue
        if not _stable_exc_cached(key, rho, c1V):
            continue
        best = cand
    return best
