"""E13-M2 / G18: the prioritary sharp non-emptiness bound ``delta^p_n(nu)`` on ``F_e``.

Coskun-Huizenga (arXiv:1907.06739 Thm 1.2 = Prop 4.15 + Cor 4.17) give an
**explicitly computable, rank-free, purely rational** function ``delta^p_n(nu)`` with:
for a character ``v = (r, nu, Delta)`` on ``F_e`` with ``Delta >= 0``, the prioritary
stack ``P_{F, H_n}(v)`` is nonempty **iff ``Delta >= delta^p_n(nu)``**.  This is a
theorem (an ``iff``), and it is the first genuine *sharpening* of the ``F_e``
non-emptiness program past E13-M1 (which only relocated the existing lower bound to
the del Pezzo base).

Honest scope (do NOT over-claim).
---------------------------------
* ``delta^p_n`` is the sharp bound for the **prioritary** problem -- the stack of
  ``F``- and ``H_n``-prioritary torsion-free sheaves (``Ext^2(V, V(-L)) = 0``), a
  strictly *weaker* condition than ``mu_H``-semistability.
* By the strong Bogomolov inequality (Remark 1.4) it is a **lower bound** for the
  ``mu``-stable Gieseker bound the package already computes
  (:func:`bridgeland_stability.dlp_hirzebruch.dlp_envelope`):
  ``delta^{mu-s}_m(nu) >= delta^p_{ceil(m)+1}(nu) >= 0``.  It is **not itself** the
  Gieseker / semistable-sheaf existence bound -- that is E13-M3.
* ``F_e`` only.  Off a Hirzebruch surface :func:`delta_prioritary` raises
  ``NotImplementedError`` (via :func:`~bridgeland_stability.dlp_hirzebruch.hirzebruch_index`).
* Exact ``Fraction`` throughout.  ``delta^p_n`` is *purely rational* -- unlike the
  ``sqrt``-bearing Gieseker ``delta^{mu-s}`` -- so there is no float anywhere in a
  returned value.  ``math.floor`` / ``math.ceil`` appear ONLY inside the twist/dual
  reduction (exact ``int`` of a ``Fraction``), never as a returned value.  Stdlib-only
  at import time.

The mathematics (arXiv:1907.06739 Sec. 2.4, Sec. 4.1-4.3), in the package basis.
--------------------------------------------------------------------------------
Total slope ``nu = eps * E + phi * F`` on ``F_e`` (``E`` the section, ``E^2 = -e``,
``E.F = 1``, ``F^2 = 0``; ``F`` the fiber).  The package stores ``NS(F_e)`` in the
basis ``(f, s) = (F, E)``, so a package NS-vector ``(v0, v1) = v0 f + v1 s`` has

    eps = v1  (the s = E coefficient),   phi = v0  (the f = F coefficient).

* **Prioritary sheaf** (Def 2.1).  ``V`` torsion-free is ``L``-prioritary if
  ``Ext^2(V, V(-L)) = 0``.  The relevant polarizations are the fiber ``F`` and
  ``H_m = E + (m + e) F`` (so ``H_m . F = 1``); ``P_F(v)`` is irreducible (Walter) and
  nonempty whenever ``Delta >= 0``.
* **Prop 4.15 -- the closed form.**  On the triangle ``T`` with vertices
  ``(eps, phi) = (-1, n-1), (0, 0), (0, -1)`` write
  ``(eps, phi) = lam1 (-1, n-1) + lam2 (0, 0) + lam3 (0, -1)``, i.e.

      lam1 = -eps,   lam3 = -((n-1) eps + phi),   lam2 = 1 - lam1 - lam3.

  The direct sum ``V = O(-E + (n-1) F)^A (+) O^B (+) O(-F)^C`` with
  ``A = r lam1, B = r lam2, C = r lam3`` has rank ``r``, slope ``nu``, and is ``F``-
  and ``H_n``-prioritary.  From ``ch(V)`` one gets
  ``Delta(V) = A/(2 r^2) (B (e + 2n - 2) + C (e + 2n))``, and the ``r``-factors
  cancel (``A = r lam1`` etc.), giving the **rank-free master formula**

      delta^p_n(nu) = max{ 1/2 lam1 (lam2 (e + 2n - 2) + lam3 (e + 2n)), 0 }   on T.

  If ``eps in Z`` then ``delta^p_n(nu) = 0`` (Def 4.11).
* **Reduction to T for arbitrary nu** (Remark 4.13).  ``delta^p_n`` is invariant under
  integer *twists* ``(eps, phi) -> (eps + a, phi + b)`` (``V -> V (x) O(aE + bF)``:
  ``Delta`` and the prioritary condition are twist-invariant) and *duals*
  ``(eps, phi) -> (-eps, -phi)``.  ``T`` and ``-T`` each have area ``1/2`` and, up to an
  integer translation, tile a ``Z^2``-fundamental domain, so every ``nu`` reduces into
  ``T``: twist ``E`` so ``eps in (-1, 0)``, twist ``F`` so ``lam3 in [0, 1)``; if the
  remaining ``lam2 < 0``, dualize once and re-normalize.  This lands in ``T`` (proved in
  ``docs/CORRECTIONS.md`` Sec. 10: the dual maps the ``lam2 < 0`` sub-triangle back into
  ``T`` with ``lam2 = lam1 + lam3 - 1 > 0``), and the value is preserved throughout.
* **Cor 4.17** -- the ``iff`` above; ``Rigor.PROVEN``.
* **Cor 4.18** -- the generic prioritary index.  For ``eps not in Z``,
  ``rho_gen(v) = floor( Delta / ((ceil(eps) - eps)(eps - floor(eps))) - e/2 + 1
  - (ceil(psi) - psi) )`` with the L-Gaeta parameter
  ``psi = phi + 1/2 e (ceil(eps) - eps) - Delta / (1 - (ceil(eps) - eps))`` and
  ``L_0 = ceil(eps) E + ceil(psi) F``; then ``P_{F, H_n}(v) != empty <=> n <= rho_gen``.

References
----------
* Coskun-Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
  arXiv:1907.06739 -- Sec. 2.4 (prioritary sheaves), Sec. 4.1-4.3 (the L-Gaeta
  resolution and the direct-sum construction), Prop 4.15, Cor 4.17, Cor 4.18,
  Remark 1.4, Remark 4.13.  See ``docs/CORRECTIONS.md`` Sec. 10 for the two-way
  exact-``Fraction`` evidence.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Sequence, Tuple, Union

from .dlp_hirzebruch import hirzebruch_index
from .exceptional_surface import SurfaceBundle
from .varieties import Surface

Number = Union[int, Fraction]

__all__ = [
    "delta_prioritary",
    "prioritary_nonempty",
    "generic_prioritary_index",
    "delta_prioritary_bundle",
]


# --------------------------------------------------------------------------
# Coordinate helpers.  Package NS basis is (f, s) = (F, E); the paper's
# (eps, phi) are the (E, F) = (s, f) coefficients of nu.  (Invariant-2 danger
# zone: a swapped mapping silently yields the wrong delta^p; the two-way
# `discriminant` pins in tests/test_prioritary.py fail loudly if swapped.)
# --------------------------------------------------------------------------
def _eps_phi(nu: Sequence[Number]) -> Tuple[Fraction, Fraction]:
    """Map a package ``(f, s)`` NS-vector ``nu`` to the paper's ``(eps, phi)``.

    ``eps = nu[1]`` (the ``s = E`` coefficient), ``phi = nu[0]`` (the ``f = F``
    coefficient).  Coerces to exact ``Fraction``; raises on a non-length-2 vector.
    """
    nu_t = tuple(Fraction(x) for x in nu)
    if len(nu_t) != 2:
        raise ValueError("delta_prioritary needs a length-2 F_e NS-vector nu = (f, s)")
    return nu_t[1], nu_t[0]


def _lambda_coords(eps: Fraction, phi: Fraction, n: int) -> Tuple[Fraction, Fraction, Fraction]:
    """The barycentric coordinates ``(lam1, lam2, lam3)`` of ``(eps, phi)`` in ``T``.

    ``T`` has vertices ``(-1, n-1), (0, 0), (0, -1)``; ``(eps, phi)`` lies in ``T`` iff
    all three are ``>= 0``.
    """
    lam1 = -eps
    lam3 = -((n - 1) * eps + phi)
    lam2 = 1 - lam1 - lam3
    return lam1, lam2, lam3


def _reduce_to_triangle(eps: Fraction, phi: Fraction, n: int) -> Tuple[Fraction, Fraction]:
    """Twist/dual-reduce ``(eps, phi)`` into the triangle ``T`` (Remark 4.13).

    Preserves ``delta^p_n`` exactly (integer twists and the dual are symmetries).
    **Precondition:** ``eps`` is not an integer (the ``eps in Z`` case returns ``0``
    upstream in :func:`delta_prioritary`).  **Postcondition:** the returned
    ``(eps', phi')`` lies in ``T`` -- all three barycentric coordinates ``>= 0``.

    ``math.floor`` / ``math.ceil`` of a ``Fraction`` are exact ``int``; they only choose
    the integer twist amounts, never a returned value.
    """
    # Twist E so that eps in (-1, 0):  eps -> eps - ceil(eps).
    eps = eps - math.ceil(eps)
    # Twist F so that lam3 in [0, 1):  phi -> phi + floor(lam3)  (lam3 -> lam3 - floor(lam3)).
    lam3 = -((n - 1) * eps + phi)
    phi = phi + math.floor(lam3)
    _, lam2, _ = _lambda_coords(eps, phi, n)
    if lam2 < 0:
        # The lam2 < 0 sub-triangle: dualize once and re-normalize (lands in T).
        eps, phi = -eps, -phi
        eps = eps - math.ceil(eps)
        lam3 = -((n - 1) * eps + phi)
        phi = phi + math.floor(lam3)
    return eps, phi


def _delta_p_on_triangle(eps: Fraction, phi: Fraction, n: int, e: int) -> Fraction:
    """The Prop 4.15 master formula ``max{1/2 lam1 (lam2 (e+2n-2) + lam3 (e+2n)), 0}``.

    Valid only on ``T`` (all ``lam_i >= 0``).  A defensive check guards against a
    reduction bug landing outside ``T`` -- it can never fire after
    :func:`_reduce_to_triangle`, but if it did the bare formula would return a wrong
    value silently.
    """
    lam1, lam2, lam3 = _lambda_coords(eps, phi, n)
    if lam1 < 0 or lam2 < 0 or lam3 < 0:
        raise RuntimeError(
            f"internal: (eps, phi) = ({eps}, {phi}) with n={n} is not in the triangle T "
            f"(lam = {lam1}, {lam2}, {lam3}); the twist/dual reduction is broken.")
    val = Fraction(1, 2) * lam1 * (lam2 * (e + 2 * n - 2) + lam3 * (e + 2 * n))
    return val if val > 0 else Fraction(0)


def _L0_and_psi(
    eps: Fraction, phi: Fraction, Delta: Fraction, e: int
) -> Tuple[Tuple[int, int], Fraction]:
    """The L-Gaeta line bundle ``L_0`` and parameter ``psi`` (Def 4.2 / Cor 4.18).

    ``psi = phi + 1/2 e (ceil(eps) - eps) - Delta / (1 - (ceil(eps) - eps))`` and
    ``L_0 = ceil(eps) E + ceil(psi) F``, returned as ``((a0, b0), psi)`` in ``(E, F)``
    coordinates (matching the paper's Figure-2 caption ``(a0, b0)``).  Requires
    ``eps not in Z``.
    """
    if eps.denominator == 1:
        raise ValueError("the L-Gaeta resolution / Cor 4.18 needs eps not in Z")
    ce = math.ceil(eps)
    frac_up = ce - eps                       # ceil(eps) - eps, in (0, 1)
    # Denominator 1 - frac_up is the Prop 4.15 proof form; for eps not in Z it equals
    # Cor 4.18's printed eps - floor(eps) (since ceil(eps) - floor(eps) = 1).
    psi = phi + Fraction(1, 2) * e * frac_up - Fraction(Delta) / (1 - frac_up)
    return (int(ce), int(math.ceil(psi))), psi


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------
def delta_prioritary(nu: Sequence[Number], n: int, surface: Surface) -> Fraction:
    """``delta^p_n(nu)`` -- the sharp prioritary non-emptiness bound (Prop 4.15).

    Parameters
    ----------
    nu : length-2 package NS-vector ``(f, s)`` (``int`` / ``Fraction``)
        The total slope ``nu = c1 / r``.
    n : int
        The polarization index of ``H_n = E + (n + e) F``.
    surface : Surface
        An ``F_e`` (any polarization -- ``surface.H`` is **not** used, mirroring the
        polarization-independence of ``Delta``; only ``n`` carries the polarization
        index and ``e = hirzebruch_index(surface)``).  Off ``F_e`` raises
        ``NotImplementedError``.

    Returns
    -------
    Fraction
        ``delta^p_n(nu) >= 0``, exact.
    """
    e = hirzebruch_index(surface)                    # NotImplementedError off F_e
    if n != int(n):
        raise ValueError("n must be an integer polarization index")
    n = int(n)
    eps, phi = _eps_phi(nu)
    if eps.denominator == 1:                          # Def 4.11: integer section slope
        return Fraction(0)
    eps_r, phi_r = _reduce_to_triangle(eps, phi, n)
    return _delta_p_on_triangle(eps_r, phi_r, n, e)


def prioritary_nonempty(
    r: int, nu: Sequence[Number], Delta: Number, n: int, surface: Surface
) -> bool:
    """Cor 4.17: ``P_{F, H_n}(v) != empty  <=>  Delta >= delta^p_n(nu)`` (``Delta >= 0``).

    A PROVEN ``iff`` verdict for the *prioritary* stack.  Validates that ``v`` is a
    genuine integral character of ``F_e``: ``c1 = r nu in NS`` and
    ``c2 = 1/2 <c1,c1> - ch2 in Z`` with ``ch2 = r (1/2 <nu,nu> - Delta)`` -- a
    non-integral class is out of the domain of Cor 4.17 (``v in K(F_e)``) and raises
    ``ValueError``.  Raises ``ValueError`` on ``Delta < 0`` (Cor 4.17 assumes the
    Bogomolov floor) and ``NotImplementedError`` off ``F_e``.
    """
    # E13 re-audit R3: never int()-truncate the rank -- Fraction(3,2) used to be
    # silently answered as r = 1, a Cor 4.17 verdict for a different character.
    r_frac = Fraction(r)
    if r_frac.denominator != 1 or r_frac < 1:
        raise ValueError(f"rank must be a positive integer; got {r!r}")
    r = int(r_frac)
    Delta = Fraction(Delta)
    if Delta < 0:
        raise ValueError("Cor 4.17 assumes the Bogomolov floor Delta >= 0")
    bound = delta_prioritary(nu, n, surface)         # validates F_e + nu length
    lat = surface.lattice
    nu_f = tuple(Fraction(x) for x in nu)
    c1 = tuple(r * x for x in nu_f)
    if any(c.denominator != 1 for c in c1):
        raise ValueError(f"c1 = r*nu = {c1} is not an integral NS class")
    ch2 = r * (Fraction(1, 2) * lat.self_pairing(nu_f) - Delta)
    c2 = Fraction(1, 2) * lat.self_pairing(c1) - ch2
    if c2.denominator != 1:
        raise ValueError(
            f"c2 = 1/2<c1,c1> - ch2 = {c2} is not an integer: (r, nu, Delta) is not a "
            "genuine K(F_e) character (Cor 4.17's domain)")
    return Delta >= bound


def generic_prioritary_index(nu: Sequence[Number], Delta: Number, surface: Surface) -> int:
    """Cor 4.18: the generic prioritary index ``rho_gen(v)`` (requires ``Delta >= 0``
    and ``eps not in Z``).

    ``P_{F, H_n}(v) != empty  <=>  n <= rho_gen(v)``.  Returns the exact integer

        rho_gen = floor( Delta / ((ceil(eps) - eps)(eps - floor(eps))) - e/2 + 1
                         - (ceil(psi) - psi) )

    with ``psi`` the L-Gaeta parameter (:func:`_L0_and_psi`).  Reproduces the paper's
    Example 4.9 / Figure 2 worked example: ``nu = 1/2 E + 1/3 F``, ``Delta = 11/10``,
    ``e = 1`` gives ``psi = -97/60``, ``L_0 = (1, -1)``, ``rho_gen = 4``.  (``psi`` and the
    ``rho_gen = 4`` conclusion are stated in Example 4.9's text; the Figure 2 caption
    itself carries only ``nu``, ``Delta``, ``e``, and ``(a0, b0) = (1, -1)``.)
    """
    e = hirzebruch_index(surface)                    # NotImplementedError off F_e
    Delta = Fraction(Delta)
    if Delta < 0:
        # E13 re-audit R3: same domain as Cor 4.17 -- the prioritary stack P_F(v) is
        # nonempty iff Delta >= 0 (Walter), so rho_gen is defined only on the Bogomolov
        # floor; below it the formula returned a meaningless index (e.g. -4 at Delta=-1).
        raise ValueError("Cor 4.18 rho_gen assumes the Bogomolov floor Delta >= 0")
    eps, phi = _eps_phi(nu)
    if eps.denominator == 1:
        raise ValueError("Cor 4.18 rho_gen needs eps not in Z (integer section slope)")
    _, psi = _L0_and_psi(eps, phi, Delta, e)
    frac_up = math.ceil(eps) - eps                   # ceil(eps) - eps
    frac_low = eps - math.floor(eps)                 # eps - floor(eps)
    val = (Delta / (frac_up * frac_low) - Fraction(e, 2) + 1
           - (math.ceil(psi) - psi))
    return int(math.floor(val))


def general_betti(r: int, c1: Sequence[Number], ch2: Number,
                  surface: Surface) -> Tuple[int, int, int]:
    """``(h^0, h^1, h^2)`` of the GENERAL ``F``-prioritary sheaf of the integral
    character ``(r, c1, ch2)`` with ``Delta >= 0`` — arXiv:1907.06739 Thm "BN"
    (= CoskunHuizengaBN Thm 3.1), implemented as stated (E15-M1e):

    * ``nu.F = -1``: ``(0, -chi, 0)``;
    * ``nu.F > -1``: ``h^2 = 0``; if ``nu.E >= -1`` at most one nonzero group
      (``h^0 = max(chi, 0)``, ``h^1 = max(-chi, 0)``); else
      ``h^0(V) = h^0(V(-E))`` inductively (each ``-E`` twist lowers ``nu.F``
      by 1, so the recursion terminates on every ``F_e``);
    * ``nu.F < -1``: ``h^0 = 0`` and ``h^2`` by Serre duality from the dual
      character (rank >= 2 for local freeness of the general sheaf).

    Package coordinates: ``nu.F`` is the ``s``-component of ``nu``; ``nu.E`` is
    ``f``-component − e·(s-component).  Cross-validated in the tests against
    ``rho_gen`` through the E15-M1e Gaeta dimension inequality.
    """
    e = hirzebruch_index(surface)
    lat = surface.lattice
    r = int(r)
    c1 = tuple(int(x) for x in c1)
    ch2 = Fraction(ch2)

    def chi_of(w) -> int:
        rr, cc, hh = w
        # chi(O, w) via RR: r*chi(O) - <c1, K>/2 + ch2  (chi(O) = 1)
        K = (-(e + 2), -2)
        val = rr - Fraction(lat.pairing(cc, K), 2) + hh
        if val.denominator != 1:
            raise ValueError(f"non-integral chi for {w!r}")
        return int(val)

    def tw(w, D):
        rr, cc, hh = w
        return (rr, (cc[0] + rr * D[0], cc[1] + rr * D[1]),
                hh + lat.pairing(cc, D) + rr * Fraction(lat.self_pairing(D), 2))

    def dual_serre(w):
        rr, cc, hh = w
        return tw((rr, (-cc[0], -cc[1]), hh), (-(e + 2), -2))

    def go(w, depth=0):
        if depth > 10000:
            raise AssertionError("thm-BN recursion failed to terminate")
        rr, cc, hh = w
        nuF = Fraction(cc[1], rr)
        nuE = Fraction(cc[0], rr) - e * Fraction(cc[1], rr)
        x = chi_of(w)
        if nuF == -1:
            return (0, -x, 0)
        if nuF > -1:
            if nuE >= -1:
                return (max(x, 0), max(-x, 0), 0)
            h0 = go(tw(w, (0, -1)), depth + 1)[0]
            return (h0, h0 - x, 0)
        h2 = go(dual_serre(w), depth + 1)[0]
        return (0, h2 - x, h2)

    return go((r, c1, ch2))


def delta_prioritary_bundle(
    nu: Sequence[Number], n: int, r: int, surface: Surface
) -> SurfaceBundle:
    """The Prop 4.15 direct-sum witness ``V`` at slope ``nu`` (``nu`` must be in ``T``).

    ``V = O(-E + (n-1) F)^A (+) O^B (+) O(-F)^C`` with ``A = r lam1, B = r lam2,
    C = r lam3``.  Requires the multiplicities ``A, B, C = r lam_i`` to be non-negative
    integers -- i.e. ``nu`` already lies in ``T`` and ``r lam_i in Z`` (raise
    ``ValueError`` otherwise).  For an in-``T`` slope,
    ``discriminant(delta_prioritary_bundle(nu, n, r, surface), surface)`` equals
    ``delta_prioritary(nu, n, surface)`` whenever the latter is ``> 0`` (the two-way
    ``Delta`` check).
    """
    e = hirzebruch_index(surface)                    # NotImplementedError off F_e
    if n != int(n):
        raise ValueError("n must be an integer polarization index")
    n = int(n)
    # E13 re-audit R3: never int()-truncate the rank -- Fraction(3,2) used to yield
    # the Prop 4.15 witness for r = 1, a different character than the caller supplied.
    r_frac = Fraction(r)
    if r_frac.denominator != 1 or r_frac < 1:
        raise ValueError(f"rank must be a positive integer; got {r!r}")
    r = int(r_frac)
    eps, phi = _eps_phi(nu)
    lam1, lam2, lam3 = _lambda_coords(eps, phi, n)
    A, B, C = r * lam1, r * lam2, r * lam3
    if any(x.denominator != 1 for x in (A, B, C)) or A < 0 or B < 0 or C < 0:
        raise ValueError(
            f"A, B, C = r*lambda = ({A}, {B}, {C}) must be non-negative integers "
            "(nu must be in the triangle T with r*lambda_i integral)")
    A, B, C = int(A), int(B), int(C)
    lat = surface.lattice
    # Summands in the (f, s) = (F, E) basis:
    #   O(-E + (n-1) F): c1 = (n-1) f + (-1) s = (n-1, -1)
    #   O:               c1 = (0, 0)
    #   O(-F):           c1 = (-1) f + 0 s = (-1, 0)
    L1 = (Fraction(n - 1), Fraction(-1))
    L3 = (Fraction(-1), Fraction(0))
    c1 = (A * (n - 1) + C * (-1), A * (-1))          # = (A(n-1) - C, -A)
    ch2 = (A * (Fraction(1, 2) * lat.self_pairing(L1))
           + C * (Fraction(1, 2) * lat.self_pairing(L3)))
    return SurfaceBundle(A + B + C, c1, ch2)          # rank A+B+C = r
