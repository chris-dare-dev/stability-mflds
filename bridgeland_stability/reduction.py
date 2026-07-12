"""E13-M1 / G18: the Coskun-Huizenga F_e -> F_{e-2} reduction map pi (arXiv:1907.06739 Sec.11.1).

``pi`` is the unimodular NS isometry  ``M = [[1,-1],[0,1]]``  acting on ``c1`` in the ``(f,s)``
basis, fixing ``r`` and ``ch2``:

    pi(r, (x, y), ch2) = (r, (x - y, y), ch2),        M^T G_{e-2} M = G_e.

Here ``(x, y)`` is ``c1 = x f + y s`` (``x`` the fiber ``f = F`` coefficient = the paper's ``b``,
``y`` the section ``s = E`` coefficient = the paper's ``a``); the paper writes the same map in the
``(E, F)`` basis as ``pi(r, aE + bF, d) = (r, aE' + (b-a)F', d)``.  ``M in SL_2(Z)`` (``det = 1``,
unimodular) and ``M^T G_{e-2} M = G_e`` -- exactly re-derived below -- so ``pi`` is an ISOMETRY of
Neron-Severi lattices ``NS(F_e) -> NS(F_{e-2})`` that fixes ``r`` and ``ch2``.

    M^T G_{e-2} M = [[1,0],[-1,1]] . [[0,1],[1,-(e-2)]] . [[1,-1],[0,1]] = [[0,1],[1,-e]] = G_e.

Every Lemma 11.3 property of the paper follows from those three facts (unimodular + isometric +
``r, ch2`` fixed), because each downstream invariant is built only from the pairing ``<.,.>``, ``r``,
``ch2``, ``K_X``, and ``chi(O_X)`` -- all of which ``pi`` preserves:

* (1) intersection pairing:  ``<pi u, pi v>_{e-2} = <u, v>_e``  (isometry).
* (2) discriminant:  ``Delta = 1/2 <nu,nu> - ch2/r``, ``nu = c1/r``; isometry + fixed ``ch2/r``.
* (3) ``pi(K_{F_e}) = K_{F_{e-2}} = (-e, -2)`` and ``pi(ch O_{F_e}) = ch O_{F_{e-2}}``
  (``K_{F_e} = (-(e+2), -2)`` maps to ``(-e, -2)``; ``ch O = (1,(0,0),0)`` is fixed).
* (4) ``chi(v)``, ``chi(v, w)`` (RR Euler form; ``chi(O_X) = 1`` on every ``F_e``), and -- ``M``
  unimodular -- integral -> integral, primitive -> primitive (``gcd(x-y, y) = gcd(x, y)``).
* (5) polarization / Hilbert:  ``A_m = -1/2 K_{F_e} + m F`` maps to ``A'_m = -1/2 K_{F_{e-2}} + m F'``,
  so ``mu_{A_m}`` and ``hilbert_P(nu)`` are preserved.
* (6) the paper's named three-term direct-sum character of Sec.11.1(6): ``pi`` is additive on Chern
  characters (``r, c1, ch2`` all add under direct sum and ``M`` is linear), so it commutes with
  ``(+)``.

Iterating ``e -> e-2`` (:func:`reduce_to_del_pezzo`) lands on the del Pezzo cases ``F_0`` (even ``e``)
or ``F_1`` (odd ``e``), where E11-M6's sharp ``DLP_{-K}`` theory lives.

Honest scope (E13-M1).  Reducing a strictly **ample** ``F_e`` (``e >= 2``) polarization can never
land on an anticanonical (sharp) ``F_0`` / ``F_1`` ray: ``pi(H)`` is proportional to ``-K_{F_{e-2}}``
iff ``H`` is proportional to ``-K_{F_e}`` (``pi`` is a bijection with ``pi(-K_{F_e}) = -K_{F_{e-2}}``),
and ``-K_{F_e}`` is not ample for ``e >= 2``.  So the reduced envelope is a **certified lower bound
equal to the direct one**, not a sharp value; the sharp ``delta_H`` off the ``-K`` ray is E13-M2/M3
(see ``docs/CORRECTIONS.md`` Sec. 9, open question O2).

Exact ``Fraction`` throughout; stdlib-only at import (only intra-package imports, all stdlib-only).

References
----------
* Coskun-Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
  arXiv:1907.06739 Sec. 11.1 and Lemma 11.3 -- the reduction map and its exactness properties.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence, Tuple, Union

from .dlp_hirzebruch import hirzebruch_index, is_ample
from .exceptional_surface import SurfaceBundle
from .nonemptiness_rational import hirzebruch_with_polarization
from .varieties import Surface

Number = Union[int, Fraction]

__all__ = ["pi_c1", "reduce", "reduce_to_del_pezzo", "REDUCTION_MATRIX"]

#: The unimodular SL_2(Z) matrix ``M`` acting on the column ``(f-coeff, s-coeff)`` of ``c1``.
#: ``pi(c1) = M c1``; ``M^T G_{e-2} M = G_e`` makes ``pi`` an NS isometry ``F_e -> F_{e-2}``.
REDUCTION_MATRIX: Tuple[Tuple[int, int], Tuple[int, int]] = ((1, -1), (0, 1))


def pi_c1(c1: Sequence[Number]) -> Tuple[Number, ...]:
    """``pi`` on a length-2 ``F_e`` NS-vector in the ``(f, s)`` basis: ``(x, y) -> (x - y, y)``.

    Linear and exact: ``int`` -> ``int``, ``Fraction`` -> ``Fraction`` (no float).  This is the
    matrix action ``M . (x, y)`` with ``M = REDUCTION_MATRIX``.
    """
    if len(c1) != 2:
        raise ValueError("pi_c1 needs a length-2 F_e NS-vector (f, s)")
    x, y = c1
    return (x - y, y)


def reduce(xi: SurfaceBundle, surface: Surface) -> Tuple[SurfaceBundle, Surface]:
    """Reduce a character on an ample ``F_e`` (``e >= 2``) to ``F_{e-2}`` via ``pi``.

    arXiv:1907.06739 Sec. 11.1.  Returns ``(xi', surface')`` with
    ``xi' = (r, pi(c1), ch2)`` and ``surface'`` the ample ``F_{e-2}`` carrying the transported
    polarization ``pi(H) = A'_m`` (Lemma 11.3(5)).  ``r`` and ``ch2`` are untouched; exact
    ``Fraction``, no float.

    Raises
    ------
    ValueError
        if ``e < 2`` (``F_e`` is already a del Pezzo -- nothing to reduce), or if ``H`` is not
        strictly ample (the transported theory needs strict ampleness; the nef-and-big factory
        polarization is refused -- use
        :func:`bridgeland_stability.nonemptiness_rational.hirzebruch_with_polarization`).
    NotImplementedError
        (via :func:`bridgeland_stability.dlp_hirzebruch.hirzebruch_index`) for a non-``F_e`` surface.

    Note.  ``pi(H)`` is automatically ample on ``F_{e-2}``: writing ``H = (a, b)`` with ``b > 0`` and
    ``a - e b > 0`` (Nakai on ``F_e``), ``pi(H) = (a - b, b)`` satisfies
    ``(a - b) - (e-2) b = (a - e b) + b > 0``, so the constructed ``F_{e-2}`` is never refused.
    """
    e = hirzebruch_index(surface)               # NotImplementedError off F_e
    if e < 2:
        raise ValueError(
            f"reduce needs e >= 2 (F_{e} is already a del Pezzo; nothing to reduce)")
    if not is_ample(surface):
        raise ValueError(
            f"{surface.name}: reduce needs a strictly ample H (use "
            "hirzebruch_with_polarization); the nef-and-big factory H is refused")
    if len(xi.c1) != 2:
        raise ValueError("xi.c1 must be a length-2 F_e NS-vector")
    xi_red = SurfaceBundle(xi.r, pi_c1(xi.c1), xi.ch2)
    H_red = tuple(int(v) for v in pi_c1(surface.H))          # pi(H): integral, ample on F_{e-2}
    return xi_red, hirzebruch_with_polarization(e - 2, H_red)


def reduce_to_del_pezzo(xi: SurfaceBundle, surface: Surface) -> Tuple[SurfaceBundle, Surface]:
    """Iterate ``pi`` until the base is ``F_0`` (even ``e``) or ``F_1`` (odd ``e``).

    A no-op for ``e in {0, 1}`` (the character/surface are returned unchanged).  Property (2)
    telescopes: :func:`bridgeland_stability.dlp_hirzebruch.discriminant` of the result equals that of
    the input exactly.  Raises like :func:`reduce` for a non-ample ``H`` (each intermediate step is
    validated), or ``NotImplementedError`` for a non-``F_e`` surface.
    """
    e = hirzebruch_index(surface)               # NotImplementedError off F_e
    while e >= 2:
        xi, surface = reduce(xi, surface)
        e -= 2
    return xi, surface
