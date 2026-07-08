"""E8-M1 / G12.1: the Neron-Severi lattice type + rank-1 backward-compat shim.

Stdlib-only at import time (fractions, dataclasses, typing) -- preserving the
zero-runtime-dependency invariant.  Introduces the data model the H-numerical
scalar (r, c=ch1.H, e) core already implicitly uses; ChernChar is NOT migrated
here (that is E8-M2).  The rank-1 shim (rank1_shim / shim_ch1) pins the
bit-for-bit encoding <ch1,H>=c, <ch1,ch1>=c^2/d that protects the pinned gate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Sequence, Tuple, Union

Number = Union[int, Fraction]


def _Q(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


@dataclass(frozen=True)
class NSLattice:
    """A Neron-Severi lattice: rank + symmetric integer Gram (intersection form).

    ``pairing(u, v)`` is the symmetric bilinear form  sum_ij u_i G_ij v_j,
    returned as an exact ``Fraction`` (vectors may be rational).  Frozen and
    hashable; the Gram is stored as a tuple-of-tuples of ``int``.
    """

    rank: int
    gram: Tuple[Tuple[int, ...], ...]

    def __post_init__(self) -> None:
        for row in self.gram:                 # fail LOUD on a non-integral Gram
            for x in row:                     # (intersection numbers are integers);
                if int(x) != x:               # never silently truncate int(1/2)->0
                    raise ValueError(
                        f"Gram entries must be integers (intersection numbers); "
                        f"got non-integral {x!r}")
        g = tuple(tuple(int(x) for x in row) for row in self.gram)
        if self.rank < 1:
            raise ValueError("NSLattice rank must be >= 1")
        if len(g) != self.rank or any(len(row) != self.rank for row in g):
            raise ValueError(f"Gram must be {self.rank}x{self.rank}")
        for i in range(self.rank):
            for j in range(i):
                if g[i][j] != g[j][i]:
                    raise ValueError("Gram matrix must be symmetric")
        object.__setattr__(self, "gram", g)

    def pairing(self, u: Sequence[Number], v: Sequence[Number]) -> Fraction:
        if len(u) != self.rank or len(v) != self.rank:
            raise ValueError("vector length must equal lattice rank")
        total = Fraction(0)
        g = self.gram
        for i in range(self.rank):
            ui = _Q(u[i])
            if not ui:
                continue
            row = g[i]
            for j in range(self.rank):
                gij = row[j]
                if gij:
                    total += ui * gij * _Q(v[j])
        return total

    def self_pairing(self, u: Sequence[Number]) -> Fraction:
        return self.pairing(u, u)

    def det(self) -> Fraction:
        """Exact determinant of the Gram (for Hodge-index signature checks)."""
        return _det(tuple(tuple(Fraction(x) for x in row) for row in self.gram))


def _det(m) -> Fraction:
    n = len(m)
    if n == 1:
        return m[0][0]
    if n == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]
    total = Fraction(0)
    for k in range(n):
        minor = tuple(tuple(row[j] for j in range(n) if j != k) for row in m[1:])
        total += ((-1) ** k) * m[0][k] * _det(minor)
    return total


# -- rank-1 backward-compat shim (pins the pinned-test gate) ---------------
RANK1_AMPLE: Tuple[int, ...] = (1,)


@lru_cache(maxsize=None)
def rank1_shim(d: int) -> NSLattice:
    """Rank-1 NS lattice of a Picard-rank-1 surface with H^2 = d.

    Basis = the ample generator H, so <H,H> = 1*d*1 = d.  This is the encoding
    the scalar (r, c=ch1.H, e) model implicitly uses.  Memoized: NSLattice is
    frozen/immutable so the shared instance is safe, and the slope/discriminant/
    wall-minor/Mukai paths call this per-iteration in hot loops (E8-M3).
    """
    return NSLattice(1, ((d,),))


def shim_ch1(c: Number, d: int) -> Tuple[Fraction, ...]:
    """ch1 as an NS-vector for a rank-1 class of H-degree c = ch1.H: ch1=(c/d)H.

    Stores the H-COORDINATE c/d (NOT the naive scalar c): with <H,H>=d this
    yields <ch1,H> = (c/d)*d = c and <ch1,ch1> = (c/d)^2 * d = c^2/d -- exactly
    the terms the shipped discriminant / bogomolov_discriminant use.  The naive
    (c,) encoding would instead give c*d and c^2*d, breaking the gate.
    """
    return (Fraction(c, d),)
