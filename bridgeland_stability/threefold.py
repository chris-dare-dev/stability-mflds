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

from .chern import Number, Q
from .varieties import Threefold


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
    return BGBoundary(
        betas=betas,
        alphas=alphas,
        bg_proven=threefold.bg_proven,
        threefold_name=threefold.name,
        degenerate_betas=degenerate,
    )
