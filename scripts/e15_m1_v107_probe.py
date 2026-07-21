"""E15-M1: the v_107 rigid-factor probe (the open case of the Sec. 11 conjecture).

Computes the generic H_m-Harder-Narasimhan factors of the paper's candidate
``v_107 = (107, 25/107 E + 76/107 F, 5724/11449)`` on ``F_3`` at a
rational chamber midpoint near ``m = 1``, and checks each factor against the
rigid bound ``Delta_i <= (1 - 1/r_i^2)/2`` (CORRECTIONS Sec. 21: every factor
of a rigid sheaf is rigid by prop-mukai, so a factor above the bound REFUTES
the existence of an exceptional bundle of character v_107 -- which, combined
with the Sec. 24 sweep, would verify the Sec. 11 conjecture through rank 130
outright).  A length-one filtration conflicts with the independently empty
stability interval and must fail closed; it is not a refutation.  All factors
passing is honestly INCONCLUSIVE.

DO NOT RE-RUN AS-IS.  The 2026-07 attempt was killed after 67.4 CPU-hours
with zero output; the post-mortem (CORRECTIONS Sec. 21) measured a ~x10^4
per-rank-doubling scaling law (rank-107 extrapolation: 10^15..10^19 s), a
17-frame-deep Fraction-bound recursion, and a large constant-factor penalty
from the rational chamber midpoint's forced ~7.3e5 denominators.  Re-attempting
requires a genuinely different high-rank algorithm, frontier telemetry with
an explicit CPU budget, and small-denominator sampling with post-hoc wall
certification.  This script is kept as the exact statement of the target
computation and its decision rule.

Usage (only after such a replacement):  python scripts/e15_m1_v107_probe.py
"""

import sys
import time
from fractions import Fraction as F

sys.path.insert(0, __file__.rsplit("scripts", 1)[0])

from bridgeland_stability.generic_hn import generic_hn_factors          # noqa: E402
from bridgeland_stability.dlp_hirzebruch import (                        # noqa: E402
    discriminant,
    is_semiexceptional,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle       # noqa: E402
from bridgeland_stability.delta_sharp import surface_with_index          # noqa: E402


def main() -> None:
    r = 107
    c1 = (76, 25)
    nu = (F(76, 107), F(25, 107))
    delta = F(5724, 11449)                        # = (1 - 1/107^2)/2, potentially exceptional
    S0 = surface_with_index(3, F(1))
    ch2 = r * (F(1, 2) * S0.lattice.self_pairing(nu) - delta)
    print("ch2 =", ch2, flush=True)               # -89/2

    g = F(1, 32 * 107 * 107)                      # the E14-M1 chamber gap at anchor 1
    m = 1 + g / 2
    S = surface_with_index(3, m)
    print(f"sample m = 1 + {g/2}", flush=True)

    t0 = time.time()
    facts = generic_hn_factors(r, c1, ch2, S)
    dt = time.time() - t0
    print(f"[{dt:.0f}s] factors:",
          "None (prioritary stack empty)" if facts is None else len(facts), flush=True)
    if facts is None:
        return
    certificate_conflict = len(facts) == 1        # fail closed; see Sec. 26.2
    refuted = False
    for (ri, c1i, ch2i) in facts:
        xi = SurfaceBundle(ri, c1i, ch2i)
        d = discriminant(xi, S)
        bound = F(1, 2) * (1 - F(1, ri * ri))
        ok = d <= bound
        refuted = refuted or not ok
        print(f"  r={ri} c1={tuple(c1i)} ch2={ch2i} Delta={d} bound={bound} "
              f"RIGID-OK={ok} semiexc={is_semiexceptional(xi, S)}", flush=True)
    if certificate_conflict:
        raise AssertionError(
            "certificate conflict: length-one semistability transfers to a "
            "generic irrational polarization, but v_107's stability interval is empty")
    print("VERDICT:", "REFUTED -- no exceptional bundle of character v_107 "
          "(run the adversarial pass before recording)" if refuted
          else "INCONCLUSIVE -- all factors satisfy the rigid bound", flush=True)


if __name__ == "__main__":
    main()
