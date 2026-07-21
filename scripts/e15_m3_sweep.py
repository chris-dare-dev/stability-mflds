"""E15-M3: the Conjecture A falsification sweep.

Scans generic HN filtrations over a bounded character x polarization grid on
F_0 / F_1 (``c1`` modulo rank, ``0 <= c2 <= 2r``, and eight right-offset
polarization samples per surface),
classifying every computed factor's semiexceptionality:

* a length >= 3 filtration with >= 2 non-semiexceptional factors is an
  IMMEDIATE counterexample to the Sec. 1.5 conjecture -- printed loudly;
* every length-2 both-non-semiexceptional pair is fed to the block-Kronecker
  decomposition search; pairs with NO witness in the bounded family are
  printed as RANKED CANDIDATES (search-bounded, not counterexamples).

Usage: python scripts/e15_m3_sweep.py [--rank-max R] [--verbose]
"""

import argparse
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("scripts", 1)[0])

from bridgeland_stability.varieties import P1xP1, hirzebruch          # noqa: E402
from bridgeland_stability.block_kronecker import (                     # noqa: E402
    block_decomposition,
    classify_generic_filtration,
)
from bridgeland_stability.delta_sharp import surface_with_index        # noqa: E402
from bridgeland_stability.dlp_hirzebruch import discriminant           # noqa: E402
from bridgeland_stability.exceptional_surface import SurfaceBundle     # noqa: E402

# Right chamber-offsets from the listed anchors, including the Sec. 18 walls.
# These straddle the anticanonical index 1 on F_0; on F_1 every sample is to
# the right of the anticanonical index 1/2.
ANCHORS = (Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2),
           Fraction(5, 2), Fraction(25, 9), Fraction(12, 7), Fraction(3))
EPS = Fraction(1, 1000)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank-max", type=int, default=10)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    n_class = 0
    lengths = {}
    n_pairs = 0
    n_pairs_decomposed = 0
    candidates = []
    violations = []

    for ep in (0, 1):
        for anchor in ANCHORS:
            m = anchor + EPS
            S = surface_with_index(ep, m)
            for r in range(2, args.rank_max + 1):
                for x in range(r):
                    for y in range(r):
                        c1 = (x, y)
                        # ch2 lattice: c2 in Z around the interesting band
                        c1c1 = S.lattice.self_pairing(c1)
                        for c2 in range(0, 2 * r + 1):
                            ch2 = Fraction(c1c1, 2) - c2
                            xi = SurfaceBundle(r, c1, ch2)
                            d = discriminant(xi, S)
                            if d < 0 or d > 2:
                                continue
                            p = classify_generic_filtration(r, c1, ch2, S)
                            n_class += 1
                            if p.length is None:
                                continue
                            lengths[p.length] = lengths.get(p.length, 0) + 1
                            if p.conjecture_a_violation:
                                violations.append((ep, str(m), r, c1, str(ch2), p))
                                print(f"VIOLATION e={ep} m={m} v=({r},{c1},{ch2}) "
                                      f"pattern={p.semiexceptional}", flush=True)
                            elif p.length == 2 and p.non_semiexceptional_count == 2:
                                n_pairs += 1
                                w = block_decomposition(p.factors[0], p.factors[1], S)
                                if w is None:
                                    candidates.append(
                                        (ep, str(m), r, c1, str(ch2), p.factors))
                                    print(f"CANDIDATE (no block witness) e={ep} m={m} "
                                          f"v=({r},{c1},{ch2}) factors={p.factors}",
                                          flush=True)
                                else:
                                    n_pairs_decomposed += 1
                                    if args.verbose:
                                        print(f"  pair decomposed: e={ep} m={m} "
                                              f"v=({r},{c1},{ch2}) l={w.l} "
                                              f"twist={w.twist}", flush=True)
            print(f"[e={ep} anchor={anchor}: cumulative {n_class} classified, "
                  f"{time.time() - t0:.0f}s]", flush=True)

    print(f"\nSWEEP COMPLETE rank_max={args.rank_max}  [{time.time() - t0:.0f}s]",
          flush=True)
    print(f"classified: {n_class}; length histogram: {dict(sorted(lengths.items()))}",
          flush=True)
    print(f"length-2 both-non-semiexceptional pairs: {n_pairs} "
          f"({n_pairs_decomposed} block-decomposed)", flush=True)
    print(f"LENGTH>=3 VIOLATIONS: {len(violations)}", flush=True)
    print(f"UNDECOMPOSED CANDIDATES: {len(candidates)}", flush=True)


if __name__ == "__main__":
    main()
