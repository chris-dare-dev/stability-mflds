"""E15-M2: extend the Sec. 11 conjecture's verified range by sweep.

Methodology (the paper's own, arXiv:1907.06739 Sec. 11, made executable):

1. enumerate the EXCEPTIONAL characters on the del Pezzo bases F_0 / F_1 up to
   rank R, one twist-normalized representative each (0 <= c1_f, c1_s < r;
   potentially exceptional + mu_{-K}-stable-exceptional, i.e. the bundle
   exists -- cor-delPezzoExceptional);
2. compute each stability interval by the E14-M3 memoized rank induction;
3. lift k = 1, 2 reduction steps to F_{e'+2k} (pi^{-1}: (x, y) -> (x + k y, y));
   the transported interval {t > 0 : t + k in I_W} is EMPTY iff hi(I_W) <= k --
   an empty lift is a would-be counterexample IF an exceptional bundle of the
   lifted character exists on F_{e'+2k};
4. dispatch each empty lift with the E15-M1 battery's prioritary condition:
   rho_gen = 1 refuses existence (the conjecture holds for that character);
   rho_gen >= 2 goes on the SURVIVOR ledger for M1-style attack.

Scope caveat (recorded in CORRECTIONS Sec. 22): the swept family is the
pi-lifts of del Pezzo exceptional characters.  A hypothetical exceptional
bundle on F_e whose pi-reduction is not an exceptional character lies outside
the family -- same implicit scope as the paper's Example computations.

Usage:  python scripts/e15_m2_sweep.py [--rank-max R]      (default R = 130)
"""

import argparse
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("scripts", 1)[0])

from bridgeland_stability.varieties import P1xP1, hirzebruch  # noqa: E402
from bridgeland_stability.dlp_hirzebruch import (              # noqa: E402
    is_potentially_exceptional,
    is_stable_exceptional,
)
from bridgeland_stability.stability_interval import stability_interval  # noqa: E402
from bridgeland_stability.prioritary import generic_prioritary_index   # noqa: E402
from bridgeland_stability.delta_sharp import surface_with_index        # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank-max", type=int, default=130)
    args = ap.parse_args()
    R = args.rank_max

    t_start = time.time()
    n_exc = 0
    n_empty = 0
    survivors = []
    dispatched = 0

    for ep in (0, 1):
        base = P1xP1 if ep == 0 else hirzebruch(1)
        S_antic = surface_with_index(ep, Fraction(2 - ep, 2))
        for r in range(2, R + 1):
            t_rank = time.time()
            found = 0
            for x in range(r):
                for y in range(r):
                    if not is_potentially_exceptional(r, (x, y), base):
                        continue
                    if not is_stable_exceptional(r, (x, y), S_antic):
                        continue
                    n_exc += 1
                    found += 1
                    iv = stability_interval(r, (x, y), base)
                    if iv.empty or iv.hi is None:
                        continue
                    for k in (1, 2):
                        if iv.hi <= k:
                            n_empty += 1
                            e = ep + 2 * k
                            lift = (x + k * y, y)
                            nu = (Fraction(lift[0], r), Fraction(lift[1], r))
                            delta = Fraction(1, 2) * (1 - Fraction(1, r * r))
                            rho = generic_prioritary_index(nu, delta, hirzebruch(e))
                            if rho >= 2:
                                survivors.append((e, r, lift, iv.hi, rho))
                                print(f"SURVIVOR e={e} r={r} c1={lift} "
                                      f"delPezzo_hi={iv.hi} rho_gen={rho}", flush=True)
                            else:
                                dispatched += 1
            if (found and time.time() - t_rank > 5) or r % 10 == 0:
                print(f"[F_{ep} rank {r}: {found} exceptional chars, "
                      f"{time.time() - t_rank:.0f}s, total {time.time() - t_start:.0f}s]",
                      flush=True)

    print(f"\nSWEEP COMPLETE rank_max={R}  [{time.time() - t_start:.0f}s]", flush=True)
    print(f"del Pezzo exceptional characters: {n_exc}", flush=True)
    print(f"empty lifts (would-be counterexamples): {n_empty}", flush=True)
    print(f"dispatched by rho_gen = 1: {dispatched}", flush=True)
    print(f"SURVIVORS (rho_gen >= 2): {len(survivors)}", flush=True)
    for row in survivors:
        print(f"  {row}", flush=True)


if __name__ == "__main__":
    main()
