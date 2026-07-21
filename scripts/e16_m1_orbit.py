"""E16-M1: the geometric mutation-orbit search for v_107 on F_3, done right.

BFS over full exceptional collections (as K-theory quadruples) with:
* the six pairwise mutations  [L_A B] = chi(A,B)[A] - [B],
  [R_B A] = chi(A,B)[B] - [A];
* the two HELIX shifts (Bondal): right (E2,E3,E4, E1 x (-K)), left
  (E4 x K, E1,E2,E3) -- moves the spec-time probe lacked;
* TWIST-NORMALIZATION: each visited collection is canonicalized by x L so
  that c1(E1) lands in [0, r1)^2 -- collapsing the twist orbit;
* the GEOMETRIC filter: every member has rank >= 1 and integral c2 (a
  character no sheaf carries is pruned);
* enforced budgets (node cap / wall cap) + telemetry.

Target: any member of rank 107 whose slope fractional part is nu(v_107) or
its dual (Delta is forced for exceptional characters, and is twist-invariant).
A hit yields the geometric path for E16-M3; sustained absence feeds the
E16-M2 constructibility interpretation.

Usage: python scripts/e16_m1_orbit.py [--node-cap N] [--wall-cap SECONDS]
"""

import argparse
import sys
import time
from collections import deque
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("scripts", 1)[0])

from bridgeland_stability.varieties import hirzebruch                   # noqa: E402
from bridgeland_stability.exceptional_surface import (                   # noqa: E402
    SurfaceBundle,
    chi as chi_pair,
)

S = hirzebruch(3)
E_INDEX = 3
LAT = S.lattice
HALF = Fraction(1, 2)
K = (-(E_INDEX + 2), -2)
TARGET_FRACS = ((Fraction(76, 107), Fraction(25, 107)),
                (Fraction(31, 107), Fraction(82, 107)))     # nu and its dual class


def line(c1):
    return (1, c1, HALF * LAT.self_pairing(c1))


def chi_k(A, B):
    return chi_pair(SurfaceBundle(*A), SurfaceBundle(*B), S)


def twist(w, D):
    r, c1, ch2 = w
    return (r, (c1[0] + r * D[0], c1[1] + r * D[1]),
            ch2 + LAT.pairing(c1, D) + r * HALF * LAT.self_pairing(D))


def integral_c2(w):
    r, c1, ch2 = w
    return (HALF * LAT.self_pairing(c1) - ch2).denominator == 1


def normalize(coll):
    r1, c11, _ = coll[0]
    D = (-(c11[0] // r1), -(c11[1] // r1))
    if D == (0, 0):
        return coll
    return tuple(twist(w, D) for w in coll)


def moves(coll):
    for i in range(3):
        A, B = coll[i], coll[i + 1]
        x = chi_k(A, B)
        LAB = (x * A[0] - B[0], (x * A[1][0] - B[1][0], x * A[1][1] - B[1][1]),
               x * A[2] - B[2])
        RBA = (x * B[0] - A[0], (x * B[1][0] - A[1][0], x * B[1][1] - A[1][1]),
               x * B[2] - A[2])
        yield coll[:i] + (LAB, A) + coll[i + 2:], (i, "L")
        yield coll[:i] + (B, RBA) + coll[i + 2:], (i, "R")
    mK = (-K[0], -K[1])
    yield coll[1:] + (twist(coll[0], mK),), ("helix", "R")
    yield (twist(coll[3], K),) + coll[:3], ("helix", "L")


def is_target(w):
    if w[0] != 107:
        return False
    fx, fy = Fraction(w[1][0], 107) % 1, Fraction(w[1][1], 107) % 1
    return (fx, fy) in TARGET_FRACS


_BUNDLE_CACHE = {}


def certified_bundle(w) -> bool:
    """A REALIZABLE node: rank 1 (line bundle), or the exceptional bundle of
    the character provably exists (nonempty E14-M3 stability interval on F_3;
    twist-invariant, so cached by the normalized slope class).  The first
    geometric hit (2026-07-21) was killed by exactly this: its rank-2 node is
    the sweep-dispatched (2,(2,1)) with rho_gen = 1 -- no bundle."""
    r = w[0]
    if r == 1:
        return True
    c1 = (w[1][0] % r, w[1][1] % r)              # twist-normalized slope class
    key = (r, c1)
    got = _BUNDLE_CACHE.get(key)
    if got is None:
        from bridgeland_stability.stability_interval import stability_interval
        try:
            got = not stability_interval(r, c1, S).empty
        except ValueError:                        # not potentially exceptional
            got = False
        _BUNDLE_CACHE[key] = got
    return got


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--node-cap", type=int, default=2_000_000)
    ap.add_argument("--wall-cap", type=float, default=1800.0)
    ap.add_argument("--rank-cap", type=int, default=300)
    ap.add_argument("--coord-cap", type=int, default=4000)
    args = ap.parse_args()

    start = normalize((line((0, 0)), line((1, 0)),
                       line((E_INDEX, 1)), line((E_INDEX + 1, 1))))
    t0 = time.time()
    parent = {start: None}
    queue = deque([start])
    nodes = 0
    max_rank = 1
    hit = None
    while queue and nodes < args.node_cap and time.time() - t0 < args.wall_cap:
        coll = queue.popleft()
        nodes += 1
        if nodes % 50_000 == 0:
            print(f"[{nodes} nodes, {len(parent)} seen, max rank {max_rank}, "
                  f"{time.time() - t0:.0f}s]", flush=True)
        for nc, move in moves(coll):
            ok = True
            for w in nc:
                if (w[0] < 1 or w[0] > args.rank_cap
                        or abs(w[1][0]) > args.coord_cap
                        or abs(w[1][1]) > args.coord_cap
                        or not integral_c2(w)):
                    ok = False
                    break
                # every non-target node must be a CERTIFIED bundle (the
                # realizability filter; the target itself is the question)
                if not is_target(w) and not certified_bundle(w):
                    ok = False
                    break
            if not ok:
                continue
            nc = normalize(nc)
            if nc in parent:
                continue
            parent[nc] = (coll, move)
            for w in nc:
                max_rank = max(max_rank, w[0])
                if is_target(w):
                    hit = (nc, w)
            if hit:
                queue.clear()
                break
            queue.append(nc)
        if hit:
            break

    print(f"\nDONE: {nodes} nodes expanded, {len(parent)} collections, "
          f"max rank {max_rank}, {time.time() - t0:.0f}s", flush=True)
    if hit is None:
        print("NO geometric-orbit hit within budget (evidence, not proof; "
              "see the E16-M2 constructibility interpretation)", flush=True)
        return
    nc, w = hit
    print(f"GEOMETRIC HIT: {w} in collection "
          f"{[(m[0], m[1], str(m[2])) for m in nc]}", flush=True)
    path = []
    cur = nc
    while parent[cur] is not None:
        prev, move = parent[cur]
        path.append(move)
        cur = prev
    path.reverse()
    print(f"path ({len(path)} moves): {path}", flush=True)


if __name__ == "__main__":
    main()
