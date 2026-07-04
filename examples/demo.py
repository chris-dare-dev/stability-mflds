"""Reproduce every validated value and write the figures.

Run from the repo root:

    python examples/demo.py            # prints the validated table
    python examples/demo.py --figures  # also writes figures/*.html and .png

Numbers printed here are the *corrected* values (see docs/CORRECTIONS.md), all
computed exactly.
"""

from __future__ import annotations

import sys
from fractions import Fraction as F

import bridgeland_stability as bs
from bridgeland_stability.threefold import ThreefoldChern, bg_boundary_curve


def section(title: str) -> None:
    print("\n" + title)
    print("-" * len(title))


def main(write_figures: bool = False) -> None:
    section("1. Exceptional bundles on P^2 (Markov ranks)")
    bundles = bs.enumerate_exceptional(0, 1, R_max=30)
    print("slopes (rank):", [(str(b.slope), b.r) for b in bundles])
    print("Markov numbers <= 40:", sorted(bs.markov_numbers(40)))
    print("rank-3 'bundle' at 1/3 genuine? ", bs.is_exceptional(bs.Bundle.from_slope(F(1, 3))),
          "(c2 =", bs.Bundle.from_slope(F(1, 3)).c2, "is not integral)")
    print("rank-5 bundle at 2/5:", bs.Bundle.from_slope(F(2, 5)),
          "Delta =", bs.Bundle.from_slope(F(2, 5)).discriminant)

    section("2. Drezet-Le Potier curve delta(mu)  [CH convention]")
    for mu in [F(0), F(1, 4), F(1, 3), F(2, 5), F(1, 2)]:
        print(f"  delta({mu}) = {bs.delta(mu, bundles)}")
    print("  (brief claimed delta(1/2)=3/4, delta(1/3)=8/9 -- both wrong)")
    res = bs.moduli_nonempty(bs.Bundle.T_minus1())
    print("  T(-1):", res["reason"])

    section("3. Bridgeland walls (exact)")
    v = bs.ChernChar(1, 0, F(-2))            # P^2[2] ideal sheaf
    w = bs.ChernChar(1, -1, F(1, 2))         # O(-1)
    wall = bs.numerical_wall(v, w, bs.P2.d)
    print(f"  P^2[2] wall: center {wall.center}, radius {wall.radius}  (ABCH: -5/2, 3/2)")
    v2 = bs.ChernChar(2, 0, F(-1, 4))
    w2 = bs.ChernChar(1, 1, F(1, 2))
    wall2 = bs.numerical_wall(v2, w2, bs.P2.d)
    print(f"  v=(2,0,-1/4), w=(1,1,1/2): center {wall2.center}, R^2 {wall2.radius_sq}")
    print("  actual (certified-necessary) walls for P^2[n]:")
    for n in (2, 3, 4, 5):
        ws, complete = bs.actual_walls_complete(bs.ChernChar(1, 0, F(-n)), bs.P2)
        gie = ws[0]
        print(f"    P^2[{n}]: {len(ws)} wall(s) (complete={complete}); "
              f"Gieseker center {gie.center}, radius {gie.radius}")

    section("4. Bogomolov-Gieseker (surface)")
    for label, ch in [("T(-1)", bs.Bundle.T_minus1().chern_char()),
                      ("O + O", bs.ChernChar(2, 0, 0)),
                      ("(1,0,1/4)", bs.ChernChar(1, 0, F(1, 4)))]:
        r = bs.check_bg_surface(ch, bs.P2)
        print(f"  {label}: Delta={r.discriminant}  -> {r.note}")

    section("5. Threefold tilt-BG boundary alpha_crit(beta)")
    nc = ThreefoldChern(2, 0, 1, 0)          # P^3 null-correlation bundle
    for beta in [F(1, 2), F(1)]:
        print(f"  beta={beta}: Q={bs.bmt_Q_at(nc, beta, bs.P3.d3)}, "
              f"alpha_crit={bs.alpha_crit(nc, beta, bs.P3.d3)}")
    print("  beta=0 degenerate (ch1^beta=0):", bs.alpha_crit(nc, F(0), bs.P3.d3))
    O5 = ThreefoldChern(1, 0, 0, 0)
    print("  quintic O: Q(beta=2/3) =", bs.bmt_Q_at(O5, F(2, 3), bs.QUINTIC.d3),
          "(O sits on the BG boundary)")

    section("6. K3 Mukai lattice")
    vO = bs.MukaiVector.from_chern(r=1, l=0, ch2=0)
    print(f"  v(O) = {tuple(vO)},  <v,v> = {bs.self_pairing(vO, 2)} (spherical),"
          f"  dim M = {bs.moduli_dim(vO, 2)}")
    vh = bs.MukaiVector(2, 0, -1)
    cls = bs.classify_wall(vh, bs.MukaiVector(1, 0, 0), d=2)
    print(f"  wall of v={tuple(vh)} via u=(1,0,0): {cls.wall_type}/{cls.subtype} (v^2={cls.v_squared})")

    if write_figures:
        section("Writing figures/")
        import os
        from bridgeland_stability import viz, walls_from_subobjects
        os.makedirs("figures", exist_ok=True)
        viz.plot_dlp_curve(mu_min=0, mu_max=1, R_max=40, html_path="figures/dlp_curve.html")
        viz.plot_dlp_curve(mu_min=0, mu_max=1, R_max=30,
                           backend="matplotlib").savefig("figures/dlp_curve.png",
                                                         dpi=110, bbox_inches="tight")
        abch = walls_from_subobjects(v, [w], bs.P2)
        viz.plot_walls(v, bs.P2, s_range=(-4, 1), t_max=3, walls=abch,
                       html_path="figures/walls_hilb2_abch.html")
        viz.plot_walls(v2, bs.P2, s_range=(-3, 3), t_max=3, rank_bound=4,
                       html_path="figures/walls_numerical.html")
        # actual nested walls for P^2[4]
        v4 = bs.ChernChar(1, 0, F(-4))
        viz.plot_walls(v4, bs.P2, s_range=(-6, 1), t_max=4,
                       walls=bs.actual_walls(v4, bs.P2),
                       html_path="figures/walls_hilb4_actual.html")
        bd = bg_boundary_curve(nc, bs.P3, beta_range=(-2, 2), N=300)
        viz.plot_threefold_bg(bd, html_path="figures/p3_bg_boundary.html")
        print("  wrote figures/{dlp_curve.html,.png, walls_*.html, p3_bg_boundary.html}")


if __name__ == "__main__":
    main(write_figures="--figures" in sys.argv)
