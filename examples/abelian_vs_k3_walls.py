"""Worked example: the SAME class on an abelian surface vs a K3, honest limits.

Run from the repo root:

    python examples/abelian_vs_k3_walls.py

Takes one Chern class ``v = (1, 0, -2)`` and one destabilizer and computes its
Bridgeland (s, t) wall on an abelian surface (``abelian_wall``) and on a K3
(``k3_wall``).  The two semicircles share a CENTER (-5/2) but the K3 radius^2 is
exactly ``2/d = 1`` larger: 17/4 (abelian) vs 21/4 (K3).  The reason is the
Mukai vector ``v(E) = ch(E) sqrt(td)``: an abelian surface has
``sqrt(td) = (1, 0, 0)`` (bare triple IS the Mukai vector, NO shift), a K3 has
``sqrt(td) = (1, 0, 1)`` so ``ch2 -> ch2 + r`` (docs/CORRECTIONS.md section 6;
Bridgeland arXiv:math/0307164; Bayer-Macri arXiv:1301.6968).

Then it attaches a genuine Bayer-Macri wall TYPE on an integral-l K3 lattice
class via ``k3_wall_classified``, and states plainly what the current model does
NOT cover.  All wall invariants are exact ``fractions.Fraction``; only the
``sqrt(td)`` annotations and the honest-limits prose are commentary.  This script
imports only the zero-dependency core (chern, varieties, walls, mukai) -- no viz,
no matplotlib/plotly.
"""

from __future__ import annotations

import sys
from fractions import Fraction as F
from pathlib import Path

# Allow ``python examples/abelian_vs_k3_walls.py`` from a checkout where the
# package is not pip-installed: prepend the repo root (this file's parent's
# parent) so ``import bridgeland_stability`` resolves.  Stdlib-only; no effect
# when the package is already importable.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.varieties import abelian_surface, K3
from bridgeland_stability.walls import abelian_wall
from bridgeland_stability.mukai import k3_wall, k3_wall_classified


def section(title: str) -> None:
    print("\n" + title)
    print("-" * len(title))


def worked_example_values() -> dict:
    """Compute and RETURN the exact ``Fraction`` wall invariants (no float, no
    stdout parsing).  The test asserts against this dict directly."""
    S_ab = abelian_surface(2)
    S_k3 = K3(2)

    # Shared class + the pinned synthetic +2/d destabilizer fixture.  Note
    # w = (1, -1, 1/2) has l = c/d = -1/2 not in Z, so it is NOT a genuine K3
    # lattice class -- it is a synthetic +2/d demonstration only.
    v = ChernChar(1, F(0), F(-2))
    w = ChernChar(1, F(-1), F(1, 2))

    ab = abelian_wall(v, w, S_ab)      # bare, NO Mukai shift
    k3 = k3_wall(v, w, S_k3.d)         # ch2 -> ch2 + r shim

    # A genuine integral-l K3 lattice class for the Bayer-Macri TYPE.
    #   v_g = ChernChar(2, 0, -3)  -> Mukai (2, 0, -1)
    #   w_g = ChernChar(2, 2, -1)  -> Mukai (2, 1,  1)  (spherical, (w,v)=0)
    v_g = ChernChar(2, F(0), F(-3))
    w_g = ChernChar(2, F(2), F(-1))
    geom_g, cls_g = k3_wall_classified(v_g, w_g, S_k3)

    return {
        "abelian_center": ab.center,
        "abelian_radius_sq": ab.radius_sq,
        "k3_center": k3.center,
        "k3_radius_sq": k3.radius_sq,
        "shift": k3.radius_sq - ab.radius_sq,
        "genuine_center": geom_g.center,
        "genuine_radius_sq": geom_g.radius_sq,
        "genuine_wall_type": cls_g.wall_type,
        "genuine_subtype": cls_g.subtype,
    }


def main() -> None:
    vals = worked_example_values()

    section("Abelian vs K3 (s, t) walls for v = (1, 0, -2), d = 2")
    print("  shared destabilizer w = (1, -1, 1/2)")
    print(f"  abelian_wall: center {vals['abelian_center']}, "
          f"radius^2 {vals['abelian_radius_sq']}   "
          "[sqrt(td) = (1, 0, 0): bare triple IS the Mukai vector, NO shift]")
    print(f"  k3_wall:      center {vals['k3_center']}, "
          f"radius^2 {vals['k3_radius_sq']}   "
          "[sqrt(td) = (1, 0, 1): ch2 -> ch2 + r shim]")
    print(f"  K3 - abelian radius^2 shift = {vals['shift']} = 2/d = 2/2 = 1 "
          "(center invariant)")
    print("  NOTE: w = (1, -1, 1/2) has l = c/d = -1/2 (not in Z), so this is a")
    print("        SYNTHETIC +2/d demonstration, not a genuine K3 lattice class.")

    section("Bayer-Macri wall TYPE on a genuine integral-l K3 lattice class")
    print("  v = ChernChar(2, 0, -3)  [Mukai (2, 0, -1)]")
    print("  w = ChernChar(2, 2, -1)  [Mukai (2, 1,  1)]  (spherical, (w,v)=0)")
    print(f"  k3_wall_classified: center {vals['genuine_center']}, "
          f"radius^2 {vals['genuine_radius_sq']} -> "
          f"{vals['genuine_wall_type']}/{vals['genuine_subtype']}")

    section("Current-model limits (honest)")
    print("  * Picard rank 1 only.  rho >= 2 abelian (E1 x E2) and rho >= 2 K3 are")
    print("    gated on the NS-lattice refactor (E8 / G12) -- NOT covered here.")
    print("  * The ch2 -> ch2 + r shim is K3-ONLY.  Never feed abelian input to")
    print("    k3_wall: an abelian surface has sqrt(td) = (1, 0, 0), so the shim")
    print("    would inject a spurious +2/d.  Use abelian_wall for abelian input.")
    print("  * actual_walls / actual_walls_complete are NOT abelian/K3 completeness")
    print("    oracles: their _is_integral_chern filter tests the d = 1 P^2 Chern-")
    print("    lattice condition only, and actual_walls layers ABCH / P^2-specific")
    print("    necessary conditions that do not transfer to d != 1.")
    print("  * The +2/d demo w = (1, -1, 1/2) (l = -1/2) has no Bayer-Macri type:")
    print("    k3_wall_classified raises on non-integral l (requires c % d == 0).")


if __name__ == "__main__":
    main()
