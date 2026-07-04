"""Algorithm 3: Bridgeland numerical walls (semicircles)."""

from fractions import Fraction as F

import pytest

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.varieties import P2
from bridgeland_stability.walls import (
    VerticalWall,
    Wall,
    actual_walls,
    actual_walls_complete,
    compute_walls,
    numerical_wall,
)


def test_basic_wall_center_and_radius():
    v = ChernChar(2, F(0), F(-1, 4))
    w = ChernChar(1, F(1), F(1, 2))
    wall = numerical_wall(v, w, P2.d)
    assert isinstance(wall, Wall)
    assert wall.center == F(5, 8)
    assert wall.radius_sq == F(9, 64)
    assert abs(wall.radius - 0.375) < 1e-12


def test_hilbert_scheme_two_points_principal_wall():
    # P^2[2]: ideal sheaf of 2 points ch = (1,0,-2); destabilized by O(-1)=(1,-1,1/2).
    # ABCH (arXiv:1203.0316, sec 9): unique wall, center -5/2, radius 3/2.
    v = ChernChar(1, F(0), F(-2))
    w = ChernChar(1, F(-1), F(1, 2))
    wall = numerical_wall(v, w, P2.d)
    assert wall.center == F(-5, 2)
    assert wall.radius_sq == F(9, 4)
    assert abs(wall.radius - 1.5) < 1e-12


def test_vertical_wall_on_equal_slopes():
    v = ChernChar(2, F(0), F(-1, 4))
    w = ChernChar(1, F(0), F(1, 2))  # same Mumford slope (both 0)
    wall = numerical_wall(v, w, P2.d)
    assert isinstance(wall, VerticalWall)
    assert wall.s_value == F(0)


def test_compute_walls_finds_the_known_wall_and_is_nested():
    v = ChernChar(2, F(0), F(-1, 4))
    walls = compute_walls(v, P2, s_range=(-3, 3), rank_bound=4)
    assert all(w.is_real for w in walls)
    assert all(-3 <= w.center <= 3 for w in walls)
    # the known wall (center 5/8, radius^2 9/64) must appear
    found = [w for w in walls if w.center == F(5, 8) and w.radius_sq == F(9, 64)]
    assert found, "expected the (5/8, 9/64) wall in the enumerated set"
    # sorted by radius descending (Gieseker wall first)
    radii = [w.radius_sq for w in walls]
    assert radii == sorted(radii, reverse=True)


def test_radius_sq_can_be_negative_meaning_no_real_wall():
    # A subobject whose center is within sqrt(Delta_v) of mu_v gives rho^2 < 0.
    v = ChernChar(2, F(0), F(-1, 4))  # mu_v = 0, Delta_v(CH)=1/8, so need |s0|>=1/(2sqrt2)
    w = ChernChar(1, F(1), F(0))  # gives a small-offset center
    wall = numerical_wall(v, w, P2.d)
    # center s0 = (2*0 - 1*(-1/4))/(2*1 - 1*0) = (1/4)/2 = 1/8; rho^2 = 1/64 - 1/8 < 0
    assert wall.center == F(1, 8)
    assert wall.radius_sq < 0
    assert not wall.is_real


# ---- actual (certified-necessary) walls ----------------------------------
def test_actual_walls_hilbert_two_points_is_unique():
    # ABCH (arXiv:1203.0316, sec 9): P^2[2] has a SINGLE actual wall (-5/2, 3/2).
    v = ChernChar(1, 0, F(-2))
    walls, complete = actual_walls_complete(v, P2)
    assert complete
    assert len(walls) == 1
    w = walls[0]
    assert w.center == F(-5, 2)
    assert w.radius_sq == F(9, 4)
    assert abs(w.radius - 1.5) < 1e-12


def test_actual_walls_exclude_spurious_numerical_wall():
    # The numerical wall from (1,-9,34) (center -4) fails the heart-ordering /
    # integrality conditions and must NOT appear among the actual walls.
    v = ChernChar(1, 0, F(-2))
    walls = actual_walls(v, P2)
    assert all(w.center != F(-4) for w in walls)


def test_actual_walls_reject_non_integral_class():
    with pytest.raises(ValueError):
        actual_walls(ChernChar(2, 0, F(-1, 4)), P2)  # c2 = 1/4, not a real object


def test_actual_walls_classes_are_integral():
    # every destabilizing sub and quotient must be a genuine Chern character (c2 in Z)
    for n in (2, 3, 4, 5):
        for w in actual_walls(ChernChar(1, 0, F(-n)), P2):
            for ch in (w.subobject, w.quotient):
                c2 = F(ch.c * ch.c, 2) - ch.e
                assert c2.denominator == 1


def test_actual_walls_gieseker_wall_formula():
    # For P^2[n] the largest (Gieseker) wall is center -(2n+1)/2, radius (2n-1)/2,
    # realized by O(-1).  Verified against the explicit ABCH structure.
    for n in (2, 3, 4, 5):
        walls = actual_walls(ChernChar(1, 0, F(-n)), P2)
        gieseker = walls[0]  # sorted largest-first
        assert gieseker.center == F(-(2 * n + 1), 2)
        assert gieseker.radius_sq == F((2 * n - 1) ** 2, 4)


def test_actual_walls_nested():
    # distinct walls for a fixed v are nested: ordering by radius also orders centers
    walls = actual_walls(ChernChar(1, 0, F(-5)), P2)
    radii = [w.radius_sq for w in walls]
    assert radii == sorted(radii, reverse=True)
    # all centers are < mu_v = 0 (walls to the left for this ideal-sheaf class)
    assert all(w.center < 0 for w in walls)
