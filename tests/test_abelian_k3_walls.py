"""Abelian- and K3-surface (s, t) wall affordances (E2 shared fixture module).

Shared fixture (re-derived by hand AND re-executed against the shipped code):
``S = abelian_surface(2)`` so ``S.d == 2``, ``S.kind == "abelian"``;
``v = (1, 0, -2)``, ``w = (1, -1, 1/2)``.  The three ``walls._minors`` (exact
``Fraction``) are

    W_rc = r c' - r' c  = 1*(-1) - 1*0   = -1
    W_re = r e' - r' e  = 1*(1/2) - 1*(-2) = 5/2
    W_ce = c e' - c' e  = 0*(1/2) - (-1)*(-2) = -2

Since ``W_rc != 0`` this is a semicircular wall with

    center  = W_re / W_rc              = (5/2)/(-1)           = -5/2
    radius^2 = center^2 - 2 W_ce/(d W_rc) = 25/4 - 2(-2)/(2*(-1)) = 25/4 - 2 = 17/4

This is the genuine abelian (s, t) semicircle: an abelian surface has
``sqrt(td) = (1, 0, 0)``, so the bare Chern triple IS the Mukai vector and NO
``ch2 -> ch2 + r`` shift is applied (docs/CORRECTIONS.md §6).  The abelian value
``radius^2 = 17/4`` is exactly ``2/d = 1`` LESS than the K3 shim value ``21/4``
(the "abelian != K3 by exactly 2/d" regression guard, exercised in E2-M2).
"""

from dataclasses import replace
from fractions import Fraction as F

import pytest

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.varieties import abelian_surface, K3, P2
from bridgeland_stability.walls import Wall, abelian_wall, numerical_wall
from bridgeland_stability.mukai import k3_wall

# Shared module-level fixtures (E2-M2 / E2-M3 append their K3 tests here).
S_ABELIAN = abelian_surface(2)
D = S_ABELIAN.d
V = ChernChar(1, F(0), F(-2))
W = ChernChar(1, F(-1), F(1, 2))


def test_abelian_wall_value():
    assert D == 2
    assert S_ABELIAN.kind == "abelian"
    wall = abelian_wall(V, W, S_ABELIAN)
    assert isinstance(wall, Wall)
    assert wall.center == F(-5, 2)
    assert wall.radius_sq == F(17, 4)
    # Returns EXACTLY numerical_wall on the bare (r, c, e): no Mukai shift.
    bare = numerical_wall(V, W, D)
    assert (wall.center, wall.radius_sq) == (bare.center, bare.radius_sq)


def test_abelian_wall_asserts_kind():
    # K3 kind rejected -- guards the K3 (ch2 -> ch2 + r) shift being applied to
    # abelian input (and vice versa).
    with pytest.raises(ValueError):
        abelian_wall(V, W, K3(2))
    # Any non-abelian surface is rejected (P^2 here).
    with pytest.raises(ValueError):
        abelian_wall(V, W, P2)


def test_k3_wall_value():
    wall = k3_wall(V, W, D)
    assert isinstance(wall, Wall)
    # Center invariant vs bare/abelian; radius^2 = 17/4 + 2/d = 21/4.
    assert wall.center == F(-5, 2)
    assert wall.radius_sq == F(21, 4)


def test_k3_minus_abelian_is_2_over_d():
    k3 = k3_wall(V, W, D)
    ab = abelian_wall(V, W, S_ABELIAN)
    assert k3.center == ab.center == F(-5, 2)
    assert k3.radius_sq - ab.radius_sq == F(2, D) == F(2, 2) == 1


def test_center_invariant():
    # On a grid of small (r, c, e) triples the K3 shim leaves the wall CENTER
    # invariant and raises radius^2 by exactly 2/d (PROVEN, Picard rank 1, G3).
    ranks, cs, es, ds = (1, 2, 3), (-2, -1, 0, 1, 2), (F(-1), F(0), F(1, 2), F(1)), (1, 2, 3)
    checked = 0
    for rv in ranks:
        for cv in cs:
            for ev in es:
                for rw in (1, 2):
                    for cw in cs:
                        for ew in es:
                            for d in ds:
                                v = ChernChar(rv, F(cv), ev)
                                w = ChernChar(rw, F(cw), ew)
                                bare = numerical_wall(v, w, d)
                                k3 = k3_wall(v, w, d)
                                # Skip W_rc == 0 (vertical wall: no radius_sq).
                                if not (isinstance(bare, Wall) and isinstance(k3, Wall)):
                                    continue
                                assert k3.center == bare.center
                                assert k3.radius_sq - bare.radius_sq == F(2, d)
                                checked += 1
    assert checked > 0


# --------------------------------------------------------------------------
# E2-M5: baseline surface Certificates + PROVEN-tag the (s,t) walls (G2/G3/G5).
# METADATA ONLY -- no center/radius_sq value changes.  Rigor gate: the scalar-
# minor wall is PROVEN only at Picard rank 1 (Bridgeland arXiv:math/0307164,
# Bayer-Macri arXiv:1301.6968; ABL arXiv:0708.2247, Yanagida-Yoshioka
# arXiv:1203.0884); at rho >= 2 c = ch1.H is only the H-projection of ch1
# (Maciocia arXiv:1202.4587) -> HEURISTIC pending the NS-lattice refactor
# (E8/G12).
# --------------------------------------------------------------------------
def test_surface_baseline_certificates_proven():
    # K3(2) and abelian_surface(2) carry a PROVEN baseline Certificate with the
    # correct literature citations and the "Picard rank 1" hypothesis.
    k3_cert = K3(2).certificate
    assert k3_cert.rigor == Rigor.PROVEN
    assert "Picard rank 1" in k3_cert.hypotheses
    assert set(("arXiv:math/0307164", "arXiv:1301.6968")).issubset(k3_cert.citations)

    ab_cert = abelian_surface(2).certificate
    assert ab_cert.rigor == Rigor.PROVEN
    assert "Picard rank 1" in ab_cert.hypotheses
    assert set(("arXiv:0708.2247", "arXiv:1203.0884")).issubset(ab_cert.citations)


def test_abelian_wall_certificate_proven():
    # The rho = 1 abelian (s,t) wall geometry is PROVEN = meet(PROVEN alg, PROVEN
    # surface).  Value UNCHANGED: center -5/2, radius^2 17/4.
    wall = abelian_wall(V, W, S_ABELIAN)
    assert (wall.center, wall.radius_sq) == (F(-5, 2), F(17, 4))
    assert wall.certificate.rigor == Rigor.PROVEN


def test_abelian_wall_rho2_downgraded_to_heuristic():
    # Defensive rho >= 2 case: bare c = ch1.H is only the H-projection of ch1, so
    # the scalar-minor wall is NOT PROVEN -- meet(HEURISTIC alg, PROVEN surface) =
    # HEURISTIC, with a note flagging the H-projection and the E8/G12 gate.  The
    # surface still carries a >= HEURISTIC (here PROVEN) cert so meet does NOT
    # collapse to UNKNOWN.  Geometry identical to rho = 1 (metadata only).
    S_rho2 = replace(abelian_surface(2), picard_rank=2)
    assert S_rho2.kind == "abelian" and S_rho2.picard_rank == 2
    wall = abelian_wall(V, W, S_rho2)
    assert (wall.center, wall.radius_sq) == (F(-5, 2), F(17, 4))  # unchanged
    assert wall.certificate.rigor == Rigor.HEURISTIC
    assert "H-project" in wall.certificate.note
    assert "G12" in wall.certificate.note
    assert "arXiv:1202.4587" in wall.certificate.citations
