"""E9-M1 / G13: the full-Neron-Severi Bridgeland wall solver ``numerical_wall_ns``.

``numerical_wall_ns`` computes the ``(s, t)`` wall from the WHOLE intersection form
-- via ``<ch1, ch1>`` and the polarization-orthogonal ``gamma^2`` -- rather than
only the H-projection ``<ch1, H>`` the scalar :func:`numerical_wall` sees.

Maciocia ``beta = b*omega + gamma`` decomposition (arXiv:1202.4587).  With
``beta = b*H + gamma`` and ``gamma`` H-orthogonal (``<gamma, H> = 0``, hence
``gamma^2 <= 0`` by the Hodge index theorem), the imaginary part of ``Z_{beta,
omega}`` is unchanged and the real part shifts by a ``gamma``-twist of ``ch2``:

    e~(x) = ch2(x) - <ch1(x), gamma> + (r_x / 2) * gamma^2   (ch2 of ch(x).e^{-gamma})

so the full-NS wall in that slice is the scalar wall on the twisted characters
``(r, c, e~)``.  The chosen convention is ``gamma := H^perp(r_v ch1_w - r_w ch1_v)``.

Two-way verified (exact Fraction hand-arithmetic AND re-executed against the shipped
``numerical_wall`` / ``NSLattice``):

* rho = 1 (H^perp = {0}, gamma = 0, e~ = e) reduces bit-for-bit to
  :func:`numerical_wall`:
  - P^2[2]  v=(1,0,-2), w=O(-1)=(1,-1,1/2)  ->  center -5/2, radius^2 9/4
    (ABCH arXiv:1203.0316 sec 9, pinned in test_walls.py);
  - K3(2)   same v, w (bare primitive, no +2/d Mukai shim) -> center -5/2, radius^2 17/4.

* rho >= 2 gain, P^1 x P^1 (NS = hyperbolic U, <(a,b),(c,e)> = a e + b c, H=(1,1), d=2).
  v = (1, ch1=(0,0), -2).  Two subobjects that SHARE the scalar rep (1,-1,-1)
  (equal <ch1,H> = -1) but differ in bidegree:
  - w1: ch1=(0,-1), ch1^2 = 0   -> gamma=(1/2,-1/2),  gamma^2=-1/2 -> (center -3/2, radius^2 0)
  - w2: ch1=(1,-2), ch1^2 = -4  -> gamma=(3/2,-3/2),  gamma^2=-9/2 -> (center -11/2, radius^2 26)
  The scalar solver conflates them: numerical_wall((1,0,-2),(1,-1,-1),2) = (-1, -1),
  radius^2 < 0 (non-real: MISSES).  numerical_wall_ns FINDS the real wall (-11/2, 26)
  and distinguishes the two bidegrees (radius^2 0 vs 26).

[RESEARCH-light].  The rho = 1 reduction and the presence/absence sign behaviour are
convention-robust; the EXACT rho > 1 radius depends on the documented
``gamma = H^perp(relative-ch1)`` normalization and is the acceptance evidence until
one value is read off a Maciocia arXiv:1202.4587 worked example.
"""

from fractions import Fraction as F

import pytest

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.varieties import P2, K3, abelian_surface, P1xP1
from bridgeland_stability.walls import (
    Wall,
    VerticalWall,
    numerical_wall,
    numerical_wall_ns,
    maciocia_wall_bound,
    enumerate_ns_walls,
    actual_walls_complete,
)


def test_ns_reduces_to_scalar():
    # (a) P^2[2] regression: ideal sheaf (1,0,-2) destabilized by O(-1)=(1,-1,1/2).
    #     ABCH arXiv:1203.0316 sec 9: center -5/2, radius 3/2 (radius^2 9/4).
    v = ChernChar(1, F(0), F(-2))
    w = ChernChar(1, F(-1), F(1, 2))
    ns = numerical_wall_ns(v, w, P2)
    assert isinstance(ns, Wall)
    assert ns.center == F(-5, 2)
    assert ns.radius_sq == F(9, 4)
    assert ns.subobject == w and ns.v == v
    # bit-for-bit identical to the scalar primitive (incl. default certificate).
    assert ns == numerical_wall(v, w, P2.d)

    # (b) K3(2), same bare (r,c,e): no +2/d Mukai shim -> center -5/2, radius^2 17/4.
    ns_k3 = numerical_wall_ns(v, w, K3(2))
    assert ns_k3.center == F(-5, 2)
    assert ns_k3.radius_sq == F(17, 4)
    assert ns_k3 == numerical_wall(v, w, K3(2).d)

    # (c) bit-for-bit grid over every Picard-rank-1 catalog surface x small (r,c,e)
    #     for BOTH v and w: proves the rho=1 reduction (Wall AND VerticalWall cases,
    #     including non-real radius^2 < 0).  H^perp = {0} => gamma = 0 => e~ = e.
    surfaces = [P2, K3(2), abelian_surface(2), K3(4)]
    ranks = (1, 2, 3)
    degrees = (F(-2), F(-1), F(0), F(1))
    ch2s = (F(-1), F(0), F(1, 2))
    for surface in surfaces:
        for rv in ranks:
            for cv in degrees:
                for ev in ch2s:
                    v = ChernChar(rv, cv, ev)
                    for rw in ranks:
                        for cw in degrees:
                            for ew in ch2s:
                                w = ChernChar(rw, cw, ew)
                                ns = numerical_wall_ns(v, w, surface)
                                sc = numerical_wall(v, w, surface.d)
                                assert type(ns) is type(sc)
                                # full dataclass equality: center/radius_sq (Wall) or
                                # s_value (VerticalWall), plus subobject/v/certificate.
                                assert ns == sc


def test_ns_finds_bidegree_wall():
    # P^1 x P^1: NS = hyperbolic U, <(a,b),(c,e)> = a*e + b*c, H = (1,1), d = 2.
    L, H, d = P1xP1.lattice, P1xP1.H, P1xP1.d
    assert L.self_pairing(H) == d == 2

    v = ChernChar(1, F(0), F(-2))            # ch1 = (0,0)
    # Two subobjects sharing the scalar rep (1, -1, -1) but differing only by bidegree.
    w = ChernChar(1, F(-1), F(-1))
    w1_ch1 = (F(0), F(-1))
    w2_ch1 = (F(1), F(-2))

    # Equal H-projection, different self-intersection (bidegree-only difference).
    assert L.pairing(w1_ch1, H) == L.pairing(w2_ch1, H) == F(-1)
    assert L.self_pairing(w1_ch1) == 0
    assert L.self_pairing(w2_ch1) == F(-4)

    # Scalar solver conflates them and MISSES: non-real (-1, -1).
    sc = numerical_wall(v, w, d)
    assert isinstance(sc, Wall)
    assert sc.center == F(-1) and sc.radius_sq == F(-1)
    assert sc.is_real is False

    # Full-NS solver distinguishes the two bidegrees and FINDS the real wall.
    ns1 = numerical_wall_ns(v, w, P1xP1, v_ch1=(0, 0), w_ch1=(0, -1))
    ns2 = numerical_wall_ns(v, w, P1xP1, v_ch1=(0, 0), w_ch1=(1, -2))
    assert ns1.center == F(-3, 2) and ns1.radius_sq == F(0)
    assert ns2.center == F(-11, 2) and ns2.radius_sq == F(26)
    assert ns2.is_real is True
    assert ns2.subobject == w and ns2.v == v

    # Distinguished by bidegree, and it FINDS what the scalar model MISSES.
    assert ns1.radius_sq != ns2.radius_sq
    assert ns2.is_real is True and sc.is_real is False


def test_ns_ch1_projection_guard():
    # The supplied ch1 must have H-projection equal to the stored scalar c.
    # <(1,1),(1,1)> = 2 on the hyperbolic U, but w.c = -1, so this is rejected.
    v = ChernChar(1, F(0), F(-2))
    w = ChernChar(1, F(-1), F(-1))
    with pytest.raises(ValueError):
        numerical_wall_ns(v, w, P1xP1, v_ch1=(0, 0), w_ch1=(1, 1))


# ---- E9-M2: PROVEN Maciocia bound + finite NS-wall enumeration (G13) ------


def test_maciocia_wall_bound_surface_signature():
    # The Surface overload is value-preserving vs. the raw-degree signature.
    # v = (1,0,-2): bog = c^2/d - 2 r e = 0 - 2*1*(-2) = 4 (any d, since c=0),
    # so (bog-1)^2/4 = 9/4 on BOTH P^2 (d=1) and P^1 x P^1 (d=2).
    v = ChernChar(1, F(0), F(-2))
    assert maciocia_wall_bound(v, P2) == maciocia_wall_bound(v, 1) == F(9, 4)
    assert maciocia_wall_bound(v, P1xP1) == maciocia_wall_bound(v, 2) == F(9, 4)


def test_enumerated_within_bound():
    # rho = 1: every enumerated wall lies inside the (proven) Maciocia bound.
    for n in (2, 3, 4, 5):
        v = ChernChar(1, F(0), F(-n))
        walls = enumerate_ns_walls(v, P2)
        bound = maciocia_wall_bound(v, P2)
        assert walls  # non-empty (the Gieseker wall at least)
        assert all(w.radius_sq <= bound for w in walls)

    # n = 2: exactly the single ABCH Gieseker wall (-5/2, 9/4), tagged PROVEN.
    v2 = ChernChar(1, F(0), F(-2))
    walls2 = enumerate_ns_walls(v2, P2)
    assert {(w.center, w.radius_sq) for w in walls2} == {(F(-5, 2), F(9, 4))}
    assert all(w.certificate.rigor == Rigor.PROVEN for w in walls2)

    # rho >= 2: the bound filter holds by construction (P^1 x P^1).
    v3 = ChernChar(1, F(0), F(-2))
    e = enumerate_ns_walls(v3, P1xP1, v_ch1=(0, 0))
    b = maciocia_wall_bound(v3, P1xP1)
    assert all(w.radius_sq <= b for w in e)


def test_actual_walls_complete_proven():
    # P^2[2]: completeness now certified from the PROVEN bounding semicircle
    # (window^2 = 144 >= bound 9/4, single wall <= bound), not from doubling.
    v = ChernChar(1, F(0), F(-2))
    walls, complete = actual_walls_complete(v, P2)
    assert complete is True
    assert len(walls) == 1
    assert walls[0].center == F(-5, 2)
    assert walls[0].radius_sq == F(9, 4)
    assert walls[0].certificate.rigor == Rigor.PROVEN
    assert walls[0].radius_sq <= maciocia_wall_bound(v, P2)  # bound-backed


def test_ns_bound_is_heuristic_rho2():
    # Honesty walk-back / [RESEARCH-light] witness: a genuine rho >= 2 NS wall can
    # EXCEED the H-projected (bog-1)^2/4 estimate, so the estimate is NOT a proven
    # bound at Picard rank >= 2.  P^1 x P^1, v = (1,0,-2), w = (1,-1,-1) bidegree (1,-2)
    # has the E9-M1-pinned wall (-11/2, 26), while the estimate is only 9/4.
    v = ChernChar(1, F(0), F(-2))
    w = ChernChar(1, F(-1), F(-1))
    ns2 = numerical_wall_ns(v, w, P1xP1, v_ch1=(0, 0), w_ch1=(1, -2))
    assert ns2.radius_sq == F(26)
    assert maciocia_wall_bound(v, P1xP1) == F(9, 4)
    assert ns2.radius_sq > maciocia_wall_bound(v, P1xP1)

    # The bound filter DROPS that real wall (radius^2 26 > 9/4): the rho >= 2 in-bound
    # set is therefore NOT complete, and every returned wall is HEURISTIC (never PROVEN).
    e = enumerate_ns_walls(v, P1xP1, v_ch1=(0, 0))
    assert all(x.radius_sq != F(26) for x in e)
    assert all(x.certificate.rigor == Rigor.HEURISTIC for x in e)


def test_outermost_gieseker_wall_not_dropped_when_certified():
    # E9-M2 soundness regression (BLOCKER/MAJOR fix).  For P^2[n] the outermost
    # (Gieseker) wall is (-(2n+1)/2, ((2n-1)/2)^2), realised by O(-1); its CENTER lies
    # at distance (2n+1)/2 = radius+1 from mu_v = 0.  With the default center_window=12
    # that center falls OUTSIDE [-12, 12] for n >= 12, so the old ``win^2 >= bound``
    # certification silently truncated it while still reporting complete=True / PROVEN.
    # The invariant now: NEVER certify complete / tag a PROVEN "complete enumeration"
    # while the outermost Gieseker wall has been dropped.
    for n in (12, 13, 20):
        v = ChernChar(1, 0, F(-n))
        gieseker = (F(-(2 * n + 1), 2), F((2 * n - 1) ** 2, 4))
        bound = maciocia_wall_bound(v, P2)
        # radius^2 of the outermost wall equals the Maciocia bound exactly.
        assert gieseker[1] == bound

        # -- actual_walls_complete: never complete=True with the outermost wall omitted.
        walls, complete = actual_walls_complete(v, P2)
        keys = {(w.center, w.radius_sq) for w in walls}
        assert not (complete and gieseker not in keys)  # the exact defect being fixed
        assert all(w.radius_sq <= bound for w in walls)  # no wall exceeds the bound

        # -- enumerate_ns_walls: never a PROVEN "complete enumeration" that drops it.
        ns = enumerate_ns_walls(v, P2)
        ns_keys = {(w.center, w.radius_sq) for w in ns}
        all_proven = ns and all(w.certificate.rigor == Rigor.PROVEN for w in ns)
        assert not (all_proven and gieseker not in ns_keys)
        assert all(w.radius_sq <= bound for w in ns)

    # Concrete P^2[12] witness from the finding: the outermost wall (-25/2, 529/4) is
    # actually PRESENT now (window sized from the bound), tagged PROVEN and complete.
    v12 = ChernChar(1, 0, F(-12))
    assert maciocia_wall_bound(v12, P2) == F(529, 4)
    ns12 = enumerate_ns_walls(v12, P2)
    ns12_keys = {(w.center, w.radius_sq) for w in ns12}
    assert (F(-25, 2), F(529, 4)) in ns12_keys
    assert all(w.certificate.rigor == Rigor.PROVEN for w in ns12)
    walls12, _ = actual_walls_complete(v12, P2)
    assert (F(-25, 2), F(529, 4)) in {(w.center, w.radius_sq) for w in walls12}
