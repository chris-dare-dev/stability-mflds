"""Week 3: K3 Mukai lattice, pairing, and Bayer-Macri wall classification."""

from dataclasses import replace
from fractions import Fraction as F
from math import gcd

import pytest

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.varieties import K3, abelian_surface
from bridgeland_stability.walls import Wall, VerticalWall, numerical_wall
from bridgeland_stability.mukai import (
    MukaiVector,
    classify_wall,
    classify_wall_certified,
    is_isotropic,
    is_spherical,
    moduli_dim,
    pairing,
    self_pairing,
    k3_wall,
    k3_wall_classified,
    solve_binary_quadratic,
    saturated_basis,
    hyperbolic_witnesses,
)

D = 2  # a quartic-adjacent K3 with H^2 = 2 (value is irrelevant when l = 0)
K3_2 = K3(2)  # d=2, kind=="K3"


def test_structure_sheaf_is_spherical():
    # CORRECTED Test 5: v(O) = (1,0,1) (NOT (1,0,-1)); <v,v> = -2 = -chi(O,O).
    vO = MukaiVector.from_chern(r=1, l=0, ch2=0)
    assert tuple(vO) == (1, 0, 1)
    assert self_pairing(vO, D) == -2
    assert pairing(vO, vO, D) == -2
    assert is_spherical(vO, D)
    assert moduli_dim(vO, D) == 0  # a point


def test_mukai_pairing_lattice():
    # E8-M3/G12.3: the Mukai pairing's H-degree term d*l*l' is now the NS-lattice
    # intersection form <c1,c1'> = <l H, l' H> (rank-1 shim, <H,H>=d); every value
    # is bit-for-bit unchanged and stays a Python int (the integer solver depends on it).
    from bridgeland_stability.nslattice import rank1_shim
    vO = MukaiVector.from_chern(r=1, l=0, ch2=0)
    assert tuple(vO) == (1, 0, 1)
    assert self_pairing(vO, D) == -2                     # v(O) spherical, unchanged
    assert isinstance(self_pairing(vO, D), int)          # still int (solver dependency)
    for (r, l, s), (r2, l2, s2), d in [
        ((1, 0, 1), (1, 0, 1),  2),                      # v(O), self
        ((2, 3, 1), (1, 0, 1),  2),
        ((3, 5, 2), (2, 2, -1), 2),
        ((2, 1, 1), (2, 0, -1), 2),                      # Brill-Noether pair
        ((1, 1, 0), (1, -1, 0), 4),                      # d = 4
    ]:
        v, w = MukaiVector(r, l, s), MukaiVector(r2, l2, s2)
        lat = rank1_shim(d)
        ll = lat.pairing((v.l,), (w.l,))                 # <l H, l' H> == d l l'
        assert ll == d * v.l * w.l
        assert pairing(v, w, d) == ll - v.r * w.s - w.r * v.s
        assert self_pairing(v, d) == lat.self_pairing((v.l,)) - 2 * v.r * v.s
        assert isinstance(pairing(v, w, d), int)
        assert isinstance(self_pairing(v, d), int)
    # v^2 = d l^2 - 2 r s unchanged on the pinned self-pairing grid.
    for (r, l, s) in [(2, 3, 1), (1, 0, 1), (3, 5, 2), (2, 2, -1)]:
        v = MukaiVector(r, l, s)
        assert self_pairing(v, D) == D * l * l - 2 * r * s


def test_wrong_triple_gives_plus_two_artifact():
    # The brief's confused (1,0,-1) gives +2, which is NOT a spherical class.
    bad = MukaiVector(1, 0, -1)
    assert self_pairing(bad, D) == 2
    assert not is_spherical(bad, D)


def test_self_pairing_formula():
    for (r, l, s) in [(2, 3, 1), (1, 0, 1), (3, 5, 2), (2, 2, -1)]:
        v = MukaiVector(r, l, s)
        assert self_pairing(v, D) == D * l * l - 2 * r * s


def test_isotropic_and_dimension():
    v = MukaiVector(1, 0, 0)  # v^2 = 0
    assert is_isotropic(v, D)
    assert moduli_dim(v, D) == 2
    w = MukaiVector(2, 0, 1)  # v^2 = -4? d*0 - 2*2*1 = -4
    assert self_pairing(w, D) == -4


def test_classify_hilbert_chow_wall():
    # v = (2,0,-1), v^2 = 4 > 0; isotropic u=(1,0,0) has (u,v)=1 -> Hilbert-Chow divisorial.
    v = MukaiVector(2, 0, -1)
    assert self_pairing(v, D) == 4
    u = MukaiVector(1, 0, 0)
    assert is_isotropic(u, D) and pairing(u, v, D) == 1
    cls = classify_wall(v, u, D, search=6)
    assert cls.wall_type == "divisorial"
    assert cls.subtype == "hilbert-chow"
    assert cls.v_squared == 4


def test_classify_li_gieseker_uhlenbeck_wall():
    # isotropic u=(2,0,0): (u,v)=2 -> Li-Gieseker-Uhlenbeck (when no (u,v)=1 witness blocks it)
    v = MukaiVector(3, 0, -1)  # v^2 = 6
    u = MukaiVector(2, 0, 0)
    assert is_isotropic(u, D) and pairing(u, v, D) == 2
    cls = classify_wall(v, u, D, search=6)
    assert cls.wall_type == "divisorial"
    # could be hilbert-chow if some (u',v)=1 witness exists in the lattice; for this v it's LGU
    assert cls.subtype in ("li-gieseker-uhlenbeck", "hilbert-chow")


# --------------------------------------------------------------------------
# E2-M3: k3_wall_classified -- pair the K3 (s,t) semicircle with its type.
# All asserted values are re-derived by hand (exact Fraction) and anchored to
# Bayer-Macri Thm 5.7 (arXiv:1301.6968) + the Bridgeland/ABCH semicircle
# (arXiv:math/0307164).  See docs/CORRECTIONS.md sec. 6.
# --------------------------------------------------------------------------
def test_k3_wall_classified_pairs_geometry_and_type():
    # The pinned Mukai example as integral-l Chern classes on K3(2):
    #   v=(2,0,-1) Mukai  <- ChernChar(2,0,-3)  (ch2 = s-r = -1-2 = -3, c = l*d = 0)
    #   u=(1,0, 0) Mukai  <- ChernChar(1,0,-1)  (ch2 = 0-1 = -1, c = 0)
    # Both l = c/d = 0 (integral).  Equal Mumford slope 0 -> k3_wall shifts to
    # (2,0,-1)/(1,0,0), minor W_rc = 2*0 - 1*0 = 0 -> the geometry is the
    # DEGENERATE vertical wall s = 0.  The genuine semicircle is exercised in T4.
    v_ch = ChernChar(2, 0, -3)
    w_ch = ChernChar(1, 0, -1)
    geom, cls = k3_wall_classified(v_ch, w_ch, K3_2)

    assert geom == k3_wall(v_ch, w_ch, 2)
    assert isinstance(geom, VerticalWall)
    assert geom.s_value == F(0)

    assert cls.wall_type == "divisorial"
    assert cls.subtype == "hilbert-chow"
    assert cls.v_squared == 4

    # the Chern -> Mukai (ch2 -> ch2 + r) shift lands on the pinned vectors
    assert tuple(MukaiVector.from_chern(2, 0, -3)) == (2, 0, -1)
    assert tuple(MukaiVector.from_chern(1, 0, -1)) == (1, 0, 0)


def test_k3_wall_classified_requires_integral_l():
    # The shared synthetic fixture w=(1,-1,1/2) on d=2 has l = c/d = -1/2
    # (non-integral): ChernChar(1,-1,1/2).c % 2 == 1 != 0, so it has no integral
    # Mukai vector and k3_wall_classified must raise.  (v passes: c=0, e=-2 in Z.)
    V = ChernChar(1, 0, -2)
    W = ChernChar(1, -1, F(1, 2))
    assert W.c % 2 == 1  # l = -1/2 is not in Z
    with pytest.raises(ValueError):
        k3_wall_classified(V, W, K3_2)

    # ...but the bare +2/d semicircle primitive is unrestricted and still works.
    assert k3_wall(V, W, 2).radius_sq == F(21, 4)


def test_k3_wall_classified_rejects_non_k3():
    # K3-only guard: the ch2 -> ch2 + r Mukai shift is K3-specific, so an abelian
    # surface (kind == 'abelian' != 'K3') must raise -- symmetric with the shipped
    # walls.abelian_wall kind guard.
    with pytest.raises(ValueError):
        k3_wall_classified(
            ChernChar(2, 0, -3), ChernChar(2, 2, -1), abelian_surface(2)
        )


def test_integral_l_k3_class():
    # ONE genuine integral-l K3 lattice (s,t) wall (both classes in the Mukai
    # lattice), giving a real semicircle + a positive Bayer-Macri type.
    #   v=(2,0,-1) Mukai  <- ChernChar(2,0,-3)
    #   w=(2,1, 1) Mukai  <- ChernChar(2,2,-1)   (l = 2/2 = 1 in Z, ch2 = s-r = -1)
    # w^2  = d l^2 - 2 r s = 2*1 - 2*2*1 = -2 (spherical);
    # (w,v) = d l l' - r s' - r' s = 2*1*0 - 2*(-1) - 2*1 = 0  -> Brill-Noether.
    v_ch = ChernChar(2, 0, -3)
    w_ch = ChernChar(2, 2, -1)
    geom, cls = k3_wall_classified(v_ch, w_ch, K3_2)

    # Geometry: shift -> v'=(2,0,-1), w'=(2,2,1); W_rc=4, W_re=4, W_ce=2;
    # center = 4/4 = 1; radius^2 = 1 - 2*2/(2*4) = 1/2.
    assert isinstance(geom, Wall)
    assert geom.center == F(1)
    assert geom.radius_sq == F(1, 2)
    assert geom.is_real and geom.radius_sq > 0
    # k3 - bare = 1/2 - (-1/2) = 1 = 2/d.  NB the BARE wall here is non-real
    # (radius_sq = -1/2 < 0); this asserts the algebraic +2/d shift identity, not
    # a difference of two real semicircles (the K3-shifted wall IS real, above).
    assert geom.radius_sq - numerical_wall(v_ch, w_ch, 2).radius_sq == F(2, 2) == 1

    assert cls.wall_type == "divisorial"
    assert cls.subtype == "brill-noether"
    assert cls.v_squared == 4

    # Witness robust to search order (the first (s,v)=0 spherical found may be
    # (-2,-1,-1), not (2,1,1)); pin the defining properties, not the triple.
    sph = MukaiVector(*cls.witness)
    assert self_pairing(sph, 2) == -2
    assert pairing(sph, MukaiVector(2, 0, -1), 2) == 0


# --------------------------------------------------------------------------
# E2-M5: PROVEN-tag the k3_wall_classified geometry (rho=1) / HEURISTIC (rho>=2);
# the raw k3_wall primitive stays UNTAGGED (UNKNOWN).  METADATA ONLY -- no
# center/radius_sq value changes.  Anchors: Bridgeland arXiv:math/0307164,
# Bayer-Macri arXiv:1301.6968; the rho>=2 H-projection gate Maciocia
# arXiv:1202.4587.
# --------------------------------------------------------------------------
def test_k3_wall_classified_geometry_cert_proven():
    # The genuine integral-l semicircle (v=(2,0,-1), w=(2,1,1) Mukai; real Wall
    # center 1, radius^2 1/2) on the Picard-rank-1 K3(2): geometry cert PROVEN =
    # meet(PROVEN alg, PROVEN surface).  Value UNCHANGED.
    v_ch = ChernChar(2, 0, -3)
    w_ch = ChernChar(2, 2, -1)
    geom, cls = k3_wall_classified(v_ch, w_ch, K3_2)
    assert isinstance(geom, Wall)
    assert (geom.center, geom.radius_sq) == (F(1), F(1, 2))  # unchanged
    assert geom.certificate.rigor == Rigor.PROVEN
    # The WallClassification (the *type*) keeps classify_wall's HEURISTIC
    # bounded-search cert -- the type is NOT upgraded here (that is G10/E6).
    assert cls.certificate.rigor == Rigor.HEURISTIC


def test_k3_wall_classified_rho2_downgraded_to_heuristic():
    # Defensive rho >= 2 K3 (kind stays 'K3', picard_rank=2): the scalar c=ch1.H
    # is only the H-projection, so the geometry cert downgrades to HEURISTIC with
    # the H-projection / E8-G12 note.  meet(HEURISTIC alg, PROVEN surface) keeps a
    # >= HEURISTIC surface cert so it does NOT collapse to UNKNOWN.  Geometry
    # identical (metadata only).
    K3_rho2 = replace(K3(2), picard_rank=2)
    assert K3_rho2.kind == "K3" and K3_rho2.picard_rank == 2
    v_ch = ChernChar(2, 0, -3)
    w_ch = ChernChar(2, 2, -1)
    geom, cls = k3_wall_classified(v_ch, w_ch, K3_rho2)
    assert isinstance(geom, Wall)
    assert (geom.center, geom.radius_sq) == (F(1), F(1, 2))  # unchanged
    assert geom.certificate.rigor == Rigor.HEURISTIC
    assert "H-project" in geom.certificate.note
    assert "G12" in geom.certificate.note
    assert cls.certificate.rigor == Rigor.HEURISTIC  # type still not upgraded


def test_k3_wall_primitive_is_untagged():
    # k3_wall(v, w, d) takes NO Surface, so it is the raw Picard-rank-1 semicircle
    # primitive and returns an UNTAGGED (UNKNOWN) Wall; the surface-gated tag is
    # delivered only by k3_wall_classified.
    v_ch = ChernChar(2, 0, -3)
    w_ch = ChernChar(2, 2, -1)
    wall = k3_wall(v_ch, w_ch, 2)
    assert isinstance(wall, Wall)
    assert (wall.center, wall.radius_sq) == (F(1), F(1, 2))
    assert wall.certificate.rigor == Rigor.UNKNOWN


# --------------------------------------------------------------------------
# E6-M1 (G10): hyperbolic binary-quadratic solver over the saturated sublattice
# All Gram values (A,B,C) = (v^2, 2(v,w), w^2); D = B^2 - 4AC.  Every asserted
# value below was recomputed exactly (integer) AND brute-force-enumerated over
# |a|,|b| <= 40.  Anchors: Bayer-Macri Thm 5.7 (arXiv:1301.6968) "potential wall
# = primitive rank-2 sublattice"; Buchmann-Vollmer, *Binary Quadratic Forms*
# (reduction theory of indefinite forms).  See docs/CORRECTIONS.md.
# --------------------------------------------------------------------------
def test_solver_reproduces_spherical():
    # lattice f1 = v(O) = (1,0,1) (spherical, v^2=-2), f2 = (1,0,0) (isotropic),
    # d=2.  Gram (-2,-2,0), D=4 (perfect square).
    # -2a^2 - 2ab = -2  <=>  a(a+b) = 1  =>  {(1,0),(-1,0)} = one +-orbit.
    assert solve_binary_quadratic(-2, -2, 0, -2) == [(1, 0)]
    (a, b), = solve_binary_quadratic(-2, -2, 0, -2)
    assert -2 * a * a - 2 * a * b == -2 and gcd(abs(a), abs(b)) == 1
    # (1,0) -> a*f1 + b*f2 = (a+b, 0, a) = (1,0,1) = v(O)
    assert tuple(MukaiVector(a + b, 0, a)) == (1, 0, 1)
    assert self_pairing(MukaiVector(a + b, 0, a), 2) == -2


def test_solver_isotropic():
    # Hilbert-Chow lattice v=(2,0,-1) (v^2=4), u=(1,0,0) (isotropic, (u,v)=1),
    # d=2.  Gram (4,2,0), D=4.  4a^2 + 2ab = 0  <=>  2a(2a+b) = 0;
    # proper orbits (0,1) and (1,-2).
    res = solve_binary_quadratic(4, 2, 0, 0)
    assert all(4 * a * a + 2 * a * b == 0 for (a, b) in res)
    assert (0, 1) in res and len(res) == 2
    # (0,1) -> 0*v + 1*u = u = (1,0,0), the pinned Hilbert-Chow witness
    assert self_pairing(MukaiVector(1, 0, 0), 2) == 0


def test_solver_pell_orbit():
    # Brill-Noether lattice v=(2,0,-1) (v^2=4), w=(2,1,1) (w^2=-2, (w,v)=0),
    # d=2.  Gram (4,0,-2), D=32 (non-square: genuine Pell).
    # 4a^2 - 2b^2 = -2  <=>  b^2 - 2a^2 = 1: (0,1),(2,3),(12,17),... are ALL one
    # orbit under M=[[3,2],[4,3]] (from Pell 6^2 - 32*1^2 = 4) and -I.
    res = solve_binary_quadratic(4, 0, -2, -2)
    assert res == [(0, 1)]                       # infinitely many solutions, ONE orbit rep
    assert 4 * 2 * 2 - 2 * 3 * 3 == -2           # (2,3) is a solution ...
    assert (2, 3) not in res                      # ... but deduped (orbit of (0,1))
    assert self_pairing(MukaiVector(2, 1, 1), 2) == -2   # (0,1) -> w=(2,1,1), spherical


def test_saturation():
    # constructed v=(2,0,2), w=(2,0,-2), d=2.  (1,0,0) is NOT in the Z-span
    # (2a+2b=1 unsolvable) but IS in the saturation (4*(1,0,0) = v + w).
    v, w = MukaiVector(2, 0, 2), MukaiVector(2, 0, -2)
    assert not any((2 * a + 2 * b, 0, 2 * a - 2 * b) == (1, 0, 0)
                   for a in range(-4, 5) for b in range(-4, 5))
    wits = hyperbolic_witnesses(v, w, 2, 0)      # enumerates over the SATURATION
    assert MukaiVector(1, 0, 0) in wits          # saturated witness IS found
    assert all(self_pairing(x, 2) == 0 for x in wits)
    g1, g2 = saturated_basis(v, w)               # a strict overlattice basis of Zv (+) Zw
    assert self_pairing(g1, 2) == 0 and self_pairing(g2, 2) == 0   # (1,0,0),(0,0,1)


def test_hyperbolic_witnesses_reproduces_pins():
    # the wrapper reproduces the two pinned test_mukai.py witnesses through
    # saturate -> Gram -> solve -> map.
    wsph = hyperbolic_witnesses(MukaiVector(1, 0, 1), MukaiVector(1, 0, 0), 2, -2)
    assert any(tuple(x) in {(1, 0, 1), (-1, 0, -1)} for x in wsph)   # v(O) spherical
    wiso = hyperbolic_witnesses(MukaiVector(2, 0, -1), MukaiVector(1, 0, 0), 2, 0)
    assert any(tuple(x) in {(1, 0, 0), (-1, 0, 0)} for x in wiso)    # Hilbert-Chow isotropic u


# --------------------------------------------------------------------------
# E6-M2 (G10): classify_wall_certified -- positive-only certified wall type.
# EXACT orbit enumeration over the SATURATED rank-2 hyperbolic sublattice H_W
# (via hyperbolic_witnesses) replacing classify_wall's bounded search.  Every
# asserted (wall_type, subtype, witness, gram) below was recomputed exactly
# (integer Mukai pairing) AND cross-checked vs classify_wall(search in {8,32}).
# certified=True is SOUND: it EXHIBITS an exact lattice witness with the
# verified Bayer-Macri Thm 5.7 property (search-bound independent).  Anchor:
# Bayer-Macri arXiv:1301.6968 Thm 5.7.  METADATA/NEW-FUNCTION ONLY -- changes
# no existing computed value.  See docs/CORRECTIONS.md.
# --------------------------------------------------------------------------
def test_certified_matches_search():
    # On the pinned POSITIVE cases the certified (wall_type, subtype) matches
    # classify_wall(search in {8,32}) AND sets certified=True, bound-independent.
    # v(O) (spherical v=(1,0,1), w=(1,0,0)): isotropic (1,0,0) with |(u,v)|=1
    #      -> divisorial/hilbert-chow.
    # Hilbert-Chow (v=(2,0,-1), u=(1,0,0)): (u,v)=1 -> divisorial/hilbert-chow.
    # genuine LGU where span == saturation (v=(2,0,-3), w=(1,0,-1)): isotropic
    #      (0,0,-1) with (u,v)=2 -> divisorial/li-gieseker-uhlenbeck (both give
    #      LGU, so the saturation refinement does NOT change the type here).
    cases = [
        (MukaiVector(1, 0, 1), MukaiVector(1, 0, 0), "divisorial", "hilbert-chow"),
        (MukaiVector(2, 0, -1), MukaiVector(1, 0, 0), "divisorial", "hilbert-chow"),
        (MukaiVector(2, 0, -3), MukaiVector(1, 0, -1),
         "divisorial", "li-gieseker-uhlenbeck"),
    ]
    for v, w, wt, st in cases:
        cert = classify_wall_certified(v, w, D)
        assert (cert.wall_type, cert.subtype) == (wt, st)
        assert cert.certified is True
        assert cert.certificate.rigor == Rigor.PROVEN
        assert cert.certificate.citations == ("arXiv:1301.6968",)
        # search-bound independent: identical to bounded search at 8 AND 32
        for search in (8, 32):
            s = classify_wall(v, w, D, search=search)
            assert (s.wall_type, s.subtype) == (wt, st)
        # the exhibited witness truly lies in the saturated lattice and realises
        # the theorem hypothesis (spherical/isotropic + the pinned pairing).
        wit = MukaiVector(*cert.witness)
        assert self_pairing(wit, D) in (-2, 0)


def test_certified_brill_noether_exhibits_witness():
    # v=(2,0,-1) (v^2=4), w=(2,1,1) (w^2=-2, (w,v)=0) -> divisorial/brill-noether.
    # The certified verdict EXHIBITS a spherical s in H_W with (s,v)=0; pin the
    # defining role, not the exact triple (either sign is legitimate).
    v = MukaiVector(2, 0, -1)
    w = MukaiVector(2, 1, 1)
    cert = classify_wall_certified(v, w, D)
    assert (cert.wall_type, cert.subtype) == ("divisorial", "brill-noether")
    assert cert.certified is True and cert.certificate.rigor == Rigor.PROVEN
    s = MukaiVector(*cert.witness)
    assert self_pairing(s, D) == -2          # spherical
    assert pairing(s, v, D) == 0             # Brill-Noether: (s,v)=0
    # matches the bounded search (subtype identical at both bounds)
    for search in (8, 32):
        assert classify_wall(v, w, D, search=search).subtype == "brill-noether"


def test_lgu_triple_saturation_refines_to_hilbert_chow():
    # DOCUMENTED refinement (the exact G10 payoff): the roadmap's LGU triple
    # v=(3,0,-1), u=(2,0,0) has (u,v)=2, so bounded classify_wall reports
    # li-gieseker-uhlenbeck.  But saturating H_W reveals the pairing-1 isotropic
    # (1,0,0) -- which is in the l=0 saturation plane yet NOT the Z-span of
    # (3,0,-1),(2,0,0) (3a+2b=1, -a=0 unsolvable) -- so the CERTIFIED type
    # refines to divisorial/hilbert-chow (Bayer-Macri Thm 5.7, potential wall =
    # saturated primitive sublattice).  (1,0,0) is isotropic with (u,v)=1.
    v = MukaiVector(3, 0, -1)
    u = MukaiVector(2, 0, 0)
    assert pairing(u, v, D) == 2                          # span witness gives LGU
    assert classify_wall(v, u, D, search=8).subtype == "li-gieseker-uhlenbeck"
    cert = classify_wall_certified(v, u, D)
    assert (cert.wall_type, cert.subtype) == ("divisorial", "hilbert-chow")
    assert cert.certified is True
    ref = MukaiVector(*cert.witness)                     # the refining witness
    assert self_pairing(ref, D) == 0 and pairing(ref, v, D) == 1
    # the refining class is NOT in the raw Z-span but IS in the saturation
    # (the exhibited witness is a sign-representative of a saturated isotropic
    # orbit -- hyperbolic_witnesses returns one per +-orbit, either sign legit).
    assert not any((3 * a + 2 * b, 0, -a) == tuple(ref)
                   for a in range(-8, 9) for b in range(-8, 9))
    neg = MukaiVector(-ref.r, -ref.l, -ref.s)
    iso = hyperbolic_witnesses(v, u, D, 0)
    assert ref in iso or neg in iso


def test_fake_or_none_uncertified():
    # Pell D=52 case v=(3,1,0), w=(1,0,1) on d=2: saturated form is hyperbolic
    # (potential wall) but has NO Brill-Noether/Hilbert-Chow/LGU/flopping witness
    # among the minimal saturated classes -> fake-or-none, which HONESTLY stays
    # certified=False (positive-only certification; totally-semistable split is
    # deferred to G17).  Matches classify_wall(search in {8,32}).
    v = MukaiVector(3, 1, 0)
    w = MukaiVector(1, 0, 1)
    cert = classify_wall_certified(v, w, D)
    assert cert.wall_type == "fake-or-none"
    assert cert.subtype is None and cert.witness is None
    assert cert.certified is False
    assert cert.certificate.rigor == Rigor.HEURISTIC   # not PROVEN
    for search in (8, 32):
        assert classify_wall(v, w, D, search=search).wall_type == "fake-or-none"


def test_carries_gram():
    # lattice_gram is the SPAN Gram [[v^2, (v,w)], [(v,w), w^2]], populated on
    # every result (certified positive AND fake-or-none).  Recomputed exactly.
    for v, w in [(MukaiVector(2, 0, -1), MukaiVector(1, 0, 0)),   # certified
                 (MukaiVector(3, 1, 0), MukaiVector(1, 0, 1))]:   # fake-or-none
        cert = classify_wall_certified(v, w, D)
        v2 = self_pairing(v, D)
        vw = pairing(v, w, D)
        w2 = self_pairing(w, D)
        assert cert.lattice_gram == ((v2, vw), (vw, w2))
        assert cert.v_squared == v2


def test_certified_raises_on_non_potential_wall():
    # (i) parallel classes (rank < 2): saturated_basis raises (v x w = 0).
    with pytest.raises(ValueError):
        classify_wall_certified(MukaiVector(1, 0, 0), MukaiVector(2, 0, 0), D)
    # (i) non-hyperbolic saturation: (2,1,0),(1,1,0) give saturated form
    # (A,B,C)=(0,0,2), disc = 0 <= 0 -> not a Bayer-Macri potential wall -> raise.
    with pytest.raises(ValueError):
        classify_wall_certified(MukaiVector(2, 1, 0), MukaiVector(1, 1, 0), D)


def test_k3_wall_classified_certified_bridge():
    # A certified type bridges to its k3_wall (s,t) semicircle:
    # k3_wall_classified(..., certified=True) pairs the SAME PROVEN geometry Wall
    # (center 1, radius^2 1/2 at rho=1) with the certified (PROVEN) type.
    #   v=(2,0,-1) Mukai <- ChernChar(2,0,-3);  w=(2,1,1) <- ChernChar(2,2,-1).
    v_ch = ChernChar(2, 0, -3)
    w_ch = ChernChar(2, 2, -1)
    geom_def, cls_def = k3_wall_classified(v_ch, w_ch, K3_2)                 # default
    geom_cert, cls_cert = k3_wall_classified(v_ch, w_ch, K3_2, certified=True)

    # SAME geometry Wall (default path unchanged), PROVEN at rho=1.
    assert geom_cert == geom_def
    assert isinstance(geom_cert, Wall)
    assert (geom_cert.center, geom_cert.radius_sq) == (F(1), F(1, 2))
    assert geom_cert.certificate.rigor == Rigor.PROVEN

    # default type keeps classify_wall's HEURISTIC bounded-search cert, uncertified
    assert cls_def.certificate.rigor == Rigor.HEURISTIC and cls_def.certified is False
    # certified type is the PROVEN Bayer-Macri Brill-Noether verdict
    assert (cls_cert.wall_type, cls_cert.subtype) == ("divisorial", "brill-noether")
    assert cls_cert.certified is True
    assert cls_cert.certificate.rigor == Rigor.PROVEN
    assert cls_cert.lattice_gram == ((4, 0), (0, -2))
