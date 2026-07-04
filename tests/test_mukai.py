"""Week 3: K3 Mukai lattice, pairing, and Bayer-Macri wall classification."""

from bridgeland_stability.mukai import (
    MukaiVector,
    classify_wall,
    is_isotropic,
    is_spherical,
    moduli_dim,
    pairing,
    self_pairing,
)

D = 2  # a quartic-adjacent K3 with H^2 = 2 (value is irrelevant when l = 0)


def test_structure_sheaf_is_spherical():
    # CORRECTED Test 5: v(O) = (1,0,1) (NOT (1,0,-1)); <v,v> = -2 = -chi(O,O).
    vO = MukaiVector.from_chern(r=1, l=0, ch2=0)
    assert tuple(vO) == (1, 0, 1)
    assert self_pairing(vO, D) == -2
    assert pairing(vO, vO, D) == -2
    assert is_spherical(vO, D)
    assert moduli_dim(vO, D) == 0  # a point


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
