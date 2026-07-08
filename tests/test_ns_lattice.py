"""E8-M1 / G12.1 tests: NSLattice type + rank-1 backward-compat shim pin.

No new mathematics.  The shim assertions <ch1,H>==c and <ch1,ch1>==c^2/d are the
bit-for-bit encoding that protects the pinned surface gate (delta(1/2)=5/8,
P^2[2] wall (-5/2, 3/2), K3 v(O)=(1,0,1)); they are exactly re-derived, and
cross-checked against the shipped chern.bogomolov_discriminant.
"""
import subprocess
import sys
from fractions import Fraction as F

import pytest

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.nslattice import (
    NSLattice, rank1_shim, shim_ch1, RANK1_AMPLE,
)
from bridgeland_stability.varieties import P2, P1xP1, K3, abelian_surface, hirzebruch


def test_shim_pins_pairings():
    # Roadmap worked example: r=1, c=1, d=2.
    d = 2
    lat = rank1_shim(d)                       # NSLattice(1, ((2,),))
    H = RANK1_AMPLE                           # (1,)
    assert lat.pairing(H, H) == d             # <H,H> = 1*2*1 = 2 = H^2
    ch1 = shim_ch1(1, d)                      # (1/2,)
    assert lat.pairing(ch1, H) == 1           # <ch1,H> == c == 1
    assert lat.pairing(ch1, ch1) == F(1, 2)   # <ch1,ch1> == c^2/d == 1/2

    # More (c, d), each exactly re-derived: <ch1,H>=c and <ch1,ch1>=c^2/d.
    for c, dd, cc in [(-1, 1, F(1)), (3, 2, F(9, 2)), (2, 4, F(1)), (2, 2, F(2))]:
        L, e1 = rank1_shim(dd), shim_ch1(c, dd)
        assert L.pairing(e1, RANK1_AMPLE) == c
        assert L.self_pairing(e1) == cc

    # Bit-for-bit cross-check vs the SHIPPED formula bogomolov_discriminant = c^2/d - 2re:
    cc_obj = ChernChar(1, 1, 0)               # d=2 -> 1/2 - 0 = 1/2
    shim_bog = rank1_shim(2).self_pairing(shim_ch1(cc_obj.c, 2)) - 2 * cc_obj.r * cc_obj.e
    assert shim_bog == cc_obj.bogomolov_discriminant(2) == F(1, 2)


def test_naive_encoding_rejected():
    # The naive (c,) + Gram=[[d]] encoding yields c*d and c^2*d (NOT c, c^2/d),
    # which would break delta(1/2)=5/8; the shim stores c/d instead.
    d, c = 2, 1
    lat = NSLattice(1, ((d,),))               # same Gram as the shim
    H, naive, shim = (1,), (c,), shim_ch1(c, d)
    assert lat.pairing(naive, H) == c * d               # = 2  (naive: c*d)
    assert lat.self_pairing(naive) == c * c * d         # = 2  (naive: c^2*d)
    assert lat.pairing(shim, H) == c                    # = 1  (correct)
    assert lat.self_pairing(shim) == F(1, 2)            # correct
    assert lat.pairing(naive, H) != lat.pairing(shim, H)  # naive != shim (d != 1)
    # Naive would corrupt bogomolov_discriminant; the shim reproduces it exactly:
    v = ChernChar(1, 1, 0)
    assert lat.self_pairing(naive) - 2 * v.r * v.e != v.bogomolov_discriminant(2)
    assert lat.self_pairing(shim) - 2 * v.r * v.e == v.bogomolov_discriminant(2)


def test_p1xp1_pairing():
    # P^1 x P^1 hyperbolic plane, basis (f, s), Gram [[0,1],[1,0]].
    lat = NSLattice(2, ((0, 1), (1, 0)))
    f, s, H = (1, 0), (0, 1), (1, 1)          # H = f + s
    assert lat.pairing(f, f) == 0             # f^2 = 0
    assert lat.pairing(s, s) == 0             # s^2 = 0
    assert lat.pairing(f, s) == 1             # f.s = 1
    assert lat.self_pairing(H) == 2           # <f+s, f+s> = 2  (reproduces d = H^2 = 2)
    # Two classes the scalar model conflates: equal <.,H> but different self-pairing.
    a, b = (1, 1), (2, 0)                     # f+s vs 2f
    assert lat.pairing(a, H) == lat.pairing(b, H) == 2
    assert lat.self_pairing(a) == 2
    assert lat.self_pairing(b) == 0
    assert lat.self_pairing(a) != lat.self_pairing(b)   # distinguished


def test_surface_carries_shim_lattice():
    # An ample H is stored on Surface; a genuinely Picard-rank-1 row with
    # ns_lattice unset derives the rank-1 shim NSLattice(1, ((d,),)) with <H,H>=d.
    # (E8-M4: P^1xP^1 / F_n now carry genuine rank-2 NS lattices, so they are no
    # longer shim rows -- see test_p1xp1_distinguishes_classes / test_fn_hodge_signature.)
    for surf in (P2, K3(2), K3(4), abelian_surface(2)):
        assert surf.H == RANK1_AMPLE                    # (1,)
        assert surf.ns_lattice is None                  # unset -> shim
        lat = surf.lattice
        assert lat == NSLattice(1, ((surf.d,),))
        assert lat.pairing(surf.H, surf.H) == surf.d    # <H,H> == H^2


def test_p1xp1_distinguishes_classes():
    # E8-M4 / G12.4: P^1 x P^1 carries the genuine rank-2 hyperbolic NS lattice
    # (basis f, s; Gram [[0,1],[1,0]]), NOT the rank-1 shim.  H = f + s reproduces
    # d = H^2 = 2, and the lattice DISTINGUISHES two classes the old scalar model
    # conflated: ch1 = f+s and ch1 = 2f share <ch1,H> = 2 but differ in ch1^2.
    surf = P1xP1
    lat = surf.lattice
    assert surf.ns_lattice is not None                       # genuine, not the shim
    assert lat == NSLattice(2, ((0, 1), (1, 0)))
    assert surf.H == (1, 1)                                  # H = f + s
    assert lat.self_pairing(surf.H) == surf.d == 2           # <H,H> = 2 reproduces d
    fps, twof = (1, 1), (2, 0)                                # f+s vs 2f
    assert lat.pairing(fps, surf.H) == lat.pairing(twof, surf.H) == 2   # equal H-degree
    assert lat.self_pairing(fps) == 2                        # <f+s, f+s> = 2
    assert lat.self_pairing(twof) == 0                       # <2f, 2f> = 0
    assert lat.self_pairing(fps) != lat.self_pairing(twof)   # DISTINGUISHED (scalar model conflated)


def test_fn_hodge_signature():
    # E8-M4 / G12.4: F_n carries the rank-2 NS lattice with Gram [[0,1],[1,-n]]
    # (fiber f, negative section s with s^2 = -n).  det(Gram) = 0*(-n) - 1*1 = -1 < 0
    # is the Hodge-index signature (1,1) for every n >= 1; H = (n,1) reproduces d = n.
    for n in (1, 2, 3, 5):
        surf = hirzebruch(n)
        lat = surf.lattice
        assert surf.ns_lattice is not None
        assert lat == NSLattice(2, ((0, 1), (1, -n)))
        assert lat.gram[1][1] == -n                          # s^2 = -n
        assert lat.det() == -1                               # 0*(-n) - 1*1
        assert lat.det() < 0                                 # Hodge-index signature (1,1)
        assert lat.self_pairing(surf.H) == surf.d == n       # <H,H> = n reproduces d
    # F_0 stays the rank-1 backward-compat shim (F_0 = P^1 x P^1 is the P1xP1 row):
    assert hirzebruch(0).ns_lattice is None


def test_nslattice_validation():
    with pytest.raises(ValueError):
        NSLattice(2, ((0, 1), (1, 0), (0, 0)))          # not 2x2
    with pytest.raises(ValueError):
        NSLattice(2, ((0, 1), (2, 0)))                  # not symmetric
    with pytest.raises(ValueError):
        NSLattice(2, ((0, 1), (1, 0))).pairing((1,), (1, 1))  # wrong length
    # 2x2 det for the Hodge-index signature (E8-M4 consumes this on F_n):
    assert NSLattice(2, ((0, 1), (1, 0))).det() == -1   # 0*0 - 1*1
    assert NSLattice(2, ((0, 1), (1, -3))).det() == -1  # F_3: 0*(-3) - 1*1 < 0


def test_nslattice_is_frozen():
    lat = NSLattice(1, ((2,),))
    with pytest.raises(Exception):
        lat.rank = 2                                     # frozen
    assert hash(lat) == hash(NSLattice(1, ((2,),)))      # frozen -> hashable


def test_nslattice_is_stdlib_only():
    # nslattice.py must not pull matplotlib/plotly at import time.
    child = r'''
import sys
class _Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in {"matplotlib", "plotly"}:
            raise ImportError("blocked: " + name)
        return None
sys.meta_path.insert(0, _Blocker())
import bridgeland_stability.nslattice
assert "matplotlib" not in sys.modules
assert "plotly" not in sys.modules
print("STDLIB_ONLY_OK")
'''
    r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "STDLIB_ONLY_OK" in r.stdout
