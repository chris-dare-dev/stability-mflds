"""E15-M2 / CORRECTIONS Sec. 26.4: the del Pezzo mutation oracle differential.

``tests/oracle/delpezzo_mutation_reference.py`` re-derives the exceptional
CLASSES on ``F_0``/``F_1`` from mutation theory, sharing no code with the
package: its own NS Gram, canonical class and Riemann-Roch Euler form, and
generation by mutation rather than by the Drezet-Le Potier rank induction the
E15-M2 sweep scans with.  These tests are the differential gate between the two.

The full-scale run behind CORRECTIONS Sec. 26.4 (rank <= 130: 376 on ``F_0`` plus
211 on ``F_1``, totalling the sweep's 587, saturated and cap-invariant, matching
the production scan member for member at rank <= 15, <= 30 and <= 40) is far too
slow for the suite; what is pinned here is the same differential at rank <= 12
plus the charter properties that make the oracle worth having.
"""

import os
from fractions import Fraction as F

from oracle.delpezzo_mutation_reference import (
    _canonical,
    _self_pairing,
    euler_pairing,
    exceptional_classes,
    line_bundle,
    standard_collection,
)

from bridgeland_stability.delta_sharp import surface_with_index
from bridgeland_stability.dlp_hirzebruch import (
    is_potentially_exceptional,
    is_stable_exceptional,
)
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi as pkg_chi
from bridgeland_stability.varieties import P1xP1, hirzebruch

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

RANK_MAX = 12
SAT = 8000

# rank <= 12, one twist-normalized representative per class, rank >= 2
PINNED = {0: 20, 1: 15}


_ENUM_MEMO = {}


def _classes(e, rank_cap=None):
    """Memoized so the walk is paid for once per (surface, cap), not per test."""
    key = (e, rank_cap)
    if key not in _ENUM_MEMO:
        _ENUM_MEMO[key] = exceptional_classes(e, RANK_MAX, rank_cap=rank_cap,
                                              sat_patience=SAT)
    return _ENUM_MEMO[key]


def _base(e):
    return P1xP1 if e == 0 else hirzebruch(1)


def _production_classes(e, rank_max):
    """The E15-M2 sweep's own enumeration filter, over ``c1`` in ``[0,r)^2``."""
    base = _base(e)
    antic = surface_with_index(e, F(2 - e, 2))
    out = set()
    for r in range(2, rank_max + 1):
        for x in range(r):
            for y in range(r):
                if not is_potentially_exceptional(r, (x, y), base):
                    continue
                if not is_stable_exceptional(r, (x, y), antic):
                    continue
                out.add((r, (x, y)))
    return out


# --------------------------------------------------------------------------- #
# 1. Charter: the oracle is import-independent and exact.                      #
# --------------------------------------------------------------------------- #
def _oracle_src():
    path = os.path.join(_TESTS_DIR, "oracle", "delpezzo_mutation_reference.py")
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def test_oracle_has_no_package_import():
    for line in _oracle_src().splitlines():
        stripped = line.strip()
        if stripped.startswith(("import ", "from ")):
            assert "bridgeland" not in stripped, stripped


def test_oracle_uses_no_float():
    src = _oracle_src()
    assert "import math" not in src
    for line in src.splitlines():
        stripped = line.strip()
        if stripped.startswith(("import ", "from ")):
            assert "math" not in stripped, stripped
    assert "float(" not in src


# --------------------------------------------------------------------------- #
# 2. The re-derived Riemann-Roch agrees with the package's chi.                #
# --------------------------------------------------------------------------- #
def test_reimplemented_euler_form_matches_the_package():
    # genuine K-theory classes only: ch2 = <c1,c1>/2 - n, n in Z
    for e in (0, 1, 2, 3):
        surf = P1xP1 if e == 0 else hirzebruch(e)
        for c1 in ((0, 0), (1, 0), (0, 1), (2, -3), (-1, 4)):
            for c12 in ((1, 1), (0, 2), (3, -1), (-2, -1)):
                for r, r2, n, n2 in ((1, 1, 0, 0), (2, 3, 1, -2), (5, 4, -3, 7)):
                    A = (r, c1, F(surf.lattice.self_pairing(c1), 2) - n)
                    B = (r2, c12, F(surf.lattice.self_pairing(c12), 2) - n2)
                    assert euler_pairing(A, B, e) == pkg_chi(
                        SurfaceBundle(*A), SurfaceBundle(*B), surf), (e, A, B)


def test_seed_collection_is_numerically_exceptional():
    # chi(E_i,E_i) = 1 and chi(E_j,E_i) = 0 for j > i -- by the oracle's own form
    for e in (0, 1):
        coll = standard_collection(e)
        for i, Ei in enumerate(coll):
            assert euler_pairing(Ei, Ei, e) == 1
            for j in range(i + 1, len(coll)):
                assert euler_pairing(coll[j], Ei, e) == 0, (e, i, j)


def test_canonical_class_has_square_eight():
    # K^2 = 8 is the degree-8 del Pezzo value, and pins the re-derived Gram:
    # a wrong E^2 or F.E would break it.  Holds on every F_e, not just 0/1.
    for e in range(0, 7):
        assert _self_pairing(_canonical(e), e) == 8, e
    # and O is exceptional under the oracle's own form
    for e in (0, 1):
        assert euler_pairing(line_bundle((0, 0), e), line_bundle((0, 0), e), e) == 1


# --------------------------------------------------------------------------- #
# 3. The differential gate against production's scan.                          #
# --------------------------------------------------------------------------- #
def test_oracle_matches_production_enumeration():
    for e in (0, 1):
        oracle, info = _classes(e)
        assert info["status"] == "SATURATED", info
        assert set(oracle) == _production_classes(e, RANK_MAX), e


def test_pinned_class_counts():
    # a joint drift of both enumerators would still be caught here
    for e in (0, 1):
        oracle, _ = _classes(e)
        assert len(oracle) == PINNED[e], (e, len(oracle))


def test_result_is_invariant_under_the_intermediate_rank_cap():
    # the full-scale analogue (caps 160 vs 220 at rank 130) is CORRECTIONS 26.4
    for e in (0, 1):
        low, _ = _classes(e)                       # default cap = 3 * RANK_MAX
        high, _ = _classes(e, rank_cap=5 * RANK_MAX)
        assert low == high, e


def test_f0_ranks_are_odd_and_f1_admits_rank_two():
    # CLAUDE.md invariant 4, reproduced by the independent generator
    f0, _ = _classes(0)
    assert all(r % 2 == 1 for r, _ in f0)
    f1, _ = _classes(1)
    assert (2, (1, 1)) in f1
