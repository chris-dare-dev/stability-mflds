"""E10-M1 / G16: lazy M2 oracle skeleton + graceful-skip guard.

E10-M1 adds NO mathematics: every assertion here is structural (import purity,
guard-raises, discovery order), not a literature value.  It confirms the oracle
mirrors the viz lazy-backend culture -- ``import bridgeland_stability`` stays
zero-dependency and the core never pulls the oracle, while calling an oracle
function with M2 absent raises a clear ``require_m2()`` guard.
"""
import subprocess
import sys
from fractions import Fraction

import pytest

from bridgeland_stability import mukai
from bridgeland_stability.exceptional import Bundle, chi as chi_rr
from bridgeland_stability.oracle import m2 as m2mod
from bridgeland_stability.oracle.m2 import (
    M2NotFoundError,
    M2Session,
    chi_via_ext,
    ext_dims,
    find_m2,
    moduli_nonempty_by_construction,
    require_m2,
    _qq_roundtrip,
    _rank1_ideal_length,
)
from bridgeland_stability.varieties import K3, P2

#: E10-M2 tests need a real Macaulay2 install; skipped (not failed) when absent
#: so the suite stays green on an M2-less host (this one) while still exercising
#: the sheaf-level Euler-pairing cross-check wherever a user provisions M2.
requires_m2 = pytest.mark.skipif(
    find_m2() is None, reason="Macaulay2 (M2) not installed"
)


def test_zero_dep_import():
    # core imports with ZERO third-party deps AND does not pull the oracle;
    # oracle subpackage itself is stdlib-only at import time (M2 spawned lazily).
    child = r'''
import sys
BLOCK = {"matplotlib","plotly","numpy","scipy","sympy","pandas","pexpect"}
class _Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in BLOCK:
            raise ImportError("blocked for zero-dep test: " + name)
        return None
sys.meta_path.insert(0, _Blocker())
import bridgeland_stability
assert "bridgeland_stability.oracle" not in sys.modules, "core must not import oracle"
import bridgeland_stability.oracle
import bridgeland_stability.oracle.m2
for name in BLOCK:
    assert name not in sys.modules, name + " imported at import time"
assert "bridgeland_stability.oracle.m2" in sys.modules
print("ZERO_DEP_OK")
'''
    r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "ZERO_DEP_OK" in r.stdout


def test_oracle_guard_raises(monkeypatch):
    # Force M2 "absent" regardless of host (mirrors viz require_mpl guard).
    monkeypatch.setattr(m2mod, "find_m2", lambda: None)
    with pytest.raises(M2NotFoundError) as ei:
        require_m2()
    assert "Macaulay2" in str(ei.value)
    # every oracle entry point guards on require_m2() first
    with pytest.raises(M2NotFoundError):
        m2mod.chi_via_ext(None, None, None)
    with pytest.raises(M2NotFoundError):
        m2mod.ext_dims(None, None, None)
    with pytest.raises(M2NotFoundError):
        m2mod.moduli_nonempty_by_construction(1, 0, 0, None)
    # the context manager also guards on entry (never spawns when M2 absent)
    with pytest.raises(M2NotFoundError):
        with M2Session():
            pass


def test_find_m2_env_override(monkeypatch):
    # discovery honors the documented BRIDGELAND_M2 override, no real binary needed;
    # proves require_m2() is not over-eager (present -> returns path, no guard).
    fake = "/nonexistent/path/to/M2"
    monkeypatch.setenv("BRIDGELAND_M2", fake)
    assert m2mod.find_m2() == fake
    assert require_m2() == fake


def test_oracle_public_surface():
    # facade re-exports the M1 surface; core __all__ is NOT polluted by the oracle.
    import bridgeland_stability as bs
    import bridgeland_stability.oracle as orc
    for name in ("require_m2", "M2Session", "M2NotFoundError",
                 "chi_via_ext", "ext_dims", "moduli_nonempty_by_construction"):
        assert name in orc.__all__ and hasattr(orc, name)
    assert "oracle" not in bs.__all__  # oracle stays out of the core surface


# ---------------------------------------------------------------------------
# E10-M2: sheaf-level Euler-pairing cross-check (skipped when M2 is absent).
#
# Every asserted value is literature-anchored AND recomputed against the shipped
# formula layer (CLAUDE.md two-way verification).  No computed core value changes
# -- these tests DIFF the M2 sheaf-level answer against exceptional.chi /
# mukai.pairing, which stay byte-for-byte as shipped.
# ---------------------------------------------------------------------------
@requires_m2
def test_chi_matches_formula():
    # chi(O, O(n)) = (n+1)(n+2)/2 = Hilbert polynomial of O_{P^2}; the M2
    # sheaf-level Euler pairing must reproduce the shipped exceptional.chi AND
    # the closed form for every n in {-3,...,5}.  Sequence: 1,0,0,1,3,6,10,15,21.
    expected = [1, 0, 0, 1, 3, 6, 10, 15, 21]
    for n, want in zip(range(-3, 6), expected):
        got = chi_via_ext(Bundle.O(0), Bundle.O(n), P2)
        assert got == want, (n, got, want)
        assert got == (n + 1) * (n + 2) // 2          # closed form
        assert got == chi_rr(Bundle.O(0), Bundle.O(n))  # shipped formula layer
    # explicit pinned anchors (the '15' red herring is chi(O,O(4)), not chi(O,O(3)))
    assert chi_via_ext(Bundle.O(0), Bundle.O(3), P2) == 10
    assert chi_via_ext(Bundle.O(0), Bundle.O(4), P2) == 15
    assert chi_via_ext(Bundle.O(0), Bundle.O(1), P2) == 3
    assert chi_via_ext(Bundle.O(0), Bundle.O(0), P2) == 1


@requires_m2
def test_ext_dims_o_minus3():
    # Ext^i_{P^2}(O, O(-3)) = H^i(P^2, O(-3)): h0=0 (negative degree), h1=0 (no
    # middle cohomology for P^2 line bundles), h2=1 by Serre duality
    # H^2(O(-3)) = H^0(O(K-(-3)))^* = H^0(O(0))^* = C  (K_{P^2}=O(-3)).
    dims = ext_dims(Bundle.O(0), Bundle.O(-3), P2)
    assert dims == (0, 0, 1)
    alt = sum((-1) ** i * d for i, d in enumerate(dims))
    assert alt == 1                                    # 0 - 0 + 1 = +1
    assert alt == chi_rr(Bundle.O(0), Bundle.O(-3))    # shipped chi(O,O(-3)) = +1
    assert chi_via_ext(Bundle.O(0), Bundle.O(-3), P2) == 1


@requires_m2
def test_ext_dims_k3_structure_sheaf():
    # H^*(K3, O_X) = (1, 0, 1): h0=1, h1=0 (K3 is regular), h2 = h0(O_X)^* = 1 by
    # Serre duality with K_X = 0.  Holds for ANY K3 (the Fermat quartic used by
    # the translator is fine; Picard rank is irrelevant to H^i(O_X)).
    X = K3(4)  # smooth quartic in P^3 has H^2 = 4, Pic = ZH
    dims = ext_dims(Bundle.O(0), Bundle.O(0), X)
    assert dims == (1, 0, 1)
    alt = sum((-1) ** i * d for i, d in enumerate(dims))  # 1 - 0 + 1 = 2
    v_O = mukai.MukaiVector(1, 0, 1)                       # v(O) = (1,0,1)
    # the -ext alternating sum equals the Mukai self-pairing (pinned <v,v> = -2)
    assert -alt == mukai.self_pairing(v_O, X.d)
    assert mukai.self_pairing(v_O, X.d) == -2
    assert -alt == mukai.pairing(v_O, v_O, X.d)


@requires_m2
def test_roundtrip():
    # Fraction -> QQ literal -> toString text -> Fraction is lossless for any
    # rational (G16: M2 rationals are QQ; Fraction("p/q") <-> p/q exact).
    assert _qq_roundtrip(Fraction(-7, 6)) == Fraction(-7, 6)
    assert _qq_roundtrip(Fraction(4)) == Fraction(4)  # whole number prints "4"


# ---------------------------------------------------------------------------
# E10-M3: constructive SUFFICIENT-ONLY non-emptiness witness (G16).
#
# The load-bearing arithmetic (l = c1^2/2 - ch2 = c2) runs on-host with no M2.
# The full sheaf-level construction is [UNVERIFIED] (M2 absent), so its value
# test is @requires_m2; the two monkeypatched tests exercise the surrounding
# gate/parse logic on-host with a canned Euler-characteristic transcript.
#
# Anchor (CORRECTIONS.md P^2[2]): ch(I_Z(c1)) = (1, c1, c1^2/2 - l).  For
# (1,0,-2): l = 0 - (-2) = 2 -> the ideal sheaf of 2 points, moduli NONEMPTY.
# ---------------------------------------------------------------------------
def test_rank1_ideal_length():
    # l = c1^2/2 - ch2 = c2 (non-negative integer) is the constructibility gate.
    assert _rank1_ideal_length(1, 0, Fraction(-2)) == 2   # P^2[2] ideal sheaf
    assert _rank1_ideal_length(1, 2, Fraction(1)) == 1
    assert _rank1_ideal_length(1, 0, Fraction(0)) == 0    # I_empty = O
    # not constructible -> None (NOT an emptiness claim):
    assert _rank1_ideal_length(1, 0, Fraction(1)) is None     # l = -1 < 0
    assert _rank1_ideal_length(1, 1, Fraction(0)) is None     # l = 1/2 not in Z
    assert _rank1_ideal_length(1, 0, Fraction(-1, 2)) is None  # l = 1/2 not in Z
    assert _rank1_ideal_length(2, 0, Fraction(-3)) is None     # rank != 1


def test_moduli_none_when_not_constructible(monkeypatch):
    # require_m2() passes (fake path) but no subprocess is spawned: the length
    # gate returns None first.  A non-constructible class -> None, never False.
    monkeypatch.setenv("BRIDGELAND_M2", "/nonexistent/path/to/M2")
    for r, c1, ch2 in [(1, 0, Fraction(1)), (2, 0, Fraction(-3)), (1, 1, Fraction(0))]:
        res = moduli_nonempty_by_construction(r, c1, ch2, P2)
        assert res is None
        assert res is not False


def test_moduli_witness_success_and_mismatch(monkeypatch):
    # require_m2() passes via fake path; _run_m2 is stubbed so no real M2 runs.
    monkeypatch.setenv("BRIDGELAND_M2", "/nonexistent/path/to/M2")
    # Canned Euler transcript for (r,c1,ch2)=(1,0,-2) via P^2 Riemann-Roch
    # chi(E(t)) = (r/2)t^2 + (c1+3r/2)t + (ch2+3/2 c1 + r):
    #   chi(0) = -2 + 0 + 1 = -1;  chi(1) = 1/2 + 3/2 - 1 = 1;
    #   chi(2) = 2 + 3 - 1 = 4.  _parse_witness must invert to (1, 0, -2).
    success = "CHI 0 -1\nCHI 1 1\nCHI 2 4\nOK\n"
    monkeypatch.setattr(m2mod, "_run_m2", lambda code, **kw: success)
    res = moduli_nonempty_by_construction(1, 0, Fraction(-2), P2)
    assert res is True
    # A mismatched transcript (reconstructs to (1,2,1)) -> None, never False.
    #   (1,2,1): chi(0)=1+3+1=5; chi(1)=1/2+(2+3/2)+5=9; chi(2)=2+7+5=14.
    mismatch = "CHI 0 5\nCHI 1 9\nCHI 2 14\nOK\n"
    monkeypatch.setattr(m2mod, "_run_m2", lambda code, **kw: mismatch)
    res2 = moduli_nonempty_by_construction(1, 0, Fraction(-2), P2)
    assert res2 is None
    assert res2 is not False


@requires_m2
def test_moduli_witness_is_sufficient_only():
    # The milestone's named test: a real M2 builds I_Z(c1) and verifies its
    # Chern character.  A known-nonempty class -> True; a construction failure
    # -> None (NOT False) -- a failed build is not a proof of emptiness (G16).
    assert moduli_nonempty_by_construction(1, 0, Fraction(-2), P2) is True  # P^2[2]
    r_neg = moduli_nonempty_by_construction(1, 0, Fraction(1), P2)   # l=-1: not built
    assert r_neg is None and r_neg is not False
    r_rk2 = moduli_nonempty_by_construction(2, 0, Fraction(-3), P2)  # rank!=1 witness
    assert r_rk2 is None and r_rk2 is not False
