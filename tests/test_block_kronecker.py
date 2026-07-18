"""E15-M3: the Conjecture A falsification harness.

Positive controls: the Sec. 8 / CORRECTIONS Sec. 14-15 Kronecker pins classify
as length-2 both-non-semiexceptional and the block-decomposition search finds
the paper's own collections with the paper's own exponents.  Plus unit tests of
the exact Z-span solver and the collection tripwire.
"""

from fractions import Fraction as F

from bridgeland_stability.block_kronecker import (
    block_decomposition,
    classify_generic_filtration,
    _collection,
    _in_z_span,
    _line,
)
from bridgeland_stability.delta_sharp import surface_with_index
from bridgeland_stability.exceptional_surface import is_exceptional_collection
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)


def test_f1_kronecker_pin_classifies_and_decomposes():
    # (13,(6,3),-13/2) at m = 121/70 (the Sec. 14 sample H = (191,70)):
    # factors (2,(0,1),-3/2), (11,(6,2),-5), both non-semiexceptional; the
    # block search finds the paper's l = 3 untwisted collection with
    # v1 = 1*E3 + 1*E4 and v2 = -2*E1 + 13*E2 (ex-KroneckerF1: a = b = 1,
    # c = 2, d = 13).
    S = surface_with_index(1, F(121, 70))
    p = classify_generic_filtration(13, (6, 3), F(-13, 2), S)
    assert p.length == 2
    assert p.semiexceptional == (False, False)
    assert p.non_semiexceptional_count == 2
    assert p.conjecture_a_violation is False
    assert p.factors == ((2, (0, 1), F(-3, 2)), (11, (6, 2), F(-5)))
    w = block_decomposition(p.factors[0], p.factors[1], S)
    assert w is not None
    assert (w.l, w.twist) == (3, (0, 0))
    assert w.coeffs_v1 == (1, 1)
    assert w.coeffs_v2 == (-2, 13)


def test_f0_kronecker_pin_classifies_and_decomposes():
    # (15,(5,3),-8) just above the m = 25/9 wall: factors (2,(-1,1),-2),
    # (13,(6,2),-6) (Sec. 14); witness l = 3, untwisted, v2 = -2*E1 + 15*E2
    # (ex-KroneckerF0: c = 2, d = 15).
    S = surface_with_index(0, F(25, 9) + F(1, 100))
    p = classify_generic_filtration(15, (5, 3), F(-8), S)
    assert p.length == 2 and p.non_semiexceptional_count == 2
    assert p.factors == ((2, (-1, 1), F(-2)), (13, (6, 2), F(-6)))
    w = block_decomposition(p.factors[0], p.factors[1], S)
    assert w is not None
    assert (w.l, w.twist) == (3, (0, 0))
    assert w.coeffs_v1 == (1, 1)
    assert w.coeffs_v2 == (-2, 15)


def test_z_span_solver_is_exact():
    E3 = _line(F1, (1, 0))
    E4 = _line(F1, (-1, 1))
    # 1*E3 + 1*E4 = (2, (0,1), ch2): solvable, integral
    assert _in_z_span((2, (0, 1), E3.ch2 + E4.ch2), E3, E4) == (1, 1, True)
    # 3*E3 - 2*E4 likewise, with a negative coefficient
    v = (1, (3 - (-2), -2), 3 * E3.ch2 - 2 * E4.ch2)
    assert _in_z_span(v, E3, E4) == (3, -2, True)
    # a ch2 mismatch is rejected even when (r, c1) solves
    assert _in_z_span((2, (0, 1), E3.ch2 + E4.ch2 + 1), E3, E4) is None
    # a c1 mismatch on the checked component is rejected
    assert _in_z_span((2, (5, 1), E3.ch2 + E4.ch2), E3, E4) is None


def test_sec8_collections_pass_the_orthogonality_tripwire():
    for surf, e in ((F0, 0), (F1, 1)):
        for l in (3, 4, 6):
            for twist in ((0, 0), (2, -1), (-3, 2)):
                E = _collection(surf, e, l, twist)
                assert is_exceptional_collection(list(E), surf)


def test_region_k_single_nonsemiexceptional_is_not_a_violation():
    # Thm 1.13-shaped filtrations (at most ONE non-semiexceptional factor) are
    # conjecture-consistent regardless of length; the violation flag needs
    # length >= 3 AND >= 2 non-semiexceptional factors.
    S = surface_with_index(1, F(121, 70))
    p = classify_generic_filtration(13, (6, 3), F(-13, 2), S)
    assert p.length == 2                          # length-2 pair: never a violation
    assert p.conjecture_a_violation is False
