"""Exact provenance pins for the E16 v_107 numerical mutation orbit."""

from fractions import Fraction as F

from scripts.e16_m1_orbit import (
    POSITIVE_PATH,
    chi_k,
    follow_path,
    is_target,
    normalize,
    standard_collection,
    twist,
)


def test_five_step_positive_path_and_euler_gram():
    final = follow_path(standard_collection(), POSITIVE_PATH)
    assert final == (
        (5, (4, 2), F(-2)),
        (6, (4, 1), F(-5, 2)),
        (107, (76, 25), F(-89, 2)),
        (4, (3, 1), F(-3, 2)),
    )
    assert [w for w in final if is_target(w)] == [final[2]]
    assert [[chi_k(a, b) for b in final] for a in final] == [
        [1, -1, -7, 1],
        [0, 1, 13, 2],
        [0, 0, 1, 27],
        [0, 0, 0, 1],
    ]


def test_normalization_is_twist_invariant_for_signed_and_rank_zero_inputs():
    start = standard_collection()
    signed = tuple((-r, (-c1[0], -c1[1]), -ch2) for r, c1, ch2 in start)
    huge = (1234, -987)
    assert normalize(tuple(twist(w, huge) for w in signed)) == normalize(signed)

    rank_zero_first = ((0, (1, 0), F(0)),) + start[1:]
    assert normalize(tuple(twist(w, huge) for w in rank_zero_first)) == normalize(
        rank_zero_first)


def test_target_recognizes_every_common_twist_and_the_dual_residue():
    target = (107, (76, 25), F(-89, 2))
    assert all(is_target(twist(target, D)) for D in ((0, 0), (7, -11), (-20, 31)))
    dual = (107, (-76, -25), F(-89, 2))
    assert is_target(dual)
