"""E15-M2: spot checks of the Conjecture B range-extension sweep.

The full sweep (``scripts/e15_m2_sweep.py``) runs offline; its frozen ledger is
recorded in ``docs/CORRECTIONS.md`` Sec. 22.  These tests pin the sweep's LOGIC
on hand-verified rows: the lift arithmetic, the empty-lift criterion
``hi(I_W) <= k``, a dispatched example (the paper's F_4 example arises exactly
as the k = 2 lift of the F_0 row (3,(1,1))), and the survivor row v_107.
"""

from fractions import Fraction as F

from bridgeland_stability.prioritary import generic_prioritary_index
from bridgeland_stability.stability_interval import stability_interval
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)
F3 = hirzebruch(3)
F4 = hirzebruch(4)


def _lift(c1, k):
    """pi^{-k} in package coordinates: (x, y) -> (x + k*y, y)."""
    return (c1[0] + k * c1[1], c1[1])


def test_the_f4_example_is_the_k2_lift_of_an_f0_table_row():
    # F_0 (3,(1,1)) has I = (1/2, 2): hi = 2 > 1 (nonempty on F_2, E14-M3
    # pinned (0,1) there) but hi <= 2 -- EMPTY on F_4.  Its k = 2 lift is
    # exactly the paper's F_4 example (3, (3,1)), dispatched by rho_gen = 1.
    iv = stability_interval(3, (1, 1), F0)
    assert iv.hi == 2
    assert not (iv.hi <= 1)                       # not empty on F_2
    assert iv.hi <= 2                             # empty on F_4
    assert _lift((1, 1), 2) == (3, 1)
    nu = (F(1), F(1, 3))
    delta = F(1, 2) * (1 - F(1, 9))
    assert delta == F(4, 9)                       # the paper's discriminant
    assert generic_prioritary_index(nu, delta, F4) == 1   # dispatched


def test_a_rank_two_empty_lift_is_dispatched():
    # F_1 (2,(1,1)) has I = (0, 1): hi <= 1, so its k = 1 lift (2,(2,1)) on
    # F_3 is an empty-interval would-be counterexample -- dispatched by
    # rho_gen = 1 (consistent with the paper's claim that the FIRST potential
    # counterexample has rank 107).
    iv = stability_interval(2, (1, 1), F1)
    assert iv.hi == 1 and iv.hi <= 1
    assert _lift((1, 1), 1) == (2, 1)
    assert generic_prioritary_index(
        (F(1), F(1, 2)), F(3, 8), F3) == 1


def test_v107_is_the_survivor_row():
    # F_1 (107,(51,25)) has hi = 13/23 <= 1 (E14-M3 / E15 spec probe), its
    # k = 1 lift is v_107 = (107, (76,25)) on F_3, and rho_gen = 2 (E15-M1
    # pin): the sweep must place exactly this row on the survivor ledger.
    assert _lift((51, 25), 1) == (76, 25)
    assert generic_prioritary_index(
        (F(76, 107), F(25, 107)), F(5724, 11449), F3) == 2


def test_lift_arithmetic_is_pi_inverse():
    # pi: (x, y) -> (x - y, y), so the k-fold lift must invert it.
    from bridgeland_stability.reduction import pi_c1
    for c1 in ((51, 25), (1, 1), (4, 7), (0, 3)):
        for k in (1, 2, 3):
            lifted = _lift(c1, k)
            back = lifted
            for _ in range(k):
                back = tuple(pi_c1(back))
            assert back == c1
