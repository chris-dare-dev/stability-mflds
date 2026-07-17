"""E14-M3: stability intervals ``I_V`` of exceptional bundles (thm-stabilityInterval).

Pins BOTH of arXiv:1907.06739's stability-interval tables bit-for-bit --
Table "stabilityInterval0" (F_0, 13 rows, ranks to 19) and Table
"stabilityInterval1" (F_1, 15 rows, ranks to 20) -- including the destabilizing
witness columns, after the paper -> package coordinate transport
(paper ``(a, b)`` = ``a E + b F`` resp. ``a F_1 + b F_2``  ->  package ``(b, a)``).
Plus: the two-route boundary differential against ``is_stable_exceptional``, the
``e >= 2`` transport (``cor-highermus``), the paper's F_4 no-stable-bundle
example, and the Sec. 11 conjecture's first candidate on F_3.
"""

from fractions import Fraction as F

import pytest

from bridgeland_stability.delta_sharp import surface_with_index
from bridgeland_stability.dlp_hirzebruch import is_stable_exceptional
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.stability_interval import StabilityInterval, stability_interval
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)
F2 = hirzebruch(2)
F3 = hirzebruch(3)
F4 = hirzebruch(4)

INF = None

# Table "stabilityInterval0" (F_0; paper c1 = a F_1 + b F_2, normalized
# 0 <= a < r/2, a <= b < r) -- rows as printed: (r, (a,b), lo, hi, W0, W1).
TABLE_F0 = [
    (1, (0, 0), 0, INF, None, None),
    (3, (1, 1), F(1, 2), 2, (1, (-1, 1)), (1, (1, -1))),
    (5, (1, 2), F(1, 2), 3, (1, (-1, 1)), (1, (1, -2))),
    (7, (1, 3), F(1, 2), 4, (1, (-1, 1)), (1, (1, -3))),
    (9, (1, 4), F(1, 2), 5, (1, (-1, 1)), (1, (1, -4))),
    (11, (1, 5), F(1, 2), 6, (1, (-1, 1)), (1, (1, -5))),
    (11, (4, 4), F(4, 7), F(7, 4), (5, (-2, 4)), (5, (4, -2))),
    (13, (1, 6), F(1, 2), 7, (1, (-1, 1)), (1, (1, -6))),
    (15, (1, 7), F(1, 2), 8, (1, (-1, 1)), (1, (1, -7))),
    (17, (1, 8), F(1, 2), 9, (1, (-1, 1)), (1, (1, -8))),
    (17, (5, 5), F(8, 9), F(9, 8), (7, (1, 3)), (7, (3, 1))),
    (19, (1, 9), F(1, 2), 10, (1, (-1, 1)), (1, (1, -9))),
    (19, (4, 7), F(8, 9), 3, (7, (1, 3)), (1, (1, -2))),
]

# Table "stabilityInterval1" (F_1; paper c1 = a E + b F, 0 <= a <= r/2,
# 0 <= b < r).
TABLE_F1 = [
    (1, (0, 0), 0, INF, None, None),
    (2, (1, 1), 0, 1, None, (1, (1, 0))),
    (4, (1, 2), 0, 2, None, (1, (1, -1))),
    (5, (2, 2), 0, F(2, 3), None, (1, (1, 0))),
    (6, (1, 3), 0, 3, None, (1, (1, -2))),
    (8, (1, 4), 0, 4, None, (1, (1, -3))),
    (10, (1, 5), 0, 5, None, (1, (1, -4))),
    (11, (3, 5), F(3, 7), 2, (6, (1, 3)), (1, (1, -1))),
    (12, (1, 6), 0, 6, None, (1, (1, -5))),
    (13, (5, 5), 0, F(5, 8), None, (1, (1, 0))),
    (14, (1, 7), 0, 7, None, (1, (1, -6))),
    (16, (1, 8), 0, 8, None, (1, (1, -7))),
    (18, (1, 9), 0, 9, None, (1, (1, -8))),
    (19, (5, 10), F(1, 9), F(9, 5), (5, (-2, 3)), (6, (5, -3))),
    (20, (1, 10), 0, 10, None, (1, (1, -9))),
]


def _pkg(c):
    """paper (a, b) -> package (f, s) = (b, a)."""
    return (c[1], c[0])


def _wpkg(w):
    return None if w is None else (w[0], _pkg(w[1]))


@pytest.mark.parametrize("row", TABLE_F0, ids=lambda r: f"F0-{r[0]}-{r[1]}")
def test_table_stability_interval_f0(row):
    r, c1, lo, hi, w0, w1 = row
    res = stability_interval(r, _pkg(c1), F0)
    assert not res.empty
    assert res.lo == F(lo)
    assert res.hi == (None if hi is None else F(hi))
    assert res.witness_lo == _wpkg(w0)
    assert res.witness_hi == _wpkg(w1)
    assert res.certificate.rigor == Rigor.PROVEN


@pytest.mark.parametrize("row", TABLE_F1, ids=lambda r: f"F1-{r[0]}-{r[1]}")
def test_table_stability_interval_f1(row):
    r, c1, lo, hi, w0, w1 = row
    res = stability_interval(r, _pkg(c1), F1)
    assert not res.empty
    assert res.lo == F(lo)
    assert res.hi == (None if hi is None else F(hi))
    assert res.witness_lo == _wpkg(w0)
    assert res.witness_hi == _wpkg(w1)


def test_anticanonical_index_is_always_inside():
    # Gorodentsev: every exceptional bundle on a del Pezzo F_e is
    # mu_{-K}-stable, so m0 = 1 - e/2 lies in every nonempty interval.
    for surf, e, rows in ((F0, 0, TABLE_F0), (F1, 1, TABLE_F1)):
        m0 = F(2 - e, 2)
        for r, c1, *_ in rows:
            assert stability_interval(r, _pkg(c1), surf).contains(m0)


def test_two_route_boundary_differential():
    # I_v = {m : Delta >= DLP^{<r}_{H_m}(nu)} (cor-DLPexcdelPezzo), so interval
    # membership must agree with is_stable_exceptional at sampled m -- inside,
    # outside, and AT the open endpoints (V is strictly semistable there).
    for surf, e, r, c1_paper, samples in (
            (F0, 0, 11, (4, 4), (F(1, 2), F(4, 7), 1, F(7, 4), 2)),
            (F1, 1, 11, (3, 5), (F(1, 3), F(3, 7), F(1, 2), 1, 2, F(5, 2))),
    ):
        c1 = _pkg(c1_paper)
        res = stability_interval(r, c1, surf)
        for m in samples:
            assert res.contains(m) == bool(
                is_stable_exceptional(r, c1, surface_with_index(e, m))), (r, c1, m)


def test_fiber_swap_symmetry_on_f0():
    # Swapping the two rulings of F_0 maps H_m to the ray of H_{1/m}: the
    # interval of the swapped character is the reciprocal interval.
    res = stability_interval(5, _pkg((1, 2)), F0)             # (1/2, 3)
    swapped = stability_interval(5, _pkg((2, 1)), F0)
    assert (swapped.lo, swapped.hi) == (1 / res.hi, 1 / res.lo)


def test_rank_one_is_everywhere_stable():
    for surf in (F0, F1, F2, F3):
        res = stability_interval(1, (0, 0), surf)
        assert (res.lo, res.hi, res.empty) == (0, None, False)


def test_non_exceptional_characters_are_refused():
    with pytest.raises(ValueError):
        stability_interval(2, (1, 1), F0)        # c2 non-integral: not exceptional
    with pytest.raises(ValueError):
        stability_interval(F(3, 2), (1, 1), F0)  # never truncate the rank
    with pytest.raises(ValueError):
        stability_interval(3, (F(1, 2), 1), F0)  # non-integral c1


def test_e_ge_2_transport_and_the_f4_example():
    # cor-highermus: I_v on F_e = {t > 0 : t + 1 in I_{pi(v)}}.  The paper's F_4
    # example (3, 1/3 E + F, 4/9) reduces to the F_0 row (3,(1,1)) with
    # I = (1/2, 2): on F_2 the window is (0, 1), and on F_4 it closes -- there
    # are NO slope-stable sheaves of that character for any polarization.
    res2 = stability_interval(3, (2, 1), F2)                  # pi of the F_4 char
    assert (res2.lo, res2.hi, res2.empty) == (0, 1, False)
    res4 = stability_interval(3, (3, 1), F4)
    assert res4.empty
    # transported right witness on F_2: pi^{-1} of the F_0 witness -- paper
    # W_1 = (1, (1,-1)) is package (-1, 1), lifting to (0, 1)
    assert res2.witness_hi == (1, (0, 1))


def test_conjecture_candidate_on_f3_has_empty_interval():
    # Sec. 11 conjecture, first potential counterexample: (107, 25/107 E +
    # 76/107 F, 5724/11449) on F_3 -- "the stability interval is empty"
    # (paper, verbatim).  Package c1 = (76, 25); reduces to (107, (51, 25)) on
    # F_1, whose interval must sit entirely at or below the transport threshold.
    res = stability_interval(107, (76, 25), F3)
    assert res.empty
