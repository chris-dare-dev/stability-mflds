"""E15-M1: existence obstructions for exceptional bundles on ``F_e``.

Pins the necessary-condition battery for the Sec. 11 conjecture of
arXiv:1907.06739: the paper's F_4 example is refuted by the prioritary index
exactly as in the paper; provably-existing bundles (E14-M3 table rows) are
NEVER refuted (the soundness controls); and the paper's first open case
``v_107`` on ``F_3`` is recorded with whatever the battery honestly decides.
"""

from fractions import Fraction as F

import pytest

from bridgeland_stability.exceptional_existence import (
    ExceptionalRefutation,
    exceptional_refutation,
)
from bridgeland_stability.prioritary import generic_prioritary_index
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)
F3 = hirzebruch(3)
F4 = hirzebruch(4)


def test_f4_example_refuted_by_the_prioritary_index():
    # The paper's own dispatch of (3, 1/3 E + F, 4/9) on F_4: rho_gen = 1 < 2
    # (an exceptional bundle is simple, hence H_2-prioritary by lem-simple, and
    # rigid, hence generic in the irreducible stack -- prop-excPrior /
    # cor-prioritaryRho).
    res = exceptional_refutation(3, (3, 1), F4)
    assert res.refuted is True
    assert res.rho_gen == 1
    assert "rho_gen = 1" in res.reason


def test_rho_gen_probes_pin_the_spec_findings():
    # Spec-time probes (E15 spec, 2026-07-16): the F_4 control matches the
    # paper bit-for-bit, and the easy route is INCONCLUSIVE for the paper's
    # first open case v_107 -- rho_gen = 2.
    assert generic_prioritary_index((F(1), F(1, 3)), F(4, 9), F4) == 1
    assert generic_prioritary_index(
        (F(76, 107), F(25, 107)), F(5724, 11449), F3) == 2


@pytest.mark.parametrize("surf,r,c1", [
    (F1, 2, (1, 1)),
    (F1, 11, (5, 3)),
    (F0, 3, (1, 1)),
    (F0, 11, (4, 4)),
])
def test_existing_bundles_are_never_refuted(surf, r, c1):
    # Soundness: these characters carry exceptional bundles (E14-M3 pinned
    # their nonempty stability intervals against the paper's tables), so every
    # necessary condition must PASS.  A refutation here falsifies the
    # derivation, not the bundle.
    res = exceptional_refutation(r, c1, surf)
    assert res.refuted is False
    assert res.rho_gen >= 2


def test_invalid_character_is_trivially_refuted():
    res = exceptional_refutation(2, (1, 1), F0)   # c2 not integral
    assert res.refuted is True
    assert "not potentially exceptional" in res.reason


def test_anchor_regime_guard():
    # The rigid-factor sampling is only justified where ceil(m) <= 2
    # (lem-simple gives H_1/H_2-prioritariness of a hypothetical V; beyond
    # that the sheaf need not sit in the filtered stack).
    with pytest.raises(ValueError):
        exceptional_refutation(11, (5, 3), F1, anchors=(F(5, 2),))
