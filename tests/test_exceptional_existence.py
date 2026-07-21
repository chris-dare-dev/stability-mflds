"""E15-M1: existence obstructions for exceptional bundles on ``F_e``.

Pins the necessary-condition battery for the Sec. 11 conjecture of
arXiv:1907.06739: the paper's F_4 example is refuted by the prioritary index
exactly as in the paper; provably-existing bundles (E14-M3 table rows) are
NEVER refuted (the soundness controls); and the paper's first open case
``v_107`` on ``F_3`` is recorded with whatever the battery honestly decides.
"""

from fractions import Fraction as F
from types import SimpleNamespace

import pytest

from bridgeland_stability.exceptional_existence import (
    ExceptionalRefutation,
    exceptional_refutation,
)
from bridgeland_stability.prioritary import generic_prioritary_index
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)
F2 = hirzebruch(2)
F3 = hirzebruch(3)
F4 = hirzebruch(4)


def test_f4_example_refuted_by_the_battery():
    # The paper's own dispatch of (3, 1/3 E + F, 4/9) on F_4.  Since E15-M1e
    # the battery's FIRST firing condition is the Gaeta star (at D' = (6,1) =
    # H_2 on F_4 -- reproducing rho_gen = 1 < 2 through shared machinery; the
    # test_rho_gen_probes still pins directly); the verdict is unchanged.
    res = exceptional_refutation(3, (3, 1), F4)
    assert res.refuted is True
    assert "Gaeta star" in res.reason
    from bridgeland_stability.exceptional_existence import gaeta_star_conditions
    assert [(D, l, r) for (D, l, r) in gaeta_star_conditions(3, (3, 1), F4)
            if l > r] == [((6, 1), 2, 0)]


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


def test_chi_box_shape_and_v107_values():
    # E15-M1d: the box on F_e has (e+3)*3 - 2 divisors; on v_107 every value
    # is deeply negative (pinned spot values) -- the family does NOT refute
    # the candidate, recorded as such.
    from bridgeland_stability.exceptional_existence import chi_box_conditions
    rows = chi_box_conditions(107, (76, 25), F3)
    assert len(rows) == 6 * 3 - 2
    vals = dict(rows)
    assert vals[(0, 1)] == -11448                 # = chi(v,v) - r^2
    assert vals[(0, 2)] == -57244
    assert vals[(4, 2)] == -11448
    assert all(x <= 0 for x in vals.values())


def test_chi_box_controls_and_the_f4_weakness_record():
    # Existing bundles satisfy every box condition (soundness); the F_4
    # example ALSO satisfies them all -- the chi-box is strictly weaker than
    # rho_gen on every known refutation case, and ships as a battery widener
    # only (the honest track record, CORRECTIONS Sec. 21).
    from bridgeland_stability.exceptional_existence import chi_box_conditions
    for surf, rr, cc in ((F1, 2, (1, 1)), (F1, 11, (5, 3)),
                         (F0, 3, (1, 1)), (F0, 11, (4, 4))):
        assert all(x <= 0 for (_, x) in chi_box_conditions(rr, cc, surf))
    assert all(x <= 0 for (_, x) in chi_box_conditions(3, (3, 1), F4))


def test_general_betti_pins():
    # thm-BN on hand-checkable line-bundle characters (F_3): h0(O) = 1;
    # h0(O(F)) = 2 (fiber sections); h2(O(K)) = 1 by Serre duality.
    from bridgeland_stability.prioritary import general_betti
    lat = F3.lattice
    assert general_betti(1, (0, 0), F(0), F3) == (1, 0, 0)
    assert general_betti(1, (1, 0), F(0), F3) == (2, 0, 0)
    K = (-5, -2)
    chK = F(1, 2) * lat.self_pairing(K)
    h = general_betti(1, K, chK, F3)
    assert h[2] == 1 and h[0] == 0
    # Rank one is not locally free when Delta > 0: I_Z(K), not a fictitious
    # dual character retaining -length(Z).  From
    # 0 -> I_Z(K) -> O(K) -> O_Z -> 0, h=(0,length(Z),1).
    for n in (1, 2, 5):
        assert general_betti(1, K, chK - n, F3) == (0, n, 1)
    K0 = (-2, -2)
    chK0 = F(1, 2) * F0.lattice.self_pairing(K0)
    for n in (1, 2, 5):
        assert general_betti(1, K0, chK0 - n, F0) == (0, n, 1)


def test_gaeta_star_machinery_cross_checks_rho_gen():
    # E15-M1e: the Gaeta exponents of v_107 are (329, 411, 7, 18) with the
    # rank identity -329+411+7+18 = 107; the H_2 box row holds (the general
    # sheaf IS H_2-prioritary, rho_gen = 2), and the H_3 delta-term Betti
    # number is 81 -- giving LHS = 18*81 = 1458 > 0, i.e. the star inequality
    # VIOLATED at H_3, reproducing rho_gen(v_107) = 2 through a second
    # computation (not independent: both share the paper's prioritary and
    # Brill--Noether chain).  H_3 is outside the
    # refutation box (-(K+H_3) is not effective), so it appears here as a
    # machinery differential only.
    from bridgeland_stability.exceptional_existence import gaeta_star_conditions
    from bridgeland_stability.prioritary import general_betti
    rows = dict(((D, (l, r)) for (D, l, r) in
                 ((D, l, r) for (D, l, r) in gaeta_star_conditions(107, (76, 25), F3))))
    assert len(rows) == 16
    assert rows[(5, 1)] == (0, 0)                 # H_2: holds
    # the H_3 delta-twist character: v(-L0 - (6,1)) with L0 = (0,1) pkg
    lat = F3.lattice
    r, c1 = 107, (76, 25)
    delta = F(1, 2) * (1 - F(1, r * r))
    nu = (F(76, 107), F(25, 107))
    ch2 = r * (F(1, 2) * lat.self_pairing(nu) - delta)
    D = (-0 - 6, -1 - 1)                          # -L0 - H_3
    c1t = (c1[0] + r * D[0], c1[1] + r * D[1])
    ch2t = ch2 + lat.pairing(c1, D) + r * F(1, 2) * lat.self_pairing(D)
    assert general_betti(r, c1t, ch2t, F3)[2] == 81


def test_gaeta_star_v107_and_controls():
    # v_107 passes the star inequality at ALL 16 box divisors (every relevant
    # twist has h^2 = 0 -- a distinct implementation check, not an independent
    # theorem chain);
    # the existing-bundle controls pass everywhere (soundness: their general
    # sheaves ARE the bundles, D-prioritary throughout the box).
    from bridgeland_stability.exceptional_existence import gaeta_star_conditions
    rows = gaeta_star_conditions(107, (76, 25), F3)
    assert all(l == r for (_, l, r) in rows)
    assert all((l, r) == (0, 0) for (D, l, r) in rows if D[1] == 2)
    for surf, rr, cc in ((F1, 11, (5, 3)), (F0, 11, (4, 4)), (F1, 2, (1, 1))):
        assert all(l == r for (_, l, r) in gaeta_star_conditions(rr, cc, surf))


def test_rigid_factor_refutation_records_the_firing_sample():
    # Exact F_2 witness: at anchor 1, g=1/1568 and m=3137/3136.  The first
    # generic HN factor has Delta=3/4 > 15/32, so the refutation fires, but the
    # evidence contract must retain the very sample that caused it.
    res = exceptional_refutation(7, (0, 2), F2)
    assert res.refuted is True
    assert "generic HN factor" in res.reason
    assert len(res.samples) == 1
    anchor, length, facts = res.samples[0]
    assert (anchor, length) == (F(1), 2)
    assert facts[0] == (4, (0, 2), F(3, 4), F(15, 32))


def test_length_one_empty_interval_is_a_certificate_conflict(monkeypatch):
    # The old branch mislabeled two contradictory PROVEN subsystems as a
    # nonexistence proof.  Force only that dormant branch and require a loud,
    # fail-closed result.  The genuine normalized sweep found no live hits.
    import bridgeland_stability.exceptional_existence as ee

    monkeypatch.setattr(ee, "is_potentially_exceptional", lambda *args: True)
    monkeypatch.setattr(ee, "chi_box_conditions", lambda *args: ())
    monkeypatch.setattr(ee, "gaeta_star_conditions", lambda *args: ())
    monkeypatch.setattr(ee, "generic_prioritary_index", lambda *args: 2)
    monkeypatch.setattr(ee, "stability_interval",
                        lambda *args: SimpleNamespace(empty=True))
    monkeypatch.setattr(ee, "_chamber_gap", lambda *args, **kwargs: F(1, 128))
    monkeypatch.setattr(ee, "surface_with_index", lambda *args: F2)
    monkeypatch.setattr(
        ee, "generic_hn_factors",
        lambda r, c1, ch2, surface: ((r, tuple(c1), ch2),))

    with pytest.raises(AssertionError, match="certificate conflict"):
        ee.exceptional_refutation(2, (1, 0), F2)


def test_anchor_regime_guard():
    # The rigid-factor sampling is only justified where ceil(m) <= 2
    # (lem-simple gives H_1/H_2-prioritariness of a hypothetical V; beyond
    # that the sheaf need not sit in the filtered stack).
    with pytest.raises(ValueError):
        exceptional_refutation(11, (5, 3), F1, anchors=(F(5, 2),))
