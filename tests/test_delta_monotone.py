"""E14-M4: cor-deltaMonotone / cor-deltaMonotoneHigher + cross-theorem sweeps.

The qualitative theorems of arXiv:1907.06739 as executable gates over the E14
machinery: monotonicity of ``delta_m^{mu-s}`` in the polarization on both sides
of the anticanonical index (``cor-deltaMonotone``), its one-sided ``e >= 2``
analogue via the reduction (``cor-deltaMonotoneHigher``), monotone transport of
per-character mu-stable existence (``cor-KstabilityEasy``), and the region-R
strictness ``DLP^{<r} < delta`` seeded by the Kronecker construction (the Sec. 8
closing corollary).  Append-only: part of the standing differential suite.
"""

from fractions import Fraction as F

from bridgeland_stability.delta_sharp import (
    delta_kronecker,
    delta_mu_stable,
    mu_stable_exists,
    surface_with_index,
)
from bridgeland_stability.dlp_hirzebruch import dlp_restricted
from bridgeland_stability.stability_interval import stability_interval
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)
F2 = hirzebruch(2)

NU_F0 = (F(1, 3), F(1, 5))       # Ex. KroneckerF0 slope, package (f, s)
NU_F1 = (F(6, 13), F(3, 13))     # Ex. KroneckerF1 slope


def test_existence_transports_toward_the_anticanonical_index():
    # cor-KstabilityEasy: mu-stable existence at m' transports to every m
    # between m0 and m'.  The F_0 Kronecker class one lattice step above the
    # wall exists at m' = 25/9 (E14-M1 pin), hence at every sampled m in
    # [1, 25/9]; the wall class itself does not exist at 25/9, hence (same
    # transport, contrapositive) at no m'' >= 25/9.
    S = lambda m: surface_with_index(0, m)
    assert mu_stable_exists(15, NU_F0, F(2, 3), S(F(25, 9))) is True
    for m in (F(2), F(3, 2), F(1)):
        assert mu_stable_exists(15, NU_F0, F(2, 3), S(m)) is True
    assert mu_stable_exists(15, NU_F0, F(3, 5), S(F(25, 9))) is False
    for m in (F(3), F(7, 2)):
        assert mu_stable_exists(15, NU_F0, F(3, 5), S(m)) is False


def test_left_side_monotonicity_on_f1():
    # cor-deltaMonotone, the 0 < m' <= m <= m0 clause, witnessed through the
    # exceptional branch: (2,(1,1)) on F_1 has I_V = (0, 1) (E14-M3 table), so
    # existence at any m' < m0 = 1/2 transports to every m in [m', 1/2].
    nu = (F(1, 2), F(1, 2))
    for m in (F(1, 8), F(1, 4), F(3, 8), F(1, 2)):
        assert mu_stable_exists(2, nu, F(3, 8), surface_with_index(1, m)) is True
    # and the interval route agrees
    iv = stability_interval(2, (1, 1), F1)
    assert iv.contains(F(1, 8)) and iv.contains(F(1, 2)) and not iv.contains(F(3, 2))


def test_sandwich_bounds_are_monotone_in_m():
    # For m >= m0, per-character existence shrinks as m grows
    # (cor-KstabilityEasy), so at a FIXED scanned rank set the certified upper
    # bound is nondecreasing in m; the DLP-based lower bound is nondecreasing by
    # prop-DLPmonotone.  Two theorem-backed sweeps on each surface.
    f0 = [delta_mu_stable(NU_F0, m, F0, max_rank=15)
          for m in (F(5, 2), F(8, 3), F(25, 9), F(14, 5))]
    f1 = [delta_mu_stable(NU_F1, m, F1, max_rank=13)
          for m in (F(3, 2), F(5, 3), F(12, 7), F(9, 5))]
    for seq in (f0, f1):
        for a, b in zip(seq, seq[1:]):
            assert a.lower <= b.lower
            assert a.upper is not None and b.upper is not None
            assert a.upper <= b.upper


def test_delta_monotone_higher_via_the_reduction():
    # cor-deltaMonotoneHigher (e >= 2, 0 < m < m'): sampled through the
    # transported Kronecker formula on F_2 -- strictly increasing here -- and
    # each value equals the F_0 value at m + 1 (the transport identity).
    ms = (F(3, 2), F(29, 18), F(16, 9), F(9, 5))
    vals = [delta_kronecker((F(8, 15), F(1, 5)), m, F2) for m in ms]
    assert all(v is not None for v in vals)
    assert all(a < b for a, b in zip(vals, vals[1:]))
    for m, v in zip(ms, vals):
        assert v == delta_kronecker(NU_F0, m + 1, F0)


def test_region_r_strictness_dlp_below_delta():
    # The Sec. 8 closing corollary ("Kronecker beats exceptional", the seed of
    # the E15 Sec. 1.5 conjecture): on the Kronecker triangle the restricted
    # DLP function sits STRICTLY below delta_m^{mu-s}.  <= is cor-deltaDLPe;
    # the strictness at these rational m is asserted as observed fact.
    for e, nu, r, pairs in (
            (0, NU_F0, 15, (F(5, 2), F(8, 3), F(25, 9), F(14, 5))),
            (1, NU_F1, 13, (F(3, 2), F(5, 3), F(12, 7), F(9, 5)))):
        surf = F0 if e == 0 else F1
        for m in pairs:
            dk = delta_kronecker(nu, m, surf)
            assert dk is not None
            assert dlp_restricted(nu, surface_with_index(e, m), r) < dk


def test_interval_endpoints_respect_delta_monotonicity():
    # Cross-theorem gate: an exceptional bundle V is mu-stable exactly on I_V,
    # and delta-monotonicity away from m0 means once V destabilizes it stays
    # destabilized: is_stable_exceptional-membership must be monotone along m
    # on each side.  Sampled through the E14-M3 interval of (11,(4,4))/F_0.
    iv = stability_interval(11, (4, 4), F0)                   # (4/7, 7/4)
    inside = (F(3, 5), F(1), F(3, 2))
    outside_right = (F(7, 4), F(2), F(3))
    outside_left = (F(1, 2), F(4, 7), F(1, 3))
    assert all(iv.contains(m) for m in inside)
    assert not any(iv.contains(m) for m in outside_right)
    assert not any(iv.contains(m) for m in outside_left)
