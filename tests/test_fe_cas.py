"""E10-M4 / G16: the F_e CAS cross-check -- Macaulay2 toric cohomology vs the RR layer.

The E13 re-audit closed the P^2 epsilon-recursion common mode with the mutation
oracle (docs/CORRECTIONS.md Sec. 13) but noted the F_e side stayed single-sourced.
This module closes that: Macaulay2's ``NormalToricVarieties`` computes the
cohomology of a toric line bundle by POLYTOPE LATTICE-POINT COMBINATORICS -- a
route entirely independent of the package's Riemann-Roch transcription -- and the
resulting Euler characteristics must match ``hilbert_P`` / the RR Euler pairing
over whole windows of line-bundle classes on ``F_0 .. F_3``.

No basis convention is trusted: the M2 transcript is SELF-DESCRIBING (it emits
the prime-divisor classes and ``-K_X`` in M2's own Cl(X) basis), and the test
FITS the unimodular identification ``T`` with the package ``(f, s)`` basis from
the data -- requiring ``T(-K_M2) = (e+2, 2)`` and full-table chi agreement --
then asserts how many identifications exist (two on ``F_0``, whose rulings are
interchangeable; exactly one for ``e >= 1``).

The live tests are ``@requires_m2`` (Macaulay2 in WSL via ``scripts/m2-wsl.cmd``;
set ``BRIDGELAND_M2`` to that shim).  The transcript parser is exercised without
M2 via a canned transcript, so the plumbing is never skip-hidden (the E13 R4
lesson).
"""

from fractions import Fraction as F
from itertools import product

import pytest

import bridgeland_stability.oracle.m2 as m2mod
from bridgeland_stability.dlp_hirzebruch import hilbert_P
from bridgeland_stability.exceptional_surface import SurfaceBundle, chi as surface_chi
from bridgeland_stability.nonemptiness_rational import hirzebruch_with_polarization
from bridgeland_stability.oracle.m2 import fe_line_bundle_cohomology, find_m2
from bridgeland_stability.varieties import P1xP1

requires_m2 = pytest.mark.skipif(
    find_m2() is None, reason="Macaulay2 (M2) not installed"
)

#: An ample F_e per index (H is irrelevant to chi; hilbert_P needs the surface).
_SURFACES = {
    0: P1xP1,
    1: hirzebruch_with_polarization(1, (3, 2)),
    2: hirzebruch_with_polarization(2, (5, 2)),
    3: hirzebruch_with_polarization(3, (7, 2)),
}


def _chi_table(table):
    return {cd: hs[0] - hs[1] + hs[2] for cd, hs in table["coh"].items()}


def _fit_identifications(e, table):
    """All unimodular T: Cl_M2 -> (f, s) with T(-K_M2) = (e+2, 2) and FULL-table
    chi agreement chi_M2(c, d) == P(T.(c, d)) (hilbert_P on the package F_e)."""
    S = _SURFACES[e]
    chi_m2 = _chi_table(table)
    aK = table["antik"]
    fits = []
    for t00, t01, t10, t11 in product(range(-4, 5), repeat=4):
        if t00 * t11 - t01 * t10 not in (1, -1):
            continue
        if (t00 * aK[0] + t01 * aK[1], t10 * aK[0] + t11 * aK[1]) != (e + 2, 2):
            continue
        ok = all(
            hilbert_P((F(t00 * c + t01 * d), F(t10 * c + t11 * d)), S) == chi
            for (c, d), chi in chi_m2.items()
        )
        if ok:
            fits.append(((t00, t01), (t10, t11)))
    return fits


# --------------------------------------------------------------------------- #
# 1. The transcript parser, exercised WITHOUT M2 (never skip-hidden).           #
# --------------------------------------------------------------------------- #
_CANNED = (
    "DEG 0 1 0\nDEG 1 -1 1\nDEG 2 1 0\nDEG 3 0 1\n"
    "ANTIK 1 2\n"
    "COH 0 0 0 1\nCOH 0 0 1 0\nCOH 0 0 2 0\n"
    "OK\n"
)


def test_parser_on_a_canned_transcript(monkeypatch):
    monkeypatch.setenv("BRIDGELAND_M2", "/nonexistent/path/to/M2")
    monkeypatch.setattr(m2mod, "_run_m2", lambda code, **kw: _CANNED)
    got = fe_line_bundle_cohomology(1, 0, 0)          # 1x1 window matches _CANNED
    assert got["degrees"] == {0: (1, 0), 1: (-1, 1), 2: (1, 0), 3: (0, 1)}
    assert got["antik"] == (1, 2)
    assert got["coh"] == {(0, 0): (1, 0, 0)}


def test_parser_rejects_a_truncated_transcript(monkeypatch):
    monkeypatch.setenv("BRIDGELAND_M2", "/nonexistent/path/to/M2")
    monkeypatch.setattr(m2mod, "_run_m2", lambda code, **kw: "DEG 0 1 0\nOK\n")
    with pytest.raises(m2mod.M2NotFoundError, match="drift"):
        fe_line_bundle_cohomology(1, 0, 0)


# --------------------------------------------------------------------------- #
# 2. The live cross-check: toric combinatorics vs the package RR layer.        #
# --------------------------------------------------------------------------- #
@requires_m2
@pytest.mark.parametrize("e", [0, 1, 2, 3])
def test_toric_chi_table_matches_the_package_rr_layer(e):
    """chi from M2 lattice-point cohomology == hilbert_P over a full window,
    under a FITTED unimodular basis identification pinned by -K.  Any
    transcription error in the package's Gram, K, chi(O), or hilbert_P -- the
    layer every F_e verdict is built on -- fails here (zero fits) against an
    independent computation.

    The fit COUNT is lattice-theoretic: it equals the number of isometries of
    (NS, Gram) fixing K.  Two on F_0 (the ruling swap f <-> s) and two on F_2
    (sigma: f -> f + s, s -> -s -- both F and E + F are isotropic there, and
    sigma^T G sigma = G with sigma(K) = K, hand-verified; the chi-table cannot
    distinguish the two isotropic rays, only effectivity can).  Exactly one
    for the odd indices e = 1, 3."""
    table = fe_line_bundle_cohomology(e, -3, 3)
    assert table["coh"][(0, 0)] == (1, 0, 0)          # h^*(O_X) pin
    fits = _fit_identifications(e, table)
    assert len(fits) == (2 if e in (0, 2) else 1), fits


@requires_m2
@pytest.mark.parametrize("e", [0, 1, 2, 3])
def test_toric_table_is_serre_self_consistent(e):
    """Package-free internal check of the M2 data itself: Serre duality
    h^2(D) = h^0(K - D) wherever both classes lie in the window."""
    table = fe_line_bundle_cohomology(e, -3, 3)
    aK = table["antik"]
    checked = 0
    for (c, d), hs in table["coh"].items():
        kd = (-aK[0] - c, -aK[1] - d)                 # K - D in M2 coordinates
        if kd in table["coh"]:
            assert hs[2] == table["coh"][kd][0], ((c, d), kd)
            checked += 1
    assert checked > 20


@requires_m2
def test_flagship_character_chi_is_cas_witnessed():
    """The M3b flagship (2,(1,1),0) on F_0 is O(1,0) (+) O(0,1): its Euler
    characteristic decomposes as chi(O(1,0)) + chi(O(0,1)) = 2 + 2 = 4, and the
    package's RR pairing agrees -- with both line-bundle values read from the
    INDEPENDENT toric table under a fitted identification."""
    table = fe_line_bundle_cohomology(0, -3, 3)
    chi_m2 = _chi_table(table)
    (t0, t1) = _fit_identifications(0, table)[0]
    # invert T (unimodular) to find the M2 classes of O(1,0) and O(0,1)
    det = t0[0] * t1[1] - t0[1] * t1[0]
    inv = ((t1[1] * det, -t0[1] * det), (-t1[0] * det, t0[0] * det))
    def m2_class(fs):
        return (inv[0][0] * fs[0] + inv[0][1] * fs[1],
                inv[1][0] * fs[0] + inv[1][1] * fs[1])
    u1, u2 = m2_class((1, 0)), m2_class((0, 1))
    assert chi_m2[u1] == 2 and chi_m2[u2] == 2
    O = SurfaceBundle(1, (0, 0), F(0))
    flagship = SurfaceBundle(2, (1, 1), F(0))
    assert surface_chi(O, flagship, P1xP1) == chi_m2[u1] + chi_m2[u2] == 4
