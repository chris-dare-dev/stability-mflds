"""E4-M2 / G9a: the PROVEN tilt-layer (nu) wall solver ``numerical_tilt_wall``.

The tilt slope ``nu_{alpha,beta}`` (E4-M1, ``threefold.nu``) depends ONLY on the
truncated triple ``(r, a1, a2) = (ch0, ch1.H^2, ch2.H)`` and NEVER on ``a3 = ch3``.
The tilt wall ``nu_{alpha,beta}(v) = nu_{alpha,beta}(w)`` is therefore EXACTLY the
surface ``(s, t)`` semicircle of ``walls.numerical_wall`` under the substitution
``(s, t) -> (beta, alpha)``, ``(r, c, e) -> (r, a1, a2)``, ``d -> d3``.  This adds
NO new mathematics -- it re-reads the existing exact primitive on a truncated
character.  Sources: Bayer-Macri-Toda arXiv:1103.5010 (tilt slope); Schmidt
arXiv:1509.04608 (threefold tilt-stability walls).

Two-way verification of the load-bearing pin (twisted cubic on P^3, d3 = 1).
``V = ThreefoldChern(1, 0, -3, 5)`` (truncated ``(1, 0, -3)``);
``W = ThreefoldChern(1, -1, 1/2, -1/6)`` = ``ch(O(-1))`` on P^3 (untwisted Chern
character; truncated to ``(1, -1, 1/2)``).

(a) via numerical_wall minors on the truncated triples::

    W_rc = r_v a1_w - r_w a1_v = 1*(-1) - 1*0     = -1
    W_re = r_v a2_w - r_w a2_v = 1*(1/2) - 1*(-3)  = 7/2
    W_ce = a1_v a2_w - a1_w a2_v = 0*(1/2) - (-1)*(-3) = -3
    center   = W_re/W_rc = (7/2)/(-1) = -7/2
    radius^2 = center^2 - 2 W_ce/(d3 W_rc) = 49/4 - 2(-3)/(1(-1)) = 49/4 - 6 = 25/4

(b) via clearing nu(v) = nu(w) directly (independent route, no minor reuse):
with ``nu(E) = (a2^b - 1/2 alpha^2 r d3)/a1^b``, ``a1^b = a1 - beta r d3``,
``a2^b = a2 - beta a1 + 1/2 beta^2 r d3``, cross-multiplying ``nu_v = nu_w`` and
collecting gives ``alpha^2 + beta^2 + 7 beta + 6 = 0``, i.e.
``(beta + 7/2)^2 + alpha^2 = 25/4`` -> center ``-7/2``, radius^2 ``25/4``.
Both routes agree exactly, and the alpha^2 / beta^2 coefficients emerge equal
(the beta-axis-centered-circle fact).
"""

from fractions import Fraction as F

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.threefold import (
    ThreefoldChern,
    ThreefoldTiltWall,
    ThreefoldVerticalTiltWall,
    numerical_tilt_wall,
    tilt_wall_coeffs,
    BridgelandWall,
    bridgeland_wall,
    bridgeland_wall_coeffs,
    enumerate_tilt_walls,
)
from bridgeland_stability.varieties import P3, QUADRIC3, BLOWUP_P3_POINT
from bridgeland_stability.walls import numerical_wall


V = ThreefoldChern(1, F(0), F(-3), F(5))
W = ThreefoldChern(1, F(-1), F(1, 2), F(-1, 6))


def test_tilt_wall_value():
    """Load-bearing pin: twisted cubic vs ch(O(-1)), d3 = 1."""
    w = numerical_tilt_wall(V, W, 1)
    assert isinstance(w, ThreefoldTiltWall)
    assert w.center == F(-7, 2)
    assert w.radius_sq == F(25, 4)
    assert w.is_real
    assert w.radius == 2.5  # float display only (sqrt(25/4))
    assert w.bridgeland_certified is False


def test_tilt_wall_ch3_independent():
    """center and radius_sq are unchanged when only a3 = ch3 varies."""
    for a3v, a3w in [(F(5), F(-1, 6)), (F(-11), F(7)), (F(1, 3), F(0))]:
        v2 = ThreefoldChern(1, F(0), F(-3), a3v)
        w2 = ThreefoldChern(1, F(-1), F(1, 2), a3w)
        wall = numerical_tilt_wall(v2, w2, 1)
        assert wall.center == F(-7, 2)
        assert wall.radius_sq == F(25, 4)


def test_tilt_wall_reduction_identity():
    """Bit-for-bit equal to numerical_wall on the truncated triple (same d3)."""
    tilt = numerical_tilt_wall(V, W, 1)
    surf = numerical_wall(
        ChernChar(V.r, V.a1, V.a2), ChernChar(W.r, W.a1, W.a2), 1
    )
    assert tilt.center == surf.center
    assert tilt.radius_sq == surf.radius_sq


def test_tilt_wall_d3_normalization():
    """beta-center is d3-invariant; radius^2 tracks d3 via -2 W_ce/(d3 W_rc)."""
    tilt = numerical_tilt_wall(V, W, 2)
    surf = numerical_wall(
        ChernChar(V.r, V.a1, V.a2), ChernChar(W.r, W.a1, W.a2), 2
    )
    assert tilt.center == surf.center == F(-7, 2)
    # 49/4 - 2(-3)/(2(-1)) = 49/4 - 3 = 37/4
    assert tilt.radius_sq == surf.radius_sq == F(37, 4)


def test_tilt_wall_beta_axis_centered():
    """Equal alpha^2 / beta^2 coefficients => beta-axis-centered semicircle."""
    a2, b2, b1, c0 = tilt_wall_coeffs(V, W, 1)
    assert a2 == b2 == F(1, 2)  # equal alpha^2 and beta^2 coefficients
    assert b1 == F(7, 2)
    assert c0 == F(3)
    wall = numerical_tilt_wall(V, W, 1)
    assert -b1 / (2 * b2) == wall.center
    assert (b1 * b1 - 4 * b2 * c0) / (4 * b2 * b2) == wall.radius_sq


def test_tilt_wall_vertical_when_equal_slopes():
    """W_rc = 0 (equal tilt-Mumford slope) yields a vertical beta = const wall."""
    v = ThreefoldChern(2, F(2), F(0), F(0))
    w = ThreefoldChern(1, F(1), F(0), F(0))
    assert tilt_wall_coeffs(v, w, 1)[:2] == (0, 0)
    wall = numerical_tilt_wall(v, w, 1)
    assert isinstance(wall, ThreefoldVerticalTiltWall)
    assert wall.beta_value == F(1)  # a1/(r d3) = 2/(2*1)
    # d3 != 1 distinguishes the correct a1/(r*d3) from a1/r (W_rc = 0 is
    # d3-independent, so the wall stays vertical): beta = 2/(2*2) = 1/2.
    wall2 = numerical_tilt_wall(v, w, 2)
    assert isinstance(wall2, ThreefoldVerticalTiltWall)
    assert wall2.beta_value == F(1, 2)


def test_tilt_wall_certificate_proven():
    """The tilt-layer reduction is PROVEN; but it is NOT a Bridgeland wall."""
    wall = numerical_tilt_wall(V, W, 1)
    assert wall.certificate.rigor == Rigor.PROVEN
    assert "arXiv:1509.04608" in wall.certificate.citations
    assert wall.bridgeland_certified is False


# ==========================================================================
# E4-M3 / G9b -- second-tilt / Bridgeland wall `bridgeland_wall` [RESEARCH-light]
#
# The exact ch3-dependent (y, z) coefficients (and hence the twisted-cubic
# center/radius numeric value) are DEFERRED pending Schmidt arXiv:1509.04608
# sec.3/sec.7; NONE of these tests pins a deferred Schmidt number.  They pin only
# STRUCTURE that is exactly re-derived: the (r,a1)-minor x coefficient (Schmidt
# Thm 3.3 / E4-M2 anchored), the quadratic-formula circle identity, the exact
# reduction to the E4-M2 tilt wall when the ch3-minors vanish, the ch3-dependence
# of (y, z), the bg_proven gate, and the deferred/conjectural certificate.
# ==========================================================================


def test_bridgeland_gate():
    """`bridgeland_certified` echoes `threefold.bg_proven`: a genuine Bridgeland
    wall only where the threefold BG inequality is a theorem."""
    assert bridgeland_wall(V, W, P3).bridgeland_certified is True
    assert bridgeland_wall(V, W, QUADRIC3).bridgeland_certified is True
    assert bridgeland_wall(V, W, BLOWUP_P3_POINT).bridgeland_certified is False


def test_bridgeland_locus_not_verified():
    # [RESEARCH-light]: do NOT conflate bridgeland_certified (the THREEFOLD's BG
    # theorem -> Bridgeland walls exist) with locus_verified (THIS computed wall's
    # shape + (y,z) coeffs).  The beta-axis circle SHAPE is itself a candidate -- the
    # cubic ch3^beta term means the true locus is not automatically a circle -- so
    # locus_verified is always False until confirmed against Schmidt sec.3/sec.7.
    bw = bridgeland_wall(V, W, P3)
    assert bw.bridgeland_certified is True     # P^3 BG is a theorem
    assert bw.locus_verified is False          # but the (y,z) circle locus is deferred
    assert bw.certificate.rigor == Rigor.CONJECTURAL
    # the certificate flags the SHAPE (not just the (y,z) coeffs) as a candidate
    assert any("SHAPE" in h for h in bw.certificate.hypotheses)


def test_bridgeland_is_circle():
    """center = -y/(2x), radius_sq = (y^2-4xz)/(4x^2) exactly (quadratic-formula
    circle identity); x = 1/2 is the (r,a1)-minor value (Schmidt Thm 3.3 anchored)."""
    x, y, z = bridgeland_wall_coeffs(V, W, 1)
    bw = bridgeland_wall(V, W, P3)
    assert isinstance(bw, BridgelandWall)
    assert bw.center == -y / (2 * x)
    assert bw.radius_sq == (y * y - 4 * x * z) / (4 * x * x)
    assert x == F(1, 2)  # (r,a1) minor: matches the E4-M2 tilt wall's alpha^2 coeff


def test_bridgeland_x_from_r_a1_minor():
    """x depends ONLY on the (r,a1) minor: invariant under (a2,a3) variations."""
    d3 = 1
    for a2v, a3v, a2w, a3w in [
        (F(-3), F(5), F(1, 2), F(-1, 6)),
        (F(0), F(0), F(7), F(11)),
        (F(-100), F(2, 3), F(5, 2), F(-9)),
    ]:
        vv = ThreefoldChern(1, F(0), a2v, a3v)
        ww = ThreefoldChern(1, F(-1), a2w, a3w)
        x = bridgeland_wall_coeffs(vv, ww, d3)[0]
        assert x == -F(d3, 2) * (vv.r * ww.a1 - ww.r * vv.a1)
        assert x == F(1, 2)


def test_bridgeland_reduces_to_tilt_when_ch3_minors_vanish():
    """With a3 = 0 (both), (x,y,z) collapse to the E4-M2 tilt wall exactly."""
    V0 = ThreefoldChern(1, F(0), F(-3), F(0))
    W0 = ThreefoldChern(1, F(-1), F(1, 2), F(0))
    bw = bridgeland_wall(V0, W0, P3)
    tw = numerical_tilt_wall(V0, W0, 1)
    assert bw.center == tw.center == F(-7, 2)  # E4-M2 / Schmidt Thm 3.3 pin
    assert bw.radius_sq == tw.radius_sq == F(25, 4)


def test_bridgeland_ch3_enters_yz():
    """a3 = ch3 enters the (y,z) coefficients but never x (why it stays a circle)."""
    c1 = bridgeland_wall_coeffs(ThreefoldChern(1, F(0), F(-3), F(5)), W, 1)
    c2 = bridgeland_wall_coeffs(ThreefoldChern(1, F(0), F(-3), F(99)), W, 1)
    assert c1[0] == c2[0]                       # x unchanged by ch3
    assert (c1[1], c1[2]) != (c2[1], c2[2])     # (y, z) move with ch3


def test_bridgeland_vertical_when_equal_slopes():
    """W_rc = 0 (equal tilt-Mumford slope) => vertical beta = a1/(r d3) wall,
    still gated by bg_proven."""
    v = ThreefoldChern(2, F(2), F(0), F(0))
    w = ThreefoldChern(1, F(1), F(0), F(0))
    bw = bridgeland_wall(v, w, P3)
    assert isinstance(bw, ThreefoldVerticalTiltWall)
    assert bw.beta_value == F(1)  # a1/(r d3) = 2/(2*1)
    assert bw.bridgeland_certified is True


def test_bridgeland_certificate_deferred():
    """The ch3-dependent (y,z) is a CANDIDATE: certificate is CONJECTURAL, cites
    Schmidt arXiv:1509.04608, and flags the deferred transcription."""
    bw = bridgeland_wall(V, W, P3)
    assert bw.certificate.rigor <= Rigor.CONJECTURAL
    assert "arXiv:1509.04608" in bw.certificate.citations
    note = bw.certificate.note.lower()
    assert "candidate" in note or "transcription" in note
    assert bridgeland_wall(V, W, BLOWUP_P3_POINT).certificate.rigor == Rigor.CONJECTURAL


def test_enumerate_tilt_walls():
    """`enumerate_tilt_walls` returns the PROVEN, ch3-independent tilt (nu) layer,
    sorted by radius_sq desc, analogous to `walls.actual_walls`."""
    walls = enumerate_tilt_walls(V, P3)
    assert len(walls) >= 1
    assert all(isinstance(w, ThreefoldTiltWall) for w in walls)
    radii = [w.radius_sq for w in walls]
    assert radii == sorted(radii, reverse=True)
    # ch3-independence: varying only a3 leaves the (center, radius_sq) set fixed.
    walls99 = enumerate_tilt_walls(ThreefoldChern(1, F(0), F(-3), F(99)), P3)
    assert [(w.center, w.radius_sq) for w in walls] == \
           [(w.center, w.radius_sq) for w in walls99]
