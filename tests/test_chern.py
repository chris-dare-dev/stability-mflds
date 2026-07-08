"""Dedicated tests for :mod:`bridgeland_stability.chern` (E1-M3).

Convention pinned here (Coskun-Huizenga / Drezet-Le Potier), matching
``chern.py`` and CLAUDE.md invariant 2:

    mu(v)    = c / (r * d)                     (Mumford H-slope)
    Delta(v) = (1/2) * mu^2 - e / (r * d)      (normalized discriminant)
    discriminant_brief(v) = 2 * Delta(v)       (doubled -- comparison only)

ROADMAP-TEXT CORRECTION.  The E1-M3 roadmap prose states
``ChernChar(1,1,0).discriminant(2) == F(1,4)``.  That value is the *brief*
(doubled) discriminant, NOT the CH ``discriminant``.  The shipped code returns
``discriminant(2) == F(1,8)`` -- exactly re-derived below:
``mu = 1/(1*2) = 1/2``, ``Delta = (1/2)*(1/2)^2 - 0/(1*2) = 1/8`` -- and
``discriminant_brief(2) == 2*Delta == F(1,4)`` recovers the roadmap number.
Per invariants 5 (code is ground truth) and 6 (values must be exactly
re-derived), and the E1 mandate that NO computed value may change, these tests
pin the code-true CH re-derivation, not the mislabelled prose.

Every asserted class is either literature-anchored or exactly re-derived:
  * ``ChernChar(1,2,1)`` = O(1) on ``d=2``; line bundles have Delta = 0
    (chern.py module docstring).
  * ``ChernChar(2,0,-1/4)`` on ``d=1`` is the pinned P^2 wall class ``v`` from
    tests/test_walls.py (``Delta_v(CH) = 1/8``).
  * ``ChernChar(2,0,-3)`` on ``d=2`` is the G6 K3-accept class
    (docs/GOALS.md, ``Delta == 3/4``).

The module imports only stdlib ``fractions`` plus ``pytest`` and
``bridgeland_stability.chern``; no viz / matplotlib is pulled in, so the
zero-runtime-dependency import guarantee is preserved.
"""

import subprocess
import sys
from fractions import Fraction
from fractions import Fraction as F

from bridgeland_stability.chern import ChernChar
from bridgeland_stability.dlp import delta
from bridgeland_stability.exceptional import enumerate_exceptional
from bridgeland_stability.nslattice import rank1_shim, RANK1_AMPLE
from bridgeland_stability.varieties import P2, K3, abelian_surface


def test_discriminant_convention():
    """CH convention  Delta = (1/2) mu^2 - e/(rd)  and  brief == 2*Delta.

    Pinned on >= 3 distinct classes.  Note the roadmap-text value F(1,4) for
    ``ChernChar(1,1,0).discriminant(2)`` is the *brief*; the CH discriminant
    is F(1,8) (exactly re-derived in this module's docstring).
    """
    v = ChernChar(1, 1, 0)
    # Exact CH re-derivation: mu = 1/2, Delta = (1/2)(1/4) = 1/8.
    assert v.discriminant(2) == F(1, 8)
    # The doubled "brief" convention recovers the roadmap-prose number F(1/4).
    assert v.discriminant_brief(2) == F(1, 4)
    assert v.discriminant_brief(2) == 2 * v.discriminant(2)

    # (class, d, expected CH discriminant) -- literature-anchored / re-derived.
    cases = [
        (ChernChar(1, 2, 1), 2, F(0)),        # O(1) on d=2: line bundle, Delta=0
        (ChernChar(2, 0, F(-1, 4)), 1, F(1, 8)),  # pinned P^2 wall class v
        (ChernChar(2, 0, -3), 2, F(3, 4)),    # G6 K3-accept class, Delta=3/4
        (ChernChar(1, 1, 0), 2, F(1, 8)),     # exact re-derivation above
    ]
    for cls, d, delta in cases:
        assert cls.discriminant(d) == delta
        assert cls.discriminant_brief(d) == 2 * cls.discriminant(d)


def test_twist_shift_invariance():
    """``twist(s, d)`` (the B = sH twist, ch . e^{-sH}) is translation-invariant.

    The B-twist LOWERS the Mumford slope by exactly s (``mu -> mu - s``,
    since ``c -> c - s r d``) and leaves the discriminant unchanged.  Asserting
    a ``+ s`` shift would be wrong; the code-true relation is ``mu - s``.
    """
    for (r, c, e) in [(2, 3, 1), (1, 0, -2), (3, -1, 4), (2, 0, -3)]:
        for s in [1, -2, F(1, 2)]:
            for d in [1, 2]:
                v = ChernChar(r, c, e)
                tw = v.twist(s, d)
                # Additive slope shift by -s (the B=sH twist lowers mu by s).
                assert tw.slope(d) == v.slope(d) - s
                # Discriminant is translation-invariant under the B-twist.
                assert tw.discriminant(d) == v.discriminant(d)


def test_central_charge_matches_fraction_parts():
    """``central_charge`` float output == exact-Fraction derivation to 1e-12.

    Re Z = -(e - s c + (s^2 - t^2)/2 r d),  Im Z = t (c - s r d).  The exact
    parts stay ``Fraction`` right up to the ``float()`` comparison boundary --
    the only place float is allowed to appear.
    """
    r, c, e = F(1), F(0), F(-2)
    for (s, t, d) in [(F(-1), F(1), 2), (F(1, 2), F(3, 2), 1)]:
        # Exact Fraction derivation (no float until the boundary below).
        re_ex = -(e - s * c + (s * s - t * t) / 2 * r * d)
        im_ex = t * (c - s * r * d)
        re, im = ChernChar(1, 0, -2).central_charge(s, t, d)
        assert abs(re - float(re_ex)) < 1e-12
        assert abs(im - float(im_ex)) < 1e-12

    # Pin the literal primary case: v=(1,0,-2) at (s,t,d)=(-1,1,2) -> (2, 2).
    re, im = ChernChar(1, 0, -2).central_charge(-1, 1, 2)
    assert abs(re - 2.0) < 1e-12
    assert abs(im - 2.0) < 1e-12

    # Second INDEPENDENT literal anchor (not a mirror of the formula):
    # v=(1,0,-2) at (s,t,d)=(1/2, 3/2, 1).
    #   Re = -(-2 - 0 + ((1/4)-(9/4))/2 * 1 * 1) = -(-2 - 1) = 3
    #   Im = (3/2) * (0 - (1/2)*1*1) = (3/2)*(-1/2) = -3/4
    re, im = ChernChar(1, 0, -2).central_charge(F(1, 2), F(3, 2), 1)
    assert abs(re - 3.0) < 1e-12
    assert abs(im - (-0.75)) < 1e-12


def test_discriminant_lattice_terms():
    """slope/discriminant/bogomolov rewritten in NS-lattice terms; values UNCHANGED.

    (A) DLP fractal-envelope regression -- delta routes through Bundle.discriminant
        (roadmap-named coarse gate); values pinned in CORRECTIONS.md sec2 / test_dlp.py.
    (B) ChernChar's discriminant / bogomolov now come from the NS-lattice pairing:
        ch1(d)=(c/d,); <ch1,H>=c; <ch1,ch1>=c^2/d on the rank-1 shim (<H,H>=d).
    """
    # (A) unchanged after the rewrite.
    bs = enumerate_exceptional(-3, 3, R_max=60)
    assert delta(F(1, 2), bs) == F(5, 8)     # CORRECTIONS.md sec2 / test_dlp.py
    assert delta(F(2, 5), bs) == F(13, 25)
    # (B) discriminant/bogomolov now from NS-lattice terms, values identical.
    # (r,  c,  e,     d, ch1_H,  ch1_sq, disc,   bog)
    cases = [
        (1,  1, F(0),  2, F(1),   F(1, 2), F(1, 8), F(1, 2)),  # mu=1/2, disc=1/2*1/4=1/8
        (2,  0, F(-3), 2, F(0),   F(0),    F(3, 4), F(12)),    # disc=0+3/4; bog=0+12
        (1,  0, F(-2), 1, F(0),   F(0),    F(2),    F(4)),     # P^2[2] v; disc=2, bog=4
        (1, -1, F(1, 2), 1, F(-1), F(1),    F(0),    F(0)),    # O(-1); disc=1/2-1/2=0
    ]
    for (r, c, e, d, ch1H, ch1sq, disc, bog) in cases:
        v = ChernChar(r, c, e)
        assert v.ch1(d) == (Fraction(c, d),)
        assert v.ch1_dot_H(d) == ch1H
        assert v.ch1_squared(d) == ch1sq
        lat = rank1_shim(d)                                  # independent recompute
        assert v.ch1_dot_H(d) == lat.pairing(v.ch1(d), RANK1_AMPLE)
        assert v.ch1_squared(d) == lat.self_pairing(v.ch1(d))
        assert v.slope(d) == Fraction(c, r * d)              # unchanged vs closed form
        assert v.discriminant(d) == disc
        assert v.bogomolov_discriminant(d) == bog
        assert v.discriminant_brief(d) == 2 * v.discriminant(d)


def test_bogomolov_binds_to_surface_lattice():
    """<ch1,ch1> - 2 r e from the ACTUAL Surface.lattice == shipped bogomolov."""
    for surf in (P2, K3(2), abelian_surface(2)):
        for (r, c, e) in [(1, 1, F(0)), (2, 0, F(-3)), (1, 0, F(-2))]:
            v = ChernChar(r, c, e)
            assert (surf.lattice.self_pairing(v.ch1(surf.d)) - 2 * r * e
                    == v.bogomolov_discriminant(surf.d))
            assert v.ch1_dot_H(surf.d) == surf.lattice.pairing(v.ch1(surf.d), surf.H)


def test_naive_encoding_not_used_by_chernchar():
    """Guards the E8 danger-zone trap: ChernChar stores c/d, never the naive scalar c.

    Naive (c,)+Gram=[[d]] would give <ch1,H>=c*d=2, <ch1,ch1>=c^2*d=2 for c=1,d=2,
    breaking delta(1/2)=5/8; ChernChar stores c/d -> 1 and 1/2.
    """
    v = ChernChar(1, 1, 0)
    assert v.ch1_dot_H(2) == 1 and v.ch1_dot_H(2) != 1 * 2
    assert v.ch1_squared(2) == F(1, 2) and v.ch1_squared(2) != 1 * 2


def test_chern_is_stdlib_only():
    """chern.py's new nslattice import stays stdlib-only (invariant 3)."""
    child = (
        "import sys\n"
        "class _B:\n"
        "    def find_spec(self, n, p=None, t=None):\n"
        "        if n.split('.')[0] in {'matplotlib','plotly'}:\n"
        "            raise ImportError('blocked: '+n)\n"
        "        return None\n"
        "sys.meta_path.insert(0, _B())\n"
        "import bridgeland_stability.chern\n"
        "assert 'matplotlib' not in sys.modules and 'plotly' not in sys.modules\n"
        "print('STDLIB_ONLY_OK')\n"
    )
    r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "STDLIB_ONLY_OK" in r.stdout
