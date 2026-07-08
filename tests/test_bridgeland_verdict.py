"""E5-M1 / G15 -- conjecture-gated Bridgeland-wall verdict + rigor propagation.

PURE ADDITIVE LABELING (G15).  This milestone adds a machine-readable verdict on
whether a PROVEN tilt (nu) wall lifts to a genuine *Bridgeland* wall on a given
threefold; it changes NO ``Fraction`` value.  The verdict is PROVEN/certified iff
``threefold.bg_proven`` and CONJECTURAL/tilt-only otherwise -- it never emits
"no stability condition".

Provenance of the pinned numbers:
  * The tilt-wall geometry pin ``center=-7/2, radius_sq=25/4`` (d3=1) is the
    E4-M2 load-bearing value, re-derived below from the three ``numerical_wall``
    minors on the truncated triple ``(r, a1, a2)`` (ch3-independent):
        W_rc = 1*(-1) - 1*0   = -1
        W_re = 1*(1/2) - 1*(-3) = 7/2
        W_ce = 0*(1/2) - (-1)*(-3) = -3
        center   = W_re/W_rc = (7/2)/(-1) = -7/2         (d3-independent)
        radius^2 = center^2 - 2 W_ce/(d3 W_rc)
                 = 49/4 - 2(-3)/(d3*(-1)) = 49/4 - 6/d3.
    d3=1 -> 25/4;  d3=2 -> 37/4;  d3=5 -> 221/20.
  * NO deferred Schmidt second-tilt (y,z) number is pinned here -- that is E5-M2.

The meet law: rigor = min(TILT_SOLVER_RIGOR, PROVEN if bg_proven else CONJECTURAL),
with TILT_SOLVER_RIGOR == Rigor.PROVEN (the tilt (nu) solver is a proven identity).
"""

from fractions import Fraction as F

from bridgeland_stability.rigor import Rigor, Certificate
from bridgeland_stability.threefold import (
    ThreefoldChern,
    ThreefoldTiltWall,
    numerical_tilt_wall,
    bg_boundary_curve,
    BridgelandWallVerdict,
    bridgeland_wall_verdict,
    is_bridgeland_certified,
    TILT_SOLVER_RIGOR,
)
from bridgeland_stability.varieties import (
    ALL_THREEFOLDS,
    P3,
    QUADRIC3,
    QUINTIC,
    BLOWUP_P3_POINT,
    Threefold,
    fano_picard_one,
    abelian_threefold,
)

# Twisted-cubic-style class and a destabilizing subobject (E4-M2 fixture).
V = ThreefoldChern(1, F(0), F(-3), F(5))
W = ThreefoldChern(1, F(-1), F(1, 2), F(-1, 6))

# Every bg_proven=True catalog row, including the two factory families.
PROVEN_ROWS = [
    P3,
    QUADRIC3,
    QUINTIC,
    fano_picard_one("Fano rho=1 (generic)", 4),
    abelian_threefold(6),
]


def test_propagation_proven():
    """For every bg_proven=True row the verdict is PROVEN and certified."""
    for tf in PROVEN_ROWS:
        assert tf.bg_proven is True
        verdict = bridgeland_wall_verdict(V, W, tf)
        assert verdict.rigor == Rigor.PROVEN
        assert verdict.certified is True
        assert is_bridgeland_certified(tf) is True


def test_propagation_blowup():
    """Bl_p(P^3): tilt-only, rigor <= CONJECTURAL, NOT certified."""
    verdict = bridgeland_wall_verdict(V, W, BLOWUP_P3_POINT)
    assert BLOWUP_P3_POINT.bg_proven is False
    assert verdict.rigor <= Rigor.CONJECTURAL
    assert verdict.rigor == Rigor.CONJECTURAL
    assert verdict.certified is False
    assert is_bridgeland_certified(BLOWUP_P3_POINT) is False


def test_meet_law():
    """rigor == min(TILT_SOLVER_RIGOR, PROVEN if bg_proven else CONJECTURAL) exactly.

    Includes the synthetic-row guard (design decision 2): a Threefold with the
    default UNKNOWN_CERTIFICATE must still yield CONJECTURAL (not UNKNOWN) when
    bg_proven is False -- the rigor is derived from the bg_proven BOOL, NOT folded
    from threefold.certificate.rigor.
    """
    assert TILT_SOLVER_RIGOR == Rigor.PROVEN
    # Catalog rows.
    for tf in list(ALL_THREEFOLDS) + PROVEN_ROWS:
        bg_rigor = Rigor.PROVEN if tf.bg_proven else Rigor.CONJECTURAL
        expected = min(TILT_SOLVER_RIGOR, bg_rigor)
        assert bridgeland_wall_verdict(V, W, tf).rigor == expected
    # Synthetic rows carrying the DEFAULT UNKNOWN_CERTIFICATE (rigor 0).
    synth_unproven = Threefold(name="synthetic X", d3=1, bg_proven=False)
    synth_proven = Threefold(name="synthetic Y", d3=3, bg_proven=True)
    assert synth_unproven.certificate.rigor == Rigor.UNKNOWN  # would poison a fold
    assert bridgeland_wall_verdict(V, W, synth_unproven).rigor == Rigor.CONJECTURAL
    assert bridgeland_wall_verdict(V, W, synth_proven).rigor == Rigor.PROVEN


def test_never_certifies_unproven():
    """NO code path returns certified=True when bg_proven is False."""
    unproven = [BLOWUP_P3_POINT, Threefold(name="synthetic Z", d3=2, bg_proven=False)]
    for tf in unproven:
        verdict = bridgeland_wall_verdict(V, W, tf)
        assert verdict.certified is False
        assert verdict.rigor < Rigor.PROVEN


def test_blowup_note():
    """Bl_p(P^3) note states stability conditions EXIST; never 'no stability condition'."""
    note = bridgeland_wall_verdict(V, W, BLOWUP_P3_POINT).note
    assert "tilt-stability wall only; Bridgeland upgrade unproven" in note
    assert "stability conditions nonetheless EXIST" in note
    assert "strong BMT boundary FAILS" in note
    assert "BMSZ arXiv:1607.08199" in note
    assert "no stability condition" not in note.lower()


def test_unproven_note_exact_phrase():
    """A non-blowup unproven threefold carries the EXACT roadmap unproven phrase,
    with NO blowup-only extra clause."""
    tf = Threefold(name="generic conjectural 3-fold", d3=1, bg_proven=False, kind="general")
    note = bridgeland_wall_verdict(V, W, tf).note
    assert note == (
        "tilt-stability wall only; Bridgeland upgrade unproven "
        "(threefold BG open here)"
    )
    assert "no stability condition" not in note.lower()


def test_certified_note_and_never_no_stability():
    """Certified rows carry the certified note; NO verdict emits 'no stability condition'."""
    for tf in list(ALL_THREEFOLDS) + PROVEN_ROWS:
        note = bridgeland_wall_verdict(V, W, tf).note
        assert "no stability condition" not in note.lower()
    p3_note = bridgeland_wall_verdict(V, W, P3).note
    assert "Bridgeland wall certified" in p3_note
    assert "Piyaratne-Toda arXiv:1504.01177" in p3_note


def test_tilt_wall_geometry_unchanged():
    """The wrapped tilt_wall is the RAW numerical_tilt_wall locus; no Fraction moves.

    Re-derived pins: center=-7/2 for all d3; radius_sq = 49/4 - 6/d3.
    """
    cases = [(P3, F(25, 4)), (QUADRIC3, F(37, 4)), (QUINTIC, F(221, 20))]
    for tf, expected_r2 in cases:
        verdict = bridgeland_wall_verdict(V, W, tf)
        raw = numerical_tilt_wall(V, W, tf.d3)
        assert isinstance(verdict.tilt_wall, ThreefoldTiltWall)
        assert verdict.tilt_wall.center == F(-7, 2)
        assert verdict.tilt_wall.radius_sq == expected_r2
        # identical to the raw solver output (byte-for-byte Fraction)
        assert verdict.tilt_wall.center == raw.center
        assert verdict.tilt_wall.radius_sq == raw.radius_sq
        # the wrapped tilt wall stays a tilt wall; the bg gate lives in .certified
        assert verdict.tilt_wall.bridgeland_certified is False


def test_citations_and_certificate_property():
    """citations dedup the tilt-cert refs + threefold refs + Piyaratne-Toda; the
    .certificate property surfaces (rigor, citations, note)."""
    verdict = bridgeland_wall_verdict(V, W, P3)
    cits = verdict.citations
    assert "arXiv:1504.01177" in cits                 # conditional properness (always)
    assert "arXiv:1103.5010" in cits                  # BMT tilt slope (tilt cert)
    assert "arXiv:1509.04608" in cits                 # Schmidt tilt walls (tilt cert)
    for ref in P3.references:                         # threefold BG refs unioned in
        assert ref in cits
    assert len(cits) == len(set(cits))                # order-preserving dedup, no dups
    cert = verdict.certificate
    assert isinstance(cert, Certificate)
    assert cert.rigor == verdict.rigor == Rigor.PROVEN
    assert cert.citations == verdict.citations
    assert cert.note == verdict.note


def test_is_bridgeland_certified_routed_through_bg_boundary():
    """is_bridgeland_certified agrees with the Certificate bg_boundary_curve attaches
    to every emitted threefold wall (PROVEN iff certified)."""
    for tf in list(ALL_THREEFOLDS) + PROVEN_ROWS:
        bg = bg_boundary_curve(V, tf)
        certified = is_bridgeland_certified(tf)
        assert certified == (bg.certificate.rigor == Rigor.PROVEN)
        assert certified == bridgeland_wall_verdict(V, W, tf).certified


def test_verdict_is_frozen_dataclass():
    """BridgelandWallVerdict is a frozen dataclass (create, don't mutate)."""
    import dataclasses

    verdict = bridgeland_wall_verdict(V, W, P3)
    assert dataclasses.is_dataclass(verdict)
    assert isinstance(verdict, BridgelandWallVerdict)
    try:
        verdict.certified = False  # type: ignore[misc]
    except dataclasses.FrozenInstanceError:
        pass
    else:  # pragma: no cover - would signal a non-frozen regression
        raise AssertionError("BridgelandWallVerdict must be frozen")


# -----------------------------------------------------------------------------
# E5-M2 / G15 -- Schmidt first-wall ground-truth cross-check (TEST-ONLY).
#
# These tests cross-check the E5-M1 verdict pipeline against the E4 PROVEN tilt
# (nu) wall solver on Schmidt's twisted-cubic object (arXiv:1509.04608), and
# contrast the SAME class on P^3 (PROVEN) vs Bl_p(P^3) (CONJECTURAL) on identical
# ``Fraction`` arithmetic.  NO new mathematics; NO computed value changes -- they
# only pin values the shipped E4/E5 code already produces.
# -----------------------------------------------------------------------------


def test_twisted_cubic_is_schmidt_equal_chern_object():
    """HRR re-derivation of ``ch(I_twisted cubic) = (1,0,-3,5)`` and Schmidt's
    equal-``ch`` plane-cubic-union-point identity (arXiv:1509.04608).

    Two-way verification (HRR on P^3, d3 = H^3 = 1): for a curve C of degree e
    (``a2 = ch2.H = e``) and arithmetic genus g, ``chi(O_C) = 1 - g`` and on P^3
    ``chi(O_C) = ch3 + 2 a2``, so ``ch3(O_C) = (1 - g) - 2e``.  Pure exact-
    ``Fraction`` identity checks; no floats, no computed value moves.
    """
    O = ThreefoldChern(1, F(0), F(0), F(0))          # ch(O_{P^3})
    # twisted cubic O_C: e=3, g=0 -> ch3 = (1-g) - 2e = 1 - 6 = -5
    O_C = ThreefoldChern(0, F(0), F(3), F(-5))
    I_twisted = ThreefoldChern(1, O.a1 - O_C.a1, O.a2 - O_C.a2, O.a3 - O_C.a3)
    assert I_twisted == V == ThreefoldChern(1, F(0), F(-3), F(5))
    # plane cubic (e=3, g=1 -> ch3 = 0 - 6 = -6) UNION a point (ch = (0,0,0,1))
    O_planecubic = ThreefoldChern(0, F(0), F(3), F(-6))
    O_union = ThreefoldChern(0, F(0), O_planecubic.a2, O_planecubic.a3 + F(1))
    I_union = ThreefoldChern(1, O.a1 - O_union.a1, O.a2 - O_union.a2, O.a3 - O_union.a3)
    assert I_union == V   # identical Chern character -> Schmidt's example object


def test_schmidt_first_wall_p3():
    """Schmidt (arXiv:1509.04608) twisted-cubic first-wall object V on P^3.

    The E5-M1 verdict reproduces the E4 PROVEN tilt (nu) wall locus exactly and
    labels it PROVEN.  The pinned (beta,alpha) locus is the exactly-re-derived
    E4-M2 value (HRR ch + ``numerical_wall`` minors: center=-7/2, radius_sq=25/4),
    cross-checked byte-for-byte equal to the raw E4 solver.

    [RESEARCH] The exact Schmidt arXiv:1509.04608 section 7 numeric transcription
    (and confirmation that O(-1) is the true destabilizer) is DEFERRED -- no
    computed value may change in Epic 1.  The pinned first-wall (beta,alpha) value
    is the exactly-re-derived E4-M2 tilt wall, NOT a raw section 7 read-off.
    """
    raw = numerical_tilt_wall(V, W, P3.d3)            # P3.d3 == 1
    verdict = bridgeland_wall_verdict(V, W, P3)
    assert isinstance(verdict.tilt_wall, ThreefoldTiltWall)
    # (beta,alpha) locus matches the E4 solver EXACTLY (as Fraction)
    assert verdict.tilt_wall.center == raw.center == F(-7, 2)
    assert verdict.tilt_wall.radius_sq == raw.radius_sq == F(25, 4)
    assert verdict.tilt_wall.radius == 2.5            # float display only (sqrt(25/4))
    # labeled PROVEN on P^3
    assert verdict.rigor == Rigor.PROVEN
    assert verdict.certified is True
    # it is the ch3-independent tilt (nu) wall, not a certified Bridgeland wall
    assert verdict.tilt_wall.bridgeland_certified is False


def test_same_class_blowup_is_conjectural():
    """The IDENTICAL class V (destabilizer W) on Bl_p(P^3): identical ``Fraction``
    tilt arithmetic (both d3=1), downgraded rigor; never emits 'no stability
    condition'.  Same theorem-vs-conjecture contrast on byte-identical geometry.
    """
    vp = bridgeland_wall_verdict(V, W, P3)
    vb = bridgeland_wall_verdict(V, W, BLOWUP_P3_POINT)
    # identical (beta,alpha) locus -- same arithmetic, different rigor label
    assert vb.tilt_wall.center == vp.tilt_wall.center == F(-7, 2)
    assert vb.tilt_wall.radius_sq == vp.tilt_wall.radius_sq == F(25, 4)
    # same class: PROVEN on P^3, CONJECTURAL / tilt-only on Bl_p(P^3)
    assert vp.rigor == Rigor.PROVEN and vp.certified is True
    assert vb.rigor == Rigor.CONJECTURAL
    assert vb.rigor <= Rigor.CONJECTURAL
    assert vb.certified is False
    # never "no stability condition"; Bl_p(P^3)-specific existence statement
    assert "no stability condition" not in vb.note.lower()
    assert "stability conditions nonetheless EXIST" in vb.note
