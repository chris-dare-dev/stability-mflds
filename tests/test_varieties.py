"""E1-M1 / G1: threefold citation-provenance guard.

Pins the corrected threefold references and enforces the curated
algebraic-geometry allowlist that catches a stray non-AG citation -- the one
confirmed shipped defect: QUADRIC3 formerly cited arXiv:1607.07182
("Interlacing Diffusions", a math.PR probability paper) and mis-attributed
Piyaratne's arXiv:1705.04011 as BMSZ.

Maintenance note: AG_ALLOWLIST is deliberately curated and must be appended
whenever a legitimately new AG reference is added.  The membership check, NOT a
format regex, is the primary guard: the probability paper's ID matches
^arXiv:\\d{4}\\.\\d{4,5}$ perfectly, so a regex-only guard would not have caught
the original defect.
"""

import re
from dataclasses import replace
from fractions import Fraction as F

import pytest

from bridgeland_stability.rigor import Rigor
from bridgeland_stability.varieties import (
    ALL_THREEFOLDS,
    AG_ALLOWLIST,
    QUADRIC3,
    BLOWUP_P3_POINT,
    fano_picard_one,
    abelian_threefold,
    koseki_product,
    require_faithful_computation,
    KOSEKI_BG_VERIFIED,
    P2,
    P1xP1,
    K3,
    abelian_surface,
    hirzebruch,
    enriques,
    bielliptic,
)
from bridgeland_stability.chern import ChernChar
from bridgeland_stability.walls import (
    compute_walls,
    actual_walls,
    numerical_wall,
    abelian_wall,
)
from bridgeland_stability.mukai import k3_wall_classified
from bridgeland_stability.dlp import moduli_nonempty
from bridgeland_stability.bg_check import check_existence_k3, check_existence_abelian
from bridgeland_stability.exceptional import Bundle

PROBABILITY_PAPER = "arXiv:1607.07182"  # "Interlacing Diffusions" (math.PR)
ARXIV_FORMAT = re.compile(r"^arXiv:\d{4}\.\d{4,5}$")


def _all_refs(threefolds):
    return sum((tf.references for tf in threefolds), [])


def test_no_probability_paper_in_catalog():
    assert PROBABILITY_PAPER not in _all_refs(ALL_THREEFOLDS)


def test_quadric_refs():
    assert QUADRIC3.references == [
        "arXiv:1510.04089",
        "arXiv:1607.08199",
        "arXiv:1509.04608",
    ]


def test_quadric_note_drops_piyaratne_bmsz_conflation():
    # The old note implied 1705.04011 == BMSZ; the ID is gone from QUADRIC3 and
    # BMSZ is now correctly attributed to 1607.08199.
    assert "1705.04011" not in QUADRIC3.note
    assert "1607.08199" in QUADRIC3.note


def test_blowup_keeps_piyaratne():
    assert "arXiv:1705.04011" in BLOWUP_P3_POINT.references
    assert "arXiv:1602.05055" in BLOWUP_P3_POINT.references


def test_blowup_note_attributes_piyaratne_and_fano():
    note = BLOWUP_P3_POINT.note
    assert "Piyaratne" in note                       # 1705.04011 -> Piyaratne
    assert "strong BMT boundary" in note             # only the strong bound fails
    assert "no stability condition" not in note.lower()  # Bl_p(P^3) is Fano
    # Pin the load-bearing POSITIVE half of the fix too: Bl_p(P^3) is Fano and
    # DOES carry Bridgeland stability conditions (BMSZ, 1607.08199).  A regression
    # dropping this clause must fail, not just one that adds "no stability condition".
    assert "Fano" in note
    assert "1607.08199" in note
    assert "stability conditions still exist" in note


def test_allowlist_membership():
    catalog = list(ALL_THREEFOLDS) + [
        fano_picard_one("Fano rho=1 (generic)", 4),
        abelian_threefold(6),
    ]
    for tf in catalog:
        for ref in tf.references:
            assert ref in AG_ALLOWLIST, f"{tf.name}: {ref!r} not in AG_ALLOWLIST"


def test_allowlist_excludes_probability_paper():
    # Membership -- not the format regex -- is the guard: the probability paper
    # matches the arXiv format perfectly yet must be rejected.
    assert ARXIV_FORMAT.match(PROBABILITY_PAPER)     # format check passes...
    assert PROBABILITY_PAPER not in AG_ALLOWLIST      # ...membership rejects it.


def test_allowlist_is_frozen_and_wellformed():
    assert isinstance(AG_ALLOWLIST, frozenset)
    assert all(ARXIV_FORMAT.match(ref) for ref in AG_ALLOWLIST)


# --------------------------------------------------------------------------
# E7-M1 / G11: Enriques + bielliptic record rows and the canonical_order field.
#
# No new mathematics: chi(O) = 1 - q + p_g is re-derived per row and cross-checked
# against the primary sources cited in each Certificate.  Enriques: Nuer-Yoshioka
# arXiv:1901.04848, Yoshioka arXiv:1607.04946 (K is 2-torsion; q=0, p_g=0).
# Bielliptic (order-2 Bagnera-de Franchis, q=1, p_g=0): arXiv:2107.13370.
# --------------------------------------------------------------------------


def test_enriques_row():
    # chi(O) = 1 - q + p_g = 1 - 0 + 0 = 1; K_X 2-torsion -> K.H = 0 but order 2.
    s = enriques()
    assert s.chi_O == 1
    assert s.canonical_order == 2
    assert s.kind == "enriques"
    assert s.d == 2               # representative even polarization H^2 = 2
    assert s.K_H == 0             # K_X numerically trivial (2-torsion)
    assert s.picard_rank == 10    # Num = U + E8(-1)
    # Numerically K-trivial (K.H = 0) yet NOT genuinely trivial (order 2):
    assert s.trivial_canonical is True
    assert s.certificate.rigor is Rigor.PROVEN


def test_bielliptic_row():
    # chi(O) = 1 - q + p_g = 1 - 1 + 0 = 0; K_X torsion (order 2 here) -> K.H = 0.
    s = bielliptic()
    assert s.chi_O == 0
    assert s.canonical_order == 2
    assert s.kind == "bielliptic"
    assert s.d == 2
    assert s.K_H == 0
    assert s.picard_rank == 2
    assert s.trivial_canonical is True
    assert s.certificate.rigor is Rigor.PROVEN


def test_distinguished_from_k3_abelian():
    # The (chi_O, canonical_order) pair separates the three numerically-K-trivial
    # families: Enriques (1,2), K3 (2,0), abelian (0,0) -- all distinct.
    enr = enriques()
    k3 = K3()
    ab = abelian_surface()
    sigs = {
        "enriques": (enr.chi_O, enr.canonical_order),
        "k3": (k3.chi_O, k3.canonical_order),
        "abelian": (ab.chi_O, ab.canonical_order),
    }
    assert sigs["enriques"] == (1, 2)
    assert sigs["k3"] == (2, 0)
    assert sigs["abelian"] == (0, 0)
    assert len(set(sigs.values())) == 3
    # All three are numerically K-trivial, so K.H alone cannot separate them:
    assert enr.trivial_canonical and k3.trivial_canonical and ab.trivial_canonical


def test_faithful_computation_supported_flag():
    # False exactly for the torsion-canonical record rows (canonical_order != 0).
    assert enriques().faithful_computation_supported is False
    assert bielliptic().faithful_computation_supported is False
    # True for every non-torsion row, INCLUDING the rho>=2 rational ones: the flag
    # is keyed on canonical torsion, not Picard rank (see E7-M1 design note).
    for s in (P2, P1xP1, K3(), abelian_surface(), hirzebruch(2)):
        assert s.faithful_computation_supported is True


def test_certificates_cite_primary_sources():
    enr = enriques().certificate
    assert enr.rigor is Rigor.PROVEN
    assert "arXiv:1901.04848" in enr.citations   # Nuer-Yoshioka
    assert "arXiv:1607.04946" in enr.citations   # Yoshioka
    bie = bielliptic().certificate
    assert bie.rigor is Rigor.PROVEN
    assert "arXiv:2107.13370" in bie.citations   # Nuer, stable sheaves on bielliptic
    # Must NOT conflate with Koseki's separate weighted-hypersurface / Calabi-Yau
    # double/triple-solid paper arXiv:2007.00044 (E7-M2 scope).  (arXiv:1510.04474
    # is Chuang-Lai, WITHDRAWN -- not Koseki -- and likewise must never appear.)
    assert "arXiv:2007.00044" not in enr.citations
    assert "arXiv:2007.00044" not in bie.citations
    assert "arXiv:1510.04474" not in enr.citations
    assert "arXiv:1510.04474" not in bie.citations


def test_polarization_validation():
    # Enriques forces an even positive H^2 (the U + E8(-1) lattice is even).
    with pytest.raises(ValueError):
        enriques(3)
    with pytest.raises(ValueError):
        enriques(0)
    # A bielliptic polarization also forces an even positive H^2: chi(O)=0 with
    # numerically-trivial K gives chi(O(D)) = D^2/2 in Z (Riemann-Roch), so H^2
    # must be even -- H^2 = 1 cannot occur.
    with pytest.raises(ValueError):
        bielliptic(0)
    with pytest.raises(ValueError):
        bielliptic(1)
    with pytest.raises(ValueError):
        bielliptic(3)
    assert bielliptic(2).d == 2


# --------------------------------------------------------------------------
# E7-M2 / G11 [RESEARCH/UNVERIFIED]: Koseki product-threefold rows and the
# surface-consuming-API faithful-computation guard.
#
# No new mathematics.  The one asserted math value, chi(O)=0, is re-derived by
# Kunneth: chi(O_{XxY}) = chi(O_X)*chi(O_Y); chi(O_{P^1}) = chi(O_{P^2}) =
# chi(O_{P^1xP^1}) = 1, and an abelian surface (h^{0,1,2}=1,2,1 -> 0) or an
# elliptic curve (chi = 1 - g = 0) contributes 0, so all three products have
# chi(O)=0.  Koseki, "Stability conditions on threefolds with nef tangent
# bundles", arXiv:1811.03267 -- distinct from Koseki's separate Calabi-Yau
# double/triple-solid (weighted-hypersurface) paper arXiv:2007.00044 (must not
# be conflated).  (arXiv:1510.04474 is NOT Koseki: it is Chuang-Lai, WITHDRAWN.)
# --------------------------------------------------------------------------
KOSEKI_PAPER = "arXiv:1811.03267"
# Koseki's product-SPECIFIC primary source (title is exactly this row class):
# "Stability conditions on product threefolds of projective spaces and Abelian
# varieties" (Bull. LMS); KOSEKI_PAPER above is the later nef-tangent generalization.
KOSEKI_PRODUCT_PAPER = "arXiv:1703.07042"
# Koseki's SEPARATE Calabi-Yau double/triple-solid (weighted-hypersurface) BG
# paper -- a genuine Koseki paper distinct from the nef-tangent-bundle paper
# above; koseki_product must NOT cite it (the G11 conflation guard).  This is
# NOT arXiv:1510.04474 (Chuang-Lai, WITHDRAWN, not Koseki) -- that earlier ID
# was a miscitation.
KOSEKI_WEIGHTED_HYPERSURFACE_PAPER = "arXiv:2007.00044"
KOSEKI_TYPES = ("P1xS", "P2xC", "P1xP1xC")


def test_koseki_product_rows():
    # Each product row carries chi(O)=0 (Kunneth, re-derived above), the caller-
    # supplied H^3, kind "koseki-product", cites the product-specific
    # arXiv:1703.07042 then the nef-tangent generalization arXiv:1811.03267, and a
    # Certificate.  d3 is passed through verbatim (polarization-dependent).
    for t in KOSEKI_TYPES:
        r = koseki_product(t, 3)
        assert r.chi_O == 0                       # Kunneth: chi(O)=0 for all three
        assert r.d3 == 3                          # caller-supplied H^3 preserved
        assert r.kind == "koseki-product"
        assert r.references == [KOSEKI_PRODUCT_PAPER, KOSEKI_PAPER]
        assert r.certificate is not None
        assert KOSEKI_PAPER in r.certificate.citations
        assert KOSEKI_PRODUCT_PAPER in r.certificate.citations
    # d3 is genuinely caller-supplied (not hard-coded):
    assert koseki_product("P1xP1xC", 6).d3 == 6
    assert koseki_product("P1xS", 12).d3 == 12


def test_koseki_bg_gate():
    # [RESEARCH/UNVERIFIED] bg_proven and the Certificate rigor are gated on
    # KOSEKI_BG_VERIFIED: until Koseki's exact hypotheses are checked against
    # arXiv:1811.03267 the row is bg_proven=False / CONJECTURAL, NEVER PROVEN.
    assert KOSEKI_BG_VERIFIED is False           # the flag ships un-verified
    for t in KOSEKI_TYPES:
        r = koseki_product(t, 3)
        # bg_proven tracks the flag exactly (True IFF the hypothesis check is done):
        assert r.bg_proven is KOSEKI_BG_VERIFIED
        assert r.bg_proven is False
        # rigor is PROVEN IFF verified; while un-verified it is CONJECTURAL, never PROVEN:
        expected = Rigor.PROVEN if KOSEKI_BG_VERIFIED else Rigor.CONJECTURAL
        assert r.certificate.rigor is expected
        assert r.certificate.rigor is Rigor.CONJECTURAL
        assert r.certificate.rigor is not Rigor.PROVEN


def test_koseki_cites_correct_paper():
    # The note AND the Certificate cite arXiv:1811.03267 and NEVER Koseki's
    # separate weighted-hypersurface paper arXiv:2007.00044 (the conflation G11
    # warns of), nor the withdrawn Chuang-Lai paper arXiv:1510.04474.
    for t in KOSEKI_TYPES:
        r = koseki_product(t, 3)
        assert KOSEKI_PAPER in r.note
        assert KOSEKI_PAPER in r.certificate.citations
        assert KOSEKI_WEIGHTED_HYPERSURFACE_PAPER not in r.note
        assert KOSEKI_WEIGHTED_HYPERSURFACE_PAPER not in r.certificate.citations
        assert KOSEKI_WEIGHTED_HYPERSURFACE_PAPER not in r.references
        # The misattributed / withdrawn ID must never resurface either:
        assert "arXiv:1510.04474" not in r.note
        assert "arXiv:1510.04474" not in r.certificate.citations
        assert "arXiv:1510.04474" not in r.references
    # Provenance guard (G1): the correct paper is allowlisted, the wrong one is not.
    assert KOSEKI_PAPER in AG_ALLOWLIST
    assert KOSEKI_WEIGHTED_HYPERSURFACE_PAPER not in AG_ALLOWLIST


def test_koseki_validation():
    # Unknown product tag and non-positive H^3 both raise ValueError.
    with pytest.raises(ValueError):
        koseki_product("P3xC", 3)                # not a Koseki product tag
    with pytest.raises(ValueError):
        koseki_product("P2xC", 0)
    with pytest.raises(ValueError):
        koseki_product("P2xC", -1)


def test_faithful_guard_raises():
    # The guard lives in the Surface-consuming APIs and fires on the torsion-
    # canonical rows (Enriques / bielliptic, faithful_computation_supported is
    # False), raising the "NS-lattice refactor (G12)" NotImplementedError.
    v = ChernChar(1, 0, -2)
    w = ChernChar(1, -1, F(1, 2))
    for surf in (enriques(), bielliptic()):
        with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
            compute_walls(v, surf)
        with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
            actual_walls(v, surf)
        with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
            moduli_nonempty(Bundle.O(0), surface=surf)
    # The bare helper is also directly callable and refuses these rows:
    with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
        require_faithful_computation(enriques())
    # ...and is a no-op (returns None, does not raise) on every faithful row:
    for s in (P2, P1xP1, K3(), abelian_surface(), hirzebruch(2)):
        assert require_faithful_computation(s) is None
    # Defensive: the same-kind wrappers guard AFTER their kind check, so a hand-
    # built torsion-canonical class of the right kind is still refused (catalog
    # rows never trigger it -- they have canonical_order=0).
    tors_ab = replace(abelian_surface(), canonical_order=2)
    with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
        abelian_wall(v, w, tors_ab)
    with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
        check_existence_abelian(ChernChar(2, 2, 0), tors_ab)
    tors_k3 = replace(K3(), canonical_order=2)
    with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
        k3_wall_classified(v, w, tors_k3)
    with pytest.raises(NotImplementedError, match="NS-lattice refactor"):
        check_existence_k3(ChernChar(1, 0, 0), tors_k3)


def test_numerical_wall_not_guarded():
    # numerical_wall(v, w, d) takes only a ChernChar and an int d -- it never sees
    # a Surface, so it CANNOT and does NOT raise the G12 guard (the caveat is
    # docs-only).  It returns the exact shared fixture wall: for v=(1,0,-2),
    # w=(1,-1,1/2) on d=2, center -5/2 and radius^2 17/4 (re-derived: W_rc=-1,
    # W_re=5/2, W_ce=-2 => center=(5/2)/(-1)=-5/2, radius_sq=25/4-2*(-2)/(2*(-1))
    # =25/4-2=17/4; identical to test_abelian_k3_walls).
    v = ChernChar(1, 0, -2)
    w = ChernChar(1, -1, F(1, 2))
    wall = numerical_wall(v, w, 2)
    assert wall.center == F(-5, 2)
    assert wall.radius_sq == F(17, 4)


def test_guard_noop_on_faithful_consumers():
    # The guard is a pure no-op on faithful surfaces: the Surface-consuming APIs
    # still compute normally (no value changes).  P^2[2] keeps its pinned wall.
    v = ChernChar(1, 0, -2)
    walls = actual_walls(v, P2)
    assert len(walls) == 1
    outer = walls[0]
    assert outer.center == F(-5, 2)              # P^2[2] Gieseker wall (ABCH)
    assert outer.radius_sq == F(9, 4)            # radius 3/2
    # moduli_nonempty on the default P^2 surface is unaffected by the guard:
    res = moduli_nonempty(Bundle.O(0))
    assert res["nonempty"] is True
