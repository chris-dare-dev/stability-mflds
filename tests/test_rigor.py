"""E1-M2 / G5: the Rigor/Certificate provenance lattice.

These tests pin the lattice laws (total order, meet = min, set-union of
hypotheses/citations) and confirm that attaching a `certificate` field to the
catalog/result objects changes NO computed value: the P^2[2] actual wall stays
exactly `(center=-5/2, radius_sq=9/4)` (radius 3/2; ABCH, arXiv:1203.0316 sec.9),
recomputed by hand below with exact `Fraction` arithmetic.

Arithmetic anchor (exact `Fraction`), P^2[2]:
  v = (1, 0, -2),  w = O(-1) = (1, -1, 1/2),  d = 1.
    W_rc = r c' - r' c = 1*(-1) - 1*0   = -1
    W_re = r e' - r' e = 1*(1/2) - 1*(-2) = 5/2
    W_ce = c e' - c' e = 0*(1/2) - (-1)*(-2) = -2
    center    s0 = W_re / W_rc = (5/2)/(-1) = -5/2
    radius_sq    = s0^2 - 2 W_ce/(d W_rc) = 25/4 - 2*(-2)/(1*-1) = 25/4 - 4 = 9/4
So the certificate never perturbs the pinned wall.
"""

from fractions import Fraction as F

import pytest

from bridgeland_stability.rigor import (
    Rigor,
    Certificate,
    meet,
    UNKNOWN_CERTIFICATE,
)
from bridgeland_stability.chern import ChernChar
from bridgeland_stability.varieties import P2, BLOWUP_P3_POINT
from bridgeland_stability.walls import actual_walls, compute_walls


# --------------------------------------------------------------------------
# 1. The lattice: total order + meet = min
# --------------------------------------------------------------------------
def test_rigor_total_order():
    assert Rigor.PROVEN == 3
    assert Rigor.CONJECTURAL == 2
    assert Rigor.HEURISTIC == 1
    assert Rigor.UNKNOWN == 0
    assert Rigor.PROVEN > Rigor.CONJECTURAL > Rigor.HEURISTIC > Rigor.UNKNOWN


def test_meet_is_min():
    # min over IntEnum members returns the smallest-valued member (the meet).
    assert min(Rigor.PROVEN, Rigor.CONJECTURAL) == Rigor.CONJECTURAL
    proven = Certificate(Rigor.PROVEN, ("hyp-a",), ("arXiv:1203.0316",), "alg")
    conj = Certificate(Rigor.CONJECTURAL, ("hyp-b",), ("arXiv:1602.05055",), "input")
    m = meet(proven, conj)
    assert m.rigor == Rigor.CONJECTURAL
    assert m.rigor == min(proven.rigor, conj.rigor)


def test_meet_unions_hypotheses_and_citations():
    a = Certificate(Rigor.PROVEN, ("h1", "h2"), ("c1",), "na")
    b = Certificate(Rigor.HEURISTIC, ("h2", "h3"), ("c1", "c2"), "nb")
    m = meet(a, b)
    # order-preserving union, first occurrence wins, no duplicates
    assert m.hypotheses == ("h1", "h2", "h3")
    assert m.citations == ("c1", "c2")
    assert m.rigor == Rigor.HEURISTIC
    assert m.note == "na; nb"


def test_meet_no_args_is_unknown():
    assert meet() == UNKNOWN_CERTIFICATE
    assert meet().rigor == Rigor.UNKNOWN


def test_meet_skips_empty_notes():
    a = Certificate(Rigor.PROVEN, (), (), "")
    b = Certificate(Rigor.PROVEN, (), (), "only")
    assert meet(a, b).note == "only"


# --------------------------------------------------------------------------
# 2. Certificate is a frozen, hashable dataclass
# --------------------------------------------------------------------------
def test_certificate_is_frozen():
    c = Certificate(Rigor.PROVEN, ("h",), ("c",), "n")
    with pytest.raises(Exception):
        c.rigor = Rigor.UNKNOWN  # type: ignore[misc]
    # frozen -> hashable -> usable as a dataclass field default
    assert hash(c) == hash(Certificate(Rigor.PROVEN, ("h",), ("c",), "n"))


def test_certificate_defaults():
    c = Certificate(Rigor.UNKNOWN)
    assert c.hypotheses == ()
    assert c.citations == ()
    assert c.note == ""
    assert UNKNOWN_CERTIFICATE == Certificate(Rigor.UNKNOWN, (), (), "")


# --------------------------------------------------------------------------
# 3. Catalog tagging
# --------------------------------------------------------------------------
def test_p2_certificate_proven():
    assert P2.certificate.rigor == Rigor.PROVEN
    assert "arXiv:1203.0316" in P2.certificate.citations


def test_blowup_certificate_conjectural():
    assert BLOWUP_P3_POINT.certificate.rigor == Rigor.CONJECTURAL
    assert BLOWUP_P3_POINT.certificate.rigor <= Rigor.CONJECTURAL


def test_blowup_certificate_note_phrasing():
    note = BLOWUP_P3_POINT.certificate.note
    assert "strong BMT boundary not rigorous" in note
    assert "no stability condition" not in note.lower()


# --------------------------------------------------------------------------
# 4. Result tagging does NOT perturb the pinned P^2[2] wall
# --------------------------------------------------------------------------
def test_actual_wall_p2_value_and_proven_tag():
    v = ChernChar(1, 0, -2)
    walls = actual_walls(v, P2)
    assert len(walls) == 1  # the single ABCH wall realized by O(-1)
    w = walls[0]
    # certificate attached but value unchanged (center -5/2, radius_sq 9/4)
    assert w.center == F(-5, 2)
    assert w.radius_sq == F(9, 4)
    assert w.radius == pytest.approx(1.5)
    assert w.certificate.rigor == Rigor.PROVEN  # meet(PROVEN alg, PROVEN surface)


def test_actual_walls_noncoprime_p2_is_heuristic_not_proven():
    # Completeness of the actual-wall set is a THEOREM only for ABCH/Coskun-Huizenga
    # -covered P^2 classes (P^2[n] and coprime gcd(rank,degree)=1).  A general
    # non-coprime rank-2 class must NOT be over-tagged PROVEN -- the provenance
    # lattice must stay honest about proven-vs-empirical completeness.
    v = ChernChar(2, 0, -4)  # gcd(2, 0) = 2  ->  non-coprime, ABCH/CH-uncertified
    walls = actual_walls(v, P2)
    assert walls, "expected actual walls for (2,0,-4) on P^2"
    assert all(w.certificate.rigor == Rigor.HEURISTIC for w in walls)
    # contrast: the coprime P^2[2] class (gcd(1,0)=1) stays PROVEN, value unchanged
    p2_2 = actual_walls(ChernChar(1, 0, -2), P2)
    assert p2_2[0].certificate.rigor == Rigor.PROVEN
    assert (p2_2[0].center, p2_2[0].radius_sq) == (F(-5, 2), F(9, 4))


def test_compute_walls_tagged_heuristic():
    v = ChernChar(1, 0, -2)
    walls = compute_walls(v, P2)
    assert walls, "compute_walls should return a non-empty dense set for P^2[2]"
    # meet(HEURISTIC alg, PROVEN P2 surface) == HEURISTIC: never over-claims
    assert all(w.certificate.rigor == Rigor.HEURISTIC for w in walls)
    # dense set is a superset of the certified actual wall (center -5/2)
    assert any(w.center == F(-5, 2) and w.radius_sq == F(9, 4) for w in walls)


# --------------------------------------------------------------------------
# 5. Backward-compatible defaults + zero-dependency export
# --------------------------------------------------------------------------
def test_default_certificate_is_unknown():
    # An untagged wall (no certificate passed) defaults to UNKNOWN.
    from bridgeland_stability.walls import numerical_wall
    w = numerical_wall(ChernChar(1, 0, -2), ChernChar(1, -1, F(1, 2)), 1)
    assert w.certificate is UNKNOWN_CERTIFICATE
    assert w.certificate.rigor == Rigor.UNKNOWN


def test_lattice_is_exported_from_package():
    import bridgeland_stability as bs
    assert bs.Rigor is Rigor
    assert bs.Certificate is Certificate
    assert bs.meet is meet
    assert bs.UNKNOWN_CERTIFICATE is UNKNOWN_CERTIFICATE
    for name in ("Rigor", "Certificate", "meet", "UNKNOWN_CERTIFICATE"):
        assert name in bs.__all__
