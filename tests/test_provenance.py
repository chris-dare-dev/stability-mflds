"""E12-M5 (A12, A13): provenance / citation / docstring regression guards.

Pure-stdlib, offline, deterministic -- no network, no arXMCP call, no Macaulay2.
These tests do NOT flip any oracle xfail (A12/A13 are provenance defects the
P^2-only reference oracle cannot adjudicate; they carry no strict-xfail tripwire
in ``tests/test_differential.py``).  They pin the corrected provenance strings so
the two wrong citations and the one false capability claim cannot silently return.

Modelled on ``tests/test_exceptional_surface.py::test_citation_provenance`` -- a
static-docstring guard against a miscited-arXiv defect class.
"""

from fractions import Fraction as F

import bridgeland_stability.exceptional_surface as es
import bridgeland_stability.nonemptiness_rational as nr
from bridgeland_stability.nonemptiness_rational import (
    HNMode,
    VerdictStatus,
    _MODE_CERT,
    moduli_nonempty,
    paper_delta_H_targets,
)
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.varieties import P2, P1xP1, K3, abelian_surface


def test_invalid_character_is_proven_empty_on_every_surface():
    # E12 code review: an invalid character (non-integral c2) is not the Chern
    # character of ANY sheaf, so M(xi) is empty on EVERY surface for EVERY
    # polarization -- a K_X-independent theorem.  The verdict must be PROVEN_EMPTY
    # uniformly.  Before the review the K3/abelian path (mode=HEURISTIC) under-
    # claimed it as UNKNOWN while P^2/F_e reported PROVEN_EMPTY.
    #
    # (2,(1),1/2): c2 = 1/2*<(1),(1)> - 1/2.  On K3(2)/abelian(2) the rank-1 shim
    # has <(1),(1)> = d = 2, so c2 = 1 - 1/2 = 1/2, non-integral.  On P^2 (d=1),
    # <(1),(1)> = 1, c2 = 1/2 - 1/2 = 0 -- integral, so use ch2=1/2 there too but
    # a genuinely non-integral case: (1,(0,),1/2) has c2 = -1/2 on P^2.
    cases = [
        (2, (1,), F(1, 2), K3(2)),
        (2, (1,), F(1, 2), abelian_surface(2)),
        (1, (0,), F(1, 2), P2),           # c2 = 0 - 1/2 = -1/2
        (2, (1, 0), F(1, 3), P1xP1),      # non-integral c2 on F_0
    ]
    for r, c1, ch2, surface in cases:
        v = moduli_nonempty(r, c1, ch2, surface)
        assert v.nonempty is False
        assert v.status is VerdictStatus.PROVEN_EMPTY, (
            f"{surface.name} {(r, c1, ch2)}: invalid character must be PROVEN_EMPTY, "
            f"got {v.status} (rigor {v.certificate.rigor.name})")


def _P(m: F) -> F:
    """Drezet-Le Potier ``P(m) = (m^2 + 3m + 2)/2`` in exact ``Fraction``."""
    return (m * m + 3 * m + 2) / F(2)


def test_delta_third_note_not_fictitious_rank3():
    # A12 (row note): delta(1/3) = 5/9 is attained by O (rank 1, slope 0) via
    # P(-1/3) - Delta_O, NOT by a "rank-3 exceptional bundle" -- an object that
    # does not exist (rank 3 is not Markov; Bundle.from_slope(1/3) has c2=5/3).
    # The two derivations agreed only by the numerical coincidence P(-1/3)=1-4/9.
    row = next(e for e in paper_delta_H_targets()
               if e.surface.is_p2 and e.r == 3 and e.c1 == (1,))
    assert row.delta_H == F(5, 9)                       # value unchanged
    assert "Delta_{rk3 exc}" not in row.note            # fictitious derivation removed
    assert "attained by O" in row.note
    assert "not Markov" in row.note
    # and the arithmetic the note asserts is actually true:
    assert _P(F(-1, 3)) - F(0) == F(5, 9)               # P(-1/3) - Delta_O = 5/9


def test_rank5_citation_not_gorodentsev_misquote():
    # A12 (citation): the rank-5 slope-2/5 exceptional bundle EXISTS by classical
    # Drezet-Le Potier (Thm A); CH arXiv:1907.06739 Cor 9.13 is delta_H = DLP_{-K}
    # on the del Pezzo F_e (e in {0,1}) -- it does not state that exceptional
    # bundles are -K-stable (that is attributed to Gorodentsev).
    row = next(e for e in paper_delta_H_targets()
               if e.surface.is_p2 and e.r == 5 and e.c1 == (2,))
    assert row.delta_H == F(13, 25)                     # value unchanged
    assert "exceptional bundles are -K-stable" not in row.citation   # debunked misquote gone
    assert "Gorodentsev" in row.citation                             # correct attribution
    assert "Ann. Sci. ENS 18 (1985)" in row.citation                # DLP existence citation
    assert "arXiv:19" in row.citation                                # keeps test_paper_table invariant
    assert _P(F(0)) - F(12, 25) == F(13, 25)            # delta(2/5) = P(0) - Delta_{r5} = 13/25


def test_paper_targets_relabelled_regression_fixture():
    # A12 (relabel): the table is a REGRESSION FIXTURE of hand-derived targets,
    # not a verbatim paper table.
    doc = (nr.paper_delta_H_targets.__doc__ or "")
    assert ("not a verbatim paper table" in doc.lower()
            or "regression fixture" in doc.lower())


def test_oracle_mode_certificate_describes_ideal_sheaf():
    # A13: _MODE_CERT[ORACLE] must describe what oracle/m2.py actually does --
    # construct a rank-1 ideal sheaf I_Z(c1) on P^2 and return True | None --
    # not claim a prioritary-sheaf HN filtration no code computes.  The dict is
    # asserted directly (M2 only gates minting, not the frozen certificate).
    cert = _MODE_CERT[HNMode.ORACLE]
    blob = " ".join(cert.hypotheses) + " " + cert.note
    assert "prioritary-sheaf HN filtration" not in blob   # false capability claim removed
    assert "ideal sheaf" in blob.lower()
    assert "I_Z" in blob
    assert ("never False" in blob or "True or None" in blob or "True|None" in blob)
    assert cert.rigor is Rigor.PROVEN                     # a real construction IS a proof of nonemptiness
    # Cross-check the string against the code it describes:
    import bridgeland_stability.oracle.m2 as m2
    assert (m2.moduli_nonempty_by_construction.__doc__
            and "NEVER" in m2.moduli_nonempty_by_construction.__doc__)


def test_canonical_arxiv_ids_resolve():
    # Retrieval regression (offline form): each of the four epic-canonical arXiv
    # IDs is cited somewhere in package source, and the two A12 debunked pairings
    # never reappear.  This is the falsifiable, network-free half of "the IDs
    # resolve"; live arXMCP resolution belongs to the adversarial source-checker
    # lens (invariant 3: the suite stays zero-dependency and offline).
    CANON = {
        "arXiv:1401.1613": ("Coskun-Huizenga-Woolf", "effective cone"),
        "arXiv:1907.06739": ("Coskun-Huizenga", "Existence of semistable sheaves"),
        "arXiv:1910.14060": ("Levine-Zhang", None),
        "arXiv:1611.02674": ("Coskun-Huizenga", "Weak Brill-Noether"),
    }
    surface = (nr.__doc__ or "") + "".join(str(getattr(nr, n)) for n in ("_CH",))
    # gather citation surface from the table + module docstrings:
    surface += "".join(e.citation + " " + e.note for e in nr.paper_delta_H_targets())
    surface += "".join(" ".join(c.citations) for c in nr._MODE_CERT.values())
    surface += (es.__doc__ or "")
    for aid in CANON:
        assert aid in surface, aid            # every canonical ID is cited in package source
    # the two A12 debunked pairings must never reappear:
    assert "birational geometry" not in surface.lower()          # 1611.02674 mis-pairing (already guarded)
    assert "exceptional bundles are -K-stable" not in surface    # 1907.06739 Cor 9.13 misquote
