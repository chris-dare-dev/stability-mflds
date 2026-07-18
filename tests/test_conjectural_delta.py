"""E15-M4: the Conjecture A-gated evaluator of ``delta_m^{mu-s}``.

The conditional computation the paper promises under its Sec. 1.5 conjecture,
with the provenance lattice carrying the conditionality: CONJECTURAL wherever
the Kronecker part decides, PROVEN exactly where the certified-sharp DLP part
dominates (the anticanonical del Pezzo ray).  Differentialed against the
E14-M1 sandwich and the E14-M2 formula.
"""

from fractions import Fraction as F

from bridgeland_stability.conjectural_delta import delta_conjectural
from bridgeland_stability.delta_sharp import delta_kronecker, delta_mu_stable
from bridgeland_stability.rigor import Rigor
from bridgeland_stability.varieties import P1xP1, hirzebruch

F0 = P1xP1
F1 = hirzebruch(1)


def test_kronecker_grid_points_reproduce_the_closed_formula():
    # On the Sec. 18 grid the Kronecker part wins and equals thm-deltaKronecker
    # exactly; the rigor is honestly CONJECTURAL.
    for surf, nu, m, expect in (
            (F0, (F(1, 3), F(1, 5)), F(25, 9), F(3, 5)),
            (F0, (F(1, 3), F(1, 5)), F(5, 2), F(26, 45)),
            (F1, (F(6, 13), F(3, 13)), F(12, 7), F(98, 169)),
            (F1, (F(6, 13), F(3, 13)), F(3, 2), F(379, 676))):
        r = delta_conjectural(nu, m, surf)
        assert r.kronecker_part == expect
        assert r.value == expect
        assert r.value == delta_kronecker(nu, m, surf)
        assert r.dlp_part < expect                    # Kronecker strictly beats DLP here
        assert r.certificate.rigor == Rigor.CONJECTURAL


def test_anticanonical_ray_is_proven_and_kronecker_never_exceeds_it():
    # cor-deltaDLP: the sharp bound on -K IS the envelope; the Kronecker orbit
    # must not contribute above it (a value above the PROVEN sharp bound would
    # falsify the formula or the machinery).
    for surf, nu, m, expect in (
            (F0, (F(1, 2), F(1, 2)), F(1), F(3, 4)),
            (F0, (0, 0), F(1), F(1)),
            (F1, (0, 0), F(1, 2), F(1))):
        r = delta_conjectural(nu, m, surf)
        assert r.value == expect and r.dlp_part == expect
        assert r.kronecker_part is None or r.kronecker_part <= r.dlp_part
        assert r.certificate.rigor == Rigor.PROVEN


def test_sandwich_brackets_the_conjectural_value():
    # E14-M1: lower <= delta <= upper must hold for the TRUE delta; under the
    # conjecture the evaluator computes it, so the sandwich must bracket the
    # conjectural value too -- two independent routes.
    r = delta_conjectural((F(1, 3), F(1, 5)), F(25, 9), F0)
    s = delta_mu_stable((F(1, 3), F(1, 5)), F(25, 9), F0, max_rank=15)
    assert s.lower <= r.value < s.upper


def test_twist_invariance_of_the_conjectural_value():
    # delta is twist-invariant; the evaluator's orbit construction must respect
    # it on a shifted slope (shift within the searched twist box).
    base = delta_conjectural((F(1, 3), F(1, 5)), F(25, 9), F0)
    shifted = delta_conjectural((F(1, 3) + 1, F(1, 5) - 2), F(25, 9), F0)
    assert shifted.value == base.value


def test_conjectural_certificate_names_its_hypotheses():
    r = delta_conjectural((F(6, 13), F(3, 13)), F(12, 7), F1)
    assert r.certificate.rigor == Rigor.CONJECTURAL
    assert any("Conjecture A" in h for h in r.certificate.hypotheses)
    assert any("orbit suffices" in h for h in r.certificate.hypotheses)
