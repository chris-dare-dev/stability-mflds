"""E13-M4: the K-theoretic mutation oracle vs the epsilon-recursion (common-mode gate).

The E13 adversarial re-audit found the E12 differential oracle import-independent
but not ALGORITHMICALLY independent: it generates exceptional slopes by the same
epsilon-mediant subdivision production uses.  This module closes that common mode
by triangulating THREE independent generators of the same finite sets:

  1. the production epsilon-recursion (Drezet-Le Potier Theoreme A; slope
     mediants seeded on the integers) -- ``bridgeland_stability.exceptional``;
  2. the K-theoretic MUTATION walk (Gorodentsev-Rudakov constructibility;
     integer linear algebra on (r, c1, chi), Markov-equation invariant) --
     ``tests/oracle/mutation_reference.py``;
  3. Springborn's MARKOV FRACTIONS (Veselov, arXiv:2501.06779, Thm 3.1: the
     slopes of the exceptional bundles are exactly the Markov fractions) --
     the purely number-theoretic mediant  p1/q1 * p2/q2 = (p1q1 + p2q2)/(q1^2
     + q2^2)  seeded on 0/1 and 1/2, transcribed INLINE below.

No pair of these shares a formula, a coordinate system, or a traversal.  The
delta-curve differential additionally runs production's certified sharp cutoff
(rank <= denominator(mu)) against the oracle's deliberately larger ``4q + 64``
cutoff, so the sharp-cutoff theorem itself (docs/CORRECTIONS.md Sec. 8) is no
longer a shared assumption.  See docs/CORRECTIONS.md Sec. 13.
"""

import os
import sys
from fractions import Fraction as F
from math import gcd

import pytest

from bridgeland_stability.dlp import delta as production_delta
from bridgeland_stability.exceptional import (
    Bundle,
    certified_rank_cutoff,
    enumerate_exceptional,
    is_exceptional,
)
from bridgeland_stability.dlp_hirzebruch import total_slope, discriminant
from bridgeland_stability.exceptional_surface import SurfaceBundle
from bridgeland_stability.hn_filtration import semistable_exists
from bridgeland_stability.varieties import P2

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _TESTS_DIR not in sys.path:
    sys.path.insert(0, _TESTS_DIR)

from oracle.corpus import CORPUS
from oracle.dlp_reference import reference_semistable_exists
from oracle.mutation_reference import (
    mutation_collections,
    mutation_delta,
    mutation_euler,
    mutation_exceptional_classes,
    mutation_exceptional_in_window,
    mutation_is_exceptional,
    mutation_is_semiexceptional,
    mutation_semistable_exists,
)

# --------------------------------------------------------------------------- #
# Literature constants.                                                        #
# --------------------------------------------------------------------------- #
#: Markov numbers <= 1000, hardcoded from OEIS A002559 (NOT computed by the
#: package): the possible ranks of exceptional bundles on P^2 (Rudakov).
MARKOV_A002559 = (1, 2, 5, 13, 29, 34, 89, 169, 194, 233, 433, 610, 985)

#: The rank-<= 610 mutation-generated base set, computed once (fast: the walk
#: visits ~64 canonical collections).
_BASE_610 = mutation_exceptional_classes(610)


def springborn_markov_fractions(q_max):
    """Markov fractions in [0, 1/2] with denominator <= ``q_max``.

    Transcribed from Veselov, arXiv:2501.06779 (after B. Springborn): seeds
    0/1 and 1/2; between adjacent fractions insert the SPRINGBORN mediant

        p1/q1 * p2/q2 = (p1 q1 + p2 q2) / (q1^2 + q2^2)

    (``Fraction`` reduces it to lowest terms, matching the paper's reduced
    form).  Denominators strictly increase down the tree, so the bound prunes
    exactly.  Every Markov fraction p/q satisfies p^2 + 1 = 0 (mod q) -- a
    NECESSARY condition asserted here, and NOT sufficient (see the impostor
    test below).
    """
    seeds = (F(0), F(1, 2))
    out = set(seeds)
    stack = [seeds]
    while stack:
        a, b = stack.pop()
        g = F(a.numerator * a.denominator + b.numerator * b.denominator,
              a.denominator ** 2 + b.denominator ** 2)
        if g.denominator > q_max:
            continue
        out.add(g)
        stack.append((a, g))
        stack.append((g, b))
    for f in out:
        assert (f.numerator ** 2 + 1) % f.denominator == 0, f
    return out


def _exceptional_ch2(r, c1):
    """The ch2 forced by chi(E,E)=1 at (r, c1): (c1^2 - r^2 + 1)/(2r)."""
    return F(c1 * c1 - r * r + 1, 2 * r)


# --------------------------------------------------------------------------- #
# 1. Charter: the second oracle is import-independent and exact.               #
# --------------------------------------------------------------------------- #
def _mutation_src():
    path = os.path.join(_TESTS_DIR, "oracle", "mutation_reference.py")
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def test_mutation_reference_has_no_package_import():
    src = _mutation_src()
    assert "bridgeland_stability" not in src.replace(
        "``bridgeland_stability.exceptional``", "")   # docstring mention only
    # the load-bearing form: no import statement can name the package
    for line in src.splitlines():
        stripped = line.strip()
        if stripped.startswith(("import ", "from ")):
            assert "bridgeland" not in stripped, stripped


def test_mutation_reference_uses_no_float():
    src = _mutation_src()
    for token in ("float(", "import math", "from math", "sqrt"):
        assert token not in src, f"forbidden token {token!r}"


def test_mutation_reference_shares_no_epsilon_machinery():
    """The independence charter, checked structurally: the mutation oracle never
    mentions the epsilon-mediant operation or its Delta_alpha rank formula."""
    src = _mutation_src()
    assert "(3 + a - b)" not in src and "(3 + alpha - beta)" not in src
    assert "exceptional_mediant" not in src


# --------------------------------------------------------------------------- #
# 2. The derived closed forms, pinned against classical values.                #
# --------------------------------------------------------------------------- #
def test_euler_form_pins_to_line_bundle_cohomology():
    """chi(O(a), O(b)) = (b-a+2)(b-a+1)/2 -- the classical chi of O(b-a) on P^2,
    checked over a grid.  This anchors the all-integer Euler form independently
    of the package's Fraction-based RR pairing."""
    def On(n):
        return (1, n, (n + 2) * (n + 1) // 2)
    for a in range(-6, 7):
        for b in range(-6, 7):
            m = b - a
            assert mutation_euler(On(a), On(b)) == (m + 2) * (m + 1) // 2, (a, b)


def test_euler_form_agrees_with_package_chi_on_generated_classes():
    """Cross-check the two independently-derived Euler forms on real classes."""
    from bridgeland_stability.exceptional import chi as package_chi
    classes = sorted(_BASE_610)[:12]
    for (r1, c1, ch21) in classes:
        for (r2, c2, ch22) in classes:
            x1 = int(ch21 + F(3 * c1, 2) + r1)
            x2 = int(ch22 + F(3 * c2, 2) + r2)
            assert mutation_euler((r1, c1, x1), (r2, c2, x2)) == package_chi(
                Bundle(r1, c1, ch21), Bundle(r2, c2, ch22))


# --------------------------------------------------------------------------- #
# 3. Internal structure of the generated set ([R], Rudakov).                   #
# --------------------------------------------------------------------------- #
def test_collections_satisfy_the_markov_equation():
    """Every visited full collection has Markov rank triple a^2+b^2+c^2 = 3abc.
    (The walk asserts this internally too; this re-checks it from the outside.)"""
    cols = mutation_collections(233)
    assert len(cols) >= 30
    for col in cols:
        a, b, c = (X[0] for X in col)
        assert a * a + b * b + c * c == 3 * a * b * c, (a, b, c)


def test_generated_ranks_are_exactly_the_markov_numbers():
    """Ranks of the generated classes = OEIS A002559 up to 610, with the right
    multiplicity in [0,1): one class each for ranks 1 and 2, two for the rest."""
    by_rank = {}
    for (r, c1, ch2) in _BASE_610:
        by_rank.setdefault(r, set()).add((c1, ch2))
    assert sorted(by_rank) == [m for m in MARKOV_A002559 if m <= 610]
    assert len(by_rank[1]) == 1 and len(by_rank[2]) == 1
    for m in sorted(by_rank):
        if m > 2:
            assert len(by_rank[m]) == 2, (m, by_rank[m])


def test_generated_classes_have_unit_self_pairing_and_primitive_c1():
    for (r, c1, ch2) in _BASE_610:
        x = int(ch2 + F(3 * c1, 2) + r)
        assert mutation_euler((r, c1, x), (r, c1, x)) == 1
        assert gcd(r, c1) == 1 or r == 1
        c2 = F(c1 * c1, 2) - ch2
        assert c2.denominator == 1
        # Delta = (1 - 1/r^2)/2 emerges from chi(X,X)=1; it is never assumed
        # by the oracle, so assert it here as the classical cross-check.
        mu = F(c1, r)
        assert mu * mu / 2 - ch2 / r == (1 - F(1, r * r)) / 2


# --------------------------------------------------------------------------- #
# 4. THE TRIANGULATION: three independent recursions, one set.                 #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("R", [13, 89, 610])
def test_mutation_set_equals_production_epsilon_set(R):
    """Generator 2 (mutations) == generator 1 (production epsilon-recursion),
    as FULL (r, c1, ch2) triples with slope in [0, 1)."""
    mset = {t for t in mutation_exceptional_classes(R)}
    prod = {(b.r, b.c1, b.ch2) for b in enumerate_exceptional(0, 1, R)
            if 0 <= F(b.c1, b.r) < 1}
    assert mset == prod


def test_mutation_window_equals_production_window():
    mset = mutation_exceptional_in_window(-3, 4, 89)
    prod = {(b.r, b.c1, b.ch2) for b in enumerate_exceptional(-3, 4, 89)}
    assert mset == prod


def test_springborn_markov_fractions_are_the_slopes():
    """Generator 3 (Springborn/Veselov Markov fractions) == the slope sets of
    generators 1 and 2 in [0, 1): the fractions in [0, 1/2] plus their mirrors
    1 - f (the affine action's other representative in the unit interval)."""
    fr = springborn_markov_fractions(610)
    expected = {f for f in fr} | {1 - f for f in fr if f != 0}
    mut_slopes = {F(c1, r) for (r, c1, ch2) in _BASE_610}
    prod_slopes = {b.slope for b in enumerate_exceptional(0, 1, 610)
                   if 0 <= b.slope < 1}
    assert mut_slopes == expected
    assert prod_slopes == expected
    # the explicit small fractions of arXiv:2501.06779:
    for f in (F(0), F(1, 2), F(2, 5), F(5, 13), F(12, 29), F(13, 34),
              F(34, 89), F(70, 169), F(75, 194), F(89, 233), F(233, 610)):
        assert f in fr


def test_impostor_survives_the_congruence_but_no_recursion():
    """(610, 133, -581/2) -- the docs/CORRECTIONS.md Sec. 8 impostor -- passes
    chi(E,E)=1, integral c2, Markov rank, AND Springborn's necessary congruence
    p^2 = -1 (mod q) (133^2 + 1 = 29 * 610): only the tree structure kills it.
    All three generators must reject it; the genuine rank-610 slopes are
    233/610 and 377/610 (Fibonacci branch F(2k-1)/F(2k+1))."""
    assert (133 ** 2 + 1) % 610 == 0                      # congruence-proof
    assert F(133, 610) not in springborn_markov_fractions(610)
    assert not mutation_is_exceptional(610, 133, F(-581, 2))
    assert not is_exceptional(Bundle(610, 133, F(-581, 2)))
    for p in (233, 377):
        ch2 = _exceptional_ch2(610, p)
        assert (610, p, ch2) in _BASE_610
        assert is_exceptional(Bundle(610, p, ch2))
    # Fibonacci branch F(2k-1)/F(2k+1): 2/5, 5/13, 13/34, 34/89, 89/233, 233/610
    fib = [1, 1]
    while fib[-1] < 700:
        fib.append(fib[-1] + fib[-2])
    pairs = [(fib[i], fib[i + 2]) for i in range(2, len(fib) - 2, 2)
             if fib[i + 2] <= 610]
    assert pairs[0] == (2, 5) and pairs[-1] == (233, 610)
    for p, q in pairs:
        assert (q, p, _exceptional_ch2(q, p)) in _BASE_610


def test_impostor_family_sweep_agrees_with_production():
    """Every 'potentially exceptional' candidate (Markov rank q, exceptional ch2,
    integral c2) up to rank 610: production membership == mutation membership.
    This is the family the E12 audit's impostors live in; the sweep must
    actually contain rejected candidates for the gate to have teeth."""
    rejected = 0
    for q in (m for m in MARKOV_A002559 if m <= 610):
        for p in range(q):
            if gcd(p, q) != 1 and q != 1:
                continue
            ch2 = _exceptional_ch2(q, p)
            c2 = F(p * p, 2) - ch2
            if c2.denominator != 1:
                continue
            in_mutation = (q, p, ch2) in _BASE_610
            in_production = is_exceptional(Bundle(q, p, ch2))
            assert in_mutation == in_production, (q, p)
            if not in_mutation:
                rejected += 1
    assert rejected > 0, "sweep contained no impostor: the gate is vacuous"


# --------------------------------------------------------------------------- #
# 5. The delta curve, cross-cutoff.                                            #
# --------------------------------------------------------------------------- #
def test_delta_pinned_literature_values():
    """CLAUDE.md invariant-4 pins, reproduced by the mutation oracle alone."""
    assert mutation_delta(F(1, 2)) == F(5, 8)
    assert mutation_delta(F(1, 3)) == F(5, 9)
    assert mutation_delta(F(1, 4)) == F(21, 32)
    assert mutation_delta(F(2, 5)) == F(13, 25)
    assert mutation_delta(0) == F(1)


def test_delta_differential_with_disjoint_cutoffs():
    """production delta (certified sharp cutoff rank <= q) == mutation delta
    (cutoff 4q + 64) over a dense low-denominator sweep + high-q spot checks.
    Agreement across DIFFERENT cutoffs also tests the sharp-cutoff theorem
    (docs/CORRECTIONS.md Sec. 8): if any bundle of rank in (q, 4q+64] were
    binding, this would fail."""
    def prod(mu):
        mu = F(mu)
        R = certified_rank_cutoff(mu)
        return production_delta(mu, enumerate_exceptional(mu - 3, mu + 3, R))

    for q in range(1, 33):
        for p in range(0, q):
            if gcd(p, q) != 1 and q != 1:
                continue
            mu = F(p, q)
            assert mutation_delta(mu) == prod(mu), mu
    # high-denominator spot checks, including near an exceptional cusp
    for mu in (F(233, 610), F(89, 233), F(55, 144), F(100, 199), F(355, 113)):
        assert mutation_delta(mu) == prod(mu), mu


# --------------------------------------------------------------------------- #
# 6. The full verdict: triple differential.                                    #
# --------------------------------------------------------------------------- #
def _production_exists(r, c1, ch2):
    xi = SurfaceBundle(r, (c1,), ch2)
    return semistable_exists(r, total_slope(xi), discriminant(xi, P2), P2)


def test_semistable_triple_differential_on_the_frozen_corpus():
    """production == dlp_reference (epsilon) == mutation_reference (K-theory)
    on the E12 falsification corpus."""
    for row in CORPUS:
        assert row.surface == "P^2", row
        a = _production_exists(row.r, row.c1, row.ch2)
        b = reference_semistable_exists(row.r, row.c1, row.ch2)
        c = mutation_semistable_exists(row.r, row.c1, row.ch2)
        assert a is b is c, (row.r, row.c1, row.ch2, a, b, c)


def test_semistable_triple_differential_over_a_box():
    """The same triple agreement over an integral-c2 box sweep (the E12/E13
    audit box: r <= 6, c1 in [-8, 8], c2 in [0, 7])."""
    for r in range(1, 7):
        for c1 in range(-8, 9):
            for c2 in range(0, 8):
                ch2 = F(c1 * c1, 2) - c2
                a = _production_exists(r, c1, ch2)
                b = reference_semistable_exists(r, c1, ch2)
                c = mutation_semistable_exists(r, c1, ch2)
                assert a is b is c, (r, c1, ch2, a, b, c)


def test_mutation_membership_smoke():
    """Direct membership checks, including twisted spellings."""
    assert mutation_is_exceptional(1, 0, F(0))                # O
    assert mutation_is_exceptional(1, 7, F(49, 2))            # O(7)
    assert mutation_is_exceptional(2, 1, F(-1, 2))            # T(-1)
    assert mutation_is_exceptional(5, 7, F(5, 2))             # slope 7/5 = 2/5 + 1
    assert not mutation_is_exceptional(3, 1, F(-7, 6))        # the brief's fake rank 3
    assert not mutation_is_exceptional(2, 1, F(0))            # wrong ch2
    assert mutation_is_semiexceptional(4, 2, F(-1))           # T(-1)^{+2}
    assert not mutation_is_semiexceptional(4, 2, F(0))
