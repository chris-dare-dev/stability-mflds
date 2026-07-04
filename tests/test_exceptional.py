"""Algorithm 1: exceptional bundles on P^2 (corrected Markov/epsilon construction)."""

from fractions import Fraction as F

import pytest

from bridgeland_stability.exceptional import (
    Bundle,
    chi,
    enumerate_exceptional,
    exceptional_mediant,
    is_exceptional,
    is_markov_number,
    markov_numbers,
)


def test_line_bundles_are_exceptional_with_zero_discriminant():
    for n in range(-3, 4):
        O = Bundle.O(n)
        assert O.r == 1
        assert O.slope == F(n)
        assert O.discriminant == 0
        assert O.c2 == 0  # a line bundle has c2 = ch1^2/2 - ch2 = n^2/2 - n^2/2 = 0
        assert is_exceptional(O)


def test_T_minus1():
    T = Bundle.T_minus1()
    assert (T.r, T.c1, T.ch2) == (2, 1, F(-1, 2))
    assert T.slope == F(1, 2)
    assert T.discriminant == F(3, 8)  # CH convention
    assert T.discriminant_brief == F(3, 4)  # brief's doubled value
    assert T.c2 == 1
    assert chi(T, T) == 1
    assert is_exceptional(T)


def test_riemann_roch_offdiagonal():
    # chi(O, O(n)) = (n+1)(n+2)/2; independent validation of the RR implementation
    O = Bundle.O(0)
    assert chi(O, Bundle.O(1)) == 3
    assert chi(O, Bundle.O(3)) == 10
    assert chi(O, O) == 1
    assert chi(Bundle.O(1), O) == 0  # chi(O(1),O) = chi(O(-1)) = 0


def test_rank5_bundle_at_two_fifths_is_genuine():
    B = Bundle.from_slope(F(2, 5))
    assert (B.r, B.c1, B.ch2) == (5, 2, F(-2))
    assert B.discriminant == F(12, 25)  # (1/2)(1 - 1/25)
    assert B.c2 == 4  # integral -> genuine bundle
    assert is_exceptional(B)


def test_rank3_pseudobundle_does_not_exist():
    # The brief's "rank-3 exceptional at slope 1/3" is fictitious: c2 = 5/3 not integral.
    fake = Bundle.from_slope(F(1, 3))
    assert (fake.r, fake.c1, fake.ch2) == (3, 1, F(-7, 6))
    assert fake.c2 == F(5, 3)  # NOT an integer
    assert not is_exceptional(fake)
    # chi(E,E)=1 alone cannot detect it (it only sees (r,c1,ch2)):
    assert chi(fake, fake) == 1
    assert not is_markov_number(3)
    assert not is_markov_number(4)


def test_epsilon_mediant_recursion():
    assert exceptional_mediant(0, 1) == F(1, 2)  # -> T(-1), rank 2
    assert exceptional_mediant(0, F(1, 2)) == F(2, 5)  # -> rank 5
    assert exceptional_mediant(F(1, 2), 1) == F(3, 5)  # -> rank 5 (symmetric)
    assert exceptional_mediant(0, F(2, 5)) == F(5, 13)  # -> rank 13
    assert exceptional_mediant(F(2, 5), F(1, 2)) == F(12, 29)  # -> rank 29


def test_markov_numbers_set():
    assert markov_numbers(40) == {1, 2, 5, 13, 29, 34}
    assert 3 not in markov_numbers(100)
    assert 4 not in markov_numbers(100)
    assert {1, 2, 5, 13, 29, 34, 89, 169, 194, 233} <= markov_numbers(250)


def test_enumeration_only_produces_markov_ranks():
    bundles = enumerate_exceptional(0, 1, R_max=30)
    slopes = {b.slope for b in bundles}
    # the simplest exceptional slopes in [0,1]
    assert F(0) in slopes and F(1) in slopes and F(1, 2) in slopes
    assert F(2, 5) in slopes and F(3, 5) in slopes and F(5, 13) in slopes
    # 1/3, 2/3, 1/4 are NOT exceptional slopes
    assert F(1, 3) not in slopes
    assert F(2, 3) not in slopes
    assert F(1, 4) not in slopes
    # every enumerated bundle is genuine, with Markov rank
    mk = markov_numbers(30)
    for b in bundles:
        assert is_exceptional(b)
        assert b.r in mk
        assert b.discriminant == F(1, 2) * (1 - F(1, b.r * b.r))


def test_enumeration_sorted_and_bounded():
    bundles = enumerate_exceptional(0, 2, R_max=20)
    assert all(b.r <= 20 for b in bundles)
    slopes = [b.slope for b in bundles]
    assert slopes == sorted(slopes)
