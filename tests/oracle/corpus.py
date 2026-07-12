"""Falsification corpus for the P^2 non-emptiness oracle (E12-M0).

APPEND-ONLY.  Each row carries, without exception:
  * the surface (always ``P^2`` -- ``dlp_reference`` is P^2-only);
  * the character (r, c_1, ch_2);
  * the expected verdict, transcribed from the theorem, NOT from the library;
  * a primary-source citation naming the numbered statement it rests on;
  * the exact Fraction arithmetic that carries the class to the verdict.

A row that cannot be justified that way does not belong here.  ``test_oracle_integrity``
asserts the corpus is self-consistent with ``dlp_reference`` and that a literal frozen
map of the rows below is unchanged; appending is free, mutating a frozen row's verdict
fails the suite and requires a ``docs/CORRECTIONS.md`` entry (enforced by pre-commit).

CITATIONS
  [CHW 2.2]  Coskun-Huizenga-Woolf, arXiv:1401.1613, Sec.2, Thm 2.2 (the P^2
             non-emptiness iff:  c_1 in Z, chi in Z, Delta >= delta(mu)).
  [CHW exc]  same, Sec.2: exceptional bundles are the stable E with Delta(E) < 1/2;
             their moduli spaces are a single reduced point.
  [CH 1.9]   Coskun-Huizenga, arXiv:1907.06739, Ex.1.9 (non-exceptional stable
             sheaves exist iff Delta >= delta(mu)).
  [CH 1.14]  same, Ex.1.14 (a semiexceptional bundle = a direct sum of copies of
             an exceptional bundle).
  [DLP]      the Drezet-Le Potier delta curve; P(m)=(m^2+3m+2)/2,
             Delta_alpha=(1-1/rho^2)/2, delta(alpha)=P(0)-Delta_alpha at an
             exceptional slope alpha of rank rho.

Verdicts on P^2 use Delta = (1/2)mu^2 - ch_2/r (H^2 = 1), mu = c_1/r,
chi = ch_2 + (3/2)c_1 + r, and delta the DLP curve.
"""

from dataclasses import dataclass
from fractions import Fraction

# Local import of the enum so the corpus speaks the oracle's vocabulary.
from .dlp_reference import Status


@dataclass(frozen=True)
class CorpusRow:
    surface: str
    r: int
    c1: int
    ch2: Fraction
    expected: Status
    defect: str = ""        # "" or "A1".."A6": the audit defect this row pins
    citation: str = ""      # numbered primary statement
    derivation: str = ""    # exact Fraction arithmetic, class -> verdict


@dataclass(frozen=True)
class DeltaRow:
    mu: Fraction
    expected_delta: Fraction
    citation: str = ""
    derivation: str = ""


# --------------------------------------------------------------------------- #
# Non-emptiness verdicts.                                                      #
# --------------------------------------------------------------------------- #
NONEMPTY_CORPUS = (
    # ---- Character-decidable defects (A1..A4). ----
    CorpusRow(
        "P^2", 4, 2, Fraction(-1), Status.NONEMPTY, "A1",
        "[CH 1.14] + [CHW 2.2]",
        "mu=1/2; Delta=(1/2)(1/2)^2-(-1)/4=1/8+1/4=3/8 < delta(1/2)=5/8, so the "
        "curve clause fails; but xi = 2*ch(T(-1)) with T(-1)=(2,1,-1/2) exceptional "
        "(Delta_{T(-1)}=3/8=(1-1/4)/2), so xi is semiexceptional => NONEMPTY.",
    ),
    CorpusRow(
        "P^2", 2, 0, Fraction(0), Status.NONEMPTY, "A1",
        "[CH 1.14] + [CHW 2.2]",
        "mu=0; Delta=0 < delta(0)=1, curve clause fails; xi = 2*ch(O) with "
        "O=(1,0,0) exceptional, so xi = O^{+2} is semiexceptional => NONEMPTY.",
    ),
    CorpusRow(
        "P^2", 10, 3, Fraction(-9, 2), Status.EMPTY, "A2",
        "[CHW 2.2] + [CHW exc]",
        "mu=3/10 (rank 10, and 10 is NOT a Markov number); chi=-9/2+9/2+10=10 in Z "
        "and c_2=9/2+9/2=9 in Z; chi(E,E)=1 since Delta=(1/2)(3/10)^2+9/20=9/200+90/200"
        "=99/200=(1-1/100)/2. But 3/10 is not in the epsilon-recursion image, so it is "
        "NOT exceptional; Delta=99/200 < delta(3/10)=119/200 and it is not "
        "semiexceptional => EMPTY.  (The library accepts chi(E,E)=1 as exceptionality.)",
    ),
    CorpusRow(
        "P^2", 610, 133, Fraction(-581, 2), Status.EMPTY, "A2",
        "[CHW 2.2] + [CHW exc]",
        "mu=133/610; 610=2*5*61 IS a Markov number and ch_2=(133^2-610^2+1)/(2*610)"
        "=-581/2 matches the exceptional formula, so chi(E,E)=1 and "
        "Delta=(1-1/610^2)/2=372099/744200. But the epsilon-slopes of denominator 610 "
        "are 233/610 and 377/610; 133/610 is an impostor, hence NOT exceptional. "
        "Delta < 1/2 <= delta(133/610) and it is not semiexceptional => EMPTY.",
    ),
    CorpusRow(
        "P^2", 1, 0, Fraction(-3, 2), Status.INVALID, "A3",
        "[CHW 2.2]",
        "chi = -3/2 + (3/2)*0 + 1 = -1/2 is not in Z (equivalently c_2=3/2 not in Z), "
        "so (1,0,-3/2) is not the Chern character of any sheaf => INVALID.",
    ),
    CorpusRow(
        "P^2", 8010, 3060, Fraction(-3421), Status.EMPTY, "A4",
        "[CHW 2.2] + [CHW exc]",
        "mu=3060/8010=34/89; Delta=(1/2)(34/89)^2+3421/8010=356489/712890 and "
        "delta(34/89)=3961/7921=356490/712890, so Delta = delta - 1/712890 < delta. "
        "The only exceptional bundle of slope 34/89 is E=(89,34,-38); 90*E=(8010,3060,-3420) "
        "has ch_2=-3420 != -3421, so xi is neither exceptional nor semiexceptional => EMPTY. "
        "(The library's R_max=60 truncation misses the rank-89 cusp and reports NONEMPTY.)",
    ),

    # ---- Canonical DLP / non-emptiness anchors (no defect). ----
    CorpusRow(
        "P^2", 2, 1, Fraction(-5, 2), Status.NONEMPTY, "",
        "[CHW 2.2]",
        "mu=1/2; Delta=(1/2)(1/4)+5/4=1/8+10/8=11/8 >= delta(1/2)=5/8 => NONEMPTY.",
    ),
    CorpusRow(
        "P^2", 2, 1, Fraction(-1, 2), Status.NONEMPTY, "",
        "[CHW exc]",
        "T(-1): Delta=1/8+1/4=3/8 < delta(1/2)=5/8, but 1/2 is an epsilon-slope of "
        "rank 2 and ch_2=(1-4+1)/4=-1/2 matches, so T(-1) is exceptional (m=1) => NONEMPTY.",
    ),
    CorpusRow(
        "P^2", 5, 2, Fraction(-2), Status.NONEMPTY, "",
        "[CHW exc]",
        "mu=2/5 (epsilon-slope, rank 5); ch_2=(4-25+1)/10=-2 matches; "
        "Delta=(1/2)(4/25)+2/5=2/25+10/25=12/25 < delta(2/5)=13/25, exceptional => NONEMPTY.",
    ),
    CorpusRow(
        "P^2", 1, 0, Fraction(0), Status.NONEMPTY, "",
        "[CHW exc]",
        "O_{P^2}: Delta=0 < delta(0)=1, but 0 is an integer (rank-1) epsilon-slope and "
        "ch_2=0 matches, so O is exceptional => NONEMPTY.",
    ),
    CorpusRow(
        "P^2", 2, 1, Fraction(1, 2), Status.EMPTY, "",
        "[CHW 2.2]",
        "mu=1/2; c_2=1/2-1/2=0 in Z and chi=1/2+3/2+2=4 in Z, so the character is valid; "
        "Delta=(1/2)(1/4)-1/4=-1/8 < delta(1/2)=5/8, and chi(E,E) != 1 "
        "(Delta != 3/8) so not exceptional and not semiexceptional => EMPTY.",
    ),
    CorpusRow(
        "P^2", 2, 1, Fraction(0), Status.INVALID, "",
        "[CHW 2.2]",
        "chi = 0 + (3/2)(1) + 2 = 7/2 not in Z (c_2 = 1/2 not in Z) => INVALID.",
    ),
    CorpusRow(
        "P^2", 1, 0, Fraction(-1), Status.NONEMPTY, "",
        "[CHW 2.2]",
        "P^2[1] ideal sheaf: Delta=(1/2)(0)-(-1)/1=1 >= delta(0)=1 (boundary is "
        "inclusive on P^2) => NONEMPTY.",
    ),
    CorpusRow(
        "P^2", 1, 0, Fraction(-5), Status.NONEMPTY, "",
        "[CHW 2.2]",
        "P^2[5]: Delta=5 >= delta(0)=1 => NONEMPTY.",
    ),
)


# --------------------------------------------------------------------------- #
# The DLP delta curve at pinned slopes (CLAUDE.md invariant 4).               #
# --------------------------------------------------------------------------- #
DELTA_CORPUS = (
    DeltaRow(
        Fraction(1, 2), Fraction(5, 8), "[DLP]",
        "alpha=1/2 (rank 2), Delta_alpha=3/8; delta=P(0)-3/8=1-3/8=5/8.",
    ),
    DeltaRow(
        Fraction(1, 3), Fraction(5, 9), "[DLP]",
        "3 is not Markov, so there is no cusp at 1/3; the controlling slope is "
        "alpha=0 (rank 1) at distance 1/3: delta=P(-1/3)-0=(1/9-1+2)/2=5/9.",
    ),
    DeltaRow(
        Fraction(1, 4), Fraction(21, 32), "[DLP]",
        "4 is not Markov; controlling slope alpha=0 at distance 1/4: "
        "delta=P(-1/4)=(1/16-3/4+2)/2=(21/16)/2=21/32.",
    ),
    DeltaRow(
        Fraction(2, 5), Fraction(13, 25), "[DLP]",
        "alpha=2/5 (rank 5), Delta_alpha=12/25; delta=P(0)-12/25=1-12/25=13/25.",
    ),
    DeltaRow(
        Fraction(0), Fraction(1), "[DLP]",
        "alpha=0 (rank 1), Delta_alpha=0; delta=P(0)-0=1.",
    ),
    DeltaRow(
        Fraction(34, 89), Fraction(3961, 7921), "[DLP]",
        "34/89 is an epsilon-slope of rank 89 (89 Markov); "
        "delta=1-Delta_{34/89}=1-(1-1/89^2)/2=1-3960/7921=3961/7921.",
    ),
)


# --------------------------------------------------------------------------- #
# Out-of-scope defects, stated not guessed.                                    #
# --------------------------------------------------------------------------- #
#
# A5 (forgeable certificate, fixed in E12-M4) and A6 (F_e ch_2 guard, fixed in
# E12-M3) are NOT P^2 class->verdict rows -- A5 is a certificate-provenance defect
# and A6 lives on P^1xP^1 (F_0), which this P^2-only oracle cannot adjudicate.  Per
# the charter they are pinned by their own xfail tests in test_differential.py with
# primary citations, not invented as P^2 rows here.
#
# A7..A13 are metadata/provenance defects (moduli_dim doubled; Surface.K_H wrong on
# F_n; hirzebruch_with_polarization truncation; line_bundle fractional divisors;
# canonical_class keyed on the Gram matrix; the M4 table mislabel; the ORACLE mode
# certificate).  They are not class->verdict facts either; the fixing milestones are
# E12-M5 (A12, A13) and E12-M6 (A7..A11).  They are recorded here, not faked into rows.
#
# The CLAUDE.md P^2[n] GIESEKER WALL (center -(2n+1)/2, radius (2n-1)/2) is a
# `walls`-module fact, not a DLP delta / non-emptiness fact; this oracle computes no
# walls, so that value belongs to a future walls oracle, not to this corpus.  The
# P^2[n] NON-EMPTINESS rows (1,0,-n) above ARE in scope and are included.
METADATA_DEFECTS = ("A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13")


# Public exports.
CORPUS = NONEMPTY_CORPUS
DELTAS = DELTA_CORPUS
