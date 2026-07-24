"""Integrity gate for the independent oracle (E12-M0).

These tests protect the properties that make ``tests/oracle`` a trustworthy gate:
  * the reference implementation never imports the package it tests;
  * the reference implementation stays exact (no float, no square root);
  * the corpus is self-consistent with the reference (both transcribe the theorem);
  * a literal frozen map of every corpus row is unchanged -- appending rows is free,
    but mutating a frozen row's verdict fails here and needs a docs/CORRECTIONS.md
    entry to land (also enforced by .githooks/pre-commit).

This file imports only the oracle, never the package under test; the
package-importing differential sweep lives in test_differential.py.
"""

import os
import sys
from fractions import Fraction

# The oracle lives at tests/oracle/; put tests/ on sys.path so `import oracle`
# resolves regardless of pytest's import mode.
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _TESTS_DIR not in sys.path:
    sys.path.insert(0, _TESTS_DIR)

from oracle.dlp_reference import Status, reference_nonempty, reference_delta
from oracle.corpus import CORPUS, DELTAS

_ORACLE_DIR = os.path.join(_TESTS_DIR, "oracle")
_REFERENCE_SRC = os.path.join(_ORACLE_DIR, "dlp_reference.py")


def _read(path):
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def test_reference_has_no_package_import():
    """The charter's load-bearing rule: the oracle names no package symbol."""
    src = _read(_REFERENCE_SRC)
    assert "bridgeland_stability" not in src


def test_reference_uses_no_float():
    """The oracle is exact-integer/Fraction only; the rank cutoff is a denominator
    bound, never a floating-point square root."""
    src = _read(_REFERENCE_SRC)
    for token in ("float(", "import math", "from math", "sqrt"):
        assert token not in src, (
            "reference must stay exact/integer-only; found forbidden token "
            f"{token!r}"
        )


def test_corpus_agrees_with_reference():
    """Every corpus row reproduces under the reference (self-consistency)."""
    for row in CORPUS:
        assert row.surface == "P^2", row
        got = reference_nonempty(row.r, row.c1, row.ch2)
        assert got is row.expected, (row, "got", got)
    for row in DELTAS:
        got = reference_delta(row.mu)
        assert got == row.expected_delta, (row, "got", got)


# --------------------------------------------------------------------------- #
# Frozen map of every corpus verdict known at freeze time (E12-M0).           #
# Copy of the corpus table, as literals.  Appending new corpus rows is free;  #
# changing a verdict below (or in the corpus) fails test_frozen_corpus_unchanged
# and requires a docs/CORRECTIONS.md entry.                                    #
# --------------------------------------------------------------------------- #
FROZEN_STATUS = {
    ("P^2", 4, 2, Fraction(-1)): Status.NONEMPTY,        # A1
    ("P^2", 2, 0, Fraction(0)): Status.NONEMPTY,         # A1
    ("P^2", 10, 3, Fraction(-9, 2)): Status.EMPTY,       # A2
    ("P^2", 610, 133, Fraction(-581, 2)): Status.EMPTY,  # A2
    ("P^2", 1, 0, Fraction(-3, 2)): Status.INVALID,      # A3
    ("P^2", 8010, 3060, Fraction(-3421)): Status.EMPTY,  # A4
    ("P^2", 2, 1, Fraction(-5, 2)): Status.NONEMPTY,
    ("P^2", 2, 1, Fraction(-1, 2)): Status.NONEMPTY,
    ("P^2", 5, 2, Fraction(-2)): Status.NONEMPTY,
    ("P^2", 1, 0, Fraction(0)): Status.NONEMPTY,
    ("P^2", 2, 1, Fraction(1, 2)): Status.EMPTY,
    ("P^2", 2, 1, Fraction(0)): Status.INVALID,
    ("P^2", 1, 0, Fraction(-1)): Status.NONEMPTY,
    ("P^2", 1, 0, Fraction(-5)): Status.NONEMPTY,
}


def test_frozen_corpus_unchanged():
    """Every frozen row is still present in the corpus with the same verdict."""
    index = {(row.surface, row.r, row.c1, row.ch2): row for row in CORPUS}
    # No two corpus rows may collide on the (surface, r, c1, ch2) key.
    assert len(index) == len(CORPUS), "duplicate (surface, r, c1, ch2) key in corpus"
    for key, status in FROZEN_STATUS.items():
        assert key in index, f"frozen corpus row vanished: {key}"
        assert index[key].expected is status, (
            f"frozen corpus row mutated: {key} expected {status} "
            f"but corpus says {index[key].expected}"
        )


def test_oracle_dir_holds_exactly_the_expected_modules():
    """Make an accidental extra file under tests/oracle/ visible."""
    names = sorted(n for n in os.listdir(_ORACLE_DIR) if n.endswith(".py"))
    assert names == ["__init__.py", "corpus.py", "delpezzo_mutation_reference.py",
                     "dlp_reference.py", "mutation_reference.py"], names
