# Independent test oracle package (E12-M0).
#
# This package is the independent gate for EPIC E12. Nothing in it may import
# the package under test, and its rows are frozen (see test_oracle_integrity.py
# and .githooks/pre-commit). It exists so that a bug in the library cannot hide
# behind a reference implementation derived from that same library.
