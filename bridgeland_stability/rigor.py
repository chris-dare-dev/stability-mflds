"""E1-M2 / G5: the Rigor/Certificate provenance lattice.

A sheaf-free taint lattice: it records WHICH theorem is invoked and whether its
stated hypotheses are met; it never proves anything.  Import-time deps are
``enum``, ``dataclasses``, ``typing`` (standard library only) -- preserving the
zero-runtime-dependency invariant.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Tuple


class Rigor(IntEnum):
    """Provenance levels as a total order (PROVEN strongest; ``min`` = meet)."""

    PROVEN = 3        # a cited theorem covers these hypotheses
    CONJECTURAL = 2   # holds modulo a named open conjecture (e.g. threefold BG)
    HEURISTIC = 1     # numerical/empirical (dense compute_walls; doubling cert)
    UNKNOWN = 0       # untagged / no claim


@dataclass(frozen=True)
class Certificate:
    """A provenance stamp: a rigor level plus the hypotheses/citations backing it."""

    rigor: Rigor
    hypotheses: Tuple[str, ...] = ()
    citations: Tuple[str, ...] = ()
    note: str = ""


# Shared immutable default: an untagged certificate.
UNKNOWN_CERTIFICATE = Certificate(Rigor.UNKNOWN, (), (), "")


def _union(seq_of_tuples) -> Tuple[str, ...]:
    """Order-preserving union (first occurrence wins) of string tuples."""
    seen = {}
    for t in seq_of_tuples:
        for x in t:
            seen.setdefault(x, None)
    return tuple(seen)


def meet(*certs: "Certificate") -> "Certificate":
    """Combine by the meet (weakest link): rigor = min; hypotheses/citations
    are set-unions; notes are joined.  No args -> UNKNOWN_CERTIFICATE."""
    if not certs:
        return UNKNOWN_CERTIFICATE
    return Certificate(
        rigor=min(c.rigor for c in certs),
        hypotheses=_union(c.hypotheses for c in certs),
        citations=_union(c.citations for c in certs),
        note="; ".join(c.note for c in certs if c.note),
    )
