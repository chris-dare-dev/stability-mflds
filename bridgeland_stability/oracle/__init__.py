"""E10 / G16: optional Macaulay2 (M2) subprocess oracle -- NEVER imported by the core.

Mirrors the ``viz/`` lazy-backend culture: importing this subpackage pulls only
the standard library; the external Macaulay2 *binary* is discovered and spawned
lazily behind :func:`require_m2`, exactly as
:func:`bridgeland_stability.viz.style.require_mpl` gates matplotlib.  Because the
core package (:mod:`bridgeland_stability`) never imports this subpackage,
``import bridgeland_stability`` stays zero-third-party-dependency; and even
``import bridgeland_stability.oracle`` is stdlib-only -- only *calling* an oracle
function requires an M2 install the user provisions.

Sheaf-level bodies (chi/ext/moduli) land in E10-M2 / E10-M3; E10-M1 ships the
lazy skeleton, discovery + :func:`require_m2` guard, and :class:`M2Session`.
"""
from __future__ import annotations

from .m2 import (
    M2NotFoundError,
    M2Session,
    chi_via_ext,
    ext_dims,
    find_m2,
    moduli_nonempty_by_construction,
    require_m2,
)

__all__ = [
    "M2NotFoundError",
    "M2Session",
    "require_m2",
    "find_m2",
    "chi_via_ext",
    "ext_dims",
    "moduli_nonempty_by_construction",
]
