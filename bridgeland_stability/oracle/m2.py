"""E10-M1 / G16: lazy Macaulay2 (M2) subprocess oracle -- NEVER imported by the core.

Importing this module pulls only the standard library (os, shutil, subprocess);
the external Macaulay2 binary is discovered and spawned lazily behind
:func:`require_m2`, exactly as :func:`bridgeland_stability.viz.style.require_mpl`
gates matplotlib.  So ``import bridgeland_stability`` (which never imports this
subpackage) and ``import bridgeland_stability.oracle[.m2]`` stay
zero-third-party-dependency; only *calling* an oracle function needs M2.

Host-path note [UNVERIFIED on Windows]: on Windows, Macaulay2 runs either
natively or under WSL; the discovered invocation (``BRIDGELAND_M2`` env override,
then ``shutil.which('M2')``) has NOT been verified on this host -- confirm the M2
call before relying on any oracle value.  (Goal G16; epic E10.)

E10-M2 (Euler-pairing cross-check) note: the *mathematics* of
:func:`chi_via_ext` / :func:`ext_dims` is **[PROVEN]** -- P^2 / K3 line-bundle
cohomology plus global ``Ext^i`` via truncated graded-module Ext (Greg Smith,
arXiv:math/9807170).  But the *exact M2 accessor* used to read the QQ-dimension
back (``rank Ext^i(F, G)`` printed by ``--script``) is **[UNVERIFIED]** against a
real Macaulay2 install: print-format drift is the named E10 danger-zone risk, so
the hand-maintained regex parser MUST be confirmed on a provisioned M2 before the
[PROVEN] tag is earned end-to-end.  On this host M2 is absent, so the cross-check
runs only where a user provisions M2 (the four E10-M2 tests skip otherwise).

E10-M3 (constructive non-emptiness) note: :func:`moduli_nonempty_by_construction`
is a **sufficient-only** witness (G16).  It CONSTRUCTS one explicit sheaf with
the requested numerical class and returns ``True`` when the build verifies, or
``None`` when the class is outside the constructible family or the build/verify
fails -- it NEVER returns ``False`` (a construction failure is not a proof of
emptiness).  The witness family is the rank-1 ideal-sheaf-of-points class on P^2,
``I_Z(c1)`` with length ``l = c1^2/2 - ch2 = c2 >= 0`` in ZZ: such a sheaf is
torsion-free of rank 1, hence mu-stable, so its existence genuinely proves the
moduli space is nonempty.  Higher rank / non-P^2 classes return ``None`` (out of
witness scope), never ``False``.  As with E10-M2 the M2 sheaf-level idiom is
**[UNVERIFIED]** on this host (M2 absent); the ``@requires_m2`` test exercises it
where a user provisions Macaulay2.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from fractions import Fraction
from typing import Optional, Tuple

__all__ = [
    "M2NotFoundError",
    "find_m2",
    "require_m2",
    "M2Session",
    "chi_via_ext",
    "ext_dims",
    "moduli_nonempty_by_construction",
]

#: Environment variables consulted (in order) to locate M2 before a PATH lookup.
_M2_ENV_VARS = ("BRIDGELAND_M2", "M2_PATH")


class M2NotFoundError(RuntimeError):
    """Raised when an oracle op needs Macaulay2 but no M2 binary is found.

    Deliberately a subprocess/environment guard (not ``ImportError``): M2 is an
    external CAS binary, not a Python wheel.  It mirrors the *intent* of
    :func:`viz.style.require_mpl` -- one clear, actionable guard -- while being
    honest that the missing dependency is a program on ``PATH``.
    """


def find_m2() -> Optional[str]:
    """Return the M2 executable path, or ``None`` if absent.  Never raises/spawns.

    Order: ``BRIDGELAND_M2`` / ``M2_PATH`` env overrides (first non-empty wins),
    then ``shutil.which('M2')`` / ``shutil.which('m2')``.
    """
    for var in _M2_ENV_VARS:
        val = os.environ.get(var)
        if val:
            return val
    return shutil.which("M2") or shutil.which("m2")


def require_m2() -> str:
    """Return the M2 path or raise a clear :class:`M2NotFoundError`.

    The oracle analogue of :func:`viz.style.require_mpl`: every oracle entry
    point calls this first, so with M2 absent the failure is a single helpful
    guard, not an opaque ``FileNotFoundError`` from ``subprocess``.
    """
    path = find_m2()
    if path is None:
        raise M2NotFoundError(
            "Macaulay2 (M2) is required for the sheaf-level oracle but was not "
            "found. Install Macaulay2 (https://macaulay2.com) and put 'M2' on "
            "PATH, or set BRIDGELAND_M2 to its path. The core package does not "
            "need M2; only oracle functions do."
        )
    return path


class M2Session:
    """Context-managed Macaulay2 ``--script`` subprocess, gated by :func:`require_m2`.

    Constructing the object is cheap and import-safe; the external process is
    spawned only inside a ``with`` block, after :func:`require_m2` (which raises
    :class:`M2NotFoundError` when M2 is absent).  The concrete script protocol is
    finalized in E10-M2; E10-M1 pins only the lazy spawn-behind-guard contract.
    """

    def __init__(self, executable: Optional[str] = None, *, timeout: float = 60.0):
        self._executable = executable
        self._timeout = timeout
        self._proc: Optional[subprocess.Popen] = None

    def __enter__(self) -> "M2Session":
        exe = self._executable or require_m2()  # raises M2NotFoundError if absent
        self._executable = exe
        self._proc = subprocess.Popen(   # noqa: S603 -- trusted, user-provisioned exe
            [exe, "--script"],           # concrete args refined in E10-M2
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        proc, self._proc = self._proc, None
        if proc is not None and proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=self._timeout)
            except subprocess.TimeoutExpired:  # pragma: no cover
                proc.kill()
        return None


# --- E10-M2 batch runner + script generators (stdlib-only at import) ---------
def _run_m2(code: str, *, executable: Optional[str] = None, timeout: float = 60.0) -> str:
    """Run one M2 script string in a fresh subprocess; return stdout. Guarded by require_m2().

    One-shot batch per query is a deliberate choice: it is the robust shape for
    the "hand-maintained text parser" the E10 danger-zone warns about, and it
    does not require finalizing an interactive stdin REPL (that is E10-M3
    territory).  :class:`M2Session` is left exactly as E10-M1 shipped.
    """
    exe = executable or require_m2()
    with tempfile.NamedTemporaryFile("w", suffix=".m2", delete=False) as fh:
        fh.write(code)
        path = fh.name
    try:
        proc = subprocess.run(  # noqa: S603 -- trusted, user-provisioned exe
            [exe, "--script", path],
            capture_output=True,
            text=True,
            encoding="utf-8",        # locale-independent QQ/ZZ round-trip (E10-M2 nit)
            errors="strict",
            timeout=timeout,
        )
    finally:
        os.unlink(path)
    if proc.returncode != 0:
        raise M2NotFoundError(  # reuse the oracle guard family; carries M2's stderr
            "Macaulay2 script failed (returncode %d):\n%s"
            % (proc.returncode, proc.stderr or proc.stdout)
        )
    return proc.stdout


def _variety_m2(X) -> str:
    """Return M2 code defining ``X = Proj(...)`` for a supported Surface.

    [UNVERIFIED idiom] -- the ``Varieties``/``Proj`` construction is standard M2
    but the exact printed output is not yet confirmed on a real install.
    """
    if getattr(X, "is_p2", False) or getattr(X, "kind", "") == "P2":
        return "S = QQ[x_0,x_1,x_2];\nX = Proj S;\n"
    if getattr(X, "kind", "") == "K3":
        # smooth Fermat quartic in P^3: a K3 (d=4).  H^i(O_X)=(1,0,1) for ANY K3,
        # so the Picard rank is irrelevant to this cohomology witness.
        return (
            "S = QQ[x_0,x_1,x_2,x_3];\n"
            "X = Proj(S/ideal(x_0^4+x_1^4+x_2^4+x_3^4));\n"
        )
    raise NotImplementedError(
        "oracle Ext cross-check supports P^2 and K3 (quartic) only; got %r" % (X,)
    )


def _line_twist(B) -> int:
    """Twist ``n`` of a rank-1 line bundle ``O(n)``; raise for higher-rank input."""
    if getattr(B, "r", None) != 1:
        raise NotImplementedError(
            "chi_via_ext/ext_dims support line bundles (rank 1) only in E10-M2"
        )
    return int(B.c1)


def _qq_roundtrip(x: Fraction) -> Fraction:
    """Send ``x`` to M2 as a ``QQ`` literal, read the echoed ``toString`` back,
    parse losslessly to :class:`~fractions.Fraction`.

    Both ``Fraction("-7/6")`` and ``Fraction("4")`` parse exactly, so a rational
    survives the ``Fraction -> QQ -> text -> Fraction`` round-trip bit-for-bit
    (G16: "M2 rationals are QQ; Fraction('p/q') <-> p/q round-trips exactly").
    """
    require_m2()
    code = 'v = (%d)/%d;\n<< "QQ " << toString v << endl;\n' % (
        x.numerator,
        x.denominator,
    )
    out = _run_m2(code)
    for line in out.splitlines():
        mo = re.match(r"^QQ\s+(\S+)\s*$", line.strip())
        if mo:
            return Fraction(mo.group(1))
    raise M2NotFoundError("could not parse QQ round-trip from M2 output:\n" + out)


# --- E10-M2 sheaf-level Euler-pairing cross-check ----------------------------
def ext_dims(E, F, X) -> Tuple[int, ...]:
    """``(ext^0, ext^1, ext^2) = dim_QQ Ext^i_X(E, F)`` via M2.

    Mathematics [PROVEN]: P^2 / K3 line-bundle cohomology + global ``Ext^i`` via
    truncated graded-module Ext (Greg Smith, arXiv:math/9807170).  The exact M2
    accessor ``rank Ext^i(F, G)`` used to read the QQ-dimension is [UNVERIFIED]
    on a real M2 install (print-format drift; see the module docstring).
    """
    require_m2()  # FIRST -- preserves the E10-M1 None-arg guard test
    a, b = _line_twist(E), _line_twist(F)
    code = _variety_m2(X) + (
        "F = OO_X(%d);\nG = OO_X(%d);\n"
        'scan(0..2, i -> ( m := Ext^i(F, G); << "EXT " << i << " " << (rank m) << endl; ));\n'
        '<< "OK" << endl;\n' % (a, b)
    )
    out = _run_m2(code)
    dims = {}
    for line in out.splitlines():
        mo = re.match(r"^EXT (\d+) (\d+)\s*$", line.strip())
        if mo:
            dims[int(mo.group(1))] = int(mo.group(2))
    return tuple(dims[i] for i in range(3))


def chi_via_ext(E, F, X) -> int:
    """``chi(E, F) = sum (-1)^i ext^i(E, F)``, M2-computed.

    A second, independent verification path diffed against the formula layer
    (:func:`bridgeland_stability.exceptional.chi` on P^2,
    :func:`bridgeland_stability.mukai.pairing` ``= -chi`` on K3).  No new core
    mathematics -- every shipped value stays byte-for-byte unchanged.
    """
    require_m2()
    dims = ext_dims(E, F, X)
    return sum((-1) ** i * d for i, d in enumerate(dims))


# --- E10-M3 constructive (sufficient-only) non-emptiness witness -------------
def _rank1_ideal_length(r: int, c1, ch2) -> Optional[int]:
    """Length l (= c2) of the ideal-sheaf witness I_Z(c1) realizing the rank-1
    class (1, c1, ch2) on P^2, or None if this class is not constructible by the
    ideal-sheaf-of-points routine.

    ch(I_Z(c1)) = (1, c1, c1^2/2 - l), so l = c1^2/2 - ch2 = c2.  Returns l only
    when r == 1, c1 is an integer, and l is a NON-NEGATIVE INTEGER (the P^2
    Chern-lattice integrality c2 in Z plus effectivity l >= 0); otherwise None.
    All arithmetic is exact Fraction -- never int '/' (which is float in Py3).
    """
    if r != 1:
        return None
    try:
        c1f = Fraction(c1)
    except (TypeError, ValueError):
        return None
    if c1f.denominator != 1:
        return None
    ell = Fraction(c1f * c1f, 2) - Fraction(ch2)   # = c2
    if ell < 0 or ell.denominator != 1:
        return None
    return int(ell)


def _moduli_witness_m2(X, c1: int, ell: int) -> str:
    """M2 code building I_Z(c1) for a length-`ell` 0-dim scheme Z on P^2, then
    echoing the constructed sheaf's (r, c1, ch2) reconstructed from three Euler
    characteristics chi(E(0)), chi(E(1)), chi(E(2)) via P^2 Riemann-Roch.

    [UNVERIFIED idiom] -- P^2/Proj construction is standard M2, but the exact
    printed output is not confirmed on a real install (see module docstring).
    Only P^2 is supported; callers gate on X being P^2 first.
    """
    if not (getattr(X, "is_p2", False) or getattr(X, "kind", "") == "P2"):
        raise NotImplementedError(
            "moduli_nonempty_by_construction witness supports P^2 only; got %r" % (X,)
        )
    return (
        "S = QQ[x_0,x_1,x_2];\nX = Proj S;\n"
        "IZ = sheaf module ideal(x_1, x_2^%d);\n"   # colength-ell scheme at [1:0:0]
        "E = IZ ** OO_X(%d);\n"
        'scan({0,1,2}, t -> ( << "CHI " << t << " " << euler(E(t)) << endl; ));\n'
        '<< "OK" << endl;\n' % (ell, c1)
    )


def _parse_witness(out: str) -> Optional[Tuple[int, int, Fraction]]:
    """Parse three 'CHI t k' lines into (r, c1, ch2) via P^2 Riemann-Roch, or None.

    With chi(E(t)) = (r/2) t^2 + (c1 + 3r/2) t + (ch2 + 3/2 c1 + r), the three
    integer Euler characteristics k_t = chi(E(t)) invert to
        r   = k_0 - 2 k_1 + k_2         (second difference = r),
        c1  = (k_1 - k_0) - 2 r         (chi(E(1)) - chi(E(0)) = c1 + 2r),
        ch2 = k_0 - 3/2 c1 - r          (chi(E(0)) = ch2 + 3/2 c1 + r).
    Returns None unless exactly the three CHI lines parsed and the reconstructed
    class is rank 1 with integral c1 (the ideal-sheaf witness family).
    """
    chi = {}
    for line in out.splitlines():
        mo = re.match(r"^CHI (\d+) (-?\d+)\s*$", line.strip())
        if mo:
            chi[int(mo.group(1))] = int(mo.group(2))
    if set(chi) != {0, 1, 2}:
        return None
    k0, k1, k2 = chi[0], chi[1], chi[2]
    r = k0 - 2 * k1 + k2                          # second difference = r
    c1 = Fraction(k1 - k0 - 2 * r)               # from chi(E(1)) - chi(E(0))
    ch2 = Fraction(k0) - Fraction(3, 2) * c1 - r  # chi(E(0)) = ch2 + 3/2 c1 + r
    if r != 1 or c1.denominator != 1:
        return None
    return (int(r), int(c1), ch2)


def moduli_nonempty_by_construction(
    r: int, c1, ch2: Fraction, X
) -> Optional[bool]:
    """SUFFICIENT-ONLY non-emptiness witness (E10-M3 / G16).

    Attempts to CONSTRUCT one Gieseker-semistable sheaf with numerical class
    (r, c1, ch2) on X in Macaulay2.  Returns:
      * True  -- a witness was built and its Chern character verified to match
                 (moduli space M_H(r,c1,ch2) is NONEMPTY);
      * None  -- the construction did not yield a matching witness (class not in
                 the constructible family, or M2 build/verify failed).
    A construction failure is NOT a proof of emptiness, so this function NEVER
    returns False (G16 risk note: 'only a sufficient witness -- label it so').

    Constructible family (E10-M3): rank-1 ideal-sheaf-of-points classes on P^2,
    I_Z(c1) with length l = c1^2/2 - ch2 = c2 >= 0 in Z.  Such a sheaf is
    torsion-free rank 1, hence mu-stable, so its existence is a genuine
    sufficient witness.  Higher rank / non-P^2 return None (out of witness
    scope), never False.  The M2 sheaf-level idiom is [UNVERIFIED] on this host
    (M2 absent); the @requires_m2 test exercises it where a user provisions M2.
    """
    require_m2()  # FIRST -- preserves the E10-M1 None-arg guard test
    ell = _rank1_ideal_length(r, c1, ch2)
    if ell is None:
        return None
    if not (getattr(X, "is_p2", False) or getattr(X, "kind", "") == "P2"):
        return None
    out = _run_m2(_moduli_witness_m2(X, int(c1), ell))
    got = _parse_witness(out)
    if got == (1, int(c1), Fraction(ch2)):
        return True
    return None
