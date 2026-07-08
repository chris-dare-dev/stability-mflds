"""E1-M5 / CH provenance: dedicated tests for the stdlib-only ``_latex`` helpers.

These pin the exact LaTeX strings emitted by ``bridgeland_stability._latex`` and
guard the import-time purity invariant (the module, and the ``viz.style`` layer
that re-exports it, must not pull matplotlib/plotly at import time).

Display-defect fix (E1-M5, Option A)
------------------------------------
The shipped ``latex_sqrt`` returned malformed LaTeX for a bare surd whose radical
argument is an integer::

    latex_sqrt(3) -> '\\frac{\\sqrt{3}}{1}'   # BUG: \\frac{...}{1}
    latex_sqrt(2) -> '\\frac{\\sqrt{2}}{1}'   # same bug

Trace (``_latex.py`` second branch): for ``latex_sqrt(3)`` we have ``n=3, d=1``;
``_isqrt_exact(3) is None`` (not a perfect square) and ``_isqrt_exact(1) == 1``,
so ``rn=None, rd=1``.  The perfect-square branch is skipped and the ``sqrt(n)/rd``
branch is entered with ``rd == 1`` -- but the missing ``rd == 1`` case produced
``\\frac{\\sqrt{3}}{1}`` instead of the bare surd.  E1-M5 adds the ``rd == 1``
special case so a denominator of 1 renders as ``\\sqrt{n}``.  This is a
string-only, display-layer change: no ``Fraction`` / computed invariant moves,
and the function's own docstring already promised ``3 -> \\sqrt{3}``.

Literature anchors
------------------
* ``\\sqrt{3}`` is the display of the pinned P^3 null-correlation bundle result
  ``alpha_crit(beta=1/2) = sqrt(3)`` (docs/CORRECTIONS.md section 5).
* ``v=(2, 0, -1/4)`` is a *generic* Chern-tuple formatting example -- it mirrors
  the ``e.g.`` in ``_latex.latex_chern``'s own docstring and exercises the
  integer / zero / proper-fraction render paths.  It is **not** the pinned
  P^2[2] wall class and anchors no specific wall: it is not even an integral
  Chern character (``c_2 = ch_1^2/2 - ch_2 = 1/4``).  The genuine pinned P^2[2]
  wall is the ideal sheaf ``ch(I_Z) = (1, 0, -2)`` destabilized by
  ``O(-1) = (1, -1, 1/2)``, giving center ``-5/2`` and radius ``3/2``
  (ABCH, arXiv:1203.0316) -- that class is not the tuple displayed here.
"""

import subprocess
import sys
from fractions import Fraction as F

from bridgeland_stability import _latex as L


def test_latex_frac():
    assert L.latex_frac(F(-4, 3)) == r"-\frac{4}{3}"   # pinned roadmap value
    assert L.latex_frac(0) == "0"
    assert L.latex_frac(-3) == "-3"                     # integers render bare
    assert L.latex_frac(F(-1, 4)) == r"-\frac{1}{4}"   # sign pulled to front
    assert L.latex_frac(F(5, 8)) == r"\frac{5}{8}"


def test_latex_sqrt():
    assert L.latex_sqrt(3) == r"\sqrt{3}"               # pinned; P^3 alpha_crit(1/2)=sqrt3
    assert L.latex_sqrt(2) == r"\sqrt{2}"
    assert L.latex_sqrt(4) == "2"                        # perfect square
    assert L.latex_sqrt(F(9, 4)) == r"\frac{3}{2}"     # perfect-square ratio
    assert L.latex_sqrt(F(3, 4)) == r"\frac{\sqrt{3}}{2}"


def test_latex_chern():
    # Generic Chern-tuple formatting example (int, zero, proper fraction),
    # matching the ``e.g.`` in _latex.latex_chern's docstring.  Note (2,0,-1/4)
    # is NOT an integral Chern character (c_2 = 1/4) and anchors no wall; this
    # is a pure string-format check.
    assert L.latex_chern(2, 0, F(-1, 4)) == r"v=(r,c,\mathrm{ch}_2)=(2,\,0,\,-\frac{1}{4})"


def test_latex_is_stdlib_only():
    child = r'''
import sys
BLOCK = {"matplotlib", "plotly"}
class _Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in BLOCK:
            raise ImportError("blocked for import-purity test: " + name)
        return None
sys.meta_path.insert(0, _Blocker())
import bridgeland_stability._latex          # stdlib-only at import time
import bridgeland_stability.viz.style       # lazy backends -> also clean
assert "matplotlib" not in sys.modules, "matplotlib imported at import time"
assert "plotly" not in sys.modules, "plotly imported at import time"
print("STDLIB_ONLY_OK")
'''
    r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "STDLIB_ONLY_OK" in r.stdout
