"""E2-M4 smoke tests for the abelian-vs-K3 worked example.

Loads the non-package example script via ``importlib`` (it lives under
``examples/``, not the installed package) and asserts the two printed radius^2
values are the pinned exact fractions, plus that the script runs to a clean exit.
All asserted numbers are already ground-truth-pinned by
``tests/test_abelian_k3_walls.py`` and ``tests/test_mukai.py``; this module just
guards that the worked example agrees with them.
"""

import importlib.util
import subprocess
import sys
from fractions import Fraction as F
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "examples" / "abelian_vs_k3_walls.py"


def _load_example():
    spec = importlib.util.spec_from_file_location("abelian_vs_k3_walls", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_worked_example_radii():
    vals = _load_example().worked_example_values()
    # Same class v = (1, 0, -2), same center on both surfaces.
    assert vals["abelian_center"] == F(-5, 2)
    assert vals["k3_center"] == F(-5, 2)
    # Abelian (bare) 17/4 vs K3 (ch2 -> ch2 + r shim) 21/4.
    assert vals["abelian_radius_sq"] == F(17, 4)
    assert vals["k3_radius_sq"] == F(21, 4)
    # The shift is exactly 2/d = 2/2 = 1.
    assert vals["k3_radius_sq"] - vals["abelian_radius_sq"] == F(2, 2) == 1
    assert vals["shift"] == 1
    # Genuine integral-l K3 lattice class (2,0,-3)/(2,2,-1): Bayer-Macri type.
    assert vals["genuine_center"] == F(1)
    assert vals["genuine_radius_sq"] == F(1, 2)
    assert vals["genuine_wall_type"] == "divisorial"
    assert vals["genuine_subtype"] == "brill-noether"


def test_worked_example_exits_zero():
    result = subprocess.run(
        [sys.executable, str(_SCRIPT)],
        cwd=str(_REPO),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "17/4" in result.stdout
    assert "21/4" in result.stdout
