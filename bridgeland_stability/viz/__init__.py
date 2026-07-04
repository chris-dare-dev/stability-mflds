"""Visualization layer (matplotlib static + plotly interactive).

Backward-compatible facade: ``from bridgeland_stability import viz`` then
``viz.plot_dlp_curve(...)`` / ``viz.plot_walls(...)`` / ``viz.plot_threefold_bg(...)``
work exactly as before.  Each plot family lives in its own submodule and shares
the house style defined in :mod:`bridgeland_stability.viz.style`.

Both backends are imported lazily (inside the plot functions, via
:func:`style.require_mpl` / :func:`style.require_plotly`) so importing this
package never requires matplotlib or plotly.

Functions return the underlying figure object; pass ``show=True`` to display,
``html_path=...`` (plotly) to write standalone interactive HTML, or
``theme="light"|"dark"`` to pick the palette.  Use :func:`save_publication` for
vector PDF+SVG export.
"""

from __future__ import annotations

from . import style
from .plot_dlp import plot_dlp_curve
from .plot_threefold import plot_threefold_bg
from .plot_walls import plot_walls
from .style import (
    INK,
    LABEL,
    OKABE_ITO,
    PALETTE,
    SYM,
    save_publication,
)

__all__ = [
    "plot_dlp_curve",
    "plot_walls",
    "plot_threefold_bg",
    "save_publication",
    "style",
    "OKABE_ITO",
    "PALETTE",
    "INK",
    "SYM",
    "LABEL",
]
