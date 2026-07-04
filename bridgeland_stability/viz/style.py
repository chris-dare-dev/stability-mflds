"""Shared style/theme/typography layer for the visualization package.

This module is the single source of truth that every plot family
(:mod:`.plot_dlp`, :mod:`.plot_walls`, :mod:`.plot_threefold`) binds to, so the
figures share one coherent, publication-ready, colorblind-safe identity.

Nothing here pulls matplotlib/plotly at import time -- the heavy libraries are
imported lazily inside :func:`require_mpl` / :func:`require_plotly` so that
``import bridgeland_stability`` (and ``import bridgeland_stability.viz``) stay
dependency-free.

What it provides
----------------
* **Palette** -- :data:`OKABE_ITO` (colorblind-safe categorical cycle),
  :data:`PALETTE` (semantic roles reused by *both* artists and legend handles),
  :data:`INK` (per-theme foreground/grid colors), and discrete/sequential color
  helpers (never a smooth ramp for a discrete quantity).
* **Typography** -- :data:`SYM`/:data:`LABEL` LaTeX constants, :func:`m` (wrap
  as mathtext), and the exact-fraction tick tools :func:`fraction_formatter` /
  :func:`set_fraction_ticks` (re-exported :func:`latex_frac` &c. from
  :mod:`bridgeland_stability._latex`).
* **Theming** -- :func:`mpl_context` (light/dark ``.mplstyle`` context),
  :func:`style_axes` (defensive grid/spine/ink application),
  :func:`apply_plotly_theme`, and :func:`emit_plotly` (HTML w/ MathJax).
* **Export** -- :func:`save_publication` (vector PDF+SVG, font-embedded).
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# Pure-stdlib LaTeX helpers (re-exported for convenience).
from .._latex import latex_chern, latex_frac, latex_frac_from_float, latex_sqrt

__all__ = [
    "require_mpl",
    "require_plotly",
    "OKABE_ITO",
    "PALETTE",
    "INK",
    "SYM",
    "LABEL",
    "m",
    "latex_variety_name",
    "latex_frac",
    "latex_frac_from_float",
    "latex_sqrt",
    "latex_chern",
    "fraction_formatter",
    "set_fraction_ticks",
    "discrete_colors",
    "rank_color_map",
    "sequential_colors",
    "proxy_handle",
    "add_colorbar",
    "add_badge",
    "style_dir",
    "style_path",
    "mpl_context",
    "style_axes",
    "theme_ink",
    "apply_plotly_theme",
    "plotly_font",
    "emit_plotly",
    "save_publication",
    "THEMES",
]

THEMES = ("light", "dark")


# --------------------------------------------------------------------------
# Lazy backend requires
# --------------------------------------------------------------------------
def require_mpl():
    """Return ``matplotlib.pyplot`` or raise a helpful ImportError."""
    try:
        import matplotlib.pyplot as plt  # noqa: F401

        return plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "matplotlib is required for static plots: "
            "pip install 'bridgeland-stability[viz]'"
        ) from exc


def require_plotly():
    """Return ``plotly.graph_objects`` or raise a helpful ImportError."""
    try:
        import plotly.graph_objects as go  # noqa: F401

        return go
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "plotly is required for interactive plots: "
            "pip install 'bridgeland-stability[viz]'"
        ) from exc


# --------------------------------------------------------------------------
# Palette  (Okabe-Ito categorical + semantic roles + per-theme ink)
# --------------------------------------------------------------------------
#: The Okabe-Ito colorblind-safe qualitative palette (8 colors).
OKABE_ITO: List[str] = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]

#: Semantic color roles, theme-independent and chosen to read on light *and*
#: dark grounds.  Artists and their legend handles MUST both pull from here so
#: the legend never lies about the plotted color.
PALETTE: Dict[str, str] = {
    "curve": "#0072B2",        # DLP delta(mu) envelope
    "fill": "#0072B2",         # nonempty-region fill (used with low alpha)
    "cusp": "#009E73",         # cusp tips lying ON the curve
    "exceptional": "#CC79A7",  # exceptional bundles, below the curve
    "boundary": "#D55E00",     # threefold alpha_crit(beta)
    "asymptote": "#E69F00",    # vertical asymptote at the degenerate beta
    "forbidden": "#D55E00",    # forbidden/degenerate region shade (low alpha)
    "floor": "#888888",        # delta = 1/2 floor guide line
    "proven": "#009E73",       # PROVEN status badge (green)
    "conjectural": "#E69F00",  # CONJECTURAL status badge (amber)
}

#: Per-theme foreground/background/grid ink.
INK: Dict[str, Dict[str, object]] = {
    "light": {
        "fig_bg": "white",
        "axes_bg": "white",
        "ink": "#222222",
        "edge": "#444444",
        "grid": "#b8b8b8",
        "grid_alpha": 0.35,
        "minor_grid_alpha": 0.18,
        "muted": "#666666",
    },
    "dark": {
        "fig_bg": "#11151c",
        "axes_bg": "#11151c",
        "ink": "#e6e6e6",
        "edge": "#8aa0bf",
        "grid": "#3a4658",
        "grid_alpha": 0.5,
        "minor_grid_alpha": 0.22,
        "muted": "#9fb0c8",
    },
}


def theme_ink(theme: str = "light") -> Dict[str, object]:
    """Return the :data:`INK` dict for ``theme`` (validates the name)."""
    if theme not in INK:
        raise ValueError(f"theme must be one of {THEMES!r}, got {theme!r}")
    return INK[theme]


# --------------------------------------------------------------------------
# LaTeX symbol constants  (mathtext: use \frac, \mathrm, \mathbb -- never \tfrac)
# --------------------------------------------------------------------------
SYM: Dict[str, str] = {
    "mu": r"\mu",
    "Delta": r"\Delta",
    "delta": r"\delta",
    "alpha": r"\alpha",
    "beta": r"\beta",
    "chi": r"\chi",
    "s": r"s",
    "t": r"t",
    "P2": r"\mathbb{P}^2",
    "P3": r"\mathbb{P}^3",
    "Rmax": r"R_{\max}",
    "alpha_crit": r"\alpha_{\mathrm{crit}}(\beta)",
    "ch1b": r"\mathrm{ch}_1^{\beta}",
}


def m(expr: str) -> str:
    """Wrap a LaTeX body in ``$...$`` for matplotlib/plotly labels."""
    return f"${expr}$"


def latex_variety_name(name: str) -> str:
    """Render a catalog variety ``name`` with LaTeX math spans for use in titles.

    Turns the ASCII catalog names into mixed text+math, e.g.
    ``"P^3" -> "$\\mathbb{P}^3$"``, ``"P^1 x P^1" -> "$\\mathbb{P}^1$ x $\\mathbb{P}^1$"``,
    ``"quadric 3-fold Q in P^4" -> "quadric 3-fold Q in $\\mathbb{P}^4$"``,
    ``"K3 (H^2=2)" -> "K3 ($H^2=2$)"``.  Plain prose passes through unchanged.
    """
    import re

    s = name
    # P^n  ->  $\mathbb{P}^n$
    s = re.sub(r"\bP\^(\d+)", lambda mo: m(rf"\mathbb{{P}}^{{{mo.group(1)}}}"), s)
    # H^2=2  ->  $H^2=2$
    s = re.sub(r"\bH\^(\d+)\s*=\s*([0-9]+)", lambda mo: m(rf"H^{{{mo.group(1)}}}={mo.group(2)}"), s)
    return s


#: Canonical, rendered axis labels (English glosses are fine; the symbols are
#: real math, never ASCII ``mu``/``Delta``/``alpha``).
LABEL: Dict[str, str] = {
    "mu": rf"slope $\mu$",
    "Delta": rf"discriminant $\Delta$",
    "s": r"$s$",
    "t": r"$t$",
    "beta": r"$\beta$",
    "alpha": r"$\alpha$",
    "P2": m(SYM["P2"]),
    "P3": m(SYM["P3"]),
}


# --------------------------------------------------------------------------
# Exact-fraction tick tooling
# --------------------------------------------------------------------------
def fraction_formatter(max_denominator: int = 16):
    """A ``FuncFormatter`` that renders each tick as its simplest fraction.

    ``0.4 -> $\\frac{2}{5}$``, ``-0.5 -> $-\\frac{1}{2}$``, ``2.0 -> $2$``.
    Use on axes where the salient locations are simple rationals.
    """
    from matplotlib.ticker import FuncFormatter

    def _fmt(value: float, _pos=None) -> str:
        if abs(value) < 1e-12:
            return "$0$"
        return m(latex_frac_from_float(value, max_denominator))

    return FuncFormatter(_fmt)


def set_fraction_ticks(ax, axis: str, values: Sequence, *, minor: bool = False) -> None:
    """Place explicit ticks at the *exact* rationals ``values`` with fraction labels.

    ``values`` may be :class:`~fractions.Fraction` or numbers; labels are the
    exact fractions (no float round-trip).  ``axis`` is ``'x'`` or ``'y'``.
    """
    from fractions import Fraction

    locs = [float(v) for v in values]
    labels = [m(latex_frac(v if isinstance(v, Fraction) else Fraction(v))) for v in values]
    target = ax.xaxis if axis == "x" else ax.yaxis
    target.set_ticks(locs, minor=minor)
    if not minor:
        target.set_ticklabels(labels)


# --------------------------------------------------------------------------
# Discrete / sequential color helpers
# --------------------------------------------------------------------------
def discrete_colors(n: int, *, prefer_qualitative: bool = True) -> List[str]:
    """``n`` distinct hex colors: Okabe-Ito for ``n <= 8`` else sampled cividis.

    For genuinely *categorical* data (e.g. integer sheaf ranks) -- never a
    smooth ramp.
    """
    if n <= 0:
        return []
    if prefer_qualitative and n <= len(OKABE_ITO):
        # Skip pure black (index 7) unless we truly need all 8.
        order = [0, 1, 2, 3, 4, 5, 6, 7]
        return [OKABE_ITO[i] for i in order[:n]]
    import matplotlib as mpl
    from matplotlib import colormaps

    cmap = colormaps["cividis"]
    return [mpl.colors.to_hex(cmap(i / max(1, n - 1))) for i in range(n)]


def rank_color_map(ranks: Iterable[int]) -> Dict[int, str]:
    """Map each *distinct* integer rank to a categorical color (sorted)."""
    uniq = sorted(set(int(r) for r in ranks))
    colors = discrete_colors(len(uniq))
    return {r: c for r, c in zip(uniq, colors)}


def sequential_colors(n: int, cmap: str = "viridis") -> List[str]:
    """``n`` hex colors sampled along a perceptually-uniform sequential ``cmap``.

    For an *ordered* family (e.g. walls ordered by radius), paired with
    :func:`add_colorbar`.
    """
    if n <= 0:
        return []
    import matplotlib as mpl
    from matplotlib import colormaps

    cm = colormaps[cmap]
    if n == 1:
        return [mpl.colors.to_hex(cm(0.5))]
    return [mpl.colors.to_hex(cm(i / (n - 1))) for i in range(n)]


# --------------------------------------------------------------------------
# Legend / colorbar / badge artists
# --------------------------------------------------------------------------
def proxy_handle(
    *,
    color: str,
    label: str,
    marker: Optional[str] = None,
    linestyle: str = "None",
    markeredgecolor: Optional[str] = None,
    linewidth: float = 2.0,
    markersize: float = 8.0,
):
    """A ``Line2D`` legend proxy whose color is *identical* to the artist's.

    Pass the SAME color object you used for the plotted artist so the legend
    swatch provably equals the plotted color (see ``tests/test_viz.py``).
    """
    from matplotlib.lines import Line2D

    return Line2D(
        [0],
        [0],
        color=color,
        marker=marker,
        linestyle=linestyle,
        markerfacecolor=color,
        markeredgecolor=markeredgecolor if markeredgecolor is not None else color,
        linewidth=linewidth,
        markersize=markersize,
        label=label,
    )


def add_colorbar(
    fig,
    ax,
    *,
    values: Sequence[float],
    cmap: str = "viridis",
    label: str = "",
    theme: str = "light",
):
    """Attach a labeled sequential colorbar spanning ``[min(values), max(values)]``.

    Returns the ``Colorbar``.  Use for *ordered continuous* scalars only
    (e.g. wall radius); for discrete categories use :func:`rank_color_map` +
    :func:`proxy_handle`.
    """
    import matplotlib as mpl
    from matplotlib import colormaps

    ink = theme_ink(theme)
    lo, hi = (min(values), max(values)) if values else (0.0, 1.0)
    norm = mpl.colors.Normalize(vmin=lo, vmax=hi)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=colormaps[cmap])
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(label, color=ink["ink"])
    cbar.ax.yaxis.set_tick_params(color=ink["ink"], labelcolor=ink["ink"])
    cbar.outline.set_edgecolor(ink["edge"])
    return cbar


def add_badge(
    ax,
    text: str,
    *,
    color: str,
    theme: str = "light",
    xy: Tuple[float, float] = (0.985, 0.97),
    ha: str = "right",
    va: str = "top",
    fontsize: float = 11.5,
):
    """A rounded status badge in axes-fraction coords (e.g. PROVEN / CONJECTURAL).

    The badge text is white-on-``color``; returns the ``Text`` artist.
    """
    return ax.text(
        xy[0],
        xy[1],
        text,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=fontsize,
        fontweight="bold",
        color="white",
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": color,
            "edgecolor": "none",
            "alpha": 0.95,
        },
        zorder=10,
    )


# --------------------------------------------------------------------------
# Matplotlib theming
# --------------------------------------------------------------------------
def style_dir() -> Path:
    """Directory holding the bundled ``.mplstyle`` sheets."""
    return Path(__file__).parent / "styles"


def style_path(theme: str = "light") -> Path:
    """Path to the ``bridgeland-{theme}.mplstyle`` sheet."""
    if theme not in INK:
        raise ValueError(f"theme must be one of {THEMES!r}, got {theme!r}")
    return style_dir() / f"bridgeland-{theme}.mplstyle"


@contextmanager
def mpl_context(theme: str = "light"):
    """Context manager applying the bundled house style for ``theme``.

    Resets to matplotlib defaults first, then layers the ``.mplstyle`` sheet, so
    styling is deterministic regardless of the caller's rc state.  Create the
    figure *inside* this block (fonts/sizes/cycle/constrained-layout are read at
    artist-creation time); :func:`style_axes` re-asserts colors defensively.
    """
    plt = require_mpl()
    with plt.style.context(["default", str(style_path(theme))]):
        yield


def style_axes(ax, theme: str = "light", *, minor: bool = False) -> None:
    """Defensively apply theme ink, grid, spines and facecolors to ``ax``.

    Safe to call after labels/title are set; it re-asserts the foreground color
    on them so the figure is theme-correct even outside :func:`mpl_context`.
    """
    ink = theme_ink(theme)
    fig = ax.figure
    fig.set_facecolor(ink["fig_bg"])
    ax.set_facecolor(ink["axes_bg"])

    # grid
    ax.grid(True, which="major", color=ink["grid"], alpha=ink["grid_alpha"], linewidth=0.6)
    if minor:
        ax.minorticks_on()
        ax.grid(
            True,
            which="minor",
            color=ink["grid"],
            alpha=ink["minor_grid_alpha"],
            linewidth=0.4,
        )
    ax.set_axisbelow(True)

    # spines
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(ink["edge"])

    # ink on ticks + labels + title
    ax.tick_params(colors=ink["ink"], which="both")
    ax.xaxis.label.set_color(ink["ink"])
    ax.yaxis.label.set_color(ink["ink"])
    if ax.get_title():
        ax.title.set_color(ink["ink"])


# --------------------------------------------------------------------------
# Plotly theming
# --------------------------------------------------------------------------
def plotly_font(theme: str = "light") -> dict:
    """Serif font dict (Computer-Modern-like) honoring the theme ink color."""
    ink = theme_ink(theme)
    return {
        "family": "CMU Serif, Latin Modern Roman, Georgia, serif",
        "color": ink["ink"],
        "size": 14,
    }


def apply_plotly_theme(fig, theme: str = "light"):
    """Apply transparent background, serif font, Okabe-Ito colorway and grid.

    Call this first; set figure-specific axis ranges / ``scaleanchor`` after.
    Returns ``fig`` for chaining.
    """
    ink = theme_ink(theme)
    fig.update_layout(
        template="plotly_white" if theme == "light" else "plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=plotly_font(theme),
        colorway=OKABE_ITO,
        title=dict(font=dict(size=18, color=ink["ink"])),
        legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(color=ink["ink"])),
        margin=dict(l=72, r=36, t=72, b=60),
    )
    grid_kw = dict(
        showgrid=True,
        gridcolor=ink["grid"],
        gridwidth=0.6,
        zeroline=False,
        linecolor=ink["edge"],
        color=ink["ink"],
    )
    fig.update_xaxes(**grid_kw)
    fig.update_yaxes(**grid_kw)
    return fig


def emit_plotly(fig, show: bool, html_path: Optional[str]) -> None:
    """Write standalone HTML (plotly + MathJax via CDN) and/or show the figure."""
    if html_path:
        fig.write_html(html_path, include_plotlyjs="cdn", include_mathjax="cdn")
    if show:  # pragma: no cover
        fig.show()


# --------------------------------------------------------------------------
# Publication export (vector, font-embedded)
# --------------------------------------------------------------------------
def save_publication(
    fig,
    path: str,
    *,
    formats: Sequence[str] = ("pdf", "svg"),
    dpi: int = 300,
    transparent: bool = True,
) -> List[Path]:
    """Export ``fig`` to vector formats with embedded (TrueType/none) fonts.

    For a matplotlib figure: writes each of ``formats`` (default PDF+SVG) with
    ``pdf.fonttype=42``, ``ps.fonttype=42``, ``svg.fonttype='none'`` (editable,
    journal-safe -- no Type-3 fonts), ``bbox_inches='tight'``.  For a plotly
    figure: uses Kaleido (raises a clear error if Kaleido is absent).

    ``path`` may carry an extension or not; one file per format is written using
    the stem.  Returns the list of written paths.
    """
    base = Path(path)
    stem = base.with_suffix("")
    written: List[Path] = []
    fmts = [f.lstrip(".").lower() for f in formats]

    if hasattr(fig, "savefig"):  # matplotlib Figure
        plt = require_mpl()
        rc = {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}
        with plt.rc_context(rc):
            for fmt in fmts:
                out = stem.with_suffix(f".{fmt}")
                fig.savefig(out, dpi=dpi, bbox_inches="tight", transparent=transparent)
                written.append(out)
        return written

    if hasattr(fig, "write_image"):  # plotly Figure
        for fmt in fmts:
            out = stem.with_suffix(f".{fmt}")
            try:
                fig.write_image(str(out))
            except Exception as exc:  # pragma: no cover - needs kaleido
                raise RuntimeError(
                    "plotly static export requires kaleido: pip install kaleido"
                ) from exc
            written.append(out)
        return written

    raise TypeError("save_publication expects a matplotlib or plotly figure")
