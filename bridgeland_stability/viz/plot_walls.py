"""Bridgeland wall plot -- nested concentric semicircular walls in (s, t).

Renders the numerical / actual walls of a Chern character ``v`` as semicircles
in the ``(s, t)`` upper half-plane (each centered on the ``s``-axis at ``s0``
with radius ``R = sqrt(R^2)``).  The figure typesets real mathematics, keeps the
circles geometrically round (equal aspect), orders the walls by radius with a
perceptually-uniform sequential colormap + labelled colorbar, derives ``t_max``
from the data so nothing is clipped, and -- for dense wall families that pinch
into an unreadable smear near the origin / near ``s = mu_v`` -- adds a zoomed
inset re-drawing the same walls.

Shares the house style in :mod:`bridgeland_stability.viz.style`; obeys exact
rational arithmetic (centers / radius^2 are :class:`fractions.Fraction`, exposed
as LaTeX fractions in plotly hovers and degraded to ``float`` only for the
geometry of the arcs).
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import List, Optional, Sequence, Tuple

from ..chern import ChernChar, Number, Q
from ..varieties import Surface
from ..walls import VerticalWall, Wall, compute_walls
from . import style


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
def _is_drawable(w) -> bool:
    """True for a real, finite-radius semicircular wall (skip vertical / nan)."""
    if isinstance(w, VerticalWall):
        return False
    if not getattr(w, "is_real", True):
        return False
    R = getattr(w, "radius", float("nan"))
    if R is None or math.isinf(R) or math.isnan(R) or R <= 0:
        return False
    return True


def _arc(w, n: int = 360) -> Tuple[List[float], List[float]]:
    """The upper semicircle ``(s, t)`` of wall ``w`` sampled with ``n+1`` points."""
    cx, R = float(w.center), w.radius
    th = [math.pi * i / n for i in range(n + 1)]
    return [cx + R * math.cos(a) for a in th], [R * math.sin(a) for a in th]


def _vertical_walls(walls: Sequence) -> List[VerticalWall]:
    return [w for w in walls if isinstance(w, VerticalWall) and w.s_value is not None]


def _inset_window(
    drawable: Sequence,
    mu_v: float,
    t_max: float,
) -> Optional[Tuple[float, float, float, float]]:
    """Data-driven near-origin / near ``s = mu_v`` zoom window ``(x0, x1, y0, y1)``.

    The dense smear is the cluster of *small* walls around ``s = mu_v``; we take
    the walls whose radius is at most a quarter of the largest radius and bound
    their leftmost / rightmost extent, then pad.  Returns ``None`` when there is
    no meaningful cluster (too few walls / everything the same size).
    """
    radii = [w.radius for w in drawable]
    if not radii:
        return None
    R_big = max(radii)
    cluster = [w for w in drawable if w.radius <= 0.30 * R_big]
    if len(cluster) < 3:
        # No dense small-wall smear; fall back to a generic near-mu_v box.
        cluster = sorted(drawable, key=lambda w: w.radius)[: max(3, len(drawable) // 3)]
    if not cluster:
        return None
    lefts = [float(w.center) - w.radius for w in cluster]
    rights = [float(w.center) + w.radius for w in cluster]
    apex = max(w.radius for w in cluster)
    x0, x1 = min(lefts), max(rights)
    # keep the window centered-ish on mu_v and never zero-width
    half = max((x1 - x0) / 2.0, apex, 0.05)
    cx = mu_v
    x0, x1 = cx - half, cx + half
    y0 = 0.0
    y1 = min(t_max, apex * 1.18 + 0.02)
    if y1 <= 0:
        y1 = t_max * 0.25
    return x0, x1, y0, y1


def _color_for_radius(R: float, R_lo: float, R_hi: float, cmap_lut):
    import matplotlib as mpl

    if R_hi <= R_lo:
        frac = 0.5
    else:
        frac = (R - R_lo) / (R_hi - R_lo)
    return mpl.colors.to_hex(cmap_lut(frac))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def plot_walls(
    v: ChernChar,
    surface: Surface,
    *,
    s_range: Tuple[Number, Number] = (-4, 4),
    t_max: Optional[float] = None,
    rank_bound: int = 5,
    walls: Optional[List] = None,
    max_walls: int = 30,
    backend: str = "plotly",
    theme: str = "light",
    inset: Optional[bool] = None,
    cmap: str = "cividis",
    show: bool = False,
    html_path: Optional[str] = None,
):
    """Plot the nested semicircular Bridgeland walls for ``v`` on ``surface``.

    Walls are ordered by radius (largest first) and colored along a
    perceptually-uniform sequential colormap (``cmap``, default ``"cividis"``)
    with a labelled "wall radius R" colorbar -- never a categorical rainbow.

    Parameters
    ----------
    v, surface
        The Chern character and its polarized surface.
    s_range
        ``(s_min, s_max)`` window on the ``s``-axis (used both to clip and, when
        ``walls`` is None, to bound :func:`~bridgeland_stability.walls.compute_walls`).
    t_max
        Optional override for the top of the ``t``-axis.  Default ``None`` makes
        it *data-driven*: ``1.1 x`` the largest wall apex so no semicircle is
        clipped.
    walls
        An explicit list of :class:`~bridgeland_stability.walls.Wall` /
        :class:`~bridgeland_stability.walls.ActualWall` /
        :class:`~bridgeland_stability.walls.VerticalWall` (e.g. from
        :func:`~bridgeland_stability.walls.actual_walls`).  When None, numerical
        walls are computed.
    inset
        Whether to draw the near-origin zoom inset (matplotlib only).  Default
        ``None`` -> on automatically when there are several drawable walls.
    cmap
        Sequential colormap name (``"cividis"`` or ``"viridis"``).

    Returns
    -------
    The matplotlib :class:`~matplotlib.figure.Figure` or the plotly
    :class:`~plotly.graph_objects.Figure`.
    """
    if walls is None:
        walls = compute_walls(v, surface, s_range=s_range, rank_bound=rank_bound)
    if max_walls is not None:
        walls = list(walls[:max_walls])

    s_min, s_max = float(s_range[0]), float(s_range[1])

    drawable = [w for w in walls if _is_drawable(w)]
    # largest radius first (sort is defensive; compute_walls/actual_walls already do)
    drawable.sort(key=lambda w: w.radius, reverse=True)
    verticals = _vertical_walls(walls)

    radii = [w.radius for w in drawable]
    R_hi = max(radii) if radii else 1.0
    R_lo = min(radii) if radii else 0.0
    if t_max is None:
        t_max = (R_hi * 1.1) if radii else 4.0

    try:
        mu_v = float(v.slope(surface.d)) if v.r != 0 else 0.0
    except Exception:
        mu_v = 0.0

    title = (
        rf"Bridgeland walls of ${v.to_latex()}$"
        rf" on {style.latex_variety_name(surface.name)}"
    )

    if backend == "plotly":
        return _plot_plotly(
            v, surface, drawable, verticals, s_min, s_max, t_max, mu_v,
            R_lo, R_hi, cmap, title, theme, show, html_path,
        )

    return _plot_mpl(
        v, surface, drawable, verticals, s_min, s_max, t_max, mu_v,
        R_lo, R_hi, cmap, title, theme, inset,
    )


# --------------------------------------------------------------------------
# matplotlib backend
# --------------------------------------------------------------------------
def _draw_arcs(ax, drawable, R_lo, R_hi, cmap_lut, lw):
    for w in drawable:
        xs, ys = _arc(w)
        c = _color_for_radius(w.radius, R_lo, R_hi, cmap_lut)
        ax.plot(xs, ys, lw=lw, color=c, solid_capstyle="round", antialiased=True)


def _plot_mpl(
    v, surface, drawable, verticals, s_min, s_max, t_max, mu_v,
    R_lo, R_hi, cmap, title, theme, inset,
):
    plt = style.require_mpl()
    from matplotlib import colormaps

    cmap_lut = colormaps[cmap]
    ink = style.theme_ink(theme)

    if inset is None:
        inset = len(drawable) >= 5

    with style.mpl_context(theme):
        fig, ax = plt.subplots(figsize=(8.6, 6.2))

        _draw_arcs(ax, drawable, R_lo, R_hi, cmap_lut, lw=1.5)

        # vertical (degenerate) walls drawn as dashed guide lines
        for w in verticals:
            sv = float(w.s_value)
            if s_min <= sv <= s_max:
                ax.plot(
                    [sv, sv], [0, t_max], ls="--", lw=1.1,
                    color=ink["muted"], alpha=0.9, zorder=2,
                )

        ax.set_xlim(s_min, s_max)
        ax.set_ylim(0, t_max)
        ax.set_aspect("equal", adjustable="datalim")

        ax.set_xlabel(style.LABEL["s"])
        ax.set_ylabel(style.LABEL["t"])
        ax.set_title(title, pad=12)

        # major + minor grid, ink, spines
        style.style_axes(ax, theme, minor=True)

        if drawable:
            style.add_colorbar(
                fig, ax,
                values=[R_lo, R_hi],
                cmap=cmap,
                label=r"wall radius $R$",
                theme=theme,
            )

        # near-origin / near s=mu_v zoom inset
        if inset and len(drawable) >= 3:
            win = _inset_window(drawable, mu_v, t_max)
            if win is not None:
                x0, x1, y0, y1 = win
                axins = ax.inset_axes([0.62, 0.55, 0.36, 0.40])
                _draw_arcs(axins, drawable, R_lo, R_hi, cmap_lut, lw=1.1)
                for w in verticals:
                    sv = float(w.s_value)
                    if x0 <= sv <= x1:
                        axins.plot(
                            [sv, sv], [y0, y1], ls="--", lw=0.9,
                            color=ink["muted"], alpha=0.9,
                        )
                axins.set_xlim(x0, x1)
                axins.set_ylim(y0, y1)
                axins.set_aspect("equal", adjustable="box")
                axins.tick_params(labelsize=7)
                style.style_axes(axins, theme, minor=False)
                ax.indicate_inset_zoom(axins, edgecolor=ink["edge"], alpha=0.85)

        # re-assert defensively (title color, facecolors) after inset work
        style.style_axes(ax, theme, minor=True)

    return fig


# --------------------------------------------------------------------------
# plotly backend
# --------------------------------------------------------------------------
def _plotly_hover(w) -> str:
    """An exact-fraction hovertemplate: center, R^2, and the destabilizer class."""
    sub = getattr(w, "subobject", None)
    parts = [
        f"center $s_0={style.latex_frac(w.center)}$",
        f"$R^2={style.latex_frac(w.radius_sq)}$",
    ]
    if sub is not None:
        parts.append(f"destabilizer ${sub.to_latex(name='w')}$")
    return "<br>".join(parts) + "<extra></extra>"


def _plot_plotly(
    v, surface, drawable, verticals, s_min, s_max, t_max, mu_v,
    R_lo, R_hi, cmap, title, theme, show, html_path,
):
    go = style.require_plotly()
    from matplotlib import colormaps

    cmap_lut = colormaps[cmap]

    fig = go.Figure()
    style.apply_plotly_theme(fig, theme)

    # exact-fraction semicircles (go.Scatter ONLY -- float64 SVG, deep-zoom-safe)
    for w in drawable:
        xs, ys = _arc(w)
        color = _color_for_radius(w.radius, R_lo, R_hi, cmap_lut)
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(color=color, width=1.6),
                name=f"R={w.radius:.3f}",
                showlegend=False,
                hovertemplate=_plotly_hover(w),
            )
        )

    # vertical degenerate walls as dashed lines
    ink = style.theme_ink(theme)
    for w in verticals:
        sv = float(w.s_value)
        if s_min <= sv <= s_max:
            fig.add_trace(
                go.Scatter(
                    x=[sv, sv],
                    y=[0, t_max],
                    mode="lines",
                    line=dict(color=ink["muted"], width=1.1, dash="dash"),
                    name=f"vertical s={style.latex_frac(w.s_value)}",
                    showlegend=False,
                    hovertemplate=(
                        f"vertical wall $s={style.latex_frac(w.s_value)}$<extra></extra>"
                    ),
                )
            )

    # a hidden marker trace carrying the colorbar for "wall radius R"
    if drawable:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                showlegend=False,
                hoverinfo="skip",
                marker=dict(
                    colorscale=cmap,
                    cmin=R_lo,
                    cmax=R_hi,
                    color=[R_lo],
                    showscale=True,
                    colorbar=dict(
                        title=dict(text="wall radius $R$", side="right"),
                        outlinecolor=ink["edge"],
                        tickcolor=ink["ink"],
                    ),
                ),
            )
        )

    fig.update_layout(
        title=title,
        xaxis_title=style.LABEL["s"],
        yaxis_title=style.LABEL["t"],
    )
    # ranges + equal aspect locked so the semicircles stay round under zoom
    fig.update_xaxes(range=[s_min, s_max])
    fig.update_yaxes(range=[0, t_max], scaleanchor="x", scaleratio=1)

    style.emit_plotly(fig, show, html_path)
    return fig
