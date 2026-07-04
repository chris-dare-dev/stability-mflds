"""Drézet–Le Potier curve plot for P^2.

Renders the fractal upper envelope ``delta(mu)`` with the nonempty moduli
region (``Delta >= delta``) shaded, the exceptional-bundle cusp tips lying ON
the curve, and the exceptional bundles themselves sitting STRICTLY BELOW it.

Design decisions tied to the "exact data" value proposition
-----------------------------------------------------------
* **Category vs. rank are decoupled.** Marker *shape* encodes the geometric
  category (``^`` triangle = cusp tip on the curve; ``x`` = exceptional bundle
  below it); marker *color* encodes the integer Markov *rank* via
  :func:`style.rank_color_map`, a DISCRETE categorical map -- never a smooth
  colorbar (ranks are a finite Markov set {1,2,5,13,29,...}, not a continuum).
* **The legend cannot lie.** Every legend swatch is built with
  :func:`style.proxy_handle` using the IDENTICAL hex that colored the artist.
* **Fractal self-similarity is made visible.** The main panel samples at a high
  ``R_max`` / density so many cusps resolve, and a zoomed inset re-plots a
  *finer* :func:`~bridgeland_stability.dlp.dlp_curve` over a narrow window
  around a chosen cusp at an even higher ``R_max`` to reveal the self-similar
  cusp cascade.
* **Ticks are exact rationals.** Salient Markov-fraction cusp slopes are shown
  as exact fractions via :func:`style.set_fraction_ticks`.

Shares the house style in :mod:`bridgeland_stability.viz.style`.
"""

from __future__ import annotations

from fractions import Fraction
from typing import List, Optional, Tuple

from ..chern import Number, Q
from ..dlp import DLPCurve, dlp_curve
from . import style

__all__ = ["plot_dlp_curve"]

# Salient Markov-fraction cusp slopes worth annotating / ticking inside [0, 1].
_SALIENT = (Fraction(0), Fraction(2, 5), Fraction(1, 2), Fraction(3, 5), Fraction(1))


def _in_window(value: float, lo: float, hi: float) -> bool:
    return lo <= value <= hi


def _salient_ticks(lo: float, hi: float) -> List[Fraction]:
    """Salient cusp slopes (plus window endpoints) that fall inside ``[lo, hi]``."""
    vals = {Q(lo), Q(hi)}
    for fr in _SALIENT:
        if lo - 1e-12 <= float(fr) <= hi + 1e-12:
            vals.add(fr)
    return sorted(vals)


def _cusp_rank(alpha: Fraction) -> int:
    """Markov rank of a cusp tip at slope ``alpha`` (= denominator in lowest terms)."""
    return alpha.denominator


def plot_dlp_curve(
    curve: Optional[DLPCurve] = None,
    *,
    mu_min: Number = 0,
    mu_max: Number = 1,
    R_max: int = 200,
    samples_per_unit: int = 1200,
    backend: str = "plotly",
    theme: str = "light",
    show: bool = False,
    html_path: Optional[str] = None,
    inset_center: Number = Fraction(2, 5),
    inset_halfwidth: Number = Fraction(1, 20),
    inset_R_max: int = 600,
    inset_samples_per_unit: int = 4000,
    show_inset: bool = True,
):
    """Plot the Drézet–Le Potier curve ``delta(mu)`` on ``P^2``.

    Parameters
    ----------
    curve : DLPCurve, optional
        A pre-built curve.  If ``None``, one is built from ``mu_min``/``mu_max``
        at the (now higher) default ``R_max`` / ``samples_per_unit`` so the
        fractal cusp cascade resolves.
    mu_min, mu_max : Number
        Slope window for the main panel (used only when ``curve is None``).
    R_max, samples_per_unit : int
        Resolution controls for the main panel.  Higher ``R_max`` resolves more
        (higher-rank, finer) cusps; ``samples_per_unit`` controls smoothness.
    backend : {"plotly", "matplotlib"}
    theme : {"light", "dark"}
    inset_center, inset_halfwidth : Number
        A cusp slope to zoom on and the half-width of the zoom window.  The
        inset re-plots a finer curve here to reveal self-similarity (matplotlib
        backend only).
    inset_R_max, inset_samples_per_unit : int
        Resolution for the zoomed inset (higher than the main panel).
    show_inset : bool
        Whether to draw the zoomed inset (matplotlib backend).

    Returns
    -------
    A matplotlib ``Figure`` or a plotly ``Figure``.
    """
    if curve is None:
        curve = dlp_curve(mu_min, mu_max, R_max, samples_per_unit)
    xs, ys = curve.as_floats()
    x_lo, x_hi = float(curve.mus[0]), float(curve.mus[-1])
    top = max(1.05, max(ys) * 1.05)

    # Exceptional bundles whose slope falls inside the window (for hover text).
    in_window = [
        b for b in curve.bundles if _in_window(float(b.slope), x_lo, x_hi)
    ]

    # Rank -> categorical color, spanning BOTH cusps and exceptional points so a
    # rank shares one color across the two shapes.
    all_ranks = sorted(
        {_cusp_rank(a) for a, _ in curve.cusps}
        | {_cusp_rank(a) for a, _ in curve.exceptional_points}
    )
    rcolors = style.rank_color_map(all_ranks)

    title = (
        rf"Drézet–Le Potier curve ${style.SYM['delta']}({style.SYM['mu']})$ on "
        rf"${style.SYM['P2']}$ (${style.SYM['Rmax']}={curve.R_max}$)"
    )

    if backend == "plotly":
        return _plot_plotly(
            curve, xs, ys, top, in_window, rcolors, title, theme, show, html_path
        )

    return _plot_mpl(
        curve,
        xs,
        ys,
        top,
        x_lo,
        x_hi,
        rcolors,
        all_ranks,
        title,
        theme,
        inset_center,
        inset_halfwidth,
        inset_R_max,
        inset_samples_per_unit,
        show_inset,
    )


# --------------------------------------------------------------------------
# Plotly backend
# --------------------------------------------------------------------------
def _plot_plotly(curve, xs, ys, top, in_window, rcolors, title, theme, show, html_path):
    go = style.require_plotly()
    fig = go.Figure()
    style.apply_plotly_theme(fig, theme)

    floor_body = rf"{style.SYM['delta']} \geq \frac{{1}}{{2}}"

    # nonempty region fill (Delta >= delta) -- drawn first, no hover.
    fig.add_trace(
        go.Scatter(
            x=xs + xs[::-1],
            y=ys + [top] * len(ys),
            fill="toself",
            fillcolor="rgba(0,114,178,0.12)",
            line=dict(width=0),
            name=rf"moduli nonempty (${style.SYM['Delta']} \geq {style.SYM['delta']}$)",
            hoverinfo="skip",
        )
    )

    # the delta(mu) curve itself (exact fractions in hover).
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="lines",
            name=rf"${style.SYM['delta']}({style.SYM['mu']})$ curve",
            line=dict(color=style.PALETTE["curve"], width=2),
            customdata=[
                [
                    style.m(style.latex_frac(mu)),
                    style.m(style.latex_frac(d)),
                ]
                for mu, d in zip(curve.mus, curve.deltas)
            ],
            hovertemplate=(
                rf"${style.SYM['mu']}$ = %{{customdata[0]}}<br>"
                rf"${style.SYM['delta']}$ = %{{customdata[1]}}<extra></extra>"
            ),
        )
    )

    # cusp tips ON the curve, colored by rank (one trace per rank => discrete
    # categorical legend, never a colorbar).
    _add_rank_traces_plotly(
        fig,
        go,
        points=curve.cusps,
        rcolors=rcolors,
        symbol="triangle-up",
        size=9,
        category=r"cusp (on curve)",
    )

    # exceptional bundles BELOW the curve, colored by rank.
    _add_rank_traces_plotly(
        fig,
        go,
        points=curve.exceptional_points,
        rcolors=rcolors,
        symbol="x",
        size=9,
        category=r"exceptional (below)",
    )

    # floor guide delta = 1/2.
    fig.add_hline(
        y=0.5,
        line=dict(color=style.PALETTE["floor"], dash="dot", width=1.4),
        annotation_text=rf"floor ${floor_body}$",
        annotation_position="bottom right",
        annotation_font=dict(color=style.theme_ink(theme)["muted"]),
    )

    fig.update_layout(
        title=title,
        xaxis_title=style.LABEL["mu"],
        yaxis_title=style.LABEL["Delta"],
        legend=dict(itemsizing="constant"),
    )
    # Exact-fraction tick labels on x at salient slopes.
    ticks = _salient_ticks(float(curve.mus[0]), float(curve.mus[-1]))
    fig.update_xaxes(
        tickmode="array",
        tickvals=[float(t) for t in ticks],
        ticktext=[style.m(style.latex_frac(t)) for t in ticks],
    )
    fig.update_yaxes(range=[0.42, top])

    style.emit_plotly(fig, show, html_path)
    return fig


def _add_rank_traces_plotly(fig, go, *, points, rcolors, symbol, size, category):
    """One go.Scatter per rank so the discrete rank legend has constant swatches."""
    by_rank: dict = {}
    for alpha, val in points:
        by_rank.setdefault(_cusp_rank(alpha), []).append((alpha, val))
    for rank in sorted(by_rank):
        pts = by_rank[rank]
        color = rcolors[rank]
        fig.add_trace(
            go.Scatter(
                x=[float(a) for a, _ in pts],
                y=[float(v) for _, v in pts],
                mode="markers",
                name=rf"{category}, rank ${rank}$",
                legendgroup=category,
                legendgrouptitle_text=category,
                marker=dict(color=color, size=size, symbol=symbol,
                            line=dict(width=0.6, color=color)),
                customdata=[
                    [
                        style.m(style.latex_frac(a)),
                        style.m(style.latex_frac(v)),
                        rank,
                    ]
                    for a, v in pts
                ],
                hovertemplate=(
                    rf"{category}<br>"
                    rf"${style.SYM['mu']}$ = %{{customdata[0]}}<br>"
                    rf"${style.SYM['delta']}$ = %{{customdata[1]}}<br>"
                    r"rank = %{customdata[2]}<extra></extra>"
                ),
            )
        )


# --------------------------------------------------------------------------
# Matplotlib backend
# --------------------------------------------------------------------------
def _scatter_by_rank(ax, points, rcolors, *, marker, s, zorder):
    """Scatter ``points`` colored by rank; return ranks actually drawn (sorted)."""
    by_rank: dict = {}
    for alpha, val in points:
        by_rank.setdefault(_cusp_rank(alpha), []).append((float(alpha), float(val)))
    for rank in sorted(by_rank):
        pts = by_rank[rank]
        color = rcolors[rank]
        ax.scatter(
            [p[0] for p in pts],
            [p[1] for p in pts],
            marker=marker,
            s=s,
            color=color,
            edgecolors=color,
            linewidths=0.6,
            zorder=zorder,
        )
    return sorted(by_rank)


def _plot_mpl(
    curve,
    xs,
    ys,
    top,
    x_lo,
    x_hi,
    rcolors,
    all_ranks,
    title,
    theme,
    inset_center,
    inset_halfwidth,
    inset_R_max,
    inset_samples_per_unit,
    show_inset,
):
    plt = style.require_mpl()
    ink = style.theme_ink(theme)
    floor_label = rf"floor ${style.SYM['delta']} \geq \frac{{1}}{{2}}$"

    with style.mpl_context(theme):
        fig, ax = plt.subplots(figsize=(9.6, 5.8))

        # nonempty region (Delta >= delta) stays shaded.
        ax.fill_between(xs, ys, top, color=style.PALETTE["fill"], alpha=0.12, zorder=1)
        # the delta(mu) envelope.
        ax.plot(xs, ys, color=style.PALETTE["curve"], lw=2.0, zorder=4)
        # floor guide.
        ax.axhline(0.5, color=style.PALETTE["floor"], ls=":", lw=1.3, zorder=2)

        # cusps ON the curve (triangle), exceptional BELOW (x) -- both by rank.
        cusp_ranks = _scatter_by_rank(
            ax, curve.cusps, rcolors, marker="^", s=46, zorder=6
        )
        exc_ranks = _scatter_by_rank(
            ax, curve.exceptional_points, rcolors, marker="x", s=52, zorder=6
        )

        ax.set_xlabel(style.LABEL["mu"])
        ax.set_ylabel(style.LABEL["Delta"])
        ax.set_title(title)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(0.42, top)

        # exact-fraction ticks at salient cusp slopes.
        style.set_fraction_ticks(ax, "x", _salient_ticks(x_lo, x_hi))

        # annotate the prominent cusps with their exact fraction.
        _annotate_cusps(ax, curve, ink, theme)

        # ---- legends (two, so the discrete-rank decode is unambiguous) ----
        # (1) shape/curve legend.
        shape_handles = [
            style.proxy_handle(
                color=style.PALETTE["curve"],
                label=rf"${style.SYM['delta']}({style.SYM['mu']})$ curve",
                linestyle="-",
            ),
            style.proxy_handle(
                color=ink["muted"],
                label=r"cusp $\blacktriangle$ (on curve)",
                marker="^",
            ),
            style.proxy_handle(
                color=ink["muted"],
                label=r"exceptional $\times$ (below)",
                marker="x",
            ),
            style.proxy_handle(
                color=style.PALETTE["floor"],
                label=floor_label,
                linestyle=":",
            ),
        ]
        leg_shape = ax.legend(
            handles=shape_handles,
            loc="upper right",
            framealpha=0.92,
            fontsize=9.5,
            title="category",
        )
        if leg_shape.get_title():
            leg_shape.get_title().set_color(ink["ink"])
        for txt in leg_shape.get_texts():
            txt.set_color(ink["ink"])
        ax.add_artist(leg_shape)

        # (2) DISCRETE rank legend -- one swatch per rank, SAME hex as artists.
        rank_handles = [
            style.proxy_handle(
                color=rcolors[r], label=rf"rank ${r}$", marker="o", markersize=8
            )
            for r in all_ranks
        ]
        leg_rank = ax.legend(
            handles=rank_handles,
            loc="upper left",
            framealpha=0.92,
            fontsize=9.0,
            ncol=2 if len(rank_handles) > 5 else 1,
            title=r"Markov rank $r$",
        )
        if leg_rank.get_title():
            leg_rank.get_title().set_color(ink["ink"])
        for txt in leg_rank.get_texts():
            txt.set_color(ink["ink"])

        # ---- zoomed inset revealing the self-similar cusp cascade ----
        if show_inset:
            _add_zoom_inset(
                fig,
                ax,
                rcolors,
                ink,
                theme,
                inset_center,
                inset_halfwidth,
                inset_R_max,
                inset_samples_per_unit,
            )

        style.style_axes(ax, theme)
        # re-assert exact-fraction ticks (style_axes does not touch locators,
        # but be defensive in case a future change does).
        style.set_fraction_ticks(ax, "x", _salient_ticks(x_lo, x_hi))

    return fig


def _annotate_cusps(ax, curve, ink, theme):
    """Label the prominent (low-rank) cusps with their exact fraction slope."""
    prominent = {Fraction(2, 5), Fraction(1, 2), Fraction(3, 5)}
    for alpha, val in curve.cusps:
        if alpha in prominent:
            ax.annotate(
                style.m(style.latex_frac(alpha)),
                xy=(float(alpha), float(val)),
                xytext=(0, 7),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=10,
                color=ink["ink"],
                zorder=7,
            )


def _add_zoom_inset(
    fig,
    ax,
    rcolors,
    ink,
    theme,
    inset_center,
    inset_halfwidth,
    inset_R_max,
    inset_samples_per_unit,
):
    """Inset that re-plots a FINER dlp_curve around a cusp to show self-similarity."""
    center = Q(inset_center)
    half = Q(inset_halfwidth)
    lo, hi = center - half, center + half

    fine = dlp_curve(lo, hi, inset_R_max, inset_samples_per_unit)
    fxs, fys = fine.as_floats()
    f_top = max(float(v) for _, v in fine.cusps) if fine.cusps else max(fys)
    f_top = max(f_top, max(fys)) * 1.02

    axins = ax.inset_axes([0.60, 0.30, 0.36, 0.40])
    axins.set_facecolor(ink["axes_bg"])
    axins.fill_between(fxs, fys, f_top, color=style.PALETTE["fill"], alpha=0.12, zorder=1)
    axins.plot(fxs, fys, color=style.PALETTE["curve"], lw=1.4, zorder=4)
    axins.axhline(0.5, color=style.PALETTE["floor"], ls=":", lw=1.0, zorder=2)

    # finer cusps colored by rank too -- the inset rank set extends the main one.
    fine_ranks = sorted({a.denominator for a, _ in fine.cusps})
    fine_rcolors = dict(rcolors)
    missing = [r for r in fine_ranks if r not in fine_rcolors]
    if missing:
        extra = style.rank_color_map(fine_ranks)
        for r in missing:
            fine_rcolors[r] = extra[r]
    _scatter_by_rank(axins, fine.cusps, fine_rcolors, marker="^", s=22, zorder=6)
    _scatter_by_rank(
        axins, fine.exceptional_points, fine_rcolors, marker="x", s=24, zorder=6
    )

    axins.set_xlim(float(lo), float(hi))
    y_floor = 0.49
    axins.set_ylim(y_floor, f_top)

    # exact-fraction tick at the cusp center inside the inset.
    style.set_fraction_ticks(axins, "x", [lo, center, hi])
    axins.tick_params(labelsize=7.5, colors=ink["ink"])
    for side in ("top", "right"):
        axins.spines[side].set_visible(True)
        axins.spines[side].set_color(ink["edge"])
    for side in ("left", "bottom"):
        axins.spines[side].set_color(ink["edge"])
    axins.grid(True, color=ink["grid"], alpha=float(ink["grid_alpha"]) * 0.7, linewidth=0.4)
    axins.set_axisbelow(True)

    axins.set_title(
        rf"zoom @ ${style.SYM['mu']}={style.latex_frac(center)}$ "
        rf"(${style.SYM['Rmax']}={fine.R_max}$)",
        fontsize=8.5,
        color=ink["ink"],
        pad=2.0,
    )

    # connector box drawn on the main axes.
    ax.indicate_inset_zoom(axins, edgecolor=ink["edge"], alpha=0.7)
