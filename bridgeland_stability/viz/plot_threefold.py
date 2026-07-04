"""Threefold tilt-BG boundary plot.

Renders ``alpha_crit(beta)`` for a threefold, breaking the curve at each
degenerate beta (where ``ch_1^beta = 0``), capping the y-axis so the
meaningful band is readable, annotating the vertical asymptote(s), shading the
degenerate/forbidden locus, and showing the BG proof status as a styled badge.
Shares the house style in :mod:`bridgeland_stability.viz.style`.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

from ..threefold import BGBoundary
from . import style


def _split_at_degenerate(
    betas: List[float],
    alphas: List[float],
    degenerate_betas: List[float],
):
    """Return ``(xs, ys)`` with ``nan`` inserted across each degenerate beta.

    ``boundary.betas`` already *excludes* the degenerate beta (sampling skips
    it), so a degenerate value lies strictly between two consecutive sampled
    betas.  Inserting ``nan`` (in both x and y) there severs the polyline so the
    plotted curve does not draw a spurious near-vertical connector across the
    asymptotic blow-up.
    """
    import math

    degs = sorted(float(d) for d in degenerate_betas)
    xs: List[float] = []
    ys: List[float] = []
    for i, (bx, by) in enumerate(zip(betas, alphas)):
        if i > 0:
            prev = betas[i - 1]
            lo, hi = (prev, bx) if prev <= bx else (bx, prev)
            # break wherever a degenerate beta sits between the two samples
            if any(lo < d < hi for d in degs):
                xs.append(math.nan)
                ys.append(math.nan)
        xs.append(float(bx))
        ys.append(float(by))
    return xs, ys


def plot_threefold_bg(
    boundary: BGBoundary,
    *,
    backend: str = "plotly",
    theme: str = "light",
    y_cap: float = 6.0,
    show: bool = False,
    html_path: Optional[str] = None,
):
    """Plot the tilt-stability BG boundary ``alpha_crit(beta)`` for a threefold.

    Parameters
    ----------
    boundary:
        Sampled :class:`~bridgeland_stability.threefold.BGBoundary`.  Its
        ``betas`` exclude each degenerate beta (where ``ch_1^beta = 0``);
        consecutive betas straddle every entry of ``degenerate_betas``.
    backend:
        ``"plotly"`` (interactive HTML, exact fractions in hover) or
        ``"matplotlib"`` (publication-quality static figure).
    theme:
        ``"light"`` or ``"dark"``.
    y_cap:
        Upper alpha limit.  The asymptote spikes toward ``+inf`` near each
        degenerate beta; capping keeps the meaningful ``alpha in [0, ~3]`` band
        legible instead of crushing it into a sliver.
    show, html_path:
        Forwarded to :func:`style.emit_plotly` for the plotly backend.

    Returns the backend figure object (matplotlib ``Figure`` or plotly
    ``Figure``).
    """
    proven = bool(boundary.bg_proven)
    status = "PROVEN" if proven else "CONJECTURAL"
    status_color = style.PALETTE["proven" if proven else "conjectural"]

    variety = style.latex_variety_name(boundary.threefold_name)
    title = f"Tilt-BG boundary {style.m(style.SYM['alpha_crit'])} on {variety}"

    ch1b_label = style.m(style.SYM["ch1b"] + " = 0")  # $\mathrm{ch}_1^{\beta} = 0$

    betas = list(boundary.betas)
    alphas = list(boundary.alphas)
    degs = sorted(float(d) for d in boundary.degenerate_betas)

    # NaN-break the polyline across every degenerate beta.
    xs, ys = _split_at_degenerate(betas, alphas, degs)

    cap = float(y_cap)

    # ------------------------------------------------------------------
    # Plotly
    # ------------------------------------------------------------------
    if backend == "plotly":
        go = style.require_plotly()
        fig = go.Figure()
        style.apply_plotly_theme(fig, theme)

        # forbidden/degenerate shaded band around each degenerate beta.
        span = 0.02
        if len(betas) >= 2:
            span = max(span, 0.4 * (max(betas) - min(betas)) / max(1, len(betas)))
        for d in degs:
            fig.add_vrect(
                x0=d - span,
                x1=d + span,
                fillcolor=style.PALETTE["forbidden"],
                opacity=0.12,
                line_width=0,
                layer="below",
            )

        # NOTE: go.Scatter only (float64 / SVG); never Scattergl.  Exact
        # rationals are not available on the float grid, so the hover carries
        # high-precision floats with the real LaTeX symbols.
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                connectgaps=False,
                line=dict(color=style.PALETTE["boundary"], width=2.4),
                name=style.m(style.SYM["alpha_crit"]),
                hovertemplate=(
                    style.m(style.SYM["beta"]) + " = %{x:.5f}<br>"
                    + style.m(style.SYM["alpha"]) + " = %{y:.5f}"
                    + "<extra>" + style.m(style.SYM["alpha_crit"]) + "</extra>"
                ),
            )
        )

        for d in degs:
            fig.add_vline(
                x=d,
                line=dict(color=style.PALETTE["asymptote"], dash="dash", width=1.6),
                annotation_text=ch1b_label,
                annotation_position="top",
                annotation_font=dict(color=style.PALETTE["asymptote"]),
            )
            # blow-up arrow near the asymptote.
            fig.add_annotation(
                x=d,
                y=cap * 0.86,
                text=style.m(r"\to \infty"),
                showarrow=True,
                arrowhead=2,
                arrowcolor=style.PALETTE["asymptote"],
                arrowwidth=1.6,
                ax=28,
                ay=22,
                font=dict(color=style.PALETTE["asymptote"], size=14),
            )

        fig.update_layout(
            title=title,
            xaxis_title=style.LABEL["beta"],
            yaxis_title=style.LABEL["alpha"],
        )
        # capped y-axis; no aspect lock (beta and alpha are different quantities).
        fig.update_yaxes(range=[0.0, cap])
        if betas:
            fig.update_xaxes(range=[min(betas), max(betas)])

        # proof-status badge as a styled annotation in paper coords.
        fig.add_annotation(
            x=0.985,
            y=0.97,
            xref="paper",
            yref="paper",
            xanchor="right",
            yanchor="top",
            text=f"<b>{status}</b>",
            showarrow=False,
            font=dict(color="white", size=13),
            bgcolor=status_color,
            bordercolor=status_color,
            borderpad=5,
            opacity=0.95,
        )

        style.emit_plotly(fig, show, html_path)
        return fig

    # ------------------------------------------------------------------
    # Matplotlib
    # ------------------------------------------------------------------
    plt = style.require_mpl()
    with style.mpl_context(theme):
        fig, ax = plt.subplots(figsize=(8.5, 5.2))

        # forbidden/degenerate shaded band around each degenerate beta.
        span = 0.02
        if len(betas) >= 2:
            span = max(span, 0.4 * (max(betas) - min(betas)) / max(1, len(betas)))
        for d in degs:
            ax.axvspan(
                d - span,
                d + span,
                color=style.PALETTE["forbidden"],
                alpha=0.12,
                lw=0,
                zorder=1,
            )

        # the NaN-broken boundary curve.
        ax.plot(
            xs,
            ys,
            color=style.PALETTE["boundary"],
            lw=2.4,
            zorder=5,
            label=style.m(style.SYM["alpha_crit"]),
        )

        for d in degs:
            ax.axvline(
                d,
                color=style.PALETTE["asymptote"],
                ls="--",
                lw=1.5,
                zorder=3,
            )
            ax.annotate(
                ch1b_label,
                xy=(d, cap),
                xytext=(0, -4),
                textcoords="offset points",
                ha="center",
                va="top",
                color=style.PALETTE["asymptote"],
                fontsize=10.5,
                zorder=6,
            )
            # blow-up arrow indicating alpha_crit -> +inf at the asymptote.
            ax.annotate(
                style.m(r"\to \infty"),
                xy=(d, cap * 0.82),
                xytext=(span + 0.18, cap * 0.62),
                textcoords="data",
                ha="left",
                va="center",
                color=style.PALETTE["asymptote"],
                fontsize=12.5,
                arrowprops=dict(
                    arrowstyle="->",
                    color=style.PALETTE["asymptote"],
                    lw=1.5,
                ),
                zorder=7,
            )

        # capped y-axis; NO aspect lock (beta and alpha are different quantities).
        ax.set_ylim(0.0, cap)
        if betas:
            ax.set_xlim(min(betas), max(betas))

        ax.set_xlabel(style.LABEL["beta"])
        ax.set_ylabel(style.LABEL["alpha"])
        ax.set_title(title)

        leg = ax.legend(loc="upper left", framealpha=0.0)
        if leg is not None:
            ink = style.theme_ink(theme)
            for txt in leg.get_texts():
                txt.set_color(ink["ink"])

        style.add_badge(ax, status, color=status_color, theme=theme)
        style.style_axes(ax, theme)

    return fig
