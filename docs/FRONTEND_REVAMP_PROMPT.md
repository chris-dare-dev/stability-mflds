# Prompt — Visualization / Frontend Revamp for `bridgeland_stability`

## For a fresh Opus 4.8 "ultracode" session

> Paste this whole document as the opening prompt of a new session. It is
> self-contained: it assumes no memory of the conversation that produced it.
> **Scope:** redesign and harden the plotting / frontend layer only. Do **not**
> change the mathematics, the exact-arithmetic core, or any algorithm result.

---

## 0. Mission

`bridgeland_stability` is a Python package (in this repo) that computes, in
**exact rational arithmetic**, three families of objects from the theory of
moduli of sheaves and Bridgeland stability, and plots them:

1. the **Drézet–Le Potier curve** `δ(μ)` on P² — a fractal upper envelope with
   sharp upward **cusps** at a dense set of Markov-fraction slopes;
2. **Bridgeland stability walls** — families of **nested concentric
   semicircles** in the `(s, t)` upper half-plane, meant to be **deep-zoomed**;
3. the **threefold tilt-BG boundary** `α_crit(β)` — a curve with a **vertical
   asymptote / degenerate point** where `ch₁ᵝ = 0`.

The plotting lives in [`bridgeland_stability/viz.py`](../bridgeland_stability/viz.py),
with a **matplotlib** backend (static PNG, the subject of most complaints) and a
**plotly** backend (interactive HTML). The package's whole value proposition is
*exact* data — so the visuals must (a) typeset real mathematics, (b) preserve
precision under deep zoom, and (c) serve **both** interactive web exploration
**and** publication-quality static figures.

Your job: bring the frontend up to that standard. The detailed critique (§2),
recommended stack (§3), and task list with acceptance criteria (§4) are below.

---

## 1. Current state & how to reproduce the figures

Public plotting API (in `viz.py`), all returning the underlying figure object:

```python
from bridgeland_stability import viz
viz.plot_dlp_curve(mu_min=0, mu_max=1, R_max=40, backend="matplotlib"|"plotly", show=, html_path=)
viz.plot_walls(v, surface, s_range=, t_max=, rank_bound=, walls=None, max_walls=30, backend=, ...)
viz.plot_threefold_bg(boundary, backend=, ...)
```

Data objects the plotters consume (all exact `fractions.Fraction` internally):
`DLPCurve` (`.mus`, `.deltas`, `.bundles`, `.cusps`, `.exceptional_points`),
`Wall` / `ActualWall` (`.center`, `.radius_sq`, `.radius`, `.subobject`,
`.quotient`), `BGBoundary` (`.betas`, `.alphas`, `.bg_proven`,
`.degenerate_betas`), `ChernChar` (`.r`, `.c`, `.e`), `Bundle` (`.r`, `.c1`,
`.ch2`, `.slope`, `.discriminant`).

Reproduce the exact figures that were reviewed:

```bash
pip install -e ".[dev]"          # plotly + matplotlib + pytest
python examples/demo.py --figures # writes figures/*.{html,png}
# the reviewed matplotlib PNGs:
python - <<'PY'
from fractions import Fraction as F
from bridgeland_stability import viz, actual_walls, ChernChar, P2, P3
from bridgeland_stability.threefold import bg_boundary_curve, ThreefoldChern
viz.plot_dlp_curve(mu_min=0, mu_max=1, R_max=30, backend="matplotlib").savefig("/tmp/01_dlp.png", dpi=120, bbox_inches="tight")
v=ChernChar(1,0,F(-4)); viz.plot_walls(v,P2,s_range=(-6,1),t_max=4,walls=actual_walls(v,P2),backend="matplotlib").savefig("/tmp/02_walls.png",dpi=120,bbox_inches="tight")
viz.plot_walls(ChernChar(2,0,F(-1,4)),P2,s_range=(-3,3),t_max=3,rank_bound=4,backend="matplotlib").savefig("/tmp/03_numerical.png",dpi=120,bbox_inches="tight")
bd=bg_boundary_curve(ThreefoldChern(2,F(0),F(1),F(0)),P3,beta_range=(-2,2),N=300); viz.plot_threefold_bg(bd,backend="matplotlib").savefig("/tmp/04_bg.png",dpi=120,bbox_inches="tight")
PY
```

You can also serve a gallery and inspect it in a browser preview: there is a
ready page at `figures/mpl_review/index.html` and a `.claude/launch.json` config
named `figure-review` (static server on port 8765).

---

## 2. The critique — UI/UX shortcomings to fix

Grouped by theme; each has a severity and the intended direction. This list was
produced from a hands-on visual review of the four rendered figures plus an
independent data-viz review; it is the spec for what "fixed" means.

### A. Mathematical typography (no LaTeX) — **high**

1. **ASCII math words instead of symbols.** Axes/titles literally read
   `slope mu`, `discriminant Delta (CH)`, `delta(mu)`, `alpha`, `beta`, `t`,
   `s`, `P^2`, `P^3`. An algebraic-geometry audience expects rendered `μ`, `Δ`,
   `α`, `β`, `δ`, `χ`, and `ℙ²`/`ℙ³`. The ASCII forms look unfinished and are
   unacceptable in a paper. → Use mathtext (`r"$\mu$"`, `r"$\Delta$"`,
   `r"$\alpha_{\mathrm{crit}}(\beta)$"`, `r"$\mathbb{P}^2$"`) everywhere; mirror
   in plotly via its MathJax LaTeX support.
2. **Python `repr` leaking into titles.** Figs 2–3 title with
   `v=ChernChar(r=2, c=0, e=-1/4)` — a raw dataclass repr that exposes internals
   and is meaningless to a reader. → Add a `ChernChar.to_latex()` / `__format__`
   producing e.g. `$v=(r,c,\mathrm{ch}_2)=(2,0,-\tfrac14)$`.
3. **Decimal ticks where exact rationals are the point.** Ticks read
   `0.0, 0.2, 0.4 …` and `-0.5, 0.5`, but the salient locations are rationals
   (cusps at Markov fractions `1/2, 2/5, 3/5`; the `Δ=1/2` line; the degenerate
   `β`). → A `FuncFormatter` rendering key positions as fractions, and/or
   annotate cusps/asymptotes with their exact fraction.

### B. Color & legends (the "heat-scale" issue) — **high**

4. **Continuous viridis colorbar for a discrete rank** (fig 1). Integer sheaf
   rank `∈ {1,2,5,13,29,…}` is mapped to a smooth viridis ramp with a continuous
   colorbar — implying a nonexistent continuum; you cannot read an exact integer
   off a marker. → Discrete/qualitative colormap with a **categorical legend**
   (one swatch per rank), or a stepped `BoundaryNorm` bar on integer boundaries,
   or annotate markers with their rank.
5. **Legend swatch ≠ plotted color** (fig 1). The legend shows one green
   triangle and one purple ✗, but the actual markers are recolored by the
   viridis rank scale, so the key lies. → Decouple *category* (marker shape,
   neutral edge) from *rank* (color); bind legend handles to the same color
   objects as the artists; unit-test that they match.
6. **Rainbow default cycle, no legend** (fig 3, numerical walls). ~12–15 walls
   in the repeating `tab10` cycle with no legend/labels — decorative, not
   informative; no curve maps to a wall/radius/destabilizer. → Order walls by a
   meaningful scalar (radius / destabilizer) and color with a sequential
   perceptually-uniform colormap + a labeled colorbar, or a compact legend.
7. **No colorblind-safe palette.** → Use **Okabe–Ito** for categorical series
   and **viridis/cividis** only for genuinely continuous fields.

### C. Theming, backgrounds, "plain design" — **high/medium**

8. **Opaque white figure background** (all figs). Jarring on a dark web UI and
   prevents clean compositing on slides/reports. → `savefig(transparent=True)`
   and/or theme-matched facecolor; ship matched **light/dark `.mplstyle`** sheets
   with ink colors that keep contrast on both; plotly
   `paper_bgcolor/plot_bgcolor='rgba(0,0,0,0)'` + per-theme template.
9. **No gridlines** (all figs). Reading `Δ` at a given `μ`, or wall radii, or
   `α` at a `β`, forces eyeballing off distant ticks. → Light major gridlines
   (`alpha≈0.3`), minor gridlines on the wall plots, subordinate to the data.
10. **Default small fonts, no DPI bump** (all figs). ~10 pt text won't survive
    downscaling into a column or read well on high-DPI screens. → A style with
    ~12–14 pt base, `savefig dpi 200–300`, `constrained_layout=True`.
11. **Overall "default matplotlib" look.** No coherent style identity. → A
    bundled house style (consider SciencePlots) so every figure is publication-
    ready by default.

### D. Geometry & framing — **high/medium**

12. **Asymptote blows out the y-axis** (fig 4). `α_crit(β)` spikes to ~75 near
    `β=0` (where `ch₁ᵝ→0`), so autoscale crushes the meaningful band (`α≈0–3`)
    into the bottom ~7 % of the panel — it reads as a flat line. → Break the
    curve with `np.nan` at the degenerate `β`; cap `ylim` (e.g. 0–6); draw the
    asymptote as an annotated dashed `axvline`; `axvspan`-shade the
    degenerate/forbidden region; arrow-annotate "→ ∞". Use `BGBoundary
    .degenerate_betas` (already computed).
13. **Largest semicircle clipped / `t_max` not data-driven** (fig 2). `t_max` is
    a fixed argument, so a wall whose apex exceeds it is sliced flat and the
    nested structure breaks at the top. → Derive default `t_max` from the
    largest wall radius (`max apex × ~1.1`).
14. **Aspect handling.** The walls plots *do* set `set_aspect("equal")`
    (matplotlib) and `scaleanchor` (plotly), so circles are round — **verify
    this is preserved** after restyling, ensure it interacts well with
    data-driven limits (`adjustable='datalim'` to keep circles round while
    filling the axes), and do **not** aspect-lock the `α_crit(β)` plot (β and α
    are different quantities).
15. **Fractal cusp structure undersampled** (fig 1). The DLP curve shows only a
    few coarse cusps with straight segments between; the dense self-similar
    structure that defines the object is invisible at default sampling. →
    Adaptively refine sampling near Markov fractions; add a **zoomed inset** on
    one cusp to show self-similarity; expose a density / `R_max` control.
16. **No deep-zoom affordance into the wall cluster** (fig 3). Walls pinch
    together near the origin into an unreadable smear — exactly the region the
    use case wants to explore. → Add an `inset_axes` / `indicate_inset_zoom`
    near-origin inset for the static figure; ensure the plotly version supports
    smooth zoom with preserved aspect (see §3 on SVG vs WebGL).
17. **Unexplained `Δ=0.5` guide line** (fig 1). A dotted horizontal line at
    `Δ=1/2` has no label or legend entry. → Label it inline/legend as the floor
    `δ ≥ 1/2`.
18. **Unstyled `-- PROVEN` status with a double hyphen** (fig 4). The
    proof-status flag is appended to the title as `-- PROVEN` (ASCII `--`). It is
    semantically important (proven vs conjectural BG inequality). → Render as a
    styled badge (green PROVEN / amber CONJECTURAL) or a separate line with a
    proper dash; drive from `BGBoundary.bg_proven`.

### E. Interactivity & surfacing exact data — **high/medium**

19. **Deep-zoom precision: use SVG plotly traces, not WebGL.** This is the
    headline interactive requirement. Plotly's `go.Scattergl` (WebGL) renders
    coordinates in **float32** (~1e-7), with documented jitter/pixelation under
    deep zoom (plotly.js #6636/#6820, plotly.py #2881/#4440). For concentric
    semicircles deep-zoomed many decades, that **visibly breaks** the geometry
    and discards the float64 fidelity the exact engine produced. → Draw all
    curves/semicircles with **`go.Scatter` (SVG, float64)**; gate `Scattergl`
    behind an explicit large-N flag (not needed here — data is sparse). Add a
    regression test that zooms to ~1e-10 of axis range.
20. **No interactive readouts of exact data** (all figs). Flat PNGs with no
    hover, no exact-rational tooltips, no zoom cues — yet the value is exact
    rationals at cusps, wall centers, and the degenerate point. → In the plotly
    variant, attach `hovertemplate`s exposing exact values **as fractions**
    (e.g. `s₀ = -5/2`, `R² = 9/4`, destabilizer `(r,c,ch₂)`), not floats.

### F. Publication export — **medium**

21. **No vector / font-embedding path.** Raster PNG pixelates in print, and
    matplotlib's default PDF Type-3 fonts are often un-editable / journal-
    rejected. → A `save_publication()` emitting **PDF + SVG** with
    `pdf.fonttype=42`, `ps.fonttype=42`, `svg.fonttype='none'`, `dpi≥300`,
    `bbox_inches='tight'`; PNG only for previews.

---

## 3. Recommended frontend / plotting stack (2026)

The reviewed conclusion: **keep the matplotlib (static) + plotly (interactive)
split, but harden both**, and add targeted complements. Rationale: exact-rational
*sparse* curves that need real LaTeX, lossless deep-zoom, and dual web+print
output uniquely reward a **float64 SVG** interactive path and a **true-LaTeX
vector** static path — and penalize every WebGL / rasterization /
grammar-of-graphics option.

**Use (core):**
- **Matplotlib 3.10+** — static core. `mathtext` by default (no TeX install
  needed); opt-in `text.usetex=True` for camera-ready. float64 vector
  PDF/SVG/EPS export.
- **SciencePlots (v2.1.x)** — one-line journal styling layer on matplotlib
  (Computer Modern, colorblind-safe palettes, sensible ticks). Note: its `usetex`
  path needs a LaTeX install (MacTeX on macOS); keep `mathtext` as fallback.
- **Plotly 6.x** — interactive core, **SVG `go.Scatter` only** (never
  `go.Scattergl`). LaTeX-in-labels via MathJax. Exports PDF/SVG/PNG via Kaleido.

**Consider (complementary, scoped):**
- **PGFPlots / TikZ (1.18.2+)** — for 1–2 *hero* camera-ready figures typeset by
  the paper's own LaTeX engine (perfect typographic match). Slow to author; keep
  matplotlib for bulk static output.
- **D3.js v7 (+ KaTeX)** — *only* if a polished deep-zoom semicircle UI becomes a
  product feature: true float64 semantic zoom with per-level re-tessellation, far
  beyond any high-level lib. Large build effort; defer unless Plotly SVG zoom
  proves insufficient.
- **Manim CE (v0.20.1)** — *only* for **animating** wall-crossing (parameter
  sweeps where semicircles deform/merge) into MP4/GIF for talks/supplementary
  media. Not a plotting library; keep out of the figure pipeline.
- **Math typesetting:** **KaTeX** (fast, for any custom web UI) and **MathJax**
  (broadest; what plotly uses).

**Skip (with reasons):** `go.Scattergl`/three.js/raw WebGL (float32 — the core
disqualifier for deep zoom); Bokeh / HoloViews / Datashader (built for
large-data rasterization, irrelevant to sparse exact curves and fights zoom
sharpness); Altair / Vega-Lite / Observable Plot (no native LaTeX;
grammar-of-graphics is a poor fit for hand-crafted cusp/asymptote curves);
ECharts (no native math, no edge over plotly); mpld3 (dormant); Asymptote
(niche/slow — PGFPlots covers it); seaborn-as-engine (styling only).

---

## 4. Task list & acceptance criteria

Implement the following. Keep changes confined to the viz layer + a small style
package; do not touch the math modules or their tests.

**Typography & data formatting**
- [ ] All axis labels, titles, tick annotations use mathtext/LaTeX symbols
      (`μ, Δ, α, β, δ, χ, ℙ²`). No ASCII `mu`/`Delta`/`alpha`/`beta`/`P^2`
      anywhere in rendered output.
- [ ] `ChernChar.to_latex()` (and `Bundle.to_latex()`); titles use it — no
      `repr` leakage.
- [ ] Fraction tick/annotation formatter for the salient rational locations.

**Color, legend, theme**
- [ ] Rank shown as discrete categorical (legend or stepped bar), never a smooth
      colorbar; legend swatches provably equal plotted colors (add a unit test).
- [ ] Okabe–Ito categorical palette centralized in one module-level dict reused
      by artists *and* legend handles; viridis/cividis only for continuous.
- [ ] Bundled light **and** dark `.mplstyle` sheets; default export
      `transparent=True`; plotly per-theme template with transparent bg.
- [ ] Gridlines (major; minor on wall plots), ~12–14 pt base font,
      `constrained_layout`, `dpi≥200`.

**Geometry & framing**
- [ ] `plot_threefold_bg`: NaN-break at degenerate `β`, capped `ylim`, annotated
      asymptote (`axvline`+label), shaded forbidden region, "→ ∞" annotation;
      proof-status **badge** (green/amber) from `bg_proven`.
- [ ] `plot_walls`: data-driven default `t_max` so no semicircle is clipped;
      verify equal-aspect circles survive restyling; near-origin zoom inset for
      dense wall families.
- [ ] `plot_dlp_curve`: adaptive cusp sampling + a self-similarity zoom inset;
      labeled `δ=1/2` floor line.

**Interactivity & export**
- [ ] Plotly: `go.Scatter` (SVG) for all curves/semicircles; `hovertemplate`s
      surfacing **exact fractions** (wall center/radius², cusp `μ`/`δ`,
      destabilizer Chern data); deep-zoom regression check.
- [ ] `viz.save_publication(fig, path)` → PDF+SVG with `pdf.fonttype=42`,
      `ps.fonttype=42`, `svg.fonttype='none'`, `bbox_inches='tight'`.

**Definition of done**
- [ ] The 4 figures regenerated and visually reviewed (use the preview
      gallery + a browser screenshot) against this critique.
- [ ] New `tests/test_viz.py`: backends import and produce a figure without a
      display; legend-color == artist-color; titles contain no `ChernChar(`;
      labels contain LaTeX (`$`); threefold y-axis is capped; `save_publication`
      writes a valid PDF and SVG.
- [ ] The existing **42 tests still pass**; the math/core modules are unchanged
      (`git diff` touches only `viz.py`, a new style package, `tests/test_viz.py`,
      and the `to_latex` helpers).
- [ ] `viz` stays optional: importing `bridgeland_stability` must not require
      matplotlib/plotly (lazy imports inside `viz`).

---

## 5. Constraints & guardrails

- **Do not change the mathematics** or any computed value. This is a pure
  presentation-layer task. The discriminant convention is Coskun–Huizenga
  (`Δ = ½μ² − ch₂/(rd)`); see `docs/CORRECTIONS.md`.
- Keep the public `viz` function signatures backward-compatible where reasonable
  (add keyword args; don't remove existing ones without need).
- `matplotlib`/`plotly` are **optional extras** (`[viz]`); guard imports.
- Verify by *rendering and looking* — regenerate the gallery, screenshot it, and
  check each critique item is addressed. Do not declare done from code alone.

## 6. Suggested orchestration (this is an ultracode session)

A reasonable workflow: (1) a short design phase deciding the style system
(palette dict, light/dark `.mplstyle`, mathtext helpers, `to_latex`); (2) a
fan-out implementing the three plot families in parallel against that shared
style; (3) an adversarial visual-review pass that re-renders and screenshots
each figure and checks it against §2/§4; loop until clean. Keep the math tests
green throughout.

---

## 7. Sources

Plotly perf & WebGL precision: plotly.com/python/performance,
github.com/plotly/plotly.py issues 2881 & 4440, plotly.js issues 6636 & 6820,
webglfundamentals.org WebGL precision issues. Matplotlib text/usetex/pgf/fonts:
matplotlib.org/stable/users/explain/text/{mathtext,usetex,pgf,fonts}.html.
SciencePlots: github.com/garrettj403/SciencePlots. PGFPlots: github.com/pgf-tikz/pgfplots.
Okabe–Ito palette references (sci-draw, conceptviz, easystats `see`). Discrete
vs continuous color & `BoundaryNorm`/`ListedColormap`: matplotlib colorbar/colors
docs. Asymptote/broken-axis handling: matplotlib `axline`, `broken_axis`
examples. Manim: docs.manim.community. D3: d3js.org.
