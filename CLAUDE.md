# CLAUDE.md — bridgeland_stability

Guidance for Claude Code (and any AI coding agent) working in this repository.

## What this is

`bridgeland_stability` is a pure-Python, **exact-arithmetic** toolkit for
research-level algebraic geometry: Drézet–Le Potier curves, Bogomolov–Gieseker
boundaries, Bridgeland stability walls in the (s, t) upper half-plane, and K3
Mukai lattices. Every invariant — slope, discriminant, wall center, radius² — is
computed with `fractions.Fraction`; floats appear only for geometry and
plotting. The package explicitly **corrects substantive errors in its originating
project brief**; [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) documents each
correction with a primary-source citation and an independent exact-`Fraction`
verification. That correctness culture is the point of the project — respect it.

## Commands

```bash
pip install -e .                    # core — pure standard library, zero runtime deps
pip install -e ".[viz]"             # + plotly>=5.18 + matplotlib>=3.7 for figures
pip install -e ".[dev]"             # + pytest (includes viz)
pytest -q                           # run the suite (42 tests)
python examples/demo.py             # print the corrected validated-value table
python examples/demo.py --figures   # write figures/*.{html,png}  (untracked build output)
# figure-review gallery (use `python` on Windows):
python -m http.server 8765 --directory figures/mpl_review
sh scripts/install-hooks.sh         # once per clone: activate .githooks/pre-commit
```

## Module map

| Module | Responsibility |
|---|---|
| `chern` | `ChernChar(r, c, e)`: H-numerical Chern character on a surface. `slope(d)`, `discriminant(d)` (CH normalization), `discriminant_brief(d)` (doubled — comparison only), `twist(s, d)`, `central_charge(s, t, d)` → `(float, float)`, `to_latex(name)`. Plus the `Q()` coerce-to-`Fraction` helper and `Number = Union[int, Fraction]`. |
| `varieties` | `Surface` / `Threefold` frozen dataclasses + catalog: `P2`, `P1xP1`, `P3`, `QUADRIC3`, `QUINTIC`, `BLOWUP_P3_POINT`; factories `K3(h2)`, `hirzebruch(n)`, `abelian_surface(h2)`, `abelian_threefold(d3)`, `fano_picard_one(name, d3)`. |
| `exceptional` | Algorithm 1. `Bundle(r, c1, ch2)` with `.discriminant` / `.discriminant_brief` properties and constructors `Bundle.O(n)`, `Bundle.T_minus1()`, `Bundle.from_slope(alpha)`. Markov/ε-recursion: `enumerate_exceptional`, `exceptional_mediant`, `markov_numbers`, `is_markov_number`, `is_exceptional`, `chi` (Riemann–Roch). |
| `dlp` | Algorithm 2. DLP curve as a fractal upper envelope: `delta(mu, bundles)`, `dlp_curve(...)` → `DLPCurve`, `control_interval_halfwidth`, `moduli_nonempty`. |
| `walls` | Algorithm 3. `numerical_wall(v, w, d)` (exact primitive), `compute_walls` (dense numerical set — *not* certified), `walls_from_subobjects`, `actual_walls` / `actual_walls_complete` (certified-necessary, finite set). Types: `Wall`, `VerticalWall`, `ActualWall` (`.center`, `.radius_sq` exact; `.radius` is `float`). |
| `bg_check` | Algorithm 4. `check_bg_surface(...)`; re-exports `check_bg_threefold`. |
| `threefold` | Algorithm 5. `ThreefoldChern(r, a1, a2, a3)`, `bmt_Q` / `bmt_Q_at` (return `Fraction`), `alpha_crit(v, beta, d3)` → `Optional[float]`, `bg_boundary_curve(...)` → `BGBoundary`. `Threefold.bg_proven` flags proven vs conjectural. |
| `mukai` | K3 Mukai lattice. `MukaiVector(r, l, s)`, `pairing`, `self_pairing`, `moduli_dim`, `is_spherical`, `is_isotropic`, `classify_wall` → `WallClassification` (Bayer–Macrì Thm 5.7). `solve_binary_quadratic`/`saturated_basis`/`hyperbolic_witnesses` (G10 lattice solver). |
| `rigor` | Provenance lattice (E1-M2/G5), stdlib-only. `Rigor(IntEnum)` total order `PROVEN=3 > CONJECTURAL=2 > HEURISTIC=1 > UNKNOWN=0` (`min` = meet); frozen `Certificate(rigor, hypotheses: tuple, citations: tuple, note)`; `meet(*certs)` = min-rigor with order-preserving unioned hypotheses/citations; shared `UNKNOWN_CERTIFICATE` default. Lets every downstream wall/verdict advertise theorem-vs-conjecture without changing any value. |
| `_latex` | Stdlib-only LaTeX helpers (`latex_frac`, `latex_sqrt`, `latex_chern`, …). No matplotlib/plotly — keep it that way. |
| `viz/` | Optional extra. Public API: `plot_dlp_curve`, `plot_walls`, `plot_threefold_bg`, `save_publication`. Style exports: `OKABE_ITO`, `PALETTE`, `INK`, `SYM`, `LABEL`. Submodules `plot_dlp`, `plot_walls`, `plot_threefold`, `style`; bundled `bridgeland-light` / `bridgeland-dark` `.mplstyle`. |

## Invariants — do not break these

### 1. Exact arithmetic: never introduce float into the exact core
All core invariants stay `fractions.Fraction`: `ChernChar.slope` /
`.discriminant`, `Wall.center` / `.radius_sq`, `delta()`, `bmt_Q()`, the
`MukaiVector` pairings. Float is permitted **only** at the outermost
display/geometry layer: `Wall.radius` and `ActualWall.radius` (properties),
`central_charge()` (returns `(float, float)`), the `bg_boundary_curve()`
sampling grid, `control_interval_halfwidth()`, `alpha_crit()`'s final
`math.sqrt`, and everything in `viz/`.

### 2. The Coskun–Huizenga discriminant convention is fixed
The standard throughout: `Δ = ½·μ² − ch₂/(r·d)`, `d = H²`. `discriminant_brief`
(= 2Δ) exists only for comparison with the original brief — never make it the
default. Every DLP, wall, and BG formula here is stated in CH normalization. Do
not change this convention.

### 3. Every mathematical change requires two-way verification
Before changing any formula, output value, or constant: (1) recompute it with
exact `fractions.Fraction` arithmetic, and (2) cite the primary source (read the
paper, not from memory). [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) is the
template. Math correctness is safety-critical here — this package exists because
it corrected a brief that got the mathematics wrong in several independent ways.

### 4. Pinned test values are ground truth
The 42 tests pin exact values checked against the literature. If code produces a
different answer, **the code is wrong** — do not edit a test to match a new
value without completing step 3 above. A subset:
- `delta(1/2) = 5/8`, `delta(1/3) = 5/9`, `delta(1/4) = 21/32`, `delta(2/5) = 13/25`
- P²[2] wall: center `−5/2`, radius `3/2` (ABCH, arXiv:1203.0316 §9)
- Gieseker (outermost) wall of P²[n]: center `−(2n+1)/2`, radius `(2n−1)/2`
- P³ null-correlation bundle (2,0,1,0): `alpha_crit(β=1/2) = √3` (not √29/4)
- K3 structure sheaf: `v(O) = (1,0,1)`, `⟨v,v⟩ = −2` (not `(1,0,−1)`)
- Ranks 3 and 4 are **not** Markov numbers — no rank-3 or rank-4 exceptional bundle exists

### 5. viz is an optional extra — keep matplotlib/plotly imports lazy
`import bridgeland_stability` must succeed with **zero** runtime dependencies.
All matplotlib/plotly imports live inside `viz` functions, gated by
`style.require_mpl()` / `style.require_plotly()` — never at module top level
anywhere in `viz/`. `_latex.py` and `viz/style.py` are stdlib-only at import
time; keep them so. After any `viz/` change, confirm the core still imports
without viz installed.

## Conventions

- **`d = H²`** for surfaces (`surface.d`); **`d3 = H³`** for threefolds
  (`threefold.d3`). Pass the correct degree — methods do not auto-detect it.
- `Bundle.discriminant` is a **property** (P² only, `d=1` implicit).
  `ChernChar.discriminant(d)` requires `d` explicitly and is surface-agnostic.
  This asymmetry is a common bug surface.
- All data objects are **frozen dataclasses** — create new instances, don't mutate.
- `actual_walls` / `actual_walls_complete` is the entry point for certified wall
  enumeration; `compute_walls` returns the dense, *uncertified* numerical set.
- `BLOWUP_P3_POINT.bg_proven = False` — the BG inequality is not proven there, so
  Algorithm 5 output is not rigorous for it.
- Tests and examples use `from fractions import Fraction as F`.

## Reference docs

- **[`docs/CORRECTIONS.md`](docs/CORRECTIONS.md)** — read before touching any
  formula or value. Every correction, with exact-`Fraction` evidence and a
  primary-source citation. **This is the only document under `docs/` that is
  tracked in git.**

### Working documents are NOT committed

Planning, research and handoff documents — `docs/ROADMAP.md`, `docs/GOALS.md`,
`docs/HANDOFF.md`, the literature surveys and gap analyses, the viz revamp spec —
live on disk (and in the Obsidian vault) but are **deliberately untracked**. So
are `.claude/` (agent tooling) and `figures/` (generated by
`python examples/demo.py --figures`).

If they exist in your working tree, read them: `docs/ROADMAP.md` is the plan,
`docs/HANDOFF.md` is the resume doc, `docs/GOALS.md` is the goal set (G1–G19).
If they do not, you are in a fresh clone and the code plus `docs/CORRECTIONS.md`
are self-contained.

`.gitignore` hides them and `.githooks/pre-commit` refuses to commit them. **Do
not `git add -f` them, do not `--no-verify` past the hook, and do not cite their
paths from source docstrings** — cite the goal or epic ID (`G3`, `E10`) instead,
which stays meaningful without the file. To add a genuinely public document: link
it from `README.md`, un-ignore it in `.gitignore`, and add it to `ALLOWED_DOCS` in
`.githooks/pre-commit`.

Do not add YAML frontmatter to a tracked document; vault metadata belongs in the
vault. Run `sh scripts/install-hooks.sh` once per clone to activate the hook.

## Contributing / security

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the workflow and the mathematical
correctness standard, and [`SECURITY.md`](SECURITY.md) for vulnerability
reporting. This is a solo, personal project (MIT-licensed) — **unrelated to any
employer or platform**.
