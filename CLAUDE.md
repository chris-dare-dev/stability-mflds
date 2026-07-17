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
pytest -q                           # run the suite (M2-gated tests skip without Macaulay2)
# opt in to the live Macaulay2 tests (M2 lives in WSL Debian; E10-M4):
#   BRIDGELAND_M2=$PWD/scripts/m2-wsl.cmd pytest -q      # bash
#   $env:BRIDGELAND_M2="$PWD\scripts\m2-wsl.cmd"; pytest -q   # PowerShell
python examples/demo.py             # print the corrected validated-value table
python examples/demo.py --figures   # write figures/*.{html,png}  (untracked build output)
# figure-review gallery (use `python` on Windows):
python -m http.server 8765 --directory figures/mpl_review
sh scripts/setup-clone.sh           # once per clone: hooks + the Obsidian strip filter
```

## Module map

| Module | Responsibility |
|---|---|
| `chern` | `ChernChar(r, c, e)`: H-numerical Chern character on a surface. `slope(d)`, `discriminant(d)` (CH normalization), `discriminant_brief(d)` (doubled — comparison only), `twist(s, d)`, `central_charge(s, t, d)` → `(float, float)`, `to_latex(name)`. Plus the `Q()` coerce-to-`Fraction` helper and `Number = Union[int, Fraction]`. |
| `varieties` | `Surface` / `Threefold` frozen dataclasses + catalog: `P2`, `P1xP1`, `P3`, `QUADRIC3`, `QUINTIC`, `BLOWUP_P3_POINT`; factories `K3(h2)`, `hirzebruch(n)`, `abelian_surface(h2)`, `abelian_threefold(d3)`, `fano_picard_one(name, d3)`. |
| `nslattice` | E8-M1/G12.1. `NSLattice` frozen dataclass (Gram matrix + exact `pairing` on NS⊗ℚ) with `rank1_shim(d)` / `shim_ch1(c, d)` pinning the rank-1 backward-compat encoding (`⟨ch1,H⟩ = c`, `⟨ch1,ch1⟩ = c²/d` bit-for-bit). Stdlib-only. |
| `exceptional` | Algorithm 1. `Bundle(r, c1, ch2)` with `.discriminant` / `.discriminant_brief` properties and constructors `Bundle.O(n)`, `Bundle.T_minus1()`, `Bundle.from_slope(alpha)`. Markov/ε-recursion: `enumerate_exceptional`, `exceptional_mediant`, `markov_numbers`, `is_markov_number`, `is_exceptional`, `chi` (Riemann–Roch). |
| `exceptional_surface` | E11-M1/M2 / G14. `SurfaceBundle(r, c1, ch2)` (`c1` an NS-vector in `Surface.lattice` coordinates) + the exact Riemann–Roch Euler pairing `chi(E, F, surface)`, `euler_gram`, `canonical_class`; per-family generators `p1xp1_collection`, `hirzebruch_collection(n)`, `del_pezzo_collection(deg)` and the necessary-only `is_exceptional_collection`. |
| `dlp` | Algorithm 2. DLP curve as a fractal upper envelope: `delta(mu, bundles)`, `dlp_curve(...)` → `DLPCurve`, `control_interval_halfwidth`, `moduli_nonempty`. **P² only.** |
| `dlp_hirzebruch` | E11-M6 / G18. The polarization-dependent CH Drézet–Le Potier envelope on `F_e`. `discriminant` (full-NS Δ — see invariant 2), `total_slope`, `hilbert_P`, `is_potentially_exceptional`, `is_stable_exceptional` (Cor. "DLPExceptional" rank induction), `is_semiexceptional`, `dlp_bundle` (`DLP_{H,V}`), `dlp_restricted` (`DLP_H^{<r}`, an exact max), `dlp_envelope` → `DLPEnvelope(value, exact, sharp, witness)`, `emptiness_bound` (strictly weaker than the envelope — only the branches that are theorems). Also `hirzebruch_index` — the single `F_e` dispatch choke point; it authenticates the surface *family* (Gram + `K` + `χ(O)` + kind), not just the lattice (CORRECTIONS §12 R1). |
| `reduction` | E13-M1 / G18. The Coskun–Huizenga `F_e → F_{e−2}` reduction `π` (arXiv:1907.06739 §11.1 / Lemma 11.3 — a unimodular NS isometry fixing `r`, `ch₂`, `Δ`, `χ`, transporting `A_m`): `pi_c1`, `reduce(xi, surface)`, `reduce_to_del_pezzo`. |
| `prioritary` | E13-M2 / G18. The prioritary sharp bound `δ^p_n(ν)` (Thm 1.2 = Prop 4.15 + Cor 4.17): `delta_prioritary(nu, n, surface)`, the iff-verdict `prioritary_nonempty` (domain `Δ ≥ 0`), `generic_prioritary_index`, `delta_prioritary_bundle`. |
| `generic_hn` | E13-M3b / G18. The generic `H_m`-Harder–Narasimhan filtration (arXiv:1907.06739 §5 `cor-algorithm`; rank induction + the linear-orthogonality solve): `generic_hn_factors` → the exact factor characters, `semistable_exists_hn`. Makes non-emptiness decidable on every strictly-ample `F_e`. |
| `hn_filtration` | E13-M3a/M3c / G18. The sharp-verdict layer, total on P² + every strictly-ample `F_e`: `hn_verdict(r, nu, Delta, surface)` → `HNVerdict` (regions earned from computed filtrations, factors exhibited), `semistable_exists`, `generic_hn_length`, `hn_region`. Region semantics: CORRECTIONS §12 R2, §14–15. |
| `nonemptiness_rational` | E11-M3 → E12 / G18a. The certified non-emptiness verdict layer: `moduli_nonempty` → `NonemptinessVerdict` (`VerdictStatus` + `Rigor` certificate), `validate_character`, `delta_H`, `discriminant_H` (the H-projected scalar — see invariant 2), `hirzebruch_with_polarization`, class-bound `SharpBoundEvidence`, `paper_delta_H_targets`. |
| `delta_sharp` | E14-M1/M2 / G18. The sharp Bogomolov function `δ_m^{μs}(ν)` (`def-deltass`): `mu_stable_exists(r, nu, Delta, surface)` — a PROVEN decision procedure for slope-stable existence at rational `m` (chamber-sample certificate; honest `None` only on `Δ = 1/2, r ≥ 2`) — and `delta_mu_stable(nu, m, surface, max_rank)` → `DeltaSharp(lower, upper, exact)` certified sandwich; `polarization_index` / `surface_with_index`. E14-M2: `delta_kronecker` / `kronecker_data` — the `thm-deltaKronecker` closed formula on the Kronecker triangle (exact quadratic-surd window tests). The inf need not be attained — see CORRECTIONS §17–18. |
| `stability_interval` | E14-M3 / G18. `stability_interval(r, c1, surface)` → `StabilityInterval(lo, hi, empty, witnesses)` — the exact stability interval `I_V` of an exceptional character on any `F_e` (`thm-stabilityInterval` memoized rank induction; `e ≥ 2` via the `cor-highermus` transport; both paper tables pinned). CORRECTIONS §19. |
| `walls` | Algorithm 3. `numerical_wall(v, w, d)` (exact primitive), `compute_walls` (dense numerical set — *not* certified), `walls_from_subobjects`, `actual_walls` / `actual_walls_complete` (certified-necessary, finite set). Types: `Wall`, `VerticalWall`, `ActualWall` (`.center`, `.radius_sq` exact; `.radius` is `float`). |
| `bg_check` | Algorithm 4. `check_bg_surface(...)`; re-exports `check_bg_threefold`. |
| `threefold` | Algorithm 5. `ThreefoldChern(r, a1, a2, a3)`, `bmt_Q` / `bmt_Q_at` (return `Fraction`), `alpha_crit(v, beta, d3)` → `Optional[float]`, `bg_boundary_curve(...)` → `BGBoundary`. `Threefold.bg_proven` flags proven vs conjectural. |
| `mukai` | K3 Mukai lattice. `MukaiVector(r, l, s)`, `pairing`, `self_pairing`, `moduli_dim`, `is_spherical`, `is_isotropic`, `classify_wall` → `WallClassification` (Bayer–Macrì Thm 5.7). `solve_binary_quadratic`/`saturated_basis`/`hyperbolic_witnesses` (G10 lattice solver). |
| `rigor` | Provenance lattice (E1-M2/G5), stdlib-only. `Rigor(IntEnum)` total order `PROVEN=3 > CONJECTURAL=2 > HEURISTIC=1 > UNKNOWN=0` (`min` = meet); frozen `Certificate(rigor, hypotheses: tuple, citations: tuple, note)`; `meet(*certs)` = min-rigor with order-preserving unioned hypotheses/citations; shared `UNKNOWN_CERTIFICATE` default. Lets every downstream wall/verdict advertise theorem-vs-conjecture without changing any value. |
| `oracle/` | E10 / G16 (+ the E12-M4 evidence mint). Optional Macaulay2 bridge, `oracle.m2`: `M2Session` / `find_m2` / `require_m2` (opt-in via `BRIDGELAND_M2`, graceful skip when absent), `ext_dims` / `chi_via_ext`, constructive `moduli_nonempty_by_construction` (sufficient-only), `mint_oracle_evidence`, `fe_line_bundle_cohomology` (the E10-M4 toric `F_e` cross-check). Core stays zero-dep without M2. |
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
Two objects, and they are **not** the same once ρ(X) ≥ 2. Do not conflate them.

- **The CH discriminant** — the one in the primary sources, and the one non-emptiness
  verdicts compare against: `Δ = ½·⟨ν,ν⟩ − ch₂/r`, with the *full-NS* total slope
  `ν = c₁/r`. It is **polarization-independent**; all polarization dependence lives in
  `δ_H`. See `dlp_hirzebruch.discriminant`.
- **The H-numerical scalar** carried by the `(r, c = ch₁·H, ch₂)` model that the wall /
  BG / P²-DLP machinery is built on: `Δ_H = ½·μ² − ch₂/(r·d)`, `μ = ⟨c₁,H⟩/(r·d)`,
  `d = H²`. See `ChernChar.discriminant(d)` and `nonemptiness_rational.discriminant_H`.

They satisfy `Δ = d · Δ_H` **iff c₁ ∥ H** — automatic at Picard rank 1, and on P²
(`d = 1`) they are *equal*, which is why every pinned P² value is unambiguous. For a
non-diagonal `c₁` on `P¹×P¹` / `F_n` they differ, and `Δ_H` is a lossy surrogate that
spuriously depends on `H`. Read `docs/CORRECTIONS.md` §7 before touching either.

`discriminant_brief` (= 2Δ_H) exists only for comparison with the original brief —
never make it the default. Do not change either convention.

### 3. Every mathematical change requires two-way verification
Before changing any formula, output value, or constant: (1) recompute it with
exact `fractions.Fraction` arithmetic, and (2) cite the primary source (read the
paper, not from memory). [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) is the
template. Math correctness is safety-critical here — this package exists because
it corrected a brief that got the mathematics wrong in several independent ways.

### 4. Pinned test values are ground truth
The suite pins exact values checked against the literature. If code produces a
different answer, **the code is wrong** — do not edit a test to match a new
value without completing step 3 above. A subset:
- `delta(1/2) = 5/8`, `delta(1/3) = 5/9`, `delta(1/4) = 21/32`, `delta(2/5) = 13/25`
- P²[2] wall: center `−5/2`, radius `3/2` (ABCH, arXiv:1203.0316 §9)
- Gieseker (outermost) wall of P²[n]: center `−(2n+1)/2`, radius `(2n−1)/2`
- P³ null-correlation bundle (2,0,1,0): `alpha_crit(β=1/2) = √3` (not √29/4)
- K3 structure sheaf: `v(O) = (1,0,1)`, `⟨v,v⟩ = −2` (not `(1,0,−1)`)
- Ranks 3 and 4 are **not** Markov numbers — no rank-3 or rank-4 exceptional bundle exists
- `F_0` exceptional-bundle ranks are **odd** (`(1,(0,0))`, `(3,(1,1))`, `(5,(1,2))`, …);
  `F_1` admits rank 2 (`(2,(1,1))`). Both tables reproduce Coskun–Huizenga
  arXiv:1907.06739 Tables 1–2 bit-for-bit (`tests/test_dlp_hirzebruch.py`)
- `DLP_{-K_{F_e}}(ν) ≥ 1/2` (the `F_e` analogue of the P² ½-clamp), and `δ_H = DLP_{-K}`
  is **sharp only** for `e ∈ {0,1}` with the anticanonical `H`; otherwise a lower bound

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
are `.claude/` and `AGENTS.md` (agent tooling) and `figures/` (generated by
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

### Obsidian frontmatter is stripped, not forbidden

This repo sits inside an Obsidian vault whose stamper writes `project:` / `type:` /
`authorship:` frontmatter into `CLAUDE.md`, `README.md` and `docs/*.md` (and rewrites
them CRLF), fired from a detached PostToolUse hook. **Leave it alone** — the working
tree is supposed to carry it. `.gitattributes` routes those paths through the
`obsidian-strip` clean filter (`scripts/strip-obsidian-frontmatter.py`), so git stores
the clean, LF version. Never hand-strip it, and never commit it.

Run **`sh scripts/setup-clone.sh`** once per clone: it sets `core.hooksPath` *and* the
filter, both of which are local git config that `git clone` does not carry. Forget it
and `.githooks/pre-commit` blocks the commit rather than letting metadata through.

After the stamper runs, `git status` may show these files as ` M` while `git diff` is
empty — stat churn, not content. Clear it with `git add --renormalize .`.

## Contributing / security

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the workflow and the mathematical
correctness standard, and [`SECURITY.md`](SECURITY.md) for vulnerability
reporting. This is a solo, personal project (MIT-licensed) — **unrelated to any
employer or platform**.
