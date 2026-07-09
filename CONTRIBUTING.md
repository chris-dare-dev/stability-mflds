# Contributing to bridgeland_stability

Thanks for your interest. `bridgeland_stability` is a **solo research tool**
maintained by one person in spare time. This guide explains how to contribute
in a way that is genuinely useful rather than creating maintenance burden.

**The single highest-priority rule:** every mathematical claim in this codebase
is verified two independent ways. If you touch mathematics, you must meet that
standard — see [The mathematical correctness standard](#the-mathematical-correctness-standard).

---

## What contributions are welcome

- **Mathematical corrections** — if a formula, value, or convention disagrees
  with the primary literature, please report or fix it.
  [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) shows the project's culture: every
  past correction is recorded with its exact computation and its primary source.
- **New algorithms** within scope (Drézet–Le Potier curves, Bogomolov–Gieseker
  boundaries, Bridgeland walls, K3 Mukai lattice).
- **Additional tests** that pin documented values against the literature.
- **Documentation** — clearer docstrings, more illustrative examples, corrected
  LaTeX.
- **Visualization** improvements that stay inside the `viz` subpackage and the
  `[viz]` optional extra (see [The pure-stdlib core invariant](#the-pure-stdlib-core-invariant)).
- **Bug reports** — open a GitHub issue with a minimal reproducible example.

Purely aesthetic refactors and out-of-scope features are unlikely to be merged;
if in doubt, open an issue first.

---

## Development setup

**Python ≥ 3.9** is required.

```bash
git clone https://github.com/chris-dare-dev/stability-mflds.git
cd stability-mflds
pip install -e ".[dev]"      # pytest + plotly + matplotlib
sh scripts/install-hooks.sh  # activate the pre-commit hook (see below)
```

The `[dev]` extra installs everything you need: `pytest>=7.4`, `plotly>=5.18`,
`matplotlib>=3.7`. The **core** package has **no runtime dependencies** — only
the Python standard library (`fractions`, `math`). For the core alone:

```bash
pip install -e .
```

On Windows without a POSIX shell, run
`powershell -ExecutionPolicy Bypass -File scripts\install-hooks.ps1` instead. Both
scripts just set `git config core.hooksPath .githooks`; that setting is local to a
clone and is not carried by `git clone`, which is why it must be run once.

---

## What this repository tracks

The library, its tests and examples, and a deliberately small set of documents:

- the standard project files — `README.md`, `CLAUDE.md`, `CONTRIBUTING.md`,
  `SECURITY.md`, `LICENSE`;
- **[`docs/CORRECTIONS.md`](docs/CORRECTIONS.md)** — the only tracked document
  under `docs/`, because the README links to it and the source cites it.

Everything else is a **working artifact and stays out of git**: planning and
research documents (`docs/ROADMAP.md`, `docs/GOALS.md`, `docs/HANDOFF.md`, the
literature surveys and gap analyses), agent tooling (`.claude/`), and generated
output (`figures/`, `*junit*.xml`). They may exist in a maintainer's working tree;
they are not part of the project's public surface, and they go stale fast.

`.gitignore` hides them; [`.githooks/pre-commit`](.githooks/pre-commit) refuses any
commit that contains them, as a backstop against `git add -f`. Regenerate the
figures with `python examples/demo.py --figures`.

**To add a document to the repository**, in this order:

1. link it from `README.md` — if the README does not link it, it does not belong;
2. un-ignore it in `.gitignore`;
3. add it to `ALLOWED_DOCS` in `.githooks/pre-commit`.

Never cite an untracked document's path from a source docstring or a test — cite
the goal or epic ID (`G3`, `E10`), which stays meaningful in a fresh clone. And do
not add YAML frontmatter to a tracked document: vault metadata belongs in the vault.

---

## Running the tests

```bash
pytest -q
```

All 42 tests should pass. They are organized by algorithm:

| File | Coverage |
|---|---|
| `tests/test_exceptional.py` | Exceptional bundles, Markov numbers, ε-recursion |
| `tests/test_dlp.py` | DLP curve δ(μ), moduli non-emptiness |
| `tests/test_walls.py` | Numerical and actual Bridgeland walls |
| `tests/test_threefold.py` | Threefold tilt-BG boundary α_crit(β) |
| `tests/test_mukai.py` | K3 Mukai lattice, pairing, Bayer–Macrì classification |

The demo script reproduces every validated value and (optionally) the figures:

```bash
python examples/demo.py              # prints the validated-value table
python examples/demo.py --figures    # also writes figures/*.{html,png}
```

**Every new mathematical value must be pinned by a test.**

---

## The mathematical correctness standard

> This is the most important section. Non-compliance is grounds for rejection
> regardless of code quality.

The originating project brief contained several substantive mathematical errors;
[`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) records each one with two independent
verifications. **All contributions that add or change mathematical content must
meet the same bar.**

### Two-way verification

Every mathematical claim — a formula, a computed value, a boundary condition, a
sign convention — must be verified both of these ways:

1. **Exact computation.** Compute the value with `fractions.Fraction` arithmetic
   and confirm it matches the expected result. "It looks right" or "I checked it
   numerically" is not sufficient.
2. **Primary-source citation.** Identify the theorem, formula, or example in the
   original literature and cite it precisely: author, paper (or arXiv id),
   section, theorem/example number. Read the source directly — not from memory,
   not from a secondary source, not from the project brief.

Each entry in [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) is the template: it
shows both the exact computation (e.g. `c₂ = 5/3` is not an integer) and the
primary citation (e.g. Coskun–Huizenga survey §4.2, Example 4.13).

### Pin new values with tests

If a change introduces or modifies a numerical value, assert it exactly:

```python
from fractions import Fraction as F

# Good — exact assertion, matches the primary source
assert delta(F(2, 5), bundles) == F(13, 25)   # cusp of the rank-5 bundle: 1 − 12/25

# Bad — tolerance test for a value that is exactly rational
assert abs(delta(F(2, 5), bundles) - 0.52) < 1e-6
```

Values already pinned by tests — δ(½)=5/8, δ(⅓)=5/9, the P²[2] wall
(center −5/2, radius 3/2), α_crit(β=½)=√3 for the P³ null-correlation bundle,
⟨v(O),v(O)⟩=−2, and that ranks 3 and 4 are *not* Markov numbers — must not be
changed without a documented correction and citations.

### What a math PR must include

- The exact `fractions.Fraction` computation (in a docstring, comment, or the
  test itself).
- At least one primary-source citation per new formula or value, in the codebase's
  established format: `Author, "Title", arXiv:NNNN.NNNNN, §X.Y, Theorem/Example Z`.
- A new or updated test asserting the value exactly.
- If it corrects an existing value, an entry in
  [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) following the existing format.

---

## The exact-arithmetic invariant

The core modules (`chern`, `varieties`, `exceptional`, `dlp`, `walls`,
`bg_check`, `threefold`, `mukai`) use `fractions.Fraction` for **all**
mathematically-rational values. Float is permitted **only** at explicit output
boundaries — e.g. `ChernChar.central_charge` (returns `(float, float)` for
plotting), `Wall.radius` (casts `radius_sq`), the sampling grid in
`bg_boundary_curve`, `control_interval_halfwidth`, and everything under
`bridgeland_stability/viz/`.

```python
# Correct — exact
return Fraction(1, 2) * mu * mu - Fraction(self.e, self.r * d)

# Wrong — float drift into an exact computation
return 0.5 * float(mu) ** 2 - float(self.e) / (self.r * d)
```

Rule of thumb: if a value is compared to a `Fraction` in any test, it must be
computed as a `Fraction`. If unsure which side of the boundary a value is on,
ask in the issue before writing code.

---

## The Coskun–Huizenga discriminant convention

The discriminant used throughout is the **Coskun–Huizenga normalization** (the
standard in the modern literature):

```
Δ = ½·μ² − ch₂/(r·d),   μ = (ch₁·H)/(r·d),   d = H²
```

so line bundles have Δ=0, an exceptional bundle of rank r has Δ=½(1−1/r²)<½, and
the DLP curve takes values in [½, 1]. **This convention is fixed.** The brief's
doubled convention `Δ_brief = 2Δ` is available as `ChernChar.discriminant_brief`
/ `Bundle.discriminant_brief` for comparison output only; all new formulas,
docstrings, tests, and docs must use the CH convention.

---

## The pure-stdlib core invariant

Everything reachable from `import bridgeland_stability` must stay **pure Python
standard library**. Two consequences:

1. **No new runtime dependencies in the core.** The only third-party libraries
   are declared as optional extras in `pyproject.toml`:

   ```toml
   [project.optional-dependencies]
   viz = ["plotly>=5.18", "matplotlib>=3.7"]
   dev = ["pytest>=7.4", "plotly>=5.18", "matplotlib>=3.7"]
   ```

   If a genuinely new extra is warranted, propose it in an issue first.

2. **Plotting imports must stay lazy/guarded.** matplotlib and plotly are
   imported only inside `viz` functions, via `style.require_mpl()` /
   `style.require_plotly()`. Never add a top-level `import matplotlib`/`import
   plotly` anywhere in `viz/`. `import bridgeland_stability` must succeed with
   only the standard library installed — verify this after any `viz/` change.

---

## Code style and conventions

No linter is enforced yet. The conventions visible in the codebase:

- `from __future__ import annotations` at the top of every module.
- Type hints on public functions and methods; `Fraction` inputs are typed
  `Number = Union[int, Fraction]` and coerced with `Q(x)`.
- Value objects are `@dataclass(frozen=True)` (`ChernChar`, `Bundle`, `Wall`,
  `Surface`, `Threefold`, `MukaiVector`, …). Prefer frozen dataclasses for new
  mathematical objects; don't mutate — create new instances.
- Docstrings explain the mathematical meaning, state the formula, and cite the
  source.
- Tests use `from fractions import Fraction as F` and assert exact `Fraction`
  equality, not float closeness, for rational values.

`.gitignore` references `.ruff_cache/`, so [Ruff](https://docs.astral.sh/ruff/)
has been used informally (`pip install ruff && ruff check .`). It is encouraged
but not yet a hard requirement.

---

## Submitting a change

1. **Open an issue first** for anything non-trivial — a new algorithm, a changed
   convention, or a correction to an existing value.
2. **Fork** and branch from `main`:
   ```bash
   git checkout -b fix/dlp-boundary-condition
   ```
3. **Write the code**, following the conventions above.
4. **Run the tests** before pushing — `pytest -q`; all must pass.
5. **Open a pull request** against `main`. The description should state what
   mathematical claim is added/changed, give the two-way verification (exact
   value + primary citation), and list any pinned test values that changed and
   why.
6. **Response time.** This is a solo project maintained in spare time; reviews
   may take days to weeks. Smaller, well-scoped PRs are reviewed fastest. A polite
   ping on the issue after two weeks is welcome.

---

## Reporting issues

Open a [GitHub issue](https://github.com/chris-dare-dev/stability-mflds/issues)
with:

- A minimal Python snippet that reproduces the problem.
- Your Python version (`python --version`) and package version
  (`python -c "import bridgeland_stability as bs; print(bs.__version__)"`).
- For a mathematical discrepancy: the value the code returns, the value you
  expect, and the primary source that gives the expected value.

Security vulnerabilities go through a separate private channel — see
[`SECURITY.md`](SECURITY.md).

---

## Contact

For serious mathematical questions — corrections to the literature or edge cases
not covered by existing tests — reach Chris Dare at **me@chrisdare.net**.
