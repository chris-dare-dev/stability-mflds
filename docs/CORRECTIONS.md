# Corrections to the project brief

This package implements the mathematics **as it actually is in the literature**,
which differs from the original project brief in several substantive ways. Every
correction below was verified two independent ways: by exact `fractions.Fraction`
computation, and against primary sources (read directly, not from memory). This
document records each discrepancy, the correct statement, and the citation.

The single most important convention choice: we use the **Coskun–Huizenga
normalized discriminant**. On P² and on any Picard-rank-1 surface this is the
H-numerical scalar

> Δ_H = ½·μ² − ch₂/(r·d),   d = H²,   μ = (ch₁·H)/(r·d)

(the brief used `Δ_brief = μ² − 2 ch₂/(r d) = 2·Δ_H`). The CH normalization is the
one in which every DLP / wall / BG formula in the literature is stated, so all
the explicit formulas in §§1–6 below are clean in it. `ChernChar.discriminant_brief`
returns the doubled value when needed.

> ⚠️ Once ρ(X) ≥ 2 this scalar is **not** the discriminant of the primary sources.
> The real one is the full-NS `Δ = ½⟨ν,ν⟩ − ch₂/r` with `ν = c₁/r`, and it is
> polarization-independent. They coincide on P² (`d = 1`) and agree as `Δ = d·Δ_H`
> whenever `c₁ ∥ H`. **See §7**, which records the error this conflation caused.

---

## 1. Algorithm 1 (exceptional bundles) — the brief is WRONG

**Brief:** child of consecutive exceptionals via Farey rank-addition,
`r_G = r_E + r_F`, giving ranks 1, 2, 3, 4, 5, … and a "rank-3 exceptional
bundle (3, 1, −7/6) at slope 1/3."

**Correct:** the ranks of exceptional bundles on P² are exactly the **Markov
numbers** {1, 2, 5, 13, 29, 34, 89, 169, …} (Rudakov; an exceptional *triple*
satisfies x²+y²+z²=3xyz). Their slopes are the **Markov fractions**: a slope
p/q (lowest terms) is exceptional iff q is a Markov number; then the bundle has
rank r = q, c₁ = p, ch₂ = (p²−q²+1)/(2q), Δ = ½(1 − 1/r²).

**Smoking gun:** the brief's "rank-3 bundle (3,1,−7/6)" has
`c₂ = ch₁²/2 − ch₂ = 1/2 + 7/6 = 5/3`, which is **not an integer** — so it is
not the Chern character of any vector bundle. `3` and `4` are not Markov
numbers; there is no rank-3 or rank-4 exceptional bundle on P². The note
`χ(E,E)=1` does *not* detect this, because χ depends only on (r, c₁, ch₂).

**Correct recursion** (Coskun–Huizenga survey §4.2): exceptional slopes are
`ε(ℤ[½])` where `ε(n)=n` and
`ε((2p+1)/2^{q+1}) = ε(p/2^q) · ε((p+1)/2^q)` with

> α·β = (α+β)/2 + (Δ_β − Δ_α)/(3 + α − β).

The exceptional bundle between O (slope 0) and T(−1) (slope ½) has **rank 5 at
slope 2/5**, not rank 3 at 1/3. Implemented in `exceptional.enumerate_exceptional`.

*Sources:* I. Coskun, J. Huizenga, "The birational geometry of the moduli
spaces of sheaves on P²" (Gökova survey), §4.2, Example 4.13; A. P. Veselov,
"Markov fractions and the slopes of the exceptional bundles on P²",
[arXiv:2501.06779](https://arxiv.org/abs/2501.06779); A. N. Rudakov, "The Markov
numbers and exceptional bundles on P²" (1988); Drézet–Le Potier, Ann. Sci. ENS
18 (1985).

---

## 2. Algorithm 2 (DLP curve) — the brief is WRONG (three ways)

**Brief:** ν(μ) is a piecewise parabola through consecutive exceptionals and
their Farey mediant; exceptional bundles lie **on** the curve (ν(α)=Δ_α); and
e.g. `ν(1/2)=3/4`, `ν(1/3)=ν(2/3)=8/9`.

**Correct** (Coskun–Huizenga, Thm 4.15 / Fig. 1): the DLP curve is the fractal
upper envelope

> δ(μ) = sup over exceptional slopes α with |μ−α|<3 of ( P(−|μ−α|) − Δ_α ),
> clamped below by 1/2,   where P(m) = ½(m²+3m+2).

Each exceptional bundle contributes an **upward cusp of height 1−Δ_α at μ=α**;
between cusps the curve dips to 1/2. The control interval is
`I_α = (α − x_α, α + x_α)` with `x_α = (3 − √(5+8Δ_α))/2`.

| | brief ν | **correct δ (CH)** | note |
|---|---|---|---|
| δ(0), δ(1) | 0 | **1** | line-bundle cusps |
| δ(1/2) | 3/4 | **5/8** | cusp of T(−1): 1 − 3/8 |
| δ(1/3) = δ(2/3) | 8/9 | **5/9** | controlled by O: P(−1/3); 1/3 is *not* exceptional |
| δ(1/4) | — | **21/32** | P(−1/4) |
| δ(2/5) | — | **13/25** | cusp of the rank-5 bundle: 1 − 12/25 |

The three independent errors:
1. **Exceptional bundles are isolated points strictly BELOW the curve** (at
   Δ_α < ½ ≤ δ), not on it. The brief's `ν(α)=Δ_α` is wrong (and its Test 1).
2. The local shape is a **two-branch cusp** of a single bundle's parabola, not
   one parabola through a mediant.
3. It uses the non-existent rank-3 mediant (see §1).

Implemented in `dlp.delta` / `dlp.dlp_curve`. The brief's value `ν(1/2)=3/4` is
in fact the *discriminant of the bundle T(−1)* (= Δ_brief), not the curve value.

*Source:* Coskun–Huizenga survey §4.3, Theorem 4.15, Figure 1.

---

## 3. Algorithm 3 (Bridgeland walls) — the brief's FORMULA is correct

The wall center/radius formulas are right (and match Coskun–Huizenga §5 and
ABCH). We re-derived a robust `(r,c,e)` form:

> W_rc = r c′ − r′ c,  W_re = r e′ − r′ e,  W_ce = c e′ − c′ e
> center s₀ = W_re/W_rc,  ρ² = s₀² − 2 W_ce/(d W_rc) = (s₀−μ_v)² − 2Δ_v.

**But the brief's Test 4 numbers are wrong.** For P²[2]:
* the ideal sheaf of 2 points has **ch(I_Z) = (1, 0, −2)** (the brief's
  `(1,0,−1/2)` is non-integral and `(1,0,−1)` is P²[1]);
* the unique wall is **center −5/2, radius 3/2** (destabilizer O(−1)=(1,−1,1/2)),
  not the brief's "center −1/2, radius 1/2".

*Verified verbatim against* Arcara–Bertram–Coskun–Huizenga,
[arXiv:1203.0316](https://arxiv.org/abs/1203.0316), §9; Maciocia,
[arXiv:1202.4587](https://arxiv.org/abs/1202.4587).

### Actual vs. numerical walls

`numerical_wall(v, w)` is the exact primitive (a single semicircle).
`compute_walls` enumerates *numerical* walls, of which there are densely many.
`actual_walls(v, surface)` is the certified refinement: it keeps only walls
`W(v, w)` whose destabilizer can really occur, by imposing the conditions that
are **necessary** for an actual wall (Coskun–Huizenga survey §6; Maciocia
[1202.4587](https://arxiv.org/abs/1202.4587)):

1. **rank reduction** — `0 ≤ rank(w) ≤ rank(v)` (the first destabilizing object
   has rank ≤ rank v);
2. **integral classes** — both `w` and `v−w` lie in the Chern-character lattice
   (`c₂ ∈ ℤ`), i.e. are classes of actual objects;
3. **Bogomolov on both pieces** — `Δ(w) ≥ 0` and `Δ(v−w) ≥ 0`;
4. **real semicircle** — `radius² > 0`;
5. **heart/phase ordering** — `Im Z(w) > 0` and `Im Z(v−w) > 0` on the wall, so
   `w` is a genuine sub-object in the tilted heart.

This set is **finite** (`actual_walls_complete` certifies stability under
doubled search bounds). For the Hilbert scheme P²[n] and the coprime / small-rank
cases covered by the ABCH–Coskun–Huizenga theorems it is exactly the set of
actual walls. **Validated:** P²[2] returns the single ABCH wall (center −5/2,
radius 3/2, destabilizer O(−1)); the Gieseker (outermost) wall of P²[n] is
center −(2n+1)/2, radius (2n−1)/2 for all n; the dense spurious numerical walls
(e.g. the (1,−9,34) semicircle at center −4) are correctly excluded.

---

## 4. Algorithm 4 (BG, surface) — correct

`Δ ≥ 0` for μ-semistable sheaves. Verified: T(−1) → Δ=3/8 (brief 3/4); O^⊕2 →
Δ=0 (equality); (1, 0, ¼) → Δ=−1/4 (brief −1/2), BG violated. The brief's own
gotcha #8 (O(1)⊕O(−1) is not μ-semistable, so BG does not apply to it) is
correct and respected.

---

## 5. Algorithm 5 (threefold BG boundary) — formula correct, brief's NUMBERS wrong

`Q = 4(ch₂ᵇ)² − 6 ch₁ᵇ ch₃ᵇ`, `α_crit(β) = √(max(0,Q))/|ch₁ᵇ|`. For the P³
null-correlation bundle v=(2,0,1,0), d₃=1:

| β | correct ch₃ᵇ | correct Q | correct α_crit | brief claim |
|---|---|---|---|---|
| 1/2 | −13/24 | 3 | **√3 ≈ 1.732** | √29/4 ≈ 1.34 (wrong) |
| 1 | **−4/3** | **0** | **0** | Q=2, √2/2 ≈ 0.707 (wrong) |

The brief's β=1 error is a **dropped rank factor** in the cubic term: it used
`ch₃ᵇ = −7/6` (as if r=1) instead of `−4/3` (r=2), giving the bogus Q=2. β=0 is
degenerate (ch₁ᵇ=0 → vertical wall). The quintic structure sheaf O=(1,0,0,0),
d₃=5 has Q≡0 (it sits exactly on the BG boundary).

**BG proven** (so Algorithm 5 is rigorous): P³ and all Fano 3-folds of Picard
rank 1 ([1103.5010](https://arxiv.org/abs/1103.5010),
[1510.04089](https://arxiv.org/abs/1510.04089)); abelian 3-folds
([1410.1585](https://arxiv.org/abs/1410.1585)); quintic
([1810.03434](https://arxiv.org/abs/1810.03434)). The **stronger** form FAILS on
Bl_p(P³) (Schmidt, [1602.05055](https://arxiv.org/abs/1602.05055)) — flagged by
`Threefold.bg_proven=False`.

---

## 6. K3 Mukai lattice — Test 5 and the wall types are garbled in the brief

* **v(O) = (1, 0, 1)** (since √td_K3 = (1,0,1), so the third coordinate is
  ch₂+r = 1), and ⟨v,v⟩ = −2 = −χ(O,O) with χ(O,O)=2. The brief's confused
  "(1,0,−1)" gives ⟨v,v⟩=+2, which is **not** a spherical class. Pairing:
  ⟨(r,l,s),(r′,l′,s′)⟩ = d·l·l′ − r s′ − r′ s; v² = d l² − 2rs; dim M(v)=v²+2.

  **Convention (canonical): the `ch₂ → ch₂ + r` Mukai shift is K3-only.**
  `MukaiVector.from_chern` / `mukai.classify_wall` apply it because √td(K3) =
  (1,0,1), raising radius² by exactly +2/d. For an **abelian** surface
  √td = (1,0,0), so the bare Chern triple `(r, c₁, ch₂)` *is already* the Mukai
  vector — never apply the shift there (it injects a spurious +2/d). This ledger
  is the **canonical statement** of the caveat; goals G2/G3 restate the same math.
* The brief's wall trichotomy "δ²=−2/0/2" is wrong: the only invariants are
  **spherical s²=−2** and **isotropic w²=0** (no "+2" type — that was only the
  wrong-sign artifact of (1,0,−1)). The correct four-case classification is
  Bayer–Macrì Thm 5.7 (Brill–Noether / Hilbert–Chow / Li–Gieseker–Uhlenbeck
  divisorial, then flopping, then fake), implemented in `mukai.classify_wall`.

*Source:* A. Bayer, E. Macrì, "MMP for moduli of sheaves on K3s via
wall-crossing", [arXiv:1301.6968](https://arxiv.org/abs/1301.6968), Thm 2.15
(dimension) and Thm 5.7 (classification).

---

## 7. The discriminant off P²: the H-projected scalar is NOT the CH discriminant

**Status:** corrected in E11-M6 / G18 (`bridgeland_stability/dlp_hirzebruch.py`).
This correction is against *the package's own earlier code*, not the brief.

Through E11-M5 the rational-surface non-emptiness layer compared

```
discriminant_H(xi, X) = 1/2 mu^2 - ch2/(r d),   mu = <c1,H>/(r d),  d = H^2
```

against `delta_H`. But the primary sources define the discriminant with the **full
Néron–Severi slope**, not its H-projection. Verbatim, Coskun–Huizenga
[arXiv:1907.06739](https://arxiv.org/abs/1907.06739) §2.1:

> the *total slope* ν and *discriminant* Δ of a Chern character **v** ∈ K(X) are defined by
> ν = c₁/r,  Δ = ½ν² − ch₂/r.

### Why it matters

* **Δ is polarization-independent.** `discriminant_H` is built out of μ_H, so it *moves
  with H*. Every Bogomolov-type statement ("Δ ≥ 0", "Δ ≥ δ_H(ν)") is a statement about the
  intrinsic Δ; the polarization dependence lives entirely in **δ_H**, never in Δ.
* The two agree exactly when **c₁ ∥ H**:  `Δ = d · discriminant_H`. That covers every
  Picard-rank-1 surface, and on P² (`d = 1`) they are *equal* — so **no P² value in this
  package changes**, and all pinned P² tests are untouched.
* For a non-diagonal c₁ at ρ(X) ≥ 2 they genuinely differ, and the surrogate is lossy.

### Exact-`Fraction` evidence

On `P^1 x P^1` (Gram `[[0,1],[1,0]]`, `H = f+s`, `d = 2`), for `xi = (2, f, 0)`:

| quantity | value |
|---|---|
| ν = c₁/r | `(1/2, 0)` |
| ⟨ν,ν⟩ | `0` |
| **Δ = ½⟨ν,ν⟩ − ch₂/r** | **`0`** |
| μ_H = ⟨c₁,H⟩/(r d) | `1/4` |
| `discriminant_H` | `1/32` |

`d · discriminant_H = 1/16 ≠ 0 = Δ` — the c₁ ∥ H identity fails, as it must, since
`c₁ = f` is not proportional to `H = f+s`.

The consequence is not cosmetic. With Δ = 0 the class lies *on* the Bogomolov boundary,
and the line bundle `O` (a μ_H-stable exceptional bundle) forces every μ_H-semistable
sheaf of this slope to satisfy `Δ ≥ P(−w) = 1/2` where `w = ν − ν(O) = (1/2,0)`. So
`M_H(2, f, 0)` is **provably empty**. The old code, comparing `1/32 ≥ 0`, reported
"non-empty (HEURISTIC)" — the wrong verdict, from the wrong invariant.

### What changed

* `dlp_hirzebruch.discriminant(xi, X)` is the CH discriminant and is what
  `moduli_nonempty` now compares against.
* `nonemptiness_rational.discriminant_H` is **retained**, documented as the H-projected
  scalar of the `(r, ch1·H, ch2)` model (it still agrees bit-for-bit with
  `ChernChar.discriminant(d)`), exactly as `discriminant_brief` is retained for
  comparison. It is no longer the basis of any verdict.
* The E11-M4 paper table dropped its `delta_H_paper / d` rescaling: targets are now
  stored in the paper's own normalization.
* The E11-M5 polarization-dependence witness was **rebuilt**. Its old class
  `xi = (2,(1,1),1/2)` on 𝔽₁ has `Δ = −1/8 < 0` under the true discriminant — it violates
  Bogomolov and is empty for *every* polarization; its apparent "polarization dependence"
  (`discriminant_H` = −1/36 vs 1/196) was an artifact of the surrogate. The replacement
  witness fixes `Δ = 3/8` and varies `H` so that **δ_H** moves (5/8 vs 7/8), flipping the
  verdict with both sides PROVEN. See `tests/test_nonemptiness.py::test_fn_polarization_dependence`.

*Sources:* Coskun–Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
[arXiv:1907.06739](https://arxiv.org/abs/1907.06739) §2.1 (definition), §5.4 (the
DLP surface), Cor. "deltaDLP" / "deltaDLPe" (sharpness), Cor. "K1/2" (the ½ floor),
Cor. "DLPExceptional" (the rank induction), Lemma "excFacts" (2) (the integrality
congruence), Tables 1–2 (regression data).

---

## 8. The independent E12 oracle: certified rank cutoff and the frozen corpus

**Status:** added in E12-M0 (`tests/oracle/dlp_reference.py`, `tests/oracle/corpus.py`,
`tests/test_differential.py`, `tests/test_oracle_integrity.py`). This is a **gate**, not a
package change — E12-M0 touches no file under `bridgeland_stability/`. The gate exists
because the E11 non-emptiness layer returned `Rigor.PROVEN` on false verdicts in *both*
directions on inputs the suite never exercised; a reference implementation *derived from*
that code could not have caught them. So the oracle transcribes the published theorem
statements directly and imports nothing from the package (asserted by
`test_oracle_integrity.py::test_reference_has_no_package_import`); it is exact-`Fraction`
only, with no float and no square root anywhere (`…::test_reference_uses_no_float`).

### The statement it transcribes (P²)

> **CHW 2.2** (arXiv:1401.1613 §2): a *positive-dimensional* `M(ξ)` exists iff
> `c₁ = rμ ∈ ℤ`, `χ = r(P(μ) − Δ) ∈ ℤ`, and `Δ ≥ δ(μ)`, with `Δ = ½μ² − ch₂/r`.
> Exceptional bundles are the stable `E` with `Δ(E) < ½`; their moduli space is a single
> reduced point.  **CH Ex. 1.9 / 1.14** (arXiv:1907.06739): non-exceptional μ-stable
> sheaves exist iff `Δ ≥ δ(μ)`, and a *semiexceptional* bundle is a direct sum of copies
> of an exceptional bundle.

Combined, the P² verdict the oracle computes is

```
M(ξ) ≠ ∅  ⟺  c₁ ∈ ℤ  ∧  χ ∈ ℤ  ∧  ( Δ ≥ δ(μ)  ∨  ξ = m·ch(E), E exceptional, m ≥ 1 ).
```

Exceptionality is **ε-membership** — Drézet–Le Potier Théorème A: a slope `α` carries an
exceptional bundle iff `α` is in the image of the ε-recursion `ε(n) = n`,
`ε((2p+1)/2^{q+1}) = ε(p/2^q)·ε((p+1)/2^q)` under
`α·β = (α+β)/2 + (Δ_β − Δ_α)/(3 + α − β)`. `χ(E,E) = 1` (equivalently `Δ = ½(1 − 1/r²)`)
and Markov rank are each only **necessary**; the oracle checks all three — slope in lowest
terms of denominator `r`, the exceptional `ch₂`, *and* ε-membership — which is what rejects
the impostors pinned below.

### The certified rank cutoff is a theorem, not a truncation

`reference_delta(μ)` enumerates only ε-slopes of denominator ≤ `denominator(μ)`. This bound
is **exact**, and — the point of the exercise — it is applied as an integer denominator
bound, never as the square root in which it is naturally derived.

An exceptional bundle of slope `α`, rank `ρ = denom(α)`, `Δ_α = (1 − 1/ρ²)/2`, contributes
to the DLP envelope the parabola `m ↦ P(−m) − Δ_α`, `m = |μ − α|`. It exceeds the ½ floor
exactly on `m < x_ρ`, where `x_ρ` is the smaller root of

```
P(−x) − Δ_α = ½   ⟺   x² − 3x + 1/ρ² = 0.
```

By Vieta the roots multiply to `1/ρ²` and sum to `3`; the larger root
`x₊ = (3 + √(9 − 4/ρ²))/2 > 5/2` (since `9 − 4/ρ² > 4` for every `ρ ≥ 1`), so — **without
evaluating the root** —

```
x_ρ = (1/ρ²)/x₊ < (1/ρ²)/(5/2) = 2/(5ρ²).
```

If `α ≠ μ`, write `μ = a/q` and `α = b/ρ` in lowest terms (`q = denom(μ)`); then
`|μ − α| = |aρ − bq|/(qρ) ≥ 1/(qρ)` because `aρ − bq` is a nonzero integer. For `α` to lift
`δ` above ½ at `μ` we need `|μ − α| < x_ρ < 2/(5ρ²)`, hence

```
1/(qρ) < 2/(5ρ²)  ⟹  5ρ < 2q  ⟹  ρ < 2q/5 < q.
```

The only remaining exceptional slope that can control `μ` is `α = μ` itself, of rank `q`.
So every ε-slope that raises `δ(μ)` above ½ has denominator ≤ `q`: enumerating denominators
≤ `denom(μ)` misses nothing. The package instead hard-codes `R_max = 60`, silently dropping
every cusp of rank > 60 — e.g. the rank-89 cusp at `μ = 34/89` — while carrying
`Rigor.PROVEN`. That is defect **A4** (fixed in E12-M1, where this derivation is re-verified
against a 400-random-μ sweep on the package side).

### Exact-`Fraction` evidence (the two character-decidable defects the oracle catches)

| ξ = (r, c₁, ch₂) | fact | exact value | consequence |
|---|---|---|---|
| **(610, 133, −581/2)** | `ch₂ = (133² − 610² + 1)/(2·610)` | `−581/2` | matches the exceptional `ch₂` |
| | `Δ = (1 − 1/610²)/2` | `372099/744200` | so `χ(E,E) = 1` (necessary only) |
| | ε-slopes of denominator 610 | `{233/610, 377/610}` | `133/610 ∉` ⟹ **not** exceptional |
| | `Δ < ½ ≤ δ(133/610)`, not semiexceptional | | **EMPTY** (defect **A2**) |
| **(8010, 3060, −3421)** | `μ = 3060/8010` | `34/89`  (denom 89 > 60) | |
| | `δ(34/89) = 1 − (1 − 1/89²)/2` | `3961/7921 = 356490/712890` | |
| | `Δ = ½·(34/89)² + 3421/8010` | `356489/712890` | `Δ = δ − 1/712890 < δ` |
| | the exceptional bundle of slope 34/89 is `(89, 34, −38)`; `90·(89,34,−38) = (8010,3060,−3420)` | `ch₂ = −3420 ≠ −3421` | not semiexceptional ⟹ **EMPTY** (defect **A4**) |

`610 = 2·5·61` and `89` are Markov numbers, so neither witness is caught by a "rank is not
Markov" heuristic; A2 requires the ε-recursion and A4 the denominator-`q` cutoff. The
pre-E12 package reports both **non-empty, `Rigor.PROVEN`**. Every corpus verdict was
recomputed exactly from the theorem and independently reproduced by a from-scratch
transcription of the ε-recursion that imports neither the package nor the oracle.

### The freeze contract

`test_oracle_integrity.py::FROZEN_STATUS` is a literal 14-row map
`(surface, r, c₁, ch₂) → Status`. **Appending corpus rows is free**; mutating a frozen
verdict fails `test_frozen_corpus_unchanged`, and `.githooks/pre-commit` refuses any commit
that stages a `tests/oracle/` change without a same-commit `docs/CORRECTIONS.md` entry —
this section is that entry. The intent is that a later milestone can only *strengthen* the
oracle by adding rows, never quietly relabel a verdict to match a regressing implementation.

*Sources:* Coskun–Huizenga–Woolf, "The effective cone of the moduli space of sheaves on the
plane", [arXiv:1401.1613](https://arxiv.org/abs/1401.1613) §2, Thm 2.2 and the
exceptional-bundle characterization; Coskun–Huizenga, "Existence of semistable sheaves on
Hirzebruch surfaces", [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Ex. 1.9 and
Ex. 1.14; Drézet–Le Potier, "Fibrés stables et fibrés exceptionnels sur P₂", Ann. Sci. ÉNS
**18** (1985) 193–244 (NUMDAM `ASENS_1985_4_18_2_193_0`), Théorème A (exceptional slopes =
image of the ε-recursion) and Théorème B (`Δ = (1/r)(c₂ − ((r−1)/2r)·c₁²)`, which expands to
`½μ² − ch₂/r`).

### E12-M1 (package side): the fix lands in the library

E12-M0 built the oracle and pinned A2/A4 as `xfail(strict=True)` tripwires. **E12-M1**
now repairs the package so those tripwires flip. Two edits, both `[PROVEN]`, touching only
`bridgeland_stability/exceptional.py` and `bridgeland_stability/nonemptiness_rational.py`
(no file under `tests/oracle/` is touched, so the pre-commit oracle-guard is not triggered).

**A2 — `is_exceptional` becomes ε-membership.** It was `χ(E,E)=1 ∧ c₂∈ℤ`, a merely
*necessary* condition met by infinitely many non-exceptional integral classes. It is now the
Drézet–Le Potier Théorème A test, with `χ=1 ∧ c₂∈ℤ` kept only as a cheap pre-filter:

```
is_exceptional(E):  χ(E,E)=1 ∧ c₂∈ℤ   (cheap necessary pre-filter)
                 ∧  μ.denominator == r                     (rank = reduced denominator)
                 ∧  ch₂ == Bundle.from_slope(μ).ch₂         (the unique exceptional ch₂)
                 ∧  is_exceptional_slope(μ, r)              (ε-recursion image membership)
```

`is_exceptional_slope(α, r_max=denom(α))` and `exceptional_slopes(μ, r_max)` are new public
helpers wrapping `enumerate_exceptional` (no new recursion → the same ε-image as `dlp.delta`).
The last two clauses are equivalent given the pre-filter (χ=1 ∧ denom==r ⟹ ch₂ matches
`from_slope`), but both are kept for legibility. This matches the oracle's
`reference_is_exceptional` (denom==r, exceptional ch₂, ε-membership) bit-for-bit. **Pinned
rejected** (`tests/test_exceptional.py::test_is_exceptional_rejects_epsilon_impostors`,
`tests/test_differential.py::test_A2_impostor_not_exceptional`): the ε-impostor table

```
{ 133/610, 477/610, 183/985, 802/985, 182/1325, 1143/1325 }   and   3/10
```

Every rank here IS a Markov number (610=2·5·61, 985=5·197, 1325=5²·53), so a "rank is Markov"
heuristic would still accept them; only ε-membership rejects them. **Pinned accepted**
(`::test_is_exceptional_accepts_genuine_epsilon_bundles`): `(2,1,−½)`, `(5,2,−2)`, `(1,0,0)`,
`(13,5,−11/2)`. The box divergence count `_FROZEN_A2` drops `6 → 0`; no genuine exceptional
bundle is newly rejected (the fix only *removes* impostors — `is_exceptional` becomes strictly
more conservative — so it creates no new missed-non-empty A1 divergence; `_FROZEN_A1` stays 99).

**A4 — the P² `delta_H` cutoff becomes `R_max = max(R_max, denominator(μ))`.** The hard-coded
`R_max = 60` silently dropped every DLP cusp of rank > 60. The certified-exact cutoff is
`denominator(μ)` (theorem proved above: every ε-slope that lifts `δ(μ)` over ½ has denominator
≤ `denom(μ)`). Exact-`Fraction` witness at `μ = 34/89` (denominator `89 > 60`):

| quantity | at `R_max = 60` (buggy) | at `R_max = 89` (fixed) |
|---|---|---|
| `δ(34/89)` | `½ = 356445/712890` (cusp missed) | `1 − (1 − 1/89²)/2 = 3961/7921 = 356490/712890` |
| `Δ(8010,3060,−3421)` | `356489/712890` | `356489/712890` |
| `Δ ≥ δ ?` | `356489 ≥ 356445` → **True** | `356489 ≥ 356490` → **False** |
| `moduli_nonempty(8010,3060,−3421,P2)` | **`Rigor.PROVEN`, `nonempty=True`** (wrong) | **`Rigor.PROVEN`, `nonempty=False`** (correct) |

`Δ = δ − 1/712890 < δ`, and the class is not (semi)exceptional (the only slope-34/89
exceptional bundle is `(89,34,−38)`, and `90·(89,34,−38) = (8010,3060,−3420)` has `ch₂=−3420 ≠
−3421`), so it is genuinely **EMPTY**. Pinned in
`tests/test_nonemptiness.py::test_delta_H_certified_rank_cutoff_at_rank_89`,
`tests/test_differential.py::test_A4_truncation_flips_to_empty`, and
`::test_A4_cutoff_now_captures_rank_89_cusp` (which now asserts the package's `δ` *equals* the
oracle's exact `δ`, the STRONGER corrected value — the one bug-documenting assertion this
milestone updates, justified by this entry).

**The 400-random-μ verification of the cutoff (package side).** The `§8` derivation is confirmed
empirically by `tests/test_exceptional.py::test_certified_cutoff_stable_under_margin` (seed
fixed): for 400 random `μ` of denominator up to 300,

```
delta(μ, enumerate_exceptional(μ−3, μ+3, denom(μ)))
      == delta(μ, enumerate_exceptional(μ−3, μ+3, denom(μ) + 200))    for all 400 μ (0 mismatches),
```

i.e. no ε-cusp of rank `> denom(μ)` ever contributes to `δ(μ)` — the empirical face of the
theorem that `R_max = denom(μ)` is exact. Since the box slopes all have denominator ≤ 20, the
`max(60, denom)` bump leaves every previously-pinned small-`μ` value byte-identical; only
high-denominator μ (like 34/89) change, and they change from *wrong* to *right*.

### E12-M2 (package side): the P² semiexceptional disjunct, character validation, and the dlp twin

E12-M1 closed A2/A4. **E12-M2** closes the remaining two *character-decidable* P² defects — **A1**
(the missing semiexceptional disjunct) and **A3** (no Chern-character validation) — plus **A4b**, the
same rank-truncation A4 fixed, but surviving in the P²-only twin `dlp.moduli_nonempty`. Three files
change: `bridgeland_stability/exceptional.py` (new `is_semiexceptional_p2`),
`bridgeland_stability/nonemptiness_rational.py` (new `validate_character`, `_is_p2_semiexceptional`,
`VerdictStatus`, boundary certificate), and `bridgeland_stability/dlp.py` (the twin's disjunct + cutoff).
No file under `tests/oracle/` is touched. The theorem, verbatim from the epic anchor
([arXiv:1401.1613](https://arxiv.org/abs/1401.1613) Thm 2.2 + [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Ex. 1.9/1.14):

```
M(ξ) ≠ ∅  ⟺  c₁ ∈ ℤ  ∧  χ ∈ ℤ  ∧  ( Δ ≥ δ(μ)  ∨  ξ = m·ch(E), E exceptional, m ≥ 1 ).
```

**A1 — the semiexceptional disjunct `ξ = m·ch(E)`, `m ≥ 1`.** A *semiexceptional* bundle is a direct
sum of copies of an exceptional bundle ([CH] Ex. 1.14); its moduli point (a Gieseker-polystable sheaf)
exists even though it sits **strictly below** the DLP curve. `is_semiexceptional_p2(r, c1, ch2)` mirrors
the oracle's `reference_is_semiexceptional` bit-for-bit: it divides the character by each `m | gcd(r, c₁)`
and tests the quotient with the (E12-M1-corrected, ε-membership) `is_exceptional`. Exact-`Fraction`
witnesses:

| ξ = (r, c₁, ch₂) | `Δ = ½μ² − ch₂/r` | `δ(μ)` | quotient | verdict |
|---|---|---|---|---|
| **(4, 2, −1)** = `2·ch(T(−1))` | `½·(½)² − (−1)/4 = 1/8 + 1/4 = 3/8` | `δ(1/2) = 5/8` | `m=2 → (2,1,−½) = T(−1)`, exceptional | `3/8 < 5/8`, rank 4 ≠ denom 2 (**not** a single exceptional bundle), yet `T(−1)^{⊕2}` is Gieseker-polystable ⟹ **NONEMPTY** |
| **(2, 0, 0)** = `O^{⊕2}` | `½·0 − 0/2 = 0` | `δ(0) = 1` | `m=2 → (1,0,0) = O`, exceptional | `0 < 1`, semiexceptional ⟹ **NONEMPTY** |

The package previously returned `Rigor.PROVEN`, `nonempty=False` for both — a PROVEN verdict *against*
a class that exists, the worst outcome the package can produce. `(4,2,−1)` and `(2,0,0)` are corpus rows
`("P^2",4,2,−1)=NONEMPTY`, `("P^2",2,0,0)=NONEMPTY`. **This corrects the pinned value in
`tests/test_dlp.py::test_moduli_empty_between_exceptional_and_curve`** (was wrongly
`nonempty=False`, "EMPTY" — defect A1 encoded as a passing test); it now pins `semiexceptional=True`,
`nonempty=True`. *Sources:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Ex. 1.14;
[arXiv:1401.1613](https://arxiv.org/abs/1401.1613) Thm 2.2.

**A3 — Chern-character validation (Thm 2.2 integrality).** Thm 2.2 requires `c₁ = rμ ∈ ℤ` **and**
`χ = r(P(μ) − Δ) ∈ ℤ`; a character failing either is not the Chern character of any sheaf, so `M(ξ)` is
empty. On P², `χ = r(P(μ) − Δ) = ch₂ + 3c₁/2 + r` reproduces the oracle's `_chi` identically (so
`χ ∈ ℤ ⟺ c₂ ∈ ℤ`). Witness `(1, 0, −3/2)`:

```
μ = 0,  Δ = ½·0 − (−3/2)/1 = 3/2,  χ = 1·(P(0) − 3/2) = 1·(1 − 3/2) = −1/2 ∉ ℤ   (c₂ = 0 − (−3/2) = 3/2 ∉ ℤ).
```

The package previously returned `Rigor.PROVEN`, `nonempty=True`. `validate_character(1,(0,),−3/2,P2)` now
returns `False` and `moduli_nonempty(1,(0,),−3/2,P2).nonempty` is `False` (`status = PROVEN_EMPTY`,
reason names "invalid Chern character"). Corpus row `("P^2",1,0,−3/2)=INVALID`. *Source:*
[arXiv:1401.1613](https://arxiv.org/abs/1401.1613) Thm 2.2 integrality clause.

**A3 off P² — the `c₂`-integrality clause (IMPROVE round).** The first cut checked integrality only via
the P² Euler polynomial `P` and *deferred* off-P² integrality to E12-M6's `K_H` repair. That left a hole
the adversarial stage missed: on 𝔽ₑ the native `_hirzebruch_verdict` ran **no** integrality check, so a
non-integral character reached its "`Δ > δ_H` sharp" and exceptional-bundle branches and was stamped
`PROVEN_NONEMPTY` — a forged PROVEN verdict *against* a class that is trivially empty (invariant 7's worst
outcome), the surviving P²/𝔽ₑ asymmetry being the tell. The repair needs **no** `K_H`: a coherent sheaf
has `c₂ = ½⟨c₁,c₁⟩ − ch₂ ∈ ℤ` (Chern classes are integral), which is `K_X`-independent and computed from
the NS self-pairing (`surface.lattice` — a rank-1 shim on P²). `validate_character` now enforces it on
**every** surface, and `_hirzebruch_verdict` routes through it first (short-circuiting to `PROVEN_EMPTY`
via `_INVALID_CHARACTER_CERT` before any exceptional branch). This is not a fragment of Thm 2.2: with
`c₁` integral, `c₂`-integrality **implies** χ-integrality by Riemann-Roch — `χ = ch₂ + ½⟨c₁,−K⟩ + r·χ(𝒪)`
and `c₁·(c₁−K)` is even on any surface (Wu) — so `(c₁,c₂)`-integrality is the *whole* integrality clause;
verified with **0 counterexamples** sweeping `1 ≤ r ≤ 3`, `|c₁ᵢ| ≤ 3`, integral `c₂` on P², P¹×P¹, 𝔽₁,
𝔽₂, 𝔽₃. Forge witnesses on P¹×P¹ (`= 𝔽₀`, NS Gram `[[0,1],[1,0]]`), each `PROVEN_NONEMPTY → PROVEN_EMPTY`:

| ξ = (r, c₁, ch₂) | `⟨c₁,c₁⟩` | `c₂ = ½⟨c₁,c₁⟩ − ch₂` | forged branch (before) |
|---|---|---|---|
| **(2, (0,0), −7/2)** | 0 | `0 + 7/2 = 7/2 ∉ ℤ` | `Δ = 7/4 > δ_H = 1` (sharp) |
| **(3, (1,1), −9/2)** | 2 | `1 + 9/2 = 11/2 ∉ ℤ` | exceptional-bundle disjunct |
| **(2, (2,2), −3/2)** | 8 | `4 + 3/2 = 11/2 ∉ ℤ` | `Δ = 7/4 > δ_H = 1` (sharp) |

Regression `tests/test_nonemptiness.py::test_fe_invalid_character_is_empty_not_forged_nonempty` pins all
three (plus a valid-character positive control on 𝔽₁). *Sources:*
[arXiv:1401.1613](https://arxiv.org/abs/1401.1613) Thm 2.2 (integrality clause);
[arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Sec. 2.1 (integral characters on 𝔽ₑ).

**A4b — the same rank truncation, surviving in `dlp.py`.** E12-M1 patched only
`nonemptiness_rational.delta_H`; the P²-only twin `dlp.moduli_nonempty` still enumerated exceptional
bundles at its default `R_max = 50`. The differential lens had never exercised it. Reusing the E12-M1 §8
witness at `μ = 34/89` (denominator `89 > 50`):

| quantity | exact value |
|---|---|
| `δ(34/89) = 1 − (1 − 1/89²)/2` | `3961/7921 = 356490/712890` |
| `Δ(8010,3060,−3421) = ½·(34/89)² + 3421/8010` | `356489/712890` |
| `Δ − δ` | `−1/712890 < 0` |

`dlp.moduli_nonempty` now applies `R_max = max(R_max, denominator(μ)) = 89` before
`enumerate_exceptional`, sees the rank-89 cusp, and returns `nonempty=False` — matching
`nonemptiness_rational.moduli_nonempty` and the oracle. The class is not (semi)exceptional (the only
slope-34/89 exceptional bundle is `(89,34,−38)`, and `90·(89,34,−38) = (8010,3060,−3420)` has
`ch₂ = −3420 ≠ −3421`). **`dlp.dlp_curve`'s `R_max` is deliberately left at 50**: it *draws a picture*,
and a truncated picture is honest where a truncated **decision** is not. Delivered as a direct passing
regression `tests/test_differential.py::test_A4b_dlp_truncation_flips_to_empty` (matching how E12-M1
landed A2/A4), and the differential sweep now compares the reference against **both** evaluators
(`test_box_dlp_moduli_nonempty_matches_reference`, and the strengthened
`test_dlp_moduli_nonempty_cross_check` — a hard equality now that A1/A2/A4b are all closed for the twin).

**Differential baseline.** With A1 closed on the `nonemptiness_rational` side, the box divergence count
`_FROZEN_A1` drops `99 → 0` in `tests/test_differential.py` (the sibling of E12-M1's `_FROZEN_A2 6→0`):
the package now matches the oracle on **every** integral character in `1 ≤ r ≤ 20`, `|c₁| ≤ 40`,
`0 ≤ c₂ ≤ 60`, so `test_box_status_divergences_are_exactly_A1_and_A2` sees `(n_a1, n_a2) = (0, 0)`. The
`xfail(strict=True)` tripwires `test_A1_semiexceptional_nonempty` and `test_A3_invalid_character_is_empty`
flip to passing.

**Status is branch-derived, and `≥` vs `>` off P².** `NonemptinessVerdict.status` is a computed
`VerdictStatus ∈ {PROVEN_NONEMPTY, PROVEN_EMPTY, UNKNOWN}` derived from `(nonempty, certificate.rigor)`,
never from the mode. On P² the CHW Thm 2.2 boundary `Δ = δ(μ)` is **inclusive** (`Δ ≥ δ` ⟹
`PROVEN_NONEMPTY`); off P² the CH Thm "deltaSurface"(1) needs a **strict** inequality, so an external
(PAPER/ORACLE) target with `Δ == δ_H` is downgraded to a non-PROVEN `_BOUNDARY_CERT` and reads `UNKNOWN`
(e.g. `(2,(0,0),−4)` on `P¹×P¹` with `δ_H = 2 = Δ`). This matches `_hirzebruch_verdict`'s own boundary
handling and removes the last place the shared `disc >= dH` path silently applied `≥` where the off-P²
theorem is strict. *Source:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Thm "deltaSurface"(1).

**Defect B (IMPROVE round) — the certified-target disjunct gap off P² (a certificate forger).** The
disjunct-on-every-surface principle above was applied to `_hirzebruch_verdict` (the native 𝔽ₑ path) but
**not** to the shared certified-target tail of `moduli_nonempty` (the documented `delta_H_target` +
`hn_source` entry point). There, off P², the two disjunct detectors were the P²-only
`_is_p2_exceptional` / `_is_p2_semiexceptional`, both hard-`False` off P² — so the branch collapsed to
`nonempty = valid ∧ (Δ ≥ δ_H)`. A **genuine** μ_H-stable exceptional bundle fed its **own** correct sharp
`δ_H` then returned `nonempty=False` with `Rigor.PROVEN` — a forged `PROVEN_EMPTY` against a class whose
moduli space is a single reduced point (invariant 7's worst outcome), and one *contradicting the same
function's native verdict for the identical class*. This is the off-P² analogue of the pinned P²
`test_paper_exceptional_coexists_with_target`, which is exactly why it is a defect, not an ambiguity.
Exact-`Fraction` witness — `ξ = (3, (1,1), −1)` on 𝔽₀ = `P¹×P¹` (NS Gram `[[0,1],[1,0]]`), `ν = c₁/r =
(1/3, 1/3)`:

```
⟨ν,ν⟩ = 2·(1/3)² = 2/9,   Δ = ½·(2/9) − (−1)/3 = 1/9 + 1/3 = 4/9,
Δ_V(3) = (1 − 1/3²)/2 = 4/9  ⟹  Δ = Δ_V  (a genuine rank-3 μ_H-stable exceptional bundle);
exceptional_ch2(3,(1,1)) = ⟨c₁,c₁⟩/(2r) − r·Δ_V = 2/6 − 3·(4/9) = 1/3 − 4/3 = −1 = ch₂  ✓.
sharp  δ_H = DLP_{−K}(1/3,1/3) = 5/9  (native envelope, exact and sharp),   Δ = 4/9 < 5/9.
```

So the bundle sits **strictly below** its own sharp envelope yet is NONEMPTY. `moduli_nonempty(3,(1,1),−1,
P¹×P¹)` (native) already returned `PROVEN_NONEMPTY`; the certified `PAPER`-target call with the class's own
`δ_H = 5/9` returned `PROVEN_EMPTY` — the forge. **Fix.** A surface-aware `_exceptional_disjunct(xi,
surface)` now feeds the shared tail on **every** surface: the P² detectors on P²; off P², the
surface-native `dlp_hirzebruch.is_semiexceptional`, which requires `Δ = Δ_V` and hence already subsumes the
`m = 1` pure-exceptional case with the correct `ch₂ = exceptional_ch2` — so it does **not** carry the A6
raw-`is_stable_exceptional` (an `(r, c₁)`-only test) bug; the `m = 1` sub-flag is re-derived under the same
`ch₂` guard, purely for the reason string. The boundary downgrade to `_BOUNDARY_CERT` is additionally
guarded by `∧ ¬(exceptional ∨ semiexceptional)`, so a class the disjunct proves non-empty is never
downgraded to `UNKNOWN` by the strict-inequality subtlety. The A6 impostor `(3,(1,1),0)` (`Δ = 1/9 ≠ Δ_V =
4/9`) is correspondingly **not** rescued on the target path — it stays `PROVEN_EMPTY` — which is why the
A6 fix (native path, E12-M3) is orthogonal to this one. Regression
`tests/test_nonemptiness.py::test_certified_target_off_p2_keeps_exceptional_disjunct` pins the witness, its
agreement with the native verdict, and the impostor guard. *Sources:*
[arXiv:1401.1613](https://arxiv.org/abs/1401.1613) Thm 2.2 (exceptional-bundle disjunct);
[arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Ex. 1.14 (semiexceptional = ⊕ copies of an
exceptional bundle).

**Defect B (IMPROVE round 3) — the certified-target *emptiness* threshold off P² (the same forger, other
side).** Round 2 restored the *non-emptiness* disjunct to the shared certified-target tail, but the tail
still read **emptiness** off the flat `Δ < δ_H`: off P², a certified mode with `nonempty=False` kept
`Rigor.PROVEN`, giving `PROVEN_EMPTY` for **every** non-(semi)exceptional class below the supplied `δ_H`.
That over-claims. On an ample 𝔽ₑ the converse of the CH existence theorem — "`Δ < δ_H ⟹ empty`" — is a
theorem only **below the certified `emptiness_bound`**, which is *strictly weaker* than the envelope: it
drops the `(ν − ν_V)·H = 0, ν ≠ ν_V` branch that the paper calls "somewhat arbitrary" (arXiv:1907.06739
Sec. 5.4, the counterexample before Cor. "K1/2"; see `dlp_hirzebruch.emptiness_bound` and CLAUDE.md's
"emptiness_bound is strictly weaker than the envelope" invariant). In the gap `emptiness_bound ≤ Δ < δ_H`
emptiness is **not** a theorem, and the package's **own** native `_hirzebruch_verdict` returns
`HEURISTIC`/`UNKNOWN` there — so the certified tail's `PROVEN_EMPTY` contradicted the same function's native
verdict for the identical class fed its **own** `δ_H` (invariant 7's worst outcome; unlike Defect A5 the
target is *not* forged — it equals the native `δ_H`, `sharp` and `exact` both true — so E12-M4's class-bound
evidence would not remove it). Exact-`Fraction` witness — `ξ = (3, (1,2), −1)` on 𝔽₀ = `P¹×P¹`,
`ν = c₁/r = (1/3, 2/3)`:

```
⟨ν,ν⟩ = 2·(1/3)(2/3) = 4/9,   Δ = ½·(4/9) − (−1)/3 = 2/9 + 1/3 = 5/9,
emptiness_bound = 2/9   (theorem-branch max),   sharp δ_H = DLP_{−K}(1/3,2/3) = 8/9,
2/9 ≤ 5/9 < 8/9  — strictly inside the non-theorem gap.
```

`moduli_nonempty(3,(1,2),−1, P¹×P¹)` (native) returns `UNKNOWN`; the `PAPER`-target call with the class's
own `δ_H = 8/9` returned `PROVEN_EMPTY` — the forge. **Fix.** The certified-target tail now mirrors
`_hirzebruch_verdict`'s emptiness gate: a surface-aware `_fe_emptiness_bound(xi, surface, rank_max)`
(`emptiness_bound` on an ample 𝔽ₑ, else `None`) downgrades the certificate to `_BOUNDARY_CERT` (→ `UNKNOWN`)
for the whole band `emptiness_bound ≤ Δ ≤ δ_H` when no exceptional disjunct fires. `Δ < emptiness_bound`
stays `PROVEN_EMPTY` (theorem — e.g. the A6 impostor `(3,(1,1),0)`, `Δ = 1/9 < 5/9`); `Δ > δ_H` strict stays
`PROVEN_NONEMPTY` (Thm "deltaSurface"(1)); the exceptional/semiexceptional disjunct is never downgraded. Off
𝔽ₑ (K3, abelian, nef-and-big 𝔽ₙ) there is no `emptiness_bound` theory, so only the boundary `Δ = δ_H` is
downgraded, exactly as before. Regression
`tests/test_nonemptiness.py::test_certified_target_off_p2_band_is_unknown_not_proven_empty` pins the band
witness (native ≡ target ≡ `UNKNOWN`) and both theorem-backed boundaries. *Sources:*
[arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Sec. 5.4 (the two emptiness-theorem branches) and Thm
"deltaSurface"(1) (strict `>` for existence).

**Defect C (IMPROVE round 4) — the target-LESS certified source off P² (the forger, other entry point).**
Every earlier round patched the `delta_H_target`-**supplied** tail. But `moduli_nonempty` also accepts a
certified `hn_source` (`ORACLE`/`PAPER`/`DLP`) with **no** `delta_H_target` — the documented E11-M5 hook
(*"HEURISTIC Bogomolov floor unless a certified hn_source is passed"*). Off P², that path substituted the
package's **own** native envelope value `δ_H = env.value` as the target and stamped
`_MODE_CERT[hn_source] = Rigor.PROVEN` **without ever consulting `env.sharp`**. When `H` is ample but not
anticanonical (or `e ≥ 2`) the envelope is only a *certified lower bound* (`env.sharp = False`), so
`Δ ≥ env.value` does **not** imply `Δ ≥` the (larger, uncomputed) sharp `δ_H^{μ-s}` and cannot certify
existence via Thm "deltaSurface"(1). The round-3 band downgrade covered only `[emptiness_bound, δ_H]` and
the `Δ = δ_H` boundary, leaving the **whole `Δ > env.value` region** over a lower-bound envelope stamped
`PROVEN` — so a class in the gap `[env.value, sharp δ_H)` is EMPTY yet reported `PROVEN_NONEMPTY`
(invariant 7's worst outcome), reachable through the ORACLE hook with one public call. The package's *own*
native `_hirzebruch_verdict` returns `HEURISTIC`/`UNKNOWN` for the identical class (it gates `PROVEN` on
`env.certified_sharp`), so the target-less source path **contradicted the same function's native verdict**.
A second, equivalent forge appeared where there is no envelope at all (`env is None`: K3, abelian, or a
nef-and-big non-ample 𝔽ₙ), where `δ_H` falls back to the Bogomolov floor `0` and `Δ ≥ 0` forged `PROVEN`.
Exact-`Fraction` witnesses:

```
𝔽₀, ample non-anticanonical H = 2f + s  (ray ≠ −K = (2,2), so env.sharp = False):
  ξ = (2, (−3,−2), −1),  ν = c₁/r = (−3/2, −1),  ⟨ν,ν⟩ = 2·(−3/2)(−1) = 3,
  Δ = ½·3 − (−1)/2 = 2;   env.value = 1/2  (a certified LOWER bound, env.certified_sharp = False).
  native → UNKNOWN;   pre-fix moduli_nonempty(…, hn_source=ORACLE) → PROVEN_NONEMPTY  (the forge).
K3(2) (env is None):  ξ = (2, (1), −1),  ν = (1/2),  ⟨ν,ν⟩ = ½,  Δ = ½·½ + ½ = 3/4;   δ_H = 0 (floor).
  native → UNKNOWN;   pre-fix hn_source=ORACLE → PROVEN_NONEMPTY  (the forge).
```

**Fix.** Off P² with **no** `delta_H_target`, `moduli_nonempty` now routes **every** class carrying a native
CH envelope (`env is not None`) through `_hirzebruch_verdict` — the honest native evaluator that gates
`PROVEN` on `env.certified_sharp` per-branch — regardless of `hn_source`; and forces the `HEURISTIC`
Bogomolov floor (never a certified source label) when there is no envelope (`env is None`). A bare certified
source certifies only the HN-length-one hypothesis, **not** a sharp `δ_H`; it reaches `PROVEN` off P² solely
where the native envelope is itself certified sharp (`e ∈ {0,1}`, `H` anticanonical) or where the caller
supplies the sharp `δ_H` as `delta_H_target` (the untouched first branch). The verdict for a target-less
certified source is now identical to the native one, term for term. Regression
`tests/test_nonemptiness.py::test_target_less_certified_source_off_p2_is_not_forged_proven` pins the
lower-bound-envelope and `env is None` witnesses (source ≡ native ≡ `UNKNOWN`) and checks the fix is not
over-broad: on the anticanonical del Pezzo 𝔽₀ = `P¹×P¹`, `(2,(0,0),−4)` with `Δ = 2 > 1 = δ_H` still reads
`PROVEN_NONEMPTY` with or without a bare source. *Sources:*
[arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Cor. "deltaDLP"/"deltaDLPe" (sharp only on the
anticanonical del Pezzo ray; a certified lower bound otherwise) and Thm "deltaSurface"(1) (strict `>` for
existence).

### E12-M3 (package side): the 𝔽ₑ `ch₂` guard

**Defect A6 — the native 𝔽ₑ exceptional shortcut ignores `ch₂`.** `_hirzebruch_verdict`
(`bridgeland_stability/nonemptiness_rational.py`) computed
`exceptional = is_stable_exceptional(xi.r, c1i, surface)` — a predicate of `(r, c₁)` **only**,
which never sees `ch₂` — and then short-circuited
`semiexceptional = exceptional or is_semiexceptional(xi, surface)`. So **any** class whose
`(r, c₁)` happens to carry a μ_H-stable exceptional bundle was stamped `exceptional = True`
regardless of its `ch₂`; the `exceptional or …` forced `semiexceptional = True`, and the
semiexceptional branch returned `PROVEN_NONEMPTY`. `is_semiexceptional` *does* guard on `ch₂`
(it requires `Δ = Δ_V`); the raw `exceptional or` was the entire bug.

Exact-`Fraction` witness — `ξ = (3, (1,1), 0)` on 𝔽₀ = `P¹×P¹` (NS Gram `[[0,1],[1,0]]`),
`ν = c₁/r = (1/3, 1/3)`:

```
⟨ν,ν⟩ = 2·(1/3)(1/3) = 2/9,       Δ = ½·(2/9) − 0/3 = 1/9,
Δ_V(3) = (1 − 1/3²)/2 = 4/9,      Δ = 1/9 ≠ Δ_V  ⟹  NOT an exceptional bundle;
exceptional_ch2(3,(1,1)) = ⟨c₁,c₁⟩/(2r) − r·Δ_V = 2/6 − 3·(4/9) = 1/3 − 4/3 = −1 ≠ 0 = ch₂.
emptiness_bound = 5/9   (native envelope, theorem branches only),   Δ = 1/9 < 5/9.
```

`(r, c₁) = (3, (1,1))` **does** carry a rank-3 μ_H-stable exceptional bundle (its own `ch₂`
is `−1`, so `Δ = Δ_V = 4/9`), but the character with `ch₂ = 0` is *not* it. Pre-fix,
`moduli_nonempty(3,(1,1),0, P¹×P¹)` returned `PROVEN_NONEMPTY`; the class is `PROVEN_EMPTY`
(`Δ = 1/9` sits strictly below the certified `emptiness_bound = 5/9`).

**Fix.** `exceptional` now additionally requires `xi.ch2 == exceptional_ch2(xi.r, c1i, surface)`
(exactly `Δ = Δ_V`). This is the **identical** `ch₂` guard already carried by the certified-target
twin `_exceptional_disjunct` (E12-M2, Defect B above), which is why the E12-M2 block records the
native-path A6 fix as *orthogonal to* the certified-target one. The `integral and …` short-circuit
keeps a non-integral `c₁` (`c1i = None`) out of `is_stable_exceptional`/`exceptional_ch2`, unchanged;
`xi.ch2` and `exceptional_ch2(…)` are both `Fraction`, so the `==` is exact (invariant 1). The
`semiexceptional` OR-arm is untouched: `V^{⊕m}` with `m > 1` has `gcd(r, c₁) ≠ 1`, so
`is_stable_exceptional`/`is_potentially_exceptional` is already `False` for it and `is_semiexceptional`
(the `Δ = Δ_V` sum-of-copies detector) still catches it. The impostor `(3,(1,1),0)` now falls through
to the certified emptiness branch: `Δ = 1/9 < emptiness_bound = 5/9 ⟹ PROVEN_EMPTY`.

**Emptiness here is a theorem, so `PROVEN` is honest (invariant 7).** `emptiness_bound` is
*strictly weaker* than the envelope — it keeps only the two branches of `DLP_{H,V}` that are
theorems about Gieseker-semistable sheaves (`0 < |(ν−ν_V)·H| ≤ −½K·H`, and `ν = ν_V` with
`Δ ≠ Δ_V`), dropping the "somewhat arbitrary" `(ν−ν_V)·H = 0, ν ≠ ν_V` branch
(arXiv:1907.06739 Sec. 5.4). `Δ` below that bound is a *proof* of emptiness, not an envelope
comparison, so the `Rigor.PROVEN` on `(3,(1,1),0)` is theorem-backed.

**Branch-disjointness.** A6 was a class satisfying **both** a PROVEN-empty branch predicate
(`Δ = 1/9 < emptiness_bound = 5/9`) and the (buggy) PROVEN-nonempty exceptional predicate; only
the source *order* of the branches decided the verdict. Post-guard the two families are provably
disjoint: `PROVEN_NONEMPTY` fires only via **(a)** the (semi)exceptional branch — where `Δ = Δ_V`
and `emptiness_bound` explicitly excludes the class's own `V`, so `emptiness_bound ≤ Δ_V = Δ` — or
**(b)** `env.certified_sharp ∧ Δ > δ_H`, where `emptiness_bound ≤ δ_H < Δ`. Both give
`Δ ≥ emptiness_bound`, so no character is ever reported `PROVEN_NONEMPTY` while sitting strictly
below its own certified `emptiness_bound`. `tests/test_nonemptiness.py::
test_hirzebruch_branch_disjointness_no_double_fire` pins this crisp invariant across a runtime-bounded
box on 𝔽₀/𝔽₁ containing the A6 witness and firing both families; a wider box (`r ≤ 6`, `|c₁ᵢ| ≤ 4`,
`c₂ ∈ [−4,4]`, ~8.7k integral characters) was verified offline with **zero** double-fires (the sweep
is `moduli_nonempty`-heavy — ~65 ms/character off P² because it enumerates the DLP envelope — so the
full ROADMAP index range is offline-only). Regressions:
`tests/test_differential.py::test_A6_F0_ch2_guard_empty` (the flipped strict-`xfail` tripwire) and
`tests/test_nonemptiness.py::test_A6_native_ch2_guard_is_proven_empty` (the full exact-arithmetic pin).
*Sources:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Lemma "excFacts"(1)
(`χ(v,v) = 1 ⟺ Δ = ½ − 1/(2r²)`, so `ch₂` is pinned by `(r, c₁)` — an exceptional character's `ch₂`
is not free), Cor. "DLPExceptional" (`Δ ≥ DLP_H^{<r}(ν)` certifies μ_H-stability) and Sec. 5.4 (the
two emptiness-theorem branches carried by `emptiness_bound`).

### E12-M4 (package side): class-bound sharp-bound evidence (A5)

**Defect A5 — a `(delta_H_target, hn_source)` pair forges `Rigor.PROVEN`.** `moduli_nonempty`
(`bridgeland_stability/nonemptiness_rational.py`) treated the pair as an *unverified caller assertion*: it
substituted `delta_H_target` for the sharp bound and stamped `_MODE_CERT[hn_source] = Rigor.PROVEN` without
ever checking the value against anything. Two live PROVEN-false forges survived E12-M2/M3:

```
P¹×P¹ ORACLE (the A5 tripwire):
  moduli_nonempty(2,(0,0),−4, P¹×P¹, delta_H_target=10⁶, hn_source=ORACLE)
  ν = c₁/r = (0,0),  Δ = ½·⟨ν,ν⟩ − ch₂/r = ½·0 − (−4)/2 = 2,   native sharp δ_H = 1.
  Returned a PROVEN verdict for the absurd bound 10⁶ (2 < 10⁶ ⟹ PROVEN "empty").

P² PAPER / P² ORACLE (roadmap crit. 5, verified 2026-07-10):
  moduli_nonempty(3,(0),−2, P², delta_H_target=0, hn_source=ORACLE)
  μ = 0,  Δ = ½·0² − (−2)/3 = 2/3,   native δ(0) = χ(O_{P²}) − Δ_O = 1 − 0 = 1.
  Returned PROVEN nonempty=True (2/3 ≥ 0) for a class that is natively PROVEN EMPTY
  (2/3 < δ(0) = 1, and rank 3 is not a Markov number so it is not exceptional).
```

The P² forge is the same "fixed one path, missed its twin" shape as A4b: the certified-target branch skips
the M2 `band_unknown` downgrade (`not surface.is_p2`) and so stamped `_MODE_CERT[mode] = PROVEN`
unconditionally on P².

**The gate.** An external sharp bound is honoured (→ PROVEN-eligible) **iff** the package can independently
certify a sharp bound for that class **and** the supplied value equals it, exactly (`Fraction` `==`,
invariant 1). A new `_certified_sharp_bound(xi, surface, R_max, rank_max)` returns the package's OWN
theorem-certified bound and `None` where no theorem gives one:

- **P²** — the Drézet–Le Potier closed form `δ(μ)` (always sharp; HN length one implicit).
- an **ample 𝔽ₑ with `env.certified_sharp`** — `e ∈ {0,1}`, `H` anticanonical, where CH Cor. "deltaDLP"
  gives `δ_H^{μ-s}(ν) = DLP_{−K}(ν)` and the truncation is certified exact → `env.value`.
- **nowhere else** (non-anticanonical ample 𝔽ₑ, `e ≥ 2`, K3, abelian, nef-and-big 𝔽ₙ) — `None`, and an
  external target is then refused (unverifiable ⟹ not trusted, invariant 7).

A `SharpBoundEvidence` (frozen) carries the class it was derived for and **two now-separate claims** (the
audit's crit. 2): the VALUE claim `sharp_bound` + its `sharp_bound_source`, and the sheaf-theoretic
`hn_length_one_source` ("the generic prioritary HN filtration has length one"). `moduli_nonempty` refuses
it unless `evidence.matches(r, c₁, ch₂, surface)` (crit. 1, class-bound) **and** `sharp_bound == native`.
The forgeable pair is wrapped into evidence internally, so the legacy signature is unchanged for every
honest caller. **Key property keeping the suite green:** every surviving certified-target call in the suite
already passes `target == native δ_H`, so `dH` is byte-identical to before and no accepted verdict moves —
the gate converts only the *mismatch* calls into `ValueError`.

**ORACLE is now a capability object (crit. 3).** A raw `hn_source=ORACLE` target is refused outright.
ORACLE-sourced `SharpBoundEvidence` carries a module-private `_ORACLE_TOKEN` that only
`bridgeland_stability.oracle.mint_oracle_evidence` holds; that mint runs only *after*
`moduli_nonempty_by_construction` actually returned `True` (a verified Gieseker-semistable witness), and
`SharpBoundEvidence.__post_init__` raises `TypeError` on any ORACLE object built without the token. The
`oracle → core` import stays one-directional (the core never imports `oracle`), and the token import is
deferred to the mint body, so `import bridgeland_stability` remains zero-dependency (invariant 3).

**Two pinned tests deliberately changed** (invariant 5 requires this entry). Both pinned the *same*
forgeable-override behaviour — "an absurd target overrides the native bound and flips the verdict while
staying certified" — which is exactly the A5 bug:

1. `tests/test_nonemptiness.py::test_certified_external_target_is_proven` — its second half fed
   `delta_H_target=5` (≠ native `1`) and asserted `w.nonempty is False and w.delta_H == 5`, i.e. it pinned
   the forge as a feature. Rewritten to assert `ValueError`. (The roadmap names this as the one place a
   pinned test is deliberately changed.)
2. The `bdry` probe of `tests/test_nonemptiness.py::test_verdict_status_is_branch_derived` — collateral: it
   reached the `Δ == δ_H` boundary by *forging* `delta_H_target=2` (≠ native `1`) and asserting `UNKNOWN`.
   With forged targets refused it cannot use a mismatched value to reach the boundary; it is replaced by a
   `ValueError` assertion, and the native `Δ == δ_H` band → `UNKNOWN` semantics remain pinned by
   `test_certified_target_off_p2_band_is_unknown_not_proven_empty` (which reaches the band with each class's
   *own* correct sharp bound, not a forged one).

The A5 strict-`xfail` tripwire `tests/test_differential.py::test_A5_forged_target_rejected` flips and is
expanded to all three forge paths (off-P² ORACLE, P² PAPER, P² ORACLE). New pins:
`test_sharp_bound_evidence_is_class_bound`, `test_sharp_bound_evidence_wrong_value_refused`,
`test_oracle_evidence_is_mint_guarded`, `test_raw_oracle_target_is_refused`,
`test_p2_forged_target_refused_both_paths`. **`tests/oracle/` is untouched** (invariant 6): A5 is a
certificate-provenance defect the P²-only reference oracle does not adjudicate. *Sources:*
[arXiv:1401.1613](https://arxiv.org/abs/1401.1613) §2 Thm 2.2 (existence needs `Δ ≥ δ(μ)`, not a
caller-asserted bound); [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Cor. "deltaDLP" (the sole
off-P² sharp-bound theorem, `e ∈ {0,1}`, `H` anticanonical).

### E12-M5 (package side): provenance repair (A12, A13)

M5 changes **no numeric value and flips no `xfail`** — A12/A13 are provenance defects the P²-only
reference oracle cannot adjudicate, so they carry no strict-`xfail` tripwire. But they are *wrong
citations* and a *false capability claim* reaching real verdict strings, which is exactly what this
ledger exists to correct. Every `δ_H` in `paper_delta_H_targets()` is byte-identical before and
after; only the surrounding prose (one block comment, two docstrings, one row note, one row
citation, one certificate, one enum comment, one module-reference bullet) changed.

**(i) The `paper_delta_H_targets()` fixture is regression-derived, not paper-tabulated (A12).** Since
E11-M6 (`§7` above) `moduli_nonempty` compares against the full-NS `Δ`, and every P² entry is
regressed against the package's own `dlp.delta` curve (`test_paper_p2_targets_match_native_dlp`)
while every `F₀` entry is regressed against `dlp_envelope`
(`test_paper_p1xp1_targets_match_native_envelope`). The `δ_H` values are **hand-derived from general
theorems** (the Drézet–Le Potier closed form `δ(μ)` on P²; `DLP_{−K}` on the del Pezzo `Fₑ`) and
checked against that machinery — the per-entry arXiv citation names the primary source for the
**existence verdict**, not for the numeric value. The block comment, the function docstring, and the
`PaperDeltaHTarget` class docstring are relabelled accordingly.

**(ii) The `δ(1/3) = 5/9` row note derived `5/9` from a nonexistent rank-3 exceptional bundle (A12).**
The old note read `δ(1/3) = χ(O) − Δ_{rk3 exc} = 1 − (1 − 1/3²)/2 = 1 − 4/9 = 5/9`. **There is no
rank-3 exceptional bundle**: rank 3 is not a Markov number, and `Bundle.from_slope(1/3)` has
`c₂ = 5/3 ∉ ℤ` (`tests/test_exceptional.py::test_rank3_pseudobundle_does_not_exist`). The true
controlling bundle is `O` (rank 1, slope 0). The two derivations agreed only by the numerical
coincidence `P(−1/3) = 1 − 4/9`.

Exact-`Fraction` evidence, `P(m) = (m² + 3m + 2)/2`, `Δ_α = (1 − 1/r_α²)/2`,
`δ(μ) = max(1/2, sup_α [P(−|μ−α|) − Δ_α])`:

```
P(−1/3) = ((1/9) + (−1) + 2)/2 = (1/9 + 1)/2 = (10/9)/2 = 5/9
Δ_O     = (1 − 1/1²)/2 = 0
δ(1/3)  = P(−1/3) − Δ_O = 5/9 − 0 = 5/9          ← attained at O (rank 1, slope 0)

competitors (all strictly smaller, so the sup is at O):
  α = 1/2 (r=2):  P(−1/6) − Δ_{r2} = 55/72 − 3/8 = 55/72 − 27/72 = 28/72 = 7/18 < 5/9
  α = 1   (r=1):  P(−2/3) − 0       = (4/9)/2 = 2/9                            < 5/9
  floor:          1/2                                                          < 5/9
```

**(iii) The `δ(2/5) = 13/25` row cited a misquote of CH Cor 9.13 (A12).** The `(5,(2),−2)` class is
the **genuine rank-5 slope-2/5 exceptional bundle** (rank 5 *is* Markov; `2/5` is in the image of
the Drézet–Le Potier `ε`-recursion), `Δ = Δ_E = 12/25`; its moduli space is a single reduced point,
so `M(5,(2),−2) ≠ ∅` via the bundle itself. The existence citation is now classical Drézet–Le Potier
1985 (Thm A). The old citation attributed to **CH `arXiv:1907.06739` Cor 9.13** the statement
"exceptional bundles are −K-stable on an anticanonically polarized del Pezzo". Cor 9.13 actually
states `δ^{μ-s}_{1−e/2}(ν) = DLP_{−K}(ν)` on the del Pezzo `Fₑ` (`e ∈ {0,1}`); the −K-stability of
exceptional bundles is a separate result the paper **attributes to Gorodentsev**, not a statement of
Cor 9.13.

```
Δ_{r5} = (1 − 1/5²)/2 = (24/25)/2 = 12/25
P(0)   = (0 + 0 + 2)/2 = 1
δ(2/5) = P(0) − Δ_{r5} = 1 − 12/25 = 13/25         ← the row's target, unchanged
Δ_E    = 12/25 < 13/25 = δ(2/5)                     ← the bundle sits strictly below the curve
```

**(iv) `_MODE_CERT[ORACLE]` claimed a prioritary-sheaf HN filtration no code computes (A13).** A minted
`ORACLE` verdict stamps `_MODE_CERT[HNMode.ORACLE]`, whose hypothesis string previously read
"HN-length-one datum supplied by an M2/OSCAR-constructed prioritary-sheaf HN filtration". But
`oracle/m2.py::moduli_nonempty_by_construction` constructs a **rank-1 ideal sheaf `I_Z(c₁)` on P²**
of length `l = c₁²/2 − ch₂ = c₂ ≥ 0` (torsion-free of rank 1, hence μ-stable), returns `True | None`
and **never `False`**, and handles **P² only** (no `Fₙ`). The certificate now describes that
sufficient-only witness. Its rigor stays `Rigor.PROVEN`: a construction genuinely proves
non-emptiness. The `HNMode.ORACLE` enum comment is corrected in the same change.

*Residual (honest scope).* The `ORACLE` certificate's citation tuple is left as
`("arXiv:1907.06739",)`, but the ideal-sheaf witness is really P² Riemann–Roch; this arXiv id is a
mild residual mismatch, flagged here rather than fixed to keep M5 minimal (the value it certifies is
unaffected).

**Retrieval provenance.** `arXiv:1401.1613` (Coskun–Huizenga–Woolf, "The effective cone of the
moduli space of sheaves on the plane") — the P² Thm 2.2 the module implements via
`validate_character` and `_is_p2_exceptional` — was absent from all package source. A module
`References` bullet now cites it, so all four epic-canonical ids
(`1401.1613`, `1907.06739`, `1910.14060`, `1611.02674`) are resolvable from source; a new offline
regression `tests/test_provenance.py::test_canonical_arxiv_ids_resolve` pins this and pins that the
two debunked pairings ("birational geometry" for `1611.02674`; "exceptional bundles are -K-stable"
for `1907.06739` Cor 9.13) never reappear.

**Tests (all new, none changed; `tests/oracle/` untouched, invariant 6).**
`tests/test_provenance.py`: `test_delta_third_note_not_fictitious_rank3`,
`test_rank5_citation_not_gorodentsev_misquote`, `test_paper_targets_relabelled_regression_fixture`,
`test_oracle_mode_certificate_describes_ideal_sheaf`, `test_canonical_arxiv_ids_resolve`.

*Sources:* [arXiv:1401.1613](https://arxiv.org/abs/1401.1613) §2 Thm 2.2 (the P² non-emptiness
criterion: integrality + `Δ ≥ δ(μ)` **or** exceptional); Drézet–Le Potier, *Ann. Sci. ENS* **18**
(1985) Thm A/B (existence of an exceptional bundle of slope `α` ⟺ `α` in the image of `ε`; the rank
cutoff and `Δ`-form); [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Cor 9.13
(`δ^{μ-s}_{1−e/2}(ν) = DLP_{−K}(ν)` on del Pezzo `Fₑ`; −K-stability attributed to Gorodentsev).

### E12 code-review fixes (three behaviour changes, all verified two-way)

A high-effort review of the assembled E12 diff (six finder angles + per-finding verify) found **no
soundness bug** — every A1–A6 fix reproduced correct — but three changes landed here:

1. **`canonical_class` restored on K3 / abelian.** An interim E12-M6 revision keyed on `surface.kind`
   and *raised* `NotImplementedError` for K3 / abelian, which silently broke the general-purpose
   Riemann–Roch `chi` / `euler_gram` there (`chi(O,O,K3) = χ(O_K3) = 2`, `= 0` on an abelian surface —
   both well-defined and previously correct). Since `K` is now a **stored** field (A8), the fix returns
   `surface.K` for every `canonical_order = 0` surface and raises only for torsion-canonical ones
   (Enriques / bielliptic). This *is* the real content of A11: never infer `K` from the Gram matrix. A
   K3 with `NS = U` (Gram `[[0,1],[1,0]]`, shared with `F₀`) now returns its true stored `K = (0,0)`,
   not the Gram-inferred `(−2,−2)`. The A11 test was reframed from "raises" to "returns the stored
   `(0,0)`" (`test_A11_canonical_class_returns_stored_K_not_gram_inferred`).
2. **Invalid character → `PROVEN_EMPTY` on every surface.** A character with non-integral
   `c₂ = ½⟨c₁,c₁⟩ − ch₂` is not the Chern character of any sheaf, so `M(ξ)` is empty for every
   polarization — a `K_X`-independent theorem. `_hirzebruch_verdict` already returned `PROVEN_EMPTY`,
   but the `moduli_nonempty` common tail (K3 / abelian / nef-and-big `Fₙ`) kept `cert = _MODE_CERT[mode]`
   and under-claimed it as `UNKNOWN` (`Rigor.HEURISTIC`). Never a false `PROVEN` — a conservative
   under-claim — but inconsistent; the tail now swaps in `_INVALID_CHARACTER_CERT`
   (`test_invalid_character_is_proven_empty_on_every_surface`).
3. **Single certified-rank-cutoff helper.** The A4/A4b cutoff `max(R_max, denom(μ))` was copy-pasted in
   `dlp.moduli_nonempty` and `nonemptiness_rational.delta_H`; A4b was exactly such a cutoff landing in one
   and not the other. Extracted to `exceptional.certified_rank_cutoff(μ, R_max)`, the single source both
   P² decision procedures now call, so they cannot drift apart again.

The two review follow-ups were then addressed (a second commit):

- **The ε-tree double-build is gone.** `is_exceptional_slope` — the pure ε-recursion-membership test
  that dominates a P² query's cost — is now `@lru_cache`-memoized, so `is_semiexceptional_p2`'s `m=1`
  re-check of the slope `is_exceptional` already tested is a cache hit rather than a second full
  enumeration. `Fraction(n)` hashes equal to `n`, so `int`/`Fraction` spellings of one slope share a
  cache entry. The `m ≥ 1` contract (bit-for-bit agreement with the frozen oracle's
  `reference_is_semiexceptional`) is unchanged — the alternative "start the loop at `m = 2`" was
  rejected precisely because it would have broken that oracle-agreement contract.
- **The two verdict engines were NOT merged — deliberately.** The finder read `_hirzebruch_verdict`
  and the `moduli_nonempty` common tail as duplicated regimes; on inspection their *cores* encode
  different theorems. The tail's non-emptiness signal is the **non-strict** `Δ ≥ δ(μ)` (CHW Thm 2.2 on
  P²); `_hirzebruch_verdict`'s is the **strict** `env.certified_sharp ∧ Δ > δ_H` (CH Thm "deltaSurface"
  (1)) plus an `emptiness_bound` band that exists only off P². A single "branch-derived builder" would
  have to reconcile `>` vs `≥` and the band per surface — re-introducing exactly the boundary bug the
  audit closed. So only the one genuinely-shared regime, the invalid-character verdict, was extracted
  (`_invalid_character_verdict`, called by both engines); the divergent theorems stay in their own
  engines, and both docstrings now record why a merge must not be attempted.

---

## 9. The F_e -> F_{e-2} reduction map `π` (E13-M1 / G18)

**Not a correction of the brief — a new exact structure**, recorded here under the same two-way
standard because it is math-load-bearing (it decides which envelope values are honest lower bounds
vs. sharp), and because E13-M1 appends an independent oracle reference (`reference_reduce_pi`), which
the freeze contract pairs with a `docs/CORRECTIONS.md` entry.

**The map.** Coskun–Huizenga (arXiv:1907.06739 §11.1) reduce the whole `F_e` non-emptiness problem to
the del Pezzo cases `F_0` / `F_1` by a linear map on Chern characters. In their `(E, F)` basis (`E` the
section, `F` the fiber),

> `π(r, aE + bF, d) = (r, aE′ + (b − a)F′, d)`   (their `d` is `ch₂`).

This package stores `NS(F_e)` in the basis `(f, s) = (F, E)` with Gram `G_e = [[0,1],[1,−e]]`
(`f² = 0`, `f·s = 1`, `s² = −e`). Writing `c₁ = x·f + y·s` — so `x` is the fiber (`= b`) and `y` the
section (`= a`) coefficient — the same map is the **`r`- and `ch₂`-fixing** NS map

> `π(r, (x, y), ch₂) = (r, (x − y, y), ch₂)`,   matrix `M = [[1, −1], [0, 1]]` on the column `(x, y)`.

`M ∈ SL₂(ℤ)` (`det M = 1·1 − (−1)·0 = 1`, **unimodular**), and — the key identity —

> `Mᵀ · G_{e−2} · M = [[1,0],[−1,1]]·[[0,1],[1,−(e−2)]]·[[1,−1],[0,1]] = [[0,1],[1,−e]] = G_e`,

so **`π` is an isometry `NS(F_e) → NS(F_{e−2})`** (verified exactly for `e = 2..6` in
`tests/test_reduction.py::test_pi_is_an_isometry_M_T_G_M_equals_G`). From unimodular + isometric +
`(r, ch₂)`-fixed, **every** Lemma 11.3 property follows, because each invariant is built only from
`⟨·,·⟩`, `r`, `ch₂`, `K_X`, and `χ(O_X)`:

- **(1) pairing** `⟨πu, πv⟩_{e−2} = ⟨u, v⟩_e` (isometry). E.g. `u=(3,1), v=(1,2), e=2`: both `= 3`.
- **(2) discriminant** `Δ = ½⟨ν,ν⟩ − ch₂/r`, `ν = c₁/r`: isometry + fixed `ch₂/r` ⟹ `Δ(πv) = Δ(v)`.
- **(3) `π(K)`, `π(O)`.** `K_{F_e} = (−(e+2), −2)`, `π(K_{F_e}) = (−(e+2)+2, −2) = (−e, −2) =
  K_{F_{e−2}}`; the `−K` ray `(e+2, 2) ↦ (e, 2)`. `ch(O) = (1,(0,0),0)` is fixed. `K² = 8` on every
  `F_e`, so `ch(O(K)) = (1, K, 4)` transports with `ch₂ = 4` unchanged.
- **(4) `χ(v)`, `χ(v,w)`; integral→integral, primitive→primitive.** `χ(O_X) = 1` on every `F_e`; the RR
  Euler form and `c₂ = ½⟨c₁,c₁⟩ − ch₂` are isometry/`r`/`ch₂`-built, hence preserved. `M` unimodular ⟹
  the integral lattice bijects onto itself (`c₂ ∈ ℤ` preserved) and `gcd(x−y, y) = gcd(x, y)` (primitive
  → primitive).
- **(5) polarization / Hilbert.** `A_m = −½K_{F_e} + mF ↦ A′_m = −½K_{F_{e−2}} + mF′`
  (`π(A_m) = ((e+2)/2 + m − 1, 1) = (e/2 + m, 1) = A′_m`); `μ_{A_m}` and `hilbert_P(ν)` are preserved.
- **(6) direct sums.** `π` is additive on Chern characters (`r, c₁, ch₂` all add under `⊕` and `M` is
  linear), so `π(A⊕B⊕C) = πA ⊕ πB ⊕ πC`. The paper's *named* §11.1(6) character is now pinned exactly
  (closing what was open item O1): `O(−E+(n−1)F)^A ⊕ O^B ⊕ O(−F)^C ↦ O(−E′+nF′)^A ⊕ O′^B ⊕ O(−F′)^C`,
  since in the `(f,s)` basis `O(−E+(n−1)F)` has `c₁=(n−1,−1) ↦ (n,−1)=O(−E′+nF′)`, `O(−F)=(−1,0)` and
  `O=(0,0)` are fixed, and the isometry forces the matching `ch₂`
  (`test_lemma_113_6_named_direct_sum_character`).

**Two worked characters (hand-computed, then confirmed by the package).**

| source (`F_e`, `H`, `d`) | `π` target (`F_{e−2}`, `H′`, `d`) | `Δ` | `χ` | `μ_H` | `hilbert_P(ν)` | `c₂` |
|---|---|---|---|---|---|---|
| `(2,(3,1),−1)` / F₂, `H=A₁=(3,1)`, `d=4` | `(2,(2,1),−1)` / F₀, `H′=A′₁=(2,1)`, `d=4` | `1` | `4` | `1/2` | `3` | `3` |
| `(2,(2,1),−1/2)` / F₃, `H=(4,1)`, `d=5` | `(2,(1,1),−1/2)` / F₁, `H′=(3,1)`, `d=5` | `3/8` `(=Δ_V(2))` | — | — | — | `1` |

Exact recompute of row 1: source `ν=(3/2,1/2)` on `G₂=[[0,1],[1,−2]]` gives `⟨ν,ν⟩ = 2·(3/2)(1/2) −
2·(1/2)² = 3/2 − 1/2 = 1`, so `Δ = ½·1 − (−1)/2 = 1`; `⟨c₁,c₁⟩ = 4`, `c₂ = ½·4 − (−1) = 3`; `⟨c₁,H⟩ = 4`,
`μ_H = 4/(2·4) = 1/2`; `⟨ν,K⟩ = −3`, `hilbert_P = 1 + ½(1 − (−3)) = 3`; `χ(O,·) = 4`. On the `F₀` image
`ν=(1,1/2)` on `G₀=[[0,1],[1,0]]`: `⟨ν,ν⟩ = 2·1·(1/2) = 1` — every invariant matches (the isometry). Row 2:
`Δ = 3/8 = ½ − 1/(2·2²) = Δ_V(2)` on both sides, `c₂ = 1`. **Telescope** `(3,(3,1),0)` / F₄ `H=(5,1)` →
F₂ (`H=(4,1)`) → F₀ (`H=(3,1)`): `Δ = 1/9` at both ends (`⟨(1,1/3),(1,1/3)⟩ = 2/3 − 4/9 = 2/9`,
`Δ = ½·2/9 = 1/9`).

**Honest scope — the reduced envelope is a lower bound equal to the direct one, not a sharp value.**
E11-M6's `DLP_H(ν)` is sharp (`= δ_H^{μ-s}`) only for `e ∈ {0,1}` with the **anticanonical** `H`.
Reducing a *strictly ample* `F_e` (`e ≥ 2`) can never land on that anticanonical ray: `π` is a bijection
with `π(−K_{F_e}) = −K_{F_{e−2}}`, so `π(H) ∝ −K_{F_{e−2}}` iff `H ∝ −K_{F_e}`, and `−K_{F_e}` is **not
ample** for `e ≥ 2` (Nakai: `−K = (e+2, 2)` has `a − e·b = 2 − e ≤ 0`). Consequently both envelopes carry
`sharp = False`, and by Lemma 11.3 the reduced value **equals** the direct one: the flagship slope gives
`dlp_envelope(ν, F₂).value = dlp_envelope(π(ν), F₀).value = 1`, and a genuine cusp slope
`ν = (2/3, 1/3)` / F₂ → `(1/3, 1/3)` / F₀ gives `10/9` on **both** sides (both non-sharp);
`emptiness_bound` transports identically (`1/2 = 1/2`). So the acceptance inequality
`lower_bound ≤ reduced` holds as an **exact equality**. Obtaining a *strictly sharper* `δ_H` off the
`−K` ray needs the sharp non-anticanonical theory (the prioritary bound `δ^p_n`); that is **open
question O2**, deferred to E13-M2/M3.

*Source:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) §11.1 and Lemma 11.3 (Coskun–Huizenga,
"Existence of semistable sheaves on Hirzebruch surfaces"). Package: `bridgeland_stability/reduction.py`
(`pi_c1`, `reduce`/`reduce_character`, `reduce_to_del_pezzo`, `REDUCTION_MATRIX`); tests in
`tests/test_reduction.py`; independent oracle `tests/oracle/dlp_reference.py::reference_reduce_pi`.


## 10. The prioritary sharp bound `δ^p_n(ν)` (E13-M2 / G18)

**Not a correction of the brief — a new exact structure**, recorded here under the two-way standard
because it is math-load-bearing (it is the first genuine *sharpening* of the `F_e` non-emptiness program
past the E13-M1 reduction, and it partially closes open question O2 of §9) and because E13-M2 appends an
independent oracle reference (`reference_delta_prioritary`), which the freeze contract pairs with a
`docs/CORRECTIONS.md` entry.

Notation on `F_e` (package basis `(f, s) = (F, E)`, Gram `[[0,1],[1,−e]]`): total slope
`ν = ε·E + φ·F`, so a package NS-vector `(v0, v1)` has **`ε = v1`** (the `s = E` coeff) and
**`φ = v0`** (the `f = F` coeff). Discriminant `Δ = ½⟨ν,ν⟩ − ch₂/r` — the full-NS
`dlp_hirzebruch.discriminant` (invariant 2), **never** the H-projected `discriminant_H`.

**Prioritary sheaves (Def 2.1).** A torsion-free `V` is `L`-prioritary if `Ext²(V, V(−L)) = 0` — weaker
than `μ_H`-semistability. The relevant polarizations are the fiber `F` and `H_m = E + (m+e)F`
(so `H_m·F = 1`); `P_{F,H_n}(v)` is the stack of `F`- and `H_n`-prioritary sheaves, `P_F(v)` irreducible
(Walter) and nonempty whenever `Δ ≥ 0`. Coskun–Huizenga (Thm 1.2 = Prop 4.15 + Cor 4.17) give an
explicitly computable `δ^p_n(ν)` with, for `v = (r, ν, Δ) ∈ K(F_e)`, `Δ ≥ 0`:

> **`P_{F,H_n}(v) ≠ ∅  ⟺  Δ ≥ δ^p_n(ν)`**   (Cor 4.17, an `iff` — `Rigor.PROVEN`).

**The rank-free master formula (Prop 4.15).** On the triangle `T` with vertices
`(ε,φ) = (−1, n−1), (0,0), (0,−1)`, write `(ε,φ) = λ₁(−1,n−1) + λ₂(0,0) + λ₃(0,−1)`, i.e.
`λ₁ = −ε`, `λ₃ = −((n−1)ε + φ)`, `λ₂ = 1 − λ₁ − λ₃`. The direct sum
`V = O(−E+(n−1)F)^A ⊕ O^B ⊕ O(−F)^C` with `A = rλ₁, B = rλ₂, C = rλ₃` has rank `r`, slope `ν`, and is
`F`- and `H_n`-prioritary. From `ch(V) = (A+B+C, −A·E + (A(n−1)−C)·F, ½A(−e−2(n−1)))` one computes,
using `E² = −e, E·F = 1, F² = 0`, `⟨c₁,c₁⟩ = −A²(e+2n−2) + 2AC` and `ch₂/r = −A(e+2n−2)/(2r)`, so

> `Δ(V) = ½⟨ν,ν⟩ − ch₂/r = A/(2r²)·(B(e+2n−2) + C(e+2n))`,   and — the `r`-factors **cancel** —
> **`δ^p_n(ν) = max{ ½·λ₁·( λ₂·(e+2n−2) + λ₃·(e+2n) ), 0 }`   on `T`.**

If `ε ∈ ℤ` then `δ^p_n(ν) = 0` (Def 4.11). For `n ≤ −e` the coefficient `e+2n ≤ 0`, so the bracket is
`≤ 0` everywhere and `δ^p_n ≡ 0` (Example 4.12).

**Reduction to `T` for arbitrary `ν` (Remark 4.13).** `δ^p_n` is invariant under integer *twists*
`(ε,φ) ↦ (ε+a, φ+b)` (`V ↦ V⊗O(aE+bF)`: both `Δ` and the prioritary condition are twist-invariant) and
*duals* `(ε,φ) ↦ (−ε,−φ)`. `T` and `−T` each have area `½`; up to an integer translation `{T, −T}` tile a
`ℤ²`-fundamental domain (the parallelogram spanned by `(1, 1−n)` and `(1, −n)`, `det = −1`). Algorithm:
twist `E` so `ε ∈ (−1,0)` (i.e. `ε ↦ ε − ⌈ε⌉`), twist `F` so `λ₃ ∈ [0,1)` (i.e. `φ ↦ φ + ⌊λ₃⌋`); if the
remaining `λ₂ = 1 − λ₁ − λ₃ < 0`, dualize once and re-normalize. **Proof it lands in `T`.** In `(λ₁, λ₃)`
coordinates the normalized region is the unit square `(0,1)×[0,1)`; its lower-left triangle `λ₁+λ₃ ≤ 1`
is `T`, the upper-right `λ₁+λ₃ > 1` is `T′`. The dual-plus-renormalize map is `ρ(λ₁,λ₃) = (1−λ₁, 1−λ₃)`
(for `λ₃ > 0`; `E`-twist sends `λ₁ ↦ 1−λ₁`, `F`-twist sends `λ₃ = (n−1)−λ₃_old ↦ frac = 1−λ₃`), so for
`(λ₁,λ₃) ∈ T′` the image has `λ₂ = 1 − (1−λ₁) − (1−λ₃) = λ₁+λ₃−1 > 0` — **in `T`**, and one dual suffices.
Every step is a symmetry, so the value is preserved. *Verified:* interior slope `(ε,φ) = (−½,−¼)` and its
twists (package `(−1/4,−1/2)`, `(−13/4,3/2)`, `(19/4,−3/2)`) and dual `(1/4,1/2)` all return
`δ^p_1 = 1/8` on `F_0` (`test_delta_p_twist_dual_invariant`).

**Two worked characters (hand-computed, then confirmed by the package's own `discriminant`).**

| slope `ν` (`ε,φ`) | `n`, `F_e` | `λ = (λ₁,λ₂,λ₃)` | `δ^p_n` (formula) | witness `V = (r, c₁, ch₂)` | `discriminant(V)` |
|---|---|---|---|---|---|
| `−½E − ¼F` | `1`, `F_0` (`e=0`) | `(½, ¼, ¼)` | `½·½·(¼·0 + ¼·2) = 1/8` | `(4, (−1,−2), 0)` | `1/8` |
| `−½E − ¼F` | `1`, `F_1` (`e=1`) | `(½, ¼, ¼)` | `½·½·(¼·1 + ¼·3) = 1/4` | `(4, (−1,−2), −1)` | `1/4` |

Exact recompute of row 1 (`F_0`, Gram `[[0,1],[1,0]]`): `V = O(−E)² ⊕ O ⊕ O(−F)`, `c₁ = 2·(0,−1) +
(0,0) + (−1,0) = (−1,−2)` in `(f,s)`, `ch₂ = 2·½⟨(0,−1),(0,−1)⟩ + ½⟨(−1,0),(−1,0)⟩ = 0`; `ν = (−¼,−½)`,
`⟨ν,ν⟩ = 2·(−¼)(−½) = ¼`, `Δ = ½·¼ − 0 = 1/8`. Row 2 (`F_1`, Gram `[[0,1],[1,−1]]`): same `c₁`,
`ch₂ = 2·½·⟨(0,−1),(0,−1)⟩ = ½·2·(−1) = −1`; `⟨ν,ν⟩ = 2·(−¼)(−½) − (−½)² = ¼ − ¼ = 0`,
`Δ = 0 − (−1)/4 = 1/4`. The monotone sequence on `F_1` at `ν = (−¼,−½)` is
`[δ^p_n]_{n=−2..4} = [0, 0, 0, 1/4, 1/2, 3/4, 1]` (`test_delta_p_monotonic_in_n`).

**Cor 4.18 (generic prioritary index), Example 4.9 / Figure 2 anchor.** For `ε ∉ ℤ`, with L-Gaeta parameter
`ψ = φ + ½e(⌈ε⌉−ε) − Δ/(1−(⌈ε⌉−ε))` and `L₀ = ⌈ε⌉E + ⌈ψ⌉F`,
`ρ_gen(v) = ⌊ Δ/((⌈ε⌉−ε)(ε−⌊ε⌋)) − e/2 + 1 − (⌈ψ⌉−ψ) ⌋`, and `P_{F,H_n}(v) ≠ ∅ ⟺ n ≤ ρ_gen`. The paper's
Example 4.9 / Figure 2 (`ν = ½E + ⅓F`, `Δ = 11/10`, `e = 1`): `ψ = ⅓ + ½·½ − (11/10)/(½) = ⅓ + ¼ − 11/5 = −97/60`,
`L₀ = ⌈½⌉E + ⌈−97/60⌉F = (1,−1)` (matching the Figure 2 caption `(a₀,b₀) = (1,−1)`), and
`ρ_gen = ⌊ (11/10)/(¼) − ½ + 1 − 37/60 ⌋ = ⌊257/60⌋ = 4` (`test_generic_prioritary_index_figure2`). `ψ` and the
`ρ_gen = 4` conclusion are from Example 4.9's text; the Figure 2 caption carries only `ν, Δ, e, (a₀,b₀)`.

**Form note (Prop 4.15 vs Cor 4.18).** The code writes the `ψ`-denominator as `1 − (⌈ε⌉ − ε)` (the Prop 4.15
proof form) while Cor 4.18's printed form is `ε − ⌊ε⌋`. For `ε ∉ ℤ` these are equal (`⌈ε⌉ − ⌊ε⌋ = 1`, so
`1 − (⌈ε⌉ − ε) = ε − ⌊ε⌋`), and `ε ∉ ℤ` is a precondition of both — so the two are interchangeable here.

**Honest scope — `δ^p` is the *prioritary* bound, NOT the Gieseker / semistable bound.** By the strong
Bogomolov inequality (Remark 1.4), `δ^{μ-s}_m(ν) ≥ δ^p_{⌈m⌉+1}(ν) ≥ 0`: `δ^p` sits **between** Bogomolov
(`Δ ≥ 0`) and the sharp `μ`-stable Gieseker bound `δ^{μ-s}` that E11-M6 computes as `dlp_envelope`. It is
a **lower** bound for `δ^{μ-s}`, *not itself* the semistable-sheaf existence bound — that is E13-M3, which
assembles the E13-M1 reduction, this `δ^p`, `is_stable_exceptional`, and the generic-HN filtration
(Thm 1.6 / §5) into the sharp `δ^{μ-s}` off the `−K` ray and finally retires O2. On both anticanonical del
Pezzo cases `dlp_envelope = δ^{μ-s}` exactly, and both use `n = 2` (`F_0`, `H=(1,1)=H_1`, `m=1`,
`⌈1⌉+1=2`; `F_1`, `H=(3,2)=H_{1/2}`, `m=½`, `⌈½⌉+1=2`); over a slope sweep
`dlp_envelope(ν).value ≥ δ^p_2(ν)` with **zero** violations, and on `F_1` the bound is *tight* at the
`ε = ±½` slopes (`δ^p_2 = 5/8 = dlp_envelope`, e.g. `ν = ±½E`), so the check is not vacuously met by the
`½` floor (`test_remark_1_4_vs_certified_sharp_envelope_F0/F1`).

*Source:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) §2.4, §4.1–4.3, Prop 4.15, Cor 4.17,
Cor 4.18, Remark 1.4, Remark 4.13 (Coskun–Huizenga, "Existence of semistable sheaves on Hirzebruch
surfaces"). Package: `bridgeland_stability/prioritary.py` (`delta_prioritary`, `prioritary_nonempty`,
`generic_prioritary_index`, `delta_prioritary_bundle`); tests in `tests/test_prioritary.py`; independent
oracle `tests/oracle/dlp_reference.py::reference_delta_prioritary` (purely rational — no square root,
unlike the sharp `δ^{μ-s}`).

## 11. The HN-length-one existence criterion + Thm 1.13 structure (E13-M3a / G18)

**Not a correction of the brief — a new exact interface**, recorded here under the two-way standard
because it is the first partial closure of open question **O2** of §9 (the sharp *semistable-sheaf*
bound `δ^{μ-s}` off the anticanonical ray) and because E13-M3a appends an independent oracle reference
(`reference_semistable_exists`), which the freeze contract pairs with a `docs/CORRECTIONS.md` entry.

**The load-bearing theorem (arXiv:1907.06739 §1.6, verbatim).**

> "there exists an `H_m`-semistable sheaf with Chern character `v` **if and only if the generic
> `H_m`-Harder-Narasimhan filtration has length 1**."

That HN-length-one datum is exactly what the E11-M3 numerical evaluator *delegated* ("genuinely
sheaf-theoretic, not pure Chern arithmetic"). On the decidable regions of a del Pezzo `F_e` (`e ∈ {0,1}`)
and on `P²` the generic HN length is **already determined by the shipped `moduli_nonempty` verdict**, so
E13-M3a is a thin, faithful reframing — it re-derives **no** envelope/verdict logic (which would
re-introduce the `>`-vs-`≥` boundary over-claim the E12 audit closed), it delegates and maps the
branch-derived `VerdictStatus`:

| `moduli_nonempty(...).status` | `semistable_exists` | region | generic HN length |
|---|---|---|---|
| `PROVEN_NONEMPTY` | `True`  | `S`     | `1` |
| `PROVEN_EMPTY`    | `False` | `EMPTY` | `None` (≥ 2; the exact value is M3b) |
| `UNKNOWN`         | `None`  | `K`     | `None` (pending M3b) |

Region **S** (length 1) is `Δ > δ_H^{μ-s}` (sharp) OR the character is (semi)exceptional — a non-empty
point below the envelope, the `F_e` analogue of the Drezet–Le Potier disjunct. Region **EMPTY** is the
certified obstruction `Δ < 0` (Bogomolov) or `Δ <` the certified `emptiness_bound`. The remaining band
`emptiness_bound ≤ Δ ≤ δ_H` (and the boundary `Δ = δ_H`) is the length-2 **Kronecker** region **K**,
where the sharp `δ_H^{μ-s}` *is* a Kronecker-module computation — **deferred to E13-M3b** and honestly
`None`, never a fabricated verdict (invariant 7).

**Thm 1.13 = Cor 7.7 structure (§7; Example 1.14).** For `e ∈ {0,1}`, `Δ ≥ 3/8`, and `H` sufficiently
close to `−K`: if there are no `H`-semistable sheaves then **at most one HN factor of the general
prioritary sheaf is not a semiexceptional bundle** (that one factor is the Kronecker module of region K).
The threshold `THM_1_13_MIN_DELTA = 3/8` and the shape `THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS = 1` are
pinned. **Two-way:** `3/8 = exceptional_discriminant(2) = ½ − 1/(2·2²) = ½ − 1/8 = 3/8` — the rank-2
exceptional discriminant `Δ_V = ½ − 1/(2r²)` at `r = 2`, the smallest `Δ_V` above the rank-1 floor `0`.

**The flagship K-region example (probe-confirmed, hand-recomputed).** `(2,(1,1),0)` on `F_0`
(`H=(1,1)`, `d=2`, Gram `[[0,1],[1,0]]`) **is** `O(1,0) ⊕ O(0,1)`: `ch(O(1,0)) = (1,(1,0),0)`,
`ch(O(0,1)) = (1,(0,1),0)` (each `ch₂ = ½⟨D,D⟩ = 0`), sum `(2,(1,1),0)`. Its slope `ν = (½,½)`,
`⟨ν,ν⟩ = 2·½·½ = ½`, so `Δ = ½·½ − 0 = 1/4`. The certified-sharp `dlp_envelope((½,½)) = DLP_{−K}(½,½) = 3/4`
and `emptiness_bound = 1/4`, so `Δ = 1/4 ∈ [1/4, 3/4]` — **region K**. Both summands individually exist
(each is a line bundle, region S), so the class is **genuinely non-empty**, yet M3a honestly returns
`None`: the sharp `δ_H^{μ-s}` in this length-2 band is the M3b Kronecker datum. This is the honest UNKNOWN,
the E13-M1 O2 / E13-M2 discipline applied to the Gieseker bound.

**Pinned anchors (all probe-confirmed; `Δ = ½⟨ν,ν⟩ − ch₂/r`, full-NS — invariant 2, never `discriminant_H`).**

| class `(r, c₁, ch₂)` | surface | `ν` | `Δ` | region | `exists` | why |
|---|---|---|---|---|---|---|
| `(1,(0),−5)` | `P²` | `0` | `5` | S | `True` | `5 ≥ δ(0) = 1` |
| `(2,(1),−½)` = `T(−1)` | `P²` | `½` | `3/8` | S | `True` | exceptional (below `δ(½)=5/8`) |
| `(5,(2),−2)` | `P²` | `2/5` | `12/25` | S | `True` | rank-5 exceptional (below `δ(2/5)=13/25`) |
| `(3,(0),−2)` | `P²` | `0` | `2/3` | EMPTY | `False` | `2/3 < δ(0)=1`, rank 3 not Markov |
| `(2,(0,0),−4)` | `F_0` | `(0,0)` | `2` | S | `True` | `2 > DLP_{−K}=1` |
| `(3,(1,1),−1)` | `F_0` | `(⅓,⅓)` | `4/9` | S | `True` | rank-3 `μ_H`-stable exceptional |
| `(2,(1,1),0)` | `F_0` | `(½,½)` | `1/4` | **K** | **`None`** | `emptiness_bound=1/4 ≤ Δ ≤ 3/4` (flagship) |
| `(3,(0,0),−1)` | `F_0` | `(0,0)` | `1/3` | EMPTY | `False` | `1/3 <` `emptiness_bound=1` |
| `(2,(0,0),½)` | `F_0` | `(0,0)` | `−1/4` | EMPTY | `False` | Bogomolov `Δ < 0` |
| `(2,(0,0),−4)` | `F_1` | `(0,0)` | `2` | S | `True` | `2 > DLP_{−K}=1` |
| `(2,(0,0),−2)` | `F_1` | `(0,0)` | `1` | **K** | **`None`** | boundary `Δ = DLP_{−K} = emptiness_bound = 1` (strict-inequality open question) |
| `(2,(1,1),0)` | `F_1` | `(½,½)` | `1/8` | EMPTY | `False` | `1/8 <` `emptiness_bound=5/8` |

(`F_1` = `H=(3,2)`, `d=8`, Gram `[[0,1],[1,−1]]`. Row `(2,(1,1),0)/F_1`: `⟨ν,ν⟩ = 2·½·½ − (½)² = ½ − ¼ = ¼`,
`Δ = ½·¼ = 1/8`. Reconstruction round-trip verified: from `(r, ν, Δ)` the criterion rebuilds
`c₁ = r·ν`, `ch₂ = r(½⟨ν,ν⟩ − Δ)`, e.g. `(5,(2/5),12/25) → c₁=(2), ch₂=−2`.)

**Honest scope.** M3a decides region **S** (HN length 1 → semistable sheaves exist) and the **certified-empty**
regions (`Δ < 0`; `Δ <` `emptiness_bound`) on `P²` and the ample anticanonical del Pezzo `F_0`/`F_1`. `P²` is
**total** (the DLP closed form is sharp everywhere — never `None`, no K region); `F_0`/`F_1` have a genuine
`None` (K / boundary) region, deferred to **E13-M3b** (the Kronecker-module invariants). `e ≥ 2` is out of
scope (a `NotImplementedError`) — **E13-M3c** assembles it via the E13-M1 reduction `π`. A K3 / abelian /
nef-and-big factory `F_n` carries no del Pezzo CH theory and is refused. On the anticanonical del Pezzo ray
the verdict's `sharp_bound` is bit-for-bit the certified-sharp `dlp_envelope.value` (a regression tying M3a
to the shipped sharp theory), and `semistable_exists` equals the `moduli_nonempty` status-map over a P²/F_0/F_1
grid (the no-fabrication guarantee).

**Erratum (E13 adversarial re-audit — corrected in §12).** Three claims above over-reached and are
retracted by §12 (R2): (1) the undecided band is now labelled **`UNCLASSIFIED`**, not `K` — an epistemic
UNKNOWN is not evidence of a Kronecker region, and the flagship itself refutes the `K` label (it is the
polystable `O(1,0) ⊕ O(0,1)`, so a semistable sheaf *exists* and the §1.6 criterion gives generic HN
length **one**; its `Δ = 1/4 < 3/8` is also outside Thm 1.13's stated range). `HNRegion.K` is reserved
for M3b, which will actually compute HN factors. Read the two **K** rows of the anchor table above as
**UNCLASSIFIED**. (2) "P² … no K region" holds only for the package's **existence boolean** (totality);
Example 1.14's S/K/R/empty generic-HN *shapes* occur on P² too — the region label describes the verdict,
never the sheaf's HN structure. (3) Thm 1.13 bounds the count of non-semiexceptional factors; it does
**not** assert an exactly length-two Kronecker filtration.

*Source:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) §1.6, §5 (generic HN filtration; Thm 1.6),
§7 (Thm 1.13 = Cor 7.7), Example 1.9, Example 1.14 (Coskun–Huizenga, "Existence of semistable sheaves on
Hirzebruch surfaces"). Package: `bridgeland_stability/hn_filtration.py` (`semistable_exists`,
`generic_hn_length`, `hn_region`, `hn_verdict`, `HNRegion`, `HNVerdict`, `THM_1_13_MIN_DELTA`,
`THM_1_13_MAX_NON_SEMIEXCEPTIONAL_FACTORS`); tests in `tests/test_hn_filtration.py`; independent oracle
`tests/oracle/dlp_reference.py::reference_semistable_exists` (imports nothing from the package; no float).

## 12. The E13 adversarial re-audit remediation (R1–R5)

The E13 re-audit (2026-07-13; report on disk as `docs/E12_E13_ADVERSARIAL_REAUDIT.md`, untracked) found
one theorem-level false `PROVEN_NONEMPTY` and four supporting defects. All five are fixed; each is
recorded here with two-way evidence. Suite: 465 → 478 items (13 new regressions), 6 Macaulay2 skips
unchanged.

### R1 — Gram-only `F_e` recognition minted a false proof on a K3 [P1]

`hirzebruch_index` identified `F_e` solely from the NS Gram `[[0,1],[1,-e]]`. A projective K3 with
`NS(X) ≅ U` carries the `F_0` Gram (and, rebased to `(f, s−f)`, the `F_2` Gram), so with ample
`H = (2,1)` and `ch = (5,(−3,1),−3)` the package returned `PROVEN_NONEMPTY` ("exceptional bundle").
**That verdict is false.** Exact evidence: `ν = (−3/5, 1/5)`, `⟨ν,ν⟩ = −6/25`,
`Δ = ½(−6/25) − (−3)/5 = 12/25 = Δ_exc(5)` — which is *why* the exceptional branch fired — but the Mukai
vector is `v = (r, c₁, r + ch₂) = (5, (−3,1), 2)` with `v² = ⟨c₁,c₁⟩ − 2rs = −6 − 20 = −26 < −2`, and
`c₁.H = −1` is coprime to `r = 5` (semistable ⇒ stable), so a stable K3 sheaf would need `v² ≥ −2`
(Mukai; equivalently RR + Serre duality force `ext¹ ≥ 0`). The moduli space is **empty**. The same
disguise (U in the `F_2` basis) slipped past the E13-M1 reduction `π` while breaking both Lemma 11.3
target identities: `π(K_X) = (0,0) ≠ (−2,−2) = K_{F_0}` and `χ_X(O,O) = 2 ≠ 1`.

**Fix.** `hirzebruch_index` — the single dispatch choke point for the whole CH `F_e` theory (`reduce`,
`delta_prioritary`, `hn_verdict`, `_hirzebruch_envelope`, `is_ample`, …) — now authenticates the surface
*family*, not just the lattice shape: `e ≥ 0`, `K == (−(e+2), −2)` (the Lemma 11.3(3) normalization),
`chi_O == 1`, and `kind` not in {P2, K3, abelian, enriques, bielliptic}. `(Gram, K, χ(O))` jointly pin
the deformation class: a smooth projective surface with this rank-2 NS, `K² = 8`, `χ(O) = 1` is a minimal
rational surface, hence Hirzebruch. The disguised K3 now falls back to the honest HEURISTIC Bogomolov
floor (`UNKNOWN`) in `moduli_nonempty` and is refused outright by every `F_e`-native API. Tests:
`test_dlp_hirzebruch.py::test_k3_with_hyperbolic_ns_is_not_F0` / `…_is_not_F2_either` /
`…_wrong_K_or_chi_O_is_refused_even_untagged` / `…_genuine_hirzebruch_surfaces_still_authenticate`,
`test_nonemptiness.py::test_k3_with_hyperbolic_ns_never_mints_a_proof`,
`test_reduction.py::test_reduce_refuses_a_disguised_k3`.

### R2 — M3a fabricated region K from epistemic `UNKNOWN` [P1]

`hn_verdict` mapped every underlying `UNKNOWN` to `HNRegion.K` ("length-2 Kronecker region"). The
flagship `(2,(1,1),0)` on `F_0` is itself the counterexample: it *is* `O(1,0) ⊕ O(0,1)`, the summands
share a reduced Hilbert polynomial, so the direct sum is polystable — a semistable sheaf **exists** and
§1.6 gives generic HN length **one**; and `Δ = 1/4 < 3/8` is outside Thm 1.13's range, so nothing
licensed a Kronecker claim. **Fix:** a new `HNRegion.UNCLASSIFIED` carries the undecided band; `K` is
reserved for M3b (never returned today); the P²-totality and Thm 1.13 prose were narrowed (§11 erratum).
`None` remains the honest existence answer — possibly a conservative under-claim, never a structural
claim. Tests: `test_hn_filtration.py::test_flagship_unclassified_is_a_genuinely_nonempty_undecided_class`,
`::test_P2_existence_verdict_is_total_over_grid`.

### R3 — integral-rank validation was bypassable [P1]

`validate_character` documented `r ∈ Z` but tested only `r < 1`, so `Fraction(3,2)` passed; and
`prioritary_nonempty`, `delta_prioritary_bundle`, `hn_verdict` silently `int()`-truncated it —
theorem-level answers/witnesses **for a different character than the caller supplied** (e.g.
`prioritary_nonempty(Fraction(3,2), …) == True` answered Cor 4.17 for `r = 1`). Also
`generic_prioritary_index(ν, −1, F₀) == −4` although Cor 4.18's prioritary stack lives on the Bogomolov
floor (`P_F(v) ≠ ∅ ⇔ Δ ≥ 0`, Walter — the same domain `prioritary_nonempty` already enforced for
Cor 4.17). **Fix:** `validate_character` rejects non-integral `r` (a coherent sheaf's rank is a positive
integer, so a fractional rank is the invalid-character `PROVEN_EMPTY` on every surface); the prioritary
APIs raise `ValueError("rank must be a positive integer…")`; `hn_verdict` passes a non-integral `r`
through to the invalid-character verdict instead of truncating; `generic_prioritary_index` requires
`Δ ≥ 0`. Tests: `test_nonemptiness.py::test_non_integral_rank_is_invalid`,
`test_prioritary.py::test_non_integral_rank_is_never_truncated`,
`test_hn_filtration.py::test_non_integral_rank_is_never_truncated`.

### R4 — the ORACLE evidence mint had no usable input shape [P2]

`mint_oracle_evidence` required the *scalar* P² `c1` for the construction gate (`_rank1_ideal_length`
does `Fraction(c1)`) but then executed `tuple(c1)` when minting — a `TypeError` on the success branch;
the vector spelling `(0,)` could never reach that branch (the length gate needed the scalar). No input
shape could mint, and the six `@requires_m2` skips hid it. **Fix:** a `_scalar_c1` normalizer accepts
both the scalar and the length-1 NS-vector spelling everywhere on the capability path; the minted
evidence stores `c1 = (scalar,)` — the spelling `SharpBoundEvidence.matches` compares against
`tuple(SurfaceBundle.c1)`. The repaired path is exercised **without** M2 via the canned-transcript
monkeypatch (the E10-M3 technique): mint → `moduli_nonempty(evidence=…)` → `PROVEN_NONEMPTY`
end-to-end. Tests: `test_oracle.py::test_mint_oracle_evidence_accepts_scalar_and_vector_c1`,
`::test_minted_oracle_evidence_is_honoured_end_to_end`,
`::test_mint_oracle_evidence_refuses_shapeless_c1`.

### R5 — three recorded inconsistencies

* **(a)** `_MODE_CERT[DLP]`'s hypothesis read "HN filtration has length one (implicit)…" and was stamped
  on **every** P² verdict — contradicting a `PROVEN_EMPTY` verdict (whose class has generic HN length
  ≥ 2). Reworded: the hypothesis certifies the *bound* (`delta_H = delta(mu)` is a theorem for every
  character), asserting nothing about the queried class's own HN filtration.
* **(b)** `dlp.moduli_nonempty(Bundle(1, 0, −3/2))` returned `integral=False, nonempty=False` **and**
  `positive_dimensional=True, moduli_dim=3`. A non-integral character carries no sheaves; both fields are
  now gated on integrality (`test_dlp.py::test_moduli_nonempty_non_integral_is_internally_consistent`).
* **(c)** `exceptional.py`'s module docstring still said a reduced slope `p/q` is exceptional **iff** `q`
  is a Markov number, contradicting the corrected ε-image-membership implementation (§8): a Markov
  denominator is necessary, never sufficient — e.g. `133/610` has Markov `q = 610 = 2·5·61` yet
  `(610, 133, −581/2)` is exactly the §8 impostor. The docstring now states the ε-image (Théorème A) as
  the test.

**What the re-audit confirmed sound** (recorded for scope): the P² DLP differential gate survived three
widened passes (≈23,000 characters, denominators past 11,983, ranks to 100) with zero package/reference
disagreements, and 337,000+ exact `F_e` probes found no E13-M1/M2 formula failure. The escape was in
surface-family *dispatch*, outside those gates — which is why R1's fix authenticates at the dispatch
choke point rather than adding another numeric gate.

*Sources:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) (Lemma 11.3, Prop 4.15, Cor 4.17/4.18,
§1.6, Thm 1.13 = Cor 7.7, Example 1.14); the `v² ≥ −2` bound for stable sheaves on a K3 (Mukai; see
Huybrechts, *Lectures on K3 Surfaces*, Ch. 10). Package files: `dlp_hirzebruch.py`, `hn_filtration.py`,
`nonemptiness_rational.py`, `prioritary.py`, `dlp.py`, `exceptional.py`, `oracle/m2.py`.

## 13. Retiring the ε-recursion common mode: the K-theoretic mutation oracle (E13-M4)

**The finding.** The E13 re-audit noted that the E12 differential oracle (`tests/oracle/dlp_reference.py`)
is import-independent but not *algorithmically* independent: its `_exceptional_slopes` is the **same**
ε-mediant interval subdivision production's `enumerate_exceptional` uses — same integer seeds, same
`α·β = (α+β)/2 + (Δ_β−Δ_α)/(3+α−β)`, same binary tree, same certified cutoff. A shared misunderstanding
of any of those would reproduce identically on both sides, and the widened differential passes (§12)
would stay blind to it.

**The hardening.** `tests/oracle/mutation_reference.py` generates the same finite sets from a **different
theorem, in different coordinates, with different arithmetic**:

- classes live in the numerical K-group `K(P²) ≅ ℤ³` as **integer** triples `(r, c₁, χ)` — no slope, no
  `Fraction` in the generator;
- generation is by **mutation of full exceptional collections** — `[L_A B] = χ(A,B)[A] − [B]`,
  `[R_B A] = χ(A,B)[B] − [A]` — starting from `(O, O(1), O(2)) = ((1,0,1),(1,1,3),(1,2,6))`; completeness
  is Gorodentsev–Rudakov's constructibility theorem, not the ε-image description;
- the Euler form collapses (RR, eliminating `ch₂ = χ − 3c₁/2 − r`) to the all-integer
  `χ(E,F) = r_E χ_F + r_F χ_E − r_E r_F − 3 r_F c_E − c_E c_F`, pinned against
  `χ(O(a),O(b)) = (b−a+2)(b−a+1)/2` over a grid; the twist is
  `(r, c, χ)⊗O(n) = (r, c+nr, χ+nc+(n(n+3)/2)r)`;
- `Δ_α` is **computed from each generated class's own `(r,c₁,ch₂)`** — the `(1−1/ρ²)/2` rank formula is
  never transcribed (it is forced by `χ(X,X)=1`, which the walk asserts on every class); the tests then
  verify it *emerges*;
- two live tripwires run on every visited collection: the **Markov equation** `a²+b²+c² = 3abc` on the
  rank triple (Rudakov), and `χ(X,X) = 1` per member.

**The triangulation** (`tests/test_mutation_oracle.py`). Three independent recursions must produce one
set, as full `(r, c₁, ch₂)` triples where applicable:

1. production's ε-recursion (Drezet–Le Potier Théorème A);
2. the mutation walk (K-theory);
3. **Springborn's Markov fractions** (Veselov, arXiv:2501.06779, Thm 3.1: the exceptional slopes are
   exactly the Markov fractions) — the purely number-theoretic mediant
   `p₁/q₁ ∗ p₂/q₂ = (p₁q₁+p₂q₂)/(q₁²+q₂²)` seeded on `0/1, 1/2` in `[0,½]`, transcribed inline in the
   test; slopes in `[0,1)` are `{f} ∪ {1−f}`.

Set equality is asserted at ranks 13/89/610 (and over the window `[−3,4]`); the rank multiset is pinned
against a **hardcoded OEIS A002559** Markov list; the Fibonacci branch `F₂ₖ₋₁/F₂ₖ₊₁` (2/5, 5/13, 13/34,
34/89, 89/233, 233/610) is pinned explicitly.

**The impostor family, swept.** Exact evidence worth recording: the §8 impostor `(610, 133, −581/2)`
satisfies not only `χ(E,E)=1`, integral `c₂`, and Markov rank, but even Springborn's **necessary**
congruence `p² ≡ −1 (mod q)`: `133² + 1 = 17690 = 29·610` (the genuine numerator has
`233² + 1 = 54290 = 89·610`). `610 = 2·5·61` is composite, so `p² ≡ −1` has four roots mod 610 — only
the tree structure separates slopes from impostors, which is exactly why a membership test must be
generative, never congruence-local. The sweep asserts production membership == mutation membership over
every `(q Markov ≤ 610, p)` candidate with exceptional `ch₂` and integral `c₂`, and asserts the sweep
rejected at least one candidate (no vacuous gate).

**Cross-cutoff δ differential.** Production evaluates `δ(μ)` with the certified sharp cutoff
`rank ≤ q = denominator(μ)` (§8); the mutation oracle deliberately uses `4q + 64`. Their agreement over a
dense `q ≤ 32` sweep plus high-`q` spot checks (including `233/610` and `355/113`) tests the sharp-cutoff
theorem itself: a binding bundle of rank in `(q, 4q+64]` would now surface as a mismatch. The full
verdict runs as a **triple differential** (production == ε-reference == mutation-reference) on the E12
frozen corpus and the audit box `r ≤ 6, c₁ ∈ [−8,8], c₂ ∈ [0,7]`.

**What remains shared, honestly.** All three sides still transcribe (a) the Euler polynomial
`P(m) = (m²+3m+2)/2` and (b) the CHW Thm 2.2 verdict statement — both are one-line paper transcriptions
anchored by pinned literature values (`δ(1/2)=5/8`, …, and the binomial pins), which is the appropriate
mitigation for formula-level risk. The **F_e** envelope machinery (congruence-enumerated, not
ε-recursive) has no second-generator analogue yet; its residual hardening candidate remains a
CAS/Macaulay2 cross-check (E10/G16 infrastructure).

*Sources:* Gorodentsev–Rudakov, "Exceptional vector bundles on projective spaces", Duke Math. J. 54
(1987), 115–130 (mutations; constructibility on P²). Rudakov, "Exceptional vector bundles on P² and
Markov numbers", Izv. Akad. Nauk SSSR Ser. Mat. 52 (1988); Engl. transl. Math. USSR-Izv. 32 (1989),
99–112 (Markov rank triples; (rank, slope) determines the bundle).
[arXiv:2501.06779](https://arxiv.org/abs/2501.06779) (Veselov, "Markov fractions and the slopes of the
exceptional bundles on P²", after B. Springborn — the third recursion and the `p²+1 ≡ 0 (mod q)`
necessary congruence). [arXiv:1401.1613](https://arxiv.org/abs/1401.1613) Thm 2.2 and
[arXiv:1907.06739](https://arxiv.org/abs/1907.06739) Ex. 1.9/1.14 (the verdict statements, as in §8).
Files: `tests/oracle/mutation_reference.py`, `tests/test_mutation_oracle.py`,
`tests/test_oracle_integrity.py`.

## 14. The generic Harder–Narasimhan filtration: total verdicts on the del Pezzo scope (E13-M3b)

**What shipped.** `bridgeland_stability/generic_hn.py` implements arXiv:1907.06739 §5's finite inductive
procedure (`thm-HNcriterion` / `cor-algorithm`) computing the characters of the generic
`H_m`-Harder–Narasimhan factors of a character `v` on any strictly-ample `F_e` (`m = a/b − e ∈ ℚ₊`,
`H_m = H/b`). Because "there exists an `H_m`-semistable sheaf iff the generic HN filtration has
length 1" (§1.6), this makes non-emptiness **decidable** — including E13-M3a's honest UNKNOWN band and
the boundary `Δ = δ_H`. The `hn_verdict` layer is now **total** on P²/ample-F₀/F₁: a PROVEN
`moduli_nonempty` status binds; an UNKNOWN status is decided by the computed filtration, with the factor
characters **exhibited** in `HNVerdict.factors` and region `K` **earned** (length ≥ 2 with a
non-semiexceptional factor), never asserted — the §12 R2 discipline preserved.

**The algorithm** (all citations are `\label`s in the paper's source, verified against the fetched
LaTeX). Standing hypothesis: the `F`-and-`H_⌈m⌉`-prioritary stack of `v` is nonempty; if it is empty,
`M_{H_m}(v)` is empty by `prop-ssPrior` (semistable ⇒ `H_n`-prioritary for every integer `n < m+2`).
Search `w₁ = gr₁` over the `lem-HNclose`/`lem-slopeQuad` windows
(`0 ≤ (ν₁−ν)·H_m ≤ 1`, `|(ν₁−ν)·F| ≤ max{1, 2/(e+2m)}`, closed supersets of the strict bounds);
`w₁ = gr₁` iff (1) `u = v − w₁` carries `H_⌈m⌉`-prioritary sheaves — decided by `cor-prioritaryDelta`
(= the shipped E13-M2 `prioritary_nonempty`) for `Δ(u) ≥ 0` — with `u`'s factors known by rank
induction, (2) `q₁ > q₂` (reduced `H_m`-Hilbert polynomials, compared exactly as the lexicographic key
`(ν·H_m, P(ν) − Δ)` — the `t²`- and common linear coefficients cancel), (3) `μ(w₁) − μ(w_k) ≤ 1`,
(4) `χ(w₁, wᵢ) = 0` for all of `u`'s factors, (5) `M_{H_m}(w₁)` nonempty (induction). No valid `w₁`
⟺ length 1 ⟺ nonempty. Rank-1 base: integral with `c₂ ≥ 0`.

**Two derivations recorded for the domain edge and the speed:**

* **`Δ < 0` refuses condition (1) for every `n ≥ 1`, all slopes.** `thm-prioritaryNecessary` forces
  `χ(v(−L₀−H_n)) ≤ 0`. For `ε ∉ ℤ`, `cor-equivalentInequality` reads it as
  `n ≤ Δ/((⌈ε⌉−ε)(ε−⌊ε⌋)) − e/2 + 1 − (⌈ψ⌉−ψ) ≤ 4Δ + 1 < 1` when `Δ < 0` (the denominator is `≤ ¼`).
  For `ε ∈ ℤ`, `def-L0` gives `⌈ε⌉−ε = 0` and the same RR expansion collapses to `χ/r = −Δ`
  (`rem-epInteger`), so the necessary inequality is exactly `Δ ≥ 0`.
* **The linear-orthogonality solve.** Condition (4) and bilinearity force `χ(w₁, u) = Σχ(w₁, wᵢ) = 0`
  — computable with NO recursion, and **linear in ch₂(w₁)** with coefficient `2(r_u − r₁)`. For
  `r₁ ≠ r_u` it pins ch₂(w₁) to at most one lattice value per `(r₁, c₁)`; the balanced case scans a
  window bounded by two further proven prunes: `Δ₁ ≤ P(ν_j − ν₁) ≤ C(e, m)` (orthogonality + Bogomolov
  over the doubled slope window) and the `lem-discBound` break (`Δ₁` is the minimal semistable
  discriminant at its slope, so the ascending scan stops at the first semistable hit). Effect measured
  on the rank-15 paper pin: 225 s → 0.14 s.

**The paper's two orthogonal-Kronecker pins, reproduced bit-for-bit** (§1.5; package `(f,s)` coords:
a paper slope `εE + φF` of rank `r` is `c₁ = (rφ, rε)`; sums and χ-orthogonality re-verified with the
package's general RR pairing in `tests/test_generic_hn.py`, at TWO distinct ε each):

| surface, `m` | `v` (package) | computed factors | paper Δ's |
|---|---|---|---|
| `F₁`, `12/7+ε` | `(13, (6,3), −13/2)` | `(2,(0,1),−3/2)`, `(11,(6,2),−5)` | `5/8`, `65/121` |
| `F₀`, `25/9+ε` | `(15, (5,3), −8)` | `(2,(−1,1),−2)`, `(13,(6,2),−6)` | `3/4`, `90/169` |

Neither factor is semiexceptional (asserted) — the Kronecker-pair shape the paper constructs to show
`δ^{μs}_m > DLP^{<r}_{H_m}` can be strict. Through `hn_verdict` the F₁ pin yields
`exists=False, generic_hn_length=2, region=K, factors=…` with a PROVEN certificate.

**The flagship flips, the boundary closes.** `(2,(1,1),0)` on `F₀` — §12 R2's counterexample —
now **decides to `exists=True` with computed length 1**, matching the polystable truth (the paper's
own example after `cor-delPezzoKss`: `O(F₁)⊕O(F₂)` is `−K`-semistable with `Δ = ¼ < DLP₋ₖ = ¾`).
The F₁ boundary anchor `(2,(0,0),−2)` (`Δ = 1 = δ_H`, the E11-M6 strict-inequality open question O2)
**decides to existence**: no `w₁` passes the iff — hand-checked rejections: `w₁ = O` fails (4) with
`χ(O, I_{Z₂}) = 1 − 2 = −1 ≠ 0`; `w₁ = (1,(0,0),−1)` fails `q₁ > q_v` (both slope 0 with
`P(0) − Δ = 0`). O2 is closed on the del Pezzo scope: the boundary is decided per-class by the
computed filtration.

**Consistency gates** (`tests/test_generic_hn.py`): the algorithm never contradicts the envelope
verdicts over integral-`c₂` grids on F₀/F₁ (two independent theorem routes — §5 vs the DLP envelope —
with the UNKNOWN band genuinely exercised); existence implies the Cor 4.17 prioritary bound; the
`cor-algorithm` uniqueness of `gr₁` is asserted under a full-sweep flag (`PARANOID_UNIQUENESS`) over a
grid; the inlined integer `χ` and `q`-key forms are cross-pinned against the package's general
`exceptional_surface.chi` / `hilbert_P` / `discriminant`; the module accepts any ample `F_e` (an `e=2`
smoke runs now; E13-M3c will differential the direct computation against the E13-M1 reduction `π`).

Suite: 498 → 516 items (18 new), 6 Macaulay2 skips unchanged, 0 failures.

*Source:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) §1.5, §1.6, §4
(`thm-prioritaryNecessary`, `def-L0`, `cor-equivalentInequality`, `rem-epInteger`, `prop-triangle`,
`cor-prioritaryDelta`), §5 (`lem-HNclose`, `lem-HNorthogonal`, `thm-HNcriterion`, `lem-slopeQuad`,
`lem-discBound`, the ℓ ≤ 4 lemma, `cor-algorithm`), §7 (the `O(F₁)⊕O(F₂)` example after
`cor-delPezzoKss`). Package: `bridgeland_stability/generic_hn.py`, `hn_filtration.py`;
tests in `tests/test_generic_hn.py`, `tests/test_hn_filtration.py`.

## 15. `e ≥ 2` unlocked: π-equivariance of the generic HN filtration (E13-M3c)

**What shipped.** The verdict layer (`hn_filtration._require_del_pezzo_scope`) now admits **every
strictly ample `F_e`**: the §14 architecture is uniform in `e` (the envelope's PROVEN branches —
Bogomolov, the exceptional disjunct, `emptiness_bound` — hold on every ample `F_e`; off the del Pezzo
base the envelope is only a certified lower bound, so more classes fall to the computed filtration).
The gate is a **π-equivariance differential** against the E13-M1 reduction.

**The transport arithmetic (re-derived, then asserted over grids).** Lemma 11.3 transports every
ingredient of `thm-HNcriterion` exactly:

- pairing/χ/Δ (isometry, `r` and `ch₂` fixed), `K` and `χ(O)` (family invariants);
- the polarization: `H = (a, b) ↦ π(H) = (a−b, b)` shifts the H-index `m = a/b − e ↦ m + 1`
  (equivalently `A_k ↦ A_k` with `H_n = A_{n−1+e/2}`), so the algorithm's prioritary index
  `⌈m⌉ ↦ ⌈m⌉ + 1` transports **automatically** with condition (1) — and `cor-equivalentInequality`'s
  bound is equivariant because `ψ' − ψ = −⌈ε⌉ ∈ ℤ` leaves the fractional part `⌈ψ⌉ − ψ` unchanged while
  `−e/2 ↦ −(e−2)/2` absorbs the index shift;
- the `lem-slopeQuad` width `e + 2m` is invariant under `(e, m) ↦ (e−2, m+1)`; the reduced-Hilbert
  q-keys `(ν·H_m, P(ν) − Δ)` are isometry-invariant (`π(K_e) = K_{e−2}`).

Hence the computed factor lists must transport **bit-for-bit**:
`factors(π(v), π(surface)) == π(factors(v, surface))`.

**Gates** (`tests/test_generic_hn.py`):

- π-equivariance over integral-`c₂` grids `F₂ → F₀` and `F₃ → F₁` (180 classes each, multi-factor
  cases required to appear) plus a `reduce_to_del_pezzo` telescoping spot-check `F₄ → F₂ → F₀`;
- verdict totality on `F₂`: `exists ∈ {True, False}`, PROVEN certificates, UNCLASSIFIED never fires;
- the envelope-consistency differential of §14 extended to `F₂`;
- **a new pinned `e ≥ 2` Kronecker datum**: the §14 F₁ paper pin transported UP by `π⁻¹(x,y) = (x+y, y)`
  — `v = (13, (9,3), −13/2)` on `F₃` with `H = (261, 70)` (the lift of `(191, 70)`).  The envelope is
  UNKNOWN there (never certified sharp off the del Pezzo base), and the verdict layer **earns region K**
  with computed length 2 and factors `(2,(1,1),−3/2)`, `(11,(8,2),−5)` — exactly `π⁻¹` of the paper's
  factors, as Lemma 11.3 demands.

One honest scope note: region `K` is earned only where the *verdict* ran the filtration (the UNKNOWN
band).  A class the envelope already proves empty (below `emptiness_bound`) reports region `EMPTY` with
`factors=None` — call `generic_hn_factors` directly for its filtration.  Suite: 516 → 522 items,
6 Macaulay2 skips unchanged, 0 failures.

*Source:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) §11 (Lemma 11.3, the reduction `π`),
§4–5 as in §14.  Package: `bridgeland_stability/hn_filtration.py` (scope), `generic_hn.py` (unchanged —
already uniform in `e`); tests in `tests/test_generic_hn.py`, `tests/test_hn_filtration.py`.

## 16. The F_e CAS cross-check + the first fully-live Macaulay2 suite run (E10-M4)

**The finding this closes.** The E13 re-audit's residual note: the P² ε-recursion common mode was
retired by the mutation oracle (§13), but the **F_e numerical layer** — the Gram matrices, `K_{F_e}`,
`χ(O)`, `hilbert_P`, and the RR Euler pairing that the envelope, the prioritary bounds, and the §14
generic-HN algorithm are all built on — remained single-sourced. Separately, the E10 M2-gated tests had
**never run on this host** ("six Macaulay2 skips hide this defect" was the audit's own R4 evidence).

**Provisioning (host).** Macaulay2 1.24.11 installed in WSL Debian (`apt install macaulay2`);
`scripts/m2-wsl.cmd` bridges Windows → WSL, translating script paths via `wslpath`, so the E10 oracle's
`BRIDGELAND_M2` discovery works from Windows pytest unchanged. Opt in per run:
`BRIDGELAND_M2=<repo>/scripts/m2-wsl.cmd pytest -q`. The `[UNVERIFIED on Windows]` /
`[UNVERIFIED idiom]` notes in `oracle/m2.py` are updated to **[VERIFIED]** — the four E10-M2 Ext tests
(P² Proj and the K3 Fermat quartic), the QQ round-trip, the E10-M3 witness construction, and the
E12-M4/R4-repaired mint all **pass live**.

**The cross-check** (`fe_line_bundle_cohomology` in `oracle/m2.py`; `tests/test_fe_cas.py`).
Macaulay2's `NormalToricVarieties` computes `h^i(F_e, O(D))` by **polytope lattice-point
combinatorics** — a route entirely independent of Riemann-Roch. The protocol is self-describing (the
transcript emits the prime-divisor classes and `−K` in M2's own Cl basis; no convention is trusted),
and the test **fits** the unimodular identification `T` with the package `(f, s)` basis from the data,
requiring `T(−K_{M2}) = (e+2, 2)` and full-table χ agreement `χ_{M2}(c,d) = P(T·(c,d))` over the
`[−3,3]²` window for `e ∈ {0,1,2,3}`. Any transcription error in the package's Gram / `K` / `χ(O)` /
`hilbert_P` yields **zero** fits.

A lattice-theoretic fact surfaced by the fit count (recorded, hand-verified): the number of
identifications equals the number of isometries of `(NS, Gram)` fixing `K` — **two** on `F_0` (the
ruling swap) and **two** on `F_2` (`σ: f ↦ f+s, s ↦ −s`; `σᵀGσ = G` and `σ(K) = K` exactly — both `F`
and `E+F` are isotropic on `F_2` and χ-data cannot distinguish the isotropic rays, only effectivity
can), **one** for `e = 1, 3`. Additional gates: `h^•(O) = (1,0,0)` pins; a package-free Serre-duality
self-consistency sweep of the M2 data (`h²(D) = h⁰(K−D)` within the window); and the M3b flagship's
character arithmetic CAS-witnessed (`χ(O(1,0)) + χ(O(0,1)) = 2 + 2 = 4` from the toric table equals the
package RR pairing on `(2,(1,1),0)`). The transcript parser also runs WITHOUT M2 via canned transcripts
— the plumbing is never skip-hidden (the R4 lesson).

**Suite:** default mode 533 items / 518 passed / **15 skips** (6 legacy + 9 new gated); with
`BRIDGELAND_M2` set: **533 / 533 / 0 skips / 0 failures** — the first fully-live run, ~266 s.
Re-run the gated tests after any Macaulay2 upgrade (print-format drift is the standing E10 risk).

*Sources:* Macaulay2 1.24.11 + the bundled `NormalToricVarieties` package (Gregory G. Smith et al.;
toric sheaf cohomology via polyhedral combinatorics — see the package documentation); the package-side
values cross-checked are those of §§7, 9–15 (arXiv:1907.06739 conventions). Files:
`bridgeland_stability/oracle/m2.py`, `scripts/m2-wsl.cmd`, `tests/test_fe_cas.py`, `CLAUDE.md`
(opt-in documented).

## 17. The sharp Bogomolov function δ_m^{μs} as a computable sandwich (E14-M1)

**What shipped.** `bridgeland_stability/delta_sharp.py`: (a) `mu_stable_exists(r, ν, Δ, surface)` — a
PROVEN decision procedure for "is there a `μ_{H_m}`-stable sheaf of character `(r, ν, Δ)`?" at any
rational `m > 0` on any strictly-ample `F_e`, honest `None` only on the single band `Δ = ½, r ≥ 2`;
(b) `delta_mu_stable(ν, m, surface, max_rank)` — the paper's headline function
`δ_m^{μs}(ν) = inf{Δ ≥ ½ : ∃ μ_{H_m}-stable sheaf of slope ν, discriminant Δ}` (`def-deltass`) as a
certified sandwich `DeltaSharp(lower, upper, exact)`.

**Two research findings that reshaped the milestone (recorded as an erratum to the E14 spec draft of
2026-07-16).** (1) *A rational `m` is never "generic".* The generic-polarization bridges
`prop-ssIMPs` (ss + `Δ > ½` ⟹ stable) and `prop-sIMPmus` (stable ⟹ μ-stable) hypothesize generic `m`;
on the rank-2 NS lattice, `(ν' − ν)·H_m = 0` has rational solutions iff `m ∈ ℚ` (in package `(f,s)`
coordinates `ξ·H_m = x + y m`), so every scannable polarization is special and the chain
`δ^p_{⌈m⌉+1} ≤ δ^{μss}_m ≤ δ^{ss}_m ≤ δ^s_m ≤ δ^{μs}_m` (paper §3.1) can be strict in every link —
an `hn_verdict` hit certifies `δ^{ss}`, NOT `δ^{μs}`. (2) *The inf need not be attained.* At the §8
Kronecker values the general sheaf AT the sharp discriminant is strictly μ-semistable
(`thm-intervalKronecker`) and the paper reaches `δ^{μs}` as a limit `Δ_{m±ε} → δ`; the drafted M1
criterion "the §8 values reproduced exactly at finite rank with `exact=True`" was therefore
**mathematically impossible** and is replaced by the sandwich/convergence pins below. (The M3a
lesson applied to a spec: an acceptance criterion can itself fabricate.)

**The decision procedure (each step a theorem; two-way evidence in `tests/test_delta_sharp.py`).**
Existence of a `μ_{H_m}`-stable sheaf of character `v` ⟺ `m ∈ I_v`, the *generic stability interval*
(slope stability is open in flat families and `P_F(v)` is irreducible — Walter — so one stable sheaf
makes the general sheaf stable). `I_v` is open, convex (slope stability passes to positive rational
combinations of ample classes; two dense-open loci of an irreducible stack meet), and contains the
anticanonical index `m₀ = 1 − e/2` whenever nonempty (`cor-KstabilityEasy`, `e ∈ {0,1}`). So for
`Δ > ½` the certifier samples the first **wall-free chamber** beside `m` on the side away from `m₀`:

* *The chamber gap.* Gieseker-ss existence of `v` is constant on `(m, m+g)` (mirror for left samples)
  with `g = 1/(32·Ymax·r²·q)`, `q = den(m)`, `Ymax = max(1, 2/(e+2m))`: every condition of the §5
  criterion (`thm-HNcriterion`/`cor-algorithm`) flips only where (i) a slope relation `ξ·H_{m'} ∈
  {0,1}` crosses, `ξ = ν_w − ν_u` a slope difference of recursion characters — coordinates in
  `(1/r²)ℤ`, `|ξ·F| ≤ 8·Ymax` (the `lem-slopeQuad` window stacked over recursion depth ≤ 4, doubled
  for pairs), giving `|m' − m| = |x + m y ∓ 1|/|y| ≥ 1/(8·Ymax·r²·q)`; (ii) the F-window boundary
  `|ξ·F| = 2/(e+2m')` crosses a candidate — same lattice bound up to a factor 2; (iii) an integer
  (the prioritary index `⌈m'⌉`) — at distance `≥ 1/q`. `g` under-runs all three.
* *Both directions.* Semistable at the rational chamber midpoint (one exact `hn_verdict` call) ⟹
  semistable at an **irrational** `m'` of the chamber (constancy) ⟹ `μ_{H_{m'}}`-stable sheaves exist
  (irrational ⟹ generic, `prop-ssIMPs` + `prop-sIMPmus`, `Δ > ½`) ⟹ `m ∈ [m₀, m'] ⊂ I_v` ⟹ exists at
  `m`. Not semistable there ⟹ no μ-stable sheaf anywhere in the open chamber (μ-stable ⟹ Gieseker-ss)
  ⟹ `I_v` misses it ⟹ `m ∉ I_v` (`I_v` open) ⟹ none at `m`.
* *`e ≥ 2`.* `cor-highermus` transports μ-stable existence bijectively along the E13-M1 reduction `π`
  (every strictly-ample `H` on `F_e` is an `A_m` in range; `r`, `ch₂`, `Δ` are π-invariant), so the
  question reduces to the del Pezzo base.
* *`Δ < ½`.* A μ-stable sheaf is simple with `Ext²(V,V) = Hom(V, V(K))* = 0` (`K·H < 0` on ample
  `F_e`), so `χ(v,v) = r²(1 − 2Δ) ≤ 1` (identity = Lemma "excFacts"(1), evaluated by the package RR
  `chi` and tripwired, never transcribed): `χ ≥ 2` refuses; `χ = 1` is exactly `cor-DLPExceptional`
  (`is_stable_exceptional`). Rank 1 is the ideal-sheaf test (`Δ = c₂ ∈ ℤ≥0`); `Δ < 0` is Bogomolov.

**The sandwich.** `lower = max(½, dlp_envelope value, δ^p_{⌈m⌉+1}(ν))` (`cor-deltaDLPe` + the §3.1
chain); `upper` = the least lattice `Δ > ½` of a scanned rank certified by the decision procedure
(per-rank first hits suffice: elementary modifications keep slope-stable sheaves slope-stable, so the
per-rank existence set is upward closed; termination ≤ one lattice step above the true value by
`thm-deltaSurface`(1) + totality). `exact = (upper == lower)`.

**Exact evidence (all pinned).**

| fact | value |
|---|---|
| `F₀`, `ν = (1/3,1/5)` pkg, `m = 25/9`: wall class `(15, ν, 3/5)` | `mu_stable_exists = False` (strictly μ-ss; inf not attained) |
| … one lattice step up `(15, ν, 2/3)` | `True`; scan `r ≤ 15`: sandwich `[19/35, 2/3] ∋ 3/5` |
| `F₁`, `ν = (6/13,3/13)` pkg, `m = 12/7`: wall class `(13, ν, 98/169)` | `False`; scan `r ≤ 13`: sandwich `[523/1014, 111/169] ∋ 98/169` |
| the paper's two `DLP^{<r}` computer values (§8) | `dlp_restricted` = `19/35` and `523/1014` **bit-for-bit** — first literature cross-check off the `−K` ray |
| anticanonical pinches | `ν = 0`: `δ = 1` exact (both `e`); `ν = (½,½)`/`F₀`: `δ = 3/4` exact, witness `(2,(1,1),−1)` |
| the E13 flagship cousin `(2,(1,1),0)`/`F₀` at `−K` | Gieseker-ss exists (`hn_verdict` True) **and** `mu_stable_exists = False` (`χ(v,v) = 2`) — the two stabilities separate on one class |
| exceptional branch vs Table "stabilityInterval1" | `(2,(1,1))`/`F₁` (`Δ = 3/8`): True at `m = ½ ∈ I_V = (0,1)`, False at `m = 3/2` |

Suite: 533 → 550 items (17 new), 15 default-mode skips unchanged, 0 failures.

*Source:* [arXiv:1907.06739](https://arxiv.org/abs/1907.06739) `def-deltass`, §3.1 (the δ-chain),
`thm-deltaGeneric`, `prop-ssIMPs`, `prop-sIMPmus`, `cor-KstabilityEasy`, `cor-deltaDLPe`,
`thm-deltaSurface`, `prop-divideSpace`, `lem-muss12`, `thm-intervalKronecker`, `thm-deltaKronecker`
(the §8 values), `cor-highermus`, Lemma "excFacts", `cor-DLPExceptional`; Walter, *Irreducibility of
moduli spaces of vector bundles on birationally ruled surfaces* (the prioritary stack). Package:
`bridgeland_stability/delta_sharp.py`; tests in `tests/test_delta_sharp.py`.

