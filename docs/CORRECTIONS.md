# Corrections to the project brief

This package implements the mathematics **as it actually is in the literature**,
which differs from the original project brief in several substantive ways. Every
correction below was verified two independent ways: by exact `fractions.Fraction`
computation, and against primary sources (read directly, not from memory). This
document records each discrepancy, the correct statement, and the citation.

The single most important convention choice: we use the **Coskun–Huizenga
normalized discriminant**

> Δ = ½·μ² − ch₂/(r·d),   d = H²,   μ = (ch₁·H)/(r·d)

(the brief used `Δ_brief = μ² − 2 ch₂/(r d) = 2·Δ`). The CH normalization is the
one in which every DLP / wall / BG formula in the literature is stated, so all
the explicit formulas below are clean in it. `ChernChar.discriminant_brief`
returns the doubled value when needed.

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
* The brief's wall trichotomy "δ²=−2/0/2" is wrong: the only invariants are
  **spherical s²=−2** and **isotropic w²=0** (no "+2" type — that was only the
  wrong-sign artifact of (1,0,−1)). The correct four-case classification is
  Bayer–Macrì Thm 5.7 (Brill–Noether / Hilbert–Chow / Li–Gieseker–Uhlenbeck
  divisorial, then flopping, then fake), implemented in `mukai.classify_wall`.

*Source:* A. Bayer, E. Macrì, "MMP for moduli of sheaves on K3s via
wall-crossing", [arXiv:1301.6968](https://arxiv.org/abs/1301.6968), Thm 2.15
(dimension) and Thm 5.7 (classification).
