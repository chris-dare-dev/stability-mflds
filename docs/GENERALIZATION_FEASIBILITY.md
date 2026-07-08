# Generalization Feasibility Study

*Feasibility of moving `bridgeland_stability` beyond the numerical (Chern-character),
Picard-rank-1 regime it currently occupies. Every code claim is grounded in the actual files
(`chern.py`, `walls.py`, `mukai.py`, `varieties.py`, `dlp.py`, `exceptional.py`, `threefold.py`).
Every mathematical claim is tagged `[PROVEN]` / `[CONJECTURED]` / `[SPECULATIVE]` / `[UNVERIFIED]`,
and the three-way distinction below is maintained throughout:*

- **(a)** computable in principle by a known algorithm;
- **(b)** already implemented in existing software (named);
- **(c)** implementable *in this codebase's exact-`Fraction`, numerical-Chern-character architecture,
  with no sheaf-level computation*.

To the author's knowledge there is **no maintained public package** covering (b) for these surface/threefold
wall and existence algorithms; existing implementations are bespoke research scripts (e.g. B. Schmidt's
tilt-stability computations). So the interesting column throughout is (c).

---

## 0. What the code actually is today (baseline)

- `chern.py` fixes a class to a **scalar triple** `(r, c = ch1·H, e = ch2)` plus `d = H²`. The central
  charge is `Z_{s,t} = −∫ e^{−(s+it)H} ch` — **no Todd twist**.
- `walls.py::numerical_wall(v,w,d)` is exact and convention-clean: minors `W_rc, W_re, W_ce`, center
  `s₀ = W_re/W_rc`, `ρ² = s₀² − 2 W_ce/(d·W_rc)`. It sees only `d` and the three scalars.
- `mukai.py` is **Picard-rank-1 K3 only**: `MukaiVector(r,l,s)` with `s = ch2+r`, pairing
  `d l l' − r s' − r' s`, and `classify_wall` (Bayer–Macrì Thm 5.7) that **searches integer combinations**
  `a v + b w` for spherical/isotropic witnesses. It produces a wall *type*, **never an (s,t) semicircle**.
  Note `MukaiVector.from_chern` hard-codes the `ch2+r` shift, which is **K3-specific** (see §1b, §1c).
- `varieties.py::Surface` stores `d, K_H, chi_O, picard_rank` — `picard_rank` is a bare **integer count,
  not a lattice**. `P1xP1` and `hirzebruch(n)` carry `picard_rank=2` but only a single degree `d`.
- `threefold.py` has `ThreefoldChern(r,a1,a2,a3)`, the BMT form `Q = 4 ch2² − 6 ch1 ch3`, and
  `alpha_crit`; there is **no `ν_{α,β}` tilt-slope function** and no wall-crossing.
- `dlp.py`/`exceptional.py` are **P²-only**, `d = 1` hardcoded, Markov recursion.

Everything below is measured against this baseline.

> **Prioritization headline (read before the roadmap).** The genuinely cheap, *in-architecture* wins are
> all **Picard-rank-1**: K3 `(s,t)` walls (a small Mukai-coordinate shim), abelian-surface `(s,t)` walls
> (**zero code change** — see §1c), Picard-1 threefold tilt walls, and higher-rank surface BG bounds. The
> rational surfaces that "beyond P²" most evokes — `P¹×P¹`, Hirzebruch `𝔽_n`, del Pezzo — are **all** gated
> behind the expensive NS-lattice refactor (§2a) **and** a partly *sheaf-theoretic* (prioritary + Harder–
> Narasimhan) algorithm whose core is not pure Chern arithmetic. Do **not** schedule rational surfaces as
> near-term.

---

## 1. Easiest generalizations (highest research value per unit effort)

### 1a. Coskun–Huizenga non-emptiness beyond P²

**On P²** the non-emptiness of `M(r,c₁,ch₂)` is *completely* solved by the Drézet–Le Potier curve
`δ(μ)` and the exceptional (Markov) bundles — a **closed form**, and it is already implemented here
(`dlp.delta`, `exceptional.enumerate_exceptional`). `[PROVEN]` (Drézet–Le Potier 1985; Coskun–Huizenga
survey). This is category (c)-complete already.

**On Hirzebruch surfaces** (`P¹×P¹ = 𝔽₀`, `𝔽_n`) the existence problem **is algorithmic**, but *not* as a
single closed-form `δ`-curve. Coskun–Huizenga (arXiv:1907.06739, *Adv. Math.* 2021) give a **decision
procedure**: pass to stacks of *prioritary* sheaves, solve existence there, then compute the
Harder–Narasimhan filtration of the generic prioritary sheaf; semistable sheaves exist **iff that HN
filtration has length one**. They also produce *sharp* Bogomolov-type bounds `Δ ≥ δ_H(ξ)` that **depend on
the polarization and the full slope**. `[PROVEN]` for Hirzebruch surfaces. This is category (a) (algorithm
exists) and **(c)-feasible only for the numerical output** — the HN/prioritary *core* of the algorithm is
sheaf-theoretic and lives outside this codebase's constraint (see §4). It is **not** a small port of
`dlp.py`.

**Precisely what changes vs P²** (this is the crux):

1. **Two slope parameters, not one.** With `ρ(X)=2`, the discrete invariant `c₁` is a *lattice vector* in
   `NS(X) = ℤf ⊕ ℤs`, so "slope" is a point in `NS(X)⊗ℚ` (two coordinates), and the existence function lives
   over a **2-dimensional** base, not the real line. The current scalar `c = ch1·H` cannot represent it.
   `[PROVEN]` (this is the structural content of 1907.06739).
2. **The polarization matters.** The sharp bound `δ_H` genuinely depends on `H`; different chambers of the
   ample cone give different existence loci. On P² there is only one ray, so this dependence is invisible.
   `[PROVEN]`.
3. **The exceptional collection changes.** The clean Markov/`ε`-recursion is a P²-specific miracle.
   Hirzebruch and del Pezzo surfaces have their **own** exceptional collections (Kuleshov–Orlov on del
   Pezzo; Rudakov), so `exceptional.py`'s Markov machinery does **not** transfer; a surface-specific
   exceptional-bundle generator is required. `[PROVEN]` that they exist. `[SPECULATIVE]` — and, on the
   evidence, **doubtful** — that a clean uniform recursion as tidy as Markov exists for each: Rudakov and
   Kuleshov–Orlov supply these **case by case, with no single Markov-style recursion**. This is a genuine
   per-surface cataloguing effort, not a reusable recursion (see node E1, §3).

**Del Pezzo surfaces.** Levine–Zhang (arXiv:1910.14060, *Ann. Inst. Fourier* 74, 2024;
*Brill–Noether and existence of semistable sheaves on del Pezzo surfaces*) settle **existence for del Pezzo
of degree ≥ 3** and **general-sheaf cohomology for degree ≥ 4**, in both cases for the **anticanonical
polarization**, using exceptional bundles. `[PROVEN]` for those cases; degree ≤ 2 del Pezzos and arbitrary
(non-anticanonical) polarizations are only partially covered. The rational-surface backbone is Coskun–
Huizenga's **weak Brill–Noether for rational surfaces** (arXiv:1611.02674) `[PROVEN]`.

**Is a `δ(μ)`-analogue implementable for P¹×P¹?** `[SPECULATIVE]` **Only after a lattice refactor, and only
for the numerical verdict.** A `δ_H(ξ)`-analogue is category (c)-feasible *provided* the codebase first (i)
represents `ch₁` as an `NS`-lattice vector, (ii) stores the intersection form, and (iii) grows a
`𝔽_n`-specific exceptional-bundle generator. Without (i)–(iii) the current scalar architecture literally
cannot express the input, so a naive port would be silently wrong. Even with (i)–(iii), the honest
deliverable is the **numerical decision verdict** (a tabulated/recomputed `δ_H` and yes/no), **not** the
sheaf-level HN/prioritary proof engine, which is out of architecture (§4).

### 1b. Are Arcara–Bertram / K3 walls already reproduced by `numerical_wall`?

This has a **clean, partly surprising answer**, and it splits by Picard rank.

**Picard rank 1 (the K3/abelian case Arcara–Bertram actually treat).** Arcara–Bertram (arXiv:0708.2247,
*JEMS* 2013) work with `H` generating `Pic(S)`, i.e. `ρ = 1`. There the `(s,t)`-slice wall between two
classes is the locus where two central charges of the form `∫ e^{−(s+it)H}(·)` are parallel — exactly the
locus `numerical_wall` computes. **So the geometry is the right shape.** `[PROVEN]` that the wall is such a
semicircle.

**But the K3 central charge carries a Todd twist that the code omits.** Bridgeland's K3 stability function
(arXiv:math/0307164) is `Z_{s,t}(E) = ⟨exp((s+it)H), v(E)⟩` with the **Mukai vector**
`v(E) = ch(E)·√td`, and `√td_{K3} = (1,0,1)` — a fact `mukai.py`'s own docstring states — so the third
coordinate is `ch₂ + r`, **not** `ch₂`. `chern.py::central_charge` uses bare `ch₂`. Therefore:

- **`numerical_wall` fed a K3's bare `(r,c,ch₂)` returns the *tilt-stability* (ch-based) wall, not the
  genuine Arcara–Bertram/Bridgeland K3 wall.** `[PROVEN]` from reading `chern.py` + `walls.py` (no `√td`
  anywhere) against the standard K3 central charge.
- **The fix is trivial and needs no new machinery:** call `numerical_wall` on the **Mukai coordinates**
  `(r, c, ch₂+r)`. `[PROVEN for ρ=1 — exact derivation, independently reproduced]` Under `e ↦ e+r`,
  `e' ↦ e'+r'`: `W_re` is **unchanged**, `W_rc` is unchanged, and `W_ce ↦ W_ce − W_rc`. Hence the
  **center is identical** and the **radius² shifts by exactly `+2/d`**. So the current output is off by a
  known, exact `2/H²` in `ρ²` — small, but wrong, and it will fail any pinned literature K3 wall. This is a
  *theorem*, not a guess: because `√td_{K3} = (1,0,1)`, feeding the Mukai triple `(r, c, ch₂+r)` to
  `chern.central_charge` reproduces Bridgeland's K3 central charge `⟨exp((s+it)H), v(E)⟩` exactly
  (Bridgeland arXiv:math/0307164).

**Does the Mukai pairing "replace the scalar minors"?** `[PROVEN for ρ=1]` Effectively yes: the correct
radius is what you get by running the *same* `(r,c,e)` minor formula on the **Mukai** third coordinate. The
scalar minors are not wrong; they are being fed the wrong invariant. Structurally, `numerical_wall`'s
`(W_rc, W_re, W_ce)` are the antisymmetrized 2×2 minors of the two classes, and on a `ρ=1` K3 the Mukai
pairing `d l l' − r s' − r' s` **is** the relevant bilinear form — so once you pass Mukai vectors the two
descriptions coincide.

**Crucially, the two K3 modules are disconnected in the current code.** `[PROVEN]` `mukai.py` never calls
`numerical_wall`, and `numerical_wall` never sees a `MukaiVector`. `mukai.py` gives the *type* of a wall
(divisorial / flopping / fake, Bayer–Macrì Thm 5.7) — which is **already correct for `ρ=1`** `[PROVEN]` —
but **no (s,t) geometry**. So today: K3 wall *classification* exists; K3 wall *semicircles in (s,t)* do
**not**, and the surface `numerical_wall` would give the wrong radius by `+2/d` if naively reused.

**Picard rank ≥ 2 (P¹×P¹, higher-`ρ` K3).** Arcara–Bertram's explicit formulas do **not** apply (they
assumed `ρ=1`). Maciocia (arXiv:1202.4587) shows the walls then require the **full Néron–Severi intersection
form**: he decomposes `β = bω + γ` with `γ ⟂ ω` and the semicircle depends on `γ² = −d ≤ 0` (Hodge index).
`[PROVEN]`. The current scalar architecture (one degree `c = ch1·H`) **cannot see the `γ`-orthogonal
component**, so on P¹×P¹ `numerical_wall` computes only the *H-projected* walls and is **incomplete** — it
misses destabilizers whose class has the same `ch1·H` but different bidegree. `[SPECULATIVE — direct
consequence of Maciocia's decomposition + reading the scalar code]`.

**Bottom line for §1b:** the cheapest high-value win for K3 is a small addition (`[SPECULATIVE]` ~20 lines) —
a `k3_wall(v, w, d)` that maps `(r,c,ch₂) ↦ (r,c,ch₂+r)` and calls the existing `numerical_wall`, plus a
pinned test against a published K3 wall. Category (c), no new data structures. This shim is **K3-only** (do
**not** reuse it for abelian surfaces — see §1c). Everything for `ρ ≥ 2` is gated behind the NS-lattice
refactor (§2a).

### 1c. Abelian surfaces: `numerical_wall` is already correct — a free win

**An abelian surface is the exception that proves the rule, and it is NOT the K3 case.** `[PROVEN]` An
abelian surface has **trivial tangent bundle**, so its Todd class is `td = (1,0,0)` and therefore
`√td = (1,0,0)`. The Mukai vector is `v(E) = ch(E)·√td = (r, c₁, ch₂)` with **no `+r` shift**. Consequently:

- **`numerical_wall` on the bare `(r,c,ch₂)` already returns the genuine abelian `(s,t)` wall** — with
  **zero code change**. By the same reduction used in §1b, feeding the *bare* Chern triple to
  `chern.central_charge` reproduces Bridgeland's abelian central charge `⟨exp((s+it)H), v(E)⟩` exactly,
  precisely because `√td = (1,0,0)`. `[PROVEN for ρ=1]`.
- **The K3 `+2/d` correction and `MukaiVector.from_chern`'s `ch₂+r` are K3-SPECIFIC and would be WRONG if
  reused for abelian surfaces.** `[PROVEN]` Applying the K3 shim to an abelian class would introduce a
  spurious `+2/d` error in `ρ²`. Any future `k3_wall` helper (§1b) must be documented as **K3-only** and
  must **not** be pointed at abelian input.

This is the single cheapest item in the entire study: the correct abelian `(s,t)` walls are *already* what
`numerical_wall` computes, so exposing them costs nothing beyond a pinned test and documentation. (Wall
*classification* on abelian surfaces uses a different Mukai-lattice signature than K3 and is out of scope
here; `mukai.classify_wall` is a K3 construction and should not be reused for abelian surfaces either.)

---

## 2. Prerequisites for each deeper generalization

### 2a. Maciocia 1202.4587 walls on general surfaces → **needs the full NS lattice**

`[PROVEN]` Maciocia's algorithm requires **the intersection form on `NS(X)`, not merely `d = H²`.** His
wall/center/radius formulas use `γ²` (the self-intersection of the polarization-orthogonal part of the
Chern character), and his boundedness result (Thm 3.11: all walls of a class lie inside one bounding
semicircle; for all but one ray the wall count is globally finite) is stated in those terms. For `ρ = 1` the
walls are *nested* and collapse to the H-projection — which is exactly why the current code is correct on P²
and on `ρ=1` K3/abelian surfaces but nowhere else.

**Data structure the codebase must add** (category (c) prerequisite, and the single highest-leverage item):

1. An **`NSLattice`**: a rank-`ρ` free ℤ-module with a **symmetric bilinear form** (an integer Gram matrix).
   For P¹×P¹: `⟨f,f⟩=0, ⟨s,s⟩=0, ⟨f,s⟩=1`; for `𝔽_n`: the standard `[[0,1],[1,−n]]`.
2. `ch₁` promoted from a **scalar** `c` to an **`NS` lattice vector**; `ch₂` stays a `Fraction`.
3. `Surface` gains the Gram matrix and the ample class `H` as a lattice vector; `slope`/`discriminant` use
   the true `ch₁²` (currently proxied by `(ch1·H)²/d`, an upper bound by Hodge index).

This is a **frozen-dataclass-friendly, exact-`Fraction`, sheaf-free** refactor — fully in-architecture. It
is the parent node of almost every `ρ>1` generalization (see §3). **Effort: high / strategic**, `[SPECULATIVE]`
on the exact magnitude but decisively **not** a "medium" change: it mutates the **core frozen `ChernChar`**
(`c` scalar → NS lattice vector), and with it `walls.py`, `mukai.py`, and **every one of the 42 pinned
tests** (whose literal `(r, c, e)` inputs would change shape). Treat it as a v2 architecture decision, not a
weekend refactor. `[PROVEN]` that it is the mathematically necessary input.

### 2b. Threefolds: a `ν_{α,β}` slope function and tilt-stability walls → **cheap for the slope; wall locus needs a check**

The tilt slope is a **known closed formula** `[PROVEN]` (standard; BMT arXiv:1103.5010, Schmidt):
```
ν_{α,β}(E) = ( H^{n-2}·ch₂^β(E) − ½ α² H^n·ch₀^β(E) ) / ( H^{n-1}·ch₁^β(E) ),   +∞ if ch₁^β·H^{n-1}=0.
```
`ThreefoldChern` already stores `(r, a1, a2, a3)` and `twist(β)`, so `ν_{α,β}` is a small pure-`Fraction`
function (`[SPECULATIVE]` ~5 lines) — category (c), no lattice, no sheaves.

The **locus where an object changes tilt-stability** (a tilt wall: `ν_{α,β}(F) = ν_{α,β}(E)` for a fixed
potential sub-`F`) is a curve in the `(α,β)` half-plane. Schmidt (arXiv:1509.04608, *Bridgeland Stability on
Threefolds — Some Wall Crossings*) works these out explicitly. **Caveat to resolve before pinning any
tests** `[SPECULATIVE]`: on a threefold `ν_{α,β}` involves a **cubic `ch₃^β` term** (which `ThreefoldChern`
carries via `a3`), so it is **not automatic** that the tilt-wall locus reduces to the same scalar-2×2-minor
conic as the surface `numerical_wall`. Before writing `threefold_tilt_wall(v, w, d3)` or pinning literature
values, **confirm from Schmidt 1509.04608 that the `(α,β)` tilt-wall locus is the expected conic despite the
cubic term**. If it is, a `threefold_tilt_wall` mirroring `numerical_wall` is category (c)-feasible; the
"mirrors `numerical_wall`" claim and any effort/line-count figures stay `[SPECULATIVE]` until that check is
done.

**What the literature gives algorithmically on threefolds beyond the BG boundary** (all `[PROVEN]` unless
noted):

- **Tilt-stability walls** (numerical, in `(α,β)`): explicit — Schmidt 1509.04608. Category (a); (c) pending
  the conic check above.
- **The tilt-BG inequality `Q ≥ 0`** is a **theorem** for P³ / all Fano `ρ=1` (Li arXiv:1510.04089), abelian
  3-folds (Bayer–Macrì–Stellari arXiv:1410.1585), the quintic (Li arXiv:1810.03434), the quadric (Schmidt),
  and — via **Koseki (arXiv:1811.03267, *Stability conditions on threefolds with nef tangent bundles*, Adv.
  Math. 372, 2020)** — **threefolds with nef tangent bundle** and the products `P¹×S`, `P²×C`, `P¹×P¹×C`
  (`S` an abelian surface, `C` an elliptic curve). `[UNVERIFIED]` — before adding any `bg_proven=True`
  catalog row for these products, **verify the exact hypotheses on `S` and `C`** against 1811.03267.
  (Koseki's separate weighted-hypersurface results — stability conditions on Calabi–Yau double/triple
  solids, e.g. hypersurfaces in `P(1,1,1,1,2)` / `P(1,1,1,1,4)` — are arXiv:2007.00044, a **different**
  paper; do not conflate the two IDs. arXiv:1510.04474 is instead Chuang–Lai's *withdrawn* CY/Fano BG
  paper — not Koseki at all.) The general
  threefold BG statement is a `[CONJECTURED]` inequality (BMT 1103.5010), and its *stronger* form is
  `[PROVEN]` false on `Bl_p(P³)` (Schmidt arXiv:1602.05055) — already flagged by `Threefold.bg_proven=False`.
- **Full Bridgeland stability + proper moduli** on a threefold follow *conditionally* on BG: Piyaratne–Toda
  (arXiv:1504.01177) prove the moduli stacks are proper of finite type **assuming** the BG conjecture. So a
  category-(c) implementation can compute the tilt-stability picture unconditionally, but any statement about
  *actual* Bridgeland walls/moduli must be tagged conjectural exactly where `bg_proven` is `False`.

**Prerequisite chain for threefold wall-crossing:** `ν_{α,β}` function → (conic check) → tilt-wall solver →
(optionally) a `bg_proven`-gated "is this a genuine Bridgeland wall" verdict. No sheaf theory needed for the
numerical layer; the derived-category content stays outside the code, as it already does on surfaces.

### 2c. K3 `(s,t)` walls → **needs only the K3-specific Mukai-coordinate shim** (see §1b). Abelian surfaces need **nothing** (§1c). No lattice for `ρ=1`.

### 2d. Higher-rank / mutation / quivers → **needs the derived category** (out of architecture; see §4).

---

## 3. Dependency graph of generalizations

Nodes are ordered so that each depends only on earlier ones. `→` means "unlocks".

```
[N0]  NSLattice + symmetric form + ch1-as-vector          (HIGH/STRATEGIC; foundation for all ρ>1 work)
[K1]  K3 (s,t) walls via K3-only Mukai-coordinate shim     (INDEPENDENT of N0; ρ=1 only)
[A1]  abelian (s,t) walls — ALREADY CORRECT (free win)     (INDEPENDENT; zero code change, see §1c)
[T1]  ν_{α,β} tilt-slope function on ThreefoldChern         (INDEPENDENT; uses existing twist)
[R0]  RigorStatus provenance tag + propagation             (INDEPENDENT; orthogonal, §5)

[N0] → [M1]  Maciocia full-NS surface walls (P1xP1, 𝔽_n, del Pezzo)
[N0] → [E1]  surface-specific exceptional-bundle generators (𝔽_n, del Pezzo) — CASE-BY-CASE, no Markov recursion
[M1]+[E1] + non-Chern HN step → [X1]  Coskun–Huizenga existence / δ_H decision procedure on rational surfaces
[T1] + conic check → [T2]  threefold tilt-stability wall solver (Schmidt-style)
[T2] + bg_proven → [T3]  conjecture-gated Bridgeland-wall verdicts on threefolds
[R0] → colours every node's output proven/conjectural/heuristic
```

Adjacency list:

- `N0`: prerequisites = none. **High-effort / strategic (v2 architecture decision)** — mutates the core
  frozen `ChernChar` and all 42 pinned tests. Unlocks `M1`, `E1` (hence `X1`).
- `K1`: prerequisites = none. Unlocks correct **K3** wall diagrams (pairs with existing
  `mukai.classify_wall`). The Mukai shim is **K3-only**.
- `A1`: prerequisites = none. **Free win** — abelian `(s,t)` walls are already what `numerical_wall`
  returns (§1c); deliverable is a pinned test + docs, no new code.
- `T1`: prerequisites = none. Unlocks `T2` (after the §2b conic check).
- `M1`: prerequisites = `N0`. Unlocks certified surface wall enumeration for `ρ>1`.
- `E1`: prerequisites = `N0`. **Case-by-case per-surface cataloguing** (Rudakov / Kuleshov–Orlov), **no
  uniform Markov-style recursion**. Unlocks `X1`.
- `X1`: prerequisites = `M1`, `E1`, **and a non-Chern-level HN/prioritary step**. See the scope caveat below.
- `T2`: prerequisites = `T1` **and** the §2b conic check. Unlocks `T3`.
- `T3`: prerequisites = `T2` + `bg_proven`. Conjecture-gated.
- `R0`: prerequisites = none. Orthogonal; should be built early because everything downstream is a mix of
  proven and conjectural.

**X1 is not a near-term deliverable.** `[SPECULATIVE]` The Coskun–Huizenga `δ_H` existence decision procedure
on rational surfaces sits behind **both** `N0` (the high-effort NS-lattice refactor) **and** `E1` (a
case-by-case per-surface exceptional-bundle cataloguing effort with no reusable recursion), and its
algorithmic core is a **Harder–Narasimhan / prioritary-sheaf** computation that is **not pure Chern
arithmetic** (only its numerical *output* is; see §4). It is therefore a **multi-month, per-surface-family
research effort**, not a days-scale port. Do not place it in any near-term tier.

**Recommended build order (value/effort).** Near-term, in-architecture, high-value (days-scale):
`A1` (free) → `K1` → `T1` → `R0`. Then the strategic refactor `N0` (high effort), which unlocks the
`ρ>1` surface work `M1`, and the threefold solver `T2` (after the conic check). The per-surface research
program `E1`/`X1` and the conjecture-gated `T3` come last and are **not** days-scale. **Note the ordering
trap:** the rational surfaces (`P¹×P¹`, `𝔽_n`, del Pezzo) that "generalize beyond P²" most evokes are the
*expensive* end of this list, gated behind `N0` **and** the partly sheaf-theoretic `E1`/`X1`; the cheap wins
are all Picard-rank-1.

---

## 4. What is mathematically impossible or open-problem-hard

Sharp separation of **open problem** (nobody knows the answer) from **just hard/large** (known but heavy)
from **out of architecture** (needs sheaves/derived category this codebase forbids).

**Open-problem-hard (do not promise these):**

- **DLP-style *sharp* non-emptiness on surfaces of general type / `p_g > 0`.** Beyond the Bogomolov
  inequality `Δ ≥ 0`, essentially **nothing sharp is known** for general-type surfaces; the existence region
  is not described by any known `δ`-curve. `[CONJECTURED/OPEN]` (survey consensus; cf. discriminant studies
  on general-type surfaces, arXiv:1812.02735). A `δ`-analogue here is **not** category (a) — there is no
  algorithm to implement. Even on K3 (`p_g = 1`), the *sharp* existence bound is subtler than on rational
  surfaces, though the K3 lattice makes the *moduli-dimension* count (`v²+2`, already in `mukai.py`) clean.
- **The threefold Bogomolov–Gieseker conjecture in general.** `[CONJECTURED]` (BMT 1103.5010), and **false**
  in its stronger form on `Bl_p(P³)` `[PROVEN]` (Schmidt 1602.05055). For any threefold **not** on the
  proven list, the BG boundary the code draws is **not rigorous** — this is a genuine open problem, correctly
  surfaced by `bg_proven`. **Important distinction** `[PROVEN]`: "the strong BMT boundary is not rigorous on
  `X`" is **not** the same as "`X` carries no stability condition." `Bl_p(P³)` **is Fano and does carry
  Bridgeland stability conditions** (Bernardara–Macrì–Schmidt–Zhao, arXiv:1607.08199); what fails there is
  only the *strong/naive* BMT inequality (Schmidt 1602.05055). Any code note for `BLOWUP_P3_POINT` should say
  "**strong BMT boundary not rigorous**", **not** "no stability condition exists" (the latter is false). The
  code must never assert rigor where the theorem is absent — but it also must not overstate the failure.
- **`actual_walls` completeness for arbitrary `v`.** The current `actual_walls`/`actual_walls_complete` give
  the **necessary** conditions and an *empirical* doubling certificate. Provable completeness of the wall set
  for a *general* class `v` (all `ρ`, all ranks) is **not** a solved problem; it is a theorem only in the
  ABCH / Coskun–Huizenga-covered cases (`P²[n]`, coprime/small-rank). `[CONJECTURED/OPEN]` in general;
  `[PROVEN]` in the covered cases the tests pin. Do **not** upgrade the docstring's honesty here.

**Just hard/large (feasible, but heavy):**

- Maciocia walls for high `ρ` — feasible after `N0`, but the wall enumeration grows with `ρ` and needs
  careful bounding (Maciocia Thm 3.11 gives the bound). Category (a)+(c), large.
- Exceptional collections on each Fano/rational surface — known case by case (Kuleshov–Orlov, Rudakov), **no
  single uniform recursion**; a real cataloguing effort. `[PROVEN]` existence, `[SPECULATIVE]` uniform
  tidiness (and, per §1a, likely *no* Markov-style recursion exists).

**Out of architecture (would violate the zero-sheaf constraint):**

- Anything requiring **Ext groups, resolutions, HN filtrations of actual objects, mutations, quivers, the
  Serre functor, or the derived category**. The prioritary-sheaf existence algorithm (§1a) is genuinely
  *sheaf-theoretic* in its proof, but its **output** (the sharp `δ_H` and the yes/no verdict) is numerical
  and can be *tabulated/recomputed* at the Chern-character level — so the **decision function** is (c), while
  a from-scratch *proof engine* is not. Be honest about which is being shipped.

---

## 5. The conjecture-mode question: a general provenance mechanism

`Threefold.bg_proven` is the right instinct but too narrow: it is one boolean on one dataclass. The general
pattern is a **rigor/provenance tag carried on every result and propagated like a taint lattice.**

**Sketch (in-architecture, `Fraction`-clean, sheaf-free):**

```python
class Rigor(IntEnum):          # a total order: PROVEN is strongest
    PROVEN      = 3            # a cited theorem covers these hypotheses
    CONJECTURAL = 2            # holds modulo a named open conjecture (e.g. threefold BG)
    HEURISTIC   = 1            # numerical/empirical (dense compute_walls; doubling certificate)
    UNKNOWN     = 0

@dataclass(frozen=True)
class Certificate:
    rigor: Rigor
    hypotheses: tuple[str, ...]     # e.g. ("rho(X)=1", "H generates Pic")
    citations: tuple[str, ...]      # arXiv IDs / DOIs
    note: str = ""
```

- **Attach a `Certificate` to variety objects and to results.** `Surface`/`Threefold` list which theorems
  apply (P² DLP: `PROVEN`; `Bl_p(P³)` strong BMT boundary: `CONJECTURAL`/`UNKNOWN` — *not* "no stability
  condition"). A `Wall`/`ActualWall`/`BGBoundary` carries the `Certificate` under which it is valid.
- **Propagate by the *minimum* (meet) over inputs.** A `PROVEN` algorithm run on a `CONJECTURAL` BG input
  yields a `CONJECTURAL` result: `result.rigor = min(algorithm_rigor, *input_rigors)`. This is exactly the
  `bg_proven` logic, generalized to a lattice so it composes across the whole pipeline.
- **Hypothesis checking is data, not proof.** The tag records *which* theorem is being invoked and whether
  its stated hypotheses (`ρ=1`, coprimality, polarization genericity, BG-proven) are met by the object; it
  does **not** attempt to prove anything. This keeps it sheaf-free and honest.
- **Rendering.** `viz`/`__repr__` shows the tag (solid line = `PROVEN`, dashed = `CONJECTURAL`, dotted =
  `HEURISTIC`), matching the existing `bg_boundary_curve` warning behaviour. `docs/CORRECTIONS.md`'s
  "verified two ways" culture becomes machine-readable.

This mechanism is category (c), orthogonal to all math work (`R0` in §3), and should be built early: it is
the engineering answer to "the same function is a theorem on P² and a conjecture on a general threefold."

---

## Per-family verdict summary

| Family | Non-emptiness / DLP | Surface walls (s,t) | Status |
|---|---|---|---|
| **P²** | closed-form `δ(μ)`, **implemented** | correct (`ρ=1`) | `[PROVEN]`, in code |
| **P¹×P¹, 𝔽_n** | algorithmic (prioritary+HN), sharp `δ_H`; numerical output only | needs full NS form | `[PROVEN]` math; needs `N0` |
| **del Pezzo (deg ≥ 3 existence; ≥ 4 cohomology)** | solved for anticanonical `H` | needs full NS form | `[PROVEN]` (Levine–Zhang) |
| **K3 (`ρ=1`)** | moduli dim `v²+2` in code; sharp existence subtle | center OK, radius off by `+2/d`; fix via **K3-only** Mukai shim | `[PROVEN]` type; `[PROVEN, ρ=1]` radius fix |
| **abelian surface (`ρ=1`)** | not treated here | **already correct** — bare `numerical_wall` IS the genuine wall; **no `+r` shift** (free win) | `[PROVEN, ρ=1]` no shift |
| **general type / `p_g>0`** | **no sharp `δ` known** | Bogomolov only | `[CONJECTURED/OPEN]` |
| **Fano `ρ=1`, abelian, quintic, quadric 3-folds; nef-tangent-bundle 3-folds** | — | tilt-BG boundary rigorous | `[PROVEN]`, `bg_proven=True` |
| **`Bl_p(P³)` and general 3-folds** | — | strong BMT boundary not rigorous (but stability conditions **do** exist) | `[CONJECTURED]`/strong-form-false |

---

## References

- Arcara, Bertram — *Bridgeland-stable moduli spaces for K-trivial surfaces*, arXiv:0708.2247 (JEMS 2013).
- Arcara, Bertram, Coskun, Huizenga (ABCH) — *The minimal model program for the Hilbert scheme of points on
  P² and Bridgeland stability*, arXiv:1203.0316.
- Bayer, Macrì — *MMP for moduli of sheaves on K3s via wall-crossing*, arXiv:1301.6968 (Thm 2.15, Thm 5.7).
- Bayer, Macrì, Toda (BMT) — *Bridgeland stability conditions on threefolds I: Bogomolov–Gieseker type
  inequalities*, arXiv:1103.5010.
- Bayer, Macrì, Stellari — *Stability conditions on abelian threefolds*, arXiv:1410.1585.
- Bernardara, Macrì, Schmidt, Zhao (BMSZ) — *Bridgeland stability conditions on Fano threefolds*,
  arXiv:1607.08199. (`Bl_p(P³)` is Fano and carries stability conditions; only the strong BMT boundary fails
  there.)
- Bridgeland — *Stability conditions on K3 surfaces*, arXiv:math/0307164. (K3 central charge
  `⟨exp((s+it)H), v(E)⟩`; `√td_{K3}=(1,0,1)`.)
- Coskun, Huizenga — *The birational geometry of the moduli spaces of sheaves on P²* (Gökova survey), §4–6.
- Coskun, Huizenga — *Weak Brill–Noether for rational surfaces*, arXiv:1611.02674.
- Coskun, Huizenga — *Existence of semistable sheaves on Hirzebruch surfaces*, arXiv:1907.06739 (Adv. Math.
  2021).
- Drézet, Le Potier — *Fibrés stables et fibrés exceptionnels sur P₂*, Ann. Sci. ENS 18 (1985), 193–243.
- Koseki — *Stability conditions on threefolds with nef tangent bundles*, arXiv:1811.03267 (Adv. Math. 372,
  2020). BG for `P¹×S`, `P²×C`, `P¹×P¹×C` (`S` abelian surface, `C` elliptic curve) and nef-tangent-bundle
  threefolds. **Verify exact hypotheses before adding any `bg_proven=True` catalog row.**
- Koseki — *Stability conditions on Calabi–Yau double/triple solids* (weighted hypersurfaces in
  `P(1,1,1,1,2)`, `P(1,1,1,1,4)`), arXiv:2007.00044. (A **distinct** paper from 1811.03267. NB:
  arXiv:1510.04474 is Chuang–Lai's *withdrawn* CY/Fano BG paper — not Koseki — and was a prior miscitation.)
- Levine, Zhang — *Brill–Noether and existence of semistable sheaves on del Pezzo surfaces*,
  arXiv:1910.14060 (Ann. Inst. Fourier 74, 2024). Existence for degree ≥ 3; general-sheaf cohomology for
  degree ≥ 4; anticanonical polarization.
- C. Li — *Stability conditions on Fano threefolds of Picard number one*, arXiv:1510.04089.
- C. Li — *On stability conditions for the quintic threefold*, arXiv:1810.03434.
- Maciocia — *Computing the walls associated to Bridgeland stability conditions on projective surfaces*,
  arXiv:1202.4587 (Asian J. Math. 18, 2014).
- Piyaratne — *Stability conditions, Bogomolov–Gieseker type inequalities and Fano 3-folds*, arXiv:1705.04011.
  (Solo Piyaratne paper — **not** Bernardara–Macrì–Schmidt–Zhao.)
- Piyaratne, Toda — *Moduli of Bridgeland semistable objects on 3-folds and Donaldson–Thomas invariants*,
  arXiv:1504.01177.
- Schmidt — *Bridgeland stability on threefolds — some wall crossings*, arXiv:1509.04608.
- Schmidt — *Counterexample to the generalized Bogomolov–Gieseker inequality for threefolds*,
  arXiv:1602.05055.

---

## Adversarial review notes

**Overall Phase-1 verdict: proceed-with-amendments.**

Substantive changes made in response to the critic:

- **Corrected the abelian-surface math error.** The draft treated abelian surfaces like K3 (a `+r` Mukai
  shim, tagged `[SPECULATIVE]`). Abelian surfaces have trivial tangent bundle, so `√td = (1,0,0)` and
  `v(E) = (r,c₁,ch₂)` with **no shift**; `numerical_wall` on the bare triple is *already* the genuine abelian
  `(s,t)` wall. Rewrote the claim as a new §1c "free win", added node `A1`, changed the table row to
  `[PROVEN, ρ=1]` no shift, and flagged that the K3 shim / `mukai.from_chern` are K3-only and wrong for
  abelian.
- **Upgraded the K3 `+2/d` derivation** from `[SPECULATIVE]` to `[PROVEN for ρ=1]`, citing Bridgeland
  arXiv:math/0307164 and noting `√td_{K3}=(1,0,1)` makes `chern.central_charge` on Mukai coords exact.
- **Fixed the threefold BG citations.** Added Koseki arXiv:1811.03267 (nef tangent bundle; `P¹×S`, `P²×C`,
  `P¹×P¹×C`) as the correct ID, recorded Koseki's distinct weighted-hypersurface / CY double–triple-solid
  paper arXiv:2007.00044 (correcting a prior miscitation to arXiv:1510.04474, which is actually
  Chuang–Lai's *withdrawn* CY/Fano paper, not Koseki), and
  corrected the reference that mislabeled arXiv:1705.04011 as Bernardara–Macrì–Schmidt–Zhao (it is
  Piyaratne; the actual BMSZ Fano-threefold paper is arXiv:1607.08199). Tagged the product-threefold rows
  `[UNVERIFIED]` pending an exact-hypothesis check before any `bg_proven=True` row.
- **Corrected the Levine–Zhang scope** to existence for del Pezzo degree ≥ 3, general-sheaf cohomology for
  degree ≥ 4, anticanonical polarization (title "on del Pezzo surfaces").
- **Reconciled the NS-lattice refactor (`N0`) to high/strategic**, not "medium": it mutates the core frozen
  `ChernChar` and all 42 pinned tests (a v2 architecture decision).
- **Moved `X1` out of the near-term tier** and relabeled it a multi-month, per-surface research effort gated
  on `N0` + `E1` + a non-Chern-level HN/prioritary step; flagged that `E1` has **no** uniform Markov-style
  recursion (Rudakov / Kuleshov–Orlov are case-by-case).
- **Added the threefold tilt-wall caveat**: the cubic `ch₃^β` term means the `(α,β)` wall locus is not
  automatically the surface conic; it must be confirmed against Schmidt arXiv:1509.04608 before pinning
  tests. Kept the `~5`/`~20`-line and effort figures tagged `[SPECULATIVE]`.
- **Adopted the finer BMT framing**: "strong BMT boundary not rigorous on `Bl_p(P³)`" is distinguished from
  "no stability condition exists" (false — `Bl_p(P³)` is Fano and carries stability conditions, BMSZ
  arXiv:1607.08199).
- **Made the prioritization ordering explicit**: the cheap in-architecture wins are all Picard-rank-1; the
  rational surfaces "beyond P²" evokes are the expensive, `N0`+sheaf-gated end.
