# Literature Gap Analysis — `bridgeland_stability`

*Phase-1 research artifact. Scope: map what is **known in the literature** about wall-and-chamber
decompositions and existence of stability conditions for each catalog variety, against what is **implemented**
in this repository, and identify where the computable/open boundary lies for a toolkit whose architecture is
**numerical Chern-character + exact `fractions.Fraction`, Picard rank one, no sheaf-level computation**.*

Every mathematical claim is flagged **[PROVEN]** (a theorem, with citation), **[CONJECTURED]** (stated as a
conjecture/expected, with citation), or **[SPECULATIVE]** (my own inference). Claims I could not pin to a
primary source carry **[UNVERIFIED]**. arXiv IDs were checked against the arXiv abstract page before citing;
the six secondary surface/threefold IDs (1901.04848, 1203.0884, 1107.5304, 1607.04946, 2501.06779, 1811.03267)
were re-verified against their abstracts during this review and are confirmed.

## 0. Three distinctions this document keeps separate

Per the project's correctness culture, I never conflate:

- **(a) computable in principle** — a known algorithm exists on paper;
- **(b) implemented in existing software** — a named package/function does it. *Finding:* for Bridgeland
  walls there is **no standard, maintained library**; published computations are ad-hoc scripts (e.g. Benjamin
  Schmidt's threefold wall computations accompanying arXiv:1509.04608). I flag **(b)** as **[SPECULATIVE]**
  wherever I cannot name a package (this is an unprovable negative — retained as a flag, not a claim), and this
  repo is in fact a rare public (b) for the numerical layer.
- **(c) implementable *here*** — expressible purely at the numerical-Chern-character level (one degree number
  `c = ch1·H` per class, exact `Fraction`), **without** ever constructing a sheaf, Ext group, resolution, or
  the derived category.

The single sharpest architectural constraint: **a class is `(r, c, e)` with one degree number.** This is
faithful only for **Picard-rank-one** varieties. Any variety whose interesting wall structure lives in a
higher-rank Néron–Severi lattice (a product of two non-isogenous elliptic curves, most K3s of Picard rank ≥ 2,
`P^1×P^1` with mixed polarizations, blow-ups) **cannot be faithfully represented** in the current model. This
is the recurring gap below.

---

## 1. Per-variety knowledge table

### 1a. Surfaces

Existence of Bridgeland stability conditions on **every** smooth projective surface is a theorem (Bridgeland
for K3, Arcara–Bertram in general, via tilting) **[PROVEN]** (Bridgeland, Duke 141 (2008), math/0307164;
Arcara–Bertram, arXiv:0708.2247). So the "existence" question is settled for all surface rows; the live
content is the **wall/chamber** structure and **non-emptiness**.

| Surface (catalog) | Wall-and-chamber / non-emptiness — known results | Implemented here | Gap |
|---|---|---|---|
| **P²** (`P2`) | Full wall-and-chamber for `P²[n]` and small-rank/coprime classes; MMP of `P²[n]` via wall-crossing; DLP non-emptiness curve exact **[PROVEN]** (ABCH arXiv:1203.0316; Coskun–Huizenga Gökova survey; Drézet–Le Potier 1985) | DLP curve `δ(μ)`, exceptional bundles (Markov), `actual_walls` certified for `P²[n]` | Higher-rank Brill–Noether cohomology; non-`P²[n]` completeness proof |
| **P¹×P¹** (`P1xP1`) | Non-emptiness algorithm and Bridgeland walls known (Picard rank 2; `F_0` case of Hirzebruch existence) **[PROVEN]** (Coskun–Huizenga arXiv:1907.06739; Rudakov exceptional collections) | Stored as `d=2` rank-1 proxy; wall geometry via `numerical_wall` | **Picard rank 2**: the (1,1) polarization loses the two-parameter NS lattice; mixed/asymmetric walls invisible **[SPECULATIVE]** (architectural inference, not theorem) |
| **K3** (`K3(h2)`) | Walls ↔ Mukai-lattice geometry; finiteness of walls; MMP via wall-crossing; `dim M_σ(v)=v²+2`; wall classification **[PROVEN]** (Bayer–Macri arXiv:1301.6968; Bridgeland Duke 2008; Yoshioka arXiv:1607.04946) | `mukai.py`: pairing, `v²`, `moduli_dim`, `classify_wall` (BM Thm 5.7); rank-1 lattice; Mukai vector `(r, c₁, ch₂+r)` via the `ch₂+r` shift (see note) | Picard rank ≥ 2 K3s; totally-semistable/fake walls detection is search-boxed, not certified |
| **Abelian surface** (`abelian_surface`) | Walls, finiteness, positive-cone/MMP fully described; rank-1 moduli explicit **[PROVEN]** (Yanagida–Yoshioka arXiv:1203.0884; Maciocia–Meachan arXiv:1107.5304; Bayer–Macri arXiv:1301.6968) | `χ(O)=0`, `K=0` recorded; same nested-semicircle `(s,t)` geometry as K3 — but **no K3 Mukai shift** (see note) | No abelian-specific Mukai module; the bare `(r,c,ch₂)` is *already* the Mukai vector; `E₁×E₂` (Picard rank ≥ 2) unrepresentable |
| **Hirzebruch** `F_n` (`hirzebruch`) | **Algorithmic** non-emptiness of semistable sheaves for *any* ample divisor **[PROVEN]** (Coskun–Huizenga arXiv:1907.06739; Brill–Noether arXiv, Nagoya 238 (2020)) | Only `d=H²` stored; wall geometry via `numerical_wall` | The Coskun–Huizenga algorithm needs prioritary-sheaf stacks + HN filtration — **more than Chern arithmetic**; only partially portable |
| **Enriques** (special) | Projective moduli for generic σ; MMP via wall-crossing; birational to Hilbert schemes for odd-rank primitive `v` **[PROVEN]** (Nuer–Yoshioka arXiv:1901.04848; Yoshioka arXiv:1607.04946) | Not in catalog (kind tag only) | `K` is 2-torsion (not numerically trivial, not ample-negative); needs a torsion-aware invariant the model lacks |
| **Bielliptic** (special) | Chern-character classification of stable sheaves + Bridgeland wall-crossing now complete **[PROVEN]** (arXiv:2107.13370, Math. Z. 2024) | Not in catalog | Picard rank up to 4 with torsion; unrepresentable in rank-1 model |

**Mukai-coordinate convention — K3 vs. abelian (project-wide).** For a **K3** surface the Mukai vector is
`v(E) = ch(E)·√td(X) = (r, c₁, ch₂ + r)`, because `√td(K3) = (1,0,1)`; this is the `ch₂ → ch₂ + r` shift in
`mukai.from_chern`, and it raises the wall radius² by `2/d`. **[PROVEN]** (standard; Mukai, *Invent. Math.* 77
(1984); Bayer–Macri arXiv:1301.6968). For an **abelian** surface `√td(X) = (1,0,0)`, so the Mukai vector is the
*bare* `(r, c₁, ch₂)` with **no shift** and **no `2/d` radius correction**. Consequently `mukai.from_chern`
(with its `ch₂ + r`) is **K3-only**; the abelian row's "same semicircle geometry" refers to the shared
`(s,t)`-plane nested-semicircle picture, **not** to the K3 coordinate shift. **[SPECULATIVE]** as an
implementation directive (the Todd-class computation is [PROVEN]; the K3-only scoping of `from_chern` is my
architectural reading, offered to keep the two families from being conflated in any future `(s,t)` wall module).

### 1b. Threefolds

Here "existence" is the deep content: it rests on a **Bogomolov–Gieseker-type inequality for tilt-stable
objects** (the BMT conjecture). The general threefold case is **open** (see §3); the catalog rows are exactly
the varieties where it is a theorem. Note the sharper framing carried through the blow-up row below: existence
of *a* stability condition and rigor of the *naive strong* tilt-BG boundary are **different** questions — the
first can hold while the second fails.

| Threefold (catalog) | Tilt-BG / existence of stability conditions | Wall-crossing known? | Implemented here | Gap |
|---|---|---|---|---|
| **P³** (`P3`, `bg_proven=True`) | BG **[PROVEN]** (Bayer–Macri–Toda arXiv:1103.5010; Macri arXiv:1207.4980) | Explicit first wall-crossings (e.g. twisted cubics) **[PROVEN]** (Schmidt arXiv:1509.04608) | `bmt_Q`, `alpha_crit`, `bg_boundary_curve` (exact `Q`, float boundary) | No `ν_{α,β}` tilt-slope function; no threefold wall-crossing/HN |
| **Quadric Q⊂P⁴** (`QUADRIC3`, `bg_proven=True`) | BG **[PROVEN]** as a Fano 3-fold of Picard rank 1 (C. Li arXiv:1510.04089; BMSZ arXiv:1607.08199) | Explicit wall-crossing examples **[PROVEN]** (Schmidt arXiv:1509.04608) | Same BG boundary machinery | **Citation bug:** code cites `1607.07182` (a *probability* paper) and `1705.04011` (Piyaratne, not BMSZ) — see §5; fix first |
| **Fano Picard-1** (`fano_picard_one`) | BG **[PROVEN]** for *all* Fano 3-folds of Picard rank 1 (C. Li arXiv:1510.04089); stability conditions on *all* Fano 3-folds (BMSZ arXiv:1607.08199) | Case-by-case (index 1/2 Kuznetsov-component work, large literature) | `alpha_crit` boundary | No index/Kuznetsov-component structure (needs derived category) |
| **Abelian 3-fold** (`abelian_threefold`) | BG (strong form) + explicit stability-condition component **[PROVEN]** (Bayer–Macri–Stellari arXiv:1410.1585) | Less developed; Fourier–Mukai symmetry known | BG boundary; `χ(O)=0` | Picard rank ≥ 3 generic; rank-1 proxy only |
| **Quintic CY3** (`QUINTIC`, `bg_proven=True`) | BG + first geometric stability on a strict CY3 **[PROVEN]** (C. Li arXiv:1810.03434); related CY weighted hypersurfaces / double–triple solids (Koseki arXiv:2007.00044) | Essentially open beyond existence | BG boundary (`Q≡0` on `O`) | Wall-crossing on CY3s largely open (§3) |
| **Bl_p P³** (`BLOWUP_P3_POINT`, `bg_proven=False`) | The *stronger/generalized* BMT inequality **FAILS** **[PROVEN]** (Schmidt arXiv:1602.05055); yet a **modified** optimal inequality holds and stability conditions **do exist** (it is Fano) **[PROVEN]** (Piyaratne arXiv:1705.04011; BMSZ arXiv:1607.08199) | — | `bg_proven=False` correctly warns Algorithm 5's naive boundary is not rigorous | **Nuance to preserve:** stability conditions *exist* here; only the *naive strong tilt-BG boundary* is wrong — this is NOT "no stability condition exists" (§5) |

---

## 2. Algorithmic references (explicit recipes translatable to code)

Papers containing procedures/formulas portable to the numerical-Chern-character layer. Each is verified against
its arXiv abstract.

1. **Maciocia — arXiv:1202.4587**, *Computing the Walls Associated to Bridgeland Stability Conditions on
   Projective Surfaces* (Asian J. Math. 18 (2014)). **[PROVEN]** Gives the semicircular-wall structure in the
   `(s,t)` plane and **finiteness/boundedness criteria** for walls on a general surface. *Procedure:* nested
   semicircles; a class `w` gives a genuine wall for `v` only inside the "potential wall" region. **This is the
   theoretical backbone of `walls.numerical_wall` / `actual_walls`.** Directly portable; already partly used.

2. **Arcara–Bertram–Coskun–Huizenga (ABCH) — arXiv:1203.0316**, *The MMP for the Hilbert Scheme of Points on P²
   and Bridgeland Stability* (Adv. Math. 235 (2013)). **[PROVEN]** For `P²[n]` gives the **complete finite wall
   set** with explicit destabilizers; the Gieseker (outermost) wall is center `−(2n+1)/2`, radius `(2n−1)/2`.
   *Formula:* the center/radius derivation this repo re-implements exactly. Fully portable (done).

3. **Coskun–Huizenga — Gökova survey**, *The birational geometry of the moduli spaces of sheaves on P²*
   (Gökova 2014); expanded in **Huizenga survey arXiv:1606.02775**. **[PROVEN]** §4: the **DLP curve** as a
   fractal cusp envelope `δ(μ)=sup(P(−|μ−α|)−Δ_α)`, control interval `x_α=(3−√(5+8Δ_α))/2`; §5–6: the
   **rank-reduction and Bogomolov-on-both necessary conditions** for an *actual* wall. Portable (this is exactly
   `dlp.delta` and `actual_walls`).

4. **Coskun–Huizenga–Woolf — arXiv:1401.1613**, *The effective cone of the moduli space of sheaves on the
   plane* (JEMS 19 (2017)). **[PROVEN]** Gives, via the Beilinson spectral sequence, the resolution of the
   general sheaf that determines the effective cone. *Portability:* **partial** — the numerical shadow (which
   exceptional bundle controls the extremal ray) is Chern-level, but the resolution itself is sheaf-level.

5. **Coskun–Huizenga — arXiv:1907.06739**, *Existence of semistable sheaves on Hirzebruch surfaces* (Adv. Math.
   381 (2021)). **[PROVEN]** An **effective non-emptiness algorithm** for any ample divisor, via prioritary
   stacks and the HN filtration of the general sheaf. *Portability:* the decision output is combinatorial but
   the intermediate HN/prioritary steps exceed pure `(r,c,e)` arithmetic — a genuine extension target, and one
   that is **not** near-term (see the §4 framing note: it is sheaf-theoretic and, for `P¹×P¹`/`F_n`, also gated
   behind the NS-lattice refactor).

6. **Bayer–Macri — arXiv:1301.6968**, *MMP for moduli of sheaves on K3s via wall-crossing* (Invent. Math. 198
   (2014)). **[PROVEN]** Thm 2.15: `dim M_σ(v)=⟨v,v⟩+2`; **Thm 5.7**: the divisorial (Brill–Noether /
   Hilbert–Chow / Li–Gieseker–Uhlenbeck) vs. flopping vs. fake classification via spherical (`s²=−2`) and
   isotropic (`w²=0`) classes in the rank-2 hyperbolic sub-lattice. Fully portable at the Mukai-lattice level
   (implemented in `mukai.classify_wall`).

7. **Bayer–Macri–Toda — arXiv:1103.5010**, *Bridgeland Stability conditions on threefolds I* (J. Alg. Geom. 23
   (2014)). **[PROVEN for the listed varieties / CONJECTURED in general]** The tilt-BG quadratic form
   `Q = 4 ch₂ᵇ² − 6 ch₁ᵇ ch₃ᵇ ≥ 0` and the critical-radius boundary. Portable (implemented in `threefold.py`).

8. **C. Li — arXiv:1510.04089** (Fano Picard-1, JEMS 21 (2019)) and **arXiv:1810.03434** (quintic, Invent.
   Math. 218 (2019)). **[PROVEN]** Provide the *stronger* explicit inequalities that make `bg_proven=True`
   rigorous for those catalog rows. The inequalities themselves are Chern-level (portable as sharper BG checks);
   the *proofs* use Clifford-type bounds on curves (not portable).

9. **Bernardara–Macri–Schmidt–Zhao — arXiv:1607.08199** (Épijournal Géom. Alg. 1 (2017)) and **Piyaratne —
   arXiv:1705.04011**. **[PROVEN]** Existence of stability conditions on **all** Fano 3-folds, incl. an
   optimal modified inequality on `Bl_p P³`. These are the *correct* IDs for the quadric/blow-up rows (§5).

10. **Schmidt — arXiv:1509.04608**, *Bridgeland Stability on Threefolds — Some Wall Crossings*. **[PROVEN]**
    The first explicit threefold wall computations (P³, quadric); the closest thing to a reference
    *implementation* of threefold wall-crossing. A model for a future `ν_{α,β}` + wall module here.

11. **Koseki — arXiv:2007.00044**, *Stability conditions on Calabi–Yau double/triple solids*, and
    **Koseki — arXiv:1811.03267**, *Stability conditions on threefolds with nef tangent bundles*. **[PROVEN]**
    The first proves BG on CY double/triple solids (weighted hypersurfaces in `P(1,1,1,1,2)`/`P(1,1,1,1,4)`);
    the second proves BG (hence existence of
    Bridgeland stability conditions) for **threefolds with nef tangent bundle** — a class that includes products
    whose factors individually have nef tangent bundle (e.g. `P¹×S`, `P²×C`, `P¹×P¹×C` with `S` abelian /
    `C` elliptic, since a product of nef tangent bundles is nef). Both are candidates for new catalog rows with
    `bg_proven=True`, **provided the per-variety hypotheses are checked first** (Picard rank of such a product is
    ≥ 2, so those rows would be `bg_proven` markers only, not faithful wall models — see §3, Néron–Severi rank
    > 1).

---

## 3. Open problems that block generalization

- **The BMT conjecture for general threefolds is OPEN, and its naive strong form is FALSE.** The
  Bogomolov–Gieseker-type inequality for tilt-stable objects (arXiv:1103.5010) is **[CONJECTURED]** in general;
  Schmidt's `Bl_p P³` counterexample **[PROVEN]** (arXiv:1602.05055) shows the *generalized/strong* form fails,
  so there is **no known uniform tilt-BG boundary** for an arbitrary polarized threefold. *Consequence for the
  code:* Algorithm 5's boundary is rigorous **only** on the explicitly-flagged `bg_proven=True` rows; extending
  the threefold catalog requires a per-variety theorem, not a general algorithm. (Note the distinction the
  blow-up row makes: the failure of the strong boundary is *not* the non-existence of stability conditions —
  those exist on `Bl_p P³` because it is Fano, BMSZ arXiv:1607.08199.)

- **Tilt-stability walls are not automatically Bridgeland walls.** The `(s,t)` / `(β,α)` semicircle a
  numerical-wall computation returns is the **tilt-stability** wall. Upgrading it to a genuine **Bridgeland**
  wall — the object of an actual wall-crossing statement on the threefold — requires the BG (BMT) inequality to
  hold, i.e. the `bg_proven=True` gate; the numerics **do not certify** that upgrade on their own. Any future
  `ν_{α,β}` enumerator therefore produces certified *Bridgeland* walls only on the `bg_proven=True` rows, and
  tilt-only candidates elsewhere. **[PROVEN]** (this is exactly the content of the arXiv:1103.5010 / BMSZ
  framework: the second tilt and the Bridgeland heart depend on the BG inequality).

- **Certified completeness of surface wall sets for arbitrary `v` is not a closed theorem.** ABCH/Coskun–Huizenga
  give the finite actual-wall set for `P²[n]` and coprime/small-rank classes **[PROVEN]**; for general `(r,c,e)`
  on general surfaces the finiteness is known (Maciocia arXiv:1202.4587) **[PROVEN]** but a *provably complete*
  enumeration with a stopping certificate is not published in closed form. The repo's `actual_walls_complete`
  gives an **empirical** doubling certificate — honest, but **[SPECULATIVE]** as a completeness proof for
  arbitrary `v`.

- **Totally-semistable / "fake" walls on K3s.** Deciding whether a numerical wall is an *actual* wall vs. a
  totally-semistable locus requires the full Bayer–Macri stratification; `mukai.classify_wall` searches a finite
  `|a|,|b|≤search` box, so a "fake-or-none" verdict is **[SPECULATIVE]** (not certified) — the theorem
  (arXiv:1301.6968) is exact, the *bounded search* is not.

- **Wall-crossing on Calabi–Yau threefolds is largely open.** Even with existence proven for the quintic
  (arXiv:1810.03434) and some CY weighted hypersurfaces / double–triple solids (arXiv:2007.00044), an explicit wall-and-chamber
  decomposition for CY3 Chern characters is **[CONJECTURED]/open** — a frontier, not a porting target.

- **Higher-rank / non-Bogomolov surface bounds.** `bg_check` asserts only `Δ≥0`. Sharper existence bounds
  (Hoppe–Gieseker on abelian/K3, the DLP curve on P², the Coskun–Huizenga Hirzebruch algorithm) are the *actual*
  non-emptiness criteria; the general sharp bound for an arbitrary surface is **open** **[CONJECTURED]** (it is
  the subject of ongoing Brill–Noether-for-sheaves work).

- **Néron–Severi rank > 1.** There is no mathematical obstruction to wall theory in higher Picard rank — the
  obstruction is purely architectural. Representing `P¹×P¹` mixed polarizations, `E₁×E₂`, or general K3s
  requires a full NS lattice vector, not one degree number. **[SPECULATIVE]** (my architectural assessment).

---

## 4. Prioritization — research value per unit effort (rank-1 / exact-`Fraction` constraint respected)

Ranked by *(research value) / (effort achievable **without** sheaf-level computation)*.

**Framing note — "beyond P²" is not one target but two.** The genuinely cheap, *in-architecture* wins are all
**Picard-rank-one**: K3/abelian `(s,t)` wall enumeration, Picard-1 threefold tilt-stability walls, and sharper
higher-rank surface BG bounds. The **rational** surfaces that "beyond P²" most evokes — `P¹×P¹`, Hirzebruch
`F_n`, del Pezzo — are *not* near-term: they are gated behind BOTH the expensive Néron–Severi lattice refactor
(item 7) AND a partly **sheaf-theoretic** non-emptiness algorithm (prioritary stacks + HN filtration,
arXiv:1907.06739), so their wall/non-emptiness content is not expressible in pure `(r,c,e)` arithmetic. The
ranking below respects this split so that a roadmap does not mis-prioritize rational surfaces as near-term.
**[SPECULATIVE]** (architectural assessment).

1. **Fix the miscited threefold references (§5).** Effort: minutes. Value: high — this is a
   correctness-culture project, and this is the **only confirmed defect in shipped code**: two catalog rows cite
   a probability paper and a mis-attributed author. **Do this first — it is the roadmap's first task.**

2. **A general tilt-slope `ν_{α,β}` function + threefold tilt-stability-wall enumerator (Picard-1).** Effort:
   medium. Value: high. The `(β,α)` semicircle geometry is the exact 3-fold analogue of the surface
   `numerical_wall` already present; Schmidt arXiv:1509.04608 gives the recipe, and it stays entirely at the
   Chern-character level. **Caveat (see §3):** the numerics return the **tilt-stability** wall; upgrading it to a
   genuine **Bridgeland** wall requires the `bg_proven=True` gate and is **not certified by the numerics alone**,
   so the enumerator must label its output tilt-only vs. Bridgeland-certified accordingly. **Highest
   value/effort after the citation fix.**

3. **K3/abelian wall enumerator in the `(s,t)` plane driven by the Mukai lattice.** Effort: low–medium. Value:
   high. `mukai.classify_wall` already encodes BM Thm 5.7; adding the actual semicircle centers/radii for a
   Mukai vector `v` closes the loop between the numerical wall (`walls.py`) and the lattice classification. Fully
   Chern/lattice-level. **Mind the K3-vs-abelian convention (§1a note):** K3 uses `v=(r,c₁,ch₂+r)` (radius²
   shifted by `2/d`), abelian uses the bare `(r,c₁,ch₂)` with no shift — `mukai.from_chern` is K3-only.
   **[PROVEN]** foundations (arXiv:1301.6968).

4. **Sharper surface non-emptiness beyond `Δ≥0`: port the Hirzebruch decision output (arXiv:1907.06739) and
   del Pezzo bounds.** Effort: medium–high (the algorithm's HN step is not pure `(r,c,e)`; and for `P¹×P¹`/`F_n`
   the *faithful* model is additionally gated behind item 7 — see the framing note). Value: medium-high — would
   make `hirzebruch(n)`/`P¹×P¹` first-class. Partial portability; flag honestly. **Not near-term.**

5. **Widen exact completeness certificates.** Effort: low. Value: medium. Replace the empirical doubling in
   `actual_walls_complete` with the Maciocia boundedness bound (arXiv:1202.4587) where a closed radius bound
   exists, upgrading `[SPECULATIVE]`→`[PROVEN]` for the classes it covers.

6. **Add Enriques/bielliptic catalog rows.** Effort: low to *record*, high to *compute*. Value: medium. The
   existence/wall theory is now proven (arXiv:1607.04946, arXiv:2107.13370, arXiv:1901.04848), but the torsion
   canonical class and Picard rank ≥ 2 make faithful representation hard. Record as catalog entries with
   `bg`-style flags; defer computation. **[SPECULATIVE]** on feasibility of full computation here.

7. **Néron–Severi rank-`ρ` refactor (`c` → NS lattice vector).** Effort: **very high** — this mutates the core
   frozen `ChernChar` (scalar `c` → NS lattice vector), `walls.py`, `mukai.py`, and every one of the 42 pinned
   tests, and it contradicts the current "one degree number" invariant. Value: high but foundational. This is
   the one change that unlocks products of curves, mixed polarizations, and higher-rank K3s. **Strategic, not
   incremental — a deliberate v2 architecture decision, not a roadmap increment.** (Do not mistake this for a
   "medium refactor": it is the gate under every rational-surface item above.)

---

## 5. Correctness note — miscited arXiv IDs in `varieties.py` (the roadmap's first task)

Because citation accuracy is this project's core value, two errors found while verifying — this is the **only
confirmed defect in shipped code** and should be the first roadmap task:

- `QUADRIC3.references = ["arXiv:1607.07182", "arXiv:1705.04011"]` — **both wrong.**
  - **arXiv:1607.07182** is *"Interlacing Diffusions"* (Assiotis–O'Connell–Warren, math.PR) — an unrelated
    **probability** paper. **[PROVEN]** (verified against the arXiv abstract).
  - **arXiv:1705.04011** is Piyaratne, *"Stability conditions, Bogomolov–Gieseker type inequalities and Fano
    3-folds"* — **not** Bernardara–Macri–Schmidt–Zhao. **[PROVEN]** (verified).
  - **Concrete fix:** set
    `QUADRIC3.references = ["arXiv:1510.04089", "arXiv:1607.08199", "arXiv:1509.04608"]`
    (C. Li Fano-Picard-1 BG theorem covering the quadric; BMSZ Fano-3-fold stability conditions; Schmidt's
    quadric wall-crossing). Update the `note` string accordingly so it no longer implies "1705.04011 = BMSZ".

- `BLOWUP_P3_POINT.references = ["arXiv:1602.05055", "arXiv:1705.04011"]` — the **IDs are correct** and should
  be **kept**: 1602.05055 (Schmidt counterexample) and 1705.04011 (Piyaratne's optimal *modified* inequality on
  `Bl_p P³`) are both legitimately relevant here.
  - **Concrete fix:** keep the two IDs, but correct the attribution in the code comment — **1705.04011 is
    Piyaratne**, not Bernardara–Macri–Schmidt–Zhao.
  - **Also add to the note:** the blow-up **is Fano**, so stability conditions **do exist** on it
    (arXiv:1607.08199) — only the *naive strong* tilt-BG boundary fails. The note should distinguish "strong
    BMT boundary not rigorous" (true, and what `bg_proven=False` encodes) from "no stability condition exists"
    (false).

*(These are documentation/citation fixes; no computed value or pinned test changes.)*

---

## References

- Bridgeland, *Stability conditions on K3 surfaces*, Duke Math. J. 141 (2008) 241–291 — arXiv:math/0307164.
- Arcara–Bertram, *Bridgeland-stable moduli spaces for K-trivial surfaces* — arXiv:0708.2247 (surface existence).
- Maciocia, *Computing the walls associated to Bridgeland stability conditions on projective surfaces*, Asian J.
  Math. 18 (2014) — arXiv:1202.4587.
- Arcara–Bertram–Coskun–Huizenga (ABCH), *The MMP for the Hilbert scheme of points on P² and Bridgeland
  stability*, Adv. Math. 235 (2013) 580–626 — arXiv:1203.0316.
- Coskun–Huizenga, *The birational geometry of the moduli spaces of sheaves on P²*, Gökova 2014; survey:
  Huizenga, *Birational geometry of moduli spaces of sheaves and Bridgeland stability* — arXiv:1606.02775.
- Coskun–Huizenga–Woolf, *The effective cone of the moduli space of sheaves on the plane*, JEMS 19 (2017) —
  arXiv:1401.1613.
- Coskun–Huizenga, *Existence of semistable sheaves on Hirzebruch surfaces*, Adv. Math. 381 (2021) —
  arXiv:1907.06739.
- Bayer–Macri, *MMP for moduli of sheaves on K3s via wall-crossing…*, Invent. Math. 198 (2014) —
  arXiv:1301.6968.
- Mukai, *Symplectic structure of the moduli space of sheaves on an abelian or K3 surface*, Invent. Math. 77
  (1984) 101–116 (Mukai lattice / `√td` convention).
- Yanagida–Yoshioka, *Bridgeland's stabilities on abelian surfaces*, Math. Z. — arXiv:1203.0884 (verified).
- Maciocia–Meachan, *Rank one Bridgeland stable moduli spaces on a principally polarized abelian surface*, IMRN
  2013 — arXiv:1107.5304 (verified).
- Yoshioka, *Some remarks on Bridgeland stability conditions on K3 and Enriques surfaces* — arXiv:1607.04946
  (verified).
- Nuer–Yoshioka, *MMP via wall-crossing for moduli spaces of stable sheaves on an Enriques surface* —
  arXiv:1901.04848 (verified).
- *Stable sheaves on bielliptic surfaces: from the classical to the modern*, Math. Z. (2024) — arXiv:2107.13370.
- Bayer–Macri–Toda (BMT), *Bridgeland stability conditions on threefolds I: Bogomolov–Gieseker type
  inequalities*, J. Alg. Geom. 23 (2014) 117–163 — arXiv:1103.5010.
- Macri, *A generalized Bogomolov–Gieseker inequality for the three-dimensional projective space* —
  arXiv:1207.4980.
- Schmidt, *Bridgeland stability on threefolds — some wall crossings* — arXiv:1509.04608.
- Schmidt, *Counterexample to the generalized Bogomolov–Gieseker inequality for threefolds* — arXiv:1602.05055.
- C. Li, *Stability conditions on Fano threefolds of Picard number one*, JEMS 21 (2019) 709–726 —
  arXiv:1510.04089.
- C. Li, *On stability conditions for the quintic threefold*, Invent. Math. 218 (2019) 301–340 —
  arXiv:1810.03434.
- Bayer–Macri–Stellari, *The space of stability conditions on abelian threefolds, and on some Calabi–Yau
  threefolds*, Invent. Math. 206 (2016) 869–933 — arXiv:1410.1585.
- Bernardara–Macri–Schmidt–Zhao, *Bridgeland stability conditions on Fano threefolds*, Épijournal Géom. Alg. 1
  (2017) — arXiv:1607.08199.
- Piyaratne, *Stability conditions, Bogomolov–Gieseker type inequalities and Fano 3-folds* — arXiv:1705.04011.
- Koseki, *Stability conditions on Calabi–Yau double/triple solids* — arXiv:2007.00044.
- Koseki, *Stability conditions on threefolds with nef tangent bundles* — arXiv:1811.03267 (verified).
- Drézet–Le Potier, *Fibrés stables et fibrés exceptionnels sur P₂*, Ann. Sci. ENS 18 (1985) 193–243.
- Rudakov, *The Markov numbers and exceptional bundles on P²* (1988); Veselov, *Markov fractions and the slopes
  of the exceptional bundles on P²* — arXiv:2501.06779 (verified against the arXiv abstract this review).

---

## Adversarial review notes

**Overall Phase-1 verdict: proceed-with-amendments.**

Substantive changes made in response to the adversarial critic:

- **Made the citation-bug fix concrete and elevated it to the roadmap's first task.** §5 now gives the exact
  edits: `QUADRIC3.references → ["arXiv:1510.04089", "arXiv:1607.08199", "arXiv:1509.04608"]`, and
  `BLOWUP_P3_POINT.references` kept as `["arXiv:1602.05055", "arXiv:1705.04011"]` with the code comment fixed to
  attribute 1705.04011 to **Piyaratne** (not BMSZ). Flagged as the only confirmed defect in shipped code.
- **Reworded §4 item 2** from "actual wall" to **"tilt-stability wall"**, added a cross-reference to §3, and
  surfaced the tacit assumption that a tilt wall equals a Bridgeland wall — it does not; the upgrade needs the
  `bg_proven` gate and is not certified by numerics alone. Added a matching new §3 bullet.
- **Verified all six secondary IDs** (1901.04848, 1203.0884, 1107.5304, 1607.04946, 2501.06779, 1811.03267)
  against their arXiv abstracts; all confirmed, so no `[UNVERIFIED]` tags were needed. Updated the Veselov,
  Yanagida–Yoshioka, Maciocia–Meachan, Yoshioka and Nuer–Yoshioka reference lines to mark them verified.
- **Unified the K3-vs-abelian Mukai-coordinate convention project-wide** (new §1a note): K3 shifts
  `ch₂ → ch₂ + r` (radius² `+= 2/d`), abelian does **not** (bare `(r,c,ch₂)` is already the Mukai vector), and
  `mukai.from_chern` is documented as **K3-only**. Corrected the abelian table row accordingly.
- **Made the "beyond P²" prioritization ordering explicit** (new §4 framing note): the cheap in-architecture
  wins are the Picard-rank-1 varieties; the rational surfaces (`P¹×P¹`, `F_n`, del Pezzo) are gated behind BOTH
  the NS-lattice refactor and a partly sheaf-theoretic algorithm, so they are not near-term.
- **Reconciled the NS-lattice-refactor effort to "very high / strategic v2 decision"** (§4 item 7), noting it
  mutates the core frozen `ChernChar`, `walls.py`, `mukai.py`, and all 42 pinned tests — the gate under every
  rational-surface item.
- **Sharpened the existence-vs-BG framing** on `Bl_p P³`: stability conditions **exist** (it is Fano); only the
  naive strong tilt-BG boundary fails. Kept `bg_proven=False` as encoding the boundary-rigor claim, not a
  non-existence claim.
- **Added Koseki arXiv:1811.03267** (nef-tangent-bundle threefolds, incl. products with abelian/elliptic
  factors) as an additional `bg_proven` candidate source, with the hypothesis caveat that such products have
  Picard rank ≥ 2 (marker rows only, not faithful wall models).
- **Retained** the `[SPECULATIVE]` tag on "no standard, maintained Bridgeland-walls package" (unprovable
  negative) and on the `P¹×P¹` "mixed/asymmetric walls invisible" claim (architectural inference, not theorem).
