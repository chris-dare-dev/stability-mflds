# Computational Homological Algebra Survey for `bridgeland_stability`

*Discovery survey. Scope: what deep homological-algebra / derived-category
software actually exists, how (and whether) to reach it from a zero-dependency
Python core, what "computing Ext" means computationally, what is genuinely
uncomputable today, and the right backend strategy for this project.*

Every feasibility/mathematical claim is tagged `[PROVEN]` (theorem in the
literature), `[CONJECTURED]` (stated as conjecture/expected), or `[SPECULATIVE]`
(my own inference); software facts I could not pin to a paper/manual are tagged
`[UNVERIFIED]`. Throughout I keep three things separate: **(a)** computable in
principle by a known algorithm, **(b)** already implemented in named software,
**(c)** implementable in *this* codebase's numerical-Chern-character / exact-`Fraction`
architecture **without any sheaf-level computation**. The current code
(`chern.py`, `walls.py`, `exceptional.py`, `dlp.py`, `bg_check.py`,
`threefold.py`, `mukai.py`, `varieties.py`) lives entirely at level **(c)**: it
never constructs a sheaf, an Ext group, a resolution, or an object of the derived
category. That is the context against which everything below is judged.

---

## 1. Existing software landscape

### Macaulay2 (M2)
M2 is the center of gravity for computational homological algebra in algebraic
geometry. Verified against the official "packages provided with Macaulay2" list:

- **`Varieties`** — the modern home of coherent-sheaf computation on affine/projective
  varieties. Computes `HH^i(CoherentSheaf)` (sheaf cohomology), Euler
  characteristics, and crucially **`Ext^i(F,G)` of coherent sheaves on a
  projective scheme**. This is genuine derived-functor computation, via the
  graded-module Ext algorithm of Greg Smith `[PROVEN]` (arXiv:math/9807170).
- **`Complexes`** — the current chain-complex / homological-algebra engine
  ("Towards computing in the derived category"). Functorial cones, resolutions,
  quasi-isomorphisms, `Ext`, `Hom` of complexes. It is the intended successor to
  the old `ChainComplex` machinery for derived-category-style work.
- **`BGG`** — the Bernstein–Gel'fand–Gel'fand correspondence: sheaf cohomology on
  projective space via free resolutions over the exterior algebra (Tate
  resolutions), `cohomologyTable`, `tateResolution`. Extended by the published
  **`MultigradedBGG`** package (arXiv:2402.12293) to the multigraded/toric case
  and differential modules.
- **`TateOnProducts`** — cohomology of sheaves on **products of projective
  spaces**, `cohomologyHashTable`, **`beilinsonMonad`** (the Beilinson spectral
  sequence / monad), `directImageComplex` (derived pushforward) (arXiv:1905.10230).
- **`NormalToricVarieties`** — `HH^i(NormalToricVariety, CoherentSheaf)`: sheaf
  cohomology on normal toric varieties, reducing to Cox-ring module cohomology
  (the local-cohomology-with-monomial-supports reduction of Eisenbud–Mustață–Stillman,
  arXiv:math/0001159) `[PROVEN]`; companion `ToricHigherDirectImages` /
  higher-direct-image work (arXiv:2505.22835).

**What M2 does NOT have as a distributed package:** there is **no** package named
`DerivedCategories.m2` and **no** package named `ExceptionalCollections.m2` in the
standard distribution (verified against the packages-provided list). Derived-category
"objects" are represented indirectly (as complexes, or via a tilting/exceptional
description), not as first-class citizens. The closest thing to *distributed*
exceptional-collection / mutation functionality is now the **Macaulay2 package
accompanying arXiv:2509.25103** (Brown–Dey–Li–Sayrafi, "Computing global Ext for
complexes"): it computes `Ext` of bounded complexes and **mutations of exceptional
collections**, and checks whether bounded complexes form an exceptional collection.
That partially softens the blunt "no ExceptionalCollections package" claim — a
usable mutation tool now ships with a paper, even if it is not a general
derived-category engine. (Its abstract does **not** claim spherical twists; any
spherical-twist capability is `[UNVERIFIED]`.) The **`QuiversToricVarieties`**
package (Prabhu-Naik, arXiv:1501.05861) is a *database* of full strong exceptional
collections of line bundles for smooth toric Fano varieties of dimension ≤ 4 (plus
functions to build a quiver of sections and test the strong-exceptional property)
— useful, but a database, not a general-purpose derived-category toolkit.
Maturity: `Varieties`, `BGG`, `TateOnProducts`, `NormalToricVarieties` are mature,
core-distributed, and battle-tested; `Complexes` is stable but still evolving.

### SageMath
Sage has a large **abstract category framework** (`Category`, `Modules`,
`HomCategory`, …) — but this is category-theory *scaffolding for its algebra
objects*, **not** derived categories of coherent sheaves. For geometry it has real
but shallow support: toric varieties (`ToricVariety`), toric divisors with
`ToricDivisor.cohomology()` (cohomology of line bundles / toric divisor sheaves via
a combinatorial chamber algorithm), Chow groups, and a limited toric-sheaf
constructor. General coherent-sheaf `Ext`, free resolutions, and sheaf cohomology
of arbitrary coherent sheaves are **absent natively**; module-level `Ext`/resolutions
are delegated to **Singular** under the hood. Net: Sage is good for toric line-bundle
cohomology and commutative-algebra module computations, but there is **no**
derived-category-of-`Coh(X)` machinery `[verified against Sage schemes/toric docs]`.

### OSCAR (Julia)
OSCAR unifies GAP + Singular + polymake + Nemo/AbstractAlgebra in Julia. Strong,
actively growing algebraic-geometry stack: toric varieties
(`projective_space`, `hirzebruch_surface`, `del_pezzo_surface`, …), coherent sheaves
on toric varieties as graded modules over the Cox ring, sheaf/line-bundle cohomology,
and general homological algebra on modules (`free_resolution`, `hom`, `ext`) backed
by Singular (arXiv:2303.08110, "Toric Geometry in OSCAR"). It does **not** provide a
turnkey derived-category-of-sheaves or Bridgeland-stability layer. Maturity: rapidly
maturing; the toric and commutative-algebra layers are solid, the higher
derived-category layer is research-stage.

### QPA + CAP / homalg (GAP)
- **QPA** ("Quivers and Path Algebras") — finite-dimensional quotients of path
  algebras, quiver representations, module categories, **`Ext`**, projective
  resolutions, almost-split sequences, tilting. QPA2 is reimplemented on CAP.
- **CAP** ("Categories, Algorithms, Programming") + **homalg** — constructive
  category theory: build abelian/additive/triangulated categories and functors
  generically; homalg does homological algebra (resolutions, derived functors,
  spectral sequences) over arbitrary computable rings.

These operate at the level of **modules over finite-dimensional algebras / quivers**,
not sheaves directly. Their relevance to *this* project is via **tilting**: when
`D^b(Coh X)` admits a full strong exceptional collection, it is equivalent to
`D^b(mod A)` for an explicit quiver-with-relations algebra `A`. For `P^n` this is
Beilinson's theorem `[PROVEN]` (`⟨O, O(1), …, O(n)⟩`, endomorphism algebra = a
Beilinson quiver). So `Ext` in `D^b(P^2)` *can* be routed through QPA — a real, if
indirect, bridge.

### Others
- **Singular** — the commutative-algebra backbone (Gröbner bases, free resolutions,
  module `Ext`). Its `sheafcoh.lib` computes sheaf cohomology on `P^n` (`sheafCoh`,
  `sheafCohBGG`) via the BGG/Tate method.
- **Magma** (proprietary/paid) — projective resolutions and `Ext` over commutative
  *and* non-commutative algebras and exterior algebras (cached, incremental
  resolutions), plus some sheaf-cohomology tooling. Strong but closed-source and
  license-gated.
- **polymake** — computes **cellular** (co)sheaf cohomology on polyhedral complexes
  (the applied-topology sense), **not** coherent-sheaf cohomology in the
  algebraic-geometry sense. Easy to conflate; they are different objects.

**Landscape summary.** Deep, *sheaf-level* homological algebra for `Coh(X)` is
essentially **Macaulay2 (+Singular) and OSCAR**; the derived-category /
exceptional-collection layer is thin everywhere and mostly research code (the
strongest distributed piece being the arXiv:2509.25103 M2 package for `Ext` of
complexes and mutations); QPA/CAP give the quiver side of any tilting equivalence.
No system ships a Bridgeland wall-crossing or stability-condition engine. There is,
in particular, **no general derived-category-of-`Coh(X)` engine in any system** —
a fact that bounds what any of the strategies below can deliver.

---

## 2. Python-interop strategies

The core must stay **zero-runtime-dependency**; any CAS is an *optional*, lazily
imported (or subprocess) bridge, mirroring how `viz/` gates matplotlib/plotly.

- **`juliacall` / PythonCall.jl → OSCAR.** `pip install juliacall`, then
  `jl.seval("using Oscar")`. This is the actively maintained successor to the older
  **PyJulia** (`pip install julia`), which the JuliaPy org now points away from.
  Friction: adds a pip dependency *and* a Julia toolchain + a multi-GB OSCAR install
  (first `using Oscar` precompiles for minutes). **Windows caveat (this host is
  Windows 11):** OSCAR depends on Singular/GAP/polymake, which are Unix-oriented;
  native-Windows OSCAR support has historically been limited and it is typically run
  under **WSL** `[UNVERIFIED — current native-Windows status may have improved and
  needs a fresh check; the "strongest optional extra" ranking below is contingent
  on this]`. Serialization: OSCAR rationals are GMP `QQFieldElem`; `Fraction ↔ QQ`
  round-trips **losslessly** as `"p/q"` text — no float ever appears.
- **Macaulay2 via subprocess / `pexpect`.** M2 has no official Python binding.
  The route that adds no *pip* dependency is to spawn `M2 --script` (or a `pexpect`
  session) and exchange **text**: send a generated `.m2` script, parse the printed
  result. M2 rationals are exact `QQ`; `Fraction("−7/6") ↔ -7/6` is a trivial,
  lossless string round-trip. Be honest about the cost, though: this is "zero
  *Python packages*", not zero dependency — it still requires an **external M2
  install** the user must provision, plus a **hand-maintained text-protocol parser**
  and process-lifecycle handling, which is a real (non-pip) maintenance burden.
  Windows: M2 runs natively via its Windows build or under WSL.
- **Sage as a Python library (`sage -python`, `from sage.all import *`).** Only works
  inside Sage's own bundled CPython; you cannot `pip install sage` into an arbitrary
  interpreter. Sage dropped native Windows support — on Windows it is **WSL/Cygwin
  only**. Heavyweight and awkward for an "optional bridge"; not recommended here.
- **`pexpect` drivers generically.** Robust for any REPL-style CAS (M2, Singular,
  Sage), text-only, no compiled bindings. Main cost is parser fragility and process
  lifecycle management.

**Verdict for a zero-dep core that only *optionally* reaches a CAS:** a **Macaulay2
subprocess bridge** (text I/O, lossless `Fraction↔QQ`, no *pip* dependency) is the
most architecturally aligned — with the caveat above that "no pip dependency" is not
"no dependency". **`juliacall`→OSCAR** is the strongest *optional extra*
(`pip install ".[cas]"`) when toric coherent-sheaf work is wanted, accepting the
Julia/OSCAR install and the Windows-via-WSL friction (which itself is `[UNVERIFIED]`
and should be re-checked before committing to that ranking). Sage is the weakest fit
on this platform. All exact-arithmetic round-tripping is a non-issue: every one of
these CAS uses GMP rationals, so `fractions.Fraction` serialized as `p/q` is exact
in both directions `[PROVEN — trivially, both sides are exact ℚ]`.

---

## 3. What "computing `Ext^i(E,F)`" means computationally, and the state of the art

`Ext^i(E,F) = Hom_{D^b(X)}(E, F[i])`. Three standard reductions, with where each is
implemented:

1. **Via sheaf cohomology of `E^∨ ⊗ F`.** When `E` is locally free,
   `Ext^i(E,F) ≅ H^i(X, E^∨ ⊗ F)` `[PROVEN]`; in general the local-to-global
   spectral sequence `H^p(Ext^q) ⇒ Ext^{p+q}` intervenes. On a projective scheme,
   `H^i(X, \mathcal F)` and global `Ext` are computed as **truncations of `Ext` of
   the associated graded modules over the homogeneous coordinate ring** — the
   algorithm of Greg Smith `[PROVEN]` (arXiv:math/9807170). On a toric variety the
   ring is the **Cox ring** and cohomology reduces to a **local-cohomology /
   Cox-module** computation `[PROVEN]` (Eisenbud–Mustață–Stillman, "Cohomology on
   toric varieties and local cohomology with monomial supports", arXiv:math/0001159).
   **Implemented:** M2 `Varieties.Ext^i(F,G)` and `NormalToricVarieties.HH^i`;
   OSCAR; Singular (module `Ext`); Magma. This is category **(b)** — it *constructs*
   sheaves/modules and takes free resolutions, which is exactly what the current
   `bridgeland_stability` core deliberately never does.
2. **Via the Beilinson spectral sequence / monad.** Resolve the diagonal / build the
   Beilinson monad to express a sheaf on `P^n` (or products) through the exceptional
   collection `{O(−i)}` (equivalently, the Tate resolution over the exterior algebra
   — the BGG correspondence), turning cohomology into linear algebra. **Implemented:**
   M2 `TateOnProducts.beilinsonMonad`, `BGG` / `MultigradedBGG`, Singular's
   BGG-based sheaf cohomology `[PROVEN + implemented]`.
3. **Via an exceptional collection + mutations.** If `D^b(X) = ⟨E_0, …, E_n⟩`, then
   `Ext` between arbitrary objects becomes computation in the endomorphism DG-algebra
   / bound quiver algebra; `Ext^•(⊕E_i, ⊕E_i)` is the quiver-with-relations data.
   **Implemented:** partially — quiver side via **QPA**; the geometry-to-quiver and
   mutation side now has a distributed M2 tool (`Ext` of bounded complexes and
   **mutations of exceptional collections**, Brown–Dey–Li–Sayrafi, arXiv:2509.25103)
   plus databases (`QuiversToricVarieties`, arXiv:1501.05861), but no single unified
   general-purpose package.

**Bottom line for this project:** at level **(a)** and **(b)**, `Ext^i` between
concrete sheaves is a *solved, implemented* problem on projective and toric
varieties (M2/OSCAR). At level **(c)** — the numerical-Chern-character architecture
here — only the **Euler pairing** `χ(E,F) = Σ(−1)^i ext^i` is accessible, because
it is a Riemann–Roch *numerical* invariant (already implemented: `exceptional.chi`
on `P^2`, `mukai.pairing = −χ` on K3). Individual `ext^i` are **not** a function of
Chern characters alone (they depend on the actual sheaves), so they are
**out of scope for the core by construction** — reachable only through a CAS bridge.

---

## 4. What is not currently computable in ANY software — separating open math from unimplemented

- **Automatic wall-and-chamber decomposition for arbitrary Chern characters (arbitrary
  surface).** The *numerical* walls are computable and **finite** in a fixed slice:
  Maciocia proves the walls are bounded/nested and gives finiteness criteria
  `[PROVEN]` (arXiv:1202.4587); Bertram's nested-wall theorem underlies this. This
  project already enumerates the numerical set (`walls.compute_walls`,
  `numerical_wall`). The hard part is certifying which numerical walls are **actual**
  destabilizing walls — that requires knowing an object of the sub-Chern-character
  *exists* and genuinely sub-objects `v`. On `P^2` the actual walls are fully
  characterized (ABCH arXiv:1203.0316; CHW), which is why `walls.actual_walls` is
  empirically complete there. For a **general** surface, a complete actual-wall
  classification is **partly open mathematics** (tied to the non-emptiness / moduli
  structure problem) and, even where theory exists (K3/abelian via
  Bridgeland/Bayer–Macrì), **not packaged** in any software. → *Both* open math *and*
  unimplemented, depending on the surface.
- **Automatic construction/verification of Bridgeland stability conditions on
  arbitrary threefolds.** The standard **sufficient** route to a stability condition
  on a threefold `X` is a generalized Bogomolov–Gieseker (BMT) inequality for
  tilt-stable objects on `X`; where that (strong) inequality holds, a stability
  condition exists. It is a **theorem for specific `X`** — `P^3`, the quadric, all
  Fano 3-folds of Picard rank 1, abelian 3-folds, the quintic `[PROVEN]` (the
  codebase's `Threefold.bg_proven=True` catalog: arXiv:1103.5010, 1207.4980,
  1510.04089, 1410.1585, 1810.03434). Crucially, the implication is *not* an
  equivalence in the other direction: existence of a stability condition can still
  hold via a **modified** inequality where the *strong* BMT form **fails**. The sharp
  example is `Bl_p(P^3)`: Schmidt showed the stronger/generalized BMT inequality
  fails there `[PROVEN]` (arXiv:1602.05055; encoded as `BLOWUP_P3_POINT.bg_proven=False`),
  yet `Bl_p(P^3)` is Fano and **does carry Bridgeland stability conditions** via a
  modified inequality `[PROVEN]` (Bernardara–Macrì–Schmidt–Zhao, arXiv:1607.08199,
  which proves existence on *all* Fano threefolds). So the correct statement is: the
  BMT inequality is the standard **sufficient** route, not an equivalence; "the
  strong BMT boundary is not rigorous on `Bl_p(P^3)`" does **not** mean "no stability
  condition exists there" (that is false). For a *general* (non-Fano, otherwise
  uncovered) threefold, whether the needed inequality holds is **OPEN**
  `[CONJECTURED]` — no software can close it because the theorem does not exist.
- **The full set of ACTUAL walls outside ABCH / Maciocia-covered cases.** As above:
  open where the destabilizer-existence theory is missing; unimplemented elsewhere.
- **BG-boundary where the inequality is unproven or not strong enough.** Which
  *stronger/modified* inequality gives the sharp boundary on varieties like
  `Bl_p(P^3)` is **unknown** `[CONJECTURED / open]`; Algorithm 5's boundary there is
  correctly flagged non-rigorous (`bg_proven=False`) — but that flag should be read
  as "strong BMT boundary not rigorous", **not** "no stability condition exists". No
  amount of engineering computes a *sharp* inequality that has not been proven.

The clean separation: *numerical wall enumeration, χ, BG on the proven catalog,
non-emptiness on `P^2`* are computable and (here) implemented; *actual-wall
certification on general surfaces* and *the sharp inequality / stability existence
on truly general threefolds* are **open mathematics**; *the packaged, turnkey
versions of the cases theory already covers beyond `P^2`* are **merely
unimplemented**.

---

## 5. Recommendation: backend strategy for this project

**Core question:** can an exact-`Fraction`, formula-level layer built on the existing
core generalize beyond `P^2` using *published closed-form / decision-procedure
results*, without a CAS — and where (if ever) does an optional CAS bridge earn its
keep?

**Answer: yes for a valuable slice — but the cheap, in-architecture wins are the
Picard-rank-1 varieties, not the rational surfaces.** It is worth being explicit
about the ordering, because "generalize beyond `P^2`" most evokes `P^1×P^1`,
Hirzebruch, and del Pezzo surfaces — and those are precisely the *expensive* ones.

**Genuinely near-term (fit the current `(r, c, e)` + lattice architecture natively,
category (c)):**

- **Wall certification on Picard-rank-1 surfaces (K3 / abelian).**
  Arcara–Bertram–Lieblich give the Picard-rank-1 wall-crossing family and its
  numerical description `[PROVEN]` (arXiv:0708.2247); Maciocia gives finiteness /
  location bounds turning the numerical wall set into a certified-finite one
  `[PROVEN]` (arXiv:1202.4587). Both are semicircle/lattice arithmetic — exactly
  `walls.numerical_wall`'s domain, and a class is still a single degree number.
  *Convention note:* the `(s, t)` semicircle geometry is identical for K3 and abelian
  surfaces, but the Mukai-lattice bookkeeping differs — on a **K3** the Mukai vector
  shifts `ch₂ → ch₂ + r` (so `radius²` picks up `+2/d`), because `√td = (1, 0, 1)`;
  on an **abelian** surface `√td = (1, 0, 0)` and there is **no** shift — the bare
  `(r, c, ch₂)` is already the correct Mukai vector. `mukai.from_chern`'s `ch₂+r`
  shift is therefore **K3-only** and must be documented as such; do not apply it to
  abelian surfaces.
- **K3 lattice machinery.** Bridgeland's stability conditions on K3s
  `[PROVEN]` (arXiv:math/0307164) and the Bayer–Macrì wall/chamber and MMP results
  `[PROVEN]` (arXiv:1301.6968) are **entirely computations in the Mukai lattice** —
  already the substrate of `mukai.py` (`pairing`, `self_pairing`, `classify_wall`).
  Extending `mukai.py` (rank-1 → the full rank-2 hyperbolic lattice analysis, totally
  semistable / fake walls) is pure lattice arithmetic, no CAS.
- **Higher-rank BG on surfaces.** Beyond the current `Δ ≥ 0`, Hoppe–Gieseker–type
  and Mukai/O'Grady/Marian–Oprea numerical bounds are closed inequalities in the
  Chern data `[PROVEN — as inequalities]` and drop straight into `bg_check.py`
  without touching the class representation.

**Gated behind a much larger refactor (do NOT mis-schedule as near-term):**

- **Non-emptiness on Hirzebruch / del Pezzo surfaces.** This is *not* a closed-form
  formula port. The relevant results are **decision procedures**, and they are tied
  to the anticanonical polarization, not an arbitrary `H`:
  - *del Pezzo:* Levine–Zhang prove **existence** of (semi)stable sheaves with
    respect to the **anti-canonical** polarization for degree ≥ 3, and compute the
    cohomology of a general sheaf for degree ≥ 4 `[PROVEN]` (arXiv:1910.14060) —
    but there is **no closed δ-curve** analogous to the `P^2` DLP curve.
  - *Hirzebruch:* the Coskun–Huizenga existence algorithm's core is a
    **prioritary-sheaf / Harder–Narasimhan** argument `[PROVEN]` (arXiv:1907.06739),
    whose HN-recursion step **exceeds pure `(r, c, e)` arithmetic** — it is partly
    sheaf-theoretic, not a numeric identity.

  Worse, all of these rational surfaces have **Picard rank ≥ 2** (`P^1×P^1`,
  Hirzebruch `F_n`, del Pezzo), so they cannot even be *represented* until the core
  is refactored from a scalar degree `c = ch₁·H` to a full **NS-lattice** Chern
  character. That refactor mutates the frozen `ChernChar` (`c` scalar → NS vector),
  `walls.py`, `mukai.py`, and **every one of the 42 pinned tests** — it is a
  **high-effort, strategic, v2-architecture decision**, not a medium tweak. Do not
  schedule del Pezzo/Hirzebruch non-emptiness as a Stage-0 deliverable; it depends on
  this refactor *and* on importing a partly sheaf-theoretic (prioritary + HN)
  decision procedure. `[SPECULATIVE — the scheduling judgement is mine; the
  underlying results are PROVEN]`

Because the near-term items are *formulas in `(r, c, e)` and lattice pairings*, they
respect the inviolable constraints (exact `Fraction`, CH convention, zero deps,
frozen dataclasses) and extend the pinned-test culture. This is where the project's
comparative advantage lies: **a CAS gives you `Ext` and moduli by construction, but
this project's value is the exact, citable, closed-form numerical layer that no CAS
packages together.** `[SPECULATIVE — the architectural judgement is mine; the
individual formulas are PROVEN]`

**Where a CAS bridge finally earns its keep** — only for genuinely **sheaf-level**
questions the numerics cannot answer:

- deciding non-emptiness / dimension of a moduli space on a surface where **no
  closed-form criterion exists** (construct sheaves, compute `Ext^1`/obstructions →
  M2 `Varieties`);
- exhibiting an **actual destabilizing sub-object** for a specific class (build the
  object, verify the inclusion and phases → M2/Singular);
- **exceptional collections / mutations** on new varieties (the arXiv:2509.25103 M2
  tool, `QuiversToricVarieties`, or QPA via a tilting equivalence);
- **ground-truth cross-checks**: recompute a `χ`, a cohomology dimension, or a wall
  destabilizer in M2 and diff against the formula layer — feeding the two-way
  verification discipline in `docs/CORRECTIONS.md`.

### Staged recommendation

- **Stage 0 — formula/lattice layer only (now, no CAS).** Add Picard-rank-1 surface
  wall certification (Maciocia + Arcara–Bertram), flesh out the K3 lattice module
  (Bridgeland + Bayer–Macrì), and add higher-rank surface BG bounds. **Unlocks:**
  mathematically correct, citable answers for the Picard-rank-1 classes the
  literature has closed-form results for, with zero new dependencies. *Explicitly
  excluded from Stage 0:* del Pezzo/Hirzebruch non-emptiness (needs the NS-lattice
  refactor and a prioritary/HN decision procedure — see above).
- **Stage 0.5 — NS-lattice refactor (strategic, v2).** The prerequisite for any
  Picard-rank-≥2 surface (the rational surfaces). High effort; touches the frozen
  `ChernChar`, `walls.py`, `mukai.py`, and all 42 pinned tests; schedule and cost it
  as an architecture decision, not a quick win.
- **Stage 1 — thin, optional Macaulay2 subprocess bridge (later).** Behind an
  optional install/import (never imported by the core), text I/O, lossless
  `Fraction↔QQ`. Used **only** for sheaf-level queries: `Ext^i`, moduli non-emptiness
  by construction, and as an independent oracle to certify formula-layer outputs.
  **Unlocks:** cases with no closed formula, plus a verification oracle that
  strengthens the correctness culture. **Cost, stated plainly:** an external M2
  install plus a hand-maintained text-protocol parser — a real maintenance burden,
  just not a pip one.
- **Stage 2 — OSCAR via `juliacall` (research escalation, separate track).** Toric
  coherent-sheaf cohomology and exceptional-collection / mutation experiments.
  **Unlocks:** toric-sheaf computation and (partial) mutation tooling — but note the
  hard ceiling: **no system provides a general derived-category-of-`Coh(X)` engine**,
  so a full derived-category layer is *not a tractable bolt-on* to an exact-`Fraction`
  core; Stage 2 buys specific, constructive sheaf/quiver computations, not
  "derived categories". Justified only if the roadmap explicitly moves past the
  numerical layer, and accepting the Julia/OSCAR install (Windows-via-WSL, itself
  `[UNVERIFIED]`) cost.

**One-line strategy:** stay formula-first and Picard-rank-1-first (it is where this
project is differentiated and where correctness is cheapest to guarantee); treat the
NS-lattice refactor as the gate to rational surfaces; add an optional Macaulay2
oracle when — and only when — a question is intrinsically sheaf-level; treat
OSCAR/derived-category work as a separate research escalation with no turnkey
derived-category engine anywhere to lean on.

---

## References

- I. Coskun, J. Huizenga, J. Woolf, *The effective cone of the moduli space of
  sheaves on the plane*, JEMS 2017 — arXiv:1401.1613.
- I. Coskun, J. Huizenga, *Existence of semistable sheaves on Hirzebruch surfaces*
  (prioritary-sheaf / HN existence algorithm) — arXiv:1907.06739.
- D. Levine, S. Zhang, *Brill–Noether and existence of semistable sheaves on del
  Pezzo surfaces* (existence w.r.t. the anti-canonical polarization, degree ≥ 3;
  general-sheaf cohomology, degree ≥ 4) — arXiv:1910.14060.
- D. Arcara, A. Bertram, I. Coskun, J. Huizenga (ABCH), *The minimal model program
  for the Hilbert scheme of points on `P^2` and Bridgeland stability*, Adv. Math.
  235 (2013) — arXiv:1203.0316.
- A. Maciocia, *Computing the walls associated to Bridgeland stability conditions on
  projective surfaces*, Asian J. Math. 18 (2014) 263–280 — arXiv:1202.4587.
- T. Bridgeland, *Stability conditions on K3 surfaces*, Duke Math. J. 141 (2008)
  241–291 — arXiv:math/0307164.
- A. Bayer, E. Macrì, *MMP for moduli of sheaves on K3s via wall-crossing*, Invent.
  Math. 198 (2014) 505–590 — arXiv:1301.6968.
- D. Arcara, A. Bertram, M. Lieblich, *Bridgeland-stable moduli spaces for K-trivial
  surfaces*, JEMS 15 (2013) 1–38 — arXiv:0708.2247.
- A. Bayer, E. Macrì, Y. Toda, *Bridgeland stability conditions on threefolds I:
  Bogomolov–Gieseker type inequalities* — arXiv:1103.5010.
- E. Macrì, *A generalized Bogomolov–Gieseker inequality for `P^3`* — arXiv:1207.4980.
- C. Li, *Stability conditions on Fano threefolds of Picard number 1* —
  arXiv:1510.04089.
- A. Bayer, E. Macrì, P. Stellari, *The space of stability conditions on abelian
  threefolds* — arXiv:1410.1585.
- C. Li, *On stability conditions for the quintic threefold* — arXiv:1810.03434.
- B. Schmidt, *Counterexample to the generalized Bogomolov–Gieseker inequality for
  threefolds* — arXiv:1602.05055.
- M. Bernardara, E. Macrì, B. Schmidt, X. Zhao, *Bridgeland stability conditions on
  Fano threefolds* (existence on all Fano threefolds via a modified inequality; the
  strong BMT form is not needed), Épijournal Géom. Algébrique 1 (2017) —
  arXiv:1607.08199.
- G. G. Smith, *Computing global extension modules for coherent sheaves on a
  projective scheme*, J. Symbolic Comput. 29 (2000) — arXiv:math/9807170.
- D. Eisenbud, M. Mustață, M. Stillman, *Cohomology on toric varieties and local
  cohomology with monomial supports*, J. Symbolic Comput. 29 (2000) 583–600 —
  arXiv:math/0001159.
- D. Eisenbud, D. Erman, F.-O. Schreyer, *Tate resolutions on products of projective
  spaces* — arXiv:1905.10230.
- M. K. Brown, D. Erman, et al., *The multigraded BGG correspondence in Macaulay2* —
  arXiv:2402.12293.
- *Toric geometry in OSCAR* — arXiv:2303.08110.
- *Computing higher direct images in Macaulay2* — arXiv:2505.22835.
- M. K. Brown, S. Dey, Y. Li, M. Sayrafi, *Computing global Ext for complexes*
  (`Ext` of bounded complexes + mutations of exceptional collections in Macaulay2;
  no spherical twists) — arXiv:2509.25103.
- N. Prabhu-Naik, *QuiversToricVarieties: a package to construct quivers of sections
  on complete toric varieties* (database of full strong exceptional collections of
  line bundles for smooth toric Fano varieties of dimension ≤ 4) — arXiv:1501.05861.
- J.-M. Drézet, J. Le Potier, *Fibrés stables et fibrés exceptionnels sur `P_2`*,
  Ann. Sci. ÉNS 18 (1985) 193–243.
- Macaulay2 package docs: `Varieties`, `Complexes`, `BGG`, `TateOnProducts`,
  `NormalToricVarieties` — macaulay2.com/doc (packages-provided list).
- QPA manual (v1.35, 2024); CAP / homalg-project; SageMath schemes/toric reference
  manual; Singular `sheafcoh.lib` (`sheafCoh`, `sheafCohBGG`) — software documentation.

---

## Adversarial review notes

**Overall Phase-1 verdict: proceed-with-amendments.**

Substantive changes made in response to the critic:

- **Fixed the Smith / EMS citation split.** `arXiv:math/9807170` is Gregory G.
  Smith's sole-author global-Ext algorithm (kept for the projective-scheme case);
  the toric Cox-module / local-cohomology reduction is now correctly attributed to
  Eisenbud–Mustață–Stillman, *Cohomology on toric varieties…*, **arXiv:math/0001159**
  (verified), removing the previous misattribution of the EMS work to the Smith ID.
- **Removed "spherical twists" from all mentions of arXiv:2509.25103.** Rescoped to
  Brown–Dey–Li–Sayrafi, *Computing global Ext for complexes* — `Ext` of bounded
  complexes + mutations of exceptional collections; any spherical-twist capability is
  now tagged `[UNVERIFIED]`.
- **De-flagged and cited `QuiversToricVarieties`** as Prabhu-Naik, **arXiv:1501.05861**
  (database of full strong exceptional collections of line bundles for smooth toric
  Fano varieties of dim ≤ 4); removed its `[UNVERIFIED]` tag.
- **Softened the threefold framing.** Changed "existence ⇔ generalized BG (BMT)
  inequality" to "BMT is the standard *sufficient* route"; added that `Bl_p(P^3)` is
  Fano and **does carry stability conditions** via a modified inequality
  (Bernardara–Macrì–Schmidt–Zhao, arXiv:1607.08199), so `bg_proven=False` means
  "strong BMT boundary not rigorous", not "no stability condition".
- **Qualified the del Pezzo/Hirzebruch non-emptiness item** as a *decision-procedure*
  port, not a closed-form formula port (Levine–Zhang arXiv:1910.14060, anticanonical,
  del Pezzo deg ≥ 3; Hirzebruch prioritary+HN, arXiv:1907.06739), noting the HN core
  exceeds `(r, c, e)` arithmetic, and moved it out of Stage 0.
- **Reconciled the effort estimate and re-ordered the roadmap.** The NS-lattice
  refactor is now marked high-effort / strategic (touches frozen `ChernChar`,
  `walls.py`, `mukai.py`, all 42 tests) and made an explicit gate: rational surfaces
  are *not* near-term; the cheap wins are Picard-rank-1 (K3/abelian walls, K3 lattice,
  higher-rank BG).
- **Unified the K3-vs-abelian Mukai convention:** K3 shifts `ch₂ → ch₂ + r`
  (`radius² += 2/d`); abelian does **not** (`√td = (1,0,0)`); documented
  `mukai.from_chern`'s shift as K3-only.
- **De-flagged the Singular `sheafcoh.lib` / `sheafCoh` / `sheafCohBGG` names**
  (confirmed) and honestly re-scoped the "zero Python packages" M2 bridge as still a
  real (non-pip) external-install + parser-maintenance burden.
- **Bounded Stage 2 / scope creep:** stated outright that no system provides a
  general derived-category-of-`Coh(X)` engine, so Stage 2 is not a tractable bolt-on;
  added that the arXiv:2509.25103 M2 package is the closest distributed
  exceptional-collection/mutation tool, partially softening the "no
  ExceptionalCollections package" claim.
