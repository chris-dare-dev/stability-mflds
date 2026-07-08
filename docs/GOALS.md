# GOALS — bridgeland_stability roadmap

This document is the final, adversarially-reviewed goal set for extending `bridgeland_stability`
beyond its current Picard-rank-1, numerical-Chern-character, exact-`Fraction` regime. It merges the
Phase-1 feasibility findings (`HOMOLOGICAL_ALGEBRA_SURVEY.md`, `GENERALIZATION_FEASIBILITY.md`,
`LITERATURE_GAP_ANALYSIS.md`) into one coherent DAG, then applies the Phase-2 adversarial critic's
per-goal dispositions and amendments. Every goal carries the project's confidence flags —
**[PROVEN]** (a cited theorem), **[CONJECTURED]** (an open conjecture), **[SPECULATIVE]** (an in-repo
inference), **[UNVERIFIED]** (a claim not yet pinned to a primary source) — and any `[RESEARCH]` tag
marks a dependence on reading a specific paper or on an open problem.

**Adversarial verdict: proceed-with-amendments.** The goal set is fundamentally sound and executable
by a solo maintainer while respecting the correctness culture. No goal hides an unsolved problem in a
way that would silently produce wrong values; the sheaf-level dependencies (G14/G16/G18) are
explicitly delegated, not concealed. Six goals received surgical corrections (none a redesign), two
goals (G17, G19) were deferred; see the appendix and the closing review notes.

## How to read this

Goals have **stable global IDs `G1..G19`** (unique across tiers; grouped by tier but never
renumbered). Each goal states: Tier, Mathematical prerequisites (with arXiv IDs), Code prerequisites
(modules + prior goal IDs), an Implementation sketch (concrete functions/signatures), Testability (at
least one exact literature-anchored or exactly-derived value — every Tier-1 goal names one), Research
value, Risks/unknowns, and Dependencies (by Goal ID; the union forms an acyclic DAG — see the
Dependency DAG section).

**Tier definitions.**
- **Tier 1 (1–3 months):** directly implementable in the *current* architecture (exact `Fraction`,
  numerical Chern-character level, no sheaf-level computation); each must be testable against a
  published or exactly-derived value.
- **Tier 2 (3–12 months):** moderate architectural extensions (full NS/Picard lattice, threefold
  tilt-wall solver, Maciocia ρ>1 walls, exceptional-collection generators) still implementable from
  KNOWN FORMULAS.
- **Tier 3 (12+ months):** significant new math/computational infrastructure (Python↔M2/OSCAR bridge
  for real Ext, sheaf-level HN/prioritary decision procedures, or a derived-category engine); MAY
  depend on open problems. Note a Tier-3 label may reflect **epistemic status** (dependence on an open
  conjecture) rather than raw effort — G15 is the clearest example, and is scheduled early.

**Global honesty invariants preserved throughout.** (a) The exact-`Fraction` core, the CH
discriminant convention `Δ = ½μ² − ch₂/(rd)`, frozen dataclasses, zero runtime deps, and the 42
pinned tests are inviolable ground truth. (b) The K3 Mukai shift `ch₂ → ch₂ + r`
(`mukai.from_chern`) is **K3-only** — never reuse it for abelian surfaces (abelian `√td = (1,0,0)`,
so the bare `(r,c,ch₂)` is already the Mukai vector). (c) `Bl_p(P³)` **is Fano and does carry
Bridgeland stability conditions** (BMSZ arXiv:1607.08199); only the *strong* BMT boundary fails there
(Schmidt arXiv:1602.05055) — no goal may ever emit "no stability condition exists" for it. (d) The
threefold BG inequality is a **theorem only on the `bg_proven=True` catalog rows** and
**[CONJECTURED]** in general. (e) No system anywhere provides a general derived-category-of-`Coh(X)`
engine — this bounds every CAS-bridge goal.

---

# Tier 1 — in current architecture, testable against published values

Ordered by value-per-effort (highest first). Every Tier-1 anchor was re-derived AND re-executed
against the shipped code before submission. Shared verified fixtures: abelian wall `(−5/2, 17/4)` vs
K3 wall `(−5/2, 21/4)` differ by exactly `2/d = 1` (G2/G3); `ν` on `(2,0,1,0)` gives `−1, −1/2, 0,
None` and `twist(1,1).a3 == −4/3` (G4); the identity `Δ = (c²/d − 2re)/(2r²d)` giving the K3 bound
`Δ ≥ (1 − 1/r²)/d ⟺ v² ≥ −2` (G6).

## G1 — Fix the miscited threefold references in `varieties.py`
- **Tier:** 1
- **Mathematical prerequisites:** None new — provenance only, all orchestrator-VERIFIED against arXiv
  abstracts **[PROVEN]**: C. Li, *Stability conditions on Fano threefolds of Picard number one*
  (arXiv:1510.04089, covers the quadric as a Fano ρ=1); Bernardara–Macrì–Schmidt–Zhao, *Bridgeland
  stability conditions on Fano threefolds* (arXiv:1607.08199 — the true "BMSZ"); Schmidt,
  *Bridgeland stability on threefolds — some wall crossings* (arXiv:1509.04608 — general threefold
  wall-crossing techniques, worked examples on P³ (twisted cubics / complete intersections), **not**
  quadric-specific; the quadric's own generalized BG inequality is Schmidt arXiv:1309.4265, subsumed
  by the Fano ρ=1 theorem 1510.04089); Piyaratne (arXiv:1705.04011, **solo** — NOT BMSZ); Schmidt counterexample
  (arXiv:1602.05055). Confirmed defects: `arXiv:1607.07182` is *"Interlacing Diffusions"* (a math.PR
  probability paper), and `arXiv:1705.04011` was mis-attributed to BMSZ.
- **Code prerequisites:** `varieties.py` only.
- **Implementation sketch:** Set `QUADRIC3.references = ["arXiv:1510.04089", "arXiv:1607.08199",
  "arXiv:1509.04608"]` and rewrite its `note` so it no longer implies `1705.04011 = BMSZ`. For
  `BLOWUP_P3_POINT`, KEEP both IDs `["arXiv:1602.05055", "arXiv:1705.04011"]` (both legitimately
  relevant) but correct the note to attribute `1705.04011` to **Piyaratne** (his optimal *modified*
  inequality), and add the nuance that `Bl_p(P³)` is Fano and **does** carry stability conditions
  (BMSZ 1607.08199) — only the strong BMT boundary fails, which is exactly what `bg_proven=False`
  encodes. Add a guard test (`tests/test_varieties.py`): iterate every catalog `Threefold` and assert
  every reference is a member of a small **curated algebraic-geometry allowlist** (a frozen set of the
  vetted AG arXiv IDs), not merely a well-formed arXiv string.
- **Testability:** `assert "arXiv:1607.07182" not in sum((tf.references for tf in ALL_THREEFOLDS),
  [])`; `assert QUADRIC3.references == ["arXiv:1510.04089","arXiv:1607.08199","arXiv:1509.04608"]`;
  `assert "arXiv:1705.04011" in BLOWUP_P3_POINT.references`; and, for every catalog reference `ref`,
  `assert ref in AG_ALLOWLIST`. **Guard-mechanism caveat (amended):** the format regex
  `^arXiv:\d{4}\.\d{4,5}$` is a **format check only** and does *not* distinguish an AG ID from a
  probability paper — `arXiv:1607.07182` matches it perfectly, so a regex-only guard would NOT have
  caught the original defect. The **membership-in-allowlist** assertion is the mechanism that actually
  catches a stray non-AG ID; the regex may be used as a cheap secondary sanity filter but must not be
  the primary check. The allowlist is maintenance-brittle (every legitimately added reference must be
  appended) — acceptable for this small catalog, but state it in the test docstring. No computed value
  or pinned test changes.
- **Research value:** This is a correctness-culture project whose differentiator is citation accuracy;
  this is the ONLY confirmed defect in shipped code. Anyone citing the catalog inherits the fix.
- **Risks and unknowns:** None to the math. The only ongoing cost is allowlist maintenance (above).
  Pure docs/citation edit; zero math-value change. **[PROVEN]** corrections.
- **Dependencies:** none. (Roadmap FIRST task.)

## G2 — Expose abelian-surface (s,t) walls (the free win)
- **Tier:** 1
- **Mathematical prerequisites:** An abelian surface has trivial tangent bundle, so `td = (1,0,0)`,
  `√td = (1,0,0)`, and the Mukai vector is the **bare** `(r, c₁, ch₂)` with **no** shift **[PROVEN]**.
  Feeding the bare Chern triple to `chern.central_charge` reproduces Bridgeland's abelian central
  charge `⟨exp((s+it)H), v(E)⟩` exactly, so `numerical_wall(v,w,d)` on the bare triple already IS the
  genuine abelian (s,t) wall **[PROVEN, ρ=1]**. Literature anchors: Yanagida–Yoshioka
  (arXiv:1203.0884); Maciocia–Meachan (arXiv:1107.5304); Bayer–Macrì (arXiv:1301.6968).
- **Code prerequisites:** `walls.numerical_wall`, `varieties.abelian_surface` (both exist). No new
  math.
- **Implementation sketch (revised — scope the free win precisely):** No new computation — the
  deliverable is a thin, well-documented public affordance plus tests. **Only the semicircle
  primitive is a genuine free win:** `numerical_wall(v, w, surface.d)` on the bare `(r,c,e)` (and a
  convenience wrapper `abelian_wall(v, w, surface)` that asserts `surface.kind == "abelian"` and
  returns exactly that, guarding against accidental application of the K3 `+r` shift) returns the
  genuine abelian (s,t) semicircle with NO Mukai shift. **Do NOT claim this extends to
  `compute_walls`/`actual_walls`.** `actual_walls`' completeness is certified only for `P²[n]`, and
  its `_is_integral_chern` filter tests `c₁²/2 − ch₂ ∈ ℤ` — the `d = 1` P² Chern lattice condition —
  so on an abelian surface (`d ≠ 1`) it applies a *wrong* integrality test and layers on the ABCH
  rank-reduction / heart-ordering necessary conditions that are P²-specific; a genuine abelian
  destabilizer can be spuriously kept or dropped, and the "these are ALL the walls" certification does
  not hold. Document loudly that `abelian_wall`/`numerical_wall` are the abelian entry points, and
  route abelian input away from the K3 shim of G3.
- **Testability:** With `d = abelian_surface(2).d == 2`, `v = ChernChar(1,0,-2)`, `w =
  ChernChar(1,-1,1/2)`: `numerical_wall(v,w,2).center == F(-5,2)` and `.radius_sq == F(17,4)`
  (RE-VERIFIED by execution: `W_rc=-1, W_re=5/2, W_ce=-2 ⟹ ρ² = 25/4 − 2 = 17/4`). This is exactly
  `2/d = 1` LESS than the K3 shim value `F(21,4)` of G3 — the "abelian ≠ K3 by exactly 2/d" assertion
  is the key regression guard against the K3 shift being wrongly applied to abelian input.
- **Research value:** Makes abelian-surface wall diagrams first-class at zero math cost; used by
  anyone studying Bridgeland walls / moduli on abelian surfaces (Yanagida–Yoshioka MMP program).
- **Risks and unknowns:** The chief risk is exactly the over-claim removed above: a future contributor
  routing abelian classes through `actual_walls` (wrong lattice filter) or the K3 shim; the regression
  test and the `abelian_wall` guard catch both. `ρ≥2` abelian surfaces (`E₁×E₂`) are NOT covered —
  gated by the Tier-2 NS-lattice refactor G12 (do not regress ρ>1 work into Tier 1).
- **Dependencies:** none. (Shares test fixtures with G3; documented alongside it.)

## G3 — K3-only `k3_wall` shim (ch₂→ch₂+r) + connect semicircle to `classify_wall`
- **Tier:** 1
- **Mathematical prerequisites:** `√td_{K3} = (1,0,1)`, so the K3 Mukai vector is `(r, c₁, ch₂+r)` and
  Bridgeland's K3 central charge (arXiv:math/0307164) is reproduced by feeding the Mukai triple to
  `chern.central_charge`. Under `e ↦ e+r`, `e' ↦ e'+r'`, the minors satisfy `W_rc` unchanged, `W_re`
  unchanged, `W_ce ↦ W_ce − W_rc`, so the wall **center is invariant** and **radius² increases by
  exactly `+2/d`** **[PROVEN, ρ=1]** (re-checked and executed). Wall-type layer: Bayer–Macrì Thm 5.7
  (arXiv:1301.6968), already in `mukai.classify_wall`.
- **Code prerequisites:** `walls.numerical_wall`, `mukai.MukaiVector`/`classify_wall` (all exist).
- **Implementation sketch:** Add `k3_wall(v: ChernChar, w: ChernChar, d: int) -> Wall` mapping
  `(r,c,e) ↦ (r,c,e+r)` for both classes and calling `numerical_wall`. Place it in `mukai.py`
  importing `numerical_wall`/`ChernChar` — acyclic, since `walls.py` does not import `mukai`; this
  keeps the shift co-located with the Mukai convention. Document loudly that it is **K3-only** and must
  not receive abelian input. Second function `k3_wall_classified(v, w, surface) -> (Wall,
  WallClassification)`: build Mukai vectors via `MukaiVector.from_chern(r, l, ch2)` (recall `l = c/d`;
  **assert `c % d == 0`** and raise a clear error otherwise), compute the semicircle via `k3_wall`,
  and the Bayer–Macrì type via `classify_wall`, returning both — closing the currently-disconnected
  gap that `mukai.py` never calls `numerical_wall`.
- **Testability:** `d=2`, `v=ChernChar(1,0,-2)`, `w=ChernChar(1,-1,1/2)`: `k3_wall(v,w,2).center ==
  F(-5,2)` (invariant vs bare) and `.radius_sq == F(21,4)` = bare `F(17,4)` + `2/d` (RE-VERIFIED by
  execution). General invariant on random small triples: `k3_wall(...).center ==
  numerical_wall(...).center` and `k3_wall(...).radius_sq − numerical_wall(...).radius_sq ==
  Fraction(2,d)`. For the classified pairing, reuse the pinned Mukai example `v=MukaiVector(2,0,-1)`,
  isotropic `u=(1,0,0)`, `(u,v)=1` → hilbert–chow divisorial (already pinned in `test_mukai.py`).
  **Fixture caveat (amended):** the shared fixture `w=(1,-1,1/2)` on `d=2` has `l = c/d = −1/2`
  (non-integral), so it is a **SYNTHETIC `+2/d` regression demonstration** (exactly-derived,
  **[PROVEN]**), *not* a genuine Picard-rank-1 K3 lattice class. Additionally pin one **integral-`l`**
  K3 class as a literature-anchored wall so the "genuine K3 (s,t) wall" claim is backed by a class
  actually in the lattice; retain the `c % d == 0` assertion in `k3_wall_classified` (the semicircle
  primitive `k3_wall` itself imposes no such restriction).
- **Research value:** Correct K3 (s,t) wall diagrams (previously off by `+2/d`) unified with wall
  *type* — the single most-used K3 output for wall-crossing / MMP studies (Bayer–Macrì).
- **Risks and unknowns:** `l = c/d` integrality: `classify_wall` requires integral `l`, so
  `k3_wall_classified` asserts/raises on non-integral `l`. `ρ≥2` K3s and totally-semistable/fake-wall
  certification beyond the search box are NOT covered (Tier 2/3, G10/deferred G17).
- **Dependencies:** none. (Contrast test shared with G2.)

## G4 — `ν_{α,β}` tilt-slope function on `ThreefoldChern`
- **Tier:** 1
- **Mathematical prerequisites:** Standard tilt slope **[PROVEN]** (BMT arXiv:1103.5010; Schmidt
  arXiv:1509.04608): `ν_{α,β}(E) = (H^{n-2}·ch₂^β − ½α²H^n·ch₀^β)/(H^{n-1}·ch₁^β)`, `+∞` when
  `H^{n-1}·ch₁^β = 0`. In `ThreefoldChern(r,a1,a2,a3)` with `d3=H³`: numerator `= a2^β − ½α²·r·d3`,
  denominator `= a1^β`, where `a1^β,a2^β` are `twist(β)` components. Note `ν` depends only on
  `(ch₀,ch₁,ch₂)` — NOT on `ch₃`.
- **Code prerequisites:** `ThreefoldChern.twist` (exists).
- **Implementation sketch:** Add `nu(v: ThreefoldChern, alpha_sq: Number, beta: Number, d3: int) ->
  Optional[Fraction]`. Take `alpha_sq = α²` as an exact `Fraction` (NOT `α`) to keep `√`-valued radii
  representable exactly. Compute `tw = v.twist(beta, d3)`; return `None` (encoding `+∞`, the
  vertical/degenerate case) when `tw.a1 == 0`; else `return (tw.a2 - Fraction(alpha_sq)/2 * v.r * d3)
  / tw.a1`. Pure `Fraction`, no float; mirror the existing `bmt_Q_at` naming.
- **Testability:** P³ null-correlation `v=ThreefoldChern(2,0,1,0)`, `d3=1`. **These four values are
  EXACTLY-DERIVED from the BMT tilt-slope formula above** (not read out of Schmidt) and RE-VERIFIED by
  execution: `nu(v, 0, 1, 1) == -1`; `nu(v, 1, 1, 1) == F(-1,2)`; `nu(v, 2, 1, 1) == 0`; `nu(v, 0, 0,
  1) is None` (β=0 gives `a1^β=0`, vertical). Consistency with the pinned twist: `v.twist(1,1).a3 ==
  F(-4,3)` (RE-VERIFIED: `a3^β = 0 − 1·1 + 0 − (1/6)·2·1 = −4/3`). **Amendment:** before any claim of
  cross-paper agreement with Schmidt arXiv:1509.04608, confirm the H-power normalization (num =
  `H^{n-2}ch₂^b − ½α²H^n ch₀^b`, den = `H^{n-1}ch₁^b`) matches that source's convention; until then
  the pinned values are self-consistent exact derivations, not a Schmidt cross-check.
- **Research value:** The missing scalar primitive for all threefold tilt-stability work; direct
  prerequisite for the Tier-2 threefold tilt-wall solver G9. Used by anyone reproducing Schmidt's
  threefold wall computations.
- **Risks and unknowns:** Sign/normalization conventions for `ν` vary across sources; the exact-derived
  pins insulate the code from that, but the Schmidt cross-check must confirm normalization (above).
  The tilt WALL locus (where two `ν` coincide) is handled in G9, not here — this goal ships only the
  slope function.
- **Dependencies:** none. (Unlocks G9.)

## G5 — `Rigor` provenance lattice + `Certificate`, propagated by minimum
- **Tier:** 1
- **Mathematical prerequisites:** None (engineering). Encodes existing epistemic distinctions: P² DLP
  / P²[n] walls **[PROVEN]**; threefold BMT boundary **[PROVEN]** on `bg_proven=True` rows,
  **[CONJECTURED]** / strong-form-false on `Bl_p(P³)` (Schmidt arXiv:1602.05055; nuance from BMSZ
  arXiv:1607.08199); `compute_walls` dense set and `actual_walls_complete`'s doubling certificate =
  **[HEURISTIC]**.
- **Code prerequisites:** touches `varieties.py`, `walls.py`, `threefold.py`, `mukai.py` at the
  attribute/repr level only (no formula changes). Build EARLY, in parallel with G1.
- **Implementation sketch:** New `rigor.py`: `class Rigor(IntEnum): PROVEN=3; CONJECTURAL=2;
  HEURISTIC=1; UNKNOWN=0` (total order so `min` = meet) and `@dataclass(frozen=True) class
  Certificate: rigor: Rigor; hypotheses: tuple[str,...]; citations: tuple[str,...]; note: str = ""`.
  Provide `meet(*certs)` returning a `Certificate` whose `rigor = min(...)` and whose
  `hypotheses`/`citations` are the unions — generalizing the single `bg_proven` boolean to a
  composable taint lattice. Attach, don't rewire: give `Surface`/`Threefold` an optional `certificate`
  field (default `UNKNOWN`); populate the catalog (P2 → PROVEN with DLP citations; `BLOWUP_P3_POINT` →
  CONJECTURAL/UNKNOWN for the *strong BMT boundary* with the "stability conditions DO exist" note).
  Let `Wall`/`ActualWall`/`BGBoundary` optionally carry a `Certificate = meet(algorithm_rigor,
  *input_certificates)`. Render the tag in `__repr__` (and expose to `viz` for solid/dashed/dotted
  styling). Hypothesis fields are DATA, not proof — sheaf-free.
- **Testability:** `min(Rigor.PROVEN, Rigor.CONJECTURAL) == Rigor.CONJECTURAL`; `meet(proven_alg,
  conjectural_input).rigor == Rigor.CONJECTURAL`; `P2.certificate.rigor == Rigor.PROVEN`;
  `BLOWUP_P3_POINT.certificate.rigor <= Rigor.CONJECTURAL` with its note distinguishing "strong BMT
  boundary not rigorous" from "no stability condition" (assert the phrase present/absent
  respectively); a P²[2] `ActualWall.certificate.rigor == Rigor.PROVEN`; a `compute_walls` output
  tagged `HEURISTIC`.
- **Research value:** Machine-readable answer to "the same function is a theorem on P² and a conjecture
  on a general threefold"; makes `docs/CORRECTIONS.md`'s "verified two ways" culture programmatic.
  Orthogonal — it colours every downstream goal in all three tiers.
- **Risks and unknowns:** Must NOT attempt to prove hypotheses — it only records which theorem is
  invoked and whether its stated hypotheses are met; keep it sheaf-free. Adding optional fields to
  frozen dataclasses is backward-compatible (defaults) and touches no pinned value. Low risk.
- **Dependencies:** none. (Soft-consumed by G6, G7, G9, G10, G11, G16, deferred G17, G18; hard-consumed
  by G15. Build early.)

## G6 — Higher-rank surface BG bounds (Yoshioka / Yanagida–Yoshioka existence) in `bg_check`
- **Tier:** 1
- **Mathematical prerequisites:** Rank-dependent sharpenings of Bogomolov `Δ≥0`, all closed
  inequalities in Chern data **[PROVEN as inequalities]**. **K3 existence** (Yoshioka): for primitive
  Mukai `v` and generic polarization, `M_H(v) ≠ ∅ ⟺ v² ≥ -2`. **Abelian existence** (Yanagida–Yoshioka
  arXiv:1203.0884): `M_H(v) ≠ ∅ ⟺ v² ≥ 0`. Both reduce to the algebraic identity (RE-VERIFIED by
  execution) `c²/d − 2re = v² + 2r²`, i.e. `Δ = (c²/d − 2re)/(2r²d)`, so the K3 bound is exactly
  `Δ ≥ (1 − 1/r²)/d ⟺ v² ≥ -2`. References: Bayer–Macrì arXiv:1301.6968; Yoshioka (moduli on
  abelian/K3).
- **Code prerequisites:** `chern.bogomolov_discriminant`, `mukai.self_pairing`, `mukai.from_chern`
  (exist); `bg_check`. Optionally consumes G5 (`Certificate`).
- **Implementation sketch (revised — add an integrality guard):** Extend `bg_check` with
  rank-and-surface-aware existence checks beyond `Δ≥0`. Add `check_existence_k3(v, surface)` returning
  the verdict `Δ ≥ (1 − 1/r²)/d` (equivalently `v² ≥ -2`), and `check_existence_abelian(...)`
  returning `v² ≥ 0` (using the **bare, un-shifted** Mukai self-pairing — do NOT apply the K3 `+r`
  shift). **Before computing `l`, assert the class lies in the Picard lattice:** `check_existence_k3`
  (and `check_existence_abelian`) must `assert v.c % surface.d == 0` and raise a clear error otherwise,
  *then* set `l = v.c // surface.d`. Without this guard, non-Picard-lattice input silently floors: e.g.
  `v = ChernChar(2,1,0)` on `K3(2)` has `c//d = 1//2 = 0`, so `from_chern(2,0,0) = (2,0,2)` gives
  `v² = -8` and the class is wrongly REJECTED — the honest behavior is that `l = 1/2` is not in the
  Picard lattice, so it must raise (RE-VERIFIED against the shipped `MukaiVector.from_chern`, which
  takes an integer `l`). Surface the threshold `Δ_min(r,d) = (1 − 1/r²)/d` as a helper so callers see
  WHY a class fails (the Mukai/O'Grady moduli-dimension bound, sharper than Bogomolov for `r ≥ 2`).
  Keep the existing `check_bg_surface` (`Δ≥0`) untouched as the rank-independent baseline.
- **Testability:** K3 structure sheaf `v=ChernChar(1,0,0)` on `K3(2)` (note ch=(1,0,0), not (1,0,1)):
  `Δ==0`, Mukai `(1,0,1)`, `v²==-2` — exactly on the r=1 existence boundary `(1 − 1/1²)/2 == 0`
  (nonempty; the rigid spherical class, matching the pinned `v(O)=(1,0,1)`, `⟨v,v⟩=-2`). Empty rank-2
  example `v=ChernChar(2,2,0)` on `d=2`: `Δ==1/8 < 3/8 == (1−1/4)/2`, `v²==-6 < -2` → REJECT
  (VERIFIED). Nonempty rank-2 `v=ChernChar(2,0,-3)`: `Δ==3/4 ≥ 3/8`, `v²==4` → ACCEPT, `moduli_dim ==
  v²+2 == 6` (VERIFIED). Abelian: `v=ChernChar(1,0,0)` → bare `v²==0`, boundary nonempty (VERIFIED).
  Guard: `check_existence_k3(ChernChar(2,1,0), K3(2))` RAISES (not a silent `v²=-8`).
- **Research value:** Turns K3/abelian rows from "Bogomolov only" toward the *actual* non-emptiness
  criterion for Picard-rank-1 K3/abelian; used by anyone deciding whether a moduli space is nonempty
  without a CAS.
- **Risks and unknowns:** Hoppe–Gieseker on `Pⁿ` is a *cohomological* (H⁰-vanishing) criterion whose
  core is NOT pure-Chern — keep it OUT of scope (ship only the K3/abelian closed inequalities; flag the
  cohomological criterion as CAS territory, G16). Yoshioka's `v²≥-2` requires **primitive** `v` and
  **generic** polarization — encode those as `Certificate.hypotheses`; **do not silently assert
  existence for imprimitive `v`.** General-surface (arbitrary ρ) sharp bounds are OPEN/**[CONJECTURED]**
  and NOT Tier 1 — restrict this goal to the Picard-rank-1 K3/abelian closed forms.
- **Dependencies:** none required; optional integration with G5.

## G7 — Exact wall-completeness certificate via the Maciocia radius bound (ρ=1 / ABCH classes)
- **Tier:** 1
- **Mathematical prerequisites:** Maciocia (arXiv:1202.4587) Thm 3.11: for a fixed class the
  semicircular walls are **nested** and lie inside one bounding semicircle; the wall set is finite with
  an explicit location bound **[PROVEN]**. ABCH (arXiv:1203.0316): for `P²[n]` the outermost (Gieseker)
  wall is center `−(2n+1)/2`, radius `(2n−1)/2` **[PROVEN]** (already pinned). Bertram's nested-wall
  theorem underlies boundedness.
- **Code prerequisites:** `walls.actual_walls` / `actual_walls_complete` (exist). Optionally consumes
  G5.
- **Implementation sketch:** Replace the EMPIRICAL doubling certificate in `actual_walls_complete` with
  a proven radius/center bound. Add `maciocia_wall_bound(v: ChernChar, d: int) -> Fraction` returning
  the Thm-3.11 bounding-semicircle radius² for `v`. Change `actual_walls_complete` to run
  `actual_walls` with a search window chosen to *cover* the Maciocia bound, then certify completeness by
  the theorem (return `complete = True` because the search provably exhausts the bounding semicircle)
  rather than by comparing against doubled bounds. Where the closed-form constant is extracted and
  verified, tag the certificate `Rigor.PROVEN` (G5); where not yet extracted, fall back to the current
  `HEURISTIC` doubling — so the function never over-claims.
- **Testability:** For `P²[n]` (`v=ChernChar(1,0,-n)`), `n∈{2,3,4,5}`: `maciocia_wall_bound(v,1) ≥
  ((2n−1)/2)²` (contains the Gieseker wall, already pinned), and `actual_walls_complete(v, P2)` returns
  the same wall set as today with `complete=True` now backed by the theorem — regression: the P²[2] set
  is exactly the single wall center `−5/2`, radius `3/2`. Assert no returned wall has `radius² >
  maciocia_wall_bound(v,d)`.
- **Research value:** Upgrades the completeness claim from **[SPECULATIVE]** (empirical doubling) to
  **[PROVEN]** for covered classes — anyone relying on "these are ALL the walls" gets a theorem.
- **Risks and unknowns:** **[RESEARCH]** Extracting the EXACT closed-form radius constant from Maciocia
  §3 and validating it against the `P²[n]` Gieseker value; until done, keep the `HEURISTIC` fallback.
  General-`v`, general-surface completeness is OPEN outside ABCH/Coskun–Huizenga-covered cases — do NOT
  upgrade the docstring honesty for arbitrary `v`. The general-surface extension escalates to Tier 2
  (G13, which needs the NS `γ²` decomposition).
- **Dependencies:** none required; optional integration with G5.

## G8 — Documentation-only worked example (K3 + abelian (s,t) walls), honest about model limits
- **Tier:** 1
- **Mathematical prerequisites:** Reuses G2 (abelian) and G3 (K3); `√td` conventions
  (arXiv:math/0307164); Bayer–Macrì types (arXiv:1301.6968). No new math.
- **Code prerequisites:** G2 and G3 shipped. Optionally references G5 tags.
- **Implementation sketch:** Add a short `docs/` worked example (or an `examples/` script following
  `examples/demo.py` style) that: (1) computes the SAME class `v=(1,0,-2)` on both `abelian_surface(2)`
  and `K3(2)`, shows the abelian wall `(−5/2, 17/4)` (via `abelian_wall`) vs the K3 wall `(−5/2,
  21/4)` (via `k3_wall`), and explains the exact `+2/d` shift and WHY (`√td = (1,0,0)` vs `(1,0,1)`);
  (2) attaches the Bayer–Macrì wall type via `k3_wall_classified`; (3) states plainly the current-model
  limits — Picard-rank-1 only, no `E₁×E₂`, no `ρ≥2` K3, `mukai.from_chern`/K3 shim are K3-only, and
  `actual_walls` is NOT an abelian/K3 completeness oracle. Keep it stdlib-only (no viz) so it runs in
  the zero-dependency core.
- **Testability:** A doctest / example-smoke test asserting the two printed radius² values are `F(17,4)`
  and `F(21,4)` and differ by exactly `F(2,2) == 1`; the script exits 0 under `python examples/…`.
- **Research value:** Onboarding / citable illustration of the single subtlest convention in the
  toolkit (K3-vs-abelian Mukai shift), preventing the exact misuse the docs warn about.
- **Risks and unknowns:** Must stay honest — explicitly NOT claim `ρ≥2` coverage nor `actual_walls`
  completeness for these surfaces. Doc-only, lowest risk; lowest value-per-effort of Tier 1 (hence
  last), but cheap.
- **Dependencies:** G2, G3 (documents both); optional G5.

---

# Tier 2 — moderate architectural extensions from known formulas

Ordered by value-per-effort. Scheduling gate: the NS-lattice refactor G12 must precede G13 and G14
regardless of listing order, and it gates the *faithful-computation* half of G11. **G12 is the
roadmap's dominant effort/schedule risk** (see Pipeline-blocking goals).

## G9 — Threefold tilt-stability wall solver in (α,β), with bg_proven-gated Bridgeland labeling
- **Tier:** 2 (RESEARCH-light — the tilt-layer conic is resolved; the Bridgeland-layer coefficient
  transcription is not yet)
- **Mathematical prerequisites (revised — corrected citation attribution):** BMT tilt-slope /
  tilt-stability (arXiv:1103.5010; Schmidt survey arXiv:1607.01262). **Tilt layer — Schmidt
  arXiv:1509.04608 Theorem 3.3 (VERIFIED against the paper):** numerical *tilt* walls are
  `x(α²+β²) + yβ + z = 0` — β-axis-centered semicircles (or vertical rays) in the (α,β) upper
  half-plane, with data `(R,C,D) = (ch₀, ch₁, ch₂)` and **no `ch₃`**. Because `ν_{α,β}` (G4) depends
  only on `(ch₀,ch₁,ch₂)`, the *tilt* wall between two fixed classes is unconditionally a conic and is
  structurally IDENTICAL to the surface `numerical_wall` on the truncated triple `(r, ch₁·H², ch₂·H)`
  with `d → d3` **[PROVEN]**. Worked example: the twisted-cubic ideal sheaf `I_C ⊂ P³` has `ch(I_C) =
  (1, 0, −3, 5)`. **Second-tilt / Bridgeland layer — this is a DIFFERENT theorem.** Thm 3.3 covers the
  tilt layer only; the second-tilt (Bridgeland) wall, where the cubic `ch₃^β` term enters the central
  charge `Z_{α,β}`, is governed by Schmidt's second-tilt / central-charge construction (Schmidt
  arXiv:1509.04608 §3 and the BMT central-charge framework arXiv:1103.5010), NOT by Thm 3.3. Do **not**
  cite Thm 3.3 for the `ch₃`-dependent semicircle — that would be a miscitation. Piyaratne–Toda
  arXiv:1504.01177 (moduli properness *assuming* BG).
- **Code prerequisites:** G4 (the `ν_{α,β}` tilt-slope function — hard, cross-tier). Existing
  `threefold.py` (`ThreefoldChern`, `.twist`, `bmt_Q`). Soft: G5 for the tilt-only-vs-Bridgeland
  label.
- **Implementation sketch:** Add `numerical_tilt_wall(v: ThreefoldChern, w: ThreefoldChern, d3: int)
  -> ThreefoldTiltWall`. Because `ν` depends only on `(r, a1, a2)`, the tilt wall
  `ν_{α,β}(v)=ν_{α,β}(w)` is computed by the *same three 2×2 minors* as `walls.numerical_wall`, applied
  to `ChernChar(r,a1,a2)` / `ChernChar(r',a1',a2')` with degree `d3`, in the `(β,α)` plane (β↔s, α↔t);
  this tilt layer can delegate to `numerical_wall`. Return a frozen `ThreefoldTiltWall(center,
  radius_sq, subobject, v, bridgeland_certified: bool)` mirroring `walls.Wall`; degenerate `W_rc=0` → a
  vertical `β=const` wall. The genuinely Tier-2 content is the **second-tilt/Bridgeland wall**
  `bridgeland_wall(v, w, threefold)` which uses the full central charge `Z_{α,β}` (Schmidt §3), so
  `a3=ch₃` enters the `y, z` coefficients of `x(α²+β²)+yβ+z=0`; the `x` coefficient still comes only
  from the `(r,a1)` minor (why it remains a circle). Compute `center = −y/(2x)`, `radius_sq =
  (y²−4xz)/(4x²)` exactly as `Fraction`. **Keep the exact `(y,z)` coefficient transcription and the
  twisted-cubic center/radius pin flagged `[RESEARCH-light]` until read off Schmidt's second-tilt
  construction** — do not attribute them to Thm 3.3. Set `bridgeland_certified = threefold.bg_proven`
  (a genuine Bridgeland wall only where BG is a theorem; Piyaratne–Toda). Finally add
  `enumerate_tilt_walls(v, threefold, bounds…)` analogous to `actual_walls`.
- **Testability:** Pinned input `ThreefoldChern(1, 0, −3, 5)` on `P3` (twisted cubic). **ch₃-
  independence of the TILT wall (exact):** `numerical_tilt_wall(v, w)` is unchanged when `a3, a3'`
  vary. **Reduction identity (exact):** `numerical_tilt_wall(v,w,d3).center` equals
  `numerical_wall(ChernChar(v.r,v.a1,v.a2), ChernChar(w.r,w.a1,w.a2), d3).center`, radius² matching up
  to the fixed `d3` normalization. **Semicircle form:** the fitted α² and β² coefficients of
  `bridgeland_wall` are equal (assert β-axis-centered). **bg_proven gate:** `bridgeland_certified` is
  `True` for `P3`/`QUADRIC3`, `False` for `BLOWUP_P3_POINT`. **[RESEARCH-light]** pin the exact
  twisted-cubic Bridgeland-wall centers/radii from Schmidt §7 once the second-tilt `(y,z)` coefficients
  are read off the paper.
- **Research value:** Threefold wall-crossing researchers (Schmidt, Macrì, Li circle); the first
  faithful threefold (α,β) walls in the toolkit on `P³`, the quadric, Fano ρ=1, abelian and quintic
  threefolds, with honest tilt-only-vs-Bridgeland labeling.
- **Risks and unknowns:** The second-tilt central-charge coefficients (`y,z` with `a3`) must be
  transcribed carefully from Schmidt's second-tilt construction — **[RESEARCH-light]**. The tilt-layer
  conic form is de-risked (Thm 3.3), but the Bridgeland-layer form must be read off the second-tilt
  source, not asserted from Thm 3.3. On `bg_proven=False` threefolds the Bridgeland upgrade is
  conjectural (correctly flagged, not asserted). No sheaf theory needed.
- **Dependencies:** Hard: G4. Soft: G5. Feeds G15.

## G10 — K3 rank-2 hyperbolic-lattice certified wall analysis (beyond the search box)
- **Tier:** 2
- **Mathematical prerequisites:** Bayer–Macrì arXiv:1301.6968 (Thm 5.7 wall types; the "potential
  wall" = a primitive rank-2 hyperbolic sublattice `H_W` containing `v`). Bridgeland
  arXiv:math/0307164. The certification replaces the bounded `|a|,|b|≤search` scan with exact
  enumeration of spherical (`x²=−2`) and isotropic (`x²=0`) classes in the **saturation** of `⟨v,w⟩` by
  solving the associated binary-quadratic (Pell-type) Diophantine equation.
- **Code prerequisites:** Existing `mukai.py` (`MukaiVector`, `pairing`, `self_pairing`,
  `classify_wall`). **No NS refactor needed** — the Picard-rank-1 K3 Mukai lattice is already the
  rank-3 `(r,l,s)` lattice with pairing `d l l' − r s' − r' s`. Soft: G3 (the `k3_wall` (s,t) shim, to
  link classification to geometry); G5 for the `certified` flag.
- **Implementation sketch (revised — saturated sublattice + positive-only certification):** Form the
  rank-2 sublattice and its Gram `[[v², (v,w)], [(v,w), w²]]`; a class `x = a v + b w` has `x² = a²v² +
  2ab(v,w) + b²w²`. Add `solve_binary_quadratic(A, B, C, target) -> list` returning primitive integer
  solutions of `A a² + B ab + C b² = target` for `target ∈ {−2, 0}`; for the hyperbolic (signature
  (1,1)) case it returns **one representative per orbit of the sublattice automorphism group** (the
  form is indefinite, so lattice points are infinite but orbits are finite). **(a) Enumerate over the
  SATURATED primitive rank-2 sublattice `H_W`** — the Bayer–Macrì potential wall, i.e. the primitive
  closure of `ℤv ⊕ ℤw` in the Mukai lattice — not merely the `ℤ`-span, or the "certified-complete"
  claim fails when a witness lives in the saturation but not the span. Then
  `classify_wall_certified(v, w, d) -> WallClassification`: (i) check `⟨v,w⟩` is genuinely rank 2 and
  hyperbolic, (ii) enumerate the minimal spherical/isotropic classes exactly over `H_W`, (iii) apply
  Thm 5.7 with a proof witness — no `search` parameter. **(b) Set `certified=True` ONLY for positive
  verdicts** (divisorial/flopping, where a witness is exhibited); a "fake-or-none" verdict must stay
  `certified=False` and defer the totally-semistable split to deferred G17. Extend `WallClassification`
  with `certified: bool` and `lattice_gram: tuple`. Bridge to `k3_wall` (G3) so a certified type
  carries its (s,t) semicircle.
- **Testability:** Pins already in `test_mukai.py`: `v(O)=(1,0,1)` spherical (`v²=−2`); Hilbert–Chow
  `v=(2,0,−1)`, isotropic `u=(1,0,0)`, `(u,v)=1`; LGU `v=(3,0,−1)`, `u=(2,0,0)`, `(u,v)=2`.
  **Certification invariant:** on all pinned *positive* cases `classify_wall_certified` returns the
  SAME `wall_type`/`subtype` as `classify_wall(..., search=large)` AND sets `certified=True`, with the
  verdict provably independent of any search bound; on a fake-or-none case it returns
  `certified=False`. **Solver test:** `solve_binary_quadratic` reproduces the spherical/isotropic
  witnesses `(a,b)` for the pinned lattices exactly, as one representative per automorphism orbit.
- **Research value:** K3 wall-crossing / MMP researchers (Bayer–Macrì program): certified wall types
  with an exhibited lattice witness, replacing a heuristic bounded search — directly strengthens
  `mukai.py`.
- **Risks and unknowns:** The hyperbolic binary-quadratic solver (finitely many orbits, infinitely
  many points) is a known but fiddly algorithm — must return orbit representatives correctly, and must
  saturate (above). **[RESEARCH]** the general fake/totally-semistable decision beyond the rank-2
  sublattice needs the full stratification (hard; escalates to deferred G17). Carry the guard that
  `mukai.from_chern`'s `ch₂+r` shift is **K3-only**.
- **Dependencies:** None hard (uses shipped `mukai.py`). Soft: G3, G5. Independent of G12. (Its solver
  is the shared engine deferred G17 would reuse.)

## G11 — Enriques / bielliptic + Koseki product-threefold catalog rows (record cheaply, defer faithful computation)
- **Tier:** 2 (the "record a row" sub-part is Tier-1-cheap; faithful computation is NS-gated → Tier 2)
- **Mathematical prerequisites:** Enriques: Nuer–Yoshioka arXiv:1901.04848; Yoshioka arXiv:1607.04946
  (`K` is 2-torsion). Bielliptic: arXiv:2107.13370 (Picard rank up to 4, with torsion). **Koseki,
  *Stability conditions on threefolds with nef tangent bundles*, arXiv:1811.03267** — BG (hence
  stability conditions) for `P¹×S`, `P²×C`, `P¹×P¹×C` with `S` abelian, `C` elliptic. **[UNVERIFIED —
  exact hypotheses on `S`, `C` MUST be checked against 1811.03267 before setting `bg_proven=True`.]**
  Distinct from Koseki's *weighted-hypersurface* / Calabi–Yau double–triple-solid paper arXiv:2007.00044
  — do not conflate (arXiv:1510.04474 is Chuang–Lai's *withdrawn* CY/Fano paper, not Koseki).
- **Code prerequisites:** Existing `varieties.py` `Surface`/`Threefold`. Soft: G1 (citation-fix task,
  same file) and G5 (`Certificate`). Faithful-computation half gated on G12.
- **Implementation sketch (revised — guard routing corrected):** Add catalog entries `enriques()`,
  `bielliptic()`, and threefold rows `koseki_product(...)` for `P¹×S`, `P²×C`, `P¹×P¹×C`. Add a
  `canonical_order: int` field to `Surface` (0 = numerically trivial K3/abelian, 2 = Enriques
  2-torsion, else the torsion order) so Enriques is distinguished from K3/abelian *without* the full NS
  refactor. These are **record-only** rows carrying correct `d`/`K_H`/`chi_O`, a `Certificate` (G5)
  with exact hypotheses + citations, and — for the threefolds — `bg_proven` set **only after** the
  Koseki hypothesis check. Add a `faithful_computation_supported` property returning `False` for these
  rows. **The guard must live in the Surface-consuming entry points** (`compute_walls`, `actual_walls`,
  the `dlp` functions, `bg_check`, and any surface-aware `k3_wall`/`abelian_wall` wrapper) — NOT in
  `numerical_wall(v, w, d)`, which receives only a `ChernChar` and an `int d`, never sees the `Surface`,
  and so *cannot* detect `surface.kind == "enriques"` (RE-VERIFIED against the shipped
  `numerical_wall` signature). Those Surface-aware APIs raise a clear "requires NS-lattice refactor
  (G12)" error when `faithful_computation_supported is False`; `numerical_wall` carries only a doc
  caveat.
- **Testability:** Enriques row: `chi_O=1`, `canonical_order=2`, `kind="enriques"` — distinguishable
  from K3 (`chi_O=2`, order 0) and abelian (`chi_O=0`, order 0). Koseki product rows: `bg_proven` is
  `True` **iff** the 1811.03267 hypotheses are met; assert the `note`/`Certificate` cites
  **arXiv:1811.03267** (not the weighted-hypersurface paper arXiv:2007.00044) and
  `faithful_computation_supported is False`. Guard test: a
  `compute_walls`/`actual_walls`/`dlp` call on these rows raises the G12-required error (and
  `numerical_wall` on the raw triple does not — the caveat lives in docs). **[RESEARCH]** the Koseki
  exact-hypothesis verification is a required literature sub-task; until done, product rows carry rigor
  `CONJECTURAL`/`UNKNOWN`, never `PROVEN`.
- **Research value:** Catalog completeness with honest, machine-readable provenance; surfaces exactly
  which known-in-the-literature varieties are not yet faithfully computable and why. The Koseki
  hypothesis check is a small, citable literature-verification contribution.
- **Risks and unknowns:** **[RESEARCH/UNVERIFIED]** the Koseki hypotheses. Faithful (not record-only)
  computation is entirely gated on G12. Low risk for the record part; the temptation to compute walls
  with the scalar model must be actively blocked (the Surface-consuming-API guard above). The Koseki
  `bg_proven=True` gate stays strictly behind the [UNVERIFIED] hypothesis check against arXiv:1811.03267.
- **Dependencies:** Soft: G1, G5. Faithful-computation half is gated on G12 (partial/soft — the primary
  record-row deliverable does NOT block on G12; see the DAG section).

## G12 — NS-lattice refactor: ch₁ scalar → NS lattice vector + symmetric intersection form
- **Tier:** 2 (HIGH-effort / strategic — a v2 architecture decision, the gate under every Picard-rank
  ≥ 2 item, and the **largest single schedule risk in the roadmap**)
- **Mathematical prerequisites:** Néron–Severi lattice with its symmetric intersection form; Hodge
  index theorem (signature `(1, ρ−1)`); Maciocia arXiv:1202.4587 `β = bω + γ` decomposition (`γ ⊥ ω`,
  `γ² ≤ 0`). No new theorem — this is the data model the existing formulas already implicitly require;
  the current scalar `c = ch₁·H` only sees the H-projection. **[PROVEN]** that it is the mathematically
  necessary input.
- **Code prerequisites:** None (foundational). **MUTATES** the frozen `ChernChar` in `chern.py`, plus
  `walls.py`, `mukai.py`, `bg_check.py`, and all 42 pinned tests.
- **Implementation sketch (revised — pin the backward-compat shim exactly):** Introduce a frozen
  `NSLattice(rank: int, gram: tuple[tuple[int,...],...])` with `pairing(u, v) -> Fraction` and a
  designated ample class `H` (integer vector); store it on `Surface`. Promote `ChernChar.c` from a
  scalar `Fraction` to `ch1`, an NS-vector (a length-ρ tuple bound to the surface's lattice); keep `e =
  ch₂` a `Fraction`. **Backward-compat shim (mandatory, and its encoding must be pinned bit-for-bit):**
  `ChernChar(r, c, e)` with scalar `c` must still work for Picard rank 1 and MUST reproduce the current
  formulas exactly, i.e. it must yield `⟨ch1, H⟩ = c` and `⟨ch1, ch1⟩ = c²/d` (so `discriminant = ½μ²
  − e/(rd)` and `bogomolov_discriminant = c²/d − 2re` are UNCHANGED — RE-VERIFIED against the shipped
  `chern.py`: `discriminant` and `bogomolov_discriminant` compute exactly these). **The naive encoding
  "`(c,)` with `Gram=[[d]]`, `H=(1,)` is WRONG and breaks the 42-test gate:** under `Gram=[[d]]` and
  `H=(1,)` you get `⟨ch1,H⟩ = c·1·d = c·d` and `⟨ch1,ch1⟩ = c²·d`, not `c` and `c²/d`, so `delta(1/2)`
  would no longer equal `5/8`. Instead store the H-coordinate `c/d` (or choose a basis where the ample
  generator has self-intersection `d` and `ch1 = (c/d)H`), so `⟨(c/d)H, H⟩ = (c/d)·d = c` and
  `⟨(c/d)H, (c/d)H⟩ = (c/d)²·d = c²/d`. Rewrite every formula that used `c·H`/`c²/d` in lattice terms:
  `ch₁·H = ⟨ch1, H⟩`, and the TRUE `ch₁² = ⟨ch1, ch1⟩` (replacing the `(c/d)²` Hodge-index proxy in
  `discriminant`/`bogomolov_discriminant`). `walls.numerical_wall` generalizes its `W_rc, W_re, W_ce`
  minors to the lattice pairing, and the ρ=1 path MUST reproduce the current minors bit-for-bit.
  `mukai.py`'s `d·l·l'` generalizes to `⟨·,·⟩` on the `l`-component.
- **Testability:** **Regression gate (ground truth): all 42 existing tests pass after migration** —
  P²[2] wall `(−5/2, 3/2)`, `delta(1/2)=5/8`, `delta(2/5)=13/25`, K3 `v(O)=(1,0,1)` `v²=−2`, P³
  `alpha_crit(β=1/2)=√3`, Gieseker wall `(−(2n+1)/2, (2n−1)/2)`. **Shim pin (new):** for the rank-1
  shim class, `⟨ch1,H⟩ == c` and `⟨ch1,ch1⟩ == c²/d` bit-for-bit (this is the assertion that protects
  the 42-test gate). New ρ=2 test — **P¹×P¹** with Gram `[[0,1],[1,0]]` (basis `f,s`), `H=f+s`,
  `⟨H,H⟩=2` reproduces `d=2`; two classes with `ch1=f+s` vs `ch1=2f` have equal `⟨ch1,H⟩` but different
  true `ch1²` (`⟨f+s,f+s⟩=2` vs `⟨2f,2f⟩=0`) and must be DISTINGUISHED (the old scalar model conflated
  them). **𝔽_n** Gram `[[0,1],[1,−n]]`: assert `det(Gram) < 0` (Hodge-index signature `(1,1)`).
- **Research value:** THE gate for all Picard-rank ≥ 2 work: rational surfaces (P¹×P¹, 𝔽_n, del
  Pezzo), products of curves (E₁×E₂), higher-rank K3s, mixed polarizations. Enables G13, G14, the
  faithful half of G11, and G18.
- **Risks and unknowns:** High effort; touches the frozen core and every test — risk of subtle
  convention drift during migration, mitigated by the 42-test regression gate and the explicit shim
  pin. **This is the dominant schedule risk in the roadmap** (upper-end Tier 2 / borderline strategic).
  **[SPECULATIVE]** on exact effort magnitude; decisively not a "medium" tweak. Must preserve the
  inviolables (exact `Fraction`, CH convention, frozen dataclasses, zero deps).
- **Dependencies:** None. Unlocks G13, G14, and the faithful half of G11; also unlocks G18.

## G13 — Maciocia ρ>1 full-NS surface walls + PROVEN boundedness certificate
- **Tier:** 2
- **Mathematical prerequisites:** Maciocia arXiv:1202.4587 — wall centers/radii via the `β = bω + γ`
  decomposition with `γ² = −d ≤ 0`; **Thm 3.11 boundedness** (all walls of a fixed class lie inside one
  bounding semicircle; wall count globally finite for all but one ray). Bertram nested-wall theorem;
  Hodge index. Arcara–Bertram–Lieblich arXiv:0708.2247 (Picard-rank-1 family, the reduction check).
- **Code prerequisites:** G12 — hard. Existing `walls.py`.
- **Implementation sketch:** With `NSLattice` available, add `numerical_wall_ns(v, w, surface) -> Wall`
  computing the semicircle from the full intersection form: center/radius depend on `⟨ch1,ch1⟩` and the
  polarization-orthogonal `γ²`, not merely `⟨ch1,H⟩`. Add `maciocia_wall_bound(v, surface) -> Fraction`
  (radius² of the outermost bounding semicircle, Thm 3.11), and `enumerate_ns_walls(v, surface)`
  enumerating the finitely many walls inside that bound. **Replace `actual_walls_complete`'s empirical
  doubling with this PROVEN bound** for the ρ>1 classes it covers — upgrading
  **[SPECULATIVE]→[PROVEN]**. (The ρ=1 boundedness upgrade is G7; this is the ρ>1 full-NS
  generalization.)
- **Testability:** **ρ=1 reduction (regression):** `numerical_wall_ns` on any Picard-rank-1 surface
  reproduces `numerical_wall` exactly — P²[2] → `(−5/2, 3/2)`; K3 walls. **ρ>1 gain:** on P¹×P¹, a
  destabilizer distinguished only by bidegree (same `⟨ch1,H⟩`) yields a wall that `numerical_wall_ns`
  finds but scalar `numerical_wall` MISSES — assert the presence/absence pair. **Boundedness:** every
  enumerated wall's `radius_sq ≤ maciocia_wall_bound(v)`; `actual_walls_complete` returns
  `complete=True` from the PROVEN bound. **[RESEARCH-light]** pin one exact ρ>1 wall value from a worked
  example in Maciocia 1202.4587 (his P¹×P¹ / abelian examples), read off the paper before pinning.
- **Research value:** Certified, finite surface wall enumeration for Picard rank ≥ 2 — products of
  curves, mixed polarizations, higher-ρ K3s — the first faithful surface walls beyond P²/ρ=1 in the
  toolkit.
- **Risks and unknowns:** Wall count grows with ρ and needs careful bounding (Thm 3.11 supplies it).
  **[RESEARCH-light]** transcribing a published ρ>1 example. Feasible but heavy; entirely gated on G12.
- **Dependencies:** Hard: G12. Feeds G18.

## G14 — Surface-specific exceptional-bundle generators for P¹×P¹ / 𝔽_n / del Pezzo (numerical layer)
- **Tier:** 2 (case-by-case per-surface cataloguing; the numerical Euler/Chern OUTPUT is
  in-architecture; the sheaf-level exceptionality *proof* is out of architecture and delegated to the
  CAS oracle G16). **Part of the single "rational-surface program" epic with G18 — one arc, not two
  independent wins.**
- **Mathematical prerequisites:** Existence of full (strong) exceptional collections: `P¹×P¹` has the
  Beilinson-type `⟨O, O(1,0), O(0,1), O(1,1)⟩` **[PROVEN]**; `𝔽_n` (Rudakov) and del Pezzo
  (Kuleshov–Orlov) have their own **[PROVEN]** — but **case by case, with NO uniform Markov-style
  recursion** (unlike P²'s ε/Markov miracle; Coskun–Huizenga weak Brill–Noether arXiv:1611.02674;
  Levine–Zhang arXiv:1910.14060). K-theory rank of a surface `= 2 + ρ`, so collection length is `4` on
  `P¹×P¹`/`𝔽_n` (ρ=2) and `12 − deg` on a del Pezzo of degree `deg` (e.g. length 9 on a cubic surface)
  **[PROVEN]**. Euler pairing via Riemann–Roch: `χ(E,F) = ∫ ch(E)^∨·ch(F)·td(X)`, exactly computable
  from `(r, c₁∈NS, ch₂)` and the intersection form — the Chern-numerical shadow of exceptionality
  **[PROVEN]**.
- **Code prerequisites:** G12 — hard (needs the NS-vector `ch₁`; `O(a,b)` bidegree lives in a rank-2 NS
  lattice). Existing `exceptional.py` (P²-only Markov machinery retained but not reused). Optionally G16
  (M2 oracle) to *verify* exceptionality at the sheaf level (`Ext^{>0}(E,E)=0`).
- **Implementation sketch:** Add `exceptional_surface.py` that, given a surface's `NSLattice` +
  intersection form + `K_X` (from G12), builds exceptional collections as `SurfaceBundle(r, c1:
  NSVector, ch2: Fraction)` objects and a Riemann–Roch `chi(E, F, surface) -> int`. Hard-code the
  *known* generators per family (no reusable recursion): `p1xp1_collection()` returns the 4-object
  list; `hirzebruch_collection(n)` the Rudakov `𝔽_n` collection; `del_pezzo_collection(deg)` the
  Kuleshov–Orlov collection. The uniform deliverable is the **Euler-form Gram matrix** `G[i][j] =
  chi(E_i, E_j, surface)` and the check `is_exceptional_collection(list, surface)` = "`G`
  upper-triangular unipotent with unit diagonal". **State explicitly that
  `is_exceptional_collection` is a NECESSARY numerical condition only:** genuine exceptionality
  (`Ext^{>0}(E_i,E_i)=0` plus Hom-vanishing) is sheaf-level and must be delegated to G16, **never
  asserted by the core**. Deliberately ship *only* the numerical layer — the concrete meaning of "proof
  engine OUT of architecture; numerical decision OUTPUT in-architecture".
- **Testability:** **P¹×P¹ Euler Gram** in basis `(O, O(1,0), O(0,1), O(1,1))` equals the exact matrix
  `[[1,2,2,4],[0,1,0,2],[0,0,1,2],[0,0,0,1]]` (from `χ(O(a,b),O(c,d)) = (c−a+1)(d−b+1)`); upper-
  triangular unipotent ⟹ `is_exceptional_collection` True (necessary condition). Diagonal:
  `chi(E,E,surface) == 1` for every generator. Off-diagonal witness: `χ(O(1,0),O(0,1)) == 0` and
  `χ(O(0,1),O(1,0)) == 0`. Length invariants: `len(p1xp1_collection()) == 4`;
  `len(del_pezzo_collection(3)) == 9` (`= 12 − 3`). (If G16 present) M2 confirms `Ext^{>0}(E_i,E_i)=0`
  for each generator — the ONLY route to certifying exceptionality.
- **Research value:** The numerical backbone of rational-surface non-emptiness — the exceptional-
  collection DATA + Euler-form table feeds G18 (Coskun–Huizenga `δ_H`); a citable exact table no CAS
  packages at the numerical level.
- **Risks and unknowns:** **[RESEARCH]** no uniform recursion — each family is hand-catalogued and
  pinned against Rudakov / Kuleshov–Orlov; correctness is per-surface literature verification. Del Pezzo
  `deg ≤ 2` / non-anticanonical polarizations only partially understood. Exceptionality itself is
  unprovable at Chern level (needs the sheaf-level oracle G16); the core must not claim to certify it.
  Gated entirely on G12.
- **Dependencies:** Hard: G12. Optional: G16 (verification). Feeds G18. (Epic-paired with G18.)

---

# Tier 3 — new math / computational infrastructure; may depend on open problems

Ordered by value-per-effort. G15 is the lightest (a labeling layer) and is Tier 3 only because its
non-`bg_proven` verdicts depend on the OPEN threefold BG conjecture — its **code is Tier-1-light and it
should be built as soon as G9 and G5 land**, not artificially delayed behind the heavier Tier-3 goals.
G18 is the heaviest. G16 (M2) is the DEFAULT CAS oracle. (G17 and G19 have been DEFERRED — see the
appendix.)

## G15 — Conjecture-gated Bridgeland-wall verdicts on threefolds
- **Tier:** 3 (epistemic, not effort — code is light; Tier 3 because correctness on non-`bg_proven`
  rows depends on the OPEN threefold BG conjecture and a *conditional* moduli theorem. **Schedule
  immediately after G9/G5.**)
- **Mathematical prerequisites:** The tilt→Bridgeland upgrade: a numerical tilt wall in `(β,α)` is a
  genuine *Bridgeland* wall only where the BMT/BG inequality `Q ≥ 0` holds. **[PROVEN for the catalog
  rows]** (BMT arXiv:1103.5010; the `bg_proven=True` rows: P³ 1103.5010/1207.4980; Fano ρ=1 Li
  1510.04089; abelian 1410.1585; quintic 1810.03434; quadric via BMSZ 1607.08199). **[CONJECTURED in
  general]**; strong form **[PROVEN]** false on `Bl_p(P³)` (Schmidt 1602.05055). Conditional moduli
  properness: Piyaratne–Toda arXiv:1504.01177 (proper of finite type **assuming** BG). Schmidt
  arXiv:1509.04608 for concrete P³ first-wall examples (twisted cubics / complete intersections).
- **Code prerequisites:** G9 (the threefold tilt-wall solver — hard). G5 (the `Rigor`/`Certificate`
  lattice — hard; G15 is its first real cross-pipeline consumer). Existing `threefold.py`,
  `varieties.Threefold.bg_proven`.
- **Implementation sketch:** A thin certification wrapper, not a new solver.
  `bridgeland_wall_verdict(v, w, threefold, d3) -> BridgelandWallVerdict` calls `numerical_tilt_wall`/
  `bridgeland_wall` (G9) to get the locus, then attaches a `Certificate` (G5) whose rigor is
  `min(TILT_SOLVER_RIGOR, PROVEN if threefold.bg_proven else CONJECTURAL)`. The dataclass
  `BridgelandWallVerdict(tilt_wall, certified: bool, rigor: Rigor, citations, note)` records
  `certified = threefold.bg_proven` and, when `False`, the note "tilt-stability wall only; Bridgeland
  upgrade unproven (threefold BG open here)". For `Bl_p(P³)` the note must ADDITIONALLY state "strong
  BMT boundary FAILS (Schmidt 1602.05055); stability conditions nonetheless EXIST (Fano, BMSZ
  1607.08199)" — i.e. NEVER emit "no stability condition". **NO code path may emit `certified=True` when
  `bg_proven` is False, nor emit "no stability condition exists".** Route
  `is_bridgeland_certified(threefold)` through `bg_boundary_curve` / any wall enumerator so every
  emitted threefold wall carries its `Certificate`. No `Fraction` value changes — this goal only adds
  the provenance label that distinguishes a theorem from a conjecture on identical arithmetic.
- **Testability:** Propagation invariant (exact, no open math): for every `bg_proven=True` row
  `bridgeland_wall_verdict(...).rigor == PROVEN` and `.certified is True`; for `BLOWUP_P3_POINT`
  `.rigor <= CONJECTURAL` and `.certified is False`. Meet law: `verdict.rigor == min(tilt_rigor, PROVEN
  if bg_proven else CONJECTURAL)`. Ground truth: a Schmidt 1509.04608 first-wall example on P³ yields
  the same `(β,α)` locus as G9 AND is labeled `PROVEN`; the identical class on `Bl_p(P³)` yields
  `CONJECTURAL`/tilt-only. Negative test: NO code path returns `certified=True` when `bg_proven is
  False`.
- **Research value:** Threefold-stability researchers get machine-checked honesty about *which* wall
  statements are theorems — operationalizing the correctness culture at the frontier where the same
  arithmetic is a theorem on P³ and an open conjecture on a general threefold.
- **Risks and unknowns:** **[RESEARCH]** Entirely gated on G9. The Piyaratne–Toda properness statement
  is *conditional* ("assuming BG"); the citation must carry that hypothesis even on `bg_proven` rows. No
  risk of regressing pinned values (pure additive labeling). Its Tier-3 label is epistemic — do not
  schedule it after the heavier Tier-3 goals; build it right after G9/G5.
- **Dependencies:** Hard: G9, G5. No cycles.

## G16 — Optional Macaulay2 subprocess oracle (sheaf-level Ext/moduli + formula-layer cross-check)
- **Tier:** 3 (significant new computational infrastructure: an external CAS bridge reaching sheaf-level
  data the numerical core forbids by construction)
- **Mathematical prerequisites:** `Ext^i(E,F)` on a projective scheme via truncated graded-module Ext
  (Greg Smith, arXiv:math/9807170) **[PROVEN]**; sheaf cohomology on projective/toric X (M2
  `Varieties`, `NormalToricVarieties`) **[PROVEN + implemented]**. Euler-pairing consistency `χ(E,F) =
  Σ(−1)^i ext^i(E,F)` is the only Chern-numerical shadow of the individual `ext^i` **[PROVEN]** — the
  exact quantity the core already computes (`exceptional.chi` on P², `mukai.pairing = −χ` on K3), a
  rigorous cross-check target. Lossless interchange: M2 rationals are `QQ`; `Fraction("p/q") ↔ p/q`
  round-trips exactly as text **[PROVEN — both sides exact ℚ]**.
- **Code prerequisites:** None in the math core. Depends only on an **external, non-pip M2 install** the
  user provisions. Optionally consumes G5 to tag oracle-verified results.
- **Implementation sketch:** Ship a never-core subpackage `bridgeland_stability/oracle/m2.py`, imported
  lazily exactly like `viz/` (so `import bridgeland_stability` stays zero-dependency). A context-managed
  `class M2Session` spawns `M2 --script` (or a `subprocess`/`pexpect` pipe), with `require_m2()`
  mirroring `viz.style.require_mpl()`. A generator layer emits `.m2` scripts from Chern data:
  `chi_via_ext(E, F, X) -> int` builds the sheaves, calls `Ext^i(F,G)` for `i=0..dim`, parses printed
  dimensions; `ext_dims(E, F, X) -> tuple[int,...]` returns the full `(ext^0,…,ext^n)`;
  `moduli_nonempty_by_construction(r,c1,ch2,X) -> Optional[bool]` attempts an explicit sheaf/monad and
  reports success/obstruction. A hand-maintained regex parser turns `QQ`/`ZZ` output back into
  `Fraction`/`int`. Two use-cases: (a) independent oracle — recompute a `χ`/`Ext^1`/cohomology
  dimension in M2 and diff against the formula layer (feeding `docs/CORRECTIONS.md`'s two-way
  discipline); (b) unlock sheaf-level questions (actual `ext^i`, moduli non-emptiness *by
  construction*). This is the single most architecturally-aligned CAS bridge; it must NOT try to
  reimplement M2 or become a general derived-category layer (none exists anywhere).
- **Testability (two illustrative values CORRECTED per the amendment):** `chi_via_ext(O(0), O(n), P2)`
  equals `exceptional.chi(...)` and the closed form `(n+1)(n+2)/2` for `n ∈ {−3,…,5}` — RE-VERIFIED
  against the shipped `exceptional.chi`: `χ(O,O)=1`, `χ(O,O(1))=3`, **`χ(O,O(3))=10`** (the shipped and
  closed-form value; the earlier draft's `15` is the `n=4` value `χ(O,O(4))`). `ext_dims(O(0), O(-3),
  P2) == (0,0,1)` (only `Ext²=H²(O(-3))≅ℂ` by Serre duality, `K_{P²}=O(-3)`), and its **alternating sum
  is `+1`** (`= (−1)^0·0 − (−1)^1·0 + (−1)^2·1 = +1`), matching the shipped `exceptional.chi(O,O(-3)) =
  +1` (RE-VERIFIED: `deg2=9/2, deg1=−3, deg0=1 ⟹ 9/2 − 9/2 + 1 = 1`). Both wrong draft values (`15`,
  `−1`) would have failed the very cross-check the oracle exists to provide. On a K3 with `Pic=ℤH`,
  `−ext`-alternating-sum equals `mukai.pairing` for a pinned pair (`⟨v(O),v(O)⟩ = −2`). Round-trip:
  `Fraction(-7,6)` sent to M2 and read back is exactly `Fraction(-7,6)`. Graceful-skip: with M2 absent,
  `import bridgeland_stability` still succeeds with zero deps and `oracle` raises a clear guard (mirrors
  the viz lazy-import culture).
- **Research value:** Anyone auditing the formula layer gets an independent, exact second opinion;
  researchers needing a *specific* `Ext^1`/obstruction or constructive non-emptiness witness get it
  without leaving Python. The enabling oracle for the rational-surface program (G18) where the
  HN/prioritary core is out of architecture.
- **Risks and unknowns:** External M2 install + hand-maintained text parser are a real (non-pip)
  maintenance burden; M2 print-format drift can break the parser. Windows: M2 runs natively or under WSL
  **[UNVERIFIED — confirm the host's M2 path before relying on it]**. `moduli_nonempty_by_construction`
  is only a *sufficient* witness (a construction failure is not a non-existence proof) — label it so. No
  risk to the core (never imported by it).
- **Dependencies:** None (standalone bridge). Optional: G5 for tagging. Enables deferred G17, and is the
  default (soft) supplier for G18's non-Chern HN step.

## G18 — Coskun–Huizenga `δ_H` non-emptiness decision procedure on rational surfaces
- **Tier:** 3 (multi-month, per-surface; its algorithmic CORE is a Harder–Narasimhan / prioritary-sheaf
  computation that is NOT pure Chern arithmetic — only the numerical *verdict* is in-architecture).
  **Together with G14 this is the single "rational-surface program" epic — one arc, not two independent
  wins.**
- **Mathematical prerequisites:** Coskun–Huizenga arXiv:1907.06739 (existence on Hirzebruch surfaces
  for *any* ample divisor): pass to prioritary-sheaf stacks, compute the HN filtration of the generic
  prioritary sheaf; semistable sheaves exist **iff that HN filtration has length one**; sharp bounds `Δ
  ≥ δ_H(ξ)` depending on the polarization and the *full* (two-coordinate) slope **[PROVEN]**.
  Levine–Zhang arXiv:1910.14060 (del Pezzo deg ≥ 3 existence, deg ≥ 4 cohomology, anticanonical)
  **[PROVEN]**; Coskun–Huizenga weak Brill–Noether arXiv:1611.02674 **[PROVEN]**. Crucial scope fact:
  unlike P²'s closed-form DLP curve, there is **no closed `δ`-curve** — the answer is a decision
  procedure with a sheaf-theoretic (HN/prioritary) core.
- **Code prerequisites:** G12 (NS-lattice refactor — two slope coordinates; the scalar `c=ch₁·H` cannot
  express the input). G13 (Maciocia full-NS walls — the sharp `δ_H` interacts with the wall structure
  via the true `ch₁²`). G14 (the surface's exceptional-bundle generators — the bounds are stated through
  the exceptional collection). A non-Chern HN/prioritary step supplied *out of architecture* via **G16
  (M2, default) OR a hand-tabulated finite class list from the papers OR deferred G19 (OSCAR)** — this
  supplier is a **SOFT** prerequisite (see the DAG amendment), so G18's in-architecture slice is NOT
  gated on an optional external install.
- **Implementation sketch:** Add `nonemptiness_rational.py` exposing `delta_H(xi, surface) -> Fraction`
  (the sharp bound for `ξ=(r,c₁∈NS,ch₂)` w.r.t. `surface.H`) and `moduli_nonempty(r, c1, ch2, surface)
  -> Verdict`. Two-layer split: (1) an **in-architecture numerical layer** that, *given* the
  HN-length-one datum, evaluates `Δ(ξ) ≥ δ_H(ξ)` exactly in `Fraction` and returns yes/no with a G5
  `Certificate`; (2) an **out-of-architecture HN/prioritary layer** — the generic-prioritary-sheaf HN
  filtration — either *delegated* to the M2/OSCAR oracle (G16/deferred G19) which constructs the
  prioritary sheaf and reads off its HN length, or *tabulated* from the paper (1907.06739 /
  1910.14060) for a fixed finite class list. `moduli_nonempty` must clearly report which mode produced
  the verdict. A verdict from the in-architecture inequality carries `rigor=PROVEN` only when the
  HN-length-one hypothesis was supplied by a certified source (oracle or paper table); otherwise
  `HEURISTIC`. The shippable in-architecture slice (evaluate `Δ ≥ δ_H` given HN-length-one) is **thin**;
  the HN/prioritary core is sheaf-theoretic and delegated. This is the sharpest test of the "numeric
  output in, sheaf engine out" discipline.
- **Testability:** Reproduce a published existence verdict: for a table of `(r,c₁,ch₂)` classes worked
  in 1907.06739 on `P¹×P¹`/`𝔽_n` (and Levine–Zhang del Pezzo deg ≥ 3, anticanonical),
  `moduli_nonempty(...)` returns the paper's yes/no. **[RESEARCH — exact `δ_H` numeric targets must be
  extracted from 1907.06739/1910.14060; these are the only genuinely paper-dependent test values in the
  roadmap.]** Sanity floor: every returned bound satisfies `δ_H(ξ) ≥ 0` and, in the P² limit /
  single-ray case, degrades to the existing DLP `delta(μ)` (regression against pinned `delta(1/2)=5/8`,
  `delta(1/3)=5/9`). Cross-check (G16): M2's `moduli_nonempty_by_construction` agrees with the
  formula-layer verdict on a shared class list. Polarization-dependence witness: on `𝔽_n`, two ample `H`
  in different chambers give *different* `δ_H` for a fixed `ξ` (the dependence P² cannot see).
- **Research value:** Would make `hirzebruch(n)`/`P¹×P¹`/del Pezzo *first-class* for non-emptiness — the
  most-requested "beyond P²" capability — delivering exact, citable numerical verdicts (with an
  oracle-backed HN certificate) that no packaged tool currently provides.
- **Risks and unknowns:** **[RESEARCH]** The HN/prioritary core is genuinely sheaf-theoretic and outside
  the zero-sheaf constraint; only its output is portable — and the in-architecture slice that ships
  without an oracle/paper table is thin. Gated on G12 **and** G13 **and** G14 — the deepest hard
  dependency chain in the roadmap; not a near-term deliverable. Extracting closed `δ_H` numeric targets
  requires reading 1907.06739 carefully (**[UNVERIFIED]** at the exact-value level). Non-anticanonical /
  low-degree del Pezzo only partially covered by theory.
- **Dependencies:** Hard: G12, G13, G14. Soft: G16 (M2 default) OR paper-tabulation OR deferred G19 —
  the non-Chern HN step; NOT a hard gate. No cycles. (G18 + G14 together = the single "rational-surface
  program" epic.)

---

## Dependency DAG

The union of all dependencies is acyclic. **Hard edges** are blocking prerequisites (they define the
DAG and the topological order); **soft edges** are recommended, non-blocking precedence (provenance
tagging, build-early, or partial/optional integration) documented for scheduling but not required for a
goal to ship its primary deliverable.

**Hard edges (Gx → Gy: Gx must precede Gy):**
- G2 → G8   (abelian free win is documented by the worked example)
- G3 → G8   (K3 shim is documented by the worked example)
- G4 → G9   (ν tilt-slope is the primitive the threefold tilt-wall solver is built on)
- G5 → G15  (the Rigor/Certificate lattice is the substrate G15 propagates)
- G9 → G15  (G15 decorates G9's tilt-wall locus with a Bridgeland certificate)
- G12 → G13 (Maciocia ρ>1 walls need the NS lattice)
- G12 → G14 (exceptional bundles on ρ>1 surfaces need the NS-vector ch₁)
- G12 → G18 (δ_H needs two NS slope coordinates)
- G13 → G18 (the sharp δ_H interacts with the full-NS wall structure)
- G14 → G18 (the δ_H bounds are stated through the exceptional collection)

That is **10 hard edges**, all pointing forward in the `G1..G19` ordering ⟹ the DAG is acyclic and the
topological order is valid.

**Soft edges (recommended, non-blocking):**
- **G16 → G18 (RECLASSIFIED HARD→SOFT):** G18's own text supplies a non-G16 path for the
  HN-length-one datum (*tabulated from the paper for a fixed finite class list*), so G18 is **not**
  gated on an optional external Macaulay2 install. The HN supplier is *G16 OR paper-tabulation OR
  deferred G19*.
- G5 → G6, G5 → G7, G5 → G9, G5 → G10, G5 → G11, G5 → G16, G5 → G18, G5 → G17(deferred)
  (Certificate tagging; build G5 EARLY so every downstream goal can attach certificates from day one)
- G1 → G11 (G11 reuses the citation-fix task and the same `varieties.py`)
- G3 → G10 (G10 links its certified classification to G3's (s,t) semicircle)
- G12 → G11 (only the *faithful-computation* half of G11; the primary record-row deliverable does NOT
  block on G12)
- G16 → G14 (optional sheaf-level verification of exceptionality, `Ext^{>0}(E,E)=0`)
- G16 → G17(deferred) (optional cross-check of a fake-wall verdict against a constructed moduli space)
- G10 → G17(deferred) (G17 reuses G10's `solve_binary_quadratic`; G17 is a possible follow-on of G10,
  NOT an independent win — see below)
- G19(deferred) → G18 (OSCAR as an alternative HN-step supplier to G16)

**Scheduling notes (from the DAG review):**
- **G15 is epistemically Tier 3 but its code is Tier-1-light.** Its only hard prerequisites are G9 and
  G5 (both buildable early); the Tier-3 label reflects dependence on the open threefold BG conjecture
  for non-`bg_proven` rows, not effort. **Schedule G15 immediately after G9/G5**, ahead of the heavier
  Tier-3 goals G16/G18 — do not artificially delay it.
- **G17 substantially overlaps G10.** Both implement the hyperbolic binary-quadratic solver and
  certified witness enumeration; G17's only net-new deliverable (the fake/totally-semistable decision)
  is exactly the part the goal admits may not close. G17 is therefore DEFERRED (appendix) and, if
  revived, should be a follow-on of G10, not scheduled as an independent win.

**Topological-ish ordering (respects tiers and value-per-effort):**

```
G1 → G2 → G3 → G4 → G5 → G6 → G7 → G8      (Tier 1)
→ G9 → G10 → G11 → G12 → G13 → G14          (Tier 2; G12 must precede G13/G14)
→ G15 → G16 → G18                            (Tier 3 active; schedule G15 right after G9/G5;
                                              G16 or paper-tabulation before G18)
[deferred: G17 (follow-on of G10), G19 (redundant with G16)]
```

Every hard edge points forward in this ordering, confirming acyclicity. Independent roots (no hard
prerequisites): G1, G2, G3, G4, G5, G6, G7, G10, G11, G12, G16. Build-early recommendation: **G1** (the
only confirmed shipped defect) and **G5** (the provenance substrate consumed everywhere) first, in
parallel with the Tier-1 math wins.

**Dedup / overlap notes.** The surface-specific exceptional-bundle generators are the single Tier-2
goal **G14** (the numerical Euler/Chern layer is in-architecture and testable via the exact P¹×P¹ Euler
Gram; the sheaf-level exceptionality *proof* is out of architecture, optionally verified by G16). G10
and deferred G17 share the hyperbolic binary-quadratic solver (`solve_binary_quadratic` /
`hyperbolic_witnesses`) — implement once in `mukai.py` (in G10) and reuse. G16 (Macaulay2) is the
DEFAULT CAS oracle; deferred G19 (OSCAR) is escalation-only and overlaps G16. The G2/G3 shared fixture
`v=(1,0,-2), w=(1,-1,1/2), d=2` (abelian `17/4` vs K3 `21/4`, differing by exactly `2/d`) is the
project-wide regression guard against wrongly applying the K3 `+r` shift to abelian input — keep it in
one shared test module (reused by G8). Note this fixture's `l = −1/2` is non-integral, so it is a
synthetic `+2/d` demonstration; G3 additionally pins an integral-`l` genuine K3 lattice class.

## Pipeline-blocking goals

These goals gate the largest downstream fan-out; sequence them deliberately.

- **G12 (NS-lattice refactor)** — the dominant effort/schedule risk. It gates the entire
  rational-surface arc: G13, G14, the faithful half of G11, and G18. Its backward-compat shim must be
  pinned so the rho=1 path reproduces `⟨ch1,H⟩=c` and `⟨ch1,ch1⟩=c²/d` bit-for-bit (the 42-test
  regression gate). Treat as a v2 architecture decision, not an increment.
- **G5 (Rigor/Certificate lattice)** — consumed (soft) by nearly every downstream goal and hard by G15;
  build it EARLY so certificates attach from day one.
- **G4 (ν tilt-slope)** — the scalar primitive that hard-gates the entire threefold arc (G9 → G15).
- **G9 (threefold tilt-wall solver)** — hard-gates G15; the first faithful threefold (α,β) walls.

## Deferred / Abandoned goals

No goals were abandoned. Two were **deferred** (kept, with full fields preserved, but moved out of the
active tiers).

### G17 — Certified K3 totally-semistable / "fake" wall stratification (beyond the bounded search) — DEFERRED
- **Deferred:** the certified hyperbolic-witness enumeration this goal needs is already delivered by G10
  `solve_binary_quadratic`; G17's only net-new content is the totally-semistable / fake-vs-no-wall
  decision, which the goal itself flags may require the full Bayer–Macrì movable-cone stratification
  (beyond the rank-2 lattice) and may stay HEURISTIC. Revisit only after G10 ships and a concrete
  fake-wall classification need arises; if revived, build it as a follow-on of G10, not an independent
  win.
- **Tier:** 3 (upgrading the bounded-box `fake-or-none` verdict to a *certified* classification needs
  the full Bayer–Macrì totally-semistable stratification, which may exceed pure lattice arithmetic)
- **Mathematical prerequisites:** Bayer–Macrì arXiv:1301.6968: Thm 5.7 (divisorial / flopping / fake)
  via spherical (`s²=−2`) and isotropic (`w²=0`) classes in `H_w = ⟨v,w⟩` **[PROVEN — already in
  `mukai.classify_wall`]**; §5–6 the *totally semistable* ("fake") stratum. A *provable* search bound:
  on the hyperbolic plane of discriminant `−D`, the `{−2,0}`-solutions form a Pell-type binary-
  quadratic-form set whose primitive representatives are bounded by an explicit function of
  `v²,w²,(v,w),d`. **[SPECULATIVE that a clean closed bound suffices; PROVEN that the lattice is 2-dim
  hyperbolic so the solution set is Pell-structured.]**
- **Code prerequisites:** Existing `mukai.py`. **Reuse G10's `solve_binary_quadratic` — do not
  re-implement.** Optionally G16 (M2 oracle) to cross-check a stratification verdict against a
  constructed moduli space. Optionally G5 to tag `PROVEN` vs `HEURISTIC`.
- **Implementation sketch:** Replace the `search=8` box in `classify_wall` with a certified enumerator.
  Add `hyperbolic_witnesses(v, w, d) -> WitnessSet` solving `⟨x,x⟩ ∈ {−2,0}` on the saturated `H_w`
  exactly (the form `a²v² + 2ab(v,w) + b²w²` is an integral binary quadratic form of discriminant `Δ_H
  = 4((v,w)² − v²w²)`; enumerate `{−2,0}`-representatives via reduction theory / Pell fundamental
  solutions, returning a *proved-complete* finite orbit list); feed it to the existing Thm-5.7 branch
  logic. Add `fake_wall_certificate(v, w, d) -> Stratification` that, when no divisorial/flopping
  witness exists, decides "totally semistable (fake wall)" vs "not a wall" using the Bayer–Macrì
  totally-semistable criterion. Output carries a G5 `Certificate`: `PROVEN` when the witness enumeration
  is certified-exhaustive and the criterion is decisive; `HEURISTIC` where the stratification needs more
  than lattice data.
- **Testability:** Certified enumeration agrees with the box on all pinned `test_mukai.py` cases but is
  now bound-independent: `classify_wall(v,w,d)` returns the SAME verdict for `search ∈ {8, 32}` AND the
  enumerator asserts completeness. Spherical anchor: `v(O_{K3})=(1,0,1)` has `⟨v,v⟩=−2`
  (`is_spherical` True); the Brill–Noether divisorial branch fires for a spherical `s` with `(s,v)=0`.
  Hyperbolic-form invariant: for a pinned `(v,w)` the enumerated `{−2}`-witness count equals the number
  of Pell fundamental-solution orbits of the binary form of discriminant `Δ_H` (a checkable integer).
  Cross-check (if G16 present): a `PROVEN` fake verdict is consistent with M2 finding `M_σ(v)` nonempty
  but with strictly-semistable generic object.
- **Research value:** K3 wall-crossing / MMP researchers get a *certified* wall type instead of a
  search-box guess — the difference between a publishable classification and a heuristic. Removes a
  standing **[SPECULATIVE]** flag on `classify_wall`'s completeness.
- **Risks and unknowns:** **[RESEARCH]** The totally-semistable stratum is the hard part: deciding "fake
  vs no wall" may need more than the 2-dim lattice (the full movable-cone stratification), in which case
  the honest deliverable is *certified witness enumeration* (which duplicates G10, upgrading
  **[SPECULATIVE]→[PROVEN]** for the divisorial/flopping decision) while the fake/no-wall split stays
  `HEURISTIC`. Do not overclaim. Picard-rank-≥2 K3s remain out of scope (need G12).
- **Dependencies:** Hard: none (builds on shipped `mukai.py` + G10's solver). Optional: G16, G5. No
  cycles.

### G19 — OSCAR-via-`juliacall` research escalation (toric sheaf cohomology + mutation experiments) — DEFERRED
- **Deferred:** redundant with G16 as a sheaf-level oracle; native-Windows OSCAR support on this host is
  **[UNVERIFIED]** and historically WSL-only, with a multi-GB Julia/OSCAR install. Build only if the
  roadmap takes on a concrete toric coherent-sheaf-cohomology or exceptional-collection/mutation
  experiment that Macaulay2 (G16) cannot serve as conveniently; do not build as a second verification
  oracle. Never position it as a route to derived categories — no system provides a general
  derived-category-of-`Coh(X)` engine.
- **Tier:** 3 (heaviest infrastructure: a Julia toolchain + multi-GB OSCAR, on a Windows host where the
  path is **[UNVERIFIED]** and typically WSL-only; a separate research track, not a bolt-on. Most
  droppable.)
- **Mathematical prerequisites:** Toric coherent-sheaf / line-bundle cohomology via Cox-ring module
  cohomology (Eisenbud–Mustață–Stillman, arXiv:math/0001159) **[PROVEN + implemented in OSCAR/M2]**;
  *Toric geometry in OSCAR* arXiv:2303.08110. Mutation experiments: OSCAR homological algebra
  (`free_resolution`, `hom`, `ext`) + the distributed M2 mutation tool (Brown–Dey–Li–Sayrafi,
  arXiv:2509.25103) and the `QuiversToricVarieties` database (Prabhu-Naik, arXiv:1501.05861). Hard
  ceiling: **no system provides a general derived-category-of-`Coh(X)` engine** **[PROVEN — established
  negative in the survey]**.
- **Code prerequisites:** None in-core. An optional `pip install ".[cas]"` extra (`juliacall`) plus an
  external Julia + OSCAR install. Overlaps in purpose with G16 (both are sheaf-level oracles).
- **Implementation sketch:** Add a second never-core oracle backend `oracle/oscar.py`, lazily importing
  `juliacall` and running `jl.seval("using Oscar")`, behind `require_oscar()` (mirrors `require_mpl`).
  Expose `toric_line_bundle_cohomology(X, divisor) -> tuple[int,...]`, `toric_sheaf_chi(...)`, and
  experimental `exceptional_collection(X)` / `mutate(collection, i)` wrappers (and/or shell out to the
  M2 2509.25103 tool). Serialize `Fraction ↔ QQFieldElem` losslessly as `"p/q"`. **LEAN DEFER:**
  justified only if the roadmap explicitly moves into toric coherent-sheaf cohomology / mutation
  experiments M2 does not cover as conveniently. For the project's actual need (a sheaf-level
  cross-check oracle), G16 is cheaper and more aligned — build G16 first; add OSCAR only if a concrete
  toric/mutation experiment demands it. Do NOT position this as a route to "derived categories" — that
  engine exists nowhere.
- **Testability:** Toric cohomology closed-form checks: on `P²`, `h⁰(O(n)) = C(n+2,2)`, `hⁱ=0` for
  `0<i<2`, `h²(O(n))=C(−n−1,2)`; on `P¹×P¹`, `h⁰(O(a,b))=(a+1)(b+1)` for `a,b≥0` — OSCAR output must
  match these exact integers and agree with G14's `chi`. `toric_sheaf_chi` on `P²` line bundles equals
  `exceptional.chi`. Mutation: mutating `⟨O,O(1),O(2)⟩` on `P²` reproduces a known exceptional
  collection (cross-checked against the M2 2509.25103 tool). Graceful-skip + lossless round-trip tests
  as in G16.
- **Research value:** Toric-geometry / exceptional-collection researchers get scripted access to toric
  sheaf cohomology and mutation experiments from Python, plus a *second independent* CAS oracle (OSCAR
  vs M2) for the highest-stakes cross-checks.
- **Risks and unknowns:** **[UNVERIFIED]** Native-Windows OSCAR support is historically limited
  (typically WSL) — re-check current status on this host before committing. Multi-GB install +
  minutes-long first `using Oscar` precompile. Mutation tooling is research-stage; no spherical-twist
  capability **[UNVERIFIED]**. Redundant with G16 for pure verification — avoid building both purely as
  oracles.
- **Dependencies:** None in-core (standalone, optional). Overlaps G16. Can supply the non-Chern HN step
  for G18 as an alternative to G16 (soft). No cycles.

---

## Adversarial review notes

**Verdict: proceed-with-amendments.**

The goal set is fundamentally sound and executable by a solo maintainer while respecting the
correctness culture. Every Tier-1 anchor was re-derived and re-run against the shipped code: G2
(abelian `17/4`), G3 (K3 `21/4`, difference `2/d = 1`), G4 (`ν = −1, −1/2, 0, None`; `twist(1,1).a3 =
−4/3`), G6 (the four existence identities), and G16 (the two illustrative Euler-pairing values) all
reproduce exactly. The DAG is acyclic with a coherent topological order, and the one confirmed shipped
defect (G1 miscitations) is correctly prioritized first.

**Amendments applied (six surgical corrections, no redesigns):**
- **G1** — clarified that the guard must assert *membership in a curated AG allowlist*, not merely
  arXiv-format regex (the format regex matches the probability-paper ID `1607.07182` perfectly); noted
  the allowlist is maintenance-brittle.
- **G2** — scoped the "free win" to `numerical_wall`/`abelian_wall` only; removed the false claim that
  `compute_walls`/`actual_walls` return genuine abelian walls (their P²[n]-only certification and the
  `_is_integral_chern` `d=1` filter do not transfer to abelian).
- **G6** — added a `c % d == 0` Picard-lattice guard before `l = c // d` (else `v=(2,1,0)` on `K3(2)`
  silently floors to a wrong `v²=−8`); encoded Yoshioka primitivity/genericity as
  `Certificate.hypotheses`.
- **G9** — corrected the citation: Schmidt Thm 3.3 justifies the TILT layer only (data `(ch₀,ch₁,ch₂)`,
  no `ch₃`); the second-tilt/Bridgeland (`ch₃`-dependent) semicircle gets its own reference from
  Schmidt's second-tilt construction and its `(y,z)` transcription stays `[RESEARCH-light]`.
- **G10** — enumerate over the SATURATED rank-2 sublattice `H_W`; emit `certified=True` only for
  positive (divisorial/flopping) witnesses; the solver returns one representative per automorphism
  orbit.
- **G11** — moved the `faithful_computation_supported` guard into the Surface-consuming APIs (not the
  surface-blind `numerical_wall(v,w,d)`); kept the Koseki `bg_proven` gate strictly behind the
  [UNVERIFIED] check against arXiv:1811.03267.
- **G12** — pinned the backward-compat shim so the rho=1 path reproduces `⟨ch1,H⟩=c` and
  `⟨ch1,ch1⟩=c²/d` bit-for-bit (the literal `(c,)`+`Gram=[[d]]` encoding is dimensionally wrong and
  breaks the 42-test gate); flagged G12 as the roadmap's largest schedule risk.
- **G16** — fixed the two wrong illustrative values: `χ(O,O(3))=10` (not 15; 15 is the `n=4` value) and
  the alternating sum of `ext_dims(O(0),O(-3),P2)=(0,0,1)` is `+1` (not −1), matching the shipped
  `exceptional.chi(O,O(-3))=+1`.
- Minor: G3 flagged its shared fixture as a synthetic `+2/d` demonstration and recommended an
  additional integral-`l` K3 pin; G4 pinned its values as exactly-derived from BMT (Schmidt agreement
  contingent on a normalization check); G14 stated `is_exceptional_collection` is a necessary numerical
  condition only; G15 flagged as Tier-1-light code to be scheduled right after G9/G5.

**DAG fixes applied:**
- Reclassified **G16 → G18 from HARD to SOFT** (the HN-length-one datum has a non-G16 path via
  paper-tabulation), so G18 is not gated on an optional external Macaulay2 install.
- Noted **G15**'s Tier-3 label is epistemic (open BG conjecture), not effort — scheduled right after
  G9/G5.
- Flagged the **G10/G17 overlap** and folded G17 into a possible follow-on of G10.

**Deferred (kept in the appendix, full fields preserved):**
- **G17** — its certifiable content duplicates G10's `solve_binary_quadratic`; only the possibly-
  unclosable fake/totally-semistable stratum is net-new. Revisit after G10 ships.
- **G19** — redundant with G16 as an oracle; native-Windows/WSL OSCAR support is [UNVERIFIED] on this
  host and the install is multi-GB. Build only for a concrete toric/mutation need M2 cannot serve.

**Nothing abandoned.** No goal hides an unsolved problem in a way that would silently produce wrong
values; the sheaf-level dependencies in G14/G16/G18 are explicitly delegated, not concealed. G12 is the
roadmap's dominant effort/schedule risk and gates the entire rational-surface arc.
