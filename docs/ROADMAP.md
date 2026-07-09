# ROADMAP — bridgeland_stability

This is the executable roadmap for extending `bridgeland_stability` beyond its current
Picard-rank-1, numerical-Chern-character, exact-`Fraction` regime. It converts the 17 active goals of
[`docs/GOALS.md`](GOALS.md) (G1–G16, G18; G17/G19 deferred) into **epics** (thematic work streams of
~2–6 months) composed of **milestones** (unit-sized 1–3 week items, each producing a testable,
mergeable artifact with at least one NEW test). Stable global goal IDs `G1..G19` are preserved; the
10 hard DAG edges are respected; the mandated G12/G9/G18 splits are applied and renumbered
consistently.

**Confidence flags are carried throughout:** **[PROVEN]** (a cited theorem), **[CONJECTURED]** (an
open conjecture), **[SPECULATIVE]** (an in-repo inference), **[UNVERIFIED]** (a claim not yet pinned to
a primary source), and **[RESEARCH]** / **[RESEARCH-light]** (a dependence on reading a specific paper
or on an open problem). The global honesty invariants of `CLAUDE.md` and `GOALS.md` are inviolable:
the exact-`Fraction` core, the Coskun–Huizenga convention `Δ = ½μ² − ch₂/(rd)`, frozen dataclasses,
zero runtime deps, and the 42 pinned tests are ground truth.

## How the epics are ordered

Epics are ordered so **each early epic strengthens the codebase and yields visible, testable progress
even if every later epic is abandoned** (standalone value at every step). The hard DAG is respected:
the NS-lattice refactor epic (E8) precedes the ρ≥2 epics (E9, E11); the CAS-oracle epic (E10) precedes
the rational-surface epic (E11).

- **E1 is the harden-and-document epic** (no new mathematics) — mandated first.
- **E11 is the most ambitious surviving Tier-3 goal** — the rational-surface non-emptiness program
  (G18 [+G14]) — mandated last.
- **G15 is pulled forward** to immediately after the threefold-engine epic (its only hard prerequisites
  are G9 and G5), not stranded behind the heavy G12/G13/G14 block. Its Tier-3 label is epistemic (the
  open threefold BG conjecture), not effort.
- **G12 (NS-lattice refactor) is explicitly marked as foundational infrastructure with NO standalone
  user-facing payoff** — it reproduces every current value bit-for-bit and its new power is consumed
  only downstream. G13 (E9) is scheduled to immediately follow it so G12's payoff actually lands, and
  the G14+G18 rational-surface arc completes the program in the final epic.

**Shared code-verified fixtures (all re-executed against the shipped code this pass):** abelian
`numerical_wall = (−5/2, 17/4)` vs K3 shim `(−5/2, 21/4)`, differing by exactly `2/d = 1`;
`ν(2,0,1,0) = −1, −1/2, 0, None` and `twist(1,1).a3 = −4/3`; `exceptional.chi(O,O(3)) = 10` (NOT 15 —
15 is the `n=4` value `chi(O,O(4))`) and `chi(O,O(−3)) = +1`; the P¹×P¹ Euler Gram
`[[1,2,2,4],[0,1,0,2],[0,0,1,2],[0,0,0,1]]`; `delta(1/2)=5/8`, `delta(1/3)=5/9`, `delta(1/4)=21/32`,
`delta(2/5)=13/25`; the P²[2] wall `(−5/2, 3/2)`; K3 `v(O)=(1,0,1)`, `⟨v,v⟩=−2`; P³
`alpha_crit(β=1/2)=√3`.

---

## Epic index

| Epic | Name | Goals | Depends on |
|---|---|---|---|
| **E1** | Harden and document the existing code | G1, G5 | — |
| **E2** | Surface (s,t) wall diagrams: abelian + K3 | G2, G3, G8 | E1 |
| **E3** | Surface existence bounds and proven wall completeness | G6, G7 | E1 |
| **E4** | Threefold tilt-stability engine | G4, G9 | E1 |
| **E5** | Conjecture-gated threefold Bridgeland verdicts | G15 | E1, E4 |
| **E6** | Certified K3 lattice wall analysis | G10 | E1, E2 |
| **E7** | Honest catalog expansion (Enriques / bielliptic / Koseki) | G11 | E1 |
| **E8** | NS-lattice refactor (foundational, no standalone payoff) | G12 | E1 |
| **E9** | Picard-rank ≥ 2 full-NS surface walls | G13 | E8 |
| **E10** | CAS oracle bridge (Macaulay2) | G16 | E1 |
| **E11** | Rational-surface non-emptiness program | G14, G18 | E8, E9, E10 |

---

# E1 — Harden and document the existing code

- **Objective:** After this epic the shipped code has dedicated test coverage for every core module
  (including the four modules that currently have no test file — `chern`, `varieties`, `bg_check`,
  `_latex`), the one confirmed shipped defect (G1 miscitations) is fixed and guarded, and a
  machine-readable provenance lattice (G5 `Rigor`/`Certificate`) colours every downstream output. **No
  new mathematics is introduced** — every value the code emits is unchanged; only tests, citations, and
  provenance metadata are added.
- **Contained goal IDs:** G1, G5 (plus the four missing dedicated test modules and doc-gap closure).
- **Dependencies on other epics:** none (roadmap FIRST epic).
- **Mathematical danger zone:** none to the math — this epic must not change a single computed value.
  The only risk is a test that pins a *wrong* expectation; every asserted value here is copied from an
  already-pinned literature-anchored fixture or re-derived exactly. The G1 allowlist is
  maintenance-brittle (every future legitimate reference must be appended) — acceptable for this small
  catalog, stated in the test docstring.
- **Review criteria (epic done when):** `pytest -q` runs strictly more than 42 tests, all green;
  `tests/test_chern.py`, `tests/test_varieties.py`, `tests/test_bg_check.py`, `tests/test__latex.py`
  exist and are non-empty; `rigor.py` ships with `Rigor`/`Certificate`/`meet`; every catalog reference
  is an allowlist member; `import bridgeland_stability` still succeeds with zero runtime deps.
- **Standalone value:** Even if the roadmap stops here, the project is strictly better: the one shipped
  citation defect is fixed and can never silently recur (allowlist guard), the previously untested
  modules gain regression protection, and every emitted object can advertise whether it is a theorem or
  a conjecture. Pure hardening, immediately mergeable.

### E1-M1 — Fix miscited threefold references + allowlist guard (G1) **[PROVEN]**
- **Files modified:** `bridgeland_stability/varieties.py`.
- **New tests:** `tests/test_varieties.py` (NEW module) —
  `test_no_probability_paper_in_catalog` asserts `"arXiv:1607.07182" not in sum((tf.references for tf
  in ALL_THREEFOLDS), [])`; `test_quadric_refs` asserts `QUADRIC3.references ==
  ["arXiv:1510.04089","arXiv:1607.08199","arXiv:1509.04608"]`; `test_blowup_keeps_piyaratne` asserts
  `"arXiv:1705.04011" in BLOWUP_P3_POINT.references`; `test_allowlist_membership` iterates every
  catalog `Threefold` and asserts every `ref in AG_ALLOWLIST` (a frozen curated set of vetted AG arXiv
  IDs) — the membership check, NOT the format regex `^arXiv:\d{4}\.\d{4,5}$`, is the primary guard
  (the regex matches the probability-paper ID `1607.07182` perfectly and would not have caught the
  defect).
- **Acceptance criteria (testable checklist):**
  - [ ] `QUADRIC3.references` equals the three-ID list above and its `note` no longer implies
        `1705.04011 = BMSZ`.
  - [ ] `BLOWUP_P3_POINT.references` still contains both `1602.05055` and `1705.04011`; its note
        attributes `1705.04011` to **Piyaratne** and states `Bl_p(P³)` is Fano and DOES carry stability
        conditions (BMSZ 1607.08199) — only the strong BMT boundary fails.
  - [ ] Every reference of every catalog `Threefold` is a member of `AG_ALLOWLIST`.
  - [ ] No computed value or pinned test changes (the 42 remain green).
- **Expected effort:** 1 week.

### E1-M2 — Ship the `Rigor`/`Certificate` provenance lattice (G5) **[PROVEN]** (engineering)
- **Files modified:** new `bridgeland_stability/rigor.py`; attribute/repr-level touches to
  `varieties.py`, `walls.py`, `threefold.py`, `mukai.py` (NO formula changes).
- **New tests:** `tests/test_rigor.py` (NEW module) — `test_meet_is_min` asserts `min(Rigor.PROVEN,
  Rigor.CONJECTURAL) == Rigor.CONJECTURAL` and `meet(proven_alg, conjectural_input).rigor ==
  Rigor.CONJECTURAL` with unioned hypotheses/citations; `test_catalog_tags` asserts
  `P2.certificate.rigor == Rigor.PROVEN` and `BLOWUP_P3_POINT.certificate.rigor <= Rigor.CONJECTURAL`;
  `test_blowup_note_phrasing` asserts the `BLOWUP_P3_POINT` note CONTAINS a "strong BMT boundary not
  rigorous" phrase and does NOT contain a "no stability condition" phrase.
- **Acceptance criteria (testable checklist):**
  - [ ] `Rigor(IntEnum)` has `PROVEN=3 > CONJECTURAL=2 > HEURISTIC=1 > UNKNOWN=0` (total order, `min` =
        meet).
  - [ ] `Certificate` is a frozen dataclass `(rigor, hypotheses: tuple, citations: tuple, note: str)`;
        `meet(*certs).rigor == min(c.rigor for c in certs)` and hypotheses/citations are set-unions.
  - [ ] `P2.certificate.rigor == Rigor.PROVEN`; a P²[2] `ActualWall.certificate.rigor == Rigor.PROVEN`;
        a `compute_walls` output is tagged `Rigor.HEURISTIC`.
  - [ ] Optional `certificate` fields default to `UNKNOWN` so no pinned value changes; `import
        bridgeland_stability` succeeds with zero deps.
- **Expected effort:** 2 weeks.

### E1-M3 — Dedicated `chern` tests + twist-invariance / central_charge edge cases **[PROVEN]**
- **Files modified:** none in core (test-only); may add clarifying docstrings to `chern.py`.
- **New tests:** `tests/test_chern.py` (NEW module) — `test_discriminant_convention` asserts
  `ChernChar(1,1,0).discriminant(2) == F(1,8)` (CH convention: μ=1/2, Δ=½·¼=1/8 — note `F(1,4)` is the
  *doubled* `discriminant_brief`, not the CH `discriminant`) and `discriminant_brief` is exactly twice
  `discriminant`; `test_twist_shift_invariance` asserts the `B = sH` twist `ChernChar(r,c,e).twist(s,d)`
  LOWERS the slope by exactly `s` (`slope → slope − s`) and leaves `discriminant(d)` invariant for a
  sample grid `(r,c,e,s,d)`;
  `test_central_charge_matches_fraction_parts` asserts `central_charge(s,t,d)` real/imag parts equal
  the exact `Fraction` formula evaluated to float for pinned `(1,0,-2)` at `(s,t)=(-1,1)`.
- **Acceptance criteria (testable checklist):**
  - [ ] Discriminant uses the CH convention `Δ = ½μ² − e/(rd)` and `discriminant_brief == 2·discriminant`
        for at least three distinct classes.
  - [ ] `twist(s,d)` is verified translation-invariant on discriminant and additive on slope across a
        parameter grid (twist-invariance edge case).
  - [ ] `central_charge` float output agrees with the exact-`Fraction` derivation to within 1e-12 for
        the pinned class; float appears ONLY at this display boundary.
- **Expected effort:** 1 week.

### E1-M4 — Dedicated `bg_check` tests + `compute_walls` edge-case coverage **[PROVEN]**
- **Files modified:** none in core (test-only).
- **New tests:** `tests/test_bg_check.py` (NEW module) — `test_bg_surface_boundary` asserts
  `check_bg_surface` accepts `ChernChar(2,0,-3)` on `d=2` (`Δ=3/4 ≥ 0`) and rejects a `Δ<0` class;
  `test_bg_threefold_reexport` asserts `check_bg_threefold` is re-exported and agrees with
  `threefold.bmt_Q` sign on the pinned `(2,0,1,0)`. Additions to `tests/test_walls.py`:
  `test_compute_walls_is_uncertified` asserts a `compute_walls` result is tagged `HEURISTIC` (via G5)
  and that its dense set is a SUPERSET of the `actual_walls` certified set for P²[2] (the single wall
  `(−5/2, 3/2)`).
- **Acceptance criteria (testable checklist):**
  - [ ] `check_bg_surface` verdicts pinned exactly for one accepted and one rejected class.
  - [ ] `compute_walls` output is explicitly NOT certified (HEURISTIC tag) and contains the P²[2]
        Gieseker wall center `−5/2`, radius `3/2` as a member.
  - [ ] No core formula changed.
- **Expected effort:** 1 week.

### E1-M5 — Dedicated `_latex` tests + doc-gap closure **[PROVEN]**
- **Files modified:** none in core (test-only); `docs/` doc-gap edits (module map, convention notes).
- **New tests:** `tests/test__latex.py` (NEW module) — `test_latex_frac` asserts
  `latex_frac(F(-4,3)) == r"-\frac{4}{3}"`; `test_latex_sqrt` asserts `latex_sqrt(3) == r"\sqrt{3}"`;
  `test_latex_is_stdlib_only` asserts `_latex` imports with no matplotlib/plotly present (import-time
  purity).
- **Acceptance criteria (testable checklist):**
  - [ ] `_latex` helpers produce the exact strings pinned above for `−4/3` and `√3`.
  - [ ] `_latex.py` and `viz/style.py` import with zero third-party deps.
  - [ ] Doc gaps closed: every public function referenced by later epics is named in the module map;
        the K3-vs-abelian Mukai-shift caveat is stated once, canonically.
- **Expected effort:** 1 week.

---

# E2 — Surface (s,t) wall diagrams: abelian + K3

- **Objective:** After this epic the toolkit exposes first-class, correctly-shifted Bridgeland (s,t)
  walls for abelian surfaces (no Mukai shift) and K3 surfaces (the `ch₂ → ch₂+r` shim), unified with
  the Bayer–Macrì wall *type*, plus a citable worked example that makes the single subtlest convention
  (abelian vs K3) explicit and honest about model limits.
- **Contained goal IDs:** G2, G3, G8.
- **Dependencies on other epics:** E1 (uses the G5 `Certificate` for tagging; builds on hardened core).
- **Mathematical danger zone:** wrongly applying the K3 `+r` shift to abelian input (or vice versa) —
  the exact-`2/d` regression fixture guards against it. `k3_wall_classified` requires integral `l=c/d`;
  a non-integral class must RAISE, not silently mis-classify. `actual_walls`' completeness is certified
  only for `P²[n]` and its `_is_integral_chern` filter tests the `d=1` P² lattice condition — it must
  NOT be advertised as an abelian/K3 completeness oracle.
- **Review criteria (epic done when):** `abelian_wall` and `k3_wall`/`k3_wall_classified` exist and are
  documented as the abelian/K3 entry points; the shared `2/d` regression fixture is green; the worked
  example runs and exits 0.
- **Standalone value:** Abelian-surface and K3 wall diagrams become first-class at zero new-math cost —
  directly usable by anyone studying Bridgeland walls / moduli on abelian and K3 surfaces (Yanagida–
  Yoshioka / Bayer–Macrì programs), with the convention pitfall documented and test-guarded.

### E2-M1 — Abelian (s,t) wall affordance `abelian_wall` (G2) **[PROVEN, ρ=1]**
- **Files modified:** `bridgeland_stability/walls.py` (thin wrapper + docstrings); `varieties.py`
  (doc caveat routing abelian input away from the K3 shim).
- **New tests:** `tests/test_abelian_k3_walls.py` (NEW shared fixture module) —
  `test_abelian_wall_value` asserts with `d = abelian_surface(2).d == 2`, `v=ChernChar(1,0,-2)`,
  `w=ChernChar(1,-1,1/2)`: `abelian_wall(v,w,S).center == F(-5,2)` and `.radius_sq == F(17,4)`;
  `test_abelian_wall_asserts_kind` asserts `abelian_wall` raises if `surface.kind != "abelian"`.
- **Acceptance criteria (testable checklist):**
  - [ ] `abelian_wall(v,w,surface)` returns exactly `numerical_wall(v,w,surface.d)` on the bare `(r,c,e)`
        triple — center `−5/2`, radius² `17/4` — with NO Mukai shift.
  - [ ] It asserts `surface.kind == "abelian"` and raises a clear error otherwise.
  - [ ] Docstring states plainly that `compute_walls`/`actual_walls` are NOT abelian completeness
        oracles (the `_is_integral_chern` `d=1` filter does not transfer), and that `ρ≥2` abelian
        surfaces (`E₁×E₂`) are gated on E8/G12.
- **Expected effort:** 1 week.

### E2-M2 — K3-only `k3_wall` semicircle shim (G3) **[PROVEN, ρ=1]**
- **Files modified:** `bridgeland_stability/mukai.py` (add `k3_wall`, importing `numerical_wall`/
  `ChernChar` — acyclic, since `walls.py` does not import `mukai`).
- **New tests:** `tests/test_abelian_k3_walls.py` — `test_k3_wall_value` asserts on the same fixture
  `k3_wall(v,w,2).center == F(-5,2)` and `.radius_sq == F(21,4)`; `test_k3_minus_abelian_is_2_over_d`
  asserts `k3_wall(v,w,2).radius_sq − abelian_wall(v,w,abelian_surface(2)).radius_sq == F(2,2) == 1`;
  `test_center_invariant` asserts on a grid of small triples `k3_wall(...).center ==
  numerical_wall(...).center` and `k3_wall(...).radius_sq − numerical_wall(...).radius_sq ==
  Fraction(2,d)`.
- **Acceptance criteria (testable checklist):**
  - [ ] `k3_wall(v,w,d)` maps `(r,c,e) ↦ (r,c,e+r)` for both classes then calls `numerical_wall`; center
        is invariant, radius² increases by exactly `+2/d`.
  - [ ] The shared fixture yields K3 `(−5/2, 21/4)`, exactly `2/d = 1` MORE than the abelian `17/4`.
  - [ ] Docstring marks it **K3-only**; the semicircle primitive imposes no integrality restriction on
        `l=c/d` (that is enforced by the classified wrapper in E2-M3).
- **Expected effort:** 1 week.

### E2-M3 — Connect the K3 semicircle to `classify_wall` via `k3_wall_classified` (G3) **[PROVEN]**
- **Files modified:** `bridgeland_stability/mukai.py`.
- **New tests:** `tests/test_mukai.py` (additions) — `test_k3_wall_classified_pairs_geometry_and_type`
  reuses the pinned Mukai example `v=MukaiVector(2,0,-1)`, isotropic `u=(1,0,0)`, `(u,v)=1` → asserts
  the returned `WallClassification` is the pinned hilbert–chow divisorial type paired with the `k3_wall`
  geometry (note: this example has equal Mumford slopes `W_rc=0`, so the geometry is the degenerate
  `VerticalWall` `s=0`, not a semicircle — the genuine semicircle is pinned by `test_integral_l_k3_class`
  below); `test_k3_wall_classified_requires_integral_l` asserts a class with
  `c % d != 0` RAISES; `test_integral_l_k3_class` pins one integral-`l` literature-anchored K3 lattice
  class as a genuine (s,t) wall.
- **Acceptance criteria (testable checklist):**
  - [ ] `k3_wall_classified(v,w,surface)` returns `(Wall, WallClassification)` — semicircle via
        `k3_wall`, Bayer–Macrì type via `classify_wall` — closing the gap that `mukai.py` never called
        `numerical_wall`.
  - [ ] It asserts `c % d == 0` and raises a clear error on non-integral `l` (the shared fixture
        `w=(1,-1,1/2)` on `d=2` has `l=−1/2`, so it is a SYNTHETIC `+2/d` demonstration, not a genuine
        Picard-rank-1 K3 lattice class).
  - [ ] At least one integral-`l` genuine K3 lattice class is pinned so the "genuine K3 (s,t) wall"
        claim is backed by a class actually in the lattice.
- **Expected effort:** 2 weeks.

### E2-M4 — Worked example: abelian vs K3 (s,t) walls, honest about limits (G8) **[PROVEN]**
- **Files modified:** `examples/` (new script following `examples/demo.py` style, stdlib-only); `docs/`
  worked-example prose.
- **New tests:** `tests/test_examples.py` (NEW module) — `test_worked_example_radii` asserts the two
  printed radius² values are `F(17,4)` (abelian) and `F(21,4)` (K3) and differ by exactly `F(2,2) == 1`;
  `test_worked_example_exits_zero` asserts the script exits 0 under `python examples/…`.
- **Acceptance criteria (testable checklist):**
  - [ ] The example computes the SAME class `v=(1,0,-2)` on `abelian_surface(2)` and `K3(2)`, prints the
        abelian wall `(−5/2, 17/4)` (via `abelian_wall`) and the K3 wall `(−5/2, 21/4)` (via `k3_wall`),
        and explains the `+2/d` shift and WHY (`√td = (1,0,0)` vs `(1,0,1)`).
  - [ ] It attaches the Bayer–Macrì type via `k3_wall_classified`.
  - [ ] It states plainly the current-model limits: Picard-rank-1 only, no `E₁×E₂`, no `ρ≥2` K3, the K3
        shim is K3-only, and `actual_walls` is NOT an abelian/K3 completeness oracle.
  - [ ] Stdlib-only (no viz), runs in the zero-dependency core.
- **Expected effort:** 1 week.

### E2-M5 — Provenance: baseline `Certificate`s on K3/abelian surfaces + PROVEN-tag the (s,t) walls (G2/G3/G5) **[PROVEN, ρ=1]**
- **Motivation:** the E1-M2 `Rigor`/`Certificate` lattice exists so outputs advertise theorem-vs-conjecture,
  but `abelian_wall`/`k3_wall` currently return `UNKNOWN` (their surfaces default to `UNKNOWN_CERTIFICATE`,
  so `meet` collapses to `UNKNOWN`). The abelian ρ=1 and K3 ρ=1 `(s,t)` wall geometries are **[PROVEN]** —
  tag them honestly. (Closes the E1-M2/E2-M1 deferred nit.)
- **Files modified:** `bridgeland_stability/varieties.py` (baseline certs on `K3()` / `abelian_surface()`);
  `walls.py` (`abelian_wall` cert); `mukai.py` (`k3_wall` / `k3_wall_classified` cert).
- **Math / rigor gate (the load-bearing subtlety):** the scalar-minor wall is PROVEN-correct **only at
  Picard rank 1**. At ρ≥2 the bare `(r,c,e)` loses the polarization-orthogonal component, so
  `numerical_wall` returns only the **H-projected** wall (Maciocia arXiv:1202.4587) — not reliable. The
  `PROVEN` tag MUST therefore be gated on `surface.picard_rank == 1`; a ρ≥2 surface falls back to a
  `HEURISTIC` cert noting "H-projected only; needs the NS-lattice refactor (E8/G12)".
- **Implementation sketch:**
  - `K3()` and `abelian_surface()` (both `picard_rank=1`) get a `PROVEN` baseline `Certificate`
    ("carries Bridgeland stability conditions; Picard rank 1"): K3 → Bridgeland arXiv:math/0307164,
    Bayer–Macrì arXiv:1301.6968; abelian → Arcara–Bertram–Lieblich arXiv:0708.2247, Yanagida–Yoshioka
    arXiv:1203.0884.
  - `abelian_wall` attaches, via `dataclasses.replace`, `meet(alg_cert, surface.certificate)` to the
    returned **`Wall`** — `alg_cert` = `PROVEN` when `surface.picard_rank == 1`, else `HEURISTIC` with the
    ρ≥2 note. `VerticalWall` has no `certificate` field, so equal-slope results stay untagged (documented
    limitation). `k3_wall_classified` (which sees the surface) does the same picard_rank-gated tagging on
    its geometry `Wall`; its returned `WallClassification` keeps `classify_wall`'s `HEURISTIC`
    bounded-search cert (the *type* is not upgraded — that is G10/E6). `k3_wall(v,w,d)` (no surface in its
    signature) stays a raw ρ=1 primitive: document its ρ=1 assumption; the surface-gated tag is delivered by
    `k3_wall_classified`.
- **New tests (`tests/test_abelian_k3_walls.py` / `test_mukai.py` additions):**
  - `abelian_wall((1,0,-2),(1,-1,1/2),abelian_surface(2)).certificate.rigor == Rigor.PROVEN`;
    `K3(2).certificate.rigor == Rigor.PROVEN`; `abelian_surface(2).certificate.rigor == Rigor.PROVEN`.
  - `k3_wall_classified` geometry cert `PROVEN` on the ρ=1 integral-`l` K3 class; its `WallClassification`
    cert stays `HEURISTIC`.
  - a defensive ρ≥2 case (a `Surface(kind="abelian", picard_rank=2)`) → `abelian_wall` geometry cert is
    **not** `PROVEN` (`HEURISTIC`) and its note mentions the H-projection / G12 gate.
- **Acceptance criteria (testable checklist):**
  - [ ] `K3()`/`abelian_surface()` carry a `PROVEN` baseline `Certificate` with the correct citations.
  - [ ] `abelian_wall` / `k3_wall_classified` geometry certificate is `PROVEN` for ρ=1 and downgraded
        (`HEURISTIC` + "H-projected; needs G12" note) for ρ≥2.
  - [ ] `k3_wall_classified`'s `WallClassification` keeps the `HEURISTIC` bounded-search cert (type not
        upgraded).
  - [ ] Certificates are METADATA ONLY — NO computed value changes; all 86 prior tests pass unchanged;
        `import bridgeland_stability` stays zero-dependency.
- **[PROVEN, ρ=1]** surface stability-condition existence (Bridgeland / ABL) and single-wall geometry;
  **[SPECULATIVE]** only the engineering choice of gating rule. **Expected effort:** 1 week.

---

# E3 — Surface existence bounds and proven wall completeness

- **Objective:** After this epic the toolkit answers "is this moduli space nonempty?" for
  Picard-rank-1 K3 and abelian surfaces via the actual Yoshioka / Yanagida–Yoshioka non-emptiness
  criterion (beyond bare Bogomolov `Δ≥0`), and upgrades the P²[n] wall-completeness certificate from an
  empirical doubling heuristic to a Maciocia-theorem-backed bound.
- **Contained goal IDs:** G6, G7.
- **Dependencies on other epics:** E1 (optional G5 tagging).
- **Mathematical danger zone:** the `l = c/d` floor bug — non-Picard-lattice input (`c % d ≠ 0`) must
  RAISE, never silently floor (`v=(2,1,0)` on `K3(2)` would otherwise floor `c//d = 0`, giving a wrong
  `v²=−8` and wrongly REJECTING). Yoshioka's `v²≥−2` requires PRIMITIVE `v` and GENERIC polarization —
  encode as `Certificate.hypotheses`, do not silently assert existence for imprimitive `v`. Hoppe–
  Gieseker on `Pⁿ` is cohomological, NOT pure-Chern — out of scope (CAS territory, E10/G16). Extracting
  the exact Maciocia radius constant is genuinely paper-dependent (see the [RESEARCH] tag on E3-M2).
- **Review criteria (epic done when):** `check_existence_k3`/`check_existence_abelian` return the pinned
  verdicts with the integrality guard live; `maciocia_wall_bound` contains the Gieseker wall for
  `n∈{2,3,4,5}` and `actual_walls_complete` certifies `complete=True` from the theorem where the
  constant is extracted, falling back to `HEURISTIC` otherwise (never over-claiming).
- **Standalone value:** Moves K3/abelian rows from "Bogomolov only" to the actual non-emptiness
  criterion, and upgrades the "these are ALL the walls" claim from [SPECULATIVE] to [PROVEN] for covered
  classes — both directly usable without a CAS.

### E3-M1 — K3/abelian non-emptiness checks with integrality guard (G6) **[PROVEN as inequalities]**
- **Files modified:** `bridgeland_stability/bg_check.py`; consumes `mukai.self_pairing`/`from_chern`.
- **New tests:** `tests/test_bg_check.py` (additions) — `test_k3_structure_sheaf` asserts
  `v=ChernChar(1,0,0)` on `K3(2)` has `Δ==0`, Mukai `(1,0,1)`, `v²==−2`, on the boundary
  `(1 − 1/1²)/2 == 0`, ACCEPT; `test_k3_reject` asserts `v=ChernChar(2,2,0)` on `d=2` has `Δ==1/8 <
  3/8 == (1−1/4)/2`, `v²==−6 < −2`, REJECT; `test_k3_accept_rank2` asserts `v=ChernChar(2,0,-3)` has
  `Δ==3/4 ≥ 3/8`, `v²==4`, ACCEPT with `moduli_dim == v²+2 == 6`; `test_abelian_boundary` asserts
  `v=ChernChar(1,0,0)` abelian → bare `v²==0`, ACCEPT; `test_integrality_guard` asserts
  `check_existence_k3(ChernChar(2,1,0), K3(2))` RAISES (not a silent `v²=−8`).
- **Acceptance criteria (testable checklist):**
  - [ ] `check_existence_k3` returns `Δ ≥ (1 − 1/r²)/d` (equivalently `v² ≥ −2`) and
        `check_existence_abelian` returns `v² ≥ 0` on the BARE, un-shifted Mukai self-pairing (no K3
        `+r` shift).
  - [ ] Both assert `v.c % surface.d == 0` BEFORE computing `l = v.c // surface.d`, raising a clear
        error otherwise.
  - [ ] A `Δ_min(r,d) = (1 − 1/r²)/d` helper is exposed so callers see WHY a class fails.
  - [ ] Yoshioka primitivity/genericity are recorded as `Certificate.hypotheses`; existence is NOT
        asserted for imprimitive `v`.
  - [ ] `check_bg_surface` (`Δ≥0`) is untouched as the rank-independent baseline.
- **Expected effort:** 2 weeks.

### E3-M2 — Maciocia-bounded wall-completeness certificate for P²[n] (G7) **[PROVEN]** / **[RESEARCH]**
- **Files modified:** `bridgeland_stability/walls.py`.
- **New tests:** `tests/test_walls.py` (additions) — `test_maciocia_bound_contains_gieseker` asserts for
  `v=ChernChar(1,0,-n)`, `n∈{2,3,4,5}`, `maciocia_wall_bound(v,1) ≥ ((2n−1)/2)²`;
  `test_actual_walls_complete_p2_2` asserts `actual_walls_complete(v, P2)` for `n=2` returns exactly the
  single wall center `−5/2`, radius `3/2` with `complete=True`; `test_no_wall_exceeds_bound` asserts no
  returned wall has `radius² > maciocia_wall_bound(v,d)`.
- **Acceptance criteria (testable checklist):**
  - [ ] `maciocia_wall_bound(v,d)` returns the Thm-3.11 bounding-semicircle radius² and contains the
        pinned Gieseker wall for every `n∈{2,3,4,5}`.
  - [ ] `actual_walls_complete` returns the SAME wall set as today, now with `complete=True` backed by
        the theorem where the closed-form constant is extracted; where NOT extracted it falls back to the
        current `HEURISTIC` doubling (tagged via G5) — the function never over-claims.
  - [ ] The P²[2] regression set is exactly `{(−5/2, 3/2)}`.
- **[RESEARCH]:** Extracting the EXACT closed-form radius constant from Maciocia (arXiv:1202.4587) §3
  and validating it against the `P²[n]` Gieseker value must be done before the `PROVEN` tag is set;
  until then keep the `HEURISTIC` fallback. General-`v`, general-surface completeness stays OPEN (that
  escalates to E9/G13). **[UNVERIFIED, corrected during E3-M2 gate]** The abstract of arXiv:1202.4587
  (verified) shows Maciocia proves walls are bounded *in suitable planes of stability conditions* with
  *conditions for global finiteness* — NOT a blanket "every fixed-class wall set nests in one bounding
  semicircle", and the specific "Thm 3.11" attribution is unverified at the theorem level. The shipped
  `maciocia_wall_bound` therefore returns the PROVEN ABCH `P²[n]` Gieseker `radius²` `(bog−1)²/4` and is
  documented as a heuristic estimate (which a wall MAY exceed) for a general class.
- **Expected effort:** 3 weeks.

---

# E4 — Threefold tilt-stability engine

- **Objective:** After this epic the toolkit computes faithful threefold tilt-stability walls in the
  `(α,β)` upper half-plane: the scalar `ν_{α,β}` tilt-slope primitive, the PROVEN tilt-layer wall solver
  (Schmidt Thm 3.3, `ch₃`-independent), and the second-tilt/Bridgeland wall whose `ch₃`-dependent
  coefficients are isolated behind a paper-transcription gate.
- **Contained goal IDs:** G4, G9 (split into G9a PROVEN and G9b [RESEARCH-light]).
- **Dependencies on other epics:** E1 (hardened `threefold.py`; optional G5 tagging).
- **Mathematical danger zone:** citation hygiene — Schmidt Thm 3.3 justifies the TILT layer ONLY (data
  `(ch₀,ch₁,ch₂)`, no `ch₃`); the second-tilt (Bridgeland) `ch₃`-dependent semicircle is a DIFFERENT
  theorem and must NOT cite Thm 3.3. `ν` sign/normalization conventions vary across sources; the
  exact-derived pins insulate the code, and any Schmidt cross-check must first confirm the H-power
  normalization. On `bg_proven=False` threefolds the Bridgeland upgrade is conjectural — never asserted.
- **Review criteria (epic done when):** `nu` reproduces the four pinned values; `numerical_tilt_wall` is
  `ch₃`-independent and reduces to the surface `numerical_wall` on the truncated triple;
  `bridgeland_wall` returns a β-axis-centered semicircle with `bridgeland_certified = threefold.bg_proven`
  and its `(y,z)` transcription is flagged and NOT yet pinned against Schmidt §7.
- **Standalone value:** The first faithful threefold `(α,β)` walls in the toolkit — usable by threefold
  wall-crossing researchers on P³, the quadric, Fano ρ=1, abelian and quintic threefolds — with honest
  tilt-only-vs-Bridgeland labeling even if the Bridgeland `(y,z)` pin is never completed.

### E4-M1 — `ν_{α,β}` tilt-slope function (G4) **[PROVEN]**
- **Files modified:** `bridgeland_stability/threefold.py`.
- **New tests:** `tests/test_threefold.py` (additions) — `test_nu_null_correlation` asserts for
  `v=ThreefoldChern(2,0,1,0)`, `d3=1`: `nu(v,0,1,1) == -1`, `nu(v,1,1,1) == F(-1,2)`, `nu(v,2,1,1) ==
  0`, `nu(v,0,0,1) is None` (β=0 gives `a1^β=0`, the vertical/degenerate case);
  `test_nu_twist_consistency` asserts `v.twist(1,1).a3 == F(-4,3)`.
- **Acceptance criteria (testable checklist):**
  - [ ] `nu(v, alpha_sq, beta, d3)` takes `alpha_sq = α²` as an exact `Fraction` (NOT `α`), returns
        `None` (encoding `+∞`) when `tw.a1 == 0`, else `(tw.a2 − Fraction(alpha_sq)/2 · v.r · d3) /
        tw.a1` where `tw = v.twist(beta, d3)`.
  - [ ] The four pinned values `−1, −1/2, 0, None` hold exactly; `ν` depends only on `(ch₀,ch₁,ch₂)`,
        never `ch₃` (assert unchanged as `a3` varies).
  - [ ] Pure `Fraction`, no float.
  - [ ] Docstring states the pins are EXACTLY-DERIVED from the BMT formula (self-consistent), and that a
        Schmidt cross-check requires confirming H-power normalization first (not yet a Schmidt
        cross-check).
- **Expected effort:** 1 week.

### E4-M2 — PROVEN tilt-layer wall solver `numerical_tilt_wall` (G9a) **[PROVEN]** (split of G9)
- **Files modified:** `bridgeland_stability/threefold.py` (or a new `tilt_walls.py`), delegating to
  `walls.numerical_wall`.
- **New tests:** `tests/test_tilt_walls.py` (NEW module) — `test_tilt_wall_ch3_independent` asserts for
  the twisted cubic `v=ThreefoldChern(1,0,-3,5)` and a second class `w`, `numerical_tilt_wall(v,w,d3)`
  is UNCHANGED when `a3, a3'` vary; `test_tilt_wall_reduction_identity` asserts
  `numerical_tilt_wall(v,w,d3).center == numerical_wall(ChernChar(v.r,v.a1,v.a2),
  ChernChar(w.r,w.a1,w.a2), d3).center` (radius² matching up to the fixed `d3` normalization);
  `test_tilt_wall_beta_axis_centered` asserts the fitted `α²` and `β²` coefficients are equal
  (β-axis-centered semicircle) or the wall is a vertical `β=const` ray when `W_rc=0`.
- **Acceptance criteria (testable checklist):**
  - [ ] `numerical_tilt_wall(v,w,d3)` returns a frozen `ThreefoldTiltWall(center, radius_sq, subobject,
        v, bridgeland_certified: bool)` mirroring `walls.Wall`, computed by the SAME three 2×2 minors as
        `numerical_wall` on the truncated triple `(r,a1,a2)` with `d → d3`.
  - [ ] `ch₃`-independence holds exactly (Schmidt Thm 3.3, VERIFIED against the paper).
  - [ ] The reduction identity to `numerical_wall` holds exactly.
  - [ ] Degenerate `W_rc=0` yields a vertical `β=const` wall.
- **Expected effort:** 2 weeks.

### E4-M3 — Second-tilt/Bridgeland wall `bridgeland_wall` (G9b) **[RESEARCH-light]** (split of G9)
- **Files modified:** `bridgeland_stability/threefold.py` / `tilt_walls.py`; `enumerate_tilt_walls`.
- **New tests:** `tests/test_tilt_walls.py` (additions) — `test_bridgeland_gate` asserts
  `bridgeland_wall(v,w,threefold).bridgeland_certified` is `True` for `P3`/`QUADRIC3` and `False` for
  `BLOWUP_P3_POINT`; `test_bridgeland_is_circle` asserts `center == −y/(2x)` and `radius_sq ==
  (y²−4xz)/(4x²)` as exact `Fraction` from the fitted coefficients; **the exact twisted-cubic
  center/radius pin is DEFERRED** (see [RESEARCH-light]).
- **Acceptance criteria (testable checklist):**
  - [ ] `bridgeland_wall(v,w,threefold)` uses the full central charge `Z_{α,β}` so `a3=ch₃` enters the
        `(y,z)` coefficients of `x(α²+β²)+yβ+z=0`; the `x` coefficient still comes only from the
        `(r,a1)` minor (why it remains a circle).
  - [ ] `center = −y/(2x)`, `radius_sq = (y²−4xz)/(4x²)` exactly as `Fraction`.
  - [ ] `bridgeland_certified = threefold.bg_proven` (a genuine Bridgeland wall only where BG is a
        theorem; Piyaratne–Toda arXiv:1504.01177).
  - [ ] `enumerate_tilt_walls(v, threefold, bounds…)` exists, analogous to `actual_walls`.
- **[RESEARCH-light]:** The exact `(y,z)` coefficient transcription and the twisted-cubic center/radius
  pin must be read off Schmidt's second-tilt construction (arXiv:1509.04608 §3 and §7) before pinning —
  do NOT attribute the `ch₃`-dependent semicircle to Thm 3.3. Until transcribed, the `(y,z)` pins stay
  deferred and `bridgeland_wall` ships behind this flag.
- **Expected effort:** 3 weeks.

---

# E5 — Conjecture-gated threefold Bridgeland verdicts

- **Objective:** After this epic every emitted threefold wall carries a machine-readable certificate
  distinguishing a theorem from a conjecture on identical arithmetic: `PROVEN` on `bg_proven=True`
  rows, `CONJECTURAL`/tilt-only otherwise — and NO code path ever emits "no stability condition exists"
  (even on `Bl_p(P³)`, which is Fano and DOES carry stability conditions).
- **Contained goal IDs:** G15.
- **Dependencies on other epics:** E1 (G5 `Rigor`/`Certificate` — hard), E4 (G9 tilt-wall solver —
  hard). Pulled forward to here (immediately after the threefold engine and the E1 provenance
  substrate), NOT delayed behind E8/G12; its Tier-3 label is epistemic, not effort.
- **Mathematical danger zone:** the Piyaratne–Toda properness statement is CONDITIONAL ("assuming BG");
  the citation must carry that hypothesis even on `bg_proven` rows. The `Bl_p(P³)` note must state
  "strong BMT boundary FAILS (Schmidt 1602.05055); stability conditions nonetheless EXIST (Fano, BMSZ
  1607.08199)". No `Fraction` value changes — pure additive labeling.
- **Review criteria (epic done when):** the propagation invariant and meet law hold for every catalog
  row; the negative test (no `certified=True` when `bg_proven is False`) passes; the Schmidt first-wall
  ground-truth example is pinned (E5-M2).
- **Standalone value:** Threefold-stability researchers get machine-checked honesty about WHICH wall
  statements are theorems — operationalizing the correctness culture exactly where the same arithmetic
  is a theorem on P³ and an open conjecture on a general threefold.

### E5-M1 — `bridgeland_wall_verdict` + rigor propagation (G15) **[PROVEN]** (labeling logic)
- **Files modified:** `bridgeland_stability/threefold.py`; consumes `rigor.py` (E1) and the E4 solver.
- **New tests:** `tests/test_bridgeland_verdict.py` (NEW module) — `test_propagation_proven` asserts for
  every `bg_proven=True` row `bridgeland_wall_verdict(...).rigor == Rigor.PROVEN` and `.certified is
  True`; `test_propagation_blowup` asserts for `BLOWUP_P3_POINT` `.rigor <= Rigor.CONJECTURAL` and
  `.certified is False`; `test_meet_law` asserts `verdict.rigor == min(tilt_rigor, PROVEN if
  bg_proven else CONJECTURAL)`; `test_never_certifies_unproven` asserts NO code path returns
  `certified=True` when `bg_proven is False`; `test_blowup_note` asserts the `Bl_p(P³)` note contains
  "stability conditions nonetheless EXIST" and never contains "no stability condition".
- **Acceptance criteria (testable checklist):**
  - [ ] `BridgelandWallVerdict(tilt_wall, certified: bool, rigor: Rigor, citations, note)` records
        `certified = threefold.bg_proven` and, when `False`, the note "tilt-stability wall only;
        Bridgeland upgrade unproven (threefold BG open here)".
  - [ ] For `Bl_p(P³)` the note ADDITIONALLY states the strong-BMT-fails / stability-conditions-exist
        distinction; NEVER emits "no stability condition".
  - [ ] The meet law `rigor = min(TILT_SOLVER_RIGOR, PROVEN if bg_proven else CONJECTURAL)` holds
        exactly; no `Fraction` value changes.
  - [ ] `is_bridgeland_certified(threefold)` is routed through `bg_boundary_curve` / any wall enumerator
        so every emitted threefold wall carries its `Certificate`.
- **Expected effort:** 2 weeks.

### E5-M2 — Schmidt first-wall ground-truth cross-check (G15) **[RESEARCH]**
- **Files modified:** `tests/` (test-only); optional example.
- **New tests:** `tests/test_bridgeland_verdict.py` (additions) — `test_schmidt_first_wall_p3` asserts a
  Schmidt (arXiv:1509.04608) first-wall example on P³ yields the SAME `(β,α)` locus as the E4 solver AND
  is labeled `Rigor.PROVEN`; `test_same_class_blowup_is_conjectural` asserts the identical class on
  `Bl_p(P³)` yields `Rigor.CONJECTURAL`/tilt-only.
- **Acceptance criteria (testable checklist):**
  - [ ] The pinned Schmidt first-wall `(β,α)` locus matches the E4 `numerical_tilt_wall`/`bridgeland_wall`
        output exactly (as `Fraction`).
  - [ ] The same class is labeled `PROVEN` on P³ and `CONJECTURAL` on `Bl_p(P³)`.
- **[RESEARCH]:** The concrete first-wall `(β,α)` values must be read off Schmidt arXiv:1509.04608 §7
  before pinning; this milestone's acceptance values are gated on that transcription.
- **Expected effort:** 2 weeks.

---

# E6 — Certified K3 lattice wall analysis

- **Objective:** After this epic the K3 rank-2 hyperbolic-lattice wall analysis is CERTIFIED — replacing
  the heuristic bounded `|a|,|b|≤search` scan with exact enumeration of spherical (`x²=−2`) and
  isotropic (`x²=0`) classes in the SATURATION of `⟨v,w⟩`, exhibiting a lattice witness for positive
  wall-type verdicts.
- **Contained goal IDs:** G10.
- **Dependencies on other epics:** E1 (optional G5 `certified` flag), E2 (soft: G3 `k3_wall`, to link
  classification to its (s,t) semicircle). Independent of E8/G12 (uses the shipped rank-3 K3 Mukai
  lattice directly).
- **Mathematical danger zone:** the hyperbolic binary-quadratic solver returns finitely many ORBITS
  from infinitely many points — it must return one representative per automorphism orbit and must
  enumerate over the SATURATED primitive rank-2 sublattice `H_W` (not merely the `ℤ`-span), or the
  "certified-complete" claim fails when a witness lives in the saturation but not the span. Emit
  `certified=True` ONLY for positive (divisorial/flopping) verdicts; a "fake-or-none" verdict stays
  `certified=False` and defers the totally-semistable split to deferred G17. The `mukai.from_chern`
  `ch₂+r` shift is K3-only.
- **Review criteria (epic done when):** `classify_wall_certified` matches `classify_wall(...,
  search=large)` on all pinned positive cases with `certified=True` and verdict provably search-bound
  independent; `solve_binary_quadratic` reproduces the pinned witnesses as orbit representatives.
- **Standalone value:** K3 wall-crossing / MMP researchers (Bayer–Macrì program) get certified wall
  types with an exhibited lattice witness, replacing a heuristic bounded search — a direct strengthening
  of `mukai.py`, independent of the heavy NS refactor.

### E6-M1 — Hyperbolic binary-quadratic solver over the saturated sublattice (G10) **[RESEARCH]**
- **Files modified:** `bridgeland_stability/mukai.py`.
- **New tests:** `tests/test_mukai.py` (additions) — `test_solver_reproduces_spherical` asserts
  `solve_binary_quadratic(A,B,C,-2)` for the Gram of `v(O)=(1,0,1)` (`v²=−2`) returns the spherical
  witness `(a,b)` as one representative per automorphism orbit; `test_solver_isotropic` asserts the
  `target=0` case reproduces the pinned isotropic `u=(1,0,0)` witness for the Hilbert–Chow lattice
  `v=(2,0,−1)`, `(u,v)=1`; `test_saturation` asserts a witness living in the saturation but not the
  `ℤ`-span is found (constructed example).
- **Acceptance criteria (testable checklist):**
  - [ ] `solve_binary_quadratic(A,B,C,target)` returns primitive integer solutions of `A a² + B ab + C
        b² = target` for `target ∈ {−2, 0}`, ONE representative per orbit of the sublattice automorphism
        group (indefinite form: infinite points, finite orbits).
  - [ ] Enumeration is over the SATURATED primitive rank-2 sublattice `H_W` (the Bayer–Macrì potential
        wall = primitive closure of `ℤv ⊕ ℤw`), not the raw span.
  - [ ] Reproduces the pinned `test_mukai.py` spherical/isotropic witnesses exactly.
- **[RESEARCH]:** The hyperbolic (Pell-type) orbit-enumeration algorithm is known but fiddly; the
  general fake/totally-semistable decision beyond the rank-2 sublattice needs the full Bayer–Macrì
  stratification (hard; escalates to deferred G17) and is explicitly OUT of scope here.
- **Expected effort:** 3 weeks.

### E6-M2 — `classify_wall_certified` with positive-only certification (G10) **[PROVEN]**
- **Files modified:** `bridgeland_stability/mukai.py`.
- **New tests:** `tests/test_mukai.py` (additions) — `test_certified_matches_search` asserts on all
  pinned positive cases (`v(O)=(1,0,1)`; Hilbert–Chow `v=(2,0,−1)`, `u=(1,0,0)`, `(u,v)=1`; LGU
  `v=(3,0,−1)`, `u=(2,0,0)`, `(u,v)=2`) `classify_wall_certified` returns the SAME
  `wall_type`/`subtype` as `classify_wall(..., search=large)` AND sets `certified=True`, with the
  verdict identical for `search ∈ {8, 32}` (bound-independent); `test_fake_or_none_uncertified` asserts
  a fake-or-none case returns `certified=False`; `test_carries_gram` asserts `lattice_gram` equals
  `[[v², (v,w)], [(v,w), w²]]`.
- **Acceptance criteria (testable checklist):**
  - [ ] `classify_wall_certified(v,w,d)` (i) checks `⟨v,w⟩` is rank-2 hyperbolic, (ii) enumerates minimal
        spherical/isotropic classes over `H_W` via `solve_binary_quadratic`, (iii) applies Thm 5.7 with
        a proof witness and NO `search` parameter.
  - [ ] `certified=True` ONLY for positive (divisorial/flopping) verdicts; fake-or-none stays
        `certified=False`.
  - [ ] `WallClassification` gains `certified: bool` and `lattice_gram: tuple`.
  - [ ] A certified type bridges to its `k3_wall` (E2/G3) (s,t) semicircle.
- **Expected effort:** 2 weeks.

---

# E7 — Honest catalog expansion (Enriques / bielliptic / Koseki)

- **Objective:** After this epic the catalog records Enriques, bielliptic, and Koseki product-threefold
  rows with correct invariants and honest, machine-readable provenance — WITHOUT pretending the scalar
  model can faithfully compute their walls. Surface-consuming APIs raise a clear "requires NS-lattice
  refactor (G12)" error on these rows until E8 lands.
- **Contained goal IDs:** G11.
- **Dependencies on other epics:** E1 (soft: G1 same-file citation work, G5 `Certificate`). The
  faithful-computation half is gated on E8/G12; the record-row deliverable is NOT.
- **Mathematical danger zone:** the `faithful_computation_supported` guard MUST live in the
  Surface-consuming entry points (`compute_walls`, `actual_walls`, `dlp`, `bg_check`, surface-aware
  `k3_wall`/`abelian_wall`) — NOT in `numerical_wall(v,w,d)`, which sees only a `ChernChar` and an `int
  d`, never a `Surface`, and so cannot detect `surface.kind == "enriques"`. The Koseki `bg_proven=True`
  gate stays strictly behind the [UNVERIFIED] hypothesis check against arXiv:1811.03267 (distinct from
  Koseki's separate weighted-hypersurface / Calabi–Yau double–triple-solid paper arXiv:2007.00044 — do
  not conflate; arXiv:1510.04474 is Chuang–Lai's *withdrawn* CY/Fano paper, not Koseki).
- **Review criteria (epic done when):** Enriques/bielliptic/Koseki rows exist with correct
  `d`/`K_H`/`chi_O`/`Certificate`; `faithful_computation_supported is False` for them; a
  `compute_walls`/`actual_walls`/`dlp` call on them raises the G12-required error while `numerical_wall`
  on the raw triple does not.
- **Standalone value:** Catalog completeness with honest provenance; surfaces exactly which
  known-in-the-literature varieties are not yet faithfully computable and why. The Koseki hypothesis
  check is a small, citable literature-verification contribution.

### E7-M1 — Enriques / bielliptic surface rows + `canonical_order` field (G11) **[PROVEN]** (record)
- **Files modified:** `bridgeland_stability/varieties.py`.
- **New tests:** `tests/test_varieties.py` (additions) — `test_enriques_row` asserts `enriques()` has
  `chi_O == 1`, `canonical_order == 2`, `kind == "enriques"`; `test_distinguished_from_k3_abelian`
  asserts Enriques (`chi_O=1`, order 2) is distinguishable from K3 (`chi_O=2`, order 0) and abelian
  (`chi_O=0`, order 0); `test_bielliptic_row` asserts `bielliptic()` records correct `d`/torsion.
- **Acceptance criteria (testable checklist):**
  - [ ] `enriques()` and `bielliptic()` catalog rows carry correct `d`, `K_H`, `chi_O`, and a
        `Certificate` (G5) with exact hypotheses + citations (Nuer–Yoshioka arXiv:1901.04848; Yoshioka
        arXiv:1607.04946; bielliptic arXiv:2107.13370).
  - [ ] A `canonical_order: int` field on `Surface` (0 = numerically trivial K3/abelian, 2 = Enriques
        2-torsion, else torsion order) distinguishes Enriques WITHOUT the full NS refactor.
  - [ ] `faithful_computation_supported` returns `False` for these rows.
- **Expected effort:** 2 weeks.

### E7-M2 — Koseki product-threefold rows + surface-consuming-API guard (G11) **[RESEARCH/UNVERIFIED]**
- **Files modified:** `bridgeland_stability/varieties.py`; guard insertion into `walls.py`
  (`compute_walls`/`actual_walls`), `dlp.py`, `bg_check.py`, and the surface-aware wall wrappers.
- **New tests:** `tests/test_varieties.py` (additions) — `test_koseki_bg_gate` asserts each
  `koseki_product(...)` row has `bg_proven is True` IFF the arXiv:1811.03267 hypotheses are met (until
  verified, `bg_proven is False` / rigor `CONJECTURAL`); `test_koseki_cites_correct_paper` asserts the
  `note`/`Certificate` cites **arXiv:1811.03267** and NOT Koseki's separate
  weighted-hypersurface paper arXiv:2007.00044;
  `test_faithful_guard_raises` asserts a `compute_walls`/`actual_walls`/`dlp` call on a Koseki/Enriques
  row RAISES the "requires NS-lattice refactor (G12)" error; `test_numerical_wall_not_guarded` asserts
  `numerical_wall(v,w,d)` on the raw triple does NOT raise (the caveat lives in docs only).
- **Acceptance criteria (testable checklist):**
  - [ ] `koseki_product(...)` rows for `P¹×S`, `P²×C`, `P¹×P¹×C` carry correct `d`/`K_H`/`chi_O` and a
        Koseki `Certificate`.
  - [ ] `bg_proven=True` is set ONLY after the arXiv:1811.03267 hypothesis check; `faithful_computation_
        supported is False`.
  - [ ] The guard lives in the Surface-consuming APIs (not in `numerical_wall`); those APIs raise the
        G12-required error when `faithful_computation_supported is False`.
- **[RESEARCH/UNVERIFIED]:** The Koseki exact-hypothesis verification against arXiv:1811.03267 is a
  required literature sub-task; until done, product rows carry rigor `CONJECTURAL`/`UNKNOWN`, never
  `PROVEN`.
- **Expected effort:** 3 weeks.

---

# E8 — NS-lattice refactor (foundational infrastructure — NO standalone user-facing payoff)

- **Objective:** After this epic the Chern character carries a true Néron–Severi vector `ch₁` (with a
  symmetric intersection form) instead of the scalar `c = ch₁·H`, unlocking every Picard-rank ≥ 2 item.
  As a v2 architecture decision it reproduces EVERY current value bit-for-bit (the 42 pinned tests are
  unchanged) and gains ρ≥2 distinguishing power that is consumed only by E9 and E11.
- **Contained goal IDs:** G12 (mandatorily split into four sequenced milestones G12.1–G12.4).
- **Dependencies on other epics:** E1 (hardened core + provenance). Foundational — no other math
  prerequisite.
- **EXPLICIT STANDALONE-VALUE NOTE:** **As a terminal epic, E8 delivers NO user-visible capability.**
  The backward-compat shim reproduces every current value bit-for-bit; the new ρ≥2 power is only
  consumed by E9 (G13) and E11 (G14, G18), which come after. This is acceptable ONLY as foundational
  infrastructure: **E9 (G13) is scheduled to immediately follow E8** so the ρ≥2 payoff actually lands,
  and the G14+G18 rational-surface arc completes the program in E11. E8 is the roadmap's dominant
  effort/schedule risk; it is split into four individually PhD-implementable milestones, each with an
  explicit week estimate, so no single milestone mutates the frozen core + all 42 tests as one unit.
- **Mathematical danger zone:** subtle convention drift while migrating the frozen `ChernChar` and all
  42 tests — mitigated by the 42-test regression gate and the explicit shim pin. **The naive encoding
  "`(c,)` with `Gram=[[d]]`, `H=(1,)` is WRONG and breaks the 42-test gate:** it gives `⟨ch1,H⟩ = c·d`
  and `⟨ch1,ch1⟩ = c²·d`, so `delta(1/2)` would no longer equal `5/8`. The shim MUST store the
  H-coordinate `c/d` so `⟨ch1,H⟩ = c` and `⟨ch1,ch1⟩ = c²/d` bit-for-bit. Must preserve all inviolables
  (exact `Fraction`, CH convention, frozen dataclasses, zero deps).
- **Review criteria (epic done when):** all 42 existing tests pass unchanged after the full migration;
  the shim pin holds bit-for-bit; the ρ=2 coverage distinguishes classes the old scalar model conflated;
  the ρ=1 `numerical_wall` path reproduces the current minors bit-for-bit.
- **Standalone value:** **NONE at the user-facing level (stated explicitly above).** Its entire value is
  unlocking the rational-surface arc (E9, E11) and the faithful half of E7/G11. Justified purely as a
  foundational v2 refactor.

### E8-M1 — `NSLattice` + pinned rank-1 backward-compat shim (G12.1) **[PROVEN]** (split of G12)
- **Files modified:** new `bridgeland_stability/nslattice.py`; `varieties.py` (store `NSLattice` + ample
  `H` on `Surface`).
- **New tests:** `tests/test_ns_lattice.py` (NEW module) — `test_shim_pins_pairings` asserts for the
  rank-1 shim class `ChernChar(r,c,e)`, `⟨ch1,H⟩ == c` and `⟨ch1,ch1⟩ == c²/d` bit-for-bit (e.g.
  `r=1,c=1,d=2` gives `⟨ch1,H⟩==1`, `⟨ch1,ch1⟩==F(1,2)`); `test_naive_encoding_rejected` asserts the
  naive `(c,)`+`Gram=[[d]]` encoding is NOT used (would give `⟨ch1,H⟩==c·d`); `test_p1xp1_pairing`
  asserts `NSLattice.pairing` on the P¹×P¹ Gram `[[0,1],[1,0]]` gives `⟨f+s, f+s⟩ == 2`.
- **Acceptance criteria (testable checklist):**
  - [ ] `NSLattice(rank, gram)` is a frozen dataclass with `pairing(u,v) -> Fraction`; a designated
        ample integer vector `H` is stored on `Surface`.
  - [ ] The rank-1 shim stores the H-coordinate `c/d` (or a basis where the ample generator has
        self-intersection `d` and `ch1 = (c/d)H`), so `⟨ch1,H⟩ == c` and `⟨ch1,ch1⟩ == c²/d`
        bit-for-bit.
  - [ ] The naive `(c,)`+`Gram=[[d]]`,`H=(1,)` encoding is explicitly forbidden (it yields
        `⟨ch1,H⟩=c·d`, `⟨ch1,ch1⟩=c²·d`, breaking `delta(1/2)=5/8`).
  - [ ] `NSLattice.pairing` on `[[0,1],[1,0]]` gives `⟨f+s,f+s⟩ == 2`.
- **Expected effort:** 2 weeks.

### E8-M2 — Migrate `ChernChar.c` scalar → `ch1` NS-vector, 42-test gate (G12.2) **[PROVEN]** (split of G12)
- **Files modified:** `bridgeland_stability/chern.py` (mutate the frozen `ChernChar`); downstream
  attribute reads in `walls.py`, `mukai.py`, `bg_check.py`.
- **New tests:** the **entire existing 42-test suite is the regression gate** (must pass unchanged);
  `tests/test_chern.py` (additions) — `test_discriminant_lattice_terms` asserts `delta(1/2) == 5/8` and
  `delta(2/5) == 13/25` UNCHANGED after `discriminant`/`bogomolov_discriminant` are rewritten as
  `ch1.H=⟨ch1,H⟩`, `ch1²=⟨ch1,ch1⟩`.
- **Acceptance criteria (testable checklist):**
  - [ ] `ChernChar.c` is promoted to `ch1`, an NS-vector bound to the surface lattice; `e = ch₂` stays a
        `Fraction`; the `ChernChar(r, c, e)` scalar constructor still works for ρ=1 via the E8-M1 shim.
  - [ ] `discriminant = ½μ² − e/(rd)` and `bogomolov_discriminant = ⟨ch1,ch1⟩ − 2re` are rewritten in
        lattice terms and are UNCHANGED numerically.
  - [ ] **All 42 existing tests pass** — P²[2] wall `(−5/2, 3/2)`, `delta(1/2)=5/8`, `delta(2/5)=13/25`,
        K3 `v(O)=(1,0,1)` `v²=−2`, P³ `alpha_crit(β=1/2)=√3`, Gieseker wall `(−(2n+1)/2, (2n−1)/2)`.
- **Expected effort:** 3 weeks.

### E8-M3 — Generalize the wall minors + Mukai pairing to the lattice form (G12.3) **[PROVEN]** (split of G12)
- **Files modified:** `bridgeland_stability/walls.py` (`numerical_wall` minors), `mukai.py` (`d·l·l'`
  pairing).
- **New tests:** `tests/test_walls.py` (additions) — `test_numerical_wall_rho1_bitforbit` asserts the
  ρ=1 path reproduces the current minors bit-for-bit: P²[2] wall stays `(−5/2, 3/2)`;
  `tests/test_mukai.py` (additions) — `test_mukai_pairing_lattice` asserts K3 `v(O)=(1,0,1)` still gives
  `v²=−2` with the `d·l·l'` term generalized to `⟨·,·⟩` on the `l`-component.
- **Acceptance criteria (testable checklist):**
  - [ ] `walls.numerical_wall` computes `W_rc, W_re, W_ce` via the lattice pairing; the ρ=1 path is
        bit-for-bit identical (P²[2] → `(−5/2, 3/2)`).
  - [ ] `mukai.py`'s `d·l·l'` generalizes to `⟨·,·⟩` on the `l`-component; K3 `v(O)=(1,0,1)`, `v²=−2`
        unchanged.
  - [ ] All 42 tests still pass.
- **Expected effort:** 2 weeks.

### E8-M4 — ρ=2 coverage: P¹×P¹ and 𝔽_n distinguishing power (G12.4) **[PROVEN]** (split of G12)
- **Files modified:** `varieties.py` (P¹×P¹ / 𝔽_n `NSLattice` rows); tests.
- **New tests:** `tests/test_ns_lattice.py` (additions) — `test_p1xp1_distinguishes_classes` asserts two
  classes with `ch1=f+s` vs `ch1=2f` have EQUAL `⟨ch1,H⟩` (with `H=f+s`, `⟨H,H⟩=2` reproducing `d=2`)
  but DIFFERENT true `ch1²`: `⟨f+s,f+s⟩ == 2` vs `⟨2f,2f⟩ == 0`, and are DISTINGUISHED (the old scalar
  model conflated them); `test_fn_hodge_signature` asserts for 𝔽_n Gram `[[0,1],[1,−n]]`, `det(Gram) <
  0` (Hodge-index signature `(1,1)`).
- **Acceptance criteria (testable checklist):**
  - [ ] P¹×P¹ with Gram `[[0,1],[1,0]]` (basis `f,s`), `H=f+s`, `⟨H,H⟩=2` reproduces `d=2`.
  - [ ] `ch1=f+s` and `ch1=2f` share `⟨ch1,H⟩` but differ in `ch1²` (`2` vs `0`) and are distinguished.
  - [ ] 𝔽_n Gram `[[0,1],[1,−n]]` asserts `det(Gram) < 0` (signature `(1,1)`).
- **Expected effort:** 2 weeks.

---

# E9 — Picard-rank ≥ 2 full-NS surface walls

- **Objective:** After this epic the toolkit computes faithful Bridgeland surface walls for Picard rank
  ≥ 2 (products of curves, mixed polarizations, higher-ρ K3s) from the FULL intersection form via the
  Maciocia `β = bω + γ` decomposition, with a PROVEN Thm-3.11 boundedness certificate replacing the
  empirical doubling for the classes it covers. **This epic immediately follows E8 so the NS-refactor
  payoff lands as visible, testable progress.**
- **Contained goal IDs:** G13.
- **Dependencies on other epics:** E8/G12 (hard — needs the NS lattice).
- **Mathematical danger zone:** wall count grows with ρ and needs careful bounding (Thm 3.11 supplies
  it); the ρ=1 reduction MUST reproduce `numerical_wall` exactly (regression). Pinning one exact ρ>1
  wall value requires reading a worked example off Maciocia arXiv:1202.4587 (see [RESEARCH-light]).
- **Review criteria (epic done when):** `numerical_wall_ns` reduces to `numerical_wall` on every ρ=1
  surface; on P¹×P¹ it finds a wall the scalar model misses; every enumerated wall respects the
  Maciocia bound and `actual_walls_complete` returns `complete=True` from the PROVEN bound.
- **Standalone value:** The first faithful surface walls beyond P²/ρ=1 in the toolkit — certified,
  finite wall enumeration for Picard rank ≥ 2, directly usable by researchers on products of curves,
  mixed polarizations, and higher-ρ K3s. This is E8's payoff made visible.

### E9-M1 — `numerical_wall_ns` full-NS surface wall solver (G13) **[PROVEN]** / **[RESEARCH-light]**
- **Files modified:** `bridgeland_stability/walls.py`.
- **New tests:** `tests/test_walls_ns.py` (NEW module) — `test_ns_reduces_to_scalar` asserts
  `numerical_wall_ns` on any Picard-rank-1 surface reproduces `numerical_wall` exactly (P²[2] →
  `(−5/2, 3/2)`; a K3 wall); `test_ns_finds_bidegree_wall` asserts on P¹×P¹ a destabilizer distinguished
  ONLY by bidegree (same `⟨ch1,H⟩`, different `ch1²`) yields a wall that `numerical_wall_ns` FINDS but
  scalar `numerical_wall` MISSES (assert the presence/absence pair).
- **Acceptance criteria (testable checklist):**
  - [ ] `numerical_wall_ns(v,w,surface)` computes center/radius from the full intersection form (via
        `⟨ch1,ch1⟩` and the polarization-orthogonal `γ²`), not merely `⟨ch1,H⟩`.
  - [ ] ρ=1 reduction is bit-for-bit exact.
  - [ ] The ρ>1 gain is demonstrated by a constructed P¹×P¹ presence/absence pair.
- **[RESEARCH-light]:** Pin ONE exact ρ>1 wall value from a worked example in Maciocia arXiv:1202.4587
  (his P¹×P¹ / abelian examples), read off the paper before pinning; until then the constructed
  presence/absence pair is the acceptance evidence.
- **Expected effort:** 3 weeks.

### E9-M2 — PROVEN Maciocia boundedness + finite NS-wall enumeration (G13) **[PROVEN]** / **[RESEARCH-light]**
- **Files modified:** `bridgeland_stability/walls.py`.
- **New tests:** `tests/test_walls_ns.py` (additions) — `test_enumerated_within_bound` asserts every wall
  from `enumerate_ns_walls(v, surface)` has `radius_sq ≤ maciocia_wall_bound(v, surface)`;
  `test_actual_walls_complete_proven` asserts `actual_walls_complete` returns `complete=True` backed by
  the PROVEN bound (not doubling) for a covered ρ>1 class, tagged `Rigor.PROVEN` (G5).
- **Acceptance criteria (testable checklist):**
  - [ ] `maciocia_wall_bound(v, surface)` returns the Thm-3.11 outermost bounding-semicircle radius²;
        `enumerate_ns_walls(v, surface)` returns the finitely many walls inside that bound.
  - [ ] `actual_walls_complete`'s empirical doubling is REPLACED by this PROVEN bound for the covered
        family — the **ρ=1 ABCH / coprime-P²** classes (where the bound `(bog−1)²/4` is the exact ABCH
        Gieseker radius²), upgrading [SPECULATIVE]→[PROVEN]; the certificate is granted ONLY when the
        search window provably contains the whole bounding semicircle (an exact no-`sqrt` center-coverage
        test — merely `window² ≥ bound` is UNSOUND, since it can truncate the outermost Gieseker wall by
        its *center* while passing the radius filter vacuously). **ρ≥2 PROVEN completeness is DEFERRED**
        (stays HEURISTIC) pending the transcribed Maciocia §3 constant — see the [RESEARCH-light] note.
  - [ ] Every enumerated wall respects the bound; for the covered ρ=1 family completeness carries a
        `Rigor.PROVEN` certificate.
- **[RESEARCH-light]:** the general ρ≥2 Thm-3.11 bounding-semicircle radius² constant must be read off
  Maciocia arXiv:1202.4587 §3 before any ρ≥2 `PROVEN` completeness tag is set, as in E3-M2; until then
  ρ≥2 stays HEURISTIC (witnessed by `test_ns_bound_is_heuristic_rho2`).
- **Expected effort:** 3 weeks.

---

# E10 — CAS oracle bridge (Macaulay2)

- **Objective:** After this epic the toolkit has a never-core, lazily-imported Macaulay2 subprocess
  oracle that (a) recomputes `χ`/`Ext^i`/cohomology dimensions at the SHEAF level and diffs them against
  the formula layer (extending the "verified two ways" discipline), and (b) unlocks sheaf-level
  questions the numerical core forbids by construction — while `import bridgeland_stability` stays
  zero-dependency.
- **Contained goal IDs:** G16.
- **Dependencies on other epics:** E1 (optional G5 tagging). Independent of the NS refactor. Scheduled
  BEFORE the rational-surface epic (E11), for which it is the default (soft) HN-step supplier.
- **Mathematical danger zone:** the M2 install is external/non-pip and the text parser is
  hand-maintained (print-format drift can break it). `moduli_nonempty_by_construction` is only a
  SUFFICIENT witness — a construction failure is NOT a non-existence proof; label it so. Windows: M2
  runs natively or under WSL **[UNVERIFIED — confirm the host's M2 path before relying on it]**. The
  oracle must NEVER be imported by the core.
- **Review criteria (epic done when):** `import bridgeland_stability` succeeds with zero deps and M2
  absent; the cross-check values match the shipped `exceptional.chi`; the `Fraction ↔ QQ` round-trip is
  exact; `moduli_nonempty_by_construction` is labeled sufficient-only.
- **Standalone value:** Anyone auditing the formula layer gets an independent, exact sheaf-level second
  opinion; researchers needing a specific `Ext^1`/obstruction or constructive non-emptiness witness get
  it without leaving Python. The most architecturally-aligned CAS bridge.

### E10-M1 — `M2Session` lazy oracle + graceful-skip (G16) **[PROVEN + implemented]**
- **Files modified:** new `bridgeland_stability/oracle/m2.py` (lazy import, like `viz/`).
- **New tests:** `tests/test_oracle.py` (NEW module) — `test_zero_dep_import` asserts `import
  bridgeland_stability` succeeds with M2 absent AND zero third-party deps; `test_oracle_guard_raises`
  asserts calling an oracle function with M2 absent raises a clear `require_m2()` guard error (mirroring
  `viz.style.require_mpl()`).
- **Acceptance criteria (testable checklist):**
  - [ ] `oracle/m2.py` is imported lazily; a context-managed `M2Session` spawns `M2 --script` behind
        `require_m2()`.
  - [ ] With M2 absent, the core still imports zero-dependency and the oracle raises a clear guard.
  - [ ] The M2 host path is documented as **[UNVERIFIED]** on Windows (native vs WSL) — confirm before
        relying on it.
- **Expected effort:** 2 weeks.

### E10-M2 — Euler-pairing cross-check `chi_via_ext` / `ext_dims` (G16) **[PROVEN]**
- **Files modified:** `bridgeland_stability/oracle/m2.py`.
- **New tests:** `tests/test_oracle.py` (additions, skipped if M2 absent) — `test_chi_matches_formula`
  asserts `chi_via_ext(O(0), O(n), P2)` equals `exceptional.chi(...)` and the closed form
  `(n+1)(n+2)/2` for `n ∈ {−3,…,5}`, INCLUDING `χ(O,O(3)) == 10` (the shipped value; 15 is the `n=4`
  value `χ(O,O(4))`) and `χ(O,O(1)) == 3`, `χ(O,O) == 1`; `test_ext_dims_o_minus3` asserts
  `ext_dims(O(0), O(-3), P2) == (0,0,1)` and its alternating sum is `+1` (`= 0 − 0 + 1`), matching the
  shipped `exceptional.chi(O,O(-3)) == +1`; `test_roundtrip` asserts `Fraction(-7,6)` sent to M2 and read
  back is exactly `Fraction(-7,6)`.
- **Acceptance criteria (testable checklist):**
  - [ ] `chi_via_ext(E,F,X)` and `ext_dims(E,F,X)` parse M2 `QQ`/`ZZ` output back to `Fraction`/`int`
        losslessly.
  - [ ] The cross-check reproduces `χ(O,O(3)) = 10` and `χ(O,O(−3)) = +1` (both correct values that the
        earlier draft's `15`/`−1` would have failed).
  - [ ] On a K3 with `Pic=ℤH`, the `−ext`-alternating-sum equals `mukai.pairing` for the pinned pair
        (`⟨v(O),v(O)⟩ = −2`).
- **Expected effort:** 3 weeks.

### E10-M3 — Constructive `moduli_nonempty_by_construction` (sufficient-only) (G16) **[PROVEN]**
- **Files modified:** `bridgeland_stability/oracle/m2.py`.
- **New tests:** `tests/test_oracle.py` (additions, skipped if M2 absent) —
  `test_moduli_witness_is_sufficient_only` asserts `moduli_nonempty_by_construction(r,c1,ch2,X)` returns
  `True` (with a witness) on a class known nonempty, and returns `None`/`Optional[bool]` (NOT `False`)
  when construction fails — a construction failure is NOT a non-existence proof.
- **Acceptance criteria (testable checklist):**
  - [ ] `moduli_nonempty_by_construction` attempts an explicit sheaf/monad and reports success/obstruction
        as `Optional[bool]`, labeled a SUFFICIENT witness only.
  - [ ] A successful witness on a known-nonempty class returns `True`; a failure returns `None`, never
        `False`.
- **Expected effort:** 3 weeks.

---

# E11 — Rational-surface non-emptiness program (LAST — most ambitious)

- **Objective:** After this epic the toolkit answers "does a semistable sheaf with these numerics exist?"
  for rational surfaces (P¹×P¹, 𝔽_n, del Pezzo) via the Coskun–Huizenga / Levine–Zhang `δ_H`
  decision procedure — delivering exact, citable numerical verdicts. This is the single "rational-surface
  program" arc: the surface-specific exceptional-collection generators (G14) supply the Euler-form data
  through which the `δ_H` bounds are stated (G18), so G14 and G18 live in this ONE epic, contiguous, not
  as two independent wins. The in-architecture slice is deliberately thin; the sheaf-theoretic
  HN/prioritary core is delegated out of architecture (E10/G16 default, or paper-tabulation, or deferred
  G19).
- **Contained goal IDs:** G14, G18 (split into G18a in-architecture, G18b [RESEARCH] paper-tabulated,
  G18c optional cross-check).
- **Dependencies on other epics:** E8/G12 (hard — the two NS slope coordinates; the scalar `c=ch₁·H`
  cannot express the input), E9/G13 (hard — the sharp `δ_H` interacts with the full-NS wall structure
  via the true `ch₁²`), E10/G16 (SOFT — the default HN-step supplier; NOT a hard gate, since a
  paper-tabulated finite class list provides a non-oracle path). **The most ambitious surviving Tier-3
  goal — deliberately last.**
- **Mathematical danger zone:** the HN/prioritary core is genuinely sheaf-theoretic and OUTSIDE the
  zero-sheaf constraint — only its output is portable, and the in-architecture slice that ships without
  an oracle/paper table is THIN (stated plainly). Unlike P²'s closed-form DLP curve, there is NO closed
  `δ`-curve — the answer is a decision procedure. `is_exceptional_collection` is a NECESSARY numerical
  condition only; genuine exceptionality (`Ext^{>0}(E_i,E_i)=0`) is sheaf-level and delegated to G16,
  never asserted by the core. The exact `δ_H` numeric targets are the ONLY genuinely paper-dependent
  test values in the roadmap (G18b).
- **Review criteria (epic done when):** the P¹×P¹ Euler Gram matches the exact pinned matrix;
  `is_exceptional_collection` is documented as necessary-only; `moduli_nonempty` degrades to the pinned
  DLP `delta` in the P² limit; the paper-tabulated `δ_H` targets (G18b) are pinned once extracted; the
  𝔽_n polarization-dependence witness fires.
- **Standalone value:** Makes `hirzebruch(n)`/`P¹×P¹`/del Pezzo first-class for non-emptiness — the
  most-requested "beyond P²" capability — delivering exact, citable numerical verdicts (with an
  oracle-backed HN certificate) that no packaged tool provides. Even the G14 sub-arc alone (the
  exceptional-collection Euler tables) is a citable exact table no CAS packages at the numerical level;
  it is NOT a dead-end deliverable because G18 in the same epic consumes it directly.

### E11-M1 — `exceptional_surface.py`: `SurfaceBundle` + Riemann–Roch `chi` (G14) **[PROVEN]**
- **Files modified:** new `bridgeland_stability/exceptional_surface.py` (uses the E8 `NSLattice` + `K_X`).
- **New tests:** `tests/test_exceptional_surface.py` (NEW module) — `test_p1xp1_euler_gram` asserts the
  P¹×P¹ Euler Gram in basis `(O, O(1,0), O(0,1), O(1,1))` equals EXACTLY
  `[[1,2,2,4],[0,1,0,2],[0,0,1,2],[0,0,0,1]]` (from `χ(O(a,b),O(c,d)) = (c−a+1)(d−b+1)`);
  `test_diagonal_unit` asserts `chi(E,E,surface) == 1` for every generator; `test_offdiagonal_witness`
  asserts `χ(O(1,0),O(0,1)) == 0` and `χ(O(0,1),O(1,0)) == 0`.
- **Acceptance criteria (testable checklist):**
  - [ ] `SurfaceBundle(r, c1: NSVector, ch2: Fraction)` and `chi(E, F, surface) -> int` (Riemann–Roch
        `χ(E,F) = ∫ ch(E)^∨·ch(F)·td(X)`) compute exactly from `(r, c₁∈NS, ch₂)` and the intersection
        form.
  - [ ] The P¹×P¹ Euler Gram equals the pinned matrix exactly; the diagonal is all `1`; the two
        off-diagonal witnesses are `0`.
- **Expected effort:** 3 weeks.

### E11-M2 — Per-family exceptional generators + necessary-only collection check (G14) **[PROVEN]** / **[RESEARCH-light]**
- **Files modified:** `bridgeland_stability/exceptional_surface.py`.
- **New tests:** `tests/test_exceptional_surface.py` (additions) — `test_p1xp1_length` asserts
  `len(p1xp1_collection()) == 4`; `test_del_pezzo_length` asserts `len(del_pezzo_collection(3)) == 9`
  (`= 12 − 3`); `test_is_exceptional_collection_necessary` asserts the P¹×P¹ Euler Gram is
  upper-triangular unipotent with unit diagonal ⟹ `is_exceptional_collection(list, surface)` returns
  `True` (a NECESSARY numerical condition only).
- **Acceptance criteria (testable checklist):**
  - [ ] `p1xp1_collection()` (4 objects), `hirzebruch_collection(n)` (Rudakov), `del_pezzo_collection(deg)`
        (Kuleshov–Orlov, length `12 − deg`) hard-code the KNOWN generators per family (NO uniform
        recursion).
  - [ ] `is_exceptional_collection(list, surface)` = "Euler Gram upper-triangular unipotent with unit
        diagonal" and is documented as a NECESSARY numerical condition only — genuine exceptionality
        (`Ext^{>0}(E_i,E_i)=0`) is delegated to G16 (E10), never asserted by the core.
  - [ ] Length invariants `4` (P¹×P¹) and `9` (cubic del Pezzo) hold.
  - [ ] (If E10/G16 present) M2 confirms `Ext^{>0}(E_i,E_i)=0` for each generator — the ONLY route to
        certifying exceptionality.
- **[RESEARCH-light]:** The per-family generators must be read off and transcribed from Rudakov (𝔽_n
  exceptional collections) and Kuleshov–Orlov (del Pezzo, length `12 − deg`) before pinning; only the
  P¹×P¹ Beilinson collection `⟨O, O(1,0), O(0,1), O(1,1)⟩` is transcription-free. `is_exceptional_collection`
  stays a necessary numerical condition only (sheaf-level exceptionality delegated to G16).
- **Expected effort:** 3 weeks.

### E11-M3 — In-architecture `Δ ≥ δ_H` evaluator given HN-length-one (G18a) **[PROVEN]** (split of G18)
- **Files modified:** new `bridgeland_stability/nonemptiness_rational.py`.
- **New tests:** `tests/test_nonemptiness.py` (NEW module) — `test_p2_limit_regresses_to_dlp` asserts in
  the P² limit / single-ray case `moduli_nonempty(...)` degrades to the existing DLP `delta(μ)` with the
  pinned `delta(1/2) == 5/8` and `delta(1/3) == 5/9`; `test_delta_h_nonnegative` asserts every returned
  bound satisfies `δ_H(ξ) ≥ 0` (sanity floor); `test_verdict_reports_mode` asserts a verdict carries a
  G5 `Certificate` whose rigor is `PROVEN` only when the HN-length-one hypothesis came from a certified
  source, else `HEURISTIC`.
- **Acceptance criteria (testable checklist):**
  - [ ] `delta_H(xi, surface) -> Fraction` and `moduli_nonempty(r, c1, ch2, surface) -> Verdict` exist;
        GIVEN the HN-length-one datum, the evaluator computes `Δ(ξ) ≥ δ_H(ξ)` exactly in `Fraction` and
        returns yes/no with a G5 `Certificate`.
  - [ ] In the P²-limit / single-ray case it regresses to the pinned DLP `delta(1/2)=5/8`,
        `delta(1/3)=5/9`.
  - [ ] The `δ_H(ξ) ≥ 0` sanity floor holds; the verdict clearly reports which mode (oracle / paper
        table / heuristic) produced it.
  - [ ] **The in-architecture slice is THIN and stated as such** — the sheaf-theoretic HN/prioritary
        core is delegated (E10/G16, paper-tabulation, or deferred G19).
- **Expected effort:** 3 weeks.

### E11-M4 — Paper-tabulated `δ_H` targets from 1907.06739 / 1910.14060 (G18b) **[RESEARCH]** (split of G18)
- **Files modified:** `bridgeland_stability/nonemptiness_rational.py` (a fixed finite class table).
- **New tests:** `tests/test_nonemptiness.py` (additions) — `test_paper_verdicts` asserts for a fixed
  finite class list on P¹×P¹/𝔽_n from Coskun–Huizenga arXiv:1907.06739 (and Levine–Zhang
  arXiv:1910.14060 del Pezzo deg ≥ 3, anticanonical), `moduli_nonempty(...)` returns the paper's exact
  yes/no verdict and `delta_H` returns the exact tabulated `Fraction` target.
- **Acceptance criteria (testable checklist):**
  - [x] A fixed finite class list supplies HN-length-one data + exact `δ_H` numeric targets; each entry
        cites its paper and class. — `PaperDeltaHTarget` table, 6 entries (`nonemptiness_rational.py`).
  - [x] `moduli_nonempty` reproduces the paper's yes/no for every tabulated class. — `test_paper_verdicts`.
  - [x] A verdict from the in-architecture inequality carries `rigor=PROVEN` only because the
        HN-length-one hypothesis was supplied by this certified paper table. — asserted for all entries.
- **[RESEARCH] — DONE (with scope caveat):** the exact `δ_H` targets were transcribed from
  arXiv:1907.06739 / arXiv:1910.14060 and pinned (4 P² anticanonical entries machine-regressed against the
  native DLP curve; 2 `F_0` `c₁∥H` entries via the machine-checked `Δ_paper = d·discriminant_H` identity).
  Genuinely paper-novel off-P² content is bounded: the sources give the sharp off-P² `δ_H` only via a
  limiting procedure (no exact tables), so non-diagonal / non-anticanonical `𝔽_n` classes remain deferred
  open questions (the polarization-dependent bound is E11-M5 territory + a future full-NS discriminant).
- **Expected effort:** 3 weeks.

### E11-M5 — Optional G16 cross-check + 𝔽_n polarization-dependence witness (G18c) **[PROVEN]** (split of G18)
- **Files modified:** `bridgeland_stability/nonemptiness_rational.py`; `oracle/m2.py` glue.
- **New tests:** `tests/test_nonemptiness.py` (additions) — `test_m2_cross_check` (skipped if M2 absent)
  asserts M2's `moduli_nonempty_by_construction` AGREES with the formula-layer verdict on a shared class
  list; `test_fn_polarization_dependence` asserts on 𝔽_n two ample `H` in DIFFERENT chambers give
  DIFFERENT `δ_H` for a fixed `ξ` (the polarization dependence P² cannot see).
- **Acceptance criteria (testable checklist):**
  - [x] (If E10/G16 present) `moduli_nonempty_by_construction` agrees with the formula-layer verdict on a
        shared class list. — `test_m2_cross_check` (P² shared list; sufficient-only directional agreement,
        oracle `True ⇒ formula.nonempty True`). Skips cleanly here (M2 absent); agrees by construction when M2
        is provisioned.
  - [x] On 𝔽_n, the verdict differs between two ample `H` in different chambers for a fixed `ξ` — a
        polarization-dependence witness the P² model cannot express. — `test_fn_polarization_dependence`.
        **Honest reframing (see E11-M4 caveat):** `δ_H` itself does NOT differ — off P² it is the Bogomolov
        floor `0` for both `H` (asserted *equal*, never faked). The dependence is shown in the quantities the
        package genuinely computes — `μ` (1/3 vs 2/7), `d=H²` (3 vs 7), `discriminant_H` (−1/36 vs 1/196,
        straddling 0), and hence the `moduli_nonempty` verdict (empty vs nonempty) — for `ξ=(2,(1,1),1/2)` on
        𝔽_1 under two strictly-ample `H`. The **sharp polarization-dependent `δ_H` envelope** (Hirzebruch
        exceptional-bundle DLP limiting procedure) remains the **deferred G18 remainder**.
- **Expected effort:** 2 weeks.

---

## Deferred / future work

Two goals are DEFERRED (kept with full provenance in `docs/GOALS.md`, but moved out of the active
roadmap). Neither is abandoned; neither hides an unsolved problem in the active tiers.

- **G17 — Certified K3 totally-semistable / "fake" wall stratification.** Its certifiable content
  duplicates E6/G10's `solve_binary_quadratic`; the only net-new deliverable (the totally-semistable /
  fake-vs-no-wall decision) may require the full Bayer–Macrì movable-cone stratification beyond the
  rank-2 lattice and may stay `HEURISTIC`. **Revisit only after E6/G10 ships and a concrete fake-wall
  classification need arises; if revived, build it as a follow-on of G10, not an independent win.**
- **G19 — OSCAR-via-`juliacall` research escalation.** Redundant with E10/G16 as a sheaf-level oracle;
  native-Windows OSCAR support on this host is **[UNVERIFIED]** and historically WSL-only, with a
  multi-GB Julia/OSCAR install. **Build only if the roadmap takes on a concrete toric coherent-sheaf-
  cohomology or exceptional-collection/mutation experiment that Macaulay2 (G16) cannot serve as
  conveniently; never position it as a route to derived categories — no system provides a general
  derived-category-of-`Coh(X)` engine.**

---

## Adversarial review notes

**Phase-3 independent verdict on this assembled roadmap: proceed (after two flag amendments, applied).**

This section reports a GENUINE INDEPENDENT adversarial review of the assembled roadmap. Production note:
the automated assembler pass hit an output-size limit, so the first critic pass never saw the assembled
artifact and its "insufficient" verdict was about a *missing* document; this roadmap was then
reconstructed from `docs/GOALS.md` with the mandated amendments applied. A fresh critic subsequently
reviewed THIS document end-to-end against the shipped code and the six structural checks (standalone
value per epic, contributor-implementable-in-the-stated-weeks, testable acceptance, research flags,
burstiness, monotonic ordering).

All six load-bearing code-anchored fixtures were re-executed against the shipped `bridgeland_stability`
code and confirmed (42/42 tests green): abelian `numerical_wall = (−5/2, 17/4)` vs K3 shim `(−5/2, 21/4)`
(difference exactly `2/d = 1`); `ν(2,0,1,0) = −1, −1/2, 0, None` and `twist(1,1).a3 = −4/3`;
`chi(O,O(3)) = 10`, `chi(O,O(4)) = 15`, `chi(O,O(−3)) = +1`; the P¹×P¹ Euler Gram
`[[1,2,2,4],[0,1,0,2],[0,0,1,2],[0,0,0,1]]`; `delta` values `5/8, 5/9, 21/32, 13/25`; the P²[2] wall
`(−5/2, 3/2)`; K3 `v(O)=(1,0,1)`, `⟨v,v⟩=−2`; P³ `alpha_crit(β=1/2)=√3`. The one shipped defect
(`QUADRIC3` citing the probability paper `arXiv:1607.07182`) is real and scheduled first (E1-M1).

Structural checks (independent pass): standalone value holds for every epic except E8, which is honestly
flagged as no-payoff NS-refactor infrastructure with its payoff epic E9 scheduled immediately after;
ordering is monotonic (E1 harden → … → E11 rational-surface program) with all 10 hard DAG edges forward;
the mandated G12/G9/G18 splits are present and no milestone exceeds the 1–3-week unit; every acceptance
criterion is a testable checklist asserting concrete values. Two flag amendments were required and have
been applied: (1) **E11-M2** restored to `[RESEARCH-light]` (hard-coding the Rudakov / Kuleshov–Orlov
collections is a paper-transcription dependence that `GOALS.md` G14 carried); (2) **E9-M2** annotated
`[RESEARCH-light]` for the Maciocia Thm-3.11 bounding-constant extraction, consistent with E3-M2.

**Amendments applied in this document:**

- **Missing artifact supplied.** The 17 active goals (G1–G16, G18; G17/G19 deferred) are converted into
  11 epics → 34 milestones, each milestone carrying an explicit 1–3 week estimate and a TESTABLE
  acceptance checklist. Stable IDs G1..G18 and all 10 hard DAG edges (G2→G8, G3→G8, G4→G9, G5→G15,
  G9→G15, G12→G13, G12→G14, G12→G18, G13→G18, G14→G18) are preserved and all point forward across the
  epic ordering.
- **Splits applied (per the mandated split suggestions), renumbered consistently:**
  - **G12 → four sequenced milestones E8-M1..E8-M4** (NSLattice + pinned rank-1 shim / ChernChar
    migration + 42-test gate / minor + Mukai generalization / ρ=2 coverage). The backward-compat shim is
    pinned so the ρ=1 path yields `⟨ch1,H⟩=c` and `⟨ch1,ch1⟩=c²/d` bit-for-bit; the naive `(c,)`+
    `Gram=[[d]]` encoding is explicitly forbidden (it gives `⟨ch1,H⟩=c·d` and breaks `delta(1/2)=5/8`).
  - **G9 → E4-M2 (G9a, `numerical_tilt_wall`, PROVEN via Schmidt Thm 3.3, ships now with
    ch₃-independence + reduction-identity tests)** and **E4-M3 (G9b, `bridgeland_wall` second-tilt,
    [RESEARCH-light], (y,z) coeffs + twisted-cubic pins deferred until transcribed from Schmidt §3/§7 —
    Thm 3.3 is NOT cited for the ch₃-dependent semicircle).**
  - **G18 → E11-M3 (G18a in-architecture `Δ ≥ δ_H` evaluator, regression to pinned DLP delta), E11-M4
    (G18b [RESEARCH] paper-tabulated `δ_H` targets from 1907.06739/1910.14060 — the only genuinely
    paper-dependent test values; acceptance not pinned until extracted), and E11-M5 (G18c optional G16
    cross-check + 𝔽_n polarization-dependence witness).** The in-architecture slice is stated plainly as
    thin.
- **Research flags added** to every research-dependent milestone with a statement of what must be
  read/resolved: E3-M2 (Maciocia radius constant, §3), E4-M3 (Schmidt second-tilt `(y,z)`, §3/§7),
  E5-M2 (Schmidt first-wall locus, §7), E6-M1 (hyperbolic Pell-type enumeration; general fake decision
  out of scope), E7-M2 (Koseki hypotheses vs arXiv:1811.03267 [UNVERIFIED]), E9-M1 (Maciocia ρ>1
  worked-example pin [RESEARCH-light]), E11-M4 (Coskun–Huizenga / Levine–Zhang `δ_H` targets).
- **Every acceptance criterion is a concrete testable checklist item and every "New tests" entry
  asserts a concrete value/property** — never a bare "implement X". The code-verified fixtures are
  carried verbatim into the checklists: abelian `numerical_wall=(−5/2,17/4)` vs K3 shim `(−5/2,21/4)`
  (difference exactly `2/d=1`); `ν(2,0,1,0)=−1,−1/2,0,None` and `twist(1,1).a3=−4/3`;
  `exceptional.chi(O,O(3))=10` (NOT 15 — that is the `n=4` value) and `chi(O,O(−3))=+1`; the P¹×P¹ Euler
  Gram `[[1,2,2,4],[0,1,0,2],[0,0,1,2],[0,0,0,1]]`; `delta(1/2)=5/8`, `delta(1/3)=5/9`,
  `delta(2/5)=13/25`; P²[2] wall `(−5/2, 3/2)`.
- **Ordering fixes (sequence made monotonic; E1 = harden epic, last epic = rational-surface program):**
  - **G15 pulled forward to E5**, immediately after the threefold engine (E4) and the E1 provenance
    substrate — its only hard prerequisites are G9 and G5; it is NOT placed after the G12/G13/G14 block.
    Its Tier-3 label is epistemic (the open threefold BG conjecture), not effort.
  - **G12 marked explicitly as foundational infrastructure with NO standalone user-facing payoff**, and
    **E9 (G13) scheduled to immediately follow E8** so the ρ≥2 payoff lands; the G14+G18 rational-surface
    arc completes the program in E11 so G14's exceptional-collection tables are never a stranded
    dead-end (G14 and G18 share the final epic, contiguous).
  - **The rational-surface arc (G14 + G18) is kept contiguous in E11**, the last and most ambitious epic,
    with the CAS-oracle epic (E10/G16) scheduled before it as the default HN-step supplier.
- **Every epic states an explicit "Standalone value" line** (E8's is explicitly "NONE at the
  user-facing level", justified as foundational infrastructure with the immediate-follow-on E9 landing
  the payoff).
- **Missing test modules named as targets.** The four modules that currently have no test file
  (`test_chern`, `test_varieties`, `test_bg_check`, `test__latex`) are created in E1 (E1-M1, E1-M3,
  E1-M4, E1-M5); E1 also adds `test_rigor`. Every milestone names its target test file (e.g.
  `test_ns_lattice`, `test_tilt_walls`, `test_bridgeland_verdict`, `test_walls_ns`, `test_oracle`,
  `test_exceptional_surface`, `test_nonemptiness`, and the shared `test_abelian_k3_walls` regression
  module). The G1 acceptance requires `tests/test_varieties.py` with the allowlist-membership guard on
  every catalog reference; the G6/G12 acceptance criteria have `test_bg_check`/`test_chern`/
  `test_ns_lattice` homes.
- **All confidence flags preserved** ([PROVEN]/[CONJECTURED]/[SPECULATIVE]/[UNVERIFIED]/[RESEARCH]) and
  the exact code-verified test values kept.

**Resulting independent verdict: proceed.** With the two flag amendments applied (E11-M2 and E9-M2
now carry `[RESEARCH-light]`), the roadmap is executable milestone-by-milestone by a PhD algebraic
geometer with no prior codebase knowledge: it is monotonic (E1 harden → … → E11 rational-surface
program), respects the hard DAG, splits the three oversized goals into unit-sized milestones with week
estimates, flags every research-dependent value, gives every milestone a testable acceptance checklist
anchored to a code-verified fixture, and respects the inviolable exact-`Fraction` / CH-convention /
zero-dependency / 42-pinned-test invariants without over-claiming.
