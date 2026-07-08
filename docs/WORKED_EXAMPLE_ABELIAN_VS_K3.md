# Worked example — abelian vs K3 (s, t) walls for the same class

This note walks the same Chern class `v = (1, 0, -2)` through its Bridgeland
(s, t) wall on **both** an abelian surface and a K3 surface, shows the
hand-arithmetic behind every exact `Fraction`, explains the `+2/d` shift and
*why* it happens, attaches a Bayer–Macrì wall **type**, and states plainly what
the current Picard-rank-1 model does **not** cover.

Run it:

```bash
python examples/abelian_vs_k3_walls.py
```

Everything below is exact `fractions.Fraction`; floats never enter the core.

## The shared fixture and its three minors

Take `S = abelian_surface(2)` and `S = K3(2)` (both `d = H² = 2`), the class
`v = (r, c, e) = (1, 0, -2)` and the destabilizer `w = (r', c', e') = (1, -1, 1/2)`.
The three `walls._minors` are

```
W_rc = r c' - r' c  = 1·(-1) - 1·0     = -1
W_re = r e' - r' e  = 1·(1/2) - 1·(-2) = 5/2
W_ce = c e' - c' e  = 0·(1/2) - (-1)·(-2) = -2
```

Since `W_rc ≠ 0` this is a semicircular wall with

```
center   = W_re / W_rc                 = (5/2)/(-1)             = -5/2
radius²  = center² - 2 W_ce/(d·W_rc)   = 25/4 - 2·(-2)/(2·(-1)) = 25/4 - 2 = 17/4
```

That is the genuine **abelian** semicircle `(center, radius²) = (-5/2, 17/4)`,
returned by `abelian_wall(v, w, abelian_surface(2))` (which is exactly
`numerical_wall` on the bare triple — no shift).

## The K3 shift and WHY it is exactly `+2/d`

The Mukai vector is `v(E) = ch(E)·√td`. The Todd-square-root differs by surface:

* **Abelian:** trivial tangent bundle ⇒ `√td = (1, 0, 0)`, so the bare Chern
  triple `(r, c, ch₂)` **is already** the Mukai vector — no shift.
* **K3:** `√td_K3 = (1, 0, 1)`, so `v(E) = (r, c₁, ch₂ + r)` — the `ch₂`
  coordinate is raised by the rank.

The K3 wall applies `e ↦ e + r`, `e' ↦ e' + r'` and recomputes:

```
v' = (1, 0, -1),  w' = (1, -1, 3/2)
W_rc = -1  (unchanged)
W_re = 3/2 - (-1) = 5/2  (unchanged)
W_ce = 0·(3/2) - (-1)·(-1) = -1   (was -2; shifted by -W_rc)
center   = -5/2                              (invariant)
radius²  = 25/4 - 2·(-1)/(2·(-1)) = 25/4 - 1 = 21/4
```

So `k3_wall(v, w, 2) = (-5/2, 21/4)`. Under `e ↦ e + r` the minors satisfy
`W_rc` unchanged, `W_re` unchanged, `W_ce ↦ W_ce - W_rc`; hence the **center is
invariant** and

```
radius²_K3 - radius²_abelian = 21/4 - 17/4 = 1 = 2/d = 2/2.
```

This is a proven identity at Picard rank 1, not a coincidence of this fixture
(see `docs/CORRECTIONS.md` §6; Bridgeland, *Stability conditions on K3
surfaces*, arXiv:math/0307164; Bayer–Macrì, *MMP for moduli of sheaves on K3s*,
arXiv:1301.6968).

**Caveat on this fixture.** `w = (1, -1, 1/2)` has `l = c/d = -1/2 ∉ ℤ`, so it is
**not** a genuine K3 Mukai-lattice class — it is a *synthetic* `+2/d`
demonstration only. The exact `17/4` vs `21/4` numbers are produced solely by
this specific `w`; any integral-`l` (even-`c`) destabilizer yields different
minors and a different radius². Because it is non-integral, `k3_wall_classified`
**raises** on it (`l = -1/2`); it has no Bayer–Macrì type.

## Attaching a Bayer–Macrì wall type (genuine integral-`l` class)

For the wall **type** we use a genuine integral-`l` K3 lattice class:

```
v_g = ChernChar(2, 0, -3)   [Mukai (2, 0, -1)]
w_g = ChernChar(2, 2, -1)   [Mukai (2, 1,  1)]  (spherical, (w,v) = 0)
```

Its (s, t) semicircle, via the same `ch₂ ↦ ch₂ + r` shim (`d = 2`):

```
v' = (2, 0, -1),  w' = (2, 2, 1)
W_rc = 2·2 - 2·0 = 4
W_re = 2·1 - 2·(-1) = 4
W_ce = 0·1 - 2·(-1) = 2
center  = 4/4 = 1
radius² = 1 - 2·2/(2·4) = 1 - 1/2 = 1/2
```

Its Bayer–Macrì type (Thm 5.7), from the Mukai vectors `v_g → (2, 0, -1)`,
`w_g → (2, 1, 1)`:

```
w_g² = d l² - 2 r s = 2·1 - 2·2·1 = -2   (spherical)
(w, v) = d l l' - r s' - r' s = 0 - 2·(-1) - 2·1 = 0
```

A spherical class `s ∈ H` with `(s, v) = 0` ⇒ **divisorial / brill-noether**. So

```
k3_wall_classified(v_g, w_g, K3(2)) = ((center 1, radius² 1/2), divisorial/brill-noether)
```

## Current-model limits (honest)

The example prints these, and they bound what the numbers mean:

* **Picard rank 1 only.** `ρ ≥ 2` abelian (`E₁ × E₂`) and `ρ ≥ 2` K3 are gated
  on the NS-lattice refactor (E8 / G12) and are not covered here.
* **The `ch₂ ↦ ch₂ + r` shim is K3-only.** Never feed abelian input to
  `k3_wall`: an abelian surface has `√td = (1, 0, 0)`, so the shim would inject a
  spurious `+2/d`. Use `abelian_wall` for abelian input, and `numerical_wall`
  for the bare, surface-agnostic primitive.
* **`actual_walls` / `actual_walls_complete` are NOT abelian/K3 completeness
  oracles.** Their `_is_integral_chern` filter tests the `d = 1` P² Chern-lattice
  condition only, and `actual_walls` layers ABCH / P²-specific necessary
  conditions that do not transfer to `d ≠ 1`.
* **The `+2/d` demo `w = (1, -1, 1/2)` (`l = -1/2`) has no Bayer–Macrì type:**
  `k3_wall_classified` requires `c % d == 0` and raises on non-integral `l`.
