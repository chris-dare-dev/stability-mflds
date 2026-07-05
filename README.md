# bridgeland_stability

Exact-arithmetic tools for three tightly-coupled objects in the theory of
moduli of sheaves on projective surfaces and threefolds:

1. the **Drézet–Le Potier curve** δ(μ) on P² and the exceptional bundles that
   shape it (the existence boundary for stable sheaves);
2. the **Bogomolov–Gieseker** positivity boundary (surfaces and threefolds);
3. **Bridgeland stability walls** in the (s, t) upper half-plane;
4. the **K3 Mukai lattice** and Bayer–Macrì wall classification.

All core computations are exact (`fractions.Fraction`); floats appear only for
geometry and plotting.

> ⚠️ This package corrects several substantive mathematical errors in the
> originating project brief (the exceptional-bundle ranks are **Markov numbers**,
> not Farey; the DLP curve is a fractal cusp-envelope with bundles *below* it,
> not on it; several test values are wrong). See **[docs/CORRECTIONS.md](docs/CORRECTIONS.md)**
> for the full list with citations.

## Install

```bash
pip install -e .          # core (pure standard library)
pip install -e ".[viz]"   # + plotly / matplotlib for figures
pip install -e ".[dev]"   # + pytest
```

Pure-stdlib core: the algorithms need only `fractions`. Python ≥ 3.9.

## Convention

Discriminant is the Coskun–Huizenga normalization (the literature standard):

```
mu    = (ch1 . H) / (r * d),   d = H^2
Delta = (1/2) mu^2 - ch2 / (r * d)
```

In it, line bundles have Δ=0, an exceptional bundle of rank r has
Δ = ½(1 − 1/r²) < ½, and the DLP curve takes values in [½, 1].
`ChernChar.discriminant_brief` returns `2·Δ` (the brief's doubled convention).

## Quickstart

```python
from fractions import Fraction as F
import bridgeland_stability as bs

# 1. Exceptional bundles on P^2 (Markov ranks 1,2,5,13,...)
bundles = bs.enumerate_exceptional(0, 1, R_max=30)
print([(b.r, b.slope) for b in bundles[:5]])     # rank 1 at 0, 13 at 5/13, 5 at 2/5, ...
bs.is_exceptional(bs.Bundle.from_slope(F(1, 3)))  # False: no rank-3 exceptional bundle

# 2. Drézet–Le Potier curve
d = bs.delta(F(1, 2), bundles)                    # 5/8  (cusp of T(-1))
bs.moduli_nonempty(bs.Bundle(2, 1, F(-5, 2)))     # positive-dimensional moduli

# 3. Bridgeland walls (exact)
v, w = bs.ChernChar(1, 0, F(-2)), bs.ChernChar(1, -1, F(1, 2))   # P^2[2], destabilizer O(-1)
wall = bs.numerical_wall(v, w, bs.P2.d)
print(wall.center, wall.radius)                   # -5/2  1.5   (matches ABCH)

# ...and the *actual* (certified-necessary) walls, finite + de-spuriated:
walls, complete = bs.actual_walls_complete(v, bs.P2)
print(len(walls), complete)                       # 1 True   (the unique ABCH wall)

# 4. Threefold BG boundary
v3 = bs.ThreefoldChern(2, 0, 1, 0)                # P^3 null-correlation bundle
bs.alpha_crit(v3, F(1, 2), bs.P3.d3)              # sqrt(3)

# 5. K3 Mukai lattice
vO = bs.MukaiVector.from_chern(r=1, l=0, ch2=0)   # (1,0,1)
bs.self_pairing(vO, d=2)                          # -2  (spherical)

# Figures (needs the [viz] extra)
from bridgeland_stability import viz
viz.plot_dlp_curve(mu_min=0, mu_max=1, R_max=40, html_path="dlp.html")
```

See [`examples/demo.py`](examples/demo.py) for a script that reproduces every
validated value and writes the figures in `figures/`.

## Package layout

| module | contents |
|---|---|
| `chern` | `ChernChar` (surface), slope/discriminant/twist/central charge |
| `varieties` | `Surface`/`Threefold` data classes + catalog (P², P¹×P¹, K3, P³, quintic, …) |
| `exceptional` | Algorithm 1: exceptional bundles (ε/Markov recursion), Riemann–Roch χ |
| `dlp` | Algorithm 2: the DLP curve δ(μ) and the non-emptiness criterion |
| `walls` | Algorithm 3: `numerical_wall`, `compute_walls`, `actual_walls` (certified), `walls_from_subobjects` |
| `bg_check` | Algorithm 4: Bogomolov–Gieseker (surface) |
| `threefold` | Algorithm 5: tilt-BG boundary α_crit(β) (proven cases flagged) |
| `mukai` | K3 Mukai vectors, pairing, Bayer–Macrì wall classification |
| `viz` | plotly (interactive) + matplotlib (static) rendering |

## Validation

`pytest` (42 tests) pins every documented value, e.g. δ(½)=5/8, δ(1/3)=5/9,
the P²[2] wall (−5/2, 3/2), α_crit(β=½)=√3 for the P³ null-correlation bundle,
⟨v(O),v(O)⟩=−2 on K3, and that ranks 3, 4 are not Markov numbers.

```bash
pytest -q
```

## References

Drézet–Le Potier (Ann. Sci. ENS 18, 1985); Coskun–Huizenga (Gökova survey);
Arcara–Bertram–Coskun–Huizenga ([1203.0316](https://arxiv.org/abs/1203.0316));
Maciocia ([1202.4587](https://arxiv.org/abs/1202.4587)); Bayer–Macrì
([1301.6968](https://arxiv.org/abs/1301.6968)); Bayer–Macrì–Toda
([1103.5010](https://arxiv.org/abs/1103.5010)); C. Li
([1510.04089](https://arxiv.org/abs/1510.04089),
[1810.03434](https://arxiv.org/abs/1810.03434)); Veselov
([2501.06779](https://arxiv.org/abs/2501.06779)). Full per-claim citations in
[docs/CORRECTIONS.md](docs/CORRECTIONS.md).
