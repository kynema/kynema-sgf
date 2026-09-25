# KLAxellSeparation on the Ishihara–Zhou ridge and hill

These scripts extend the KLAxell cases of `make_cases.py` and `README.md` to the
`KLAxellSeparation` model. They rely on those two files and on the
`TabulatedProfile` inflow boundary condition (kynema/kynema-sgf#2031).

| File | What it does |
|---|---|
| `separation_cases.py` | Writes the 2D ridge and 3D hill cases for the `limiter` and `all` setups |
| `plot_ridge_overlay.py` | Near-surface U along the ridge plus U and k profiles for KLAxell, `limiter` and `all` |
| `les_agreement.py` | Near-surface U against the LES and the tunnel, and the reversed-flow extent |
| `ref_les_liu2019.py` | LES reference values, digitized |

## Setups

| Folder name | Settings |
|---|---|
| `klaxell` (from `make_cases.py`) | `turbulence.model = KLAxell` |
| `limiter` | `KLAxellSeparation` with `pressure_gradient_sensor`, `realizable_cmu` and `KLAxellSeparation_coeffs.gate_relaxation_time = 10` |
| `all` | `limiter` plus `production_cap`, `destruction_boost`, `curvature_correction` (default `richardson`) and `implicit_dissipation` |

```bash
python3 make_cases.py
python3 separation_cases.py
python3 plot_ridge_overlay.py rough overlay_rough_4m.png 4
python3 les_agreement.py 4
```

`plot_ridge_overlay.py` and `les_agreement.py` read finished runs from
`cases/ridge2d/<surface>_<setup>_<grid>`. The grid argument is `4`, `2` or `8x1`.

## LES reference

`ref_les_liu2019.py` holds values read from the figures of Liu, Diao & Ishihara,
"Study of the flow fields over simplified topographies with different roughness
conditions using large eddy simulations", *Renewable Energy* 136 (2019) 968–992:

- **Station values:** U at 5 mm (smooth) and 7 mm (rough) above the surface from
  their Fig. 8c.
- **Reversed-flow extent:** the same heights, from their separation-bubble
  outlines in Fig. 11.

Their LES models the rough surface as a 5 mm canopy and was validated against
earlier measurements of the same ridge. Reading the figures is accurate to about
0.3–0.5 m/s. U/Uref is converted with the 2025 approach-flow U at 160 mm.

## Results on the 2D ridge, 4 m grid

Reversed flow at 5 m (smooth) and 7 m (rough) above the surface:

| | Smooth: separation | Smooth: reattachment | Rough: separation | Rough: reattachment |
|---|---|---|---|---|
| Tunnel | 0–50 m | 150–200 m | about 50 m | 250–300 m |
| LES | 45 m | 180 m | 37 m | 239 m |
| KLAxell | 66 m | 118 m | 58 m | 126 m |
| `limiter` | 62 m | 138 m | 58 m | 146 m |
| `all` | 62 m | 138 m | 58 m | 150 m |

Station rms of U against the tunnel (m/s):

| | Smooth | Rough |
|---|---|---|
| KLAxell | 0.34 | 0.55 |
| `limiter` | 0.27 | 0.51 |
| `all` | 0.29 | 0.51 |

- **The limiter carries the improvement:** it lengthens the bubble and brings U
  closer to both the tunnel and the LES. The other treatments change little on
  this case.
- **Reattachment is still early.** Near the surface, 150–250 m behind the crest,
  all setups are 0.5–0.9 m/s faster than the tunnel and the LES.
- **Rough crest:** the lowest point (7 mm, just above the grass) is overpredicted
  by every setup and by the canopy LES alike.
