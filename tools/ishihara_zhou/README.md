# Ishihara–Zhou ridge and hill cases

Scripts that set up RANS runs of the wind-tunnel dataset of Ishihara & Zhou,
"Wind tunnel measurement dataset of turbulent flow over two-dimensional ridges and
three-dimensional hills with smooth and rough surfaces", *Data in Brief* 63 (2025)
112260 ([doi:10.1016/j.dib.2025.112260](https://doi.org/10.1016/j.dib.2025.112260)).
The cases use the `TabulatedProfile` inflow boundary condition and KLAxell.

## Files

| File | What it does |
|---|---|
| `make_cases.py` | Writes every case folder from the measured approach flows |
| `run_sequential.sh` | Runs case folders one after another, each in place |
| `compare_ridge.py` | Compares a 2D ridge or approach-flow run with the tunnel data |

## Data

The measurements are CC BY 4.0 supplementary files of the paper and are not
stored in this repository. Download them into `data/`:

```bash
mkdir -p data
for n in 1 2 3 4 5 6; do curl -L -o data/mmc$n.xls https://ars.els-cdn.com/content/image/1-s2.0-S2352340925009813-mmc$n.xls; done
for n in 7 8; do curl -L -o data/mmc$n.xlsx https://ars.els-cdn.com/content/image/1-s2.0-S2352340925009813-mmc$n.xlsx; done
```

| File | Contents |
|---|---|
| mmc1.xls / mmc2.xls | Approach flow, smooth / rough surface (used to build the inflow) |
| mmc3.xls / mmc4.xls | 2D ridge, smooth / rough, vertical plane y = 0 |
| mmc5.xls / mmc6.xls | 3D hill, smooth / rough, vertical plane y = 0 |
| mmc7.xlsx / mmc8.xlsx | 3D hill horizontal planes: smooth z/h = 0.125 and 1.0, rough z/h = 0.25 and 1.0 |

## Generating the cases

```bash
python3 -m pip install xlrd
python3 make_cases.py
```

This writes `cases/`, which git ignores. Each folder holds `case.inp`,
`terrain.amrwind`, `inflow_profile.txt` (`z u v T tke`) and `rans_1d.info`
(`z u v w tke`). The two profile files use different column orders.

| Cases | Grid (dx, dz) | Cells |
|---|---|---|
| `approach/<surface>_<grid>`, `ridge2d/<surface>_klaxell_<grid>` | 4 m, 2 m, 8 m × 1 m | 0.54 M, 2.15 M, 1.08 M |
| `hill3d/<surface>_klaxell_<grid>` | 8 m, 4 m | 4.84 M, 38.7 M |

Surfaces are `smooth` (z0 = 0.01 m, u* = 0.21 m/s) and `rough` (z0 = 0.3 m,
u* = 0.32 m/s, 5 mm artificial grass). Cases are at full scale (tunnel ×1000):
the ridge or hill is 40 m high with a 100 m half-width.

- **Domain:** x from 0 to 2400 m with the crest or hill center at x = 800 m;
  z from 0 to 896 m under a slip top, standing in for the tunnel ceiling.
- **Boundaries:** mass inflow from the measured profile, pressure outflow,
  wall model on the floor, TerrainDrag with uniform z0.
- **Physics:** no Coriolis force, geostrophic forcing or buoyancy.
- **3D hill:** y is periodic over 1152 m. The tunnel is 1.1 m wide with side
  walls, which the cases do not include.
- **Anisotropic grids** (dx ≠ dz) use 16 pre- and post-smoothing sweeps in the
  MAC and nodal projections so the multigrid solves stay stable.

## Running

```bash
EXE=/path/to/kynema_sgf NP=4 ./run_sequential.sh cases/approach/rough_4m cases/ridge2d/rough_klaxell_4m
```

`ARGS` passes extra ParmParse overrides to every run. Run the approach cases
first to confirm that the inflow holds over the fetch.

## Comparing with the tunnel

```bash
python3 compare_ridge.py cases/approach/rough_4m rough empty
python3 compare_ridge.py cases/ridge2d/rough_klaxell_4m rough ridge rough_klaxell_4m.png
```

This needs `yt`, `xlrd`, `numpy` and `matplotlib`. For each measurement station
it prints the rms U and k differences and the lowest-point velocity, then the
reversed-flow extent at the first cell above the surface.

On the rough surface the lowest measured points (7 mm) sit just above the 5 mm
grass, so a wall-function roughness model is not expected to match them closely.

## Regression test

`test/test_files/terrain_ishihara_zhou_ridge` is the rough ridge on a coarse
16 m grid. It runs the tabulated inflow of velocity, temperature and TKE together
with TerrainDrag and KLAxell.
