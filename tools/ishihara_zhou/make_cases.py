"""Generate RANS cases for the Ishihara-Zhou wind-tunnel ridge and hill.

Usage: python3 make_cases.py [output directory, default ./cases]

Data: Ishihara & Zhou, "Wind tunnel measurement dataset of turbulent flow over
two-dimensional ridges and three-dimensional hills with smooth and rough
surfaces", Data in Brief 63 (2025) 112260, CC BY 4.0. Download the approach
flow files mmc1.xls (smooth) and mmc2.xls (rough) into ./data (see README.md).

Wind tunnel: test section 7.0 x 1.1 x 0.9 m, boundary layer ~280 mm deep under
a 5.9 m/s free stream. 2D ridge z = h cos^2(pi x / 2L) and 3D hill
z = h cos^2(pi r / 2L) for |x|, r <= L, h = 40 mm, L = 100 mm (maximum slope
about 32 degrees). Cases are written at full scale (x1000): h = 40 m,
L = 100 m, tunnel height 900 m; velocities are not scaled.

Surfaces (both geometries share the approach flows): smooth z0 = 0.01 m,
u* = 0.21 m/s (mmc1); rough z0 = 0.3 m, u* = 0.32 m/s, 5 mm artificial grass
(mmc2). The inflow is the measured approach profile: U, and
k = (su^2 + sv^2 + sw^2) / 2. Below the lowest measurement the log law with the
paper's u* and z0 is used (capped at the lowest measured U, k held); above the
highest point the free stream is held.

Cases, one folder each with case.inp, terrain.amrwind, inflow_profile.txt
(z u v T tke) and rans_1d.info (z u v w tke):
  approach/<surface>_<grid>          flat thin domain: drift of the inflow over
                                     the fetch, compare with mmc1 / mmc2
  ridge2d/<surface>_<model>_<grid>   2D ridge, 4 periodic cells in y (mmc3 / mmc4)
  hill3d/<surface>_<model>_<grid>    3D hill, periodic in y over 1152 m, close
                                     to the 1.1 m tunnel width (mmc5-mmc8)

Grids (dx = dy, dz): 2D 4x4, 2x2 and 8x1 m; 3D 8x8 and 4x4 m. A grid is named
"4m" when dx = dz and "8x1m" otherwise. Anisotropic grids get 16 pre- and
post-smoothing sweeps in the MAC and nodal projections so MLMG stays stable.
2D terrain files are sampled every min(dx, dz), which TerrainDrag interpolates
bilinearly; 3D terrain files every dx.

Domain: x 0-2400 m (crest or hill center at x = 800 m), z 0-896 m with a slip
top (tunnel ceiling). Mass inflow at xlo from the TabulatedProfile boundary
condition, pressure outflow at xhi, wall model on the floor, TerrainDrag with
uniform roughness z0, KLAxell. No Coriolis, geostrophic forcing or buoyancy;
ABL.meso_sponge_start is above the domain top, so KLAxell uses its neutral
length scale everywhere and the TKE sponge never acts.

Other turbulence models can reuse the setup through generate(out, models),
where models maps a folder name to the input keys that select the model.
"""

import math
import pathlib
import sys

import xlrd

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
HEIGHT, HALF_LENGTH, XC = 40.0, 100.0, 800.0
LX, LZ = 2400.0, 896.0
NY_2D = 4
LY_3D, YC_3D = 1152.0, 576.0
KAPPA, T0, RHO = 0.41, 300.0, 1.225
SURFACES = {
    "smooth": {"file": "mmc1.xls", "z0": 0.01, "ustar": 0.21},
    "rough": {"file": "mmc2.xls", "z0": 0.3, "ustar": 0.32},
}
GRIDS_2D = ((4.0, 4.0), (2.0, 2.0), (8.0, 1.0))
GRIDS = {"approach": GRIDS_2D, "ridge2d": GRIDS_2D, "hill3d": ((8.0, 8.0), (4.0, 4.0))}
STOP = {"approach": 900.0, "ridge2d": 1500.0, "hill3d": 1500.0}
ANISOTROPIC_SMOOTHING = 16
MODELS = {"klaxell": {"turbulence.model": "KLAxell"}}


def approach_profile(name):
    """Measured approach profile as rows (z [m], U [m/s], k [m^2/s^2])."""
    sheet = xlrd.open_workbook(DATA / SURFACES[name]["file"]).sheet_by_index(0)
    rows = []
    for r in range(sheet.nrows):
        z = sheet.cell_value(r, 0)
        if isinstance(z, float):
            u, su, sv, sw = (sheet.cell_value(r, c) for c in range(1, 5))
            rows.append((z, u, 0.5 * (su * su + sv * sv + sw * sw)))
    rows.sort()
    z0, ustar = SURFACES[name]["z0"], SURFACES[name]["ustar"]
    zmin, umin, kmin = rows[0]
    below = []
    for z in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
        if z < zmin:
            u = min(max(ustar / KAPPA * math.log(z / z0), 0.0), umin) if z > z0 else 0.0
            below.append((z, u, kmin))
    ztop, utop, ktop = rows[-1]
    return below + rows + [(LZ, utop, ktop)]


def terrain_height(kind, x, y):
    if kind == "approach":
        return 0.0
    s = abs(x - XC) if kind == "ridge2d" else math.hypot(x - XC, y - YC_3D)
    return HEIGHT * math.cos(0.5 * math.pi * s / HALF_LENGTH) ** 2 if s <= HALF_LENGTH else 0.0


def blocking(n):
    """Largest of 8, 4, 2, 1 that divides n cells."""
    b = 8
    while n % b:
        b //= 2
    return b


def grid_name(dx, dz):
    return f"{dx:g}m" if dx == dz else f"{dx:g}x{dz:g}m"


def write_case(out, kind, surface, model, model_keys, dx, dz):
    grid = grid_name(dx, dz)
    name = f"{surface}_{grid}" if kind == "approach" else f"{surface}_{model}_{grid}"
    dest = out / kind / name
    dest.mkdir(parents=True, exist_ok=True)
    profile = approach_profile(surface)
    with open(dest / "inflow_profile.txt", "w") as f:
        f.write("# z u v temperature tke\n")
        f.writelines(f"{z:.3f} {u:.5f} 0.0 {T0:.1f} {k:.6f}\n" for z, u, k in profile)
    with open(dest / "rans_1d.info", "w") as f:
        f.writelines(f"{z:.3f} {u:.5f} 0.0 0.0 {k:.6f}\n" for z, u, k in profile)

    three_d = kind == "hill3d"
    ly = LY_3D if three_d else NY_2D * dx
    nx, ny, nz = int(LX / dx), int(ly / dx), int(LZ / dz)
    spacing = dx if three_d else min(dx, dz)
    xs = [spacing * i for i in range(int(LX / spacing) + 1)]
    ys = [dx * j for j in range(ny + 1)] if three_d else [0.0, ly]
    with open(dest / "terrain.amrwind", "w") as f:
        f.write(f"{len(xs)}\n{len(ys)}\n")
        f.writelines(f"{x:.3f}\n" for x in xs)
        f.writelines(f"{y:.3f}\n" for y in ys)
        for x in xs:
            f.writelines(f"{terrain_height(kind, x, y):.6f}\n" for y in ys)

    z0 = SURFACES[surface]["z0"]
    stop = STOP[kind]
    keys = {
        "geometry.prob_lo": "0.0 0.0 0.0",
        "geometry.prob_hi": f"{LX:g} {ly:g} {LZ:g}",
        "geometry.is_periodic": "0 1 0",
        "amr.n_cell": f"{nx} {ny} {nz}",
        "amr.max_level": "0",
        "amr.max_grid_size": "128",
        "amr.blocking_factor_x": f"{blocking(nx)}",
        "amr.blocking_factor_y": f"{blocking(ny) if three_d else 1}",
        "amr.blocking_factor_z": f"{blocking(nz)}",
        "time.stop_time": f"{stop:g}",
        "time.max_step": "-1",
        "time.initial_dt": "0.1",
        "time.cfl": "0.5",
        "time.init_shrink": "0.1",
        "time.regrid_interval": "-1",
        "time.plot_interval": "-1",
        # 3D plot files are large: write only the final state
        "time.plot_time_interval": f"{stop:g}" if three_d else f"{stop / 2:g}",
        "time.checkpoint_interval": "-1",
        "io.int_outputs": "terrain_blank terrain_drag",
        "incflo.physics": "ABL TerrainDrag",
        "incflo.density": f"{RHO}",
        "incflo.gravity": "0.0 0.0 -9.81",
        "incflo.velocity": f"{profile[-1][1]:.3f} 0.0 0.0",
        "incflo.verbose": "0",
        "incflo.initial_iterations": "8",
        "incflo.do_initial_proj": "true",
        "incflo.constant_density": "true",
        "incflo.use_godunov": "true",
        "incflo.godunov_type": '"ppm"',
        "incflo.diffusion_type": "2",
        "transport.model": "ConstTransport",
        "transport.viscosity": "1.0e-5",
        "transport.laminar_prandtl": "0.7",
        "transport.turbulent_prandtl": "0.3333",
        "transport.reference_temperature": f"{T0:g}",
        "turbulence.model": "KLAxell",
        "TKE.source_terms": "KransAxell",
        "ABL.kappa": f"{KAPPA}",
        "ABL.normal_direction": "2",
        "ABL.surface_roughness_z0": f"{z0:g}",
        "ABL.surface_temp_flux": "0.0",
        "ABL.wall_shear_stress_type": "local",
        "ABL.perturb_velocity": "false",
        "ABL.perturb_temperature": "false",
        "ABL.initial_wind_profile": "true",
        "ABL.rans_1dprofile_file": '"rans_1d.info"',
        "ABL.meso_sponge_start": "5000.0",
        "ABL.temperature_heights": f"0.0 {LZ:g}",
        "ABL.temperature_values": f"{T0:g} {T0:g}",
        "ABL.stats_output_frequency": "100",
        # No DragTempForcing: the cases are neutral without buoyancy, so the
        # temperature only has to stay uniform, and its explicit terrain-cell
        # relaxation (rate 10 / dz) is unstable on the 1 m vertical grid
        "ICNS.source_terms": "DragForcing",
        "TerrainDrag.terrain_file": '"terrain.amrwind"',
        "TerrainDrag.uniform_roughness": f"{z0:g}",
        "TabulatedProfile.filename": "inflow_profile.txt",
        "xlo.type": '"mass_inflow"',
        "xlo.density": f"{RHO}",
        "xlo.velocity.inflow_type": "TabulatedProfile",
        "xlo.temperature.inflow_type": "TabulatedProfile",
        "xlo.tke.inflow_type": "TabulatedProfile",
        "xhi.type": '"pressure_outflow"',
        "zlo.type": '"wall_model"',
        "zlo.tke_type": '"fixed_gradient"',
        "zlo.tke": "0.0",
        "zhi.type": '"slip_wall"',
        "zhi.temperature_type": '"fixed_gradient"',
        "zhi.temperature": "0.0",
        "zhi.tke_type": '"fixed_gradient"',
        "zhi.tke": "0.0",
    }
    keys.update(model_keys)
    # Absolute tolerances: 1e-6 for the diffusion solves; the projections keep
    # 1e-4, since on the 8 m x 1 m grid the nodal projection stalls near 2e-5
    for solver in ("mac_proj", "nodal_proj", "diffusion", "temperature_diffusion", "tke_diffusion"):
        keys[f"{solver}.mg_rtol"] = "-1"
        keys[f"{solver}.mg_atol"] = "1e-4" if solver.endswith("_proj") else "1e-6"
    if dx != dz:
        for solver in ("mac_proj", "nodal_proj"):
            keys[f"{solver}.num_pre_smooth"] = f"{ANISOTROPIC_SMOOTHING}"
            keys[f"{solver}.num_post_smooth"] = f"{ANISOTROPIC_SMOOTHING}"
    geometry = {"approach": "flat approach flow", "ridge2d": "2D ridge", "hill3d": "3D hill"}[kind]
    header = (
        f"# Ishihara-Zhou {geometry}, {surface} surface, model {model}, dx = {dx:g} m, "
        f"dz = {dz:g} m, full scale (x1000)\n"
    )
    (dest / "case.inp").write_text(header + "".join(f"{k:<46}= {v}\n" for k, v in keys.items()))
    cells = nx * ny * nz
    print(f"{kind}/{name}: n_cell {nx} x {ny} x {nz} = {cells / 1e6:.2f} M, stop {stop:g} s")


def generate(out, models, approach=True):
    """Write every case for the given models ({folder name: input keys})."""
    if approach:
        for surface in SURFACES:
            for dx, dz in GRIDS["approach"]:
                write_case(out, "approach", surface, "klaxell", MODELS["klaxell"], dx, dz)
    for kind in ("ridge2d", "hill3d"):
        for surface in SURFACES:
            for model, model_keys in models.items():
                for dx, dz in GRIDS[kind]:
                    write_case(out, kind, surface, model, model_keys, dx, dz)


if __name__ == "__main__":
    generate(pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "cases", MODELS)
