"""Compare Ishihara-Zhou ridge runs with the wind-tunnel data.

Usage: python3 compare_ridge.py <run dir> <smooth|rough> <empty|ridge> [png]

Full scale is the tunnel x1000, so tunnel mm equal model m. The measured x is
relative to the crest (x = 800 m in the runs) and z is height above the tunnel
floor; velocities are not scaled. The model k is compared with
k = (su^2 + sv^2 + sw^2) / 2 from the data.

empty: profiles at x = 150, 400 and 650 m (the ridge's -150 mm station)
against the approach profile, as the inflow drift over the fetch.
ridge: profiles at every station against the ridge data, the reversed-flow
extent at the first cell above the surface, and the reattachment point.
"""

import pathlib
import sys

import numpy as np
import xlrd
import yt

yt.set_log_level(50)

DATA = pathlib.Path(__file__).resolve().parent / "data"
XC = 800.0
FILES = {
    ("smooth", "approach"): "mmc1.xls",
    ("rough", "approach"): "mmc2.xls",
    ("smooth", "ridge"): "mmc3.xls",
    ("rough", "ridge"): "mmc4.xls",
}


def read_sheet(name):
    sheet = xlrd.open_workbook(DATA / name).sheet_by_index(0)
    header = [str(sheet.cell_value(1, c)) for c in range(sheet.ncols)]
    rows = [
        [sheet.cell_value(r, c) for c in range(sheet.ncols)]
        for r in range(2, sheet.nrows)
        if isinstance(sheet.cell_value(r, 0), float)
    ]
    return header, np.array(rows)


def last_plot(run):
    plots = sorted(pathlib.Path(run).glob("plt*"))
    if not plots:
        sys.exit(f"no plot files in {run}")
    return plots[-1]


def load(plot):
    ds = yt.load(str(plot))
    dims = ds.domain_dimensions
    cg = ds.covering_grid(0, ds.domain_left_edge, dims)
    j = dims[1] // 2
    out = {n: np.asarray(cg[("boxlib", n)])[:, j, :] for n in ("velocityx", "velocityz", "tke")}
    out["terrain_height"] = np.asarray(cg[("boxlib", "terrain_height")])[:, j, :]
    lo = ds.domain_left_edge.d
    dx = (ds.domain_right_edge.d - lo) / dims
    out["x"] = lo[0] + (np.arange(dims[0]) + 0.5) * dx[0]
    out["z"] = lo[2] + (np.arange(dims[2]) + 0.5) * dx[2]
    out["time"] = float(ds.current_time)
    return out


def column(f, x):
    i = int(np.argmin(np.abs(f["x"] - x)))
    zs = f["terrain_height"][i, 0]
    fluid = f["z"] > zs
    return i, zs, f["z"][fluid], {n: f[n][i, fluid] for n in ("velocityx", "velocityz", "tke")}


def interp(zm, zd, v):
    return np.interp(zd, zm, v, left=np.nan, right=np.nan)


def main():
    run, surface, variant = sys.argv[1], sys.argv[2], sys.argv[3]
    png = sys.argv[4] if len(sys.argv) > 4 else None
    f = load(last_plot(run))
    print(f"run {run}, t = {f['time']:.0f} s")

    if variant == "empty":
        _, d = read_sheet(FILES[(surface, "approach")])
        zd, ud = d[:, 0], d[:, 1]
        kd = 0.5 * (d[:, 2] ** 2 + d[:, 3] ** 2 + d[:, 4] ** 2)
        for x in (150.0, 400.0, 650.0):
            _, _, zm, v = column(f, x)
            um, km = interp(zm, zd, v["velocityx"]), interp(zm, zd, v["tke"])
            print(f"x = {x:.0f} m vs approach profile:")
            for z, a, b, c, e in zip(zd, ud, um, kd, km):
                print(f"  z {z:5.0f}  U data {a:.3f} model {b:.3f} ({b / a - 1:+.1%})  k data {c:.4f} model {e:.4f}")
            ok = ~np.isnan(um)
            print(f"  U rms rel error {np.sqrt(np.mean((um[ok] / ud[ok] - 1) ** 2)):.1%}")
        return

    _, d = read_sheet(FILES[(surface, "ridge")])
    stations = sorted(set(d[:, 0]))
    rms_u, rms_k = [], []
    for xs in stations:
        sel = d[:, 0] == xs
        zd, ud, wd = d[sel, 1], d[sel, 2], d[sel, 3]
        kd = 0.5 * (d[sel, 4] ** 2 + d[sel, 5] ** 2 + d[sel, 6] ** 2)
        _, zsurf, zm, v = column(f, XC + xs)
        um = interp(zm, zd, v["velocityx"])
        km = interp(zm, zd, v["tke"])
        ok = ~np.isnan(um)
        eu = np.sqrt(np.mean((um[ok] - ud[ok]) ** 2))
        ek = np.sqrt(np.mean((km[ok] - kd[ok]) ** 2))
        rms_u.append(eu)
        rms_k.append(ek)
        near = int(np.argmin(zd))
        print(
            f"x {xs:5.0f} mm (surface {zsurf:5.1f} m): U rms {eu:.3f} m/s, k rms {ek:.4f}; "
            f"lowest point z {zd[near]:.0f}: U data {ud[near]:+.3f} model {um[near]:+.3f}"
        )
    print(f"mean U rms over stations {np.mean(rms_u):.3f} m/s, mean k rms {np.mean(rms_k):.4f}")

    # Reversed flow in the first fluid cell above the surface, lee side
    first_u = []
    for i, x in enumerate(f["x"]):
        zs = f["terrain_height"][i, 0]
        k = int(np.argmax(f["z"] > zs))
        first_u.append(f["velocityx"][i, k])
    first_u = np.array(first_u)
    lee = (f["x"] > XC) & (f["x"] < XC + 800)
    rev = lee & (first_u < 0)
    if rev.any():
        xr = f["x"][rev] - XC
        print(f"model reversed flow at the first cell: x = {xr.min():.0f} to {xr.max():.0f} m from the crest "
              f"(reattachment ~{xr.max():.0f} m = {xr.max() / 40:.1f} h)")
    else:
        print("model: no reversed flow at the first cell")
    lowest = [(xs, d[(d[:, 0] == xs), 2][np.argmin(d[(d[:, 0] == xs), 1])]) for xs in stations]
    neg = [xs for xs, u in lowest if u < 0 and xs > 0]
    pos_after = [xs for xs, u in lowest if u >= 0 and neg and xs > max(neg)]
    if neg:
        print(f"data reversed flow at the lowest point: x = {min(neg):.0f} to {max(neg):.0f} mm; "
              f"first positive after: {min(pos_after) if pos_after else 'none'} mm")

    if png:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, len(stations), figsize=(2.0 * len(stations), 7), sharey="row")
        for a, xs in enumerate(stations):
            sel = d[:, 0] == xs
            zd = d[sel, 1]
            kd = 0.5 * (d[sel, 4] ** 2 + d[sel, 5] ** 2 + d[sel, 6] ** 2)
            _, zsurf, zm, v = column(f, XC + xs)
            axes[0, a].plot(d[sel, 2], zd, "ko", ms=3, label="tunnel")
            axes[0, a].plot(v["velocityx"], zm, "b-", label="KLAxell")
            axes[0, a].axvline(0, color="0.7", lw=0.5)
            axes[0, a].set_title(f"x = {xs:.0f} mm")
            axes[1, a].plot(kd, zd, "ko", ms=3)
            axes[1, a].plot(v["tke"], zm, "b-")
            for r in (0, 1):
                axes[r, a].set_ylim(0, 300)
        axes[0, 0].set_ylabel("z (m)  [U, m/s]")
        axes[1, 0].set_ylabel("z (m)  [k, m2/s2]")
        axes[0, 0].legend(fontsize=7)
        fig.suptitle(f"{surface} ridge, {pathlib.Path(run).name}, t = {f['time']:.0f} s")
        fig.tight_layout()
        fig.savefig(png, dpi=110)
        print(f"wrote {png}")


if __name__ == "__main__":
    main()
