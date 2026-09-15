"""Ishihara-Zhou 2D ridge: KLAxell and KLAxellSeparation against the tunnel.

Usage: python3 plot_ridge_overlay.py <smooth|rough> <output png> [dx in m, default 4]

Runs are read from cases/ridge2d/<surface>_<model>_<dx>m.

Top panel: U at the lowest measured height above the local surface (5 mm
smooth, 7 mm rough, i.e. 5 / 7 m at full scale) along x, which shows the
reversed-flow region and the reattachment; the model is sampled at the same
height above its terrain. Lower panels: U and k profiles at every station.
Runs that do not exist yet are skipped. Full scale = tunnel x1000.
"""

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_ridge import FILES, XC, column, last_plot, load, read_sheet  # noqa: E402
from ref_les_liu2019 import H, REVERSED_XH, STATIONS_XH, U_OVER_UREF, UREF  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
# Reference categorical slots 1-3 (validated all-pairs, light mode)
MODELS = (
    ("KLAxell", "{s}_klaxell_{dx}m", "#2a78d6"),
    ("KLAxellSeparation, limiter + 10 s gate", "{s}_limiter_{dx}m", "#eb6834"),
    ("KLAxellSeparation, all features", "{s}_all_{dx}m", "#1baf7a"),
)
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 8,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK2,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.titlecolor": INK,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
    }
)


def style(ax):
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main():
    surface, png = sys.argv[1], sys.argv[2]
    dx = sys.argv[3] if len(sys.argv) > 3 else "4"
    _, d = read_sheet(FILES[(surface, "ridge")])
    stations = sorted(set(d[:, 0]))
    hmin = 5.0 if surface == "smooth" else 7.0

    runs = []
    for label, pattern, color in MODELS:
        run = HERE / "cases" / "ridge2d" / pattern.format(s=surface, dx=dx)
        log = run / "run.log"
        # Only finished runs: a running case already has its t = 0 and
        # half-time plot files
        if log.exists() and "Time spent in Evolve" in log.read_text(errors="ignore"):
            runs.append((label, load(last_plot(run)), color))
        else:
            print(f"skipping {run.name}: not finished")

    n = len(stations)
    fig = plt.figure(figsize=(1.55 * n, 9.5))
    gs = fig.add_gridspec(3, n, height_ratios=(1.0, 1.3, 1.3), hspace=0.45, wspace=0.12)

    ax0 = fig.add_subplot(gs[0, :])
    xd, ud = [], []
    for xs in stations:
        sel = d[:, 0] == xs
        i = int(np.argmin(d[sel, 1]))
        xd.append(xs)
        ud.append(d[sel, 2][i])
    for label, f, color in runs:
        xm = f["x"] - XC
        keep = (xm >= min(stations) - 50) & (xm <= max(stations) + 50)
        um = []
        for i in np.where(keep)[0]:
            zs = f["terrain_height"][i, 0]
            fluid = f["z"] > zs
            um.append(np.interp(zs + hmin, f["z"][fluid], f["velocityx"][i, fluid]))
        ax0.plot(xm[keep], um, color=color, lw=2.0, label=label)
    # LES reference (Liu et al. 2019), digitized: station values and the
    # reversed-flow extent at the same height, drawn under the model lines
    x0r, x1r = H * REVERSED_XH[surface][0], H * REVERSED_XH[surface][1]
    ax0.plot([x0r, x1r], [0.0, 0.0], color=INK2, lw=4.0, solid_capstyle="butt", zorder=1.8,
             label="LES reversed flow (Liu et al. 2019)")
    ax0.plot(xd, ud, "o", ms=8, mfc=SURFACE, mec=INK, mew=1.2, label="Tunnel")
    # Smaller and on top, so a diamond stays visible inside a tunnel circle
    ax0.plot(H * np.array(STATIONS_XH), UREF[surface] * np.array(U_OVER_UREF[surface]), "D", ms=4.5,
             mfc=INK2, mec=INK2, mew=0.0, zorder=3, label="LES, Liu et al. 2019 (digitized)")
    ax0.axhline(0.0, color=AXIS, lw=1.0)
    ax0.axvspan(-100, 100, color=GRID, alpha=0.5, lw=0)
    # Tunnel reattachment bracket: from the last lee station with reversed flow
    # at the lowest point to the first station after it with forward flow
    lee_neg = [x for x, u in zip(xd, ud) if x > 0 and u < 0]
    if lee_neg:
        after = [x for x, u in zip(xd, ud) if x > max(lee_neg) and u >= 0]
        if after:
            x0, x1 = max(lee_neg), min(after)
            ax0.axvspan(x0, x1, facecolor="none", edgecolor=MUTED, hatch="///", lw=0.0)
            ax0.text((x0 + x1) / 2, 0.98, "tunnel\nreattachment", ha="center", va="top",
                     transform=ax0.get_xaxis_transform(), color=INK2, fontsize=8)
    ax0.text(0, ax0.get_ylim()[1] if False else 0.02, "ridge", ha="center", va="bottom",
             transform=ax0.get_xaxis_transform(), color=MUTED)
    ax0.set_xlabel("x from the crest (m, full scale; tunnel mm)")
    ax0.set_ylabel(f"U at {hmin:g} m above surface (m/s)")
    ax0.set_title(f"{surface.capitalize()} ridge: near-surface U along the ridge (reversed flow below 0)",
                  loc="left", fontsize=10)
    style(ax0)
    # Above the panel, right of the title, so it hides no data
    ax0.legend(frameon=False, ncol=3, fontsize=8, loc="lower right", bbox_to_anchor=(1.0, 1.01))

    # One k scale for every station: the largest k in the data or the runs
    kmax = float(np.max(0.5 * (d[:, 4] ** 2 + d[:, 5] ** 2 + d[:, 6] ** 2)))
    for _, f, _ in runs:
        for xs in stations:
            _, _, zm, v = column(f, XC + xs)
            kmax = max(kmax, float(np.max(v["tke"][zm <= 300])))

    for a, xs in enumerate(stations):
        sel = d[:, 0] == xs
        zd = d[sel, 1]
        kd = 0.5 * (d[sel, 4] ** 2 + d[sel, 5] ** 2 + d[sel, 6] ** 2)
        axu = fig.add_subplot(gs[1, a])
        axk = fig.add_subplot(gs[2, a])
        for label, f, color in runs:
            _, _, zm, v = column(f, XC + xs)
            top = zm <= 300
            axu.plot(v["velocityx"][top], zm[top], color=color, lw=1.6)
            axk.plot(v["tke"][top], zm[top], color=color, lw=1.6)
        axu.plot(d[sel, 2], zd, "o", ms=4, mfc=SURFACE, mec=INK, mew=0.9)
        axk.plot(kd, zd, "o", ms=4, mfc=SURFACE, mec=INK, mew=0.9)
        axu.axvline(0, color=AXIS, lw=0.8)
        axu.set_title(f"x = {xs:.0f} mm", fontsize=8)
        for ax in (axu, axk):
            ax.set_ylim(0, 300)
            style(ax)
            if a > 0:
                ax.set_yticklabels([])
        axu.set_xlim(-1.5, 6.5)
        axk.set_xlim(0, 1.08 * kmax)
    fig.axes[1].set_ylabel("z (m)\nU (m/s)")
    fig.axes[2].set_ylabel("z (m)\nk (m²/s²)")
    fig.savefig(png, dpi=120, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
