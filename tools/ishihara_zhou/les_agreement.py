"""Near-surface agreement of the ridge runs with the digitized LES and the tunnel.

Usage: python3 les_agreement.py [dx in m, default 4]

For each surface and model (runs in cases/ridge2d/<surface>_<model>_<dx>m): U at
5 m (smooth) / 7 m (rough) above the surface at the LES stations, rms differences
against the LES and the tunnel (all stations and without the crest), and the lee
reversed-flow extent at that height.
"""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_ridge import FILES, XC, last_plot, load, read_sheet  # noqa: E402
from ref_les_liu2019 import H, REVERSED_XH, STATIONS_XH, U_OVER_UREF, UREF  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DX = sys.argv[1] if len(sys.argv) > 1 else "4"
RUNS = (("KLAxell", "klaxell"), ("limiter", "limiter"), ("all", "all"))


def u_at_height(f, hmin):
    out = []
    for i in range(len(f["x"])):
        zs = f["terrain_height"][i, 0]
        fluid = f["z"] > zs
        out.append(np.interp(zs + hmin, f["z"][fluid], f["velocityx"][i, fluid]))
    return np.array(out)


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


for surface, hmin in (("smooth", 5.0), ("rough", 7.0)):
    xs = H * np.array(STATIONS_XH)
    les = UREF[surface] * np.array(U_OVER_UREF[surface])
    _, d = read_sheet(FILES[(surface, "ridge")])
    tun = np.array([d[d[:, 0] == x, 2][np.argmin(d[d[:, 0] == x, 1])] for x in xs])
    nocrest = xs != 0.0
    print(f"== {surface} at {hmin:g} m; stations {xs.tolist()}")
    print(f"  LES    {np.round(les, 2).tolist()}")
    print(f"  tunnel {np.round(tun, 2).tolist()}")
    print(f"  LES vs tunnel rms {rms(les, tun):.2f} (no crest {rms(les[nocrest], tun[nocrest]):.2f})")
    print(f"  LES reversed flow {H * REVERSED_XH[surface][0]:.0f}-{H * REVERSED_XH[surface][1]:.0f} m")
    for label, pattern in RUNS:
        f = load(last_plot(HERE / "cases" / "ridge2d" / f"{surface}_{pattern}_{DX}m"))
        xm = f["x"] - XC
        um = u_at_height(f, hmin)
        mod = np.interp(xs, xm, um)
        lee = (xm > 0) & (xm < 600)
        rev = xm[lee & (um < 0)]
        extent = f"{rev.min():.0f}-{rev.max():.0f} m" if rev.size else "none"
        print(f"  {label:8s} {np.round(mod, 2).tolist()}")
        print(f"           vs LES rms {rms(mod, les):.2f} (no crest {rms(mod[nocrest], les[nocrest]):.2f}), "
              f"vs tunnel rms {rms(mod, tun):.2f} (no crest {rms(mod[nocrest], tun[nocrest]):.2f}), "
              f"reversed flow {extent}")
