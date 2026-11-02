"""Plot how the RANSConvergence monitor behaved on a KLAxell run.

Reads two outputs from a run directory:

  post_processing/<monitor label>NNNNN.txt   the monitor's diagnostics file
  post_processing/<sampler label>NNNNN.txt   ASCII probe samples of velocity
                                             and tke at the same points

and writes a three-panel figure: the horizontal speed at each point, the
spread-over-tolerance history the monitor tests together with its hold periods
and exponential fit, and the sample-to-sample change that a simpler criterion
would have used.

The probe samples must come from a Sampling post-processor with
output_format = ascii, a single ProbeSampler reading the same point file as
the monitor, and fields = velocity tke. ASCII sample files carry a step index
but no time, so the fixed timestep of the run is required.

Example:

  python3 plot_rans_convergence.py run_dir --dt 5.0 -o convergence.png
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_monitor(path):
    """Read the monitor's diagnostics file into a dict of columns."""
    data = np.atleast_2d(np.loadtxt(path, skiprows=1))
    # Rows written while the window is still filling are dropped: a window
    # holding a couple of samples has a small spread whatever the flow is
    # doing, so those rows say nothing about convergence
    data = data[data[:, 2] > 0.5]
    return {
        "time": data[:, 0],
        "vel_spread": data[:, 6],
        "vel_tol": data[:, 7],
        "tke_spread": data[:, 9],
        "tke_tol": data[:, 10],
        "hold_elapsed": data[:, 11],
    }


def read_probes(post_dir, label, dt):
    """Read the ASCII probe files into times and per-point arrays.

    Each row of a probe file is x y z, some particle identifiers, and then the
    sampled fields in the order they were requested: u, v, w, tke. Only the
    trailing four columns are relied on, since the identifier block is an
    implementation detail of the particle container.
    """
    pattern = re.compile(re.escape(label) + r"(\d+)\.txt$")
    times, rows = [], []
    for fname in glob.glob(os.path.join(post_dir, label + "*.txt")):
        match = pattern.search(os.path.basename(fname))
        if match is None:
            continue  # e.g. the _info.txt files
        vals = []
        with open(fname) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 7:
                    continue  # header lines
                vals.append([float(p) for p in parts[-4:]])
        if vals:
            times.append(int(match.group(1)) * dt)
            rows.append(vals)
    if not rows:
        sys.exit(f"No ASCII probe files named {label}NNNNN.txt in {post_dir}")

    order = np.argsort(times)
    # Shape: (nsamples, npoints, 4) for u, v, w, tke
    return np.array(times)[order], np.array(rows)[order]


def fit_decay(times, spreads, threshold=1.0):
    """Least squares of ln(s) against t, the same fit the monitor performs."""
    mask = spreads > 0.0
    if mask.sum() < 5:
        return None
    logs = np.log(spreads[mask])
    slope, intercept = np.polyfit(times[mask], logs, 1)
    rate, amplitude = -slope, np.exp(intercept)
    if rate <= 0.0:
        return None
    resid = logs - (intercept + slope * times[mask])
    ss_tot = np.sum((logs - np.mean(logs)) ** 2)
    rsq = 1.0 - (np.sum(resid**2) / ss_tot) if ss_tot > 0 else 0.0
    return {
        "rate": rate,
        "amplitude": amplitude,
        "rsq": rsq,
        "t_cross": np.log(amplitude / threshold) / rate,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Plot RANSConvergence diagnostics from a run directory"
    )
    parser.add_argument("run_dir", help="directory holding post_processing/")
    parser.add_argument(
        "--dt", type=float, required=True, help="fixed timestep of the run [s]"
    )
    parser.add_argument(
        "--monitor-label", default="convergence", help="RANSConvergence label"
    )
    parser.add_argument(
        "--sampler-label", default="sampling", help="Sampling label"
    )
    parser.add_argument(
        "-o", "--output", default="rans_convergence_monitor.png"
    )
    args = parser.parse_args()

    post_dir = os.path.join(args.run_dir, "post_processing")
    monitor_files = sorted(
        glob.glob(os.path.join(post_dir, args.monitor_label + "[0-9]*.txt"))
    )
    if not monitor_files:
        sys.exit(f"No {args.monitor_label}NNNNN.txt found in {post_dir}")
    mon = read_monitor(monitor_files[0])
    if mon["time"].size == 0:
        sys.exit("The monitor never filled its window; nothing to plot")
    ptimes, pvals = read_probes(post_dir, args.sampler_label, args.dt)

    hours = 3600.0
    speed = np.hypot(pvals[:, :, 0], pvals[:, :, 1])
    npts = speed.shape[1]

    # What the criterion actually tests: spread divided by that point's own
    # tolerance, so convergence is the single value one for both quantities
    vel_ratio = mon["vel_spread"] / mon["vel_tol"]
    tke_ratio = mon["tke_spread"] / mon["tke_tol"]

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 10.5), sharex=True)

    # Panel 1: the signal the monitor watches
    ax = axes[0]
    for i in range(npts):
        ax.plot(ptimes / hours, speed[:, i], lw=1.2, label=f"point {i}")
    ax.set_ylabel("horizontal speed  $\\sqrt{u^2+v^2}$  [m/s]")
    ax.set_title("Horizontal speed at the monitor points", fontsize=11)
    ax.legend(fontsize=9, loc="lower right", ncol=min(npts, 6))
    ax.grid(alpha=0.3)

    # Panel 2: the envelope test, the hold, and the extrapolation
    ax = axes[1]
    ax.semilogy(mon["time"] / hours, vel_ratio, lw=1.4, label="speed envelope")
    ax.semilogy(mon["time"] / hours, tke_ratio, lw=1.4, label="tke envelope")
    ax.axhline(
        1.0, color="k", ls="--", lw=1.2, label="tolerance (converged below)"
    )
    fit = fit_decay(mon["time"], vel_ratio)
    if fit is not None:
        tspan = np.linspace(
            mon["time"][0], max(fit["t_cross"], mon["time"][-1]), 200
        )
        ax.semilogy(
            tspan / hours,
            fit["amplitude"] * np.exp(-fit["rate"] * tspan),
            ":",
            color="C0",
            lw=1.6,
            label=(
                f"exponential fit, $R^2$={fit['rsq']:.2f}, "
                f"crossing at {fit['t_cross'] / hours:.1f} h"
            ),
        )
    holding = mon["hold_elapsed"] > 0.0
    if holding.any():
        ax.fill_between(
            mon["time"] / hours,
            *ax.get_ylim(),
            where=holding,
            color="C2",
            alpha=0.18,
            label="all points within tolerance, holding",
        )
    ax.set_ylabel("spread / tolerance")
    ax.set_title(
        "What the monitor tests: peak-to-trough spread over a trailing window",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3, which="both")

    # Panel 3: the sample-to-sample change a simpler criterion would use
    ax = axes[2]
    worst = int(np.argmax(np.ptp(speed, axis=0)))
    naive = np.abs(np.diff(speed[:, worst]))
    ax.semilogy(
        ptimes[1:] / hours,
        np.maximum(naive, 1e-12),
        lw=1.0,
        color="C3",
        label=f"|change between samples|, point {worst}",
    )
    ax.semilogy(
        mon["time"] / hours,
        mon["vel_spread"],
        lw=1.6,
        color="C0",
        label="window spread, worst point",
    )
    ax.axhline(
        mon["vel_tol"][0], color="k", ls="--", lw=1.2, label="speed tolerance"
    )
    ax.set_ylabel("speed change [m/s]")
    ax.set_xlabel("time [hours]")
    ax.set_title(
        "Sample-to-sample change against the window spread", fontsize=11
    )
    ax.legend(fontsize=9, loc="lower left")
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(args.output, dpi=160)
    print("wrote", args.output)


if __name__ == "__main__":
    main()
