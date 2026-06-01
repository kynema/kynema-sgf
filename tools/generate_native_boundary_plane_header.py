"""\
A tool to generate native boundary planes header files
------------------------------------------------------

This script is only for generating the header files for boundary plane
data that does not already have them. Most likely they were generated
with Kynema-SGF version 3.2.0 and prior.

"""

import argparse
import glob
import pathlib
import itertools
import amrex.space3d as amr
import amrex_utils as au
import native_boundary_plane as nbp
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="A tool to generate native boundary planes header files"
    )
    parser.add_argument(
        "-f",
        "--fdir",
        help="Kynema-SGF directory with boundary data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-i",
        "--iname",
        help="Kynema-SGF input file that generated the boundary data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Force overwrite of existing header files",
    )
    args = parser.parse_args()

    amr.initialize([])
    spacedim = amr.Config.spacedim

    spath = pathlib.Path(args.fdir)
    tname = "time.dat"
    time_file = spath / tname
    times = pd.read_csv(time_file, sep="\\s+", names=["step", "time"], header=None)
    pfx = "bndry_output"
    lvl_pfx = "Level_"

    for fname in sorted(glob.glob(f"{args.fdir}/{pfx}" + "*")):
        print(f"Generating Header files for data in {fname}")
        fpath = pathlib.Path(fname)
        step = int(fpath.name.replace(pfx, ""))
        time = (times.time[times.step == step]).values[0]

        finest_level = len(sorted(glob.glob(f"{fname}/{lvl_pfx}" + "*"))) - 1
        nlevels = finest_level + 1
        mf_h_names = sorted(glob.glob(f"{fname}/{lvl_pfx}0/" + "*" + "_H"))

        fields = []
        oris = []
        for mf_h_name in mf_h_names:
            fields.append(pathlib.Path(mf_h_name).name.split("_")[0])
            oris.append(int(pathlib.Path(mf_h_name).name.split("_")[1]))
        fields = list(set(fields))
        oris = list(set(oris))

        for field, ori in itertools.product(fields, oris):
            hname = fpath / f"Header_{ori}_{field}"
            if hname.exists() and not args.overwrite:
                print(f"{hname} exists already. Skipping.")
                continue

            mfs = []
            for ilev in range(nlevels):
                mf_h_name = fpath / f"{lvl_pfx}{ilev}" / f"{field}_{ori}_H"
                mfs.append(amr.VisMF.Read(str(mf_h_name).replace("_H", "")))
            ncomp = mfs[0].num_comp

            bp = nbp.NativeBoundaryPlane(field, ncomp, ori, fpath)
            bp.define_from_mfs(args.iname, step, time, mfs)
            bp.plt.write_header(hname)


if __name__ == "__main__":
    main()
