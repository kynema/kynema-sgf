import argparse
import glob
import pathlib
import shutil
import amrex.space3d as amr
from amrex_plotfile import AmrexPlotFile
import amrex_utils as au
import native_boundary_plane as nbp
import numpy as np


def parse_fields_from_header_name(hname):
    return hname.split("_")[-1]


def parse_orientations_from_header_name(hname):
    return int(hname.split("_")[-2])


def main():
    parser = argparse.ArgumentParser(description="A tool to refine boundary planes")
    parser.add_argument(
        "-f",
        "--fdir",
        help="Kynema-SGF directory with boundary data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-r",
        "--refinement_ratio",
        help="Refinement ratio",
        default=2,
        type=int,
    )
    args = parser.parse_args()

    amr.initialize([])

    for fname in sorted(glob.glob(f"{args.fdir}/bndry_output" + "*")):
        print(f"Refining {fname}")
        headers = [pathlib.Path(x) for x in glob.glob(f"{fname}/Header_*")]
        fpath = pathlib.Path(fname)
        wname = pathlib.Path(fpath.parent.name) / pathlib.Path(fpath.name)

        for hname in headers:
            field = parse_fields_from_header_name(hname.name)
            ori = parse_orientations_from_header_name(hname.name)
            ncomp = au.ncomp_from_field(field)
            bp = nbp.NativeBoundaryPlane(field, ncomp, ori, wname)

            plt = AmrexPlotFile(hname)
            mfs = plt()

            bp.define_from_plotfile(hname)
            bp.refine(args.refinement_ratio)
            au.interpolate(plt, bp.plt)
            bp.write()

    # Copy the time file
    tname = "time.dat"
    spath = pathlib.Path(args.fdir)
    src_time_file = spath / tname
    dst_time_file = pathlib.Path(spath.name) / tname
    shutil.copyfile(src_time_file, dst_time_file)


if __name__ == "__main__":
    main()
