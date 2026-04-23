"""\
A tool to generate native boundary planes
-----------------------------------------

This script is for generating boundary plane data.

"""

import argparse
import pathlib
import numpy as np
import amrex.space3d as amr
import native_boundary_plane as nbp


class Gaussian:
    """Two dimensional (y, z) Gaussian functor."""

    def __init__(self, sigmay, sigmaz):
        self.sigmay = sigmay
        self.sigmaz = sigmaz
        self.yc = 500
        self.zc = 500

    def __call__(self, xg, yg, zg, time):
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape
        return (
            np.exp(
                -(
                    (yg - self.yc) ** 2 / (2 * self.sigmay**2)
                    + (zg - self.zc) ** 2 / (2 * self.sigmaz**2)
                )
            )
            * time
        )


class Constant:
    """Constant functor."""

    def __init__(self, constant):
        self.constant = constant

    def __call__(self, xg, yg, zg, time):
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape
        return self.constant * np.ones(xg.shape)


def main():
    parser = argparse.ArgumentParser(
        description="A tool to generate native boundary planes"
    )
    parser.add_argument(
        "-i",
        "--iname",
        help="Kynema-SGF input file for boundary data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Force overwrite of existing files",
    )
    args = parser.parse_args()

    amr.initialize([])

    bdir = pathlib.Path("bndry_file")
    if bdir.exists() and not args.overwrite:
        raise Exception(
            f"{bdir} exists in this directory. Skipping. Use -o to overwrite."
        )

    field = "velocity"
    ori = 0
    ncomp = 3

    nsteps = 100
    steps = [x for x in range(nsteps + 1)]
    times = np.linspace(0, 5, nsteps + 1)
    functors = {
        "velocityx": Gaussian(100, 100),
        "velocityy": Gaussian(100, 100),
        "velocityz": Gaussian(100, 100),
    }

    for step, time in zip(steps, times):
        odir = bdir / f"bndry_output{step:05d}"
        print(f"Generating {odir}")
        bp = nbp.NativeBoundaryPlane(field, ncomp, ori, odir)
        bp.define_from_file(args.iname, step, time)
        bp.evaluate(functors)
        bp.write()

    np.savetxt(bdir / "time.dat", np.c_[steps, times], fmt="%.17g")


if __name__ == "__main__":
    main()
