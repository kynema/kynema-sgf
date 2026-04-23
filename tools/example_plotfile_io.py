"""Example script for reading and writing a Kynema-SGF plot file."""

import argparse
import pathlib
import amrex.space3d as amr
from amrex_plotfile import AmrexPlotFile
import numpy as np


class Gaussian:
    """Three dimensional Gaussian functor."""

    def __init__(self, sigmax, sigmay, sigmaz):
        self.sigmax = sigmax
        self.sigmay = sigmay
        self.sigmaz = sigmaz
        self.xc = 500
        self.yc = 500
        self.zc = 500

    def __call__(self, xg, yg, zg, time):
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape
        return np.exp(
            -(
                (xg - self.xc) ** 2 / (2 * self.sigmax**2)
                + (yg - self.yc) ** 2 / (2 * self.sigmay**2)
                + (zg - self.zc) ** 2 / (2 * self.sigmaz**2)
            )
        )


def main():
    parser = argparse.ArgumentParser(
        description="An example tool to manipulate plot file data"
    )
    parser.add_argument(
        "-f",
        "--fname",
        help="Kynema-SGF plot file",
        required=True,
        type=str,
    )
    args = parser.parse_args()

    amr.initialize([])

    hname = pathlib.Path(args.fname) / "Header"
    wname = pathlib.Path(hname.parent.name)
    plt = AmrexPlotFile(hname)
    mfs = plt()
    plt.evaluate(
        {
            "velocityx": Gaussian(100, 200, 300),
            "velocityy": Gaussian(200, 300, 100),
            "velocityz": Gaussian(300, 100, 200),
        }
    )
    plt.write(wname)


if __name__ == "__main__":
    main()
