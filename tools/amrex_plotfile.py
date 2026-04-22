# -*- coding: utf-8 -*-

"""\
AMReX Plot File Processor
-------------------------

This module provides a python class :class:`AmrexPlotFile` that can parse
the plot files written out by Kynema-SGF.

"""

import os
import pathlib
import numpy as np
import amrex.space3d as amr


class AmrexPlotFile:
    """AmrexPlotFile reader

    Loads AMReX plot file data for a single timestep
    """

    def __init__(self, hname):
        """
        Args:
            hname (path like): header file name
        """

        # Disallow OMP/MPI for now
        if amr.Config.have_mpi:
            nprocs = amr.ParallelDescriptor.NProcs()
            assert nprocs == 1
        if amr.Config.have_omp:
            nthreads = int(os.environ["OMP_NUM_THREADS"])
            assert nthreads == 1

        self.fname = pathlib.Path(hname)

    def __call__(self):
        """Parse the header file"""
        if not hasattr(self, "mfs"):
            self.parse_header()
            self.load_data()
        return self.mfs

    def parse_header(self):
        """Parse the header"""
        assert self.fname.exists()
        with open(self.fname, "r") as f:
            self.version = f.readline()
            self.ncomp = int(f.readline())
            self.names = []
            for _ in range(self.ncomp):
                self.names.append(f.readline().rstrip())
            self.spacedim = int(f.readline())
            assert self.spacedim == 3
            self.time = float(f.readline())
            self.finest_level = int(f.readline())
            self.nlevels = self.finest_level + 1

            self.prob_lo = [float(x) for x in f.readline().split()]
            self.prob_hi = [float(x) for x in f.readline().split()]
            self.ref_ratio = [int(x) for x in f.readline().split()]

            line = [
                int(x)
                for x in f.readline()
                .rstrip()
                .replace("(", " ")
                .replace(")", " ")
                .replace(",", " ")
                .split()
            ]
            self.prob_domain = []
            for ilev in range(self.nlevels):
                offset = ilev * 9
                small_end = amr.IntVect(line[0 + offset : 3 + offset])
                big_end = amr.IntVect(line[3 + offset : 6 + offset])
                t = amr.IntVect(line[6 + offset : 9 + offset])
                self.prob_domain.append(amr.Box(small_end, big_end))

            self.level_steps = [int(x) for x in f.readline().rstrip().split()]

            self.cell_sizes = []
            for ilev in range(self.nlevels):
                line = [float(x) for x in f.readline().rstrip().split()]
                self.cell_sizes.append(line)

            self.coord_system = int(f.readline())
            self.bwidth = int(f.readline())

            self.mf_names = []
            self.level_steps = []
            self.ngrids = []
            self.glohis = [[] for _ in range(self.nlevels)]
            for ilev in range(self.nlevels):
                line = f.readline().split()
                lev = int(line[0])
                self.ngrids.append(int(line[1]))
                gtime = float(line[2])
                self.level_steps.append(int(f.readline()))
                nglohis = []
                for igrid in range(self.ngrids[ilev]):
                    glos = []
                    ghis = []
                    for idim in range(self.spacedim):
                        glo, ghi = [float(x) for x in f.readline().split()]
                        glos.append(glo)
                        ghis.append(ghi)
                    nglohis.append([glos, ghis])
                self.glohis[ilev] = nglohis
                assert len(self.glohis[ilev]) == self.ngrids[ilev]
                self.mf_names.append(f.readline().rstrip())

    def load_data(self):
        """Read data into memory"""
        self.mfs = []
        for mf_name in self.mf_names:
            self.mfs.append(amr.VisMF.Read(str(self.fname.parent / mf_name)))

    def write(self, fname, overwrite=False):
        """
        Args:
            fname (path): folder name for writing
        """
        hname = pathlib.Path(fname) / self.fname.name
        if not overwrite:
            assert hname != self.fname
        self.write_header(hname)
        self.write_data(fname)

    def write_header(self, fname):
        """Write header file

        Args:
            fname (path): header file name for writing
        """
        fname.parent.mkdir(exist_ok=True, parents=True)
        with open(fname, "w") as f:
            f.write(self.version)
            f.write(str(self.ncomp) + "\n")
            f.write("\n".join(self.names))
            f.write("\n")
            f.write(f"{self.spacedim}\n")
            f.write(f"{self.time:.17g}\n")
            f.write(f"{self.finest_level}\n")
            f.write(" ".join([f"{x:.17g}" for x in self.prob_lo]))
            f.write(" \n")
            f.write(" ".join([f"{x:.17g}" for x in self.prob_hi]))
            f.write(" \n")
            f.write(" ".join([f"{x:.17g}" for x in self.ref_ratio]))
            f.write(" \n")
            line = "("
            for ilev in range(self.nlevels):
                line += "("
                for idim in range(self.spacedim):
                    end = "," if idim != self.spacedim - 1 else ")"
                    line += str(self.prob_domain[ilev].small_end[idim]) + end
                line += " ("
                for idim in range(self.spacedim):
                    end = "," if idim != self.spacedim - 1 else ")"
                    line += str(self.prob_domain[ilev].big_end[idim]) + end
                line += " ("
                for idim in range(self.spacedim):
                    end = "," if idim != self.spacedim - 1 else ")"
                    line += str(0) + end
                end = " (" if ilev != self.finest_level else ""
                line += ")" + end
            f.write(f"{line} \n")

            f.write(" ".join([str(x) for x in self.level_steps]))
            f.write(" \n")

            for ilev in range(self.nlevels):
                f.write(" ".join([str(x) for x in self.cell_sizes[ilev]]))
                f.write(" \n")

            f.write(f"{self.coord_system}\n")
            f.write(f"{self.bwidth}\n")

            for ilev in range(self.nlevels):
                f.write(f"{ilev} {self.ngrids[ilev]} {self.time:.17g}\n")
                f.write(f"{self.level_steps[ilev]}\n")
                for igrid in range(self.ngrids[ilev]):
                    for idim in range(self.spacedim):
                        glo = self.glohis[ilev][igrid][0][idim]
                        ghi = self.glohis[ilev][igrid][1][idim]
                        f.write(f"{glo:.17g} {ghi:.17g}\n")
                f.write(f"{self.mf_names[ilev]}\n")

    def write_data(self, fname):
        """Write data file

        Args:
            fname (path): header file name for writing
        """
        assert fname.exists()
        assert hasattr(self, "mfs")
        for mf_name, mf in zip(self.mf_names, self.mfs):
            fdir = pathlib.Path(fname / mf_name).parent
            fdir.mkdir(exist_ok=True, parents=True)
            amr.VisMF.Write(mf, str(fname / mf_name))

    def define(
        self,
        mfs,
        vnames,
        mf_names,
        step,
        time,
        prob_lo,
        prob_hi,
        ref_ratio,
        prob_domain,
        glohis,
    ):
        assert len(mfs) == len(mf_names)
        assert len(mfs) == len(prob_domain)
        assert len(mfs) - 1 <= len(ref_ratio)
        self.version = "HyperCLaw-V1.1\n"
        self.mfs = mfs
        self.ncomp = self.mfs[0].num_comp
        self.names = vnames
        self.spacedim = amr.Config.spacedim
        self.time = time
        self.finest_level = len(self.mfs) - 1
        self.nlevels = self.finest_level + 1
        self.prob_lo = prob_lo
        self.prob_hi = prob_hi
        self.ref_ratio = ref_ratio
        self.prob_domain = prob_domain
        self.level_steps = [step] * self.nlevels
        self.cell_sizes = [
            [
                (self.prob_hi[i] - self.prob_lo[i]) / (prob_domain[0].size[i])
                for i in range(self.spacedim)
            ]
        ]
        for ilev in range(1, self.nlevels):
            self.cell_sizes.append(
                [x / (self.ref_ratio[ilev - 1] ** ilev) for x in self.cell_sizes[0]]
            )
        self.coord_system = 0
        self.bwidth = 0
        self.ngrids = [mf.box_array().size for mf in self.mfs]
        self.glohis = glohis
        self.mf_names = mf_names

    def evaluate(self, functors):
        """Evaluate dictionary of functors on the plot file MultiFabs

        Args:
            functors (dict): dictionary of functors to be evaluated of format {field_name: functor}
        """
        for name, functor in functors.items():
            icomp = self.names.index(name)
            for ilev in range(self.nlevels):
                for igrid in range(self.ngrids[ilev]):
                    data = self.mfs[ilev].to_xp()[igrid]
                    xg, yg, zg = self.coordinates(ilev, igrid)
                    data[:, :, :, icomp] = functor(xg, yg, zg, self.time)

    def coordinates(self, ilev, igrid):
        los = self.glohis[ilev][igrid][0]
        his = self.glohis[ilev][igrid][1]
        nxs = self.mfs[ilev].box_array()[igrid].size
        dxs = self.cell_sizes[ilev]
        xs = []
        for lo, hi, nx, dx in zip(los, his, nxs, dxs):
            xs.append(np.linspace(lo + 0.5 * dx, hi - 0.5 * dx, nx))
        xg, yg, zg = np.meshgrid(xs[0], xs[1], xs[2], indexing="ij")
        return xg, yg, zg

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.fname.stem)
