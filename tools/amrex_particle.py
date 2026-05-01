# -*- coding: utf-8 -*-

"""\
AMReX Particle Data Processor
-----------------------------

This module provides a python class :class:`AmrexParticleFile` that can parse
the sampling data written out by Kynema-SGF using the native AMReX binary format.

"""

from pathlib import Path
import numpy as np
import pandas as pd
import yaml

class AmrexParticleFile:
    """AmrexParticleFile reader

    Loads AMReX Particle data in binary format for a single timestep
    """

    @classmethod
    def load(cls, time_index, label="sampling", root_dir="post_processing", suffix="/particles"):
        """
        Args:
            time_index (int): Time index to load
            label (str): Label used for this set of probes
            root_dir (path): Path to the post_processing directory
        """
        fpath = Path(root_dir) / ("%s%05d%s"%(label, time_index, suffix))
        return cls(fpath)

    def __init__(self, pdir):
        """
        Args:
            pdir (path): Directory path
        """
        self.pdir = Path(pdir)
        if not self.pdir.exists():
            raise ValueError(f'Path {self.pdir} does not exist.')

    def __call__(self):
        """Parse the header and load all binary files

        Returns:
            DataFrame: A pandas dataframe with data for all the particles
        """
        if not hasattr(self, "df"):
            self.parse_header()
            self.parse_info()
            self.load_binary_data()
        return self.df

    def parse_header(self):
        """Parse the header and load metadata"""
        with open(self.pdir / "Header", 'r') as fh:
            # Skip first line
            fh.readline()
            self.ndim = int(fh.readline().strip())
            self.num_reals = int(fh.readline().strip())
            self.real_var_names = [fh.readline().strip() for _ in
                                   range(self.num_reals)]
            self.num_ints = int(fh.readline().strip())
            self.int_var_names = [fh.readline().strip() for _ in
                                  range(self.num_ints)]
            # skip flag
            fh.readline()
            self.num_particles = int(fh.readline().strip())
            self.next_id = int(fh.readline().strip())
            self.finest_level = int(fh.readline().strip())

            self.num_grids = np.empty((self.finest_level+1,), dtype=int)
            self.grid_info = []
            for lev in range(self.finest_level+1):
                self.num_grids[lev] = int(fh.readline().strip())

            for lev in range(self.finest_level+1):
                ginfo = []
                for _ in range(self.num_grids[lev]):
                    ginfo.append(
                        [int(ix) for ix in fh.readline().strip().split()])
                self.grid_info.append(ginfo)

    def parse_info(self):
        """Parse the sampling info file"""
        with open(self.pdir.parent / "sampling_info.yaml", 'r') as fh:
            try:
                self.info = yaml.safe_load(fh)
            except yaml.YAMLError as exc:
                print(exc)

    def load_binary_data(self):
        """Read binary data into memory"""
        self.real_data = np.empty((self.num_particles,
                                   self.ndim+self.num_reals), dtype=float)
        self.int_data = np.empty((self.num_particles,
                                  self.num_ints), dtype=int)

        idata = self.int_data
        rdata = self.real_data
        nints = self.num_ints + 2
        nreals = self.num_reals + self.ndim
        for lev, ginfo in enumerate(self.grid_info):
            bfiles = {}
            for idx, npts, offset in ginfo:
                if npts < 1:
                    continue
                fname = self.bin_file_name(lev, idx)
                fstr = str(fname)
                if fstr not in bfiles:
                    bfiles[fstr] = open(fname, 'rb')
                fh = bfiles[fstr]
                fh.seek(offset)
                ivals = np.fromfile(fh, dtype=np.int32, count=nints*npts)
                rvals = np.fromfile(fh, dtype=float, count=nreals*npts)

                for i, ii in enumerate(range(0, nints * npts, nints)):
                    pidx = ivals[ii + 2]
                    idata[pidx, 0] = pidx
                    idata[pidx, 1] = ivals[ii + 3]
                    idata[pidx, 2] = ivals[ii + 4]

                    offset = i * nreals
                    for jj in range(nreals):
                        rdata[pidx, jj] = rvals[offset + jj]
            for fh in bfiles.values():
                fh.close()

        idict = {key: idata[:, i] for i, key in enumerate(self.int_var_names)}
        rvar_names = "xco yco zco".split() + self.real_var_names
        rdict = {key: rdata[:, i] for i, key in enumerate(rvar_names)}
        idict.update(rdict)
        self.df = pd.DataFrame(idict, index=idata[:,0])

    def bin_file_name(self, lev, idx):
        """Return the name of the binary file"""
        lvl_name = "Level_%d"%lev
        fname = "DATA_%05d"%idx
        return (self.pdir / lvl_name / fname)

    def __repr__(self):
        return "<%s: %s>"%(
            self.__class__.__name__, self.pdir.stem)
