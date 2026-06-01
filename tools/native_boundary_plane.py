# -*- coding: utf-8 -*-

"""\
Kynema-SGF Native Boundary Plane
------------------------------

This module provides a python class :class:`NativeBoundaryPlane` that can be used to define, manipulate, and output native boundary planes for Kynema-SGF

"""
from amrex_plotfile import AmrexPlotFile
import amrex.space3d as amr
import amrex_utils as au


def parse_input_file(iname):
    spacedim = amr.Config.spacedim
    pp = amr.ParmParse("")
    pp.addfile(iname)
    pp_geom = amr.ParmParse("geometry")
    prob_lo = [
        pp_geom.get_real("prob_lo", 0),
        pp_geom.get_real("prob_lo", 1),
        pp_geom.get_real("prob_lo", 2),
    ]
    prob_hi = [
        pp_geom.get_real("prob_hi", 0),
        pp_geom.get_real("prob_hi", 1),
        pp_geom.get_real("prob_hi", 2),
    ]
    pp_amr = amr.ParmParse("amr")
    n_cell = [
        pp_amr.get_int("n_cell", 0),
        pp_amr.get_int("n_cell", 1),
        pp_amr.get_int("n_cell", 2),
    ]
    finest_level = pp_amr.get_int("max_level")
    nlevels = finest_level + 1
    ref_ratio = []
    with open(iname, "r") as f:
        for line in f:
            if "amr.ref_ratio" in line:
                line = line.split()
                ref_ratio = [int(x) for x in line[2:]]
    if not ref_ratio:
        ref_ratio = [2] * (nlevels - 1)
    cell_sizes = []
    for ilev in range(nlevels):
        cs = [(prob_hi[i] - prob_lo[i]) / (n_cell[i]) for i in range(spacedim)]
        if ilev > 0:
            cs = [x / (ref_ratio[ilev - 1] ** ilev) for x in cs]
        cell_sizes.append(cs)
    return prob_lo, prob_hi, n_cell, nlevels, ref_ratio, cell_sizes


def adjust_plane(ori, prob_lo, prob_hi, n_cell):
    normal = au.normal_from_ori(ori)
    dx_normal = (prob_hi[normal] - prob_lo[normal]) / (n_cell[normal])
    plo_normal = prob_lo[normal]
    phi_normal = prob_hi[normal]
    n_cell[normal] = 2
    prob_lo[normal] = (
        plo_normal - dx_normal if au.is_low(ori) else phi_normal - dx_normal
    )
    prob_hi[normal] = (
        plo_normal + dx_normal if au.is_low(ori) else phi_normal + dx_normal
    )

    return prob_lo, prob_hi, n_cell


class NativeBoundaryPlane:
    def __init__(self, field, ncomp, ori, odir):
        self.field = field
        self.ncomp = ncomp
        self.ori = ori
        self.odir = odir
        self.hname = self.odir / f"Header_{ori}_{field}"

    def define_from_file(self, iname, step, time):
        spacedim = amr.Config.spacedim
        prob_lo, prob_hi, n_cell, nlevels, ref_ratio, cell_sizes = parse_input_file(
            iname
        )

        normal = au.normal_from_ori(self.ori)
        perp = au.perpendicular_from_ori(self.ori)
        domain_lo = amr.IntVect(0)
        domain_hi = amr.IntVect(n_cell[0] - 1, n_cell[1] - 1, n_cell[2] - 1)
        domain_lo[normal] = -1 if au.is_low(self.ori) else domain_hi[normal]
        domain_hi[normal] = 0 if au.is_low(self.ori) else domain_hi[normal] + 1
        domain = amr.Box(domain_lo, domain_hi)

        prob_lo, prob_hi, n_cell = adjust_plane(self.ori, prob_lo, prob_hi, n_cell)

        prob_domain = [domain]
        bas = [amr.BoxArray(domain)]
        for ilev in range(1, nlevels):
            small_end = prob_domain[ilev - 1].small_end
            big_end = prob_domain[ilev - 1].big_end
            for p in perp:
                small_end[p] *= ref_ratio[ilev - 1]
                big_end[p] = (big_end[p] + 1) * ref_ratio[ilev - 1] - 1
            if not au.is_low(self.ori):
                small_end[normal] = small_end[normal] * ref_ratio[ilev - 1] + 1
                big_end[normal] = small_end[normal] + 1

            prob_domain.append(amr.Box(small_end, big_end))
            bas.append(amr.BoxArray(prob_domain[ilev]))

        mfs = []
        n_grow = 0
        for ilev in range(nlevels):
            dm = amr.DistributionMapping(bas[ilev])
            mf = amr.MultiFab(bas[ilev], dm, self.ncomp, n_grow)
            mf.set_val(0.0)
            mfs.append(mf)

        glohis = [[] for _ in range(nlevels)]
        for ilev in range(nlevels):
            nglohis = []
            for igrid in range(mfs[ilev].box_array().size):
                glos = []
                ghis = []
                for idim in range(spacedim):
                    if idim == normal:
                        offset = (prob_lo[idim] + prob_hi[idim]) / 2
                        glo = offset - cell_sizes[ilev][idim]
                        ghi = offset + cell_sizes[ilev][idim]
                    else:
                        glo = prob_lo[idim]
                        ghi = prob_hi[idim]
                    glos.append(glo)
                    ghis.append(ghi)
                nglohis.append([glos, ghis])
            glohis[ilev] = nglohis

        self.plt = AmrexPlotFile(self.hname)
        lvl_pfx = "Level_"
        vnames = au.variable_names(self.field, self.ncomp)
        mf_names = [f"{lvl_pfx}{x}/{self.field}_{self.ori}" for x in range(nlevels)]
        self.plt.define(
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
        )

    def define_from_plotfile(self, hname):
        self.plt = AmrexPlotFile(hname)
        mfs = self.plt()
        assert mfs[0].num_comp == self.ncomp

    def define_from_mfs(self, iname, step, time, mfs):
        nlevels = len(mfs)
        spacedim = amr.Config.spacedim
        prob_lo, prob_hi, n_cell, nlevels, ref_ratio, cell_sizes = parse_input_file(
            iname
        )
        prob_lo, prob_hi, n_cell = adjust_plane(self.ori, prob_lo, prob_hi, n_cell)
        assert nlevels == len(mfs)

        bas = [mf.box_array() for mf in mfs]
        prob_domain = []
        for ilev in range(nlevels):
            prob_domain.append(amr.Box(bas[ilev][0].small_end, bas[ilev][0].big_end))

        ngrids = [ba.size for ba in bas]
        for ilev in range(nlevels):
            assert ngrids[ilev] == 1

        normal = au.normal_from_ori(self.ori)
        glohis = [[] for _ in range(nlevels)]
        for ilev in range(nlevels):
            nglohis = []
            for igrid in range(mfs[ilev].box_array().size):
                glos = []
                ghis = []
                for idim in range(spacedim):
                    if idim == normal:
                        offset = (prob_lo[idim] + prob_hi[idim]) / 2
                        glo = offset - cell_sizes[ilev][idim]
                        ghi = offset + cell_sizes[ilev][idim]
                    else:
                        glo = prob_lo[idim]
                        ghi = prob_hi[idim]
                    glos.append(glo)
                    ghis.append(ghi)
                nglohis.append([glos, ghis])
            glohis[ilev] = nglohis

        self.plt = AmrexPlotFile(self.hname)
        lvl_pfx = "Level_"
        mf_names = [f"{lvl_pfx}{x}/{self.field}_{self.ori}" for x in range(nlevels)]
        vnames = au.variable_names(self.field, self.ncomp)
        self.plt.define(
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
        )

    def evaluate(self, functors):
        self.plt.evaluate(functors)

    def write(self):
        overwrite = True
        self.plt.write(self.odir, overwrite)

    def refine(self, refinement_ratio):
        """Update data in a boundary plane for integer refinement."""
        normal = au.normal_from_ori(self.ori)
        perp = au.perpendicular_from_ori(self.ori)

        for ilev in range(self.plt.nlevels):
            self.plt.cell_sizes[ilev] = [
                x / refinement_ratio for x in self.plt.cell_sizes[ilev]
            ]

        plo = self.plt.prob_lo[normal]
        phi = self.plt.prob_hi[normal]
        offset = (phi + plo) / 2
        self.plt.prob_lo[normal] = offset - self.plt.cell_sizes[0][normal]
        self.plt.prob_hi[normal] = offset + self.plt.cell_sizes[0][normal]

        for ilev in range(self.plt.nlevels):
            small_endi = self.plt.prob_domain[ilev].small_end
            big_endi = self.plt.prob_domain[ilev].big_end
            for p in perp:
                small_endi[p] *= refinement_ratio
                big_endi[p] = (big_endi[p] + 1) * refinement_ratio - 1
            self.plt.prob_domain[ilev] = amr.Box(small_endi, big_endi)

            dx = self.plt.cell_sizes[ilev][normal]
            for igrid in range(self.plt.ngrids[ilev]):
                glo = self.plt.glohis[ilev][igrid][0][normal]
                ghi = self.plt.glohis[ilev][igrid][1][normal]
                offset = (ghi + glo) / 2
                self.plt.glohis[ilev][igrid][0][normal] = offset - dx
                self.plt.glohis[ilev][igrid][1][normal] = offset + dx

            bac = self.plt.mfs[ilev].box_array()
            rr = amr.IntVect(refinement_ratio)
            rr[normal] = 1
            baf = bac.refine(rr)

            dm = self.plt.mfs[ilev].dm()
            self.plt.mfs[ilev] = amr.MultiFab(baf, dm, self.plt.ncomp, 0)
            self.plt.mfs[ilev].set_val(0.0)
