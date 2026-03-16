#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

std::string getFileRoot(const std::string& infile)
{
    const auto& tokens = Tokenize(infile, std::string("/"));
    return tokens[tokens.size() - 1];
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::ParmParse pp;

        // Open plotfile headers and create an amrData objects pointing into it
        // Reference file: (should be finer than the finest level of the
        // multilevel file)
        std::string plotFileName1;
        pp.get("referencefile", plotFileName1);
        // Multilevel file
        std::string plotFileName2;
        pp.get("multilevelfile", plotFileName2);

        PlotFileData pf1(plotFileName1);
        PlotFileData pf2(plotFileName2);
        int finest_level = pf2.finestLevel();
        int Nlev = finest_level + 1;

        // Refinement ratio between levels of multilevelfile
        int ratio = 2;
        pp.query("ratio", ratio);

        // Ratio of finest level to reference
        int referenceRefinement = 2;
        pp.query("referenceratio", referenceRefinement);

        const auto& dx_1 = pf1.cellSize(0);
        const auto& dx_2 = pf2.cellSize(Nlev - 1);
        amrex::Print() << "dx_1 = " << dx_1 << "; dx_2 = " << dx_2 << std::endl;
        bool not_match = AMREX_D_TERM(
            dx_1[0] * referenceRefinement != dx_2[0],
            || dx_1[1] * referenceRefinement != dx_2[1],
            || dx_1[2] * referenceRefinement != dx_2[2]);
        if (not_match) {
            amrex::Abort(
                "ERROR: grid dx at finest level does not match ratio to "
                "reference dx");
        }

        int nComp, sComp;
        const auto& varNames = pf1.varNames();
        amrex::Vector<std::string> varnames(1);
        if (pp.countval("varname") != 0) {
            std::string varname;
            pp.get("varname", varname);
            nComp = 1;
            sComp = -1;
            for (int j = 0; j < varNames.size() && sComp == -1; ++j) {
                if (varNames[j] == varname) sComp = j;
            }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                sComp != -1, "variable not found in plotfile");
            varnames[0] = varname;
        } else {
            nComp = varNames.size();
            sComp = 0;
            varnames = varNames;
        }

        for (int j = 0; j < varnames.size(); ++j) {
            amrex::Print() << "varnames[" << j << "] = " << varnames[j]
                           << std::endl;
        }

        amrex::Vector<amrex::MultiFab*> data(Nlev);
        int coordSys = pf1.coordSys();
        const auto problo = pf1.probLo();
        const auto probhi = pf1.probHi();
        amrex::RealBox rb(problo, probhi);
        amrex::Array<int, AMREX_SPACEDIM> is_per = {{0, 0, 1}};
        amrex::Vector<int> levelSteps(Nlev);
        amrex::Vector<amrex::IntVect> refRatio(Nlev - 1);
        amrex::Vector<amrex::Geometry> geoms(Nlev);
        amrex::Box fpd = pf1.probDomain(0);
        amrex::Geometry fgeom(fpd, rb, coordSys, is_per);
        auto mf1 = (nComp == 1 ? pf1.get(0, varnames[0]) : pf1.get(0));
        const auto& fba = mf1.boxArray();
        const auto& dm = mf1.DistributionMap();

        amrex::Real L2error = 0.0;
        amrex::Real L2exact = 0.0;
        amrex::Real Linferror = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            int r = std::pow(ratio, finest_level - lev) * referenceRefinement;
            amrex::Print() << "coarsening reference solution (level " << Nlev
                           << ") to level " << lev
                           << "; ratio from reference = " << r << std::endl;

            amrex::Box cpd = amrex::Box(fpd).coarsen(r);
            geoms[lev].define(cpd, rb, coordSys, is_per);
            levelSteps[lev] = pf1.levelStep(0);
            if (lev < finest_level) {
                refRatio[lev] =
                    amrex::IntVect(AMREX_D_DECL(ratio, ratio, ratio));
                amrex::Print() << "refRatio[" << lev << "] = " << refRatio[lev]
                               << std::endl;
            }

            amrex::BoxArray cba = amrex::BoxArray(fba).coarsen(r);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                amrex::BoxArray(cba).refine(r) == fba,
                "BoxArray doesn't coarsen correctly");

            data[lev] = new amrex::MultiFab(cba, dm, nComp, 0);
            average_down(mf1, *data[lev], fgeom, geoms[lev], 0, nComp, r);

            auto mf2_lev =
                (nComp == 1 ? pf2.get(lev, varnames[0]) : pf2.get(lev));

            amrex::MultiFab mf_b;
            mf_b.define(mf2_lev.boxArray(), mf2_lev.DistributionMap(), 1, 0);
            mf_b.ParallelCopy(*data[lev]);

            // create mask to avoid counting the error in the refined area
            // multiple times

            amrex::iMultiFab level_mask;
            if (lev < finest_level) {
                auto mf2_levp1 =
                    (nComp == 1 ? pf2.get(lev + 1, varnames[0])
                                : pf2.get(lev + 1));
                level_mask = makeFineMask(
                    mf2_lev.boxArray(), mf2_lev.DistributionMap(),
                    mf2_levp1.boxArray(), amrex::IntVect(2), 1, 0);
            } else {
                level_mask.define(
                    mf2_lev.boxArray(), mf2_lev.DistributionMap(), 1, 0,
                    amrex::MFInfo());
                level_mask.setVal(1);
            }

            // compare with refined solution

            const auto& dx = pf2.cellSize(lev);
            for (amrex::MFIter mfi(mf2_lev); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
                auto mask = level_mask.array(mfi);
                auto var1 = mf_b.array(mfi);
                auto var2 = mf2_lev.array(mfi);
                amrex::Loop(
                    vbx, [=, &L2error, &L2exact,
                          &Linferror](int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                        const amrex::Real diff =
                            mask(i, j, k) *
                            (var2(i, j, k, 0) - var1(i, j, k, 0));
                        L2error += std::pow(diff, 2);
                        L2exact +=
                            std::pow(mask(i, j, k) * var1(i, j, k, 0), 2);

                        if (std::abs(diff) > Linferror) {
                            Linferror = std::abs(diff);
                        }
                    });
            }
        }

        amrex::ParallelDescriptor::ReduceRealSum(L2error);
        amrex::ParallelDescriptor::ReduceRealSum(L2exact);
        amrex::ParallelDescriptor::ReduceRealMax(Linferror);

        L2error = std::sqrt(L2error);
        L2exact = std::sqrt(L2exact);
        amrex::Print() << "L2error = " << L2error << std::endl;
        amrex::Print() << "L2exact = " << L2exact << std::endl;
        amrex::Print() << "relL2error = " << L2error / L2exact << std::endl;
        amrex::Print() << "Linferror = " << Linferror << std::endl;

        // Create multilevel plotfile from coarsened reference solution

        // Vector<const MultiFab*> dat(Nlev);
        // for (int i=0; i<Nlev; ++i) {
        //  dat[i] = data[i];
        //}
        // WriteMultiLevelPlotfile(getFileRoot(plotFileName1)+"_crse",Nlev,dat,varnames,geoms,pf1.time(),levelSteps,refRatio);
    }
    amrex::Finalize();
    return 0;
}
