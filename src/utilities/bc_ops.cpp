#include "amr-wind/utilities/bc_ops.H"

std::pair<bool, bool>
amr_wind::utils::has_extdir(amrex::BCRec const* bcrec, int ncomp, int dir)
{
    std::pair<bool, bool> r{false, false};
    for (int n = 0; n < ncomp; ++n) {
        r.first = r.first or bcrec[n].lo(dir) == amrex::BCType::ext_dir;
        r.second = r.second or bcrec[n].hi(dir) == amrex::BCType::ext_dir;
    }
    return r;
}

std::pair<bool, bool>
amr_wind::utils::has_extdir_or_ho(amrex::BCRec const* bcrec, int ncomp, int dir)
{
    std::pair<bool, bool> r{false, false};
    for (int n = 0; n < ncomp; ++n) {
        r.first = r.first or (bcrec[n].lo(dir) == amrex::BCType::ext_dir) or
                  (bcrec[n].lo(dir) == amrex::BCType::hoextrap);
        r.second = r.second or (bcrec[n].hi(dir) == amrex::BCType::ext_dir) or
                   (bcrec[n].hi(dir) == amrex::BCType::hoextrap);
    }
    return r;
}
