#ifndef REFINECHECKPT_H
#define REFINECHECKPT_H

#include "src/incflo.H"

namespace amr_wind {
namespace tools {

class RefineCheckpt : public incflo
{
public:
    RefineCheckpt();

    virtual ~RefineCheckpt() = default;

    virtual void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

    void read_chkpt_file();

    void refine_chkpt_file();

    void run_utility();

private:
    amrex::Vector<amrex::BoxArray> m_orig_ba;
};

} // namespace tools
} // namespace amr_wind

#endif /* REFINECHECKPT_H */
