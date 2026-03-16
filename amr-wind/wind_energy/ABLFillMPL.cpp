#include <utility>

#include "amr-wind/wind_energy/ABLFillMPL.H"

namespace amr_wind {

ABLFillMPL::ABLFillMPL(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time,
    const ABLModulatedPowerLaw& abl_mpl)
    : FieldFillPatchOps<FieldBCDirichlet>(
          field, mesh, time, FieldInterpolator::CellConsLinear)
    , m_abl_mpl(abl_mpl)
{}

ABLFillMPL::~ABLFillMPL() = default;

void ABLFillMPL::fillpatch(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillpatch(
        lev, time, mfab, nghost, fstate);

    if (m_field.base_name() == "velocity") {
        m_abl_mpl.set_velocity(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "temperature") {
        m_abl_mpl.set_temperature(lev, time, m_field, mfab);
    }
}

void ABLFillMPL::fillpatch_from_coarse(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillpatch_from_coarse(
        lev, time, mfab, nghost, fstate);

    if (m_field.base_name() == "velocity") {
        m_abl_mpl.set_velocity(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "temperature") {
        m_abl_mpl.set_temperature(lev, time, m_field, mfab);
    }
}

void ABLFillMPL::fillphysbc(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    const FieldState fstate)
{
    FieldFillPatchOps<FieldBCDirichlet>::fillphysbc(
        lev, time, mfab, nghost, fstate);

    if (m_field.base_name() == "velocity") {
        m_abl_mpl.set_velocity(lev, time, m_field, mfab);
    } else if (m_field.base_name() == "temperature") {
        m_abl_mpl.set_temperature(lev, time, m_field, mfab);
    }
}

void ABLFillMPL::fillpatch_sibling_fields(
    const int lev,
    const amrex::Real time,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& ffabs,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& cfabs,
    const amrex::IntVect& nghost,
    const amrex::Vector<amrex::BCRec>& bcrec,
    const amrex::Vector<amrex::BCRec>& /* unused */,
    const FieldState fstate)
{
    if (m_field.base_name() == "velocity") {
        // For an ABL MPL, we first just foextrap the mac velocities
        amrex::Vector<amrex::BCRec> lbcrec(m_field.num_comp());
        const auto& ibctype = m_field.bc_type();
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const auto side = ori.faceDir();
            const auto bct = ibctype[ori];
            const int dir = ori.coordDir();
            for (int i = 0; i < m_field.num_comp(); ++i) {
                if ((bct == BC::mass_inflow) ||
                    (bct == BC::mass_inflow_outflow)) {
                    if (side == amrex::Orientation::low) {
                        lbcrec[i].setLo(dir, amrex::BCType::foextrap);
                    } else {
                        lbcrec[i].setHi(dir, amrex::BCType::foextrap);
                    }
                } else {
                    if (side == amrex::Orientation::low) {
                        lbcrec[i].setLo(dir, bcrec[i].lo(dir));
                    } else {
                        lbcrec[i].setHi(dir, bcrec[i].hi(dir));
                    }
                }
            }
        }

        FieldFillPatchOps<FieldBCDirichlet>::fillpatch_sibling_fields(
            lev, time, mfabs, ffabs, cfabs, nghost, lbcrec, lbcrec, fstate);

        for (int i = 0; std::cmp_less(i, mfabs.size()); i++) {
            m_abl_mpl.set_velocity(lev, time, m_field, *mfabs[i], 0, i);
        }
    }
}

} // namespace amr_wind
