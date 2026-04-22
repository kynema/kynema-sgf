#ifndef VOF_BCIFACE_H
#define VOF_BCIFACE_H

#include "amr-wind/boundary_conditions/BCInterface.H"

namespace amr_wind {

/** BC interface for VOF
 *  \ingroup vof
 */

class BCVOF : public BCIface
{
public:
    explicit BCVOF(Field& field) : BCIface(field) {}

protected:
    // Only the process of reading BCs needs to be modified
    void set_bcrec() override
    {
        const auto& ibctype = m_field.bc_type();

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const auto side = ori.faceDir();
            const auto bct = ibctype[ori];
            const int dir = ori.coordDir();

            switch (bct) {
            case BC::periodic:
                if (side == amrex::Orientation::low) {
                    set_bcrec_lo(dir, amrex::BCType::int_dir);
                } else {
                    set_bcrec_hi(dir, amrex::BCType::int_dir);
                }
                break;

            // For the wall cases, the fillpatch BC is effectively a
            // contact-line model
            case BC::pressure_outflow:
            case BC::zero_gradient:
            case BC::no_slip_wall:
            case BC::symmetric_wall:
            case BC::slip_wall:
            case BC::wall_model:
                if (side == amrex::Orientation::low) {
                    set_bcrec_lo(dir, amrex::BCType::foextrap);
                } else {
                    set_bcrec_hi(dir, amrex::BCType::foextrap);
                }
                break;

                // Never do hoextrap because it is unbounded

            case BC::wave_generation:
            case BC::mass_inflow:
            case BC::mass_inflow_outflow:
                if (side == amrex::Orientation::low) {
                    set_bcrec_lo(dir, amrex::BCType::ext_dir);
                } else {
                    set_bcrec_hi(dir, amrex::BCType::ext_dir);
                }
                break;

            // Fixed gradient BC is not allowed for VOF
            case BC::fixed_gradient:
            default:
                amrex::Abort(
                    "Invalid BC type for VOF encountered: fixed_gradient");
            }
        }
    }
    // Copied from BC interface for scalars
    void read_values() override
    {
        const auto& fname = m_field.name();
        const auto& bctype = m_field.bc_type();
        auto& bcval = m_field.bc_values();
        const int ndim = m_field.num_comp();
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const auto& bcid = bcnames[ori];
            const auto bct = bctype[ori];

            amrex::ParmParse pp(bcid);
            if (bct == BC::mass_inflow || bct == BC::mass_inflow_outflow) {
                pp.getarr(fname, bcval[ori], 0, ndim);
            } else {
                pp.queryarr(fname, bcval[ori], 0, ndim);
            }
        }
    }
};

class BCDensityMod : public BCIface
{
public:
    explicit BCDensityMod(Field& field) : BCIface(field) {}

protected:
    void set_bcrec() override
    {
        const auto& ibctype = m_field.bc_type();

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const auto side = ori.faceDir();
            const auto bct = ibctype[ori];
            const int dir = ori.coordDir();

            switch (bct) {
            case BC::periodic:
                if (side == amrex::Orientation::low) {
                    set_bcrec_lo(dir, amrex::BCType::int_dir);
                } else {
                    set_bcrec_hi(dir, amrex::BCType::int_dir);
                }
                break;

            // Both extrapolations are consistent between vof and density
            case BC::pressure_outflow:
            case BC::zero_gradient:
            case BC::no_slip_wall:
            case BC::symmetric_wall:
            case BC::slip_wall:
            case BC::wall_model:
                if (side == amrex::Orientation::low) {
                    set_bcrec_lo(dir, amrex::BCType::foextrap);
                } else {
                    set_bcrec_hi(dir, amrex::BCType::foextrap);
                }
                break;

                // Never do hoextrap because it is unbounded

            case BC::wave_generation:
            case BC::mass_inflow:
            case BC::mass_inflow_outflow:
                if (side == amrex::Orientation::low) {
                    set_bcrec_lo(dir, amrex::BCType::ext_dir);
                } else {
                    set_bcrec_hi(dir, amrex::BCType::ext_dir);
                }
                break;

            case BC::fixed_gradient:
            default:
                amrex::Abort(
                    "Invalid BC type for density (vof-compatible) encountered: "
                    "fixed_gradient");
            }
        }
    }

    void read_values() override
    {
        const auto& fname = m_field.name();
        const auto& bctype = m_field.bc_type();
        auto& bcval = m_field.bc_values();
        const int ndim = m_field.num_comp();
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const auto& bcid = bcnames[ori];
            const auto bct = bctype[ori];

            amrex::ParmParse pp(bcid);
            if (bct == BC::mass_inflow || bct == BC::mass_inflow_outflow) {
                pp.getarr(fname, bcval[ori], 0, ndim);
            } else {
                pp.queryarr(fname, bcval[ori], 0, ndim);
            }
        }
    }
};

} // namespace amr_wind

#endif /* VOF_BCIFACE_H */
