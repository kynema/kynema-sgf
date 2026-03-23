#include "amr-wind/boundary_conditions/BCInterface.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/boundary_conditions/FixedGradientBC.H"
#include "amr-wind/boundary_conditions/MassInflowOutflowBC.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

BCIface::BCIface(Field& field) : m_field(field) {}

void BCIface::operator()(const amrex::Real value)
{
    amrex::Print() << "Initializing boundary conditions for " << m_field.name()
                   << '\n';
    set_default_value(value);
    read_bctype();
    set_bcrec();
    read_values();
    set_bcfuncs();
    m_field.copy_bc_to_device();
}

void BCIface::set_default_value(const amrex::Real value)
{
    auto& bcval = m_field.bc_values();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        for (int i = 0; i < m_field.num_comp(); ++i) {
            bcval[ori][i] = value;
        }
    }
}

void BCIface::read_bctype()
{
    const std::string key = m_field.name() + "_type";
    auto& ibctype = m_field.bc_type();
    const auto& geom = m_field.repo().mesh().Geom(0);
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        amrex::ParmParse pp(bcid);
        std::string bcstr = "null";
        pp.query("type", bcstr);
        pp.query(key, bcstr);
        bcstr = amrex::toLower(bcstr);

        // Protect against copy/paste errors where user intended to add a BC but
        // forgot to turn off periodic in that direction, or vice versa.
        if (geom.isPeriodic(ori.coordDir()) &&
            ((bcstr != "null") && (bcstr != "periodic"))) {
            amrex::Abort(
                "Setting " + bcstr + " BC on a periodic face " + bcid +
                " is not allowed");
        }

        // Process and quit early if this face is periodic
        if (geom.isPeriodic(ori.coordDir())) {
            ibctype[ori] = BC::periodic;
            continue;
        }

        if ((bcstr == "pressure_outflow") || (bcstr == "po")) {
            ibctype[ori] = BC::pressure_outflow;
        } else if ((bcstr == "mass_inflow") || (bcstr == "mi")) {
            ibctype[ori] = BC::mass_inflow;
        } else if ((bcstr == "mass_inflow_outflow") || (bcstr == "mio")) {
            ibctype[ori] = BC::mass_inflow_outflow;
        } else if ((bcstr == "no_slip_wall") || (bcstr == "nsw")) {
            ibctype[ori] = BC::no_slip_wall;
        } else if ((bcstr == "slip_wall") || (bcstr == "sw")) {
            ibctype[ori] = BC::slip_wall;
        } else if ((bcstr == "wall_model") || (bcstr == "wm")) {
            ibctype[ori] = BC::wall_model;
        } else if ((bcstr == "zero_gradient") || (bcstr == "zg")) {
            ibctype[ori] = BC::zero_gradient;
        } else if ((bcstr == "fixed_gradient") || (bcstr == "fg")) {
            ibctype[ori] = BC::fixed_gradient;
        } else if ((bcstr == "wave_generation") || (bcstr == "wg")) {
            ibctype[ori] = BC::wave_generation;
        } else if ((bcstr == "symmetric_wall") || (bcstr == "symw")) {
            ibctype[ori] = BC::symmetric_wall;
        } else {
            ibctype[ori] = BC::undefined;
        }

        if (ibctype[ori] == BC::undefined) {
            amrex::Abort(
                "Invalid BC specification for non-periodic boundary = " +
                bcstr);
        }
    }
}

void BCIface::set_bcfuncs()
{
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto bct = ibctype[ori];

        if (bct == BC::fixed_gradient) {
            m_field.register_custom_bc<FixedGradientBC>(ori);
        }

        if (((m_field.name() == "velocity") ||
             (m_field.name() == "temperature")) &&
            (bct == BC::mass_inflow_outflow)) {

            m_field.register_custom_bc<MassInflowOutflowBC>(ori);
        }
    }
}

amrex::Array<const std::string, 3> BCIface::get_dirichlet_udfs()
{
    const auto& fname = m_field.name();
    const auto& bctype = m_field.bc_type();
    const std::string inflow_key = fname + ".inflow_type";
    const std::string inflow_outflow_key = fname + ".inflow_outflow_type";
    const std::string wall_key = fname + ".wall_type";
    std::string inflow_udf{"ConstDirichlet"};
    std::string inflow_outflow_udf{"ConstDirichlet"};
    std::string wall_udf{"ConstDirichlet"};
    bool has_inflow_udf = false;
    bool has_inflow_outflow_udf = false;
    bool has_wall_udf = false;

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        const auto bct = bctype[ori];
        amrex::ParmParse pp(bcid);

        if (bct == BC::mass_inflow) {
            if (pp.contains(inflow_key)) {
                std::string val;
                pp.get(inflow_key, val);

                if (has_inflow_udf && (inflow_udf != val)) {
                    amrex::Abort(
                        "BC: Inflow UDF must be same for all inflow faces");
                } else {
                    inflow_udf = val;
                    has_inflow_udf = true;
                }
            }
        }

        if (bct == BC::mass_inflow_outflow) {
            if (pp.contains(inflow_outflow_key)) {
                std::string val;
                pp.get(inflow_outflow_key, val);

                if (has_inflow_outflow_udf && (inflow_outflow_udf != val)) {
                    amrex::Abort(
                        "BC: Inflow-outflow UDF must be same for all "
                        "inflow-outflow faces");
                } else {
                    inflow_outflow_udf = val;
                    has_inflow_outflow_udf = true;
                }
            }
        }

        if (bct == BC::slip_wall) {
            if (pp.contains(wall_key)) {
                std::string val;
                pp.get(wall_key, val);

                if (has_wall_udf && (wall_udf != val)) {
                    amrex::Abort(
                        "BC: Wall UDF must be same for all wall faces");
                } else {
                    wall_udf = val;
                    has_wall_udf = true;
                }
            }
        }
    }

    return {inflow_udf, inflow_outflow_udf, wall_udf};
}

void BCVelocity::set_bcrec()
{
    const auto& ibctype = m_field.bc_type();
    auto& bcrec = m_field.bcrec();

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

        case BC::pressure_outflow:
        case BC::zero_gradient:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            }
            break;

        case BC::wave_generation:
        case BC::mass_inflow:
        case BC::no_slip_wall:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::ext_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::ext_dir);
            }
            break;

        case BC::mass_inflow_outflow:
            m_field.set_inout_bndry();
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::direction_dependent);
            } else {
                set_bcrec_hi(dir, amrex::BCType::direction_dependent);
            }
            break;

        case BC::slip_wall:
        case BC::wall_model:
            if (side == amrex::Orientation::low) {
                // Tangential directions use higher-order extrapolation
                set_bcrec_lo(dir, amrex::BCType::hoextrap);
                // Normal direction uses dirichlet (ext-dir)
                bcrec[dir].setLo(dir, amrex::BCType::ext_dir);
            } else {
                // Tangential directions use higher-order extrapolation
                set_bcrec_hi(dir, amrex::BCType::hoextrap);
                // Normal direction uses dirichlet (ext-dir)
                bcrec[dir].setHi(dir, amrex::BCType::ext_dir);
            }
            break;
        case BC::symmetric_wall:
            if (side == amrex::Orientation::low) {
                // Tangential directions use first-order extrapolation
                set_bcrec_lo(dir, amrex::BCType::foextrap);
                // Normal direction uses dirichlet (ext-dir)
                bcrec[dir].setLo(dir, amrex::BCType::ext_dir);
            } else {
                // Tangential directions use first-order extrapolation
                set_bcrec_hi(dir, amrex::BCType::foextrap);
                // Normal direction uses dirichlet (ext-dir)
                bcrec[dir].setHi(dir, amrex::BCType::ext_dir);
            }
            break;

        default:
            amrex::Abort("Invalid incflo BC type encountered");
        }
    }
}

void BCVelocity::read_values()
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
        pp.queryarr(fname, bcval[ori], 0, ndim);
        if (bct == BC::no_slip_wall) {
            // Set normal component to zero
            bcval[ori][ori.coordDir()] = 0.0_rt;
        }
    }
}

void BCScalar::set_bcrec()
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

        case BC::pressure_outflow:
        case BC::zero_gradient:
        case BC::symmetric_wall:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            }
            break;

        case BC::wave_generation:
        case BC::mass_inflow:
        case BC::no_slip_wall:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::ext_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::ext_dir);
            }
            break;

        case BC::mass_inflow_outflow:
            m_field.set_inout_bndry();
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::direction_dependent);
            } else {
                set_bcrec_hi(dir, amrex::BCType::direction_dependent);
            }
            break;

        case BC::slip_wall:
        case BC::wall_model:
        case BC::fixed_gradient:
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::hoextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::hoextrap);
            }
            break;

        default:
            amrex::Abort("Invalid incflo BC type encountered");
        }
    }
}

void BCScalar::read_values()
{
    const auto& fname = m_field.name();
    const auto& bctype = m_field.bc_type();
    auto& bcval = m_field.bc_values();
    const int ndim = m_field.num_comp();
    const auto udfs = get_dirichlet_udfs();
    const auto const_dirichlet_inflow = udfs[0] == "ConstDirichlet";
    const auto const_dirichlet_inflow_outflow = udfs[1] == "ConstDirichlet";

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto& bcid = bcnames[ori];
        const auto bct = bctype[ori];

        amrex::ParmParse pp(bcid);
        if (((bct == BC::mass_inflow) && (const_dirichlet_inflow)) ||
            ((bct == BC::mass_inflow_outflow) &&
             (const_dirichlet_inflow_outflow))) {
            pp.getarr(fname, bcval[ori], 0, ndim);
        } else {
            pp.queryarr(fname, bcval[ori], 0, ndim);
        }
    }
}

void BCPressure::read_values()
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
        if ((bct == BC::pressure_outflow)) {
            pp.queryarr(fname, bcval[ori], 0, ndim);
        }
    }
}

void BCSrcTerm::set_bcrec()
{
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        if (bct == BC::periodic) {
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            }
        } else {
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::foextrap);
            } else {
                set_bcrec_hi(dir, amrex::BCType::foextrap);
            }
        }
    }
}

void BCFillPatchExtrap::set_bcrec()
{
    const auto& ibctype = m_field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        const auto side = ori.faceDir();
        const auto bct = ibctype[ori];
        const int dir = ori.coordDir();

        if (bct == BC::periodic) {
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, amrex::BCType::int_dir);
            } else {
                set_bcrec_hi(dir, amrex::BCType::int_dir);
            }
        } else {
            if (side == amrex::Orientation::low) {
                set_bcrec_lo(dir, m_extrap_type);
            } else {
                set_bcrec_hi(dir, m_extrap_type);
            }
        }
    }
}

} // namespace amr_wind
