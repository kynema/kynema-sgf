#include "src/physics/DampingLayer.H"
#include "src/utilities/math_ops.H"
#include "src/utilities/constants.H"
#include "src/utilities/IOManager.H"
#include "src/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::damping_layer {

DampingLayer::DampingLayer(CFDSim& sim) : m_repo(sim.repo()), m_mesh(sim.mesh())
{
    amrex::Vector<std::string> labels;
    amrex::ParmParse pp(identifier());
    pp.getarr("fields", labels);

    amrex::Vector<amrex::Array<amrex::Real, 6>> layers_thickness;
    amrex::Vector<amrex::Array<amrex::Real, 6>> layers_blending_fraction;
    amrex::Vector<amrex::Array<amrex::Real, 4>> layers_min_height;
    amrex::Vector<amrex::Array<BlendingFunctionType, 6>>
        layers_blending_function_type;

    for (const auto& lbl : labels) {
        const std::string key = identifier() + "." + lbl;

        amrex::Array<amrex::Real, 6> bc_thickness;
        amrex::Array<amrex::Real, 6> bc_blending_fraction;
        amrex::Array<amrex::Real, 4> bc_min_height;
        amrex::Array<BlendingFunctionType, 6> bc_blending_function_type;
        int bc_index = 0;
        for (const auto& name : m_bc_names) {
            // Get arguments specific to this boundary
            const std::string key_bc = key + "." + name;
            amrex::ParmParse pp_bc(key_bc);
            amrex::Real thickness = -1.0_rt;
            pp_bc.query("thickness", thickness);
            amrex::Real blending_fraction = 0.0_rt;
            pp_bc.query("blending_fraction", blending_fraction);
            amrex::Real min_height = constants::LOW_NUM;
            pp_bc.query("minimum_height", min_height);
            std::string blending_function_str = "cosine";
            pp_bc.query("blending_function_type", blending_function_str);
            BlendingFunctionType blending_function_type =
                string_to_blending_function_type(blending_function_str);
            // Abort statement if height is specified for a z boundary
            if (name == "zlo" || name == "zhi") {
                if (pp_bc.contains("minimum_height")) {
                    amrex::Abort(
                        "DampingLayer: minimum_height is not supported for z "
                        "boundaries.");
                }
            }
            // Create field to go with this boundary damping layer if thickness
            // > 0
            if (thickness > 0.0_rt) {
                const auto field_name = "damping_layer_" + lbl + "_" + name;
                m_repo.declare_field(field_name, 1, 1, 1);
                sim.io_manager().register_io_var(field_name);
                // Record damping layer parameters for field creation (next
                // step)
                bc_thickness[bc_index] = thickness;
                bc_blending_fraction[bc_index] = blending_fraction;
                bc_blending_function_type[bc_index] = blending_function_type;
                if (name != "zlo" && name != "zhi") {
                    bc_min_height[bc_index] = min_height;
                }
            }
            ++bc_index;
        }

        layers_thickness.emplace_back(bc_thickness);
        layers_blending_fraction.emplace_back(bc_blending_fraction);
        layers_min_height.emplace_back(bc_min_height);
        layers_blending_function_type.emplace_back(bc_blending_function_type);
    }

    const int nfields = static_cast<int>(layers_thickness.size());
    m_field_names.resize(nfields);
    m_layers_thickness.resize(nfields);
    m_layers_blending_fraction.resize(nfields);
    m_layers_min_height.resize(nfields);
    m_layers_blending_function_type.resize(nfields);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, labels.begin(), labels.end(),
        m_field_names.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, layers_thickness.begin(),
        layers_thickness.end(), m_layers_thickness.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, layers_blending_fraction.begin(),
        layers_blending_fraction.end(), m_layers_blending_fraction.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, layers_min_height.begin(),
        layers_min_height.end(), m_layers_min_height.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, layers_blending_function_type.begin(),
        layers_blending_function_type.end(),
        m_layers_blending_function_type.begin());
}

void DampingLayer::initialize_fields(int level, const amrex::Geometry& geom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    const int nfields = static_cast<int>(m_field_names.size());
    const std::string* field_names_ptr = m_field_names.data();
    const amrex::Array<amrex::Real, 6>* bc_thickness_ptr =
        m_layers_thickness.data();
    const amrex::Array<amrex::Real, 6>* bc_blending_fraction_ptr =
        m_layers_blending_fraction.data();
    const amrex::Array<amrex::Real, 4>* bc_min_height_ptr =
        m_layers_min_height.data();
    const amrex::Array<BlendingFunctionType, 6>* bc_blending_function_type_ptr =
        m_layers_blending_function_type.data();

    for (int field_idx = 0; field_idx < nfields; ++field_idx) {
        const std::string lbl = field_names_ptr[field_idx];

        const amrex::Array<amrex::Real, 6>& bc_thickness =
            bc_thickness_ptr[field_idx];
        const amrex::Array<amrex::Real, 6>& bc_blending_fraction =
            bc_blending_fraction_ptr[field_idx];
        const amrex::Array<amrex::Real, 4>& bc_min_height =
            bc_min_height_ptr[field_idx];
        const amrex::Array<BlendingFunctionType, 6>& bc_blending_function_type =
            bc_blending_function_type_ptr[field_idx];

        // Field pointer to the damping layer field for this boundary condition
        Field* damping_layer_ptr{nullptr};

        for (int bc_idx = 0; bc_idx < 6; ++bc_idx) {
            if (bc_thickness[bc_idx] > 0.0_rt) {
                // Form damping layer field name
                const auto field_name =
                    "damping_layer_" + lbl + "_" + m_bc_names[bc_idx];
                damping_layer_ptr = &m_repo.get_field(field_name);

                // Get other damping layer parameters for this boundary
                const amrex::Real blending_fraction =
                    bc_blending_fraction[bc_idx];
                const amrex::Real min_height =
                    (bc_idx < 4) ? bc_min_height[bc_idx] : constants::LOW_NUM;
                const BlendingFunctionType blending_function_type =
                    bc_blending_function_type[bc_idx];

                auto& damping_layer_mfab = (*damping_layer_ptr)(level);
                auto& damping_layer_arrs = damping_layer_mfab.arrays();

                amrex::ParallelFor(
                    damping_layer_mfab,
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k) noexcept {
                        const amrex::Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
                        const amrex::Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
                        const amrex::Real z = prob_lo[2] + (k + 0.5_rt) * dx[2];

                        // Determine the distance from the boundary
                        amrex::Real distance_from_boundary = 0.0_rt;
                        switch (bc_idx) {
                        case 0: // xlo
                            distance_from_boundary = x - prob_lo[0];
                            break;
                        case 1: // xhi
                            distance_from_boundary = prob_hi[0] - x;
                            break;
                        case 2: // ylo
                            distance_from_boundary = y - prob_lo[1];
                            break;
                        case 3: // yhi
                            distance_from_boundary = prob_hi[1] - y;
                            break;
                        case 4: // zlo
                            distance_from_boundary = z - prob_lo[2];
                            break;
                        case 5: // zhi
                            distance_from_boundary = prob_hi[2] - z;
                            break;
                        }

                        amrex::Real damping_coeff = damping_calc(
                            distance_from_boundary, bc_thickness[bc_idx],
                            blending_fraction, blending_function_type);

                        if (min_height > constants::LOW_NUM) {
                            const amrex::Real distance_from_zhi =
                                prob_hi[2] - z;
                            const amrex::Real vertical_thickness =
                                prob_hi[2] - min_height;
                            const amrex::Real vertical_damping_coeff =
                                damping_calc(
                                    distance_from_zhi, vertical_thickness,
                                    blending_fraction, blending_function_type);
                            damping_coeff =
                                std::min(damping_coeff, vertical_damping_coeff);
                        }
                        // Set the damping coefficient in the field
                        damping_layer_arrs[nbx](i, j, k, 0) = damping_coeff;
                    });
                amrex::Gpu::streamSynchronize();
            }
        }
    }
}

void DampingLayer::post_regrid_actions()
{
    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_mesh.Geom(lev));
    }
}

} // namespace kynema_sgf::damping_layer
