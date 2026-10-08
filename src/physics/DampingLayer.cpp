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
    amrex::ParmParse pp(identifier());
    pp.getarr("fields", m_field_names);

    for (const auto& lbl : m_field_names) {
        const std::string key = identifier() + "." + lbl;

        amrex::Array<amrex::Real, 6> bc_thickness;
        amrex::Array<amrex::Real, 6> bc_blending_fraction;
        amrex::Array<amrex::Real, 4> bc_min_height;
        amrex::Array<amrex::Real, 4> bc_vertical_blending_thickness;
        amrex::Array<BlendingFunctionType, 4>
            bc_vertical_blending_function_type;
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
            amrex::Real vert_blend_thickness = -1.0_rt;
            if (pp_bc.contains("minimum_height")) {
                pp_bc.get("vertical_blending_thickness", vert_blend_thickness);
            }
            std::string vert_blend_function_str = "cosine";
            pp_bc.query(
                "vertical_blending_function_type", vert_blend_function_str);
            BlendingFunctionType vert_blend_function_type =
                string_to_blending_function_type(vert_blend_function_str);
            // Abort statement if height is specified for a z boundary
            if (name == "zlo" || name == "zhi") {
                if (pp_bc.contains("minimum_height")) {
                    amrex::Abort(
                        "DampingLayer: minimum_height is not supported for z "
                        "boundaries.");
                }
                if (pp_bc.contains("vertical_blending_thickness")) {
                    amrex::Abort(
                        "DampingLayer: vertical_blending_thickness is not "
                        "supported for z boundaries.");
                }
                if (pp_bc.contains("vertical_blending_function_type")) {
                    amrex::Abort(
                        "DampingLayer: vertical_blending_function_type is not "
                        "supported for z boundaries.");
                }
            }
            // Create field to go with this boundary damping layer
            bc_thickness[bc_index] = thickness;
            if (thickness > 0.0_rt) {
                const auto field_name = "damping_layer_" + lbl + "_" + name;
                m_repo.declare_field(field_name, 1, 1, 1);
                sim.io_manager().register_io_var(field_name);
                // Record damping layer parameters for field creation
                bc_blending_fraction[bc_index] = blending_fraction;
                bc_blending_function_type[bc_index] = blending_function_type;
                if (name != "zlo" && name != "zhi") {
                    bc_min_height[bc_index] = min_height;
                    bc_vertical_blending_thickness[bc_index] =
                        vert_blend_thickness;
                    bc_vertical_blending_function_type[bc_index] =
                        vert_blend_function_type;
                }
            }
            ++bc_index;
        }

        m_layers_thickness.emplace_back(bc_thickness);
        m_layers_blending_fraction.emplace_back(bc_blending_fraction);
        m_layers_min_height.emplace_back(bc_min_height);
        m_layers_vertical_blending_thickness.emplace_back(
            bc_vertical_blending_thickness);
        m_layers_vertical_blending_function_type.emplace_back(
            bc_vertical_blending_function_type);
        m_layers_blending_function_type.emplace_back(bc_blending_function_type);
    }
}

void DampingLayer::initialize_fields(int level, const amrex::Geometry& geom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    const int nfields = static_cast<int>(m_field_names.size());
    const std::string* field_names_ptr = m_field_names.data();

    for (int field_idx = 0; field_idx < nfields; ++field_idx) {
        const std::string lbl = field_names_ptr[field_idx];

        const amrex::Array<amrex::Real, 6>& bc_thickness =
            m_layers_thickness[field_idx];
        const amrex::Array<amrex::Real, 6>& bc_blending_fraction =
            m_layers_blending_fraction[field_idx];
        const amrex::Array<amrex::Real, 4>& bc_min_height =
            m_layers_min_height[field_idx];
        const amrex::Array<amrex::Real, 4>& bc_vertical_blending_thickness =
            m_layers_vertical_blending_thickness[field_idx];
        const amrex::Array<BlendingFunctionType, 4>&
            bc_vertical_blending_function_type =
                m_layers_vertical_blending_function_type[field_idx];
        const amrex::Array<BlendingFunctionType, 6>& bc_blending_function_type =
            m_layers_blending_function_type[field_idx];

        // Field pointer to the damping layer field for this boundary condition
        Field* damping_layer_ptr{nullptr};

        for (uint bc_idx = 0; bc_idx < 6; ++bc_idx) {
            const amrex::Real thickness = bc_thickness[bc_idx];
            if (thickness > 0.0_rt) {
                // Form damping layer field name
                const auto field_name =
                    "damping_layer_" + lbl + "_" + m_bc_names[bc_idx];
                damping_layer_ptr = &m_repo.get_field(field_name);

                // Get other damping layer parameters for this boundary
                const amrex::Real blending_fraction =
                    bc_blending_fraction[bc_idx];
                const amrex::Real min_height =
                    (bc_idx < 4) ? bc_min_height[bc_idx] : constants::LOW_NUM;
                const amrex::Real vertical_blending_thickness =
                    (bc_idx < 4) ? bc_vertical_blending_thickness[bc_idx]
                                 : -1.0_rt;
                const BlendingFunctionType vertical_blending_function_type =
                    (bc_idx < 4) ? bc_vertical_blending_function_type[bc_idx]
                                 : BlendingFunctionType::Cosine;
                const BlendingFunctionType blending_function_type =
                    bc_blending_function_type[bc_idx];
                const auto coordinate_idx = static_cast<int>(bc_idx / 2);
                const bool high_side = (bc_idx % 2) == 1;

                auto& damping_layer_mfab = (*damping_layer_ptr)(level);
                auto damping_layer_arrs = damping_layer_mfab.arrays();

                amrex::ParallelFor(
                    damping_layer_mfab,
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k) noexcept {
                        const amrex::Real x =
                            prob_lo[0] + ((i + 0.5_rt) * dx[0]);
                        const amrex::Real y =
                            prob_lo[1] + ((j + 0.5_rt) * dx[1]);
                        const amrex::Real z =
                            prob_lo[2] + ((k + 0.5_rt) * dx[2]);

                        const amrex::Real coordinate = (coordinate_idx == 0) ? x
                                                       : (coordinate_idx == 1)
                                                           ? y
                                                           : z;
                        const amrex::Real distance_from_boundary =
                            high_side ? prob_hi[coordinate_idx] - coordinate
                                      : coordinate - prob_lo[coordinate_idx];

                        amrex::Real damping_coeff = damping_calc(
                            distance_from_boundary, thickness,
                            blending_fraction, blending_function_type);

                        if (min_height > constants::LOW_NUM) {
                            const amrex::Real distance_from_zhi =
                                prob_hi[2] - z;
                            const amrex::Real vertical_thickness =
                                prob_hi[2] - min_height;
                            const amrex::Real vertical_blending_fraction =
                                vertical_blending_thickness /
                                vertical_thickness;
                            const amrex::Real vertical_damping_coeff =
                                damping_calc(
                                    distance_from_zhi, vertical_thickness,
                                    vertical_blending_fraction,
                                    vertical_blending_function_type);
                            damping_coeff = amrex::min<amrex::Real>(
                                damping_coeff, vertical_damping_coeff);
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
