#ifndef DAMPINGLAYER_H
#define DAMPINGLAYER_H

#include "src/core/Physics.H"
#include "src/core/Field.H"
#include "src/CFDSim.H"
#include "AMReX_Gpu.H"

namespace kynema_sgf::damping_layer {

enum class BlendingFunctionType : std::uint8_t {
    Linear = 0,
    Quadratic = 1,
    Exponential = 2,
    Cosine = 3
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE BlendingFunctionType
string_to_blending_function_type(const std::string& op_str)
{
    std::string op = amrex::toLower(op_str);
    if (op == "linear") {
        return BlendingFunctionType::Linear;
    }
    if (op == "quadratic") {
        return BlendingFunctionType::Quadratic;
    }
    if (op == "exponential") {
        return BlendingFunctionType::Exponential;
    }
    if (op == "cosine") {
        return BlendingFunctionType::Cosine;
    }
    amrex::Abort(
        "Invalid blending function type: " + op_str +
        ". Must be one of: linear, quadratic, exponential, cosine.");
    // Unreachable, but avoids compiler warning
    return BlendingFunctionType::Cosine;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real blending_function(
    amrex::Real position_nondim, BlendingFunctionType function_type)
{
    switch (function_type) {
    case BlendingFunctionType::Linear:
        return 1.0_rt - position_nondim;
    case BlendingFunctionType::Quadratic:
        return 1.0_rt - position_nondim * position_nondim;
    case BlendingFunctionType::Exponential: // exponential difference and ratio
        return 1.0_rt - (std::expm1(std::pow(position_nondim, 3.5_rt)) /
                         std::expm1(1.0_rt));
    case BlendingFunctionType::Cosine:
    default:
        return 0.5_rt +
               0.5_rt *
                   std::cos(std::numbers::pi_v<amrex::Real> * position_nondim);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real damping_calc(
    amrex::Real distance_from_boundary,
    amrex::Real layer_thickness,
    amrex::Real blending_fraction,
    BlendingFunctionType function_type)
{
    const amrex::Real length_full_damp =
        layer_thickness * (1.0_rt - blending_fraction);
    const amrex::Real blend_position_nondim =
        (distance_from_boundary - length_full_damp) /
        (layer_thickness * blending_fraction);
    if (blend_position_nondim <= 0.0_rt) {
        return 1.0_rt;
    } else if (blend_position_nondim <= 1.0_rt) {
        return blending_function(blend_position_nondim, function_type);
    } else {
        return 0.0_rt;
    }
}

class DampingLayer : public Physics::Register<DampingLayer>
{
public:
    static std::string identifier() { return "DampingLayer"; }

    explicit DampingLayer(CFDSim& sim);

    ~DampingLayer() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override {}

    void post_init_actions() override {}

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    const amrex::Vector<const std::string> m_bc_names = {"xlo", "xhi", "ylo",
                                                         "yhi", "zlo", "zhi"};

    //! Quantities for constructing damping layers
    amrex::Gpu::DeviceVector<std::string> m_field_names;
    amrex::Gpu::DeviceVector<amrex::Array<amrex::Real, 6>> m_layers_thickness;
    amrex::Gpu::DeviceVector<amrex::Array<amrex::Real, 6>>
        m_layers_blending_fraction;
    amrex::Gpu::DeviceVector<amrex::Array<amrex::Real, 4>> m_layers_min_height;
    amrex::Gpu::DeviceVector<amrex::Array<BlendingFunctionType, 6>>
        m_layers_blending_function_type;
};
} // namespace kynema_sgf::damping_layer

#endif
