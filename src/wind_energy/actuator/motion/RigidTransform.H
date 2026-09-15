#ifndef ACTUATOR_RIGID_TRANSFORM_H
#define ACTUATOR_RIGID_TRANSFORM_H

#include "src/core/vs/vector_space.H"

namespace kynema_sgf::actuator::motion {

/** Position and orientation mapping body-frame coordinates to the CFD frame. */
struct RigidTransform
{
    vs::Vector position{vs::Vector::zero()};
    vs::Tensor orientation{vs::Tensor::identity()};

    [[nodiscard]] vs::Vector apply_point(const vs::Vector& point) const
    {
        return position + (orientation & point);
    }

    [[nodiscard]] vs::Vector apply_vector(const vs::Vector& vector) const
    {
        return orientation & vector;
    }
};

} // namespace kynema_sgf::actuator::motion

#endif
