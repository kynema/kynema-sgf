#ifndef OCEAN_WAVES_OPS_H
#define OCEAN_WAVES_OPS_H

#include "src/ocean_waves/OceanWavesTypes.H"
#include "AMReX_Vector.H"
#include "src/CFDSim.H"

namespace kynema_sgf::ocean_waves::ops {

template <typename WaveTheoryTrait, typename = void>
struct ReadInputsOp;

template <typename WaveTheoryTrait, typename = void>
struct InitDataOp;

template <typename WaveTheoryTrait, typename = void>
struct UpdateTargetFieldsOp;

template <typename WaveTheoryTrait, typename = void>
struct ApplyRelaxZonesOp;

template <typename WaveTheoryTrait, typename = void>
struct UpdateTargetFieldsSecondaryOp;

template <typename WaveTheoryTrait, typename = void>
struct ProcessOutputsOp;

} // namespace kynema_sgf::ocean_waves::ops

#endif /* OCEAN_WAVES_OPS_H */
