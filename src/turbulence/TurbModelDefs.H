#ifndef TURBMODELDEFS_H
#define TURBMODELDEFS_H

#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/transport_models/ConstTransport.H"
#include "amr-wind/transport_models/TwoPhaseTransport.H"

namespace amr_wind {

#define INSTANTIATE_TURBULENCE_MODEL(Model)                                    \
    template struct ::amr_wind::turbulence::TurbulenceModel::Register<         \
        ::amr_wind::turbulence::Model<transport::ConstTransport>>;             \
    template struct ::amr_wind::turbulence::TurbulenceModel::Register<         \
        ::amr_wind::turbulence::Model<transport::TwoPhaseTransport>>

} // namespace amr_wind

#endif /* TURBMODELDEFS_H */
