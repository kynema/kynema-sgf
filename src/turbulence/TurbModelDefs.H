#ifndef TURBMODELDEFS_H
#define TURBMODELDEFS_H

#include "src/turbulence/TurbulenceModel.H"
#include "src/transport_models/ConstTransport.H"
#include "src/transport_models/TwoPhaseTransport.H"

namespace kynema_sgf {

#define INSTANTIATE_TURBULENCE_MODEL(Model)                                    \
    template struct ::kynema_sgf::turbulence::TurbulenceModel::Register<         \
        ::kynema_sgf::turbulence::Model<transport::ConstTransport>>;             \
    template struct ::kynema_sgf::turbulence::TurbulenceModel::Register<         \
        ::kynema_sgf::turbulence::Model<transport::TwoPhaseTransport>>

} // namespace kynema_sgf

#endif /* TURBMODELDEFS_H */
