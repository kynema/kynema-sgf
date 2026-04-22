#include "src/equation_systems/temperature/temperature.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"

namespace kynema_sgf::pde {

template class PDESystem<Temperature, fvm::Godunov>;
template class PDESystem<Temperature, fvm::MOL>;

} // namespace kynema_sgf::pde
