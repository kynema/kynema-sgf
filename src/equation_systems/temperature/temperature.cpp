#include "src/equation_systems/temperature/temperature.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"

namespace amr_wind::pde {

template class PDESystem<Temperature, fvm::Godunov>;
template class PDESystem<Temperature, fvm::MOL>;

} // namespace amr_wind::pde
