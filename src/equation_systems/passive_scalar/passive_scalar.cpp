#include "src/equation_systems/passive_scalar/passive_scalar.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"

namespace amr_wind::pde {

template class PDESystem<PassiveScalar, fvm::Godunov>;
template class PDESystem<PassiveScalar, fvm::MOL>;

} // namespace amr_wind::pde
