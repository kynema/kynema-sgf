#include "src/equation_systems/passive_scalar/passive_scalar.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"

namespace kynema_sgf::pde {

template class PDESystem<PassiveScalar, fvm::Godunov>;
template class PDESystem<PassiveScalar, fvm::MOL>;

} // namespace kynema_sgf::pde
