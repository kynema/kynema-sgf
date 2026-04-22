#include "amr-wind/equation_systems/passive_scalar/passive_scalar.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"

namespace amr_wind::pde {

template class PDESystem<PassiveScalar, fvm::Godunov>;
template class PDESystem<PassiveScalar, fvm::MOL>;

} // namespace amr_wind::pde
