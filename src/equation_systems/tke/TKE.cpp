#include "src/equation_systems/tke/TKE.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"
#include "src/equation_systems/tke/tke_ops.H"

namespace kynema_sgf::pde {

template class PDESystem<TKE, fvm::Godunov>;
template class PDESystem<TKE, fvm::MOL>;

} // namespace kynema_sgf::pde
