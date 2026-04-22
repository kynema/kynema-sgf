#include "src/equation_systems/vof/vof.H"
#include "src/equation_systems/vof/vof_advection.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"
#include "src/equation_systems/vof/vof_ops.H"
#include "src/equation_systems/vof/vof_bcop.H"

namespace kynema_sgf::pde {

template class PDESystem<VOF, fvm::Godunov>;
template class PDESystem<VOF, fvm::MOL>;

} // namespace kynema_sgf::pde
