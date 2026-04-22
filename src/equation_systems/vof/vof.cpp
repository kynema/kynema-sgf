#include "src/equation_systems/vof/vof.H"
#include "src/equation_systems/vof/vof_advection.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"
#include "src/equation_systems/vof/vof_ops.H"
#include "src/equation_systems/vof/vof_bcop.H"

namespace amr_wind::pde {

template class PDESystem<VOF, fvm::Godunov>;
template class PDESystem<VOF, fvm::MOL>;

} // namespace amr_wind::pde
