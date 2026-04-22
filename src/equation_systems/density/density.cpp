#include "src/equation_systems/density/density.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"
#include "src/equation_systems/density/density_ops.H"

namespace amr_wind::pde {

template class PDESystem<Density, fvm::Godunov>;
template class PDESystem<Density, fvm::MOL>;

} // namespace amr_wind::pde
