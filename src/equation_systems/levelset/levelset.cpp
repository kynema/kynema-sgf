#include "src/equation_systems/levelset/levelset.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"
#include "src/equation_systems/levelset/levelset_ops.H"

namespace amr_wind::pde {

template class PDESystem<Levelset, fvm::Godunov>;
template class PDESystem<Levelset, fvm::MOL>;

} // namespace amr_wind::pde
