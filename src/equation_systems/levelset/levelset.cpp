#include "amr-wind/equation_systems/levelset/levelset.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"
#include "amr-wind/equation_systems/levelset/levelset_ops.H"

namespace amr_wind::pde {

template class PDESystem<Levelset, fvm::Godunov>;
template class PDESystem<Levelset, fvm::MOL>;

} // namespace amr_wind::pde
