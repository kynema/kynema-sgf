#include "src/equation_systems/levelset/levelset.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"
#include "src/equation_systems/levelset/levelset_ops.H"

namespace kynema_sgf::pde {

template class PDESystem<Levelset, fvm::Godunov>;
template class PDESystem<Levelset, fvm::MOL>;

} // namespace kynema_sgf::pde
