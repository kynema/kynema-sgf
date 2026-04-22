#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"
#include "amr-wind/equation_systems/tke/tke_ops.H"

namespace amr_wind::pde {

template class PDESystem<TKE, fvm::Godunov>;
template class PDESystem<TKE, fvm::MOL>;

} // namespace amr_wind::pde
