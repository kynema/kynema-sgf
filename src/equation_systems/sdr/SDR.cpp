#include "amr-wind/equation_systems/sdr/SDR.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"
#include "amr-wind/equation_systems/sdr/sdr_ops.H"

namespace amr_wind::pde {

template class PDESystem<SDR, fvm::Godunov>;
template class PDESystem<SDR, fvm::MOL>;

} // namespace amr_wind::pde
