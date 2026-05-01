#include "src/equation_systems/sdr/SDR.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/BCOps.H"
#include "src/equation_systems/sdr/sdr_ops.H"

namespace kynema_sgf::pde {

template class PDESystem<SDR, fvm::Godunov>;
template class PDESystem<SDR, fvm::MOL>;

} // namespace kynema_sgf::pde
