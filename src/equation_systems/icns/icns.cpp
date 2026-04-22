#include "src/equation_systems/icns/icns.H"
#include "src/equation_systems/icns/icns_ops.H"
#include "src/equation_systems/icns/icns_advection.H"
#include "src/equation_systems/icns/icns_diffusion.H"
#include "src/equation_systems/icns/icns_bcop.H"

namespace kynema_sgf::pde {

template class PDESystem<ICNS, ::kynema_sgf::fvm::Godunov>;
template class PDESystem<ICNS, ::kynema_sgf::fvm::MOL>;

} // namespace kynema_sgf::pde
