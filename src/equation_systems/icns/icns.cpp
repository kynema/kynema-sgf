#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "amr-wind/equation_systems/icns/icns_advection.H"
#include "amr-wind/equation_systems/icns/icns_diffusion.H"
#include "amr-wind/equation_systems/icns/icns_bcop.H"

namespace amr_wind::pde {

template class PDESystem<ICNS, ::amr_wind::fvm::Godunov>;
template class PDESystem<ICNS, ::amr_wind::fvm::MOL>;

} // namespace amr_wind::pde
