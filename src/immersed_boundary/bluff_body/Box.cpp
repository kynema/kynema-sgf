#include "amr-wind/immersed_boundary/bluff_body/Box.H"
#include "amr-wind/immersed_boundary/bluff_body/box_ops.H"
#include "amr-wind/immersed_boundary/IBModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ib {

template class IBModel<Box>;

} // namespace amr_wind::ib
