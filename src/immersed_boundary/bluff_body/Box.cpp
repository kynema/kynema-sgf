#include "src/immersed_boundary/bluff_body/Box.H"
#include "src/immersed_boundary/bluff_body/box_ops.H"
#include "src/immersed_boundary/IBModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::ib {

template class IBModel<Box>;

} // namespace kynema_sgf::ib
