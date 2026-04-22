#include "amr-wind/core/IntScratchField.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_Gpu.H"
#include "AMReX_IArrayBox.H"
#include "AMReX_Geometry.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_iMultiFab.H"

namespace amr_wind {

void IntScratchField::setVal(int value)
{
    BL_PROFILE("amr-wind::IntScratchField::setVal 1");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value);
    }
}

} // namespace amr_wind
