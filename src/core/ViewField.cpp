#include "src/core/ViewField.H"
#include "src/core/Field.H"
#include "src/core/ScratchField.H"
#include "src/core/FieldRepo.H"

namespace amr_wind {

template <typename T>
ViewField<T>::ViewField(T& src, const int scomp, const int ncomp)
    : m_src(src), m_scomp(scomp), m_ncomp(ncomp)
{
    int nlevels = src.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        m_data.emplace_back(
            amrex::MultiFab(src(lev), amrex::make_alias, m_scomp, m_ncomp));
    }
}

template class ViewField<Field>;
template class ViewField<ScratchField>;

} // namespace amr_wind
