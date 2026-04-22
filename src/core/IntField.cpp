#include <utility>

#include "amr-wind/core/IntField.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {

IntField::IntField(
    FieldRepo& repo,
    std::string name,
    const unsigned fid,
    const int ncomp,
    const int ngrow,
    const FieldLoc floc)
    : m_repo(repo)
    , m_name(std::move(name))
    , m_id(fid)
    , m_ncomp(ncomp)
    , m_ngrow(ngrow)
    , m_floc(floc)
{}

amrex::iMultiFab& IntField::operator()(int lev)
{
    AMREX_ASSERT(lev < m_repo.num_active_levels());
    return m_repo.get_int_fab(m_id, lev);
}

const amrex::iMultiFab& IntField::operator()(int lev) const
{
    AMREX_ASSERT(lev < m_repo.num_active_levels());
    return m_repo.get_int_fab(m_id, lev);
}

amrex::Vector<amrex::iMultiFab*> IntField::vec_ptrs()
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<amrex::iMultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(&m_repo.get_int_fab(m_id, lev));
    }
    return ret;
}

amrex::Vector<const amrex::iMultiFab*> IntField::vec_const_ptrs() const
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<const amrex::iMultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(
            static_cast<const amrex::iMultiFab*>(
                &m_repo.get_int_fab(m_id, lev)));
    }
    return ret;
}

void IntField::setVal(int value)
{
    BL_PROFILE("amr-wind::IntField::setVal 1");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value);
    }
}

void IntField::setVal(int value, int start_comp, int num_comp, int nghost)
{
    BL_PROFILE("amr-wind::IntField::setVal 2");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value, start_comp, num_comp, nghost);
    }
}

void IntField::setVal(const amrex::Vector<int>& values, int nghost)
{
    BL_PROFILE("amr-wind::IntField::setVal 3");
    AMREX_ASSERT(num_comp() == static_cast<int>(values.size()));

    // Update 1 component at a time
    const int ncomp = 1;
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        auto& mf = operator()(lev);
        for (int ic = 0; ic < num_comp(); ++ic) {
            int value = values[ic];
            mf.setVal(value, ic, ncomp, nghost);
        }
    }
}

} // namespace amr_wind
