#include <utility>

#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/boundary_conditions/BCInterface.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

FieldInfo::FieldInfo(
    std::string basename,
    const int ncomp,
    const int ngrow,
    const int nstates,
    const FieldLoc floc)
    : m_basename(std::move(basename))
    , m_ncomp(ncomp)
    , m_ngrow(ngrow)
    , m_nstates(nstates)
    , m_floc(floc)
    , m_bc_values(
          static_cast<long>(AMREX_SPACEDIM) * 2,
          amrex::Vector<amrex::Real>(ncomp, 0.0_rt))
    , m_bc_values_dview(static_cast<long>(ncomp) * AMREX_SPACEDIM * 2)
    , m_bcrec(ncomp)
    , m_bcrec_d(ncomp)
    , m_states(FieldInfo::max_field_states, nullptr)
{
    for (int i = 0; i < AMREX_SPACEDIM * 2; ++i) {
        m_bc_type[i] = BC::undefined;
        m_bc_values_d[i] = nullptr;
    }
}

FieldInfo::~FieldInfo() = default;

bool FieldInfo::bc_initialized()
{
    bool has_undefined = false;
    for (int i = 0; i < AMREX_SPACEDIM * 2; ++i) {
        if (m_bc_type[i] == BC::undefined) {
            has_undefined = true;
        }
    }

    // Check that BC has been initialized properly
    bool has_bogus = false;
    for (int dir = 0; dir < m_ncomp; ++dir) {
        const auto* bcrec = m_bcrec[dir].vect();
        for (int i = 0; i < AMREX_SPACEDIM * 2; ++i) {
            if (bcrec[i] == amrex::BCType::bogus) {
                has_bogus = true;
            }
        }
    }

    return ((!has_undefined) && (!has_bogus));
}

void FieldInfo::copy_bc_to_device()
{
    if (!bc_initialized()) {
        amrex::Abort("Invalid BC type encountered");
    }

    amrex::Vector<amrex::Real> h_data(
        static_cast<long>(m_ncomp) * AMREX_SPACEDIM * 2);

    // Copy data to a flat array for transfer to device
    {
        amrex::Real* hp = h_data.data();
        for (const auto& v : m_bc_values) {
            for (const auto& x : v) {
                *(hp++) = x;
            }
        }
    }

    // Transfer BC values to device
    amrex::Real* ptr = m_bc_values_dview.data();
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_data.begin(), h_data.end(),
        m_bc_values_dview.begin());

    for (int i = 0; i < AMREX_SPACEDIM * 2; ++i) {
        m_bc_values_d[i] = ptr;
        ptr += m_ncomp;
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_bcrec.begin(), m_bcrec.end(),
        m_bcrec_d.begin());

    m_bc_copied_to_device = true;
}

Field::Field(
    FieldRepo& repo,
    std::string name,
    std::shared_ptr<FieldInfo> info,
    const unsigned fid,
    const FieldState state)
    : m_repo(repo)
    , m_name(std::move(name))
    , m_info(std::move(info))
    , m_id(fid)
    , m_state(state)
{}

Field::~Field() = default;

Field& Field::state(const FieldState fstate)
{
    const auto& fstates = m_info->m_states;
    const int fint = static_cast<int>(fstate);
    AMREX_ASSERT(fstates[fint] != nullptr);
    return *fstates[fint];
}

const Field& Field::state(const FieldState fstate) const
{
    const auto& fstates = m_info->m_states;
    const int fint = static_cast<int>(fstate);
    AMREX_ASSERT(fstates[fint] != nullptr);
    return *fstates[fint];
}

amrex::MultiFab& Field::operator()(int lev)
{
    BL_ASSERT(lev < m_repo.num_active_levels());
    return m_repo.get_multifab(m_id, lev);
}

const amrex::MultiFab& Field::operator()(int lev) const
{
    BL_ASSERT(lev < m_repo.num_active_levels());
    return m_repo.get_multifab(m_id, lev);
}

amrex::Vector<amrex::MultiFab*> Field::vec_ptrs()
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<amrex::MultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(&m_repo.get_multifab(m_id, lev));
    }
    return ret;
}

amrex::Vector<const amrex::MultiFab*> Field::vec_const_ptrs() const
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<const amrex::MultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(
            static_cast<const amrex::MultiFab*>(
                &m_repo.get_multifab(m_id, lev)));
    }
    return ret;
}

void Field::fillpatch(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost)
{
    BL_PROFILE("amr-wind::Field::fillpatch 2");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    auto& fop = *(m_info->m_fillpatch_op);

    fop.fillpatch(lev, time, mfab, nghost, field_state());
}

void Field::fillpatch_from_coarse(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost)
{
    BL_PROFILE("amr-wind::Field::fillpatch_from_coarse");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    auto& fop = *(m_info->m_fillpatch_op);

    fop.fillpatch_from_coarse(lev, time, mfab, nghost, field_state());
}

void Field::fillpatch(const amrex::Real time, const amrex::IntVect ng)
{
    BL_PROFILE("amr-wind::Field::fillpatch");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    auto& fop = *(m_info->m_fillpatch_op);
    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        fop.fillpatch(
            lev, time, m_repo.get_multifab(m_id, lev), ng, field_state());
    }
}

void Field::fillpatch(const amrex::Real time) { fillpatch(time, num_grow()); }

void Field::fillpatch_sibling_fields(
    const amrex::Real time,
    const amrex::IntVect ng,
    amrex::Array<Field*, AMREX_SPACEDIM>& fields) const
{
    BL_PROFILE("amr-wind::Field::fillpatch array");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    BL_ASSERT(m_info->m_ncomp == static_cast<int>(fields.size()));
    auto& fop = *(m_info->m_fillpatch_op);
    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mfabs = {AMREX_D_DECL(
            &((*fields[0])(lev)), &((*fields[1])(lev)), &((*fields[2])(lev)))};
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> cfabs;
        if (lev > 0) {
            for (int i = 0; std::cmp_less(i, fields.size()); i++) {
                cfabs[i] = &((*fields[i])(lev - 1));
            }
        }

        fop.fillpatch_sibling_fields(
            lev, time, mfabs, mfabs, cfabs, ng, m_info->m_bcrec,
            m_info->m_bcrec, field_state());
    }
}

void Field::fillphysbc(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& ng)
{
    BL_PROFILE("amr-wind::Field::fillphysbc");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    auto& fop = *(m_info->m_fillpatch_op);
    fop.fillphysbc(lev, time, mfab, ng, field_state());
}

void Field::fillphysbc(const amrex::Real time, const amrex::IntVect ng)
{
    BL_PROFILE("amr-wind::Field::fillphysbc");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    auto& fop = *(m_info->m_fillpatch_op);
    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        fop.fillphysbc(
            lev, time, m_repo.get_multifab(m_id, lev), ng, field_state());
    }
}

void Field::fillphysbc(const amrex::Real time) { fillphysbc(time, num_grow()); }

void Field::apply_bc_funcs(const FieldState rho_state)
{
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    for (const auto& func : m_info->m_bc_func) {
        (*func)(*this, rho_state);
    }
}

void Field::set_inflow(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& ng)
{
    BL_PROFILE("amr-wind::Field::set_inflow");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    auto& fop = *(m_info->m_fillpatch_op);
    fop.set_inflow(lev, time, mfab, ng, field_state());
}

void Field::set_inflow_sibling_fields(
    const int lev,
    const amrex::Real time,
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mfabs)
{
    BL_PROFILE("amr-wind::Field::set_inflow_sibling_fields");
    BL_ASSERT(m_info->m_fillpatch_op);
    BL_ASSERT(m_info->bc_initialized() && m_info->m_bc_copied_to_device);
    auto& fop = *(m_info->m_fillpatch_op);
    fop.set_inflow_sibling_fields(lev, time, mfabs);
}

void Field::advance_states()
{
    BL_PROFILE("amr-wind::Field::advance_states");
    if (num_time_states() < 2) {
        return;
    }

    for (int i = num_time_states() - 1; i > 0; --i) {
        const auto sold = static_cast<FieldState>(i);
        const auto snew = static_cast<FieldState>(i - 1);
        auto& old_field = state(sold);
        const auto& new_field = state(snew);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            amrex::MultiFab::Copy(
                old_field(lev), new_field(lev), 0, 0, num_comp(), num_grow());
        }
    }
}

void Field::copy_state(FieldState to_state, FieldState from_state)
{
    BL_PROFILE("amr-wind::Field::copy_state");
    auto& to_field = state(to_state);
    const auto& from_field = state(from_state);

    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        amrex::MultiFab::Copy(
            to_field(lev), from_field(lev), 0, 0, num_comp(), num_grow());
    }
}

Field& Field::create_state(const FieldState fstate)
{
    const int sid = static_cast<int>(fstate);
    if (m_info->m_states[sid] == nullptr) {
        m_repo.create_state(*this, fstate);
    }

    return state(fstate);
}

void Field::setVal(amrex::Real value)
{
    BL_PROFILE("amr-wind::Field::setVal 1");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value);
    }
}

void Field::setVal(amrex::Real value, int start_comp, int num_comp, int nghost)
{
    BL_PROFILE("amr-wind::Field::setVal 2");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value, start_comp, num_comp, nghost);
    }
}

void Field::setVal(const amrex::Vector<amrex::Real>& values, int nghost)
{
    BL_PROFILE("amr-wind::Field::setVal 3");
    AMREX_ASSERT(num_comp() == static_cast<int>(values.size()));

    // Update 1 component at a time
    const int ncomp = 1;
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        auto& mf = operator()(lev);
        for (int ic = 0; ic < num_comp(); ++ic) {
            amrex::Real value = values[ic];
            mf.setVal(value, ic, ncomp, nghost);
        }
    }
}

void Field::set_default_fillpatch_bc(
    const SimTime& time, amrex::BCType::mathematicalBndryTypes bctype)
{
    if (!m_info->bc_initialized()) {
        BCFillPatchExtrap bc_op(*this, bctype);
        bc_op();
    }

    if (!m_info->m_fillpatch_op) {
        register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
            repo().mesh(), time);
    }
}

void Field::to_uniform_space()
{
    if (m_info->m_ncomp < AMREX_SPACEDIM) {
        amrex::Abort("Trying to transform a non-vector field:" + m_name);
    }
    if (m_mesh_mapped) {
        amrex::Print() << "WARNING: Field already in uniform mesh space: "
                       << m_name << '\n';
        return;
    }

    const auto& mesh_fac = m_repo.get_mesh_mapping_field(m_info->m_floc);
    const auto& mesh_detJ = m_repo.get_mesh_mapping_det_j(m_info->m_floc);

    // scale velocity to accommodate for mesh mapping -> U^bar = U * J/fac
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        const auto& fac = mesh_fac(lev).const_arrays();
        const auto& detJ = mesh_detJ(lev).const_arrays();
        const auto& field = operator()(lev).arrays();
        amrex::ParallelFor(
            mesh_fac(lev), num_grow(), operator()(lev).nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                field[nbx](i, j, k, n) *=
                    detJ[nbx](i, j, k) / fac[nbx](i, j, k, n);
            });
    }
    amrex::Gpu::streamSynchronize();
    m_mesh_mapped = true;
}

void Field::to_stretched_space()
{
    if (m_info->m_ncomp < AMREX_SPACEDIM) {
        amrex::Abort("Trying to transform a non-vector field:" + m_name);
    }
    if (!m_mesh_mapped) {
        amrex::Print() << "WARNING: Field already in stretched mesh space: "
                       << m_name << '\n';
        return;
    }

    const auto& mesh_fac = m_repo.get_mesh_mapping_field(m_info->m_floc);
    const auto& mesh_detJ = m_repo.get_mesh_mapping_det_j(m_info->m_floc);

    // scale field back to stretched mesh -> U = U^bar * fac/J
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {

        const auto& fac = mesh_fac(lev).const_arrays();
        const auto& detJ = mesh_detJ(lev).const_arrays();
        const auto& field = operator()(lev).arrays();
        amrex::ParallelFor(
            mesh_fac(lev), num_grow(), operator()(lev).nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                field[nbx](i, j, k, n) *=
                    fac[nbx](i, j, k, n) / detJ[nbx](i, j, k);
            });
    }
    amrex::Gpu::streamSynchronize();
    m_mesh_mapped = false;
}

} // namespace amr_wind
