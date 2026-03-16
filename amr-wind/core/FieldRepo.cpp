#include <memory>

#include "amr-wind/core/FieldRepo.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

LevelDataHolder::LevelDataHolder()
    : m_factory(new amrex::FArrayBoxFactory())
    , m_int_fact(new amrex::DefaultFabFactory<amrex::IArrayBox>())
{}

void FieldRepo::make_new_level_from_scratch(
    int lev,
    amrex::Real /* time */,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("amr-wind::FieldRepo::make_new_level_from_scratch");
    m_leveldata[lev] = std::make_unique<LevelDataHolder>();

    allocate_field_data(
        ba, dm, *m_leveldata[lev], *(m_leveldata[lev]->m_factory));
    allocate_field_data(
        ba, dm, *m_leveldata[lev], *(m_leveldata[lev]->m_int_fact));

    m_is_initialized = true;
}

void FieldRepo::make_new_level_from_coarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("amr-wind::FieldRepo::make_level_from_coarse");
    std::unique_ptr<LevelDataHolder> ldata(new LevelDataHolder());

    allocate_field_data(ba, dm, *ldata, *(ldata->m_factory));
    allocate_field_data(ba, dm, *ldata, *(ldata->m_int_fact));

    for (auto& field : m_field_vec) {
        if (!field->fillpatch_on_regrid()) {
            continue;
        }

        field->fillpatch_from_coarse(lev, time, ldata->m_mfabs[field->id()], 0);
    }

    m_leveldata[lev] = std::move(ldata);
    m_is_initialized = true;
}

void FieldRepo::remake_level(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("amr-wind::FieldRepo::remake_level");
    std::unique_ptr<LevelDataHolder> ldata(new LevelDataHolder());

    allocate_field_data(ba, dm, *ldata, *(ldata->m_factory));
    allocate_field_data(ba, dm, *ldata, *(ldata->m_int_fact));

    for (auto& field : m_field_vec) {
        if (!field->fillpatch_on_regrid()) {
            continue;
        }

        field->fillpatch(lev, time, ldata->m_mfabs[field->id()], 0);
    }

    m_leveldata[lev] = std::move(ldata);
    m_is_initialized = true;
}

void FieldRepo::clear_level(int lev)
{
    BL_PROFILE("amr-wind::FieldRepo::clear_level");
    m_leveldata[lev].reset();
}

Field& FieldRepo::declare_field(
    const std::string& name,
    const int ncomp,
    const int ngrow,
    const int nstates,
    const FieldLoc floc)
{
    BL_PROFILE("amr-wind::FieldRepo::declare_field");
    // If the field is already registered check and return the fields
    {
        auto found = m_fid_map.find(name);
        if (found != m_fid_map.end()) {
            auto& field = *m_field_vec[found->second];

            if ((ncomp != field.num_comp()) ||
                (field.num_time_states() != nstates) ||
                (floc != field.field_location())) {
                amrex::Abort(
                    "Attempt to reregister field with inconsistent "
                    "parameters: " +
                    name);
            }
            return field;
        }
    }

    // Only allow states to be at times and not half steps
    if (nstates > (static_cast<int>(FieldState::NM1) + 1)) {
        amrex::Abort("Invalid number of states specified for field: " + name);
    }

    if (!field_impl::is_valid_field_name(name)) {
        amrex::Abort("Attempt to use reserved field name: " + name);
    }

    // Create the field data structures
    std::shared_ptr<FieldInfo> finfo(
        new FieldInfo(name, ncomp, ngrow, nstates, floc));
    for (int i = 0; i < nstates; ++i) {
        const auto fstate = static_cast<FieldState>(i);
        const std::string fname =
            field_impl::field_name_with_state(name, fstate);
        const unsigned fid = m_field_vec.size();

        // Create new field instance
        std::unique_ptr<Field> field(
            new Field(*this, fname, finfo, fid, fstate));
        // If declare field is called after mesh has been initialized create
        // field multifabs
        if (m_is_initialized) {
            allocate_field_data(*field);
        }

        // Add reference to states lookup
        finfo->m_states[i] = field.get();
        // Store the field instance
        m_field_vec.emplace_back(std::move(field));
        // Name to ID lookup map
        m_fid_map[fname] = fid;
    }

    // We want the New state returned when we have multiple states involved
    return *m_field_vec[m_fid_map[name]];
}

Field&
FieldRepo::get_field(const std::string& name, const FieldState fstate) const
{
    BL_PROFILE("amr-wind::FieldRepo::get_field");
    const auto fname = field_impl::field_name_with_state(name, fstate);
    const auto found = m_fid_map.find(fname);
    if (found == m_fid_map.end()) { // NOLINT(bugprone-branch-clone)
        amrex::Abort("Cannot find field: " + name);
        exit(1); // To appease the compiler
    } else if (found->second < static_cast<unsigned>(m_field_vec.size())) {
        return *m_field_vec[found->second];
    } else {
        amrex::Abort("Cannot find field: " + name);
        exit(1); // To appease the compiler
    }
}

Field& FieldRepo::get_mesh_mapping_field(FieldLoc floc) const
{
    Field* fac = nullptr;
    switch (floc) {
    case FieldLoc::CELL:
        fac = &(get_field("mesh_scaling_factor_cc"));
        break;
    case FieldLoc::NODE:
        fac = &(get_field("mesh_scaling_factor_nd"));
        break;
    case FieldLoc::XFACE:
        fac = &(get_field("mesh_scaling_factor_xf"));
        break;
    case FieldLoc::YFACE:
        fac = &(get_field("mesh_scaling_factor_yf"));
        break;
    case FieldLoc::ZFACE:
        fac = &(get_field("mesh_scaling_factor_zf"));
        break;
    default:
        amrex::Abort("Invalid field location");
    }
    return *fac;
}

Field& FieldRepo::get_mesh_mapping_det_j(FieldLoc floc) const
{
    Field* detJ = nullptr;
    switch (floc) {
    case FieldLoc::CELL:
        detJ = &(get_field("mesh_scaling_detJ_cc"));
        break;
    case FieldLoc::NODE:
        detJ = &(get_field("mesh_scaling_detJ_nd"));
        break;
    case FieldLoc::XFACE:
        detJ = &(get_field("mesh_scaling_detJ_xf"));
        break;
    case FieldLoc::YFACE:
        detJ = &(get_field("mesh_scaling_detJ_yf"));
        break;
    case FieldLoc::ZFACE:
        detJ = &(get_field("mesh_scaling_detJ_zf"));
        break;
    default:
        amrex::Abort("Invalid field location");
    }
    return *detJ;
}

bool FieldRepo::field_exists(
    const std::string& name, const FieldState fstate) const
{
    const auto fname = field_impl::field_name_with_state(name, fstate);
    const auto found = m_fid_map.find(fname);
    return (found != m_fid_map.end());
}

IntField& FieldRepo::declare_int_field(
    const std::string& name,
    const int ncomp,
    const int ngrow,
    const int nstates,
    const FieldLoc floc)
{
    BL_PROFILE("amr-wind::FieldRepo::declare_int_field");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nstates == 1, "Multiple states not supported for integer fields");

    // If the field is already registered check and return the fields
    {
        auto found = m_int_fid_map.find(name);
        if (found != m_int_fid_map.end()) {
            auto& field = *m_int_field_vec[found->second];

            if ((ncomp != field.num_comp()) ||
                (floc != field.field_location())) {
                amrex::Abort(
                    "Attempt to reregister field with inconsistent "
                    "parameters: " +
                    name);
            }
            return field;
        }
    }

    if (!field_impl::is_valid_field_name(name)) {
        amrex::Abort("Attempt to use reserved field name: " + name);
    }

    {
        const FieldState fstate = FieldState::New;
        const std::string fname =
            field_impl::field_name_with_state(name, fstate);
        const int fid = static_cast<int>(m_int_field_vec.size());

        std::unique_ptr<IntField> field(
            new IntField(*this, fname, fid, ncomp, ngrow, floc));

        if (m_is_initialized) {
            allocate_field_data(*field);
        }

        m_int_field_vec.emplace_back(std::move(field));
        m_int_fid_map[fname] = fid;
    }

    return *m_int_field_vec[m_int_fid_map[name]];
}

IntField&
FieldRepo::get_int_field(const std::string& name, const FieldState fstate) const
{
    BL_PROFILE("amr-wind::FieldRepo::get_int_field");
    AMREX_ALWAYS_ASSERT(fstate == FieldState::New);
    const auto fname = field_impl::field_name_with_state(name, fstate);
    const auto found = m_int_fid_map.find(fname);
    if (found == m_int_fid_map.end()) { // NOLINT(bugprone-branch-clone)
        amrex::Abort("Cannot find field: " + name);
        exit(1); // To appease the compiler
    } else if (found->second < static_cast<unsigned>(m_int_field_vec.size())) {
        return *m_int_field_vec[found->second];
    } else {
        amrex::Abort("Cannot find field: " + name);
        exit(1); // To appease the compiler
    }
}

bool FieldRepo::int_field_exists(
    const std::string& name, const FieldState fstate) const
{
    AMREX_ALWAYS_ASSERT(fstate == FieldState::New);
    const auto fname = field_impl::field_name_with_state(name, fstate);
    const auto found = m_int_fid_map.find(fname);
    return (found != m_int_fid_map.end());
}

std::unique_ptr<ScratchField> FieldRepo::create_scratch_field(
    const std::string& name,
    const int ncomp,
    const int nghost,
    const FieldLoc floc) const
{
    BL_PROFILE("amr-wind::FieldRepo::create_scratch_field");
    if (!m_is_initialized) {
        amrex::Abort(
            "Scratch field creation is not permitted before mesh is "
            "initialized");
    }
    std::unique_ptr<ScratchField> field(
        new ScratchField(*this, name, ncomp, nghost, floc));

    for (int lev = 0; lev <= m_mesh.finestLevel(); ++lev) {
        const auto ba =
            amrex::convert(m_mesh.boxArray(lev), field_impl::index_type(floc));

        field->m_data.emplace_back(
            ba, m_mesh.DistributionMap(lev), ncomp, nghost, amrex::MFInfo(),
            *(m_leveldata[lev]->m_factory));
    }
    return field;
}

std::unique_ptr<ScratchField> FieldRepo::create_scratch_field(
    const int ncomp, const int nghost, const FieldLoc floc) const
{
    return create_scratch_field("scratch_field", ncomp, nghost, floc);
}

std::unique_ptr<ScratchField> FieldRepo::create_scratch_field_on_host(
    const std::string& name,
    const int ncomp,
    const int nghost,
    const FieldLoc floc) const
{
    BL_PROFILE("amr-wind::FieldRepo::create_scratch_field_on_host");
    if (!m_is_initialized) {
        amrex::Abort(
            "Scratch field creation is not permitted before mesh is "
            "initialized");
    }

    std::unique_ptr<ScratchField> field(
        new ScratchField(*this, name, ncomp, nghost, floc));

    for (int lev = 0; lev <= m_mesh.finestLevel(); ++lev) {
        const auto ba =
            amrex::convert(m_mesh.boxArray(lev), field_impl::index_type(floc));

        field->m_data.emplace_back(
            ba, m_mesh.DistributionMap(lev), ncomp, nghost,
            amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()),
            *(m_leveldata[lev]->m_factory));
    }
    return field;
}

std::unique_ptr<ScratchField> FieldRepo::create_scratch_field_on_host(
    const int ncomp, const int nghost, const FieldLoc floc) const
{
    return create_scratch_field_on_host(
        "scratch_field_host", ncomp, nghost, floc);
}
std::unique_ptr<IntScratchField> FieldRepo::create_int_scratch_field_on_host(
    const std::string& name,
    const int ncomp,
    const int nghost,
    const FieldLoc floc) const
{
    BL_PROFILE("amr-wind::FieldRepo::create_int_scratch_field_on_host");
    if (!m_is_initialized) {
        amrex::Abort(
            "Integer scratch field creation is not permitted before mesh is "
            "initialized");
    }

    std::unique_ptr<IntScratchField> field(
        new IntScratchField(*this, name, ncomp, nghost, floc));

    for (int lev = 0; lev <= m_mesh.finestLevel(); ++lev) {
        const auto ba =
            amrex::convert(m_mesh.boxArray(lev), field_impl::index_type(floc));

        field->m_data.emplace_back(
            ba, m_mesh.DistributionMap(lev), ncomp, nghost,
            amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()),
            *(m_leveldata[lev]->m_int_fact));
    }
    return field;
}
std::unique_ptr<IntScratchField> FieldRepo::create_int_scratch_field_on_host(
    const int ncomp, const int nghost, const FieldLoc floc) const
{
    return create_int_scratch_field_on_host(
        "int_scratch_field_host", ncomp, nghost, floc);
}
void FieldRepo::advance_states()
{
    for (auto& it : m_field_vec) {
        if (it->field_state() != FieldState::New) {
            continue;
        }
        it->advance_states();
    }
}

void FieldRepo::allocate_field_data(
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    LevelDataHolder& level_data,
    const amrex::FabFactory<amrex::FArrayBox>& factory)
{
    auto& mfab_vec = level_data.m_mfabs;

    for (auto& field : m_field_vec) {
        auto ba1 =
            amrex::convert(ba, field_impl::index_type(field->field_location()));

        mfab_vec.emplace_back(
            ba1, dm, field->num_comp(), field->num_grow(), amrex::MFInfo(),
            factory);

        mfab_vec.back().setVal(0.0_rt);
    }
}

void FieldRepo::allocate_field_data(
    int lev,
    const Field& field,
    LevelDataHolder& level_data,
    const amrex::FabFactory<amrex::FArrayBox>& factory)
{
    auto& mfab_vec = level_data.m_mfabs;
    AMREX_ASSERT(mfab_vec.size() == field.id());
    const auto ba = amrex::convert(
        m_mesh.boxArray(lev), field_impl::index_type(field.field_location()));

    mfab_vec.emplace_back(
        ba, m_mesh.DistributionMap(lev), field.num_comp(), field.num_grow(),
        amrex::MFInfo(), factory);

    mfab_vec.back().setVal(0.0_rt);
}

void FieldRepo::allocate_field_data(Field& field)
{
    for (int lev = 0; lev <= m_mesh.finestLevel(); ++lev) {
        allocate_field_data(
            lev, field, *m_leveldata[lev], *(m_leveldata[lev]->m_factory));
    }
}

void FieldRepo::allocate_field_data(
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    LevelDataHolder& level_data,
    const amrex::FabFactory<amrex::IArrayBox>& factory)
{
    auto& fab_vec = level_data.m_int_fabs;

    for (auto& field : m_int_field_vec) {
        auto ba1 =
            amrex::convert(ba, field_impl::index_type(field->field_location()));

        fab_vec.emplace_back(
            ba1, dm, field->num_comp(), field->num_grow(), amrex::MFInfo(),
            factory);
    }
}

void FieldRepo::allocate_field_data(
    int lev, const IntField& field, LevelDataHolder& level_data)
{
    auto& fab_vec = level_data.m_int_fabs;
    AMREX_ASSERT(fab_vec.size() == field.id());

    const auto ba = amrex::convert(
        m_mesh.boxArray(lev), field_impl::index_type(field.field_location()));

    fab_vec.emplace_back(
        ba, m_mesh.DistributionMap(lev), field.num_comp(), field.num_grow(),
        amrex::MFInfo(), *level_data.m_int_fact);
}

void FieldRepo::allocate_field_data(const IntField& field)
{
    for (int lev = 0; lev <= m_mesh.finestLevel(); ++lev) {
        allocate_field_data(lev, field, *m_leveldata[lev]);
    }
}

Field& FieldRepo::create_state(Field& infield, const FieldState fstate)
{
    BL_PROFILE("amr-wind::FieldRepo::create_state");
    AMREX_ASSERT((fstate == FieldState::NPH));
    AMREX_ASSERT(!field_exists(infield.base_name(), fstate));

    auto& finfo = infield.m_info;
    const int i = static_cast<int>(fstate);
    const std::string fname =
        field_impl::field_name_with_state(infield.base_name(), fstate);
    const unsigned fid = m_field_vec.size();

    // Create the half state
    std::unique_ptr<Field> field(new Field(*this, fname, finfo, fid, fstate));
    if (m_is_initialized) {
        allocate_field_data(*field);
    }

    // Add reference to states lookup
    finfo->m_states[i] = field.get();
    m_field_vec.emplace_back(std::move(field));
    m_fid_map[fname] = fid;

    return *m_field_vec.back();
}

} // namespace amr_wind
