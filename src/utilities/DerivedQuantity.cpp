#include <algorithm>
#include <string>
#include "amr-wind/utilities/DerivedQuantity.H"
#include "amr-wind/utilities/io_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
namespace {

std::string strip_spaces(const std::string& inp)
{
    std::string str(inp);
    str.erase(
        std::remove( // NOLINT(modernize-use-ranges)
            str.begin(), str.end(), ' '),
        str.end());
    return str;
}

std::pair<std::string, std::vector<std::string>>
parse_derived_qty(const std::string& key)
{
    auto popen = key.find('(');
    // If this is not a function type field then return it
    if (popen == std::string::npos) {
        return {key, {}};
    }

    auto pclose = key.find(')');
    if (pclose == std::string::npos) {
        amrex::Abort(
            "Error encountered when parsing derived field name: " + key);
    }

    // Get the function name
    auto dname = strip_spaces(key.substr(0, popen));

    // Collect all arguments required for the derived quantity function
    auto fargs = strip_spaces(key.substr(popen + 1, pclose - popen - 1));
    if (fargs.back() != ',') {
        fargs += ",";
    }
    std::vector<std::string> args;
    size_t start = 0;
    size_t pos;
    while ((pos = fargs.find(',', start)) != std::string::npos) {
        args.push_back(fargs.substr(start, pos - start));
        start = pos + 1;
    }

    return {dname, args};
}

} // namespace

void DerivedQty::var_names(amrex::Vector<std::string>& plt_var_names)
{
    ioutils::add_var_names(plt_var_names, this->name(), this->num_comp());
}

DerivedQtyMgr::DerivedQtyMgr(const FieldRepo& repo) : m_repo(repo) {}

DerivedQty& DerivedQtyMgr::create(const std::string& key)
{
    auto qty_name = strip_spaces(key);

    // If this quantity is already registered return early
    if (contains(qty_name)) {
        return *m_derived_vec[m_obj_map[qty_name]];
    }

    auto tokens = parse_derived_qty(qty_name);
    m_derived_vec.emplace_back(
        DerivedQty::create(tokens.first, m_repo, tokens.second));
    m_obj_map[qty_name] = static_cast<int>(m_derived_vec.size()) - 1;

    return *m_derived_vec.back();
}

void DerivedQtyMgr::create(const amrex::Vector<std::string>& keys)
{
    for (const auto& qty : keys) {
        create(qty);
    }
}

void DerivedQtyMgr::operator()(ScratchField& fld, const int scomp) const
{
    AMREX_ALWAYS_ASSERT((scomp + num_comp()) <= fld.num_comp());

    int icomp = scomp;
    for (const auto& qty : m_derived_vec) {
        (*qty)(fld, icomp);
        icomp += qty->num_comp();
    }
    fld.fillpatch(0.0_rt);
}

int DerivedQtyMgr::num_comp() const
{
    return std::accumulate(
        m_derived_vec.begin(), m_derived_vec.end(), 0,
        [](const int init, const std::unique_ptr<DerivedQty>& qty) {
            return init + qty->num_comp();
        });
}

bool DerivedQtyMgr::contains(const std::string& key) const
{
    auto it = m_obj_map.find(key);
    return (it != m_obj_map.end());
}

void DerivedQtyMgr::var_names(amrex::Vector<std::string>& plt_var_names) const
{
    for (const auto& qty : m_derived_vec) {
        qty->var_names(plt_var_names);
    }
}

void DerivedQtyMgr::filter(const amrex::Vector<std::string>& erase)
{
    const std::set<std::string> set_erase(erase.begin(), erase.end());
    filter(set_erase);
}

void DerivedQtyMgr::filter(const std::set<std::string>& erase)
{
    // first erase from the unordered map
    amrex::Vector<decltype(m_obj_map)::key_type> keys_to_erase;
    for (int i = 0; i < m_derived_vec.size(); i++) {
        if (erase.contains(m_derived_vec[i]->name())) {
            auto it = std::ranges::find_if(
                m_obj_map, [&i](auto&& p) { return p.second == i; });
            keys_to_erase.emplace_back(it->first);
        }
    }
    for (auto&& key : keys_to_erase) {
        m_obj_map.erase(key);
    }

    // then erase from the derived vec
    m_derived_vec.erase(
        std::remove_if( //NOLINT(modernize-use-ranges)
            m_derived_vec.begin(), m_derived_vec.end(),
            [=](const auto& qty) {
                return erase.find(qty->name()) != erase.end();
            }),
        m_derived_vec.end());
}

} // namespace amr_wind
