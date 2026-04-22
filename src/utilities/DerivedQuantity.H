#ifndef DERIVEDQUANTITY_H
#define DERIVEDQUANTITY_H

#include <string>
#include <unordered_map>
#include <vector>

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/CollMgr.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

class DerivedQty
    : public Factory<DerivedQty, const FieldRepo&, std::vector<std::string>&>
{
public:
    static std::string base_identifier() { return "DerivedQty"; }

    ~DerivedQty() override = default;

    [[nodiscard]] virtual std::string name() const = 0;

    [[nodiscard]] virtual int num_comp() const = 0;

    virtual void operator()(ScratchField& fld, int scomp = 0) const = 0;

    virtual void var_names(amrex::Vector<std::string>& /*plt_var_names*/);
};

class DerivedQtyMgr
{
public:
    using TypePtr = std::unique_ptr<DerivedQty>;
    using TypeVector = amrex::Vector<TypePtr>;

    explicit DerivedQtyMgr(const FieldRepo& repo);
    ~DerivedQtyMgr() = default;

    void operator()(ScratchField& fld, int scomp = 0) const;

    void create(const amrex::Vector<std::string>& keys);

    //! Add a derived quantity
    DerivedQty& create(const std::string& key);

    //! Return the total number of components across all derived quantities
    [[nodiscard]] int num_comp() const;

    [[nodiscard]] bool contains(const std::string& key) const;

    //! Populate a vector of variable names (for output)
    void var_names(amrex::Vector<std::string>& /*plt_var_names*/) const;

    void filter(const std::set<std::string>& /*erase*/);

    void filter(const amrex::Vector<std::string>& /*erase*/);

private:
    const FieldRepo& m_repo;

    TypeVector m_derived_vec;

    std::unordered_map<std::string, int> m_obj_map;
};

} // namespace amr_wind

#endif /* DERIVEDQUANTITY_H */
