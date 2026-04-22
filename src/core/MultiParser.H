#ifndef MULTIPARSER_H
#define MULTIPARSER_H

#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::utils {

/** Utility to parse inputs from multiple namespaces.
 *
 *  This class is a thin-wrapper around AMReX ParmParse utility to query user
 *  inputs from two different namespaces (default and specialized). This allows
 *  the user to consolidate common options within a default section and then
 *  override specific values in the specialization section.
 *
 */
class MultiParser
{
public:
    MultiParser(const std::string& default_prefix, const std::string& prefix)
        : m_pp_default(default_prefix), m_pp(prefix)
    {}

    //! Return the ParmParse instance for the default namespace
    [[nodiscard]] const amrex::ParmParse& default_params() const
    {
        return m_pp_default;
    }

    //! Return the ParmParse instance for the specialized namespace
    [[nodiscard]] const amrex::ParmParse& params() const { return m_pp; }

    //! Check if the keyword is present in either namespace
    [[nodiscard]] bool contains(const std::string& name) const
    {
        return m_pp.contains(name) || m_pp_default.contains(name);
    }

    /** Get a vector from the input file
     *
     *  \param [in] name Keyword to search for
     *  \param [out] value Vector value
     */
    void get(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        getarr(name, val);
        AMREX_ALWAYS_ASSERT(val.size() == AMREX_SPACEDIM);
        value.x() = val[0];
        value.y() = val[1];
        value.z() = val[2];
    }

    /** Query and return if a vector exists in the input file
     *
     *  \param [in] name Keyword to search for
     *  \param [out] value Vector value
     */
    void query(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        queryarr(name, val);
        if (!val.empty()) {
            AMREX_ALWAYS_ASSERT(val.size() == AMREX_SPACEDIM);
            value.x() = val[0];
            value.y() = val[1];
            value.z() = val[2];
        }
    }

    /** Get a vector from the input file, but could be specified uniform instead
     *
     *  \param [in] name Keyword to search for
     *  \param [out] value Vector value
     */
    void get_either(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        getarr(name, val);
        AMREX_ALWAYS_ASSERT(val.size() == AMREX_SPACEDIM || val.size() == 1);
        if (val.size() == 1) {
            value.x() = val[0];
            value.y() = val[0];
            value.z() = val[0];
        } else {
            value.x() = val[0];
            value.y() = val[1];
            value.z() = val[2];
        }
    }

    /** Query and return if a scalar (uniform) or vector exists in input file
     *
     *  \param [in] name Keyword to search for
     *  \param [out] value Vector value
     */
    void query_either(const std::string& name, vs::Vector& value) const
    {
        amrex::Vector<vs::Vector::value_type> val;
        queryarr(name, val);
        if (!val.empty()) {
            AMREX_ALWAYS_ASSERT(
                val.size() == AMREX_SPACEDIM || val.size() == 1);
            if (val.size() == 1) {
                value.x() = val[0];
                value.y() = val[0];
                value.z() = val[0];
            } else {
                value.x() = val[0];
                value.y() = val[1];
                value.z() = val[2];
            }
        }
    }

    //! Get the value for the keyword entry from either namespace
    template <typename T>
    void get(const std::string& name, T& value) const
    {
        if (m_pp.contains(name)) {
            m_pp.get(name, value);
        } else {
            m_pp_default.get(name, value);
        }
    }

    //! Get a vector of values for the given keyword entry from either namespace
    template <typename T>
    void getarr(const std::string& name, T& value) const
    {
        if (m_pp.contains(name)) {
            m_pp.getarr(name, value);
        } else {
            m_pp_default.getarr(name, value);
        }
    }

    //! Query the value for the keyword entry from either namespace
    template <typename T>
    void query(const std::string& name, T& value) const
    {
        m_pp_default.query(name, value);
        m_pp.query(name, value);
    }

    //! Query a vector of values for the given keyword entry from either
    //! namespace
    template <typename T>
    void queryarr(const std::string& name, T& value) const
    {
        m_pp_default.queryarr(name, value);
        m_pp.queryarr(name, value);
    }

private:
    amrex::ParmParse m_pp_default;
    amrex::ParmParse m_pp;
};

} // namespace amr_wind::utils

#endif /* MULTIPARSER_H */
