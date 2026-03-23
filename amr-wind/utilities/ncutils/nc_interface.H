#ifndef NC_INTERFACE_H
#define NC_INTERFACE_H

#ifdef AMR_WIND_USE_NETCDF
#include <string>
#include <unordered_map>
#include <vector>

#include "netcdf.h"
#include "netcdf_par.h"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ncutils {

//! Wrapper around NetCDF data types
struct NCDType
{
    static constexpr nc_type Int = NC_INT;
#ifdef AMREX_USE_FLOAT
    static constexpr nc_type Real = NC_FLOAT;
#else
    static constexpr nc_type Real = NC_DOUBLE;
#endif
};

//! Representation of NetCDF dimension
struct NCDim
{
    //! File/Group Identifier
    const int ncid;

    //! Dimension ID used with NetCDF API
    const int dimid;

    //! Name of this dimension
    std::string name() const;

    //! Length of this dimension
    size_t len() const;
};

//! Representation of a NetCDF variable
struct NCVar
{
    //! File/Group identifier
    const int ncid;

    //! Variable ID used with NetCDF API
    const int varid;

    //! Name of this variable
    std::string name() const;

    //! Number of array dimensions for this variable
    int ndim() const;

    //! Shape of the array (size in each array dimension)
    std::vector<size_t> shape() const;

    //! Write out the entire variable
    void put(const amrex::Real* /*ptr*/) const;
    void put(const float* /*ptr*/) const;
    void put(const int* /*ptr*/) const;

    //! Write out a slice of data
    void
    put(const amrex::Real* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/) const;

    //! Write out a slice of data with with strides (see hyperslab definition in
    //! NetCDF)
    void
    put(const amrex::Real* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/,
        const std::vector<ptrdiff_t>& /*stride*/) const;
    //! Write out a slice of data
    void
    put(const float* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/) const;

    //! Write out a slice of data with with strides (see hyperslab definition in
    //! NetCDF)
    void
    put(const float* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/,
        const std::vector<ptrdiff_t>& /*stride*/) const;
    void
    put(const int* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/) const;
    void
    put(const int* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/,
        const std::vector<ptrdiff_t>& /*stride*/) const;

    //! Read the entire variable from file
    void get(amrex::Real* /*ptr*/) const;
    void get(float* /*ptr*/) const;
    void get(int* /*ptr*/) const;

    //! Read a chunk of data from the file
    void
    get(amrex::Real* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/) const;

    //! Read a chunk of data with strides
    void
    get(amrex::Real* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/,
        const std::vector<ptrdiff_t>& /*stride*/) const;

    //! Read a chunk of data from the file
    void
    get(float* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/) const;

    //! Read a chunk of data with strides
    void
    get(float* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/,
        const std::vector<ptrdiff_t>& /*stride*/) const;

    void
    get(int* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/) const;
    void
    get(int* /*dptr*/,
        const std::vector<size_t>& /*start*/,
        const std::vector<size_t>& /*count*/,
        const std::vector<ptrdiff_t>& /*stride*/) const;

    bool has_attr(const std::string& name) const;
    void put_attr(const std::string& name, const std::string& value) const;
    void put_attr(
        const std::string& name, const std::vector<amrex::Real>& value) const;
    void
    put_attr(const std::string& name, const std::vector<float>& value) const;
    void put_attr(const std::string& name, const std::vector<int>& value) const;

    std::string get_attr(const std::string& name) const;
    void
    get_attr(const std::string& name, std::vector<amrex::Real>& value) const;
    void get_attr(const std::string& name, std::vector<float>& value) const;
    void get_attr(const std::string& name, std::vector<int>& value) const;
    void par_access(const int cmode) const;
};

//! Representation of a NetCDF group
class NCGroup
{
public:
    //! Identifier used with NetCDF API calls
    const int ncid;

    //! Is this the NetCDF file (root group)?
    bool is_root() const { return (m_parent != nullptr); }

    //! Return parent group, for root it just returns itself
    const NCGroup& parent() const { return (is_root() ? *this : *m_parent); }

    //! Return the root (file instance)
    const NCGroup& root() const
    {
        return (is_root() ? *this : m_parent->root());
    }

    //! Name of this group
    std::string name() const;

    //! Full name for this group
    std::string full_name() const;

    //! Number of sub-groups within this group
    int num_groups() const;

    //! Number of dimensions in this group
    int num_dimensions() const;

    //! Number of variables within this group
    int num_variables() const;

    //! Number of attributes within this group
    int num_attributes() const;

    //! Check if a group exists
    bool has_group(const std::string& /*name*/) const;

    //! Check if a dimension exists by name
    bool has_dim(const std::string& /*name*/) const;

    //! Check if a variable exists by name
    bool has_var(const std::string& /*name*/) const;

    //! Check if an attribute exists
    bool has_attr(const std::string& /*name*/) const;

    /** Get the group by name
     *
     *  Throws error if the group doesn't exist, use `has_group` to check
     */
    NCGroup group(const std::string& /*name*/) const;

    //! Get the dimension instance by name
    NCDim dim(const std::string& /*name*/) const;

    //! Get the variable instance by name
    NCVar var(const std::string& /*name*/) const;

    //! Define new group
    NCGroup def_group(const std::string& /*name*/) const;

    //! Define new dimension
    NCDim def_dim(const std::string& /*name*/, const size_t len) const;

    //! Define a scalar variable, i.e., 0-dimensional array
    NCVar def_scalar(const std::string& name, const nc_type dtype) const;

    //! Define an array
    NCVar def_array(
        const std::string& name,
        const nc_type dtype,
        const std::vector<std::string>& /*dnames*/) const;

    //! Define a variable (wrapper for def_array)
    NCVar def_var(
        const std::string& name,
        const nc_type dtype,
        const std::vector<std::string>& dnames) const
    {
        return def_array(name, dtype, dnames);
    }

    void put_attr(const std::string& name, const std::string& value) const;
    void put_attr(
        const std::string& name, const std::vector<amrex::Real>& value) const;
    void
    put_attr(const std::string& name, const std::vector<float>& value) const;
    void put_attr(const std::string& name, const std::vector<int>& value) const;

    std::string get_attr(const std::string& name) const;
    void
    get_attr(const std::string& name, std::vector<amrex::Real>& value) const;
    void get_attr(const std::string& name, std::vector<float>& value) const;
    void get_attr(const std::string& name, std::vector<int>& value) const;

    //! Return a list of all groups defined in this group
    std::vector<NCGroup> all_groups() const;

    //! Return a list of all dimensions defined in this group
    std::vector<NCDim> all_dims() const;

    //! Return a list of all variables defined in this group
    std::vector<NCVar> all_vars() const;

    //! Enter definition mode (not needed for NetCDF4 format)
    void enter_def_mode() const;

    //! Exit definition mode
    void exit_def_mode() const;

protected:
    explicit NCGroup(const int id) : ncid(id) {}
    NCGroup(const int id, const NCGroup* par) : ncid(id), m_parent(par) {}

    const NCGroup* m_parent{nullptr};
};

/** Representation of a NetCDF file
 *
 *  Provide wrappes to create and open file
 */
class NCFile : public NCGroup
{
public:
    static NCFile
    create(const std::string& name, const int cmode = NC_CLOBBER | NC_NETCDF4);

    static NCFile open(const std::string& name, const int cmode = NC_NOWRITE);

    static NCFile create_par(
        const std::string& name,
        const int cmode = NC_CLOBBER | NC_NETCDF4 | NC_MPIIO,
        MPI_Comm comm = MPI_COMM_WORLD,
        MPI_Info info = MPI_INFO_NULL);

    static NCFile open_par(
        const std::string& name,
        const int cmode = NC_NOWRITE,
        MPI_Comm comm = MPI_COMM_WORLD,
        MPI_Info info = MPI_INFO_NULL);

    ~NCFile();

    void close();

protected:
    explicit NCFile(const int id) : NCGroup(id), m_is_open{true} {}

    bool m_is_open{false};
};

} // namespace ncutils

#else

namespace ncutils {

struct NCDim
{
    const int ncid{-1};
    const int dimid{-1};
};

struct NCVar
{
    const int ncid{-1};
    const int varid{-1};
};

class NCGroup
{
public:
    const int ncid{-1};
};

class NCFile : public NCGroup
{};

} // namespace ncutils

#endif

#endif /* NC_INTERFACE_H */
