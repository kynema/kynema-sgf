#include <cstdio>

#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

#define abort_func amrex::Abort

namespace ncutils {

namespace {

std::array<char, NC_MAX_NAME + 1> recname;

void check_nc_error(int ierr)
{
    if (ierr != NC_NOERR) {
        printf("\n%s\n\n", nc_strerror(ierr));
        abort_func("Encountered NetCDF error; aborting");
    }
}
} // namespace

std::string NCDim::name() const
{
    check_nc_error(nc_inq_dimname(ncid, dimid, recname.begin()));
    return {recname.begin()};
}

size_t NCDim::len() const
{
    size_t dlen;
    check_nc_error(nc_inq_dimlen(ncid, dimid, &dlen));
    return dlen;
}

std::string NCVar::name() const
{
    check_nc_error(nc_inq_varname(ncid, varid, recname.begin()));
    return {recname.begin()};
}

int NCVar::ndim() const
{
    int ndims;
    check_nc_error(nc_inq_varndims(ncid, varid, &ndims));
    return ndims;
}

std::vector<size_t> NCVar::shape() const
{
    int ndims = ndim();
    std::vector<int> dimids(ndims);
    std::vector<size_t> vshape(ndims);

    for (int i = 0; i < ndims; ++i) {
        check_nc_error(nc_inq_vardimid(ncid, varid, dimids.data()));
    }

    for (int i = 0; i < ndims; ++i) {
        check_nc_error(nc_inq_dimlen(ncid, dimids[i], &vshape[i]));
    }

    return vshape;
}

void NCVar::put(const amrex::Real* ptr) const
{
    check_nc_error(nc_put_var_double(ncid, varid, ptr));
}

void NCVar::put(const float* ptr) const
{
    check_nc_error(nc_put_var_float(ncid, varid, ptr));
}

void NCVar::put(const int* ptr) const
{
    check_nc_error(nc_put_var_int(ncid, varid, ptr));
}

void NCVar::put(
    const amrex::Real* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count) const
{
    check_nc_error(
        nc_put_vara_double(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::put(
    const amrex::Real* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count,
    const std::vector<ptrdiff_t>& stride) const
{
    check_nc_error(nc_put_vars_double(
        ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put(
    const float* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count) const
{
    check_nc_error(
        nc_put_vara_float(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::put(
    const float* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count,
    const std::vector<ptrdiff_t>& stride) const
{
    check_nc_error(nc_put_vars_float(
        ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put(
    const int* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count) const
{
    check_nc_error(
        nc_put_vara_int(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::put(
    const int* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count,
    const std::vector<ptrdiff_t>& stride) const
{
    check_nc_error(nc_put_vars_int(
        ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get(amrex::Real* ptr) const
{
    check_nc_error(nc_get_var_double(ncid, varid, ptr));
}

void NCVar::get(float* ptr) const
{
    check_nc_error(nc_get_var_float(ncid, varid, ptr));
}

void NCVar::get(int* ptr) const
{
    check_nc_error(nc_get_var_int(ncid, varid, ptr));
}

void NCVar::get(
    amrex::Real* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count) const
{
    check_nc_error(
        nc_get_vara_double(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get(
    amrex::Real* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count,
    const std::vector<ptrdiff_t>& stride) const
{
    check_nc_error(nc_get_vars_double(
        ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get(
    float* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count) const
{
    check_nc_error(
        nc_get_vara_float(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get(
    float* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count,
    const std::vector<ptrdiff_t>& stride) const
{
    check_nc_error(nc_get_vars_float(
        ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get(
    int* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count) const
{
    check_nc_error(
        nc_get_vara_int(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get(
    int* dptr,
    const std::vector<size_t>& start,
    const std::vector<size_t>& count,
    const std::vector<ptrdiff_t>& stride) const
{
    check_nc_error(nc_get_vars_int(
        ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

bool NCVar::has_attr(const std::string& name) const
{
    int ierr;
    size_t lenp;
    ierr = nc_inq_att(ncid, varid, name.data(), nullptr, &lenp);
    return (ierr == NC_NOERR);
}

void NCVar::put_attr(const std::string& name, const std::string& value) const
{
    check_nc_error(
        nc_put_att_text(ncid, varid, name.data(), value.size(), value.data()));
}

void NCVar::put_attr(
    const std::string& name, const std::vector<amrex::Real>& value) const
{
    check_nc_error(nc_put_att_double(
        ncid, varid, name.data(), NC_DOUBLE, value.size(), value.data()));
}

void NCVar::put_attr(
    const std::string& name, const std::vector<float>& value) const
{
    check_nc_error(nc_put_att_float(
        ncid, varid, name.data(), NC_FLOAT, value.size(), value.data()));
}

void NCVar::put_attr(
    const std::string& name, const std::vector<int>& value) const
{
    check_nc_error(nc_put_att_int(
        ncid, varid, name.data(), NC_INT, value.size(), value.data()));
}

std::string NCVar::get_attr(const std::string& name) const
{
    size_t lenp;
    std::vector<char> aval;
    check_nc_error(nc_inq_attlen(ncid, varid, name.data(), &lenp));
    aval.resize(lenp);
    check_nc_error(nc_get_att_text(ncid, varid, name.data(), aval.data()));
    return std::string{aval.begin(), aval.end()};
}

void NCVar::get_attr(
    const std::string& name, std::vector<amrex::Real>& values) const
{
    size_t lenp;
    check_nc_error(nc_inq_attlen(ncid, varid, name.data(), &lenp));
    values.resize(lenp);
    check_nc_error(nc_get_att_double(ncid, varid, name.data(), values.data()));
}

void NCVar::get_attr(const std::string& name, std::vector<float>& values) const
{
    size_t lenp;
    check_nc_error(nc_inq_attlen(ncid, varid, name.data(), &lenp));
    values.resize(lenp);
    check_nc_error(nc_get_att_float(ncid, varid, name.data(), values.data()));
}

void NCVar::get_attr(const std::string& name, std::vector<int>& values) const
{
    size_t lenp;
    check_nc_error(nc_inq_attlen(ncid, varid, name.data(), &lenp));
    values.resize(lenp);
    check_nc_error(nc_get_att_int(ncid, varid, name.data(), values.data()));
}

void NCVar::par_access(const int cmode) const
{
    check_nc_error(nc_var_par_access(ncid, varid, cmode));
}

std::string NCGroup::name() const
{
    size_t nlen;
    std::vector<char> grpname;
    check_nc_error(nc_inq_grpname_len(ncid, &nlen));
    grpname.resize(nlen + 1);
    check_nc_error(nc_inq_grpname(ncid, grpname.data()));
    return std::string{grpname.begin(), grpname.end()};
}

std::string NCGroup::full_name() const
{
    size_t nlen;
    std::vector<char> grpname;
    check_nc_error(nc_inq_grpname_full(ncid, &nlen, nullptr));
    grpname.resize(nlen);
    check_nc_error(nc_inq_grpname_full(ncid, &nlen, grpname.data()));
    return std::string{grpname.begin(), grpname.end()};
}

NCGroup NCGroup::def_group(const std::string& name) const
{
    int newid;
    check_nc_error(nc_def_grp(ncid, name.data(), &newid));
    return {newid, this};
}

NCGroup NCGroup::group(const std::string& name) const
{
    int newid;
    check_nc_error(nc_inq_ncid(ncid, name.data(), &newid));
    return {newid, this};
}

NCDim NCGroup::dim(const std::string& name) const
{
    int newid;
    check_nc_error(nc_inq_dimid(ncid, name.data(), &newid));
    return NCDim{ncid, newid};
}

NCDim NCGroup::def_dim(const std::string& name, const size_t len) const
{
    int newid;
    check_nc_error(nc_def_dim(ncid, name.data(), len, &newid));
    return NCDim{ncid, newid};
}

NCVar NCGroup::def_scalar(const std::string& name, const nc_type dtype) const
{
    int newid;
    check_nc_error(nc_def_var(ncid, name.data(), dtype, 0, nullptr, &newid));
    return NCVar{ncid, newid};
}

NCVar NCGroup::def_array(
    const std::string& name,
    const nc_type dtype,
    const std::vector<std::string>& dnames) const
{
    int newid;
    auto ndims = dnames.size();
    std::vector<int> dimids(ndims);
    for (int i = 0; i < static_cast<int>(ndims); ++i) {
        dimids[i] = dim(dnames[i]).dimid;
    }

    check_nc_error(nc_def_var(
        ncid, name.data(), dtype, static_cast<int>(ndims), dimids.data(),
        &newid));
    return NCVar{ncid, newid};
}

NCVar NCGroup::var(const std::string& name) const
{
    int varid;
    check_nc_error(nc_inq_varid(ncid, name.data(), &varid));
    return NCVar{ncid, varid};
}

int NCGroup::num_groups() const
{
    int ngrps;
    check_nc_error(nc_inq_grps(ncid, &ngrps, nullptr));
    return ngrps;
}

int NCGroup::num_dimensions() const
{
    int ndims;
    check_nc_error(nc_inq(ncid, &ndims, nullptr, nullptr, nullptr));
    return ndims;
}

int NCGroup::num_attributes() const
{
    int nattrs;
    check_nc_error(nc_inq(ncid, nullptr, nullptr, &nattrs, nullptr));
    return nattrs;
}

int NCGroup::num_variables() const
{
    int nvars;
    check_nc_error(nc_inq(ncid, nullptr, &nvars, nullptr, nullptr));
    return nvars;
}

bool NCGroup::has_group(const std::string& name) const
{
    int ierr = nc_inq_ncid(ncid, name.data(), nullptr);
    return (ierr == NC_NOERR);
}

bool NCGroup::has_dim(const std::string& name) const
{
    int ierr = nc_inq_dimid(ncid, name.data(), nullptr);
    return (ierr == NC_NOERR);
}

bool NCGroup::has_var(const std::string& name) const
{
    int rh_id;
    int ierr = nc_inq_varid(ncid, name.data(), &rh_id);
    return (ierr == NC_NOERR);
}

bool NCGroup::has_attr(const std::string& name) const
{
    int ierr;
    size_t lenp;
    ierr = nc_inq_att(ncid, NC_GLOBAL, name.data(), nullptr, &lenp);
    return (ierr == NC_NOERR);
}

void NCGroup::put_attr(const std::string& name, const std::string& value) const
{
    check_nc_error(nc_put_att_text(
        ncid, NC_GLOBAL, name.data(), value.size(), value.data()));
}

void NCGroup::put_attr(
    const std::string& name, const std::vector<amrex::Real>& value) const
{
    check_nc_error(nc_put_att_double(
        ncid, NC_GLOBAL, name.data(), NC_DOUBLE, value.size(), value.data()));
}

void NCGroup::put_attr(
    const std::string& name, const std::vector<float>& value) const
{
    check_nc_error(nc_put_att_float(
        ncid, NC_GLOBAL, name.data(), NC_FLOAT, value.size(), value.data()));
}

void NCGroup::put_attr(
    const std::string& name, const std::vector<int>& value) const
{
    check_nc_error(nc_put_att_int(
        ncid, NC_GLOBAL, name.data(), NC_INT, value.size(), value.data()));
}

std::string NCGroup::get_attr(const std::string& name) const
{
    size_t lenp;
    std::vector<char> aval;
    check_nc_error(nc_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    aval.resize(lenp);
    check_nc_error(nc_get_att_text(ncid, NC_GLOBAL, name.data(), aval.data()));
    return std::string{aval.begin(), aval.end()};
}

void NCGroup::get_attr(
    const std::string& name, std::vector<amrex::Real>& values) const
{
    size_t lenp;
    check_nc_error(nc_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    values.resize(lenp);
    check_nc_error(
        nc_get_att_double(ncid, NC_GLOBAL, name.data(), values.data()));
}

void NCGroup::get_attr(
    const std::string& name, std::vector<float>& values) const
{
    size_t lenp;
    check_nc_error(nc_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    values.resize(lenp);
    check_nc_error(
        nc_get_att_float(ncid, NC_GLOBAL, name.data(), values.data()));
}

void NCGroup::get_attr(const std::string& name, std::vector<int>& values) const
{
    size_t lenp;
    check_nc_error(nc_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    values.resize(lenp);
    check_nc_error(nc_get_att_int(ncid, NC_GLOBAL, name.data(), values.data()));
}

std::vector<NCGroup> NCGroup::all_groups() const
{
    std::vector<NCGroup> grps;
    int ngrps = num_groups();

    // Empty list of groups return early without error
    if (ngrps < 1) {
        return grps;
    }

    std::vector<int> gids(ngrps);
    check_nc_error(nc_inq_grps(ncid, &ngrps, gids.data()));
    grps.reserve(ngrps);
    for (int i = 0; i < ngrps; ++i) {
        grps.emplace_back(NCGroup(gids[i], this));
    }
    return grps;
}

std::vector<NCDim> NCGroup::all_dims() const
{
    std::vector<NCDim> adims;
    int ndims = num_dimensions();
    adims.reserve(ndims);
    for (int i = 0; i < ndims; ++i) {
        adims.emplace_back(NCDim{ncid, i});
    }
    return adims;
}

std::vector<NCVar> NCGroup::all_vars() const
{
    std::vector<NCVar> avars;
    int nvars = num_variables();
    avars.reserve(nvars);
    for (int i = 0; i < nvars; ++i) {
        avars.emplace_back(NCVar{ncid, i});
    }
    return avars;
}

void NCGroup::enter_def_mode() const
{
    int ierr;
    ierr = nc_redef(ncid);

    // Ignore already in define mode error
    if (ierr == NC_EINDEFINE) {
        return;
    }
    // Handle all other errors
    check_nc_error(ierr);
}

void NCGroup::exit_def_mode() const { check_nc_error(nc_enddef(ncid)); }

NCFile NCFile::create(const std::string& name, const int cmode)
{
    int ncid;
    check_nc_error(nc_create(name.data(), cmode, &ncid));
    return NCFile(ncid);
}

NCFile NCFile::open(const std::string& name, const int cmode)
{
    int ncid;
    check_nc_error(nc_open(name.data(), cmode, &ncid));
    return NCFile(ncid);
}

NCFile NCFile::create_par(
    const std::string& name, const int cmode, MPI_Comm comm, MPI_Info info)
{
    int ncid;
    check_nc_error(nc_create_par(name.data(), cmode, comm, info, &ncid));
    return NCFile(ncid);
}

NCFile NCFile::open_par(
    const std::string& name, const int cmode, MPI_Comm comm, MPI_Info info)
{
    int ncid;
    check_nc_error(nc_open_par(name.data(), cmode, comm, info, &ncid));
    return NCFile(ncid);
}

NCFile::~NCFile()
{
    if (m_is_open) {
        check_nc_error(nc_close(ncid));
    }
}

void NCFile::close()
{
    m_is_open = false;
    check_nc_error(nc_close(ncid));
}

} // namespace ncutils
