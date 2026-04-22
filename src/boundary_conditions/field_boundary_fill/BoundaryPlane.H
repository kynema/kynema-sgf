#ifndef BOUNDARYPLANE_H
#define BOUNDARYPLANE_H

#include <cstdint>
#include "src/core/Field.H"
#include "src/CFDSim.H"
#include "src/utilities/ncutils/nc_interface.H"
#include "src/utilities/io_utils.H"
#include "src/boundary_conditions/field_boundary_fill/ReadERFFunction.H"
#include "AMReX_Gpu.H"
#include "AMReX_BndryRegister.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
class MultiBlockContainer;
namespace kynema_sgf {

enum struct io_mode : std::uint8_t { output, input, undefined };

struct BufferData
{
    amrex::Gpu::ManagedVector<amrex::Real> data;
    amrex::Vector<size_t> start{0, 0, 0, 0};
    amrex::Vector<size_t> count{0, 0, 0, 0};
};

/** Collection of data structures and operations for reading data
 *  \ingroup we_abl
 *
 *  This class contains the inlet data structures and operations to
 *  read and interpolate inflow data.
 */
class InletData
{
    using PlaneVector = amrex::Vector<amrex::FArrayBox>;

public:
    InletData() = default;

    void resize(int /*size*/);

    void define_plane(amrex::Orientation /*ori*/);

    void define_level_data(
        amrex::Orientation /*ori*/, const amrex::Box& /*bx*/, size_t /*nc*/);

#ifdef KYNEMA_SGF_USE_NETCDF
    void read_data(
        ncutils::NCGroup& /*grp*/,
        const amrex::Orientation /*ori*/,
        const int /*lev*/,
        const Field* /*fld*/,
        const amrex::Real /*time*/,
        const amrex::Vector<amrex::Real>& /*times*/);
#endif

    void read_data_native(
        amrex::OrientationIter oit,
        amrex::BndryRegister& bndry_n,
        amrex::BndryRegister& bndry_np1,
        int lev,
        const Field* /*fld*/,
        amrex::Real time,
        const amrex::Vector<amrex::Real>& /*times*/);

    void interpolate(amrex::Real /*time*/);
    [[nodiscard]] bool is_populated(amrex::Orientation /*ori*/) const;
    [[nodiscard]] const amrex::FArrayBox&
    interpolate_data(const amrex::Orientation ori, const int lev) const
    {
        return (*m_data_interp[ori])[lev];
    }

    int& component(const int id) { return m_components[id]; }

    [[nodiscard]] int component(const int id) const
    {
        return m_components.at(id);
    }

    [[nodiscard]] int nlevels(const amrex::Orientation ori) const
    {
        return static_cast<int>((*m_data_interp[ori]).size());
    }

    [[nodiscard]] amrex::Real tn() const { return m_tn; }
    [[nodiscard]] amrex::Real tnp1() const { return m_tnp1; }
    [[nodiscard]] amrex::Real tinterp() const { return m_tinterp; }

private:
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_n;
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_np1;
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_interp;

    //! Time for plane at n
    amrex::Real m_tn{-1.0_rt};

    //! Time for plane at n + 1
    amrex::Real m_tnp1{-1.0_rt};

    //! Time for plane at interpolation
    amrex::Real m_tinterp{-1.0_rt};

    //! Map of `{variableId : component}`
    std::unordered_map<int, int> m_components;
};

/** Interface for ABL boundary plane I/O
 *  \ingroup we_abl
 *
 *  This class performs the necessary file operations to read and
 *  write boundary planes.
 *
 *  \sa PlaneFillInflow
 */
class BoundaryPlane : public kynema_sgf::FieldBoundary::Register<BoundaryPlane>
{
    static_assert(AMREX_SPACEDIM == 3, "ABL requires 3 dimensional mesh");

public:
    static std::string identifier() { return "BoundaryPlane"; }

    explicit BoundaryPlane(CFDSim& /*sim*/);

    void post_init_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

    void pre_predictor_work() override;

    void set_velocity(
        int lev,
        amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab,
        int dcomp = 0,
        int orig_comp = 0) const override;

    // Non-inherited public methods

    void pre_advance_inner_calls();

    void initialize_data();

    void write_header();

    void write_bndry_native_header(const std::string& chkname);

    void write_file();

    void read_header();

    [[nodiscard]] amrex::Vector<amrex::BoxArray> read_bndry_native_boxarrays(
        const std::string& chkname, const Field& field) const;

    void read_file(bool /* nph_target_time*/);

    void populate_data(
        int /*lev*/,
        amrex::Real /*time*/,
        const Field& /*fld*/,
        amrex::MultiFab& /*mfab*/,
        int dcomp = 0,
        int orig_comp = 0) const;

#ifdef KYNEMA_SGF_USE_NETCDF
    static void impl_buffer_field(
        const amrex::Box& /*bx*/,
        const int /*n1*/,
        const int /*nc*/,
        const amrex::GpuArray<int, 2>& /*perp*/,
        const amrex::IntVect& /*v_offset*/,
        const amrex::Array4<const amrex::Real>& /*fld*/,
        amrex::Gpu::ManagedVector<amrex::Real>& /*buffer*/);
#endif

    [[nodiscard]] bool box_intersects_boundary(
        const amrex::Box& /*bx*/,
        int /*lev*/,
        amrex::Orientation /*ori*/) const;

    [[nodiscard]] bool is_data_newer_than(const amrex::Real time) const
    {
        return (
            m_in_data.tinterp() - time >
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }

    [[nodiscard]] io_mode mode() const { return m_io_mode; }

    MultiBlockContainer* mbc()
    {
        // Make sure that the double pointer m_mbc
        // does not point to a nullptr
        // before passing *m_mbc for use
        if ((*m_mbc) == nullptr) {
            amrex::Abort("Passing a null pointer for use!");
        }
        return *m_mbc;
    }

private:
    const kynema_sgf::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    // pointer to pointer : when ABL boundary plane gets initialized, kynema-sgf
    // CFDsim does not yet have the up to date pointer to the mbc or the
    // read_erf function, hence the use of double pointers
    MultiBlockContainer** m_mbc;
    ReadERFFunction* m_read_erf{nullptr};

#ifdef KYNEMA_SGF_USE_NETCDF
    void write_data(
        const ncutils::NCGroup& grp,
        const amrex::Orientation /*ori*/,
        const int /*lev*/,
        const Field* /*fld*/);
#endif
    [[nodiscard]] int boundary_native_file_levels() const;

    std::string m_title{"ABL boundary planes"};

    //! Normal direction for the boundary plane
    int m_write_frequency{1};

    //! Plane names
    amrex::Vector<std::string> m_plane_names{"xlo", "ylo", "zlo",
                                             "xhi", "yhi", "zhi"};

    //! IO boundary planes
    amrex::Vector<std::string> m_planes;

    //! Start outputting after this time
    amrex::Real m_out_start_time{0.0_rt};

#ifdef KYNEMA_SGF_USE_NETCDF
    //! NetCDF time output counter
    size_t m_out_counter{0};
#endif

    //! File name for IO
    std::string m_filename;

    //! File name for Native time file
    std::string m_time_file;

    //! Variables for IO
    amrex::Vector<std::string> m_var_names;

    //! List of fields for IO
    amrex::Vector<Field*> m_fields;

    //! NetCDF input level time
    amrex::Vector<amrex::Real> m_in_times;
    amrex::Vector<int> m_in_timesteps;

    //! Inlet data
    InletData m_in_data;

    //! IO mode
    io_mode m_io_mode{io_mode::undefined};

    //! Flag indicating if the simulation is overset -- needed for complicated
    //! timing of pre_advance_work
    bool m_has_overset{false};

    //! controls extents on native bndry output
    const int m_in_rad = 1;
    const int m_out_rad = 1;
    const int m_extent_rad = 0;

    //! output format for bndry output
    std::string m_out_fmt{"native"};
};

} // namespace kynema_sgf

#endif /* BOUNDARYPLANE_H */
