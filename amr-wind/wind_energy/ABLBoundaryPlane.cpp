#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"
#include "amr-wind/wind_energy/ABLFillInflow.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/constants.H"
#include <AMReX_PlotFileUtil.H>
#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

namespace {

//! Return offset vector
AMREX_FORCE_INLINE amrex::IntVect offset(const int face_dir, const int normal)
{
    amrex::IntVect offset(amrex::IntVect::TheDimensionVector(normal));
    if (face_dir == 1) {
        for (auto& o : offset) {
            o *= -1;
        }
    }
    return offset;
}

#ifdef AMR_WIND_USE_NETCDF
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
plane_idx(const int i, const int j, const int k, const int perp, const int lo)
{
    return (static_cast<int>(perp == 0) * i + static_cast<int>(perp == 1) * j +
            static_cast<int>(perp == 2) * k) -
           lo;
}

AMREX_FORCE_INLINE std::string level_name(int lev)
{
    return "level_" + std::to_string(lev);
}
#endif

} // namespace

void InletData::resize(const int size)
{
    m_data_n.resize(size);
    m_data_np1.resize(size);
    m_data_interp.resize(size);
}

void InletData::define_plane(const amrex::Orientation ori)
{
    m_data_n[ori] = std::make_unique<PlaneVector>();
    m_data_np1[ori] = std::make_unique<PlaneVector>();
    m_data_interp[ori] = std::make_unique<PlaneVector>();
}

void InletData::define_level_data(
    const amrex::Orientation ori, const amrex::Box& bx, const size_t nc)
{
    if (!this->is_populated(ori)) {
        return;
    }
    m_data_n[ori]->push_back(amrex::FArrayBox(bx, static_cast<int>(nc)));
    m_data_np1[ori]->push_back(amrex::FArrayBox(bx, static_cast<int>(nc)));
    m_data_interp[ori]->push_back(amrex::FArrayBox(bx, static_cast<int>(nc)));
}

#ifdef AMR_WIND_USE_NETCDF
void InletData::read_data(
    ncutils::NCGroup& grp,
    const amrex::Orientation ori,
    const int lev,
    const Field* fld,
    const amrex::Real time,
    const amrex::Vector<amrex::Real>& times)
{
    const size_t nc = fld->num_comp();
    const int nstart = m_components[static_cast<int>(fld->id())];

    const int idx = utils::closest_index(times, time, constants::LOOSE_TOL);
    const int idxp1 = idx + 1;
    m_tn = times[idx];
    m_tnp1 = times[idxp1];
    if (!((m_tn <= time + constants::LOOSE_TOL) &&
          (time <= m_tnp1 + constants::LOOSE_TOL))) {
        amrex::Abort(
            "ABLBoundaryPlane.cpp InletData::read_data() check failed\n"
            "Left time quantities should be <= right time quantities. Indices "
            "supplied for debugging.\n"
            "m_tn = " +
            std::to_string(m_tn) + ", time + LOOSE_TOL = " +
            std::to_string(time + constants::LOOSE_TOL) +
            "\n"
            "time = " +
            std::to_string(time) + ", m_tnp1 + LOOSE_TOL = " +
            std::to_string(m_tnp1 + constants::LOOSE_TOL) +
            "\n"
            "idx = " +
            std::to_string(idx) + ", idxp1 = " + std::to_string(idxp1));
    }

    const int normal = ori.coordDir();
    const amrex::GpuArray<int, 2> perp = utils::perpendicular_idx(normal);

    const auto& bx = (*m_data_n[ori])[lev].box();
    const auto& lo = bx.loVect();
    const size_t n0 = bx.length(perp[0]);
    const size_t n1 = bx.length(perp[1]);

    // start counting at zero because of netcdf indexing
    amrex::Vector<size_t> start{
        static_cast<size_t>(idx), static_cast<size_t>(0),
        static_cast<size_t>(0), 0};
    amrex::Vector<size_t> count{1, n0, n1, nc};
    amrex::Vector<amrex::Real> buffer(n0 * n1 * nc);
    grp.var(fld->name()).get(buffer.data(), start, count);

    amrex::FArrayBox h_datn(
        bx, (*m_data_n[ori])[lev].nComp(), amrex::The_Pinned_Arena());
    const auto& h_datn_arr = h_datn.array();
    auto* d_buffer = buffer.dataPtr();
    amrex::LoopOnCpu(bx, static_cast<int>(nc), [=](int i, int j, int k, int n) {
        const int i0 = plane_idx(i, j, k, perp[0], lo[perp[0]]);
        const int i1 = plane_idx(i, j, k, perp[1], lo[perp[1]]);
        h_datn_arr(i, j, k, n + nstart) = d_buffer[((i0 * n1) + i1) * nc + n];
    });

    start[0] = static_cast<size_t>(idxp1);
    grp.var(fld->name()).get(buffer.data(), start, count);

    amrex::FArrayBox h_datnp1(
        bx, (*m_data_np1[ori])[lev].nComp(), amrex::The_Pinned_Arena());
    const auto& h_datnp1_arr = h_datnp1.array();
    amrex::LoopOnCpu(bx, static_cast<int>(nc), [=](int i, int j, int k, int n) {
        const int i0 = plane_idx(i, j, k, perp[0], lo[perp[0]]);
        const int i1 = plane_idx(i, j, k, perp[1], lo[perp[1]]);
        h_datnp1_arr(i, j, k, n + nstart) = d_buffer[((i0 * n1) + i1) * nc + n];
    });

    const auto nelems = bx.numPts() * nc;
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, h_datn.dataPtr(nstart),
        h_datn.dataPtr(nstart) + nelems, (*m_data_n[ori])[lev].dataPtr(nstart));
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, h_datnp1.dataPtr(nstart),
        h_datnp1.dataPtr(nstart) + nelems,
        (*m_data_np1[ori])[lev].dataPtr(nstart));
    amrex::Gpu::streamSynchronize();
}

#endif

void InletData::read_data_native(
    const amrex::OrientationIter oit,
    amrex::BndryRegister& bndry_n,
    amrex::BndryRegister& bndry_np1,
    const int lev,
    const Field* fld,
    const amrex::Real time,
    const amrex::Vector<amrex::Real>& times)
{
    const size_t nc = fld->num_comp();
    const int nstart =
        static_cast<int>(m_components[static_cast<int>(fld->id())]);

    const int idx = utils::closest_index(times, time, constants::LOOSE_TOL);
    const int idxp1 = idx + 1;

    m_tn = times[idx];
    m_tnp1 = times[idxp1];

    auto ori = oit();

    if (!(m_tn <= time + constants::LOOSE_TOL) ||
        !(time <= m_tnp1 + constants::LOOSE_TOL)) {
        amrex::Abort(
            "ABLBoundaryPlane.cpp InletData::read_data_native() check "
            "failed\n"
            "Left time quantities should be <= right time quantities. Indices "
            "supplied for debugging.\n"
            "m_tn = " +
            std::to_string(m_tn) + ", time + LOOSE_TOL = " +
            std::to_string(time + constants::LOOSE_TOL) +
            "\n"
            "time = " +
            std::to_string(time) + ", m_tnp1 + LOOSE_TOL = " +
            std::to_string(m_tnp1 + constants::LOOSE_TOL) +
            "\n"
            "idx = " +
            std::to_string(idx) + ", idxp1 = " + std::to_string(idxp1));
    }
    AMREX_ALWAYS_ASSERT(fld->num_comp() == bndry_n[ori].nComp());
    AMREX_ASSERT(bndry_n[ori].boxArray() == bndry_np1[ori].boxArray());

    const int normal = ori.coordDir();
    const auto& bbx = (*m_data_n[ori])[lev].box();
    const amrex::IntVect v_offset = offset(ori.faceDir(), normal);

    amrex::MultiFab bndry(
        bndry_n[ori].boxArray(), bndry_n[ori].DistributionMap(),
        bndry_n[ori].nComp(), 0, amrex::MFInfo());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
    for (amrex::MFIter mfi(bndry); mfi.isValid(); ++mfi) {

        const auto& vbx = mfi.validbox();
        const auto& bndry_n_arr = bndry_n[ori].array(mfi);
        const auto& bndry_arr = bndry.array(mfi);

        const auto& bx = bbx & vbx;
        if (bx.isEmpty()) {
            continue;
        }

        amrex::ParallelFor(
            bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                bndry_arr(i, j, k, n) =
                    0.5_rt *
                    (bndry_n_arr(i, j, k, n) +
                     bndry_n_arr(
                         i + v_offset[0], j + v_offset[1], k + v_offset[2], n));
            });
    }

    bndry.copyTo((*m_data_n[ori])[lev], 0, nstart, static_cast<int>(nc));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
    for (amrex::MFIter mfi(bndry); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& bndry_np1_arr = bndry_np1[ori].array(mfi);
        const auto& bndry_arr = bndry.array(mfi);

        const auto& bx = bbx & vbx;
        if (bx.isEmpty()) {
            continue;
        }

        amrex::ParallelFor(
            bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                bndry_arr(i, j, k, n) =
                    0.5_rt *
                    (bndry_np1_arr(i, j, k, n) +
                     bndry_np1_arr(
                         i + v_offset[0], j + v_offset[1], k + v_offset[2], n));
            });
    }

    bndry.copyTo((*m_data_np1[ori])[lev], 0, nstart, static_cast<int>(nc));
}

void InletData::interpolate(const amrex::Real time)
{
    m_tinterp = time;
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (!this->is_populated(ori)) {
            continue;
        }

        const int lnlevels = static_cast<int>(m_data_n[ori]->size());
        for (int lev = 0; lev < lnlevels; ++lev) {

            const auto& datn = (*m_data_n[ori])[lev];
            const auto& datnp1 = (*m_data_np1[ori])[lev];
            auto& dati = (*m_data_interp[ori])[lev];
            dati.linInterp<amrex::RunOn::Device>(
                datn, 0, datnp1, 0, m_tn, m_tnp1, m_tinterp, datn.box(), 0,
                dati.nComp());
        }
    }
}

bool InletData::is_populated(amrex::Orientation ori) const
{
    return m_data_n[ori] != nullptr;
}

ABLBoundaryPlane::ABLBoundaryPlane(CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_mbc(sim.mbc())
    , m_read_erf(sim.get_read_erf())
{
    amrex::ParmParse pp("ABL");
    int pp_io_mode = -1;
    pp.query("bndry_io_mode", pp_io_mode);
    switch (pp_io_mode) {
    case 0:
        m_io_mode = io_mode::output;
        m_is_initialized = true;
        break;
    case 1:
        m_io_mode = io_mode::input;
        m_is_initialized = true;
        break;
    default:
        m_io_mode = io_mode::undefined;
        m_is_initialized = false;
        return;
    }

    pp.query("bndry_write_frequency", m_write_frequency);
    pp.queryarr("bndry_planes", m_planes);
    pp.query("bndry_output_start_time", m_out_start_time);
    pp.queryarr("bndry_var_names", m_var_names);
    pp.get("bndry_file", m_filename);
    pp.query("bndry_output_format", m_out_fmt);

#ifndef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {
        amrex::Print()
            << "Warning: boundary output format using netcdf must link netcdf "
               "library, changing output to native format"
            << '\n';
        m_out_fmt = "native";
    }
#endif

    if (!((m_out_fmt == "native") || (m_out_fmt == "netcdf") ||
          (m_out_fmt == "erf-multiblock"))) {
        amrex::Print() << "Warning: boundary output format not recognized, "
                          "changing to native format"
                       << '\n';
        m_out_fmt = "native";
    }

    // only used for native format
    m_time_file = m_filename + "/time.dat";
}

void ABLBoundaryPlane::post_init_actions()
{
    if (!m_is_initialized) {
        return;
    }
    initialize_data();
    write_header();
    write_file();
    read_header();
    read_file(false);
}

void ABLBoundaryPlane::pre_advance_work()
{
    if (!m_is_initialized) {
        return;
    }
    read_file(true);
}

void ABLBoundaryPlane::pre_predictor_work()
{
    if (!m_is_initialized) {
        return;
    }
    read_file(false);
}

void ABLBoundaryPlane::post_advance_work()
{
    if (!m_is_initialized) {
        return;
    }
    write_file();
}

void ABLBoundaryPlane::initialize_data()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::initialize_data");
    for (const auto& fname : m_var_names) {
        if (m_repo.field_exists(fname)) {
            auto& fld = m_repo.get_field(fname);
            if (m_io_mode == io_mode::input) {
                fld.register_fill_patch_op<ABLFillInflow>(
                    m_mesh, m_time, *this);
            }
            m_fields.emplace_back(&fld);
        } else {
            amrex::Abort(
                "ABLBoundaryPlane: invalid variable requested: " + fname);
        }
    }
    if ((m_io_mode == io_mode::output) && (m_out_fmt == "erf-multiblock")) {
        amrex::Abort(
            "ABLBoundaryPlane: can't output data in erf-multiblock mode");
    }
}

void ABLBoundaryPlane::write_header()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::write_header");
    if (m_io_mode != io_mode::output) {
        return;
    }

#ifdef AMR_WIND_USE_NETCDF

    if (m_out_fmt == "netcdf") {
        amrex::Print() << "Creating output NetCDF file: " << m_filename
                       << std::endl;

        auto ncf = ncutils::NCFile::create_par(
            m_filename, NC_CLOBBER | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        ncf.enter_def_mode();
        ncf.def_dim("sdim", 1);
        ncf.def_dim("pdim", 2);
        ncf.def_dim("vdim", 3);
        ncf.def_dim("nt", NC_UNLIMITED);
        ncf.def_var("time", NC_DOUBLE, {"nt"});

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const std::string plane = m_plane_names[ori];

            if (std::find(m_planes.begin(), m_planes.end(), plane) ==
                m_planes.end()) {
                continue;
            }

            auto plane_grp = ncf.def_group(plane);

            const int normal = ori.coordDir();
            auto v_normal = plane_grp.def_scalar("normal", NC_INT);
            v_normal.put(&normal);

            const int face_dir = ori.faceDir();
            auto v_face = plane_grp.def_scalar("side", NC_INT);
            v_face.put(&face_dir);

            const auto perp =
                utils::perpendicular_idx<amrex::Vector<int>>(normal);
            auto v_perp = plane_grp.def_var("perpendicular", NC_INT, {"pdim"});
            v_perp.put(perp.data());

            const int nlevels = m_repo.num_active_levels();
            for (int lev = 0; lev < nlevels; ++lev) {

                // Only do this if the output plane intersects with data on this
                // level
                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
                if (!box_intersects_boundary(minBox, lev, ori)) {
                    break;
                }

                auto lev_grp = plane_grp.def_group(level_name(lev));
                lev_grp.def_dim("nx", minBox.length(0));
                lev_grp.def_dim("ny", minBox.length(1));
                lev_grp.def_dim("nz", minBox.length(2));

                lev_grp.def_var("lengths", NC_DOUBLE, {"pdim"});
                lev_grp.def_var("lo", NC_DOUBLE, {"pdim"});
                lev_grp.def_var("hi", NC_DOUBLE, {"pdim"});
                lev_grp.def_var("dx", NC_DOUBLE, {"pdim"});

                const amrex::Vector<std::string> dirs{"nx", "ny", "nz"};
                for (auto* fld : m_fields) {
                    const std::string name = fld->name();
                    if (fld->num_comp() == 1) {
                        lev_grp.def_var(
                            name, NC_DOUBLE,
                            {"nt", dirs[perp[0]], dirs[perp[1]]});
                    } else if (fld->num_comp() == AMREX_SPACEDIM) {
                        lev_grp.def_var(
                            name, NC_DOUBLE,
                            {"nt", dirs[perp[0]], dirs[perp[1]], "vdim"});
                    }
                }
            }
        }
        ncf.put_attr("title", m_title);
        ncf.exit_def_mode();

        // Populate coordinates
        for (auto& plane_grp : ncf.all_groups()) {
            int normal;
            plane_grp.var("normal").get(&normal);
            const amrex::GpuArray<int, 2> perp =
                utils::perpendicular_idx(normal);

            const int nlevels = plane_grp.num_groups();
            for (int lev = 0; lev < nlevels; ++lev) {
                auto lev_grp = plane_grp.group(level_name(lev));

                const auto& dx = m_mesh.Geom(lev).CellSizeArray();
                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
                const auto& lo = minBox.loVect();
                const auto& hi = minBox.hiVect();
                const amrex::Vector<amrex::Real> pdx{dx[perp[0]], dx[perp[1]]};
                const amrex::Vector<amrex::Real> los{
                    lo[perp[0]] * dx[perp[0]], lo[perp[1]] * dx[perp[1]]};
                const amrex::Vector<amrex::Real> his{
                    (hi[perp[0]] + 1) * dx[perp[0]],
                    (hi[perp[1]] + 1) * dx[perp[1]]};
                const amrex::Vector<amrex::Real> lengths{
                    minBox.length(perp[0]) * dx[perp[0]],
                    minBox.length(perp[1]) * dx[perp[1]]};

                lev_grp.var("lengths").put(lengths.data());
                lev_grp.var("lo").put(los.data());
                lev_grp.var("hi").put(his.data());
                lev_grp.var("dx").put(pdx.data());
            }
        }

        amrex::Print() << "NetCDF file created successfully: " << m_filename
                       << std::endl;
    }

#endif

    if (amrex::ParallelDescriptor::IOProcessor() && m_out_fmt == "native") {
        // generate time file
        amrex::UtilCreateCleanDirectory(m_filename, false);
        std::ofstream oftime(m_time_file, std::ios::out);
        oftime.close();
    }
}

void ABLBoundaryPlane::write_bndry_native_header(const std::string& chkname)
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::write_bndry_native_header");
    if (m_io_mode != io_mode::output) {
        return;
    }

#ifndef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {
        amrex::Abort("This is only used in the native format pathway");
    }
#endif

    if (amrex::ParallelDescriptor::IOProcessor()) {
        const amrex::Real time = m_time.new_time();
        const int nlevels = m_repo.num_active_levels();

        for (auto* fld : m_fields) {
            auto& field = *fld;

            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();
                const std::string plane = m_plane_names[ori];

                if (std::ranges::find(m_planes, plane) == m_planes.end()) {
                    continue;
                }

                amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
                const std::string hdr_name(
                    chkname + "/Header_" + std::to_string(ori) + "_" +
                    field.name());
                std::ofstream hdr;
                hdr.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
                hdr.open(
                    hdr_name.c_str(), std::ofstream::out |
                                          std::ofstream::trunc |
                                          std::ofstream::binary);
                if (!hdr.good()) {
                    amrex::FileOpenFailed(hdr_name);
                }

                amrex::Vector<amrex::Geometry> bndry_geoms(nlevels);
                amrex::Vector<amrex::BoxArray> bndry_bas(nlevels);
                amrex::Vector<int> level_steps(nlevels, m_time.time_index());
                const int normal = ori.coordDir();
                for (int lev = 0; lev < nlevels; ++lev) {
                    const auto& geom = field.repo().mesh().Geom(lev);
                    const auto& dx = geom.CellSizeArray();
                    auto bndry_dom = geom.Domain();
                    auto bndry_prob = geom.ProbDomain();
                    if (ori.isLow()) {
                        const int lo = bndry_dom.smallEnd(normal);
                        const auto plo = geom.ProbLo(normal);
                        bndry_dom.setSmall(normal, lo - m_out_rad);
                        bndry_dom.setBig(normal, lo);
                        bndry_prob.setLo(
                            normal, plo - (m_out_rad * dx[normal]));
                        bndry_prob.setHi(normal, plo + dx[normal]);
                    } else {
                        const int hi = bndry_dom.bigEnd(normal);
                        const auto phi = geom.ProbHi(normal);
                        bndry_dom.setSmall(normal, hi);
                        bndry_dom.setBig(normal, hi + m_out_rad);
                        bndry_prob.setLo(normal, phi - dx[normal]);
                        bndry_prob.setHi(
                            normal, phi + (m_out_rad * dx[normal]));
                    }
                    bndry_geoms[lev] = amrex::Geometry(bndry_dom, &bndry_prob);
                    amrex::Box minBox = m_mesh.boxArray(lev).minimalBox();
                    if (ori.isLow()) {
                        minBox.setSmall(normal, bndry_dom.smallEnd(normal));
                    } else {
                        minBox.setBig(normal, bndry_dom.bigEnd(normal));
                    }
                    bndry_bas[lev] = amrex::BoxArray(bndry_dom & minBox);
                }

                amrex::Vector<std::string> var_names;
                ioutils::add_var_names(
                    var_names, field.name(), field.num_comp());
                amrex::WriteGenericPlotfileHeader(
                    hdr, nlevels, bndry_bas, var_names, bndry_geoms, time,
                    level_steps, m_mesh.refRatio(), "HyperCLaw-V1.1", "Level_",
                    field.name() + "_" + std::to_string(ori));
            }
        }
    }
}

void ABLBoundaryPlane::write_file()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::write_file");
    const amrex::Real time = m_time.new_time();
    const int t_step = m_time.time_index();

    // Only output data if at the desired timestep
    if ((t_step % m_write_frequency != 0) || ((m_io_mode != io_mode::output)) ||
        (time < m_out_start_time - constants::LOOSE_TOL)) {
        return;
    }

    for (auto* fld : m_fields) {
        fld->fillpatch(m_time.current_time());
    }

#ifdef AMR_WIND_USE_NETCDF

    if (m_out_fmt == "netcdf") {
        amrex::Print() << "\nWriting NetCDF file " << m_filename << " at time "
                       << time << std::endl;

        auto ncf = ncutils::NCFile::open_par(
            m_filename, NC_WRITE | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        auto v_time = ncf.var("time");
        v_time.par_access(NC_COLLECTIVE);
        const size_t nt = ncf.dim("nt").len();
        v_time.put(&time, {nt}, {1});

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const std::string plane = m_plane_names[ori];

            if (std::find(m_planes.begin(), m_planes.end(), plane) ==
                m_planes.end()) {
                continue;
            }

            const int nlevels = ncf.group(plane).num_groups();
            for (auto* fld : m_fields) {
                for (int lev = 0; lev < nlevels; ++lev) {
                    auto grp = ncf.group(plane).group(level_name(lev));
                    write_data(grp, ori, lev, fld);
                }
            }
        }

        m_out_counter++;
    }

#endif

    if (m_out_fmt == "native") {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream oftime(m_time_file, std::ios::out | std::ios::app);
            oftime << t_step << ' ' << std::setprecision(17) << time << '\n';
            oftime.close();
        }

        const std::string chkname =
            m_filename + amrex::Concatenate("/bndry_output", t_step);

        amrex::Print() << "Writing ABL boundary checkpoint file " << chkname
                       << " at time " << time << '\n';

        const int nlevels = m_repo.num_active_levels();
        const std::string level_prefix = "Level_";
        amrex::PreBuildDirectorHierarchy(chkname, level_prefix, nlevels, true);

        write_bndry_native_header(chkname);

        for (int lev = 0; lev < nlevels; ++lev) {
            for (auto* fld : m_fields) {
                auto& field = *fld;

                const auto& geom = field.repo().mesh().Geom();

                // note: by using one box we end up using 1
                // processor to hold all boundaries
                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
                amrex::BoxArray ba(minBox);
                amrex::DistributionMapping dm{ba};

                amrex::BndryRegister bndry(
                    ba, dm, m_in_rad, m_out_rad, m_extent_rad,
                    field.num_comp());

                bndry.setVal(1.0e13_rt);

                bndry.copyFrom(
                    field(lev), 0, 0, 0, field.num_comp(),
                    geom[lev].periodicity());

                std::string filename = amrex::MultiFabFileFullPrefix(
                    lev, chkname, level_prefix, field.name());

                // print individual faces
                for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                    auto ori = oit();
                    const std::string plane = m_plane_names[ori];

                    if (std::ranges::find(m_planes, plane) == m_planes.end()) {
                        continue;
                    }

                    std::string facename =
                        amrex::Concatenate(filename + '_', ori, 1);
                    bndry[ori].write(facename);
                }
            }
        }
    }
}

void ABLBoundaryPlane::read_header()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::read_header");
    if (m_io_mode != io_mode::input) {
        return;
    }

    // TODO: overallocate this for now
    m_in_data.resize(2 * AMREX_SPACEDIM);

#ifdef AMR_WIND_USE_NETCDF

    if (m_out_fmt == "netcdf") {
        amrex::Print() << "Reading input NetCDF file: " << m_filename
                       << std::endl;
        auto ncf = ncutils::NCFile::open_par(
            m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        // Store the input file times and reset to start at 0
        const size_t nt = ncf.dim("nt").len();
        m_in_times.resize(nt);
        ncf.var("time").get(m_in_times.data());

        // Sanity check the input file time
        if (!(m_in_times[0] <= m_time.current_time() + constants::LOOSE_TOL)) {
            amrex::Abort(
                "ABLBoundaryPlane.cpp ABLBoundaryPlane::read_header() check "
                "failed\n"
                "Left time quantities should be <= right time quantities.\n"
                "m_in_times[0] = " +
                std::to_string(m_in_times[0]) +
                ", current_time + LOOSE_TOL = " +
                std::to_string(m_time.current_time() + constants::LOOSE_TOL));
        }

        for (auto& plane_grp : ncf.all_groups()) {
            int normal, face_dir;
            plane_grp.var("normal").get(&normal);
            plane_grp.var("side").get(&face_dir);
            const amrex::GpuArray<int, 2> perp =
                utils::perpendicular_idx(normal);
            const amrex::Orientation ori(
                normal, amrex::Orientation::Side(face_dir));

            m_in_data.define_plane(ori);

            const int nlevels = plane_grp.num_groups();
            for (int lev = 0; lev < nlevels; ++lev) {
                auto lev_grp = plane_grp.group(level_name(lev));

                // sanity checks to ensure grid-to-grid matching
                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
                const auto& lo = minBox.loVect();
                const auto& hi = minBox.hiVect();
                const auto& dx = m_mesh.Geom(lev).CellSizeArray();
                const amrex::Vector<amrex::Real> pdx{dx[perp[0]], dx[perp[1]]};
                const amrex::Vector<amrex::Real> los{
                    lo[perp[0]] * pdx[0], lo[perp[1]] * pdx[1]};
                const amrex::Vector<amrex::Real> his{
                    (hi[perp[0]] + 1) * pdx[0], (hi[perp[1]] + 1) * pdx[1]};
                const amrex::Vector<amrex::Real> lengths{
                    minBox.length(perp[0]) * pdx[0],
                    minBox.length(perp[1]) * pdx[1]};

                amrex::Vector<amrex::Real> nc_dat{0, 0};
                lev_grp.var("lengths").get(nc_dat.data());
                AMREX_ALWAYS_ASSERT(nc_dat == lengths);
                lev_grp.var("lo").get(nc_dat.data());
                AMREX_ALWAYS_ASSERT(nc_dat == los);
                lev_grp.var("hi").get(nc_dat.data());
                AMREX_ALWAYS_ASSERT(nc_dat == his);
                lev_grp.var("dx").get(nc_dat.data());
                AMREX_ALWAYS_ASSERT(nc_dat == pdx);

                // Create the data structures for the input data
                amrex::IntVect plo(lo);
                amrex::IntVect phi(hi);
                plo[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
                phi[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
                const amrex::Box pbx(plo, phi);
                size_t nc = 0;
                for (auto* fld : m_fields) {
                    m_in_data.component(static_cast<int>(fld->id())) =
                        static_cast<int>(nc);
                    nc += fld->num_comp();
                }
                m_in_data.define_level_data(ori, pbx, nc);
            }
        }

        amrex::Print() << "NetCDF file read successfully: " << m_filename
                       << std::endl;
    }
#endif

    if (m_out_fmt == "native") {

        int time_file_length = 0;

        if (amrex::ParallelDescriptor::IOProcessor()) {

            std::string line;
            std::ifstream time_file(m_time_file);
            if (!time_file.good()) {
                amrex::Abort("Cannot find time file: " + m_time_file);
            }
            while (std::getline(time_file, line)) {
                ++time_file_length;
            }

            time_file.close();
        }

        amrex::ParallelDescriptor::Bcast(
            &time_file_length, 1,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());

        m_in_times.resize(time_file_length);
        m_in_timesteps.resize(time_file_length);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ifstream time_file(m_time_file);
            for (int i = 0; i < time_file_length; ++i) {
                time_file >> m_in_timesteps[i] >> m_in_times[i];
            }
            time_file.close();
        }

        amrex::ParallelDescriptor::Bcast(
            m_in_timesteps.data(), time_file_length,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());

        amrex::ParallelDescriptor::Bcast(
            m_in_times.data(), time_file_length,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());

        int nc = 0;
        for (auto* fld : m_fields) {
            m_in_data.component(static_cast<int>(fld->id())) = nc;
            nc += fld->num_comp();
        }

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();

            if (std::ranges::all_of(m_fields, [ori](const auto* fld) {
                    return (
                        (fld->bc_type()[ori] != BC::mass_inflow) &&
                        (fld->bc_type()[ori] != BC::mass_inflow_outflow));
                })) {
                continue;
            }

            m_in_data.define_plane(ori);

            const int nlevels = boundary_native_file_levels();
            for (int lev = 0; lev < nlevels; ++lev) {

                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();

                amrex::IntVect plo(minBox.loVect());
                amrex::IntVect phi(minBox.hiVect());
                const int normal = ori.coordDir();
                plo[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
                phi[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
                const amrex::Box pbx(plo, phi);
                m_in_data.define_level_data(ori, pbx, nc);
            }
        }
    } else if (m_out_fmt == "erf-multiblock") {

        m_in_times.push_back(-1.0e13_rt); // create space for storing time at
                                          // erf old and new timestep
        m_in_times.push_back(-1.0e13_rt);

        int nc = 0;
        for (auto* fld : m_fields) {
            m_in_data.component(static_cast<int>(fld->id())) = nc;
            nc += fld->num_comp();
        }

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();

            if (std::ranges::all_of(m_fields, [ori](const auto* fld) {
                    return (
                        (fld->bc_type()[ori] != BC::mass_inflow) &&
                        (fld->bc_type()[ori] != BC::mass_inflow_outflow));
                })) {
                continue;
            }

            m_in_data.define_plane(ori);

            // restrict to level 0 for now for multiblock
            const int lev = 0;

            const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();

            amrex::IntVect plo(minBox.loVect());
            amrex::IntVect phi(minBox.hiVect());
            const int normal = ori.coordDir();
            plo[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            phi[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            const amrex::Box pbx(plo, phi);
            m_in_data.define_level_data(ori, pbx, nc);
        }
    }
}

amrex::Vector<amrex::BoxArray> ABLBoundaryPlane::read_bndry_native_boxarrays(
    const std::string& chkname, const Field& field) const
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::read_bndry_native_boxarrays");
    AMREX_ALWAYS_ASSERT(m_io_mode == io_mode::input);

#ifndef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {
        amrex::Abort("This is only used in the native format pathway");
    }
#endif

    const int max_bndry_levels = boundary_native_file_levels();
    amrex::Vector<amrex::BoxArray> bndry_bas(max_bndry_levels);

    // Check if there are any header files
    bool hdr_exists = false;
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();

        if ((field.bc_type()[ori] != BC::mass_inflow) &&
            (field.bc_type()[ori] != BC::mass_inflow_outflow)) {
            continue;
        }

        const std::string hdr_name(
            chkname + "/Header_" + std::to_string(ori) + "_" + field.name());
        if (amrex::FileSystem::Exists(hdr_name)) {
            hdr_exists = true;
            break;
        }
    }

    // If there are no header files, assume the entire domain is in the BoxArray
    if (!hdr_exists) {
        for (int ilev = 0; ilev < max_bndry_levels; ++ilev) {
            amrex::Box domain = m_mesh.Geom(ilev).Domain();
            bndry_bas[ilev] = amrex::BoxArray{domain};
        }
        return bndry_bas;
    }

    amrex::Vector<amrex::Vector<amrex::Box>> bndry_boxes(max_bndry_levels);
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();

        if ((field.bc_type()[ori] != BC::mass_inflow) &&
            (field.bc_type()[ori] != BC::mass_inflow_outflow)) {
            continue;
        }

        const std::string hdr_name(
            chkname + "/Header_" + std::to_string(ori) + "_" + field.name());

        if (!amrex::FileSystem::Exists(hdr_name)) {
            continue;
        }

        amrex::Vector<char> file_char_ptr;
        amrex::ParallelDescriptor::ReadAndBcastFile(hdr_name, file_char_ptr);
        std::string file_char_ptr_string(file_char_ptr.dataPtr());
        std::istringstream is(file_char_ptr_string, std::istringstream::in);

        std::string line;

        // Title line
        is >> line;

        // Number of components
        int ncomp = -1;
        is >> ncomp;
        ioutils::goto_next_line(is);
        AMREX_ALWAYS_ASSERT(ncomp == field.num_comp());

        // Skip the names
        for (int nc = 0; nc < ncomp; nc++) {
            ioutils::goto_next_line(is);
        }

        int spacedim, finest_level;
        amrex::Real time;
        is >> spacedim >> time >> finest_level;
        const int nlevels = finest_level + 1;
        AMREX_ALWAYS_ASSERT(nlevels <= max_bndry_levels);
        AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == spacedim);
        AMREX_ALWAYS_ASSERT(finest_level >= 0);

        amrex::Array<amrex::Real, AMREX_SPACEDIM> prob_lo{
            {AMREX_D_DECL(0., 0., 0.)}};
        amrex::Array<amrex::Real, AMREX_SPACEDIM> prob_hi{
            {AMREX_D_DECL(1., 1., 1.)}};
        amrex::Array<amrex::Real, AMREX_SPACEDIM> prob_size{
            {AMREX_D_DECL(1., 1., 1.)}};
        for (int i = 0; i < spacedim; ++i) {
            is >> prob_lo[i];
        }
        for (int i = 0; i < spacedim; ++i) {
            is >> prob_hi[i];
            prob_size[i] = prob_hi[i] - prob_lo[i];
        }

        const int normal = ori.coordDir();
        const amrex::GpuArray<int, 2> perp = utils::perpendicular_idx(normal);
        AMREX_ALWAYS_ASSERT(
            constants::is_close(
                prob_lo[perp[0]], m_mesh.Geom(0).ProbLo(perp[0])));
        AMREX_ALWAYS_ASSERT(
            constants::is_close(
                prob_lo[perp[1]], m_mesh.Geom(0).ProbLo(perp[1])));
        AMREX_ALWAYS_ASSERT(
            constants::is_close(
                prob_hi[perp[0]], m_mesh.Geom(0).ProbHi(perp[0])));
        AMREX_ALWAYS_ASSERT(
            constants::is_close(
                prob_hi[perp[1]], m_mesh.Geom(0).ProbHi(perp[1])));

        amrex::Vector<int> ref_ratio;
        ref_ratio.resize(nlevels, 0);
        for (int i = 0; i < finest_level; ++i) {
            is >> ref_ratio[i];
        }
        ioutils::goto_next_line(is);

        amrex::Vector<amrex::Box> prob_domain(nlevels);
        for (int i = 0; i < nlevels; ++i) {
            is >> prob_domain[i];
            AMREX_ALWAYS_ASSERT(
                prob_domain[i].smallEnd(perp[0]) ==
                m_mesh.Geom(i).Domain().smallEnd(perp[0]));
            AMREX_ALWAYS_ASSERT(
                prob_domain[i].smallEnd(perp[1]) ==
                m_mesh.Geom(i).Domain().smallEnd(perp[1]));
            AMREX_ALWAYS_ASSERT(
                prob_domain[i].bigEnd(perp[0]) ==
                m_mesh.Geom(i).Domain().bigEnd(perp[0]));
            AMREX_ALWAYS_ASSERT(
                prob_domain[i].bigEnd(perp[1]) ==
                m_mesh.Geom(i).Domain().bigEnd(perp[1]));
        }

        amrex::Vector<int> level_steps(nlevels);
        for (int i = 0; i < nlevels; ++i) {
            is >> level_steps[i];
        }

        amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>> cell_size(
            nlevels, amrex::Array<amrex::Real, AMREX_SPACEDIM>{
                         {AMREX_D_DECL(1., 1., 1.)}});
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            for (int idim = 0; idim < spacedim; ++idim) {
                is >> cell_size[ilev][idim];
            }
        }

        int m_coordsys;
        is >> m_coordsys;
        int bwidth;
        is >> bwidth;

        for (int ilev = 0; ilev < nlevels; ++ilev) {
            int levtmp, ngrids, levsteptmp;
            amrex::Real gtime;
            is >> levtmp >> ngrids >> gtime;
            is >> levsteptmp;
            amrex::Array<amrex::Real, 3> glo = {0.0_rt};
            amrex::Array<amrex::Real, 3> ghi = {0.0_rt};
            AMREX_ASSERT(ngrids == 1);
            for (int igrid = 0; igrid < ngrids; ++igrid) {
                for (int idim = 0; idim < spacedim; ++idim) {
                    is >> glo[idim] >> ghi[idim];
                }
            }
            std::string relname;
            is >> relname;
            std::string mf_name = chkname + "/" + relname;
            const auto vismf = std::make_unique<amrex::VisMF>(mf_name);
            auto ba = vismf->boxArray();
            if (ori.isLow()) {
                ba.growLo(normal, -1);
            } else {
                ba.growHi(normal, -1);
            }
            AMREX_ALWAYS_ASSERT(ba.size() == 1);
            bndry_boxes[ilev].push_back(ba[0]);
        }
    }

    for (int ilev = 0; ilev < bndry_boxes.size(); ilev++) {
        amrex::BoxArray ba(
            bndry_boxes[ilev].data(),
            static_cast<int>(bndry_boxes[ilev].size()));
        bndry_bas[ilev] = amrex::BoxArray(ba.minimalBox());
    }
    return bndry_bas;
}

void ABLBoundaryPlane::read_file(const bool nph_target_time)
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::read_file");
    if (m_io_mode != io_mode::input) {
        return;
    }

    // populate planes and interpolate
    const amrex::Real time =
        nph_target_time ? m_time.current_time() + (0.5_rt * m_time.delta_t())
                        : m_time.new_time();

    if (m_out_fmt == "erf-multiblock") {
        ReadERFFunction read_erf = *m_read_erf;
        if (read_erf != nullptr) {
            read_erf(time, m_in_times, m_in_data, m_fields, mbc());
        } else {
            amrex::Abort("read_erf function is undefined.");
        }
        return;
    }

    if (!(m_in_times[0] <= time + constants::LOOSE_TOL) ||
        !(time < m_in_times.back() + constants::LOOSE_TOL)) {
        amrex::Abort(
            "ABLBoundaryPlane.cpp ABLBoundaryPlane::read_file() check 1"
            "failed\n"
            "Left time quantities should be <= or < right time quantities.\n"
            "m_in_times[0] = " +
            std::to_string(m_in_times[0]) + ", time + LOOSE_TOL = " +
            std::to_string(time + constants::LOOSE_TOL) +
            "\n"
            "time = " +
            std::to_string(time) + ", m_in_times.back() + LOOSE_TOL = " +
            std::to_string(m_in_times.back() + constants::LOOSE_TOL));
    }

    // return early if current data files can still be interpolated in time
    if ((m_in_data.tn() <= time) && (time < m_in_data.tnp1())) {
        m_in_data.interpolate(time);
        return;
    }

#ifdef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {

        auto ncf = ncutils::NCFile::open_par(
            m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            if (!m_in_data.is_populated(ori)) {
                continue;
            }

            const std::string plane = m_plane_names[ori];
            const int nlevels = ncf.group(plane).num_groups();
            for (auto* fld : m_fields) {
                for (int lev = 0; lev < nlevels; ++lev) {
                    auto grp = ncf.group(plane).group(level_name(lev));
                    m_in_data.read_data(grp, ori, lev, fld, time, m_in_times);
                }
            }
        }
    }

#endif

    if (m_out_fmt == "native") {

        const int index =
            utils::closest_index(m_in_times, time, constants::LOOSE_TOL);
        const int t_step1 = m_in_timesteps[index];
        const int t_step2 = m_in_timesteps[index + 1];

        if (!(m_in_times[index] <= time + constants::LOOSE_TOL) ||
            !(time <= m_in_times[index + 1] + constants::LOOSE_TOL)) {
            amrex::Abort(
                "ABLBoundaryPlane.cpp ABLBoundaryPlane::read_file() check 2"
                "failed\n"
                "Left time quantities should be <= right time quantities. "
                "Indices supplied for debugging.\n"
                "m_in_times[index] = " +
                std::to_string(m_in_times[index]) + ", time + LOOSE_TOL = " +
                std::to_string(time + constants::LOOSE_TOL) +
                "\n"
                "time = " +
                std::to_string(time) +
                ", m_in_times[index + 1] + LOOSE_TOL = " +
                std::to_string(m_in_times[index + 1] + constants::LOOSE_TOL) +
                "\n"
                "index = " +
                std::to_string(index) +
                ", index + 1 = " + std::to_string(index + 1));
        }

        const std::string chkname1 =
            m_filename + amrex::Concatenate("/bndry_output", t_step1);
        const std::string chkname2 =
            m_filename + amrex::Concatenate("/bndry_output", t_step2);

        const std::string level_prefix = "Level_";

        const int nlevels = boundary_native_file_levels();
        const auto bndry_bas =
            read_bndry_native_boxarrays(chkname1, *(m_fields[0]));
        for (int lev = 0; lev < nlevels; ++lev) {
            for (auto* fld : m_fields) {
                auto& field = *fld;

                const auto& ba = bndry_bas[lev];
                amrex::DistributionMapping dm{ba};

                amrex::BndryRegister bndry1(
                    ba, dm, m_in_rad, m_out_rad, m_extent_rad,
                    field.num_comp());
                amrex::BndryRegister bndry2(
                    ba, dm, m_in_rad, m_out_rad, m_extent_rad,
                    field.num_comp());

                bndry1.setVal(1.0e13_rt);
                bndry2.setVal(1.0e13_rt);

                std::string filename1 = amrex::MultiFabFileFullPrefix(
                    lev, chkname1, level_prefix, field.name());
                std::string filename2 = amrex::MultiFabFileFullPrefix(
                    lev, chkname2, level_prefix, field.name());

                for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                    auto ori = oit();

                    if ((!m_in_data.is_populated(ori)) ||
                        ((field.bc_type()[ori] != BC::mass_inflow) &&
                         (field.bc_type()[ori] != BC::mass_inflow_outflow))) {
                        continue;
                    }

                    std::string facename1 =
                        amrex::Concatenate(filename1 + '_', ori, 1);
                    std::string facename2 =
                        amrex::Concatenate(filename2 + '_', ori, 1);

                    bndry1[ori].read(facename1);
                    bndry2[ori].read(facename2);

                    m_in_data.read_data_native(
                        oit, bndry1, bndry2, lev, fld, time, m_in_times);
                }
            }
        }
    }

    m_in_data.interpolate(time);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void ABLBoundaryPlane::populate_data(
    const int lev,
    const amrex::Real time,
    Field& fld,
    amrex::MultiFab& mfab,
    const int dcomp,
    const int orig_comp) const
{

    BL_PROFILE("amr-wind::ABLBoundaryPlane::populate_data");

    if (m_io_mode != io_mode::input) {
        return;
    }

    if (!(m_in_data.tn() <= time + constants::LOOSE_TOL) &&
        !(time <= m_in_data.tnp1() + constants::LOOSE_TOL)) {
        amrex::Abort(
            "ABLBoundaryPlane.cpp ABLBoundaryPlane::populate_data() check 1"
            "failed\n"
            "Left time quantities should be <= right time quantities\n"
            "m_in_data.tn() = " +
            std::to_string(m_in_data.tn()) + ", time + LOOSE_TOL = " +
            std::to_string(time + constants::LOOSE_TOL) +
            "\n"
            "time = " +
            std::to_string(time) + ", m_in_data.tnp1() + LOOSE_TOL = " +
            std::to_string(m_in_data.tnp1() + constants::LOOSE_TOL));
    }
    if (!(std::abs(time - m_in_data.tinterp()) < constants::LOOSE_TOL)) {
        amrex::Abort(
            "ABLBoundaryPlane.cpp ABLBoundaryPlane::populate_data() check 2"
            "failed\n"
            "Left time quantities should be < right time quantities. "
            "Additional quantities supplied on second line for debugging.\n"
            "std::abs(time - m_in_data.tinterp()) = " +
            std::to_string(std::abs(time - m_in_data.tinterp())) +
            ", LOOSE_TOL = " + std::to_string(constants::LOOSE_TOL) +
            "\n"
            "time = " +
            std::to_string(time) +
            ", m_in_data.tinterp() = " + std::to_string(m_in_data.tinterp()));
    }

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((!m_in_data.is_populated(ori)) ||
            ((fld.bc_type()[ori] != BC::mass_inflow) &&
             (fld.bc_type()[ori] != BC::mass_inflow_outflow))) {
            continue;
        }

        // Only proceed with data population if fine levels touch the boundary
        if (lev > 0) {
            const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
            if (!box_intersects_boundary(minBox, lev, ori)) {
                continue;
            }
        }

        // Ensure inflow data exists at this level
        if (lev >= m_in_data.nlevels(ori)) {
            amrex::Abort("No inflow data at this level.");
        }

        const size_t nc = mfab.nComp();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {

            auto sbx = mfi.growntilebox(1);
            const auto& src = m_in_data.interpolate_data(ori, lev);
            auto shift_to_cc = amrex::IntVect(0);
            const auto& bx = utils::face_aware_boundary_box_intersection(
                shift_to_cc, sbx, src.box(), ori);
            if (bx.isEmpty()) {
                continue;
            }

            const auto& dest = mfab.array(mfi);
            const auto& src_arr = src.array();
            const int nstart = m_in_data.component(static_cast<int>(fld.id()));
            amrex::ParallelFor(
                bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    dest(i, j, k, n + dcomp) = src_arr(
                        i + shift_to_cc[0], j + shift_to_cc[1],
                        k + shift_to_cc[2], n + nstart + orig_comp);
                });
        }
    }

    const auto& geom = fld.repo().mesh().Geom();
    mfab.EnforcePeriodicity(
        0, mfab.nComp(), amrex::IntVect(1), geom[lev].periodicity());
}

#ifdef AMR_WIND_USE_NETCDF
void ABLBoundaryPlane::write_data(
    const ncutils::NCGroup& grp,
    const amrex::Orientation ori,
    const int lev,
    const Field* fld)
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::write_data");
    // Plane info
    const int normal = ori.coordDir();
    const amrex::GpuArray<int, 2> perp = utils::perpendicular_idx(normal);
    const amrex::IntVect v_offset = offset(ori.faceDir(), normal);

    // Field info
    const std::string name = fld->name();
    const size_t nc = fld->num_comp();

    // Domain info
    const amrex::Box& domain = m_mesh.Geom(lev).Domain();
    const auto& dlo = domain.loVect();
    const auto& dhi = domain.hiVect();

    AMREX_ALWAYS_ASSERT(dlo[0] == 0 && dlo[1] == 0 && dlo[2] == 0);

    grp.var(name).par_access(NC_COLLECTIVE);

    // TODO optimization
    // - move buffer outside this function, probably best as a member
    // - place in object to access as ori/lev/fld
    // - sizing and start/counts should be done only on init and regrid
    const auto n_buffers = m_mesh.boxArray(lev).size();
    amrex::Vector<BufferData> buffers(n_buffers);

    // Compute the minimal offset from the edge of the domain (in case
    // the refinement zones don't coincide with the low edge)
    amrex::IntVect min_lo(std::numeric_limits<int>::max());
    min_lo[normal] = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi((*fld)(lev), false); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.tilebox();
        const auto& blo = bx.loVect();
        const auto& bhi = bx.hiVect();

        if ((blo[normal] == dlo[normal] && ori.isLow()) ||
            (bhi[normal] == dhi[normal] && ori.isHigh())) {
            min_lo[perp[0]] = amrex::min(min_lo[perp[0]], blo[perp[0]]);
            min_lo[perp[1]] = amrex::min(min_lo[perp[1]], blo[perp[1]]);
        }
    }
    amrex::ParallelDescriptor::ReduceIntMin(min_lo.begin(), min_lo.size());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi((*fld)(lev), false); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.tilebox();
        const auto& blo = bx.loVect();
        const auto& bhi = bx.hiVect();

        if (blo[normal] == dlo[normal] && ori.isLow()) {
            amrex::IntVect lo(blo);
            amrex::IntVect hi(bhi);
            lo[normal] = dlo[normal];
            hi[normal] = dlo[normal];
            const amrex::Box lbx(lo, hi);

            const size_t n0 = hi[perp[0]] - lo[perp[0]] + 1;
            const size_t n1 = hi[perp[1]] - lo[perp[1]] + 1;

            auto& buffer = buffers[mfi.index()];
            buffer.data.resize(n0 * n1 * nc);

            auto const& fld_arr = (*fld)(lev).array(mfi);
            impl_buffer_field(
                lbx, static_cast<int>(n1), static_cast<int>(nc), perp, v_offset,
                fld_arr, buffer.data);
            amrex::Gpu::streamSynchronize();

            buffer.start = {
                m_out_counter,
                static_cast<size_t>(lo[perp[0]] - min_lo[perp[0]]),
                static_cast<size_t>(lo[perp[1]] - min_lo[perp[1]]), 0};
            buffer.count = {1, n0, n1, nc};
        } else if (bhi[normal] == dhi[normal] && ori.isHigh()) {
            amrex::IntVect lo(blo);
            amrex::IntVect hi(bhi);
            // shift by one to reuse impl_buffer_field
            lo[normal] = dhi[normal] + 1;
            hi[normal] = dhi[normal] + 1;
            const amrex::Box lbx(lo, hi);

            const size_t n0 = hi[perp[0]] - lo[perp[0]] + 1;
            const size_t n1 = hi[perp[1]] - lo[perp[1]] + 1;

            auto& buffer = buffers[mfi.index()];
            buffer.data.resize(n0 * n1 * nc);

            auto const& fld_arr = (*fld)(lev).array(mfi);
            impl_buffer_field(
                lbx, static_cast<int>(n1), static_cast<int>(nc), perp, v_offset,
                fld_arr, buffer.data);
            amrex::Gpu::streamSynchronize();

            buffer.start = {
                m_out_counter,
                static_cast<size_t>(lo[perp[0]] - min_lo[perp[0]]),
                static_cast<size_t>(lo[perp[1]] - min_lo[perp[1]]), 0};
            buffer.count = {1, n0, n1, nc};
        }
    }

    for (const auto& buffer : buffers) {
        grp.var(name).put(buffer.data.dataPtr(), buffer.start, buffer.count);
    }
}

void ABLBoundaryPlane::impl_buffer_field(
    const amrex::Box& bx,
    const int n1,
    const int nc,
    const amrex::GpuArray<int, 2>& perp,
    const amrex::IntVect& v_offset,
    const amrex::Array4<const amrex::Real>& fld,
    amrex::Gpu::ManagedVector<amrex::Real>& buffer)
{
    auto* d_buffer = buffer.dataPtr();
    const auto lo = bx.loVect3d();
    amrex::ParallelFor(
        bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
            const int i0 = plane_idx(i, j, k, perp[0], lo[perp[0]]);
            const int i1 = plane_idx(i, j, k, perp[1], lo[perp[1]]);
            d_buffer[((i0 * n1) + i1) * nc + n] =
                0.5_rt *
                (fld(i, j, k, n) +
                 fld(i - v_offset[0], j - v_offset[1], k - v_offset[2], n));
        });
}
#endif

// Count the number of levels defined by the native boundary files
int ABLBoundaryPlane::boundary_native_file_levels() const
{
    int nlevels = 0;
    const std::string chkname =
        m_filename + amrex::Concatenate("/bndry_output", m_in_timesteps[0]);
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        const std::string levname = amrex::LevelFullPath(lev, chkname);
        if (amrex::FileExists(levname)) {
            nlevels = lev + 1;
        } else {
            break;
        }
    }
    return nlevels;
}

//! True if box intersects the boundary
bool ABLBoundaryPlane::box_intersects_boundary(
    const amrex::Box& bx, const int lev, const amrex::Orientation ori) const
{
    const amrex::Box& domBox = m_mesh.Geom(lev).Domain();
    const int normal = ori.coordDir();
    amrex::IntVect plo(domBox.loVect());
    amrex::IntVect phi(domBox.hiVect());
    plo[normal] = ori.isHigh() ? domBox.hiVect()[normal] : 0;
    phi[normal] = ori.isHigh() ? domBox.hiVect()[normal] : 0;
    const amrex::Box pbx(plo, phi);
    const auto& intersection = bx & pbx;
    return !intersection.isEmpty();
}

} // namespace amr_wind
