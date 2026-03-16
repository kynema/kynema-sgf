#include "amr-wind/wind_energy/ABLStats.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLMesoForcingMom.H"
#include "amr-wind/equation_systems/temperature/source_terms/ABLMesoForcingTemp.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/tke/TKE.H"

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_ValLocPair.H"

using namespace amrex::literals;

namespace amr_wind {

ABLStats::ABLStats(
    CFDSim& sim, const ABLWallFunction& abl_wall_func, const int dir)
    : m_sim(sim)
    , m_abl_wall_func(abl_wall_func)
    , m_temperature(sim.repo().get_field("temperature"))
    , m_mueff(sim.pde_manager().icns().fields().mueff)
    , m_pa_vel(sim, dir)
    , m_pa_temp(m_temperature, sim.time(), dir)
    , m_pa_vel_fine(sim, dir)
    , m_pa_temp_fine(m_temperature, sim.time(), dir)
    , m_pa_mueff(m_mueff, sim.time(), dir)
    , m_pa_tt(m_pa_temp, m_pa_temp)
    , m_pa_tu(m_pa_vel, m_pa_temp)
    , m_pa_uu(m_pa_vel, m_pa_vel)
    , m_pa_uuu(m_pa_vel, m_pa_vel, m_pa_vel)
{}

ABLStats::~ABLStats() = default;

void ABLStats::post_init_actions()
{
    initialize();
    calc_averages();
}

void ABLStats::initialize()
{
    BL_PROFILE("amr-wind::ABLStats::initialize");

    {
        amrex::ParmParse pp("ABL");
        pp.query("stats_output_frequency", m_out_freq);
        pp.query("stats_output_format", m_out_fmt);
        pp.query("normal_direction", m_normal_dir);
        AMREX_ASSERT((0 <= m_normal_dir) && (m_normal_dir < AMREX_SPACEDIM));
        pp.query("kappa", m_kappa);
        pp.query("stats_do_energy_budget", m_do_energy_budget);
    }

    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> gravity{0.0_rt, 0.0_rt, -9.81_rt};
        pp.queryarr("gravity", gravity);
        m_gravity = utils::vec_mag(gravity.data());
    }

    // Get normal direction and associated stuff
    const auto& geom = (this->m_sim.repo()).mesh().Geom()[0];
    amrex::Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    if (m_normal_dir == 0) {
        m_ncells_h1 = dhi.y - dlo.y + 1;
        m_ncells_h2 = dhi.z - dlo.z + 1;
    } else if (m_normal_dir == 1) {
        m_ncells_h1 = dhi.x - dlo.x + 1;
        m_ncells_h2 = dhi.z - dlo.z + 1;
    } else if (m_normal_dir == 2) {
        m_ncells_h1 = dhi.x - dlo.x + 1;
        m_ncells_h2 = dhi.y - dlo.y + 1;
    } else {
        amrex::Abort("m_normal_dir is not 0, 1, or 2");
        exit(1); // to help the static analyzer
    }
    m_dn = geom.CellSize()[m_normal_dir];

    if (m_out_freq > 0) {
        if (m_out_fmt == "netcdf") {
            prepare_netcdf_file();
        } else {
            prepare_ascii_file();
        }
    }
}

void ABLStats::calc_averages()
{
    BL_PROFILE("amr-wind::ABLStats::calc_averages");
    m_pa_vel();
    m_pa_temp();
    m_pa_vel_fine();
    m_pa_temp_fine();
    m_pa_mueff();
}

//! Calculate sfs stress averages
void ABLStats::calc_sfs_stress_avgs(
    ScratchField& sfs_stress, ScratchField& t_sfs_stress)
{

    BL_PROFILE("amr-wind::ABLStats::calc_sfs_stress_avgs");

    const auto& repo = m_sim.repo();

    const auto& m_vel = repo.get_field("velocity");
    auto gradVel = repo.create_scratch_field(9);
    fvm::gradient(*gradVel, m_vel);

    const auto& alphaeff = repo.get_field(pde_impl::mueff_name("temperature"));
    auto gradT = repo.create_scratch_field(3);
    fvm::gradient(*gradT, m_temperature);

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& mueff_arrs = m_mueff(lev).const_arrays();
        const auto& alphaeff_arrs = alphaeff(lev).const_arrays();
        const auto& gradVel_arrs = (*gradVel)(lev).const_arrays();
        const auto& gradT_arrs = (*gradT)(lev).const_arrays();
        const auto& sfs_arrs = sfs_stress(lev).arrays();
        const auto& t_sfs_arrs = t_sfs_stress(lev).arrays();
        amrex::ParallelFor(
            m_mueff(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                sfs_arrs[nbx](i, j, k, 0) =
                    -mueff_arrs[nbx](i, j, k) * (gradVel_arrs[nbx](i, j, k, 1) +
                                                 gradVel_arrs[nbx](i, j, k, 3));
                sfs_arrs[nbx](i, j, k, 1) =
                    -mueff_arrs[nbx](i, j, k) * (gradVel_arrs[nbx](i, j, k, 2) +
                                                 gradVel_arrs[nbx](i, j, k, 6));
                sfs_arrs[nbx](i, j, k, 2) =
                    -mueff_arrs[nbx](i, j, k) * (gradVel_arrs[nbx](i, j, k, 5) +
                                                 gradVel_arrs[nbx](i, j, k, 7));

                for (int icomp = 0; icomp < AMREX_SPACEDIM; icomp++) {
                    t_sfs_arrs[nbx](i, j, k, icomp) =
                        -alphaeff_arrs[nbx](i, j, k) *
                        gradT_arrs[nbx](i, j, k, icomp);
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! Calculate sfs stress averages
void ABLStats::calc_tke_diffusion(
    ScratchField& diffusion,
    const Field& buoy_prod,
    const Field& shear_prod,
    const Field& dissipation,
    const amrex::Real dt)
{

    BL_PROFILE("amr-wind::ABLStats::calc_tke_diffusion");

    // Get tke fields
    Field& tke = m_sim.repo().get_field("tke");
    Field& tke_old = tke.state(amr_wind::FieldState::Old);
    // Get conv_term from tke eq
    auto& pde_mgr = m_sim.pde_manager();
    // Check for presence of tke-godunov
    std::string tke_pde_name = amr_wind::pde::TKE::pde_name();
    if (!pde_mgr.has_pde(tke_pde_name)) {
        amrex::Abort(
            "ABL Stats Failure: " + tke_pde_name +
            " PDE not present. Energy budget relies on TKE equation and "
            "Godunov "
            "assumptions.");
    }
    std::string tke_scheme_name = amr_wind::fvm::Godunov::scheme_name();
    if (pde_mgr.scheme() != tke_scheme_name) {
        amrex::Abort(
            "ABL Stats Failure: " + amr_wind::fvm::Godunov::scheme_name() +
            " not being used. Energy budget relies on tke equation and Godunov "
            "assumptions");
    }
    auto& tke_eqn = pde_mgr(tke_pde_name + "-" + tke_scheme_name);
    Field& conv_term = tke_eqn.fields().conv_term;

    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& diffusion_arrs = diffusion(lev).arrays();
        const auto& buoy_prod_arrs = buoy_prod(lev).const_arrays();
        const auto& shear_prod_arrs = shear_prod(lev).const_arrays();
        const auto& dissipation_arrs = dissipation(lev).const_arrays();
        const auto& tke_arrs = tke(lev).const_arrays();
        const auto& tke_old_arrs = tke_old(lev).const_arrays();
        const auto& conv_arrs = conv_term(lev).const_arrays();

        amrex::ParallelFor(
            diffusion(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                diffusion_arrs[nbx](i, j, k) =
                    ((tke_arrs[nbx](i, j, k) - tke_old_arrs[nbx](i, j, k)) /
                     dt) -
                    conv_arrs[nbx](i, j, k) - shear_prod_arrs[nbx](i, j, k) -
                    buoy_prod_arrs[nbx](i, j, k) +
                    dissipation_arrs[nbx](i, j, k);
            });
    }
    amrex::Gpu::streamSynchronize();
}

void ABLStats::post_advance_work()
{
    BL_PROFILE("amr-wind::ABLStats::post_advance_work");

    // Always compute mean velocity/temperature profiles
    calc_averages();

    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) {
        return;
    }

    compute_zi();

    m_pa_tt();
    m_pa_tu();
    m_pa_uu();
    m_pa_uuu();

    process_output();
}

void ABLStats::compute_zi()
{
    BL_PROFILE("amr-wind::ABLStats::compute_zi");

    auto gradT = (this->m_sim.repo())
                     .create_scratch_field(3, m_temperature.num_grow()[0]);
    fvm::gradient(*gradT, m_temperature);

    // Only compute zi using coarsest level
    const int lev = 0;
    const int dir = m_normal_dir;
    const auto& geom = (this->m_sim.repo()).mesh().Geom(lev);
    auto const& domain_box = geom.Domain();

    amrex::Array<bool, AMREX_SPACEDIM> decomp{AMREX_D_DECL(true, true, true)};
    decomp[dir] = false; // no domain decompose in the dir direction.
    auto new_ba = amrex::decompose(
        domain_box, amrex::ParallelDescriptor::NProcs(), decomp);

    amrex::Vector<int> pmap(new_ba.size());
    std::iota(pmap.begin(), pmap.end(), 0);
    amrex::DistributionMapping new_dm(std::move(pmap));

    amrex::MultiFab new_mf(new_ba, new_dm, 1, 0);
    new_mf.ParallelCopy((*gradT)(lev), dir, 0, 1);

    amrex::Real zi_sum = 0;
    int myproc = amrex::ParallelDescriptor::MyProc();
    if (myproc < new_mf.size()) {
        auto const& a = new_mf.const_array(myproc);
        amrex::Box fabbox(a);
#ifdef AMREX_USE_GPU
        amrex::BoxND<AMREX_SPACEDIM - 1> box2d;
        {
            const auto dlo = fabbox.smallEnd();
            const auto dhi = fabbox.bigEnd();
            AMREX_ALWAYS_ASSERT(dir == 2);
            box2d = amrex::BoxND<2>(
                amrex::IntVectND<2>(dlo[0], dlo[1]),
                amrex::IntVectND<2>(dhi[0], dhi[1]));
        }
        const int lendir = domain_box.length(dir);
        const int nblocks = box2d.numPts();
        constexpr int nthreads = 128;
        const amrex::BoxIndexerND<2> box_indexer(box2d);
        amrex::Gpu::DeviceVector<int> tmp(nblocks);
        auto* ptmp = tmp.data();
#ifdef AMREX_USE_SYCL
        constexpr std::size_t shared_mem_bytes =
            sizeof(unsigned long long) * amrex::Gpu::Device::warp_size;
        amrex::launch<nthreads>(
            nblocks, shared_mem_bytes, amrex::Gpu::gpuStream(),
            [=] AMREX_GPU_DEVICE(amrex::Gpu::Handler const& gh) {
                const amrex::Dim1 blockIdx{gh.blockIdx()};
                const amrex::Dim1 threadIdx{gh.threadIdx()};
#else
        amrex::launch<nthreads>(
            nblocks, amrex::Gpu::gpuStream(), [=] AMREX_GPU_DEVICE() {
#endif
                const auto iv2d = box_indexer.intVect(blockIdx.x);
                amrex::IntVect iv;
                int i2d = 0;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    if (idim != dir) {
                        iv[idim] = iv2d[i2d++];
                    }
                }
                amrex::KeyValuePair<amrex::Real, int> r{
                    std::numeric_limits<amrex::Real>::lowest(), 0};
                for (int k = threadIdx.x; k < lendir; k += nthreads) {
                    iv[dir] = k;
                    auto v = a(iv);
                    if (v > r.first()) {
                        r.first() = v;
                        r.second() = k;
                    }
                }
#ifdef AMREX_USE_SYCL
                r = amrex::Gpu::blockReduceMax(r, gh);
#else
                r = amrex::Gpu::blockReduceMax<nthreads>(r);
#endif
                if (threadIdx.x == 0) {
                    ptmp[blockIdx.x] = r.second();
                }
            });

        const auto dnval = m_dn;
        zi_sum = amrex::Reduce::Sum<amrex::Real>(
            nblocks, [=] AMREX_GPU_DEVICE(int iblock) {
                return (ptmp[iblock] + amrex::Real(0.5_rt)) * dnval;
            });
#else
        const auto lo = amrex::lbound(fabbox);
        const auto hi = amrex::ubound(fabbox);
        AMREX_ALWAYS_ASSERT(dir == 2);
#ifdef AMREX_USE_OMP
#pragma omp parallel for collapse(2) reduction(+ : zi_sum)
#endif
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                amrex::Real vmax = std::numeric_limits<amrex::Real>::lowest();
                int idxmax = 0;
                for (int k = lo.z; k <= hi.z; ++k) {
                    if (a(i, j, k) > vmax) {
                        vmax = a(i, j, k);
                        idxmax = i;
                    }
                }
                zi_sum += (idxmax + amrex::Real(0.5_rt)) * m_dn;
            }
        }
#endif
    }

    amrex::ParallelReduce::Sum(
        zi_sum, amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());

    amrex::Long npts = 1;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim != dir) {
            npts *= domain_box.length(idim);
        }
    }
    m_zi = zi_sum / static_cast<amrex::Real>(npts);
}

void ABLStats::process_output()
{
    BL_PROFILE("amr-wind::ABLStats::process_output");

    if (m_out_freq > 0) {
        if (m_out_fmt == "ascii") {
            write_ascii();
        } else if (m_out_fmt == "netcdf") {
            write_netcdf();
        } else {
            amrex::Abort("ABLStats: Invalid output format encountered");
        }
    }
}

void ABLStats::write_ascii()
{
    BL_PROFILE("amr-wind::ABLStats::write_ascii");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const auto& time = m_sim.time();
    m_pa_vel.output_line_average_ascii(
        post_dir + "/plane_average_velocity.txt", time.time_index(),
        time.current_time());
    m_pa_temp.output_line_average_ascii(
        post_dir + "/plane_average_temperature.txt", time.time_index(),
        time.current_time());
    m_pa_vel_fine.output_line_average_ascii(
        post_dir + "/plane_average_velocity_fine.txt", time.time_index(),
        time.current_time());
    m_pa_temp_fine.output_line_average_ascii(
        post_dir + "/plane_average_temperature_fine.txt", time.time_index(),
        time.current_time());
    m_pa_mueff.output_line_average_ascii(
        post_dir + "/plane_average_velocity_mueff.txt", time.time_index(),
        time.current_time());
    m_pa_tt.output_line_average_ascii(
        post_dir + "/second_moment_temperature_temperature.txt",
        time.time_index(), time.current_time());
    m_pa_tu.output_line_average_ascii(
        post_dir + "/second_moment_temperature_velocity.txt", time.time_index(),
        time.current_time());
    m_pa_uu.output_line_average_ascii(
        post_dir + "/second_moment_velocity_velocity.txt", time.time_index(),
        time.current_time());
    m_pa_uuu.output_line_average_ascii(
        post_dir + "/third_moment_velocity_velocity_velocity.txt",
        time.time_index(), time.current_time());

    // Only I/O processor handles this file I/O
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    amrex::RealArray abl_forcing = {0.0_rt};
    if (m_abl_forcing != nullptr) {
        abl_forcing = m_abl_forcing->abl_forcing();
    }

    amrex::Real wstar = 0.0_rt;
    const auto ref_theta = m_sim.transport_model().reference_temperature();
    auto Q = m_abl_wall_func.mo().surf_temp_flux;
    if (Q > std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt) {
        wstar = std::cbrt(m_gravity * Q * m_zi / ref_theta);
    }
    auto L = m_abl_wall_func.mo().obukhov_len;

    std::ofstream outfile;
    outfile.precision(4);
    outfile.open(
        m_ascii_file_name.c_str(), std::ios_base::out | std::ios_base::app);
    // clang-format off
    outfile << time.new_time() << ", "
            << Q << ", "
            << m_abl_wall_func.mo().surf_temp << ", "
            << m_abl_wall_func.utau() << ", "
            << wstar << ", "
            << L << ", "
            << m_zi << ", "
            << abl_forcing[0] << ", "
            << abl_forcing[1] << '\n';
    // clang-format on
    outfile.close();
}

void ABLStats::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::ABLStats::prepare_ascii_file");
    amrex::Print() << "WARNING: ABLStats: ASCII output will impact performance"
                   << '\n';

    // Only I/O processor handles this file I/O
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate("abl_statistics", m_sim.time().time_index());

    m_ascii_file_name = post_dir + "/" + sname + ".txt";

    std::ofstream outfile;
    outfile.open(m_ascii_file_name.c_str(), std::ios_base::out);
    outfile << "Time,   Q, Tsurf, ustar,   wstar,   L,   zi, abl_forcing_x, "
               "abl_forcing_y"
            << '\n';
    outfile.close();
}

void ABLStats::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::ABLStats::prepare_netcdf_file");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate("abl_statistics", m_sim.time().time_index());

    m_ncfile_name = post_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind ABL statistics output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);

    ncf.def_var("time", NC_DOUBLE, {nt_name});
    ncf.def_var("Q", NC_DOUBLE, {nt_name});
    ncf.def_var("Tsurf", NC_DOUBLE, {nt_name});
    ncf.def_var("ustar", NC_DOUBLE, {nt_name});
    ncf.def_var("wstar", NC_DOUBLE, {nt_name});
    ncf.def_var("L", NC_DOUBLE, {nt_name});
    ncf.def_var("zi", NC_DOUBLE, {nt_name});
    ncf.def_var("abl_forcing_x", NC_DOUBLE, {nt_name});
    ncf.def_var("abl_forcing_y", NC_DOUBLE, {nt_name});

    auto grp = ncf.def_group("mean_profiles");
    size_t n_levels = m_pa_vel.ncell_line();
    const std::string nlevels_name = "nlevels";
    grp.def_dim("nlevels", n_levels);
    const std::vector<std::string> two_dim{nt_name, nlevels_name};
    grp.def_var("h", NC_DOUBLE, {nlevels_name});
    grp.def_var("u", NC_DOUBLE, two_dim);
    grp.def_var("v", NC_DOUBLE, two_dim);
    grp.def_var("w", NC_DOUBLE, two_dim);
    grp.def_var("hvelmag", NC_DOUBLE, two_dim);
    grp.def_var("theta", NC_DOUBLE, two_dim);
    amrex::ParmParse pp("ABL");
    if (pp.contains("mesoscale_forcing") || pp.contains("WRFforcing")) {
        grp.def_var("abl_meso_forcing_mom_x", NC_DOUBLE, two_dim);
        grp.def_var("abl_meso_forcing_mom_y", NC_DOUBLE, two_dim);
        grp.def_var("abl_meso_forcing_theta", NC_DOUBLE, two_dim);
    }
    grp.def_var("mueff", NC_DOUBLE, two_dim);
    grp.def_var("theta'theta'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'theta'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'theta'_r", NC_DOUBLE, two_dim);
    grp.def_var("w'theta'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'u'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'v'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'w'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'v'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'w'_r", NC_DOUBLE, two_dim);
    grp.def_var("w'w'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'u'u'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'v'v'_r", NC_DOUBLE, two_dim);
    grp.def_var("w'w'w'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'theta'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("v'theta'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("w'theta'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("u'v'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("u'w'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("v'w'_sfs", NC_DOUBLE, two_dim);
    if (m_sim.repo().field_exists("tke")) {
        grp.def_var("k_sgs", NC_DOUBLE, two_dim);
    }

    // Energy budget
    if (m_do_energy_budget) {
        if (!m_sim.repo().field_exists("buoy_prod")) {
            amrex::Abort(
                "ABL Stats Failure: buoy_prod field not present, indicating "
                "OneEqKsgs turbulence model not being used. Energy budget "
                "currently only applies to this turbulence model.");
        }
        grp.def_var("tke_buoy", NC_DOUBLE, two_dim);
        grp.def_var("tke_shear", NC_DOUBLE, two_dim);
        grp.def_var("tke_dissip", NC_DOUBLE, two_dim);
        grp.def_var("tke_diff", NC_DOUBLE, two_dim);
    }

    ncf.exit_def_mode();

    {
        const std::vector<size_t> start{0};
        std::vector<size_t> count{n_levels};
        auto h = grp.var("h");
        h.put(m_pa_vel.line_centroids().data(), start, count);
    }

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void ABLStats::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::ABLStats::write_netcdf");

    // First calculate sfs stress averages
    auto sfs_stress = m_sim.repo().create_scratch_field("sfs_stress", 3);
    auto t_sfs_stress = m_sim.repo().create_scratch_field("tsfs_stress", 3);
    calc_sfs_stress_avgs(*sfs_stress, *t_sfs_stress);
    ScratchFieldPlaneAveraging pa_sfs(*sfs_stress, m_sim.time(), m_normal_dir);
    pa_sfs();
    ScratchFieldPlaneAveraging pa_tsfs(
        *t_sfs_stress, m_sim.time(), m_normal_dir);
    pa_tsfs();

    if (m_sim.repo().field_exists("tke")) {
        auto& m_ksgs = m_sim.repo().get_field("tke");
        FieldPlaneAveraging pa_ksgs(m_ksgs, m_sim.time(), m_normal_dir);
        pa_ksgs();
    }

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of the next timestep
    const size_t nt = ncf.dim(nt_name).len();
    {
        auto time = m_sim.time().new_time();
        ncf.var("time").put(&time, {nt}, {1});
        auto ustar = m_abl_wall_func.utau();
        ncf.var("ustar").put(&ustar, {nt}, {1});
        amrex::Real wstar = 0.0_rt;
        auto Q = m_abl_wall_func.mo().surf_temp_flux;
        const auto ref_theta = m_sim.transport_model().reference_temperature();
        ncf.var("Q").put(&Q, {nt}, {1});
        auto Tsurf = m_abl_wall_func.mo().surf_temp;
        ncf.var("Tsurf").put(&Tsurf, {nt}, {1});
        if (Q > std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt) {
            wstar = std::cbrt(m_gravity * Q * m_zi / ref_theta);
        }
        ncf.var("wstar").put(&wstar, {nt}, {1});
        amrex::Real L = m_abl_wall_func.mo().obukhov_len;
        ncf.var("L").put(&L, {nt}, {1});
        ncf.var("zi").put(&m_zi, {nt}, {1});

        amrex::RealArray abl_forcing = {0.0_rt};
        if (m_abl_forcing != nullptr) {
            abl_forcing = m_abl_forcing->abl_forcing();
        }
        ncf.var("abl_forcing_x").put(abl_forcing.data(), {nt}, {1});
        ncf.var("abl_forcing_y").put(&abl_forcing[1], {nt}, {1});

        auto grp = ncf.group("mean_profiles");
        size_t n_levels = m_pa_vel.ncell_line();
        amrex::Vector<amrex::Real> l_vec(n_levels);
        std::vector<size_t> start{nt, 0};
        std::vector<size_t> count{1, n_levels};

        {
            amrex::Vector<std::string> var_names{"u", "v", "w"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                m_pa_vel.line_average(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            auto var = grp.var("hvelmag");
            var.put(m_pa_vel.line_hvelmag_average().data(), start, count);
        }

        {
            auto var = grp.var("theta");
            var.put(m_pa_temp.line_average().data(), start, count);
        }

        {
            if (m_abl_meso_mom_forcing != nullptr) {
                amrex::Vector<amrex::Real> meso_forcing_mom_u(n_levels);
                amrex::Vector<amrex::Real> meso_forcing_mom_v(n_levels);

                meso_forcing_mom_u = m_abl_meso_mom_forcing->mom_u_error();
                auto var_x = grp.var("abl_meso_forcing_mom_x");
                var_x.put(meso_forcing_mom_u.data(), start, count);

                meso_forcing_mom_v = m_abl_meso_mom_forcing->mom_v_error();
                auto var_y = grp.var("abl_meso_forcing_mom_y");
                var_y.put(meso_forcing_mom_v.data(), start, count);
            }
            if (m_abl_meso_temp_forcing != nullptr) {
                amrex::Vector<amrex::Real> meso_forcing_theta(n_levels);

                meso_forcing_theta = m_abl_meso_temp_forcing->theta_error();
                auto var_theta = grp.var("abl_meso_forcing_theta");
                var_theta.put(meso_forcing_theta.data(), start, count);
            }
        }

        {
            auto var = grp.var("mueff");
            var.put(m_pa_mueff.line_average().data(), start, count);
        }

        {
            auto var = grp.var("theta'theta'_r");
            m_pa_tt.line_moment(0, l_vec);
            var.put(l_vec.data(), start, count);
        }

        {
            amrex::Vector<std::string> var_names{
                "u'theta'_r", "v'theta'_r", "w'theta'_r"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                m_pa_tu.line_moment(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{"u'u'_r", "u'v'_r", "u'w'_r",
                                                 "v'v'_r", "v'w'_r", "w'w'_r"};
            amrex::Vector<int> var_comp{0, 1, 2, 4, 5, 8};
            for (int i = 0; i < var_comp.size(); i++) {
                m_pa_uu.line_moment(var_comp[i], l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{
                "u'u'u'_r", "v'v'v'_r", "w'w'w'_r"};
            amrex::Vector<int> var_comp{0, 13, 26};
            for (int i = 0; i < var_comp.size(); i++) {
                m_pa_uuu.line_moment(var_comp[i], l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{
                "u'v'_sfs", "u'w'_sfs", "v'w'_sfs"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                pa_sfs.line_average(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{
                "u'theta'_sfs", "v'theta'_sfs", "w'theta'_sfs"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                pa_tsfs.line_average(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        if (m_do_energy_budget) {
            // TKE terms
            auto& tke_buoy_prod = m_sim.repo().get_field("buoy_prod");
            auto& tke_shear_prod = m_sim.repo().get_field("shear_prod");
            auto& tke_dissip = m_sim.repo().get_field("dissipation");
            // Scratch field for diffusion
            auto tke_diffusion =
                m_sim.repo().create_scratch_field("tke_diffusion", 1);
            // Solve for diffusion term using other terms
            calc_tke_diffusion(
                *tke_diffusion, tke_buoy_prod, tke_shear_prod, tke_dissip,
                m_sim.time().delta_t());
            {
                FieldPlaneAveraging pa_tke_buoy_prod(
                    tke_dissip, m_sim.time(), m_normal_dir);
                pa_tke_buoy_prod();
                auto var = grp.var("tke_buoy");
                var.put(pa_tke_buoy_prod.line_average().data(), start, count);
            }
            {
                FieldPlaneAveraging pa_tke_shear_prod(
                    tke_shear_prod, m_sim.time(), m_normal_dir);
                pa_tke_shear_prod();
                auto var = grp.var("tke_shear");
                var.put(pa_tke_shear_prod.line_average().data(), start, count);
            }
            {
                FieldPlaneAveraging pa_tke_dissip(
                    tke_dissip, m_sim.time(), m_normal_dir);
                pa_tke_dissip();
                auto var = grp.var("tke_dissip");
                var.put(pa_tke_dissip.line_average().data(), start, count);
            }
            {
                ScratchFieldPlaneAveraging pa_tke_diff(
                    *tke_diffusion, m_sim.time(), m_normal_dir);
                pa_tke_diff();
                auto var = grp.var("tke_diff");
                var.put(pa_tke_diff.line_average().data(), start, count);
            }
        }
    }

    ncf.close();
#endif
}

} // namespace amr_wind
