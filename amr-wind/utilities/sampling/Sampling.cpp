#include <algorithm>
#include <memory>
#include <utility>

#include "amr-wind/utilities/sampling/Sampling.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling {

Sampling::Sampling(CFDSim& sim, std::string label)
    : m_sim(sim)
    , m_derived_mgr(new DerivedQtyMgr(m_sim.repo()))
    , m_label(std::move(label))
{}

Sampling::~Sampling() = default;

void Sampling::initialize()
{
    BL_PROFILE("amr-wind::Sampling::initialize");

    // Labels for the different sampler types
    amrex::Vector<std::string> labels;
    // Fields to be sampled - requested by user
    amrex::Vector<std::string> field_names;
    // Int Fields to be sampled - requested by user
    amrex::Vector<std::string> int_field_names;
    // Derived fields to be sampled - requested by user
    amrex::Vector<std::string> derived_field_names;

    {
        amrex::ParmParse pp(m_label);
        pp.getarr("labels", labels);
        ioutils::assert_with_message(
            ioutils::all_distinct(labels),
            "Duplicates in " + m_label + ".labels");
        pp.getarr("fields", field_names);
        ioutils::assert_with_message(
            ioutils::all_distinct(field_names),
            "Duplicates in " + m_label + ".fields");
        pp.queryarr("int_fields", int_field_names);
        ioutils::assert_with_message(
            ioutils::all_distinct(int_field_names),
            "Duplicates in " + m_label + ".int_fields");
        pp.queryarr("derived_fields", derived_field_names);
        pp.query("output_format", m_out_fmt);
        pp.query("restart_sample", m_restart_sample);
        populate_output_parameters(pp);
    }

    // Process field information
    m_ncomp = 0;
    auto& repo = m_sim.repo();
    for (const auto& fname : field_names) {
        if (!repo.field_exists(fname)) {
            amrex::Print()
                << "WARNING: Sampling: Non-existent field requested: " << fname
                << ". This is a mistake or the requested field is a int field "
                   "or a derived "
                   "field and should be added to the int_fields/derived_fields "
                   "parameter"
                << '\n';
            continue;
        }

        auto& fld = repo.get_field(fname);
        m_ncomp += fld.num_comp();
        m_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }

    // Process field information
    m_nicomp = 0;
    for (const auto& fname : int_field_names) {
        if (!repo.int_field_exists(fname)) {
            amrex::Print()
                << "WARNING: Sampling: Non-existent int_field requested: "
                << fname
                << ". This is a mistake or the requested int_field is a "
                   "derived "
                   "field and should be added to the derived_fields parameter"
                << '\n';
            continue;
        }

        auto& fld = repo.get_int_field(fname);
        m_nicomp += fld.num_comp();
        m_int_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }

    // Process derived field information
    if (!derived_field_names.empty()) {
        m_derived_mgr->create(derived_field_names);
        m_derived_mgr->filter(field_names);
        m_ndcomp = m_derived_mgr->num_comp();
        m_derived_mgr->var_names(m_var_names);
    }

    // Load different probe types, default probe type is line
    int idx = 0;
    m_total_particles = 0;
    for (const auto& lbl : labels) {
        const std::string key = m_label + "." + lbl;
        amrex::ParmParse pp1(key);
        std::string stype = "LineSampler";

        pp1.query("type", stype);
        auto obj = SamplerBase::create(stype, m_sim);
        obj->label() = lbl;
        obj->sampletype() = stype;
        obj->id() = idx++;
        obj->initialize(key);

        m_total_particles += obj->num_points();
        m_samplers.emplace_back(std::move(obj));
    }

    update_container();

#ifdef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {
        prepare_netcdf_file();
        m_sample_buf.assign(m_total_particles * m_var_names.size(), 0.0_rt);
    }
#endif

    if (m_restart_sample) {
        sampling_workflow();
        sampling_post();
    }

    // Check
    for (const auto& obj : m_samplers) {
        if ((obj->do_convert_velocity_los()) && (m_out_fmt != "netcdf")) {
            amrex::Abort("Velocity line of sight capability requires NetCDF");
        }
    }
}

void Sampling::update_container()
{
    BL_PROFILE("amr-wind::Sampling::update_container");

    // Initialize the particle container based on user inputs
    m_scontainer = std::make_unique<SamplingContainer>(m_sim.mesh());

    m_scontainer->setup_container(m_ncomp + m_nicomp + m_ndcomp);

    m_scontainer->initialize_particles(m_samplers);

    // Redistribute particles to appropriate boxes/MPI ranks
    m_scontainer->Redistribute();

    m_scontainer->num_sampling_particles() =
        static_cast<int>(m_total_particles);
}

void Sampling::update_sampling_locations()
{
    BL_PROFILE("amr-wind::Sampling::update_sampling_locations");

    amrex::Vector<bool> updated_position;
    for (const auto& obj : m_samplers) {
        const bool updated_pos = obj->update_sampling_locations();
        updated_position.push_back(updated_pos);
    }

    if (std::ranges::any_of(
            updated_position, [](const auto& v) { return v; })) {
        update_container();
    }
}

void Sampling::output_actions()
{
    BL_PROFILE("amr-wind::Sampling::output_actions");

    sampling_workflow();

    process_output();

    sampling_post();
}

void Sampling::sampling_workflow()
{

    BL_PROFILE("amr-wind::Sampling::sampling_workflow");

    update_sampling_locations();

    m_scontainer->interpolate_fields(m_fields, 0);

    m_scontainer->interpolate_fields(m_int_fields, m_ncomp);

    m_scontainer->interpolate_derived_fields(
        *m_derived_mgr, m_sim.repo(), m_ncomp + m_nicomp);

    fill_buffer();

    convert_velocity_lineofsight();

    create_output_buffer();
}

void Sampling::sampling_post()
{

    BL_PROFILE("amr-wind::Sampling::sampling_post");

    for (const auto& obj : m_samplers) {
        obj->post_sample_actions();
    }

#ifdef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {
        m_output_buf.clear();
    }
#endif
}

void Sampling::post_regrid_actions()
{
    BL_PROFILE("amr-wind::Sampling::post_regrid_actions");

    for (const auto& obj : m_samplers) {
        obj->post_regrid_actions();
    }

    m_scontainer->Redistribute();
}

void Sampling::convert_velocity_lineofsight()
{
    BL_PROFILE("amr-wind::Sampling::convert_velocity_lineofsight");

    if (m_out_fmt != "netcdf") {
        return;
    }

#ifdef AMR_WIND_USE_NETCDF
    amrex::Vector<int> vel_map(AMREX_SPACEDIM, 0);
    const amrex::Vector<std::string> vnames = {
        "velocityx", "velocityy", "velocityz"};
    for (int n = 0; n < vnames.size(); n++) {
        auto vit = std::find(m_var_names.begin(), m_var_names.end(), vnames[n]);
        if (vit != m_var_names.end()) {
            vel_map[n] = static_cast<int>(vit - m_var_names.begin());
        } else {
            amrex::Abort("Can't find " + vnames[n]);
        }
    }

    long soffset = 0;
    for (const auto& obj : m_samplers) {
        long sample_size =
            obj->num_points(); // sample locs for individual sampler

        long scan_size =
            (obj->do_subsampling_interp()) ? sample_size / 2 : sample_size;

        std::vector<std::vector<amrex::Real>> temp_vel(
            scan_size, std::vector<amrex::Real>(AMREX_SPACEDIM));
        std::vector<std::vector<amrex::Real>> temp_vel_next(
            scan_size, std::vector<amrex::Real>(AMREX_SPACEDIM));

        if (obj->do_convert_velocity_los()) {
            for (int iv = 0; iv < AMREX_SPACEDIM; ++iv) {
                long vel_off = vel_map[iv];

                long offset =
                    vel_off * m_scontainer->num_sampling_particles() + soffset;
                for (int j = 0; j < scan_size; ++j) {
                    temp_vel[j][iv] = m_sample_buf[offset + j];
                    if (obj->do_subsampling_interp()) {
                        temp_vel_next[j][iv] =
                            m_sample_buf[scan_size + offset + j];
                    }
                }
            }
            obj->calc_lineofsight_velocity(temp_vel, 0);
            if (obj->do_subsampling_interp()) {
                obj->calc_lineofsight_velocity(temp_vel_next, 1);
            }
        }
        soffset += sample_size;
    }
#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void Sampling::create_output_buffer()
{
    BL_PROFILE("amr-wind::Sampling::create_output_buffer");

    if (m_out_fmt != "netcdf") {
        return;
    }

#ifdef AMR_WIND_USE_NETCDF
    const long nvars = m_var_names.size();
    for (int iv = 0; iv < nvars; ++iv) {
        long offset = iv * m_scontainer->num_sampling_particles();
        for (const auto& obj : m_samplers) {
            long sample_size = obj->num_points();
            if (obj->do_data_modification()) {
                // Run data through specific sampler's mod method
                const std::vector<amrex::Real> temp_sb_mod(
                    &m_sample_buf[offset], &m_sample_buf[offset + sample_size]);
                std::vector<amrex::Real> mod_result =
                    obj->modify_sample_data(temp_sb_mod, m_var_names[iv]);
                m_output_buf.insert(
                    m_output_buf.end(), mod_result.begin(), mod_result.end());
                offset += sample_size;
            } else {
                // Directly put m_sample_buf in m_output_buf
                std::vector<amrex::Real> temp_sb(
                    &m_sample_buf[offset], &m_sample_buf[offset + sample_size]);
                m_output_buf.insert(
                    m_output_buf.end(), temp_sb.begin(), temp_sb.end());
                offset += sample_size;
            }
        }
    }

    m_netcdf_output_particles = m_output_buf.size() / nvars;
#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void Sampling::fill_buffer()
{
    BL_PROFILE("amr-wind::Sampling::fill_buffer");
    if (m_out_fmt == "netcdf") {
#ifdef AMR_WIND_USE_NETCDF
        m_scontainer->populate_buffer(m_sample_buf);
#else
        amrex::Abort(
            "NetCDF support was not enabled during build time. Please "
            "recompile or "
            "use native format");
#endif
    }
}

void Sampling::process_output()
{
    BL_PROFILE("amr-wind::Sampling::process_output");
    if (m_out_fmt == "native") {
        impl_write_native();
    } else if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("Sampling: Invalid output format encountered");
    }
}

void Sampling::impl_write_native()
{
    BL_PROFILE("amr-wind::Sampling::write_native");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sampling_name =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    const std::string name(post_dir + "/" + sampling_name);

    m_scontainer->WritePlotFile(
        name, "particles", m_var_names, int_var_names(),
        [=] AMREX_GPU_DEVICE(const SamplingContainer::SuperParticleType& p) {
            return p.id() > 0;
        });

    const std::string info_name = name + "/sampling_info.yaml";
    write_info_file(info_name);

    const std::string header_name = name + "/Header";
    write_header_file(header_name);
}

void Sampling::write_ascii()
{
    BL_PROFILE("amr-wind::Sampling::write_ascii");
    amrex::Print()
        << "WARNING: Sampling: ASCII output will negatively impact performance"
        << '\n';

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    const std::string fname = post_dir + "/" + sname + ".txt";
    m_scontainer->WriteAsciiFile(fname);

    const std::string info_name = post_dir + "/" + sname + "_info.txt";
    write_info_file(info_name);
}

void Sampling::write_info_file(const std::string& fname)
{
    BL_PROFILE("amr-wind::Sampling::write_info_file");

    // Only I/O processor writes the info file
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    if ((m_out_fmt != "native") && (m_out_fmt != "ascii")) {
        amrex::Abort(
            "write_info_file is implemented only for native and ascii formats");
    }

    std::ofstream fh(fname.c_str(), std::ios::out);
    if (!fh.good()) {
        amrex::FileOpenFailed(fname);
    }

    // YAML formatting
    fh << "time: " << std::setprecision(6) << m_sim.time().new_time() << '\n';
    fh << "samplers:" << '\n';
    for (int i = 0; i < m_samplers.size(); ++i) {
        fh << " - index: " << i << '\n';
        fh << "   label: " << m_samplers[i]->label() << '\n';
        fh << "   type: " << m_samplers[i]->sampletype() << '\n';
    }

    fh.close();
}

void Sampling::write_header_file(const std::string& fname)
{
    BL_PROFILE("amr-wind::Sampling::write_header_file");

    // Only I/O processor writes the info file
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    if ((m_out_fmt != "native")) {
        amrex::Abort(
            "write_header_file is implemented only for native formats");
    }

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
    std::ofstream hdr;
    hdr.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    hdr.open(
        fname.c_str(),
        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    if (!hdr.good()) {
        amrex::FileOpenFailed(fname);
    }

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& bas = m_sim.mesh().boxArray();
    const auto& geoms = m_sim.mesh().Geom();
    const auto& ref_ratio = m_sim.mesh().refRatio();
    const amrex::Vector<int> level_steps(nlevels, m_sim.time().time_index());
    auto vnames = m_var_names;
    const auto inames = int_var_names();
    vnames.insert(vnames.end(), inames.begin(), inames.end());
    amrex::WriteGenericPlotfileHeader(
        hdr, nlevels, bas, vnames, geoms, m_sim.time().new_time(), level_steps,
        ref_ratio, "HyperCLaw-V1.1", "Level_", "Cell");
}

void Sampling::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    m_ncfile_name = post_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string npart_name = "num_points";
    const std::vector<std::string> two_dim{nt_name, npart_name};
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind data sampling output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_DOUBLE, {nt_name});

    // Define groups for each sampler
    for (const auto& obj : m_samplers) {
        auto grp = ncf.def_group(obj->label());
        grp.def_dim(npart_name, obj->num_output_points());
        obj->define_netcdf_metadata(grp);
        grp.def_var("coordinates", NC_DOUBLE, {npart_name, "ndim"});

        // Create variables in each sampler
        // Removing velocity components when LOS velocity is output
        for (const std::string& vname : m_var_names) {
            if (!obj->do_convert_velocity_los()) {
                grp.def_var(vname, NC_DOUBLE, two_dim);
            } else {
                if (vname.find("velocity") == std::string::npos) {
                    grp.def_var(vname, NC_DOUBLE, two_dim);
                }
            }
        }

        if (obj->do_convert_velocity_los()) {
            grp.def_var("los_velocity", NC_DOUBLE, two_dim);
        }
    }
    ncf.exit_def_mode();

    {
        const std::vector<size_t> start{0, 0};
        std::vector<size_t> count{0, AMREX_SPACEDIM};
        for (const auto& obj : m_samplers) {
            auto grp = ncf.group(obj->label());
            obj->populate_netcdf_metadata(grp);
            SampleLocType sample_locs;
            obj->output_locations(sample_locs);
            auto xyz = grp.var("coordinates");
            count[0] = obj->num_output_points();
            const auto& locs = sample_locs.locations();
            xyz.put(locs[0].begin(), start, count);
        }
    }

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void Sampling::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    amrex::Print()
        << "WARNING: Sampling: netcdf output will negatively impact performance"
        << std::endl;
    auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of the next timestep
    const size_t nt = ncf.dim(nt_name).len();
    {
        auto time = m_sim.time().new_time();
        ncf.var("time").put(&time, {nt}, {1});
    }

    for (const auto& obj : m_samplers) {
        auto grp = ncf.group(obj->label());
        obj->output_netcdf_data(grp, nt);
    }

    std::vector<size_t> start{nt, 0};
    std::vector<size_t> count{1, 0};

    // Standard sampler output from input deck
    const auto nvars = m_var_names.size();
    for (int iv = 0; iv < nvars; ++iv) {
        std::string vname = m_var_names[iv];
        start[1] = 0;
        count[1] = 0;
        auto offset = iv * num_netcdf_output_particles();
        for (const auto& obj : m_samplers) {
            auto grp = ncf.group(obj->label());
            count[1] = obj->num_output_points();

            if (!obj->do_convert_velocity_los()) {
                auto var = grp.var(vname);
                var.put(&m_output_buf[offset], start, count);
            } else {
                if (vname.find("velocity") == std::string::npos) {
                    auto var = grp.var(vname);
                    var.put(&m_output_buf[offset], start, count);
                }
            }
            offset += static_cast<int>(count[1]);
        }
    }

    // Custom sampler output from sampler function
    // Output of los_velocity goes here in addition to other custom output
    for (const auto& obj : m_samplers) {
        auto grp = ncf.group(obj->label());
        const bool custom_output =
            obj->output_netcdf_field(m_output_buf, grp, nt);
        AMREX_ALWAYS_ASSERT(custom_output);
    }

    ncf.close();
#endif
}

} // namespace amr_wind::sampling
