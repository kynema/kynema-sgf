#include "amr-wind/utilities/sampling/FieldNorms.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"
#include <utility>

using namespace amrex::literals;

namespace amr_wind::field_norms {

FieldNorms::FieldNorms(CFDSim& sim, std::string label)
    : m_sim(sim), m_label(std::move(label))
{}

FieldNorms::~FieldNorms() = default;

void FieldNorms::initialize()
{
    BL_PROFILE("amr-wind::FieldNorms::initialize");

    {
        amrex::ParmParse pp(m_label);
        populate_output_parameters(pp);
        pp.query("mask_redundant_grids", m_use_mask);
        pp.query("use_vector_magnitude", m_use_vector_magnitude);
        std::string norm_type(std::to_string(m_norm_type));
        pp.query("norm_type", norm_type);
        if (norm_type == "2") {
            m_norm_type = 2;
        } else if (norm_type == "1") {
            m_norm_type = 1;
        } else if (amrex::toLower(norm_type) == "infinity") {
            m_norm_type = -1;
        } else {
            amrex::Abort(
                "FieldNorms: norm_type not recognized. Valid options are 1, 2, "
                "and infinity.");
        }
    }

    const auto& io_mng = m_sim.io_manager();
    for (const auto& fld : io_mng.plot_fields()) {
        if (m_use_vector_magnitude) {
            m_var_names.push_back(fld->name());
        } else {
            ioutils::add_var_names(m_var_names, fld->name(), fld->num_comp());
        }
    }

    m_fnorms.resize(m_var_names.size(), 0.0_rt);

    prepare_ascii_file();
}

amrex::Real FieldNorms::get_norm(
    const amr_wind::Field& field,
    const int comp,
    const int ncomp,
    const int norm_type,
    const bool use_mask)
{
    amrex::Real nrm = 0.0_rt;

    AMREX_ASSERT(comp >= 0 && comp < field.num_comp());

    const int finest_level = field.repo().num_active_levels() - 1;
    const auto& geom = field.repo().mesh().Geom();
    constexpr int nghost = 0;

    const auto& mesh = field.repo().mesh();

    for (int lev = 0; lev <= finest_level; lev++) {
        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        auto index_type = field(lev).boxArray().ixType();
        int it_sum = 0;
        int node_dir = 0;
        for (int ix = 0; ix < AMREX_SPACEDIM; ++ix) {
            it_sum += index_type[ix];
            node_dir = index_type[ix] == 1 ? ix : node_dir;
        }

        amrex::MultiFab level_mask;
        if (use_mask) {
            if (lev < finest_level) {
                level_mask = makeFineMask(
                    mesh.boxArray(lev), mesh.DistributionMap(lev),
                    mesh.boxArray(lev + 1), mesh.refRatio(lev), 1.0_rt, 0.0_rt);
            } else {
                level_mask.define(
                    mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                    amrex::MFInfo());
                level_mask.setVal(1.);
            }
        } else {
            // Always on
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1.);
        }

        if (norm_type > 0) {
            nrm += amrex::ReduceSum(
                level_mask, field(lev), nghost,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& mask_arr,
                    amrex::Array4<amrex::Real const> const& field_arr)
                    -> amrex::Real {
                    amrex::Real nrm_fab = 0.0_rt;

                    const auto& fbx =
                        it_sum == 1
                            ? amrex::surroundingNodes(bx, node_dir)
                            : (it_sum > 1 ? amrex::surroundingNodes(bx) : bx);

                    amrex::Loop(fbx, [=, &nrm_fab](int i, int j, int k) {
                        const amrex::IntVect iv{i, j, k};
                        amrex::IntVect iv_cc = iv;
                        // adjust volume for different data locations
                        amrex::Real data_vol = cell_vol;
                        for (int nn = 0; nn < AMREX_SPACEDIM; ++nn) {
                            const bool at_node_dir =
                                index_type[nn] == 1 &&
                                (iv[nn] == fbx.bigEnd(nn) ||
                                 iv[nn] == fbx.smallEnd(nn));
                            data_vol *= at_node_dir ? 0.5_rt : 1.0_rt;
                            // limit mask array indices to cell-centered
                            iv_cc[nn] = amrex::min(
                                bx.bigEnd(nn),
                                amrex::max(bx.smallEnd(nn), iv_cc[nn]));
                        }
                        amrex::Real fval = std::abs(field_arr(i, j, k, comp));
                        // Calculate magnitude if requested
                        for (int n = 1; n < ncomp; ++n) {
                            fval *= fval;
                            fval +=
                                field_arr(i, j, k, n) * field_arr(i, j, k, n);
                            fval = std::sqrt(fval);
                        }
                        fval *= norm_type == 2 ? fval : 1.0_rt;
                        nrm_fab += data_vol * fval * mask_arr(iv_cc);
                    });
                    return nrm_fab;
                });
        } else {
            // Can directly use the field box for the infinity norm
            const amrex::Real nrm_lev = amrex::ReduceMax(
                field(lev), level_mask, nghost,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& field_arr,
                    amrex::Array4<amrex::Real const> const& mask_arr)
                    -> amrex::Real {
                    amrex::Real nrm_fab = 0.0_rt;

                    amrex::Loop(bx, [=, &nrm_fab](int i, int j, int k) {
                        amrex::Real fval = std::abs(field_arr(i, j, k, comp));
                        // Calculate magnitude if requested
                        for (int n = 1; n < ncomp; ++n) {
                            fval *= fval;
                            fval +=
                                field_arr(i, j, k, n) * field_arr(i, j, k, n);
                            fval = std::sqrt(fval);
                        }
                        if (std::abs(mask_arr(i, j, k) - 1.0_rt) <
                            constants::TIGHT_TOL) {
                            nrm_fab = amrex::max<amrex::Real>(nrm_fab, fval);
                        }
                    });
                    return nrm_fab;
                });
            nrm = amrex::max(nrm_lev, nrm);
        }
    }

    if (norm_type > 0) {
        amrex::ParallelDescriptor::ReduceRealSum(nrm);
        const amrex::Real total_volume = geom[0].ProbDomain().volume();
        nrm /= total_volume;
        if (norm_type == 2) {
            nrm = std::sqrt(nrm);
        }
    } else {
        amrex::ParallelDescriptor::ReduceRealMax(nrm);
    }

    return nrm;
}

void FieldNorms::process_field_norms()
{
    int ind = 0;
    for (const auto& fld : m_sim.io_manager().plot_fields()) {
        if (m_use_vector_magnitude) {
            m_fnorms[ind++] =
                get_norm(*fld, 0, fld->num_comp(), m_norm_type, m_use_mask);
        } else {
            for (int comp = 0; comp < fld->num_comp(); ++comp) {
                m_fnorms[ind++] =
                    get_norm(*fld, comp, 1, m_norm_type, m_use_mask);
            }
        }
    }
}

void FieldNorms::output_actions()
{
    BL_PROFILE("amr-wind::FieldNorms::output_actions");

    process_field_norms();
    write_ascii();
}

void FieldNorms::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::FieldNorms::prepare_ascii_file");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str());
        f << "time_step time";
        for (const auto& m_var_name : m_var_names) {
            f << ' ' << m_var_name;
        }
        f << '\n';
        f.close();
    }
}

void FieldNorms::write_ascii()
{
    BL_PROFILE("amr-wind::FieldNorms::write_ascii");

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str(), std::ios_base::app);
        f << m_sim.time().time_index() << std::scientific
          << std::setprecision(m_precision) << std::setw(m_width)
          << m_sim.time().new_time();
        for (amrex::Real m_fnorm : m_fnorms) {
            f << std::setw(m_width) << m_fnorm;
        }
        f << '\n';
        f.close();
    }
}

} // namespace amr_wind::field_norms
