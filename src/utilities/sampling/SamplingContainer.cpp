
#include <algorithm>

#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/core/Field.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling {

void SamplingContainer::setup_container(
    const int num_real_components, const int num_int_components)
{
    BL_PROFILE("amr-wind::SamplingContainer::setup");
    const bool communicate_comp = true;
    for (int i = 0; i < num_real_components; ++i) {
        AddRealComp(static_cast<int>(communicate_comp));
    }
    for (int i = 0; i < num_int_components; ++i) {
        AddIntComp(static_cast<int>(communicate_comp));
    }

    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
        for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
            DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
        }
    }
}

void SamplingContainer::initialize_particles(
    const amrex::Vector<std::unique_ptr<SamplerBase>>& samplers)
{
    BL_PROFILE("amr-wind::SamplingContainer::initialize");

    const int lev = 0;
    const auto iproc = amrex::ParallelDescriptor::MyProc();

    m_total_particles = 0;
    for (const auto& probes : samplers) {
        m_total_particles += probes->num_points();
    }

    const int nprobes = static_cast<int>(samplers.size());
#ifdef AMREX_DEBUG
    const auto& dxinv = m_mesh.Geom(lev).InvCellSizeArray();
    const auto& plo = m_mesh.Geom(lev).ProbLoArray();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        int np_box = 0;
        amrex::Vector<SampleLocType> vec_sample_locs(nprobes);
        for (int iprobe = 0; iprobe < nprobes; iprobe++) {
            const auto& probe = samplers[iprobe];

            SampleLocType sample_locs;
            probe->sampling_locations(sample_locs, {box});
            const auto& locs = sample_locs.locations();
            const int npts = static_cast<int>(locs.size());
            np_box += static_cast<int>(npts);
            vec_sample_locs[iprobe] = sample_locs;
        }

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& ptile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
        ptile.resize(np_box);

        if (np_box == 0) {
            continue;
        }

        int offset = 0;
        for (int iprobe = 0; iprobe < nprobes; iprobe++) {
            const auto& probe = samplers[iprobe];
            const auto& sample_locs = vec_sample_locs[iprobe];
            const auto& locs = sample_locs.locations();
            const int npts = static_cast<int>(locs.size());
            if (npts == 0) {
                continue;
            }

            const int uid_offset = std::accumulate(
                samplers.begin(), samplers.begin() + iprobe, 0,
                [&](int sum, const auto& s) { return sum + s->num_points(); });
            const auto probe_id = probe->id();
            amrex::Gpu::DeviceVector<amrex::RealVect> dlocs(npts);
            amrex::Gpu::copy(
                amrex::Gpu::hostToDevice, locs.begin(), locs.end(),
                dlocs.begin());
            const auto* p_dlocs = dlocs.data();
            const auto& ids = sample_locs.ids();
            amrex::Gpu::DeviceVector<amrex::Long> dids(npts);
            amrex::Gpu::copy(
                amrex::Gpu::hostToDevice, ids.begin(), ids.end(), dids.begin());
            const auto* p_dids = dids.data();

            const amrex::Long nextid = ParticleType::NextID();
            ParticleType::NextID(nextid + npts);

            auto* pstruct = ptile.GetArrayOfStructs()().data();
            amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE(const int ip) {
                const amrex::RealVect loc(AMREX_D_DECL(
                    p_dlocs[ip][0], p_dlocs[ip][1], p_dlocs[ip][2]));
#ifdef AMREX_DEBUG
                const amrex::IntVect div(AMREX_D_DECL(
                    static_cast<int>(
                        amrex::Math::floor((loc[0] - plo[0]) * dxinv[0])),
                    static_cast<int>(
                        amrex::Math::floor((loc[1] - plo[1]) * dxinv[1])),
                    static_cast<int>(
                        amrex::Math::floor((loc[2] - plo[2]) * dxinv[2]))));
                AMREX_ASSERT(box.contains(div));
#endif
                auto& pp = pstruct[offset + ip];
                pp.id() = nextid + ip;
                pp.cpu() = iproc;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    pp.pos(idim) = loc[idim];
                }
                pp.idata(IIx::uid) = static_cast<int>(p_dids[ip] + uid_offset);
                pp.idata(IIx::sid) = probe_id;
                pp.idata(IIx::nid) = static_cast<int>(p_dids[ip]);
            });
            offset += npts;
        }
    }
    // Skip this check if there is a DTUSpinnerSampler (may have out of domain
    // particles)
    if (std::ranges::all_of(samplers, [](const auto& probe) {
            return probe->sampletype() != "DTUSpinnerSampler";
        })) {
        AMREX_ALWAYS_ASSERT(m_total_particles == TotalNumberOfParticles(false));
    }
}

void SamplingContainer::interpolate_derived_fields(
    const DerivedQtyMgr& derived_mgr, const FieldRepo& repo, const int scomp)
{
    BL_PROFILE("amr-wind::SamplingContainer::interpolate_derived_fields");

    auto outfield = repo.create_scratch_field(derived_mgr.num_comp(), 1);
    derived_mgr(*outfield, 0);

    const int nlevels = m_mesh.finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const auto farr = (*outfield)(lev).const_array(pti);
            interpolate(
                pti, farr, lev, outfield->field_location(),
                outfield->num_comp(), scomp);
        }
    }
}

void SamplingContainer::populate_buffer(std::vector<amrex::Real>& buf)
{
    BL_PROFILE("amr-wind::SamplingContainer::populate_buffer");

    amrex::Gpu::DeviceVector<amrex::Real> dbuf(buf.size(), 0.0_rt);
    auto* dbuf_ptr = dbuf.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (int fid = 0; fid < NumRuntimeRealComps(); ++fid) {
            const long offset = fid * num_sampling_particles();
            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto* pstruct = pti.GetArrayOfStructs()().data();
                auto* parr = pti.GetStructOfArrays().GetRealData(fid).data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) {
                    auto& pp = pstruct[ip];
                    const int pidx = pp.idata(IIx::uid);
                    const long ii = offset + pidx;
                    dbuf_ptr[ii] = parr[ip];
                });
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dbuf.begin(), dbuf.end(), buf.begin());
    amrex::ParallelDescriptor::ReduceRealSum(
        buf.data(), static_cast<int>(buf.size()),
        amrex::ParallelDescriptor::IOProcessorNumber());
}

} // namespace amr_wind::sampling
