#ifndef SAMPLINGCONTAINER_H
#define SAMPLINGCONTAINER_H

#include <memory>
#include <cstdint>

#include "AMReX_AmrParticles.H"
#include "amr-wind/utilities/DerivedQuantity.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class Field;

namespace sampling {

class SamplerBase;

static constexpr int SNStructReal = 0;
static constexpr int SNStructInt = 3;
static constexpr int SNArrayReal = 0;
static constexpr int SNArrayInt = 0;

/** Indexing for the integer variables defined on the AMReX particle
 *  \ingroup sampling
 */
struct IIx
{
    enum Indices : std::uint8_t {
        uid = 0, ///< Unique identifier for this particle
        sid,     ///< Identifier of the set this particle belongs to
        nid      ///< Index within the set for this particle
    };
};

/** Specialization of the AMReX ParticleContainer object for sampling data
 *  \ingroup sampling
 *
 *  Represents the data probes as a collection of particles that are uniquely
 *  identified by three integer identifiers:
 *
 *   - a universal identifier (UID) that is unique to each particle,
 *   - a set identifier that represents a group (e.g., line, plane) of particles
 *   - a node identifier that represents indexing of this particle within the
 * group
 *
 *  The quantities of interest are themselves stored as runtime components based
 *  on user-defined inputs.
 *
 *
 *  Notes:
 *
 *   - The implementation uses linear interpolation in three directions to
 *     determine the data at a given probe location
 *
 *   - For non-nodal fields, the current implementation requires at-least one
 *     ghost cell to allow linear interpolation.
 *
 *   - Interpolation near domain boundaries does not currently handle `hoextrap`
 */
class SamplingContainer
    : public amrex::AmrParticleContainer<
          SNStructReal,
          SNStructInt,
          SNArrayReal,
          SNArrayInt>
{
public:
    explicit SamplingContainer(amrex::AmrCore& mesh)
        : amrex::AmrParticleContainer<
              SNStructReal,
              SNStructInt,
              SNArrayReal,
              SNArrayInt>(&mesh)
        , m_mesh(mesh)
    {}

    /** Define the container and allocate memory for runtime components
     */
    void setup_container(int num_real_components, int num_int_components = 0);

    /** Create particle information for all the sampling locations
     */
    void initialize_particles(
        const amrex::Vector<std::unique_ptr<SamplerBase>>& /*samplers*/);

    //! Perform field interpolation to sampling locations
    template <typename FType>
    void interpolate_fields(const amrex::Vector<FType>& fields, const int scomp)
    {
        BL_PROFILE("amr-wind::SamplingContainer::interpolate_fields");

        const int nlevels = m_mesh.finestLevel() + 1;

        for (int lev = 0; lev < nlevels; ++lev) {
            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                int scomp_curr = scomp;
                for (const auto* fld : fields) {
                    AMREX_ALWAYS_ASSERT(fld->num_grow() > amrex::IntVect{0});
                    const auto farr = (*fld)(lev).const_array(pti);
                    interpolate(
                        pti, farr, lev, fld->field_location(), fld->num_comp(),
                        scomp_curr);

                    scomp_curr += fld->num_comp();
                }
            }
        }
    }

    //! Perform derived field interpolation to sampling locations
    void interpolate_derived_fields(
        const DerivedQtyMgr& derived_mgr, const FieldRepo& repo, int scomp);

    //! Populate the buffer with data for all the particles
    void populate_buffer(std::vector<amrex::Real>& buf);

    long num_sampling_particles() const { return m_total_particles; }

    long& num_sampling_particles() { return m_total_particles; }

    /** Interpolate a field to the sampling locations
     *
     *  \param np Number of particles in the container
     *  \param ic Component of the field to be interpolated
     *  \param pvec Vector containing particle info
     *  \param pavec Array information for the real component data
     *  \param farr Array of field data for this multifab
     *  \param problo Problem low values
     *  \param dxi Inverse cell size array
     *  \param dx Cell size array
     *  \param offset Offsets for cell/node/face fields
     */
    template <typename FType>
    void sample_field(
        const int np,
        const int ic,
        SamplingContainer::ParticleVector& pvec,
        SamplingContainer::RealVector& pavec,
        const FType& farr,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& offset)
    {
        BL_PROFILE("amr-wind::SamplingContainer::sample_impl");

        auto* pstruct = pvec.data();
        auto* parr = pavec.data();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
            auto& p = pstruct[ip];
            // Determine offsets within the containing cell
            const amrex::Real x =
                (p.pos(0) - problo[0] - offset[0] * dx[0]) * dxi[0];
            const amrex::Real y =
                (p.pos(1) - problo[1] - offset[1] * dx[1]) * dxi[1];
            const amrex::Real z =
                (p.pos(2) - problo[2] - offset[2] * dx[2]) * dxi[2];

            // Index of the low corner
            const int i = static_cast<int>(std::floor(x));
            const int j = static_cast<int>(std::floor(y));
            const int k = static_cast<int>(std::floor(z));

            // Interpolation weights in each direction (linear basis)
            const amrex::Real wx_hi = (x - i);
            const amrex::Real wy_hi = (y - j);
            const amrex::Real wz_hi = (z - k);

            const amrex::Real wx_lo = 1.0_rt - wx_hi;
            const amrex::Real wy_lo = 1.0_rt - wy_hi;
            const amrex::Real wz_lo = 1.0_rt - wz_hi;

            parr[ip] = (wx_lo * wy_lo * wz_lo * farr(i, j, k, ic)) +
                       (wx_lo * wy_lo * wz_hi * farr(i, j, k + 1, ic)) +
                       (wx_lo * wy_hi * wz_lo * farr(i, j + 1, k, ic)) +
                       (wx_lo * wy_hi * wz_hi * farr(i, j + 1, k + 1, ic)) +
                       (wx_hi * wy_lo * wz_lo * farr(i + 1, j, k, ic)) +
                       (wx_hi * wy_lo * wz_hi * farr(i + 1, j, k + 1, ic)) +
                       (wx_hi * wy_hi * wz_lo * farr(i + 1, j + 1, k, ic)) +
                       (wx_hi * wy_hi * wz_hi * farr(i + 1, j + 1, k + 1, ic));
        });
    }

private:
    //! Interpolate from an array4 onto particles
    template <typename FType>
    void interpolate(
        const ParIterType& pti,
        const FType& farr,
        const int lev,
        const FieldLoc floc,
        const int ncomp,
        const int scomp)
    {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();
        const int np = pti.numParticles();
        auto& pvec = pti.GetArrayOfStructs()();
        int fidx = scomp;
        for (int ic = 0; ic < ncomp; ++ic) {
            auto& parr = pti.GetStructOfArrays().GetRealData(fidx++);

            switch (floc) {
            case FieldLoc::NODE: {
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                    0.0_rt, 0.0_rt, 0.0_rt};
                sample_field(np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                break;
            }

            case FieldLoc::CELL: {
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                    0.5_rt, 0.5_rt, 0.5_rt};
                sample_field(np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                break;
            }

            case FieldLoc::XFACE: {
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                    0.0_rt, 0.5_rt, 0.5_rt};
                sample_field(np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                break;
            }

            case FieldLoc::YFACE: {
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                    0.5_rt, 0.0_rt, 0.5_rt};
                sample_field(np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                break;
            }

            case FieldLoc::ZFACE: {
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                    0.5_rt, 0.5_rt, 0.0_rt};
                sample_field(np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                break;
            }
            }
        }
    }

    const amrex::AmrCore& m_mesh;

    long m_total_particles{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* SAMPLINGCONTAINER_H */
