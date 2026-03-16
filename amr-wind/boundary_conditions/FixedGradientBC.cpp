#include "amr-wind/boundary_conditions/FixedGradientBC.H"

namespace amr_wind {
namespace {
AMREX_FORCE_INLINE
amrex::Box lower_boundary_faces(const amrex::Box& b, int dir)
{
    amrex::IntVect lo(b.smallEnd());
    amrex::IntVect hi(b.bigEnd());
    int sm = lo[dir];
    lo.setVal(dir, sm - 1);
    hi.setVal(dir, sm - 1);
    amrex::IndexType bxtype(b.ixType());
    bxtype.set(dir);
    return {lo, hi, bxtype};
}
} // namespace

FixedGradientBC::FixedGradientBC(Field& field, amrex::Orientation ori)
    : m_field(field), m_ori(ori)
{}

void FixedGradientBC::operator()(Field& field, const FieldState /*rho_state*/)
{
    const auto& repo = m_field.repo();
    const auto bcvals = field.bc_values_device();
    const int ncomp = field.num_comp();
    const int idx = static_cast<int>(m_ori);
    const int idim = m_ori.coordDir();
    const auto islow = m_ori.isLow();
    const auto ishigh = m_ori.isHigh();

    const int nlevels = field.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& domain = repo.mesh().Geom(lev).Domain();

        amrex::MFItInfo mfi_info{};
        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(field(lev), mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            const auto& bc_a = field(lev).array(mfi);

            if (islow && (bx.smallEnd(idim) == domain.smallEnd(idim))) {
                amrex::ParallelFor(
                    lower_boundary_faces(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int n = 0; n < ncomp; ++n) {
                            bc_a(i, j, k, n) = bcvals[idx][n];
                        }
                    });
            }

            if (ishigh && (bx.bigEnd(idim) == domain.bigEnd(idim))) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int n = 0; n < ncomp; ++n) {
                            bc_a(i, j, k, n) = bcvals[idx][n];
                        }
                    });
            }
        }
    }
}

} // namespace amr_wind
