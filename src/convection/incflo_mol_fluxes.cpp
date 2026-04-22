#include "amr-wind/convection/MOL.H"
#include "amr-wind/utilities/bc_ops.H"
#include "AMReX_Geometry.H"
#include "AMReX_Slopes_K.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

void mol::compute_convective_rate(
    amrex::Box const& bx,
    int ncomp,
    amrex::Array4<amrex::Real> const& dUdt,
    amrex::Array4<amrex::Real const> const& fx,
    amrex::Array4<amrex::Real const> const& fy,
    amrex::Array4<amrex::Real const> const& fz,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxi)
{
    BL_PROFILE("amr-wind::mol::compute_convective_rate");
    const auto dxinv = dxi;
    amrex::ParallelFor(
        bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
            dUdt(i, j, k, n) =
                (dxinv[0] * (fx(i, j, k, n) - fx(i + 1, j, k, n))) +
                (dxinv[1] * (fy(i, j, k, n) - fy(i, j + 1, k, n))) +
                (dxinv[2] * (fz(i, j, k, n) - fz(i, j, k + 1, n)));
        });
}

void mol::compute_convective_fluxes(
    int lev,
    amrex::Box const& bx,
    int ncomp,
    amrex::Array4<amrex::Real> const& fx,
    amrex::Array4<amrex::Real> const& fy,
    amrex::Array4<amrex::Real> const& fz,
    amrex::Array4<amrex::Real const> const& q,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::BCRec const* h_bcrec,
    amrex::BCRec const* d_bcrec,
    amrex::Vector<amrex::Geometry> geom)
{
    BL_PROFILE("amr-wind::mol::compute_convective_fluxes");
    constexpr amrex::Real small_vel = std::numeric_limits<float>::epsilon();

    const amrex::Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    amrex::Box const& xbx = amrex::surroundingNodes(bx, 0);
    amrex::Box const& ybx = amrex::surroundingNodes(bx, 1);
    amrex::Box const& zbx = amrex::surroundingNodes(bx, 2);

    // At an ext_dir boundary, the boundary value is on the face, not cell
    // center.
    auto extdir_lohi = amr_wind::utils::has_extdir_or_ho(
        h_bcrec, ncomp, static_cast<int>(amrex::Direction::x));
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_ilo >= xbx.smallEnd(0) - 1) or
        (has_extdir_hi and domain_ihi <= xbx.bigEnd(0))) {
        amrex::ParallelFor(
            xbx, ncomp,
            [d_bcrec, q, domain_ilo, domain_ihi, umac,
             fx] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                bool extdir_or_ho_ilo =
                    (d_bcrec[n].lo(0) == amrex::BCType::ext_dir) ||
                    (d_bcrec[n].lo(0) == amrex::BCType::hoextrap);
                bool extdir_or_ho_ihi =
                    (d_bcrec[n].hi(0) == amrex::BCType::ext_dir) ||
                    (d_bcrec[n].hi(0) == amrex::BCType::hoextrap);
                amrex::Real qs;
                if (i <= domain_ilo &&
                    (d_bcrec[n].lo(0) == amrex::BCType::ext_dir)) {
                    qs = q(domain_ilo - 1, j, k, n);
                } else if (
                    i >= domain_ihi + 1 &&
                    (d_bcrec[n].hi(0) == amrex::BCType::ext_dir)) {
                    qs = q(domain_ihi + 1, j, k, n);
                } else {
                    amrex::Real qpls =
                        q(i, j, k, n) -
                        (0.5_rt * amrex_calc_xslope_extdir(
                                      i, j, k, n, 2, q, extdir_or_ho_ilo,
                                      extdir_or_ho_ihi, domain_ilo,
                                      domain_ihi));
                    amrex::Real qmns =
                        q(i - 1, j, k, n) +
                        (0.5_rt * amrex_calc_xslope_extdir(
                                      i - 1, j, k, n, 2, q, extdir_or_ho_ilo,
                                      extdir_or_ho_ihi, domain_ilo,
                                      domain_ihi));
                    if (umac(i, j, k) > small_vel) {
                        qs = qmns;
                    } else if (umac(i, j, k) < -small_vel) {
                        qs = qpls;
                    } else {
                        qs = 0.5_rt * (qmns + qpls);
                    }
                }
                fx(i, j, k, n) = qs * umac(i, j, k);
            });
    } else {
        amrex::ParallelFor(
            xbx, ncomp,
            [q, umac, fx] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                amrex::Real qpls =
                    q(i, j, k, n) -
                    (0.5_rt * amrex_calc_xslope(i, j, k, n, 2, q));
                amrex::Real qmns =
                    q(i - 1, j, k, n) +
                    (0.5_rt * amrex_calc_xslope(i - 1, j, k, n, 2, q));
                amrex::Real qs;
                if (umac(i, j, k) > small_vel) {
                    qs = qmns;
                } else if (umac(i, j, k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5_rt * (qmns + qpls);
                }
                fx(i, j, k, n) = qs * umac(i, j, k);
            });
    }

    extdir_lohi = amr_wind::utils::has_extdir_or_ho(
        h_bcrec, ncomp, static_cast<int>(amrex::Direction::y));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_jlo >= ybx.smallEnd(1) - 1) or
        (has_extdir_hi and domain_jhi <= ybx.bigEnd(1))) {
        amrex::ParallelFor(
            ybx, ncomp,
            [d_bcrec, q, domain_jlo, domain_jhi, vmac,
             fy] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                bool extdir_or_ho_jlo =
                    (d_bcrec[n].lo(1) == amrex::BCType::ext_dir) ||
                    (d_bcrec[n].lo(1) == amrex::BCType::hoextrap);
                bool extdir_or_ho_jhi =
                    (d_bcrec[n].hi(1) == amrex::BCType::ext_dir) ||
                    (d_bcrec[n].hi(1) == amrex::BCType::hoextrap);
                amrex::Real qs;
                if (j <= domain_jlo &&
                    (d_bcrec[n].lo(1) == amrex::BCType::ext_dir)) {
                    qs = q(i, domain_jlo - 1, k, n);
                } else if (
                    j >= domain_jhi + 1 &&
                    (d_bcrec[n].hi(1) == amrex::BCType::ext_dir)) {
                    qs = q(i, domain_jhi + 1, k, n);
                } else {
                    amrex::Real qpls =
                        q(i, j, k, n) -
                        (0.5_rt * amrex_calc_yslope_extdir(
                                      i, j, k, n, 2, q, extdir_or_ho_jlo,
                                      extdir_or_ho_jhi, domain_jlo,
                                      domain_jhi));
                    amrex::Real qmns =
                        q(i, j - 1, k, n) +
                        (0.5_rt * amrex_calc_yslope_extdir(
                                      i, j - 1, k, n, 2, q, extdir_or_ho_jlo,
                                      extdir_or_ho_jhi, domain_jlo,
                                      domain_jhi));
                    if (vmac(i, j, k) > small_vel) {
                        qs = qmns;
                    } else if (vmac(i, j, k) < -small_vel) {
                        qs = qpls;
                    } else {
                        qs = 0.5_rt * (qmns + qpls);
                    }
                }
                fy(i, j, k, n) = qs * vmac(i, j, k);
            });
    } else {
        amrex::ParallelFor(
            ybx, ncomp,
            [q, vmac, fy] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                amrex::Real qpls =
                    q(i, j, k, n) -
                    (0.5_rt * amrex_calc_yslope(i, j, k, n, 2, q));
                amrex::Real qmns =
                    q(i, j - 1, k, n) +
                    (0.5_rt * amrex_calc_yslope(i, j - 1, k, n, 2, q));
                amrex::Real qs;
                if (vmac(i, j, k) > small_vel) {
                    qs = qmns;
                } else if (vmac(i, j, k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5_rt * (qmns + qpls);
                }
                fy(i, j, k, n) = qs * vmac(i, j, k);
            });
    }

    extdir_lohi = amr_wind::utils::has_extdir_or_ho(
        h_bcrec, ncomp, static_cast<int>(amrex::Direction::z));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_klo >= zbx.smallEnd(2) - 1) or
        (has_extdir_hi and domain_khi <= zbx.bigEnd(2))) {
        amrex::ParallelFor(
            zbx, ncomp,
            [d_bcrec, q, domain_klo, domain_khi, wmac,
             fz] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                bool extdir_or_ho_klo =
                    (d_bcrec[n].lo(2) == amrex::BCType::ext_dir) ||
                    (d_bcrec[n].lo(2) == amrex::BCType::hoextrap);
                bool extdir_or_ho_khi =
                    (d_bcrec[n].hi(2) == amrex::BCType::ext_dir) ||
                    (d_bcrec[n].hi(2) == amrex::BCType::hoextrap);
                amrex::Real qs;
                if (k <= domain_klo &&
                    (d_bcrec[n].lo(2) == amrex::BCType::ext_dir)) {
                    qs = q(i, j, domain_klo - 1, n);
                } else if (
                    k >= domain_khi + 1 &&
                    (d_bcrec[n].hi(2) == amrex::BCType::ext_dir)) {
                    qs = q(i, j, domain_khi + 1, n);
                } else {
                    amrex::Real qpls =
                        q(i, j, k, n) -
                        (0.5_rt * amrex_calc_zslope_extdir(
                                      i, j, k, n, 2, q, extdir_or_ho_klo,
                                      extdir_or_ho_khi, domain_klo,
                                      domain_khi));
                    amrex::Real qmns =
                        q(i, j, k - 1, n) +
                        (0.5_rt * amrex_calc_zslope_extdir(
                                      i, j, k - 1, n, 2, q, extdir_or_ho_klo,
                                      extdir_or_ho_khi, domain_klo,
                                      domain_khi));
                    if (wmac(i, j, k) > small_vel) {
                        qs = qmns;
                    } else if (wmac(i, j, k) < -small_vel) {
                        qs = qpls;
                    } else {
                        qs = 0.5_rt * (qmns + qpls);
                    }
                }
                fz(i, j, k, n) = qs * wmac(i, j, k);
            });
    } else {
        amrex::ParallelFor(
            zbx, ncomp,
            [q, wmac, fz] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                amrex::Real qpls =
                    q(i, j, k, n) -
                    (0.5_rt * amrex_calc_zslope(i, j, k, n, 2, q));
                amrex::Real qmns =
                    q(i, j, k - 1, n) +
                    (0.5_rt * amrex_calc_zslope(i, j, k - 1, n, 2, q));
                amrex::Real qs;
                if (wmac(i, j, k) > small_vel) {
                    qs = qmns;
                } else if (wmac(i, j, k) < -small_vel) {
                    qs = qpls;
                } else {
                    qs = 0.5_rt * (qmns + qpls);
                }
                fz(i, j, k, n) = qs * wmac(i, j, k);
            });
    }
}
