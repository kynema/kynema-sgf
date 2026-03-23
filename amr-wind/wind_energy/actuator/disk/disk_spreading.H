#ifndef DISK_SPREADING_H
#define DISK_SPREADING_H

#include "AMReX_REAL.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/wind_energy/actuator/disk/UniformCt.H"
#include "amr-wind/core/FieldRepo.H"

using namespace amrex::literals;

namespace amr_wind::actuator::ops {

/**
 * \brief  A collection of spreading functions
 * This class allows for polymorphic spreading functions.
 * The concept is based on Impossibly Fast C++ Delegates proposed here:
 * https://www.codeproject.com/Articles/11015/The-Impossibly-Fast-C-Delegates
 *
 * @tparam T
 */
template <typename T>
class SpreadingFunction
{
public:
    using OwnerType = T;

    void operator()(
        const T& actObj,
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom)
    {
        (this->*m_function)(actObj, lev, mfi, geom);
    }

    SpreadingFunction(const SpreadingFunction&) = delete;
    void operator=(const SpreadingFunction&) = delete;

    void (SpreadingFunction::*m_function)(
        const T& actObj,
        const int,
        const amrex::MFIter&,
        const amrex::Geometry&);

    void uniform_gaussian_spreading(
        const T& actObj,
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom)
    {
        const auto& bx = mfi.tilebox();

        const auto bxa = amr_wind::utils::realbox_to_box(
            actObj.m_data.info().bound_box, geom);
        const auto& bxi = bx & bxa;
        if (bxi.isEmpty()) {
            return;
        }

        const auto& sarr = actObj.m_act_src(lev).array(mfi);
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();

        const auto& data = actObj.m_data.meta();

        const vs::Vector epsilon = vs::Vector::one() * data.epsilon;
        const vs::Vector m_normal(data.normal_vec);
        const auto* pos = actObj.m_pos.data();
        const auto* force = actObj.m_force.data();
        const int npts = data.num_force_pts;
        const int nForceTheta = data.num_force_theta_pts;
        const auto dTheta = ::amr_wind::utils::two_pi() / nForceTheta;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const vs::Vector cc{
                problo[0] + ((i + 0.5_rt) * dx[0]),
                problo[1] + ((j + 0.5_rt) * dx[1]),
                problo[2] + ((k + 0.5_rt) * dx[2]),
            };

            amrex::RealArray src_force = {0.0_rt};
            for (int ip = 0; ip < npts; ++ip) {
                const auto& pforce = force[ip] / nForceTheta;
                const auto pLoc = pos[ip];

                for (int it = 0; it < nForceTheta; ++it) {
                    const amrex::Real angle =
                        ::amr_wind::utils::degrees(it * dTheta);
                    const auto rotMatrix = vs::quaternion(m_normal, angle);
                    const auto diskPoint = pLoc & rotMatrix;
                    const auto distance = diskPoint - cc;
                    const auto projection_weight =
                        utils::gaussian3d(distance, epsilon);

                    src_force[0] += projection_weight * pforce.x();
                    src_force[1] += projection_weight * pforce.y();
                    src_force[2] += projection_weight * pforce.z();
                }
            }

            sarr(i, j, k, 0) += src_force[0];
            sarr(i, j, k, 1) += src_force[1];
            sarr(i, j, k, 2) += src_force[2];
        });
    }

    void linear_basis_spreading(
        const T& actObj,
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom)
    {
        const auto& bx = mfi.tilebox();

        const auto bxa = amr_wind::utils::realbox_to_box(
            actObj.m_data.info().bound_box, geom);
        const auto& bxi = bx & bxa;
        if (bxi.isEmpty()) {
            return;
        }

        const auto& sarr = actObj.m_act_src(lev).array(mfi);
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();

        const auto& data = actObj.m_data.meta();

        const amrex::Real dR = data.dr;
        const amrex::Real epsilon = data.epsilon;
        const vs::Vector m_normal(data.normal_vec);
        const vs::Vector m_origin(data.center);
        const auto* pos = actObj.m_pos.data();
        const auto* force = actObj.m_force.data();
        const int npts = data.num_force_pts;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const vs::Vector cc{
                problo[0] + ((i + 0.5_rt) * dx[0]),
                problo[1] + ((j + 0.5_rt) * dx[1]),
                problo[2] + ((k + 0.5_rt) * dx[2]),
            };

            amrex::RealArray src_force = {0.0_rt};
            for (int ip = 0; ip < npts; ++ip) {
                const auto R =
                    utils::delta_pnts_cyl(m_origin, m_normal, m_origin, pos[ip])
                        .x();
                const auto dist_on_disk =
                    utils::delta_pnts_cyl(m_origin, m_normal, cc, pos[ip]);
                const auto& pforce = force[ip];

                const amrex::Real weight_R =
                    utils::linear_basis_1d(dist_on_disk.x(), dR);
                const amrex::Real weight_T =
                    1.0_rt / (::amr_wind::utils::two_pi() * R);
                const amrex::Real weight_N =
                    utils::gaussian1d(dist_on_disk.z(), epsilon);
                const auto projection_weight = weight_R * weight_T * weight_N;

                src_force[0] += projection_weight * pforce.x();
                src_force[1] += projection_weight * pforce.y();
                src_force[2] += projection_weight * pforce.z();
            }

            sarr(i, j, k, 0) += src_force[0];
            sarr(i, j, k, 1) += src_force[1];
            sarr(i, j, k, 2) += src_force[2];
        });
    }
    void linear_basis_in_theta(
        const T& actObj,
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom)
    {
        const auto& bx = mfi.tilebox();

        const auto bxa = amr_wind::utils::realbox_to_box(
            actObj.m_data.info().bound_box, geom);
        const auto& bxi = bx & bxa;
        if (bxi.isEmpty()) {
            return;
        }

        const auto& sarr = actObj.m_act_src(lev).array(mfi);
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();

        const auto& data = actObj.m_data.meta();

        const amrex::Real dR = data.dr;
        const amrex::Real dTheta =
            ::amr_wind::utils::two_pi() / data.num_vel_pts_t;
        const amrex::Real epsilon = data.epsilon;
        const vs::Vector m_normal(data.normal_vec);
        const vs::Vector m_origin(data.center);
        const auto* pos = actObj.m_pos.data();
        const auto* force = actObj.m_force.data();
        const int npts = data.num_force_pts;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const vs::Vector cc{
                problo[0] + ((i + 0.5_rt) * dx[0]),
                problo[1] + ((j + 0.5_rt) * dx[1]),
                problo[2] + ((k + 0.5_rt) * dx[2]),
            };

            amrex::RealArray src_force = {0.0_rt};
            for (int ip = 0; ip < npts; ++ip) {
                const auto radius =
                    utils::delta_pnts_cyl(m_origin, m_normal, m_origin, pos[ip])
                        .x();
                const auto dArc = radius * dTheta;
                const auto dist_on_disk =
                    utils::delta_pnts_cyl(m_origin, m_normal, cc, pos[ip]);
                const amrex::Real arclength = dist_on_disk.y() * radius;
                const auto& pforce = force[ip];

                const amrex::Real weight_R =
                    utils::linear_basis_1d(dist_on_disk.x(), dR);
                const amrex::Real weight_T =
                    utils::linear_basis_1d(arclength, dArc);
                const amrex::Real weight_N =
                    utils::gaussian1d(dist_on_disk.z(), epsilon);
                const auto projection_weight = weight_R * weight_T * weight_N;

                src_force[0] += projection_weight * pforce.x();
                src_force[1] += projection_weight * pforce.y();
                src_force[2] += projection_weight * pforce.z();
            }

            sarr(i, j, k, 0) += src_force[0];
            sarr(i, j, k, 1) += src_force[1];
            sarr(i, j, k, 2) += src_force[2];
        });
    }

    SpreadingFunction() : m_function(&SpreadingFunction::linear_basis_spreading)
    {}
    void initialize(const std::string& key)
    {
        if (std::is_same_v<UniformCt, typename OwnerType::TraitType>) {
            if (key == "UniformGaussian") {
                m_function = &SpreadingFunction::uniform_gaussian_spreading;
            } else if (key == "LinearBasis") {
                m_function = &SpreadingFunction::linear_basis_spreading;
            } else {
                amrex::Abort("Invalid spreading type");
            }
        } else {
            m_function = &SpreadingFunction::linear_basis_in_theta;
        }
    }
};
} // namespace amr_wind::actuator::ops
#endif /* DISK_SPREADING_H */
