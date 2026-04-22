#include "amr-wind/physics/udfs/UDF.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/physics/udfs/LinearProfile.H"
#include "amr-wind/physics/udfs/PowerLawProfile.H"
#include "amr-wind/physics/udfs/BurggrafLid.H"
#include "amr-wind/physics/udfs/Rankine.H"
#include "amr-wind/physics/udfs/CustomVelocity.H"
#include "amr-wind/physics/udfs/CustomScalar.H"
#include "amr-wind/physics/udfs/TwoLayer.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::udf {

ConstValue::ConstValue(Field& fld) : m_field(fld)
{
    amrex::ParmParse pp_old("incflo");
    amrex::ParmParse pp(identifier() + "." + m_field.name());
    pp_old.queryarr(m_field.name(), m_value);
    pp.queryarr("value", m_value);

    if (fld.num_comp() != m_value.size()) {
        amrex::Abort("UDF: Invalid value for field: " + m_field.name());
    }
}

void ConstValue::operator()(int level, const amrex::Geometry& /*geom*/)
{
    auto& mfab = m_field(level);
    for (int i = 0; i < m_field.num_comp(); ++i) {
        mfab.setVal(m_value[i], i, 1);
    }
}

template <typename T>
UDFImpl<T>::UDFImpl(Field& fld) : m_field(fld), m_op(fld)
{}

template <typename T>
void UDFImpl<T>::operator()(int level, const amrex::Geometry& geom)
{
    auto& mfab = m_field(level);
    const auto& geomData = geom.data();

    const amrex::Real time = 0.0_rt;
    const auto ncomp = m_field.num_comp();
    const auto& dop = m_op.device_instance();
    const auto& marrs = mfab.arrays();

    amrex::ParallelFor(
        mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            for (int n = 0; n < ncomp; ++n) {
                dop({i, j, k}, marrs[nbx], geomData, time, {}, n, 0, 0);
            }
        });
    amrex::Gpu::streamSynchronize();
}

template class UDFImpl<LinearProfile>;
template class UDFImpl<PowerLawProfile>;
template class UDFImpl<BurggrafLid>;
template class UDFImpl<Rankine>;
template class UDFImpl<CustomVelocity>;
template class UDFImpl<CustomScalar>;
template class UDFImpl<TwoLayer>;

} // namespace amr_wind::udf
