#ifndef BCINTERFACE_H
#define BCINTERFACE_H

#include "amr-wind/core/Field.H"
#include "amr-wind/incflo_enums.H"

#include "AMReX_BC_TYPES.H"

#include <utility>
#include "AMReX_REAL.H"

using namespace amrex::literals;

/**
 *  \defgroup field_bc Boundary conditions management
 *  Boundary conditions management
 *
 *  This group documents the boundary condition management interface within
 *  AMR-Wind. Boundary conditions are registered on a Field instance and are
 *  used by the Field::fillpatch and Field::fillphysbc calls to perform the
 *  actual work. The classes documented in this group provide the interface to
 *  read user inputs and initialize the BC data structures in FieldInfo, moves
 *  them to the device so that they are ready for use during the solution
 *  process, and fill the ghost cells with appropriate values that require
 *  custom computations and not simple extrapolations.
 *
 *  \ingroup fields
 */

namespace amr_wind {
namespace {
amrex::Vector<std::string> bcnames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};
}

/** Abstract class that defines the interface for processing BCs
 *  \ingroup field_bc
 *
 *  This class defines the infrastructure to read and initialize the boundary
 *  condition information from the input file and populate the necessary data in
 *  FieldInfo instance for a given Field.
 */
class BCIface
{
public:
    using IncfloBC = amrex::GpuArray<BC, AMREX_SPACEDIM * 2>;

    BCIface(Field& field);

    virtual ~BCIface() = default;

    //! Operator that performs init actions and syncs the BC data to device
    virtual void operator()(amrex::Real value = 0.0_rt);

    //! User-defined functions for Dirichlet-type boundaries
    amrex::Array<const std::string, 3> get_dirichlet_udfs();

protected:
    //! Setup AMReX mathematical BC types
    virtual void set_bcrec() = 0;

    //! Read input data for specific fields based on BC types
    virtual void read_values() = 0;

    //! Parse custom BCs specified for this field
    virtual void read_bctype();

    //! BC functions for Neumann type boundaries
    virtual void set_bcfuncs();

    //! Set default BC values for the field
    void set_default_value(amrex::Real value);

    //! Set AMReX mathematical boundary types for the lower boundaries
    void set_bcrec_lo(int dir, amrex::BCType::mathematicalBndryTypes bcrec)
    {
        auto& fbcrec = m_field.bcrec();
        for (int i = 0; i < m_field.num_comp(); ++i) {
            fbcrec[i].setLo(dir, bcrec);
        }
    }

    //! Set AMReX mathematical boundary types for the upper boundaries
    void set_bcrec_hi(int dir, amrex::BCType::mathematicalBndryTypes bcrec)
    {
        auto& fbcrec = m_field.bcrec();
        for (int i = 0; i < m_field.num_comp(); ++i) {
            fbcrec[i].setHi(dir, bcrec);
        }
    }

    //! Field instance where BC is being set
    Field& m_field;
};

/** Boundary condition initialization for velocity (ICNS momentum) field
 *  \ingroup field_bc
 */
class BCVelocity : public BCIface
{
public:
    explicit BCVelocity(Field& field) : BCIface(field) {}

protected:
    void set_bcrec() override;
    void read_values() override;
};

/** Boundary condition initialization for scalar transport variables
 *  \ingroup field_bc
 */
class BCScalar : public BCIface
{
public:
    explicit BCScalar(Field& field) : BCIface(field) {}

protected:
    void set_bcrec() override;
    void read_values() override;
};

/** Specialization of BC initialization for pressure
 *
 *  For pressure there is no user value read from file. The values are assumed
 *  to be zero.
 */
class BCPressure : public BCScalar
{
public:
    explicit BCPressure(Field& field) : BCScalar(field) {}

protected:
    void read_values() override;
};

/** Specialization for source terms
 *  \ingroup field_bc
 *
 *  Sets a default first-order extrapolation for non-periodic BCs and a no-op
 *  for fillpatch operations.
 */
class BCSrcTerm : public BCIface
{
public:
    explicit BCSrcTerm(Field& field) : BCIface(field) {}

protected:
    void set_bcrec() override;
    void read_values() override
    {
        // No values specified for source terms. empty method to satisfy
        // inheritance
    }
};

/** Default BC that can be applied on any field
 *  \ingroup field_bc
 *
 *  This operator is used to initialize defaults for any field that might
 *  require valid BC data. It sets higher-order extrapolation for mathematical
 *  BC.
 */
class BCFillPatchExtrap : public BCIface
{
public:
    explicit BCFillPatchExtrap(
        Field& field,
        amrex::BCType::mathematicalBndryTypes bctype = amrex::BCType::hoextrap)
        : BCIface(field), m_extrap_type(bctype)
    {}

protected:
    void set_bcrec() override;
    void read_values() override
    {
        // No values specified for source terms. empty method to satisfy
        // inheritance
    }

    amrex::BCType::mathematicalBndryTypes m_extrap_type;
};

} // namespace amr_wind

#endif /* BCINTERFACE_H */
