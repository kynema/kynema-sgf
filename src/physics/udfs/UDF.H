#ifndef UDF_H
#define UDF_H

#include "src/core/Factory.H"
#include "AMReX_Gpu.H"

namespace kynema_sgf {
class Field;

namespace udf {

class UDF : public Factory<UDF, Field&>
{
public:
    static std::string base_identifier() { return "UDF"; };

    virtual void operator()(int level, const amrex::Geometry& geom) = 0;
};

class ConstValue : public UDF::Register<ConstValue>
{
public:
    static std::string identifier() { return "ConstValue"; }

    explicit ConstValue(Field& /*fld*/);

    void operator()(int level, const amrex::Geometry& geom) override;

private:
    Field& m_field;
    amrex::Vector<amrex::Real> m_value;
};

template <typename T>
class UDFImpl : public UDF::Register<UDFImpl<T>>
{
public:
    static std::string identifier() { return T::identifier(); }

    explicit UDFImpl(Field& /*fld*/);

    void operator()(int level, const amrex::Geometry& geom) override;

private:
    Field& m_field;
    T m_op;
};

} // namespace udf
} // namespace kynema_sgf

#endif /* UDF_H */
