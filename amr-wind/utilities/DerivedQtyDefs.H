#ifndef DERIVEDQTYDEFS_H
#define DERIVEDQTYDEFS_H

#include "amr-wind/utilities/DerivedQuantity.H"

namespace amr_wind::derived {
struct VorticityMag : public DerivedQty::Register<VorticityMag>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "mag_vorticity"; }

    VorticityMag(const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override { return identifier(); }

    [[nodiscard]] int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field& m_vel;
};

struct QCriterion : public DerivedQty::Register<QCriterion>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "q_criterion"; }

    QCriterion(const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override { return identifier(); }

    [[nodiscard]] int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field& m_vel;
};

struct QCriterionNondim : public DerivedQty::Register<QCriterionNondim>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "q_criterion_nondim"; }

    QCriterionNondim(
        const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override { return identifier(); }

    [[nodiscard]] int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field& m_vel;
};

struct StrainRateMag : public DerivedQty::Register<StrainRateMag>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "mag_strainrate"; }

    StrainRateMag(const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override { return identifier(); }

    [[nodiscard]] int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field& m_vel;
};

struct Gradient : public DerivedQty::Register<Gradient>
{
    static std::string identifier() { return "grad"; }

    Gradient(const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override
    {
        return "grad_" + m_phi->name();
    }

    [[nodiscard]] int num_comp() const override
    {
        return AMREX_SPACEDIM * m_phi->num_comp();
    }

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field* m_phi;
};

struct Divergence : public DerivedQty::Register<Divergence>
{
    static std::string identifier() { return "div"; }

    Divergence(const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override
    {
        return "div_" + m_phi->name();
    }

    [[nodiscard]] int num_comp() const override { return 1; }

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field* m_phi;
};

struct Laplacian : public DerivedQty::Register<Laplacian>
{
    static constexpr int ncomp = 1;
    static std::string identifier() { return "laplacian"; }

    Laplacian(const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override
    {
        return "laplacian_" + m_phi->name();
    }

    [[nodiscard]] int num_comp() const override { return ncomp; }

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field* m_phi;
};

struct FieldComponents : public DerivedQty::Register<FieldComponents>
{
    static std::string identifier() { return "components"; }

    FieldComponents(
        const FieldRepo& repo, const std::vector<std::string>& args);

    [[nodiscard]] std::string name() const override { return m_fld->name(); }

    [[nodiscard]] int num_comp() const override { return m_ncomp; }

    void var_names(amrex::Vector<std::string>& plt_var_names) override;

    void operator()(ScratchField& fld, int scomp = 0) const override;

private:
    const Field* m_fld;
    amrex::Vector<int> m_comp;
    int m_ncomp{0};
};

} // namespace amr_wind::derived

#endif /* DERIVEDQTYDEFS_H */
