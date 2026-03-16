#include "mms_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "amr-wind/physics/mms/MMS.H"
#include "amr-wind/physics/mms/MMSForcing.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

TEST_F(MMSMeshTest, mms_forcing)
{
#if defined(AMREX_USE_HIP)
    GTEST_SKIP();
#else
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Initialize parameters
    utils::populate_mms_params();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    for (auto& pp : sim().physics()) {
        pp->post_init_actions();
    }

    auto fields = ICNSFields(sim())(sim().time());
    auto& src_term = fields.src_term;
    amr_wind::pde::icns::mms::MMSForcing mmsforcing(sim());

    const amrex::Array<amrex::Real, AMREX_SPACEDIM> min_golds = {
        -2.1397143441391857_rt, -2.5061563892200622_rt, -2.6756003260809429_rt};
    const amrex::Array<amrex::Real, AMREX_SPACEDIM> max_golds = {
        2.0381534755116628_rt, 2.2014865191023762_rt, 2.4125363807493985_rt};
    src_term.setVal(0.0_rt);

    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        mmsforcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, min_golds[i], tol);
        EXPECT_NEAR(max_val, max_golds[i], tol);
    }
#endif
}
} // namespace amr_wind_tests
