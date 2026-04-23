#ifndef AMREXTEST_H
#define AMREXTEST_H

#include "gtest/gtest.h"
#include "AMReX_ParmParse.H"

namespace kynema_sgf_tests {

/** Base class for Kynema-SGF unit tests
 *
 *  Resets the ParmParse global dictionary for each test in this test fixture
 */
class AmrexTest : public ::testing::Test
{
public:
    AmrexTest() = default;

    ~AmrexTest() override = default;

    void SetUp() override
    {
        amrex::ParmParse pp("utest");
        pp.query("keep_parameters", m_keep_parameters);

        if (!m_keep_parameters) {
            amrex::ParmParse::Initialize(0, nullptr, nullptr);
        }
    }

    void TearDown() override
    {
        if (!m_keep_parameters) {
            amrex::ParmParse::Finalize();
        }
    }

protected:
    //! Flag indicating whether the ParmParse finalize method is called in
    //! teardown
    bool m_keep_parameters{false};
};

} // namespace kynema_sgf_tests

#endif /* AMREXTEST_H */
