/** \file utest_main.cpp
 *  Entry point for unit tests
 */

#include "gtest/gtest.h"
#include "ks_test_utils/AmrexTestEnv.H"

//! Global instance of the environment (for access in tests)
amr_wind_tests::AmrexTestEnv* utest_env = nullptr;

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    utest_env = new amr_wind_tests::AmrexTestEnv(argc, argv);
    ::testing::AddGlobalTestEnvironment(utest_env);

    return RUN_ALL_TESTS();
}
