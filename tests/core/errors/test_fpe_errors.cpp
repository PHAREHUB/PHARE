

#include "core/errors.hpp"

#include <cmath>

#include "gtest/gtest.h"

namespace PHARE
{

TEST(FloatingPointException, occurs)
{
    EXPECT_EQ(FPEWatcher::depth, 0);
    {
        PHARE_FPE_SCOPE;
        EXPECT_EQ(FPEWatcher::depth, 1);
        ASSERT_EXIT(
            {
                double volatile result = 0;
                auto const cause_fpe = [](double const a = 1, double const b = 0) { return a / b; };
                result               = cause_fpe();
                (void)result;
                // never should go here!
                exit(EXIT_SUCCESS);
            },
            testing::KilledBySignal(SIGFPE), "");
    }
    EXPECT_EQ(FPEWatcher::depth, 0);
}

TEST(FloatingPointException, doesnt_occur)
{
    EXPECT_EQ(FPEWatcher::depth, 0);
    {
        PHARE_FPE_SCOPE;
    }
    EXPECT_EQ(FPEWatcher::depth, 0);
    ASSERT_EXIT(
        {
            double volatile result = 0;
            auto const cause_fpe   = [](double const a = 1, double const b = 0) { return a / b; };
            result                 = cause_fpe();
            (void)result;
            exit(EXIT_SUCCESS);
        },
        testing::ExitedWithCode(EXIT_SUCCESS), "");
    EXPECT_EQ(FPEWatcher::depth, 0);
}


} // namespace PHARE

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
