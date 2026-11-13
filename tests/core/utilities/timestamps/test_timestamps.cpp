#include "core/utilities/timestamps.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

using namespace PHARE::core;

TEST(KahanTimeStamper, accumulatesConstantDt)
{
    KahanTimeStamper stamper{0.1};

    EXPECT_DOUBLE_EQ(stamper += 0.1, 0.1);
    EXPECT_DOUBLE_EQ(stamper += 0.1, 0.2);
    EXPECT_DOUBLE_EQ(stamper += 0.1, 0.3);
}

TEST(KahanTimeStamper, staysContinuousWhenDtChanges)
{
    KahanTimeStamper stamper{0.1};

    EXPECT_DOUBLE_EQ(stamper += 0.1, 0.1);
    EXPECT_DOUBLE_EQ(stamper += 0.1, 0.2);

    // dt grows: the accumulated time must keep advancing from where it left off, not reset
    EXPECT_DOUBLE_EQ(stamper += 0.2, 0.4);
    EXPECT_DOUBLE_EQ(stamper += 0.2, 0.6);

    // dt shrinks again
    EXPECT_DOUBLE_EQ(stamper += 0.05, 0.65);
    EXPECT_DOUBLE_EQ(stamper += 0.05, 0.7);
}

TEST(KahanTimeStamper, seedsFromANonZeroInitTime)
{
    // this is the restart case: TimeStamperFactory always seeds init_time at 0 and relies on the
    // first call resetting dt_, but the class itself supports a nonzero seed
    KahanTimeStamper stamper{0.1, /*init_time=*/5.0};

    EXPECT_DOUBLE_EQ(stamper += 0.1, 5.1);
    EXPECT_DOUBLE_EQ(stamper += 0.1, 5.2);
}

TEST(DtBetween, fromZero)
{
    EXPECT_EQ(0.0 + dt_between(0.0, 0.1), 0.1);
}

TEST(DtBetween, nearbyTimes)
{
    double const t = 3.0 * 0.1;
    EXPECT_EQ(t + dt_between(t, 0.4), 0.4);
}

TEST(DtBetween, farApartTimes)
{
    double const t = 1e-20;
    EXPECT_EQ(t + dt_between(t, 0.1), 0.1);
}

TEST(DtBetween, randomPairsLandExactlyWithinAFactorTwo)
{
    std::mt19937_64 gen{42};
    std::uniform_real_distribution<double> u{0.0, 1e3};
    for (int i = 0; i < 100000; ++i)
    {
        auto const target = u(gen);
        auto const t      = std::uniform_real_distribution<double>{target / 2, target}(gen);
        if (!(target > t))
            continue;
        EXPECT_EQ(t + dt_between(t, target), target);
    }
}

TEST(DtBetween, landsExactlyOnARoundedSum)
{
    std::mt19937_64 gen{7};
    std::uniform_real_distribution<double> u{0.0, 1e3};
    for (int i = 0; i < 100000; ++i)
    {
        auto const t      = u(gen) * 1e-3;
        auto const target = t + u(gen);
        if (!(target > t))
            continue;
        EXPECT_EQ(t + dt_between(t, target), target);
    }
}

TEST(DtBetween, unreachableEndOvershootsByOneUlp)
{
    double const t      = 19.30026167262983;
    double const target = 297.62996909823011;
    auto const landed   = t + dt_between(t, target);
    EXPECT_EQ(landed, std::nextafter(target, std::numeric_limits<double>::infinity()));
}

TEST(DtBetween, randomPairsNeverLandBelowEnd)
{
    std::mt19937_64 gen{42};
    std::uniform_real_distribution<double> u{0.0, 1e3};
    auto constexpr inf = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 100000; ++i)
    {
        auto const a = u(gen);
        auto const b = u(gen);
        if (a == b)
            continue;
        auto const t      = std::min(a, b);
        auto const target = std::max(a, b);
        auto const landed = t + dt_between(t, target);
        EXPECT_GE(landed, target);
        EXPECT_LE(landed, std::nextafter(target, inf));
    }
}

TEST(DtBetween, throwsWhenEndIsNotAfterStart)
{
    EXPECT_THROW(dt_between(2.0, 2.0), std::runtime_error);
}


TEST(ReachesTarget, acceptsTargetAndOneUlpAbove)
{
    auto constexpr inf = std::numeric_limits<double>::infinity();
    EXPECT_TRUE(reaches_target(0.4, 0.4));
    EXPECT_TRUE(reaches_target(std::nextafter(0.4, inf), 0.4));
}

TEST(ReachesTarget, rejectsBelowOrFurtherAbove)
{
    auto constexpr inf = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(reaches_target(std::nextafter(0.4, 0.0), 0.4));
    EXPECT_FALSE(reaches_target(std::nextafter(std::nextafter(0.4, inf), inf), 0.4));
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
