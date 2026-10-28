#include "amr/load_balancing/coarsest_rebalance_check.hpp"

#include "gtest/gtest.h"

#include <vector>
#include <cstddef>

using namespace PHARE::amr;

namespace
{
LoadBalancerDetails details(bool const on_init, std::size_t const backoff_multiplier = 2)
{
    return {.active                            = true,
            .automatic                         = true,
            .on_init                           = on_init,
            .mode                              = "nppc",
            .tolerance                         = 0.05,
            .next_rebalance_backoff_multiplier = backoff_multiplier,
            .next_rebalance                    = 10,
            .max_next_rebalance                = 50};
}

struct Steps
{
    std::vector<std::size_t> checked;
    std::vector<std::size_t> rebalanced;
};

// advances nSteps time steps, every check sees the workloads given for that check, the last
// ones are repeated once exhausted
Steps advance(CoarsestRebalanceCheck& check, std::size_t const nSteps,
              std::vector<std::vector<double>> const& workloads_per_check)
{
    Steps steps;
    for (std::size_t step = 0; step < nSteps; ++step)
    {
        bool const rebalanced = check([&] {
            auto const i = std::min(steps.checked.size(), workloads_per_check.size() - 1);
            steps.checked.push_back(step);
            return workloads_per_check[i];
        });
        if (rebalanced)
            steps.rebalanced.push_back(step);
    }
    return steps;
}

std::vector<double> const balanced{1, 1, 1};
std::vector<double> const imbalanced{2, 1, 1};
} // namespace


TEST(CoarsestRebalanceCheck, balancedChecksBackOffFromNextRebalance)
{
    CoarsestRebalanceCheck check{details(/*on_init=*/true)};
    auto const steps = advance(check, 200, {balanced});

    // intervals 10, 20, 40, then capped at max_next_rebalance=50
    EXPECT_EQ(steps.checked, (std::vector<std::size_t>{0, 10, 30, 70, 120, 170}));
    EXPECT_TRUE(steps.rebalanced.empty());
}

TEST(CoarsestRebalanceCheck, firstCheckIsAtNextRebalanceWithoutOnInit)
{
    CoarsestRebalanceCheck check{details(/*on_init=*/false)};
    auto const steps = advance(check, 40, {balanced});

    // then the same intervals as with on_init: 10, 20, ...
    EXPECT_EQ(steps.checked, (std::vector<std::size_t>{10, 20}));
}

TEST(CoarsestRebalanceCheck, imbalancedRebalancesEveryNextRebalanceSteps)
{
    CoarsestRebalanceCheck check{details(/*on_init=*/true)};
    auto const steps = advance(check, 31, {imbalanced});

    EXPECT_EQ(steps.checked, (std::vector<std::size_t>{0, 10, 20, 30}));
    EXPECT_EQ(steps.rebalanced, steps.checked);
}

TEST(CoarsestRebalanceCheck, backOffRestartsFromNextRebalanceAfterARebalance)
{
    CoarsestRebalanceCheck check{details(/*on_init=*/true)};
    auto const steps = advance(check, 100, {balanced, balanced, imbalanced, balanced});

    // 0 and 10 balanced so the intervals grew 10, 20, 30 rebalances so they restart at 10
    EXPECT_EQ(steps.checked, (std::vector<std::size_t>{0, 10, 30, 40, 50, 70}));
    EXPECT_EQ(steps.rebalanced, (std::vector<std::size_t>{30}));
}

TEST(CoarsestRebalanceCheck, zeroBackOffMultiplierKeepsTheInterval)
{
    CoarsestRebalanceCheck check{details(/*on_init=*/true, /*backoff_multiplier=*/0)};
    auto const steps = advance(check, 31, {balanced});

    EXPECT_EQ(steps.checked, (std::vector<std::size_t>{0, 10, 20, 30}));
}

TEST(CoarsestRebalanceCheck, toleranceIsRelativeToTheAverageLoad)
{
    // max is 4% above the average: within tolerance, although the least loaded rank is ~8%
    // below the most loaded one
    EXPECT_FALSE(CoarsestRebalanceCheck::is_imbalanced({1.04, 1.0, 0.96}, 0.05));
    EXPECT_TRUE(CoarsestRebalanceCheck::is_imbalanced({1.06, 1.0, 0.94}, 0.05));

    CoarsestRebalanceCheck check{details(/*on_init=*/true)};
    EXPECT_TRUE(advance(check, 200, {{1.04, 1.0, 0.96}}).rebalanced.empty());
}

TEST(CoarsestRebalanceCheck, noLoadIsNotImbalanced)
{
    EXPECT_FALSE(CoarsestRebalanceCheck::is_imbalanced({}, 0.05));
    EXPECT_FALSE(CoarsestRebalanceCheck::is_imbalanced({0, 0}, 0.05));
    EXPECT_TRUE(CoarsestRebalanceCheck::is_imbalanced({1, 0}, 0.05));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
