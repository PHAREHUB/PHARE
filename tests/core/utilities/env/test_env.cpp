

#include "core/env.hpp"
#include "tests/core/utilities/test_mpi_utils.hpp"

#include "gtest/gtest.h"
#include <cstdlib>

namespace PHARE
{

TEST(Env, env_setting_works)
{
    {
        auto& env = Env::reinit(); // init

        EXPECT_EQ(env.PHARE_LOG(), std::nullopt);
    }

    setenv("PHARE_LOG", "NONE", true);
    auto& env = Env::reinit(); // init

    EXPECT_EQ(env.PHARE_LOG(), "NONE");
    EXPECT_EQ(env.PHARE_SCOPE_TIMING(), "0"); // default
}


TEST(Env, logging_works_rank_files)
{
    setenv("PHARE_LOG", "RANK_FILES", true);
    auto& env = Env::reinit(); // init
    EXPECT_EQ(env.PHARE_LOG("PHARE_LOG_FILE"),
              ".log/" + std::to_string(core::mpi::rank()) + ".out");
}

} // namespace PHARE

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::core::mpi::Lifecycle mpi_init(argc, argv);
    return RUN_ALL_TESTS();
}
