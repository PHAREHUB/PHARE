// PHARE_LOG_SCOPE recording into phlop, whatever log level PHARE itself is built with
//  built only with -DwithPhlop=ON, see tests/core/utilities/CMakeLists.txt

#undef PHARE_LOG_LEVEL
#undef PHARE_LOG_SCOPE_FLAGS_1
#undef PHARE_LOG_SCOPE_FLAGS_2
#undef PHARE_LOG_SCOPE_FLAGS_3
#define PHARE_LOG_LEVEL 3
#define PHARE_LOG_SCOPE_FLAGS_1 time, mem
#define PHARE_LOG_SCOPE_FLAGS_2 time
#define PHARE_LOG_SCOPE_FLAGS_3 time

#include "core/logger.hpp"

#include "gtest/gtest.h"

#include <string>
#include <cstdlib>
#include <filesystem>

static_assert(PHARE_HAVE_PHLOP, "phlop headers not found");

// no mpi.h for core only targets: the serial timer, the mpi timer is for mpirun runs
#if __has_include("mpi.h")
#error "core tests must not see mpi.h, check this target's link libraries"
#endif
static_assert(_PHLOP_TIMER_EXPORTED_ == 0, "expected phlop's serial timer");

namespace
{
auto& timer()
{
    return phlop::scope_timer();
}

std::string const timer_file = "phlop_scopes.txt";

void level3()
{
    PHARE_LOG_SCOPE(3, "level3");
}
void level2()
{
    PHARE_LOG_SCOPE(2, "level2");
    level3();
}
void level1()
{
    PHARE_LOG_SCOPE(1, "level1");
    level2();
}

} // namespace


struct PhlopScopesTest : public ::testing::Test
{
    PhlopScopesTest() { std::filesystem::remove(timer_file); }
    ~PhlopScopesTest() override { phlop::shutdown_scope_timer(); }
};


TEST_F(PhlopScopesTest, inactiveTimerRecordsNothing)
{
    level1();
    EXPECT_TRUE(timer().traces.empty());
}


TEST_F(PhlopScopesTest, eachLevelIsAScopeWithItsFlags)
{
    timer().file_name(timer_file).init();
    level1();

    ASSERT_EQ(timer().traces.size(), 1u);
    auto const& l1 = *timer().traces[0];
    EXPECT_EQ(l1.self->k, "level1");
    EXPECT_TRUE(l1.has_mem); // PHARE_LOG_SCOPE_FLAGS_1

    ASSERT_EQ(l1.childs.size(), 1u);
    auto const& l2 = *l1.childs[0];
    EXPECT_EQ(l2.self->k, "level2");
    EXPECT_FALSE(l2.has_mem);

    ASSERT_EQ(l2.childs.size(), 1u);
    EXPECT_EQ(l2.childs[0]->self->k, "level3");
}


TEST_F(PhlopScopesTest, shutdownWritesTheTimerFile)
{
    timer().file_name(timer_file).init();
    level1();
    phlop::shutdown_scope_timer();

    EXPECT_FALSE(timer().active);
    EXPECT_TRUE(timer().traces.empty());
    ASSERT_TRUE(std::filesystem::exists(timer_file));
    EXPECT_GT(std::filesystem::file_size(timer_file), 0u);
}


TEST_F(PhlopScopesTest, resetWritesAndStartsAgain) // as SamraiLifeCycle::reset() does
{
    timer().file_name(timer_file).init();
    level1();
    phlop::reset_scope_timer();

    EXPECT_TRUE(std::filesystem::exists(timer_file));
    EXPECT_TRUE(timer().active);
    EXPECT_TRUE(timer().traces.empty());

    level1();
    EXPECT_EQ(timer().traces.size(), 1u);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::setenv("PHLOP_SCOPE_MEM", "1", 1); // read by init()
    return RUN_ALL_TESTS();
}
