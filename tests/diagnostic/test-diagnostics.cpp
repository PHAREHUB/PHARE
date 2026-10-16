
#include "phare_mpi.hpp"

#include "test_diagnostics.ipp"

// PHARE_TEST_DIM/PHARE_TEST_INTERP set per binary by CMake

template<typename SimulatorParam>
struct SimulatorDiagnosticsTest : public ::testing::Test
{
};

using DiagnosticsSimulators = SimulatorPermutations<PHARE_TEST_DIM, PHARE_TEST_INTERP>::type;
TYPED_TEST_SUITE(SimulatorDiagnosticsTest, DiagnosticsSimulators);

static std::string const job_file = "job_" + std::to_string(PHARE_TEST_DIM) + "d";
static std::string out_dir        = "phare_outputs/diags_" + std::to_string(PHARE_TEST_DIM) + "d/";

TYPED_TEST(SimulatorDiagnosticsTest, fluid)
{
    fluid_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(SimulatorDiagnosticsTest, particles)
{
    particles_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(SimulatorDiagnosticsTest, electromag)
{
    electromag_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(SimulatorDiagnosticsTest, allFromPython)
{
    allFromPython_test(TypeParam{job_file}, out_dir);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    out_dir += std::to_string(mpi::size()) + "/"; // concurrent tests
    return RUN_ALL_TESTS();
}
