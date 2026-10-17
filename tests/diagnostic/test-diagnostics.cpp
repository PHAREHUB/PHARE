
#include "phare_mpi.hpp"

#include <cstdlib>
#include <filesystem>

#include "test_diagnostics.ipp"

// PHARE_TEST_DIM/PHARE_TEST_INTERP set per binary by CMake

template<typename SimulatorParam>
struct SimulatorDiagnosticsTest : public ::testing::Test
{
};

using DiagnosticsSimulators = SimulatorPermutations<PHARE_TEST_DIM, PHARE_TEST_INTERP>::type;
TYPED_TEST_SUITE(SimulatorDiagnosticsTest, DiagnosticsSimulators);

static std::string const job_file = "job_" + std::to_string(PHARE_TEST_DIM) + "d";
static std::string out_dir; // set in main, unique per test binary

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
    out_dir = "phare_outputs/" + std::filesystem::path{argv[0]}.filename().string() + "/"
              + std::to_string(mpi::size()) + "/";      // concurrent tests
    setenv("PHARE_DIAGS_TEST_DIR", out_dir.c_str(), 1); // read by job_*d.py
    return RUN_ALL_TESTS();
}
