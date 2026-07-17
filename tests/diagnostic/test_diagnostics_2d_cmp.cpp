
#include "core/def/phare_mpi.hpp"

#include "tests/simulator/per_test.hpp"
#include "test_diagnostics.ipp"

static std::string const job_file = "job_2d";
static std::string const out_dir  = "phare_outputs/diags_2d/";


TEST(Simulator1dTest, allFromPython1)
{
    allFromPython_test(SimulatorTestParam<SimOpts{2, 1, AoSTS}>{job_file + "_cmp"},
                       out_dir + "cmp");
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
