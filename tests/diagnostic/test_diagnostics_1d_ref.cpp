
#include "core/def/phare_mpi.hpp"

#include "tests/simulator/per_test.hpp"
#include "test_diagnostics.ipp"

static std::string const job_file = "job_1d";
static std::string const out_dir  = "phare_outputs/diags_1d/";


TEST(Simulator1dTest, allFromPython0)
{
    allFromPython_test(SimulatorTestParam<SimOpts{1, 1}>{job_file + "_ref"}, out_dir + "ref");
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
