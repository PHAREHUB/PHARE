
#include "core/def/phare_mpi.hpp"

#include "test_diagnostics.ipp"

static std::string const job_file = "job_1d";
static std::string out_dir        = "phare_outputs/diags_1d/";

TYPED_TEST(Simulator1dTest, fluid)
{
    fluid_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(Simulator1dTest, particles)
{
    particles_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(Simulator1dTest, electromag)
{
    electromag_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(Simulator1dTest, allFromPython)
{
    allFromPython_test(TypeParam{job_file}, out_dir);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    out_dir += std::to_string(core::mpi::size()) + "/"; // concurrent tests
    return RUN_ALL_TESTS();
}
